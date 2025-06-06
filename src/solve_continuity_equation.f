
       subroutine solve_continuity_equation
       
       include 'pact.common'       
       integer lrw_dlsodes,liw_dlsodes
       parameter(lrw_dlsodes=20+9*nmaxzs+
     #           2*nmaxnzjac+2*nmaxzs+(nmaxnzjac+9*nmaxzs)/2)
       parameter(liw_dlsodes=31+nmaxzs+nmaxnzjac)

       integer i,j,ii,nzs,neg,nnz,ihou,imin,isec,itask,istate,
     # iwork_dlsodes(liw_dlsodes),ia(nmaxzs+1),ja(nmaxnzjac),nsoft2,
     # nminsoft2
       real*8 rwork_dlsodes(lrw_dlsodes),t,o_t,tout,o_tout,tmax,
     # y(nmaxzs),rtol,atol(nmaxzs),tstep,tmin,tmid,var,steady1,steady2,
     # minabun,tin,o_abun(nmaxspec,nmaxz),o_varabun,minvar
       external f_dlsodes,jac_dlsodes

       if(self_temp)                    write(*,1000)
       if(.not.self_temp)               write(*,1010)
       if(diffusion)                    write(*,1020)
       if(.not.diffusion)               write(*,1030)
       if(photochem)                    write(*,1040)
       if(.not.photochem)               write(*,1050)

c      ...Initialize variables
       nzs=nz*nspec                              ! number of layers x number of species
       t=0.0d0                                   ! initial time
       tout=1.0d-10                              ! first output time
       tmin=day                                  ! minimum time to integrate
       tmid=100.0d0*day                          ! middle time to evaluate stop integration
       tmax=13.75d9*year                         ! maximum time: age of the Universe
       tstep=10.0d0**0.2d0                       ! time step factor after tmin
       minabun=1.0d-20                           ! lower abundance cutoff to evaluate steady state
       steady1=1.0d-2                            ! hard relative tolerance to consider steady state
       steady2=1.0d-1                            ! soft relative tolerance to consider steady state
       nsoft2=0                                  ! number of consecutive steps with var < soft relative tolerance
       nminsoft2=3                               ! minimum number of consecutive steps with var < soft relative tolerance to consider steady state
       itime=0                                   ! counter of integration times to write out results
       varabun=huge                              ! relative abundance variation
       minvar=huge                               ! relative abundance variation
c       rtol=1.0d-8                               ! relative tolerance for DLSODES
       rtol=1.0d-3                               ! relative tolerance for DLSODES
       itask=1                                   ! itask for DLSODES
       write_supp=.false.                        ! write supplementary output files
       initial=.false.                           ! initialization at time = 0 finished
       call conservation(0)
       write(*,*)

c Big loop of time (use DLSODES to integrate the system of ODE)
       do while(t.lt.tmax.and..not.steady)

c      ...store times and abundances as "old" values, assign abundances to array y
          o_t=t
          o_tout=tout
          do j=1,nz
             do i=1,nspec
                ii=(j-1)*nspec+i
                y(ii)=abun(i,j)
                o_abun(i,j)=abun(i,j)
             enddo
          enddo
          o_varabun=varabun
          do i=1,nzs
c             atol(i)=max(1.0d-20,1.0d-15*y(i))
             atol(i)=1.0d-99
          enddo

c      ...get rate constants, photorates, and diffusion terms
          call get_reaction_rate_coefficients
          call get_diffusion_terms
          call get_photorates

c      ...update temperature
          if(self_temp.and.t.ne.0.0d0)then
             call solve_opacity
             call solve_temperature
             call adjust_vertical_structure
             call get_analytical_eddy
          endif

c      ...reset current initial time to zero
          tin=t
          t=0.0d0
          tout=tout-tin

c      ...small loop of time (normally just one step if DLSODES integrates up to tout)
          do while (t.lt.tout)
             istate=1                            ! call DLSODES
             iwork_dlsodes(:)=0
             rwork_dlsodes(:)=0.0d0
             iwork_dlsodes(6)=10000
             call compute_iaja(nmaxzs+1,nmaxnzjac,ia,ja,nnz)
             iwork_dlsodes(30+1:30+nzs+1)=ia(1:nzs+1)
             iwork_dlsodes(31+nzs+1:31+nzs+nnz)=ja(1:nnz)
             call DLSODES(f_dlsodes,nzs,y,t,tout,2,rtol,atol,itask,
     #       istate,1,rwork_dlsodes,lrw_dlsodes,iwork_dlsodes,
     #       liw_dlsodes,jac_dlsodes,21)
             neg=0
             do j=1,nz                           ! update abundances
                do i=1,nspec
                   ii=(j-1)*nspec+i
                   abun(i,j)=y(ii)
                   if(abun(i,j).lt.-atol(i))then
                      write(*,1060)j,z(j)/1.0d5,spec(i),abun(i,j)
                      neg=neg+1
                      abun(i,j)=abunmin
                   endif
                enddo
             enddo
c             if(neg.gt.0)stop 'E- Negative abundance found'
c      ...end of small loop of time
          enddo

c      ...return to real times
          t=tin+t
          tout=tin+tout
          time=t

c      ...verify conservation of elements
          call conservation(2)

c      ...update average particle mass in planet's atmosphere: apmass
          call compute_average_particle_mass

c      ...verify whether steady state has been reached
          if(t.ge.tmin)then
             varabun=tiny
             write_supp=.false.
             do j=1,nz
                do i=1,nspec
                   if(abun(i,j).lt.minabun)cycle
                   var=dabs((abun(i,j)-o_abun(i,j))/o_abun(i,j))
                   if(var.gt.varabun)varabun=var
                enddo
             enddo
             if(varabun.lt.steady1)then
                steady=.true.
                write_supp=.true.
             else
                nsoft2=nsoft2+1
                if(varabun.ge.steady2)nsoft2=0
                if(nsoft2.ge.nminsoft2)then
                   steady=.true.
                   write_supp=.true.
                else
                   if(t.ge.tmid)then
                      if(varabun.lt.minvar)then
                         write_supp=.true.
                      else
                         write_supp=.false.
                      endif
                      if(varabun.lt.minvar)minvar=varabun
                   endif
                endif
             endif

             itime=itime+1
             call write_results
          else
             varabun=huge
          endif

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          write(*,1070)t,ihou,imin,isec,varabun

c      ...get next tout
          tout=tout*tstep                        ! if t > tmin
          if(t.lt.tmin)tout=o_tout*ftstep        ! if t < tmin
          if(t.lt.tmin.and.tout.gt.tmin)tout=tmin
          if(tout.gt.tmax)tout=tmax

c End of big loop of time
       enddo

       if(steady)then
          write(*,1080)t
       else
          write(*,1090)t
       endif

       return

1000   format(1x,' I- Self-consistent temperature included')
1010   format(1x,' I- Self-consistent temperature not included')
1020   format(1x,' I- Diffusion included')
1030   format(1x,' I- Diffusion not included')
1040   format(1x,' I- Photochemistry included')
1050   format(1x,' I- Photochemistry not included')
1060   format(1x,' E- In layer ',i3,' z=',1pg10.2,' km ',/,
     #        1x,'    molar fraction of ',a12,' is ',1pg10.2)
1070   format(2x,'time[s]=',1pg10.3,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s',3x,
     #           'err=',1pg10.3) 
1080   format(1x,' I- Steady state reached at t=',1pg8.1,'s')
1090   format(1x,' I- Integration done up to t=',1pg8.1,'s')

       end




c_______________________________________________________________________

       subroutine f_dlsodes(neq,t,y,ydot)
c      user subroutine for DLSODES to be used in evolution calculation
c      it evaluates the derivatives of the abundances
       
       include 'pact.common'     
       integer i,j,ii,nzs,neq
       double precision t,y(*),ydot(*)

       nzs=nz*nspec
       i=int(t)                        ! to avoid warning
       if(neq.ne.nzs) stop ' E- Dimension error in f_dlsodes'
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             abun(i,j)=y(ii)
          enddo
       enddo

       call evaluate_derivatives

       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             ydot(ii)=f_conteq(i,j)/density(j)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine jac_dlsodes(neq,t,y,jth,ian,jan,pdj)
c      user subroutine for DLSODES to be used in evolution calculation
c      it evaluates the jacobian of the system of ODEs
       
       include 'pact.common'       
       integer i,j,l,g,nzs,gwk,id,ii,jth,neq,iwk,jwk
       double precision t,y(*),ian(*),jan(*),pdj(*)
       real*8 rwk,iswk

       nzs=nz*nspec
       i=int(t)                        ! to avoid warning
       i=int(ian(1))                   ! to avoid warning
       i=int(jan(1))                   ! to avoid warning
       if(neq.ne.nzs) stop ' E- Dimension error in jac_dlsodes'
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             abun(i,j)=y(ii)
          enddo
       enddo

       jwk=(jth-1)/nspec+1       
       iwk=jth-(jwk-1)*nspec
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             if(jwk.eq.j)then

c      ...evaluate jacobian of formation/destruction: d[dn_i^j/dt]/dy_iwk^j = d[dn_i^j/dt]/dn_iwk^j * n^j 
                df_conteq(i,j,iwk)=0.0d0
                do l=1,nterm(i)
                   id=idterm(i,l)
                   iswk=0.0d0
                   gwk=0
                   do g=1,nreag(id)
                      if(idreag(id,g).eq.iwk)then
                         iswk=iswk+1.0d0
                         gwk=g
                      endif
                   enddo
                   if(iswk.gt.0.0d0)then
                      rwk=1.0d0
                      do g=1,nreag(id)
                         if(g.ne.gwk)
     #                   rwk=rwk*abun(idreag(id,g),j)*
     #                   density(j)
                      enddo
                      df_conteq(i,j,iwk)=df_conteq(i,j,iwk)+density(j)*
     #                sterm(i,l)*fterm(i,l)*krate(id,j)*rwk*iswk
                   endif
                enddo
c     ...add diffusion contribution to jacobian: d[dn_i^j/dt]/dy_iwk^j 
                if(diffusion)then
                   if(i.eq.iwk)then
                      df_conteq(i,j,iwk)=df_conteq(i,j,iwk)+df_dif(i,j)
                   endif
                endif

c     ...add photochemical contribution to jacobian: d[dn_i^j/dt]/dy_iwk^j 
                if(photochem)then
                   do l=1,nphototerm(i)
                      id=idphototerm(i,l)
                      if(idphotoreag(id).eq.iwk)then
                         df_conteq(i,j,iwk)=df_conteq(i,j,iwk)+
     #                   density(j)*sphototerm(i,l)*fphototerm(i,l)*
     #                   jphot(id,j)
                      endif
                   enddo
                endif
c     ...assign to pdj 
                pdj(ii)=df_conteq(i,j,iwk)/density(j)
             elseif(jwk.eq.(j-1).and.iwk.eq.i)then
                pdj(ii)=df_jm_dif(i,j)/density(j)
             elseif(jwk.eq.(j+1).and.iwk.eq.i)then
                pdj(ii)=df_jp_dif(i,j)/density(j)
             else
                pdj(ii)=0.0d0
             endif
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine compute_iaja(nia,nja,ia,ja,nnz)
c      ...evaluate ia,ja for DLSODES

       include 'pact.common'       
       integer nzs,nia,nja,ia(nia),ja(nja),nnz,i,j,icol,jcol,irow,jrow

       nzs=nz*nspec
       call evaluate_jacobian
       nnz=1
       do jcol=1,nz
          do icol=1,nspec
             j=(jcol-1)*nspec+icol
             ia(j)=nnz
             do jrow=1,nz
                do irow=1,nspec
                   i=(jrow-1)*nspec+irow
                   if(jcol.eq.jrow)then
                      if(df_conteq(irow,jrow,icol).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   endif
                   if(.not.diffusion)cycle
                   if(jcol.eq.(jrow-1).and.icol.eq.irow)then
                      if(df_jm_dif(irow,jrow).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   elseif(jcol.eq.(jrow+1).and.icol.eq.irow)then
                      if(df_jp_dif(irow,jrow).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       ia(nzs+1)=nnz

       return
       end




c_______________________________________________________________________

       subroutine conservation(ido)
c      ...compute column densities of atoms as
c         the total number of atoms in the atmosphere divided by planet surface area
c         that is, the spherically averaged element column density

       include 'pact.common'
       real*8 constol                            ! tolerance for conservation
       parameter(constol=0.05d0)                 ! 5 %

       integer i,j,k,ido
       real*8 dz,rm,rp,err,cdelem(nmaxelem),cdtot

       cdelem(:)=0.0d0
       cdtot=0.0d0
       do j=1,nz
          if(j.eq.1)then
             dz=z(j)-z(j+1)
          elseif(j.eq.nz)then
             dz=z(j-1)-z(j)
          else
             dz=dabs(0.5d0*(z(j-1)+z(j+1)))
          endif
          rm=rplanet+z(j)-0.5d0*dz
          rp=rplanet+z(j)+0.5d0*dz
          do k=1,nelem
             do i=1,nspec
                if(nat(i,k).eq.0)cycle
                cdelem(k)=cdelem(k)+nat(i,k)*abun(i,j)*density(j)*
     #          ((rp/rplanet)**3.0d0-(rm/rplanet)**3.0d0)
             enddo
          enddo
          cdelem(k)=cdelem(k)*rplanet/3.0d0
          cdtot=cdtot+density(j)*
     #    ((rp/rplanet)**3.0d0-(rm/rplanet)**3.0d0)*rplanet/3.0d0
       enddo

       if(ido.eq.0)then                          ! compute initial column densities
          do k=1,nelem
             cdelem0(k)=cdelem(k)
          enddo
          cdtot0=cdtot
       elseif(ido.eq.1)then                      ! simply compute current column densities and return
          return
       elseif(ido.eq.2)then                      ! compare with initial column densities
          do k=1,nelem
             err=dabs(cdelem(k)/cdtot-cdelem0(k)/cdtot0)/
     #       (cdelem0(k)/cdtot0)
             if(err.gt.constol)then
                write(*,1000)trim(elem(k)),err
             endif
          enddo
       else
          stop ' E- Wrong input integer in conservation'
       endif

       return

1000   format(2x,'W- Conservation ',a,' : ',1pg8.2)

       end




c_______________________________________________________________________

       subroutine conservation_layer(ido)
c      ...compute sum of mole fraction of each atom in any form, for each layer

       include 'pact.common'       
       integer i,j,k,ido
       real*8 err,elab(nmaxelem,nmaxz)

       do k=1,nelem
          do j=1,nz
             elab(k,j)=0.0d0
             do i=1,nspec
                if(nat(i,k).eq.0)cycle
                elab(k,j)=elab(k,j)+nat(i,k)*abun(i,j)
             enddo
          enddo
       enddo

       if(ido.eq.0)then                          ! compute initial values
          do k=1,nelem
             do j=1,nz
                elab0(k,j)=elab(k,j)
             enddo
          enddo
       elseif(ido.eq.1)then                      ! simply compute current column densities and return
          return
       elseif(ido.eq.2)then                      ! compare with initial column densities
          i=0                                    ! ...and if different stop
          do j=1,nz
             do k=1,nelem
                err=dabs((elab(k,j)-elab0(k,j))/elab(k,j))
                if(err.gt.5.0d-3)then
                   write(*,*)' E- Element: ',elem(k)
                   write(*,*)'    layer:   ',j,z(j)/1.0d5,' km '
                   write(*,*)'    current-to-initial abundances sum =',
     #             elab(k,j)/elab0(k,j)
                   i=i+1
                endif
             enddo
          enddo
          if(i.gt.0)stop
       else
          stop ' E- Wrong input integer in conservation_layer'
       endif

       return
       end
