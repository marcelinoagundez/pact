
       subroutine solve_solid_body_rotation
c      ...To reduce CPU time the evaluation of tk, density, rate constants, ...
c         is done just after each time step (not within DLSODES),
c         therefore, the physical conditions (tk, density) do not change within DLSODES

       include 'pact.common'       
       integer lrw_dlsodes,liw_dlsodes
       parameter(lrw_dlsodes=20+9*nmaxzs+
     #           2*nmaxnzjac+2*nmaxzs+(nmaxnzjac+9*nmaxzs)/2)
       parameter(liw_dlsodes=31+nmaxzs+nmaxnzjac)

       integer i,j,ii,nzs,neg,nnz,ihou,imin,isec,itask,istate,
     # iwork_dlsodes(liw_dlsodes),ia(nmaxzs+1),ja(nmaxnzjac)
       real*8 rwork_dlsodes(lrw_dlsodes),t,o_t,tout,o_tout,tmax,
     # y(nmaxzs),rtol,atol(nmaxzs),tstep,o_abun(nmaxspec,nmaxz)
       external f_dlsodes_sbr,jac_dlsodes_sbr

c      ...write info on screen
       if(diffusion)                    write(*,1000)
       if(.not.diffusion)               write(*,1010)
       if(photochem)                    write(*,1020)
       if(.not.photochem)               write(*,1030)
                                        write(*,1040)

c      ...get input parameters related to the solid body rotation
       call read_input_solidbodyrotation

c      ...Initialize variables related to time integration
       nzs=nz*nspec                              ! number of layers x number of species
       tstep=period/dble(nlong-1)                ! time step [s]
       tmax=nperiod*period                       ! maximum time
       t=0.0d0                                   ! initial time
       tout=tstep                                ! first output time
c       rtol=1.0d-8                               ! relative tolerance for DLSODES
       rtol=1.0d-3                               ! relative tolerance for DLSODES
       itask=1                                   ! itask for DLSODES
       call conservation(0)
       write(*,*)

c      ...get initial physical conditions and write out header and initial results
       time=t
       call eval_phys_solidbodyrotation
       call write_results_solidbodyrotation
       ostat_sbr=1

c Big loop of time (use DLSODES to integrate the system of ODES)
       do while(t.lt.tmax)

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
          do i=1,nzs
c             atol(i)=max(1.0d-20,1.0d-15*y(i))
             atol(i)=1.0d-99
          enddo

c      ...get rate constants, photorates, and diffusion terms (physical conditions already evaluated)
          time=t
c          call eval_phys_solidbodyrotation
          call get_reaction_rate_coefficients
          call get_diffusion_terms
          call get_photorates

c      ...reset current initial time to zero
          tin_sbr=t
          t=0.0d0
          tout=tout-tin_sbr

c      ...small loop of time (normally just one step if DLSODES integrates up to tout)
          do while (t.lt.tout)
             istate=1                            ! call DLSODES
             iwork_dlsodes(:)=0
             rwork_dlsodes(:)=0.0d0
             iwork_dlsodes(6)=10000
             iwork_dlsodes(7)=2
             call compute_iaja_sbr(nmaxzs+1,nmaxnzjac,ia,ja,nnz)
             iwork_dlsodes(30+1:30+nzs+1)=ia(1:nzs+1)
             iwork_dlsodes(31+nzs+1:31+nzs+nnz)=ja(1:nnz)
             call DLSODES(f_dlsodes_sbr,nzs,y,t,tout,2,rtol,atol,
     #       itask,istate,1,rwork_dlsodes,lrw_dlsodes,iwork_dlsodes,
     #       liw_dlsodes,jac_dlsodes_sbr,21)
             neg=0
             do j=1,nz                           ! update abundances
                do i=1,nspec
                   ii=(j-1)*nspec+i
                   abun(i,j)=y(ii)
                   if(abun(i,j).lt.-atol(i))then
                      write(*,1070)j,z(j)/1.0d5,spec(i),abun(i,j)
                      neg=neg+1
                      abun(i,j)=abunmin
                   endif
                enddo
             enddo
c             if(neg.gt.0)stop 'E- Negative abundance found'
c      ...end of small loop of time
          enddo

c      ...return to real times
          t=tin_sbr+t
          tout=tin_sbr+tout

c      ...verify conservation of elements
          call conservation(2)

c      ...update average particle mass in planet's atmosphere: apmass
          call compute_average_particle_mass

c      ...get current physical conditions and write out results
          time=t
          call eval_phys_solidbodyrotation
          call write_results_solidbodyrotation

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          write(*,1050)t,t/period,ihou,imin,isec

c      ...get next tout
          tout=tout+tstep
          if(tout.gt.tmax)tout=tmax

c End of big loop of time
       enddo

c      ...finish
       ostat_sbr=2
       call write_results_solidbodyrotation
       write(*,1060)t

       return

1000   format(1x,' I- Diffusion included')
1010   format(1x,' I- Diffusion not included')
1020   format(1x,' I- Photochemistry included')
1030   format(1x,' I- Photochemistry not included')
1040   format(1x,' I- Solid body rotating atmosphere')

1050   format(2x,'time[s]=',1pg10.3,3x,
     #           't/period=',1pg10.3,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s')
1060   format(1x,' I- Integration done succesfully up to ',1pg8.1,'s')
1070   format(1x,' E- In layer ',i3,' z=',1pg10.2,' km ',/,
     #        1x,'    molar fraction of ',a12,' is ',1pg10.2)

       end




c_______________________________________________________________________

       subroutine f_dlsodes_sbr(neq,t,y,ydot)
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

       subroutine jac_dlsodes_sbr(neq,t,y,jth,ian,jan,pdj)
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

       subroutine compute_iaja_sbr(nia,nja,ia,ja,nnz)
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

       subroutine eval_phys_solidbodyrotation
c      ...evaluate for current time - longitude
c                                   - tk for each layer
c                                   - density for each layer
c                                   - species scale height for each layer
c                                   - height for each layer under hydrostatic equilibrium
c      ...we assume that pressure remains constant within each layer

       include 'pact.common'       
       integer i,j,l
       real*8 slope,g,h,m,norm

c     ...compute current longitude 
       longitude=(time-idint(time/period)*period)/period*2.0d0*pi
       if(longitude.lt.0.0d0.or.longitude.gt.2.0d0*pi)
     # stop' E- longitude outside (0 - 2pi) range'

c     ...locate current longitude in longitude grid 
       l=1
       do while(long(l).lt.longitude.and.l.lt.nlong)
          l=l+1
       enddo
       l=l-1
       if(l.lt.1)l=1

c     ...compute current vertical structure of tk, density, species scale height 
       do j=1,nz
          slope=(tk_grid(l+1,j)-tk_grid(l,j))/(long(l+1)-long(l))
          tk(j)=tk_grid(l,j)+slope*(longitude-long(l))
          density(j)=pressure(j)/kboltz/tk(j)
       enddo

c     ...compute z vertical structure to fulfill hydrostatic equilibrium 
       do j=nz-1,1,-1                  ! bottom height kept constant at all longitudes
          g=ggravity*mplanet/(rplanet+0.5d0*(z(j+1)+z(j)))**2.0d0
          m=0.0d0
          norm=0.0d0
          do i=1,nspec
             m=m+mspec(i)*abun(i,j+1)
             norm=norm+abun(i,j+1)
          enddo
          m=m/norm
          h=kboltz*0.5d0*(tk(j+1)+tk(j))/(m*g)
          z(j)=z(j+1)-h*dlog(pressure(j)/pressure(j+1))
       enddo

       return
       end




c_______________________________________________________________________

       subroutine read_input_solidbodyrotation
c      ...read input parameters related to the solid body rotation

       include 'pact.common'       
       integer i,j,l,jj,ll,nxl,nxp
       real*8 rwk,slope,tk1,tk2
       real*8, allocatable, dimension(:)   :: xl,xp
       real*8, allocatable, dimension(:,:) :: xtk

c      ...generate longitude grid
       rwk=2.0d0*pi/dble(nlong-1)
       do l=1,nlong
          long(l)=(l-1)*rwk
          if(l.eq.1)    long(l)=0.0d0
          if(l.eq.nlong)long(l)=2.0d0*pi
       enddo

c      ...read longitude-pressure-temperature structure
       open(unit=8,file=lptfile,status='old')
       read(8,*)windspeed                        ! equatorial wind speed
       windspeed=windspeed*1.0d5                 ! km/s -> cm/s
       period=2.0d0*pi*rplanet/windspeed         ! rotation period [s]
       read(8,*)nxl
       read(8,*)nxp
       allocate(xl(nxl),xp(nxp),xtk(nxl,nxp))
       read(8,*)
       read(8,*)(xp(j),j=1,nxp)
       do j=1,nxp
          xp(j)=xp(j)*1.0d6                      ! bar -> dyn cm-2
          if(j.gt.1.and.xp(j).ge.xp(j-1))
     #    stop' E- Bad pressure ordering in lpt file'
       enddo
       read(8,*)
       do i=1,nxl
          read(8,*)xl(i),(xtk(i,j),j=1,nxp)
          xl(i)=xl(i)*pi                         ! pi -> rad
          if(i.gt.1.and.xl(i).le.xl(i-1))
     #    stop' E- Bad longitude ordering in lpt file'
       enddo
       close(8)
       if(xl(1).ge.1.0d-2)stop' E- First lpt longitude not 0'
       if(dabs(xl(nxl)/pi-2.0d0)/2.0d0.ge.1.0d-2)
     # stop' E- Last lpt longitude not 2pi'

c      ...interpolate to current z grid
       do j=1,nz
          jj=1
          do while(xp(jj).gt.pressure(j).and.jj.lt.nxp)
             jj=jj+1
          enddo
          jj=jj-1
          if(jj.lt.1)jj=1
          do l=1,nlong
             ll=1
             do while(xl(ll).lt.long(l).and.ll.lt.nxl)
                ll=ll+1
             enddo
             ll=ll-1
             if(ll.lt.1)ll=1
             slope=(xtk(ll+1,jj)-xtk(ll,jj))/(xl(ll+1)-xl(ll))
             tk1=xtk(ll,jj)+slope*(long(l)-xl(ll))
             slope=(xtk(ll+1,jj+1)-xtk(ll,jj+1))/(xl(ll+1)-xl(ll))
             tk2=xtk(ll,jj+1)+slope*(long(l)-xl(ll))
             slope=(tk2-tk1)/(dlog10(xp(jj+1))-dlog10(xp(jj)))
             tk_grid(l,j)=tk1+
     #       slope*(dlog10(pressure(j))-dlog10(xp(jj)))
          enddo
       enddo
       deallocate(xl,xp,xtk)

       return
       end




c_______________________________________________________________________

       subroutine write_results_solidbodyrotation
c      ...write abundances at current time (longitude)

       include 'pact.common'
       integer i,j
       real*8 sum
       character*(nmaxcharfile) outfile
       character*3 difftxt,phottxt

c If first call: open files and write header
       if(ostat_sbr.eq.0)then

c      ...prepare character strings to write in header
          difftxt='No '
          if(diffusion)difftxt='Yes'
          phottxt='No '
          if(photochem)phottxt='Yes'

c      ...write out abundances
          outfile(1:nmaxcharfile)=''
          outfile(1:ncharmodel)=namemodel(1:ncharmodel)
          outfile(ncharmodel+1:ncharmodel+8)='_out.dat'
          open(unit=11,file=outfile,status='unknown')
          write(11,1000)specfile,lptfile,reacfile,starfile,
     #    dist/au,rstar/rsun,rplanet/rearth,mplanet/mearth,
     #    windspeed/1.0d5,period,nperiod,difftxt,phottxt
          write(11,1010)(i,i=1,nspec+4)
          write(11,1020)(spec(i),i=1,nspec)
          do j=nz,1,-1
             write(11,1030)2.0d0*time/period,z(j)/1.0d5,
     #       pressure(j)/1.0d6,tk(j),
     #       (dlog10(max(abunmin,abun(i,j))),i=1,nspec)
          enddo

c If second or next call: write abundances
       elseif(ostat_sbr.eq.1)then

c      ...write out abundances
          do j=nz,1,-1
             sum=1.0d0                           ! do not normalize abundances
c             sum=0.0d0                           ! normalize abundances
c             do i=1,nspec
c                sum=sum+abun(i,j)
c             enddo
             write(11,1030)2.0d0*time/period,z(j)/1.0d5,
     #       pressure(j)/1.0d6,tk(j),
     #       (dlog10(max(abunmin,abun(i,j)/sum)),i=1,nspec)
          enddo

c If last call: close files
       elseif(ostat_sbr.eq.2)then
          close(11)

       else
          stop' E- Wrong status index to write results'

c End of logical structure
       endif

       return

1000   format(
     # '!',/,
     # '! 1D solid body rotating atmosphere chemical model',
     #                                ' with parameters:',/,
     # '! species file : ',a,/,
     # '! (l,p,T) temperature structure file : ',a,/,
     # '! reaction rate coefficient file : ',a,/,
     # '! stellar spectrum file : ',a,/,
     # '! Star-planet distance                 =',1pg12.4,' AU',/,
     # '! Star radius                          =',1pg12.4,' Rsun',/,
     # '! Planet radius                        =',1pg12.4,' REarth',/,
     # '! Planet mass                          =',1pg12.4,' MEarth',/,
     # '! Equatorial wind speed at 1 bar       =',1pg12.4,' km s-1',/,
     # '! Rotation period                      =',1pg12.4,' s',/,
     # '! Number of rotation periods           =',i3,/,
     # '! Diffusion included ?                 : ',a,/,
     # '! Photochemistry included ?            : ',a,/,
     # '! .......................................................',/,
     # '! Mixing ratios : ')
1010   format('!',4(7x,i2,7x),2x,999(3x,i3,6x))
1020   format('!','  longitude[pi] ','   height[km]   ',
     # ' pressure[bar]  ',' temperature[K] ',7x,999(a12))
1030   format(3x,0pf12.3,4x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,
     # 999(2x,0pf10.6))

       end
