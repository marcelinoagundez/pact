
       subroutine get_reaction_rate_coefficients
c      Get rate coefficients of chemical reactions
       
       use pact_data
       include 'pact.common'       
       integer i,j,k,l,id,dv,im
       real*8 tkin,g,h,s,cp,dg,keq,k0,kinf,fc,m,mscale,meff,ktop,ktop2,
     # ktop3a,ktop3b,ktop0,ktopinf,kfalloff,tkinfc
       real*8, allocatable, dimension(:) :: k0_s,kinf_s,fc_s,m_s
       logical topped

       keq=0.0d0             ! to avoid warning
       ktop=0.0d0            ! ...
       ktop0=0.0d0           ! ...
       ktopinf=0.0d0         ! ...

       ktop2=2.0d-9
       ktop3a=-29.5d0
       ktop3b=2.15d0

c Loop over height 
       do j=1,nz
          tkin=tk(j)

c Loop over reactions to evaluate rate coefficients
          allocate(k0_s(nmaxmreac),kinf_s(nmaxmreac))
          allocate(fc_s(nmaxmreac),m_s(nmaxmreac))
          do i=1,nreac
             id=i
             if(rtype(i).eq.0)id=idreac_rev(i)

c         evaluate maximum allowed rate coefficient depending on molecularity of reaction
             if(rtype(id).eq.1)then
                if(nreag(i).eq.1)then
                   ktop=huge
                elseif(nreag(i).eq.2)then
                   ktop=ktop2
                elseif(nreag(i).eq.3)then
                   ktop=10.0d0**(ktop3a+ktop3b*natc(i))
                else
                   stop' E- odd molecularity in bimolecular reaction'
                endif
             elseif(rtype(id).eq.2)then
                if(nreag(i).eq.1)then
                   ktop0=ktop2
                   ktopinf=huge
                elseif(nreag(i).eq.2)then
                   ktop0=10.0d0**(ktop3a+ktop3b*natc(i))
                   ktopinf=ktop2
                else
                   stop' E- odd molecularity in p-dependent reaction'
                endif
             elseif(rtype(id).eq.3)then
                if(nreag(i).eq.2)then
                   ktop0=ktop2
                   ktopinf=huge
                else
                   stop' E- odd molecularity in special p-reaction'
                endif
             else
                stop' E- Unknown reaction type'
             endif

c         evaluate equilibrium constant if reverse reaction
             if(rtype(i).eq.0)then
                dg=0.0d0             ! DeltaG/RT
                do l=1,nreag(i)
                   call evaluate_therm(idreag(i,l),tkin,g,h,s,cp)
                   dg=dg-g
                enddo
                do l=1,nprod(i)
                   call evaluate_therm(idprod(i,l),tkin,g,h,s,cp)
                   dg=dg+g
                enddo
                dv=nprod(i)-nreag(i)
                keq=(1.0d6/kboltz/tkin)**dv*dexp(-dg)
             endif

c         evaluate rate coefficient for bimolecular reaction (forward or reverse)
             if(rtype(id).eq.1)then
                if(rtype(i).eq.1)then
                   krate(i,j)=alpha(i)*(tkin/300.0d0)**beta(i)*
     #             dexp(-gamma(i)/tkin)
                elseif(rtype(i).eq.0)then
                   krate(i,j)=keq*krate(id,j)
                endif
                if(krate(i,j).gt.ktop)then
                   if(rtype(i).eq.0)krate(id,j)=krate(id,j)*    ! if reverse, correct also forward reaction
     #             ktop/krate(i,j)
                   krate(i,j)=ktop
                endif

c         evaluate rate coefficient for pressure dependent reaction (forward or reverse)
             elseif(rtype(id).eq.2.or.rtype(id).eq.3)then
                im=idmreac(id)
                if(rtype(i).eq.2.or.rtype(i).eq.3)then
                   k0=alpha(i)*(tkin/300.)**beta(i)*dexp(-gamma(i)/tkin)
                   kinf=alphainf(im)*(tkin/300.)**betainf(im)*
     #             dexp(-gammainf(im)/tkin)
                   tkinfc=tkin
                   if(bfcent(im).ne.0.0d0)then
                      if(tkminfc(im).eq.0.0d0.or.tkmaxfc(im).eq.0.0d0)
     #                stop' E- Missing Fc temp range'
                      if(tkin.lt.tkminfc(im))tkinfc=tkminfc(im)
                      if(tkin.gt.tkmaxfc(im))tkinfc=tkmaxfc(im)
                   endif
                   fc=afcent(im)+bfcent(im)*tkinfc
                   m=density(j)
                   if(mcoll(im).ne.'')then                 ! if collider efficiency information is available
                      if(trim(mcoll(im)).ne.'N2')then      ! ...first, scale k0 to M=N2 if needed
                         mscale=0.0d0
                         do l=1,nmaxmcoll
                            if(trim(mcoll(im)).ne.trim(mcollider(l)))
     #                      cycle
                            mscale=1.0d0/mefficiency(l)
                            exit
                         enddo
                         if(mscale.ne.0.0d0)then
                            k0=k0*mscale
                         else
                            stop' E- Unknown M collider'
                         endif
                      endif
                      m=0.0d0                              ! ...second, evaluate [M] taking into account known collider efficiencies
                      do k=1,nspec
                         meff=1.0d0
                         do l=1,nmaxmcoll
                            if(trim(spec(k)).ne.trim(mcollider(l)))cycle
                            meff=mefficiency(l)
                            exit
                         enddo
                         m=m+meff*abun(k,j)*density(j)
                      enddo
                   endif
                   k0_s(im)=k0                   ! store k0, kinf, Fc, and [M] for next iteration in loop,
                   kinf_s(im)=kinf               ! when rate coefficient of reverse is to be computed 
                   fc_s(im)=fc
                   m_s(im)=m
                elseif(rtype(i).eq.0)then
                   k0=keq*k0_s(im)
                   kinf=keq*kinf_s(im)
                   fc=fc_s(im)
                   m=m_s(im)
                endif

                topped=.false.
                if(k0.gt.ktop0)then              ! if k0 or kinf exceed maximum allowed values ...
                   if(rtype(i).eq.0)k0_s(im)=k0_s(im)*ktop0/k0
                   k0=ktop0
                   topped=.true.
                endif
                if(kinf.gt.ktopinf)then
                   if(rtype(i).eq.0)kinf_s(im)=kinf_s(im)*ktopinf/kinf
                   kinf=ktopinf
                   topped=.true.
                endif
                if(rtype(i).eq.0.and.topped)then ! if reverse and k0/kin exceeded top values, correct also forward
                   krate(id,j)=kfalloff(
     #             k0_s(im),kinf_s(im),fc_s(im),m_s(im))
                   if(rtype(id).eq.3)krate(id,j)=krate(id,j)/m_s(im)
                endif
                krate(i,j)=kfalloff(k0,kinf,fc,m)
                if(rtype(id).eq.3)krate(i,j)=krate(i,j)/m        ! special pressure dependent reaction, see Eq(13) of Jasper et al (2007, J Phys Chem A, 111, 3932)

             else
                stop' E- Unknown reaction type'

             endif

c End of loop over reactions to evaluate rate constants 
          enddo
          deallocate(k0_s,kinf_s,fc_s,m_s)

c End of loop over height 
       enddo

       return

       end




c_______________________________________________________________________

       subroutine get_diffusion_terms
c      Get diffusion coefficients, and other diffusion stuff
       
       use pact_data
       include 'pact.common'       
       real*8 p_a,p_b,p_c,p_d,p_e,p_f,p_g,p_h
       parameter(p_a=1.06036d0,p_b=0.15610d0,p_c=0.19300d0,
     # p_d=0.47635d0,p_e=1.03587d0,p_f=1.52996d0,p_g=1.76474d0,
     # p_h=3.89411d0)

       integer i,j,k,l,h
       real*8 norm,s,lj_lena,lj_lenb,lj_enea,lj_eneb,omega,musqrt,
     # lj_len2,ta,atd,r,dz,rm,rp,dzm,dzp,nm,np,tkm,tkp,hspem,hspep,
     # ddifm,ddifp,kdifm,kdifp,dtkm,dtkp,fdm,fdp,hatm(nmaxz),hatmm,hatmp

       logical lj_fnda,lj_fndb

       if(.not.diffusion)return

c      ...to avoid warning
       lj_enea=0.0d0              
       lj_eneb=0.0d0              

c      ...evaluate coefficients of molecular diffusion  
       do j=1,nz
          do i=1,nspec                            !......main species a
             lj_fnda=.false.
             do l=1,nmaxlj
                if(trim(spec(i)).eq.trim(lj_spec(l)))then
                   lj_fnda=.true.
                   lj_lena=lj_length(l)*1.0d-8    ! A -> cm
                   lj_enea=lj_energy(l)
                   exit
                endif
             enddo
             if(.not.lj_fnda)lj_lena=4.0d-8       ! assumed characteristic length (4 A)
             s=0.0d0
             norm=0.0d0
             do k=1,nspec                         !......second collider b
                lj_fndb=.false.
                do l=1,nmaxlj
                   if(trim(spec(k)).eq.trim(lj_spec(l)))then
                      lj_fndb=.true.
                      lj_lenb=lj_length(l)*1.0d-8 ! A -> cm
                      lj_eneb=lj_energy(l)
                      exit
                   endif
                enddo
                if(.not.lj_fndb)lj_lenb=4.0d-8    ! assumed characteristic length (4 A)
                if(lj_fnda.and.lj_fndb.and.
     #          lj_enea.gt.0.0.and.lj_eneb.gt.0.0)then
                   ta=tk(j)/dsqrt(lj_enea*lj_eneb)
                   omega=p_a/(ta**p_b)+p_c/(p_d*ta)+
     #             p_e/(dexp(p_f*ta))+p_g/(dexp(p_h*ta))
                else
                   omega=1.0d0
                endif
                musqrt=dsqrt(mspec(i)*mspec(k)/(mspec(i)+mspec(k)))
                lj_len2=(0.5d0*(lj_lena+lj_lenb))**2.0d0
                s=s+abun(k,j)*musqrt*lj_len2*omega
                norm=norm+abun(k,j)
             enddo
             s=s*density(j)/norm
             ddif(i,j)=3.0d0/16.0d0*dsqrt(2.0d0*kboltz*tk(j)/pi)/s
          enddo
       enddo

c      ...evaluate atmospheric scale height in each layer
       do j=1,nz
          hatm(j)=kboltz*tk(j)/apmass/ggrav(j)
       enddo

c      ...compute diffusion terms to evaluate function f in continuity equation
       do j=1,nz
          r=rplanet+z(j)
          if(j.eq.1)then
             dz=z(j)-z(j+1)
          elseif(j.eq.nz)then
             dz=z(j-1)-z(j)
          else
             dz=0.5d0*(z(j-1)-z(j+1))
          endif

          h=max(2,j)                             ! evaluate at j-1/2 (if j=1, then values are not used)
          rm=rplanet+0.5d0*(z(h-1)+z(h))
          dzm=z(h-1)-z(h)
          nm=0.5d0*(density(h-1)+density(h))
          tkm=0.5d0*(tk(h-1)+tk(h))
          dtkm=(tk(h-1)-tk(h))/(z(h-1)-z(h))
          kdifm=0.5d0*(kdif(h-1)+kdif(h))
          hatmm=0.5d0*(hatm(h-1)+hatm(h))
          h=min(j,nz-1)                          ! evaluate at j+1/2 (if j=nz, then values are not used)
          rp=rplanet+0.5d0*(z(h)+z(h+1))
          dzp=z(h)-z(h+1)
          np=0.5d0*(density(h)+density(h+1))
          tkp=0.5d0*(tk(h)+tk(h+1))
          dtkp=(tk(h)-tk(h+1))/(z(h)-z(h+1))
          kdifp=0.5d0*(kdif(h)+kdif(h+1))
          hatmp=0.5d0*(hatm(h)+hatm(h+1))

          do i=1,nspec
             if(trim(spec(i)).eq.'H2'.or.
     #          trim(spec(i)).eq.'H' .or.
     #          trim(spec(i)).eq.'He')then
                atd=-0.25d0                      ! thermal diffusion factor set to -0.25 for the light species H,H2,He (Bauer 1973, Physics of Planetary Ionospheres, p 10)
             else
                atd=0.0d0                        ! thermal diffusion factor set to zero for any other species
             endif
             h=max(2,j)                          ! evaluate at j-1/2 (if j=1, then values are not used)
             hspem=0.5d0*(kboltz*tk(h-1)/mspec(i)/ggrav(h-1)+
     #                    kboltz*tk(h)  /mspec(i)/ggrav(h))
             ddifm=0.5d0*(ddif(i,h-1)+ddif(i,h))
             fdm=1.0d0/hspem-1.0d0/hatmm+atd*dtkm/tkm
             h=min(j,nz-1)                       ! evaluate at j+1/2 (if j=nz, then values are not used)
             hspep=0.5d0*(kboltz*tk(h)  /mspec(i)/ggrav(h)+
     #                    kboltz*tk(h+1)/mspec(i)/ggrav(h+1))
             ddifp=0.5d0*(ddif(i,h)+ddif(i,h+1))
             fdp=1.0d0/hspep-1.0d0/hatmp+atd*dtkp/tkp

             df_dif(i,j)=0.0d0                   ! evaluate d[f_i^j]/d(y_i^j)
             if(j.gt.1)df_dif(i,j)=df_dif(i,j)+
     #       (rm/r)**2.0d0*nm/dz*
     #       (ddifm*(0.5d0*fdm-1.0d0/dzm)-kdifm/dzm)
             if(j.lt.nz)df_dif(i,j)=df_dif(i,j)-
     #       (rp/r)**2.0d0*np/dz*
     #       (ddifp*(0.5d0*fdp+1.0d0/dzp)+kdifp/dzp)

             df_jm_dif(i,j)=0.0d0                ! evaluate d[f_i^j]/d(y_i^j-1)
             if(j.gt.1)df_jm_dif(i,j)=(rm/r)**2.0d0*nm/dz*
     #       (ddifm*(0.5d0*fdm+1.0d0/dzm)+kdifm/dzm)

             df_jp_dif(i,j)=0.0d0                ! evaluate d[f_i^j]/d(y_i^j+1)
             if(j.lt.nz)df_jp_dif(i,j)=-(rp/r)**2.0d0*np/dz*
     #       (ddifp*(0.5d0*fdp-1.0d0/dzp)-kdifp/dzp)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine get_photorates
c      Get photorates doing the ultraviolet radiative transfer
c      with absorption due to photodissociations and Rayleigh scattering by all species
c      We use a two-stream iterative algorithm assuming that 1/2 of scattered photons go up and 1/2 down
c      Algorithm from Isaksen et al 1977 (see F. Selsis thesis) slightly modified
       
       use pact_data
       include 'pact.common'       
       integer nmaxit
       parameter(nmaxit=1000)
       real*8 varmin
       parameter(varmin=1.0d-2)                  ! 1 %

       integer i,j,k,l,id,it
       real*8 csrayleigh(nmaxspec),ta,ts,n,dz,sigma,ftop(nmaxwlth),down,
     # up,fx,o_fx,var,w1,w2,s1,s2,f1,f2,evening,morning,r,path,plusangle
       real*8, allocatable, dimension(:,:) :: tau_abs,tau_sca,
     #                                        fdow,o_fdow,fupw,o_fupw
       real*8, allocatable, dimension(:)   :: polz

       if(.not.photochem)return

       if(sbrotation)then                        ! evaluate zenith angle
          plusangle=datan((rstar-rplanet)/dist)+0.5d0*pi/180.0d0
          evening=0.5d0*pi+plusangle             ! at limbs include the effect of the star apparent size
          morning=1.5d0*pi-plusangle             ! ...plus half degree due to atmospheric refraction (as in the Earth)
          if(longitude.lt.evening.or.longitude.gt.morning)then
             if(longitude.lt.evening)zenithangle=longitude
             if(longitude.gt.morning)zenithangle=2.0d0*pi-longitude
             zenithangle=min(zenithangle,0.5d0*pi-plusangle)
          else
             jphot(1:nphotoreac,1:nz)=0.0d0
             return
          endif
       endif

c      ...evaluate Rayleigh scattering cross section (part independent of wavelength)
       allocate (polz(nspec))
       do i=1,nspec                              ! evaluate polarizabilities
          polz(i)=0.0d0
          do k=1,nmaxpo
             if(trim(spec(i)).ne.trim(po_spec(k)))cycle
             polz(i)=polary(k)
             exit
          enddo
       enddo
       do i=1,nspec                              ! evaluate Rayleigh scattering cross section
          csrayleigh(i)=128.0d0/3.0d0*pi**5.0d0*(polz(i)*1.0d-24)**2.0d0
       enddo
       deallocate (polz)

c      ...evaluate flux above top layer assuming isothermal atmosphere and constant mole fractions
       dz=kboltz*tk(1)/apmass/ggrav(1)
       r=rplanet+z(1)+0.5d0*dz
c       path=dz/dcos(zenithangle)                !...plane parallel atmosphere
       path=r*(dsqrt(dcos(zenithangle)**2.0d0+   !...spherical atmosphere; see http://en.wikipedia.org/wiki/Air_mass_(astronomy)
     # (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
       do l=1,nwlth
          ta=0.0d0                               ! absorption
          do k=1,nphotospec
             id=idphotospec(k)
             n=abun(id,1)*density(1)
             ta=ta+csphotabs(k,l)*n*path
          enddo
          ts=0.0d0                               ! scattering
          do i=1,nspec
             sigma=csrayleigh(i)/(wlth(l)*1.0d-7)**4.0d0
             n=abun(i,1)*density(1)
             ts=ts+sigma*n*path
          enddo
          ftop(l)=fphot(l)*(dexp(-ta-ts)+0.5d0*(1.0d0-dexp(-ts)))
       enddo

c      ...compute opacity due to absorption and scattering for each layer
c                                      - tau(1)    is opacity between top of top layer and z(1) 
c                                      - tau(j)    is opacity between z(j-1) and z(j) 
c                                      - tau(nz+1) is opacity between z(nz) and bottom of bottom layer 
       allocate(tau_abs(nz+1,nwlth),tau_sca(nz+1,nwlth))
       do j=1,nz+1
          r=rplanet+z(min(nz,j))
          if(j.eq.1)then
             dz=0.5d0*(z(j)-z(j+1))
          elseif(j.eq.nz+1)then
             dz=0.5d0*(z(j-2)-z(j-1))
             r=rplanet+z(j-1)-0.5d0*dz
          else
             dz=z(j-1)-z(j)
          endif
c          path=dz/dcos(zenithangle)                       !...plane parallel atmosphere
          path=r*(dsqrt(dcos(zenithangle)**2.0d0+          !...spherical atmosphere
     #    (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
          do l=1,nwlth
             tau_abs(j,l)=0.0d0                  ! absorption
             do k=1,nphotospec
                id=idphotospec(k)
                if(j.eq.1)then
                   n=abun(id,j)*density(j)
                elseif(j.eq.nz+1)then
                   n=abun(id,j-1)*density(j-1)
                else
                   n=0.5d0*
     #             (abun(id,j-1)*density(j-1)+abun(id,j)*density(j))
                endif
                tau_abs(j,l)=tau_abs(j,l)+csphotabs(k,l)*n*path
             enddo
             tau_sca(j,l)=0.0d0                  ! scattering
             do i=1,nspec
                if(j.eq.1)then
                   n=abun(i,j)*density(j)
                elseif(j.eq.nz+1)then
                   n=abun(i,j-1)*density(j-1)
                else
                   n=0.5d0*
     #             (abun(i,j-1)*density(j-1)+abun(i,j)*density(j))
                endif
                sigma=csrayleigh(i)/(wlth(l)*1.0d-7)**4.0d0
                tau_sca(j,l)=tau_sca(j,l)+sigma*n*path
             enddo
          enddo
       enddo

c      ...compute flux at each layer with a two-stream iterative algorithm
       allocate(fdow(0:nz+1,nwlth),o_fdow(0:nz+1,nwlth))
       allocate(fupw(0:nz+1,nwlth),o_fupw(0:nz+1,nwlth))
       fdow(:,:)=0.0d0
       o_fdow(:,:)=0.0d0
       fupw(:,:)=0.0d0
       o_fupw(:,:)=0.0d0

c Loop over order of Isaksen algorithm
       do it=1,nmaxit

          var=0.0d0
          do l=1,nwlth
             do j=0,nz+1                         ! evaluate downward flux
                if(j.eq.0)then                   ! ...downward flux at top of top layer
                   fdow(j,l)=ftop(l)
                   cycle
                endif
                ta=tau_abs(j,l)
                ts=tau_sca(j,l)
                if(it.eq.1)then                  ! ...first iteration (zero order estimation)
                   fdow(j,l)=fdow(j-1,l)*dexp(-ta-ts)
                   cycle
                endif
                if(j.eq.1)then
                   down=ftop(l)
                else
                   down=o_fdow(j-1,l)
                endif
                up=o_fupw(j,l)
                fdow(j,l)=(0.5d0*(1.0d0-dexp(-ts))*(down+up)+
     #          fdow(j-1,l))*dexp(-ta-ts)
             enddo
             do j=nz+1,0,-1                      ! evaluate upward flux
                if(j.eq.nz+1)then                ! ...upward flux at bottom of bottom layer (planet surface)
                   if(it.eq.1)fupw(j,l)=fdow(j,l)*albedoplanet
                   if(it.gt.1)fupw(j,l)=o_fdow(j,l)*albedoplanet
                   cycle
                endif
                ta=tau_abs(j+1,l)
                ts=tau_sca(j+1,l)
                if(it.eq.1)then                  ! ...first iteration (zero order estimation)
                   fupw(j,l)=fupw(j+1,l)*dexp(-ta-ts)
                   cycle
                endif
                if(j.eq.0)then
                   down=ftop(l)
                else
                   down=o_fdow(j,l)
                endif
                up=o_fupw(j+1,l)
                fupw(j,l)=(0.5d0*(1.0d0-dexp(-ts))*(down+up)+
     #          fupw(j+1,l))*dexp(-ta-ts)
             enddo
             do j=0,nz+1                         ! evaluate convergence criterion
                if(it.eq.1)exit
                fx=fdow(j,l)+fupw(j,l)
                o_fx=o_fdow(j,l)+o_fupw(j,l)
                if(o_fx.eq.0)cycle
                if(dabs(fx-o_fx)/o_fx.gt.var)var=dabs(fx-o_fx)/o_fx
             enddo
             do j=0,nz+1                         ! store old flux values
                o_fdow(j,l)=fdow(j,l)
                o_fupw(j,l)=fupw(j,l)
             enddo
          enddo
          if(it.gt.1.and.var.lt.varmin)exit

c End of loop order of Isaksen algorithm
       enddo

       if(it.eq.nmaxit)write(*,1000)var

c      ...compute rate of photo reactions at each height
       do j=1,nz
          do k=1,nphotoreac
             jphot(k,j)=0.0d0
             do l=1,nwlth-1
                w1=wlth(l)
                w2=wlth(l+1)
                s1=csphot(k,l)
                s2=csphot(k,l+1)
                f1=fdow(j,l)+fupw(j,l)
                f2=fdow(j,l+1)+fupw(j,l+1)
                jphot(k,j)=jphot(k,j)+0.5d0*(s1*f1+s2*f2)*(w2-w1)
             enddo
          enddo
       enddo

c      ...locate tau=1 height as a function of wavelength
c         opacity is defined with respect to the non-attenuated incoming stellar flux
c         and with respect to the downward flux at each layer
c       do l=1,nwlth
c          do j=1,nz
c             tau=dlog(fphot(l)/(fdow(j,l)))
c             if(tau.gt.1.0d0)exit
c          enddo
c          if(j.eq.1)then
c             p_tauuv1(l)=pressure(1)
c             cycle
c          endif
c          if(tau.le.1.0d0)then
c             p_tauuv1(l)=pressure(nz)
c             cycle
c          endif
c          tu=dlog(fphot(l)/(fdow(j-1,l)))
c          tl=dlog(fphot(l)/(fdow(j,l)))
c          pu=pressure(j-1)
c          pl=pressure(j)
c          slope=(pu-pl)/(tu-tl)
c          p_tauuv1(l)=pl+slope*(1.0d0-tl)
c       enddo

       deallocate(tau_abs,tau_sca,fdow,o_fdow,fupw,o_fupw)

       return

1000   format(1x,
     # ' W- Rayleigh scattering two-stream algorithm not converged',/,
     # '    Reached accuracy = ',1pg10.3)

       end




c_______________________________________________________________________

       subroutine evaluate_therm(id,tkin,g,h,s,cp)
c      Evaluate H/RT, S/R, and Cp/R for requested species and temperature

       include 'pact.common'       
       integer id,j,k
       real*8 tkin,a(9),g,h,s,cp

       k=1
       do while(tkin.ge.temp_therm(id,k+1).and.k.lt.ntemp_therm(id))
          k=k+1
       enddo
       a(1:natherm(id))=atherm(id,1:natherm(id),k)

c      ...7-term polynomials fitted to G^o/RT = H^o/RT - S^o/R (standard Gibbs energy)
       if(type_therm(id).eq.1)then
          h=0.0d0
          s=0.0d0
          cp=0.0d0
          g=0.0d0
          do j=1,7
             g=g+a(j)*tkin**(j-3)
          enddo

c      ...7-term NASA polynomials (see NASA RM-4513, p 9)
       elseif(type_therm(id).eq.7)then
          h=a(1)+a(2)*tkin/2.0d0+a(3)*tkin**2.0d0/3.0d0+
     #    a(4)*tkin**3.0d0/4.0d0+a(5)*tkin**4.0d0/5.0d0+a(6)/tkin
          s=a(1)*dlog(tkin)+a(2)*tkin+a(3)*tkin**2.0d0/2.0d0+
     #    a(4)*tkin**3.0d0/3.0d0+a(5)*tkin**4.0d0/4.0d0+a(7)
          cp=a(1)+a(2)*tkin+a(3)*tkin**2.0d0+a(4)*tkin**3.0d0+
     #    a(5)*tkin**4.0d0
          g=h-s

c      ...9-term NASA polynomials (see NASA TP-2001-210959/REV1, p 2)
       elseif(type_therm(id).eq.9)then
          h=-a(1)*tkin**(-2.0d0)+a(2)*dlog(tkin)/tkin+a(3)+
     #    a(4)*tkin/2.0d0+a(5)*tkin**2.0d0/3.0d0+
     #    a(6)*tkin**3.0d0/4.0d0+a(7)*tkin**4.0d0/5.0d0+a(8)/tkin
          s=-a(1)*tkin**(-2.0d0)/2.0d0-a(2)/tkin+a(3)*dlog(tkin)+
     #    a(4)*tkin+a(5)*tkin**2.0d0/2.0d0+a(6)*tkin**3.0d0/3.0d0+
     #    a(7)*tkin**4.0d0/4.0d0+a(9)
          cp=a(1)*tkin**(-2.0d0)+a(2)/tkin+a(3)+a(4)*tkin+
     #    a(5)*tkin**2.0d0+a(6)*tkin**3.0d0+a(7)*tkin**4.0d0
          g=h-s
       endif

       return
       end




c_______________________________________________________________________

       real*8 function kfalloff(k0,kinf,fc,m)
c      computes the rate coefficient for a pressure dependent reaction:
c      input:   k0:   low pressure rate coefficient
c               kinf: high pressure rate coefficient
c               Fc:   Troe fall-off factor
c               m:    [M] taking into account enhanced efficiencies if available
c      returns: kfalloff: rate coefficient

       implicit none
       real*8 k0,kinf,fc,m,c,n,d,pr,f

       fc=max(1.0d-300,fc)
       c=-0.4d0-0.67d0*dlog10(fc)
       n=0.75d0-1.27d0*dlog10(fc)
       d=0.14d0
       pr=k0*m/kinf
       f=dlog10(fc)/(1.0d0+((dlog10(pr)+c)/(n-d*(dlog10(pr)+c)))**2.0d0)
       f=10.0d0**f
       kfalloff=kinf*(pr/(1.0d0+pr))*f

       end
