
       subroutine write_results

       include 'pact.common'
       integer i,j,k,l,n,ilen
       real*8 sum,lammin,lammax
       character*(nmaxcharfile) outfile
       character*(nmaxlenspec*nmaxspec_ktab) txt_ktab
       character*(nmaxcharfile) txt_tsol
       character*3 txt_tsel,txt_diff,txt_phot,txt_deta,txt_stea

c      ...build character strings to write in header
       txt_ktab(1:nmaxlenspec*nmaxspec_ktab)=''
       ilen=1
       do i=1,nspec_ktab
          j=len_trim(spec_ktab(i))
          write(txt_ktab(ilen:ilen+j),'(a)')trim(spec_ktab(i))
          ilen=ilen+j+2
       enddo
       k=ncharmodel

       if(solvitemp)then
          txt_tsol='computed'
       else
          txt_tsol=trim(initfile)
       endif

       txt_tsel(1:3)='No '
       if(self_temp)txt_tsel(1:3)='Yes'

       txt_diff(1:3)='No '
       if(diffusion)txt_diff(1:3)='Yes'

       txt_phot(1:3)='No '
       if(photochem)txt_phot(1:3)='Yes'

       txt_deta(1:3)='No '
       if(detached_convect)txt_deta(1:3)='Yes'

       txt_stea(1:3)='No '
       if(steady)txt_stea(1:3)='Yes'

       lammin=1.0d4/wavemax
       lammax=1.0d4/wavemin
       if(wavemin.eq.0.0d0.and.wavemax.eq.0.0d0)then
          lammin=0.0d0
          lammax=0.0d0
       endif

c      ...write output temperature and composition at initial and intermediate times
       outfile(1:nmaxcharfile)=''
       outfile(1:k+11)=namemodel(1:k)//'_out_00.dat'
       if(itime.lt.10)write(outfile(k+7:k+7),'(i1)')itime
       if(itime.ge.10)write(outfile(k+6:k+7),'(i2)')itime
       open(unit=3,file=outfile,status='unknown')
       write(3,1000)rplanet/rearth,mplanet/mearth,tintplanet,
     # albedoplanet,emplanet,fhrd,zenithangle*360.0d0/2.0d0/pi,
     # rstar/rsun,dist/au,trim(starfile),tstar,nz,pressure(nz)/1.0d6,
     # pressure(1)/1.0d6,trim(specfile),trim(thermfile),lammin,lammax,
     # respow,nspec_ktab,trim(txt_ktab),tplanet_th,tplanet,
     # pressure(jadiab)/1.0d6,txt_deta,conv,trim(reacfile),
     # trim(txt_tsol),trim(txt_tsel),trim(txt_diff),trim(txt_phot),time,
     # varabun,trim(txt_stea)
       write(3,1010)(i,i=1,nspec+3)
       write(3,1020)(spec(i),i=1,nspec)
       do j=nz,1,-1
          sum=0.0d0                              ! normalize abundances
          do i=1,nspec
             sum=sum+abun(i,j)
          enddo
          write(3,1030)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #    (dlog10(max(abunmin,abun(i,j)/sum)),i=1,nspec)
       enddo
       close(3)

       if(.not.write_supp)return

c      ...write output temperature and composition at steady state or best steady-state guess
       if(.not.initial)then
          outfile(1:nmaxcharfile)=''
          outfile(1:k+11)=namemodel(1:k)//'_steady.dat'
          open(unit=3,file=outfile,status='unknown')
          write(3,1000)rplanet/rearth,mplanet/mearth,tintplanet,
     #    albedoplanet,emplanet,fhrd,zenithangle*360.0d0/2.0d0/pi,
     #    rstar/rsun,dist/au,trim(starfile),tstar,nz,pressure(nz)/1.0d6,
     #    pressure(1)/1.0d6,trim(specfile),trim(thermfile),lammin,
     #    lammax,respow,nspec_ktab,trim(txt_ktab),tplanet_th,tplanet,
     #    pressure(jadiab)/1.0d6,txt_deta,conv,trim(reacfile),
     #    trim(txt_tsol),trim(txt_tsel),trim(txt_diff),trim(txt_phot),
     #    time,varabun,trim(txt_stea)
          write(3,1010)(i,i=1,nspec+3)
          write(3,1020)(spec(i),i=1,nspec)
          do j=nz,1,-1
             sum=0.0d0                              ! normalize abundances
             do i=1,nspec
                sum=sum+abun(i,j)
             enddo
             write(3,1030)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #       (dlog10(max(abunmin,abun(i,j)/sum)),i=1,nspec)
          enddo
          close(3)
       endif

c      ...write location where opacity at UV wavelengths is 1 as a function of wavelength
       if(.not.initial.and.photochem)then
          call compute_species_tau_uv
          outfile(1:nmaxcharfile)=''
          outfile(1:k+16)=namemodel(1:k)//'_supp_tauuv1.dat'
          open(unit=3,file=outfile,status='unknown')
          write(3,1040)time,(trim(spec(idphotospec(i))),i=1,nphotospec)
          do l=1,nwlth
             write(3,1050)wlth(l),
     #       (p_tauuv1(i,l)/1.0d6,i=nphotospec+1,nphotospec+3),
     #       (p_tauuv1(i,l)/1.0d6,i=1,nphotospec)
          enddo
          close(3)
       endif

c      ...write eddy coefficient profile
       if(.not.initial.and.diffusion)then
          outfile(1:nmaxcharfile)=''
          outfile(1:k+13)=namemodel(1:k)//'_supp_kzz.dat'
          open(unit=3,file=outfile,status='unknown')
          write(3,1060)time
          do j=nz,1,-1
             write(3,1070)pressure(j)/1.0d6,kdif(j)
          enddo
          close(3)
       endif

       if(     initial.and..not.solvitemp)return
       if(.not.initial.and..not.self_temp)return

c      ...write bolometric flux at each layer interface
       outfile(1:nmaxcharfile)=''
       outfile(1:k+14)=namemodel(1:k)//'_supp_flux.dat'
       open(unit=3,file=outfile,status='unknown')
       write(3,1080)
       do n=nz+1,1,-1
          write(3,1090)plev(n)/1.0d6,fnet(n),fup(n),fdown(n)
       enddo
       close(3)

c      ...write location where opacity at IR wavelengths is 1 for each species and wavenumber bin
       call compute_species_tau_ir
       outfile(1:nmaxcharfile)=''
       outfile(1:k+16)=namemodel(1:k)//'_supp_tauir1.dat'
       open(unit=3,file=outfile,status='unknown')
       write(3,1100)time,(trim(spec_ktab(i)),i=1,nspec_ktab)
       do l=1,nwave
          write(3,1110)wave(l),
     #    (p_tauir1(i,l)/1.0d6,i=nspec_ktab+1,nspec_ktab+3),
     #    (p_tauir1(i,l)/1.0d6,i=1,nspec_ktab)
       enddo
       close(3)

c      ...write radiative/convective character of layers
       outfile(1:nmaxcharfile)=''
       outfile(1:k+17)=namemodel(1:k)//'_supp_convect.dat'
       open(unit=3,file=outfile,status='unknown')
       write(3,1120)time
       do j=nz,1,-1
          i=0
          if(convect(j))i=1
          write(3,1130)pressure(j)/1.0d6,i
       enddo
       close(3)

c      ...write cooling/heating rates
c       call cool_heat_rates

       return

1000   format(
     # '!',/,
     # '! 1D self-consistent atmosphere model ',/,
     # '! Radiative-convective temperature ',
     # 'and non-equilibrium chemistry',/,
     # '!',/,
     # '! Planet radius                        = ',1pg10.3,' REarth',/,
     # '! Planet mass                          = ',1pg10.3,' MEarth',/,
     # '! Planet internal temperature          = ',0pf10.2,' K',/,
     # '! Planet surface albedo                = ',0pf10.2,/,
     # '! Planet surface emissivity            = ',0pf10.2,/,
     # '! Heat redistribution factor           = ',1pg10.3,/,
     # '! Zenith angle                         = ',1pg10.2,' deg',/,
     # '! Star radius                          = ',1pg10.3,' Rsun',/,
     # '! Star-planet distance                 = ',1pg10.3,' au',/,
     # '! Star spectrum file                   : ',a,/,
     # '! Star effective temperature           = ',0pf10.2,' K',/,
     # '! Number of layers                     = ',i3,/,
     # '! Pressure range                       = ',
     #    1pg8.2,'-',1pg8.2,' bar',/,
     # '! Species file                         : ',a,/,
     # '! NASA polynomials file                : ',a,/,
     # '!',/,
     # '! Radiative convective model parameters:',/,
     # '!',/,
     # '! Wavelength grid                      = ',
     #    1pg8.2,'-',1pg8.2,' um',/,
     # '! Resolution (R=lambda/Dlambda)        = ',1pg10.2,/,
     # '! Number of opacity species            = ',i3,/,
     # '! ',a,/,
     # '! Planet effective temperature (theor) = ',0pf10.2,' K',/,
     # '! Planet effective temperature         = ',0pf10.2,' K',/,
     # '! Radiative-convective transition      = ',1pg10.2,' bar',/,
     # '! Detached convective zone?            : ',a,/,
     # '! Temperature computation convergence  = ',1pg10.2,/,
     # '!',/,
     # '! Non-equilibrium chemistry parameters:',/,
     # '!',/,
     # '! Reaction file                        : ',a,/,
     # '! Initial temperature                  : ',a,/,
     # '! Self-consistent temperature included?: ',a,/,
     # '! Diffusion included?                  : ',a,/,
     # '! Photochemistry included?             : ',a,/,
     # '! Continuity equation integrated up to =',1pg12.4,' s',/,
     # '!   with mole fraction relative error  =',1pg12.4,/,
     # '! Steady state reached?                : ',a,/,
     # '!',/,
     # '! .......................................................',/,
     # '! z,p,T, and mixing ratios (decimal logarithm):',/,
     # '!')
1010   format('!',3(7x,i2,7x),4x,999(6x,i3,6x))
1020   format('    height[km]     pressure[bar]    temperature[K] ',
     # 8x,999(a15))
1030   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,999(5x,0pf10.6))
1040   format(
     # '!',/,
     # '! pressure (bar) where UV optical depth is 1 at time:',
     # 1pg12.4,' s',/,
     # '!',/,
     # 1x,'  wlth(nm)   ','   absorption   ','   scattering   ',
     # '      total     ',999(1x,a15))
1050   format(2x,1pg12.4,3x,999(4x,1pg12.4))
1060   format(
     # '!',/,
     # '! eddy diffusion profile at time:',
     # 1pg12.4,' s',/,
     # '!',/,
     # '    pressure[bar]    Kzz(cm2_s-1)')
1070   format(3x,1pg12.4,6x,1pg12.4)
1080   format('! flux at layer interfaces',/,
     # '    pressure[bar]    flux_net[erg_s-1_cm-2]    ',
     # ' flux_up[erg_s-1_cm-2]     flux_down[erg_s-1_cm-2] ')
1090   format(3x,1pg12.4,3(12x,1pg12.4))
1100   format(
     # '!',/,
     # '! pressure (bar) where IR optical depth is 1 at time:',
     # 1pg12.4,' s',/,
     # '!',/,
     # 1x,'wavenumber(cm-1)','   absorption   ','   scattering   ',
     # '      total     ',999(1x,a15))
1110   format(2x,0pf12.5,3x,999(4x,1pg12.4))
1120   format(
     # '!',/,
     # '! radiative/convective character of layers at time:',
     # 1pg12.4,' s',/,
     # '!',/,
     # '    pressure[bar]    radiative=0/convective=1')
1130   format(3x,1pg12.4,15x,i1)

       end




c_______________________________________________________________________

       subroutine compute_species_tau_ir

       include 'pact.common'
       integer i,j,k,l,ii,ip1,it1
       real*8 xg(nmaxg),wg(nmaxg),kp1t1,kp1t2,kp2t1,kp2t2,t,u,
     # ktot(nmaxg),transm,pwk,twk,tau,o_tau,pu,pl,tu,tl,slope

       o_tau=0.0d0           ! to avoid warning

c      ...get Gaussian-Legendre quadrature points and weights
       call gauleg(0.0d0,0.99d0,xg(1:20),wg(1:20),20)      ! 20 points in the range g = 0    -0.99
       call gauleg(0.99d0,0.999d0,xg(21:30),wg(21:30),10)  ! 10 points in the range g = 0.99 -0.999
       call gauleg(0.999d0,1.0d0,xg(31:36),wg(31:36),6)    !  6 points in the range g = 0.999-1

c      ...loop over species
       do i=1,nspec_ktab

c      ...loop over wavenumber grid
          do l=1,nwave

c      ...loop over atmosphere layers
             tau=0.0d0
             do j=1,nz
                pwk=pressure(j)
                if(pwk.lt.p_ktab(1))      pwk=p_ktab(1)
                if(pwk.gt.p_ktab(np_ktab))pwk=p_ktab(np_ktab)
                twk=tk(j)
                if(twk.lt.t_ktab(1))      twk=t_ktab(1)
                if(twk.gt.t_ktab(nt_ktab))twk=t_ktab(nt_ktab)
      
                ii=1                                   ! identify neighboring pressures in k-tables
                do while(pwk.gt.p_ktab(ii).and.ii.lt.np_ktab)
                   ii=ii+1
                enddo
                ip1=min(max(ii-1,1),np_ktab-1)
      
                ii=1                                   ! identify neighboring temperatures in k-tables
                do while(twk.gt.t_ktab(ii).and.ii.lt.nt_ktab)
                   ii=ii+1
                enddo
                it1=min(max(ii-1,1),nt_ktab-1)

                ktot(:)=0.0d0
                do k=1,nmaxg                     ! get k(g) at neighboring p,T
                   kp1t1=kcoeff(i,ip1,it1,l,k)
                   kp1t2=kcoeff(i,ip1,it1+1,l,k)
                   kp2t1=kcoeff(i,ip1+1,it1,l,k)
                   kp2t2=kcoeff(i,ip1+1,it1+1,l,k)
c                                                ! bilinear interpolation
                   t=(dlog10(pwk)-dlog10(p_ktab(ip1)))/
     #                (dlog10(p_ktab(ip1+1))-dlog10(p_ktab(ip1)))
                   u=(twk-t_ktab(it1))/(t_ktab(it1+1)-t_ktab(it1))
                   ktot(k)=(1.0d0-t)*(1.0d0-u)*kp1t1+
     #                             t*(1.0d0-u)*kp2t1+
     #                                     t*u*kp2t2+
     #                             (1.0d0-t)*u*kp1t2
c                                                ! convert from 'log k' to 'k * abun'
                   ktot(k)=10.0d0**ktot(k)*abun(idspec_ktab(i),j)
                enddo

c      ...compute transmission and add contribution to absorption coefficient [cm-1]
                transm=0.0d0
                do k=1,nmaxg
                   transm=transm+
     #             dexp(-ktot(k)*density(j)*deltaz(j))*wg(k)
                enddo
                transm=max(transm,dexp(-700.0d0))   ! maximum optical depth set to 700 for numerical reasons
                o_tau=tau
                tau=tau-dlog(transm)
                if(tau.gt.1.0d0)exit
             enddo
             if(j.eq.1)then
                p_tauir1(i,l)=pressure(1)
                cycle
             endif
             if(tau.le.1.0d0)then
                p_tauir1(i,l)=pressure(nz)
                cycle
             endif
             pu=pressure(j-1)
             pl=pressure(j)
             tu=o_tau
             tl=tau
             slope=(pu-pl)/(tu-tl)
             p_tauir1(i,l)=pl+slope*(1.0d0-tl)
          enddo
       enddo

c      ...get pressure level where absorption, scattering, and total opacity is 1
       do i=nspec_ktab+1,nspec_ktab+3
          do l=1,nwave
             tau=0.0d0
             do j=1,nz
                if(i.eq.nspec_ktab+1)tau=tau+kabs(l,j)*deltaz(j)
                if(i.eq.nspec_ktab+2)tau=tau+kscat(l,j)*deltaz(j)
                if(i.eq.nspec_ktab+3)
     #          tau=tau+(kabs(l,j)+kscat(l,j))*deltaz(j)
                if(tau.gt.1.0d0)exit
             enddo
             if(j.eq.1)then
                p_tauir1(i,l)=pressure(1)
                cycle
             endif
             if(tau.le.1.0d0)then
                p_tauir1(i,l)=pressure(nz)
                cycle
             endif
             pu=pressure(j-1)
             pl=pressure(j)
             tu=o_tau
             tl=tau
             slope=(pu-pl)/(tu-tl)
             p_tauir1(i,l)=pl+slope*(1.0d0-tl)
          enddo
       enddo

       return

       end




c_______________________________________________________________________

       subroutine compute_species_tau_uv

       use pact_data
       include 'pact.common'       
       integer nmaxit
       parameter(nmaxit=1000)
       real*8 varmin
       parameter(varmin=1.0d-2)                  ! 1 %

       integer i,j,k,l,ii,id,it
       real*8 csrayleigh(nmaxspec),ta,ts,n,dz,sigma,ftop(nmaxwlth),down,
     # up,fx,o_fx,var,tau,tu,tl,pu,pl,slope,r,path
       real*8, allocatable, dimension(:,:) :: tau_abs,tau_sca,
     #                                        fdow,o_fdow,fupw,o_fupw
       real*8, allocatable, dimension(:)   :: polz

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
       allocate(fdow(0:nz+1,nwlth),o_fdow(0:nz+1,nwlth))
       allocate(fupw(0:nz+1,nwlth),o_fupw(0:nz+1,nwlth))

       do ii=1,nphotospec+3

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
c             path=dz/dcos(zenithangle)                       !...plane parallel atmosphere
             path=r*(dsqrt(dcos(zenithangle)**2.0d0+          !...spherical atmosphere
     #       (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
             do l=1,nwlth
                tau_abs(j,l)=0.0d0                  ! absorption
                do k=1,nphotospec
                   if(ii.eq.nphotospec+2)cycle
                   if(ii.le.nphotospec.and.k.ne.ii)cycle
                   id=idphotospec(k)
                   if(j.eq.1)then
                      n=abun(id,j)*density(j)
                   elseif(j.eq.nz+1)then
                      n=abun(id,j-1)*density(j-1)
                   else
                      n=0.5d0*
     #                (abun(id,j-1)*density(j-1)+abun(id,j)*density(j))
                   endif
                   tau_abs(j,l)=tau_abs(j,l)+csphotabs(k,l)*n*path
                enddo
                tau_sca(j,l)=0.0d0                  ! scattering
                do i=1,nspec
                   if(ii.lt.nphotospec+2)cycle
                   if(j.eq.1)then
                      n=abun(i,j)*density(j)
                   elseif(j.eq.nz+1)then
                      n=abun(i,j-1)*density(j-1)
                   else
                      n=0.5d0*
     #                (abun(i,j-1)*density(j-1)+abun(i,j)*density(j))
                   endif
                   sigma=csrayleigh(i)/(wlth(l)*1.0d-7)**4.0d0
                   tau_sca(j,l)=tau_sca(j,l)+sigma*n*path
                enddo
             enddo
          enddo

c      ...compute flux at each layer with a two-stream iterative algorithm
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
     #             fdow(j-1,l))*dexp(-ta-ts)
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
     #             fupw(j+1,l))*dexp(-ta-ts)
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

c      ...locate tau=1 height as a function of wavelength
c         opacity is defined with respect to the non-attenuated incoming stellar flux
c         and with respect to the downward flux at each layer
          do l=1,nwlth
             do j=1,nz
                tau=dlog(fphot(l)/(fdow(j,l)))
                if(tau.gt.1.0d0)exit
             enddo
             if(j.eq.1)then
                p_tauuv1(ii,l)=pressure(1)
                cycle
             endif
             if(tau.le.1.0d0)then
                p_tauuv1(ii,l)=pressure(nz)
                cycle
             endif
             tu=dlog(fphot(l)/(fdow(j-1,l)))
             tl=dlog(fphot(l)/(fdow(j,l)))
             pu=pressure(j-1)
             pl=pressure(j)
             slope=(pu-pl)/(tu-tl)
             p_tauuv1(ii,l)=pl+slope*(1.0d0-tl)
          enddo
       enddo

       deallocate(tau_abs,tau_sca,fdow,o_fdow,fupw,o_fupw)

       return

1000   format(1x,
     # ' W- Rayleigh scattering two-stream algorithm not converged',/,
     # '    Reached accuracy = ',1pg10.3)

       end




c_______________________________________________________________________

       subroutine cool_heat_rates

       include 'pact.common'       
       integer i,j,l,id,n
       real*8 tkin,sum,dh,g,h,s,cp,dn,chrate(nmaxz),cpmean
       character*(nmaxcharfile) outfile

c Loop over height
       do j=1,nz
          sum=0.0d0                              ! normalize abundances
          do i=1,nspec
             sum=sum+abun(i,j)
          enddo
          tkin=tk(j)
          chrate(j)=0.0d0

c ---  Loop over chemical reactions: add chemical contribution to cooling/heating rate
          do i=1,nreac
             dh=0.0d0 ! DeltaH/RT                !...evaluate reaction enthalpy
             do l=1,nreag(i)
                id=idreag(i,l)
                call evaluate_therm(id,tkin,g,h,s,cp)
                dh=dh-h
             enddo
             do l=1,nprod(i)
                id=idprod(i,l)
                call evaluate_therm(id,tkin,g,h,s,cp)
                dh=dh+h
             enddo

             dn=krate(i,j)                       !...evaluate reaction rate
             do l=1,nreag(i)
                dn=dn*max(abunmin,abun(idreag(i,l),j)/sum)*density(j)
             enddo   

             chrate(j)=chrate(j)-                ! add contribution: cooling(-), heating(+)
     #       dn*dh*kboltz*tkin                   ! [=] erg cm-3 s-1
          enddo

c ---  Loop over photochemical reactions:  add photochemical contribution to cooling/heating rate
c      A break of a molecule by photodissociation does not cool the atmosphere
c      since the heat taken for breaking the bond(s) comes from UV photons
c          do i=1,nphotoreac
c             dh=0.0d0 ! DeltaH/RT                !...evaluate photoreaction enthalpy
c             id=idphotoreag(i)
c             call evaluate_therm(id,tkin,g,h,s,cp)
c             dh=dh-h
c             do l=1,nphotoprod(i)
c                id=idphotoprod(i,l)
c                call evaluate_therm(id,tkin,g,h,s,cp)
c                dh=dh+h
c             enddo
c
c             dn=jphot(i,j)*                      !...evaluate photoreaction rate
c     #       max(abunmin,abun(idreag(i,l),j)/sum)
c
c             chrate(j)=chrate(j)-                ! add contribution: cooling(-), heating(+)
c     #       dn*dh*kboltz*tkin                   ! [=] erg cm-3 s-1
c          enddo

c      ...convert cooling/heating rate to units of [K s-1]
          cpmean=0.0d0 ! Cp/R                    ! evaluate mean heat capacity
          do i=1,nspec
             call evaluate_therm(i,tkin,g,h,s,cp)
             cpmean=cpmean+cp*max(abunmin,abun(i,j)/sum)
          enddo
          chrate(j)=chrate(j)/
     #    (cpmean*kboltz*density(j))             ! [=] K s-1
       enddo

c      ...write out cooling/heating rate
       n=ncharmodel
       outfile(1:nmaxcharfile)=''
       outfile(1:n)=namemodel(1:n)
       outfile(n+1:n+11)='_chrate.dat'
       open(unit=9,file=outfile,status='unknown')
       write(9,1000)time
       do j=nz,1,-1
          write(9,1010)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),chrate(j)
       enddo
       close(9)

1000   format('! cooling/heating rate at time:',1pg12.4,' s',/,
     #        '! cooling (-) and heating (+)',/,
     #        '!','   height[km]   ',
     #            ' pressure[bar]  ',
     #            ' temperature[K] ',
     #            ' cool/heat rate[K s-1]')
1010   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,1pg12.4)

       end
