
       subroutine compute_initial_temperature

       include 'pact.common'
       integer j,ihou,imin,isec
       real*8 eps
       real*8, allocatable, dimension(:) :: tk_old

       eps=1.0d0
       allocate(tk_old(nz))

c      ...get initial temperature profile
       call get_analytical_temperature
       call adjust_vertical_structure

c      ...load opacity data and star spectrum
       call load_cia
       call read_ktable
       call get_star_spectrum

c      ...solve temperature under chemical equilibrium or not
       do bigit=1,nmaxbigit
          if(iabun.eq.0)call compute_chemical_equilibrium
          call solve_opacity
          tk_old(1:nz)=tk(1:nz)
          call solve_temperature
          call adjust_vertical_structure

          dtkmax=0.0d0
          do j=1,nz
             dtkmax=max(dtkmax,dabs(tk(j)-tk_old(j)))
          enddo
          write(*,1000)bigit,dtkmax
          if(bigit.gt.1.and.dtkmax.lt.eps)exit
       enddo

       if(bigit.eq.nmaxbigit)write(*,1010)bigit

c      ...write results at time=0
       call write_results
       call get_cputime(ihou,imin,isec)
       write(*,1020)ihou,imin,isec

       deallocate(tk_old)

       return

1000   format(/,2x,'I- Loop ',i2,' - Temperature change =',0pf10.2,' K')
1010   format(2x,'W- Temperature convergence not reached in ',i2,
     # ' iterations')
1020   format(2x,'I- CPU=',i3,'h 'i2,'m ',i2,'s')

       end




c_______________________________________________________________________

       subroutine get_analytical_temperature
c      ...analytical solution for irradiated exoplanet atmosphere
c         Eq (27) of Guillot (2010), A&A, 520, A27

       include 'pact.common'
       integer j
       real*8 tirr,mu0,kvis,kir,gam,tau,rho,g,h
       real*8, allocatable, dimension(:)   :: tklev

       allocate(tklev(nz+1))

       tirr=tstar*dsqrt(rstar/dist)
       mu0=dcos(zenithangle)
c                            ! ...values from Guillot (2010), typical of HD209458b
       kvis=4.0d-3           ! mass absorption coefficient at visible wavelengths  [=] cm2 g-1
c       kvis=1.0d-2           ! mass absorption coefficient at visible wavelengths  [=] cm2 g-1
       kir =1.0d-2           ! mass absorption coefficient at infrared wavelengths [=] cm2 g-1
       gam=kvis/kir
       if(iabun.eq.0)then    ! compute average particle mass if chemical equilibrium used for initial composition
          apmass=2.3d0*amu   ! ...typical of H2/He atmosphere; of no importance here; (p,T) profile not sensitive to it
       else                  ! compute average particle mass if prescribed initial composition
          call compute_average_particle_mass
       endif

       tau=0.0d0
       do j=1,nz+1

c      ...Guillot (2010)
          tklev(j)=3.0d0/4.0d0*tintplanet**4.0d0*(2.0d0/3.0d0+tau)+
     #       3.0d0/4.0d0*tirr**4.0d0*mu0*
     #       (2.0d0/3.0d0+mu0/gam+
     #       (gam/3.0d0/mu0-mu0/gam)*dexp(-gam*tau/mu0))
          tklev(j)=tklev(j)**0.25d0

c      ...Hansen (2008)
c          tklev(j)=3.0d0/4.0d0*tintplanet**4.0d0*(2.0d0/3.0d0+tau)+
c     #       tirr**4.0d0*mu0*
c     #       (1.0d0+3.0d0/2.0d0*(mu0/gam)**2.0d0+
c     #       -3.0d0/2.0d0*(mu0/gam)**3.0d0*dlog(1.0d0+gam/mu0)
c     #       -3.0d0/4.0d0*mu0/gam*dexp(-gam*tau/mu0))
c          tklev(j)=tklev(j)**0.25d0

          if(j.eq.nz+1)exit
          rho=pressure(j)/kboltz/tklev(j)        ! use tk(upper interface) as proxy of tk(layer)
          g=ggravity*mplanet/rplanet**2.0d0
          h=kboltz*tklev(j)/apmass/g             ! use tk(upper interface) as proxy of tk(layer)
          deltaz(j)=-h*dlog(plev(j)/plev(j+1))
          tau=tau+kir*apmass*rho*deltaz(j)
       enddo

       do j=1,nz
          tk(j)=0.5d0*(tklev(j)+tklev(j+1))
       enddo

       deallocate(tklev)

       return

       end
