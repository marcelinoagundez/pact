
       subroutine solve_radiative_transfer_stellar

       include 'pact.common'

       integer i,j,l,n
       real*8 mu0,sqrt3,mu1,f,fac,facp,facm,emlt,ssfc
       real*8, allocatable, dimension(:)   :: tauv,tau,w0,g,tauc,
     # direct,gam1,gam2,gam3,gam4,lam,gam,cp0,cm0,cpb,cmb,e1,e2,e3,e4,
     # a,b,d,e,y,y1,y2,flux_u,flux_d

       allocate(tauv(nz))
       allocate(tau(nz),w0(nz),g(nz),tauc(nz+1),direct(nz+1))
       allocate(gam1(nz),gam2(nz),gam3(nz),gam4(nz),lam(nz),gam(nz))
       allocate(cp0(nz),cm0(nz),cpb(nz),cmb(nz))
       allocate(e1(nz),e2(nz),e3(nz),e4(nz))
       allocate(a(2*nz),b(2*nz),d(2*nz),e(2*nz),y(2*nz))
       allocate(y1(nz),y2(nz),flux_u(nz+1),flux_d(nz+1))

       mu0=dcos(zenithangle)
       sqrt3=dsqrt(3.0d0)
c       mu1=0.5d0                    ! Eddington
       mu1=1.0d0/sqrt3              ! Quadrature

c      ...if no star
       if(trim(starfile).eq.'None'.and.tstar.lt.1.0d-3)then
          fup_st(:,:)=0.0d0
          fdo_st(:,:)=0.0d0
          return
       endif

c      ..Loop on wavenumber grid
       do i=nwave,1,-1
          flux_u(:)=0.0d0
          flux_d(:)=0.0d0

c      ...get opacity
          do j=1,nz
             tauv(j)=(kabs(i,j)+kscat(i,j))*deltaz(j)
          enddo

c      ...calculate coefficients gamma's, lambda, and Gamma
          do n=1,nz
             j=n

c            ...delta scalings (Joseph et al 1976, J Atm Sci, 33, 2452)
             f=gscat(i,j)*gscat(i,j)
             tau(n)=(1.0d0-walb(i,j)*f)*tauv(j)
             w0(n)=(1.0d0-f)*walb(i,j)/(1.0d0-walb(i,j)*f)
             g(n)=gscat(i,j)/(1.0d0+gscat(i,j))

c            ...Eddington
c             gam1(n)=(7.0d0-w0(n)*(4.0d0+3.0d0*g(n)))/4.0d0
c             gam2(n)=-(1.0d0-w0(n)*(4.0d0-3.0d0*g(n)))/4.0d0
c             gam3(n)=(2.0d0-3.0d0*g(n)*mu0)/4.0d0
c             gam4(n)=1.0d0-gam3(n)

c            ...Quadrature
             gam1(n)=sqrt3*(2.0d0-w0(n)*(1.0d0+g(n)))/2.0d0
             gam2(n)=sqrt3*w0(n)*(1.0d0-g(n))/2.0d0
             gam3(n)=(1.0d0-sqrt3*g(n)*mu0)/2.0d0
             gam4(n)=1.0d0-gam3(n)

             lam(n)=dsqrt(dabs(gam1(n)*gam1(n)-gam2(n)*gam2(n)))
             gam(n)=(gam1(n)-lam(n))/gam2(n)
             if(gam2(n).eq.0.0d0)gam(n)=gam2(n)/(gam1(n)+lam(n)) ! avoid possible floating point
          enddo

c      ...calculate coefficients e's
          do n=1,nz
             emlt=dexp(-lam(n)*tau(n))
             e1(n)=1.0d0+gam(n)*emlt
             e2(n)=1.0d0-gam(n)*emlt
             e3(n)=gam(n)+emlt
             e4(n)=gam(n)-emlt
          enddo

c      ...calculate coefficients A, B, and D of left-hand side of matrix equation
          l=1                          ! l=1 (top of the atmosphere)
          a(l)=0.0d0
          b(l)=e1(1)
          d(l)=-e2(1)

          do n=1,nz-1                  ! l odd   ! mistake in Toon et al (1989) for running index l=2n-1 in Eq (40)
             l=2*n+1                             ! seen in code ClimaCO2 (delta2str.f)     
             a(l)=e2(n)*e3(n)-e4(n)*e1(n)
             b(l)=e1(n)*e1(n+1)-e3(n)*e3(n+1)
             d(l)=e3(n)*e4(n+1)-e1(n)*e2(n+1)
          enddo

          do n=1,nz-1                  ! l even
             l=2*n
             a(l)=e2(n+1)*e1(n)-e3(n)*e4(n+1)
             b(l)=e2(n)*e2(n+1)-e4(n)*e4(n+1)
             d(l)=e1(n+1)*e4(n+1)-e2(n+1)*e3(n+1)
          enddo

          l=2*nz                       ! l=2*nz (bottom of the atmosphere)
          a(l)=e1(nz)-albedoplanet*e3(nz)
          b(l)=e2(nz)-albedoplanet*e4(nz)
          d(l)=0.0d0

c      ...calculate accumulative optical depth (tauc) and flux from stellar beam (direct) at each level
          do n=1,nz+1                            ! tauc(n) is optical depth above layer n
             if(n.eq.1)then                      ! direct(n) is stellar beam flux above layer n
                tauc(n)=0.0d0
             else
                tauc(n)=tauc(n-1)+tau(n-1)
             endif
             direct(n)=mu0*fhrd*fstar(i)*dexp(-tauc(n)/mu0)
          enddo

c      ...calculate coefficients C's
          do n=1,nz
             fac=w0(n)*fhrd*fstar(i)/(lam(n)*lam(n)-1.0d0/(mu0*mu0))
             facp=fac*((gam1(n)-1.0d0/mu0)*gam3(n)+gam4(n)*gam2(n))
             facm=fac*((gam1(n)+1.0d0/mu0)*gam4(n)+gam2(n)*gam3(n))
             cp0(n)=facp*dexp(-tauc(n)/mu0)
             cm0(n)=facm*dexp(-tauc(n)/mu0)
             cpb(n)=facp*dexp(-(tauc(n)+tau(n))/mu0)
             cmb(n)=facm*dexp(-(tauc(n)+tau(n))/mu0)
          enddo
          ssfc=albedoplanet*direct(nz+1)

c      ...calculate coefficients E of right-hand side of matrix equation
          l=1                          ! l=1 (top of the atmosphere)
          e(l)=-cm0(1)

          do n=1,nz-1                  ! l odd   ! mistake in Toon et al (1989) for running index l=2n-1 in Eq (40)
             l=2*n+1                             ! seen in code ClimaCO2 (delta2str.f)     
             e(l)=e3(n)*(cp0(n+1)-cpb(n))+e1(n)*(cmb(n)-cm0(n+1))
          enddo

          do n=1,nz-1                  ! l even                      ! mistake in Toon et al (1989) for E_l, last line in Eq (42)
             l=2*n                                                   ! see p 451 of Atmospheric Evolution of Inhabited and Lifeless Worlds (David C. Catling, James F. Kasting)
             e(l)=e2(n+1)*(cp0(n+1)-cpb(n))+e4(n+1)*(cmb(n)-cm0(n+1))! seen in codes ClimaCO2 (delta2str.f) and Streamer (toon.f)
          enddo

          l=2*nz                       ! l=2*nz (bottom of the atmosphere)
          e(l)=ssfc-cpb(nz)+albedoplanet*cmb(nz)

c      ...solve tridiagonal matrix
          call tridag(a,b,d,e,y,2*nz)
          do n=1,nz
             l=2*n-1         ! l odd  -> Y1
             y1(n)=y(l)
             l=2*n           ! l even -> Y2
             y2(n)=y(l)
          enddo

c      ...calculate upward/downward fluxes
          flux_u(1)  =y1(1)*e3(1)-y2(1)*e4(1)+cp0(1)
          flux_d(1)=direct(1)
          do n=1,nz
             flux_u(n+1)=y1(n)*e1(n)+y2(n)*e2(n)+cpb(n)
             flux_d(n+1)=y1(n)*e3(n)+y2(n)*e4(n)+cmb(n)+direct(n+1)
          enddo

c      ...store upward and downward fluxes
          do n=1,nz+1
             fup_st(i,n)=flux_u(n)
             fdo_st(i,n)=flux_d(n)
          enddo
       enddo
c      ..End of loop on wavenumber grid

c      ...deallocate variables
       deallocate(tauv,tau,w0,g,tauc,direct,gam1,gam2,gam3,gam4,lam,
     # gam,cp0,cm0,cpb,cmb,e1,e2,e3,e4,a,b,d,e,y,y1,y2,flux_u,flux_d)

       return

       end
