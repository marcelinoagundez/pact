
       subroutine solve_radiative_transfer_infrared

       include 'pact.common'

       integer i,j,l,n,j0,jmax,it
       real*8 mu0,mu1,slope,zl,nu,fac,emlt,ssfc,mu,dmu,emtm,eps,err,
     # maxerr
       real*8 bnuplanck
       real*8, allocatable, dimension(:)   :: tkl,tauv,tau,w0,g,
     # gam1,gam2,lam,gam,cp0,cm0,cpb,cmb,e1,e2,e3,e4,b0,b1,a,b,d,e,
     # y,y1,y2,flux_u,flux_d,gg,hh,jj,kk,alp1,alp2,sig1,sig2,iup,idown,
     # sumup,sumdown,old_flux_u,old_flux_d

       allocate(tauv(nz),tkl(nz+1),b0(nz),b1(nz))
       allocate(tau(nz),w0(nz),g(nz))
       allocate(gam1(nz),gam2(nz),lam(nz),gam(nz))
       allocate(cp0(nz),cm0(nz),cpb(nz),cmb(nz))
       allocate(e1(nz),e2(nz),e3(nz),e4(nz))
       allocate(a(2*nz),b(2*nz),d(2*nz),e(2*nz),y(2*nz))
       allocate(y1(nz),y2(nz),flux_u(nz+1),flux_d(nz+1))
       allocate(gg(nz),hh(nz),jj(nz),kk(nz))
       allocate(alp1(nz),alp2(nz),sig1(nz),sig2(nz))
       allocate(iup(nz+1),idown(nz+1))
       allocate(sumup(nz+1),sumdown(nz+1))
       allocate(old_flux_u(nz+1),old_flux_d(nz+1))

       mu0=dcos(zenithangle)
       mu1=0.5d0                                 ! Hemispheric mean
c       jmax=15                                   ! parameters for integration over mu
       jmax=10                                   ! parameters for integration over mu
c       eps=1.0d-3                                ! ...
       eps=1.0d-2                                ! ...

c      ...calculate temperature at levels between layers
       do n=1,nz+1                               ! tkl(n) is temperature at top of layer n
          j=n
          j0=j-1
          if(j0.eq.0)j0=1
          if(j0.eq.nz)j0=nz-1
          zl=z(j)+0.5d0*deltaz(j)                ! top of layer
          if(j.eq.nz+1)zl=z(nz)-0.5d0*deltaz(nz) ! bottom of atmosphere
          slope=(tk(j0+1)-tk(j0))/(z(j0+1)-z(j0))
          tkl(n)=tk(j0)+slope*(zl-z(j0))
       enddo

c      ..Loop on wavenumber grid
       do i=nwave,1,-1
          nu=wave(i)*clight
          flux_u(:)=0.0d0
          flux_d(:)=0.0d0
          old_flux_u(:)=0.0d0
          old_flux_d(:)=0.0d0

c      ...get opacity
          do j=1,nz
             tauv(j)=(kabs(i,j)+kscat(i,j))*deltaz(j)
          enddo

c      ...calculate coefficients gamma's, lambda, and Gamma
          do n=1,nz
             j=n

c            ...delta scalings (no scaling for the infrared)
             tau(n)=tauv(j)
             w0(n)=walb(i,j)
             g(n)=gscat(i,j)

c            ...Hemispheric mean
             gam1(n)=2.0d0-w0(n)*(1.0d0+g(n))
             gam2(n)=w0(n)*(1.0d0-g(n))

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

c      ...calculate coefficients C's
          do n=1,nz
             b0(n)=bnuplanck(tkl(n),nu)          ! Planck function at top of layer n
             b1(n)=(bnuplanck(tkl(n+1),nu)-b0(n))/tau(n)
             fac=2.0d0*pi*mu1
             cp0(n)=fac*(b0(n)+b1(n)*(1.0d0/(gam1(n)+gam2(n))))
             cpb(n)=fac*(b0(n)+b1(n)*(tau(n)+1.0d0/(gam1(n)+gam2(n))))
             cm0(n)=fac*(b0(n)+b1(n)*(-1.0d0/(gam1(n)+gam2(n))))
             cmb(n)=fac*(b0(n)+b1(n)*(tau(n)-1.0d0/(gam1(n)+gam2(n))))
          enddo
          ssfc=emplanet*pi*bnuplanck(tkl(nz+1),nu)

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

c      ...skip source function technique to speed-up in 1st iteration
          if(bigit.eq.1)then
             flux_u(1)=y1(1)*e3(1)-y2(1)*e4(1)+cp0(1) ! calculate upward/downward fluxes
             flux_d(1)=0.0d0
             do n=1,nz
                flux_u(n+1)=y1(n)*e1(n)+y2(n)*e2(n)+cpb(n)
                flux_d(n+1)=y1(n)*e3(n)+y2(n)*e4(n)+cmb(n)
             enddo

             do n=1,nz+1                              ! store upward and downward fluxes
                fup_ir(i,n)=flux_u(n)
                fdo_ir(i,n)=flux_d(n)
             enddo
             cycle
          endif

c      ...calculate coefficients G, H, J, K, alpha's, and sigma's
          do n=1,nz
             gg(n)=(y1(n)+y2(n))*(1.0d0/mu1-lam(n))
             hh(n)=(y1(n)-y2(n))*gam(n)*(lam(n)+1.0d0/mu1)
             jj(n)=(y1(n)+y2(n))*gam(n)*(lam(n)+1.0d0/mu1)
             kk(n)=(y1(n)-y2(n))*(1.0d0/mu1-lam(n))
             alp1(n)=2.0d0*pi*(b0(n)+b1(n)*
     #          (1.0d0/(gam1(n)+gam2(n))-mu1))
             alp2(n)=2.0d0*pi*b1(n)
             sig1(n)=2.0d0*pi*(b0(n)-b1(n)*
     #          (1.0d0/(gam1(n)+gam2(n))-mu1))
             sig2(n)=2.0d0*pi*b1(n)
          enddo

c      ...Calculate upward and downward fluxes using two-stream source function technique
          do j=1,jmax                  ! integrate I*mu over mu to compute fluxes
             it=2**(j-2)               ! use extended trapezoidal rule
             dmu=1.0d0/(dble(it))      ! algorithm qtrap/trapzd of Numerical Recipes implemented here
             mu=0.5d0*dmu
             if(j.eq.1)then            ! on first refinement (j=1) skip mu=0 and do only mu=1
                it=1
                mu=1.0d0
             endif
             sumup(1:nz+1)=0.0d0
             sumdown(1:nz+1)=0.0d0
             do l=1,it
                idown(1)=0.0d0                             ! evaluate downward intensities
                do n=1,nz                                  ! sum done in steps to avoid possible floating points
                   emtm=dexp(-tau(n)/mu)
                   idown(n+1)=idown(n)*emtm
                   if(jj(n).ne.0.0d0)idown(n+1)=idown(n+1)+
     #             jj(n)/(lam(n)*mu+1.0d0)*
     #             (1.0d0-dexp(-tau(n)*(lam(n)+1.0d0/mu)))
                   if(kk(n).ne.0.0d0)idown(n+1)=idown(n+1)+
     #             kk(n)/(lam(n)*mu-1.0d0)*(emtm-dexp(-tau(n)*lam(n)))
                   idown(n+1)=idown(n+1)+
     #             sig1(n)*(1.0d0-emtm)+sig2(n)*(mu*emtm+tau(n)-mu)
                enddo

                iup(nz+1)=ssfc/mu1+albedoplanet*idown(nz+1)! evaluate upward intensities
                do n=nz,1,-1                               ! sum done in steps to avoid possible floating points
                   emtm=dexp(-tau(n)/mu)
                   iup(n)=iup(n+1)*emtm
                   if(gg(n).ne.0.0d0)iup(n)=iup(n)+
     #             gg(n)/(lam(n)*mu-1.0d0)*(emtm-dexp(-tau(n)*lam(n)))
                   if(hh(n).ne.0.0d0)iup(n)=iup(n)+
     #             hh(n)/(lam(n)*mu+1.0d0)*
     #             (1.0d0-dexp(-tau(n)*(lam(n)+1.0d0/mu)))
                   iup(n)=iup(n)+
     #             alp1(n)*(1.0d0-emtm)+alp2(n)*(mu-(tau(n)+mu)*emtm)
                enddo

                do n=1,nz+1
                   sumup(n)=sumup(n)+mu*iup(n)
                   sumdown(n)=sumdown(n)+mu*idown(n)
                enddo
                mu=mu+dmu
             enddo

             do n=1,nz+1               ! evaluate integrals for first refinement (j=1)
                if(j.eq.1)then
                   flux_u(n)=0.5d0*sumup(n)
                   flux_d(n)=0.5d0*sumdown(n)
                else                   ! replace integrals by refined j^th value
                   flux_u(n)=0.5d0*(flux_u(n)+sumup(n)/dble(it))
                   flux_d(n)=0.5d0*(flux_d(n)+sumdown(n)/dble(it))
                endif
             enddo

             if(j.gt.5)then            ! evaluate convergence in numerical integration
                maxerr=0.0d0
                do n=1,nz+1
                   if(old_flux_u(n).ne.0.0d0)then
                      err=dabs((flux_u(n)-old_flux_u(n))/old_flux_u(n))
                      if(err.gt.maxerr)maxerr=err
                   endif
                   if(old_flux_d(n).ne.0.0d0)then
                      err=dabs((flux_d(n)-old_flux_d(n))/old_flux_d(n))
                      if(err.gt.maxerr)maxerr=err
                   endif
                enddo
                if(maxerr.lt.eps)exit
             endif
             do n=1,nz+1
                old_flux_u(n)=flux_u(n)
                old_flux_d(n)=flux_d(n)
             enddo
c             if(j.eq.jmax)write(*,*)' W- Too many steps in qtrap'
          enddo

c      ...store upward and downward fluxes
          do n=1,nz+1
             fup_ir(i,n)=flux_u(n)
             fdo_ir(i,n)=flux_d(n)
          enddo
       enddo
c      ..End of loop on wavenumber grid

c      ...deallocate variables and close binary file with kabs data
       deallocate(tkl,tauv,tau,w0,g,gam1,gam2,lam,gam,
     # cp0,cm0,cpb,cmb,e1,e2,e3,e4,b0,b1,a,b,d,e,y,y1,y2,flux_u,flux_d,
     # gg,hh,jj,kk,alp1,alp2,sig1,sig2,iup,idown,sumup,sumdown,
     # old_flux_u,old_flux_d)

       return

       end
