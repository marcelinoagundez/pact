
       subroutine compute_correlated_k_distribution

       include 'pact.common'
       integer i,j,k,l,ip1,it1,ki,kj,ko,idko,ks,is
       real*8 xg(nmaxg),wg(nmaxg),kp1t1,kp1t2,kp2t1,kp2t2,t,u,
     # ktot(nmaxg),knew(nmaxg),kcomb(nmaxg*nmaxg),wcomb(nmaxg*nmaxg),
     # kcombmin,ksort(nmaxg*nmaxg),wsort(nmaxg*nmaxg),g1,g2,gs1,gs2,
     # transm,pwk,twk
       logical kused(nmaxg*nmaxg),skip(nmaxspec_ktab)

c      ...get Gaussian-Legendre quadrature points and weights
       call gauleg(0.0d0,0.99d0,xg(1:20),wg(1:20),20)      ! 20 points in the range g = 0    -0.99
       call gauleg(0.99d0,0.999d0,xg(21:30),wg(21:30),10)  ! 10 points in the range g = 0.99 -0.999
       call gauleg(0.999d0,1.0d0,xg(31:36),wg(31:36),6)    !  6 points in the range g = 0.999-1

c      ...skip species with a mole fraction below 1 ppb at all heights
       do i=1,nspec_ktab
          skip(i)=.true.
          do j=1,nz
             if(abun(idspec_ktab(i),j).ge.1.0d-9)then
                skip(i)=.false.
                exit
             endif
          enddo
       enddo

c      ...loop over atmosphere layers
       do j=1,nz
          pwk=pressure(j)
          if(pwk.lt.p_ktab(1))      pwk=p_ktab(1)
          if(pwk.gt.p_ktab(np_ktab))pwk=p_ktab(np_ktab)
          twk=tk(j)
          if(twk.lt.t_ktab(1))      twk=t_ktab(1)
          if(twk.gt.t_ktab(nt_ktab))twk=t_ktab(nt_ktab)

          i=1                                    ! identify neighboring pressures in k-tables
          do while(pwk.gt.p_ktab(i).and.i.lt.np_ktab)
             i=i+1
          enddo
          ip1=min(max(i-1,1),np_ktab-1)

          i=1                                    ! identify neighboring temperatures in k-tables
          do while(twk.gt.t_ktab(i).and.i.lt.nt_ktab)
             i=i+1
          enddo
          it1=min(max(i-1,1),nt_ktab-1)

c      ...loop over wavenumber grid
          do l=1,nwave

c      ...compute k(g) for species i
             is=0
             do i=1,nspec_ktab
                if(skip(i))cycle
                is=is+1
                do k=1,nmaxg                     ! get k(g) at neighboring p,T
                   kp1t1=kcoeff(i,ip1,it1,l,k)
                   kp1t2=kcoeff(i,ip1,it1+1,l,k)
                   kp2t1=kcoeff(i,ip1+1,it1,l,k)
                   kp2t2=kcoeff(i,ip1+1,it1+1,l,k)
c                                                ! bilinear interpolation
                   t=(dlog10(pwk)-dlog10(p_ktab(ip1)))/
     #                (dlog10(p_ktab(ip1+1))-dlog10(p_ktab(ip1)))
                   u=(twk-t_ktab(it1))/(t_ktab(it1+1)-t_ktab(it1))
                   knew(k)=(1.0d0-t)*(1.0d0-u)*kp1t1+
     #                             t*(1.0d0-u)*kp2t1+
     #                                     t*u*kp2t2+
     #                             (1.0d0-t)*u*kp1t2
c                                                ! convert from 'log k' to 'k * abun'
                   knew(k)=10.0d0**knew(k)*abun(idspec_ktab(i),j)
                enddo

                if(is.eq.1)then                  ! if first species, assign k to ktot and go for next species
                   ktot(:)=knew(:)
                   cycle
                else                             ! if cross section is zero in the wavenumber bin, skip combining k
                   if(knew(nmaxg).lt.sectmin)cycle
                endif

c      ...combine k(g) of different species
                do ki=1,nmaxg                    ! combine k(g) of species i with previous k(g)
                   do kj=1,nmaxg
                      k=(ki-1)*nmaxg+kj
                      kcomb(k)=(ktot(ki)+knew(kj))
                      wcomb(k)=wg(ki)*wg(kj)
                   enddo
                enddo

                kused(:)=.false.                 ! sort combined k(g)
                do ks=1,nmaxg*nmaxg
                   kcombmin=huge
                   idko=0
                   do ko=1,nmaxg*nmaxg
                      if(kused(ko))cycle
                      if(kcomb(ko).lt.kcombmin)then
                         kcombmin=kcomb(ko)
                         idko=ko
                      endif
                   enddo
                   if(idko.eq.0)stop' E- Error sorting combined k(g)'
                   ksort(ks)=kcomb(idko)
                   wsort(ks)=wcomb(idko)
                   kused(idko)=.true.
                enddo

                ktot(:)=0.0d0
                g2=0.0d0
                do k=1,nmaxg                     ! rebin to nominal g-grid
                   g1=g2
                   g2=g1+wg(k)
                   gs2=0.0d0
                   do ks=1,nmaxg*nmaxg
                      gs1=gs2
                      gs2=gs1+wsort(ks)
                      if(gs1.gt.g2)cycle
                      if(gs2.lt.g1)cycle
                      if(gs1.ge.g1.and.gs2.le.g2)then
                         ktot(k)=ktot(k)+ksort(ks)*wsort(ks)/wg(k)
                      elseif(gs1.lt.g1.and.gs2.le.g2)then
                         ktot(k)=ktot(k)+ksort(ks)*(gs2-g1)/wg(k)
                      elseif(gs1.ge.g1.and.gs2.gt.g2)then
                         ktot(k)=ktot(k)+ksort(ks)*(g2-gs1)/wg(k)
                      elseif(gs1.lt.g1.and.gs2.gt.g2)then
                         ktot(k)=ktot(k)+ksort(ks)
                      else
                         stop' E- Error rebinning k to nominal g-grid'
                      endif
                   enddo
                enddo
             enddo

c      ...compute transmission and add contribution to absorption coefficient [cm-1]
             transm=0.0d0
             do k=1,nmaxg
                transm=transm+dexp(-ktot(k)*density(j)*deltaz(j))*wg(k)
             enddo
             transm=max(transm,dexp(-700.0d0))   ! maximum optical depth set to 700 for numerical reasons
             kabs(l,j)=kabs(l,j)-dlog(transm)/deltaz(j)
          enddo
       enddo

       return

       end
