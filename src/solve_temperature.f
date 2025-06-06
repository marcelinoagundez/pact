
       subroutine solve_temperature

       include 'pact.common'
       integer i,j,n,it,nup,nlo,nmaxit,nadapt
       real*8 eps,dnu,fint
       real*8, allocatable, dimension(:) :: dfnet,tk_old,tk_o,fpre,dtk,
     # tks

       write(*,1000)
       allocate(dfnet(nz),tk_old(nz),tk_o(nz),fpre(nz),dtk(nz),tks(nz))

       nmaxit=20000
c       nadapt=6
       nadapt=20
       eps=1.0d-7
       fint=sigmasb*tintplanet**4.0d0

       tk_old(1:nz)=tk(1:nz)
       call compute_heat_capacity      ! ...get heat capacity at each layer

c      ...time integration toward radiative equilibrium
       do it=1,nmaxit

c      ...solve radiative transfer
          call solve_radiative_transfer_stellar
          call solve_radiative_transfer_infrared

c      ...calculate upward, downward, and net bolometric fluxes [erg s-1 cm-2] at each atmosphere level
          fup(:)=0.0d0
          fdown(:)=0.0d0
          do i=1,nwave
             dnu=dwave(i)*clight
             do n=1,nz+1
                fup(n)=fup(n)+(fup_st(i,n)+fup_ir(i,n))*dnu
                fdown(n)=fdown(n)+(fdo_st(i,n)+fdo_ir(i,n))*dnu
             enddo
          enddo
          fnet(:)=0.0d0
          do n=1,nz
             fnet(n)=fup(n)-fdown(n)
          enddo
          fnet(nz+1)=fint

          do j=1,nz
             nup=j
             nlo=j+1
             dfnet(j)=fnet(nup)-fnet(nlo)

c      ...add smoothing temperature term to net fluxes (see Malik et al 2019, AJ, 157, 170)
c         switched off by default as it may cause the solution to deviate from radiative equilibrium
c             if(j.eq.1.or.j.eq.nz)then
c                fs=0.0d0
c             else
c                tmid=0.5d0*(tk(j-1)+tk(j+1))
c                fs=(tk(j)-tmid)**7.0d0
c             endif
c             dfnet(j)=dfnet(j)+fs

c      ...estimate time step and temperature change (see Malik et al 2017, AJ, 153, 56; modified as seen in HELIOS code)
             if(it.eq.1)fpre(j)=1.0d0
             dtk(j)=fpre(j)*pressure(j)/(plev(nup)-plev(nlo))*
     #          dfnet(j)/(dabs(dfnet(j))**0.9d0)
             if(dfnet(j).eq.0.0d0)dtk(j)=0.0d0

c      ...estimate time step and temperature change (in HELIOS code, identical results, to be tested)
c             dt(j)=1.0d3
c             if(dfnet(j).ne.0.0d0)
c     #          dt(j)=fpre(j)*pressure(nz)/(dabs(dfnet(j))**0.9d0)
c             dtk(j)=dfnet(j)/(plev(nz-1)-plev(nz))*dt(j)
c             if(dabs(dtk(j)).gt.5.0d2)
c     #          dtk(j)=5.0d2*(dfnet(j)/dabs(dfnet(j)))

c      ...estimate time step and temperature change (see Malik et al 2017, AJ, 153, 56)
c             trad=cplayer(j)*pressure(j)/(sigmasb*ggrav(j)*tk(j)**3.0d0)
c             fpre(j)=1.0d5/(dabs(dfnet(j))**0.9d0)
c             dt(j)=fpre(j)*trad
c             hrate(j)=ggrav(j)/cplayer(j)*
c             (dfnet(j)/(plev(nup)-plev(nlo)))
c             dtk(j)=hrate(j)*dt(j)

c      ...adaptative time step
             if(mod(it,nadapt).eq.0)then
                tk_o(j)=tk(j)
             elseif(mod(it,nadapt).eq.(nadapt-1))then
                if(dabs(tk(j)-tk_o(j)).lt.
     #          (dble(nadapt)/2.0d0*dabs(dtk(j))))then
                   fpre(j)=fpre(j)/1.5d0
                else
                   fpre(j)=fpre(j)*1.1d0
                endif
             endif
          enddo

c      ...update temperature
          do j=1,nz
             tk(j)=tk(j)+dtk(j)
          enddo

c      ...convergence criterion (see Malik et al 2017, AJ, 153, 56)
          conv=0.0d0
          do j=1,nz
             conv=max(conv,dabs(dfnet(j))/(sigmasb*tk(j)**4.0d0))
          enddo
c          write(*,1040,advance='no')it,conv,eps
          if(conv.lt.eps)exit
       enddo
       it=min(it,nmaxit)
       write(*,1030)it,conv,eps

c      ...if it did not converge, take previous iteration
       if(conv.ge.eps)then
          write(*,1010)
          if(time.gt.0.0d0.and.conv.gt.10.0d0*eps)then
             write(*,1020)
             tk(1:nz)=tk_old(1:nz)
             deallocate(dfnet,tk_old,tk_o,fpre,dtk,tks)
             return
          endif
       endif

c      ...smooth temperature to avoid discontinuities
       do j=2,nz-1
          tks(j)=0.25d0*tk(j-1)+0.5d0*tk(j)+0.25d0*tk(j+1)
       enddo
       do j=2,nz-1
          tk(j)=tks(j)
       enddo

c      ...calculate effective temperature of planet
       tplanet=(fup(1)/sigmasb)**0.25d0

c      ...do convective adjustment
       call convective_adjustment

       deallocate(dfnet,tk_old,tk_o,fpre,dtk,tks)

1000   format(2x,'I- Solving temperature')
1010   format(2x,'W- Maximum number of temperature iterations reached')
1020   format(2x,'W- Adopting temperature of previous time')
1030   format(2x,'I- Temperature iteration ',i5,
     # '    conv =',1pg12.4,'  -->  < ',1pe7.1)
c1040   format(70('\b'),2x,'I- Temperature iteration ',i5,
c     # '    conv =',1pg12.4,'  -->  < ',1pe7.1)

       return

       end




c_______________________________________________________________________

       subroutine compute_heat_capacity

       include 'pact.common'       

       real*8 r
       parameter(r=8.314472d7)                   ! Ideal gas constant           [erg mol-1 K-1]

       integer i,j,k
       real*8 t,a(9),c,norm,m

       do j=1,nz
          t=tk(j)
          cplayer(j)=0.0d0
          norm=0.0d0
          do i=1,nspec
             m=mspec(i)/amu                      ! [g mol-1]
             k=1
             do while(t.ge.temp_therm(i,k+1).and.k.lt.ntemp_therm(i))
                k=k+1
             enddo
             a(:)=atherm(i,:,k)

c      ...7-term NASA polynomials (see NASA RM-4513, p 9)
             if(type_therm(i).eq.7)then
                c=a(1)+a(2)*t+a(3)*t**2.0d0+a(4)*t**3.0d0+a(5)*t**4.0d0
c      ...9-term NASA polynomials (see NASA TP-2001-210959/REV1, p 2)
             elseif(type_therm(i).eq.9)then
                c=a(1)*t**(-2.0d0)+a(2)/t+a(3)+a(4)*t+a(5)*t**2.0d0+
     #          a(6)*t**3.0d0+a(7)*t**4.0d0
             else
                write(*,*)' E- Cannot compute Cp for: ',spec(i)
                stop
             endif

             cplayer(j)=cplayer(j)+c*r/m*(abun(i,j)*m)! Cp/R from NASA polynomials is unitless
             norm=norm+abun(i,j)*m                    !    we multiply by R and by the molecular mass [in g/mol] 
          enddo                                       !    to convert Cp to units of [erg g-1 K-1]
          cplayer(j)=cplayer(j)/norm                  ! To compute Cp for each atmospheric layer we weight
       enddo                                          !    Cp of each species by the mass per unit volumen of each species

       return
       end




c_______________________________________________________________________

       subroutine convective_adjustment

       include 'pact.common'
       integer j
       real*8 gtk,gadiab,cpmid

c      ...convective adjustment
       detached_convect=.false.
       convect(1)=.false.
       do j=2,nz
          gtk=tk(j)/tk(j-1)
          cpmid=0.5d0*(cplayer(j)+cplayer(j-1))
          gadiab=(pressure(j)/pressure(j-1))**(kboltz/apmass/cpmid)
          if(gtk.gt.gadiab)then
             tk(j)=tk(j-1)*gadiab
             convect(j)=.true.
             if(.not.convect(j-1))jadiab=j
          else
             convect(j)=.false.
             jadiab=nz
          endif
          if(.not.convect(j).and.convect(j-1))detached_convect=.true.
       enddo

       return
       end
