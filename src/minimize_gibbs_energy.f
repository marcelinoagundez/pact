
       subroutine minimize_gibbs_energy
       
       include 'pact.common'       
       include 'chemical_equilibrium.common'       
       integer i,j,iwk,it,nmaxit,nminit,ition,nmaxition
       real*8 x(nmaxspec),dlnng(nmaxspec),rwk,lambda,lambda1,lambda2,
     # err,maxerr,dpiion,o_pilag(nmaxeqt-1)

c      ...assign initial guesses for pi Lagrange multipliers, Delta(n_condensed), Delta(n)
       nmaxit=100                    ! maximum number of Newton-Raphson iterations
       nminit=5                      ! minimum number of Newton-Raphson iterations
       nmaxition=100                 ! maximum number of ion iterations
       neqt=nelem                    ! element conservation equations
       if(ion)neqt=neqt+1            ! charge conservation equation
       npilag=neqt
       neqt=neqt+ncondensed          ! condensed species equations
       neqt=neqt+1                   ! total mole number conservation equation
       do i=1,npilag
          pilag(i)=0.0d0
       enddo
       do i=1,ncondensed
          dnc(i)=0.0d0
       enddo
       dlnn=0.0d0

c      ...assign initial values to x array
       it=1
       converge_gibbs=.false.
       do i=1,npilag
          x(i)=pilag(i)
       enddo
       do i=1,ncondensed
          x(npilag+i)=dnc(i)
       enddo
       x(neqt)=dlnn

       do while (it.le.nmaxit.and..not.converge_gibbs)

c      ...compute chemical potential of all species
          call compute_chemical_potential(pce,tce)

c      ...call Newton-Raphson routine, 1 iteration and back
          call mnewt(1,x,neqt,tiny,tiny)

c      ...update values of pi lagrange multipliers and Delta Ln(n)
          do i=1,npilag
             pilag(i)=x(i)
          enddo
          do i=1,ncondensed
             dnc(i)=x(npilag+i)
          enddo
          dlnn=x(neqt)

c      ...compute correction to abundances: Delta Ln[n_gas(j)]
          do j=1,ngas
             dlnng(j)=0.0d0
             do i=1,npilag
                dlnng(j)=dlnng(j)+nat(j,i)*pilag(i)
             enddo
             dlnng(j)=dlnng(j)+dlnn-mupot(j)
          enddo

c      ...limit corrections to n, n_j(gas)
          if(dabs(dlnn).gt.0.4d0)dlnn=0.4d0*dlnn/dabs(dlnn)
          do j=1,ngas
             if(dabs(dlnng(j)).gt.2.0d0)
     #          dlnng(j)=2.0d0*(dlnng(j)/dabs(dlnng(j)))
          enddo

c      ...compute control factor: lambda
          lambda1=dabs(5.0d0*dlnn)
          lambda2=huge
          do j=1,ngas
             if((ncea(j)/nceatot).gt.1.0d-8)then
                if(dabs(dlnng(j)).gt.lambda1)lambda1=dabs(dlnng(j))
             else
                if(dlnng(j).ge.0.0d0)then
                   rwk=dabs((-dlog(ncea(j)/nceatot-dlog(1.0d4))/
     #             (dlnng(j)-dlnn)))
                   if(rwk.lt.lambda2)lambda2=rwk
                endif
             endif
          enddo
          lambda1=2.0d0/lambda1
          lambda=dmin1(1.0d0,lambda1,lambda2)

c      ...update estimates of n and n(j).
          do j=1,nspec
             if(.not.condensed(j))
     #          ncea(j)=ncea(j)*dexp(lambda*dlnng(j))
             if(condensed(j))ncea(j)=ncea(j)+lambda*dnc(j-ngas)
          enddo
          nceatot=nceatot*dexp(lambda*dlnn)! apply correction to n

c      ...have we converged ?
          iwk=0                                  ! convergence criteria on abundance variations
          rwk=0.0d0
          maxerr=0.0d0
          do j=1,nspec
             rwk=rwk+ncea(j)
          enddo
          do j=1,nspec+1
             if(j.le.nspec)then
                if(.not.condensed(j))err=ncea(j)*dabs(dlnng(j))/rwk
                if(condensed(j))err=dabs(dnc(j-ngas))/rwk
             else
                err=nceatot*dabs(dlnn)/rwk
             endif
             if(err.le.0.5d-5)iwk=iwk+1
             if(err.gt.maxerr)maxerr=err
          enddo
          if(iwk.eq.nspec+1)converge_gibbs=.true.
          do i=1,nelem                           ! convergence criteria on element conservation
             if(b0(i).le.1.0d-6)cycle
             rwk=0.0d0
             do j=1,nspec
                rwk=rwk+nat(j,i)*ncea(j)
             enddo
             rwk=dabs(b0(i)-rwk)/b0max
             if(rwk.gt.1.0d-6)converge_gibbs=.false.
          enddo
          if(b0min.lt.1.0d-6.and.it.gt.1)then    ! convergence criteria if trace elements
             do i=1,nelem
                rwk=dabs((o_pilag(i)-pilag(i))/pilag(i))
                if(rwk.gt.1.0d-3)converge_gibbs=.false.
             enddo
          endif
          if(it.lt.nminit)converge_gibbs=.false.
          o_pilag(1:npilag)=pilag(1:npilag)
          it=it+1
       enddo
       it=it-1
c       if(.not.converge_gibbs)stop' E - Convergence not reached'
       if(.not.converge_gibbs)return

c      ...verify convergence for ions
       if(ion)then
          dpiion=huge
          ition=1
          do while(dabs(dpiion).gt.1.0d-4.and.ition.le.nmaxition)
             dpiion=0.0d0
             rwk=0.0d0
             do j=1,ngas
                dpiion=dpiion+nat(j,npilag)*ncea(j)
                rwk=rwk+nat(j,npilag)*nat(j,npilag)*ncea(j)
             enddo
             dpiion=-dpiion/rwk
             do j=1,ngas
                if(charge(j).eq.0)cycle
                ncea(j)=ncea(j)*dexp(nat(j,npilag)*dpiion)
             enddo
             ition=ition+1
          enddo
          ition=ition-1
          if(ition.ge.nmaxition)stop' E - Ion-convergence not reached'
       endif

       return

       end




c_______________________________________________________________________

       subroutine compute_chemical_potential(p,t)

       include 'pact.common'       
       include 'chemical_equilibrium.common'       
       integer i,j,k
       real*8 p,t,a(9),h,s,g

       tflag(1:nspec)=.false.
       do i=1,nspec
          k=1
          do while(t.ge.temp_therm(i,k+1).and.k.lt.ntemp_therm(i))
             k=k+1
          enddo
          do j=1,natherm(i)
             a(j)=atherm(i,j,k)
          enddo

c      ...7-term polynomials fitted to G^o/RT = H^o/RT - S^o/R (standard Gibbs energy)
          if(type_therm(i).eq.1)then
             g=0.0d0
             do j=1,7
                g=g+a(j)*t**(j-3)
             enddo

c      ...7-term NASA polynomials (see NASA RM-4513, p 9)
          elseif(type_therm(i).eq.7)then
             h=a(1)+a(2)*t/2.0d0+a(3)*t**2.0d0/3.0d0+
     #       a(4)*t**3.0d0/4.0d0+a(5)*t**4.0d0/5.0d0+a(6)/t
             s=a(1)*dlog(t)+a(2)*t+a(3)*t**2.0d0/2.0d0+
     #       a(4)*t**3.0d0/3.0d0+a(5)*t**4.0d0/4.0d0+a(7)
             g=h-s

c      ...9-term NASA polynomials (see NASA TP-2001-210959/REV1, p 2)
          elseif(type_therm(i).eq.9)then
             h=-a(1)*t**(-2.0d0)+a(2)*dlog(t)/t+a(3)+
     #       a(4)*t/2.0d0+a(5)*t**2.0d0/3.0d0+
     #       a(6)*t**3.0d0/4.0d0+a(7)*t**4.0d0/5.0d0+a(8)/t
             s=-a(1)*t**(-2.0d0)/2.0d0-a(2)/t+
     #       a(3)*dlog(t)+a(4)*t+a(5)*t**2.0d0/2.0d0+
     #       a(6)*t**3.0d0/3.0d0+a(7)*t**4.0d0/4.0d0+a(9)
             g=h-s

c      ...unknown data type
          else
             write(*,*)' E- unknown therm data type for: ',spec(i)
             stop
          endif

          if(t.lt.temp_therm(i,1).or.
     #       t.gt.temp_therm(i,ntemp_therm(i)+1))
     #       tflag(i)=.true.
          if(condensed(i))then
             mupot(i)=g
          else
             mupot(i)=g+dlog(ncea(i)/nceatot)+dlog(p)
          endif
       enddo

       return
       end
