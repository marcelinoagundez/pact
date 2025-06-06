
c      From Fortran Numerical Recipes, Chap. 9.6
c      ...real variables changed to real*8
c      ...included call to mprove and save of fvec and fjac

       SUBROUTINE mnewt(ntrial,x,n,tolx,tolf)
       INTEGER n,ntrial,NP
       REAL*8 tolf,tolx,x(n)
       PARAMETER (NP=500)               ! Up to NP variables.
c      USES lubksb,ludcmp,usrfun
c      Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson steps to
c      improve the root. Stop if the root converges in either summed absolute variable increments
c      tolx or summed absolute function values tolf.
       INTEGER i,j,k,indx(NP)
       REAL*8 d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)
       real*8 fvec_save(NP),fjac_save(NP,NP)
       if(n.gt.NP)then
          write(*,*)' E- Too large size in mnewt'
          stop
       endif
       do 14 k=1,ntrial
          call usrfun(x,n,NP,fvec,fjac) ! User subroutine supplies function values at x in fvec
          do i=1,n                      ! *** save fvec and fjac
             fvec_save(i)=-fvec(i)
             do j=1,n
                fjac_save(i,j)=fjac(i,j)
             enddo
          enddo
          errf=0.                       !  and Jacobian matrix in fjac.
          do 11 i=1,n                   ! Check function convergence.
             errf=errf+abs(fvec(i))
11        enddo
          if(errf.le.tolf)return
          do 12 i=1,n                   ! Right-hand side of linear equations.
             p(i)=-fvec(i)
12        enddo
          call ludcmp(fjac,n,NP,indx,d) ! Solve linear equations using LU decomposition.
          call lubksb(fjac,n,NP,indx,p)
          call mprove(fjac_save,fjac,n,NP,indx,fvec_save,p)
          errx=0.                       ! Check root convergence.
          do 13 i=1,n                   ! Update solution.
             errx=errx+abs(p(i))
             x(i)=x(i)+p(i)
13        enddo
          if(errx.le.tolx)return
14     enddo
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 2.3
c      ...real variables changed to real*8

       SUBROUTINE ludcmp(a,n,np,indx,d)
       INTEGER n,np,indx(n),NMAX
       REAL*8 d,a(np,np),TINY
       PARAMETER (NMAX=500,TINY=1.0e-20)!Largest expected n, and a small number.
c      Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
c      the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
c      arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
c      row permutation effected by the partial pivoting; d is output as +-1 depending on whether
c      the number of row interchanges was even or odd, respectively. This routine is used in
c      combination with lubksb to solve linear equations or invert a matrix.
       INTEGER i,imax,j,k
       REAL*8 aamax,dum,sum,vv(NMAX)  ! vv stores the implicit scaling of each row.
       imax=0                         ! ad hoc initialization
       d=1.                         ! No row interchanges yet.
       do 12 i=1,n                  ! Loop over rows to get the implicit scaling
          aamax=0.                  ! information.
          do 11 j=1,n
             if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11        enddo
cc          if (aamax.eq.0.) pause 'singular matrix in ludcmp' ! No nonzero largest element.
          if (aamax.eq.0.) stop 'singular matrix in ludcmp'  ! No nonzero largest element.
          vv(i)=1./aamax                                     ! Save the scaling.
12     enddo
       do 19 j=1,n                  ! This is the loop over columns of Crout's method.
          do 14 i=1,j-1             ! This is equation (2.3.12) except for i = j.
             sum=a(i,j)
             do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13           enddo
             a(i,j)=sum
14        enddo
          aamax=0.                  ! Initialize for the search for largest pivot element.
          do 16 i=j,n               ! This is i = j of equation (2.3.12) and i = j+1: ::N
             sum=a(i,j)             ! of equation (2.3.13).
             do 15 k=1,j-1
                sum=sum-a(i,k)*a(k,j)
15           enddo
             a(i,j)=sum
             dum=vv(i)*abs(sum)     ! Figure of merit for the pivot.
             if (dum.ge.aamax) then ! Is it better than the best so far?
                imax=i
                aamax=dum
             endif
16        enddo
          if (j.ne.imax)then        ! Do we need to interchange rows?
             do 17 k=1,n            ! Yes, do so...
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
17           enddo
             d=-d                   ! ...and change the parity of d.
             vv(imax)=vv(j)         ! Also interchange the scale factor.
          endif
          indx(j)=imax
          if(a(j,j).eq.0.)a(j,j)=TINY
c      If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
c      For some applications on singular matrices, it is desirable to substitute TINY
c      for zero.
          if(j.ne.n)then            ! Now, Finally, divide by the pivot element.
             dum=1./a(j,j)
             do 18 i=j+1,n
                a(i,j)=a(i,j)*dum
18           enddo
          endif
19     enddo                        ! Go back for the next column in the reduction.
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 2.3
c      ...real variables changed to real*8

       SUBROUTINE lubksb(a,n,np,indx,b)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(n)
c      Solves the set of n linear equations A . X = B. Here a is input, not as the matrix A but
c      rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
c      permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
c      and returns with the solution vector X. a, n, np, and indx are not modified by this routine
c      and can be left in place for successive calls with different right-hand sides b. This routine
c      takes into account the possibility that b will begin with many zero elements, so it is efficient
c      for use in matrix inversion.
       INTEGER i,ii,j,ll
       REAL*8 sum
       ii=0                  ! When ii is set to a positive value, it will become the index
       do 12 i=1,n           ! of the first nonvanishing element of b. We now do
          ll=indx(i)         ! the forward substitution, equation (2.3.6). The only new
          sum=b(ll)          ! wrinkle is to unscramble the permutation as we go.
          b(ll)=b(i)
          if (ii.ne.0)then
             do 11 j=ii,i-1
                sum=sum-a(i,j)*b(j)
11           enddo
          else if (sum.ne.0.) then
             ii=i            ! A nonzero element was encountered, so from now on we will
          endif              ! have to do the sums in the loop above.
          b(i)=sum
12     enddo
       do 14 i=n,1,-1        ! Now we do the backsubstitution, equation (2.3.7).
          sum=b(i)
          do 13 j=i+1,n
             sum=sum-a(i,j)*b(j)
13        enddo
          b(i)=sum/a(i,i)    ! Store a component of the solution vector X.
14     enddo
       return                ! All done!
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 2.5
c      ...real variables changed to real*8

       SUBROUTINE mprove(a,alud,n,np,indx,b,x)
       INTEGER n,np,indx(n),NMAX
       REAL*8 a(np,np),alud(np,np),b(n),x(n)
       PARAMETER (NMAX=500)          ! Maximum anticipated value of n.
c      USES lubksb
c      Improves a solution vector x(1:n) of the linear set of equations A . X = B. The matrix
c      a(1:n,1:n), and the vectors b(1:n) and x(1:n) are input, as is the dimension n. Also
c      input is alud, the LU decomposition of a as returned by ludcmp, and the vector indx also
c      returned by that routine. On output, only x(1:n) is modified, to an improved set of values.
       INTEGER i,j
       REAL*8 r(NMAX)
       DOUBLE PRECISION sdp
       do 12 i=1,n                   ! Calculate the right-hand side, accumulating the residual
          sdp=-b(i)                  ! in double precision.
          do 11 j=1,n
             sdp=sdp+dble(a(i,j))*dble(x(j))
11        enddo
          r(i)=sdp
12     enddo
       call lubksb(alud,n,np,indx,r) ! Solve for the error term,
       do 13 i=1,n                   ! and subtract it from the old solution.
          x(i)=x(i)-r(i)
13     enddo
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 2.4
c      ...real variables changed to real*8

       SUBROUTINE tridag(a,b,c,r,u,n)
       INTEGER n,NMAX
       REAL*8 a(n),b(n),c(n),r(n),u(n)
       PARAMETER (NMAX=500)
c      Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1)
c      a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified
c      Parameter: NMAX is the maximum expected value of n
       INTEGER j
       REAL*8 bet,gam(NMAX)     ! One vector of workspace, gam, is needed
       if (b(1).eq.0.) stop 'tridag: rewrite equations'
c      If this happens then you should rewrite your equations as a set of order N-1,
c      with u_2 trivially eliminated
       bet=b(1)
       u(1)=r(1)/bet
       do 11 j=2,n                               ! Decomposition and forward substitution
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if (bet.eq.0.) stop 'tridag failed'   ! Algorithm fails; see Chap. 2.4
          u(j)=(r(j)-a(j)*u(j-1))/bet
11     enddo
       do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
12     enddo
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 5.8
c      ...all real variables in double precision (real*8)

       SUBROUTINE chebft(a,b,c,n,func)
       INTEGER n,NMAX
       REAL*8 a,b,c(n),func,pi
       EXTERNAL func
       PARAMETER (NMAX=50,pi=3.141592653589793d0)        
c      Chebyshev fit: Given a function func, lower and upper limits of the interval [a,b], and
c      a maximum degree n, this routine computes the n coefficients c_k such that
c      func(x) ~ [Sum_{k=1}^{n} c_k T_{k-1}(y)] - c1/2, where y and x are related by (5.8.10):
c         y = (x-0.5(b+a)) / (0.5(b-a))
c      This routine is to be used with moderately large n (e.g., 30 or 50), the array of c's subsequently
c      to be truncated at the smaller value m such that c_{m+1} and subsequent elements are negligible.
       INTEGER j,k
       REAL*8 bma,bpa,fac,y,f(NMAX),sum
       bma=0.5d0*(b-a)
       bpa=0.5d0*(b+a)
       do 11 k=1,n                      ! We evaluate the function at the n points required by (5.8.7)
          y=dcos(pi*(k-0.5d0)/n)
          f(k)=func(y*bma+bpa)
11     enddo
       fac=2.0d0/n
       do 13 j=1,n
          sum=0.0d0                     ! We will accumulate the sum in double precision
          do 12 k=1,n
             sum=sum+f(k)*dcos((pi*(j-1))*((k-0.5d0)/n))
12        enddo
          c(j)=fac*sum
13     enddo
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 5.8
c      ...all real variables in double precision (real*8)

       REAL*8 FUNCTION chebev(a,b,c,m,x)
       INTEGER m
       REAL*8 a,b,x,c(m)
c      Chebyshev evaluation: All arguments are input. c(1:m) is an array of Chebyshev coefficients,
c      the first m elements of c output from chebft (which must have been called with the same
c      a and b). The Chebyshev polynomial [Sum_{k=1}^{m} c_k T_{k-1}(y)] - c1/2 is evaluated at a
c      point y = [x-(b+a)/2]/[(b-a)/2], and the result is returned as the function value.
       INTEGER j
       REAL*8 d,dd,sv,y,y2
       if((x-a)*(x-b).gt.0.0d0) stop 'x not in range in chebev'
       d=0.0d0
       dd=0.0d0
       y=(2.0d0*x-a-b)/(b-a)        ! Change of variable
       y2=2.0d0*y
       do 11 j=m,2,-1               ! Clenshaw's recurrence
          sv=d
          d=y2*d-dd+c(j)
          dd=sv
11     enddo
       chebev=y*d-dd+0.5d0*c(1)     ! Last step is different
       return
       END




c_______________________________________________________________________

c      From Fortran Numerical Recipes, Chap. 4.5
c      ...all real variables in double precision (real*8)
c      ...PI number added as parameter

       SUBROUTINE gauleg(x1,x2,x,w,n)
       INTEGER n
       REAL*8 x1,x2,x(n),w(n)
       REAL*8 EPS
       PARAMETER (EPS=3.d-14)             ! EPS is the relative precision
       REAL*8 PI
       PARAMETER(PI=3.141592653589793d0)
c      Given the lower and upper limits of integration x1 and x2, an given n, this routine
c      returns arrays x(1:n) and w(1:n) of length n, containing the abscissas and weights
c      of the Gauss-Legendre n-point quadrature formula
       INTEGER i,j,m
       REAL*8 p1,p2,p3,pp,xl,xm,z,z1
       m=(n+1)/2                          ! The roots are symmetric in the interval, so we
       xm=0.5d0*(x2+x1)                   !    only have to find half of them
       xl=0.5d0*(x2-x1)
       do 12 i=1,m                        ! Loop over the desired roots
          z=dcos(PI*(i-.25d0)/(n+.5d0))
c      Starting with the above approximation to the ith root, we enter the main loop of
c      refinement by Newton's method
1         continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        enddo
c      p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by
c      a standard relation involving also p2, the polynomial of one lower order
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp                      ! Newton's method
          if(dabs(z-z1).gt.EPS)goto 1
          x(i)=xm-xl*z                    ! Scale the root to the desired interval,
          x(n+1-i)=xm+xl*z                !    and put in its symmetric counterpart
          w(i)=2.d0*xl/((1.d0-z*z)*pp*pp) ! compute the weight
          w(n+1-i)=w(i)                   !    and its symmetric counterpart
12     enddo
       return
       END
