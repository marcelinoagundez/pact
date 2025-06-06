
       subroutine compute_rayleigh

       use pact_data
       include 'pact.common'
       integer i,j,k,l
       real*8 sigma
       real*8, allocatable, dimension(:)  :: polz
       logical missing_polarizabilities
       logical, allocatable, dimension(:) :: misspolz

c      ...evaluate polarizabilities
       allocate (polz(nspec),misspolz(nspec))
       do i=1,nspec
          polz(i)=0.0d0
          do k=1,nmaxpo
             if(trim(spec(i)).ne.trim(po_spec(k)))cycle
             polz(i)=polary(k)
             exit
          enddo
       enddo

c      ...is there any abundant (> 1e-4) species with no polarizability data?
       missing_polarizabilities=.false.
       misspolz(:)=.false.
       do i=1,nspec
          do j=1,nz
             if(abun(i,j).gt.1.0d-4.and.polz(i).eq.0.0d0)
     #       misspolz(i)=.true.
          enddo
          if(misspolz(i))missing_polarizabilities=.true.
       enddo
       if(missing_polarizabilities)then
          write(*,1000)
          do i=1,nspec
             if(.not.misspolz(i))cycle
             write(*,*)trim(spec(i))
          enddo
       endif

c      ...compute Rayleigh scattering coefficient
       do i=1,nspec
          if(polz(i).eq.0.0d0)cycle
          do l=1,nwave
             sigma=128.0d0/3.0d0*pi**5.0d0*
     #       (polz(i)*1.0d-24)**2.0d0*wave(l)**4.0d0
             do j=1,nz
                kscat(l,j)=kscat(l,j)+sigma*abun(i,j)*density(j)
c                gscat(l,j)=0.0d0    ! gscat is initialized to zero; to be changed if clouds are included
             enddo
          enddo
       enddo
       deallocate (polz,misspolz)

       return

1000   format(1x,' W- Recommended to include polarizability for:')

       end
