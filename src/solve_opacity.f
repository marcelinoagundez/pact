
c      Compute absorption and scattering coefficients (cm-1) for each spectral bin and atmosphere layer
       subroutine solve_opacity

       include 'pact.common'
       integer j,l

       write(*,1000)
       kabs(:,:)=0.0d0                           ! initialize to zero
       kscat(:,:)=0.0d0
       walb(:,:)=0.0d0

c.........test: grey atmosphere k = 0.01 cm-2 g-1
c       do l=1,nwave
c          do j=1,nz
c             kabs(l,j)=1.0d-2*apmass*density(j)
c          enddo
c       enddo
c       return
c.........

c      ...compute contribution of CIA (collision-induced absorption) to absorption coefficient
       call compute_cia

c      ...compute contribution of species to absorption coefficient: correlated k-distribution
       call compute_correlated_k_distribution

c      ...compute contribution of Rayleigh scattering to scattering coefficient
       call compute_rayleigh

c      ...compute single scattering albedo
       do l=1,nwave
          do j=1,nz
             walb(l,j)=kscat(l,j)/(kabs(l,j)+kscat(l,j))
             if((1.0d0-walb(l,j)).lt.1.0d-6)walb(l,j)=1.0d0-1.0d-6 ! avoid pure scattering, singular matrix in two-stream algorithm
          enddo
       enddo

       return

1000   format(2x,'I- Computing opacities')

       end
