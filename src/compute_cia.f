
       subroutine compute_cia

       include 'pact.common'
       integer nmax
       parameter(nmax=30000)
       integer i,j,k,l,id1,id2,nt,nwav,k1,k2
       real*8 t,valt1,valt2,add,x1,x2,y1,y2,slope
       real*8, allocatable, dimension(:)   :: temp,wav,valt
       real*8, allocatable, dimension(:,:) :: val
       character*300 datafile
       character*5000 line,o_line
       logical read_temp

c      ...loop over number of CIA data sets
       do i=1,ncia
          id1=idspec1_cia(i)
          id2=idspec2_cia(i)

          nt=ntemp_ciafile(i)
          allocate(temp(nt),wav(nmax),val(nmax,nt),valt(nmax))
          datafile(1:300)=''
          datafile=trim(data_path)//trim(ciafile(i))
          open(unit=2,file=datafile,status='old')
          read_temp=.false.

c      ...read line by line: temperatures and CIA data
          j=0
          do
             o_line(1:5000)=line(1:5000)
             line(1:5000)=''
             read(2,'(a5000)',err=100,end=10)line
             if(line(1:1).eq.'!')cycle
             if(.not.read_temp)then              ! read temperatures
                read(o_line(2:5000),*,err=110)(temp(k),k=1,nt)
                read_temp=.true.
             else
                j=j+1
                if(j.gt.nmax)        goto 120
                read(line(1:5000),*,err=130)wav(j),(val(j,k),k=1,nt)
             endif
          enddo
10        close(2)
          nwav=j
          do k=2,nt
             if(temp(k).le.temp(k-1))goto 140
          enddo
          do j=2,nwav
             if(wav(j).le.wav(j-1))  goto 150
          enddo

c      ...extract k absorption coefficient at every height and spectral interval
          do j=1,nz
             t=tk(j)
             if(t.lt.temp(1)) t=temp(1)
             if(t.gt.temp(nt))t=temp(nt)

             if(nt.eq.1)then
                do l=1,nwav
                   valt(l)=dlog10(val(l,1))
                enddo
             else
                k=1                              ! interpolate logarithm(CIA) to current tk
                do while(t.ge.temp(k).and.k.lt.nt)
                   k=k+1
                enddo
                k=min(max(k-1,1),nt-1)
                do l=1,nwav
                   slope=(dlog10(val(l,k+1))-dlog10(val(l,k)))/
     #             (temp(k+1)-temp(k))
                   valt(l)=dlog10(val(l,k))+slope*(t-temp(k))
                enddo
             endif

             do l=1,nwave
                k=1                              ! interpolate logarithm(CIA) to lower edge of wave(l) bin
                do while(wave_l(l).ge.wav(k).and.k.lt.nwav)
                   k=k+1
                enddo
                k1=min(max(k-1,1),nwav-1)
                slope=(valt(k1+1)-valt(k1))/
     #             (dlog10(wav(k1+1))-dlog10(wav(k1)))
                valt1=valt(k1)+
     #             slope*(dlog10(wave_l(l))-dlog10(wav(k1)))
                k=k1                             ! interpolate logarithm(CIA) to upper edge of wave(l) bin
                do while(wave_u(l).ge.wav(k).and.k.lt.nwav)
                   k=k+1
                enddo
                k2=min(max(k-1,1),nwav-1)
                slope=(valt(k2+1)-valt(k2))/
     #             (dlog10(wav(k2+1))-dlog10(wav(k2)))
                valt2=valt(k2)+
     #             slope*(dlog10(wave_u(l))-dlog10(wav(k2)))

                add=0.0d0                        ! compute k absorption contribution
                do k=k1,k2                       ! to current spectral interval and height
                   x1=max(wav(k),  wave_l(l))
                   x2=min(wav(k+1),wave_u(l))
                   y1=valt(k)
                   y2=valt(k+1)
                   if(k.eq.k1)y1=valt1
                   if(k.eq.k2)y2=valt2
                   add=add+0.5d0*(y1+y2)*(x2-x1)
                enddo
                add=10.0d0**(add/dwave(l))

                if(wave_l(l).lt.wav(1))   add=10.0d0**valt(1)    ! wavenumber outside
                if(wave_u(l).gt.wav(nwav))add=10.0d0**valt(nwav) ! tabulated range

                add=add*abun(id1,j)*abun(id2,j)* ! convert from cm-1 amagat-2 to cm-1
     #          (density(j)/amagat)**2.0d0
                kabs(l,j)=kabs(l,j)+add          ! add contribution to absorption coefficient
             enddo
          enddo

          deallocate(temp,wav,val,valt)
       enddo

       return

100    write(*,*)' E- Error reading CIA file:'
       write(*,*)datafile
       stop
110    write(*,*)' E- Error reading temperatures in CIA file:'
       write(*,*)datafile
       stop
120    write(*,*)' E- Too many wavenumber points in CIA file:'
       write(*,*)datafile
       stop
130    write(*,*)' E- Error reading CIA data in file:'
       write(*,*)datafile
       stop
140    write(*,*)' E- Bad temperature ordering in CIA file:'
       write(*,*)datafile
       stop
150    write(*,*)' E- Bad wavenumber ordering in CIA file:'
       write(*,*)datafile
       stop

       end
