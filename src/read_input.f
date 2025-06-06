
       subroutine read_input
       
       include 'pact.common'       
       integer i,j,k,l

       write(*,1000)inpfile
       i=nmaxcharfile
       do while (inpfile(i:i).ne.'.')
          i=i-1
       enddo
       ncharmodel=i-1
       namemodel(1:ncharmodel)=inpfile(1:ncharmodel)

c      ...read input file
       open(unit=2,err=100,file=inpfile,status='old')
       read(2,*)rplanet,pplanet
       rplanet=rplanet*rearth                    ! R(Earth) -> cm
       pplanet=pplanet*1.0d6                     ! bar      -> dyn cm-2
       read(2,*)mplanet
       mplanet=mplanet*mearth                    ! M(Earth) -> g
       read(2,*)tintplanet
       read(2,*)albedoplanet,emplanet
       if(albedoplanet.lt.0.0d0.or.albedoplanet.gt.1.0d0)
     # stop ' E- Invalid planet surface albedo'
       read(2,*)fhrd
       read(2,*)zenithangle
       zenithangle=zenithangle*2.0d0*pi/360.0d0  ! deg      -> rad
       read(2,*)rstar
       rstar=rstar*rsun                          ! R(Sun)   -> cm
       read(2,*)dist
       dist=dist*au                              ! AU       -> cm
       read(2,*)starfile,tstar
       read(2,*)pmax,pmin,nz
       pmin=pmin*1.0d6                           ! bar      -> dyn cm-2
       pmax=pmax*1.0d6                           ! bar      -> dyn cm-2
       if(nz.gt.nmaxz)stop' E- Too many heights requested'
       read(2,*)specfile,iabun
       if(iabun.ne.0.and.iabun.ne.1)stop' E- Bad value for iabun'
       read(2,*)thermfile
       read(2,*)reacfile
       read(2,*)eddyfile
       read(2,*)i,j,k,l
       if((i.ne.0.and.i.ne.1).or.
     #    (j.ne.0.and.j.ne.1).or.
     #    (k.ne.0.and.k.ne.1).or.
     #    (l.ne.0.and.l.ne.1))then
          write(*,*)' E- Compute initial temperature'
          write(*,*)'    Self-consistent temperature'
          write(*,*)'    Diffusion'
          write(*,*)'    Photochemistry'
          write(*,*)' must be set to 0/1'
          stop
       endif
       if(i.eq.1)solvitemp=.true.
       if(j.eq.1)self_temp=.true.
       if(k.eq.1)diffusion=.true.
       if(l.eq.1)photochem=.true.
       if(sbrotation)then                        ! solid body rotation activated; to be implemented with temperature computation
          read(2,*)lptfile                       ! longitude-pressure-temperature file
          read(2,*)nlong,nperiod                 ! No longitudes per period, No rotation periods
          nlong=nlong+1                          ! to account for the first substellar longitude
          if(nlong.gt.nmaxlong)stop' E- Many longitudes per period'
       endif
       read(2,*)initfile
       if(solvitemp.and.trim(initfile).ne.'None')
     # stop' E- Initial (z,p,T,composition) file must be set to None'
       if(.not.solvitemp.and.trim(initfile).eq.'None')
     # stop' E- Provide name of initial (z,p,T,composition) file'
       read(2,*)data_path
       read(2,*)ktabledir
       read(2,*)chemdir
       read(2,*)photodir
       read(2,*)stardir
       close(2)

       tplanet_th=                               ! evaluate theoretical equilibrium temperature of planet
     # ((1.0d0-albedoplanet)/emplanet)**0.25d0*
     # dsqrt(rstar/2.0d0/dist)*tstar

       return

100    write(*,*)' E- Error opening model input file: ',inpfile
       stop

1000   format(2x,'I- Running model ',a)

       end




c_______________________________________________________________________

       subroutine get_pressure_grid
c      build a pressure grid, equispatial in log10(p), from top to bottom of atmosphere

       include 'pact.common'       
       integer j,j1,j2
       real*8 dp

c      ...pressure at midpoint of each layer
       dp=(dlog10(pmax)-dlog10(pmin))/dble(nz-1)
       do j=1,nz
          pressure(j)=10.0d0**(dlog10(pmin)+dp*(j-1))
          if(j.eq.1) pressure(j)=pmin
          if(j.eq.nz)pressure(j)=pmax
       enddo

c      ...pressure at each level (interfaces between layers)
       do j=1,nz+1
          j1=min(max(j-1,1),nz-1)
          j2=j1+1
          dp=0.5d0*(dlog10(pressure(j1))-dlog10(pressure(j2)))
          plev(j)=pressure(j)*10.0d0**dp
          if(j.eq.nz+1)plev(j)=pressure(nz)/10.0d0**dp
       enddo

       return
       end




c_______________________________________________________________________

       subroutine get_spectral_grid
c      build a spectral grid between wavemin and wavemax with the requested spectral resolution
c      the working variable is always the wavenumber in cm-1

       include 'pact.common'
       integer i,j
       real*8 lambda
       character*(nmaxcharstring) datafile

       datafile(1:nmaxcharstring)=''
       datafile=trim(data_path)//trim(ktabledir)//'ktable_spgrid.dat'
       open(unit=2,file=datafile,status='old')
       read(2,*)
       do i=1,nwave
          read(2,*)j,lambda,wave(i),wave_l(i),wave_u(i)
       enddo
       close(2)
       do i=1,nwave
          dwave(i)=wave_u(i)-wave_l(i)
       enddo
       write(*,1000)nwave,1.0d4/wavemax,1.0d4/wavemin,respow

       return

1000   format(
     # 2x,'I- ',i7,' channels (',0pf6.2,'-',0pf6.1,' um),',
     # ' R(lambda/Dlambda) =',0pf8.3)

       end




c_______________________________________________________________________

       subroutine get_star_spectrum
c      ...read stellar spectrum and resample to spectral grid

       include 'pact.common'       
       integer i,j,jl,ju,nxy
       real*8, allocatable, dimension(:)   :: xr,yr,x,y
       real*8 slope,yl,yu,x1,x2,y1,y2,fbol,nu
       real*8 bnuplanck
       character*256 line

c      ...use blackbody for stellar spectrum
       if(trim(starfile).eq.'None')then
          do i=1,nwave
             nu=wave(i)*clight
             fstar(i)=bnuplanck(tstar,nu)*       ! convert intensity [erg s-1 cm-2 Hz-1 sr-1]
     #       pi*(rstar/dist)**2.0d0              ! to flux at planet [erg s-1 cm-2 Hz-1]
          enddo
          return
       endif

c      ...read stellar spectrum from starfile
       allocate(xr(nmaxxy),yr(nmaxxy),x(nmaxxy),y(nmaxxy))
       write(*,1000)starfile
       open(unit=2,file=trim(data_path)//trim(stardir)//trim(starfile),
     # status='old')
       i=0
       do
          line(1:256)=''
          read(2,'(a256)',err=100,end=10)line
          if(line(1:1).eq.'!')cycle

          i=i+1
          if(i.gt.nmaxxy)stop ' E- Too many points in star spectrum'
          read(line(1:256),*,err=100)xr(i),yr(i)           ! wavelength [um], intensity [erg s-1 cm-2 Hz-1 sr-1]
          xr(i)=1.0d4/xr(i)                                ! wavenumber     [=] cm-1
          yr(i)=yr(i)*pi*(rstar/dist)**2.0d0               ! flux at planet [=] erg s-1 cm-2 Hz-1
       enddo
10     close(2)
       nxy=i
       do i=1,nxy
          x(i)=xr(nxy-i+1)
          y(i)=yr(nxy-i+1)
       enddo
       if(x(1).gt.wavemin.or.x(nxy).lt.wavemax)            ! missing star data in some wavenumber bin
     # stop ' E- Stellar spectrum not wide enough'

c      ...compute star effective temperature
       fbol=0.0d0
       do i=1,nxy-1
          fbol=fbol+0.5d0*(y(i+1)+y(i))*(dist/rstar)**2.0d0*
     #    (x(i+1)-x(i))*clight
       enddo
       tstar=(fbol/sigmasb)**0.25d0

c      ...resample stellar intensity to the wavenumber grid
       do i=1,nwave
          j=1
          do while(x(j).lt.wave_l(i).and.j.lt.nxy)
             j=j+1
          enddo
          j=max(j-1,1)
          jl=j
          do while(x(j).lt.wave_u(i).and.j.lt.nxy)
             j=j+1
          enddo
          ju=min(max(j-1,1),nxy-1)
          slope=(y(jl+1)-y(jl))/(x(jl+1)-x(jl))
          yl=y(jl)+slope*(wave_l(i)-x(jl))
          slope=(y(ju+1)-y(ju))/(x(ju+1)-x(ju))
          yu=y(ju)+slope*(wave_u(i)-x(ju))

          fstar(i)=0.0d0
          if(x(jl).ge.wave_u(i).or.x(ju+1).le.wave_l(i))cycle ! no star data in spectral bin
          do j=jl,ju                                          ! should not happen
             x1=x(j)
             y1=y(j)
             x2=x(j+1)
             y2=y(j+1)
             if(j.eq.jl.and.x(jl).lt.wave_l(i))then
                x1=wave_l(i)
                y1=yl
             endif
             if(j.eq.ju.and.x(ju+1).gt.wave_u(i))then
                x2=wave_u(i)
                y2=yu
             endif
             fstar(i)=fstar(i)+0.5d0*(y1+y2)*(x2-x1)
          enddo
          fstar(i)=fstar(i)/dwave(i)
       enddo

       deallocate(xr,yr,x,y)

       return

100    write(*,*)' E- Error reading stellar spectrum file: ',starfile
       stop

1000   format(1x,' I- Reading stellar spectrum from ',a45)

       end




c_______________________________________________________________________

       subroutine adjust_vertical_structure
c      compute vertical z structure to fullfil hydrostatic equilibrium

       include 'pact.common'
       integer j,j0
       real*8 tk0,p0,o_z,o_p,tkmid,g,h,slope

c      ...compute density of particles
       do j=1,nz
          density(j)=pressure(j)/kboltz/tk(j)
       enddo

c      ...locate reference level (at which pressure equals reference pressure)
       p0=pplanet
       j=1
       do while(pressure(j).lt.p0)
          j=j+1
       enddo
       j0=j-1
       if(j0.lt.1) j0=1
       if(j0.ge.nz)j0=nz-1
       slope=(tk(j0+1)-tk(j0))/
     # (dlog10(pressure(j0+1))-dlog10(pressure(j0)))
       tk0=tk(j0)+slope*(dlog10(p0)-dlog10(pressure(j0)))

c      ...compute z above reference level
       if(apmass.eq.0.0d0)goto 100
       o_z=0.0d0
       o_p=p0
       tkmid=0.5d0*(tk0+tk(j0))
       do j=j0,1,-1
          g=ggravity*mplanet/(rplanet+o_z)**2.0d0
          h=kboltz*tkmid/apmass/g
          z(j)=o_z-h*dlog(pressure(j)/o_p)
          if(j.eq.1)exit
          o_z=z(j)
          o_p=pressure(j)
          tkmid=0.5d0*(tk(j)+tk(j-1))
       enddo

c      ...compute z below reference level
       o_z=0.0d0
       o_p=p0
       tkmid=0.5d0*(tk0+tk(j0+1))
       do j=j0+1,nz
          g=ggravity*mplanet/(rplanet+o_z)**2.0d0
          h=kboltz*tkmid/apmass/g
          z(j)=o_z-h*dlog(pressure(j)/o_p)
          if(j.eq.nz)exit
          o_z=z(j)
          o_p=pressure(j)
          tkmid=0.5d0*(tk(j)+tk(j+1))
       enddo

c      ...compute vertical length of each layer
       do j=1,nz
          if(j.eq.1)then
             deltaz(j)=z(j)-z(j+1)
          elseif(j.eq.nz)then
             deltaz(j)=z(j-1)-z(j)
          else
             deltaz(j)=0.5d0*(z(j-1)-z(j))+0.5d0*(z(j)-z(j+1))
          endif
       enddo

c      ...compute gravity
       do j=1,nz
          ggrav(j)=ggravity*mplanet/(rplanet+z(j))**2.0d0
       enddo

       return

100    write(*,*)' E- Average particle mass is not provided'
       stop

       end




c_______________________________________________________________________

       subroutine compute_average_particle_mass
c      compute average particle mass in planet's atmosphere

       include 'pact.common'
       integer i,j
       real*8 rwk

c      ...compute average particle mass in planet's atmosphere
       apmass=0.0d0
       do j=1,nz
          rwk=0.0d0
          do i=1,nspec
             rwk=rwk+abun(i,j)*mspec(i)
          enddo
          apmass=apmass+rwk
       enddo
       apmass=apmass/dble(nz)

       return

       end




c_______________________________________________________________________

       subroutine read_species
c      ...read species included in the model
       
       use pact_data
       include 'pact.common'       
       integer i,j,k,l,telfab
       real*8 wkab,read_ab(nmaxz)
       character*(nmaxlenspec+1) string
       character*(nmaxcharfile) string_ab
       character*100 line
       logical fex

c      ...open species file
       open(unit=2,err=100,file=specfile,status='old')
       write(*,1000)specfile

c     ...read number of elements
       do 
          line(1:100)=''
          read(2,'(a100)',end=10)line
          if(line(1:1).eq.'!')cycle
          if(iabun.eq.0)read(line,*,err=110)nelem,telfab
          if(iabun.eq.1)read(line,*,err=110)nelem
          exit
       enddo
       if(nelem.gt.nmaxelem)stop' E- Too many elements'

c     ...read elements (and elemental abundances ?)
       i=0
       do while(i.lt.nelem)
          line(1:100)=''
          read(2,'(a100)',end=10)line
          if(line(1:1).eq.'!')cycle
          i=i+1
          if(iabun.eq.0)read(line,*,err=110)elem(i),elfab(i)
          if(iabun.eq.1)read(line,*,err=110)elem(i)
          do j=1,nmaxelem
             if(trim(elem(i)).ne.trim(element(j)))cycle
             mat(i)=atomic_weight(j)
             exit
          enddo
          if(mat(i).eq.0.0d0)goto 120
       enddo
       if(iabun.eq.0)then
          if(telfab.eq.2)then                    ! convert logarithmic to relative abundance
             do i=1,nelem
                elfab(i)=10.0d0**(elfab(i)-12.0d0)
             enddo
          endif
       endif
       write(*,1010)
       do i=1,nelem
          write(*,1020)trim(elem(i))
       enddo

c      ...read species
       nspec=0
       ncondensed=0
       do
          line(1:100)=''
          read(2,'(a100)',end=10)line
          if(line(1:1).eq.'!')cycle

c      ...read species name and properties
          nspec=nspec+1          
          if(nspec.gt.nmaxspec) stop ' E- Too many species'
          string=''
          if(iabun.eq.0)read(line,*,err=110)string
          if(iabun.eq.1)then
             read(line,*,err=110)string,string_ab
             read(string_ab,*,iostat=k)wkab
             if(k.eq.0)then
                abun(nspec,1:nz)=wkab
             else
                inquire(file=string_ab,exist=fex)
                if(.not.fex)then
                   write(*,*)trim(string_ab)
                   stop' E- Above abundance profile file does not exist'
                endif
                call read_abundance_profile(string_ab,read_ab)
                abun(nspec,1:nz)=read_ab(1:nz)
             endif
          endif
          if(len(trim(string)).gt.nmaxlenspec)goto 130
          spec(nspec)=trim(string)
          l=len(trim(string))

          condensed(nspec)=.false.
          if(spec(nspec)(l-2:l).eq.'(s)')condensed(nspec)=.true. ! solid species
          if(spec(nspec)(l-2:l).eq.'(l)')condensed(nspec)=.true. ! liquid species
          if(spec(nspec)(l-2:l).eq.'(c)')condensed(nspec)=.true. ! generic condensed species
          if(condensed(nspec))ncondensed=ncondensed+1
          if(.not.condensed(nspec).and.ncondensed.gt.0)
     #    stop ' E- Condensed species must be at end of file .spec'
       enddo
10     close(2)
       ngas=nspec-ncondensed
       if(nspec+ncondensed.gt.nmaxspec)stop' E- Need larger nmaxspec'

c      ...read NASA polynomials
       call read_therm

c      ...verify consistency between total number of atoms and those of each element
       do i=1,nspec
          k=0
          do j=1,nelem
             k=k+nat(i,j)
          enddo
          if(k.ne.nattot(i))goto 140
       enddo

c      ...evaluate species mass
       do i=1,nspec
          mspec(i)=0.0d0
          do k=1,nelem
             mspec(i)=mspec(i)+nat(i,k)*mat(k)
          enddo
          mspec(i)=mspec(i)*amu        ! amu -> g
       enddo

c      ...duplicated species?
       do i=2,nspec
          do j=1,i-1
             if(trim(spec(i)).eq.trim(spec(j)))goto 150
          enddo
       enddo

c      ...verify species charge
       j=0
       k=0
       ion=.false.
       do i=1,nspec
          if(charge(i).ne.0)ion=.true.
          if(charge(i).gt.0)j=j+1
          if(charge(i).lt.0)k=k+1
          if(charge(i).ne.0.and.condensed(i))goto 160
       enddo
       if(ion.and.(j.eq.0.or.k.eq.0))stop' E- Missing charged species'

c      ...identify electron, if present
       do i=1,nspec
          if(nattot(i).eq.0.and.charge(i).eq.-1)id_electr=i
       enddo
       if(ion.and.id_electr.eq.0)stop' E- electron is missing'

c     ...verify that number of atoms is non zero, except for electron
       do i=1,nspec
          if(ion.and.i.eq.id_electr)cycle
          if(nattot(i).eq.0)goto 170
       enddo

c      ...if species abundances are read from file, normalize them
       if(iabun.eq.1)then
          do j=1,nz
             wkab=0.0d0
             do i=1,nspec
                wkab=wkab+abun(i,j)
             enddo
             do i=1,nspec
                abun(i,j)=abun(i,j)/wkab
             enddo
          enddo
       endif

       write(*,1030)nspec

       return

100    write(*,*)' E- Error opening species file ',trim(specfile)
       stop
110    write(*,*)' E- Error reading species file at line:'
       write(*,*)' ',line
       stop
120    write(*,*)' E- Element ',trim(elem(i)),' not recognized'
       stop
130    write(*,*)' E- Too lenghty species: ',trim(spec(nspec))
       stop
140    write(*,*)' E- Inconsistent number of atoms for ',trim(spec(i))
       stop
150    write(*,*)' E- Species ',trim(spec(i)),' included twice'
       stop
160    write(*,*)' E- Charge in condensed species ',trim(spec(i))
       stop
170    write(*,*)' E- Number of atoms zero for: ',trim(spec(i))
       stop

1000   format(1x,' I- Reading species from file ',a45)
1010   format(1x,' I- Elements included:')
1020   format(1x,'       ',a)
1030   format(1x,' I- Data read for ',i4,' species')

       end




c_______________________________________________________________________

       subroutine read_therm
c      ...read therm data

       include 'pact.common'       
       integer i,j,k,read_nattot,read_nelem,read_charge,
     # read_type,read_nat(nmaxelemxs),read_ntemp,read_na
       real*8 read_a(9,nmaxtemp_therm),read_temp(nmaxtemp_therm+1)
       character*256 line
       character*(nmaxlenspec+1) string
       character*(nmaxlenspec) read_spec
       character*2 read_txt(2*nmaxelemxs),read_elem(nmaxelemxs)
       logical found(nmaxspec)

       write(*,1000)trim(thermfile)
       open(unit=2,file=trim(data_path)//trim(chemdir)//trim(thermfile),
     # status='old')
       found(:)=.false.

c Big loop over therm species
       do
          line(1:256)=''                         ! read parameters from therm file
          read(2,'(a256)',end=10)line
          if(line(1:1).eq.'!')cycle

          string=''
          read(line,*,err=100,end=10)            ! read therm parameters
     #    string,read_nattot,read_nelem,read_charge,read_type
          if(len(trim(string)).gt.nmaxlenspec)then
             write(*,*)' E- Many characters for (therm): ',string
             stop
          endif
          read_spec=trim(string)
          read(2,*,err=110)(read_txt(k),k=1,2*read_nelem)
          do k=1,2*read_nelem
             if(((-1)**k).lt.0)read_elem((k+1)/2)(1:2)=read_txt(k)(1:2)
             if(((-1)**k).gt.0)read(read_txt(k)(1:2),*)read_nat(k/2)
          enddo
          read(2,*,err=120)read_ntemp,(read_temp(k),k=1,read_ntemp+1)
          if(read_ntemp.gt.nmaxtemp_therm)
     #    stop' E- Too many temperature ranges'
          if(read_type.eq.1)read_na=7
          if(read_type.eq.7)read_na=7
          if(read_type.eq.9)read_na=9
          if(read_type.ne.1.and.read_type.ne.7.and.read_type.ne.9)then
             write(*,*)' E- unknown therm data type for: ',read_spec
             stop
          endif
          do k=1,read_ntemp
             read(2,*,err=120)(read_a(j,k),j=1,read_na)          
          enddo

          do i=1,nspec                                     ! verify whether species is included or not
             if(trim(read_spec).ne.trim(spec(i)))cycle     ! ... and if included, save parameters:
             if(found(i))goto 130
             found(i)=.true.
             nattot(i)=read_nattot
             nelemspec(i)=read_nelem
             charge(i)=read_charge
             do k=1,read_nelem
                do j=1,nelem
                   if(trim(read_elem(k)).ne.trim(elem(j)))cycle
                   nat(i,j)=read_nat(k)
                   exit
                enddo
             enddo
             type_therm(i)=read_type
             natherm(i)=read_na
             ntemp_therm(i)=read_ntemp
             do k=1,ntemp_therm(i)+1
                temp_therm(i,k)=read_temp(k)
             enddo
             do k=1,ntemp_therm(i)
                do j=1,natherm(i)
                   atherm(i,j,k)=read_a(j,k)
                enddo
             enddo
             exit
          enddo

          read_spec(1:nmaxlenspec)=''            ! reset 'read_...' variables
          read_nattot=0
          read_nelem=0
          read_charge=0
          read_type=0
          read_na=0
          read_txt(:)(1:2)=''
          read_elem(:)(1:2)=''
          read_nat(:)=0
          read_ntemp=0
          read_temp(:)=0.0d0
          read_a(:,:)=0.0d0

c End of big loop over therm species
       enddo
10     close(2)

c     ...verify that all species have therm data
       do i=1,nspec
          if(nattot(i).eq.0)goto 140
       enddo

       return

100    write(*,*)' E- Error reading therm file at line:'
       write(*,*)'    ',line
       stop
110    write(*,*)' E- Error reading therm file elements for: ',
     # trim(read_spec)
       stop
120    write(*,*)' E- Error reading therm data for: ',trim(read_spec)
       stop
130    write(*,*)' E- Duplicated therm data for species ',trim(spec(i))
       stop
140    write(*,*)' E- No therm data for species: ',trim(spec(i))
       stop

1000   format(1x,' I- Reading therm data from ',a)

       end




c_______________________________________________________________________

       subroutine read_reactions
c      ...read reactions included in the model

       use pact_data
       include 'pact.common'       
       integer i,j,k,l,m,ir,ip,ilen,nrf,nrd,balr,balp,nchunk,chunk(10),
     # nr,np
       character*1024 line
       character*(nmaxcharfile) usedreacfile
       character*1 lerror
       logical pass,inclu,rused(nmaxreag),pused(nmaxprod),found

       ilen=0                ! to avoid warning
       ip=0                  ! ...

c     ...open file to write reactions used
       usedreacfile(1:nmaxcharfile)=''
       usedreacfile(1:ncharmodel)=namemodel(1:ncharmodel)
       usedreacfile(ncharmodel+1:ncharmodel+14)='_reac_used.dat'
       open(unit=3,file=usedreacfile,status='unknown')

c      ...read reaction file
       i=0
       m=0
       write(*,1000)reacfile
       open(unit=2,file=trim(data_path)//trim(chemdir)//trim(reacfile),
     # status='old')
       do
          line(1:1024)=''
          read(2,'(a1024)',end=10)line
          if(line(1:1).eq.'!'.or.line(1:1024).eq.' ')cycle
          i=i+1
          if(i.gt.nmaxreac)stop ' E- Too many reactions'

c      ...reset variables, safety measure; needed for example if already initialized for a reaction not included
          reac(i)=''
          nreag(i)=0
          nprod(i)=0
          reag(i,:)=''
          prod(i,:)=''
          idreag(i,:)=0
          idprod(i,:)=0
          lerror=''

          nchunk=0
          chunk(:)=0
          do j=1,1024
             if(line(j:j).eq.':')then
                nchunk=nchunk+1
                chunk(nchunk)=j
                if(nchunk.eq.10)exit
             endif
          enddo
          if(nchunk.lt.10)nchunk=4
          l=chunk(1)-1
          if(l.gt.nmaxlenreac)stop' E- Too long reaction'
          reac(i)=line(1:chunk(1)-1)

c      ...decode reaction into reagents and products
          ir=1
          pass=.false.
          do j=1,l
             if(reac(i)(j:j).eq.' ')then
                cycle           
             elseif(reac(i)(j-1:j+1).eq.' + ')then
                if(.not.pass)ir=ir+1
                if(pass)     ip=ip+1
             elseif(reac(i)(j-1:j+1).eq.' = ')then
                pass=.true.
                ip=1
             else
                if(j.eq.1.or.reac(i)(j-1:j-1).eq.' ')then
                   ilen=1
                   if(.not.pass)then
                      reag(i,ir)(:)=''
                      reag(i,ir)(ilen:ilen)=reac(i)(j:j)
                   else
                      prod(i,ip)(:)=''
                      prod(i,ip)(ilen:ilen)=reac(i)(j:j)
                   endif
                else
                   ilen=ilen+1
                   if(ilen.gt.nmaxlenspec)stop' E- Too lengthy species'
                   if(.not.pass)then
                      reag(i,ir)(ilen:ilen)=reac(i)(j:j)
                   else
                      prod(i,ip)(ilen:ilen)=reac(i)(j:j)
                   endif
                endif
             endif
          enddo
          nreag(i)=ir
          nprod(i)=ip
          if(nreag(i).le.0.or.nprod(i).le.0)stop ' E- No reag/prod'
          if(nreag(i).gt.nmaxreag)stop ' E- Too many reagents'
          if(nprod(i).gt.nmaxprod)stop ' E- Too many products'

c      ...check whether species involved are actually considered in the model
          inclu=.true.
          do j=1,nreag(i)
             if(trim(reag(i,j)).eq.'M')cycle
             if(trim(reag(i,j)).eq.'(M)')cycle
             do k=1,nspec
                if(trim(reag(i,j)).ne.trim(spec(k)))cycle
                idreag(i,j)=k
                exit
             enddo
             if(idreag(i,j).eq.0)inclu=.false.
          enddo
          do j=1,nprod(i)
             if(trim(prod(i,j)).eq.'M')cycle
             if(trim(prod(i,j)).eq.'(M)')cycle
             do k=1,nspec
                if(trim(prod(i,j)).ne.trim(spec(k)))cycle
                idprod(i,j)=k
                exit
             enddo
             if(idprod(i,j).eq.0)inclu=.false.
          enddo
          if(inclu)then
             write(3,'(a)')trim(line)                
          else
             i=i-1
             cycle
          endif

c      ...do verification for pressure-dependent reactions and remove collider M from stoichiometry
          if(nchunk.eq.4)then
             rtype(i)=1                                    ! bimolecular reaction
          elseif(nchunk.eq.10)then
             if(trim(reag(i,nreag(i))).eq.'M'.and.
     #          trim(prod(i,nprod(i))).eq.'M')then
                rtype(i)=2                                 ! pressure dependent reaction
                nreag(i)=nreag(i)-1
                nprod(i)=nprod(i)-1
             elseif(trim(reag(i,nreag(i))).eq.'(M)'.and.
     #              trim(prod(i,nprod(i))).eq.'(M)')then
                nreag(i)=nreag(i)-1
                nprod(i)=nprod(i)-1
                rtype(i)=3                                 ! special pressure dependent reaction
             else
                stop' E- Unknown reaction type'
             endif
          endif

c      ...check the reaction mass and charge balance
          do k=0,nelem
             balr=0
             do j=1,nreag(i)
                if(k.eq.0)balr=balr+charge(idreag(i,j))
                if(k.gt.0)balr=balr+nat(idreag(i,j),k)
             enddo
             balp=0
             do j=1,nprod(i)
                if(k.eq.0)balp=balp+charge(idprod(i,j))
                if(k.gt.0)balp=balp+nat(idprod(i,j),k)
             enddo
             if(balr.ne.balp)goto 110
          enddo             

c      ...read reaction rate coefficient data from reaction file
          if(rtype(i).eq.1)then
             read(line(chunk(1)+1:chunk(2)-1),*,err=100)
     #       alpha(i),beta(i),gamma(i)
             read(line(chunk(2)+1:chunk(3)-1),*,err=100)
     #       tkmin(i),tkmax(i)
             read(line(chunk(3)+1:chunk(4)-1),*,err=100)lerror
             comment(i)=line(chunk(4)+1:1024)

          else
             m=m+1
             if(m.gt.nmaxmreac)stop ' E- Too many p-dependent reactions'
             idmreac(i)=m
             if(line(chunk(1)+1:chunk(2)-1).ne.'')
     #          read(line(chunk(1)+1:chunk(2)-1),*,err=100)
     #          alpha(i),beta(i),gamma(i)
             if(line(chunk(2)+1:chunk(3)-1).ne.'')
     #          read(line(chunk(2)+1:chunk(3)-1),*,err=100)
     #          alphainf(m),betainf(m),gammainf(m)
             if(line(chunk(3)+1:chunk(4)-1).ne.'')
     #          read(line(chunk(3)+1:chunk(4)-1),*,err=100)afcent(m)
             if(line(chunk(4)+1:chunk(5)-1).ne.'')
     #          read(line(chunk(4)+1:chunk(5)-1),*,err=100)bfcent(m)
             if(line(chunk(5)+1:chunk(6)-1).ne.'')
     #          read(line(chunk(5)+1:chunk(6)-1),*,err=100)mcoll(m)
             if(line(chunk(6)+1:chunk(7)-1).ne.'')
     #          read(line(chunk(6)+1:chunk(7)-1),*,err=100)
     #          tkmin(i),tkmax(i)
             if(line(chunk(7)+1:chunk(8)-1).ne.'')
     #          read(line(chunk(7)+1:chunk(8)-1),*,err=100)
     #          tkmininf(m),tkmaxinf(m)
             if(line(chunk(8)+1:chunk(9)-1).ne.'')
     #          read(line(chunk(8)+1:chunk(9)-1),*,err=100)
     #          tkminfc(m),tkmaxfc(m)
             read(line(chunk(9)+1:chunk(10)-1),*,err=100)lerror
             comment(i)=line(chunk(10)+1:1024)

             if(line(chunk(1)+1:chunk(2)-1).eq.''.and.
     #          line(chunk(2)+1:chunk(3)-1).eq.'')goto 120
          endif

c      ...evaluate number of C atoms involved
          l=0
          do j=1,nreag(i)
             if(idreag(i,j).eq.0)cycle
             do k=1,nelem
                if(trim(elem(k)).ne.'C')cycle
                l=l+nat(idreag(i,j),k)
                exit
             enddo
          enddo
          natc(i)=l
          l=0
          do j=1,nprod(i)
             if(idprod(i,j).eq.0)cycle
             do k=1,nelem
                if(trim(elem(k)).ne.'C')cycle
                l=l+nat(idprod(i,j),k)
                exit
             enddo
          enddo
          if(l.ne.natc(i))stop' E- Carbon not balanced'

c      ...adopt crude estimates if missing information in pressure dependent reaction
          if(rtype(i).gt.1)then
             if(alpha(i).eq.0.0d0)then
                alpha(i)=alphainf(m)*10.0d0**(-20.8d0+1.6d0*natc(i))
                beta(i)=betainf(m)-2.0d0
             endif
             if(alphainf(m).eq.0.0d0)then
                alphainf(m)=alpha(i)/10.0d0**(-20.8d0+1.6d0*natc(i))
                betainf(m)=beta(i)+2.0d0
             endif
             if(afcent(m).eq.0.0d0)afcent(m)=0.5d0
             if(mcoll(m).ne.'')then
                found=.false.
                do l=1,nmaxmcoll
                   if(trim(mcoll(m)).eq.trim(mcollider(l)))found=.true.
                enddo
                if(.not.found)then
                   write(*,*)trim(reac(i))
                   stop' E- Unknown M collider in above reaction'
                endif
             endif
          endif

          if(lerror.eq.'A')then       ! factor 2
             kerror(i)=10.0d0**0.30d0
          elseif(lerror.eq.'B')then   ! factor 4.5
             kerror(i)=10.0d0**0.65d0
          elseif(lerror.eq.'C')then   ! factor 10
             kerror(i)=10.0d0**1.00d0
          else
             write(*,*)trim(reac(i))
             stop' E- Unknown error code in above reaction'
          endif

c     ...check for reactions with no net change, e.g., dummy reaction He = He
          if(nreag(i).eq.nprod(i))then
             nr=0
             rused(:)=.false.
             do j=1,nreag(i)
                found=.false.
                do k=1,nprod(i)
                   if(idreag(i,j).ne.idprod(i,k))cycle
                   if(found)cycle
                   if(rused(k))cycle
                   found=.true.
                   rused(k)=.true.
                enddo
                if(found)nr=nr+1
             enddo
             if(nr.eq.nreag(i))cycle
          endif

c     ...include reverse reaction
          i=i+1
          idreac_rev(i)=i-1
          do j=1,nmaxlenreac
             if(reac(i-1)(j:j).eq.'=')k=j
             if(reac(i-1)(j:j).ne.' ')l=j
          enddo
          reac(i)(1:l-k+1)=reac(i-1)(k+1:l+1)
          reac(i)(l-k+2:l-k+2)='='
          reac(i)(l-k+3:l+1)=reac(i-1)(1:k-1)
          nreag(i)=nprod(i-1)
          nprod(i)=nreag(i-1)
          if(nreag(i).gt.nmaxreag)goto 130
          if(nprod(i).gt.nmaxprod)goto 130
          do j=1,nreag(i)
             reag(i,j)=prod(i-1,j)
             idreag(i,j)=idprod(i-1,j)
          enddo
          do j=1,nprod(i)
             prod(i,j)=reag(i-1,j)
             idprod(i,j)=idreag(i-1,j)
          enddo
          rtype(i)=0
          natc(i)=natc(i-1)
          kerror(i)=kerror(i-1)

       enddo
10     close(2)
       nreac=i
       nmreac=m
       close(3) ! close ...reactions_used.dat file

c      ...duplicated reactions?
       do i=2,nreac
          do j=1,i-1
             if(nreag(i).ne.nreag(j))cycle
             if(nprod(i).ne.nprod(j))cycle
             nr=0
             rused(:)=.false.
             do k=1,nreag(i)
                found=.false.
                do l=1,nreag(j)
                   if(idreag(i,k).ne.idreag(j,l))cycle
                   if(found)cycle
                   if(rused(l))cycle
                   found=.true.
                   rused(l)=.true.
                enddo
                if(found)nr=nr+1
             enddo
             if(nr.ne.nreag(i))cycle
             np=0
             pused(:)=.false.
             do k=1,nprod(i)
                found=.false.
                do l=1,nprod(j)
                   if(idprod(i,k).ne.idprod(j,l))cycle
                   if(found)cycle
                   if(pused(l))cycle
                   found=.true.
                   pused(l)=.true.
                enddo
                if(found)np=np+1
             enddo
             if(np.eq.nprod(i))goto 140
          enddo
       enddo

c      ...check that all species have at least one formation/destruction reaction
       do i=1,nspec
          nrd=0
          nrf=0
          do j=1,nreac
             do k=1,nreag(j)
                if(idreag(j,k).eq.i)nrd=nrd+1
             enddo
             do k=1,nprod(j)
                if(idprod(j,k).eq.i)nrf=nrf+1
             enddo
          enddo
          if(nrd.eq.0)goto 150
          if(nrf.eq.0)goto 160
          if(nrd.eq.1)
     #    write(*,*)' W- Only 1 destruction path for ',trim(spec(i))
          if(nrf.eq.1)
     #    write(*,*)' W- Only 1 formation path for ',trim(spec(i))
       enddo

c      ...build rate equations, i.e. dn_i/dt=...
       do i=1,nspec
          nterm(i)=0
          do j=1,nreac
             do k=1,nreag(j)             ! add sink terms
                if(idreag(j,k).eq.i)then
                   nterm(i)=nterm(i)+1
                   idterm(i,nterm(i))=j
                   sterm(i,nterm(i))=-1
                   fterm(i,nterm(i))=1
                   if(k.eq.nreag(j))goto 20 
                   do l=k+1,nreag(j)
                      if(idreag(j,l).eq.i)
     #                fterm(i,nterm(i))=fterm(i,nterm(i))+1
                   enddo
                   goto 20
                endif
             enddo

20           do k=1,nprod(j)             ! add source terms
                if(idprod(j,k).eq.i)then
                   nterm(i)=nterm(i)+1
                   idterm(i,nterm(i))=j
                   sterm(i,nterm(i))=1
                   fterm(i,nterm(i))=1
                   if(k.eq.nprod(j))exit 
                   do l=k+1,nprod(j)
                      if(idprod(j,l).eq.i)
     #                fterm(i,nterm(i))=fterm(i,nterm(i))+1
                   enddo
                   exit
                endif
             enddo
          enddo
          if(nterm(i).gt.nmaxrs)goto 170
       enddo

       write(*,1010)nreac

       return

100    write(*,*)trim(line)
       stop' E- Error reading reactions file at above line'
110    write(*,*)trim(reac(i))
       stop' E- Balance not fulfilled for above reaction'
120    write(*,*)trim(reac(i))
       stop' E- Missing info in above pressure dependent reaction'
130    write(*,*)trim(reac(i))
       stop' E- Too many reagents/products in above reaction'
140    write(*,*)reac(i)
       write(*,*)reac(j)
       stop' E- Duplicated reactions'
150    write(*,*)' E- No destruction path for species:',trim(spec(i))
       stop
160    write(*,*)' E- No formation path for species:',trim(spec(i))
       stop
170    write(*,*)' E- Too many reactions involving ',trim(spec(i))
       stop

1000   format(1x,' I- Reading reactions from ',a45)
1010   format(1x,' I- ',i4,' reactions included in the model')

       end




c_______________________________________________________________________

       subroutine read_eddy
c      ...read and create grid for eddy diffusion coefficient profile

       include 'pact.common'       
       integer nmaxeddy_read
       parameter(nmaxeddy_read=10000)
       integer i,j,xnz,j1,j2,js,jj
       real*8 slope
       real*8, allocatable, dimension(:) :: xpr,xeddy
       character*100 line

       if(.not.diffusion)return

       allocate(xpr(nmaxeddy_read))
       allocate(xeddy(nmaxeddy_read))

c      ...approximate Kzz profile based on Eq (1) of Moses et al (2022)
       if(trim(eddyfile).eq.'None')then
          call get_analytical_eddy
          return
       endif

c      ...read eddy diffusion coefficient profile
       open(unit=2,file=eddyfile,status='old')
       i=0
       do
          line(1:100)=''
          read(2,'(a100)',err=100,end=10)line
          if(line(1:1).eq.'!')cycle

          i=i+1
          if(i.gt.nmaxeddy_read)stop ' E- Eddy grid too fine'
          read(line(1:100),*,err=100)xpr(i),xeddy(i)       ! pressure [bar], K_eddy [cm2 s-1]
          xpr(i)=xpr(i)*1.0d6                              ! bar -> dyn cm-2
       enddo
10     close(2)
       xnz=i

c      ...interpolate to (z,p,T) grid
       if(xpr(1).lt.xpr(xnz))then                          ! input (p,eddy) ordered from top to bottom
          j1=1
          j2=xnz
          js=1
          if(xpr(1).gt.pmin)       goto 110
          if(xpr(xnz).lt.pmax)     goto 120
          do j=1,xnz-1
             if(xpr(j).ge.xpr(j+1))goto 130
          enddo
       else                                                ! input (p,eddy) ordered from bottom to top
          j1=xnz
          j2=1
          js=-1
          if(xpr(xnz).gt.pmin)     goto 110
          if(xpr(1).lt.pmax)       goto 120
          do j=1,xnz-1
             if(xpr(j).le.xpr(j+1))goto 130
          enddo
       endif
       do j=1,nz                                           ! create K_eddy grid
          do jj=j1,j2,js
             if(xpr(jj).gt.pressure(j))exit
          enddo
          slope=(dlog10(xeddy(jj))-dlog10(xeddy(jj-js)))/
     #    (dlog10(xpr(jj))-dlog10(xpr(jj-js)))
          kdif(j)=xeddy(jj-js)*(pressure(j)/xpr(jj-js))**slope
       enddo

       deallocate(xpr)
       deallocate(xeddy)

       return

100    write(*,*)' E- Error reading eddy diffusion file: ',eddyfile
       stop
110    stop ' E- Requested pressure at top outside input eddy grid'
120    stop ' E- Requested pressure at bottom outside input eddy grid'
130    write(*,*)' E- Discontinuity in eddy grid at j=',j
       stop

       end




c_______________________________________________________________________

       subroutine get_analytical_eddy
c      ...approximate Kzz profile based on Eq (1) of Moses et al (2022)

       include 'pact.common'       
       integer j
       real*8 hatm,kconvect,fint,h1,h2,slope

       if(.not.diffusion)return
       if(trim(eddyfile).ne.'None')return
       if(jadiab.lt.1.or.jadiab.gt.nz)
     # stop' E- Unknown radiative-convective transition'

c      ...Kzz in convective atmosphere from free convection, Eq (5) in Ackerman & Marley (2001)
c      with the mixing length set equal to the atmospheric scale height (L = H), which is adequate for the convective fregion
       j=jadiab
       hatm=kboltz*tk(j)/apmass/ggrav(j)
       fint=sigmasb*tintplanet**4.0d0
       kconvect=hatm/3.0d0*
     # (kboltz*fint/apmass/apmass/density(j)/cplayer(j))**(1.0d0/3.0d0)

       j=0                                                 ! evaluate atmospheric scale height at 1 mbar (1000 dyn cm-2) level
       do while(j.lt.nz)
          j=j+1
          if(pressure(j).gt.1.0d3)exit
       enddo
       j=min(max(j-1,1),nz-1)
       h1=kboltz*tk(j)/apmass/ggrav(j)
       h2=kboltz*tk(j+1)/apmass/ggrav(j+1)
       slope=(h2-h1)/(dlog10(pressure(j+1))-dlog10(pressure(j)))
       hatm=h1+slope*(dlog10(1.0d3)-dlog10(pressure(j)))
       do j=1,nz
          kdif(j)=5.0d8*(pressure(j)/1.0d6)**(-0.5d0)*     ! radiative atmosphere
     #    (hatm/6.2d7)*(tplanet_th/1450.0d0)**4.0d0
          if(kdif(j).gt.1.0d11) kdif(j)=1.0d11             ! top atmosphere if Kzz > 1e11
          if(jadiab.eq.nz)cycle
          if(j.ge.jadiab)       kdif(j)=kconvect           ! bottom convective atmosphere
       enddo

       return

       end




c_______________________________________________________________________

       subroutine read_photon
c      ...create wavelength grid and read stellar spectrum and photo cross sections

       include 'pact.common'       
       integer i,j,k,l,ii,s,nxy,ir,ip,ilen,wklenreag(2),balr,balp,npfnd,
     # idk,idl,nphotochannel
       real*8 dwlth
       real*8, allocatable, dimension(:)   :: x,y
       real*8, allocatable, dimension(:,:) :: qy
       character*1024 line
       character*(nmaxcharstring) datafile
       character*(nmaxcharfile) phfile
       character*(nmaxlenspec) wkreag(nmaxphotochannel)
       character*50 wkphotoreac(nmaxphotochannel)
       character*1 lerror
       logical pass,pfnd(nmaxprod),inclu

       if(.not.photochem)return

       ilen=0                                    ! to avoid warning
       nwlth=1000                                ! number of wavelengths
       wlthmin=1.0d0                             ! low wavelegth cutoff [nm]
       wlthmax=750.0d0                           ! up  wavelegth cutoff [nm]
       if(nwlth.gt.nmaxwlth)stop 'E- Wavelength grid is too fine'

c      ...create wavelength grid
       dwlth=(wlthmax-wlthmin)/dble(nwlth-1)
       do i=1,nwlth
          wlth(i)=wlthmin+(i-1)*dwlth
       enddo


c Stellar spectrum.........................................

c      ...read stellar spectrum from starfile
       allocate(x(nmaxxy),y(nmaxxy))
       write(*,1000)starfile
       open(unit=2,file=trim(data_path)//trim(stardir)//trim(starfile),
     # status='old')
       i=0
       do
          line(1:1024)=''
          read(2,'(a1024)',err=100,end=10)line
          if(line(1:1).eq.'!')cycle

          i=i+1
          if(i.gt.nmaxxy)stop ' E- Stellar spectrum too fine'
          read(line(1:1024),*,err=100)x(i),y(i)            ! wavelength [um], intensity [erg s-1 cm-2 Hz-1 sr-1]
          x(i)=x(i)*1.0d3                                  ! wavelength [=] nm
          y(i)=y(i)*1.0d7*clight/x(i)/x(i)                 ! I [=] erg s-1 cm-2 nm-1 sr-1
          y(i)=y(i)*pi*(rstar/dist)**2.0d0                 ! I [=] erg s-1 cm-2 nm-1
          y(i)=y(i)*x(i)*1.0d-7/hplanck/clight             ! I [=] photon s-1 cm-2 nm-1
       enddo
10     close(2)
       nxy=i

c      ...resample stellar flux to the wavelength grid
       call resample_to_wlth(x,y,nxy,fphot)


c Photo cross sections.....................................

c      ...read photo cross section info file
       allocate(qy(nmaxphotochannel,nmaxxy))
       write(*,1010)
       datafile(1:nmaxcharstring)=''
       datafile=trim(data_path)//trim(photodir)//'photoreac.ipho'
       open(unit=2,file=datafile,status='old')
       s=0                                       ! number of photospecies
       k=0                                       ! number of photoreactions
       do
          line(1:1024)=''
          read(2,'(a1024)',err=110,end=20)line
          if(line(1:1).eq.'!')cycle

c      ...open and read photo cross section file
          phfile(1:nmaxcharfile)=''
          read(line(1:1024),*,err=110)phfile
          open(unit=3,
     #    file=trim(data_path)//trim(photodir)//trim(phfile),
     #    status='old')
          j=-1
          do
             line(1:1024)=''
             read(3,'(a1024)',err=120,end=30)line
             if(line(1:1).eq.'!')cycle

c     ...read photo reaction(s)
             if(j.eq.-1)then
                lerror=''
                read(line,*,err=120)nphotochannel,lerror
                if(nphotochannel.gt.nmaxphotochannel)      goto 130
                do ii=1,nphotochannel
                   wkphotoreac(ii)(1:50)=''
                   read(3,'(a50)',err=120)wkphotoreac(ii)
                enddo
                j=0
                cycle
             endif

c      ...read photo cross section: wavelength[nm], cross section[cm2], quantum yield(s)
             j=j+1
             if(j.gt.nmaxxy)                               goto 140
             read(line(1:1024),*,err=120)x(j),y(j),
     #       (qy(ii,j),ii=1,nphotochannel)
             if(j.gt.1.and.x(j).le.x(j-1))                 goto 150
          enddo
30        close(3)
          nxy=j
          if(x(1).ge.wlthmax.or.x(nxy).le.wlthmin)         goto 160

c >>>  Absorption cross sections of species

c      ...verify that reagent is included and store absorption cross section
          do ii=1,nphotochannel
             wkreag(ii)(1:nmaxlenspec)=''
             read(wkphotoreac(ii),*,err=120)wkreag(ii)
             if(ii.gt.1)then
                if(trim(wkreag(ii)).ne.trim(wkreag(ii-1))) goto 170
             endif
          enddo
          inclu=.false.
          do i=1,nspec
             if(trim(wkreag(1)).eq.trim(spec(i)))then
                inclu=.true.
                s=s+1
                if(s.gt.nmaxspec)stop ' E- Many photo species'
                idphotospec(s)=i
                call resample_to_wlth(x,y,nxy,csphotabs(s,:))
                exit
             endif
          enddo
          if(.not.inclu)then
             write(*,1020)wkreag(1)
             cycle
          endif

c >>>  Cross sections of photo reactions (photochemical channels)

          do ii=1,nphotochannel
             k=k+1
             if(k.gt.nmaxphotoreac)stop ' E- Many photo processes'
             photoreac(k)(1:50)=''
             photoreac(k)(1:50)=wkphotoreac(ii)(1:50)

c      ...decode photo reaction into reagents and products
             ir=1
             ip=1
             pass=.false.
             wkreag(1)(1:nmaxlenspec)=''
             wkreag(2)(1:nmaxlenspec)=''
             photoreag(k)(1:nmaxlenspec)=''
             do i=1,nmaxprod
                photoprod(k,i)(1:nmaxlenspec)=''
             enddo
             do j=1,50
                if(photoreac(k)(j:j).eq.' ')then
                   cycle           
                elseif(photoreac(k)(j-1:j+1).eq.' + ')then
                   if(.not.pass)ir=ir+1
                   if(pass)     ip=ip+1
                elseif(photoreac(k)(j-1:j+1).eq.' = ')then
                   pass=.true.
                else
                   ilen=ilen+1
                   if(j.eq.1.or.photoreac(k)(j-1:j-1).eq.' ')ilen=1
                   if(.not.pass)then
                      wkreag(ir)(ilen:ilen)=photoreac(k)(j:j)
                      wklenreag(ir)=ilen
                   else
                      photoprod(k,ip)(ilen:ilen)=photoreac(k)(j:j)
                   endif
                endif
             enddo
             photoreag(k)(1:nmaxlenspec)=wkreag(1)(1:nmaxlenspec)
             nphotoprod(k)=ip

c      ...check for possible errors in the photo reaction syntax
             if(ir.ne.2.or.ip.le.0)                        goto 180
             if(wklenreag(2).ne.6)                         goto 180
             if(wkreag(2)(1:6).ne.'photon')                goto 180
             if(ip.gt.nmaxprod)                            goto 180

c      ...verify that species involved are included in the model
             inclu=.true.
             do i=1,nspec
                if(trim(photoreag(k)).eq.trim(spec(i)))then
                   idphotoreag(k)=i
                   exit
                endif
             enddo
             if(idphotoreag(k).eq.0)inclu=.false.
             do j=1,ip
                do i=1,nspec
                   if(trim(photoprod(k,j)).eq.trim(spec(i)))then
                      idphotoprod(k,j)=i
                      exit
                   endif
                enddo
                if(idphotoprod(k,j).eq.0)then
                   if(trim(photoprod(k,j)).ne.'absorption')
     #             write(*,1030)trim(wkphotoreac(ii)),
     #             trim(photoprod(k,j))
                   inclu=.false.
                endif
             enddo
             if(.not.inclu)then
                k=k-1
                cycle
             endif

c      ...check the photo reaction mass and charge balance
             do i=0,nelem
                balr=0
                if(i.eq.0)balr=charge(idphotoreag(k))
                if(i.gt.0)balr=nat(idphotoreag(k),i)
                balp=0
                do j=1,ip
                   if(i.eq.0)balp=balp+charge(idphotoprod(k,j))
                   if(i.gt.0)balp=balp+nat(idphotoprod(k,j),i)
                enddo
                if(balr.ne.balp)                           goto 190
             enddo             

c      ...check that the photo reaction did not appear before
             if(k.gt.1)then
                do j=1,k-1
                   if(idphotoreag(k).ne.idphotoreag(j))cycle
                   if(nphotoprod(k).ne.nphotoprod(j))cycle
                   pfnd(1:nmaxprod)=.false.
                   npfnd=0
                   do i=1,nphotoprod(k)
                      idk=idphotoprod(k,i)
                      do l=1,nphotoprod(j)
                         idl=idphotoprod(j,l)
                         if(idk.eq.idl.and..not.pfnd(l))then
                            pfnd(l)=.true.
                            npfnd=npfnd+1
                            exit
                         endif
                      enddo
                   enddo
                   if(npfnd.eq.nphotoprod(k))              goto 200
                enddo
             endif

c      ...convert quantum yield to cross section and resample to the wavelength grid
             do j=1,nxy
                qy(ii,j)=qy(ii,j)*y(j)
             enddo
             call resample_to_wlth(x,qy(ii,:),nxy,csphot(k,:))

c      ...assign error
             if(lerror.eq.'A')then       ! factor 2
                jerror(k)=10.0d0**0.30d0
             elseif(lerror.eq.'B')then   ! factor 4.5
                jerror(k)=10.0d0**0.65d0
             elseif(lerror.eq.'C')then   ! factor 10
                jerror(k)=10.0d0**1.00d0
             else
                write(*,*)trim(photoreac(k))
                stop' E- Unknown error code in above photoreaction'
             endif

          enddo
       enddo
20     close(2)
       nphotospec=s
       nphotoreac=k
       write(*,1040)nphotoreac
       deallocate(x,y,qy)

c      ...build rate equations, i.e. dn_i/dt=...
       do i=1,nspec
          nphototerm(i)=0
          do j=1,nphotoreac
             if(idphotoreag(j).eq.i)then ! add sink terms
                nphototerm(i)=nphototerm(i)+1
                idphototerm(i,nphototerm(i))=j
                sphototerm(i,nphototerm(i))=-1
                fphototerm(i,nphototerm(i))=1
             endif

             do k=1,nphotoprod(j)        ! add source terms
                if(idphotoprod(j,k).eq.i)then
                   nphototerm(i)=nphototerm(i)+1
                   idphototerm(i,nphototerm(i))=j
                   sphototerm(i,nphototerm(i))=1
                   fphototerm(i,nphototerm(i))=1
                   if(k.eq.nphotoprod(j))exit 
                   do l=k+1,nphotoprod(j)
                      if(idphotoprod(j,l).eq.i)
     #                fphototerm(i,nphototerm(i))=
     #                fphototerm(i,nphototerm(i))+1
                   enddo
                   exit
                endif
             enddo
          enddo
          if(nphototerm(i).gt.nmaxps)then
             write(*,*)' E- Species ',spec(i)
             stop      '    involved in too many photo reactions'
          endif
       enddo

c      ...ensure that all molecules have some photodissociation channel
       do i=1,nspec
          if(nattot(i).lt.2)cycle
          inclu=.false.
          do k=1,nphotoreac
             if(idphotoreag(k).eq.i)inclu=.true.
          enddo
          if(.not.inclu)                                   goto 210
       enddo

       return

100    write(*,*)' E- Error reading stellar spectrum file: ',
     # trim(starfile)
       stop
110    stop' E- Error reading file photoreac.ipho'
120    write(*,*)' E- Error reading photo cross section file: ',
     # trim(phfile)
       stop
130    write(*,*)' E- Many points in photo cross section file: ',
     # trim(phfile)
       stop
140    write(*,*)' E- Many photochemical channels in: ',trim(phfile)
       stop
150    write(*,*)' E- Bad wavelength ordering in file: ',trim(phfile)
       write(*,*)'    at wavelength ',x(j),' nm'
       stop
160    write(*,*)' E- File: ',trim(phfile),' out of wavelength range'
       stop
170    write(*,*)' E- Different reagents in : ',trim(phfile)
       stop
180    write(*,*)' E- In photo reaction: ',trim(photoreac(k))
       stop      '    There is an error in the syntax'
190    write(*,*)' E- In photo reaction: ',trim(photoreac(k))
       stop      '    Balance not fullfilled'
200    write(*,*)' E- : ',trim(photoreac(j))
       write(*,*)'    : ',trim(photoreac(k))
       stop      '    The above photo reactions are repeated'
210    write(*,*)' E- Molecule ',trim(spec(i)),
     # ' has no photodissociation channel'
       stop

1000   format(1x,' I- Reading stellar spectrum from ',a45)
1010   format(1x,' I- Reading photo cross section file photoreac.ipho')
1020   format(1x,' W- Photo-species ',a10,' not included in the model')
1030   format(1x,' W- In photoreaction: ',a,/,
     #        1x,'    product ',a,' is not included in the model')
1040   format(1x,' I- ',i4,' photo reactions included in the model')

       end




c_______________________________________________________________________

       subroutine resample_to_wlth(x,y,nxy,yout)
c      ...resample y(x) to wavelenth grid 

       include 'pact.common'       
       integer nxy,i,j,jl,ju
       real*8 x(*),y(*),yout(nmaxwlth),wl,wu,slope,x1,x2,y1,y2

       do i=1,nwlth
          yout(i)=0.0d0
          if(i.eq.1)then
             wl=wlth(i)
             wu=wlth(i)+0.5d0*(wlth(i+1)-wlth(i))
          elseif(i.eq.nwlth)then
             wl=wlth(i)-0.5d0*(wlth(i)-wlth(i-1))
             wu=wlth(i)
          else
             wl=wlth(i)-0.5d0*(wlth(i)-wlth(i-1))
             wu=wlth(i)+0.5d0*(wlth(i+1)-wlth(i))
          endif
          if(wu.le.x(1).or.wl.ge.x(nxy))cycle

          j=1
          do while(x(j).lt.wl.and.j.lt.nxy)
             j=j+1
          enddo
          j=max(j-1,1)
          jl=j
          do while(x(j).lt.wu.and.j.lt.nxy)
             j=j+1
          enddo
          ju=min(max(j-1,1),nxy-1)

          do j=jl,ju
             x1=max(x(j),wl)
             x2=min(x(j+1),wu)
             slope=(y(j+1)-y(j))/(x(j+1)-x(j))
             y1=y(j)+slope*(x1-x(j))
             y2=y(j)+slope*(x2-x(j))
             yout(i)=yout(i)+0.5d0*(y1+y2)*(x2-x1)
          enddo
          yout(i)=yout(i)/(wu-wl)
       enddo

       return
       end




c_______________________________________________________________________

       subroutine read_ktable
c      read k-table data

       include 'pact.common'
       integer i,j,k,l,m,nspec_ktab0
       real*8 wav
       character*(nmaxcharstring) datafile
       character*12 kfile
       character*(nmaxlenspec) spec_ktab0(nmaxspec_ktab)

       write(*,1000)

c      ...read data from ktable.info file
       datafile(1:nmaxcharstring)=''
       datafile=trim(data_path)//trim(ktabledir)//'ktable.info'
       open(unit=2,file=datafile,status='old')
       read(2,*)np_ktab
       if(np_ktab.gt.nmaxp_ktab)stop' E- Too many k-table pressures'
       read(2,*)(p_ktab(i),i=1,np_ktab)
       p_ktab(:)=p_ktab(:)*1.0d6                 ! bar      -> dyn cm-2
       read(2,*)nt_ktab
       if(nt_ktab.gt.nmaxt_ktab)stop' E- Too many k-table temperatures'
       read(2,*)(t_ktab(i),i=1,nt_ktab)
       read(2,*)wavemin,wavemax
       read(2,*)
       read(2,*)
       read(2,*)respow
       read(2,*)nwave
       if(nwave.gt.nmaxwave)stop' E- Too many wavenumber intervals'
       read(2,*)nspec_ktab0
       if(nspec_ktab0.gt.nmaxspec_ktab)
     # stop' E- Too many k-table species'
       do i=1,nspec_ktab0
          read(2,*)spec_ktab0(i)
       enddo
       close(2)
       do i=1,np_ktab-1
          if(p_ktab(i+1).le.p_ktab(i))
     #    stop' E- Bad k-table pressures ordering'
       enddo
       do i=1,nt_ktab-1
          if(t_ktab(i+1).le.t_ktab(i))
     #    stop' E- Bad k-table temperatures ordering'
       enddo

c      ...build spectral grid
       call get_spectral_grid

c      ...identify ktable species among included species
       nspec_ktab=0
       idspec_ktab(:)=0
       do i=1,nspec_ktab0
          do k=1,nspec
             if(trim(spec_ktab0(i)).ne.trim(spec(k)))cycle
             nspec_ktab=nspec_ktab+1
             spec_ktab(nspec_ktab)=spec_ktab0(i)
             idspec_ktab(nspec_ktab)=k
             exit
          enddo
       enddo

c      ...read k-coefficidents
       do i=1,nspec_ktab
          do j=1,np_ktab
             kfile='/p00_t00.dat'
             if(j.lt.10)write(kfile(4:4),'(i1)')j
             if(j.ge.10)write(kfile(3:4),'(i2)')j
             do k=1,nt_ktab
                kfile(7:8)='00'
                if(k.lt.10)write(kfile(8:8),'(i1)')k
                if(k.ge.10)write(kfile(7:8),'(i2)')k
                datafile(1:nmaxcharstring)=''
                datafile=trim(data_path)//trim(ktabledir)//
     #          trim(spec_ktab(i))//trim(kfile)
                open(unit=2,file=datafile,err=100,status='old')
                do l=1,nwave
                   read(2,*,err=110)wav,(kcoeff(i,j,k,l,m),m=1,nmaxg)
                   if(wav.ne.wave(l))goto 120
                enddo
                close(2)
             enddo
          enddo
       enddo

       write(*,1010)
       do i=1,nspec_ktab
          write(*,1020)trim(spec_ktab(i))
       enddo

       return

1000   format(2x,'I- Reading k-distribution tables')
1010   format(2x,'I- Opacity species included:')
1020   format(2x,'      ',a)

       return

100    write(*,*)' E- Error opening k-coefficients file'
       write(*,*)kfile
       stop
110    write(*,*)' E- Error reading k-coefficients file'
       write(*,*)kfile
       stop
120    write(*,*)' E- Spectral misalignment in file'
       write(*,*)kfile
       stop

       end




c_______________________________________________________________________

       subroutine read_abundance_profile(string_ab,read_ab)

       include 'pact.common'       
       integer i,j,nxy
       real*8 read_ab(nmaxz),slope
       real*8, allocatable, dimension(:)   :: x,y,xx,yy
       character*(nmaxcharfile) string_ab
       character*100 line

       allocate(x(nmaxxy),y(nmaxxy),xx(nmaxxy),yy(nmaxxy))
       open(unit=3,file=string_ab,status='old')
       i=0
       do
          line(1:100)=''
          read(3,'(a100)',end=10)line
          if(line(1:1).eq.'!')cycle
          i=i+1
          if(i.gt.nmaxxy)stop' E- Too many abundance profile points'
          read(line,*,err=100)xx(i),yy(i)
       enddo
10     close(3)
       nxy=i

       j=0
       do i=2,nxy
          if(xx(i).gt.xx(i-1))then
             j=j+1
          elseif(xx(i).lt.xx(i-1))then
             j=j-1
          endif
       enddo
       if(j.eq.nxy-1)then
          do i=1,nxy
             x(i)=xx(i)*1.0d6                    ! bar -> dyn cm-2
             y(i)=yy(i)
          enddo
       elseif(j.eq.1-nxy)then
          do i=1,nxy
             x(i)=xx(nxy-i+1)*1.0d6              ! bar -> dyn cm-2
             y(i)=yy(nxy-i+1)
          enddo
       else
          write(*,*)trim(string_ab)
          stop' E- Ordering problem in above file'
       endif
       deallocate(xx,yy)

       do j=1,nz
          i=1
          do while(pressure(j).gt.x(i).and.i.lt.nxy)
             i=i+1
          enddo
          i=min(max(2,i),nxy)
          slope=(dlog10(y(i))-dlog10(y(i-1)))/
     #          (dlog10(x(i))-dlog10(x(i-1)))
          read_ab(j)=10.0d0**(dlog10(y(i-1))+
     #    slope*(dlog10(pressure(j))-dlog10(x(i-1))))
       enddo
       deallocate(x,y)

       return

100    write(*,*)trim(string_ab)
       stop' E- Read error in above file'

       end




c_______________________________________________________________________

       subroutine get_initial_temperature_composition

       include 'pact.common'       
       integer i,j
       real*8 zx,px
       character*10000 line
       character*20 dum(3)
       character*(nmaxlenspec), allocatable, dimension(:) :: gspec

       allocate(gspec(nspec))
       open(unit=2,file=initfile,status='old')
       j=nz+2
       do
          line(1:10000)=''
          read(2,'(a10000)',end=10)line
          if(line(1:1).eq.'!')cycle
          j=j-1
          if(j.eq.nz+1)then
             read(line,*,err=100)dum(1),dum(2),dum(3),
     #       (gspec(i),i=1,nspec)
             cycle
          endif
          if(j.lt.1)stop' E- Too many initial z points'
          read(line,*,err=100)zx,px,tk(j),(abun(i,j),i=1,nspec)
          px=px*1.0d6
          if(dabs(px-pressure(j))/pressure(j).gt.1.0d-3)goto 110
          do i=1,nspec
             abun(i,j)=10.0d0**abun(i,j)
          enddo
       enddo
10     close(2)
       do i=1,nspec
          if(trim(gspec(i)).ne.trim(spec(i)))
     #    stop' E- Species names different in initial and species files'
       enddo
       deallocate(gspec)

       call compute_average_particle_mass
       call adjust_vertical_structure

c      ...do convective adjustment
       call compute_heat_capacity
       call convective_adjustment

c      ...load opacity data and star spectrum in case it will be needed
       if(self_temp)then
          call load_cia
          call read_ktable
          call get_star_spectrum
       endif

c      ...write results at time=0
       call write_results

       return

100    stop' E- Read error in initial file'
110    stop' E- Pressure read error in initial file'

       end
