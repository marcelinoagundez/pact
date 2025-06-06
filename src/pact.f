
       include 'pact.common'
       integer i

       write(*,1000)

       open(unit=1,err=100,file='pact.inp',status='old')
       do i=1,nmaxmodel
          call reset_variables
          read(1,*,err=110,end=10)inpfile
          if(inpfile(1:3).eq.'END'.or.inpfile(1:3).eq.'end')exit

c      ...read input parameters
          call read_input

c      ...get pressure grid
          call get_pressure_grid

c      ...read species and NASA polynomials
          call read_species

c      ...compute/get temperature (under chemical equilibrium or not)
          if(solvitemp)then
             call compute_initial_temperature
          else
             call get_initial_temperature_composition
          endif

c      ...read data to solve non-equilibrium chemical composition
          call read_reactions
          call read_eddy
          call read_photon

c      ...solve non-equilibrium chemical composition
          if(.not.sbrotation)then
             call solve_continuity_equation
          else
             call solve_solid_body_rotation
          endif
       enddo
10     close(1)

       write(*,1010)
       stop

100    write(*,*)' E- Error opening pact.inp'
       stop
110    write(*,*)' E- Error reading pact.inp'
       stop

1000   format(/,
     # '---------------------------------------------------------',/,
     # '                      PACT program                       ',/,
     # '    1D Planetary Atmosphere Chemistry and Temperature    ',/,
     # '                    version March 2025                   ',/,
     # '---------------------------------------------------------',/)
1010   format(/,
     # '---------------------------------------------------------',/,
     # '                     Calculations done                   ',/,
     # '                                                         ',/,
     # '                       Exiting PACT                      ',/,
     # '---------------------------------------------------------',/)

       end




c_______________________________________________________________________

       subroutine reset_variables
       
       include 'pact.common'       

c      ...reset CPU time
       call cpu_time(cpu0)

c      ...model parameters
       ion=.false.
       solvitemp=.false.
       self_temp=.false.
       diffusion=.false.
       photochem=.false.
       write_supp=.true.
       initial=.true.

c      ...physical parameters
       nz=0
       mplanet=0.0d0
       rplanet=0.0d0
       pplanet=0.0d0
       tintplanet=0.0d0
       albedoplanet=0.0d0
       emplanet=0.0d0
       fhrd=0.0d0
       tplanet=0.0d0
       tplanet_th=0.0d0
       dist=0.0d0
       zenithangle=0.0d0
       rstar=0.0d0
       tstar=0.0d0
       pmax=0.0d0
       pmin=0.0d0
       pressure(:)=0.0d0
       z(:)=0.0d0
       tk(:)=0.0d0
       density(:)=0.0d0
       deltaz(:)=0.0d0
       ggrav(:)=0.0d0
       plev(:)=0.0d0
       apmass=0.0d0
       cplayer(:)=0.0d0
       convect(:)=.false.
       detached_convect=.false.

c      ...elements parameters
       nelem=0
       elfab(:)=0.0d0
       mat(:)=0.0d0
       elem(:)(1:2)=''

c      ...species parameters
       nspec=0
       spec(:)(1:nmaxlenspec)=''
       charge(:)=0
       nattot(:)=0
       nat(:,:)=0
       mspec(:)=0.0d0
c       polary(:)=0.0d0
       nelemspec(:)=0
       abun(:,:)=0.0d0
       id_electr=0
       ngas=0
       ncondensed=0
       condensed(:)=.false.
       iabun=0

c      ...NASA polynomial parameters
       type_therm(:)=0
       natherm(:)=0
       ntemp_therm(:)=0
       atherm(:,:,:)=0.0d0
       temp_therm(:,:)=0.0d0

c      ...Radiative transfer parameters
       wavemin=0.0d0
       wavemax=0.0d0
       respow=0.0d0
       nwave=0
       wave(:)=0.0d0
       wave_l(:)=0.0d0
       wave_u(:)=0.0d0
       dwave(:)=0.0d0
       fstar(:)=0.0d0
       kabs(:,:)=0.0d0
       kscat(:,:)=0.0d0
       walb(:,:)=0.0d0
       gscat(:,:)=0.0d0
       fup_st(:,:)=0.0d0
       fdo_st(:,:)=0.0d0
       fup_ir(:,:)=0.0d0
       fdo_ir(:,:)=0.0d0
       fup(:)=0.0d0
       fdown(:)=0.0d0
       fnet(:)=0.0d0

c      ...CIA parameters
       ncia=0
       ciafile(:)(1:nmaxcharstring)=''
       ntemp_ciafile(:)=0
       idspec1_cia(:)=0
       idspec2_cia(:)=0

c      ...k-table parameters
       np_ktab=0
       p_ktab(:)=0.0d0
       nt_ktab=0
       t_ktab(:)=0.0d0
       nspec_ktab=0
       spec_ktab(:)(1:nmaxlenspec)=''
       idspec_ktab(:)=0
       kcoeff(:,:,:,:,:)=0.0d0

c      ...temperature computation
       bigit=0
       dtkmax=0.0d0
       conv=0.0d0

c      ...reaction parameters
       nreac=0
       reac(:)(1:nmaxlenreac)=''
       alpha(:)=0.0d0
       beta(:)=0.0d0
       gamma(:)=0.0d0
       rtype(:)=0
       tkmin(:)=0.0d0
       tkmax(:)=0.0d0
       kerror(:)=0.0d0
       comment(:)=''
       nreag(:)=0
       reag(:,:)(1:nmaxlenspec)=''
       idreag(:,:)=0
       nprod(:)=0
       prod(:,:)(1:nmaxlenspec)=''
       idprod(:,:)=0
       natc(:)=0
       idreac_rev(:)=0
       krate(:,:)=0.0d0

c      ...collisional reaction parameters
       nmreac=0
       idmreac(:)=0
       alphainf(:)=0.0d0
       betainf(:)=0.0d0
       gammainf(:)=0.0d0
       afcent(:)=0.0d0
       bfcent(:)=0.0d0
       mcoll(:)=''
       tkmininf(:)=0.0d0
       tkmaxinf(:)=0.0d0
       tkminfc(:)=0.0d0
       tkmaxfc(:)=0.0d0

c      ...chemical production and loss parameters
       nterm(:)=0
       idterm(:,:)=0
       sterm(:,:)=0
       fterm(:,:)=0
       nphototerm(:)=0
       idphototerm(:,:)=0
       sphototerm(:,:)=0
       fphototerm(:,:)=0

c      ...diffusion parameters
       kdif(:)=0.0d0
       ddif(:,:)=0.0d0

c      ...photo reaction parameters
       nwlth=0
       wlthmin=0.0d0
       wlthmax=0.0d0
       wlth(:)=0.0d0
       fphot(:)=0.0d0
       nphotoreac=0
       nphotospec=0
       idphotospec(:)=0
       photoreac(:)(1:50)=''
       photoreag(:)(1:nmaxlenspec)=''
       idphotoreag(:)=0
       nphotoprod(:)=0
       photoprod(:,:)(1:nmaxlenspec)=''
       idphotoprod(:,:)=0
       csphot(:,:)=0.0d0
       csphotabs(:,:)=0.0d0
       jphot(:,:)=0.0d0
       p_tauuv1(:,:)=0.0d0
       p_tauir1(:,:)=0.0d0
       jerror(:)=0.0d0

c      ...continuity equation parameters
       time=0.0d0
       steady=.false.
       varabun=0.0d0
       itime=0
       f_conteq(:,:)=0.0d0
       df_conteq(:,:,:)=0.0d0
       df_dif(:,:)=0.0d0
       df_jm_dif(:,:)=0.0d0
       df_jp_dif(:,:)=0.0d0
       cdelem0(:)=0.0d0
       cdtot0=0.0d0
       elab0(:,:)=0.0d0

c      ...solid body rotation data
       sbrotation=.false.
       windspeed=0.0d0
       period=0.0d0
       nperiod=0
       nlong=0
       long(:)=0.0d0
       tk_grid(:,:)=0.0d0
       longitude=0.0d0
       tin_sbr=0.0d0
       ostat_sbr=0

c      ...input files
       ncharmodel=0
       inpfile(1:nmaxcharfile)=''
       namemodel(1:nmaxcharfile)=''
       specfile(1:nmaxcharfile)=''
       reacfile(1:nmaxcharfile)=''
       starfile(1:nmaxcharfile)=''
       thermfile(1:nmaxcharfile)=''
       eddyfile(1:nmaxcharfile)=''
       initfile(1:nmaxcharfile)=''
       lptfile(1:nmaxcharfile)=''
       data_path(1:nmaxcharstring)=''
       ktabledir(1:nmaxcharstring)=''
       chemdir(1:nmaxcharstring)=''
       photodir(1:nmaxcharstring)=''
       stardir(1:nmaxcharstring)=''

       return
       end




c_______________________________________________________________________

       subroutine get_cputime(ihou,imin,isec)
       
       include 'pact.common'       
       integer cputime,ihou,imin,isec

       call cpu_time(cpu)
       cputime=int(cpu-cpu0)
       ihou=cputime/3600
       imin=cputime/60-ihou*60
       isec=mod(cputime,60)

       return
       end




c_______________________________________________________________________

       real*8 function bnuplanck(temp,nu)
c      computes the Planck's law:
c      input:   temp: temperature [K]
c               nu:   frequency   [Hz]
c      returns: bnuplanck: specific intensity [erg s-1 cm-2 Hz-1 sr-1]

       implicit none
       real*8 hplanck,kboltz,clight
       parameter(hplanck=6.62607015d-27) ! Planck constant              [erg s]
       parameter(kboltz=1.380649d-16)    ! Boltzmann constant           [erg/K]
       parameter(clight=2.99792458d10)   ! speed of light               [cm/s]
       real*8 temp,nu,hnukt

       hnukt=hplanck*nu/kboltz/temp
       if(hnukt.gt.250)hnukt=250.d0
       bnuplanck=(2*hplanck*nu*nu*nu)/(clight*clight)/
     #           (dexp(hnukt)-1.0d0)

       end
