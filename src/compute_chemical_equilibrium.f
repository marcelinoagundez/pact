
       subroutine compute_chemical_equilibrium
       
       include 'pact.common'       
       include 'chemical_equilibrium.common'       
       integer i,j,k,l,js,iz,iwk,nspec_save,ncondensed_save,
     # idcondensed(nmaxcondensed),ibad,ntflag(nmaxspec)
       real*8 rwk,dg,dgmin,nceatot_save
       logical usd(nmaxcondensed),bad(nmaxcondensed)

       write(*,1000)
       iwk=0                           ! to avoid warning
       ntflag(:)=0

c      ...compute terms b0(i) of initial elemental abundances (+charge)
       rwk=0.0d0
       do i=1,nelem
          rwk=rwk+elfab(i)*mat(i)*amu
       enddo
       nceatot=0.0d0
       b0max=0.0d0
       b0min=huge
       do i=1,nelem
          b0(i)=elfab(i)/navg/rwk      ! [=] mole [g of mixture]-1
          nceatot=nceatot+b0(i)        ! n_tot(gas) for a fully atomic gas
          if(b0(i).gt.b0max)b0max=b0(i)
          if(b0(i).lt.b0min)b0min=b0(i)
       enddo
       nceatot_save=nceatot
       if(ion)then                     ! charge
          b0(nelem+1)=0.0d0
          do j=1,nspec                 ! add charge to nat(i,j) variable
             nat(j,nelem+1)=charge(j)
          enddo
       endif

c      ...compute the maximum abundance reachable by each species
       do j=1,nspec
          nceamax(j)=huge
          do i=1,nelem
             iwk=nat(j,i)
             if(iwk.ne.0.and.nceamax(j).gt.(b0(i)/dble(iwk)))
     #       nceamax(j)=b0(i)/dble(iwk)
          enddo
       enddo
       if(ion.and.id_electr.ne.0)nceamax(id_electr)=nceatot

c      ...assign initial guess for species abundances
       do j=1,nspec
          ncea(j)=nceatot/(dble(ngas))
          if(condensed(j))ncea(j)=0.0d0
          if(ncea(j).gt.nceamax(j))ncea(j)=nceamax(j)
       enddo

c Loop over (p,T) grid
       do iz=1,nz                                ! from top to bottom
          pce=pressure(iz)/1.0d6                 ! dyn cm-2 -> bar
          tce=tk(iz)

c ...if only gas phase species
          if(ncondensed.eq.0)then
             call minimize_gibbs_energy
             if(.not.converge_gibbs)then
                if(iz.eq.1)then
                   nceatot=nceatot_save
                   call initialize_ncea
                else
                   stop' E - Convergence not reached'
                endif
             endif
             do i=1,nspec                        ! save abundances
                abun(i,iz)=max(ncea(i)/nceatot,abunmin)
             enddo
             cycle
          endif

c ...if condensed species are included
          nspec_save=nspec                       ! save condensed species data on js=nspec+1,...
          ncondensed_save=ncondensed
          call compute_chemical_potential(pce,tce)
          do i=1,ncondensed_save
             j=ngas+i
             js=ngas+ncondensed_save+i
             call assign_species(j,js)
             ncea(js)=0.0d0
          enddo
          nspec=ngas                             ! solve gas phase chemical equilibrium
          ncondensed=0
          call minimize_gibbs_energy
          if(.not.converge_gibbs)stop' E - Convergence not reached'

c      ...loop: one calculation per new condensed species
          usd(1:nmaxcondensed)=.false.           ! add (one by one) condensed species and
          bad(1:nmaxcondensed)=.false.           ! solve gas+condensed chemical equilibrium
          do l=1,ncondensed_save

             dgmin=huge                          ! look for new condensed species that reduces Gibbs energy
             do i=1,ncondensed_save
                js=ngas+ncondensed_save+i
                if(usd(i).or.bad(i))cycle
                if(tflag(js))cycle
                rwk=0.0d0
                do k=1,npilag
                   rwk=rwk+nat(js,k)*pilag(k)
                enddo
                dg=mupot(js)-rwk
                if(dg.lt.dgmin)then
                   dgmin=dg
                   iwk=i
                endif
             enddo
             if(dgmin.ge.0.0d0)exit

             usd(iwk)=.true.                      ! if condensed species reduces Gibbs energy
10           ncondensed=0                         ! include it and solve gas+condensed chemical equilibrium
             do i=1,ncondensed_save
                if(.not.usd(i).or.bad(i))cycle
                ncondensed=ncondensed+1
                idcondensed(ncondensed)=ngas+ncondensed_save+i
                j=ngas+ncondensed
                js=idcondensed(ncondensed)
                call assign_species(js,j)
                ncea(j)=ncea(js)
             enddo
             nspec=ngas+ncondensed
             call minimize_gibbs_energy
             if(.not.converge_gibbs)stop' E - Convergence not reached'
             do i=1,ncondensed                    ! assign computed abundances of condensed species
                ncea(idcondensed(i))=ncea(ngas+i) ! to js=nspec+1,...
             enddo
             ibad=0
             do i=1,ncondensed
                if(ncea(idcondensed(i)).ge.0.0d0)cycle
                ibad=ibad+1
                bad(idcondensed(i)-ngas-ncondensed_save)=.true.
                ncea(idcondensed(i))=0.0d0
             enddo
             if(ibad.gt.0)goto 10

c      ...end of loop: one calculation per new condensed species
          enddo

          do i=1,ncondensed_save                 ! assign back condensed species data and abundances
             j=ngas+i
             js=ngas+ncondensed_save+i
             call assign_species(js,j)
             ncea(j)=ncea(js)
          enddo
          nspec=nspec_save
          ncondensed=ncondensed_save

          do i=1,nspec                           ! save abundances
             abun(i,iz)=max(ncea(i)/nceatot,abunmin)
             if(tflag(i))ntflag(i)=ntflag(i)+1
          enddo

c End of loop over (p,T) grid
       enddo

       do i=1,nspec
          if(ntflag(i).eq.0)cycle
          write(*,1010)trim(spec(i))
       enddo

c      ...compute average particle mass in planet's atmosphere: apmass
       call compute_average_particle_mass

       return

1000   format(2x,'I- Solving chemical equilibrium')
1010   format(2x,'W- NASA polynomials extrapolated in T for: ',a)

       end




c_______________________________________________________________________

       subroutine assign_species(ji,jo)

       include 'pact.common'       
       include 'chemical_equilibrium.common'       
       integer ji,jo

c      ...Intrinsic data
       spec(jo)(1:nmaxlenspec)=spec(ji)(1:nmaxlenspec)
       condensed(jo)=condensed(ji)
       nattot(jo)=nattot(ji)
       nat(jo,:)=nat(ji,:)
       charge(jo)=charge(ji)
       type_therm(jo)=type_therm(ji)
       natherm(jo)=natherm(ji)
       ntemp_therm(jo)=ntemp_therm(ji)
       temp_therm(jo,:)=temp_therm(ji,:)
       atherm(jo,:,:)=atherm(ji,:,:)

c      ...T-dependent data (also p-dependent for gas species)
       tflag(jo)=tflag(ji)
       mupot(jo)=mupot(ji)

       return
       end




c_______________________________________________________________________

       subroutine initialize_ncea

       include 'pact.common'       
       include 'chemical_equilibrium.common'       
       integer j

c      ...assign initial guess for species abundances
       do j=1,nspec
          ncea(j)=nceatot/(dble(ngas))
          if(ncea(j).gt.nceamax(j))ncea(j)=nceamax(j)
       enddo

c      ...initialize ncea values at 0.1 mbar and 2000 K
       pce=1.0d-4                      ! [bar]
       tce=2000.0d0                    ! [K]
       call minimize_gibbs_energy

c      ...minimize Gibbs energy at first (p,T) grid point
       pce=pressure(1)/1.0d6           ! dyn cm-2 -> bar
       tce=tk(1)
       call minimize_gibbs_energy
       if(.not.converge_gibbs)stop' E- Gibbs minimization not converged'

       return
       end
