
       subroutine load_cia

       include 'pact.common'
       integer i,j,k
       character*(nmaxlenspec) name
       logical inclu

       i=0

c      ...H2-H2 =======================================================
       i=i+1
c      ...CIA H2-H2: Borysow 2002, A&A, 390, 779
c                    Borysow et al 2001, JQSRT, 68, 235
       ciafile(i)='cia/CIA_H2-H2_60-7000K_10-25000cm-1.dat'
       ntemp_ciafile(i)=20
c      ...CIA H2-H2: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
c       ciafile(i)='cia/CIA_H2-H2_HITRAN.dat'
c       ntemp_ciafile(i)=113
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='H2'
          if(k.eq.2)name(1:2)='H2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...H2-He =======================================================
       i=i+1
c      ...CIA H2-He: Borysow et al 1989, ApJ, 336, 495
c                    Borysow & Frommhold 1989, ApJ, 341, 549
c                    Borysow et al 1997, A&A, 324, 185
       ciafile(i)='cia/CIA_H2-He_100-7000K_10-25000cm-1.dat'
       ntemp_ciafile(i)=16
c      ...CIA H2-He: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
c       ciafile(i)='cia/CIA_H2-He_HITRAN.dat'
c       ntemp_ciafile(i)=334
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='H2'
          if(k.eq.2)name(1:2)='He'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...H2-H ========================================================
       i=i+1
c      ...CIA H2-H: HITRAN
c                   Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_H2-H_HITRAN.dat'
       ntemp_ciafile(i)=4
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='H2'
          if(k.eq.2)name(1:2)='H'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...He-H ========================================================
       i=i+1
c      ...CIA He-H: HITRAN
c                   Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_He-H_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='He'
          if(k.eq.2)name(1:2)='H'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...H2-CH4 ======================================================
       i=i+1
c      ...CIA H2-CH4: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_H2-CH4_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='H2'
          if(k.eq.2)name(1:2)='CH4'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...N2-CH4 ======================================================
       i=i+1
c      ...CIA N2-CH4: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_N2-CH4_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='N2'
          if(k.eq.2)name(1:2)='CH4'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...N2-H2 =======================================================
       i=i+1
c      ...CIA N2-H2: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_N2-H2_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='N2'
          if(k.eq.2)name(1:2)='H2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...N2-H2O ======================================================
       i=i+1
c      ...CIA N2-H2O: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_N2-H2O_HITRAN.dat'
       ntemp_ciafile(i)=11
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='N2'
          if(k.eq.2)name(1:2)='H2O'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...N2-He =======================================================
       i=i+1
c      ...CIA N2-He: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_N2-He_HITRAN.dat'
       ntemp_ciafile(i)=1
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='N2'
          if(k.eq.2)name(1:2)='He'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...N2-N2 =======================================================
       i=i+1
c      ...CIA N2-N2: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_N2-N2_HITRAN.dat'
       ntemp_ciafile(i)=45
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='N2'
          if(k.eq.2)name(1:2)='N2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...O2-CO2 ======================================================
       i=i+1
c      ...CIA O2-CO2: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_O2-CO2_HITRAN.dat'
       ntemp_ciafile(i)=1
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='O2'
          if(k.eq.2)name(1:2)='CO2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...O2-N2 =======================================================
       i=i+1
c      ...CIA O2-N2: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_O2-N2_HITRAN.dat'
       ntemp_ciafile(i)=17
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='O2'
          if(k.eq.2)name(1:2)='N2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...O2-O2 =======================================================
       i=i+1
c      ...CIA O2-O2: HITRAN
c                    Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_O2-O2_HITRAN.dat'
       ntemp_ciafile(i)=21
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='O2'
          if(k.eq.2)name(1:2)='O2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CO2-Ar ======================================================
       i=i+1
c      ...CIA CO2-Ar: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CO2-Ar_HITRAN.dat'
       ntemp_ciafile(i)=21
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CO2'
          if(k.eq.2)name(1:2)='Ar'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CO2-CH4 =====================================================
       i=i+1
c      ...CIA CO2-CH4: HITRAN
c                      Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CO2-CH4_HITRAN.dat'
       ntemp_ciafile(i)=4
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CO2'
          if(k.eq.2)name(1:2)='CH4'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CO2-CO2 =====================================================
       i=i+1
c      ...CIA CO2-CO2: HITRAN
c                      Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CO2-CO2_HITRAN.dat'
       ntemp_ciafile(i)=14
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CO2'
          if(k.eq.2)name(1:2)='CO2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CO2-H2 ======================================================
       i=i+1
c      ...CIA CO2-H2: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CO2-H2_HITRAN.dat'
       ntemp_ciafile(i)=4
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CO2'
          if(k.eq.2)name(1:2)='H2'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CO2-He ======================================================
       i=i+1
c      ...CIA CO2-He: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CO2-He_HITRAN.dat'
       ntemp_ciafile(i)=1
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CO2'
          if(k.eq.2)name(1:2)='He'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CH4-Ar ======================================================
       i=i+1
c      ...CIA CH4-Ar: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CH4-Ar_HITRAN.dat'
       ntemp_ciafile(i)=5
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CH4'
          if(k.eq.2)name(1:2)='Ar'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CH4-CH4 =====================================================
       i=i+1
c      ...CIA CH4-CH4: HITRAN
c                      Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CH4-CH4_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CH4'
          if(k.eq.2)name(1:2)='CH4'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

c      ...CH4-He ======================================================
       i=i+1
c      ...CIA CH4-He: HITRAN
c                     Karman et al 2019, Icarus, 160, 175
       ciafile(i)='cia/CIA_CH4-He_HITRAN.dat'
       ntemp_ciafile(i)=10
       do k=1,2
          name(1:nmaxlenspec)=''
          if(k.eq.1)name(1:2)='CH4'
          if(k.eq.2)name(1:2)='He'
          inclu=.false.
          do j=1,nspec
             if(trim(spec(j)).eq.trim(name))then
                inclu=.true.
                exit
             endif
          enddo
          if(.not.inclu)then
             i=i-1
             exit
          endif
          if(k.eq.1)idspec1_cia(i)=j
          if(k.eq.2)idspec2_cia(i)=j
       enddo

       ncia=i
       if(ncia.gt.nmaxcia)stop' E- Too many CIA data sets'

       return

       end
