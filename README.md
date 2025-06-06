### PACT (Planetary Atmosphere Chemistry and Temperature) #

  Fortran code that computes self-consistently the 1D vertical distribution of
   - temperature                         (through radiative-convective equilibrium)
   - disequilibrium chemical composition (including thermochemical kinetics, photochemistry, and vertical mixing)
  in a planetary atmosphere

### If you use PACT for your work, I would appreciate if you cite the paper
   - Agúndez 2025, A&A, in press

### Inquiries can be sent to marcelino.agundez@csic.es

### Compilation:

  Download and uncompress the zip file. This will create two directories:

   - `src`
   - `examples`

  You need to have gfortran installed in your system. Compilation has been tested with gfortran.
  Other fortran compilers such as ifort may work but they have not been tested.
  To compile simply type:

`$ cd src`

`$ make clean`       (only needed if the code has been compiled previously)

`$ make`

  This should have created a binary file named `pact`.
  It is convenient to move the binary file to the bin folder of your system where binaries are located.
  This allows you to run the code simply typing `pact` in any folder without the need to duplicate binaries.

### Usage:

  Download the directory `data` as a compressed file from [here](https://saco.csic.es/s/TsR67qdyPziEzXL), uncompress it, and put it in the `pact-main` folder.
  In the directory `data` you have sub-directories with different type of data.
  Most of these data are common to different PACT models and thus are stored once to avoid file duplication.

   - `chemistry` : reaction network and thermochemical data
   - `cia`       : collision-induced-absorption data
   - `ktable`    : pre-computed k tables (correlated k distribution) to compute temperature
   - `photoreac` : UV cross section data to compute photochemistry
   - `stars`     : stellar spectra

  In the directory `examples` you have various sub-directories with the names of different planets.
  These sub-directories contain specific input files for each particular PACT model.
  To run a PACT model of a given planet, move to corresponding sub-directory and simply type:

`$ pact`             (assuming you have the binary in the bin folder of your system)

