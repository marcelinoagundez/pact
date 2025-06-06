# Makefile for program pact
#
# Once compiled, it is convenient to move the executable to a default OS X path, e.g.:
# > sudo su
# > ...enter password
# > cp pact /usr/local/bin
#
#===============================================================================

 FC      = gfortran
# The option -fallow-argument-mismatch allows to solve a problem with DLSODES documented here:
# https://fortran-lang.discourse.group/t/how-to-fix-the-errors-in-compiling-odepack-in-debug-mode/2597
 FFLAGS  = -O2 -fno-second-underscore -Wall -fallow-argument-mismatch -fbackslash

 BINDIR  = .
 EXEC    = pact

 OBJECTS = pact_data.o \
           pact.o \
           read_input.o \
           compute_chemical_equilibrium.o \
           minimize_gibbs_energy.o \
           mnewt_aux.o \
           load_cia.o \
           compute_initial_temperature.o \
           solve_opacity.o \
           compute_cia.o \
           compute_correlated_k_distribution.o \
           compute_rayleigh.o \
           solve_radiative_transfer_stellar.o \
           solve_radiative_transfer_infrared.o \
           solve_temperature.o \
           solve_continuity_equation.o \
           solve_solid_body_rotation.o \
           get_coefficients.o \
           evaluate_derjac.o \
           write_results.o \
           numrecipes.o \
           dlsodes.o

 COMMON  = pact.common \
           chemical_equilibrium.common \
           dimension.common

all: $(EXEC)
	echo "             Done!"

$(EXEC) : $(OBJECTS) $(COMMON)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(BINDIR)/$(EXEC)

clean:
	rm -f -r *.o *.mod $(EXEC)
