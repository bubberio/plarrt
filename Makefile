# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
#configuration
FC := gfortran #ifort
LD := $(FC)
RM := rm -f
FLAG := -J src/mod
OPTS= 
#########################################################
# executable names and corresponding sources
# The goal is to make this text file as short as possible
SRCph = constants.f90 SHcomputer.f90 geometry.f90 trinomials.f90 plarrt_poly_handler.f90

SRCprop= constants.f90 SHcomputer.f90 geometry.f90 trinomials.f90 keplerian.f90 integration.f90 prop_lite.f90 spherigrav.f90 polygrav.f90 dynamics.f90 plarrt_prop.f90

SRCcomp= constants.f90 SHcomputer.f90 geometry.f90 trinomials.f90 keplerian.f90 integration.f90 prop_lite.f90 spherigrav.f90 polygrav.f90 dynamics.f90 plarrt_compare.f90

SRCfli= constants.f90 SHcomputer.f90 geometry.f90  trinomials.f90 SHcomputer.f90 keplerian.f90 integration.f90 prop_lite.f90 var_lite.f90 spherigrav.f90 polygrav.f90 dynamics.f90 plarrt_fli.f90
#########################################################
EXEPOLYHAN=plarrt_poly.x

EXEPROP=plarrt_prop.x

EXECOMP=plarrt_compare.x

EXEFLI=plarrt_fli.x
#########################################################
# list of object files with /obj prefix
OBJ_NOPREFIXPH=$(SRCph:.f90=.o)
OBJPH=$(addprefix obj/,$(OBJ_NOPREFIXPH))

OBJ_NOPREFIXPROP=$(SRCprop:.f90=.o)
OBJPROP=$(addprefix obj/,$(OBJ_NOPREFIXPROP))

OBJ_NOPREFIXCOMP=$(SRCcomp:.f90=.o)
OBJCOMP=$(addprefix obj/,$(OBJ_NOPREFIXCOMP))

OBJ_NOPREFIXFLI=$(SRCfli:.f90=.o)
OBJFLI=$(addprefix obj/,$(OBJ_NOPREFIXFLI))
################################################################
# creation of executable programs
ph: $(EXEPOLYHAN)
$(EXEPOLYHAN): $(OBJPH)
	$(FC) $(FLAG) -fopenmp $^ -o $@ $(OPTS)

prop: $(EXEPROP)
$(EXEPROP): $(OBJPROP)
	$(FC) $(FLAG) $^ -o $@ $(OPTS)

comp: $(EXECOMP)
$(EXECOMP): $(OBJCOMP)
	$(FC) $(FLAG) $^ -o $@ $(OPTS)

fli: $(EXEFLI)
$(EXEFLI): $(OBJFLI)
	$(FC) $(FLAG) $^ -o $@ $(OPTS)

#########################################################
# compilation of sources

obj/plarrt_fli.o: src/plarrt_fli.f90 obj/constants.o obj/geometry.o obj/trinomials.o obj/keplerian.o obj/integration.o obj/prop_lite.o obj/var_lite.o obj/SHcomputer.o
	$(FC) $(FLAG) -c $< -o $@

obj/plarrt_compare.o: src/plarrt_compare.f90 obj/constants.o obj/geometry.o obj/trinomials.o obj/keplerian.o obj/integration.o obj/prop_lite.o 
	$(FC) $(FLAG) -c $< -o $@

obj/plarrt_prop.o: src/plarrt_prop.f90 obj/constants.o obj/geometry.o obj/trinomials.o obj/keplerian.o obj/integration.o obj/prop_lite.o 
	$(FC) $(FLAG) -c $< -o $@

obj/plarrt_poly_handler.o: src/plarrt_poly_handler.f90 obj/constants.o obj/geometry.o obj/trinomials.o obj/SHcomputer.o
	$(FC) $(FLAG) -fopenmp -c $< -o $@

obj/var_lite.o: src/var_lite.f90 obj/constants.o obj/SHcomputer.o obj/geometry.o obj/trinomials.o obj/dynamics.o obj/polygrav.o obj/keplerian.o obj/integration.o obj/spherigrav.o obj/polygrav.o
	$(FC) $(FLAG) -c $< -o $@

obj/prop_lite.o: src/prop_lite.f90 obj/constants.o obj/SHcomputer.o obj/geometry.o obj/trinomials.o obj/dynamics.o obj/polygrav.o obj/keplerian.o obj/integration.o obj/spherigrav.o obj/polygrav.o
	$(FC) $(FLAG) -c $< -o $@

obj/dynamics.o: src/dynamics.f90 obj/constants.o obj/geometry.o obj/spherigrav.o obj/polygrav.o
	$(FC) $(FLAG) -c $< -o $@

obj/polygrav.o: src/polygrav.f90 obj/constants.o obj/trinomials.o obj/geometry.o
	$(FC) $(FLAG) -c $< -o $@

obj/spherigrav.o: src/spherigrav.f90 obj/constants.o obj/geometry.o
	$(FC) $(FLAG) -c $< -o $@

obj/integration.o: src/bottom/integration.f90 obj/constants.o
	$(FC) $(FLAG) -c $< -o $@

obj/keplerian.o: src/keplerian.f90 obj/constants.o obj/geometry.o 
	$(FC) $(FLAG) -c $< -o $@

obj/SHcomputer.o: src/SHcomputer.f90 obj/trinomials.o obj/constants.o obj/geometry.o
	$(FC) $(FLAG) -fopenmp -c $< -o $@

obj/trinomials.o: src/bottom/trinomials.f90
	$(FC) $(FLAG) -c $< -o $@

obj/geometry.o: src/bottom/geometry.f90
	$(FC)  $(FLAG) -c $< -o $@

obj/constants.o: src/bottom/constants.f90
	$(FC) $(FLAG) -c $< -o $@

#########################################################
# suppression of compilation/output files

clean:
	rm -f *.x obj/*.o src/mod/*.mod *.mod

cleanall:
	rm -f *.x obj/*.o mod/*.mod out/*

#########################################################