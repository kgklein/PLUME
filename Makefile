######################################################################## OVERVIEW
#  Makefile for PLUME / JET-PLUME:
#  (Judging Energy Transfer in a) Plasma in a Linear Uniform Magnetic Environment
#
#   VERSION 1.1
#
#  Version notes:
#
#  LAST UPDATE:  2025/08/18
###############################################################################


# SYSTEM=IFORT
#IFORT licenses are hard to obtain, let's default to GFORT
SYSTEM=GFORT

PACK = Makefile \
	src/*.f90 \
	src/*.inc \
	src/*.c \
	third_party/faddeeva/* \
	*.md \
	inputs/example/*.in \
	include/

MODDIR := include

#FLAGS=
#NOTE: For Intel Fortran compiler ifort
ifeq ($(SYSTEM),IFORT)
	FLAGS=  -O3 -r8 -double-size 128
	COMP=   ifort
	COLD_FLAGS= -fp-model precise
endif

ifeq ($(SYSTEM),GFORT)
#NOTE	: For gfortran
#gfortran orders roots a bit differently...
	FLAGS=  -O3 -DDOUBLE -fdefault-real-8 -funroll-loops -ffast-math  -I$(MODDIR) -J$(MODDIR) #
	COMP= gfortran
	COLD_FLAGS= -fno-fast-math
endif
CC ?= cc
FADDEEVA_CFLAGS = -O3 -std=c99 -fno-fast-math
LIBS= -lm

LFMOD=	faddeeva_bridge.o faddeeva.o plasma_z_faddeeva.o nrtype.o nrutil_trim.o vars.o plasma_z_legacy.o plasma_z_backend.o cold_susceptibility.o auxiliary_normalization.o kinetic_core.o harmonic_selection.o bessel_cache.o adaptive_kinetic.o cold_reference.o species_moments.o electromagnetic_fields.o heating_response.o heating_energy.o heating_diagnostics.o scan_support.o functions.o bessel.o disprels.o fpc.o

LFX=  	plume.o

VPATH= src:include:/usr/include/

###############################################################################
all: heat

heat: $(MODDIR) $(LFMOD) $(LFX)
	$(COMP) -o plume.e $(FLAGS) $(LFMOD) $(LFX) $(LIBS)

###############################################################################

###############################################################################

tidyup:
	mv *.o include/
	mv *.mod include/

clean:
	rm -f $(LFMOD) $(LFX)
	rm -f include/*.o
	rm -f include/*.mod
	rm -f plume*.e

tar: 
	tar -cvf  pack_plume_`date +'%y%m%d'`.tar ${PACK}; gzip pack_plume_`date +'%y%m%d'`.tar

#########Rules
%.o : %.f90 | $(MODDIR)
	$(COMP) -c $(FLAGS) $<

#########Dependencies
nrutil_trim.o: nrtype.o

vars.o:        nrtype.o

# Preserve IEEE validation only in the standalone cold kernel. Kinetic flags
# remain unchanged; in particular, -ffast-math must not remove finite checks.
cold_susceptibility.o: src/cold_susceptibility.f90 vars.o nrtype.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

bessel.o:      nrtype.o nrutil_trim.o

disprels.o:    functions.o heating_diagnostics.o electromagnetic_fields.o scan_support.o plasma_z_backend.o vars.o bessel.o cold_susceptibility.o cold_reference.o src/kinetic_species.inc src/kinetic_collisional_species.inc src/kinetic_harmonics.inc src/electromagnetic_electric.inc src/kinetic_parallel_response.inc src/kinetic_parallel_longitudinal.inc

functions.o:   heating_diagnostics.o vars.o adaptive_kinetic.o scan_support.o

fpc.o:         disprels.o vars.o

plume.o:       vars.o functions.o fpc.o

# Compile the explicit auxiliary core separately; preserve IEEE failure checks.
auxiliary_normalization.o: src/auxiliary_normalization.f90 vars.o nrtype.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

kinetic_core.o: src/kinetic_core.f90 vars.o nrtype.o src/kinetic_species.inc src/kinetic_harmonics.inc | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

cold_reference.o: auxiliary_normalization.o kinetic_core.o adaptive_kinetic.o vars.o

# Selection/cache checks require IEEE semantics, including in optimized builds.
harmonic_selection.o: src/harmonic_selection.f90 | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<
bessel_cache.o: src/bessel_cache.f90 harmonic_selection.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<
adaptive_kinetic.o: src/adaptive_kinetic.f90 vars.o kinetic_core.o bessel_cache.o src/kinetic_species.inc src/kinetic_collisional_species.inc | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

# Johnson's dual C/C++ source is compiled as ordinary C99; upstream is unedited.
faddeeva.o: third_party/faddeeva/Faddeeva.cc third_party/faddeeva/Faddeeva.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(FADDEEVA_CFLAGS) -x c -c $< -o $@
plasma_z_faddeeva.o: src/plasma_z_faddeeva.f90 | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

faddeeva_bridge.o: src/faddeeva_bridge.c third_party/faddeeva/Faddeeva.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(FADDEEVA_CFLAGS) -Ithird_party/faddeeva -c $< -o $@

$(MODDIR):
	mkdir -p $@

plasma_z_backend.o: src/plasma_z_backend.f90 vars.o plasma_z_faddeeva.o plasma_z_legacy.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

scan_support.o: src/scan_support.f90 vars.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

electromagnetic_fields.o: src/electromagnetic_fields.f90 src/electromagnetic_electric.inc vars.o species_moments.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

species_moments.o: src/species_moments.f90 vars.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<

heating_response.o: src/heating_response.f90 cold_reference.o species_moments.o cold_susceptibility.o src/kinetic_parallel_response.inc src/kinetic_parallel_longitudinal.inc | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<
heating_energy.o: src/heating_energy.f90 heating_response.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<
heating_diagnostics.o: src/heating_diagnostics.f90 heating_energy.o electromagnetic_fields.o | $(MODDIR)
	$(COMP) -c $(FLAGS) $(COLD_FLAGS) $<
