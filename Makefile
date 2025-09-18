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
	*.md \
	inputs/example/*.in \
	include/

MODDIR := include

#FLAGS=
#NOTE: For Intel Fortran compiler ifort
ifeq ($(SYSTEM),IFORT)
	FLAGS=  -O3 -r8 -double-size 128
	COMP=   ifort
endif

ifeq ($(SYSTEM),GFORT)
#NOTE	: For gfortran
#gfortran orders roots a bit differently...
	FLAGS=  -O3 -DDOUBLE -fdefault-real-8 -funroll-loops -ffast-math  -I$(MODDIR) -J$(MODDIR) #
	COMP= gfortran
endif
LIBS=

LFMOD=	nrtype.o nrutil_trim.o vars.o functions.o bessel.o disprels.o fpc.o

LFX=  	plume.o

VPATH= src:include:/usr/include/

###############################################################################
all: heat

heat: $(MODDIR) $(LFMOD) $(LFX)
	$(COMP) -o plume.e $(FLAGS) $(LIBS) $(LFMOD) $(LFX)

###############################################################################

###############################################################################

tidyup:
	mv *.o include/
	mv *.mod include/

clean:
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

bessel.o:      nrtype.o nrutil_trim.o

disprels.o:    vars.o bessel.o

functions.o:   vars.o

fpc.o:         disprels.o vars.o

plume.o:       vars.o functions.o fpc.o