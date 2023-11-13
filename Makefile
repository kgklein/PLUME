######################################################################## OVERVIEW
#  Makefile for PLUME / JET-PLUME:
#  (Judging Energy Transfer in a) Plasma in a Linear Uniform Magnetic Environment
#
#################################################################################


# SYSTEM=IFORT
SYSTEM=GFORT

PACK = Makefile \
	src/*.f90 \
	*.in \
	README* \
	data/gfort/ \
	include/

#FLAGS=
#NOTE: For Intel Fortran compiler ifort
ifeq ($(SYSTEM),IFORT)
	FLAGS=  -O3 -r8 -double-size 128
	COMP=   ifort
endif

ifeq ($(SYSTEM),GFORT)
#NOTE	: For gfortran
#gfortran orders roots a bit differently...
	FLAGS=  -O3 -DDOUBLE -fdefault-real-8
	COMP= gfortran
endif
LIBS=

LFMOD=	nrtype.o nrutil_trim.o vars.o functions.o bessel.o disprels.o fpc.o

LFX=  	plume.o

VPATH= src:include:/usr/include/

###############################################################################
all: heat

heat: $(LFMOD) $(LFX)
	$(COMP) -o plume.e $(FLAGS) $(LIBS) $(LFMOD) $(LFX)
	mv *.o include/
	mv *.mod include/

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
%.o : %.f90
	$(COMP) -c $(FLAGS) $<

#########Dependencies
plume.o:	vars.o functions.o fpc.o

functions.o:    vars.o

disprels.o:	vars.o bessel.o

bessel.o:	nrtype.o nrutil_trim.o
