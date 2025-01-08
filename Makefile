#################################################################### OVERVIEW
#  Makefile for PLUME:
#  Plasma in a Linear Uniform Magnetic Environment
#
#   VERSION 1.0
#
#  Version notes:
#
#  LAST UPDATE:  2024/01/08
###############################################################################

#SYSTEM=IFORT
SYSTEM=GFORT

PACK = Makefile \
	src/*.f90 \
	*.md
	inputs/example/*.in
	include/

#FLAGS=
#NOTE: For Intel Fortran compiler ifort
ifeq ($(SYSTEM),IFORT)
	FLAGS=  -O3 -r8 -double-size 128 #this works!
	COMP=   ifort
endif

ifeq ($(SYSTEM),GFORT)
	FLAGS=  -O3 -DDOUBLE -fdefault-real-8 -funroll-loops -ffast-math #
	COMP= gfortran
endif
LIBS=	

LFMOD=	nrtype.o nrutil_trim.o vars.o functions.o bessel.o disprels.o

LFX=  	plume.o 

VPATH= src:include:/usr/include/

###############################################################################
all: heat

heat: $(LFMOD) $(LFX)
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
%.o : %.f90
	$(COMP) -c $(FLAGS) $<

#########Dependencies
plume.o:	vars.o functions.o

functions.o:    vars.o     

disprels.o:	vars.o bessel.o

bessel.o:	nrtype.o nrutil_trim.o
