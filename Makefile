#################################################################### OVERVIEW
#  Makefile for PLUME:
#  Plasma in a Linear Uniform Magnetic Environment
#
#   VERSION 0.0
#
#  Version notes:
#
#  LAST UPDATE:  2023/01/23
###############################################################################

 SYSTEM=IFORT
#SYSTEM=GFORT

PACK = Makefile \
	src/*.f90 \
	*.in \
	README* \
	data/gfort/ \
	include/

#FLAGS=
#NOTE: For Intel Fortran compiler ifort
ifeq ($(SYSTEM),IFORT)
	FLAGS=  -O3 -r8 -double-size 128 #this works!
	COMP=   ifort
endif

ifeq ($(SYSTEM),GFORT)
#NOTE	: For gfortran
#gfortran orders roots a bit differently...
	FLAGS=  -O3 -DDOUBLE -fdefault-real-8 #this works!
	COMP= gfortran
endif
LIBS=

#LFMOD=	nrtype.o nrutil_trim.o bessels.o funcs.o complex_root.o \
	disprels.o

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
#	cp pack_plume_`date +'%y%m%d'`.tar.gz ~/Dropbox/Codes/PLUME

#########Rules
%.o : %.f90
	$(COMP) -c $(FLAGS) $<

#########Dependencies
plume.o:	vars.o functions.o fpc.o

functions.o:    vars.o

disprels.o:	vars.o bessel.o

bessel.o:	nrtype.o nrutil_trim.o
