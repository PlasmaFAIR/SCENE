F90 = gfortran

NAME=	scene
CMD=    scene
#
SRCS= \
param.f90 balpar.f90  balloon.f90 usrcal.f90 \
coldat.f90 scene.f90 init.f90 equil.f90 profs.f90 flxsur.f90  \
flxav.f90 hirsig.f90 hirsh.f90 splines.f signeo.f90 \
torcur.f90 extj.f90 output.f90 graphs.f90 ffdgen.f90 elite_data.f90 \
baleq.f90 test_profs.f90 helena.f90 bandat.f90 setnt.f90 \
chease_dat.f90 gfile.f90 flxorb.f90 lynton.f90 fastbs.f90 caxe_out.f90 \
bandit.f90 mercier.f90 gs2.f90 astradat.f90 tokameq.f90 \
getdata.f90  beamconsts.f beams_mfp.f \
calcbeams.f coulomb.f \
eisplit.f erf.f fastions.f frate.f gausswts.f \
getrho.f hofr.f hunt.f qsimp.f sigfit.f sinteg.f \
rinteg.f zinteg.f ssum.f sig_olson.f \
nbicur.f90 geqdsk.f90 outcdf.f90 toms790.f90

OBJS= \
param.o balpar.o balloon.o \
coldat.o scene.o init.o equil.o profs.o flxsur.o \
flxav.o hirsig.o hirsh.o splines.o signeo.o \
torcur.o extj.o output.o ffdgen.o usrcal.o elite_data.o \
#baleq.o test_profs.o helena.o bandat.o setnt.o chease_dat.o \
gfile.o flxorb.o lynton.o fastbs.o caxe_out.o bandit.o  \


OBJS_1= $(SRCS:.f90=.o)
OBJS=$(OBJS_1:.f=.o)
INPUT=
INCLUDE=  -I$(YORK_NETCDF_DIR)/include 
LIBRARIES =  \
 -llapack \
 -lblas \
 -L$(YORK_GHOST_DIR)/lib -lghost -lgrid_xghost \
 -L$(YORK_NETCDF_DIR)/lib -lnetcdff -lnetcdf \
 -L../lib -lsurfit \
# -lgcvfull

#        -lghost -lgrid_xghost -L/usr/lib/libg2c.so.0 \
# /hwdisks/sfw/nag/ifort/fllux21dcl/lib/libnag_nag.a -lpthread  \


#  -lghost -lgrid -lU77 -lnag_absoft -lf95 -lm -lc \






 # /fuslwf/home/kthcmr/GS2D/gs2dlib.a

DEBUG     = -std=f2008 -Wall -Wextra -fbounds-check #-heap-arrays  -g -traceback  -check bounds -check arg_temp_created -warn interface
OPTIMIZE  = #-O # -O3 
LISTING   = # -listing
SPECIAL   = # -p/fuslwf/home/kthcmr/GS2D
DOUBLE    = #-r8 -d8 #-qautodbl=dblpad
LIN       = # -YEXT_NAMES=LCS -B108
COMPILE   = $(F90) -c $(LIN) $(LISTING) $(SPECIAL) $(DOUBLE) $(OPTIMIZE) $(DEBUG)  $(INCLUDE)
LOAD      = $(F90) $(SPECIAL) -o $(CMD) $(LIBS)
.SUFFIXES:

.SUFFIXES: .o .f .f90

all: $(CMD)
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"

$(CMD):         $(OBJS) 
	$(LOAD) $(OBJS) $(LIBRARIES)
	chmod 755 $(CMD)

.f.o:	
	$(COMPILE) $<

.f90.o:	
	$(COMPILE) $<
#-------------------------------------------------------------------------
clean:
	rm -f $(CMD) *.a *.o *.mod





