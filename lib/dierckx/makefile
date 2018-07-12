F90 = gfortran

SRCS= \
fporde.f fprati.f fpback.f fprota.f fpgivs.f \
fpdisc.f fpbspl.f fprank.f fpsurf.f surfit.f \
fpbisp.f bispev.f

OBJS= \
fporde.o fprati.o fpback.o fprota.o fpgivs.o \
fpdisc.o fpbspl.o fprank.o fpsurf.o surfit.o \
fpbisp.o bispev.o

EXE= libsurfit.a

DEBUG = #-heap-arrays -g -traceback -check bounds -check arg_temp_created -warn interface
#OPTIMIZE = # -O3
COMPILE = $(F90) -c $(OPTIMIZE) $(DEBUG)
LOAD = ar rc  


.SUFFIXES:

.SUFFIXES: .o .f .f90

all: $(EXE)

$(EXE): $(OBJS)
	$(LOAD) $@ $(OBJS)
	cp $(EXE) ../

.f.o:
	$(COMPILE) $<

libregrid.a: $(OBJS)

libs: libsurfit.a

#---------------------------------------------------------
clean:
	rm -f *.a *.o
