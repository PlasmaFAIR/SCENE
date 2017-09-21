SUBS = nbdriver.o beamconsts.o beams_mfp.o calcbeams.o coulomb.o \
	eisplit.o erf.o fastions.o frate.o gausswts.o \
	getrho.o hofr.o hunt.o qsimp.o sigfit.o sinteg.o \
	rinteg.o zinteg.o ssum.o sig_olson.o

SRCS = nbdriver.f beamconsts.f beams_mfp.f calcbeams.f coulomb.f \
	eisplit.f erf.f fastions.f frate.f gausswts.f \
	getrho.f hofr.f hunt.f qsimp.f sigfit.f sinteg.f \
	rinteg.f zinteg.f ssum.f sig_olson.f


LIBRARIES =  \
 -L/hwdisks/home/bbp501/libscene/ \
 -lBLAS \
 -llinpack

EXE = libnbeams.a

COMPILE = ifort -c 

LOAD = ar rc 


#testtpss : testtpss.o ${SUBS}
#	f77 -o testtpss testtpss.o ${LIBRARIES} ${SUBS} 

#testptpss : testptpss.o ${SUBS}
#	f77 -o testptpss testptpss.o  ${LIBRARIES} ${SUBS}

#inteqn : inteqn.o mktpar.o mkxys.o ${SUBS}
#	f77 -o inteqn inteqn.o mktpar.o mkxys.o ${LIBRARIES} ${SUBS}

all: $(EXE)

$(EXE):		$(SUBS)
	$(LOAD) $(EXE) $(SUBS)
	cp $(EXE) ..

.f.o:
	$(COMPILE) $(LIBRARIES) $<

clean:
	rm -f *.o *.out
