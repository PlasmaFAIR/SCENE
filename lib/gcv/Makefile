SUBS = dcfcr.o dcfcr1.o dcpmut.o dcrtz.o dctsx.o ddcom.o ddiag.o\
	dftkf.o dgcv.o dgcv1.o dgemv.o dmakek.o dmaket.o dpdcr.o\
	dpmse.o dpred.o dprmut.o dptpss.o dreps.o drsap.o dset.o\
	dsetup.o dsgdc.o dsgdc1.o dsnsm.o dsuy.o dtpss.o dtsvdc.o\
	duni.o dvl.o dvlop.o dvmin.o dzdc.o fact.o mkpoly.o prmut.o

SRCS = dcfcr.f dcfcr1.f dcpmut.f dcrtz.f dctsx.f ddcom.f ddiag.f\
	dftkf.f dgcv.f dgcv1.f dgemv.f dmakek.f dmaket.f dpdcr.f\
	dpmse.f dpred.f dprmut.f dptpss.f dreps.f drsap.f dset.f\
	dsetup.f dsgdc.f dsgdc1.f dsnsm.f dsuy.f dtpss.f dtsvdc.f\
	duni.f dvl.f dvlop.f dvmin.f dzdc.f fact.f mkpoly.f prmut.f

LIBRARIES =  \
 -L/hwdisks/home/bbp501/libscene/ \
 -lBLAS \
 -llinpack

EXE = libgcv.a

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
