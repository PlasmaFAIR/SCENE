SUBS =  ddot.o dcopy.o  drotg.o  drot.o dswap.o dnrm2.o daxpy.o\
	dasum.o dqrsl.o dchdc.o dsvdc.o dqrdc.o dtrco.o dtrsl.o \
	dcfcr.o dcfcr1.o dcpmut.o dcrtz.o dctsx.o ddcom.o ddiag.o\
	dftkf.o dgcv.o dgcv1.o dgemv.o dmakek.o dmaket.o dpdcr.o\
	dpmse.o dpred.o dprmut.o dptpss.o dreps.o drsap.o dset.o\
	dsetup.o dsgdc.o dsgdc1.o dsnsm.o dsuy.o dtpss.o dtsvdc.o\
	duni.o dvl.o dvlop.o dvmin.o dzdc.o fact.o mkpoly.o prmut.o

SRCS =  ddot.f dcopy.f drotg.f  drot.f dswap.f dnrm2.f daxpy.f\
	dasum.f dqrsl.f dchdc.f dsvdc.f dqrdc.f dtrco.f dtrsl.f\
	dcfcr.f dcfcr1.f dcpmut.f dcrtz.f dctsx.f ddcom.f ddiag.f\
	dftkf.f dgcv.f dgcv1.f dgemv.f dmakek.f dmaket.f dpdcr.f\
	dpmse.f dpred.f dprmut.f dptpss.f dreps.f drsap.f dset.f\
	dsetup.f dsgdc.f dsgdc1.f dsnsm.f dsuy.f dtpss.f dtsvdc.f\
	duni.f dvl.f dvlop.f dvmin.f dzdc.f fact.f mkpoly.f prmut.f\


#LIBRARIES =  \
# -L/hwdisks/home/bbp501/libscene/ \
# -lBLAS 


EXE = libgcvfull.a

COMPILE = ifort -c 

LOAD = ifort -o

LIBLOAD = ar rc


#testtpss : testtpss.o ${SUBS}
#	f77 -o testtpss testtpss.o ${LIBRARIES} ${SUBS} 

#testptpss : testptpss.o ${SUBS}
#	f77 -o testptpss testptpss.o  ${LIBRARIES} ${SUBS}

#inteqn : inteqn.o mktpar.o mkxys.o ${SUBS}
#	f77 -o inteqn inteqn.o mktpar.o mkxys.o ${LIBRARIES} ${SUBS}

all: $(EXE)

$(EXE):		$(SUBS)
	$(LIBLOAD) $(EXE) $(SUBS) $(LIBRARIES)
	cp $(EXE) ..



.f.o:
	$(COMPILE) $< 

clean:
	rm -f *.o *.out
