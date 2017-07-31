SUBS =  dqrsl.o dchdc.o dsvdc.o dqrdc.o dtrco.o dtrsl.o



SRCS =  dqrsl.f dchdc.f dsvdc.f dqrdc.f dtrco.f dtrsl.f


EXE = liblinpack.a

LIBRARIES = \
 -L/hwdisks/home/bbp501/libscene/ \
 -lBLAS


COMPILE = ifort -c
LOAD =  ar rc 

all:	$(EXE)


$(EXE): 		$(SUBS)
	$(LOAD) $(EXE)  $(SUBS) 
	cp $(EXE) ../

.f.o:
	$(COMPILE) $(LIBRARIES) $<

clean:
	rm -f *.o
