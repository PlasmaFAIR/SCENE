all:
	make -C trunk

clean:
	make -C trunk clean
	make -C lib/dierckx clean
	rm -f lib/libsurfit.a
