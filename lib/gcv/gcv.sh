#!/bin/sh 
#    This is a shell archive.
#    Run the following text with /bin/sh to extract.

YASASTART=`pwd`

echo yasa: extracting Makefile
cat - << \Funky!Stuff! > Makefile
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

testtpss : testtpss.o ${SUBS}
	f77 -o testtpss testtpss.o ${SUBS} -llinpack

testptpss : testptpss.o ${SUBS}
	f77 -o testptpss testptpss.o ${SUBS} -llinpack

inteqn : inteqn.o mktpar.o mkxys.o ${SUBS}
	f77 -o inteqn inteqn.o mktpar.o mkxys.o ${SUBS} -llinpack
Funky!Stuff!
echo yasa: extracting dcfcr.f
cat - << \Funky!Stuff! > dcfcr.f
      subroutine dcfcr (fg,ldfg,nnull,qr,ldqr,npar,pmh,qraux,
     * sgpvt,capz,ldcapz,svals,npsing,v,ldv,nlamht,w1,z,coef,penlty,
     * work,info)
      integer ldfg,nnull,ldqr,npar,pmh,sgpvt(npar),ldcapz,npsing,ldv,
     * info
      double precision fg(ldfg,nnull),qr(ldqr,pmh),
     * qraux(pmh),capz(ldcapz,*),svals(npsing),v(ldv,npsing),nlamht,
     * w1(nnull),z(npsing),coef(npar),penlty,work(nnull)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nnull		number of columns in g
c   qr(ldqr,pmh)	information on the Householder transformations 
c			that define q and r
c   ldqr		leading dimension of qr as declared
c			in the calling program
c   npar		number of rows in q
c   pmh			number of columns in r
c   qraux(pmh)		auxiliary information on the
c			qr Householder transformations
c   sgpvt(npar)		permuted indices from the pivoted Cholesky 
c			decomposition of the matrix which defines the 
c			semi-norm 
c   capz(ldcapz,pmh)	first part of the rotated design matrix
c   ldcapz		leading dimension of capz as declared
c			in the calling program
c   svals(npsing)	singular values 
c   npsing		number of positive singular values
c   v(ldv,npsing)	right singular vectors corresponding to svals
c   ldv			leading dimension of v as declared
c			in the calling program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [(D**2 +nlamht)**-1]Dz
c   coef(npar)		estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c			  2 : error in dtrsl, r is singular
c Work arrays
c   work(nnull)		double precision work vector
c
c   Subprograms Called Directly:
c	Linpack - dqrsl dtrsl dtrco
c	Blas    - ddot dgemv
c	Other   - dprmut 
c
c $Header: dcfcr.f,v 2.100.1.1 86/10/07 12:47:09 lindstrom Exp $
c
      integer i
      double precision dummy(1),machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			form gamma-hat and penlty
      do 20 i = 1,npsing
         z(i) = z(i)*svals(i)/(svals(i)**2 + nlamht)
   20 continue
      call dgemv('N',pmh,npsing,1.0d0,v,ldv,z,1,0.0d0,coef,1)
      penlty = ddot(pmh,coef,1,coef,1)
c			form delta
      call dcopy(nnull,w1,1,coef(pmh+1),1)
      call dgemv('N',nnull,pmh,-1.0d0,capz,ldcapz,coef,1,1.0d0,
     *   coef(pmh+1),1)
c			check condition number of g
      if (nnull .ne. 0) then
          call dtrco(fg,ldfg,nnull,rcond,work,1)
          if (rcond .le. machpr*100) then
             info = 1
             return
          endif
          call dtrsl (fg,ldfg,nnull,coef(pmh+1),01,info)
      endif
c			form theta from gamma and delta
      call dtrsl (qr,ldqr,pmh,coef,11,info)
      if (info .ne. 0) then
         info = 2
         return
      endif
      call dqrsl (qr,ldqr,npar,pmh,qraux,coef,coef,dummy,dummy,dummy,
     * dummy,10000,info)
      call dprmut (coef,npar,sgpvt,1)
      return
      end
Funky!Stuff!
echo yasa: extracting dcfcr1.f
cat - << \Funky!Stuff! > dcfcr1.f
      subroutine dcfcr1(fg,ldfg,ncts1,fgaux,u,ldu,f1kf2,ldfkf,nuobs,
     * svals,npsing,nlamht,w1,z,coef,penlty,work,info)
      integer ldfg,ldfkf,ncts1,ldu,nuobs,npsing,info
      double precision fg(ldfg,ncts1),fgaux(ncts1),u(ldu,*),
     * f1kf2(ldfkf,*),svals(npsing),nlamht,w1(ncts1),z(npsing),coef(*),
     * penlty,work(*)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,ncts1)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   ncts1		number of columns in g
c   fgaux(ncts1)	auxiliary information on the fg Householder
c			transformations
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   f1kf2(ldfkf,nuobs-ncts1) f1 k f2
c   ldfkf		leading dimension of f1kf2 as declared
c			in the calling program
c   nuobs		number of rows in fg
c   svals(npsing)	singular values of f2'k f2
c   npsing		number of positive singular
c   nlamht		nobs*(lambda hat)
c   w1(ncts1)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [ (D**2 +nlamht)**-1 ] D z
c   coef(nuobs+ncts1)	estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c
c Work Arrays:
c   work(nuobs-ncts1)  	double precision work vector
c
c   Subprograms Used:
c      Linpack - dqrsl dtrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: dcfcr1.f,v 2.100.1.1 86/10/07 12:47:31 lindstrom Exp $
c
      integer i,j,nmnct,nctp1,locinf
      double precision dummy,machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      nmnct = nuobs - ncts1
      nctp1 = ncts1 + 1
c			form g and penalty
      do 20 j = 1,npsing
         z(j) = z(j)*svals(j)/(svals(j)**2 + nlamht)
   20 continue
      penlty = ddot(npsing,z,1,z,1)
c			z now contains g
c			compute xi 
      do 30 i = 1,npsing
         coef(i) = z(i)/svals(i)
   30 continue
      call dgemv('N',nmnct,npsing,1.0d0,u,ldu,coef,1,0.0d0,
     *  work,1)
      do 40 j = 1,ncts1
         coef(ncts1+j) = 0.0d0
   40 continue
      call dcopy(nmnct,work,1,coef(2*ncts1+1),1)
      call dqrsl(fg,ldfg,nuobs,ncts1,fgaux,coef(nctp1),coef(nctp1),
     * dummy,dummy,dummy,dummy,10000,locinf)
c			compute beta
      call dcopy(ncts1,w1,1,coef,1)
      call dgemv('N',ncts1,nmnct,-1.0d0,f1kf2,ldfkf,work,1,1.0d0,
     *  coef,1)
c			check condition number of g
      call dtrco(fg,ldfg,ncts1,rcond,work,1)
      if (rcond .le. machpr*100) then
         info = 1
         return
      endif
      call dtrsl (fg,ldfg,ncts1,coef,01,info)
      return
      end
Funky!Stuff!
echo yasa: extracting dcpmut.f
cat - << \Funky!Stuff! > dcpmut.f
      subroutine dcpmut (x,ldx,nobs,npar,jpvt,job)
      integer ldx,nobs,npar,jpvt(npar),job
      double precision x(ldx,npar)
c
c Purpose: permute the columns of the matrix x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(ldx,npar)		matrix whose columns are to be permuted
c   ldx			leading dimension of x as declared
c			in the calling program
c   nobs		number of rows of x used
c   npar		number of columns of x
c   jpvt(npar)		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			else backward permutation
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(ldx,npar)		matrix with columns permuted
c
c Subprograms Called Directly
c     Blas	- dswap
c
c  Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dcpmut.f,v 2.100.1.1 86/10/07 12:47:45 lindstrom Exp $
c
      integer i,j,k
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               call dswap (nobs,x(1,j),1,x(1,k),1)
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0) then
c		backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               call dswap (nobs,x(1,i),1,x(1,j),1)
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
Funky!Stuff!
echo yasa: extracting dcrtz.f
cat - << \Funky!Stuff! > dcrtz.f
      subroutine dcrtz (x,ldx,nobs,npar,pmh,qr,ldqr,qraux,sgpvt,work,
     * info)
      integer ldx,nobs,npar,pmh,ldqr,sgpvt(npar),info
      double precision x(ldx,npar),qr(ldqr,pmh),qraux(pmh),work(npar)
c
c Purpose: create z from x using permutations and householder 
c	   transformation from sigma
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   qr(ldqr,pmh)	information on the qr decomposition of the 
c			Cholesky factor of the semi-norm
c   ldqr    		leading dimension of qr as declared in the 
c			calling	program
c   qraux(pmh)		auxiliary information on the qr decomposition 
c			of the factor of the semi-norm
c   sgpvt(npar) 	permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c
c On Exit:
c   x(ldx,npar)		z (x after householder transformation)
c   info    		error indicator
c			  0 : successful completion
c			  1 : error in dtrsl, r is singular
c
c Work Arrays:
c   work(npar)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl dtrsl 
c	Blas    - dcopy 
c	Other   - dcpmut 
c
c Subprograms Called Indirectly:
c	Blas    - dswap
c
c $Header: dcrtz.f,v 2.100.1.1 86/10/07 12:47:59 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c
      info = 0
c			permute columns of x according to sgpvt
      call dcpmut (x,ldx,nobs,npar,sgpvt,0)
c			apply Householder transformations to rows of
c			permuted x to generate z1 and z2
      do 20 i = 1,nobs 
         call dcopy (npar,x(i,1),ldx,work,1)
         call dqrsl (qr,ldqr,npar,pmh,qraux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dtrsl (qr,ldqr,pmh,work,01,info)
         if (info .ne. 0) then
            info = 1
            return
         endif
         call dcopy (npar,work,1,x(i,1),ldx)
   20 continue
      return
      end
Funky!Stuff!
echo yasa: extracting dctsx.f
cat - << \Funky!Stuff! > dctsx.f
      subroutine dctsx(ts1,ldts1,ncts1,s2,lds2,ncov2,nobs,tbsb1,ldtbsb,
     * nb,sigma,ldsigm,x,ldx,npar,fgaux,work)
      integer ldts1,ncts1,lds2,ncov2,nobs,ldtbsb,nb,ldsigm,ldx,npar
      double precision ts1(ldts1,ncts1),s2(lds2,*),tbsb1(ldtbsb,ncts1),
     * sigma(ldsigm,*),x(ldx,*),fgaux(ncts1),work(nb)
c
c Purpose: compute x and sigma from t,s and k.
c
c On Entry:
c   ts1(nobs,ncts1)	[t:s1]
c   ldts1		leading dimension of ts1 as declared in the
c			calling	program 
c   ncts1	 	number of columns in [t:s1]
c   s2(lds2,ncov2) 	columns contain the covariates which do not 
c			duplicate the replication pattern of des
c   lds2		leading dimension of s2 as declared in the
c			calling	program 
c   ncov2		number of columns in s2
c   nobs		number of observations
c   tbsb1(ldtbsb,ncts1)	unique rows of t and s1 (or [t:s1] created from
c			basis functions)
c   ldtbsb		leading dimension of tbsb1 as declared in 
c			calling program
c   nb			number of basis functions
c   sigma		ku (ku must start in the ncov2+1 row and
c			the ncov2+1 column)
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   x(ldx,npar)		k (k must start in the 1st row and column)
c   ldx			leading dimension of x as
c			declared in the calling program
c
c On Exit:
c   tbsb1(ldtbsb,ncts1)	qr decomposition of tbsb1
c   sigma(ldsigm,npar)  symmetric matrix that defines the semi-norm
c				[0:   0  ]
c				[0:f2'ku f2]
c   x(ldx,npar)		design matrix [t:s1:s2:kf2]
c   npar		number of parameters (nb+ncov2)
c   fgaux(ncts1)	auxiliary vector for qr decomposition of tbsb1
c
c Work Arrays:
c   work(nb)		double precision work vector
c
c Subprograms Called Directly
c	Linpack - dqrdc dqrsl
c	Blas    - dcopy
c	Other   - dset dftkf
c
c Subprograms Called Indirectly
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dctsx.f,v 2.100.1.1 86/10/07 12:48:06 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c
      npar = nb + ncov2
      call dqrdc(tbsb1,ldtbsb,nb,ncts1,fgaux,0,0.0d0,0)
c			calculate k f2  put in last nb -ncts1 columns 
c			of x
      do 10 i=1,nobs 
         call dcopy(nb,x(i,1),ldx,work,1)
         call dqrsl(tbsb1,ldtbsb,nb,ncts1,fgaux,work,dummy,work,dummy,
     *    dummy,dummy,01000,locinf)
         call dcopy(nb-ncts1,work(ncts1+1),1,x(i,ncov2+ncts1+1),ldx)
   10 continue
c			copy [t:s1] into first ncts1 columns of x
      do 20 i=1,ncts1
         call dcopy(nobs,ts1(1,i),1,x(1,i),1)
   20 continue
c			copy s2 into next ncov2 columns of x
      do 30 i=1,ncov2
         call dcopy(nobs,s2(1,i),1,x(1,ncts1+i),1)
   30 continue
c			calculate sigma
      call dftkf(tbsb1,ldtbsb,nb,ncts1,fgaux,sigma(ncov2+1,ncov2+1),
     * ldsigm,work)
      do 40 i = 1,ncts1+ncov2
         call dset(npar,0.0d0,sigma(1,i),1)
   40 continue
      do 50 i = ncts1+ncov2+1,npar
         call dset(ncts1+ncov2,0.0d0,sigma(1,i),1)
   50 continue
      return
      end
Funky!Stuff!
echo yasa: extracting ddcom.f
cat - << \Funky!Stuff! > ddcom.f
      subroutine ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,
     * npsing,svals,sgpvt,dcaux,ldcaux,normk,work,lwa,iwork,liwa,
     * job,info)
      integer ldx,ldsigm,nobs,npar,nnull,npsing,sgpvt(npar),
     * ldcaux,lwa,liwa,iwork(liwa),job,info
      double precision x(ldx,npar),sigma(ldsigm,npar),tau,svals(*),
     * dcaux(ldcaux),normk,work(lwa)
c
c Purpose: decompose sigma and the design matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program. Must be at least max(npar,nobs)
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested
c   ldcaux		length of auxiliary vector as declared in the 
c			calling program. 
c			Must be at least (npar-nnull)**2+2*npar-nnull
c   job 		job nonzero use truncation 
c			job = 0 no truncation
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c			must be	passed to dgcv intact
c   sigma(ldsigm,npar)	overwritten with the qr decomposition of the 
c			Cholesky factor of sigma
c			must be passed intact to dgcv
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   npsing       	number of singular values calculated if info = 0
c			otherwise npsing contains nonzero info from 
c			dsvdc.
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info from dsvdc is nonzero then svals 
c			is as it was returned from dsvdc.
c   sgpvt(npar)		permuted indices from the Cholesky decomposition
c			with pivots of sigma
c			must be passed intact to dgcv 
c   dcaux(ldcaux)	auxiliary vector, must be passed intact to dgcv 
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : dimension error or tau < 0
c			   2 : ldcaux (length of dcaux) is too small
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   10 < info < 20 : 10 + nonzero info returned
c					    from dsgdc
c			   20 < info < 30 : 20 + nonzero info returned 
c					    from dcrtz
c			   30 < info < 40 : 30 + nonzero info returned 
c					    from dzdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program Must be at least 
c			(npar-nnull)*(nobs-nnull+1)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling
c			program 
c			must be at least npar-nnull
c
c Subprograms Called Directly:
c	Gcvpack - dsgdc dcrtz dzdc
c	Blas	- dcopy
c
c Subprograms Called Indirectly:
c	Gcvpack - dtsvdc
c	Linpack - dtrco dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - dcopy dswap ddot
c	Other   - dprmut dcpmut dset
c
c $Header: ddcom.f,v 2.100.1.1 86/10/07 12:48:26 lindstrom Exp $
c
      integer pmh,pp1,pmhp1,wsize,nmh,sinfo
c
      sinfo = 0
      info = 0
c			check dimensions
      nmh = nobs - nnull
      pmh = npar - nnull
      pmhp1 = pmh + 1
      pp1 = npar + 1
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0) .or. (ldx .lt. nobs).or. (ldx .lt. npar)
     * .or. (tau .lt. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 + 2*npar -nnull
      if (ldcaux .lt. wsize) then
         info = 2
         return
      endif
      wsize = pmh*nmh + pmh + nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      if (liwa .lt. pmh) then
         info = 4
         return
      endif
c			decompose sigma
      call dsgdc (sigma,ldsigm,npar,nnull,dcaux,sgpvt,info)
      if (info.gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) then
	  sinfo = info
	  pmh = npar - nnull
          pmhp1 = pmh + 1
      endif
c			create z
      call dcrtz(x,ldx,nobs,npar,pmh,sigma,ldsigm,dcaux,sgpvt,work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			create j and decompose j2
      call dzdc(x,ldx,nobs,npar,pmh,tau,dcaux(pmhp1),dcaux(pp1),
     * npsing,dcaux(npar+pmhp1),pmh,normk,work,lwa,iwork,liwa,job,info)
c			copy svals
      call dcopy(pmh,dcaux(pp1),1,svals,1)
      if (info .gt. 0) then
         info = info + 30
         return
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
Funky!Stuff!
echo yasa: extracting ddiag.f
cat - << \Funky!Stuff! > ddiag.f
      subroutine ddiag(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu,
     * nlamht,adiag,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,adiag(nobs),work(*)
c
c Purpose: determine the diagonal of hat matrix for nobs*lamhat
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder
c			transformations
c   svals(npsing)	singular values 
c   npsing		number of positive singular values 
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   nlamht		nobs*lambda hat
c
c On Exit:
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - ddot dgemv
c	Other   - dset
c
c $Header: ddiag.f,v 2.100.1.1 86/10/07 13:28:57 lindstrom Exp $
c
      integer i,j,hp1,locinf,nmh,np1
      double precision dummy(1)
      double precision ddot
c
      np1 = nobs + 1
      hp1 = nnull + 1
      nmh = nobs - nnull
c			form adiag
      do 20 i = 1,nobs 
         call dset(nobs,0.0d0,work,1)
         work(i)=1.0d0
         call dqrsl(fg,ldfg,nobs,nnull,fgaux,work,dummy,work,dummy, 
     *    dummy,dummy,01000,locinf)
         adiag(i)=ddot(nnull,work,1,work,1)
	 call dgemv('T',nmh,npsing,1.0d0,u,ldu,work(hp1),1,0.0d0,
     *    work(np1),1)
         do 10 j=1,npsing
            work(nobs+j)=work(nobs+j)*svals(j)/dsqrt(svals(j)**2+nlamht)
   10    continue
         adiag(i)=adiag(i) + ddot(npsing,work(np1),1,work(np1),1)
   20 continue
      return
      end
Funky!Stuff!
echo yasa: extracting dftkf.f
cat - << \Funky!Stuff! > dftkf.f
      subroutine dftkf(fg,ldfg,nrf,ncg,fgaux,kk,ldkk,work)
      integer ldfg,nrf,ncg,ldkk
      double precision fg(ldfg,ncg),fgaux(ncg),kk(ldkk,nrf),work(nrf)
c
c Purpose: create f'k f.
c
c On Entry:
c   fg(ldfg,ncg)	qr decomposition of [t:s1]
c   ldfg		leading dimension of fg as declared in the
c   			calling program 
c   nrf 		number of rows in f
c   ncg			number of columns in g
c   fgaux(ncg)		auxiliary information on the qr decomposition
c			of [t:s1]
c   kk(ldkk,nrf) 	k
c   ldkk		leading dimension of kk as declared in the
c   			calling program 
c
c On Exit:
c   kk(ldkk,nrf)  	f'k f
c
c Work Array:
c   work(nrf)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dftkf.f,v 2.100.1.1 86/10/07 12:48:38 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c	  		calculate k f, store in kk
      do 10 i=1,nrf 
         call dcopy(nrf,kk(i,1),ldkk,work,1)
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dcopy(nrf,work,1,kk(i,1),ldkk)
   10 continue
c	  		calculate f'k f
      do 20 i=1,nrf 
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,kk(1,i),dummy,kk(1,i),dummy,
     *    dummy,dummy,01000,locinf)
   20 continue
      return
      end
Funky!Stuff!
echo yasa: extracting dgcv.f
cat - << \Funky!Stuff! > dgcv.f
      subroutine dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,sgpvt,dcaux,
     * ldcaux,adiag,lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,lwa,
     * job,info)
      integer ldx,ldsigm,nobs,npar,nnull,sgpvt(npar),ldcaux,ntbl,ldtbl,
     * lwa,info,job
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * dcaux(ldcaux),adiag(nobs),lamlim(2),dout(4),coef(npar),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized
c 	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix as returned by ddcom
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	Cholesky factor of sigma as returned by ddcom
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   sgpvt(npar)		permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c   dcaux(ldcaux)	auxiliary vector from ddcom
c   ldcaux		length of ldcaux. Must be at least
c			(npar-nnull)**2+2*npar-nnull
c   adiag(nobs)	 	"true" y values on entry if predictive mse if 
c			requested 
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(npar)		coefficient estimates
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2) 
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : ldcaux (length of dcaux) is too small
c			  4 : lwa (length of work) is too small
c			  5 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned
c					  from dvlop
c			 20 < info < 30 : 20 +nonzero info returned
c					  from dcfcr
c
c Work Arrays:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program must be at least (npar-nnull)+nobs
c
c Subprograms Called Directly:
c	Gcvpack - drsap dvlop dpmse dcfcr dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvpack - dvmin dvl
c	Linpack - dqrsl dtrsl
c	Blas    - dcopy ddot dgemv
c	Other   - dprmut dset
c
c $Header: dgcv.f,v 2.100.1.1 86/10/07 12:49:16 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer mp1,mnh,pmh,minnp,pp1,pmhp1,wsize,npsing,i,jpmse,jlaml,
     * jadiag,sinfo
c
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      sinfo = 0
      info = 0
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      wsize = pmh**2 + 2*npar - nnull
      if (ldcaux .lt. wsize) then
         info = 3
         return
      endif
      wsize = pmh + nobs
      if (lwa .lt. wsize) then
         info = 4
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 5
         return
      endif
      pmhp1 = pmh + 1
      pp1 = npar + 1
c			calculate npsing
      minnp = min(nobs-nnull,pmh)
      do 30 i=1,minnp
         if (dcaux(npar+i)**2 .gt. 0) npsing = i
   30 continue
c			apply rotations to y
      mp1 = nnull + 1
      mnh = nobs - nnull
      call drsap (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),x(mp1,1),ldx,
     * mnh,npsing,y,ssqw2,addend,work)
c			minimize V(lambda)
      call dvlop (y(mp1),dcaux(pp1),nobs,nnull,npsing,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(x(1,pmhp1),ldx,nobs,nobs,nnull,dcaux(pmhp1),
     *    dcaux(pp1),npsing,x(mp1,1),ldx,y,y(mp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr (x(1,pmhp1),ldx,nnull,sigma,ldsigm,npar,
     * pmh,dcaux,sgpvt,x,ldx,dcaux(pp1),npsing,dcaux(npar+pmhp1),pmh, 
     * nlamht,y,y(mp1),coef,dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			calculate predicted values 
      call dpdcr (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),dcaux(pp1),
     * npsing,x(mp1,1),ldx,nlamht,y,y(mp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),
     * 	  dcaux(pp1),npsing,x(mp1,1),ldx,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
Funky!Stuff!
echo yasa: extracting dgcv1.f
cat - << \Funky!Stuff! > dgcv1.f
      subroutine dgcv1(fkf,ldfkf,y,nuobs,nobs,fg,ldfg,ncts1,fgaux,svals,
     * adiag,lamlim,ssqrep,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,lwa,
     * job,info)
      integer ldfkf,nuobs,nobs,ldfg,ncts1,ntbl,ldtbl,lwa,job,info
      double precision fkf(ldfkf,nuobs),y(nuobs),fg(ldfg,ncts1),
     * fgaux(ncts1),svals(*),adiag(nuobs),lamlim(2),ssqrep,dout(4),
     * coef(*),tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a semi-norm
c	thin plate spline model.
c
c On Entry:
c   fkf(ldfkf,nuobs) 	intermediate results as created by dsgdc1
c   ldfkf		leading dimension of fkf as declared in the 
c			calling program
c   y(nuobs)		B1'y 
c   nuobs		number of rows in fg
c   nobs		number of observations
c   fg(ldfg,ncts1)	qr decomposition of [t:s1]	
c   ldfg		leading dimension of fg as
c			declared in the calling program
c   ncts1		number of columns in [t:s1]	
c   fgaux(ncts1)	auxiliary information on the decomposition of 
c			[t:s1]
c   svals(nuobs-ncts1) 	positive singular values of the Cholesky factor
c			of f2'k f2
c   adiag(nuobs)	B1'(true y) if predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ssqrep		sum of squares for replication
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nuobs)		B1'(predicted values)
c   adiag(nuobs)	diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(nuobs+ncts1) 	estimated coefficients
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2)
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : lwa (length of work) is too small
c			  4 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned 
c					  from dvlop
c			 20 < info < 30 : 20 + nonzero info returned
c					  from dcfcr1
c
c Working Storage:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program
c			must be at least nuobs-ncts1+nobs
c
c Subprograms Called Directly:
c       Gcvlib - drsap dvlop dpmse dcfcr1 dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvlib  - dvl vmin
c	Linpack - dqrsl dtrsl
c	Blas    - ddot dcopy dgemv
c
c $Header: dgcv1.f,v 2.100.1.1 86/10/07 12:49:43 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer npsing,i,jpmse,jlaml,jadiag,nmnct,nctp1,sinfo
c
c
      sinfo = 0
      info = 0
      nmnct = nuobs - ncts1
      nctp1=ncts1+1
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      if ((nuobs.le.0).or.(ncts1 .le. 0).or.(nmnct .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      if (lwa .lt. nobs+nuobs - nmnct) then
         info = 3
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 4
         return
      endif
c			calculate npsing
      do 30 i=1,nmnct
         if (svals(i)**2 .gt. 0.0d0) npsing = i
   30 continue
c			apply rotations to y
      call drsap(fg,ldfg,nuobs,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,nmnct,
     * npsing,y,ssqw2,addend,work)
      addend = addend + ssqrep
      ssqw2 = ssqw2 + ssqrep
c			minimize V(lambda)
      call dvlop (y(nctp1),svals,nobs,ncts1,npsing,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0)  sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(fg,ldfg,nuobs,nobs,ncts1,fgaux,svals,npsing,
     *    fkf(nctp1,nctp1),ldfkf,y,y(nctp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr1(fg,ldfg,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,
     * fkf(1,nctp1),ldfkf,nuobs,svals,npsing,nlamht,y,y(nctp1),coef,
     * dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
      call dpdcr (fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     * fkf(nctp1,nctp1),ldfkf,nlamht,y,y(nctp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag(fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     *     fkf(nctp1,nctp1),ldfkf,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
Funky!Stuff!
echo yasa: extracting dgemv.f
cat - << \Funky!Stuff! > dgemv.f
Caveat receptor.  (Jack) dongarra@anl-mcs, (Eric Grosse) research!ehg
Compliments of netlib   Sun Jul  6 09:34:18 CDT 1986
C
C***********************************************************************
C
C     File of the  DOUBLE PRECISION  Level 2 BLAS routines:
C
C      DGEMV, DGBMV, DSYMV, DSBMV, DSPMV, DTRMV, DTBMV, DTPMV,
C      DGER , DSYR , DSPR ,
C      DSYR2, DSPR2,
C      DTRSV, DTBSV, DTPSV.
C
C     See:
C
C        Dongarra J. J., Du Croz J. J., Hammarling S. and Hanson R. J..
C        A proposal for an extended set of Fortran Basic Linear Algebra
C        Subprograms. Technical Memorandum No.41 (revision 1),
C        Mathematics and Computer Science Division, Argone National
C        Laboratory, 9700 South Cass Avenue, Argonne, Illinois 60439,
C        USA or NAG Technical Report TR4/85, Nuemrical Algorithms Group
C        Inc., 1101 31st Street, Suite 100, Downers Grove, Illinois
C        60515-1263, USA.
C
C***********************************************************************
C
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
      CHARACTER*1        TRANS
      INTEGER            M, N, LDA, INCX, INCY
      DOUBLE PRECISION   ALPHA, A( LDA, * ), X( * ), BETA, Y( * )
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y
*.
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           m.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Note that TRANS, M, N and LDA must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*     OK = ( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ).OR.
*    $       ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
*    $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )
*    $     .AND.
*    $     ( M.GE.0 )
*    $     .AND.
*    $     ( N.GE.0 )
*    $     .AND.
*    $     ( LDA.GE.M )
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
*
      INTEGER            I     , IX    , IY    , J     , JX    , JY
      INTEGER            KX    , KY    , LENX  , LENY
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
      DOUBLE PRECISION   TEMP
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.
     $    ( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         IF( BETA.NE.ONE )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         END IF
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( LENX - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( LENY - 1 )*INCY
         END IF
         IF( BETA.NE.ONE )THEN
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP
  100       CONTINUE
         ELSE
            JY = KY
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( IX )
                  IX    = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGEMV .
*
      END
Funky!Stuff!
echo yasa: extracting dmakek.f
cat - << \Funky!Stuff! > dmakek.f
      subroutine dmakek(m,n,dim,des,lddes,nb,desb,lddesb,kk,ldkk)
      integer m,n,dim,lddes,nb,lddesb,ldkk
      double precision des(lddes,dim),desb(lddesb,dim),kk(ldkk,nb)
c
c Purpose: create the k matrix.
c   
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			dimension of the space to be splined 
c   des(lddes,dim)	variables to be splined		
c   lddes		leading dimension of des as declared in the 
c			calling	program
c   nb			number of rows in desb 
c   desb(lddesb,dim)	positions of unique design points or basis 
c			functions
c   lddesb		leading dimension of desb as declared in the 
c			calling	program
c   ldkk		leading dimension of kk as declared in the
c			calling	program
c On Exit:
c   kk(ldkk,nb)		k matrix
c
c Subprograms Called:
c	Other   - fact
c
c $Header: dmakek.f,v 2.100.1.1 86/10/07 12:50:38 lindstrom Exp $
c
      integer i,j,k,fact
      double precision tauij,expo,theta,t,pi
c
c			t to be used in computation of theta
      pi = 4.0d0*atan(1.0d0)
      t = 2 ** (2*m) * pi**(dim/2.0d0) * fact(m-1)
c
c 			exponent for tauij
      expo = m - (dim / 2.0d0)
      if (dim .eq. 2*(dim/2)) then
c			1+dim odd
         theta = 1.0 / (0.5 * t * fact (m-dim/2))
         if ((2*m+dim) .eq. 4*((2*m+dim)/4)) theta = -theta
         do 30 i=1,n 
            do 20 j=1,nb 
               tauij = 0
               do 10 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   10          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo * 0.5 * log(tauij)
               endif
   20       continue
   30    continue
      else
c			1+dim even
c			compute theta
c			compute gamma(dim/2 - m)
         j = (1 - (dim-2*m)) / 2
         theta = sqrt(pi)
         do 40 i=1,j
	      theta = -theta / (i - 0.5d0)
   40    continue
	 theta = theta / t

         do 70 i=1,n 
            do 60 j=1,nb 
               tauij = 0
               do 50 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   50          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo
               endif
   60       continue
   70    continue
      endif
      end
Funky!Stuff!
echo yasa: extracting dmaket.f
cat - << \Funky!Stuff! > dmaket.f
      subroutine dmaket(m,n,dim,des,lddes,s1,lds1,ncov1,npoly,t,ldt,
     * wptr,info)
      integer m,n,dim,lddes,lds1,ncov1,npoly,ldt,wptr(dim),info
      double precision des(lddes,dim),s1(lds1,*),t(ldt,*)
c
c Purpose: create t matrix and append s1 to it.
c
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			number of columns in des
c   des(lddes,dim)	variables to be splined
c   lddes		leading dimension of des as declared in the 
c			calling program
c   s1(lds1,ncov1)	covariates which duplicate the replication 
c			structure of des
c   lds1		leading dimension of s1 as declared in the 
c			calling program
c   ncov1		number of columns in s1
c   ldt			leading dimension of t as declared in the 
c			calling program
c
c On Exit:
c   npoly		dimension of polynomial part of spline
c   t(ldt,npoly+ncov1)	[t:s1]
c   info 		error indication
c   			   0 : successful completion
c		 	   1 : error in creation of t
c Work Arrays:
c   wptr(dim)		integer work vector
c
c Subprograms Called Directly:
c	Blas  - dcopy
c	Other - mkpoly
c
c $Header: dmaket.f,v 2.100.1.3 86/11/21 11:32:19 lindstrom Exp $
c
      integer i,j,k,tt,nt,bptr,eptr
      integer mkpoly
c
      info = 0
      npoly = mkpoly(m,dim)
      call dset(n,1.0d0,t(1,1),1)
      nt = 1
      if (npoly .gt. 1) then
          do 10 j=1,dim 
             nt = j + 1
             wptr(j) = nt
             call dcopy(n,des(1,j),1,t(1,nt),1)
   10     continue
c
c     get cross products of x's in null space for m>2
c
c     WARNING: do NOT change next do loop unless you fully understand:
c              This first gets x1*x1, x1*x2, x1*x3, then
c              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3
c              wptr(1) is always at the beginning of the current
c	       level of cross products, hence the end of the
c	       previous level which is used for the next.
c	       wptr(j) is at the start of xj * (previous level)
c
          do 50 k=2,m-1 
             do 40 j=1,dim 
                bptr = wptr(j)
                wptr(j) = nt + 1
                eptr = wptr(1) - 1
                do 30 tt=bptr,eptr 
                   nt = nt + 1
                   do 20 i=1,n
                      t(i,nt) = des(i,j) * t(i,tt)
   20              continue
   30           continue
   40        continue
   50     continue
          if (nt .ne. npoly) then
	      info = 1
	      return
          endif
      endif
c			append s1 to t
      do 60 i = 1,ncov1
         call dcopy(n,s1(1,i),1,t(1,nt+i),1)
   60 continue
      end
Funky!Stuff!
echo yasa: extracting dpdcr.f
cat - << \Funky!Stuff! > dpdcr.f
      subroutine dpdcr(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu, 
     * nlamht,w1,g,pred,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,w1(nnull),g(npsing),pred(nobs),
     * work(*)
c
c Purpose: determine the predicted responses for a given value of 
c	nobs*lamhat and vectors g and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder
c			transformations
c   svals(npsing)	singular values 
c   npsing		number of positive singular values 
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   g(npsing)		(D**2 + nlamht*I)*-1 Dz
c
c On Exit:
c   pred(nobs)		predicted responses
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - dcopy dgemv
c
c $Header: dpdcr.f,v 2.100.1.1 86/10/07 12:51:16 lindstrom Exp $
c
      integer i,locinf,nmh,np1
      double precision dummy(1)
c
c
      np1 = nobs + 1
      nmh = nobs - nnull
c			form the response vector
      call dcopy (nnull,w1,1,pred,1)
      call dcopy (npsing,g,1,work(np1),1)
      do 10 i = 1,npsing
         work(nobs+i) = work(nobs+i)*svals(i)
   10 continue
      call dgemv('N',nmh,npsing,1.0d0,u,ldu,work(np1),1,0.0d0,
     *  pred(nnull+1),1)
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,pred,pred,dummy,dummy,dummy,
     * dummy,10000,locinf)
      return
      end
Funky!Stuff!
echo yasa: extracting dpmse.f
cat - << \Funky!Stuff! > dpmse.f
      subroutine dpmse(fg,ldfg,nuobs,nobs,nnull,fgaux,svals,npsing,u,
     * ldu,w1,z,ntbl,adiag,tbl,ldtbl,auxtbl,work)
      integer ldfg,nuobs,nobs,nnull,npsing,ldu,ntbl,ldtbl
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),w1(nnull),z(npsing),adiag(nuobs),tbl(ldtbl,3),
     * auxtbl(3,3),work(npsing)
c
c Purpose: determine the predictive mean squared error for each lambda 
c	value in tbl.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared
c			in the calling program
c   nuobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder 
c			transformations 
c   svals(npsing)	singular values 
c   npsing		number of singular values
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu			leading dimension of u as declared in the 
c			calling program
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c   ntbl		number of rows in tbl
c   adiag(nuobs)	"true" y values 
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   auxtbl(3,3)		auxiliary table
c			auxtbl(1,1) contains log10(nobs*lamhat) where
c			lamhat is the gcv estimate of lambda
c
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  3     R(lambda) 
c   auxtbl(3,3)		auxiliary table
c			3rd column contains:
c			    [R(lamhat) , R(0), R(infinity)]'
c
c Work Arrays:
c   work(npsing)	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dgemv
c
c $Header: dpmse.f,v 2.100.1.1 86/10/07 12:51:29 lindstrom Exp $
c
      integer i,nmh,k,locinf
      double precision dummy,nlam,wrk1,addtru,wrk
      double precision ddot
c
      nmh = nuobs - nnull
      addtru = 0.0d0
      call dqrsl(fg,ldfg,nuobs,nnull,fgaux,adiag,dummy,adiag,dummy, 
     * dummy,dummy,01000,locinf)
c			the first nnull positions of adiag now contain 
c			w1 true the last nuobs-nnull positions contain 
c			w2 true
      do 10 i = 1,nnull
         addtru=addtru + (w1(i)-adiag(i))**2
   10 continue
      addtru = addtru + ddot(nmh,adiag(nnull+1),1,adiag(nnull+1),1)
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,adiag(nnull+1),1,0.0d0,
     *  work,1)
      addtru = addtru - ddot(npsing,work,1,work,1)
c			addtru contains ||w1 - (w1 true)||**2 +
c			||w2 true||**2 - ||z true||**2
c			work contains z true
c
c			compute predictive mse for each lambda in tbl
      do 30 k = 1,ntbl 
         nlam = 10**tbl(k,1)
         wrk=0.0d0
         do 20 i=1,npsing 
            wrk1=(svals(i)**2)/(svals(i)**2+nlam)
            wrk = wrk + (work(i)-z(i)*wrk1)**2
   20    continue
         tbl(k,3)=(addtru+wrk)/nobs
   30 continue
c			add pred. mse for lambda hat to auxtbl
      wrk=0.0d0
      nlam=10**auxtbl(1,1)
      do 40 i=1,npsing 
         wrk1=(svals(i)**2)/(svals(i)**2+nlam)
         wrk = wrk + (work(i)-z(i)*wrk1)**2
   40 continue
      auxtbl(1,3)=(addtru+wrk)/nobs
c			add pmse for lambda = 0
      wrk=0.0d0
      do 50 i=1,npsing 
         wrk = wrk + (work(i)-z(i))**2
   50 continue
      auxtbl(2,3)=(addtru+wrk)/nobs
c			add pmse for lambda = infinity
      auxtbl(3,3)=(addtru+ddot(npsing,work,1,work,1))/nobs
      return
      end
Funky!Stuff!
echo yasa: extracting dpred.f
cat - << \Funky!Stuff! > dpred.f
      subroutine dpred(pdes,ldpdes,npred,dim,m,desb,lddesb,ndesb,ps,
     * ldps,ncov1,ncov2,coef,npar,pred,work,lwa,iwork,info) 
      integer ldpdes,npred,dim,m,lddesb,ndesb,ldps,ncov1,ncov2,npar,lwa,
     * iwork(dim),info
      double precision pdes(ldpdes,dim),desb(lddesb,dim),ps(ldps,*),
     * coef(npar),pred(npred),work(lwa)
c
c   Purpose: determine predicted values at the locations in pdes and ps
c
c  On Entry:
c   pdes(ldpdes,dim) 	prediction design for splined variables
c   ldpdes		leading dimension of pdes as declared in the 
c			calling	program 
c   npred		number of rows in pdes
c   desb(lddesb,dim) 	locations for the basis functions
c			(returned from dtpss and dptpss in the 
c			variable des)
c   lddesb		leading dimension of desb as declared in the
c			calling	program 
c   ndesb		number of rows in desb
c   dim			number of columns in desb
c   m			order of the derivatives in the penalty
c   ps(ldps,ncov1+ncov2) prediction covariates corresponding to pdes
c   ldps		leading dimension of ps as declared in the
c			calling	program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of pdes
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of pdes
c   coef(npar)		coefficient estimates  [delta':xi']'
c   npar		ndesb + (m+dim-1 choose dim) + ncov1 + ncov2
c
c On Exit:
c   pred(npred)		predicted values
c   info		error indicator
c			  0 : successful completion
c			  1 : dimension error
c			  2 : error in npar,ncov1,ncov2,m or dim
c			  3 : lwa too small
c			  4 : error in dmaket
c			
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work vector
c			must be at least npred*(nct+ndesb)
c			where nct = (m+dim-1 choose dim)
c   iwork(dim)		integer work vector
c
c Subprograms Called Directly:
c    Gcvpack - dmaket dmakek
c    Blas    - dgemv
c
c Subprograms Called Indirectly:
c    Blas    - dcopy
c    Other   - fact mkpoly
c
c $Header: dpred.f,v 2.100.1.1 86/10/07 12:51:40 lindstrom Exp $
c
      double precision dummy
      integer nct,p1,p1p1,npoly
      integer mkpoly 
c
      nct = mkpoly(m,dim) 
      if ((ndesb .le. 0) .or. (nct .le. 0) .or. (m .le. 0) .or. 
     * (dim .le. 0) .or. 2*m - dim .le. 0) then
	 info = 1
	 return
      endif
      if (npar .ne. ndesb + nct + ncov1 + ncov2) then
         info = 2
         return
      endif
      if (lwa .lt. npred*(nct+ndesb)) then
	 info = 3
	 return
      endif
c			first npred*nct positions of work contain t
      p1 = npred*nct
c			next npred*ndesb positions of work contain k
      p1p1 = p1 + 1
c
      call dmaket(m,npred,dim,pdes,ldpdes,dummy,1,0,npoly,work(1),npred,
     * iwork,info)
      if (info .ne. 0) then
         info = 4
         return
      endif
      call dmakek(m,npred,dim,pdes,ldpdes,ndesb,desb,lddesb,work(p1p1),
     * npred)
c			compute predicted values
      call dgemv('N',npred,nct,1.0d0,work,npred,coef,1,0.0d0,
     * pred,1)
      call dgemv('N',npred,ncov1+ncov2,1.0d0,ps,ldps,coef(nct+1),1,
     * 1.0d0,pred,1)
      call dgemv('N',npred,ndesb,1.0d0,work(p1+1),npred,
     * coef(nct+ncov1+ncov2+1),1,1.0d0,pred,1)
      return
      end
Funky!Stuff!
echo yasa: extracting dprmut.f
cat - << \Funky!Stuff! > dprmut.f
      subroutine dprmut (x,npar,jpvt,job)
      integer npar,jpvt(npar),job
      double precision x(npar)
c
c Purpose: permute the elements of the array x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(npar)		array to be permuted
c   npar		size of x (and jpvt)
c   jpvt		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			if job is nonzero backward permutation 
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(npar)		array with permuted entries
c
c   Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dprmut.f,v 2.100.1.1 86/10/07 12:51:58 lindstrom Exp $
c
      integer i,j,k
      double precision t
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               t = x(j)
               x(j) = x(k)
               x(k) = t
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0 ) then
c			backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               t = x(i)
               x(i) = x(j)
               x(j) = t
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
Funky!Stuff!
echo yasa: extracting dptpss.f
cat - << \Funky!Stuff! > dptpss.f
      subroutine dptpss(des,lddes,nobs,dim,m,s,lds,ncov1,ncov2,y,ntbl,
     * adiag,lamlim,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,
     * lwa,iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov1,ncov2,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),adiag(nobs),
     * lamlim(2),dout(4),coef(*),svals(*),tbl(ldtbl,3),
     * auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a partial thin
c	plate spline model.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov1+ncov2)	design for the covariates
c			first ncov1 columns contain covariates which
c			  duplicate the replication structure of des
c			next ncov2 columns contain covariates which
c			  do not duplicate the replication structure of
c			  des
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of des
c   y(nobs)		response vector
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   adiag(nobs)	 	"true" y values on entry if predictive mse is 
c			requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat 
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c On Exit:
c   des(lddes,dim)	unique rows of des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss      residual sum of squares
c			4  tr(I-A)  trace of I - A
c   iout(4)		contains:
c			1  npsing   number of positive singular values
c				    if info indicates nonzero info 
c				    from dsvdc then npsing contains
c				    info as it was returned from dsvdc
c			2  npar	    number of parameters 
c				    (npar = nuobs + nnull)
c			3  nnull    size of the null space of sigma
c				    (m+dim-1 choose dim)+ncov1+ncov2
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least 
c			nuobs+nnull
c   svals(npsing)	singular values, svals must have a dimension,
c			of at least nuobs-nnull.
c			if info indicates nonzero info in dsvdc then 
c			svals is as returned from dsvdc.
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			   0 : successful completion
c			  -1 : log10(nobs*lamhat) <= lamlim(1) 
c			       (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : error in dreps, the first ncov1 columns
c			       of s do not duplicate the replication 
c			       structure of des
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   5 : error in dmaket
c			   6 : sigma is rank deficient
c			  1000< info : 1000 + nonzero info returned from
c				       dsnsm
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least lwa1 + lwa2 where 
c			  lwa1 = (nnull-ncov2)*(nobs+nuobs+1)
c				 +npar*(nobs+npar)
c			  lwa2 = (npar-nnull)*(npar-2*nnull+2+nobs)
c			 	 +npar+nobs
c			
c			
c   iwork(liwa)		integer work vector
c   liwa		length of the iwork as declared in the calling 
c			program
c			must be at least 3*nobs - (nnull - ncov2)
c
c Subprograms Called Directly:
c 	Gcvpack - dreps dmaket duni dmakek dctsx dsnsm
c	Linpack - dqrdc dqrsl 
c	Blas    - dcopy
c	Other 	- dprmut dset prmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz ddcom dgcv dsgdc dtsvdc drsap ddiag
c		  dvlop dvlop dpmse dcfcr dpdcr dvmin dvl dzdc
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset dftkf fact mkpoly
c
c $Header: dptpss.f,v 2.100.1.2 86/11/21 12:22:18 lindstrom Exp $
c
      double precision dummy,dout2(5)
      integer i,ncts1,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,p4,p4p1,p5,p5p1,
     * ip1,ip1p1,ip2,ip2p1,nuobs,lwa2,liwa2,npoly,nnull,snnpar,q1,
     * wsize,lwa1,sinfo,ldx,iout2(4)
      integer mkpoly
c
c
      info = 0
      sinfo = 0
      ncts1 = mkpoly(m,dim) + ncov1
      nnull = ncts1 + ncov2
      iout(3) = nnull
      jadiag = mod(job,10)
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m - dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is an integer work vector
      ip2p1 = ip2 + 1
      liwa2 = liwa - 2*nobs
      call dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,dummy,iwork,nuobs,
     * iwork(ip1p1),0,info)
      iout(2) = nuobs+ncts1+ncov2
      iout(4) = nuobs
      if (info .ne. 0) then
         info = 2
         return
      endif
c			undo permutation of des, s and xrep
      do 10 i = 1,dim
          call dprmut(des(1,i),nobs,iwork,1)
   10 continue
      do 20 i = 1,ncov1+ncov2
          call dprmut(s(1,i),nobs,iwork,1)
   20 continue
      call prmut(iwork(ip1p1),nobs,iwork,1)
c
c			check size of work vectors
      snnpar = nuobs + ncov2
      lwa1 = ncts1*(nobs+nuobs+1)+snnpar*(nobs+snnpar)
      lwa2 = (snnpar-nnull)*(snnpar-2*nnull+2+nobs)+snnpar+nobs
      wsize = lwa1+lwa2
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 3*nobs - ncts1
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for work
c
c	     name     runs from  to		
c	     ----     ---------  --
c	     [t:s1]      1	 p1	  p1 = nobs*ncts1	
c	     [tu:s1u]    p1p1    p2       p2 = p1 + nuobs*ncts1
c	      x	         p2p1    p3       p3 = p2 + max(nobs,snnpar)
c								*snnpar
c	      sigma      p3p1    p4	  p4 = p3 + snnpar**2
c	      fgaux      p4p1    p5	  p5 = p4 + ncts1
c	      working    p5p1    p5+(snnpar-nnull)*(snnpar-2+nnull+
c					2+nobs)+snnpar+nobs
c
c		after the call to duni ku runs from q1 to p4
c
      p1 = nobs*ncts1
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts1
      p2p1 = p2 + 1
      ldx=max(nobs,snnpar)
      p3 = p2 + ldx*snnpar
      p3p1 = p3 + 1
      q1 = p3 + ncov2*snnpar+ncov2+1
      p4 = p3 + snnpar**2
      p4p1 = p4 + 1
      p5 = p4 + ncts1
      p5p1 = p5 + 1
      lwa2 = lwa - ( ncts1*(nobs+nuobs+1) + snnpar*(nobs+snnpar))
c			make [t:s1] and k
      if (m .ne. 1) then
          call dmaket(m,nobs,dim,des,lddes,s,lds,ncov1,npoly,work,nobs,
     *     iwork(ip2p1),info)
          if (info .ne. 0) then
             info = 5
             return
          endif
c			put unique rows of des into des
          call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put k into x
          call dmakek(m,nobs,dim,work(1+nobs),nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
      else
c			put unique rows of des in 1st nuobs 
c			positions of work
          call duni(des,lddes,nobs,dim,iwork(ip1p1),work,1)
c			put k into x
          call dmakek(m,nobs,dim,work,nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
c			copy unique rows of des from work to des 
          call dcopy(nuobs,work,1,des,1)
c			set t = 1
	  call dset(nobs,1.0d0,work,1)
      endif
c			put unique rows of k into (ncov2+1,ncov2+1)
c			position of sigma
      call duni(work(p2p1),ldx,nobs,nuobs,iwork(ip1p1),work(q1),snnpar)
c			compute unique rows of [t:s1] 
      call duni(work(1),nobs,nobs,ncts1,iwork(ip1p1),work(p1p1),nuobs)
c			make x and sigma
      call dctsx(work(1),nobs,ncts1,s(1,ncov1+1),lds,ncov2,nobs,
     * work(p1p1),nuobs, nuobs,work(p3p1),snnpar,work(p2p1),ldx,snnpar,
     * work(p4p1),work(p5p1))
c
      call dsnsm(work(p2p1),ldx,y,work(p3p1),snnpar,nobs,snnpar,nnull,
     * adiag,1.0d0,lamlim,ntbl,dout2,iout2,coef,svals,tbl,ldtbl,auxtbl,
     * iwork(ip2p1),liwa2,work(p5p1),lwa2,job,info)
      call dcopy(4,dout2,1,dout,1)
      iout(3) = nnull
      iout(1) = iout2(1)
      if (info .eq. -3) then
	  info = 6
	  return
      endif
      if (info .gt. 0) then
	 info = info + 1000
	 return
      endif
      if (info .lt. 0) sinfo = info
c			form psi = f2 zeta
      call dset(ncts1,0.0d0,work(p5p1),1)
      call dcopy(nuobs-ncts1,coef(ncts1+ncov2+1),1,work(p5+ncts1+1),1)
      call dqrsl(work(p1p1),nuobs,nuobs,ncts1,work(p4p1),work(p5p1),
     * coef(ncts1+ncov2+1),dummy,dummy,dummy,dummy,10000,info)
      if (sinfo .lt. 0) info = sinfo
      return
      end
Funky!Stuff!
echo yasa: extracting dreps.f
cat - << \Funky!Stuff! > dreps.f
      subroutine dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,c1,order,
     * nuobs,xrep,job,info)
      integer lddes,nobs,dim,lds,ncov1,ncov2,order(nobs),nuobs,
     * xrep(nobs),info,job
      double precision des(lddes,dim),s(lds,*),c1(*)
c
c Purpose: sort des and s, compute degrees of freedom for replication
c	and c1.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c			should be entered in lexicographical order 
c			(smallest to largest) if possible for efficient
c			computing
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   s(lds,ncov1+ncov2) 	design for the covariates
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate
c			replication structure of des
c   job			if job is nonzero then c1 is computed
c   			if job = 0 then c1 is not referenced
c
c On Exit:
c   des(lddes,dim) 	des sorted lexicographically
c   s(lds,ncov1+ncov2) 	s, sorted to correspond to des
c   c1(nuobs)		if job is nonzero then c1(i) the square root of
c			the number of replicates of the ith sorted 
c			design point 
c   order(nobs)		order of the sorted des
c   nuobs		number of unique rows in des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   info		error indicator
c			   0 : successful completion
c			   1 : ncov1 is incorrect
c
c
c $Header: dreps.f,v 2.100.1.1 86/10/07 12:55:24 lindstrom Exp $
c
      integer sw,oldsw,itemp,i,j,k,cont,dfrep
      double precision temp,diff,one,machpr,denom,wmin,wmax
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      do 20 i = 1,nobs
          order(i) = i
 	  xrep(i) = 0
   20 continue
      if (job .ne. 0) call dset(nobs,0.0d0,c1,1)
c			sort des and s 
      sw = nobs - 1
   30 if (sw .le. 0) goto 90 
          oldsw = sw
          sw = 0
          do 80 i = 1,oldsw 
	      cont = 1
	      k = 1
   40         if (cont .eq. 0) goto 70
                  if (k .le. dim) then
	              diff = des(i,k) - des(i+1,k)
                  else
	              diff = s(i,k-dim) - s(i+1,k-dim)
                  endif
		  if (diff .lt. 0.0d0) then
		      if (k .gt. dim) info = 1
		      cont = 0
		  else if (diff .gt. 0.0d0) then
		      if (k .gt. dim) info = 1
c		      switch the order of i and i+1
		      itemp = order(i)
		      order(i)=order(i+1)
		      order(i+1) = itemp
		      itemp = xrep(i)
		      xrep(i)= xrep(i+1)
		      xrep(i+1)= itemp
		      do 50 j = 1,dim
			  temp = des(i,j)
			  des(i,j) = des(i+1,j)
			  des(i+1,j) = temp
   50 		      continue
  		      do 60 j = 1,ncov1+ncov2
			  temp = s(i,j)
			  s(i,j) = s(i+1,j)
			  s(i+1,j) = temp
   60 		      continue
		      sw = i
		      cont = 0
                  else if (k .eq. dim + ncov1) then
		      xrep(i + 1) = 1
		      cont = 0
                  else 
	    	      k = k + 1
                  endif
              goto 40
   70         continue
   80     continue
      goto 30
   90 continue
c			compute range of design
      denom=0.0d0
      do 120 j=1,dim
          wmin = des(1,j)
	  wmax = des(1,j)
          do 110 i=1,nobs
	      if (des(i,j) .lt. wmin) wmin = des(i,j)
	      if (des(i,j) .gt. wmax) wmax = des(i,j)
  110     continue
	  denom = denom + (wmax-wmin)**2
  120 continue
	    
c			check for design points too close together
      do 140 i=1,nobs-1
	  if (xrep(i+1) .eq. 0) then
	     diff = 0.0d0
	     do 130 j=1,dim
	        diff = diff + (des(i,j)-des(i+1,j))**2
  130	     continue
	     if (abs(diff)/denom .lt. 100*machpr) xrep(i+1)=1
	  endif 
  140 continue
c			compute dfrep and c1
      dfrep = 0
      j = 0
       do 150 i = 1,nobs
	   j = j + 1 - xrep(i)
	   if (job .ne. 0) c1(j) = xrep(i)*c1(j) + 1.0d0
	   dfrep = dfrep + xrep(i)
  150 continue
      nuobs = nobs - dfrep
      if (job .eq. 0 ) return
      do 160 i = 1,nuobs
	  c1(i) = sqrt(c1(i))
  160 continue
      return 
      end
Funky!Stuff!
echo yasa: extracting drsap.f
cat - << \Funky!Stuff! > drsap.f
      subroutine drsap(fg,ldfg,nobs,nnull,fgaux,u,ldu,nmh,npsing,z,
     * ssqw2,addend,work)
      integer ldfg,nobs,nnull,ldu,nmh,npsing
      double precision fg(ldfg,nnull),fgaux(nnull),u(ldu,npsing),
     * z(nobs),ssqw2,addend,work(npsing)
c
c Purpose: apply Householder transformations to a response vector and
c	collect its inner product with u and the addend which are used 
c	to define the generalized cross validation function with a
c	semi-norm.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling program
c   nobs		number of rows in fg
c   nnull		number of columns in fg
c   fgaux(nnull)	auxiliary information on the fg	Householder 
c			transformations
c   u(ldu,npsing)		left singular vectors 
c   ldu			leading dimension of u as declared in the
c			calling program
c   nmh	     		number of rows in u. nmh = nobs - nnull
c   npsing		number of columns in u (maximum of npar - nnull)
c   z(nobs)		response vector
c
c On Exit:
c   z(nobs)		the first nnull positions contain w1 and the 
c			next npsing positions contain u'w2
c   ssqw2		the squared length of w2
c   addend		the squared length of z minus the squared length
c			of u'w2
c
c Work Arrays:
c   work(npsing) 	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: drsap.f,v 2.100.1.1 86/10/07 12:55:34 lindstrom Exp $
c
      integer locinf,hp1
      double precision dummy(1)
      double precision ddot
c			apply Householder transformations
c			which define f
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,z,dummy,z,dummy,dummy,dummy,
     * 01000,locinf)
c			w1 in first nnull positions of z,w2 in
c			last nmh
      hp1 = nnull + 1
      addend = ddot(nmh,z(hp1),1,z(hp1),1)
      ssqw2=addend
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,z(hp1),1,0.0d0,work,1)
c			u'w2 in positions nnull+1 to
c			nnull+npsing of z
      call dcopy (npsing,work,1,z(hp1),1)
      addend = addend - ddot(npsing,z(hp1),1,z(hp1),1)
      return
      end
Funky!Stuff!
echo yasa: extracting dset.f
cat - << \Funky!Stuff! > dset.f
      subroutine  dset(n,da,dx,incx)
      integer n,incx
      double precision da,dx(*)
c
c Purpose : set vector dx to constant da. Unrolled loops are used for 
c	increment equal to one.
c
c On Entry:
c   n			length of dx
c   da			any constant
c   incx		increment for dx
c
c On Exit:
c   dx(n)		vector with all n entries set to da
c
c $Header: dset.f,v 2.100.1.1 86/10/07 12:55:42 lindstrom Exp $
c
      integer i,m,mp1,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da
        dx(i + 1) = da
        dx(i + 2) = da
        dx(i + 3) = da
        dx(i + 4) = da
   50 continue
      return
      end
Funky!Stuff!
echo yasa: extracting dsetup.f
cat - << \Funky!Stuff! > dsetup.f
      subroutine dsetup(des,lddes,su1,ldsu1,dim,m,ncov1,nuobs,c1,
     * tusu1,ldtu,ncts1,fgaux,ku,ldku,work,iwork,info)
      integer lddes,ldsu1,dim,m,ncov1,nuobs,ldtu,ncts1,ldku,
     * iwork(dim),info
      double precision des(lddes,dim),su1(ldsu1,*),c1(nuobs),
     * tusu1(ldtu,ncts1),fgaux(ncts1),ku(ldku,nuobs),work(ncts1)
c
c Purpose: set up [tu:su1] as f g and ku as f'c1 ku c1'f.
c
c On Entry:
c   des(lddes,dim)  	variables to be splined (unique rows)
c   lddes		leading dimension of des as declared in the 
c			calling	program 
c   su1(ldsu1,ncov1)	covariates (unique rows)
c   ldsu1		leading dimension of su1 as declared in the
c			calling	program 
c   dim 		dimension of the variables to be splined
c   m			order of the derivatives in the penalty
c   ncov1		number of covariates
c   nuobs		number of unique rows in des
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   ldtu		leading dimension of tusu1 as declared in the 
c			calling	program 
c   ldku		leading dimension of ku as declared in the 
c			calling program 
c
c On Exit:
c   tusu1(ldtu,ncts1)	the qr decomposition of [tu:su1]
c   ncts1		number of columns in [tu:su1] = npoly + ncov1
c   fgaux(ncts1)	the auxiliary info on the qr decomposition of 
c			[tu:su1]
c   ku(p,p)  		f'ku f
c   info		error indicator
c			   0 : successful completion
c			   1 : error in dmaket
c
c Work Arrays:
c   work(ncts1)		double precision work vector
c   iwork(dim)		integer work vector
c
c Subroutines called directly
c	Gcvpack - dmaket dmakek
c	Linpack - dqrdc
c	Other   - dftkf mkpoly
c
c Subroutines called indirectly
c	Blas    - dcopy
c	Gcvpack - dqrsl
c	Other   - fact mkpoly
c
c $Header: dsetup.f,v 2.100.1.1 86/10/07 12:55:48 lindstrom Exp $
c
      integer npoly,i,j
      integer mkpoly
c
      info = 0
      npoly=mkpoly(m,dim)
      ncts1=npoly+ncov1
c			make [tu:su1] and ku
      call dmaket(m,nuobs,dim,des,lddes,su1,ldsu1,ncov1,npoly,tusu1,
     * ldtu,iwork,info)
      if (info .ne. 0) then
	 return
      endif
      call dmakek(m,nuobs,dim,des,lddes,nuobs,des,lddes,ku,ldku)
      if (c1(1) .ne. 0) then
         do 30 i = 1,nuobs 
            do 10 j = 1,npoly+ncov1
               tusu1(i,j) = tusu1(i,j) * c1(i)
   10       continue
            do 20 j = 1,nuobs
               ku(i,j) = ku(i,j) * c1(i) * c1(j)
   20       continue
   30    continue
      endif
c			decompose [tu:su1] into fg
      call dqrdc(tusu1,ldtu,nuobs,ncts1,fgaux,0,0.0d0,0)
c      			calculate f'ku f
      call dftkf(tusu1,ldtu,nuobs,ncts1,fgaux,ku,ldku,work)
      return
      end
Funky!Stuff!
echo yasa: extracting dsgdc.f
cat - << \Funky!Stuff! > dsgdc.f
      subroutine dsgdc (sigma,ldsigm,npar,nnull,qraux,sgpvt,info)
      integer ldsigm,npar,nnull,sgpvt(npar),info
      double precision sigma(ldsigm,npar),qraux(npar)
c
c Purpose: decompose the semi-norm smoothing matrix into a QR 
c	decomposition of the transpose of the Cholesky factor.
c
c On Entry:
c   sigma(ldsigm,npar)	symmetric matrix which defines the semi-norm
c   ldsigm		leading dimension of sigma as declared in the
c			calling program
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c
c On Exit:
c   sigma(ldsigm,npar)	overwritten with the QR decomposition
c			of the Cholesky factor of sigma
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   qraux(npar)		auxiliary information for the QR decomposition
c   sgpvt(npar) 	permuted indices from the Cholesky decomposition
c			with pivots of sigma
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : nnull is too large
c
c Subprograms Called Directly:
c	Linpack - dchdc dqrdc  
c	Blas    - dcopy
c	Other   - dset
c
c $Header: dsgdc.f,v 2.100.1.1 86/10/07 12:56:05 lindstrom Exp $
c
      integer locinf,i,j,nsm,idummy
      double precision dummy,machpr,one
c
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      locinf = 0
c			Cholesky decomposition of sigma
      do 20 j = 1,npar
         sgpvt(j) = 0
   20 continue
      call dchdc (sigma,ldsigm,npar,qraux,sgpvt,1,locinf)
      do 30 i=1,locinf 
         if (((sigma(i,i)/sigma(1,1))**2) .gt. machpr) nsm = i
   30 continue
      if (nsm .gt. npar - nnull) then
         info = 1
         return
      endif
      if (nsm .lt. npar - nnull) then
         nnull = npar - nsm
         info = -3
      endif
c			copy transpose of Cholesky factor to sigma
      do 40 i = 1,nsm 
         call dcopy (npar,sigma(i,1),ldsigm,sigma(1,i),1)
         j = i - 1
         call dset (j,0.0d0,sigma(1,i),1)
   40 continue
c			QR decomposition of Cholesky transpose
      call dqrdc (sigma,ldsigm,npar,nsm,qraux,idummy,dummy,0)
      return
      end
Funky!Stuff!
echo yasa: extracting dsgdc1.f
cat - << \Funky!Stuff! > dsgdc1.f
      subroutine dsgdc1(f2kf2,ldfkf,p,svals,iwork,work,lwa,info)
      integer ldfkf,p,iwork(p),lwa,info
      double precision f2kf2(ldfkf,p),svals(*),work(lwa)
c 
c Purpose: form the singular value decomposition of the Cholesky factor 
c 	of f2'k f2.
c
c On Entry:
c   f2kf2(ldfkf,p)	f2'k f2
c   ldfkf		leading dimension of f2'k f2 as declared in the 
c			calling	program 
c   p			number of rows and columns in f2'k f2
c
c On Exit:
c   f2kf2(p,p)  	overwritten with singular value decomposition 
c 			of Cholesky factor of f2'k f2
c   svals(p) 		the singular values of the Cholesky factor of 
c			f2'k f2 if info = 0.
c			if info = 3 then svals is as it was returned 
c			from dsvdc.
c   info 	   	error indicator
c			  0 : successful completion
c			  1 : lwa too small
c			  2 : f2'k f2 is not of full rank 
c			  3 : error in dsvdc
c   p			if info = 3 p contains info as it was returned 
c			from dsvdc (otherwise unchanged)
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling
c			program
c			must be at least 2*p
c   iwork(p)		integer work vector
c
c Subprograms Called:
c	Linpack - dchdc dsvdc
c	Blas    - dcopy
c	Other   - dset dprmut
c
c $Header: dsgdc1.f,v 2.100.1.1 86/10/07 12:56:12 lindstrom Exp $
c
      integer i,j,pp1,locinf,k
      double precision dummy,one,machpr
c
      info = 0
      if (lwa .lt. 2*p) then
	  info = 1
	  return
      endif
      pp1 = p + 1
      call dset(p,0.0d0,svals,1)
c
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			Cholesky decomposition of f2'k f2
      do 20 j = 1,p
         iwork(j) = 0
   20 continue
      call dchdc (f2kf2,ldfkf,p,work,iwork,1,locinf)
      do 30 i=1,locinf 
         if ((f2kf2(i,i)/f2kf2(1,1))**2 .gt. machpr) k = i
   30 continue
      if (k .lt. p) then
	 info = 2
	 return
      endif
c    		copy f2kf2' into f2kf2 
c    		svd of f2 k f2' = udv' return only u
      do 40 j = 1,p 
         call dcopy(1+p-j,f2kf2(j,j),ldfkf,f2kf2(j,j),1)
         call dset(j-1,0.0d0,f2kf2(1,j),1)
   40 continue
      call dsvdc(f2kf2,ldfkf,p,p,svals,work,f2kf2,ldfkf,dummy,1,
     * work(pp1),20,info)
      if (info .ne. 0) then
	  p = info
	  info = 3
	  return
      endif
      do 50 j=1,p
         call dprmut(f2kf2(1,j),p,iwork,1)
   50 continue
      return
      end
Funky!Stuff!

