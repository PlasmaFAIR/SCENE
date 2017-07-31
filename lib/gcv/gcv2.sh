#!/bin/sh
#    This is a shell archive.
#    Run the following text with /bin/sh to extract.

YASASTART=`pwd`

echo yasa: extracting dsnsm.f
cat - << \Funky!Stuff! > dsnsm.f
      subroutine dsnsm (x,ldx,y,sigma,ldsigm,nobs,npar,nnull,adiag,
     * tau,lamlim,ntbl,dout,iout,coef,svals,tbl,ldtbl,auxtbl,
     * iwork,liwa,work,lwa,job,info)
      integer ldx,ldsigm,nobs,npar,nnull,ntbl,iout(3),ldtbl,liwa,
     * iwork(liwa),lwa,job,info
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * adiag(nobs),tau,lamlim(2),dout(5),coef(npar),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized 
c	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as declared
c			in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   adiag(nobs)		"true" y values on entry if computation of 
c			predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested (try tau = 1 
c			to start then try 10 and 100)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abcd
c			if a is nonzero then truncation is used
c			if b is nonzero then predictive mse is computed
c			   using adiag as true y
c			if c is nonzero then user input limits on search
c			   for lambda hat are used
c			if d is nonzero then the diagonal of the hat 
c			   matrix is calculated
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   y(nobs)		predicted values
c   sigma(ldsigm,npar)	overwritten with the QR decomposition of the 
c			Cholesky factor of sigma
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1 lamhat   generalized cross validation 
c				   estimate of the smoothing parameter
c			2 penlty   smoothing penalty
c			3 rss	   residual sum of squares
c			4 tr(I-A)  trace of I - A
c			5 truncation ratio = 1/(1+(normk/(nobs*lamhat)))
c				   where normk = norm(R - R sub k)**2
c   iout(3)		contains:
c			1  npsing   number of positive singular
c				    values
c				    if info indicates nonzero info in 
c				    dsvdc then iout(1) contains info as
c				    it was returned from dsvdc
c			2  npar	    number of parameters
c			3  nnull    size of the null space of sigma
c   coef(npar)		coefficient estimates
c   svals(npar-nnull)	first npsing entries contain singular values 
c			of the matrix j2 
c			if info indicates nonzero info in dsvdc then 
c			svals is as it was returned from dsvdc
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
c			  -3 : nnull is too small (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			  -1 : log10(nobs*lamhat) <= lamlim(1)
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : lwa (length of work) is too small
c			   3 : liwa (length of iwork) is too small
c			   4 : error in ntbl or tau
c			  100< info <200 : 100 + nonzero info returned
c					   from ddcom
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least 
c			(npar-nnull)*(npar-2*nnull+2+nobs)+npar+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			must be at least 2*npar - nnull 
c
c Subprograms Called Directly:
c	Gcvpack - ddcom dgcv
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz dsgdc dcfcr drsap dvlop dtsvdc
c		  dpmse dvmin dvl dzdc dpdcr ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset 
c
c $Header: dsnsm.f,v 2.100.1.2 86/11/19 09:24:39 lindstrom Exp $
c
      integer bigp,pmh,pp1,lwa2,wsize,ldcaux,job1,job2,npsing,sinfo
      sinfo = 0
      info = 0
      iout(2) = npar
c			check dimensions
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or.(pmh .le. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 +2*npar - nnull + pmh*(nobs - nnull) + pmh + nobs
      if (lwa .lt. wsize) then
         info = 2
         return
      endif
      wsize = npar + pmh
      if (liwa .lt. wsize) then
         info = 3
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. (tau .lt. 0)) then
         info = 4
         return
      endif
c 			first p positions of iwork contain sgpvt
c 			first bigp-1  positions of work contain dcaux 
c			iwork vector starts at pp1
      pp1 = npar + 1
c 	work vector starts at bigp
      bigp = pmh**2 + 2*npar - nnull + 1
      ldcaux = pmh**2 + 2*npar - nnull
      lwa2 = lwa - ldcaux
      job1=job/1000
      job2=mod(job,1000)
c			decompose sigma and design matrix
      call ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,npsing,
     * svals,iwork,work,ldcaux,dout(5),work(bigp),lwa2,iwork(pp1),pmh,
     * job1,info)
      iout(1) = npsing
      iout(3) = nnull
      if (info .gt. 0) then
         info = info + 100
         return
      endif
      if (info .lt. 0) sinfo = info
c			compute lambda hat and other parameters
      call dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,iwork,work,ldcaux,
     * adiag, lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work(bigp),lwa2,
     * job2,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
      dout(5) = 1.0d0/(1.0d0+dout(5)/nobs/dout(1))
      if (sinfo .lt. 0) info = sinfo
      return
      end
Funky!Stuff!
echo yasa: extracting dsuy.f
cat - << \Funky!Stuff! > dsuy.f
      subroutine dsuy(y,nobs,nuobs,ytrue,c1,order,xrep,ssqrep,
     * job)
      integer nobs,nuobs,order(nobs),xrep(nobs),job
      double precision y(nobs),ytrue(nobs),c1(nuobs),ssqrep
c
c Purpose: compute B1'y, B1'ytrue and ssq for replication.
c
c On Entry:
c   y(nobs)  		response vector
c   nobs		number of observations
c   nuobs		number of unique design points
c   ytrue(nobs)		"true" response, if job is nonzero 1
c			not referenced if job = 0
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   order(nobs)		order of sorted des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   job 		job is nonzero if B1'ytrue should be calculated 
c			job = 0 otherwise
c
c On Exit:
c   y(nuobs)  		B1'y
c   ytrue(nuobs)	B1'ytrue if job is nonzero
c   ssqrep		sum of squares for replication
c
c $Header: dsuy.f,v 2.100.1.1 86/10/07 12:57:32 lindstrom Exp $
c
      integer first,i,j
      double precision accum
c
      accum = 0.0d0
      ssqrep = 0.0d0
      call dprmut (y,nobs,order,0)
      if (job .ne. 0) call dprmut (ytrue,nobs,order,0)
c			compute ssq for replication
      first = 0
      do 20 i = 1,nobs
	  if (xrep(i) .eq. 1) then
	      accum = accum + y(i)
	  else if (first .eq. i - 1) then
	      first = i
	      accum = y(i)
	  else
	      accum = accum/(i-first)
	      do 10 j = first,i-1
	          ssqrep = (y(j)-accum)**2 + ssqrep
   10         continue
	      first = i
	      accum = y(i)
          endif
   20 continue
      if (xrep(nobs) .eq. 1) then
          accum = accum/(nobs + 1 - first)
          do 30 j = first,nobs
              ssqrep = (y(j)-accum)**2 + ssqrep
   30     continue
      endif
c			compute B1'y and B1'ytrue
      j = 0
      do 40 i = 1,nobs
 	  if (xrep(i) .eq. 0) then
 	      if (j .ne. 0) then
 		  y(j) = y(j) / c1(j)
 		  if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
              endif
 	      j = j + 1
	      y(j) = y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(i)
	  else
	      y(j) = y(j) + y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(j) + ytrue(i)
 	  endif
   40 continue
      y(j) = y(j) / c1(j)
      if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
      return 
      end
Funky!Stuff!
echo yasa: extracting dtpss.f
cat - << \Funky!Stuff! > dtpss.f
      subroutine dtpss(des,lddes,nobs,dim,m,s,lds,ncov,y,ntbl,adiag,
     * lamlim,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,lwa,
     * iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),
     * adiag(nobs),lamlim(2),dout(5),coef(*),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c 	smoothing parameter and fit model parameters for a thin plate
c 	smoothing spline.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in calling
c   			program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov) 	design for the covariates. The covariates
c			must duplicate the replication structure of des.
c			See dptpss to handle covariates which do not.
c   lds			leading dimension of s as declared in calling
c   			program 
c   ncov		number of covariates 
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
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abdc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c			if d is nonzero then there are replicates in the
c			   design
c
c On Exit:
c   des(lddes,dim)	sorted unique rows of des if job indicates that
c			there are replicates otherwise not changed
c   s(lds,ncov)		unique rows of s sorted to correspond to des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   			5  ssqrep   sum of squares for replication
c   iout(4)		contains:
c			1  npsing   number of positive singular
c				    values (npsing = nuobs - ncts).
c				    if info indicates nonzero info in 
c				    dsvdc then npsing contains info as 
c				    it was returned from dsvdc.
c			2  npar	    number of parameters
c				    (npar = nuobs + ncts)
c			3  ncts     dimension of the polynomial space 
c				    plus ncov	
c				    ((m+dim-1 choose dim) + ncov)
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least nuobs+
c			ncts
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info indicates nonzero info from dsvdc then 
c			svals is as it was returned from dsvdc.
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
c			  2 : error in dreps, covariates do not 
c			      duplicate the replication structure of des
c			  3 : lwa (length of work) is too small
c			  4 : liwa (length of iwork) is too small
c			  10 < info < 20 : 10 + nonzero info returned 
c					   from dsetup
c			  100< info <200 : 100 + nonzero info returned
c					   from dsgdc1
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv1
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			Must be at least nuobs(2+ncts+nuobs)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			Must be at least 2*nobs + nuobs - ncts
c
c Subprograms Called Directly:
c       Gcvpack - dreps duni dsuy dsetup dsgdc1 dgcv1 
c       Other   - dprmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcfcr1 drsap dvlop dsvtc dpdcr dpmse 
c		  dvmin dvl dmaket dmakek ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - ddot dcopy dgemv
c	Other 	- dprmut dset dftkf fact mkpoly 
c
c $Header: dtpss.f,v 2.100.1.1 86/10/07 12:57:52 lindstrom Exp $
c
      integer ncts,jrep,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,ip1,ip1p1,ip2,
     * ip2p1,nuobs,p4,p4p1,j,i,lwa2,wsize,q1,q1p1,npsing,
     * jobgcv
      double precision ssqrep
      integer mkpoly
c
c
      info = 0
      ncts = mkpoly(m,dim) + ncov
      iout(3) = ncts
      jrep = mod(job,10)
      jadiag = mod(job,100)/10
      jobgcv = job/10
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m-dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork vector
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep    
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is a integer work vector
      ip2p1 = ip2 + 1
      if (jrep .ne. 0) then
         call dreps(des,lddes,nobs,dim,s,lds,ncov,0,work,iwork,nuobs,
     *    iwork(ip1p1),1,info)
         if (info .ne. 0) then
            info = 2
            return
         endif
c			put unique row of des into des
         call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put unique row of s into s
         call duni(s,lds,nobs,ncov,iwork(ip1p1),s,lds)
      endif
      if (jrep .eq. 0) then
         work(1)= 0
         nuobs = nobs
         ssqrep = 0.0d0
      endif
      iout(2) = nuobs + ncts
      iout(4) = nuobs
c			check size of work vectors
      wsize = nuobs*(2+ncts+nuobs)+nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 2*nobs + nuobs - ncts
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for c1,[tu:su1],ku,fgaux in work
c
c		c1  	 runs from 1    to p1,   p1 = nuobs
c		[tu:su1] runs from p1p1 to p2,   p2 = p1 + nuobs*ncts
c		ku	 runs from p2p1 to p3,   p3 = p2 + nuobs**2
c		fgaux    runs from p3p1 to p4,   p4 = p3 + ncts
c		the rest of work is a work vector

c  			after the call to dsetup 
c		f2'k f2  runs from q1p1 to p3, q1 = p2 + nuobs*ncts+ncts
c
      p1 = nuobs
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts
      p2p1 = p2 + 1
      q1 = p2 + nuobs*ncts+ncts
      q1p1 = q1 + 1
      p3 = p2 + nuobs**2
      p3p1 = p3 + 1
      p4 = p3 + ncts
      p4p1 = p4 + 1
      lwa2 = lwa - (nuobs*(1+ncts+nuobs) + ncts)
c			set up structures needed for dsgdc1 and dgcv1
      call dsetup(des,lddes,s,lds,dim,m,ncov,nuobs,work,work(p1p1),
     * nuobs,ncts,work(p3p1),work(p2p1),nuobs,work(p4p1),iwork(ip2p1),
     * info)
      if (info .ne. 0) then
         info = info + 10
         return
      endif

c			decompose f2' k f2 
      npsing = nuobs - ncts
      call dsgdc1(work(q1p1),nuobs,npsing,svals,iwork(ip2p1),
     * work(p4p1),lwa2,info)
      iout(1) = npsing
      if (info .gt. 0) then
         info = info + 100
         return
      endif
c			setup y
      if (jrep .ne. 0) then
         call dsuy(y,nobs,nuobs,adiag,work,iwork,iwork(ip1p1),
     *    ssqrep,jadiag)
      endif
      dout(5) = ssqrep
c			compute lambda hat and other parameters
      call dgcv1(work(p2p1),nuobs,y,nuobs,nobs,work(p1p1),nuobs,ncts,
     * work(p3p1),svals,adiag,lamlim,ssqrep,ntbl,dout,coef,tbl,ldtbl,
     * auxtbl,work(p4p1),lwa2,jobgcv,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
c			if there are replicates then rescale the coef. 
c			vector,	the predicted values, and diagonal of A
      if (nuobs .ne. nobs) then
         do 10 i = 1,nuobs
            coef(ncts+i) = coef(ncts+i) * work(i)
   10    continue
         j = nuobs
         do 20 i = nobs,1,-1 
            y(i)=y(j)/work(j)
            if (jadiag .ne. 0) then
               adiag(i) = adiag(j)/work(j)**2
            endif
            if (iwork(ip1 + i) .eq. 0) then
               j = j - 1
            endif
   20    continue
      endif
c			undo the permutation of the predicted values and
c			adiag if necessary
      if (jrep .ne. 0) then
         call dprmut (y,nobs,iwork,1)
         if (jadiag .ne. 0) call dprmut (adiag,nobs,iwork,1)
      endif
      return
      end
Funky!Stuff!
echo yasa: extracting dtsvdc.f
cat - << \Funky!Stuff! > dtsvdc.f
      subroutine dtsvdc(x,ldx,n,p,minrat,k,s,u,ldu,v,ldv,normk,
     * work,lwa,iwork,liwa,job,info)
      integer ldx,n,p,k,ldu,ldv,lwa,iwork(liwa),liwa,job,info
      double precision x(ldx,p),minrat,s(p),u(ldu,*),v(ldv,p),
     * normk,work(lwa)
c
c Purpose: form the singular value decomposition of a truncated matrix 
c	obtained by first taking a qr decomposition with pivoting.
c
c On Entry:
c   x(ldx,p)   		matrix to be decomposed, x is destroyed
c   ldx     		leading dimension of x in the calling program
c   n       		number of rows in x
c   p       		number of columns in x   (p <= n)
c   minrat  		minimum ratio to determine truncation
c   ldu     		leading dimension of u as declared in the 
c			calling program
c   ldv     		leading dimension of v as declared in the 
c			calling program
c   job     		controls the computation of the singular vectors
c               	it has the decimal expansion ab with the meaning
c                            a = 0 no left singular vectors
c                            a = 1 all n left singular vectors
c                            a > 1 first p left singular vectors
c                            b = 0 no right singular vectors
c                            b > 0 all p right singular vectors
c
c On Exit:
c   k       		number of positive singular values 
c   s(p)		an approximation to the first k singular values
c			in decreasing order
c   u(ldu,m)   		first k singular vectors if job > 1 or all n 
c			left singular vectors if joba = 1 otherwise it 
c			is not accessed it may be identified with x in 
c			the subroutine call if joba > 1
c   v(ldv,p)   		the first k right singular vectors if job b > 0
c               	otherwise it is not accessed
c   normk		norm of the k by k lower right sub matrix of r
c   info    		error indicator
c		   	   0 : successful completion
c		  	   info > 0 : info was returned from dsvdc
c
c Work Arrays:
c   iwork(liwa) 	integer work vector, holds the pivot vector from
c			dqrdc
c   liwa		must be at least p
c   work(lwa)		work vector
c   lwa			must be at least p**2 + p + n if joba > 0, 
c			otherwise it must be at least 2*p
c
c Subprograms Called Directly
c	Linpack - dsvdc dqrdc dqrsl
c   	Blas    - ddot dcopy
c	Other   - dset dprmut
c
c $Header: dtsvdc.f,v 2.100.1.1 86/10/07 12:58:56 lindstrom Exp $
c
      integer i,j,jj,pp1,ppnp1,nmk,sjob,locinf
      double precision trxpx,accum,dummy(1),mintr
      double precision ddot
c
      info = 0
      call dset(p,0.0d0,s,1)
      pp1 = p+1
      ppnp1 = p+n+1
c                       calculate trace of x' x
      trxpx = 0.0d0
      do 10 j = 1,p
         trxpx = trxpx+ddot(n,x(1,j),1,x(1,j),1)
   10 continue
      mintr = trxpx * minrat
c                       qr decomposition of x
      do 20 j = 1,p
         iwork(j) = 0
   20 continue
      call dqrdc(x,ldx,n,p,work,iwork,work(pp1),1)
c                       calculate ratios for the truncated matrix
      accum = 0.0d0
      k = p
      do 30 i = p,1,-1
         accum = accum+ddot(pp1-i,x(i,i),ldx,x(i,i),ldx)
         if (accum .lt. mintr) then
	     k = i
	     normk = accum
         endif
   30 continue
      if (job/10.le.0) then
c                       no left singular vectors
c                       copy rk' (transpose of the first k rows of r
c			from qr decomposition) into x
         do 40 j = 1,k 
            call dcopy(p,x(j,1),ldx,x(1,j),1)
            call dset(j-1,0.0d0,x(1,j),1)
   40    continue
c                       svd of rk'
         sjob = 0
         if (mod(job,10).ne.0) sjob = 20
         call dsvdc(x,ldx,p,k,s,work,v,ldv,dummy,0,work(pp1),sjob,info)
         if (info.ne.0) return
         if (mod(job,10).ne.0) then
            do 50 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
   50       continue
         endif
      else
c                       u or u1 is to be created
         i = n+1
c                       copy rk' to work
         do 60 j = 1,k 
            i = i+p
            call dcopy(p,x(j,1),ldx,work(i),1)
            call dset(j-1,0.0d0,work(i),1)
   60    continue
c                       create u2 if requested
         if (job/10.eq.1) then
            nmk = n-k
            do 70 i = 1,nmk 
               j = n+1-i
               call dset(n,0.0d0,work(pp1),1)
               work(p+j) = 1.0d0
               call dqrsl(x,ldx,n,min(j,p),work,work(pp1),work(pp1),
     *          dummy,dummy,dummy,dummy,10000,locinf)
               call dcopy(n,work(pp1),1,u(1,j),1)
   70       continue
         endif
         do 80 j = 1,k 
            jj = k+1-j
            call dset(n,0.0d0,work(pp1),1)
            i = p+jj
            work(i) = 1.0d0
            call dqrsl(x,ldx,n,jj,work,work(pp1),work(pp1),dummy,dummy,
     *       dummy,dummy,10000,locinf)
            call dcopy(n,work(pp1),1,u(1,jj),1)
   80    continue
c			  svd of rk'
         sjob = 1
         if (mod(job,10).ne.0) then
            sjob = 21
         endif
         call dsvdc(work(ppnp1),p,p,k,s,work,v,ldv,work(ppnp1),p,
     *    work(pp1),sjob,info)
         if (info.ne.0) return
         do 100 i = 1,n 
            call dcopy(k,u(i,1),ldu,work,1)
            jj = n+1
            do 90 j = 1,k 
               jj = jj+p
               u(i,j) = ddot(k,work,1,work(jj),1)
   90       continue
  100    continue
         if (mod(job,10).ne.0) then
c		 undo pivots on right singular vectors
            do 110 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
  110       continue
         endif
      endif
      end
Funky!Stuff!
echo yasa: extracting duni.f
cat - << \Funky!Stuff! > duni.f
      subroutine duni(x,ldx,nobs,ncx,xrep,xu,ldxu)
      integer ldx,nobs,ncx,ldxu,xrep(nobs)
      double precision x(ldx,*),xu(ldxu,*)
c
c Purpose: compute xu.
c
c On Entry:
c   x(ldx,ncx)		a matrix to be reduced to unique rows
c   ldx			leading dimension of x as declared in the
c			calling	program 
c   nobs		number of observations
c   ncx			number of columns in x
c   xrep(nobs)		xrep(i) contains 1 if ith row of x is a 
c			replicate row, 0 if not
c   ldxu		leading dimension of xu as declared in the 
c			calling	program 
c On Exit:
c   xu(ldxu,ncx) 	unique rows of x
c			may be identified with x in the calling sequence
c
c $Header: duni.f,v 2.100.1.1 86/10/07 12:59:15 lindstrom Exp $
c
      integer i,j,k
c
      j = 0
      do  20 i = 1,nobs
 	  if (xrep(i) .eq. 0) then
 	      j = j + 1
 	      do 10 k = 1,ncx
 	          xu(j,k) = x(i,k)
   10         continue
 	  endif
   20 continue
      return 
      end
Funky!Stuff!
echo yasa: extracting dvl.f
cat - << \Funky!Stuff! > dvl.f
      double precision function dvl(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block 
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
   10 continue
      rss=numrtr
      tria=denom
      dvl=dble(n)*numrtr/denom**2
      return
      end
Funky!Stuff!
echo yasa: extracting dvlop.f
cat - << \Funky!Stuff! > dvlop.f
      subroutine dvlop(z,svals,nobs,nnull,npsing,inadd,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout,job,info)
      integer nobs,nnull,npsing,ntbl,ldtbl,job,info
      double precision z(npsing),svals(npsing),inadd,ssqw2,lamlim(2),
     * nlamht,tbl(ldtbl,3),auxtbl(3,3),dout(2)
c
c Purpose: determine the optimal lambda for the generalized cross 
c	validation function given singular values and the data vector 
c	in canonical coordinates.
c
c On Entry:
c   z(npsing)		data vector in canonical coordinates
c   svals(npsing)	singular values 
c   nobs		number of observations
c   nnull		dimension of the null space of sigma
c   npsing		number of positive elements of svals 
c   inadd		constant term in expression for V
c   ssqw2		squared length of w2
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then nlamht
c			is set to 10**lamlim(1)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job      		if job is nonzero then user input limits on 
c			lambda hat search are used
c
c On Exit:
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   nlamht		nobs*(lambda hat) where lambda hat is the gcv 
c			estimate of lambda
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lambda hat), V(lambda hat)
c			2nd row contains:
c			    0, V(0) 
c			3rd row contains:
c			    0, V(infinity) 
c   dout(2)		contains:
c			1  rss
c			2  tr(I-A)
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nlamht) <= lamlim(1) (not fatal)
c			 -2 : log10(nlamht) >= lamlim(2) (not fatal)
c			  1 : svals(1) = 0.0d0
c			  2 : npsing is incorrect
c			  3 : lamlim(1) > lamlim(2)
c
c Subprograms Called Directly:
c	Gcvpack - dvmin 
c
c Subprograms Called Indirectly:
c	Gcvpack - dvl
c
c $Header: dvlop.f,v 2.100.1.1 86/10/07 12:59:38 lindstrom Exp $
c
      integer i,k
      double precision vlamht,w
      double precision dvmin
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision addend,rss,tria,machpr,one
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      n=nobs
      h=nnull
      addend = inadd
      if (svals(1) .eq. 0.0d0) then
         info = 1
         return
      endif
      k = 0
      do 20 i = 1,npsing
         if (svals(i) .gt. 0) then
            k = i
         endif
   20 continue
      if (k .ne. npsing) then
	 info = 2
    	 return
      endif
      if (job .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 3
         return
      endif
      if (job .eq. 0) then
         lamlim(2) = 2.0d0*dlog10(svals(1))+2.0d0
         lamlim(1) = 2.0d0*dlog10(svals(npsing))-2.0d0
      endif
      nlamht = dvmin (lamlim(1),lamlim(2),svals,z,npsing,ntbl,tbl,
     * ldtbl,vlamht,info)
      dout(1) = rss
      dout(2) = tria
c			compute auxtbl
      auxtbl(1,1)=nlamht
      auxtbl(1,2)=vlamht
c			lambda = 0
      auxtbl(2,1)=0.0d0
      auxtbl(2,2)=0.0d0
      if ((nobs-nnull) .ne. npsing) then
         auxtbl(2,2)=inadd*(nobs)/(nobs-nnull-npsing)**2
      endif
      if ((nobs-nnull) .eq. npsing) then
         w=0.0d0
         do 30 i=npsing,1,-1
c           w=w+(z(i)*svals(npsing)**2/(svals(i)**2))**2
            w=w+(z(i)*(svals(npsing)/svals(i))**2)**2
   30    continue
         auxtbl(2,2)=nobs*w
         w=0.0d0
         do 40 i=npsing,1,-1
            w=w+(svals(npsing)/svals(i))**2
   40    continue
         auxtbl(2,2)=auxtbl(2,2)/(w**2)
      endif
c			lambda = infinity
      auxtbl(3,1)=0.0d0
      auxtbl(3,2)=ssqw2/(nobs - nnull)
      nlamht = 10**nlamht
      return
      end
Funky!Stuff!
echo yasa: extracting dvmin.f
cat - << \Funky!Stuff! > dvmin.f
      double precision function dvmin(lower,upper,svals,z,npsing,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl
c
c $Header: dvmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl
c				null interval
      if (lower .eq. upper) then
	 dvmin = lower
	 info = -1
	 vlamht = dvl(lower,svals,z,npsing)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl(x,svals,z,npsing)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl(lower,svals,z,npsing)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl(tbl(j,1),svals,z,npsing)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl(c,svals,z,npsing)
      vd = dvl(d,svals,z,npsing)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl(d,svals,z,npsing)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl(c,svals,z,npsing)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl(x,svals,z,npsing)
      dvmin = x
      return
      end
Funky!Stuff!
echo yasa: extracting dzdc.f
cat - << \Funky!Stuff! > dzdc.f
      subroutine dzdc(x,ldx,nobs,npar,pmh,tau,fgaux,svals,npsing,
     * v,ldv,normk,work,lwa,iwork,liwa,job,info)
      integer ldx,nobs,npar,pmh,npsing,ldv,lwa,iwork(liwa),liwa,info,job
      double precision x(ldx,npar),tau,fgaux(*),svals(*),v(ldv,pmh),
     * normk,work(lwa)
c
c Purpose: create and decompose j and decompose j2
c
c On Entry:
c   x(ldx,npar)		z = design matrix after sigma rotations 
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   tau			multiplier controlling the amount of truncation 
c			if truncation is requested
c   ldv    		leading dimension of v as declared in the 
c			calling	program
c   job			controls use of truncated singular value 
c			decomposition
c			   if job is nonzero then dtsvdc is used
c			   if job = 0 then dsvdc is called directly
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   fgaux(npar-pmh)	auxiliary information on the qr decomposition 
c			of the matrix z2
c   svals(npsing)      	singular values of the matrix j2 if info = 0.
c			if info = 4,5,6 or 7 then svals is as 
c			it was returned from dsvdc.
c   npsing       	if info = 0 then npsing contains the number 
c			of singular values calculated 
c			if info = 4,5,6 or 7 then npsing contains info
c			as it was returned from dsvdc.
c   v(ldv,pmh) 		right singular vectors of the matrix j2
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info    		error indicator
c			  0 : successful completion
c			  1 : lwa (length of work) is too small
c			  2 : tau < 0
c			  3 : transpose of j2 is necessary 
c				(npar .gt. nobs) and npar .gt. ldx
c			  4 : error in dtsvdc (using j2')
c			  5 : error in dsvdc (using j2')
c			  6 : error in dtsvdc
c			  7 : error in dsvdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa     		length of work as declared in the calling 
c			program	Must be at least 
c			pmh*(nobs - (npar-pmh) ) + pmh + nobs
c   iwork(liwa)		integer work vector
c   liwa		length of integer work vector, must be at least
c			pmh
c
c Subprograms Called Directly:
c	Gcvpack - dtsvdc 
c	Linpack - dqrdc dqrsl dsvdc
c	Blas    - dcopy dswap
c
c Subprograms Called Indirectly:
c	Linpack - dqrdc dsvdc dqrsl
c	Blas    - ddot dcopy
c	Other   - dset dprmut dswap
c
c $Header: dzdc.f,v 2.100.1.1 86/10/07 13:00:03 lindstrom Exp $
c
      double precision dummy,machpr,one,minrat
      integer i,j,idummy,nnull,hp1,nmh,pmhp1,wsize,locinf
c
      info = 0
c
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
      minrat = machpr*tau
c			check dimensions
      wsize = pmh**2 + pmh + nobs
      if (lwa .lt. wsize) then
         info = 1
         return
      endif
      if (tau .lt. 0) then
	 info = 2
	 return
      endif
c			z1 is the first pmh columns of x
c			z2 is next h columns
      pmhp1 = pmh + 1
      nnull = npar - pmh
      hp1 = nnull + 1
      nmh = nobs - nnull
c			decompose z2 into f and g
      call dqrdc (x(1,pmhp1),ldx,nobs,nnull,fgaux,idummy,dummy,0)
c			create j as f' z1 
      do 30 j = 1,pmh
         call dqrsl (x(1,pmhp1),ldx,nobs,nnull,fgaux,x(1,j),dummy,
     *    x(1,j),dummy,dummy,dummy,01000,locinf)
   30 continue
c			decompose j2 (last npar - nnull rows of j)
c			as udv'
      if (pmh .gt. nmh) then
c			transpose the j2 matrix
         if (npar .gt. ldx) then
            info = 3
            return
         endif
         do 40 i = 0, nmh-1
            call dcopy(pmh,x(hp1+i,1),ldx,work(i*pmh),1)
   40    continue
         do 50 i = 0, nmh-1
            call dcopy (pmh, work(i*pmh),1,x(hp1,i+1),1)
   50    continue
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,pmh,nmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 4
	       return
            endif
         else
            call dsvdc (x(hp1,1),ldx,pmh,nmh,svals,work(pmhp1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 5
	       return
            endif
            npsing=nmh
	    normk=0.0d0
         endif
         do 60 i = 1, npsing
            call dswap(pmh,x(hp1,i),1,v(1,i),1)
   60    continue
      else
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,nmh,pmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 6
	       return
            endif
         else
            call dsvdc (x(hp1,1),ldx,nmh,pmh,svals,work(nmh+1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 7
	       return
            endif
            npsing=pmh
	    normk = 0.0d0
         endif
      endif
      return
      end
Funky!Stuff!
echo yasa: extracting fact.f
cat - << \Funky!Stuff! > fact.f
      integer function fact(i)
      integer i
c 
c Purpose: quick factorial function for the bspline routine
c	returns zero for negative i.
c
c On Entry:
c   i			a non-negative integer
c On Exit:
c   fact		i factorial
c
c $Header: fact.f,v 2.100.1.1 86/10/07 13:00:27 lindstrom Exp $
c
      integer j
      fact = 0
      if (i .ge. 0) fact = 1
      if (i .le. 1) return
      do 10 j = 2,i
	 fact = fact*j
   10	continue
      return
      end
Funky!Stuff!
echo yasa: extracting mkpoly.f
cat - << \Funky!Stuff! > mkpoly.f
      integer function mkpoly(m,dim)
      integer m,dim
c
c  Purpose: compute the binomial coefficient of m + dim - 1 choose dim.
c  	This is the dimension of the space of polynomials which are in
c  	the null space of the smoothing penalty. Uses Andy Jaworski's
c	binomial coefficient algorithm that only requires integer
c	arithmetic.
c
c  On Entry:
c   m			order of derivatives in the penalty
c   dim	 		dimension of the variables to be splined
c
c  On Exit:
c   mkploy		(m + dim - 1) choose dim
c
c $Header: mkpoly.f,v 2.100.1.1 86/10/07 13:00:35 lindstrom Exp $
c
      integer i,j,k,k1,kcoef,n
c 			compute binomial coefficient 
c			m + dim - 1 choose dim
      n = m + dim - 1
      k1 = dim
      if (k1 .gt. n .or. k1 .lt. 0) then
         mkpoly = 0
         return
      endif
      k = k1
      if ((n - k1) .lt. k) then
         k = n - k1
      endif
      kcoef = 1
      j = n - k
      do 10 i = 1, k
         j = j + 1
         kcoef = (kcoef * j) / i
   10 continue
      mkpoly = kcoef
      return
      end
Funky!Stuff!
echo yasa: extracting prmut.f
cat - << \Funky!Stuff! > prmut.f
      subroutine prmut (x,npar,jpvt,job)
      integer npar,x(npar),jpvt(npar),job
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
c $Header: prmut.f,v 2.100.1.1 86/10/07 13:00:47 lindstrom Exp $
c
      integer i,j,k,t
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
echo yasa: extracting mktpar.f
cat - << \Funky!Stuff! > mktpar.f
      subroutine mktpar(bign,nobs,testno,bandwd,tpar,t,f0,
     *   k0,info)
      double precision bandwd,tpar(4,*),t(nobs),f0(nobs),k0(nobs)
      integer bign,nobs,testno,info
c
c  Purpose: To construct the "true" parameter vectors for an
c	integral equation test case.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   testno		test number (determines the spread of the peaks)
c   bandwd		band width for the filter function
c   t(nobs)		the location of the observations
c
c  On Exit:
c   tpar(4,bign/2+1)	values of the "true" parameters
c   f0(nobs)		the "true" values of f at locations t
c   k0(nobs)		the "true" values of k at locations t
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mktpar.f,v 2.100.1.2 86/10/13 14:09:36 lindstrom Exp $ 
      double precision rt2pi,twopi,sn,cs,fo,f00,ko,k00,fdiff,
     *	iovern
      integer j,i
c     
      twopi = 8.0d0 * atan(1.0d0)
      rt2pi = sqrt(twopi)
      info = 1
      if ((testno .le. 0) .or. (testno .gt. 4)) return
      info = 2
      if ((bign .le. 1) .or. (nobs .le. 1)) return
      info = 3
      if (int(bign/2) .ne. dble(bign)/2.0d0) return
      info = 0
      fdiff = f00(0,testno,rt2pi) - f00(1,testno,rt2pi)
      call dset(4*(bign/2+1),0.0d0,tpar,1)
      do 20 j = 0, bign/2 
	  do 10 i = 1,bign
	      iovern = dble(i)/dble(bign)
	      sn = sin(twopi*j*iovern)
	      cs = cos(twopi*j*iovern)
 	      fo = f00(iovern,testno,rt2pi)+( iovern-0.5)*fdiff
	      tpar(1,j+1)= tpar(1,j+1)+cs*fo
	      tpar(2,j+1)= tpar(2,j+1)+sn*fo
	      ko = k00(iovern,bandwd,rt2pi)
	      tpar(3,j+1)= tpar(3,j+1)+cs*ko
	      tpar(4,j+1)= tpar(4,j+1)+sn*ko
   10     continue
   20 continue
      call dscal(4*(bign/2+1),1.0d0/dble(bign),tpar,1)
      do 30 j = 1,nobs
          f0(j) = f00(t(j),testno,rt2pi)+(t(j)-0.5)*fdiff
          k0(j) = k00(t(j),bandwd,rt2pi)
   30 continue
      return
      end
c
      double precision function f00(x,testno,rt2pi)
      double precision x,rt2pi
      integer testno
c
      double precision s1, s2, mu(4)
      data mu(1), mu(2), mu(3), mu(4) / 0.2, 0.15, 0.1, 0.05 /,
     *   s1, s2 / 0.015, 0.045 /
c
      f00 = (exp(-((x-0.3)/s1)**2/2.0d0)/s1+
     *    2.0d0*exp(-((x-0.3-mu(testno))/s2)**2/2.0d0)/s2)/(3.0d0*rt2pi)
      return
      end
c
      double precision function k00(x,bandwd,rt2pi)
      double precision x,bandwd,rt2pi
c
      k00 = (exp(-(x/bandwd)**2/2.0d0)+
     *    exp(-((1-x)/bandwd)**2/2.0d0))/(rt2pi*bandwd)
      return
      end
Funky!Stuff!
echo yasa: extracting mkxys.f
cat - << \Funky!Stuff! > mkxys.f
      subroutine mkxys(bign,nobs,p,m,t,tpar,x,ldx,truey,
     *  sigma,ldsig,ftrig,ktrig,info)
      double precision t(nobs),tpar(4,*),x(ldx,*),truey(nobs),
     *  sigma(ldsig,*),ftrig(nobs),ktrig(nobs)
      integer bign,nobs,p,m,ldx,ldsig,info
 
c  Purpose: to create the design matrix, the observation vector and the
c	sigma matrix for the integral equation example.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   p			the number of parameters
c   m			order of the derivatives in the penalty
c   t(nobs)		locations for observations
c   tpar		values of the "true" parameters
c   ldx			leading dimension of x
c   ldsig 		leading dimension of sigma
c
c  On Exit:
c   x(nobs,p)		design matrix for dsnsm
c   truey(nobs)		"true" values for y (no noise added)
c   sigma(ldsig,p)	sigma matrix
c   ftrig(nobs)         trigonometric interpolant to f0
c   ktrig(nobs)         trigonometric interpolant to k0
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mkxys.f,v 2.100.1.1 86/10/07 13:30:50 lindstrom Exp $ 
      double precision twopi,sn,cs,tmp
      integer i,j,kb,r,flag
      
      twopi = 8.0d0 * atan(1.0d0)
      info = 1
      if (bign .le. 1 .or. nobs .le. 1) return
      info = 2
      if (int(bign/2) .ne. dble(bign)/2.0d0) return
      info = 3
      flag = 0
      if (p .eq. bign) then
	  flag = 1
      else
          if (int(p/2) .eq. dble(p)/2.0d0) return
      endif
      info = 0
      call dset(nobs,tpar(3,1),x(1,1),1)
      if (flag .eq. 1) then
	 kb = bign/2 + 1
	 r = bign/2 - 1
      else
         kb = (p - 1)/2 + 1
	 r = (p - 1)/2 
      endif
      do 20 i = 1, nobs
	  cs = cos( twopi*bign*t(i)/2.0d0)
	  truey(i) = tpar(1,1)*tpar(3,1)+ tpar(1,bign/2+1)*
     *	    tpar(3,bign/2+1)*cs*2.0d0
	  ftrig(i)= tpar(1,1) + tpar(1,bign/2+1)*cs/2.0d0
	  ktrig(i)= tpar(3,1) + tpar(3,bign/2+1)*cs/2.0d0
	  do 10 j = 1,bign/2-1
	      sn = sin(twopi*j*t(i))
	      cs = cos(twopi*j*t(i))
	      truey(i) = truey(i) + 2*((tpar(1,j+1)*tpar(3,j+1)-
     *	       tpar(2,j+1)*tpar(4,j+1))*cs+
     *	       (tpar(1,j+1)*tpar(4,j+1)+tpar(2,j+1)*tpar(3,j+1))*sn)
	      ftrig(i) = ftrig(i) + 2.0d0*tpar(1,j+1)*cs+
     *		    2.0d0*tpar(2,j+1)*sn
	      ktrig(i) = ktrig(i) + 2.0d0*tpar(3,j+1)*cs+
     *		    2.0d0*tpar(4,j+1)*sn
	      if (j .le. r) then
	          x(i,1+j) = 2.0d0*(tpar(3,j+1)*cs+tpar(4,j+1)*sn)
	          x(i,kb+j) = 2.0d0*(tpar(3,j+1)*sn-tpar(4,j+1)*cs)
	      endif
   10     continue
      if (flag .eq. 1) then
	  x(i,kb) = tpar(3,bign/2+1)*cos(twopi*bign*t(i)/2.0d0)/2.0d0
      endif
   20 continue
c			calculate sigma
      do 30 i = 1,p
	  call dset(p,0.0d0,sigma(1,i),1)
   30 continue
      do 40 i = 1,r
	  tmp = (twopi*i)**(2*m)
	  sigma(1+i,1+i) = tmp
	  sigma(kb+i,kb+i) = tmp
   40 continue
      if (flag .eq. 1) then
	  sigma(kb,kb) = ((twopi*dble(bign/2.0d0))**(2*m))/2.0d0
      endif
      return
      end
Funky!Stuff!
echo yasa: extracting inteqn.f
cat - << \Funky!Stuff! > inteqn.f
c program inteqn
c
c Purpose: To exercise the truncated singular value capabilities of the
c   dsnsm driver with the deconvolution problem as described in 
c	Wahba, G. (1982) ``Constrained Regularization for Ill Posed
c	Linear Operator Equations, with Applications in Meteorology and
c	Medicine,'' in Statistical decision Theory and Related Topics
c	III, Vol. 2, eds. S.S.  Gupta, and J.O. Berger, Academic Press, 
c       383-418.
c   Sample input and output for one simulated example are provided in
c   the files inteqn.in and inteqn.out. Small variations in the output 
c   can be caused by different machines.
c
c   This driver can also be modified as indicated in the code to produce
c   different simulated data to test various aspects of the truncation.
c   In particular, changing the bandwidth of the kernel (i.e. the 
c   parameter "bandwd") changes the conditioning of the regularization 
c   problem.
c
c Variables which determine the example
c    nobs	number of observations
c    bign	dimension of the hilbert space (must be even)
c    npar	dimension of the subspace in which f is estimated
c		(npar <= bign, npar odd unless npar = bign)
c    bandwd 	the standard deviation of the kernel
c    testno	indicates which test function to use
c    ntbl	passed to dsnsm
c    job	job for dsnsm
c    tau	passed to dsnsm
c    sigst 	sigst is the standard deviation of the noise
c
c		The next 3 variables are not used unless the
c		user modifies this routine and supplies 
c		random number generators.
c
c    jobt	if job is nonzero then random t vector is used
c		else t is evenly spaced.
c    seedy	seed for generation of random error for y
c    seedt	seed for generation of random t vector
c
c  $Header: inteqn.f,v 2.100.1.2 86/10/12 16:36:51 lindstrom Exp $
c
      integer maxobs, maxtbl, maxpar, nnull, lwa, liwa
      parameter ( maxobs = 200 , maxtbl = 200, maxpar = 200, nnull=1,
     *  lwa = (maxpar-nnull)*(maxpar-2*nnull+2+maxobs)+maxpar+maxobs,
     *  liwa = 2*maxpar-nnull )
c			variables passed to dsnsm
      double precision x(maxobs,maxpar),y(maxobs),sigma(maxpar,maxpar),
     *  adiag(maxobs), lamlim(2), dout(5), coef(maxpar), svals(maxobs),
     *  tbl(maxtbl,3), auxtbl(3,3), work(lwa), tau
      integer nobs, ntbl, iout(3), iwork(liwa), info
c			variables which determine the example
      double precision bandwd, sigst, tpar(4,maxpar+1),t(maxobs),
     * ftrig(maxobs), f0(maxobs), ktrig(maxobs), k0(maxobs),temp,
     * ytrue(maxobs), eps(maxobs), fpred(maxobs), twopi
      integer testno, i, j, job, npar, bign, m, kb, r
c
c	The next 2 lines are commented out because the random numbers
c	are read from the input file
c     integer jobt, seedy, seedt
c     real rnor,uni
c			external routines 
      double precision dasum
c
      write(*,*) 'Enter nobs, bign, npar, bandwd,'
      write(*,*) '  testno, ntbl, job, tau and sigst'
      read(*,*) nobs
      read(*,*) bign
      read(*,*) npar
      read(*,*) bandwd
      read(*,*) testno
      read(*,*) ntbl
      read(*,*) job
      read(*,*) tau
      read(*,*) sigst
c
c	The next 4 lines are commented out because the random numbers
c	are read from the input file
c     write(*,*) 'job for t, seed for y, and seed for t'
c     read(*,*) jobt
c     read(*,*) seedy
c     read(*,*) seedt
      if (mod(job/10,10) .ne. 0) read(*,*) lamlim(1), lamlim(2)
      write(*,*) 'Nobs =', nobs, ' bign =', bign, ' npar =', npar, 
     *    ' bandwd =', bandwd, ' testno =', testno,'ntbl =',ntbl,
     *    'job =',job, ' tau =', tau,' sigst =', sigst
c
c	The next line is commented out because the random numbers are
c	read from the input file
c     write(*,*) ' jobt= ',jobt,' seedy =', seedy, ' seedt =', seedt
c
c	if a random number generator is used remove the next 3 lines
      do 10 i = 1,nobs
	  read(*,*) eps(i),t(i)
   10 continue
      if (nobs .gt. maxobs) then
	  write(*,*) 'nobs cannot exceed the maxobs of', maxobs
	  goto 999
      endif
      if (npar .gt. maxpar) then
	  write(*,*) 'npar cannot exceed the maxpar of', maxpar
	  goto 999
      endif
      if (ntbl .gt. maxtbl) then
	  write(*,*) 'ntbl cannot exceed the maxtbl of', maxtbl
	  goto 999
      endif
c
c	The next 11 lines are commented out because the random numbers 
c	are read from the input file
c     y(1) = rnor(seedy)
c     t(1) = uni(seedt)
c     if (jobt .eq. 0) then
c         do 20 i = 1, nobs
c             t(i) = dble(i) / dble(nobs)
c  20     continue
c     else 
c         do 30 i = 1,nobs
c             t(i) = uni(0)
c 30      continue
c     endif
      call mktpar(bign,nobs,testno,bandwd,tpar,t,f0,k0,info)
      if (info .ne. 0) then
	  write(*,*) 'info from mktpar is', info
	  goto 999
      endif
c     write(*,*) 'true parameters'
c     write(*,'(1p,7g11.3)') (tpar(1,i), i = 1, bign/2+1)
c     write(*,'(1p,7g11.3)') (tpar(2,i), i = 1, bign/2+1)
c     write(*,'(1p,7g11.3)') (tpar(3,i), i = 1, bign/2+1)
c     write(*,'(1p,7g11.3)') (tpar(4,i), i = 1, bign/2+1)
      m = 2
      call mkxys(bign,nobs,npar,m,t,tpar,x,maxobs,adiag,
     *  sigma,maxpar,ftrig,ktrig,info)
      if (info .ne. 0) then
	  write(*,*) 'info from mkxys is', info
	  goto 999
      endif
      call dcopy(nobs,adiag,1,ytrue,1)
      temp = 0.0d0
      do 40 i = 1,nobs
	   temp = temp + abs( f0(i)-ftrig(i) )
c
c	The next line is replaced by the one after it because 
c	the random numbers are read from the input file
c	   y(i) = ytrue(i)+sigst*rnor(0)
 	   y(i) = ytrue(i)+sigst*eps(i)
   40 continue
      write(*,*) 'sum abs(f0 - ftrig)  = ',temp
      temp = 0.0d0
      do 50 i = 1,nobs
	   temp = temp+ abs( k0(i)-ktrig(i) )
   50 continue
      write(*,*) 'sum abs(k0 - ktrig)  = ',temp
c     write(*,*) 'Location of responses'
c     write(*,'(1p,7g11.3)') (t(i), i = 1, nobs)
      write(*,*) 'Observed responses'
      write(*,'(1p,7g11.3)') (y(i), i = 1, nobs)
      call dsnsm(x, maxobs, y, sigma, maxpar, nobs, npar, nnull,
     *  adiag, tau, lamlim, ntbl, dout, iout, coef, svals, tbl, 
     *  maxtbl, auxtbl, iwork, liwa, work, lwa, job, info)
      if (info .gt. 0) then
	  write(*,*) 'Non-zero info from dsnsm of ', info
      endif
      write(*,*) ' info = ', info
      twopi = 8.0d0 * atan(1.0d0)
      if (npar .eq. bign) then
	 kb = bign/2 + 1
	 r = bign/2 - 1
      else
         kb = (npar - 1)/2 + 1
	 r = (npar-1)/2
      endif
      call dset(nobs,coef(1),fpred,1)
      do 70 i = 1, nobs
	  if (npar .eq. bign) then
 	      fpred(i)= fpred(i) + 
     *	      coef(bign/2+1)*cos(twopi*bign*t(i)/2.0d0)/2
          endif
	  do 60 j = 1,r
	      fpred(i) = fpred(i) + 
     *		    2.0d0*coef(1+j)*cos(twopi*j*t(i))+
     *		    2.0d0*coef(kb+j)*sin(twopi*j*t(i))
   60     continue
   70 continue
      write(*,*) 'ftrig'
      write(*,'(1p,7g11.3)') (ftrig(i), i = 1, nobs)
      write(*,*) 'fpred'
      write(*,'(1p,7g11.3)') (fpred(i), i = 1, nobs)
      temp = 0.0d0
      do 80 i = 1,nobs
	   temp = temp+(ytrue(i)-y(i))**2
   80 continue
      if (abs(temp/nobs-auxtbl(1,3)) .gt. 1.0d-8) then
	  write(*,*) 'error in R',temp/nobs, auxtbl(1,3)
          write(*,*) 'Predicted responses'
          write(*,'(1p,7g11.3)') (y(i), i = 1, nobs)
      endif
      write(*,*) 'Diagonal of hat matrix'
      write(*,'(1p,7g11.3)') (adiag(i), i = 1, nobs)
      write(*,*) 'lamlim =',lamlim(1),lamlim(2)
      write(*,*) 'Trace of I - A', nobs - dasum(nobs,adiag,1)
      write(*,*) 'dout:'
      write(*,*) '	lamhat           ',dout(1)
      write(*,*) '	penlty           ',dout(2)
      write(*,*) '	rss              ',dout(3)
      write(*,*) ' 	sqrt(rss/nobs)   ',sqrt(dout(3)/nobs)
      write(*,*) '	tr(I-A)          ',dout(4)
      write(*,*) '	truncation ratio ',dout(5)
      write(*,*) 'iout:'
      write(*,*) '	npsing           ',iout(1)
      write(*,*) '	npar             ',iout(2)
      write(*,*) '	nnull            ',iout(3)
      write(*,*) 'Coefficient estimates'
      write(*,'(1p,7g11.3)') (coef(i), i = 1, npar)
      write(*,*) 'Singular values'
      write(*,'(1p,7g11.3)') (svals(i), i = 1, iout(1))
c     write(*,*) 'Table'
c     do 90 i = 1, ntbl
c          write(*,'(1p,3g15.4)') (tbl(i,j), j = 1, 3)
c  90 continue
      write(*,*) 'Auxtbl'
      do 100 i = 1, 3
	  write(*,'(1p,3g15.4)') (auxtbl(i,j), j = 1, 3)
  100 continue
  999 continue
      stop
      end
Funky!Stuff!
echo yasa: extracting inteqn.in
cat - << \Funky!Stuff! > inteqn.in
100			#nobs
100			#bign
100			#p
.043			#bandwidth
4			#testno
150			#ntbl
1101			#job
50.00			#factor
0.05			#sigst
  0.369426  0.664586		#errors for y	tvec
 -0.161741   4.74875e-02
  0.110734  0.888196		#These were created using the CMLIB
 -0.892269  0.758726		#random number generators rnor and 
  0.563741   4.12114e-02	#uni with seeds 2345 and 234567
   1.45182  0.539204		#respectively. 
  -1.31448  0.333283
   1.88693  0.316610
 -0.376364  0.110962
  -1.19051  0.802525
  0.119080  0.710972
  0.355376  0.456799
   6.85814e-02  0.394802
 -0.260287  0.200948
  -1.11644   8.98874e-02
 -0.134926  0.720376
 -0.820781  0.323931
   1.64792  0.730215
  0.355515  0.153460
 -0.627039  0.201692
 -0.288621  0.961650
  0.809240  0.282720
 -0.337654  0.191011
 -0.380919  0.820178
   1.52130e-02  0.885082
  0.619401  0.850688
  -7.60663e-02  0.480195
   3.71923e-02  0.480040
   7.20295e-03  0.363379
  0.871532  0.490280
   1.47435  0.649740
 -0.226480  0.390308
  -1.04995  0.759664
 -0.194730   3.94473e-02
 -0.123651  0.760065
  0.258064  0.496280
 -0.399952  0.188616
 -0.458049  0.798014
  0.663119  0.756727
   1.40577  0.569053
 -0.214525  0.676102
  -1.17632  0.303534
  -2.29485  0.947326
  0.946033  0.276532
  -1.91029   8.90136e-02
  0.187062  0.312723
   1.11663  0.813254
  0.179122  0.297587
 -0.279380  0.886225
  -2.96400e-02  0.329350
  -5.08641e-02  0.273276
 -0.580411   5.31892e-02
  -2.51686  0.801307
  -1.16734  0.697609
  0.530081  0.531336
 -0.586954  0.516549
  -1.19345  0.484136
  -1.00249  0.125205
  0.454596  0.394075
 -0.788416  0.584009
  0.295470  0.240016
 -0.467823  0.395122
  -1.33313  0.812482
  -1.78970  0.580821
  0.352242  0.286423
 -0.237447  0.353792
 -0.960140   6.57725e-02
  -1.83060  0.539205
  0.391730  0.527632
 -0.490526  0.485116
  -1.12822  0.656183
 -0.170927  0.534437
  -3.92629e-02   2.26565e-02
  -1.43795   4.34956e-02
 -0.400290  0.359911
  0.238791  0.262108
   1.29887  0.950427
  -1.28274  0.782640
 -0.297425  0.648373
   1.13850  0.547429
  0.637391  0.681287
 -0.485349  0.664004
  0.546681  0.428849
   1.75646  0.582601
 -0.300003   8.22389e-03
  0.592774  0.153656
  -8.13912e-03  0.178889
  0.102742  0.772666
  0.974032   4.81642e-02
 -0.340513  0.985567
   1.62357  0.110160
   1.14206  0.818978
 -0.684942  0.510558
   2.40373   9.77370e-02
   1.34218  0.202927
  0.539613  0.461787
 -0.929959  0.271549
 -0.707846  0.829271
 -0.859077  0.433733
   1.90618  0.774079
Funky!Stuff!
echo yasa: extracting inteqn.out
cat - << \Funky!Stuff! > inteqn.out
  Enter nobs, bign, npar, bandwd,
    testno, ntbl, job, tau and sigst
  Nobs =  100   bign =  100   npar =  100   bandwd =   4.3000000000000d-02
   testno =  4  ntbl =  150  job =  1101   tau =   50.000000000000   sigst =
   5.0000000000000d-02
  sum abs(f0 - ftrig)  =    7.1095891316750d-04
  sum abs(k0 - ktrig)  =    1.4738158355217d-13
  Observed responses
  1.848e-02 -8.055e-03  5.537e-03 -4.461e-02  2.821e-02  0.115       6.29    
   6.53     -1.561e-02 -5.953e-02  5.954e-03   1.01       3.64      0.504    
 -5.506e-02 -6.746e-03   6.42      8.240e-02  6.347e-02  0.503     -1.443e-02
   5.14      0.313     -1.905e-02  7.607e-04  3.097e-02  0.477      0.485    
   5.28      0.381      7.376e-02   3.86     -5.250e-02 -9.719e-03 -6.183e-03
  0.283      0.275     -2.290e-02  3.316e-02  7.902e-02 -1.072e-02   6.09    
 -0.115       4.73     -9.480e-02   6.39      5.583e-02   5.92     -1.397e-02
   6.42       4.45     -2.897e-02 -0.126     -5.837e-02  8.783e-02  8.979e-02
  0.360     -4.199e-02   3.69     -3.578e-02   2.14       3.59     -6.666e-02
 -8.508e-02   5.35       5.71     -4.788e-02 -4.944e-02  9.239e-02  0.381    
 -5.639e-02  4.443e-02 -1.959e-03 -7.187e-02   5.43       3.65      6.494e-02
 -6.414e-02 -1.483e-02  8.485e-02  3.187e-02 -2.425e-02   2.00      9.179e-02
 -1.500e-02  7.586e-02  0.182      5.137e-03  4.874e-02 -1.703e-02  8.422e-02
  5.710e-02  0.119      0.121      0.630      0.884       4.29     -3.539e-02
   1.73      9.531e-02
   info =   0
  ftrig
  3.612e-06 -2.868e-06  2.778e-06 -5.110e-06 -7.894e-06  8.495e-04   6.27    
   9.29      1.303e-05 -3.576e-06  3.964e-06  0.354       3.60      2.449e-02
  4.941e-06 -4.698e-06   7.48      4.866e-06  4.427e-04  2.587e-02 -6.965e-06
   6.50      1.152e-02 -5.047e-06 -2.484e-06  5.234e-06  8.993e-02  9.083e-02
   5.66      4.586e-02  5.400e-06   3.96     -5.121e-06 -3.457e-06 -4.995e-06
  2.998e-02  9.519e-03 -4.110e-06 -3.714e-06  4.948e-05 -5.371e-06   12.1    
 -1.684e-07   4.17      1.600e-06   10.4      2.959e-06   11.8     -5.723e-07
   6.63       3.19      1.056e-05 -4.619e-06 -5.837e-06  1.760e-03  6.256e-03
  6.955e-02  8.983e-06   3.66      1.353e-05  0.301       3.58      3.835e-06
  8.534e-06   8.06       5.90     -8.497e-06  8.494e-04  2.455e-03  6.516e-02
 -6.108e-06  1.321e-03 -9.154e-06 -1.014e-05   5.77       1.24      5.623e-06
 -3.005e-06  6.769e-06  4.009e-04 -3.152e-06  2.565e-06   1.27      1.141e-05
  4.733e-07  4.510e-04  4.286e-03  2.748e-06 -7.292e-07 -5.205e-06  9.668e-06
 -4.584e-06  1.017e-02  4.885e-06  2.830e-02  0.270       2.76      4.689e-06
   1.05      6.834e-07
  fpred
  8.800e-03 -1.526e-02 -2.870e-02  0.187     -1.734e-02  0.214       8.28    
   9.31     -9.239e-02 -0.330     -7.312e-02  0.416       3.18     -0.286    
 -6.266e-02  6.610e-02   9.01      0.202      0.236     -0.309      6.149e-02
   6.68      1.938e-03 -0.136      1.502e-02  0.288     -0.183     -0.181    
   5.15     -0.231      0.149       3.37      0.174     -1.797e-02  0.169    
 -0.209      6.311e-02 -0.336      0.212      1.029e-02 -0.127       9.08    
 -2.159e-02   5.63     -6.058e-02   9.35     -0.239       8.63     -1.204e-03
   8.62       5.05     -1.434e-02 -0.334     -0.204      0.179      1.998e-02
 -0.216     -4.740e-02   3.21     -0.130      0.142       3.17     -0.249    
 -0.108       7.26       6.11     -1.921e-02  0.214      0.148     -0.221    
  0.101      0.198     -1.831e-02 -1.652e-02   5.48       3.06      1.216e-03
 -0.208      0.156      0.205     -0.175      1.579e-02   1.73     -0.121    
  6.124e-04  0.238      0.249     -4.148e-02 -1.508e-02  6.248e-02 -9.279e-02
 -0.156     -5.927e-02 -8.035e-02 -0.345      0.233       4.74      2.248e-02
   1.49     -6.666e-02
  Diagonal of hat matrix
  0.151      0.134      0.312      0.146      0.135      0.143      0.161    
  0.143      0.192      0.124      0.235      0.206      0.190      0.184    
  0.188      0.231      0.148      0.220      0.291      0.187      0.284    
  0.131      0.162      0.147      0.296      0.351      0.133      0.133    
  0.198      0.131      0.239      0.185      0.144      0.138      0.143    
  0.134      0.164      0.129      0.150      0.190      0.172      0.138    
  0.357      0.140      0.189      0.142      0.126      0.134      0.301    
  0.155      0.149      0.143      0.125      0.230      0.134      0.133    
  0.131      0.237      0.189      0.271      0.347      0.190      0.125    
  0.245      0.129      0.195      0.176      0.143      0.131      0.130    
  0.181      0.137      0.222      0.133      0.198      0.206      0.329    
  0.136      0.255      0.156      0.190      0.151      0.272      0.259    
  0.332      0.291      0.191      0.133      0.135      0.379      0.190    
  0.142      0.135      0.182      0.193      0.184      0.155      0.204    
  0.273      0.133    
  lamlim =  -18.566397372003   1.0942422965694
  Trace of I - A   81.211241602127
  dout:
  	lamhat              1.8215098946450d-10
  	penlty              1282664.0942576
  	rss                0.21618022670516
   	sqrt(rss/nobs)      4.6495185417972d-02
  	tr(I-A)             81.211241602127
  	truncation ratio   0.99999999739516
  iout:
  	npsing             35
  	npar               100
  	nnull              1
  Coefficient estimates
  0.999     -0.482     -0.444      0.699     -0.188     -0.299      0.214    
  0.164     -0.141     -2.782e-02  2.391e-02 -4.295e-03 -4.208e-03  3.679e-04
 -4.277e-04 -1.176e-04 -1.102e-05 -2.175e-06  1.464e-07 -2.450e-07  5.380e-08
  2.575e-09 -2.347e-10  2.022e-10  1.068e-11  1.590e-12  2.051e-13 -8.216e-15
  6.862e-15  6.606e-16 -8.605e-19 -8.152e-18 -1.268e-18 -7.801e-19 -1.722e-18
  1.481e-18 -1.371e-18 -1.154e-18  2.337e-18 -1.343e-18  7.317e-20 -1.569e-18
  1.350e-18  2.090e-18 -1.529e-18  7.475e-19 -7.889e-19 -2.184e-20 -1.412e-19
 -2.785e-19  4.409e-21  0.836     -0.730     -4.605e-02  0.501     -0.242    
 -0.172      0.186      8.072e-02 -0.183      1.354e-02  3.479e-02 -4.278e-03
  4.978e-04 -2.551e-06 -3.051e-05 -1.559e-05 -1.276e-06 -7.273e-07 -3.107e-07
 -1.421e-08  2.259e-09 -2.560e-10 -1.122e-10 -8.050e-12  7.009e-13 -6.547e-14
  2.524e-14 -4.014e-16  4.221e-16  2.853e-17 -6.096e-19  1.556e-18 -1.927e-18
 -2.888e-18 -2.847e-18 -1.352e-18 -6.088e-19 -1.859e-18 -9.337e-19 -4.946e-19
  4.291e-19  1.519e-18  1.088e-18  1.535e-18 -2.280e-19  7.213e-19  1.591e-19
  7.549e-20 -2.629e-19
  Singular values
  0.352      0.337      8.378e-02  6.893e-02  2.929e-02  2.635e-02  1.278e-02
  1.129e-02  5.832e-03  5.199e-03  2.539e-03  2.520e-03  1.216e-03  1.024e-03
  5.428e-04  4.402e-04  2.203e-04  1.784e-04  8.951e-05  6.713e-05  3.139e-05
  2.697e-05  1.124e-05  9.047e-06  3.773e-06  2.963e-06  1.125e-06  9.363e-07
  3.479e-07  2.590e-07  9.231e-08  7.202e-08  2.367e-08  1.729e-08  5.210e-09
  Auxtbl
     -7.740         3.2778e-03     2.0745e-04
         0.         4.6951e-03     3.4099e-04
         0.          4.370          4.320    
Funky!Stuff!
echo yasa: extracting testptpss.f
cat - << \Funky!Stuff! > testptpss.f
c program testptpss
c
c  Purpose: Test the gcvpack driver dptpss.f
c
c  $Header: testptpss.f,v 2.100.1.2 86/10/30 16:00:27 lindstrom Exp $
c
      integer maxobs, maxuni, maxpar, maxtbl, maxnul, lwa, liwa
      parameter ( maxobs = 200 , maxtbl = 200, maxnul = 10,
     *  maxuni = 150 ,maxpar = maxuni+maxnul,
     *  lwa = maxnul*(maxobs+maxuni+1)+maxpar*(maxobs+maxpar)+
     *        maxpar*(maxpar+2+maxobs)+maxpar+maxobs,
     *  liwa = 3*maxobs )
c
      integer nobs,i,j,info,ntbl,ncov1,ncov2,job,m,dim,
     * iwork(liwa),iout(4)
      double precision s(maxobs,1),s2(maxobs,1),lamlim(2),
     * des(maxobs,4),des2(maxobs,4),y(maxobs),adiag(maxobs),
     * ytrue(maxobs),tbl(maxtbl,3),coef(maxpar),auxtbl(3,3),
     * svals(maxobs),dout(4),pred(maxobs),work(lwa),pderr,r,trA,diag
c
      double precision dasum
c
      write(*,*) 'Enter nobs, dim, ncov1, ncov2, m, ntbl, job'
      read(*,*) nobs
      read(*,*) dim
      read(*,*) ncov1
      read(*,*) ncov2
      read(*,*) m
      read(*,*) ntbl
      read(*,*) job
      if (mod(job/10,10) .ne. 0) read(*,*) lamlim(1), lamlim(2)
      write(*,*) 'nobs =', nobs, 'dim =',dim,'ncov1 =',ncov1,
     *    'ncov2 =',ncov2,'m =', m, 'ntbl =',ntbl,'job =',job
      if (mod(job/10,10) .ne. 0) then
	  write(*,*) 'lamlim =',lamlim(1), lamlim(2)
      endif
      if (nobs .gt. maxobs) then
	  write(*,*) 'nobs cannot exceed the maxobs of', maxobs
	  goto 999
      endif
      if (ntbl .gt. maxtbl) then
	  write(*,*) 'ntbl cannot exceed the maxtbl of', maxtbl
	  goto 999
      endif

      do 10  i=1,nobs 
          read(*,*) (des(i,j),j=1,dim),(s(i,j),j=1,ncov1+ncov2)
          read(*,*) ytrue(i),y(i)
   10 continue
      do 20 i = 1,dim
          call dcopy(nobs,des(1,i),1,des2(1,i),1)
   20 continue
      do 30 i = 1,ncov1+ncov2
      call dcopy(nobs,s(1,i),1,s2(1,i),1)
   30 continue
      call dcopy(nobs,ytrue,1,adiag,1)

      call dptpss(des,maxobs,nobs,dim,m,s,maxobs,ncov1,ncov2,y,ntbl,
     * adiag,lamlim,dout,iout,coef,svals,tbl,maxtbl,auxtbl,
     * work,lwa,iwork,liwa,job,info)
      if (info .ne. 0) write(*,*) 'dptpss info',info
      call dpred(des2,maxobs,nobs,dim,m,des,maxobs,iout(4),s2,maxobs,
     *  ncov1,0,coef,iout(2),pred,work,lwa,iwork,info)
      if (info .ne. 0) write(*,*) 'dpred info',info

      write(*,*) 'lamlim = ',lamlim(1),lamlim(2)
      write(*,*) 'dout:'
      write(*,*) '	lamhat           ',dout(1)
      write(*,*) '	penlty           ',dout(2)
      write(*,*) '	rss              ',dout(3)
      write(*,*) ' 	sqrt(rss/nobs)   ',sqrt(dout(3)/nobs)
      write(*,*) '	tr(I-A)          ',dout(4)
      write(*,*) 'iout:'
      write(*,*) '	npsing           ',iout(1)
      write(*,*) '	npar             ',iout(2)
      write(*,*) '	nnull            ',iout(3)
      write(*,*) '	nuobs            ',iout(4)
      write(*,*) 'auxtbl'
      do 40 i = 1,3 
	   write(*,*) (auxtbl(i,j),j=1,3)
   40 continue

      write(*,*) 'Coefficient estimates',(coef(i),i=1,iout(2))
      write(*,*) 'Singular values'
      write(*,'(1p,7g11.3)') (svals(i), i = 1, iout(1))
      R=0.0d0
      do 50 i=1,nobs 
          R=R+dble((ytrue(i)-y(i))**2)
	  pderr = pred(i)-y(i)
   50 continue
      R=R/dble(nobs)
      diag=dasum(nobs,adiag,1)
      trA=dble(nobs)-dout(4)
      if (abs(trA-diag) .gt. 1.0d-8) write(*,*) 'trA',trA,'diag',diag
      if (abs(R - auxtbl(1,3)) .gt. 1.0d-8) then
	   write(*,*) 'R=',R,'auxtblR=',auxtbl(1,3)
      endif
      if (abs(pderr) .gt. 1.0d-8) write(*,*) 'pderr = ',pderr

      stop
  999 continue
      end
Funky!Stuff!
echo yasa: extracting testptpss.in
cat - << \Funky!Stuff! > testptpss.in
50		#nobs
2		#dim
1		#ncov1
0		#ncov2
2		#m
200		#ntbl
101		#job
   -0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1554483570e+02
   -0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1576312613e+02
   -0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1867397827e+02
   -0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1849722168e+02
    0.e+00   -0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1966086310e+02
    0.e+00   -0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1980231311e+02
    0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1859838649e+02
    0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1851904737e+02
    0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1586842815e+02
    0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1603913832e+02
   -0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1092383867e+02
   -0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1114066546e+02
   -0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1481392847e+02
   -0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1482830425e+02
    0.e+00   -0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1656449698e+02
    0.e+00   -0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1644307297e+02
    0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1490792284e+02
    0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1505653924e+02
    0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1091956264e+02
    0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1094227538e+02
   -0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9614920104e+01
   -0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9646480938e+01
   -0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1403133439e+02
   -0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1403122345e+02
    0.e+00    0.e+00    0.e+00
    0.1591550715e+02    0.1577400253e+02
    0.e+00    0.e+00    0.e+00
    0.1591550715e+02    0.1600412514e+02
    0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1399627680e+02
    0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1402826553e+02
    0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9557001644e+01
    0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9584670472e+01
   -0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1120625177e+02
   -0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1108651907e+02
   -0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1483723493e+02
   -0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1499369172e+02
    0.e+00    0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1655494349e+02
    0.e+00    0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1651294369e+02
    0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1498448603e+02
    0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1471816070e+02
    0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1114575565e+02
    0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1117168689e+02
   -0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1582595514e+02
   -0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1596022497e+02
   -0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1864014953e+02
   -0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1856095997e+02
    0.e+00    0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1954375504e+02
    0.e+00    0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1980902641e+02
    0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1856884576e+02
    0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1861010439e+02
    0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1586586951e+02
    0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1590136745e+02
Funky!Stuff!
echo yasa: extracting testptpss.out
cat - << \Funky!Stuff! > testptpss.out
  Enter nobs, dim, ncov1, ncov2, m, ntbl, job
  nobs =  50  dim =  2  ncov1 =  1  ncov2 =  0  m =  2  ntbl =  200  job =  101
  lamlim =   -4.1843433265259   1.6093312869689
  dout:
  	lamhat              3.2679884650307d-05
  	penlty              853.57331502786
  	rss                0.25791662762274
   	sqrt(rss/nobs)      7.1821532651809d-02
  	tr(I-A)             26.775087178935
  iout:
  	npsing             21
  	npar               29
  	nnull              4
  	nuobs              25
  auxtbl
  -2.7867494807235   1.7988193836964d-02   5.8882753901201d-03
  0.   1.9378305182360d-02   5.8804882675552d-03
  0.   4.6729241085750   4.2954700571057
  Coefficient estimates   17.624151123995   9.2540344779902d-02
   1.0302805621053d-02   6.5633168108272  -11.356469787138   16.508680645501
   2.9818553538791  -3.5837240165699   8.2441126928124  -25.483696004923
   2.5476735637523 -0.29386952458562   26.760819459264  -48.727114493771
  -27.793851311316   41.817471784635   33.017395104093   38.257112513123
  -26.476296420657  -28.609071016442   7.1373012868979   8.9368709155005
   2.1849847711197  -22.883007849984  -4.0453833228645   8.9634672527690
  -2.7571023482932   10.840610441942  -6.1887696887454
  Singular values
  0.638      0.517      0.312      0.312      0.223      0.223      0.196    
  0.176      0.162      0.134      0.131      0.131      0.129      0.115    
  0.115      0.104      0.102      9.972e-02  8.867e-02  8.867e-02  8.088e-02
Funky!Stuff!
echo yasa: extracting testtpss.f
cat - << \Funky!Stuff! > testtpss.f

c program testtpss
c
c  Purpose: Test the gcvpack driver dtpss.f
c
c  $Header: testtpss.f,v 2.100.1.2 86/10/30 15:59:54 lindstrom Exp $
c
      integer maxobs, maxuni, maxpar, maxtbl, mxncts, lwa, liwa
      parameter ( maxobs = 200 , maxtbl = 200, mxncts = 10,
     *  maxuni = 150 ,maxpar = maxuni+mxncts,
     *  lwa = maxuni*(2+mxncts+maxuni)+maxobs,
     *  liwa = 2*maxobs + maxuni)
c
      integer nobs,i,j,info,ntbl,ncov1,job,m,dim,iwork(liwa),iout(4)
      double precision s(maxobs,1),s2(maxobs,1),lamlim(2),
     * des(maxobs,4),des2(maxobs,4),y(maxobs),adiag(maxobs),
     * ytrue(maxobs),tbl(maxtbl,3),coef(maxpar),auxtbl(3,3),
     * svals(maxobs),dout(5),pred(maxobs),work(lwa),pderr,r,trA,diag
c
      double precision dasum
c
      write(*,*) 'Enter nobs, dim, ncov1, m, ntbl, job'
      read(*,*) nobs
      read(*,*) dim
      read(*,*) ncov1
      read(*,*) m
      read(*,*) ntbl
      read(*,*) job
      if (mod(job/100,10) .ne. 0) read(*,*) lamlim(1), lamlim(2)
      write(*,*) 'nobs =', nobs, 'dim =',dim,'ncov1 =',ncov1,
     *    'm =', m, 'ntbl =',ntbl,'job =',job
      if (mod(job/100,10) .ne. 0) then
	  write(*,*) 'lamlim =',lamlim(1), lamlim(2)
      endif
      if (nobs .gt. maxobs) then
	  write(*,*) 'nobs cannot exceed the maxobs of', maxobs
	  goto 999
      endif
      if (ntbl .gt. maxtbl) then
	  write(*,*) 'ntbl cannot exceed the maxtbl of', maxtbl
	  goto 999
      endif

      do 10  i=1,nobs 
          read(*,*) (des(i,j),j=1,dim),(s(i,j),j=1,ncov1)
          read(*,*) ytrue(i),y(i)
   10 continue
      do 20 i = 1,dim
          call dcopy(nobs,des(1,i),1,des2(1,i),1)
   20 continue
      do 30 i = 1,ncov1
      call dcopy(nobs,s(1,i),1,s2(1,i),1)
   30 continue
      call dcopy(nobs,ytrue,1,adiag,1)

      call dtpss(des,maxobs,nobs,dim,m,s,maxobs,ncov1,y,ntbl,adiag,
     * lamlim,dout,iout,coef,svals,tbl,maxtbl,auxtbl,work,lwa,iwork,
     * liwa,job,info)
      if (info .ne. 0) write(*,*) 'dtpss info',info
      call dpred(des2,maxobs,nobs,dim,m,des,maxobs,iout(4),s2,maxobs,
     *  ncov1,0,coef,iout(2),pred,work,lwa,iwork,info)
      if (info .ne. 0) write(*,*) 'dpred info',info

      write(*,*) 'lamlim = ',lamlim(1),lamlim(2)
      write(*,*) 'dout:'
      write(*,*) '	lamhat           ',dout(1)
      write(*,*) '	penlty           ',dout(2)
      write(*,*) '	rss              ',dout(3)
      write(*,*) ' 	sqrt(rss/nobs)   ',sqrt(dout(3)/nobs)
      write(*,*) '	tr(I-A)          ',dout(4)
      write(*,*) '      ssqrep           ',dout(5)
      write(*,*) 'iout:'
      write(*,*) '	npsing           ',iout(1)
      write(*,*) '	npar             ',iout(2)
      write(*,*) '	nnull            ',iout(3)
      write(*,*) '	nuobs            ',iout(4)
      write(*,*) 'auxtbl'
      do 40 i = 1,3 
	   write(*,*) (auxtbl(i,j),j=1,3)
   40 continue

      write(*,*) 'Coefficient estimates',(coef(i),i=1,iout(2))
      write(*,*) 'Singular values'
      write(*,'(1p,7g11.3)') (svals(i), i = 1, iout(1))
      R=0.0d0
      do 50 i=1,nobs 
          R=R+dble((ytrue(i)-y(i))**2)
	  pderr = pred(i)-y(i)
   50 continue
      R=R/dble(nobs)
      diag= dasum(nobs,adiag,1)
      trA=dble(nobs)-dout(4)
      if (abs(trA-diag) .gt. 1.0d-8) write(*,*) 'trA',trA,'diag',diag
      if (abs(R - auxtbl(1,3)) .gt. 1.0d-8) then
	   write(*,*) 'R=',R,'auxtblR=',auxtbl(1,3)
      endif
      if (abs(pderr) .gt. 1.0d-8) write(*,*) 'pderr = ',pderr

      stop
  999 continue
      end
Funky!Stuff!
echo yasa: extracting testtpss.in
cat - << \Funky!Stuff! > testtpss.in
50		#nobs
2		#dim
1		#ncov1
2		#m
200		#ntbl
1011		#job
   -0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1554483570e+02
   -0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1576312613e+02
   -0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1867397827e+02
   -0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1849722168e+02
    0.e+00   -0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1966086310e+02
    0.e+00   -0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1980231311e+02
    0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1859838649e+02
    0.5000000000e+00   -0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1851904737e+02
    0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1586842815e+02
    0.1000000000e+01   -0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1603913832e+02
   -0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1092383867e+02
   -0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1114066546e+02
   -0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1481392847e+02
   -0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1482830425e+02
    0.e+00   -0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1656449698e+02
    0.e+00   -0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1644307297e+02
    0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1490792284e+02
    0.5000000000e+00   -0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1505653924e+02
    0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1091956264e+02
    0.1000000000e+01   -0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1094227538e+02
   -0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9614920104e+01
   -0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9646480938e+01
   -0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1403133439e+02
   -0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1403122345e+02
    0.e+00    0.e+00    0.e+00
    0.1591550715e+02    0.1577400253e+02
    0.e+00    0.e+00    0.e+00
    0.1591550715e+02    0.1600412514e+02
    0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1399627680e+02
    0.5000000000e+00    0.e+00    0.e+00
    0.1404538577e+02    0.1402826553e+02
    0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9557001644e+01
    0.1000000000e+01    0.e+00    0.e+00
    0.9653243053e+01    0.9584670472e+01
   -0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1120625177e+02
   -0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1108651907e+02
   -0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1483723493e+02
   -0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1499369172e+02
    0.e+00    0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1655494349e+02
    0.e+00    0.5000000000e+00    0.2500000000e+00
    0.1654538577e+02    0.1651294369e+02
    0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1498448603e+02
    0.5000000000e+00    0.5000000000e+00    0.2500000000e+00
    0.1489500943e+02    0.1471816070e+02
    0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1114575565e+02
    0.1000000000e+01    0.5000000000e+00    0.2500000000e+00
    0.1101895709e+02    0.1117168689e+02
   -0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1582595514e+02
   -0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1596022497e+02
   -0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1864014953e+02
   -0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1856095997e+02
    0.e+00    0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1954375504e+02
    0.e+00    0.1000000000e+01    0.1000000000e+01
    0.1965324305e+02    0.1980902641e+02
    0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1856884576e+02
    0.5000000000e+00    0.1000000000e+01    0.1000000000e+01
    0.1851895709e+02    0.1861010439e+02
    0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1586586951e+02
    0.1000000000e+01    0.1000000000e+01    0.1000000000e+01
    0.1585498788e+02    0.1590136745e+02
Funky!Stuff!
echo yasa: extracting testtpss.out
cat - << \Funky!Stuff! > testtpss.out
  Enter nobs, dim, ncov1, m, ntbl, job
  nobs =  50  dim =  2  ncov1 =  1  m =  2  ntbl =  200  job =  1011
  lamlim =   -4.1843433265259   1.6093312869689
  dout:
  	lamhat              3.2679884666661d-05
  	penlty              853.57331501868
  	rss                0.25791662763774
   	sqrt(rss/nobs)      7.1821532653897d-02
  	tr(I-A)             26.775087179714
        ssqrep             0.24222881477953
  iout:
  	npsing             21
  	npar               29
  	nnull              4
  	nuobs              25
  auxtbl
  -2.7867494805062   1.7988193836964d-02   5.8882753902953d-03
  0.   1.9378305182360d-02   5.8804882675556d-03
  0.   4.6729241085750   4.2954700571057
  Coefficient estimates   17.624151123986   9.2540344777235d-02
   1.0302805622788d-02   6.5633168108123  -11.356469786953  -25.483696004900
  -27.793851310920  -28.609071016127  -4.0453833230143   16.508680644902
   2.5476735641137   41.817471783039   7.1373012871226   8.9634672524902
   2.9818553539789 -0.29386952348460   33.017395104218   8.9368709157430
  -2.7571023477229  -3.5837240157703   26.760819457853   38.257112512072
   2.1849847715526   10.840610441479   8.2441126918716  -48.727114492199
  -26.476296420734  -22.883007849908  -6.1887696887032
  Singular values
  0.638      0.517      0.312      0.312      0.223      0.223      0.196    
  0.176      0.162      0.134      0.131      0.131      0.129      0.115    
  0.115      0.104      0.102      9.972e-02  8.867e-02  8.867e-02  8.088e-02
Funky!Stuff!

