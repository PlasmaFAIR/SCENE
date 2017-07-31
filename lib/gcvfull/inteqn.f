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
