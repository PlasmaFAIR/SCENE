cprogram testptpss
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
