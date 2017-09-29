program scene
  !     *************
  !
  use param 
  implicit none
  integer i,j,icur,ierr,niter
  double precision errcur,bpold,bnold,bpnew
  !

  nbet=0
  ! Perform initialisations
  write(6,*) 'Starting SCENE'
  call init

  ! Initialise error function array for use in bootstrap routine, hirsig
  ! In call bootstrap routine hirsig, don;'t forget you need to load up
  !  error function first...and don't call it errf!
  !  Junk that needs to be taken account of before bootstrap, etc
  !  Calculate error function array
  call erfun

!!$C  Initialise ICUR parameter...if ICUR=-1 trapped particle
!!$C  effects are not included in the plasma conductivity.
!!$      ICUR=0
  !      if (igr.gt.0) then
  !	call paper(1)
  !	call paper(1)
  !       call devoff
  !	call filnam('graphs.grd')
  !        call filon
  !end if

  !  initialise number of iterations
5 niter=0
  ! Call SCENE equilibrium routines
30 niter=niter+1



  if (ipr.eq.0) write(nw,*)' on the ',niter,' iteration'
  write(6,*) 'Printing poloidal/toroidal graphs'
  call equil(niter)

  write(6,*)' done equil'
  ! Convert equilibrium to flux surface coordinates
  call flxsur
  write(6,*)' done flxsur'
  ! Perform a number of useful averages over flux surfaces
  call flxav
  ! Calculate alpha-particle bootstrap current
  call fastbs
  write(6,*)' done flxav'

  j_nb = 0.
  if (nbi.eq.1) call nbicur()
  
  print*, 'done nbicur'
  !  Calculate currents (and read in externally applied current
  !  profile if itot=0)
  call torcur(icur)

  !do i=1,nr
  !   print*, i, nbph(i,nsym)
  !end do

  !  Calculate total current contributions....
  call xarea(bsph,totbs)
  call xarea(absph,alfbs)
  call xarea(psph,totps)
  call xarea(diph,totdi)
  call xarea(exph,totex)
  call xarea(nbph,totnb)
  print*, 'Total NB ',totnb
  call xarea(exph2,totex2)
  call xarea(gradj,totgs)

  if (itot.eq.0) then
     !  First calculate scaling factor on external current to give
     !  required total current (this will be the loop voltage in the
     !  case of a neoclassical ohmic profile)
     !vloop=(totgs-totbs-totps-totdi-totex2-totnb)/(totex)
     vloop=(totgs-totbs-totps-totdi-totex2)/(totex)

     print*, ' '
     print*, 'Vloop = ', vloop
     print*, ' '
     if (abs(neo).eq.1) then

        !vnobs=(totgs-totps-totdi-totex2-totnb)/(totex)
        vnobs=(totgs-totps-totdi-totex2)/(totex)
        call xarea(spit,spiti)
        vspit=(totgs-totps-totdi-totex2)/spiti
     end if
     !  Label exph as total driven current...
     do i=1,nr
        do j=1,nz
           exph(i,j)=exph(i,j)*vloop+exph2(i,j)
        end do
     end do
     totex=totex*vloop+totex2
     call getdata
     print*, totgs, totbs, totps, totdi, totex2, totnb, totex
     
     !  calculate error in ff'
     call ffdgen(1,icur,errcur)
     write(6,*)' done ffdgen'
     if (errcur.gt.errffd) then
        !  need to iterate further on ff'
        if (niter.lt.npass) then
           !  up-date ff'
           call ffdgen(0,icur,errcur)
           !  up-date n,T mesh if appropriate
           write(6,*)' in setnt'
           if (igr.eq.5) call setnt
           write(6,*)' out setnt'
           icont=-3
           goto 30
        else
           write(nw,*)'warning***ffd has not converged to specified'
           write(nw,*)'accuracy. accuracy=',errcur,' requested=',errffd
        end if
     end if
  end if
  !  equilibrium converged to specified accuracy...
  !  generate output
  call output
  !  If we want to converge in betan, start iteration again
  if (betan.gt.0.) then
     nbet=nbet+1
     write(6,*)' betan=',3.5*betexp/betlim
     if (abs(betan-3.5*betexp/betlim)/betan.gt.0.01) then
        if (nbet.eq.1) then
           bpold=bpol
           bnold=betlim/(3.5*betexp)
           bpol=bpol*betan*betlim/(3.5*betexp)
           write(6,*)' nbet=',nbet,' betan=',3.5*betexp/betlim, &
                ' target=',betan
        else
           bpnew=(betan*(bpold-bpol)+bnold*bpol-bpold*3.5*betexp/betlim)  &
                /(bnold-3.5*betexp/betlim)
           bpold=bpol
           bnold=3.5*betexp/betlim
           bpol=bpnew
           write(6,*)' nbet=',nbet,' betan=',3.5*betexp/betlim, &
                ' target=',betan
        end if
        !          icont=1
        goto 5
     end if
  end if
  if (ipswtch.eq.3) call setnt
  write(6,*)' done setnt'
  !  Plot out some useful figs
  !call nbicur2()
  if (igr.eq.0) call getdata
  !  Plot stability plots if igr set to 3
  !      if (igr.eq.3) call stab
  !      if (igr.eq.3) call epsplot
  !  Call user interface routine
  call usrcal
  !
  !
  !
!!$c  gst stores values of ff'
!!$      do 20 k=1,ncon+1
!!$	gst(k)=0.
!!$ 20   continue
  !
  !
  !C Call user ouput interface
  !      CALL USRCAL(U,IXOUT)
  !      if (igr.gt.0) call grend
end program scene
