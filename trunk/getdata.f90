subroutine getdata

  ! Puts data in the needed in the output file


  use param
  implicit none

  double precision :: rr, zz, rat
  double precision :: fprof, tempe, tempi, psi, bp, B_tor, B_tot
  double precision :: press, densi, dense, safety, ffp, fsi, pp, bsmean
  double precision :: J_tot, J_ext, J_bs, J_di, J_ps, J_nbi, J_ext2, jnb, jbs
  double precision :: flux_r, flux_z, te

  double precision, dimension(ncon) :: rhoflux, Rflux
  double precision :: rhomax, jtot

  integer :: i, j, nh, con
  nh = 49


  !Open file
  open(unit=nh, file=runname(1:lrunname)//'_overall.dat', &
       status='unknown', iostat=ios)
  if (ios.ne.0) then
     write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_overall.dat'
     stop
  endif

  !Writes following data along Z=0
  write(nh,*) 'R(m), ne(m^-3), ni(m^_3), B_T (T) , B_P (T),  B(T), Sfac, T(keV), Pressure'

  zz = 0.

  do i=1,nr

     !If point is inside boundary
     if (ixout(i,nsym) .eq. 1) then

        psi = umax-u(i,nsym)
        rr = r(i)

        B_tor = sngl(fprof(psi,2)/rr)
        B_tot = sngl(sqrt(B_tor**2 + (bp(rr,zz)**2)) )

        con = 1

        !Find closest flux surface and interpolate safety factor
        do j=1,ncon-1
           if (psiv(j) .lt. psi) exit

           con = j
        end do

        !if (con .eq. ncon-1) con = con-1

        rat = (psi-psiv(con))/(psiv(con+1) - psiv(con) )


        !Linearly interpolate q from flux surface grid to mesh grid
        safety = sngl( sfac(con) + rat*(sfac(con+1) - sfac(con)))


        write(nh,14) rr,dense(psi,0),densi(psi,1,0), B_tor,bp(rr,zz),B_tot, safety, tempe(psi,0)/1000., &
             press(psi,0)

     !If point is on plasma boundary
     else if (ixout(i,nsym) .eq. -1) then

        psi = umax


        safety = sngl(sfac(1))

        !Interpolate onto actual flux surface of bdry
        rat = (u(i,nsym))/(u(i,nsym) - u(i-1,nsym) )

        rr = r(i) - rat*(r(i) - r(i-1))


        B_tor = sngl(fprof(psi,2)/rr)
        B_tot = sngl(sqrt(B_tor**2 + (bp(rr,zz)**2)) )

        !Linearly interpolate q from flux surface grid to mesh grid


        write(nh,14) rr,dense(psi,0),densi(psi,1,0), B_tor,bp(rr,zz),B_tot, safety, tempe(psi,0)/1000., &
             press(psi,0)

     end if



  end do

14 format(ES14.6e2,ES14.6e2, ES14.6e2, ES14.6e2,ES14.6e2, ES14.6e2, ES14.6e2,ES14.6e2, ES14.6e2 )
   close(nh)

   !write(6,*) 'Overall data written'




   !Write Flux surface points
   nh = nh+2

   open(unit=nh, file=runname(1:lrunname)//'_fluxsur.dat', &
        status='unknown', iostat=ios)
   if (ios.ne.0) then
      write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_fluxsur.dat'
      stop
   end if

   write(nh,*) npts, ncon
   do i=1,ncon
      do j= 1,npts
         flux_r = rpts(i,j)
         flux_z = zpts(i,j)


         write(nh,15) flux_r, flux_z

      end do

      write(nh,*)

   end do
15 format(ES14.5e3, ES14.5e3)

   close(nh)

   !write(6,*) 'Flux surface data written'




   ! Write Current profiles
   nh = nh+2

   open(unit=nh, file=runname(1:lrunname)//'_currents.dat', &
        status='unknown', iostat=ios)
   if (ios.ne.0) then
      write(6,*) 'problem create/opening ',runname(1:lrunname)//'_currents.dat'
      stop
   end if

   write(nh,*) 'R(m), J_tot (kA m^-2), J_ext, J_ext2, J_bs, J_di, J_ps, J_nbi'

   do i=1,nr
      if (ixout(i,nsym) .eq. 0) cycle

      !Different current values
      J_tot = sngl(gradj(i,nsym)/1000.)
      J_ext = sngl(exph(i,nsym)/1000.)
      J_ext2 = sngl(exph2(i,nsym)/1000.)
      J_bs = sngl(bsph(i,nsym)/1000.)
      J_di = sngl(diph(i,nsym)/1000.)
      J_ps = sngl(psph(i,nsym)/1000.)
      J_nbi = sngl(nbph(i,nsym)/1000.)
      rr = r(i)


      write(nh,16) rr,J_tot,J_ext,J_bs,J_di,J_ps,J_nbi, J_ext2

   end do
16 format(8e14.6)


   close(nh)

   nh = nh+2



   nh = nh + 2


   !J profile with psi values
   open(unit=nh, file=runname(1:lrunname)//'_currentprof.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening ',runname(1:lrunname)//'_currentprof.dat'
      stop
   end if

   !Used to find magnetic axis


   rhoflux =( maxval(rpts, dim=2) - minval(rpts,dim=2)) /2

   Rflux = ( maxval(rpts, dim=2) + minval(rpts,dim=2)) /2

   rhoflux(ncon) = 0.
   Rflux(ncon) = r0
   rhomax = maxval(rhoflux)
   rhoflux = rhoflux/rhomax

   do con=ncon,1,-1

      psi = psiv(con)
      rr = sngl(Rflux(con))

      fsi = sngl(fprof(psi,2))
      ffp = -sngl(fprof(psi,1))
      pp = -sngl(press(psi,1))
      bsmean = sngl(bsqav(con))
      
      jtot = +fsi*pp/bsmean + ffp/(fsi*mu0)
      jbs = bsj(con)/sqrt(bsmean)
      jnb = J_nb(con)
      write(nh,18) psi, jtot, jnb, jbs


   end do

!    psi = umax - u(nr,nsym)
!    ffp = sngl(fprof(psi,1))
!    pp = sngl(press(psi,1))
!   write(nh,18) psi, ffp, pp

18 format(6e14.6)

   close(nh)

  ! write(6,*) 'current prof data written'


   nh=nh+2
   open(unit=nh, file='profile_ne.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening profile_ne.dat'
      stop
   end if

   write(nh, *) ncon

   do i=ncon,1,-1
      write(nh, 20) psiv(i), dense(psiv(i), 0)
   end do

20 format(2e14.6)
   close(nh)

   ! write(6,*) 'locust ne data written'



  nh=nh+2
   open(unit=nh, file='profile_ni.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening profile_ni.dat'
      stop
   end if

   write(nh, *) ncon

   do i=ncon,1,-1
      write(nh, 22) psiv(i), (densi(psiv(i),j,0), j=1,nimp+1)
   end do
22 format(7e13.6)
   close(nh)

   ! write(6,*) 'locust ni data written'




      nh=nh+2
   open(unit=nh, file='profile_Te.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening profile_Te.dat'
      stop
   end if

   write(nh, *) ncon

   do i=ncon,1,-1
      write(nh, 20) psiv(i), tempe(psiv(i), 0)
   end do


   close(nh)

 !  write(6,*) 'locust Te data written'


      nh=nh+2
   open(unit=nh, file='profile_Ti.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening profile_Ti.dat'
      stop
   end if

   write(nh, *) ncon

   do i=ncon,1,-1
      write(nh, 20) psiv(i), tempi(psiv(i),1, 0)
   end do


   close(nh)

 !  write(6,*) 'locust Ti data written'






      nh=nh+2
   open(unit=nh, file='psi.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening psi.dat'
      stop
   end if

   write(nh, *) nr,nz

   write(nh,2020) ((umax-u(i,j), i=1,nr), j=1,nz)

   close(nh)
2020 format(5e16.9)

 !  write(6,*) 'psi data written'



   call popcon()

   call jetto()
   call rfdat()
   call omfit()

   call vessel()
end subroutine getdata


subroutine popcon()

  use param
  implicit none


  double precision :: fprof
  double precision :: densi, dense

  integer :: nh


  nh = 63

  open(unit=nh, file='popcon.dat', &
       status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'_popcon.dat'
     stop
  end if


  write(nh,*) 'Parameters'
  write(nh,510) 'bpol', bpol, ' '
  write(nh,510) 'rcen', rcen, 'm'
  write(nh,510) 'a minor', amin, 'm'
  write(nh,510) 'eps ', tokeps, ' '
  write(nh,510) 'elon', elon ,' '
  write(nh,510) 'tri ', tri, ' '


  write(nh,510) 'Cen e temp', sngl(te0/1000.), 'keV'
  write(nh,510) 'Cen i temp', sngl(ti0/1000.), 'keV'
  write(nh,510) 'Vol av. e temp', avt/1000.,'keV'
  write(nh,510) 'Vol av. i temp', avti/1000.,'keV'


  write(nh,500) 'Cen e den', sngl(dense(0.0,0)), 'm^-3'
  write(nh,500) 'Cen i den', sngl(densi(0.0,1,0)), 'm^-3'

  write(nh,500) 'Vol av e den', sngl(avel), 'm^-3'
  write(nh,500) 'Lin av e den', sngl(nebar*1.0d19), 'm^-3'

  write(nh,500) 'Vol av i den', sngl(avio*1.0d19), 'm^-3'

  write(nh,*) ' '
  write(nh,510) 'Z', sngl(zm), ' '
  write(nh,510) 'Zeff', sngl(zeffav), ' '

  write(nh,*) ' '
  write(nh,510) 'B tor (mag)', sngl(fprof(0.0,2)/r0), 'T'
  write(nh,510) 'B tor (geo)', sngl(fprof(0.0,2)/rcen), 'T'
  write(nh,510) 'Vac B tor (geo)', sngl(mu0*rodi/(2.*pi*rcen)), 'T'

  write(nh,*) ' '
  write(nh,510) 'Tor. tot cur', sngl(cur/1000000.), 'MA'
  write(nh,510) 'Tor. bs cur', sngl(totbs/1000000.), 'MA'
  write(nh,510) 'Tor. nb cur', sngl(totnb/1000000.), 'MA'
  write(nh,510) 'Tor. ps cur', sngl(totps/1000000.), 'MA'
  write(nh,510) 'Tor. dia cur', sngl(totdi/1000000.), 'MA'
  write(nh,510) 'Tor. ext cur', sngl((totex+totex2)/1000000.), 'MA'

  write(nh, *) ' '

  write(nh,500) 'Rod Current', sngl(rodi/1000000.), 'MA'


  write(nh,510) 'Beta', sngl(beta), '%'
  write(nh,510) 'Norm. Beta', sngl(3.5*betexp/betlim), ' '
  write(nh,510) 'Pol. Beta', sngl(betap), ' '

  write(nh,510) 'q0' , sngl(sfac(ncon)), ' '
  write(nh,510) 'qa', sngl(sfac(1)), '  '
  write(nh,510) 'qmin',sngl(minval(sfac,1)), ' '

  write(nh,510) 'Pfus',pfus*1.0d-6*5, 'MW'
  write(nh,510) 'Beam fus',bmfus, 'MW'
  write(nh,510) 'Aux Pow',paux, 'MW'
  write(nh,510) 'NBP eff.', bmfus/(paux), ' '
  write(nh,510) 'NBCD eff.', totnb/(1.0d6*paux), ' '
  write(nh,510) 'Shine through', bmshine, 'MW'
  write(nh,510) 'Q',sngl( ((pfus*1.0d-6*5)+bmfus)/paux), ' '

  write(nh,510) 'N_gw', sngl(negw), ' '

  write(nh,510) 'H_IPB98(y1)', sngl(hipb98y1), ' '
  write(nh,510) 'H_IPB98(y2)', sngl(hipb98y2), ' '
  write(nh,510) 'H_Petty08',   sngl(petty), ' '
  write(nh,510) 'Tau_e', sngl(taue*1.0e3), 'ms'
  write(nh,510) 'Tau_h', sngl(tauh), ' '

  write(nh,510) 'Area', sngl(area), 'm^2'
  write(nh,510) 'Volume', sngl(vol), 'm^3'
  write(nh,510) 'li(3)', sngl(rli3), ' '
  write(nh,510) 'li(2)', sngl(rli2), ' '
  write(nh,510) 'Energy', sngl(conft/1.0e6), 'MJ'

500 format(' ',A12, ' , ', E11.4 , ' , ', A5)
510 format(' ',A12, ' , ', F11.4 ,' , ',A5)

  close(nh)

end subroutine popcon


subroutine jetto()
  
  use param
  implicit none
  
  integer :: con, nh, i, l, nflux
  double precision :: dense, densi, tempe, tempi, psi, mass,rmax,cs
  real :: ne, ni, te, ti, zeff, q, zni, vtor, angf, mach
  double precision, dimension(ncon) :: phi

  nh = 43

  open(unit=nh, file=runname(1:lrunname)//'.jetto', &
       status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.jetto'
     stop
  end if

  call rhotor(phi)
 
  write(nh,*) 'rho_tor, e den, i den, e temp (keV), i temp (keV), q, Zeff'
  do con = ncon, 1,-1

  
     psi = psiv(con)
     ne = sngl(dense(psi,0))
     ni = sngl(densi(psi,1,0))
     te = sngl(tempe(psi,0))
     ti = sngl(tempi(psi,1,0))
     q = sngl(sfac(con))
     mass = ni*mp*zmai
     rmax = maxval(rpts(con,:))
     if (con.eq.ncon) rmax=r0
     
     vtor = sngl(nbmom(con)*taue/mass)
     angf = vtor/sngl(rmax)

     cs = sqrt(eq*ti/(zmai*mp))

     mach = vtor/sngl(cs)

     print*, psi, vtor, angf, mach
  
     !Effective Z
     zeff=sngl(zm)
     if (imp.eq.1) then
        if (ne.gt.0.) then
           zeff=0.
           do l=1,nimp+1
              zni=sngl(densi(psi,l,0))
              zeff=zeff+(zni*iz(l)**2)/ne
           end do

        end if
     end if   
     !print*, phi(con), psi
     write(nh,500) sngl(phi(con)), ne, ni, te, ti,  q, zeff, angf

  end do
500 format(8e13.6)
  close(nh)


end subroutine jetto

subroutine omfit()
  
  use param
  implicit none
  
  integer :: con, nh, i, l, nflux
  double precision :: dense, densi, tempe, tempi, psi
  real :: ne, ni, te, ti, zeff, nHe, nIm, zni
  double precision, dimension(ncon) :: phi

  nh = 43

  open(unit=nh, file=runname(1:lrunname)//'.omfit', &
       status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.omfit'
     stop
  end if

  call rhotor(phi)
 
  write(nh,*) 'psi, e den (m^-3), i den (m^-3), e temp (keV), i temp (keV), He den (m^-3), Imp den (m^-3)'
  do con = ncon, 1,-1

  
     psi = psiv(con)
     ne = sngl(dense(psi,0))
     ni = sngl(densi(psi,1,0))
     te = sngl(tempe(psi,0))
     ti = sngl(tempi(psi,1,0))
     nHe = sngl(densi(psi,2,0))

     !Effective Z
     zeff=sngl(zm)
     if (imp.eq.1) then
        if (ne.gt.0.) then
           zeff=0.
           do l=1,nimp+1
              zni=sngl(densi(psi,l,0))
              zeff=zeff+(zni*iz(l)**2)/ne
           end do

        end if
     end if   
     
     if (nimp .eq. 2) then
        nIm = sngl(densi(psi,3,0))
     else if (nimp .eq. 3) then
        nIm = sngl( densi(psi,3,0)+densi(psi,4,0))
     else
        nIm = 0.
     end if
     
     
  
     !print*, phi(con), psi
     write(nh,500) sngl(psi), ne, ni, te, ti, nHe, nIm, zeff

  end do
500 format(8e13.6)
  close(nh)


end subroutine omfit

subroutine tglf_input()

  use param
  implicit none

  integer :: con

  double precision :: psi, rhoc, shafr, dshafr, rmaj, rmin
  double precision :: kappa, dkappa, b_unit, p_prime, dpsidr
  double precision :: elong, shift, press
  

  !TGLF Parameters
  double precision :: rmin_loc, rmaj_loc, q_loc, q_prime_loc, p_prime_loc
  double precision :: drmajdx_loc, kappa_loc, s_kappa_loc, delta_loc
  double precision :: s_delta_loc

  double precision, dimension(ncon) :: dpdrs, rhos

  call dpsidrho(dpdrs,rhos)
  
  do con=1,ncon

     psi = psiv(con)

     !Derivative of psi against r
     dpsidr = dpdrs(con)
     
     psi=psiv(con)

     ! Shafranov shift/elongation and derivative in r
     shafr = shift(con,0)
     dshafr = shift(con,1)*dpsidr
     kappa = elong(con,0)
     dkappa = elong(con,1)*dpsidr

     !Major radius 
     rmaj = rcen + shafr
     rmin = maxval(rpts(con,:)) - rmaj

     p_prime = press(psi, 1)*dpsidr

     rhoc = psi/umax
     
     rmin_loc = rmin/amin
     rmaj_loc = rmaj/amin
     q_loc = sfac(con)
     q_prime_loc = q_loc**2 * qp(con) / rmin_loc**2


     drmajdx_loc = dshafr

     kappa_loc = kappa
     s_kappa_loc = rmin*dkappa/kappa

     b_unit = qp(con)*dpsidr/rmin

     p_prime_loc = q_loc * amin**2 * p_prime_loc / (rmin*b_unit**2)
     

  end do
  


end subroutine tglf_input

   

subroutine rfdat()

  use param
  implicit none

  integer :: i, j, nx, ny, nh
  double precision :: fprof
  double precision, dimension(nr,nz) :: psicoord, btcoord, psiNcoord

  nh = 53

  open(unit=nh, file=runname(1:lrunname)//'.rfdat', &
       status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.rfdat'
     stop
  end if


  psicoord = umax - ucoord
  psiNcoord = psicoord/umax

  do i=1,nr
     btcoord(i,:) = 1./rcoord(i)
     do j=1,nz

        !Set f to edge value outside plasma
        if (psiNcoord(i,j) .gt. 1) then

           btcoord(i,j) = btcoord(i,j) * fprof(umax,2)
        else
  
           btcoord(i,j) = btcoord(i,j)* fprof(psicoord(i,j),2)
        end if
        
     end do
  end do
  
  write(nh,2022) nr, nz

  write(nh,*) 'R  data (m)'
  write(nh,2020) (rcoord(i), i=1,nr)
  
  write(nh,*) 'Z  data (m)'
  write(nh,2020) (zcoord(j), j=1,nz)

  write(nh,*) 'psiN  data (T m^2)'
  write(nh,2020) ((psiNcoord(i,j), i=1,nr), j=1,nz)

  write(nh,*) 'Br  data (T)'
  write(nh,2020) ((brcoord(i,j), i=1,nr), j=1,nz)

  write(nh,*) 'Bz  data (T)'
  write(nh,2020) ((bzcoord(i,j), i=1,nr), j=1,nz)
  
  write(nh,*) 'Bt  data (T)'
  write(nh,2020) ((btcoord(i,j), i=1,nr), j=1,nz)

2020 format(5e12.4)
2022 format(2i5)

  close(nh)

end subroutine rfdat 


subroutine vessel()
  ! Set a vessel boundary by adding a set distance to the plasma
  ! boundary in both R and Z
  
  use param
  implicit none

  integer :: i, nh, ind
  double precision, dimension(npts) :: rves, zves, rbdy, zbdy
  double precision :: dwall
  
  nh = 57

  open(unit=nh, file=runname(1:lrunname)//'.vessel', &
       status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.vessel'
     stop
  end if

  rbdy = rpts(1,:)
  zbdy = zpts(1,:)

  dwall = 0.3
  
  ind = maxloc(zbdy,1)
  
  do i=1,npts

     if (zbdy(i).gt.0) then
        zves(i) = zbdy(i)+dwall
     else if (zbdy(i).lt.0) then
        zves(i) = zbdy(i)-dwall
     else
        zves(i) = zbdy(i)
     end if

     if (rbdy(i)-rbdy(ind) .gt. 0) then
        rves(i) = rbdy(i) + dwall
     else if (rbdy(i)-rbdy(ind) .lt. 0) then
        rves(i) = rbdy(i) - dwall
     else
        rves(i) = rbdy(i)
     end if


  end do

  write(nh, *) npts
  write(nh, 2020) (rves(i),zves(i), i=1,npts)


2020 format(2e14.6)
    close(nh)

end subroutine vessel

        

     
