subroutine getdata

  ! Puts data in the needed in the output file


  use param
  implicit none

  double precision :: rr, zz, rat
  double precision :: fprof, tempe, tempi, psi,psi_n, bp, B_tor, B_pol, B_tot
  double precision :: press, densi, dense, safety, ffp, fsi, pp, rho
  double precision :: J_tot, J_ext, J_bs, J_di, J_ps, J_nbi, J_ext2, jnb
  character(8) :: date
  character(10) :: time
  double precision :: flux_r, flux_z, te

  double precision, dimension(ncon) :: rhoflux, Rflux
  double precision :: rhomax, jtot

  real, dimension(ncon) :: rhos, dpdrs

  integer :: nh, i, j, con, flag


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
      jtot = +rr*pp + ffp/(mu0*rr)
      jnb = J_nb(con)*fsi/rr

      te = tempe(psi,0)
      write(nh,18) rhoflux(con), epsv(con), jtot, jnb, te


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
      write(nh, 20) psiv(i), densi(psiv(i),1,0)
   end do

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


end subroutine getdata


subroutine popcon()

  use param
  implicit none


  double precision :: rr, zz, rat
  double precision :: fprof, tempe, tempi, psi,psi_n, bp, B_tor, B_pol, B_tot
  double precision :: press, densi, dense, safety, ffp, pp

  integer :: nh, i, j, con, flag


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
  write(nh,510) 'Shine through', bmshine, ' '
  write(nh,510) 'Q',sngl( ((pfus*1.0d-6*5)+bmfus)/paux), ' '

  write(nh,510) 'N_gw', sngl(negw), ' '

  write(nh,510) 'H_IPB98(y1)', sngl(hipb98y1), ' '
  write(nh,510) 'H_IPB98(y2)', sngl(hipb98y2), ' '
  write(nh,510) 'Tau_e', sngl(taue*1.0e3), 'ms'
  write(nh,510) 'Tau_h', sngl(tauh), ' '

  write(nh,510) 'Area', sngl(area), 'm^2'
  write(nh,510) 'Volume', sngl(vol), 'm^3'
  write(nh,510) 'li(3)', sngl(rli3), ' '
  write(nh,510) 'li(2)', sngl(rli2), ' '
  write(nh,510) 'Energy', sngl(conft/1.0e6), 'MJ'

500 format(' ',A12, ' , ', E11.4 , ' , ', A5)
510 format(' ',A12, ' , ', F11.4 ,' , ',A5)
520 format(' ',A11, ' , ' ,I11, ' , ', A5)

  close(nh)

end subroutine popcon
