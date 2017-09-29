subroutine getdata

  ! Puts data in the needed in the output file


  use param
  implicit none

  double precision :: rr, zz, rat
  double precision :: fprof, tempe, tempi, psi,psi_n, bp, B_tor, B_pol, B_tot
  double precision :: press, densi, dense, safety, ffp, pp
  double precision :: J_tot, J_ext, J_bs, J_di, J_ps, J_nbi, J_ext2
  character(8) :: date
  character(10) :: time
  double precision :: flux_r, flux_z

  
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

  do i=1,nr
     if (ixout(i,nsym) .le. 0) cycle
     psi = umax-u(i,nsym)
     rr = r(i)
     zz = 0.
     B_tor = sngl(fprof(psi,2)/rr)
     B_tot = sngl(sqrt(B_tor**2 + (bp(rr,zz)**2)) )

     con = 1

     !ncon=1 is the outermost flux surface
     do j=1,ncon-1

        if (psiv(j) .lt. psi) exit

        con = j
     end do

      
     if (con .eq. ncon-1) con = con-1
        
     rat = (psi-psiv(con))/(psiv(con+1) - psiv(con) )


     !Linearly interpolate q from flux surface grid to mesh grid
     safety = sngl( sfac(con) + rat*(sfac(con+1) - sfac(con)))


     write(nh,14) rr,dense(psi,0),densi(psi,1,0), B_tor,bp(rr,zz),B_tot, safety, tempe(psi,0)/1000., &
          press(psi,0)

  end do
  
14 format(ES14.6e2,ES14.6e2, ES14.6e2 )
   close(nh)

   write(6,*) 'Overall data written'




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

   write(6,*) 'Flux surface data written'



   
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
      if (ixout(i,nsym) .le. 0) cycle

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

   !J profile with psi values
   open(unit=nh, file=runname(1:lrunname)//'_currentprof.dat', &
        status='unknown', iostat=ios)
   if (ios .ne. 0) then
      write(6,*) 'problem opening ',runname(1:lrunname)//'_currentprof.dat'
      stop
   end if

   !Used to find magnetic axis

   
   do i=1,nr
      if (ixout(i,nsym) .le. 0) cycle 
      !Go to magnetic axis and then work outwards
      if (r(i) .lt. r0) cycle

      psi = umax - u(i,nsym)
      psi_n = psi/umax

      ffp = -sngl(fprof(psi,1))
      pp = -sngl(press(psi,1))
      write(nh,18) psi_n, ffp, pp


   end do
    
!    psi = umax - u(nr,nsym)
!    ffp = sngl(fprof(psi,1))
!    pp = sngl(press(psi,1))
!   write(nh,18) psi, ffp, pp
   
18 format(6e14.6)


   close(nh)

   write(6,*) 'current prof data written'


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
   
      

   write(6,*) 'locust ne data written'


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

   write(6,*) 'locust Te data written'

   
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

   write(6,*) 'locust Ti data written'


      
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
   
   write(6,*) 'psi data written'
end subroutine getdata

   
   
    
