subroutine getdata

  ! Puts data in the needed in the outoput file


  use param
  implicit none

  double precision :: rr, zz, rat
  double precision :: fprof, tempe, dense, psi, bp, B_tor, B_pol, B_tot
  double precision :: press, densi, safety
  double precision :: J_tot, J_ext, J_bs, J_di, J_ps
 
  double precision :: flux_r, flux_z

  
  integer :: nh, i, j, con


  nh = 49


  
  open(unit=nh, file=runname(1:lrunname)//'_overall.dat', &
       status='unknown', iostat=ios)
  if (ios.ne.0) then
     write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_overall.dat'
     stop
  endif

  write(nh,*) 'R(m), ne(m^-3), ni(m^_3), B_T (T) , B_P (T),  B(T), Sfac, T(keV), Pressure'

  do i=1,nr
     if (ixout(i,nsym) .le. 0) cycle
     psi = umax-u(i,nsym)
     rr = r(i)
     zz = 0.
     B_tor = sngl(fprof(psi,2)/rr**2)
     B_tot = sngl(sqrt(B_tor**2 + (bp(rr,zz)**2)) )

     con = 1
     
     do j=1,ncon-1

        if (psiv(j) .lt. psi) exit

        con = j
     end do

      
     if (con .eq. ncon-1) con = con-1
        
     rat = (psi-psiv(con))/(psiv(con+1) - psiv(con) )


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

   write(nh,*) 'R(m), J_tot (kA m^-2), J_ext, J_bs, J_di, J_ps'
   
   do i=1,nr
      if (ixout(i,nsym) .le. 0) cycle

      
      J_tot = sngl(gradj(i,nsym)/1000.)
      J_ext = sngl(exph(i,nsym)/1000.)
      J_bs = sngl(bsph(i,nsym)/1000.)
      J_di = sngl(diph(i,nsym)/1000.)
      J_ps = sngl(psph(i,nsym)/1000.)

      rr = r(i)


      write(nh,16) rr,J_tot,J_ext,J_bs,J_di,J_ps

   end do
16 format(6e14.6)







     

end subroutine getdata

   
   
    
