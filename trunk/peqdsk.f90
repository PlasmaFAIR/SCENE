
module peqdsk_output
  implicit none
contains
subroutine peqdsk
  ! Writes out a peqdsk file that can be read in by GACODE profile_gen tool
  ! Contains information about kinetic profiles

  use param
  implicit none

  integer :: con, nh
  double precision :: dense, densi, tempe, tempi, psi
  real :: ne, ni, te, ti, psin
  real :: dne, dni, dte, dti

  nh = 71

  open(unit=nh, file=runname(1:lrunname)//'.peqdsk', &
       status='unknown', iostat=ios)

  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.peqdsk'
     stop
  end if

  write(nh, *) ncon, 'psinorm ne(10^20/m^3) dne/dpsiN'

  do con=ncon, 1, -1

     psi = psiv(con)
     psin = sngl(psi/umax)
     ne = sngl(dense(psi, 0)) * 1.e-20
     dne = sngl(dense(psi, 1) * umax) * 1.e-20

     write(nh, *) psin, ne, dne

  end do

  write(nh, *) ncon, 'psinorm te(keV) dte/dpsiN'

  do con=ncon, 1, -1

     psi = psiv(con)
     psin = sngl(psi/umax)
     te = sngl(tempe(psi, 0)) * 1.e-3
     dte = sngl(tempe(psi, 1) * umax) * 1.e-3

     write(nh, *) psin, te, dte

  end do

  write(nh, *) ncon, 'psinorm ni(10^20/m^3) dni/dpsiN'

  do con=ncon, 1, -1

     psi = psiv(con)
     psin = sngl(psi/umax)
     ni = sngl(densi(psi, 1, 0)) * 1.e-20
     dni = sngl(densi(psi, 1, 1) * umax) * 1.0e-20

     write(nh, *) psin, ni, dni

  end do

  write(nh, *) ncon, 'psinorm ti(keV) dti/dpsiN'

  do con=ncon, 1, -1

     psi = psiv(con)
     psin = sngl(psi/umax)
     ti = sngl(tempi(psi, 1, 0)) * 1.0e-3
     dti = sngl(tempi(psi, 1, 1) * umax) * 1.0e-3

     write(nh, *) psin, ti, dti

  end do

  close(nh)


  end subroutine peqdsk
  
end module peqdsk_output
