subroutine geqdsk

  use param
  implicit none

  integer :: i, j, nh, na, con, ij
  double precision :: rmax, rmin, dimr, zmax, zmin, dimz
  double precision :: fprof, press, psi, rat, diff, sqrt
  double precision :: Bv0, s
  integer :: ndat, nlim
  double precision, dimension(:), allocatable :: zbdy, rbdy, rlim, zlim, psii
  double precision, dimension(:), allocatable :: safety, f, ff, p ,pp
  double precision, dimension(nr,nz) :: psi1
  character(8) :: date
  character(10) :: time

  !Writes GEQDSK file

  nh = 49

  open(unit=nh, file=runname(1:lrunname)//'.geqdsk', &
        status='unknown', iostat=ios)
  if (ios .ne. 0) then
     write(6,*) 'problem opening ',runname(1:lrunname)//'.geqdsk'
     stop
  end if

  call DATE_AND_TIME(date, time)


  !First line header, main thing needed in nr and nz
  write(nh,2000) 'SCENE ', date, ' : ', time, 'RUN: ', runname, 0, nr, nz
  rmax = maxval(r)
  rmin = minval(r)
  dimr = rmax - rmin
  zmax = maxval(z)
  zmin = minval(z)
  dimz = zmax-zmin

  !Dimensions of R,Z grid
  write(nh,2020) dimr, dimz, rcen, rmin, 0.

  Bv0 = mu0*rodi/(2.*pi*rcen)

  !Psi values and vac B field
  write(nh,2020) r0, 0., 0., umax,  Bv0
  write(nh,2020) cur, 0.,0., r0, 0.
  write(nh,2020) 0.,0.,umax,0., 0.

  !Gets Psi along midplane (sets Psi outside plasma
  !to be umax and then writes F along midplane

  allocate(psii(nr), f(nr), ff(nr), p(nr), pp(nr))

  do i=1,nr
     psii(i) = umax - u(i,nsym)
     if (ixout(i,nsym) .gt. 0) then
        f(i) =  fprof(psii(i),2)
        ff(i) = fprof(psii(i),0)
        p(i) = press(psii(i),0)
        pp(i) = press(psii(i),1)
     else
        f(i) = 0.
        ff(i) = 0.
        p(i) = 0.
        pp(i) = 0.
     end if

  end do

  !Writes f
  write(nh,2020) ( f(i), i=1,nr)


  !Writes Pressure p
  write(nh,2020) (p(i), i=1,nr)


  write(6,*) 'Pressure written to geqdsk'

  !Writes ff'
  write(nh,2020) (ff(i), i=1,nr)


  write(6,*) 'ffp written to geqdsk'

  !Writes p'
  write(nh,2020) (pp(i), i=1,nr)

  !call gcv()

  !Extrapolates Psi to edge of grid
  call extrappsi(rmin, rmax, zmin, zmax,psi1)
  write(6,*) 'Smoothing function'

  !Writes Psi
  diff=0.
  ij=0
  do i=1,nr
     do j=1,nz
        ij=ij+1
        if (ixout(i,j).eq.1) then
           diff = diff + ((umax - u(i,j)) - psi1(i,j))**2
           psi1(i,j) = umax - u(i,j)
        end if

     end do
  end do
  write(6,*) 'Avg diff between Psi and fit is :',  sqrt(diff/ij)
  write(nh,2020)  ((psi1(i,j), i=1,nr), j=1,nz)
  write(6,*) 'psi written to geqdsk'
  allocate(safety(nr))
  !Writes Safety factor (on R,Z grid)
  do i=1,nr
     psi = psii(i)
     con = 1

     !ncon=1 is the outermost flux surface
     do j=1,ncon-1

        if (psiv(j) .lt. psi) exit

        con = j
     end do


     if (con .eq. ncon-1) con = con-1

     rat = (psi-psiv(con))/(psiv(con+1) - psiv(con) )

     !Linearly interpolate q from flux surface grid to mesh grid
     safety(i) = sngl( sfac(con) + rat*(sfac(con+1) - sfac(con)))


  end do

  write(nh,2020) (safety(i), i=1,nr)
  write(6,*) 'safety factor written to geqdsk'

  !Writes boundary points
  na=30

  open(unit=na,file='bdy.txt', &
       status='unknown',iostat=ios)
  if(ios.ne.0) then
     write(6,*) 'problem opening bdy.txt for boundary R,Z'
     stop
  endif
  read(na,*)ndat
  nlim = 1
  write(nh, 2022) ndat, nlim

  allocate(rbdy(ndat),zbdy(ndat))


  do i=1,ndat
     read(na,*)rbdy(i),zbdy(i)


  end do

  close(na)

  write(nh,2020) (rbdy(i), zbdy(i), i=1,ndat)

  write(6,*) 'bdy written to geqdsk'

!!!!!!!!
  !Writes Limiter values

  rlim = 1.

  zlim = 0.


  write(nh,2020) rlim, zlim

 ! write(nh,'(5e16.9)') ((umax-u(i,j), i=1,nr), j=1,nz)

  deallocate(psii)
  deallocate(safety)
2000 format (6a8, 3i4)
2020 format (5e16.9)
2022 format (2i5)

  close(nh)

end subroutine geqdsk


subroutine gcv()
  !Subroutines uses gcv for thin plate spline fitting of
  !Psi data

  use param
  implicit none

  integer :: i,j, ij
  double precision :: dble
  double precision, dimension(:), allocatable :: psi_in,y, adiag, dout, coef
  double precision, dimension(:), allocatable :: svals, work
  integer, dimension(4) :: iout
  double precision, dimension(:,:), allocatable :: rz, cov,tbl,auxtbl
  double precision, dimension(2) :: lamlim

  integer, dimension(:), allocatable :: iwork
  integer :: nobs,pts, m, ntbl, dim, info, lds, ncov,ncts,ldtbl,lwa,liwa, job
!  double precision
  !no. of points in the plasma/boundary
  pts=0

  do i=1,nr
     do j=1,nz
        if(ixout(i,j).eq.1) pts = pts + 1
     end do
  end do

  allocate(rz(pts,2), psi_in(pts))


  ij=0
  do i=1,nr
     do j=1,nz
        if (ixout(i,j).eq.1) then
           ij=ij+1
           psi_in(ij) = umax - u(i,j)
           rz(ij,1) = r(i)
           rz(ij,2) = z(j)

        end if
     end do
  end do


  !no. of ovbservations points
  nobs=pts
  lds = pts
  dim=2
  m = 2 !If you change m, make sure you change ncts
  ncov = 0
  ntbl = 80
  ncts = 3 !choose(dim+m-1, dim)
  ldtbl = ntbl
  lwa = nobs*(2+ncts+nobs) + nobs
  liwa = 2*nobs + nobs - ncts

!   job         integer with decimal expansion abdc
!		if a is nonzero then predictive mse is computed
!		   using adiag as true y
!		if b is nonzero then user input limits on search
!		   for lambda hat are used
!		if c is nonzero then adiag will be calculated
  !		if d is nonzero then there are replicates in the design
  !                (no duplicate (r,z) )

  job = 0000


  allocate(y(nobs), adiag(pts), dout(5), svals(pts), coef(pts+ncts))
  allocate(auxtbl(3,3), tbl(ldtbl,3), work(lwa), iwork(liwa) )

  adiag=0
  lamlim = 0
  y = psi_in


  call dtpss(psi_in, pts,nobs,dim,m,cov,lds, ncov,y,ntbl,dout, iout, coef, &
       tbl,ldtbl,auxtbl, work, lwa, iwork, liwa, job,  info)

end subroutine gcv



subroutine extrappsi(rmin, rmax, zmin,zmax,psi_out)
  !Subroutine to fit a surface to Psi and then fill in
  !values of Psi outside plasma boundary


  use param
  implicit none

  integer ::i,j, ij, m
  integer :: iopt, kr, kz, nrest, nzest, nr1, nz1
  integer :: lwrk1, lwrk2, kwrk, ier
  integer :: nmax, uw, vw,ww, km, nest,br,bz, b1,b2
  double precision ::  rmin, rmax, zmin, zmax, dble
  real sm, fp, ep, sngl, rl, ru, zl, zu
  double precision, dimension(nr,nz), intent(out) :: psi_out

  real, dimension(:), allocatable :: tr, tz, c, wrk1, wrk2, wrk3
  integer, dimension(:), allocatable :: wrki,iwrk
  real, dimension(:), allocatable :: r_in, z_in, w_in, u_in
  real, dimension(:), allocatable :: r_o,z_o, psi_extrap

  !iopt =-1 =>least sq, 0 (1) smoothing (restart)
  iopt=0

  !kr and kz choose order (3=cubic)
  kr = 3
  kz=3


  rl=sngl(rmin)
  ru=sngl(rmax)
  zl=sngl(zmin)
  zu=sngl(zmax)

  !no. of points in the plasma/boundary
  m=0

  do i=1,nr
     do j=1,nz
        if(ixout(i,j).eq.1) m = m + 1
     end do
  end do


 !smoothing factor
  sm = (m - sqrt(2.*m))/30000
  !No. of knots0
  nrest=int(kr+sqrt(m/2.))
  nzest=int(kz+sqrt(m/2.))

  !nrest=2*kr+2 !+2
  !nzest=2*kz+2 !+2

  !precision (rank)
  ep = 1e-8

  uw = nrest-kr-1
  vw = nzest-kz-1
  km = max(kr,kz)+1
  nest = max(nrest,nzest)
  br = kr*vw +kz+1
  bz = kz*uw +kr+1
  ww = max(uw,vw)


  nr1 = nrest
  nz1 = nzest
  nmax = max(nr1,nz1)

  if (br.le.bz) then
     b1 = br
     b2 = b1+vw-kz
  else
     b1 = bz
     b2 = b1+uw-kr
  end if



  lwrk1 =  uw*vw *(2+b1+b2) + 2*(uw+vw+km*(m+nest)+nest-kr-kz) +b2 + 1
  !lwrk1 = (7*uw*vw + 25*ww) * (ww+1) + 2 * (uw+vw+4*plas)+23*ww+56

  lwrk2 = uw*vw*(b2+1)+b2

  kwrk = m+(nrest-2*kr-1)*(nzest-2*kz-1) + 1

  allocate( u_in(m), r_in(m), z_in(m), w_in(m))
  allocate(tr(nmax), tz(nmax), c((uw*vw)), wrk1(lwrk1), wrk2(lwrk2), iwrk(kwrk))
  allocate(r_o(nr), z_o(nz))
  allocate(psi_extrap((nr*nz)))

  !tr(kr+2) = sngl(r(int(nr/3)))
  !tz(kz+2) = sngl(z(int(nz/3)))

  !tr(nmax-kr-1) = sngl(r(int(2*nr/3)))
  !tz(nmax-kz-1) = sngl(r(int(2*nr/3)))

  !tz(
  ij=0
  do i=1,nr
     do j=1,nz

        if (ixout(i,j).eq.1) then
           ij=ij+1
           u_in(ij) = sngl(umax - u(i,j))
           r_in(ij) = sngl(r(i))
           z_in(ij) = sngl(z(j))
           w_in(ij) = 1e0
        end if
     end do
  end do

  call surfit(iopt, m, r_in, z_in, u_in, w_in, rl, ru, zl, zu, kr, kz, &
       sm, nrest, nzest,nmax,ep, nr1,tr, nz1, tz, c, fp, wrk1,lwrk1,wrk2, &
       lwrk2, iwrk, kwrk, ier)

  !write(6,*) nr1, tr(1:nr1)
  !write(6,*) nz1, tz(1:nz1)
  !write(6,*) c(1:nr1+nz1)


  do i=1,nr
     r_o(i) = sngl(r(i))
  end do

  do i=1,nz
     z_o(i) = sngl(z(i))
  end do

  kwrk = nr+nz
  lwrk1 = nr*(kr+1) + nz*(kz+1)

  allocate(wrki(kwrk), wrk3(lwrk1))



  call bispev(tr,nr1,tz,nz1,c,kr,kz,r_o,nr,z_o,nz,psi_extrap, &
       wrk3,lwrk1,wrki, kwrk,ier)


  ij=0
  do i=1,nr
     do j=1,nz
        ij = ij+1
        psi_out(i,j) = dble(psi_extrap(ij))
     end do
  end do
  write(6,*) 'Extrap psi: ier = ',ier, ' s is : ',sm, ' fp is ',fp, ' for m points: ',m

  write(6,*) 'Bispev in Extrappsi ier = ',ier


  deallocate( u_in, r_in, z_in, w_in)
  deallocate(tr, tz, c, wrk1, wrk2)
  deallocate(wrki,wrk3)


end subroutine extrappsi
