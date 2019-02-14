subroutine geqdsk
!!! Writes out a geqdsk file
  ! In order to include a limiter the standard SCENE R,Z grid is
  ! extended by one grid space on all sides psi (nr,nz) -> psi(nr+2,nz+2)
  ! Uses surfit.f from dierckx to fit a surface inside the boundary to
  ! extend psi out to the edge of the box
  use param
  implicit none

  integer :: i, j, nh, na, con, ij, ind, nr2, nz2
  double precision :: rmax, rmin, dimr, zmax, zmin, dimz
  double precision :: fprof, press, psi, rat, diff
  double precision :: Bv0, psibdy, dpsi, jtor
  integer :: ndat, nlim
  double precision, dimension(:), allocatable :: zbdy, rbdy, rlim, zlim, psii
  double precision, dimension(:), allocatable :: safety, f, ff, p ,pp
  double precision, dimension(:,:), allocatable :: psirz, psi1
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

  !Increase box size for limiter
  nr2 = nr+2
  nz2 = nz+2

  allocate(psirz(nr,nz), psi1(nr,nz))
  
  !First line header, main thing needed in nr and nz (add one for limiter)
  write(nh,2000) 'SCENE ', date, ' : ', time, 'RUN: ', runname, 0, nr2, nz2

  !Change box to include limiter surface defined later
  rmax = maxval(r)+dr
  rmin = minval(r)-dr
  dimr = rmax - rmin
  zmax = maxval(z)+dr
  zmin = minval(z)-dr
  dimz = zmax-zmin

  !Dimensions of R,Z grid
  !print*, dimr, dimz, rcen, rmin
  write(nh,2020) dimr, dimz, rcen, rmin, 0.

  Bv0 = mu0*rodi/(2.*pi*rcen)
  psibdy=umax
  !Psi values and vac B field
  write(nh,2020) r0, 0., 0., psibdy,  Bv0
  write(nh,2020) cur, 0.,0., r0, 0.
  write(nh,2020) 0.,0.,psibdy,0., 0.

  !Gets Psi along midplane (sets Psi outside plasma
  !to be umax and then writes F along midplane

  allocate(psii(nr2), f(nr2), ff(nr2), p(nr2), pp(nr2))

  dpsi=umax/(nr2-1)

  do i=1,nr2
     psii(i) = (i-1)*dpsi
     !print*, psii(i)

     f(i) =  fprof(psii(i),2)
     ff(i) = fprof(psii(i),1)
     p(i) = press(psii(i),0)
     pp(i) = press(psii(i),1)
  end do


  !Writes f
  write(nh,2020) ( f(i), i=1,nr2)


  !Writes Pressure p
  write(nh,2020) (p(i), i=1,nr2)


  write(6,*) 'Pressure written to geqdsk'

  !Writes ff'
  write(nh,2020) (ff(i), i=1,nr2)


  write(6,*) 'ffp written to geqdsk'

  !Writes p'
  write(nh,2020) (pp(i), i=1,nr2)

  jtor=0.0
  call extrap2()
  psirz= umax*1.2
  do i=1,nr
     do j=1,nz
        if (ixout(i,j) .ne. 0) then
           psirz(i,j) = (umax-u(i,j))
           !jtor = jtor + r(i)*press(psirz(i,j),1) + fprof(psirz(i,j),1)/(r(i)*mu0)
        end if

     end do
  end do
  !print*, jtor*dr*dz


  !Extrapolates Psi to edge of grid
  call extrappsi(rmin, rmax, zmin, zmax,psi1)
  write(6,*) 'Smoothing function'
  
  allocate(psigeq(nr2,nz2))
  psigeq = 0.0
  
  !Writes Psi
  !Uses extrapolated values outside bdy,
  !USes actual value inside bdy
  diff=0.
  ij=0
  do i=1,nr
     do j=1,nz
        !Set to extrap val
        psigeq(i+1,j+1) = psi1(i,j)
        
        if (ixout(i,j).ne.0) then
           
           ij=ij+1
           !Set psi inside bdy
           psigeq(i+1,j+1) = umax -u(i,j)
           diff = diff + ( umax - u(i,j) - psi1(i,j) )**2/(umax-u(i,j))**2
           
           !if (psi1(i,j) .lt. 0.01) print*, i,j,umax-u(i,j), psi1(i,j)
           !psi1(i,j) = umax - u(i,j)
        end if

     end do
  end do
  write(6,*) 'Avg fractional diff between Psi and fit is :',  sqrt(diff/ij)

  ! Have a psi that has one more layer grid point on each side


  !Inside points same as before
  psigeq(2:nr+1,2:nz+1) = psi1

  !Extrapolate out to the side first (not incl extra z layers)
  !Left side
  psigeq(1,2:nz2-1) = 2*psigeq(2,2:nz2-1) -psigeq(3,2:nz2-1)

  !Right
  psigeq(nr2,2:nz2-1) = 2*psigeq(nr2-1,2:nz2-1) - psigeq(nr2-2,2:nz2-1)

  !Extrapolate out to top and bottom (incl extra r layers)
  !Top
  psigeq(:,1) = 2*psigeq(:,2) - psigeq(:,3)

  !Bottom
  psigeq(:,nz2) = 2*psigeq(:,nz2-1) - psigeq(:,nz2-2)

  print*, 'Writing psi values to eqdsk, extrapolated out to edge'
  write(nh,2020)  ((psigeq(i,j), i=1,nr2), j=1,nz2)
  write(6,*) 'psi written to geqdsk'

  allocate(safety(nr2))
  !Writes Safety factor (on R,Z grid)
  do i=1,nr2
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

     !print*, i,p(i),pp(i),f(i),ff(i), safety(i)
  end do

  write(nh,2020) (safety(i), i=1,nr2)
  write(6,*) 'safety factor written to geqdsk'

  !Writes boundary points

  nlim=2
  if (ibdry .eq. 2) then
     na=30
     open(unit=na,file='bdy.txt', &
          status='unknown',iostat=ios)
     if(ios.ne.0) then
        write(6,*) 'problem opening bdy.txt for boundary R,Z'
        stop
     endif
     read(na,*)ndat


     allocate(rbdy(ndat),zbdy(ndat))


     do i=1,ndat
        read(na,*)rbdy(i),zbdy(i)


     end do

     close(na)
  else

     ndat=npts
     allocate(rbdy(ndat),zbdy(ndat))
     rbdy(:) = rpts(1,:)
     zbdy(:) = zpts(1,:)

  end if

  ! excl first point in lim, no repeats
  nlim = ndat-1
  !nlim = 1
  write(nh, 2022) ndat, nlim

  write(nh,2020) (rbdy(i), zbdy(i), i=1,ndat)

  write(6,*) 'bdy written to geqdsk'


  !Writes Limiter values
  !Change box range if you change limiter
  allocate(rlim(nlim), zlim(nlim))
  !zlim = zbdy*1.
  !rlim = rbdy
  ind = maxloc(zbdy,1)
  do i=2,npts
     if (rbdy(i)-rbdy(ind) .lt. 0.) then
        rlim(i-1) = rbdy(i)-dr
     else
        rlim(i-1) = rbdy(i)+dr
     end if

     if (zbdy(i) .gt. 0) then
        zlim(i-1) = zbdy(i)+dz
     else
        zlim(i-1) = zbdy(i)-dz
     end if
  end do

    
  !rlim = (/0.9*rmin,1.1*rmax/)

  !zlim = (/1.1*zmin,1.1*zmax/)


  write(nh,2020) (rlim(i), zlim(i), i=1,nlim)

 ! write(nh,'(5e16.9)') ((umax-u(i,j), i=1,nr), j=1,nz)

  deallocate(psii)
  deallocate(safety)
2000 format (6a8, 3i4)
2020 format (5e16.9)
2022 format (2i5)

  close(nh)

end subroutine geqdsk


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
  kr=3
  kz=3


  rl=sngl(rmin)
  ru=sngl(rmax)
  zl=sngl(zmin)
  zu=sngl(zmax)

  !no. of points in the plasma/boundary
  m=0

  do i=1,nr
     do j=1,nz
        if(ixout(i,j).ne.0) m = m + 1
     end do
  end do


 !smoothing factor
  sm = (m - sqrt(2.*m))/1000000
  !sm=0.
  !No. of knots0
  nrest=int(kr+sqrt(m/2.))
  nzest=int(kz+sqrt(m/2.))
!  nrest=nr
!  nzest=nz

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



  lwrk1 =  uw*vw *(2+b1+b2) + 2*(uw+vw+km*(m+nest)+nest-kr-kz) +b2 + 2
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

  !tr(kr+2:nr1-kr-1)=sngl(r(kr+2:nr-kr-1:2))
  !print*, tr(kr+2:nr-kr-1)
  !tz(kz+2:nz1-kz-1)=sngl(z(kz+2:nz-kz-1:2))
  ij=0
  do i=1,nr
     do j=1,nz

        if (ixout(i,j).ne.0) then
           ij=ij+1
           u_in(ij) = sngl(umax - u(i,j))
           r_in(ij) = sngl(r(i))
           z_in(ij) = sngl(z(j))
           w_in(ij) = 1e0
        end if
     end do
  end do

  print*, 2*kr+2, nr1,nrest
  print*, 2*kz+2, nz1,nzest
  print*, m, ep, nmax

  !call surfit(iopt, m, r_in, z_in, u_in, w_in, rl, ru, zl, zu, kr, kz, &
  !     sm, nrest, nzest,nmax,ep, nr1,tr, nz1, tz, c, fp, wrk1,lwrk1,wrk2, &
  !     lwrk2, iwrk, kwrk, ier)

  !print*, tr
  !write(6,*) nr1, tr(1:nr1)
  !write(6,*) nz1, tz(1:nz1)
  !write(6,*) c(1:nr1+nz1)


  !print*, 'Surfit ier:',ier

  !iopt=-1
  call surfit(iopt, m, r_in, z_in, u_in, w_in, rl, ru, zl, zu, kr, kz, &
       sm, nrest, nzest,nmax,ep, nr1,tr, nz1, tz, c, fp, wrk1,lwrk1,wrk2, &
       lwrk2, iwrk, kwrk, ier)

  print*, 'Surfit round 2 ier:',ier

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
