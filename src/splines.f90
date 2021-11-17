module splines
  implicit none
contains

!> Setup routine for bicubic spline of fun[x,y].
!> output of this routine is nx*ny spline coefficients
!> stored in coef(kx,ny)
  subroutine spline2d(fun, x, y, nx, ny, kx, coef)

    implicit none
    double precision fun(*), x(*), y(*), coef(*)
    integer kx, nx, ny
    integer ind, j
    do j = 1, ny
      ind = (j - 1) * kx + 1
      call spline(x, fun(ind), nx, -1.d30, -1.d30, coef(ind))
    end do
  end subroutine spline2d

!> Evaluates bicubic spline: determines fun_new[x_new,y_new].
!> after spline2d has been called to determine coef
!> in the calling program fun_new,run_old, and coef are 2d arrays:
!> fun_old(kx_old,ky_old),fun_new(kx_new,ky_new),coef(kx_old,ky_old)
  subroutine spline2dt(fun_new, x_new, y_new, kx_new, nx_new, ny_new, &
                       fun_old, x_old, y_old, kx_old, nx_old, ny_old, coef)

    implicit none
    double precision fun_new(*), x_new(*), y_new(*)
    double precision fun_old(*), x_old(*), y_old(*)
    double precision coef(*)
    integer kx_new, nx_new, ny_new
    integer kx_old, nx_old, ny_old
    integer i, j, ky_old, ind
    parameter(ky_old=1500)
    double precision ftemp(ky_old), ctemp(ky_old)
    if (ny_old .gt. ky_old) then
      write (6, '("dimensioning problem in spline2dt")')
      stop
    end if
    do i = 1, nx_new
      do j = 1, ny_old
        ind = (j - 1) * kx_old + 1
        call splint(x_old, fun_old(ind), coef(ind), nx_old, x_new(i), &
                    ftemp(j))
      end do
      call spline(y_old, ftemp, ny_old, -1.d30, -1.d30, ctemp)
      do j = 1, ny_new
        call splint(y_old, ftemp, ctemp, ny_old, y_new(j), &
                    fun_new((j - 1) * kx_new + i))
      end do
    end do
  end subroutine spline2dt

!> Use 1d cubic spline on yold[xold] to produce ynew[xnew].
!> y2old(1:nold) is a work array
!> ynew(1:nnew) is the output
  subroutine spline1d(ynew, xnew, nnew, yold, xold, nold, y2old)

    implicit none
    double precision ynew(*), yold(*), xnew(*), xold(*)
    double precision y2old(*)
    double precision yp1, ypn
    integer nnew, nold, i
    yp1 = -1.d30
    ypn = -1.d30
    call spline(xold, yold, nold, yp1, ypn, y2old)
    do i = 1, nnew
      call splint(xold, yold, y2old, nold, xnew(i), ynew(i))
    end do
  end subroutine spline1d

!> Spline routine based upon numerical recipes.
!> This is the setup routine which needs to be called only once
!> splines y as a function of x--both arrays have n elements
!> yp1 and ypn are boundary conditions on the spline
!> yp1=y'[x] at x=x[1]
!> if yp1>=1.e30 then y''[x]=0 at x=x[1] is used
!> if yp1<=-1.e30 then y'[x[1]] is calculated from first four points
!> ypn=y'[x] at x=x[n]
!> if ypn>=1.e30 then y''[x]=0 at x=x[n] is used
!> if ypn<=-1.e30 then y'[x[n]] is calculated from last four points
!> y2[1:n] is calculated array of the second derivatives of the
!> interpolating function at the x[i]
  subroutine spline(x, y, n, yp1, ypn, y2)
    implicit none
    integer nmax, n
    parameter(nmax=7000)
    double precision x(*), y(*), y2(*), yp1, ypn
    double precision u(nmax)
    double precision yp1t, ypnt, sig, p, qn, un
    integer k, i
    if (n .gt. nmax) then
      write (6, *) 'spline;  dimensional error n=', n, ' nmax=', nmax
      stop
    end if
    if (yp1 .gt. .99d30) then
      y2(1) = 0.d0
      u(1) = 0.d0
    else if (yp1 .lt. -.99d30) then
      yp1t = (3 * x(1)**2 + x(2) * x(3) + &
              x(2) * x(4) + x(3) * x(4) - 2 * x(1) * (x(2) + x(3) + x(4))) * &
             y(1) / ((x(1) - x(2)) * (x(1) - x(3)) * (x(1) - x(4))) + &
             (-x(1)**2 + x(1) * x(3) + x(1) * x(4) - x(3) * x(4)) * y(2) / &
             ((x(1) - x(2)) * (x(2) - x(3)) * (x(2) - x(4))) + &
             (x(1)**2 - x(1) * x(2) - x(1) * x(4) + x(2) * x(4)) * y(3) / &
             ((x(1) - x(3)) * (x(2) - x(3)) * (x(3) - x(4))) + &
             (-x(1)**2 + x(1) * x(2) + x(1) * x(3) - x(2) * x(3)) * y(4) / &
             ((x(1) - x(4)) * (x(2) - x(4)) * (x(3) - x(4)))
      y2(1) = -0.5d0
      u(1) = (3.d0 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1t)
    else
      y2(1) = -0.5d0
      u(1) = (3.d0 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
    end if
    do i = 2, n - 1
      sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
      p = sig * y2(i - 1) + 2.d0
      y2(i) = (sig - 1.d0) / p
      u(i) = (6.d0 * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1)) &
                      / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do
    if (ypn .gt. .99d30) then
      qn = 0.d0
      un = 0.d0
    else if (ypn .lt. -.99d30) then
      ypnt = (-(x(-2 + n) * x(-1 + n)) + x(-2 + n) * x(n) + x(-1 + n) * x(n) - x(n)**2) * &
             y(-3 + n) / &
             ((-x(-3 + n) + x(-2 + n)) * (-x(-3 + n) + x(-1 + n)) * &
              (-x(-3 + n) + x(n))) + &
             (x(-3 + n) * x(-1 + n) - x(-3 + n) * x(n) - x(-1 + n) * x(n) + x(n)**2) * &
             y(-2 + n) / &
             ((-x(-3 + n) + x(-2 + n)) * (-x(-2 + n) + x(-1 + n)) * &
              (-x(-2 + n) + x(n))) + &
             (-(x(-3 + n) * x(-2 + n)) + x(-3 + n) * x(n) + x(-2 + n) * x(n) - x(n)**2) * &
             y(-1 + n) / &
             ((-x(-3 + n) + x(-1 + n)) * (-x(-2 + n) + x(-1 + n)) * &
              (-x(-1 + n) + x(n))) + &
             (x(-3 + n) * x(-2 + n) + x(-3 + n) * x(-1 + n) + x(-2 + n) * x(-1 + n) - &
              2 * (x(-3 + n) + x(-2 + n) + x(-1 + n)) * x(n) + 3 * x(n)**2) * y(n) / &
             ((-x(-3 + n) + x(n)) * (-x(-2 + n) + x(n)) * (-x(-1 + n) + x(n)))
      qn = 0.5d0
      un = (3.d0 / (x(n) - x(n - 1))) * (ypnt - (y(n) - y(n - 1)) / (x(n) - x(n - 1)))
    else
      qn = 0.5d0
      un = (3.d0 / (x(n) - x(n - 1))) * (ypn - (y(n) - y(n - 1)) / (x(n) - x(n - 1)))
    end if
    y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1.d0)
    do k = n - 1, 1, -1
      y2(k) = y2(k) * y2(k + 1) + u(k)
    end do
  end subroutine spline

!> Cubic spline evaluator--spline must be called first to evaluate y2a.
!> ya is a function of xa--both are arrays of length n
!> ya2[1:n] contains spline coefficients calculated in spline
!> x is the argument of y[x] where y is to be evaluated
!> y=y[x] is the returned value

  subroutine splint(xa, ya, y2a, n, x, y)

    implicit none
    double precision xa(*), ya(*), y2a(*)
    double precision x, y
    double precision a, b, h
    integer n, klo, khi, k
    klo = 1
    khi = n
1   if (khi - klo .gt. 1) then
      k = (khi + klo) / 2
      if (xa(k) .gt. x) then
        khi = k
      else
        klo = k
      end if
      goto 1
    end if
    h = xa(khi) - xa(klo)
    if (abs(h) .le. 1e-8) then
      write (6, *) 'bad xa input.'
      stop
    end if
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.d0
  end subroutine splint

!> zspline integrates the cubic spline of ya[xa]
!> it assumes that spline has already been called to evaluate y2a
!> xa[1:n], ya[1:n], y2a[1:n]
!> In Mathematica notation z[i]=zg+Integrate[y[x],{x,x[ng],x[i]}]
!> where both zg and ng are input quantities,
!> za is calculated here.
!> za can then be used in zsplint to determine z at a
!> specific x
  subroutine zspline(xa, ya, y2a, n, zg, ng, za)
    implicit none
    double precision xa(*), ya(*), y2a(*), zg, za(*)
    integer n, ng
    integer j
    double precision const
    if (ng .lt. 0 .or. ng .gt. n) stop 'zspline: wrong ng'

    za(1) = 0.d0
    do j = 2, n
      za(j) = za(j - 1) + 0.5d0 * (xa(j) - xa(j - 1)) * (ya(j) + ya(j - 1)) - &
              (xa(j) - xa(j - 1))**3 * (y2a(j) + y2a(j - 1)) / 24.d0
    end do

    const = zg - za(ng)
    do j = 1, n
      za(j) = za(j) + const
    end do
  end subroutine zspline

!> Evaluate cubic spline to determine function (y),
!> derivative (yp) and integral (z) at location x.
!> first spline must be called to obtain y2a and
!> zspline must be called to obtain za
  subroutine zsplint(xa, ya, y2a, za, n, x, y, yp, z)
    implicit none
    double precision xa(*), ya(*), y2a(*), za(*)
    integer n
    double precision x, y, yp, z
    integer klo, khi, k
    double precision h, a, b
    klo = 1
    khi = n
1   if (khi - klo .gt. 1) then
      k = (khi + klo) / 2
      if (xa(k) .gt. x) then
        khi = k
      else
        klo = k
      end if
      goto 1
    end if
    h = xa(khi) - xa(klo)
    if (abs(h) .lt. 1e-8) then
      write (6, *) 'bad xa input.'
      stop
    end if
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.d0
    yp = (ya(khi) - ya(klo)) / h - (3.d0 * a**2 - 1.d0) / 6.d0 * h * y2a(klo) &
         + (3.d0 * b**2 - 1.d0) / 6.d0 * h * y2a(khi)
    z = za(klo) + 0.5d0 * h * ya(klo) * (1.d0 - a**2) + 0.5d0 * ya(khi) * h * b**2 + &
        y2a(klo) * h**3 * (2.d0 * a**2 - a**4 - 1.d0) / 24.d0 - &
        y2a(khi) * h**3 * (2.d0 * b**2 - b**4) / 24.d0
  end subroutine zsplint
end module splines
