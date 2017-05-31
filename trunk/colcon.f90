      subroutine colcon(psi,r,z,nr,nz,psicon,flxr,flxz,npts)
!     ******************************************************
!
      use vector
      use contour
!
      implicit none
      type(vec2dl)::con
      double precision psicon
      integer nr,nz,npts,i
      double precision psi(nr,nz),r(nr),z(nz)
      double precision flxr(nr*nz),flxz(nr*nz)
!
      con=conget(psi,r,z,psicon,init=.true.)
      npts=con%npt
      do i=1,npts
        flxr(i)=con%x(i)
        flxz(i)=con%y(i)
      end do
      continue
      return
   end subroutine colcon
