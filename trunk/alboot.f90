      subroutine pdgen(f1d,f2d,falnum,vc,v0,amass,pdal,fl0,fl1,fl2)
!     *************************************************************
!
!   Calculates the psi-derivative of the fusion alpha-particle pressure (PDAL)
!   and velocity space integrals required for effective trapped alpha
!   fraction (FL1 and FL2)
!
      implicit none
      double precision rat,v0,vc,vmin,dv,v,pdal
      double precision f1,f2,f3,f4,fl0,fl1,fl2
      double precision f1d,f2d,falnum,amass
      integer nv,i

      rat=v0/vc
      vmin=v0
      if (vc.lt.vmin) vmin=vc
      dv=vmin/500
      nv=v0/dv
      dv=v0/(nv-1)
      f1=0.0d0
      f2=0.0d0
      f3=0.0d0
      f4=0.0d0
      do i=1,nv
        v=(i-1)*dv+dv/2.
        f1=f1+dv*v**4/(v**3+vc**3)
        f2=f2+dv*v**7/(v**3+vc**3)**2
        f3=f3+dv*v**4/(v**3+vc**3)**2
        f4=f4+dv*v**4
      end do
      fl0=f4*rat**3/(3.*v0**5)
      fl1=f1/(3.*v0**2)
      fl2=f3*v0/(3.*rat**3)
      pdal=7.0e-27*(f1d*f1+f2d*f2)*falnum*amass
      return
   end subroutine pdgen

      
