module tokameq_output
  implicit none
contains
      subroutine tokameq
!     ******************
!
!  Code top calculate profile shape input for TOKAMEQ code
!
      use param
      use profiles_mod, only : press, fprof
      implicit none
      integer nh,k,j,i
      double precision psi,fpsi,jphi
!
      nh=111
      open(unit=nh,file=runname(1:lrunname)//'.tokameq', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.tokameq'
         stop
      endif
      write(nh,*)' normalised psi: 0 at edge, 1 on axis'
      do k=1,ncon
        write(nh,*)(psiv(k)-psiv(1))/(psiv(ncon)-psiv(1))
      end do
      write(nh,*)' pressure (Nm**-2)'
      do k=1,ncon
        psi=psiv(k)
        write(nh,*)press(psi,0)
      end do
      write(nh,*)' Delta=5*R*Bphi-Irod/10**6 (MA)'
      do k=1,ncon
        psi=psiv(k)
        fpsi=fprof(psi,2)
        write(nh,*)5.*fpsi-rodi*1.0d-6
      end do
      write(nh,*)' dp/da'
      do k=1,ncon
        psi=psiv(k)
        write(nh,*)press(psi,1)*(psiv(ncon)-psiv(1))
      end do
      write(nh,*)' I dI/da'
      do k=1,ncon
        psi=psiv(k)
        write(nh,*)fprof(psi,1)*(psiv(ncon)-psiv(1))
      end do
      write(nh,*)' R (m), Jphi (MAm**-2)'
      j=nsym
      do i=1,nr
        psi=umax-u(i,j)
        if (psi.lt.umax) then
          jphi=-(r(i)*press(psi,1)+fprof(psi,1)/(r(i)*mu0))*1.0d-6
        else
          jphi=0.
        end if
        write(nh,*)r(i),jphi
      end do
      write(nh,*)' (R,Z) of boundary'
      do j=1,npts
        write(nh,*)rpts(1,j),zpts(1,j)
      end do
      close(nh)
      return
    end subroutine tokameq

end module tokameq_output
