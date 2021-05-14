module caxe_output
  implicit none
contains
      subroutine caxout
!     *****************
!
!  Generates various paramters for CAXE code for translation for KINX
! MHD stability code
!
      use param
      use profiles_mod, only : fprof, press
      implicit none
      integer i,nh
      double precision psi,pp,pf,qr,cu,bv0
!
      nh=117
      open(unit=nh,file=runname(1:lrunname)//'.caxe', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.caxe'
         stop
      endif
      write(nh,*)npts
      do i=1,npts
       write(nh,'(1x,2(1pe12.4))') rpts(1,i),zpts(1,i)
      end do
      write(nh,*)ncon
      write(nh,*)' '
      bv0=mu0*rodi/(2.*pi*r0)
      do i=ncon,1,-1
        psi=psiv(i)
        pp=mu0*press(psi,1)/bv0
        pf=fprof(psi,1)/bv0
        qr=0.
        cu=-pf
        write(nh,' (1x,5(1pe12.4))')psi,pp,pf,qr,cu
      end do
      return
 end subroutine caxout
end module caxe_output
