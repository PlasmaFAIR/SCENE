      subroutine lynton
!     *****************
!
!  Generates various paramters for Lynton for fast particle physics studies.
!
      use param
      implicit none
      integer i,nh,j
      double precision nelec,nion,psi,ss,mdens,tion
      double precision dense,densi,tempi
!
      nh=229
!!$      open(unit=nh,file=runname(1:lrunname)//'.lynton', &
!!$           status='unknown',iostat=ios)
!!$      if(ios.ne.0) then
!!$         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.lynton'
!!$         stop
!!$      endif
      write(nh,*)' The following data is psi_norm, mass density, Ti and ni'
      do i=ncon,1,-1
        psi=psiv(i)
        ss=psi/umax
        nelec=dense(psi,0)
        nion=densi(psi,1,0)
        tion=tempi(psi,1,0)
        mdens=me*dense(psi,0)
        write(6,*)' **** mdens=',mdens
        do j=1,nimp+1
          mdens=mdens+zmas(j)*mp*densi(psi,j,0)
          write(6,*)' mdens=',mdens
        end do
        write(nh,10)ss,mdens,tion,nion
      end do
10    format(4e14.5)
      close(nh)
  end subroutine lynton
