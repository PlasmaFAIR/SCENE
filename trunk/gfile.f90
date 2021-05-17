module gfile_output
  implicit none
contains
      subroutine gfile
!     ****************
!
!  Prepares a G-file output (mimics EFIT output)
!
!    H R Wilson 4/7/02
!
      use param
      use profiles_mod, only : fprof, press
      implicit none
      integer :: lun,idum,itime,i,j,ij
      character(len=34) :: adum1
      character(len=14) :: adum2
      double precision :: xdim,zdim,b0vac_r,rzero,zgeo
      double precision :: rmag,zmag,psimag,psi_lcfs,b0vac
      double precision :: curnt,xdum,dpsi,rat,psi
      double precision :: arr(nr)
!
!  2 dummy integers
      itime=0
      idum=0
!  Dummy real
      xdum=0.
      adum1='Output from SCENE'
      adum2='    '//runname(1:lrunname)
      lun=72
      open(unit=lun,file=runname(1:lrunname)//'.gfile', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.gfile'
         stop
      endif
! nr is number of R mesh points
! nz is number of Z mesh points
    write(lun,500) adum1,itime,adum2,nr,nz
      xdim=r(nr)-r(1)
      zdim=z(nz)-z(1)
      b0vac_r=rcen
      rzero=r(1)
      zgeo=0.
!  xdim is maximum R on grid - minimum R on grid, in metres
!  zdim is maximum Z on grid - minimum Z on grid, in metres
!  B0vac_r is the geometric axis in metres
!  Zgeo is the Z-position of the geometric axis (set to zero, cos up-down sym)
    write(lun,510)xdim,zdim,b0vac_r,rzero,zgeo
      rmag=r0
      zmag=0.
      psimag=0.
      psi_lcfs=umax
      b0vac=mu0*rodi/(2.*pi*b0vac_r)
! r0 is the magnetic axis major radius (metres)
! Zmag is Z position of magnetic axis (set to zero 'cos it's up-down sym)
! psimag is thje value of psi at the magnetic axis
! psi_lcfs is the value of psi at the plasma edge. Psi here is in Wb/(2pi)
! so that BR=(1/R) (dpsi/dZ), BZ=-(1/R) (dpsi/dR)
! b0vac -s vacuum field (in Tesla) at b0vac_r (geometric axis)
    write(lun,510)r0,zmag,psimag,psi_lcfs,b0vac
      curnt=cur*1.0e-6
! curnt is the plasma current in MA
    write(lun,510)curnt,xdum,xdum,xdum,xdum
    write(lun,510)xdum,xdum,xdum,xdum,xdum
      dpsi=umax/(nr-1.)
      do i=1,nr
        psi=(i-1.)*dpsi
        arr(i)=fprof(psi,2)
      end do
!  f(psi) profile, note f=R Bphi in Tm (axis to edge on uniform psi mesh)
    write(lun,510)(arr(i),i=1,nr)
      do i=1,nr
        psi=(i-1.)*dpsi
        arr(i)=press(psi,0)
      end do
!  p(psi) profile, the pressure profile in Nm^{-2} (axis to edge)
    write(lun,510)(arr(i),i=1,nr)
      do i=1,nr
        psi=(i-1.)*dpsi
        arr(i)=fprof(psi,1)
      end do
! ffprime profile   f df/dpsi, psi as defined above (increasing towards edge)
! axis to edge
    write(lun,510)(arr(i),i=1,nr)
      do i=1,nr
        psi=(i-1.)*dpsi
        arr(i)=press(psi,1)
      end do
! pprime profile   dp/dpsi
    write(lun,510)(arr(i),i=1,nr)
! psi mesh The psi-mesh on the points R(i), Z(j) (note that u=umax-psi, so
! this has to be shifted to correcpond to the definition used above
    write(lun,510)(((umax-u(i,j)),i=1,nr),j=1,nz)
! q-profile
      do i=1,nr
        psi=(i-1.)*dpsi
        ij=1
        do 10 j=1,ncon
          if (psi.gt.psiv(j)) goto 10
          ij=j
 10     continue
        if (ij.eq.ncon) ij=ncon-1
        rat=(psi-psiv(ij))/(psiv(ij+1)-psiv(ij))
        if ((rat.gt.1.).or.(rat.lt.0.)) then
          if ((rat.gt.-1.0d-10).and.(rat.lt.0.)) then
            rat=0.
          else if ((rat.lt.1.000000001).and.(rat.gt.1.0d0)) then
            rat=1.
          else
          write(6,*)' rat outside range in gfil, rat=',rat
          write(6,*)' psi=',psi,' psi(ij)=',psiv(ij),' psi(ij+1)=',psiv(ij+1)
          stop
          end if
        end if
        arr(i)=sfac(ij)+rat*(sfac(ij+1)-sfac(ij))
      end do
! q-profile from axis to edge on uniform mesh
    write(lun,510)(arr(i),i=1,nr)
! no. poloidal mesh points of surface
    write(lun,520)npts,idum
! poloidal mesh points on surface
    write(lun,510)(rpts(1,i),zpts(1,i),i=1,npts)
 500  format(a34,i4,a14,2i4)
 510  format(5e16.9)
 520  format(2i5)
      close(lun)
      lun=68
      open(unit=lun,file=runname(1:lrunname)//'.surface', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.surface'
         stop
      endif
     write(lun,*)npts
     write(lun,511)(rpts(1,i),zpts(1,i),i=1,npts)
 511  format(2e16.9)
      close(lun)
      return
    end subroutine gfile
end module gfile_output  
