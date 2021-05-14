module flxorb_output
  implicit none
contains
      subroutine flxorb
!     -----------------
!  writes data for Ken's orbit calculation code

!
      use equilibrium, only : bp
      use param
      use profiles_mod, only : fprof
      use toms790, only : CSHEP2, CS2VAL
      implicit none
      integer i,j,k
      double precision rr,fsi,psi,btor
      double precision rwall,zwall,rinner,rtst,ztst
      double precision, dimension(:), allocatable:: rk,zk
      double precision, dimension(:,:), allocatable:: bpk,uk

!
      double precision, dimension(:), allocatable::  x,y,f
      double precision, dimension(:,:), allocatable:: grads
      double precision px,py,pf
      integer, dimension(:), allocatable:: triang
      integer npt_nag,ifail,nrk,nzk
!
      integer, allocatable:: LCELL(:,:), LNEXT(:)
      double precision:: XMIN, YMIN, DX, DY, RMAX
      double precision, allocatable:: RW(:), A(:,:)
      integer:: NCC, NRR, NWW
!
      open(unit=51,file=runname(1:lrunname)//'.flxorb', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.flxorb'
         stop
      endif
!  Need to add some grid points for power plant case to ensure
!  it includes the wall
!  major radius of outboard wall
      rwall=6.3
      rwall=r(nr)
!  vertical height of upper wall
      zwall=8.0
      zwall=z(nz)
      rinner=r(1)
      nrk=int(1+(rwall-rinner)/step)
      nzk=int(zwall/step)
      rtst=rinner+(nrk-1)*step
      if (rtst.lt.rwall) nrk=nrk+1
      ztst=(nzk-1)*step
      if (ztst.lt.zwall) nzk=nzk+1
      ztst=(nzk-1)*step
      nzk=2*nzk-1
      allocate(rk(nrk),zk(nzk),bpk(nrk,nzk),uk(nrk,nzk))
      do i=1,nrk
        rk(i)=rinner+(i-1)*step
      end do
      if (rk(nrk).lt.rwall) then
        write(6,*)' error in setting nrk for Kens code'
        write(6,*)' rk=',rk(nrk),' nrk=',nrk,' rwall=',rwall
        stop
      end if
      do j=1,nzk
        zk(j)=-zwall+(j-1)*step
      end do
      if (zk(nzk).lt.zwall) then
        write(6,*)' error in setting nzk for Kens code'
        write(6,*)' zk=',zk(nzk),' nzk=',nzk,' zwall=',zwall
        stop
      end if
      write(51,*)'NR       ','NZ       '
      write(51,*)nrk,nzk
      write(51,*)'Rin      ','Rout     '
      write(51,*)rk(1),rk(nrk)
      write(51,*)'R        '
      do  i=1,nrk
        write(51,*)rk(i)
      end do
      write(51,*)'Z        '
      do 410 i=1,nzk
        write(51,*)zk(i)
 410  continue
      write(51,*)'Bpol     '
!-------------------------------------------------------------------
!  Evaluate bpol by interpolating/extrapolating known values (Nag routine)
!
!  Extrapolation routine to fill out the psi mesh on the R-Z grid
!  from the points calculated from the G-S solver
!
!
!  Load up points inside and just outside plasma
      npt_nag=0
      do i=1,nr
        do j=1,nz
          if (ixout(i,j).ne.0) npt_nag=npt_nag+1
        end do
      end do
      allocate( x(npt_nag),y(npt_nag),f(npt_nag),grads(2,npt_nag),    &
               triang(7*npt_nag) )
      k=0
      do  i=1,nr
        do  j=1,nz
          if (ixout(i,j).ne.0) then
            k=k+1
            x(k)=r(i)
            y(k)=z(j)
            f(k)=bp(r(i),z(j))
          end if
        end do
      end do
!      ifail=0
      ! Setup for interpolation
      NCC=MIN(17,npt_nag-1)
      NWW=MIN(30,npt_nag-1)
      NRR=ceiling( sqrt(dble(npt_nag)/dble(3.0)) )
      allocate(LCELL(NRR,NRR))
      allocate(LNEXT(npt_nag))
      allocate(RW(npt_nag))
      allocate(A(9,npt_nag))
!      call e01saf(npt_nag,x,y,f,triang,grads,ifail)
      call CSHEP2(npt_nag,x,y,f,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
      do i=1,nrk
        do j=1,nzk
          px=rk(i)
          py=zk(j)
          ifail=1
!          call e01sbf(npt_nag,x,y,f,triang,grads,px,py,pf,ifail)
          pf=CS2VAL(px,py,npt_nag,x,y,f,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
          write(51,*)pf
        end do
      end do
      deallocate( x,y,f,grads,triang )
!
!---------------------------------------------------------------------
! repeat, loading up uk now
      npt_nag=0
      do i=1,nr
        do j=1,nz
          if (ixout(i,j).ne.0) npt_nag=npt_nag+1
        end do
      end do
      allocate( x(npt_nag),y(npt_nag),f(npt_nag),grads(2,npt_nag),    &
               triang(7*npt_nag) )
      k=0
      do  i=1,nr
        do  j=1,nz
          if (ixout(i,j).ne.0) then
            k=k+1
            x(k)=r(i)
            y(k)=z(j)
            f(k)=u(i,j)
          end if
        end do
      end do
      ifail=0
      ! Setup for interpolation
!      call e01saf(npt_nag,x,y,f,triang,grads,ifail)
      call CSHEP2(npt_nag,x,y,f,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
      do i=1,nrk
        do j=1,nzk
          px=rk(i)
          py=zk(j)
          ifail=1
!          call e01sbf(npt_nag,x,y,f,triang,grads,px,py,pf,ifail)
          pf=CS2VAL(px,py,npt_nag,x,y,f,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
          uk(i,j)=pf
        end do
      end do
      deallocate( x,y,f,grads,triang )



!!$!-----------------------------------------------------------------------
!!$      do 430 i=1,nr
!!$        do 420 j=1,nz
!!$          if (ixout(i,j).eq.0) then
!!$            bth=0.
!!$          else
!!$            rr=r(i)
!!$            zz=z(j)
!!$            bth=bp(rr,zz)
!!$          end if
!!$          write(51,*)bth
!!$ 420    continue
!!$ 430  continue
      write(51,*)'Bphi     '
      do 450 i=1,nrk
        rr=rk(i)
        do 440 j=1,nzk
          if (uk(i,j).le.0) then
            psi=umax
            fsi=fprof(psi,2)
          else
            psi=umax-uk(i,j)
            fsi=fprof(psi,2)
          end if
          btor=fsi/rr
!          write(6,*)'i=',i,' j=',j,' psi=',psi,' fsi=',fsi,' bphi=',bphi
          write(51,*)btor
!          write(6,*)' done write'
 440    continue
 450  continue
      write(51,*)'psi*2*pi '
      write(6,*)'psi*2*pi '
      do i=1,nrk
        do j=1,nzk
!!$          if (ixout(i,j).eq.0) then
!!$            psi=0.
!!$          else
!!$            psi=u(i,j)
!!$          end if
          psi=uk(i,j)
          write(51,*)psi*2.*pi
        end do
      end do
      write(51,*)' npol=',npts
      write(51,*)' nflux=',ncon-1
      do i=ncon-1,1,-1
        write(51,*)' psi*2*pi=',2.*psiv(i)*pi
        write(51,*)' R values'
        do j=1,npts
          write(51,*)rpts(i,j)
        end do
        write(51,*)' Z values'
        do j=1,npts
          write(51,*)zpts(i,j)
        end do
!        write(51,*)' Bp values'
!        do j=1,npts
!          write(51,*)bppts(i,j)
!        end do
      end do
      return
   end subroutine flxorb
end module flxorb_output
