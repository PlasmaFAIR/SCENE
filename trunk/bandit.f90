module bandit_output
  implicit none
contains
      subroutine bandit
!     *****************
!
! writes out flux surface data file for use in bandit
!
      use param
      use profiles_mod, only : fprof
      use toms790, only : CSHEP2, CS2VAL
      implicit none
      double precision :: br(nr,nz),bz(nr,nz),bt(nr,nz)
      double precision, dimension(:), allocatable ::  x,y,f
      double precision, dimension(:,:), allocatable :: grads
      double precision :: px,py,pf,fsi,psi,rr
      integer, dimension(:), allocatable :: triang
      integer :: npt,k,i,j,ifail
      integer, allocatable :: LCELL(:,:), LNEXT(:)
      double precision :: XMIN, YMIN, DX, DY, RMAX
      double precision, allocatable :: RW(:), A(:,:)
      integer :: NCC, NRR, NWW
!
!
!  Extrapolation routine to fill out the psi mesh on the R-Z grid
!  from the points calculated from the G-S solver
!
!
!  Load up points inside and just outside plasma
      npt=0
      do i=1,nr
        do j=1,nz
          if (idout(i,j).ne.0) npt=npt+1
        end do
      end do
      allocate( x(npt),y(npt),f(npt),grads(2,npt),triang(7*npt) )
! BR
      k=0
      do  i=1,nr
        do  j=1,nz
          if (idout(i,j).ne.0) then
            k=k+1
            x(k)=rcoord(i)
            y(k)=zcoord(j)
            f(k)=brcoord(i,j)
          end if
        end do
      end do
!      ifail=0
      ! Setup for interpolation
      NCC=MIN(17,npt-1)
      NWW=MIN(30,npt-1)
      NRR=ceiling( sqrt(dble(npt)/dble(3.0)) )
      allocate(LCELL(NRR,NRR))
      allocate(LNEXT(npt))
      allocate(RW(npt))
      allocate(A(9,npt))
      call CSHEP2(npt,x,y,f,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
!      call e01saf(npt,x,y,f,triang,grads,ifail)
      do i=1,nr
        do j=1,nz
          px=r(i)
          py=z(j)
          ifail=1
!          call e01sbf(npt,x,y,f,triang,grads,px,py,pf,ifail)
          pf=CS2VAL(px,py,npt,x,y,f,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
          br(i,j)=pf
        end do
      end do
! BZ
      k=0
      do  i=1,nr
        do  j=1,nz
          if (idout(i,j).ne.0) then
            k=k+1
            x(k)=rcoord(i)
            y(k)=zcoord(j)
            f(k)=bzcoord(i,j)
          end if
        end do
      end do
!      ifail=0
!      call e01saf(npt,x,y,f,triang,grads,ifail)
      call CSHEP2(npt,x,y,f,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      do i=1,nr
        do j=1,nz
          px=r(i)
          py=z(j)
          ifail=1
!          call e01sbf(npt,x,y,f,triang,grads,px,py,pf,ifail)
          pf=CS2VAL(px,py,npt,x,y,f,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
          bz(i,j)=pf
        end do
      end do
      deallocate( x,y,f,grads,triang )
! BT
     do i=1,nr
       do j=1,nz
         rr=r(i)
         psi=umax-u(i,j)
         if (ixout(i,j).le.0) psi=umax
         fsi=fprof(psi,2)
         bt(i,j)=fsi/rr
       end do
     end do
!
!!$!**********************************************************************

      open(unit=51,file=runname(1:lrunname)//'.bandit', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.bandit'
         stop
      endif
      write(51,*)'NR, NZ'
      write(51,*)nr,nz
      write(51,*)'Rmin,Rmax'
      write(51,*)r1,r2
      write(51,*)'R'
      write(51,101)(r(i),i=1,nr)
      write(51,*)'Z'
      write(51,101)(z(i),i=1,nz)
      write(51,*)' BR'
      write(51,101) ((br(i,j),i=1,nr),j=1,nz)
      write(51,*)' BT'
      write(51,101) ((bt(i,j),i=1,nr),j=1,nz)
      write(51,*)' BZ'
      write(51,101) ((bz(i,j),i=1,nr),j=1,nz)
 101  format(8(1pe15.8))
      close(51)
!!$      call pspace(0.05,0.95,0.05,0.95)
!!$      zmin=z(1)
!!$      zmax=z(nz)
!!$      xmin=r(1)
!!$      xmax=r(nr)
!!$      call map(xmin,xmax,zmin,zmax)
!!$      do k=1,ncon
!!$        do i=1,npts
!!$          xp(i)=rban(i,k)
!!$          zp(i)=zban(i,k)
!!$        end do
!!$        call ptjoin(xp,zp,1,npts,0)
!!$      end do
!!$      call frame
      return
      end subroutine bandit
!
!---------------------------------------------------------------------
!
      subroutine curplt
!     *****************
!
! Plots externally applied current within a flux surface as a function of
! r/a
!
      use ext_current_mod, only : extj
      use param
      use profiles_mod, only : dense, tempe
      implicit none
      integer :: icur,k,i,nrm,kg
      double precision :: psi,eps,extapp,extapp2,drm,rat,psiup,de,te,curf,rv
      double precision :: curflx(ncon),rsm(ncon),curden(ncon)


      curflx(ncon)=0.
      do k=ncon,1,-1
        psi=psiv(k)
        eps=epsv(k)
        rsm(k)=(rpts(k,npts/2+1)-rpts(k,1))/(2.*rcen*tokeps)
!        write(6,*)' k=',k,' psi=',psi,' eps=',eps,' rsm=',rsm(k)
!  externally applied current
        if (itot.eq.0) then
          icur=1
          if (neo.lt.0) icur=-1     ! switch off trapping effects
!  external current profile calculated in extj, returned in extapp
      call extj(eps,psi,extapp,extapp2,icur)
          curden(k)=extapp*vloop+extapp2
        else
          write(6,*)' ERROR***this routine is not yet developed for itot=1'
          write(6,*)' It only generates stuff for BANDIT, so switch off call to curplt if not needed'
          stop
        end if
        if (k.lt.ncon) then
          curflx(k)=curflx(k+1)+pi*(psiv(k)-psiv(k+1))*(           &
                   sfac(k+1)*curden(k+1)+sfac(k)*curden(k))
        end if
 !       write(6,*)' r=',rsm(k),' current=',curflx(k)
      end do
!  Put on an equal-spaced r-grid
      nrm=41
      drm=rsm(1)/(nrm-1)
      open(76,file='jntprof.dat')
      do i=1,nrm
        rv=(i-1)*drm
        kg=1
        do 10 k=1,ncon
          if (rv.gt.rsm(k)) goto 10
          kg=k
 10     continue
        if (kg.eq.ncon) kg=kg-1
        rat=(rv-rsm(kg))/(rsm(kg+1)-rsm(kg))
        curf=(curflx(kg)+rat*(curflx(kg+1)-curflx(kg)))*1.0d-6
        psi=psiv(kg)
        psiup=psiv(kg+1)
        de=dense(psi,0)+rat*(dense(psiup,0)-dense(psi,0))
        te=tempe(psi,0)+rat*(tempe(psiup,0)-tempe(psi,0))
        write(76,20)rv,curf,de,te
      end do
 20   format('r/a=',e12.4,' current=',e12.4,' MA, n=',e12.4,' m**-3, Te=',e12.4,' keV')
      close(76)
end subroutine curplt
end module bandit_output
