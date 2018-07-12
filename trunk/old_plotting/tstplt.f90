      subroutine tstplt(ncon,x,y,dx,dy)
!     ************************
!
      implicit none
      integer ncon,i
      double precision x(ncon),y(ncon)
      double precision dx,dy
      real*4 sx(ncon),sy(ncon)
      real*4 sxmax,sxmin,symax,symin,sdx,sdy
!
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon-1
        call ctrmag(18)
        sx(i)=sngl(x(i))
        sy(i)=sngl(y(i))
        if (i.eq.1) then
          sxmax=sx(1)
          sxmin=sx(1)
          symax=sy(1)
          symin=sy(1)
        else
          if (sxmax.lt.sx(i)) sxmax=sx(i)
          if (sxmin.gt.sx(i)) sxmin=sx(i)
          if (symax.lt.sy(i)) symax=sy(i)
          if (symin.gt.sy(i)) symin=sy(i)
        end if
      end do
      if (symax.gt.1.0d9) symax=1.5*y(2)
!      sxmin=-3.14
!      sxmax=3.145
!      symin=0.
!      symax=2.
      write(6,*)' sxmax=',sxmax,' sxmin=',sxmin
      write(6,*)' symax=',symax,' symin=',symin
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,2,ncon-1,1)
      call thick(1)
      call frame
  end subroutine tstplt
!
!*********************************************************
!
      subroutine tstplt2(ncon,x,y,y2)
!     ************************
!
      implicit none
      integer ncon,i,imin,imax
      double precision x(ncon),y(ncon),y2(ncon)
      real*4 sx(ncon),sy(ncon),sy2(ncon)
      real*4 sxmax,sxmin,symax,symin
!
      imin=1
      imax=ncon
      do i=imin,imax
        sx(i)=sngl(x(i))
        sy(i)=sngl(y(i))
        sy2(i)=sngl(y2(i))
        if (i.eq.imin) then
          sxmax=sx(i)
          sxmin=sx(i)
          symax=sy(i)
          symin=sy(i)
          if (symax.lt.sy2(i)) symax=sy2(i)
          if (symin.gt.sy2(i)) symin=sy2(1)
        else
          if (sxmax.lt.sx(i)) sxmax=sx(i)
          if (sxmin.gt.sx(i)) sxmin=sx(i)
          if (symax.lt.sy(i)) symax=sy(i)
          if (symin.gt.sy(i)) symin=sy(i)
          if (symax.lt.sy2(i)) symax=sy2(i)
          if (symin.gt.sy2(i)) symin=sy2(i)
        end if
      end do
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
      call ctrmag(12)
      call axessi(0.,0.)
      call ptjoin(sx,sy,imin,imax,1)
      call broken(10,5,10,5)
      call lincol(3)
      call ptjoin(sx,sy2,imin,imax,1)
      call lincol(0)
      call full
      call frame
  end subroutine tstplt2
!
!*********************************************************
!
      subroutine tstplt3(ncon,ncon2,x,x2,y,y2)
!     ************************
!
      implicit none
      integer ncon,i,imin,imax,ncon2
      double precision x(ncon),y(ncon),y2(ncon2),x2(ncon2)
      real*4 sx(ncon),sy(ncon),sy2(ncon2),sx2(ncon2)
      real*4 sxmax,sxmin,symax,symin
!
      imin=1
      imax=ncon
      do i=imin,imax
        sx(i)=sngl(x(i))
        sy(i)=sngl(y(i))
        if (i.eq.imin) then
          sxmax=sx(i)
          sxmin=sx(i)
          symax=sy(i)
          symin=sy(i)
          if (symax.lt.sy2(i)) symax=sy2(i)
        else
          if (sxmax.lt.sx(i)) sxmax=sx(i)
          if (sxmin.gt.sx(i)) sxmin=sx(i)
          if (symax.lt.sy(i)) symax=sy(i)
          if (symin.gt.sy(i)) symin=sy(i)
        end if
      end do
      imax=ncon2
      do i=imin,imax
        sx2(i)=sngl(x2(i))
        sy2(i)=sngl(y2(i))
        if (sxmax.lt.sx2(i)) sxmax=sx2(i)
        if (sxmin.gt.sx2(i)) sxmin=sx2(i)
        if (symax.lt.sy2(i)) symax=sy2(i)
        if (symin.gt.sy2(i)) symin=sy2(i)
      end do
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
      call ctrmag(12)
      call axessi(0.,0.)
      call ptjoin(sx,sy,1,ncon,1)
      call broken(10,5,10,5)
      call lincol(3)
      call ptjoin(sx2,sy2,1,ncon2,1)
        call lincol(0)
      call full
      !call frame
  end subroutine tstplt3
!
!***************************************************************
!
      subroutine tstplt5(ncon,x,y0,y1,y2,y3,y4,y5)
!     ************************
!
      implicit none
      integer ncon,i,imin,imax
      double precision x(ncon),y0(ncon),y1(ncon),y2(ncon),y3(ncon),y4(ncon),y5(ncon)
      real*4 sx(ncon),sy0(ncon), sy1(ncon),sy2(ncon),sy3(ncon),sy4(ncon)
      real*4 sxmax,sxmin,symax,symin
!
      imin=1
      imax=ncon
      do i=imin,imax
        sx(i)=sngl(x(i))
        sy0(i)=sngl(y0(i))
        sy1(i)=sngl(y1(i))
        sy2(i)=sngl(y2(i))
        sy3(i)=sngl(y3(i))
        sy4(i)=sngl(y4(i))
        if (i.eq.imin) then
          sxmax=sx(i)
          sxmin=sx(i)
          symax=sy0(i)
          symin=sy0(i)
          if (symax.lt.sy1(i)) symax=sy1(i)
          if (symin.gt.sy1(i)) symin=sy1(1)
          if (symax.lt.sy2(i)) symax=sy2(i)
          if (symin.gt.sy2(i)) symin=sy2(1)
          if (symax.lt.sy3(i)) symax=sy3(i)
          if (symin.gt.sy3(i)) symin=sy3(1)
          if (symax.lt.sy4(i)) symax=sy4(i)
          if (symin.gt.sy4(i)) symin=sy4(1)
        else
          if (sxmax.lt.sx(i)) sxmax=sx(i)
          if (sxmin.gt.sx(i)) sxmin=sx(i)
          if (symax.lt.sy0(i)) symax=sy0(i)
          if (symin.gt.sy0(i)) symin=sy0(i)
          if (symax.lt.sy1(i)) symax=sy1(i)
          if (symin.gt.sy1(i)) symin=sy1(i)
          if (symax.lt.sy2(i)) symax=sy2(i)
          if (symin.gt.sy2(i)) symin=sy2(i)
          if (symax.lt.sy3(i)) symax=sy3(i)
          if (symin.gt.sy3(i)) symin=sy3(i)
          if (symax.lt.sy4(i)) symax=sy4(i)
          if (symin.gt.sy4(i)) symin=sy4(i)
        end if
      end do
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
      call ctrmag(12)
      call axessi(0.,0.)
      call broken(10,5,10,5)
      call ptjoin(sx,sy0,imin,imax,1)
      call broken(5,10,5,10)
      call lincol(3)
      call ptjoin(sx,sy1,imin,imax,1)
      call broken(20,5,20,5)
      call lincol(2)
      call ptjoin(sx,sy2,imin,imax,1)
      call broken(30,5,30,5)
      call lincol(4)
      call ptjoin(sx,sy3,imin,imax,1)
      call broken(30,10,10,10)
      call lincol(5)
      call ptjoin(sx,sy4,imin,imax,1)
      call lincol(0)
      call full
      do i=imin,imax
        sy0(i)=sngl(y5(i))
        if (i.eq.imin) then
          symax=sy0(i)
          symin=sy0(i)
        else
          if (symax.lt.sy0(i)) symax=sy0(i)
          if (symin.gt.sy0(i)) symin=sy0(i)
        end if
      end do
      call map(sxmin,sxmax,symin,symax)
      call ptjoin(sx,sy0,imin,imax,1)
      call frame
  end subroutine tstplt5
