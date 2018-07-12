      subroutine iaea_plots
!     *********************
!
      use param
      implicit none
      double precision x(ncon),y(ncon)
      double precision press,fprof,bp,dense,tempe
      double precision psi,dx,dy,extapp,extapp2
      integer k,i,icur,jj
      real*4 sx(ncon),sy(ncon)
      real*4 sxmax,sxmin,symax,symin,sdx,sdy,sypos,sxpos
!
!  density
      call filnam('temp.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        y(k)=tempe(psi,0)*1.0d-3
      end do
      dx=0.1
      dy=5.0
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'Te')
      call typecs(' keV')
      call ctrori(0.)
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  density
      call filnam('dens.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        y(k)=dense(psi,0)*1.0d-19
      end do
      dx=0.1
      dy=0.2
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'ne')
      call typecs('X10')
      call supfix
      call typecs('19')
      call normal
      call typecs(' m')
      call supfix
      call typecs('-3')
      call normal
      call ctrori(0.)
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  Bootstrap current density
      call filnam('jbs.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        y(k)=bsj(k)*1.0d-6
      end do
      dx=0.1
      dy=0.1
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'<Jbs.B>/<B')
      call supfix
      call typecs('2')
      call normal
      call typecs('>')
      call supfix
      call typecs('1/2')
      call normal
      call typecs(' MA m')
      call supfix
      call typecs('-2')
      call normal
      call ctrori(0.)
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  Total parallel current density
      call filnam('jtot.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        y(k)=(-fprof(psi,2)*press(psi,1)/sqrt(bsqav(k))- &
             fprof(psi,1)*sqrt(bsqav(k))/(mu0*fprof(psi,2)))*1.0d-6
      end do
      dx=0.1
      dy=0.1
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'<J.B>/<B')
      call supfix
      call typecs('2')
      call normal
      call typecs('>')
      call supfix
      call typecs('1/2')
      call normal
      call typecs(' MA m')
      call supfix
      call typecs('-2')
      call normal
      call ctrori(0.)
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  Driven parallel current density
      call filnam('jdriv.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        icur=1
!  external current profile calculated in extj, returned in extapp
        call extj(epsv(k),psi,extapp,extapp2,icur)
        y(k)=(vloop*extapp+extapp2)*sqrt(bsqav(k))*1.0d-6
      end do
      dx=0.1
      dy=0.1
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'<J')
      call suffix
      call typecs('aux')
      call normal
      call typecs('.B>/<B')
      call supfix
      call typecs('2')
      call normal
      call typecs('> MA m')
      call supfix
      call typecs('-2')
      call normal
      call ctrori(0.)
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  Safety factor
      call filnam('qprof.grd')
      do k=1,ncon
        x(k)=epsv(k)/epsv(1)
        psi=psiv(k)
        y(k)=sfac(k)
      end do
      dx=0.1
      dy=2.0
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,ncon
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
      symin=0.
      if (symax.gt.1.0d9) symax=1.5*y(2)
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,ncon,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call ctrfnt(2)
      call plotcs(sxpos,sypos,'e')
      call ctrfnt(0)
!      call ctrori(90.)
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'q')
      call ctrmag(18)
      call frame
! ******************Done figure***************
!  Total parallel current density
      call filnam('btot.grd')
      jj=0
      do 50 i=1,nr
        if (ixout(i,nsym).le.0) goto 50
        jj=jj+1
        x(jj)=r(i)
        psi=umax-u(i,nsym)
        y(jj)=sqrt((fprof(psi,2)/r(i))**2+(bp(r(i),0.0d0))**2)
 50   continue
      dx=1.0
      dy=1.0
!   Ghost graphics********************
      sdx=sngl(dx)
      sdy=sngl(dy)
      do i=1,jj
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
      symin=0.
      call pspace(0.2,0.9,0.2,0.9)
      call map(sxmin,sxmax,symin,symax)
      call axorig(sxmin,symin)
!      call ctrmag(12)
      call axessi(sdx,sdy)
      call thick(2)
      call ptjoin(sx,sy,1,jj,1)
      call thick(1)
      sypos=symin-0.15*(symax-symin)
      sxpos=sxmin+0.8*(sxmax-sxmin)
      call ctrmag(25)
      call plotcs(sxpos,sypos,'R (m)')
      sypos=symin+0.5*(symax-symin)
      sxpos=sxmin-0.15*(sxmax-sxmin)
      call plotcs(sxpos,sypos,'B (T)')
      call ctrmag(18)
      call frame
! ******************Done figure***************
      return
 end subroutine iaea_plots
 
