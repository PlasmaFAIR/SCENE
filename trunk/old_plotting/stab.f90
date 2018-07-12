      subroutine stab
!     ***************
!
!   Evaluates and plots a number of useful stability parameters
!
      use param
      implicit none
      integer k,ihalf
      double precision tempe,dense,tempi,densi,fprof
      double precision ne,ned,ni,nid,te,ti,ted,tid
      double precision omstre(ncon),omstri(ncon),omegaa(ncon)
      double precision psi,fsi
      real*4 omin,omax,xmin,xmax,zero
      real*4 yv,xv1,xv2,yv1,yv2
      real*4 xp(ncon),yp1(ncon),yp2(ncon)
!
      zero=0.
      ihalf=npts/2+1
!  Calculate omega*i and omega*e
      do k=1,ncon
        psi=psiv(k)
        ne=dense(psi,0)
        ned=dense(psi,1)
        ni=densi(psi,1,0)
        nid=densi(psi,1,1)
        te=tempe(psi,0)
        ted=tempe(psi,1)
        ti=tempi(psi,1,0)
        tid=tempi(psi,1,1)
        omstre(k)=-te*(ned/ne+ted/te)
        omstri(k)=-(ti/zm)*(nid/ni+tid/ti)
        fsi=fprof(psi,2)
        if (k.lt.ncon) then
          omegaa(k)=(fsi/rpts(k,1)+fsi/rpts(k,ihalf))/      &
                    ((rpts(k,1)+rpts(k,ihalf))*sqrt(mu0*ni*zmai*mp))
        else
          omegaa(k)=fsi/(r0**2*sqrt(mu0*ni*zmai*mp))
        end if
        if (k.eq.1) then
          omin=sngl(omstre(k))
          omax=sngl(omstre(k))
          if (dble(omin).gt.omstri(k)) omin=sngl(omstri(k))
          if (dble(omax).lt.omstri(k)) omax=sngl(omstri(k))
        else
          if (dble(omin).gt.omstre(k)) omin=sngl(omstre(k))
          if (dble(omax).lt.omstre(k)) omax=sngl(omstre(k))
          if (dble(omin).gt.omstri(k)) omin=sngl(omstri(k))
          if (dble(omax).lt.omstri(k)) omax=sngl(omstri(k))
        end if
        xp(k)=sngl(epsv(k))
        yp1(k)=sngl(omstre(k))
        yp2(k)=sngl(omstri(k))
      end do
      omin=0.
      xmin=xp(ncon)
      xmax=xp(1)
      xv1=0.11
      xv2=0.45
      yv1=0.47
      yv2=0.81
      call pspace(xv1,xv2,yv1,yv2)
      call border
      call map(xmin,xmax,omin,omax)
      call axorig(xmin,omin)
      call ctrmag(12)
      call axessi(zero,zero)
      call ptjoin(xp,yp1,1,ncon,1)
      call lincol(2)
      call broken(10,5,10,5)
      call ptjoin(xp,yp2,1,ncon,1)
      yv=omin+0.87*(omax-omin)
      xv1=xmin+0.02*(xmax-xmin)
      xv2=xmin+0.35*(xmax-xmin)
      call ctrmag(20)
      call positn(xv1,yv)
      call join(xv2,yv)
      call full
      call lincol(0)
      call ctrfnt(2)
      call plotcs(xv2,yv,' w')
      call ctrfnt(0)
      call suffix
      call typecs('*pi')
      call normal
      yv=omin+0.95*(omax-omin)
      call positn(xv1,yv)
      call join(xv2,yv)
      call lincol(0)
      call ctrfnt(2)
      call plotcs(xv2,yv,' w')
      call ctrfnt(0)
      call suffix
      call typecs('*pe')
      call normal
      do k=1,ncon
        if (k.eq.1) then
          omin=sngl(omegaa(k))
          omax=sngl(omegaa(k))
        else
          if (dble(omin).gt.omegaa(k)) omin=sngl(omegaa(k))
          if (dble(omax).lt.omegaa(k)) omax=sngl(omegaa(k))
        end if
        xp(k)=sngl(epsv(k))
        yp1(k)=sngl(omegaa(k))
      end do
      omin=0.
      xv1=0.57
      xv2=0.91
      yv1=0.47
      yv2=0.81
      call pspace(xv1,xv2,yv1,yv2)
      call border
      call map(xmin,xmax,omin,omax)
      call axorig(xmin,omin)
      call ctrmag(12)
      call axessi(zero,zero)
      call ptjoin(xp,yp1,1,ncon,1)
      xv1=xmin+0.02*(xmax-xmin)
      xv2=xmin+0.35*(xmax-xmin)
      yv=omin+0.95*(omax-omin)
      call ctrmag(20)
      call positn(xv1,yv)
      call join(xv2,yv)
      call lincol(0)
      call ctrfnt(2)
      call plotcs(xv2,yv,' w')
      call ctrfnt(0)
      call suffix
      call typecs('A')
      call normal
      call frame
      call filnam ('diaplot.grd')
      do k=1,ncon
        yp1(k)=6*omstri(k)/(sqrt(2.4828942e-4)*omegaa(1))
        if (k.eq.1) then
          omax=yp1(1)
        else
          if (omax.lt.yp1(k)) omax=yp1(k)
        end if
      end do
      omin=0.
      xv1=0.11
      xv2=0.45
      yv1=0.47
      yv2=0.81
      call pspace(xv1,xv2,yv1,yv2)
      call border
      call map(xmin,xmax,omin,omax)
      call axorig(xmin,omin)
      call ctrmag(12)
      call axessi(zero,zero)
      call ptjoin(xp,yp1,1,ncon,1)
      call frame
      return
  end subroutine stab
!
!************************************************************************
!
      subroutine epsplot
!     ***************
!
!   Evaluates and plots a number of useful stability parameters
!
      use param
      implicit none
      integer k,ihalf,nc,l,lmax,i,itst
      double precision tempe,dense,tempi,densi,fprof,press
      double precision ne,ned,ni,nid,te,ti,ted,tid
      double precision omstre(ncon),omstri(ncon),omegaa(ncon)
      double precision psi,fsi
      real*4 nmin,nmax,xmin,xmax,zero,tmin,tmax
      real*4 yv,xv1,xv2,yv1,yv2
      real*4 xp(ncon),yp1(ncon),yp2(ncon)
      real*4 zr(nr),zz(nz),spsi,xc(npts),yc(npts)
      real*4 zlo,zup,rlo,rup,elong
!
      call rotate(90.)
      call filnam('flxcon.grd')
      elong=sngl(2.*zs/(r2-r1))
      do i=1,nr
        zr(i)=sngl(r(i))
      end do
      do i=1,nz
        zz(i)=sngl(z(i))
      end do
      do 20 nc=ncon,1,-1
        psi=psiv(nc)
        spsi=sngl(psi)
        itst=10*((nc-1)/10)-(nc-1)
        if (itst.eq.0) then
          zlo=0.16
          zup=0.75
          rup=0.81
          rlo=rup-(zup-zlo)/elong
          call pspace(rlo,rup,zlo,zup)
          zlo=0.
          call map(zr(1),zr(nr),zz(1),zz(nz))
          call border
          call axorig(zr(1),zz(1))
          call ctrmag(16)
          call axessi(1.0,1.0)
          call ctrmag(12)
          lmax=npts
          do  l=1,lmax
            xc(l)=sngl(rpts(nc,l))
            yc(l)=sngl(zpts(nc,l))
          end do
          if (nc.eq.1) then
            call thick(3)
            call lincol(2)
          end if
          call ptjoin(xc,yc,1,lmax,-1)
 15       call thick(1)
          call lincol(0)
	end if
 20   continue
      
!
!-------------------------------------------------------------------
!
      call filnam('pplt.grd')
      call thick(0)
      zero=0.
      do k=1,ncon
        psi=psiv(k)
        yp1(k)=sngl(press(psi,0)/press(0.0d0,0))
        yp2(k)=sngl(-fprof(psi,2)*press(psi,1)/sqrt(bsqav(k))-fprof(psi,1)  &
                   *sqrt(bsqav(k))/(mu0*fprof(psi,2)))
        xp(k)=sngl(psiv(k)/psiv(1))
        xp(k)=sngl(epsv(k)/tokeps)
        write(6,*)' xp=',xp(k),' yp2=',yp2(k)
        if (k.eq.1) then
          nmax=yp1(k)
          tmax=yp2(k)
        else
          if (nmax.lt.yp1(k)) nmax=yp1(k)
          if (tmax.lt.yp2(k)) tmax=yp2(k)
        end if
      end do
      tmax=tmax/yp2(ncon)
      do k=1,ncon
        yp2(k)=yp2(k)/yp2(ncon)
      end do
      nmin=0.
      tmin=0.
      xmin=xp(ncon)
      xmax=xp(1)
      call pspace(0.2,0.6,0.2,0.6)
      call border
      call thick(1)
      call map(xmin,xmax,nmin,nmax)
      call axorig(xmin,nmin)
      call ctrmag(14)
      call axessi(0.2,0.2)
      call lincol(2)
      call ptjoin(xp,yp1,1,ncon,1)
      call lincol(0)
      yv=nmin+0.7*(nmax-nmin)
      xv1=xmin-0.25*(xmax-xmin)
      call ctrmag(20)
      call plotcs(xv1,yv,' p')
      yv=nmin-0.17*(nmax-nmin)
      xv1=xmin+0.5*(xmax-xmin)
      call ctrmag(20)
!      call ctrfnt(2)
      call plotcs(xv1,yv,' r/a')
 !     call ctrfnt(0)
 !     call suffix
 !     call typecs('a')
 !     call normal
      call frame
!-----------------------------------------------
      call filnam('jplot.grd')
      call rotate(90.)
      call pspace(0.2,0.6,0.2,0.6)
      call border
      call map(xmin,xmax,tmin,tmax)
      call axorig(xmin,tmin)
      call ctrmag(14)
      call axessi(0.2,0.2)
      call lincol(2)
      call thick(1)
      call ptjoin(xp,yp2,1,ncon,1)
      call thick(0)
      call lincol(0)
      yv=nmin+0.5*(tmax-tmin)
      xv1=xmin-0.2*(xmax-xmin)
      call ctrmag(20)
      call ctrori(90.)
      call plotcs(xv1,yv,' <J.B>/<B')
      call supfix
      call typecs('2')
      call normal
      call typecs('>')
      call supfix
      call typecs('1/2')
      call normal
      call ctrori(0.)
      yv=nmin-0.17*(nmax-nmin)
      xv1=xmin+0.5*(xmax-xmin)
      call ctrmag(20)
!      call ctrfnt(2)
      call plotcs(xv1,yv,' r/a')
!      call ctrfnt(0)
!      call suffix
!      call typecs('a')
      call normal
      call frame
!-------------------------------------------------------------------
!
      call filnam('nplt.grd')
      call thick(0)
      zero=0.
      do k=1,ncon
        psi=psiv(k)
        yp1(k)=sngl(dense(psi,0)/dense(0.0d0,0))
        yp2(k)=sngl(tempe(psi,0)/tempe(0.0d0,0))
        xp(k)=sngl(psiv(k)/psiv(1))
        xp(k)=sngl(epsv(k)/tokeps)
        if (k.eq.1) then
          nmax=yp1(k)
          tmax=yp2(k)
        else
          if (nmax.lt.yp1(k)) nmax=yp1(k)
          if (tmax.lt.yp2(k)) tmax=yp2(k)
        end if
      end do
      nmin=0.
      tmin=0.
      xmin=xp(ncon)
      xmax=xp(1)
      call pspace(0.2,0.6,0.2,0.6)
      call border
      call map(xmin,xmax,nmin,nmax)
      call axorig(xmin,nmin)
      call ctrmag(14)
      call axessi(0.2,0.2)
      call lincol(2)
      call thick(1)
      call ptjoin(xp,yp1,1,ncon,1)
      call thick(0)
      call lincol(0)
      yv=nmin+0.7*(nmax-nmin)
      xv1=xmin-0.25*(xmax-xmin)
      call ctrmag(20)
      call plotcs(xv1,yv,' n')
      yv=nmin-0.17*(nmax-nmin)
      xv1=xmin+0.5*(xmax-xmin)
      call ctrmag(20)
!      call ctrfnt(2)
      call plotcs(xv1,yv,' r/a')
!      call ctrfnt(0)
!      call suffix
!      call typecs('a')
      call normal
      call frame
!-----------------------Temperature profile------------------
      call pspace(0.2,0.6,0.2,0.6)
      call border
      call map(xmin,xmax,tmin,tmax)
      call axorig(xmin,tmin)
      call ctrmag(14)
      call axessi(0.2,0.2)
      call lincol(2)
      call thick(1)
      call ptjoin(xp,yp2,1,ncon,1)
      call thick(0)
      call lincol(0)
      yv=nmin+0.7*(nmax-nmin)
      xv1=xmin-0.25*(xmax-xmin)
      call ctrmag(20)
      call plotcs(xv1,yv,' T')
      yv=nmin-0.17*(nmax-nmin)
      xv1=xmin+0.5*(xmax-xmin)
      call ctrmag(20)
!      call ctrfnt(2)
      call plotcs(xv1,yv,' r/a')
!      call ctrfnt(0)
!      call suffix
!      call typecs('a')
      call normal
      call frame
!-----------------------------------------------
      call rotate(0.)
      return
      end subroutine epsplot





