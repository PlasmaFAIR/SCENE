      subroutine flxsur
!     *****************
!
!----------------------------------------------------------------------
!
!  Converts the equilibrium defined on the R-Z mesh to flux surface
!  coordinates. The poloidal flux is zero on axis, and rises to umax
!  at the edge
!
!------------------------------------------------------------------------
!
!
     use param
     implicit none
     integer i
     double precision deps,psi,emin
!
!  First allocate arrays
      if (icont.gt.-3) then
        allocate( psiv(ncon), epsv(ncon), circumf(ncon) )
        allocate( rpts(ncon,npts),zpts(ncon,npts),bppts(ncon,npts) )
      end if
!  Store flux surfaces equal spaces in eps
      emin=2.*dr/rcen
      deps=(tokeps-emin)/(ncon-2)
      do i=1,ncon-1
        epsv(i)=tokeps-(i-1)*deps
        call psical(epsv(i),psi)
        psiv(i)=psi
      end do
      epsv(ncon)=0.
      psiv(ncon)=0.
! Find R and Z coordinates of flux surfaces;
      do i=1,ncon-1
        psi=psiv(i)
        call flxcon(psi,i)
      end do
  end subroutine flxsur
!
!---------------------------------------------------------------------------
!
!*****************************************************************
!
      subroutine argen(uval,epsi,k)
!     *****************************
!
!  Calculates the inverse aspect ratio of the flux surface,
!  uval and returns it in epsi. (Note U is the defn of poloidal flux
!  used in the Grad-Shafranov solver...ie max on axis, zero at edge)
!
      use param
      implicit none
      double precision uval,epsi
      integer k
      integer jpass,l,i,ixmid(nr+1)
      double precision un,ul,prod
      double precision zer(2)
      double precision umid(nr+1),rmid(nr+1)
!
      if (uval.eq.umax) then
        epsi=0.
        return
      end if
      if (uval.eq.0.0d0) then
        epsi=tokeps
        return
      end if
      if ((uval.ge.umax).or.(uval.le.0.)) then
        write(6,*)' ERROR***uval is out of range in argen...uval='  &
                  ,uval,' umax=',umax
        stop
      end if
      jpass=0
      do i=1,nr+1
        if (r(i-jpass).lt.r0) then
          umid(i)=u(i,nsym)
          rmid(i)=r(i)
          ixmid(i)=ixout(i,nsym)
        else
          if (jpass.eq.0) then
            umid(i)=umax
            rmid(i)=r0
            ixmid(i)=1
            jpass=1
          else
            umid(i)=u(i-1,nsym)
            rmid(i)=r(i-1)
            ixmid(i)=ixout(i-1,nsym)
          end if
        end if
      end do
      jpass=0
 5    un=-uval
      l=0
      do 10 i=1,nr+1
!  find the two zeros of psi-u....
        ul=un
        if (ixmid(i).eq.0) goto 10
        un=umid(i)-uval
        prod=un*ul
        if (prod.lt.0.) then
          l=l+1
          if (l.lt.3) zer(l)=rmid(i-1)+ul*dr/(ul-un)
        else if (prod.eq.0.) then
          l=l+1
          if (l.lt.3) zer(l)=r(i)
          un=-ul
        end if
 10   continue
!  check that there are two zeros...
      if (l.eq.2) then
!  found two zeros
        if (k.eq.0) then
          epsi=abs(zer(1)-zer(2))/(zer(1)+zer(2))
        else
          epsi=abs(zer(1)-zer(2))/(2.*zer(2))
        end if
      else
        write(nw,*)' error***wrong number of zeros in argen'
        write(nw,*)' number=',l,' U value=',uval,' umax=',umax
        write(nw,*)' k=',k,' zer(1)=',zer(1)
        do 30 i=1,nr+1
          write(nw,*)' ir=',i,' rmid-',rmid(i),' umid=',umid(i)
 30     continue
        stop
      end if
    end subroutine argen
!
!*****************************************************************
!
      subroutine flxcon(psi,k)
!     ************************
!
!  calculates the path of the flux surface labelled by k.
!  these are loaded into flxr,flxz and icon labels the number of points
!
      use param
      implicit none
!
      integer is,ic,izz,ilor,iupr,ir,icon,igot,jpt,i,kc,jc,k
      integer irem
      double precision psi,dps,un,uval,ul,rr,zz,prod,rgot,zgot,bp,curint
      double precision flxr(nr*nz),flxz(nr*nz)
      double precision rl(nr*nz),zl(nr*nz),leng(nr*nz),work(nr*nz)
      double precision dl,leq(npts)
      double precision x1,x2,x3,d1,d2,ydd0
!
      un=0.
      uval=umax-psi
      is=1
      ic=0
      izz=nsym-1
 5    izz=izz+is
      if (is.eq.-1) then
        if (izz.eq.(nsym-1)) goto 30
      end if
      zz=z(izz)
      igot=0
      if (is.eq.1) then
        ilor=1
        iupr=nr
      else
        ilor=nr
        iupr=1
      end if
      jpt=0
      do 10 i=ilor,iupr,is
        if (igot.eq.1) goto 10
!  find the two zeros of psi-u....
        ul=un
        rr=r(i)
        if (ixout(i,izz).eq.0) then
          un=u(i,izz) -uval
          goto 10
        end if
!        if (jpt.eq.0) then
!          jpt=1
!          un=u(i,izz)-uval
!          goto 10
!        end if
        un=u(i,izz)-uval
        prod=un*ul
        if (prod.lt.0.) then
          rgot=r(i-is)+is*abs(ul*dr/(un-ul))
          ic=ic+1
          flxr(ic)=rgot
          flxz(ic)=zz
          ir=i-1
          igot=1
        else if (prod.eq.0.) then
          rgot=r(i)
          ic=ic+1
          flxr(ic)=rgot
          flxz(ic)=zz
          ir=i
          igot=1
        end if
 10   continue
      if (igot.eq.0) then
!  need to scan across top/bottom of flux surface
  15    ir=ir+is
        rr=r(ir)
        igot=0
        jpt=0
        do 20 i=nz,1,-1
          if (igot.eq.1) goto 20
!  find the two zeros of psi-u....
          ul=un
          zz=z(i)
          if (ixout(ir,i).eq.0) then
            un=u(ir,i)-uval
            goto 20
          end if
!          if (jpt.eq.0) then
!            jpt=1
!            un=u(i,izz)-uval
!            goto 20
!          end if
          un=u(ir,i)-uval
          prod=un*ul
          if (prod.lt.0.) then
            zgot=z(i+1)-abs(ul*dr/(un-ul))
            ic=ic+1
            flxr(ic)=rr
            flxz(ic)=zgot
            izz=i+1
            igot=1
          else if (prod.eq.0.) then
            zgot=z(i)
            ic=ic+1
            flxr(ic)=rr
            flxz(ic)=zgot
            izz=i
            igot=1
          end if
 20     continue
        if (igot.ne.0) goto 15
      else
        goto 5
      end if
      if ((flxz(ic)/flxz(ic-1)).lt.1.0d-10) then
!  last point too close to zero, so throw away, and go to r-scan
        ic=ic-1
        ir=ir-1
        zgot=flxz(ic)
        do 25 i=nz,1,-1
          if (zgot.gt.z(i)) goto 25
          izz=i
 25     continue
      end if
      is=-is
      goto 5
 30   continue
!  use up/down symmetry to create lower half of flux surface
      kc=ic
      jc=ic
      if (ic.gt.nz*nr/2) then
        write(nw,*)'error***mesh is too fine'
        write(nw,*)' no points for U=',uval,' is ',ic
        write(nw,*)'max allowed is nr*nz/2=',nr*nz/2
        stop
      end if
      do i=ic-1,2,-1
        jc=jc+1
        kc=kc-1
        flxr(jc)=flxr(kc)
        flxz(jc)=-flxz(kc)
      end do
      icon=jc+1
      flxr(icon)=flxr(1)
      flxz(icon)=flxz(1)
      leng(1)=0.
      do i=2,icon
        leng(i)=leng(i-1)+sqrt((flxr(i)-flxr(i-1))**2+   &
                               (flxz(i)-flxz(i-1))**2)
      end do
!  Check two adjacent points are not the same
 35   irem=0
      do i=2,icon
        if (irem.eq.0) then
          if (leng(i).eq.leng(i-1)) then
!  remove duplicated point
            irem=i
            icon=icon-1
          end if
        else
          leng(i-1)=leng(i)
          flxr(i-1)=flxr(i)
          flxz(i-1)=flxz(i)
        end if
      end do
      if (irem.gt.0) goto 35
      circumf(k)=leng(icon)
      dl=circumf(k)/(npts-1.)
      leq(1)=0.
      do i=2,npts
        leq(i)=leq(i-1)+dl
      end do
      x1=leng(icon-1)-leng(icon)
      x2=0.
      x3=leng(2)
      d1=(flxr(icon-1)-flxr(icon))/(x1-x2)
      d2=(flxr(icon)-flxr(2))/(x2-x3)
      ydd0=0.
      call spline1d(rl,leq,npts,flxr,leng,icon,work,ydd0)
      d1=(flxz(icon-1)-flxz(icon))/(x1-x2)
      d2=(flxz(icon)-flxz(2))/(x2-x3)
      ydd0=(d1-d2)*(x2-x1)/(x1-x3)+d1
      call spline1d(zl,leq,npts,flxz,leng,icon,work,ydd0)
      do i=1,npts
        rpts(k,i)=rl(i)
        zpts(k,i)=zl(i)
        rr=rl(i)
        zz=zl(i)
        bppts(k,i)=bp(rr,zz)
      end do
      curint=0.
      if (k.eq.1) then
        do i=1,npts-1
          curint=curint+dl*bppts(k,i)/mu0
        end do
        write(6,90)100.*abs((cur-curint)/cur)
 90     format('NOTE***percentage error in current is ',f7.5,'%')
        write(6,*)'If this is too large, consider using continuation'
      end if
   end subroutine flxcon
!
!*****************************************************************
!
      subroutine psical(eps,psi)
!     **************************
!
!  calculates the value of psi for the flux surface of inverse
!  aspect ratio eps and returns it in psi.
!
      use param
      implicit none
      double precision eps,psi,err,u1,eps1,diff1,u2,eps2,diff2,u3,eps3
      double precision tsterr
      integer i3,ntrip
!
!  error to be tolerated in psi
      err=errit
      if (eps.eq.0.) then
        psi=0.
        return
      end if
      if (eps.eq.tokeps) then
        psi=umax
        return
      end if
!  approx psi as a quadratic in eps for a first guess...
      u1=umax*(1.-eps/tokeps)**2
!  find the value of eps this corresponds to
      call argen(u1,eps1,0)
      diff1=(eps1-eps)/tokeps
      if (diff1.gt.0.) then
!  linear interpolate between u1 and umax
         u2=eps*u1/eps1+(eps1-eps)*umax/eps1
      else
!  linear interpolate between u1 and zero
        u2=(tokeps-eps)*u1/(tokeps-eps1)
      end if
      if ((u2.gt.umax).or.(u2.lt.0.)) then
        write(nw,*)' error in psical***diff1=',diff1
        write(nw,*)' u1=',u1,' eps1=',eps1,' eps=',eps
        stop
      end if
      ntrip=0
 5    call argen(u2,eps2,0)
      diff2=(eps2-eps)/tokeps
      if ((diff2.gt.0.).and.(diff1.gt.0.)) then
! bump up u2 until eps lies between eps1 and eps2
        u2=u2+0.1*(umax-u2)
        ntrip=ntrip+1
        if (ntrip.gt.100) then
        write(nw,*)' error in psical***cannot span U'
        stop
        end if
        goto 5
      end if
      ntrip=0
      i3=0
!  use linear interpolation between two first guesses...
 10   u3=u1+(u1-u2)*(eps-eps1)/(eps1-eps2)
      ntrip=ntrip+1
      if ((u3.gt.umax).or.(u3.lt.0.)) then
        if (ntrip.lt.100) then
          if (u3.gt.umax) then
            u3=umax
            i3=1
          end if
          if (u3.lt.0.) u3=0.
        else
          write(nw,*)' error lin inter of psical***u3=',u3
          stop
        end if
      end if
      if (i3.eq.1) then
        eps3=0.
        i3=0
      else
        call argen(u3,eps3,0)
      end if
! iterate parameters
      eps1=eps2
      eps2=eps3
      u1=u2
      u2=u3
!  check error
      tsterr=abs(eps-eps3)/(eps+eps3)
      if (tsterr.gt.err) goto 10
      psi=umax-u3
   end subroutine psical
!
!**********************************************************************
