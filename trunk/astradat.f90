module astradat_output
  implicit none
contains
      subroutine astradat
!     *******************
!
      use param
      use profiles_mod, only : dense, densi, tempe, tempi
      implicit none
!
      double precision b0,phi,rho(ncon),dpsi,te(ncon),pden(ncon),ne,   &
                       dion,arg,ti,psi
      integer k,nwr,i
!
      nwr=61
      open(unit=nwr,file=runname(1:lrunname)//'.astra', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.astra'
         stop
      endif
      b0=mu0*rodi/(2.*pi*rcen)
! Generate integral of q dpsi to evaluate rho(psi)
      phi=0.d0
      rho(ncon)=0.
      do k=ncon-1,1,-1
        dpsi=psiv(k)-psiv(k+1)
        phi=phi+(sfac(k+1)+sfac(k))*dpsi
        rho(k)=sqrt(2.*pi*rcen*phi/(mu0*rodi))
      end do
      write(nwr,*)' No. flux surfaces=',ncon
      write(nwr,*)'Rho:'
      do k=1,ncon
        write(nwr,*)rho(k)
      end do
      write(nwr,*)'1/q'
      do k=1,ncon
        write(nwr,*)1./sfac(k)
      end do
      write(nwr,*)' no. ion species='
      write(nwr,*)nimp+2
      do i=1,nimp+2
         if (i.eq.1) then
           write(nwr,*)'D density, charge=',iz(1)
         else if (i.eq.2) then
           write(nwr,*)'T density, charge=',iz(1)
         else
           write(nwr,*)'impurity density, charge=',iz(i-1)
         end if
         do k=1,ncon
           psi=psiv(k)
           if (i.eq.1) then
             write(nwr,*)densi(psi,1,0)/2.
           else if (i.eq.2) then
             write(nwr,*)densi(psi,1,0)/2.
           else
             write(nwr,*)densi(psi,i-1,0)
           end if
         end do
      end do
      write(nwr,*)'Te'
      do k=1,ncon
        psi=psiv(k)
        write(nwr,*)tempe(psi,0)
      end do
      write(nwr,*)'Ti'
      do k=1,ncon
        psi=psiv(k)
        write(nwr,*)tempi(psi,1,0)
      end do
      write(nwr,*)'power density'
      do k=1,ncon
        psi=psiv(k)
        te(k)=tempe(psi,0)
    if (te(k).gt.1.d-10) then
          ti=tempi(psi,1,0)
          ne=dense(psi,0)
          dion=1.0d-19*densi(psi,1,0)
!   include inconsistent dilution factor if no impurities...
          if (imp.eq.0) dion=dil*ne*1.0d-19/zm
          arg=-0.476*(abs(log(1.45d-5*ti)))**2.25
    else
      ne=0.
          dion=0.
          ti=0.
          arg=1.
    end if
        pden(k)=1.27d4*dion**2*exp(arg)
        write(nwr,*)pden(k)
      end do
      close(nwr)
      call astout(ncon,pden,te)
      return
   end subroutine astradat
!
!----------------------------------------------------------------------------

      subroutine astout(nma,pden,te)
!     *****************
!
      use flux_average, only : flxint
      use param
      use profiles_mod, only : fprof
      implicit none
      integer nma
      double precision rho(nma),vr(nma),g33(nma),ipol(nma),droda(nma),  &
                g22(nma),g11(nma),shif(nma),kap(nma),trian(nma)
      double precision pden(ncon),te(ncon)
!
      integer k,i,im,ipeak,nwr
      double precision v1(npts),v2(npts),v3(npts),v4(npts)
      double precision th(npts)
      double precision phi,sinarg
      double precision drhodpsi(nma),asmall(nma)
      double precision bth,rr,psin,ant,dvdpsi
      double precision ro1,ro2,ro3,x1,x2,x3,aa,bb,cc
      double precision b0,psi,zpeak,rpeak,dpsi
      double precision av0,av1,av2,dth,del0,del1,ascur
      double precision intval,intvol
!
      nwr=62
      open(unit=nwr,file=runname(1:lrunname)//'_2.astra', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_2.astra'
         stop
      endif
      b0=mu0*rodi/(2.*pi*rcen)
! Generate integral of q dpsi to evaluate rho(psi)
      phi=0.d0
      rho(ncon)=0.
      asmall(ncon)=0.
      do k=ncon-1,1,-1
        dpsi=psiv(k)-psiv(k+1)
        phi=phi+(sfac(k+1)+sfac(k))*dpsi
        rho(k)=sqrt(2.*pi*rcen*phi/(mu0*rodi))
        asmall(k)=(rpts(k,npts/2+1)-rpts(k,1))/2.
      end do
! Now we want drho/dpsi...
      do k=1,ncon
        im=k
        psin=psiv(k)
        if (im.eq.1) im=2
        if (im.eq.ncon) im=ncon-1
        ro1=rho(im-1)
        ro2=rho(im)
        ro3=rho(im+1)
        x1=psiv(im-1)
        x2=psiv(im)
        x3=psiv(im+1)
        aa=((ro1-ro2)/(x1-x2)-(ro2-ro3)/(x2-x3))/(x1-x3)
        bb=(ro1-ro2)/(x1-x2)-aa*(x1+x2)
        cc=ro1-aa*x1**2-bb*x1
        drhodpsi(k)=2.*aa*psin+bb
        x1=asmall(im-1)
        x2=asmall(im)
        x3=asmall(im+1)
        aa=((ro1-ro2)/(x1-x2)-(ro2-ro3)/(x2-x3))/(x1-x3)
        bb=(ro1-ro2)/(x1-x2)-aa*(x1+x2)
        cc=ro1-aa*x1**2-bb*x1
        droda(k)=2.*aa*asmall(k)+bb
      end do
      write(6,*)' Edge dpsi/drho=',1./drhodpsi(1)
! k=1 is edge, ncon is centre
      do k=1,ncon
        psi=psiv(k)
        zpeak=0.
        do i=1,npts
          bth=bppts(k,i)
          rr=rpts(k,i)
          v1(i)=1./bth
          v2(i)=((rcen/rr)**2)/bth
          v3(i)=bth
          v4(i)=bth*rr**2
          if (i.gt.1) then
            if (zpts(k,i).gt.zpeak) then
              zpeak=zpts(k,i)
              rpeak=rpts(k,i)
              ipeak=i
            end if
          end if
        end do
        call flxint(v1,k,ant)
        dvdpsi=2.*pi*ant
        if (k.lt.ncon) then
          vr(k)=dvdpsi/drhodpsi(k)
        else
          vr(k)=0.
        end if
        call flxint(v2,k,ant)
        if (k.lt.ncon) then
          g33(k)=2.*pi*ant/dvdpsi
        else
          g33(k)=(rcen/r0)**2
        end if
        ipol(k)=fprof(psi,2)/(rcen*b0)
        call flxint(v3,k,ant)
        g22(k)=ant*drhodpsi(k)/(2.*pi)
        g22(k)=g22(k)*rcen/(ipol(k)*sclcur)
        call flxint(v4,k,ant)
        g11(k)=2.*pi*ant*drhodpsi(k)
!  Now need shape parameters
        if (k.lt.ncon) then
          shif(k)=0.5*(rpts(k,1)+rpts(k,npts/2+1))-rcen
        else
          shif(k)=r0-rcen
        end if
        ro1=zpts(k,ipeak-1)
        ro2=zpts(k,ipeak)
        ro3=zpts(k,ipeak+1)
        x1=rpts(k,ipeak-1)
        x2=rpts(k,ipeak)
        x3=rpts(k,ipeak+1)
        aa=((ro1-ro2)/(x1-x2)-(ro2-ro3)/(x2-x3))/(x1-x3)
        bb=(ro1-ro2)/(x1-x2)-aa*(x1+x2)
        cc=ro1-aa*x1**2-bb*x1
        rpeak=-bb/(2.*aa)
        zpeak=aa*rpeak**2+bb*rpeak+cc
        if (k.lt.ncon) then
          kap(k)=zpeak/asmall(k)
        else
          kap(k)=kap(k-1)
        end if
        do i=1,npts/2+1
          sinarg=zpts(k,i)/(asmall(k)*kap(k))
          if (sinarg.gt.1.) sinarg=1.
          th(i)=asin(sinarg)
          if (th(i).lt.th(i-1)) th(i)=pi-th(i)
        end do
        av0=0.
        av1=0.
        av2=0.
        do i=1,npts/2
          dth=th(i+1)-th(i)
          av0=av0+0.5*(rpts(k,i+1)+rpts(k,i))*dth
          av1=av1+0.5*(cos(th(i+1))*rpts(k,i+1)+cos(th(i))*rpts(k,i))*dth
          av2=av2+0.5*(cos(2.*th(i+1))*rpts(k,i+1)+cos(2.*th(i))*rpts(k,i))*dth
        end do
        del0=(0.5*(rpts(k,1)+rpts(k,npts/2+1))-av0/pi)*2./asmall(k)
        del1=4.*av2/(asmall(k)*pi)
        if (k.lt.ncon) then
          trian(k)=(rcen+shif(k)-rpeak)/asmall(k)
          trian(k)=0.5*(del0+del1)
        else
          trian(k)=0.
        end if
      end do
      ascur=2.*pi*g22(1)*ipol(1)/(mu0*drhodpsi(1)*rcen)
!  Write the outputs....
      write(nwr,*)' no. flux surfaces=',ncon
      write(nwr,*)'rho'
      do k=1,ncon
        write(nwr,*)rho(k)
      end do
!!$      write(nwr,*)'apsi'
!!$      do k=1,ncon
!!$        psi=psiv(k)
!!$        write(nwr,*)-fprof(psi,1)/(mu0*rcen)
!!$      end do
!!$      write(nwr,*)'bpsi'
!!$      do k=1,ncon
!!$        psi=psiv(k)
!!$        write(nwr,*)-rcen*press(psi,1)
!!$      end do
      write(nwr,*)'vr'
      do k=1,ncon
        write(nwr,*)vr(k)
      end do
      write(nwr,*)'g33'
      do k=1,ncon
        write(nwr,*)g33(k)
      end do
      write(nwr,*)'ipol'
      do k=1,ncon
        write(nwr,*)ipol(k)
      end do
      write(nwr,*)'asmall'
      do k=1,ncon
        write(nwr,*)asmall(k)
      end do
      write(nwr,*)'g22'
      do k=1,ncon
        write(nwr,*)g22(k)
      end do
      write(nwr,*)'g11'
      do k=1,ncon
        write(nwr,*)g11(k)
      end do
      write(nwr,*)'shif'
      do k=1,ncon
        write(nwr,*)shif(k)
      end do
      write(nwr,*)'kap'
      do k=1,ncon
        write(nwr,*)kap(k)
      end do
      write(nwr,*)'trian'
      trian(ncon)=0.
      do k=1,ncon
        write(nwr,*)trian(k)
      end do
      write(nwr,*)' no. R,Z points on LCFS=',npts
      do i=1,npts
        write(nwr,*)rpts(1,i),zpts(1,i)
      end do
!   Diagnostic integrals....
      intvol=0.
      do k=1,ncon-1
        intvol=intvol+0.5*(vr(k)+vr(k+1))*(rho(k)-rho(k+1))
      end do
      write(6,*)' Volume=',intvol
      intval=0.
      do k=1,ncon-1
        intval=intval+0.5*(pden(k)*vr(k)+pden(k+1)*vr(k+1))*(rho(k)-rho(k+1))
      end do
      write(6,*)' Fus power=',intval
      intval=0.
      do k=1,ncon-1
        intval=intval+0.5*(te(k)*vr(k)+te(k+1)*vr(k+1))*(rho(k)-rho(k+1))
      end do
      write(6,*)' Fus power=',intval/intvol
      close(nwr)
      return
 end subroutine astout
end module astradat_output
