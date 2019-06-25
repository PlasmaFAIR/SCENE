      function press(psi,ip)
!     **********************
!
!   returns pressure (ip=0) or its psi derivative (ip=1) or its
!   second derivative (ip=2)
!   at a given value of psi. (defined so psi increases with minor radius)
!
      use param
      implicit none
      integer ip,id,im,ik
      double precision psin
      double precision press, psi,xps,ptem,efac
      double precision te,tempe,ti,tempi,ted,tid,tedd,tidd
      double precision pt,zni,densi,znid,pd,znidd
      double precision p1,p2,p3,x1,x2,x3,aa,bb,cc,dpsi

      xps=1.-psi/umax
      if (xps.lt.0.) then
        if (xps.gt.-1.0d-8) then
          xps=0.
        else
          write(6,*)' xps<0 in press, psi=',psi,' umax=',umax
          write(6,*)' ip=',ip
          stop
        end if
      end if
      if (xps.gt.1.) then
        if (xps.lt.1.00000001d0) then
          xps=1.
        else
          write(6,*)' xps>1 in press, psi=',psi,' umax=',umax
          write(6,*)' xps=',xps
          stop
        end if
      end if
! if ipswtch=-1 then
      if (ipswtch.eq.-1) then
        if (imp.ne.0) then
          write(6,*)' Should not use impurities with ipsw=-1'
          stop
        end if
        psin=psi/umax
        dpsi=1./nterp
        im=int(psin/dpsi)
        if (im.le.1) im=2
        if (im.gt.nterp-1) im=nterp-1
        p1=p_m(im-1)
        p2=p_m(im)
        p3=p_m(im+1)
        x1=psi_m(im-1)
        x2=psi_m(im)
        x3=psi_m(im+1)
        aa=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
        bb=(p1-p2)/(x1-x2)-aa*(x1+x2)
        cc=p1-aa*x1**2-bb*x1
        if (ip.eq.0) then
          ptem=(aa*psin**2+bb*psin+cc)
        else if (ip.eq.1) then
          ptem=(2.*aa*psin+bb)/umax
        else
          ptem=2.*aa/umax**2
        end if
        press=umax*scl*pfac*bpol*ptem/pscl
        return
      end if
      if (ipswtch.eq.-30) then
        psin=psi/umax
        dpsi=1./nterp
        im=int(psin/dpsi)
        im=1
        do 50 ik=1,nterp
          if (psin.gt.psi_m(ik)) goto 50
          im=ik
 50     continue
        if (im.le.1) im=2
        if (im.gt.nterp-1) im=nterp-1
        p1=p_m(im-1)
        p2=p_m(im)
        p3=p_m(im+1)
        x1=psi_m(im-1)
        x2=psi_m(im)
        x3=psi_m(im+1)
        aa=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
        bb=(p1-p2)/(x1-x2)-aa*(x1+x2)
        cc=p1-aa*x1**2-bb*x1
        if (ip.eq.0) then
          ptem=(aa*psin**2+bb*psin+cc)
        else if (ip.eq.1) then
          ptem=(2.*aa*psin+bb)/umax
        else
          ptem=2.*aa/umax**2
        end if
        press=umax*scl*pfac*bpol*ptem/pscl
        return
!!$          write(6,*)' in p, psi=',psi,' ip=',ip
!!$! Use T and n from TS data
!!$          if (ip.eq.0) then
!!$! return pressure
!!$            te=tempe(psi,0)
!!$            pt=0.
!!$            do id=1,nimp+1
!!$              ti=tempi(psi,id,0)
!!$              zni=densi(psi,id,0)
!!$              pt=pt+bk*(iz(id)*te+ti)*zni
!!$            end do
!!$          else if (ip.eq.1) then
!!$! return psi-derivative of pressure
!!$            te=tempe(psi,0)
!!$            ted=tempe(psi,1)
!!$            pd=0.
!!$            do id=1,nimp+1
!!$              ti=tempi(psi,id,0)
!!$              zni=densi(psi,id,0)
!!$              tid=tempi(psi,id,1)
!!$              znid=densi(psi,id,1)
!!$              pd=pd+bk*(iz(id)*ted+tid)*zni+bk*(iz(id)*te+ti)*znid
!!$            end do
!!$            pt=pd
!!$          else if (ip.eq.2) then
!!$! return second psi-derivative of pressure
!!$            te=tempe(psi,0)
!!$            ted=tempe(psi,1)
!!$            tedd=tempe(psi,2)
!!$            pd=0.
!!$            do id=1,nimp+1
!!$              ti=tempi(psi,id,0)
!!$              zni=densi(psi,id,0)
!!$              tid=tempi(psi,id,1)
!!$              znid=densi(psi,id,1)
!!$              tidd=tempi(psi,id,2)
!!$              znidd=densi(psi,id,2)
!!$              pd=pd+bk*(iz(id)*tedd+tidd)*zni    &
!!$                   +2.*bk*(iz(id)*ted+tid)*znid     &
!!$                   +bk*(iz(id)*te+ti)*znidd
!!$            end do
!!$            pt=pd
!!$          end if
!!$          press=pt
!!$          write(6,*)' out p'
!!$        return
      end if
      if (ipswtch.eq.1) then
        if (imp.eq.0) then
          if (ip.eq.0) then
!  return total pressure
            ptem=pa+(1.-pa)*xps**ppow
          else if (ip.eq.1) then
!  return the derivative
            ptem=-ppow*(1.-pa)*xps**(ppow-1.)/umax
          else
!  return the second derivative
            ptem=(ppow*(ppow-1.)*(1.-pa)/umax**2)*xps**(ppow-2.)
          end if
          press=umax*scl*pfac*bpol*ptem/pscl
          return
        else
          if (ip.eq.0) then
! return pressure
            te=tempe(psi,0)
            pt=0.
            do id=1,nimp+1
              ti=tempi(psi,id,0)
              zni=densi(psi,id,0)
              pt=pt+bk*(iz(id)*te+ti)*zni
            end do
          else if (ip.eq.1) then
! return psi-derivative of pressure
            te=tempe(psi,0)
            ted=tempe(psi,1)
            pd=0.
            do id=1,nimp+1
              ti=tempi(psi,id,0)
              zni=densi(psi,id,0)
              tid=tempi(psi,id,1)
              znid=densi(psi,id,1)
              pd=pd+bk*(iz(id)*ted+tid)*zni+bk*(iz(id)*te+ti)*znid
            end do
            pt=pd
          else if (ip.eq.2) then
! return second psi-derivative of pressure
            te=tempe(psi,0)
            ted=tempe(psi,1)
            tedd=tempe(psi,2)
            pd=0.
            do id=1,nimp+1
              ti=tempi(psi,id,0)
              zni=densi(psi,id,0)
              tid=tempi(psi,id,1)
              znid=densi(psi,id,1)
              tidd=tempi(psi,id,2)
              znidd=densi(psi,id,2)
              pd=pd+bk*(iz(id)*tedd+tidd)*zni    &
                   +2.*bk*(iz(id)*ted+tid)*znid     &
                   +bk*(iz(id)*te+ti)*znidd
            end do
            pt=pd
          end if
          press=pt
        end if
        return
      end if
!***************************************************************
      if (ipswtch.eq.5) then
        if (imp.ne.0) then
          write(6,*)' cannot use ipswtch=5 option with impurities'
          stop
        end if
        if (ip.eq.0) then
!  return total pressure
          ptem=-(1.2*psi/umax+(39./4.)*(psi/umax)**4-(40.2/5.)*(psi/umax)**5 &
              -1.2-(39/4.)+(40.2/5.))
        else if (ip.eq.1) then
!  return the derivative
          ptem=-(1.2+39.*(psi/umax)**3-40.2*(psi/umax)**4)/umax
        else
!  return the second derivative
          ptem=-(3.*39.*(psi/umax)**2-4.*40.*(psi/umax)**3)/umax**2
        end if
        press=umax*scl*pfac*bpol*ptem
        return
      end if
!  ipswtch=0
      if (imp.eq.0) then
        efac=exp(pedg*xps)
        if (ip.eq.0) then
!  return total pressure
          ptem=pa+(pped-pa)*(efac**2-1.)/(efac**2+1.)
          if (ppow.gt.1) ptem=ptem+pn*((1.+xps)**ppow-(1.+ppow*xps))
        else if (ip.eq.1) then
!  return the derivative
          ptem=-4.*(pedg/umax)*(pped-pa)/(efac+1./efac)**2
          if (ppow.gt.1.) ptem=ptem-(ppow/umax)*pn*((1.+xps)**(ppow-1.)-1.)
        else
!  return the second derivative
          ptem=-2.*(pedg/umax)**2*(pped-pa)*  &
                ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
          if (ppow.gt.1.) &
             ptem=ptem+(pn*ppow*(ppow-1.)/umax**2)*(1.+xps)**(ppow-2.)
        end if
!        press=umax*bpol*scl*ptem/ppow
        press=umax*scl*pfac*bpol*ptem/pscl
        return
     end if

     if (ipswtch.eq.15) then
        !Use Lao's et al method for p'
        !ptem




     else


        if (ip.eq.0) then
! return pressure
          te=tempe(psi,0)
          pt=0.
          do id=1,nimp+1
            ti=tempi(psi,id,0)
            zni=densi(psi,id,0)
            pt=pt+bk*(iz(id)*te+ti)*zni
          end do
        else if (ip.eq.1) then
! return psi-derivative of pressure
          te=tempe(psi,0)
          ted=tempe(psi,1)
          pd=0.
          do id=1,nimp+1
            ti=tempi(psi,id,0)
            zni=densi(psi,id,0)
            tid=tempi(psi,id,1)
            znid=densi(psi,id,1)
            pd=pd+bk*(iz(id)*ted+tid)*zni+bk*(iz(id)*te+ti)*znid
          end do
          pt=pd
        else if (ip.eq.2) then
! return second psi-derivative of pressure
          te=tempe(psi,0)
          ted=tempe(psi,1)
          tedd=tempe(psi,2)
          pd=0.
          do id=1,nimp+1
            ti=tempi(psi,id,0)
            zni=densi(psi,id,0)
            tid=tempi(psi,id,1)
            znid=densi(psi,id,1)
            tidd=tempi(psi,id,2)
            znidd=densi(psi,id,2)
            pd=pd+bk*(iz(id)*tedd+tidd)*zni    &
                 +2.*bk*(iz(id)*ted+tid)*znid     &
                 +bk*(iz(id)*te+ti)*znidd
          end do
          pt=pd
        end if
!        scl=s*bpol*umax/p0
        press=pt
     end if



  end function press
!
!**********************************************************************
!
       function pressi(psi,id,ip)
!
!    For ip=0/1/2, returns pressure/first pressure derivative/second
!    pressure derivative at a given psi, defined such that psi
!    increases with the minor radius.
!
       use param
       implicit none
       double precision psin
       double precision pressi,press,psi,xps,ptem,efac
       double precision p1,p2,p3,x1,x2,x3,aa,bb,cc,dpsi
       integer ip,id,im,ik

       xps=1.-psi/umax
       if (xps.lt.0.) then
          if (xps.gt.-1.0d-8) then
             xps=0
          else
             Write(6,*)'xps<0 in press, psi=',psi,' umax=',umax
             Write(6,*)' ip=',ip
             Stop
          end if
       end if
       if (ipswtch.eq.-1) then
          if (imp.ne.0) then
             Write(6,*)' Should not use impurities with ipswtch=-1.'
             Stop
          end if
          pressi=press(psi,ip)
          return
       end if

       if (ipswtch.eq.-30) then
          psin=psi/umax
          dpsi=1./nterp
          im=int(psin/dpsi)
          im=1
          do 60 ik=1,nterp
             if (psin.gt.psi_m(ik)) goto 60
             im=ik
60        continue
          if (im.le.1) im=2
          if (im.gt.nterp-1) im=nterp-1
          p1=p_m(im-1)
          p2=p_m(im)
          p3=p_m(im+1)
          x1=psi_m(im-1)
          x2=psi_m(im)
          x3=psi_m(im+1)
          aa=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
          bb=(p1-p2)/(x1-x2)-aa*(x1+x2)
          cc=p1-aa*x1**2-bb*x1
          if (ip.eq.0) then
             ptem=(aa*psin**2+bb*psin+cc)
          else if (ip.eq.1) then
             ptem=(2.*aa*psin+bb)/umax
          else
             ptem=2.*aa/umax**2
          end if
          pressi=umax*scl*pfac*bpol*ptem/pscl
          return
       end if

       if (ipswtch.eq.1) then
          if (imp.eq.0) then
             pressi=press(psi,ip)
          else
             if (ip.eq.0) then
                ptem=pa+(1.-pa)*xps**ppow
             else if (ip.eq.1) then
                ptem=-ppow*(1.-pa)*xps**(ppow-1.)/umax
             else
                ptem=(ppow*(ppow-1.)*(1.-pa)/umax**2)*xps**(ppow-2.)
             end if
             pressi=umax*scl*pfac*bpol*ptem/pscl
             return
          end if
       end if

       if (ipswtch.eq.5) then
          if (imp.ne.0) then
             Write(6,*)' Cannot use ipswtch=5 with impurities.'
             Stop
          else
             pressi=press(psi,ip)
             return
          end if
       end if

       if (ipswtch.eq.0) then
          if (imp.eq.0) then
             pressi=press(psi,ip)
             return
          else
             efac=exp(pedg*xps)
             if (ip.eq.0) then
                ptem=pa+(pped-pa)*(efac**2-1.)/(efac**2+1.)
                if (ppow.gt.1) ptem=ptem+pn*((1.+xps)**ppow-(1.+ppow*xps))
             else if (ip.eq.0) then
                ptem=-4.*(pedg/umax)*(pped-pa)/(efac+1./efac)**2
                if (ppow.gt.1.) ptem=ptem-(ppow/umax)*pn*((1.+xps)**(ppow-1.)-1.)
             else
                ptem=-2.*(pedg/umax)**2*(pped-pa)*((efac**2-1.)/(efac**2+1.))/(0.5*(efac+1./efac))**2
                if (ppow.gt.1.) then
                   ptem=ptem+(pn*ppow*(ppow-1.)/umax**2)*(1.+xps)**(ppow-2.)
                end if
             end if
             pressi=umax*scl*pfac*bpol*ptem/pscl
             return
          end if
       end if

 end function pressi
!
!**********************************************************************
!
       function dense(psi,i)
!      *********************
!
!  Returns electron density if I=0 and its PSI derivative if I=1
!  The value of SCL is calculated during the equilibrium calculation
!  to give the required total current. Should only be called if
!  impurities are present.
!
      use param
      implicit none
      double precision psi,dense,densi,tempe,tempi,press
      double precision ne,ned,te,ti,ted,tid
      integer i,j
!
      if (imp.eq.0) then
        ti=tempi(psi,1,0)
        te=tempe(psi,0)
        if ((ti.lt.1e-8).and.(te.lt.1e-8)) then
          ne=0.
        else
          ne=press(psi,0)/(zm*ti+te)
        end if
        if (i.eq.0) then
          dense=ne/bk
          return
        end if
        tid=tempi(psi,1,1)
        ted=tempe(psi,1)
        ned=(press(psi,1)-(zm*tid+ted)*ne)/(zm*ti+te)
        if (i.eq.1) then
          dense=ned/bk
          return
        else
       write(6,*)'ERROR*** need to code up higher derivatives of e density'
          stop
        end if
      end if
      dense=0.
      if (i.eq.0) then
        do j=1,nimp+1
          dense=dense+iz(j)*densi(psi,j,0)
!          if (psi.lt.0.000000001) write(6,*)' j=',j,' dense', dense
        end do
      else if (i.eq.1) then
        do j=1,nimp+1
          dense=dense+iz(j)*densi(psi,j,1)
        end do
      else
        do j=1,nimp+1
          dense=dense+iz(j)*densi(psi,j,2)
        end do
      end if
    end function dense
!
!**********************************************************************
!
      function densi(psi,id,i)
!     ************************
!
!  Returns ion density if i=0 and its PSI derivative if i=1, and
!  its second psi-derivative if i=2.
!  The value of SCL is calculated during the equilibrium calculation
!  to give the required total current. ID labels the ion
!  species with ID=1 the main species ion. Should only be
!  called if impurities are present.
!
      use param
      implicit none
      double precision densi,psi,efac
      double precision ni,nid,nidd,ti,te,xps,ted,tid
      double precision press,tempi,tempe
      double precision nta,ntped,ntedg,nin,ntpow,nz0
      double precision psin,dpsi,n1,n2,n3,aa,bb,cc
      double precision ps1,ps2,ps3
      integer i,id,ip,ik
!
      if (imp.eq.0) then
        ti=tempi(psi,1,0)
        te=tempe(psi,0)
        ni=zm*press(psi,0)/(zm*ti+te)
        if (i.eq.0) then
          densi=ni/bk
          return
        end if
        tid=tempi(psi,1,1)
        ted=tempe(psi,1)
        nid=(press(psi,1)-(zm*tid+ted)*ni/zm)/(zm*ti+te)
        if (i.eq.1) then
          densi=nid/bk
          return
        else
       write(6,*)'ERROR*** need to code up higher derivatives of ion density'
          stop
        end if
     end if
      xps=1.-psi/umax
      ntpow=znpow(id)
      nz0=zn0(id)
      if (ipswtch.eq.1) then
        if (i.eq.0) then
!  return density
          ni=nz0*xps**ntpow
          densi=scl*ni*umax*pfac*bpol/pscl
          return
        else if (i.eq.1) then
!  return the derivative
          nid=(-ntpow/umax)*nz0*xps**(ntpow-1.)
          densi=scl*nid*umax*pfac*bpol/pscl
          return
        else
!  return the second derivative
          nidd=(ntpow*(ntpow-1.)/umax**2)*nz0*xps**(ntpow-2.)
          densi=scl*nidd*umax*pfac*bpol/pscl
        end if
        return
      end if
      if (ipswtch.eq.3) then
!  Use input data for profile, scaled by n0 factor
        psin=psi/umax
        dpsi=1./nterp
        ip=1+int(psin/dpsi)
        ip=1
        do 50 ik=1,nterp
          if (psin.gt.psi_m(ik)) goto 50
          ip=ik
 50     continue
        if (ip.gt.nterp-1) ip=nterp-1
        if (ip.eq.1) ip=ip+1
        n1=ne_m(ip-1)
        ps1=psi_m(ip-1)
        n2=ne_m(ip)
        ps2=psi_m(ip)
        n3=ne_m(ip+1)
        ps3=psi_m(ip+1)
        aa=((n1-n2)/(ps1-ps2)-(n2-n3)/(ps2-ps3))/(ps1-ps3)
        bb=(n1-n2)/(ps1-ps2)-aa*(ps1+ps2)
        cc=n2-aa*psi_m(ip)**2-bb*psi_m(ip)
        if (i.eq.0) then
! return density
          ni=nz0*(aa*psin**2+bb*psin+cc)
          densi=scl*ni*umax*pfac*bpol/pscl
!          write(56,*)' densi=',densi,' n1=',n1,' n2=',n2,' n3=',n3
          return
        else if (i.eq.1) then
! return derivative
          nid=nz0*(2.*aa*psin+bb)/umax
          densi=scl*nid*umax*pfac*bpol/pscl
          return
        else
! return second derivative
          nidd=2.*nz0*aa/umax**2
          densi=scl*nidd*umax*pfac*bpol/pscl
        end if
        return
      end if
!--------------------------------------------------------------
      nta=zna(id)
      ntped=znped(id)
      ntedg=znedg(id)
      ntpow=znpow(id)
      if (ntpow.gt.1.) then
      nin=(zn0(id)-nta-(ntped-nta)*(exp(ntedg)-exp(-ntedg)) &
          /(exp(ntedg)+exp(-ntedg)))/(2**ntpow-1.-ntpow)
      end if
      efac=exp(ntedg*xps)
      if ((ipswtch.eq.0).or.(ipswtch.eq.2).or.   &
           (ipswtch.eq.5).or.(ipswtch.eq.7).or.(ipswtch.eq.8)) then
        if (i.eq.0) then
!  return density
          ni=nta+(ntped-nta)*(efac**2-1.)/(efac**2+1.)
          if (ntpow.gt.1.) ni=ni+nin*((1.+xps)**ntpow-(1.+ntpow*xps))
          densi=scl*ni*umax*pfac*bpol/pscl
          return
        else if (i.eq.1) then
!  return the derivative
          nid=-4.*(ntedg/umax)*(ntped-nta)/  &
               (efac+1./efac)**2
          if (ntpow.gt.1.) nid=nid-(ntpow/umax)*nin*((1.+xps)**(ntpow-1.)-1.)
          densi=scl*nid*umax*pfac*bpol/pscl
          return
        else
!  return the second derivative
          nidd=-2.*(ntedg/umax)**2*(ntped-nta)*  &
                 ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
          if (ntpow.gt.1.) &
             nidd=nidd+(nin*ntpow*(ntpow-1.)/umax**2)*(1.+xps)**(ntpow-2.)
          densi=scl*nidd*umax*pfac*bpol/pscl
        end if
      end if
      if (ipswtch.eq.9) then
        if (i.eq.0) then
!  return density
          ni=nta+(ntped-nta)*(efac**2-1.)/(efac**2+1.)
          if (id.eq.1) ni=ni+naa*xps**2+nbb*xps**3
          densi=scl*ni*umax*pfac*bpol/pscl
          return
        else if (i.eq.1) then
!  return the derivative
          nid=-4.*(ntedg/umax)*(ntped-nta)/  &
               (efac+1./efac)**2
          if (id.eq.1) nid=nid-2.*naa*xps/umax-3.*nbb*xps**2/umax
          densi=scl*nid*umax*pfac*bpol/pscl
          return
        else
!  return the second derivative
          nidd=-2.*(ntedg/umax)**2*(ntped-nta)*  &
                 ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
          if (ntpow.gt.1.) then
             if (id.eq.1)   &
             nidd=nidd+2.*naa/umax**2+6.*nbb*xps/umax**2

          end if
          densi=scl*nidd*umax*pfac*bpol/pscl
        end if
      end if
      if (ipswtch.eq.6) then

        nin=(zn0(id)-nta-(ntped-nta)*(exp(ntedg)-exp(-ntedg))/ &
                              (exp(ntedg)+exp(-ntedg)))
        if (i.eq.0) then
!  return density
           ni=nta+(ntped-nta)*(efac**2-1.)/ &
                (efac**2+1.)
          ni=ni+nin*((1.+ane)*xps**ntpow-ane*xps**npo1)
          densi=scl*ni*umax*pfac*bpol/pscl
          return
        else if (i.eq.1) then
!  return the derivative
          if ((ntpow.lt.1.).or.(npo1.lt.1.)) then
            write(6,*)' Cannot calculate derivative at psi=psi_edge'
            write(6,*)' for ntpow or npo1<1'
            stop
          end if
          nid=-4.*(ntedg/umax)*(ntped-nta)/(efac+1./efac)**2
          nid=nid-(nin/umax)*(ntpow*(1.+ane)*xps**(ntpow-1.)         &
               -npo1*ane*xps**(npo1-1.))
          densi=scl*nid*umax*pfac*bpol/pscl
          return
        else
!  return the second derivative
          if (xps.le.0.) then
            if ((ntpow.lt.2.).or.(npo1.lt.2.)) then
              write(6,*)' Cannot calculate 2nd derivative at psi=psi_edge'
              write(6,*)' for ntpow or npo1<2*****using xps=0.00001'
            end if
            xps=0.00001
          end if
          nidd=-2.*(ntedg/umax)**2*(ntped-nta)*  &
                ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
          nidd=nidd+(nin/umax**2)*(ntpow*(ntpow-1.)*(1.+ane)*xps**(ntpow-2.)  &
               -npo1*(npo1-1.)*ane*xps**(npo1-2.))
          densi=scl*nidd*umax*pfac*bpol/pscl
        end if
      end if
      if (ipswtch.eq.7) then
        if (npo1.lt.2.) then
          write(6,*)' Error, must set npo1>/=2'
          stop
        end if
!        write(6,*)' nz0=',nz0,' nta=',nta,' xps=',xps
!        write(6,*)' ane=',ane,' ntped=',ntped,' npo1=',npo1,' umax=',umax
        if (i.eq.0) then
! return density
         ni=(nz0-nta)*xps**ntpow+ntped*((xps-ane)**2)**(npo1/2.)+nta
         densi=scl*ni*umax*pfac*bpol/pscl
         return
        else if (i.eq.1) then
! return derivative
         if (abs(ntpow-1.).lt.1e-8) then
           nid=-((nz0-nta)+npo1*ntped*((xps-ane)**2)**((npo1-1.)/2.))/umax
         else if (abs(ntpow).lt.1e-8) then
           nid=-npo1*ntped*((xps-ane)**2)**((npo1-1.)/2.)/umax
         else
           if (ntpow.lt.2.) then
             write(6,*)' Error, if npow<2, must be integer'
             stop
           end if
           nid=-(ntpow*(nz0-nta)*xps**(ntpow-1.)            &
               +(xps-ane)*npo1*ntped*((xps-ane)**2)**((npo1-2.)/2.))/umax
         end if
         densi=scl*nid*umax*pfac*bpol/pscl
         return
        else if (i.eq.2) then
! return 2nd derivative
         if ((abs(ntpow-1).lt.1e-8).or.(abs(ntpow).lt.1e-8)) then
           nidd=npo1*(npo1-1.)*ntped*((xps-ane)**2)**((npo1-2.)/2.)/umax**2
         else
           if (ntpow.lt.2.) then
             write(6,*)' Error, if npow<2, must be integer'
             stop
           end if
           nidd=(ntpow*(ntpow-1.)*(nz0-nta)*xps**(ntpow-2.)            &
               +npo1*(npo1-1.)*ntped*((xps-ane)**2)**((npo1-2.)/2.))/umax**2
         end if
         densi=scl*nidd*umax*pfac*bpol/pscl
         return
        end if
      end if
  end function densi
!
!**********************************************************************
      function fprof(psi,id)
!     **********************
!
!  if abs(id)=1 returns g(psi) (=ff'), abs(id)=2 returns f,
!  abs(id)=3 returns f', abs (id)=4 returns (ff')'
!   if ipass=0 uses input profile parameterisation
!   if ipass>0 uses psi mesh set up from previous eqbm. generation.
!
      use param
      implicit none
      integer id,ij,j,im,i0,ip
      double precision psi,fprof
      double precision xps,ffp,fsq,f,g0,fff,ffpp,ff1,ff2
      double precision rat,ds,fup,ffpint,psin
      double precision x1,x2,x3,g1,g2,g3,aa,bb,cc,dpsi
!
      xps=1.-psi/umax
      psin=psi/umax
      if (ipass.eq.0) then
!  use input profile parameterisation
        g0=-mu0*scl*r0*r0*(1.-bpol)
! if ipswtch=-1 then use mesh profile for ffprime
        if (ipswtch.eq.-1) then
          psin=psi/umax
          dpsi=1./nterp
          im=int(psin/dpsi)
          if (im.le.1) im=2
          if (im.gt.nterp-1) im=nterp-1
          g1=ffp_m(im-1)
          g2=ffp_m(im)
          g3=ffp_m(im+1)
          x1=psi_m(im-1)
          x2=psi_m(im)
          x3=psi_m(im+1)
          aa=((g1-g2)/(x1-x2)-(g2-g3)/(x2-x3))/(x1-x3)
          bb=(g1-g2)/(x1-x2)-aa*(x1+x2)
          cc=g1-aa*x1**2-bb*x1
          ffp=g0*(aa*psin**2+bb*psin+cc)
          if (id.eq.1) then
             fprof=ffp
            return
          else if (id.lt.4) then
            fsq=0.
            if (nterp.gt.(im+2)) then
              do j=nterp,im+2,-1
                fsq=fsq-0.5*(ffp_m(j)+ffp_m(j-1))*dpsi
              end do
            else
              fsq=0.
            end if
            fsq=fsq+((1./3.)*aa*psin**3+0.5*bb*psin**2+cc*psin)- &
               ((1./3.)*aa*psi_m(im+1)**3+0.5*bb*psi_m(im+1)**2+cc*psi_m(im+1))
            fsq=2.*fsq*g0*umax+const
            if (fsq.lt.0.) then
              write(6,*)' ERROR for fsq in fprof...check, fsq=',fsq
              stop
            end if
            fprof=sqrt(fsq)
            if (id.eq.3) fprof=ffp/sqrt(fsq)
            return
          else
            fprof=(g0/umax)*(2.*aa*psin+bb)
            return
          end if
        end if
        if (ipswtch.eq.1) then
          ffp=g0*(xps**fpow+af1*psin**fpow1+af2*psin**fpow2)
          if (id.eq.1) then
           fprof=ffp
           return
          end if
          fff=fpow+1.
          ff1=fpow1+1
          ff2=fpow2+1.
          fsq=(2.*umax*g0)*((-1./fff)*xps**fff    &
                            +(af1/ff1)*psin**ff1+ &
                            (af2/ff2)*psin**ff2)
          fsq=fsq+const
          if (fsq.lt.0.) then
            write(6,*)'ERROR****need to reduce bpol as F^2<0 in fprof'
            stop
          end if
          f=sqrt(fsq)
          if (id.eq.2) then
            fprof=f
            return
          end if
          if (id.eq.3) then
            fprof=ffp/f
            return
          end if
          ffpp=(g0/umax)*(-fpow*xps**(fpow-1.)          &
                          +af1*fpow1*psin**(fpow1-1.)   &
                          +af2*fpow2*psin**(fpow2-1.))
          fprof=ffpp
          return
        end if
!************************************************************
        if (ipswtch.eq.5) then
          ffp=g0*xps**fpow
          if (id.eq.1) then
            fprof=ffp
            return
          end if
          fff=fpow+1
          fsq=(2.*umax*g0/fff)*xps**fff
          fsq=fsq+const
          if (fsq.lt.0.) then
            write(6,*)'ERROR****need to reduce bpol as F^2<0 in fprof'
            stop
          end if
          f=sqrt(fsq)
          if (id.eq.2) then
            fprof=f
            return
          end if
          if (id.eq.3) then
            fprof=ffp/f
            return
          end if
          ffpp=-(g0/umax)*fpow*xps**(fpow-1.)
          fprof=ffpp
          return
       end if

       !Test function from Lao et al
       if (ipswtch.eq.15) then

          ffp = g0 * ( af1*psin**1 + af2*psin**2 + &
               af3*psin**3 - (af1+af2+af3)*psin**4)
          if (id.eq.1) then
             fprof=ffp
             return
          end if

          fsq = umax*g0* (af1/2.*psin**2 + af2/3.* &
               psin**3 + af3/4.*psin**4 - (af1+af2+af3) &
               /5. *psin**5)

          fsq = fsq + const
          f = sqrt(fsq)
          if (id.eq.2) then
             fprof=f
             return
          end if
          if (id.eq.3) then
             fprof = ffp/f
             return
          end if

          ffpp = g0/umax * (af1 + 2*af2*psin + 3*af3*psin**2 - &
               4*(af1+af2+af3)*psin**3)

          fprof=ffpp
          return
       end if



        ffp=g0*((1.+xps)**fpow-1.)
        if (id.eq.1) then
          fprof=ffp
          return
        end if
        fff=fpow+1
        fsq=-(2.*umax*g0/fff)*((1.+xps)**fff-fff*xps-1.)
        fsq=fsq+const
        if (fsq.lt.0.) then
          write(6,*)'ERROR****need to reduce bpol as F^2<0 in fprof'
          stop
        end if
        f=sqrt(fsq)
        if (id.eq.2) then
          fprof=f
          return
        end if
        if (id.eq.3) then
          fprof=ffp/f
          return
        end if
        ffpp=-(g0/umax)*fpow*(1.+xps)**(fpow-1.)
        fprof=ffpp
     else
!---------------------------------------------------------------
!  use discrete form of ff'
        ij=1
        xps=psi/umax
        if ((xps.gt.1).or.(xps.lt.0.)) then
          if (abs(xps).lt.1.0d-8) then
             xps=1.0d-13
          else if (abs(xps-1.).lt.1.0d-8) then
             xps=1.
          else
            write(6,*)' ERROR***xps is outside allowed range in fprof'
            write(6,*)' xps=',xps,' psi=',psi,' umax=',umax
            stop
          end if
        end if
        if (xps.lt.1.0d-13) then
          ij=ncon
        else
          if (xps.lt.psiold(ncon/2)) ij=ncon/2
          do while (psiold(ij).ge.xps)
            ij=ij+1
            if (ij.gt.ncon) then
              write(6,*)' error in calculating xps'
              write(6,*)' psiold(ncon)=',psiold(ncon),' xps=',xps
            end if
          end do
        end if
!!$	do 10 j=1,ncon
!!$	  if (psiold(j).le.xps) goto 10
!!$	  ij=j
!!$ 10     continue
!        if (ij.eq.ncon) ij=ncon-1
        if (ij.eq.1) ij=2
        rat=(xps-psiold(ij-1))/(psiold(ij)-psiold(ij-1))
        if ((rat.gt.1.).or.(rat.lt.0.)) then
          write(6,*)' rat=',rat,' xps=',xps,' ij=',ij
          write(6,*)' psiold(ij)=',psiold(ij),' psiold(ij-1)=',psiold(ij-1)
       end if

        ffp=gst(ij-1)+rat*(gst(ij)-gst(ij-1))
        ffp=scl*ffp
!!$        if (xps.gt.0.999) then
!!$          write(6,*)' ******* ij=',ij,' scl=',scl,' ffp=',ffp,' rat=',rat
!!$          write(6,*)' gst(ij)=',gst(ij)
!!$        end if
        if (id.eq.1) then
           fprof=ffp

          return
        end if
        if (id.eq.4) then
!  need derivative of ffp
          im=ij-1
          i0=ij
          ip=ij+1
          if (im.lt.1) then
            im=im+1
            i0=i0+1
            ip=ip+1
          else if (ip.gt.ncon) then
            im=im-1
            i0=i0-1
            ip=ip-1
          end if
          x1=psiv(im)
          x2=psiv(i0)
          x3=psiv(ip)
          g1=gst(im)
          g2=gst(i0)
          g3=gst(ip)
          aa=((g1-g2)/(x1-x2)-(g2-g3)/(x2-x3))/(x1-x3)
          bb=(g1-g2)/(x1-x2)-aa*(x1+x2)
          fprof=scl*(2.*aa*psi+bb)
          return
       end if
!  integrate up to calculate f
        ffpint=0.
        if (ij.gt.2) then
          do j=1,ij-2
            ds=psiv(j+1)-psiv(j)
            ffpint=ffpint+0.5*(gst(j+1)+gst(j))*ds
          end do
        end if
        if (ij.ne.1) then
          ds=psiv(ij)-psiv(ij-1)
          fup=ffpint+0.5*(gst(ij)+gst(ij-1))*ds
   ffpint=ffpint+rat*(fup-ffpint)

          fsq=2.*scl*ffpint+const
          if (fsq.lt.0.) then
      write(nw,*)'input error1*** f**2<0 in fprof'
      write(nw,*) fsq,2.*scl*ffpint, const
            stop
          end if
          f=sqrt(fsq)
        else
          fsq=const
          if (fsq.lt.0.) then
            write(nw,*)'input error2*** f**2<0 in fprof'
            stop
          end if
          f=sqrt(fsq)
        end if
        if (id.eq.2) then
          fprof=f
          return
        else
          fprof=ffp/f
          return
       end if
    end if
  end function fprof
!
!**********************************************************************
!
      function tempe(psi,i)
!     *********************
!
!  Function defines the electron temperature profile as a function
!  of the flux.
!  T0E is the central temperature. I=0 gives the temperature and I=1
!  is the derivative w.r.t. Psi
!
!
      use param
      implicit none
!
      integer i,ip,ik
      double precision tempe,psi,xps,tem,efac,tfac
      double precision psin,aa,bb,cc,dpsi,t1,t2,t3
      double precision ps1,ps2,ps3
!
      xps=1.-psi/umax
      efac=exp(teedg*xps)
      if (xps.lt.0.) then
        if (xps.gt.-1.0d-7) then
          xps=0.
        else
          write(6,*)' xps<0 in tempe, psi=',psi,' umax=',umax
          stop
        end if
      end if
      if (xps.gt.1.) then
        if (xps.lt.1.000001) then
          xps=1.
        else
          write(6,*)' xps>1 in tempe, psi=',psi,' umax=',umax
          stop
        end if
      end if
      if ((ipswtch.eq.1)) then
        if (i.eq.0) then
!  return electron temperature
          tem=tea+(te0-tea)*xps**tpoe
        else if (i.eq.1) then
!  return the derivative
          tem=-((te0-tea)*tpoe/umax)*xps**(tpoe-1.)
        else
!  return the second derivative
          tem=((te0-tea)*tpoe*(tpoe-1.)/umax**2)*xps**(tpoe-2.)
        end if
        tempe=tem
        return
      end if
      if (ipswtch.eq.2) then
        tfac=2.**tpoe-(1.+tpoe)
        if (i.eq.0) then
!  return electron temperature
          tem=tea+(teped-tea)*(efac**2-1.)/ &
                  (efac**2+1.)
          if (tpoe.gt.1) tem=tem+ten*tfac*(xps**tpoe)
        else if (i.eq.1) then
!  return the derivative
          tem=-4.*(teedg/umax)*(teped-tea)/(efac+1./efac)**2
          if (tpoe.gt.1) tem=tem-(tpoe*tfac/umax)*ten*xps**(tpoe-1.)
        else
!  return the second derivative
          tem=-2.*(teedg/umax)**2*(teped-tea)*  &
                  ((efac**2-1.)/(efac**2+1.))/ &
                  (0.5*(efac+1./efac))**2
          if (tpoe.gt.1.)    &
              tem=tem+(ten*tfac*tpoe*(tpoe-1.)/umax**2)*xps**(tpoe-2.)
        end if
        tempe=tem
        return
      end if
      if (ipswtch.eq.3) then
!  Use discrete input for profile
        psin=psi/umax
        dpsi=1./nterp
        ip=1+int(psin/dpsi)
        ip=1
        do 50 ik=1,nterp
          if (psin.gt.psi_m(ik)) goto 50
          ip=ik
 50     continue
        if (ip.gt.nterp-1) ip=nterp-1
        if (ip.eq.1) ip=ip+1
        t1=te_m(ip-1)
        ps1=psi_m(ip-1)
        t2=te_m(ip)
        ps2=psi_m(ip)
        t3=te_m(ip+1)
        ps3=psi_m(ip+1)
        aa=((t1-t2)/(ps1-ps2)-(t2-t3)/(ps2-ps3))/(ps1-ps3)
        bb=(t1-t2)/(ps1-ps2)-aa*(ps1+ps2)
        cc=t2-aa*psi_m(ip)**2-bb*psi_m(ip)
        if (i.eq.0) then
! return electron temperature
          tem=aa*psin**2+bb*psin+cc
        else if (i.eq.1) then
! return derivative
          tem=(2.*aa*psin+bb)/umax
        else
! return second derivative
          tem=2.*aa/umax**2
        end if
        tempe=tem
        return
      end if
      if (ipswtch.eq.5) then
        ate=1.5
        tpo1=2.0
        if (i.eq.0) then
          tempe=te0*((1.+ate)*xps**tpoe-ate*xps**tpo1)+tea
        else if (i.eq.1) then
          tempe=-(te0/umax)*(tpoe*(1.+ate)*xps**(tpoe-1.)         &
               -tpo1*ate*xps**(tpo1-1.))
        else
          tempe=(te0/umax**2)*(tpoe*(tpoe-1.)*(1.+ate)*xps**(tpoe-2.)     &
               -tpo1*(tpo1-1.)**ate*xps**(tpo1-2.))
        end if
        return
      end if
      if ((ipswtch.eq.6).or.(ipswtch.eq.8)) then
!!$        ate=2.1
!!$        tpo1=3.0
        if (tpoe.lt.2.) then
          write(6,*)' ERROR***TPOE must be >2, TPOE=',tpoe
          stop
        end if
        if (tpo1.lt.2.) then
          write(6,*)' ERROR***TPO1 must be >2, TPO1=',tpo1
          stop
        end if
        if (i.eq.0) then
           tem=tea+(teped-tea)*(efac**2-1.)/ &
                (efac**2+1.)
          tempe=tem+ten*((1.+ate)*xps**tpoe-ate*xps**tpo1)
        else if (i.eq.1) then
          tem=-4.*(teedg/umax)*(teped-tea)/(efac+1./efac)**2
          tempe=tem-(ten/umax)*(tpoe*(1.+ate)*xps**(tpoe-1.)         &
               -tpo1*ate*xps**(tpo1-1.))
        else
          tem=-2.*(teedg/umax)**2*(teped-tea)*  &
                ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
          tempe=tem+(ten/umax**2)*(tpoe*(tpoe-1.)*(1.+ate)*xps**(tpoe-2.)  &
               -tpo1*(tpo1-1.)*ate*xps**(tpo1-2.))
        end if
        return
      end if
      if (i.eq.0) then
!  return electron temperature
        tem=tea+(teped-tea)*(efac**2-1.)/ &
                (efac**2+1.)
        if (tpoe.gt.1) tem=tem+ten*((1.+xps)**tpoe-(1.+tpoe*xps))
      else if (i.eq.1) then
!  return the derivative
        tem=-4.*(teedg/umax)*(teped-tea)/(efac+1./efac)**2
        if (tpoe.gt.1) tem=tem-(tpoe/umax)*ten*((1.+xps)**(tpoe-1.)-1.)
      else
!  return the second derivative
        tem=-2.*(teedg/umax)**2*(teped-tea)*  &
                ((efac**2-1.)/(efac**2+1.))/ &
                (0.5*(efac+1./efac))**2
        if (tpoe.gt.1.)    &
            tem=tem+(ten*tpoe*(tpoe-1.)/umax**2)*(1.+xps)**(tpoe-2.)
      end if
      tempe=tem
   end function tempe
!
!**********************************************************************
!
      function tempi(psi,id,i)
!     ************************
!
!  Function defines the ion temperature profile as a function
!  of the flux.
!  I=0 gives the temperature and I=1 is the derivative w.r.t. Psi
!  ID labels the impurity species (if IMP=1) with ID=1 labeling the
!  main ion species and others in the order in which they appear
!  in the input data file (after FINI command)
!
      use param
      implicit none
!
      integer i,id
      double precision tempe
      double precision tempi,psi,xps,tem
      double precision t0,ta,tped,tpow,tedg,tn
!
      xps=1.-psi/umax
      if (xps.lt.0.) then
        if (xps.gt.-1.0d-7) then
          xps=0.
        else
          write(6,*)' xps<0 in tempi, psi=',psi,' umax=',umax
          stop
        end if
      end if
      if (xps.gt.1.) then
        if (xps.lt.1.000001) then
          xps=1.
        else
          write(6,*)' xps>1 in tempi, psi=',psi,' umax=',umax
          stop
        end if
      end if
      ta=zta(id)
      t0=zt0(id)
      tpow=ztpow(id)
      if ((ipswtch.eq.1)) then
        if (i.eq.0) then
!  return ion temperature
          tem=ta+(t0-ta)*xps**tpow
        else if (i.eq.1) then
!  return the derivative
          tem=-((t0-ta)*tpow/umax)*xps**(tpow-1.)
        else
!  return the second derivative
          tem=((t0-ta)*tpow*(tpow-1.)/umax**2)*xps**(tpow-2.)
        end if
        tempi=tem
        return
      end if
      tped=ztped(id)
      tedg=ztedg(id)
      if (ipswtch.eq.2) then
        if (tpow.le.1.) then
          tn=0.
        else
          tn=(t0-ta-(tped-ta)*(exp(tedg)-exp(-tedg))/ &
                              (exp(tedg)+exp(-tedg)))
        end if
        if (i.eq.0) then
!  return ion temperature
          tem=ta+(tped-ta)*(exp(tedg*xps)-exp(-tedg*xps))/ &
                  (exp(tedg*xps)+exp(-tedg*xps)) &
                  +tn*xps**tpow
        else if (i.eq.1) then
!  return the derivative
          tem=-4.*(tedg/umax)*(tped-ta)/(exp(tedg*xps)+exp(-tedg*xps))**2 &
                 -(tpow/umax)*tn*xps**(tpow-1.)
        else
!  return the second derivative
!  put a small regularizing factor in to stop divergences at x=0
          tem=-8.*(tedg/umax)**2*(tped-ta)*  &
                  ((exp(tedg*xps)-exp(-tedg*xps))/ &
                  (exp(tedg*xps)+exp(-tedg*xps))**3) &
                 +(tn*tpow*(tpow-1.)/umax**2)*(xps+1.0d-10)**(tpow-2.)
        end if
        tempi=tem
        return
      end if
      if ((ipswtch.eq.3).or.(ipswtch.eq.5).or.(ipswtch.eq.6).or.(ipswtch.eq.8)) then
!  Use scaled electron temperature profiles
        tem=t0*tempe(psi,i)/te0
        tempi=tem
        return
      end if
      if (tpow.le.1.) then
        tn=0.
      else
        tn=(t0-ta-(tped-ta)*(exp(tedg)-exp(-tedg))/ &
                              (exp(tedg)+exp(-tedg))) &
            /(2**tpow-1.-tpow)
      end if
      if (i.eq.0) then
         !  return ion temperature
        tem=ta+(tped-ta)*(exp(tedg*xps)-exp(-tedg*xps))/ &
                (exp(tedg*xps)+exp(-tedg*xps)) &
                +tn*((1.+xps)**tpow-(1.+tpow*xps))
      else if (i.eq.1) then
!  return the derivative
        tem=-4.*(tedg/umax)*(tped-ta)/(exp(tedg*xps)+exp(-tedg*xps))**2 &
               -(tpow/umax)*tn*((1.+xps)**(tpow-1.)-1.)
      else
         !  return the second derivative
        tem=-8.*(tedg/umax)**2*(tped-ta)*  &
                ((exp(tedg*xps)-exp(-tedg*xps))/ &
                (exp(tedg*xps)+exp(-tedg*xps))**3) &
               +(tn*tpow*(tpow-1.)/umax**2)*(1.+xps)**(tpow-2.)
      end if
      tempi=tem
!!$      if (psi.lt.0.) then
!!$	write(nw,*)'error***negative psi in tempi'
!!$	write(nw,*)'psi=',psi
!!$	stop
!!$      end if
!!$      psih=psi/umax
!!$      if (imp.eq.0) id=1
!!$      if (id.eq.1) then
!!$	t0=t0i
!!$	tpow=tpoi
!!$      else
!!$	t0=zt(id)
!!$	tpow=zat(id)
!!$      end if
!!$      if (i.eq.0) then
!!$        if (id.eq.1) then
!!$	  tempi=t0i*((1.+ati)*psih**tpoi-ati*psih**tpoi1)
!!$        else
!!$	  tempi=t0*(psi/umax)**tpow
!!$        end if
!!$      else
!!$        if (id.eq.1) then
!!$	  tempi=(t0i/umax)*(tpoi*(1.+ati)*psih**(tpoi-1.)
!!$     &        -tpoi1*ati*psih**(tpoi1-1.))
!!$        else
!!$	  tempi=tpow*(t0/umax)*(psi/umax)**(tpow-1.)
!!$        end if
!!$      end if
   end function tempi
!
!**********************************************************************





   recursive function elong(con, id) result(kap)

     ! Calculate elongation (or derivative with psi if id=1)

     use param
     implicit none

     integer :: con, id

     integer :: i
     double precision :: rmin, rmax,zmax
     double precision :: rat, kap

     double precision :: elong_l, elong_u, psi_l, psi_u, elongp, elongf

     rmin = 1.e6
     rmax = 0.
     zmax = 0.
     !find max r,z and min r
     do i =1,npts

        if (rpts(con,i) .lt. rmin) then
           rmin = rpts(con,i)

        else if (rpts(con,i) .gt. rmax) then
           rmax = rpts(con,i)
        end if

        if (zpts(con,i) .gt. zmax) then
           zmax = zpts(con,i)
        end if


     end do



     !Elongation
     if (id .eq. 0) then
        elongf = 2.* zmax/(rmax - rmin)

        if (con.eq.ncon) then
           !linear extrapolate to innermost flux surface

           elong_u = elong(con-1,0)
           elongp = elong(con,1)

           elongf = elong_u + elongp*(psiv(con) - psiv(con-1))


        end if

     !Derivative
     else if (id .eq. 1) then

        !forward difference if first contour
        if (con .eq. 1) then

           elong_u = elong(con+1,0)
           elong_l = elong(con,0)

           psi_u = psiv(con+1)
           psi_l = psiv(con)
           elongf = (elong_u - elong_l)/(psi_u - psi_l)

        ! backwards difference if second last contour
        else if (con .eq. ncon-1) then

           elong_u = elong(con, 0)
           elong_l = elong(con-1,0)

           psi_u = psiv(con)
           psi_l = psiv(con-1)
           elongf = (elong_u - elong_l)/(psi_u - psi_l)



        !Extrapolate second to last derivate to get last flux surface
        else if (con .eq. ncon) then

           !Notice using first derivative
           elong_u = elong(con-1, 1)
           elong_l = elong(con-2,1)

           psi_u = psiv(con-1)
           psi_l = psiv(con-2)

           !linear extrapolation
           rat = (elong_u - elong_l)/(psi_u - psi_l)

           elongf = elong(con-1,1) + rat*(psiv(con)-psiv(con-1))

        !centred difference
        else

           elong_l = elong(con-1,0)
           elong_u = elong(con+1,0)

           psi_l = psiv(con-1)
           psi_u = psiv(con+1)
           elongf = (elong_u - elong_l)/(psi_u - psi_l)
        end if


     end if

     kap = elongf
   end function elong





   recursive function shift(con, id) result(del)

     !Calculates shift from r0 by looking at average r

     use param
     implicit none

     integer :: id, con
     double precision :: shiftf, rat, del
     double precision :: shift_l, shift_u, psi_u,psi_l

     if (id .eq. 0) then

        !Average r value
        !shift = (sum(rpts(con, :))/max(1, size(rpts(con,:)))) - r0
        shiftf = (maxval(rpts(con,:)) + minval(rpts(con,:)))/2 - rcen

        !If last fluxsurface, extrapolate
        if (con .eq. ncon) then
           shiftf = shift(con-1,0) + shift(con-1,1)*(psiv(con) -psiv(con-1))

        end if

     else if (id .eq. 1) then

        !forward difference for first point
        if (con .eq. 1) then

           shift_l =  shift(con,0)
           shift_u =  shift(con+1,0)

           psi_l = psiv(con)
           psi_u = psiv(con+1)
           shiftf = (shift_u - shift_l)/(psi_u - psi_l)

        !backward different for the last point
        else if (con .eq. ncon-1) then

           shift_l =  shift(con-1,0)
           shift_u =  shift(con,0)

           psi_l = psiv(con-1)
           psi_u = psiv(con)
           shiftf = (shift_u - shift_l)/(psi_u - psi_l)

        !Use same derivative as before (linear extrap)
        else if (con .eq. ncon) then

           shift_l =  shift(con-2,1)
           shift_u =  shift(con-1,1)

           psi_l = psiv(con-2)
           psi_u = psiv(con-1)

           rat = (shift_u-shift_l)/(psi_u -psi_l)

           shiftf = shift_u + rat*(psiv(con)-psiv(con-1))
        ! centred difference
        else

           shift_l =  shift(con-1,0)
           shift_u =  shift(con+1,0)

           psi_l = psiv(con-1)
           psi_u = psiv(con+1)
           shiftf = (shift_u - shift_l)/(psi_u - psi_l)

        end if

     end if

     del = shiftf
   end function shift


   subroutine dpsidrho (dpdr, rho)

     use param
     implicit none

     integer :: i
     double precision, dimension(ncon) :: rho,dpdr
     double precision :: rhomax


     rho =( maxval(rpts, dim=2) - minval(rpts,dim=2)) /2


     rho(ncon) = 0.
     rhomax = maxval(rho)


     !normalise rho

     rho = rho/rhomax
     do i =1,ncon
        !print*, i, rho(i)
        if (i .eq. 1) then

           dpdr(i) = (psiv(i+1)-psiv(i)) / (rho(i+1) - rho(i))

        else if (i .eq. ncon) then

           dpdr(i) = (psiv(i) - psiv(i-1)) / (rho(i) - rho(i-1))

        else

           dpdr(i) = (psiv(i+1) - psiv(i-1)) / (rho(i+1) - rho(i-1))

        end if


     end do
   
   end subroutine dpsidrho


   ! Volume of each flux surface
   subroutine dVdrho(vols, areas, volsp)

     use param
     implicit none


     double precision, dimension(ncon)  :: voldiff, flxvol, areacon, volcon, flxarea
     integer ::i, j,k


     double precision :: flxvol_l, flxvol_u, psi_u, psi_l

     double precision :: r_l, r_r, z_l, z_r, shift
     real, dimension(ncon) :: vols, areas, volsp


     flxarea = 0.
     ! for each contour calculate area with trapeze rule
     do i= ncon,1,-1

        do j = 1, npts-1

           r_l = rpts(i,j)
           r_r = rpts(i,j+1)

           z_l = zpts(i,j)
           z_r = zpts(i,j+1)

           !write(nw,*) 'Co-ordinates (',i,',',j,')',  r_l, r_r, z_l, z_r
           flxarea(i) = flxarea(i) + (abs(z_l+z_r) * abs(r_l-r_r) /2)
           !write(nw,*) flxarea
        end do

        !area time 2pi * R  (R = r0 + shift)
        flxvol(i) = flxarea(i) * 2. * pi * (shift(i,0) + rcen)

        !Area/Vol of each flux surface from i-1 to i
        if (i .eq.ncon) then
           areacon(i) = flxarea(i)
           volcon(i) = flxvol(i)

        else
           areacon(i) = flxarea(i) - flxarea(i+1)
           volcon(i) = flxvol(i) - flxvol(i+1)

        end if

       ! write(nw,*) 'vol of flux surface ', i, flxvol(i)

     end do

     vol = sum(volcon)

     ! Calculates gradient dV/drho
     do k =1,ncon
        if (k .eq. 1) then

           flxvol_l = flxvol(k)
           flxvol_u = flxvol(k+1)

           psi_l = psiv(k)
           psi_u = psiv(k+1)

        else if (k .eq. ncon) then

           flxvol_l = flxvol(k-1)
           flxvol_u = flxvol(k)

           psi_l = psiv(k-1)
           psi_u = psiv(k)

        else

           flxvol_l = flxvol(k-1)
           flxvol_u = flxvol(k+1)

           psi_l = psiv(k-1)
           psi_u = psiv(k+1)
        end if


        voldiff(k) = (flxvol_u - flxvol_l)/(psi_u - psi_l)

     end do

     !change vol/area to between i-1/2 to i+1/2 and flips array

     do i=2,ncon

        vols(ncon-i+1) = sngl(volcon(i)+volcon(i-1))/2.
        areas(ncon-i+1) = sngl(areacon(i)+areacon(i-1))/2.
        volsp(ncon-i+1) = sngl(voldiff(i)+voldiff(i-1))/2.

     end do

     print*,'Volcon sum is ',sum(volcon)
     
     volsp(ncon) = real(voldiff(1))
     vols(ncon) = real(volcon(1)/2.)
     areas(ncon) = real(areas(1)/2.)

     print*,'Vols sum is ',sum(vols)

   end subroutine dVdrho


   subroutine lam(lambdas, beam, ecomp)
     !From 'Neutral beam heating applications and development'
     !M M Menon, Oak ridge national labs
     ! E_beam must be in keV
     use param
     implicit none

     double precision, dimension(ncon) :: lambdas

     double precision :: psi, ne, dense, rat
     integer :: i, beam, ecomp

     do i=2,ncon
        psi = psiv(i)

        ne = dense(psi,0)

        lambdas(i) = 2.8e17 * E_b(beam)/(ecomp *ne)

     end do

     rat = (lambdas(2) - lambdas(3))/(psiv(2) - psiv(3))

     lambdas(1) = lambdas(2) + rat*(psiv(1) - psiv(2))

   end subroutine lam


   subroutine rhotor(phi)

     use param
     implicit none

     integer :: con
     double precision, dimension(ncon) :: phi

     phi=0.0
     do con=ncon,1,-1

        if (con.eq.ncon) then
           phi(con) = sfac(con)*psiv(con)
        else
           phi(con) = phi(con+1) + sfac(con)*(psiv(con)-psiv(con+1))


        end if
     end do

   end subroutine rhotor
