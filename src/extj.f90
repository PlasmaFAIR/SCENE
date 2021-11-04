module ext_current_mod
  implicit none
contains
    subroutine extj(eps,psi,extapp,extapp2,icur)
!       *******************************************
!
!  User--supplied function. returns the value of <j.b>/<b**2>
!  in extapp: extapp2 fills in the missing bootstrap current within
!  1-ps=psic: extapp is a user-supplied function which is scaled
!  by Vloop to get the total plasma current correct
!  The profile which you specify will be constant on a flux
!  surface and J is assumed to be along B so that the total
!  current will be J = <J.B> B /<B**2>
!  The current density is parameterised as a power of psi with
!  ffp parameters
!  The scale is set so as to obtain the required total current as set
!  in cur.
!   If neo=1 then no current profile need be specified. The code
!   assumes ohmic current drive with a neoclassical conductivity
!   (as given by hirshman, hawryluk, birge (nf 17 (1977) 611)).
!   if icur=-1 on entry then neoclassical trapping effects are
!   switched off.
! Mod 22/3/01 byr Howard Wilson:
!      extapp2 used in neoclassical option to add aux current drive
!
      use param
      use profiles_mod, only : dense, densi, fprof, press, tempe
      implicit none
      double precision, intent(in) :: eps,psi
      double precision, intent(out) :: extapp, extapp2      
      integer, intent(in) :: icur
!
      integer :: l,ik,k
      double precision :: xmax,fsi,pt,te,ne,zeff,zni
      double precision :: sigsp,coolog
      double precision :: epsfac,ratio,pro,rt,botcol,topcol,colt,vthe,fac
      double precision :: rat,rrinv,b_average,bstrap,hhb,rle,sigfac
      double precision :: psiq,bstrapq,bavq,rsqint,ravint,rla
      double precision :: t1,t2,t3,t4,scale
!
      extapp2=0.
      xmax=0.
      if (abs(neo).eq.1) then
         !  f(psi)
         fsi=fprof(psi,2)
         !  total pressure
         pt=press(psi,0)
         !  temperature
         te=tempe(psi,0)
        te=bk*te
        !  electron density
        ne=dense(psi,0)
        !  zeff...
        zeff=zm
        if (imp.eq.1) then
           if (ne.gt.0.) then
              zeff=0.
              do l=1,nimp+1
                 zni=densi(psi,l,0)
                 zeff=zeff+(zni*iz(l)**2)/ne
              end do
           else
              zni=densi(psi,1,0)
              if (zni.gt.0.) then
                 write(nw,*)'error***problem in extj, ne=0'
                 write(nw,*)'cannot evaluate zeff'
                 write(nw,*)'psi=',psi,' ne=',ne,' ni=',zni,' umax=',umax
              stop
           else
              !  zero density so no current-> can arbitrarily set zeff=zm
              zeff=zm
           end if
        end if
     end if
     !  coulomb logarithm
     if (abs(pt).lt.1e-8) then
        sigsp=0.
     else
        coolog=log(sqrt(ne*1.0d-6)*bk/te)
        coolog=24.-coolog
        !  electron--electron collision time
        epsfac=(4.*pi*eps0/eq)**2
        ratio=te/eq
        pro=ne*eq
        rt=sqrt(me)*sqrt(te)
        botcol=sqrt(2.*pi)*pro*coolog
        topcol=rt*epsfac*ratio
        colt=(3./4.)*topcol/botcol
        vthe=sqrt(2.*te/me)
        !  spitzer conductivity
        fac=1.98*(ne*eq)*(eq/me)
        sigsp=fac*colt
     end if
     !  extrapolate between flux surfaces:
     ik=0
     if (eps.lt.0.) then
        write(nw,*)'error***eps<0'
        write(nw,*)'problem in extj'
        stop
     end if
     do 5 k=1,ncon
        if (eps.ge.epsv(k)) goto 5
        ik=k
5       continue
        if (ik.eq.0) then
           ik=1
           rat=1.
        else if (ik.lt.ncon) then
           rat=(eps-epsv(ik+1))/(epsv(ik)-epsv(ik+1))
        else
           write(nw,*)'error***cannot interpolate a psi value'
           write(nw,*)'problem in extj'
           stop
        end if
        rrinv=rsqinv(ik+1)-rat*(rsqinv(ik+1)-rsqinv(ik))
        ravint=rav(ik+1)-rat*(rav(ik+1)-rav(ik))
        rsqint=rsqav(ik+1)-rat*(rsqav(ik+1)-rsqav(ik))
        b_average=bsqav(ik+1)-rat*(bsqav(ik+1)-bsqav(ik))
        !  hirshman, hawryluk, birge coefficients
        hhb=sighhb(ik+1)-rat*(sighhb(ik+1)-sighhb(ik))
        if (icur.ne.-1) then
           !  neoclassical conductivity (hirshman, hawryluk, birge)
          sigfac=hhb*sigsp/1.98
       else
          !  classical conductivity
          !  impurity correction factor (hirshman, hawryluk, birge)
          rle=3.4*(1.13+zeff)/(zeff*(2.67+zeff))
          sigfac=sigsp*rle/1.98
       end if
       extapp=sigfac*fsi*rrinv/(2.*pi*b_average)
       !  If ibv set to one, then use electric field profile induced by Bv-ramp
       !  Assumes constant Bv (in R); vloop returned is the required Bv-dot
       if (ibv.eq.1) then
!!$           if (rpk.gt.0.) then
!!$           rla=-0.5*(1.-r1**2*rrinv)                                  &
!!$              +(rsqint/4.-2.*rpk*ravint/3.+0.5*rpk**2)/(rpk-rm)**2   &
!!$              -(0.25*r1**4-2.*rpk*r1**3/3.+0.5*(rpk*r1)**2)*rrinv/(rpk-rm)**2
!!$           else
!!$             rla=-1.
!!$           end if
          rla=1.-rrinv*r1**2
          write(6,*)' rla=',rla
          extapp=extapp*rla*pi/rrinv
       end if
       !  aux current drive if neo=1
       rat=psi/umax
       extapp2=1000.*af1*(1.-(rat-af2)**2)**powj+1000.*af0*(1.-rat)**fpow2
       return
    end if
    !  evaluate <b.grad phi>
    fsi=fprof(psi,2)
    !  extrapolate between flux surfaces:
    ik=0
    do 20 k=1,ncon
       if (psi.ge.psiv(k)) goto 20
       ik=k
20     continue
       if (ik.eq.0) then
          ik=1
          rat=1.
       else if (ik.lt.ncon) then
          rat=(psi-psiv(ik+1))/(psiv(ik)-psiv(ik+1))
       else
          write(nw,*)'error***cannot interpolate a psi value'
          write(nw,*)'problem in extj'
          stop
       end if
       rrinv=rsqinv(ik+1)-rat*(rsqinv(ik+1)-rsqinv(ik))
       bstrap=bsj(ik+1)-rat*(bsj(ik+1)-bsj(ik))
       b_average=bsqav(ik+1)-rat*(bsqav(ik+1)-bsqav(ik))
       !  profile specified in terms of psi
       rat=psi/umax
       t1=(1.-rat)**powj
       t2=((1.-rat)**(fpow1))*((psi/umax)**fpow2)
       scale=((fpow1/(fpow1+fpow2))**(fpow1))*((fpow2/(fpow2+fpow1))**fpow2)
       t2=t2/scale
       t3=(((1.-psi/umax)/0.9)**(powj/10))     &
            *(((psi/umax)/0.1)**(0.111*powj/10.))
       t4=((1.-rat)**(fpow3))*((psi/umax)**fpow4)
       scale=((fpow3/(fpow3+fpow4))**(fpow3))*((fpow4/(fpow4+fpow3))**fpow4)
       t4=t4/scale
       extapp=af0*t1+af1*t2+af2*t3+af3*t4
       if (rat.lt.psic) then
          !  ensure total current profile will be flat near to axis
          psiq=umax*psic
          ik=0
          do 30 k=1,ncon
             if (psiq.ge.psiv(k)) goto 30
          ik=k
30        continue
          if (ik.eq.0) then
             ik=1
             rat=1.
          else if (ik.lt.ncon) then
             rat=(psiq-psiv(ik+1))/(psiv(ik)-psiv(ik+1))
          else
          write(nw,*)'error***cannot interpolate a psi value'
          write(nw,*)'problem in extj'
          stop
       end if
       bstrapq=bsj(ik+1)-rat*(bsj(ik+1)-bsj(ik))
       bavq=bsqav(ik+1)-rat*(bsqav(ik+1)-bsqav(ik))
       extapp2=-(bstrap*sqrt(b_average)-bstrapq*sqrt(bavq))/b_average
    else
       extapp2=0.
    end if
  end subroutine extj
end module ext_current_mod
