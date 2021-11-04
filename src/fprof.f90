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
          im=psin/dpsi
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
       if (ipswtch.eq.6) then

          ffp = g0 * ( af1*psin**1 + af2*psin**2 + &
               af3*psin**3 - (af1+af2+af3)**4
          if (id.eq.1) then
             fprof=ffp
             return
          end if

          fsq = umax*g0* (af1/2.*psin**2 + af2/3.* &
               psin**3 + af3/4.*psin**4 - (af1+af2+af3) &
               /5. *psin**5

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
               4*(af1+af2+af3)*psin**3

          fprof=ffpp
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
!********************************************************************8

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
