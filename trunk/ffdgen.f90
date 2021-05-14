      subroutine ffdgen(i,icur,errcur)
!     ********************************
!
!  If i equals zero ffdgen calculates the new value of
!  ff' required and stores it in gst. If i=1 then ffdgen
!  calculates the maximum error between the input and
!  output ff' values.
!
!!$      common/arpram/emin,emax
!!$      common/avges/sfac(50),bsqav(50),rsqav(50),rinv(50),
!!$     &rsqinv(50),rav(50),rnorm(50),ftrap(50),bsj(50),tnue(50),tnui(50)
!!$     &,cnue(50),cnui(50),sighhb(50)
!!$      common/consts/pi,rmu0,bk,eps0
!!$      common/totcur/totbs,totps,totdi,totex,totgs,vloop,vnobs,vspit
!!$      common/consis/ipass,gst(400),errcur
!!$      common/flxvar/ncon,eps(50),sist(50),flxr(50,400),
!!$     &flxz(50,400),icon(50)
!!$      common/oldmsh/epold(50),siold(50)
!!$      common/curpar/s,const,r0,ppow,fpow,umax,bpol,neo,cur,zm,zeffav
!!$     &,itot,powj,scl,nco,p0,bpfac,zmai,aj
!!$      common/mesh/r(400),z(400),dr,dz,nr,nz,nsym
!!$      dimension ixout(400,400),gstk(50)
      use param
      implicit none
      integer k,i,icur,kmax
      double precision psi,eps
      double precision gmax,gold,errcur,gtst,err,errel,err1
      double precision fprof,press
      double precision extapp,extapp2,bavg
      double precision f,ffd,fd,pd,jedge,fval
      double precision fcen,bavcen,t3
!
      gmax=0.
      do k=1,ncon
        psi=psiv(k)
        gold=fprof(psi,1)
        if (abs(gold).gt.gmax) gmax=abs(gold)
      end do
      errcur=0.
      errel=0.
      do 70 k=1,ncon-1
        psi=psiv(k)
        eps=epsv(k)
        call extj(eps,psi,extapp,extapp2,icur)
        extapp=extapp*vloop+extapp2

        f=fprof(psi,2)
        ffd=fprof(psi,1)
        fd=ffd/f
        pd=press(psi,1)
        gold=ffd
        bavg=bsqav(k)
!  match flux surface averages (same for both par/tor currents)...

        !Flag on whether to include nb in ffdgen calc
        if (nbi .eq. 2) then
           gtst=-mu0*(f*extapp+f*J_nb(k)+f*bsj(k)/sqrt(bavg)+f*f*pd/bavg)/scl
        else
           gtst=-mu0*(f*extapp+f*bsj(k)/sqrt(bavg)+f*f*pd/bavg)/scl
        end if


        !print*, extapp, j_nb(k), bsj(k)/sqrt(bav), f*pd/bav
        if (k.eq.1) then
        !write(6,*)'********* in ffdgen, k=1'
        !write(6,*)' f=',f,' pd=',pd,' fd=',fd
        !write(6,*)' extapp=',extapp,' bsj=',bsj(1)/sqrt(bav)
        !write(6,*)' gtst*scl/f=',gtst*scl/f,' scl=',scl
        end if
        if (k.gt.1) then
!  calculate error in ff' relative to max(ff')
          err=abs((gtst*scl-gold))/gmax
          if (errcur.lt.err) errcur=err
          err=abs((gtst*scl-gold)/(gtst*scl+gold))
          if (errel.lt.err) then
            errel=err
            kmax=k
          end if
        else
          err1=abs((gtst*scl-gold)/(gtst*scl+gold))
        end if
        if (i.eq.1) goto 70
! if generating new equilibrium update ff'
        gst(k)=gtst
 70   continue
      psi=psiv(1)
      fval=fprof(psi,2)
      pd=press(psi,1)
      fd=fprof(psi,1)/fval
      jedge=-(fval*pd/sqrt(bsqav(1))+fd*sqrt(bsqav(1))/mu0)
      !write(6,*)' New f-dash=',fprof(psi,1)/fprof(psi,2),' new pd=',pd
      !write(6,*)' New fval=',fval,' New scl=',scl
      !write(6,*)' ******* j.B=',jedge
      !write(6,*)' k=',kmax,' max relative error in g=',errel,' g=',gtst*scl
      !write(6,*)' ERROR in edge ffd=',err1
      !write(6,*)' g ror test=',errcur
      if (i.eq.1) return
!  g at magnetic axis
      fcen=fprof(0.0d0,2)
      pd=press(0.0d0,1)
      bavcen=(fcen/r0)**2
      t3=fcen*pd/bavcen
      call extj(0.0d0,0.0d0,extapp,extapp2,icur)
      gst(ncon)=-mu0*fcen*(extapp*vloop+extapp2+t3)/scl
!  set ipass=1 to use mesh ff' in fprof on next iteration
      ipass=1
!  store old values
      do k=1,ncon
        psiold(k)=psiv(k)/umax
        write(77,*)' psi=',psiold(k),' g=',gst(k)
      end do
   end subroutine ffdgen
