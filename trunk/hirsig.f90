module hirsig_mod
  implicit none

  private
  public :: conlen, erfun, hirsig

contains
      subroutine coltim(k)
!     *****************
!
! calculates the collision times between the particle species.
!
      use param
      use coldat
      implicit none
!
      double precision coolog,vthi,vthj,rat,chdn,fac
      integer i,j,k
!
      coolog=20.
      do i=1,nimp+2
        do j=1,nimp+2
          if (abs(ts(i)).lt.1e-8) then
!  If zero edge temperature, bootstrap current will be zero, so input
!  dummy variables
            xab(i,j)=1.
            coltau(i,j)=0.
            return
          end if
          vthi=sqrt(2.*bk*ts(i)/sm(i))
          vthj=sqrt(2.*bk*ts(j)/sm(j))
          xab(i,j)=vthj/vthi
          rat=4.*pi*eps0*4.*pi*eps0/ch(i)
          chdn=dn(j)*ch(i)
          fac=sm(i)/ch(j)
          coltau(i,j)=(3./(16.*sqrt(pi)))*(fac**2)*(vthi**3)   &
          *rat/(chdn*coolog)
          if (i.eq.1) then
            if (j.eq.1) nuee(k)=1./coltau(i,j)
            if (j.eq.2) nuei(k)=1./coltau(i,j)
            if (k.eq.ncon) then
             write(6,*)' j=',j,' coltau(i,j)=',coltau(i,j)
            end if
          end if
        end do
      end do
   end subroutine coltim
!
!******************************************************************
!
      subroutine conlen(k,fsi,conl,dotav)
!     ***********************************
!
!  calculates the connection length (Rq in large aspect ratio)
!
      use param
      use scene_errors
      implicit none

      integer k,kmax,i,iskip,ip,ks
      double precision fsi,conl,dotav,dl
      double precision t1,t2,tk
      double precision bigth(npts),dbdl(npts)
      double precision bsq,b1,b2,b3,dbint
      double precision sumk,rk1,rk2
      double precision errcon,errgot
!
      kmax=100
      errcon=1.0d-4
      bigth(1)=0.
      dl=circumf(k)/(npts-1.)
      t1=sqrt((fsi/(rpts(k,1)*bppts(k,1)))**2+1)
      b1=sqrt((fsi/rpts(k,npts-1))**2+bppts(k,npts-1)**2)
      b2=sqrt((fsi/rpts(k,1))**2+bppts(k,1)**2)
      b3=sqrt((fsi/rpts(k,2))**2+bppts(k,2)**2)
      dbdl(1)=(b3-b1)/(2.*dl)
      dbint=0.
      do i=2,npts
        ip=i+1
        if (ip.gt.npts) ip=2
        t2=sqrt((fsi/(rpts(k,i)*bppts(k,i)))**2+1.)
        bigth(i)=bigth(i-1)+0.5*(t1+t2)*dl
        t1=t2
        b1=b2
        b2=b3
        b3=sqrt((fsi/rpts(k,ip))**2+bppts(k,ip)**2)
        dbdl(i)=(b3-b1)/(2.*dl)
        bsq=(fsi/rpts(k,i))**2+bppts(k,i)**2
        dbint=dbint+dbdl(i)**2*bppts(k,i)*dl/bsq
      end do
      dotav=dbint/rnorm(k)
      do i=1,npts
        bigth(i)=2.*pi*bigth(i)/bigth(npts)
      end do
      iskip=0
      sumk=0.
      do 10 ks=1,kmax
        if (iskip.eq.1) goto 10
        rk1=0.
        rk2=0.
        do i=1,npts-1
          bsq=(fsi/rpts(k,i))**2+bppts(k,i)**2
          rk1=rk1+sin(ks*bigth(i))*dbdl(i)*dl/bsq
          rk2=rk2+sin(ks*bigth(i))*dbdl(i)*dl/(bsq*sqrt(bsq))
        end do
        tk=rk1*rk2/ks
        if (ks.gt.1) then
          if (abs(sumk).lt.1e-8) then
             write(6,*)' ERROR***sumk=0 in conlen'
             call error_msg('ERROR***sumk=0 in conlen', 1)
            stop
          end if
          errgot=abs(tk/sumk)
          if (errgot.lt.errcon) iskip=1
        end if
        sumk=sumk+tk
        if (ks.eq.kmax) then
         write(6,*)' error***plasma bdry too complictaed for'
         write(6,*)' connection length calculation'
         stop
        end if
 10   continue
      conl=sumk*bsqav(k)**2/(pi*dbint)
    end subroutine conlen
!
!******************************************************************
!
      subroutine erfun
!     ****************
!
!   calculates and stores error function
!
      implicit none
      double precision dys,ymax,phi
      double precision y,fac,pi
      integer ny,i
      common/phifun/dys,ymax,phi(2000)
!
      pi=4.*atan(1.)
      ymax=3.
      ny=2000
      dys=ymax/(ny-1.)
      fac=2./sqrt(pi)
      phi(1)=0.
      do i=1,ny-1
        y=(i-0.5)*dys
        phi(i+1)=phi(i)+fac*exp(-y*y)*dys
      end do
  end subroutine erfun
!
!******************************************************************
!!$c
!!$      subroutine foot(u,ixout)
!!$c     ************************
!!$c
!!$c   footprint for the ecrh deposition in psi space
!!$c
!!$      common/footp/rf,zf,drf,dzf,psires,sprd,tec,iso
!!$      common/mesh/r(400),z(400),dr,dz,nr,nz,nsym
!!$      dimension u(400,400),ixout(400,400)
!!$c
!!$      ir=0
!!$      do 10 i=1,nr
!!$        if (ixout(i,nsym).le.0) goto 10
!!$        if (r(i).gt.rf) goto 10
!!$        ir=i
!!$ 10   continue
!!$      if ((ir.eq.0).or.(ir.eq.nr)) then
!!$        write(6,*)' error****ecrh resonance outside r of plasma'
!!$        stop
!!$      end if
!!$      iz=0
!!$      do 20 i=1,nz
!!$        if (z(i).gt.zf) goto 20
!!$        iz=i
!!$ 20   continue
!!$      if ((iz.eq.0).or.(iz.eq.nz)) then
!!$        write(6,*)' error****ecrh resonance outside z of plasma'
!!$        stop
!!$      end if
!!$      if (ir.gt.0) then
!!$        if (ixout(ir,iz).le.0) then
!!$          write(6,*)' error***ecrh resonance outside plasma'
!!$          stop
!!$        end if
!!$      end if
!!$      psires=u(ir,iz)
!!$      psimin=psires
!!$      psimax=psires
!!$      rlo=rf-drf/2.
!!$      rup=rf+drf/2.
!!$      zlo=zf-dzf/2.
!!$      zup=zf+dzf/2.
!!$      irl=0
!!$      do 30 i=1,nr
!!$        if (r(i).gt.rlo) goto 30
!!$        irl=i
!!$ 30   continue
!!$      izl=0
!!$      do 40 i=1,nz
!!$        if (z(i).gt.zlo) goto 40
!!$        izl=i
!!$ 40   continue
!!$      iru=0
!!$      do 50 i=1,nr
!!$        if (r(i).gt.rup) goto 50
!!$        iru=i
!!$ 50   continue
!!$      izu=0
!!$      do 60 i=1,nz
!!$        if (z(i).gt.zup) goto 60
!!$        izu=i
!!$ 60   continue
!!$      if ((irl.eq.0).or.(izl.eq.0)) psimin=0.
!!$      do 80 i=irl,iru
!!$        do 70 j=izl,izu
!!$          if (ixout(i,j).le.0) then
!!$            psimin=0.
!!$            goto 70
!!$          end if
!!$          if (psimin.gt.u(i,j)) psimin=u(i,j)
!!$          if (psimax.lt.u(i,j)) psimax=u(i,j)
!!$ 70     continue
!!$ 80   continue
!!$      sprd=psimax-psimin
!!$      return
!!$      end
!!$c
!******************************************************************
!
      subroutine fric
!     ***************
!
      use param
      use coldat
      implicit none
!
      double precision sum00(nimp+2),sum01(nimp+2),sum10(nimp+2), &
                       sum11(nimp+2)
      double precision rmat00,rmat01,rmat11
      double precision rnat00,rnat01,rnat10,rnat11
      integer i,j
!
      do i=1,nspec
        sum00(i)=0.
        sum10(i)=0.
        sum01(i)=0.
        sum11(i)=0.
        do j=1,nspec
          rmat00=-(1.+sm(i)/sm(j))/(1.+xab(i,j)**2)**1.5
          sum00(i)=sum00(i)+dn(i)*sm(i)*rmat00/coltau(i,j)
          rmat01=-1.5*(1.+sm(i)/sm(j))/(1.+xab(i,j)**2)**2.5
          sum01(i)=sum01(i)+dn(i)*sm(i)*rmat01/coltau(i,j)
          rmat11=-(3.25+4.*xab(i,j)**2+7.5*xab(i,j)**4)        &
              /(1.+xab(i,j)**2)**2.5
          sum11(i)=sum11(i)+dn(i)*sm(i)*rmat11/coltau(i,j)
        end do
        sum10(i)=sum01(i)
      end do
      do i=1,nspec
        do  j=1,nspec
          rnat00=(1.+sm(i)/sm(j))/(1.+xab(i,j)**2)**1.5
          rl11(i,j)=dn(i)*sm(i)*rnat00/coltau(i,j)
          if (i.eq.j) rl11(i,j)=rl11(i,j)+sum00(i)
          rnat10=1.5*(1.+sm(i)/sm(j))/(1.+xab(i,j)**2)**2.5
          rl21(i,j)=dn(i)*sm(i)*rnat10/coltau(i,j)
          if (i.eq.j) rl21(i,j)=rl21(i,j)+sum10(i)
          rnat01=(ts(i)/(ts(j)*xab(i,j)))*1.5*            &
                 (1.+sm(j)/sm(i))/(1.+xab(j,i)**2)**2.5
          rl12(i,j)=dn(i)*sm(i)*rnat01/coltau(i,j)
          if (i.eq.j) rl12(i,j)=rl12(i,j)+sum01(i)
          rnat11=6.75*(ts(i)/ts(j))*xab(i,j)**2/     &
                 (1.+xab(i,j)**2)**2.5
          rl22(i,j)=dn(i)*sm(i)*rnat11/coltau(i,j)
          if (i.eq.j) rl22(i,j)=rl22(i,j)+sum11(i)
        end do
      end do
  end subroutine fric
!
!******************************************************************
!
      function gov(y,id)
!     ******************
!
!  Calculates Chandrasakhar function (id=1) or error function (id=0)
!
      implicit none
      double precision dys,ymax,phi
      double precision gov
      double precision y,ylo,yup,phioy,rat,pi,fac,phid,gg
      integer id,k
      common/phifun/dys,ymax,phi(2000)
!
      pi=4.*atan(1.)
      if (y.lt.0.) then
        write(6,*)'error***only positive arguments to error fn, y=',y
        stop
      end if
      k=int((y+dys)/dys)
      if (k.eq.0) k=1
      ylo=(k-1)*dys
      yup=k*dys
      if (yup.gt.ymax) then
        phioy=1.-exp(-y*y)/(y*sqrt(pi))
      else
        rat=(y-ylo)/(yup-ylo)
        phioy=phi(k)+rat*(phi(k+1)-phi(k))
      end if
      if (id.eq.0) then
        gov=phioy
        return
      end if
      fac=2./sqrt(pi)
      phid=fac*exp(-y*y)
      if (y.le.1.0d-3) then
        gg=fac*y/3.
      else
        gg=(phioy-y*phid)/(2.*y*y)
      end if
      gov=gg
  end function gov
!
!******************************************************************
!
      subroutine hirsig(psi,fc,fsi,conl,dotav,bsq,bstrap,aspin,k)
!     ***************************************************************
!
!  calculates the bootstrap current using the full hirshman-sigmar
!  formalism (nf 21 (1981) 1079)
!
      use param
      use coldat
      use profiles_mod, only : tempe, dense, tempi, densi
      implicit none

      double precision psi,fc,fsi,conl,dotav,bsq,bstrap,aspin
      double precision ft
      double precision pda,v1a,v2a
      double precision, dimension(:), allocatable:: tsd,dnd
      double precision vvec(2*nimp+4),upar(nimp+2)
      double precision tmat(2*nimp+4,2*nimp+4),rlmat(2*nimp+4,2*nimp+4)
      integer i,j,id,narr,k
!
      allocate(ch(nimp+2),ts(nimp+2),tsd(nimp+2),sm(nimp+2),  &
               dn(nimp+2),dnd(nimp+2) )
      allocate(coltau(nimp+2,nimp+2),xab(nimp+2,nimp+2),vis(nimp+2,3) )
      allocate( rl11(nimp+2,nimp+2),rl12(nimp+2,nimp+2),   &
               rl21(nimp+2,nimp+2),rl22(nimp+2,nimp+2) )


!
!  load up species properties into unifying arrays (e=1, i=2, i=
!  3,4,...nimp+2, in same order as impure arrays above
      nspec=nimp+2
      mp=1.6726e-27
      ch(1)=-eq
      ts(1)=tempe(psi,0)
      tsd(1)=tempe(psi,1)
      sm(1)=me
      dn(1)=dense(psi,0)
      dnd(1)=dense(psi,1)
      do i=1,nimp+1
        ch(i+1)=iz(i)*eq
        ts(i+1)=tempi(psi,i,0)
        tsd(i+1)=tempi(psi,i,1)
        dn(i+1)=densi(psi,i,0)
        dnd(i+1)=densi(psi,i,1)
        sm(i+1)=zmas(i)*mp
      end do
!  load up collision time mesh (coltau) and thermal velocity ratio
!  mesh xab(i,j)=vthj/vthi
      call coltim(k)
!  trapped particle fraction
      ft=1.-fc
!  calculate viscosity coefficients...
      call viscos(ft,conl,dotav,bsq,aspin,k)
!  calculate friction coefficients
      call fric
!  load up matrices
      do i=1,2*nspec
        id=i
        if (i.le.nspec) then
          pda=bk*(ts(i)*dnd(i)+tsd(i)*dn(i))
          v1a=-fsi*pda/(ch(i)*dn(i))
          v2a=-fsi*bk*tsd(i)/ch(i)
          vvec(i)=vis(i,1)*v1a+vis(i,2)*v2a
        else
          id=i-nspec
          pda=bk*(ts(id)*dnd(id)+tsd(id)*dn(id))
          v1a=-fsi*pda/(ch(id)*dn(id))
          v2a=-fsi*bk*tsd(id)/ch(id)
          vvec(i)=vis(id,2)*v1a+vis(id,3)*v2a
        end if
        id=i-nspec
!!$        if (iso.eq.1) then
!!$          if (i.eq.1) then
!!$!  if iso=1, then use the anisotropic electron pressure model...
!!$            gau=tec*exp(-((psi-psires)/sprd)**2)
!!$            tpar=ts(1)*(1.-2.*gau/3.)
!!$            tperp=ts(1)*(1.+gau/3.)
!!$            tdpar=tsd(1)*(1.-2.*gau/3.)
!!$     &      +4.*ts(1)*((psi-psires)/sprd**2)*gau/3.
!!$            tdperp=tsd(1)*(1.+gau/3.)
!!$     &      -2.*ts(1)*((psi-psires)/sprd**2)*gau/3.
!!$            pda=bk*(tperp*dnd(1)+tdperp*dn(1))
!!$            pa=bk*tperp*dn(1)
!!$            pda=pda+0.5*(tdperp/tperp-tdpar/tpar)*pa
!!$            v1a=fsi*pda/(ch(1)*dn(1))
!!$            v2a=fsi*bk*tdperp/ch(1)
!!$          end if
!!$          if (i.eq.1) then
!!$            vvec(i)=vis(i,1)*v1a+vis(i,2)*v2a
!!$          else if (i.eq.(nspec+1)) then
!!$            vvec(i)=vis(id,2)*v1a+vis(id,3)*v2a
!!$          end if
!!$          vvec(i)=sqrt(tperp/tpar)*vvec(i)
!!$        end if
        do j=1,2*nspec
          tmat(i,j)=0.
          if ((i.le.nspec).and.(j.le.nspec)) then
            rlmat(i,j)=-rl11(i,j)
            if (i.eq.j) tmat(i,j)=vis(i,1)
          else if ((i.gt.nspec).and.(j.le.nspec)) then
            rlmat(i,j)=rl21(i-nspec,j)
            if ((i-nspec).eq.j) tmat(i,j)=vis(i-nspec,2)
          else if ((i.le.nspec).and.(j.gt.nspec)) then
            rlmat(i,j)=rl12(i,j-nspec)
            if (i.eq.(j-nspec)) tmat(i,j)=vis(i,2)
          else
            rlmat(i,j)=-rl22(i-nspec,j-nspec)
            if (i.eq.j) tmat(i,j)=vis(i-nspec,3)
          end if
        end do
      end do
      narr=2*nspec
      call matinv(rlmat,tmat,vvec,upar,narr)
      bstrap=0.
      do 40 i=1,nspec
        bstrap=bstrap+ch(i)*dn(i)*upar(i)
 40   continue
      bstrap=bstrap/sqrt(bsq)
      deallocate(ch,ts,tsd,sm,  &
               dn,dnd )
      deallocate(coltau,xab,vis )
      deallocate(rl11,rl12,rl21,rl22)
    end subroutine hirsig
!
!******************************************************************
!
      subroutine matinv(rlmat,tmat,vvec,upar,narr)
!     ********************************************
!
      use param
      implicit none
      integer i,j,narr,ndim, ifail
      double precision aa(2*nimp+4,2*nimp+4),bb(2*nimp+4),cc(2*nimp+4)
      double precision rlmat(2*nimp+4,2*nimp+4),tmat(2*nimp+4,2*nimp+4),  &
                       vvec(2*nimp+4),upar(nimp+2)
      integer, allocatable:: IPIV(:)
!
      ndim=2*nimp+4
      do i=1,narr
        bb(i)=vvec(i)
        do j=1,narr
          aa(i,j)=rlmat(i,j)+tmat(i,j)
        end do
      end do
      allocate(IPIV(ndim))
      ifail=0
!      call f04arf(aa,ndim,bb,narr,cc,wkspce,ifail)
      cc(:) = bb(:)
!      write(6,*) 'calling DGESV'
      call DGESV( ndim, 1, aa, ndim, IPIV, cc, ndim, ifail )
      if (ifail.ne.0) then
       write(6,*)' error detected in lapack routine DGESV, ifail='  &
          ,ifail
       stop
      end if
      do i=1,narr/2
        upar(i)=cc(i)
      end do
   end subroutine matinv
!
!******************************************************************
!
      function phig(y)
!     ****************
!
      implicit none
      double precision dys,ymax,phi
      double precision y,phig,ylo,yup,fac,pi,phid,rat,gg
      integer k
      common/phifun/dys,ymax,phi(2000)
!
      pi=4.*atan(1.)
      if (y.lt.0.) then
        write(6,*)'error***only positive arguments to error fn, y=',y
        stop
      end if
      k=int((y+dys)/dys)
      if (k.eq.0) k=1
      ylo=(k-1)*dys
      yup=k*dys
      if (yup.gt.ymax) then
        phig=1.-exp(-y*y)/(y*sqrt(pi))
      else
        rat=(y-ylo)/(yup-ylo)
        phig=phi(k)+rat*(phi(k+1)-phi(k))
      end if
      fac=2./sqrt(pi)
      phid=fac*exp(-y*y)
      if (y.le.1.0e-3) then
        gg=fac*y/3.
      else
        gg=(phig-y*phid)/(2.*y*y)
      end if
      phig=phig-gg
   end function phig
!
!*****************************************************************
!
      subroutine viscos(ft,conl,dotav,bsq,aspin,k)
!     ******************************************
!
      use param
      use coldat
      implicit none
!
      double precision vd(nimp+2,500),vt(nimp+2,500),vtot(nimp+2,500)
      double precision rk11(nimp+2),rk12(nimp+2),rk22(nimp+2)
      double precision fc,vmax,dv
      double precision vtha,omta,vsttau,fac1,fac2
      double precision ft,conl,dotav,bsq,aspin
      double precision y,vv,vtij,vdij,arg
      integer nv,i,n,j,k
!

!
      fc=1.-ft
      vmax=5.
!  calculate 90-degree scattering freq as a function of v/vthi
      nv=500
      dv=vmax/(nv-1.)
      do i=1,nspec
!  collisionality * tau_ii...
        vtha=sqrt(2.*bk*ts(i)/sm(i))
        omta=vtha/conl
        vsttau=ft*omta*bsq/(dotav*vtha**2)
        vsttau=vsttau*8./(3.*pi)
        do n=1,nv-1
          vv=n*dv-dv/2.
          vd(i,n)=0.
          vt(i,n)=0.
          do j=1,nspec
            y=vv/xab(i,j)
!            if (y.lt.0.)
            vdij=phig(y)/(coltau(i,j)*vv**3)
            vd(i,n)=vd(i,n)+vdij
            vtij=(((gov(y,0)-3.*gov(y,1))/vv**3)+4.*(ts(i)/ts(j)  &
                 +(1./xab(i,j)**2))*gov(y,1)/vv)/coltau(i,j)
            vt(i,n)=vt(i,n)+vtij
          end do
          vd(i,n)=(3.*sqrt(pi)/4.)*vd(i,n)
          vt(i,n)=(3.*sqrt(pi)/4.)*vt(i,n)
          if (nco.le.1) then
!  collisional suppression if nco=1
            fac1=1.+nco*vsttau*vd(i,n)/vv
            fac2=1.+nco*5.*pi*vt(i,n)/(8.*vv*omta)
          else
!  shaing's result that a=1 experiences no collisional suppression
            fac1=1.+(1.-aspin**2)*vsttau*vd(i,n)/vv
            fac2=1.+5.*pi*vt(i,n)/(8.*vv*omta)
          end if
          vtot(i,n)=vd(i,n)/(fac1*fac2)
        end do
      end do
!  perform velocity integrals to calculate k's
      do i=1,nspec
        rk11(i)=0.
        rk12(i)=0.
        rk22(i)=0.
        do n=1,nv-1
          vv=n*dv-dv/2.
          arg=(vv**4)*exp(-vv**2)*vtot(i,n)*dv
          rk11(i)=rk11(i)+arg
          rk12(i)=rk12(i)+arg*vv**2
          rk22(i)=rk22(i)+arg*vv**4
        end do
        rk11(i)=ft*rk11(i)*8./(3.*fc*sqrt(pi))
        rk12(i)=ft*rk12(i)*8./(3.*fc*sqrt(pi))
        rk22(i)=ft*rk22(i)*8./(3.*fc*sqrt(pi))
        vis(i,1)=rk11(i)*dn(i)*sm(i)
        vis(i,2)=(rk12(i)-2.5*rk11(i))*dn(i)*sm(i)
        vis(i,3)=(rk22(i)-5.*rk12(i)+6.25*rk11(i))*dn(i)*sm(i)
        if (i.eq.1) then
          mu1(k)=rk11(i)
          mu2(k)=(rk12(i)-2.5*rk11(i))
          mu3(k)=(rk22(i)-5.*rk12(i)+6.25*rk11(i))
        end if
      end do
!!$      i=1
!!$      write(6,*)' mu1=',vis(1,1)*(fc/ft)/(dn(i)*sm(i)),  &
!!$                ' mu2=',vis(1,2)*(fc/ft)/(dn(i)*sm(i)),  &
!!$                ' mu3=',vis(1,3)*(fc/ft)/(dn(i)*sm(i))
  end subroutine viscos

end module hirsig_mod
