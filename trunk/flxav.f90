module flux_average
  implicit none
contains
      subroutine flxav
!     ****************
!
!  Calculates the various flux surface averages which are required.
!  Note that this routine also calculates the flux surface average
!  of the bootstrap current; also evaluates G-S current, pressure
!  and temperature on each flux surface.
!
      use hirsig_mod, only : conlen, hirsig
      use hirsh_mod, only : hirsh
      use param
      use profiles_mod, only : dense, densi, fprof, press, tempe, tempi
      use signeo_mod, only : signeo
      implicit none
      integer :: k,i,ip,im,l
      double precision ::  ffp,pd
      double precision ::  ant,bmod,bphi,bsq,bth,dl
      double precision ::  eps,erribm,fsi,psi,pt,rint,rnor,root,rr,zz
      double precision ::  drdl,d2rdl2,dzdl,d2zdl2
      double precision ::  psi1,psi2,psi3,nu,nu1,nu2,ffpp,pdd
      double precision ::  rcinv(npts),sinu(npts),cosu(npts),dpsi1dl(npts)
      double precision :: d2psi1dl2(npts),dnudpsi(npts),d2nudpsi2(npts)
      double precision :: surfa(npts),bpsur(npts)
      double precision :: bsqar(npts),bphiar(npts), rsqar(npts),rrar(npts)
      double precision :: riar(npts), risar(npts),rar(npts),avblir(npts)
      double precision :: sfloc(npts),colop(npts),vpparr(npts)
      double precision :: bmax(ncon),bmin(ncon)
      double precision :: te,ti,tau,coolog,bfac,epsfac,ratio,pro,rt
      double precision :: topcol,botcol,colte,vthe,colti,vthi,zni
      double precision :: binv,dst,bigint,bigint2,rla,btot,rlag,bp
      double precision :: bla(npts),bla2(npts)
      double precision :: fc,x,ne,pe,rj0,bstrap,conl,dotav
      double precision :: fq(ncon),fqd(ncon)
      double precision :: zeff,qq,rnust,ft,sigfac
      double precision :: ant1,ant2,ant3
      double precision :: t1(npts),t2(npts),t3(npts)
!      double precision ppar(ncon),pppar(ncon),ppan(ncon),pppan(ncon)
      integer :: nst
!
      double precision :: q1,q2,q3,x1,x2,x3,av,bv,cv
!      double precision p1,p2,p3
      double precision :: nupplt(npts),chiplt(npts),y1(npts),y2(npts),y3(npts)
      double precision :: y4(npts),y5(npts),y6(npts)
      double precision :: plotx(npts),ploty1(ncon,npts),ploty2(ncon,npts)
      double precision :: ploty3(npts),ploty4(npts)
      double precision :: dbtdrho,dbpdrho

      if (icont.gt.-3) then
        allocate( sfac(ncon),bsqav(ncon), bphiav(ncon),  &
                rsqav(ncon),rinv(ncon),rsqinv(ncon),rav(ncon),rnorm(ncon),  &
                avbli(ncon),qp(ncon),qpp(ncon),cnue(ncon),tnue(ncon), &
                cnui(ncon),tnui(ncon),bdl(ncon),bav(ncon), ftrap(ncon),ftrapd(ncon),  &
                bsj(ncon),sighhb(ncon),mu1(ncon),mu2(ncon),mu3(ncon),      &
                nuee(ncon),nuei(ncon),vpp(ncon)  )
                allocate(ajdotb(ncon),fastp(ncon))
      end if
!
!  calculate variables to be flux surface averaged
!  k=1 corresponds to boundary, k=ncon is the magnetic axis
      do  k=1,ncon-1
        bmax(k)=0.
        psi=psiv(k)
        eps=epsv(k)
        fsi=fprof(psi,2)
        ffp=fprof(psi,1)
        ffpp=fprof(psi,4)
        pt=press(psi,0)
        pd=press(psi,1)
        pdd=press(psi,2)
        dl=circumf(k)/(npts-1)
        do i=2,npts-1
          dl=sqrt((rpts(k,i)-rpts(k,i-1))**2+(zpts(k,i)-zpts(k,i-1))**2)
        end do
        do i=1,npts
!  calculate metric required for calculating q' and q''
!  (see Miller et al Phys Plasmas 5 (1998) 973
          if (i.eq.1) then
            ip=2
            im=npts-1
          else if (i.eq.npts) then
            ip=2
            im=npts-1
          else
            ip=i+1
            im=i-1
          end if
          drdl=(rpts(k,ip)-rpts(k,im))/(2.*dl)
          d2rdl2=(rpts(k,ip)-2.*rpts(k,i)+rpts(k,im))/dl**2
          dzdl=(zpts(k,ip)-zpts(k,im))/(2.*dl)
          d2zdl2=(zpts(k,ip)-2.*zpts(k,i)+zpts(k,im))/dl**2
          dpsi1dl(i)=(rpts(k,ip)*bppts(k,ip)-rpts(k,im)*bppts(k,im))/(2.*dl)
          d2psi1dl2(i)=(rpts(k,ip)*bppts(k,ip)-2.*rpts(k,i)*bppts(k,i) &
                        +rpts(k,im)*bppts(k,im))/dl**2
          rcinv(i)=-dzdl*d2rdl2+drdl*d2zdl2
          cosu(i)=drdl
          sinu(i)=-dzdl
        end do
        do  i=1,npts
          if (i.eq.1) then
            ip=2
            im=npts-1
          else if (i.eq.npts) then
            ip=2
            im=npts-1
          else
            ip=i+1
            im=i-1
          end if
          drdl=(rpts(k,ip)-rpts(k,im))/(2.*dl)
          d2rdl2=(rpts(k,ip)-2.*rpts(k,i)+rpts(k,im))/dl**2
          dzdl=(zpts(k,ip)-zpts(k,im))/(2.*dl)
          d2zdl2=(zpts(k,ip)-2.*zpts(k,i)+zpts(k,im))/dl**2
          rr=rpts(k,i)
          zz=zpts(k,i)
          if (k.eq.1) surfa(i)=2.*pi*rr
          bphi=fsi/rr
          bth=bppts(k,i)
          bsq=bphi*bphi+bth*bth
!  parameters required got q' and q''...
          psi1=rr*bth
          psi2=0.5*((rr*rcinv(i)+sinu(i))*bth-mu0*rr*rr*pd-ffp)
          psi3=(1./6.)*((-2.*bth*sinu(i)+4.*psi2)*(rcinv(i)+sinu(i)/rr)- &
               d2psi1dl2(i)+cosu(i)*dpsi1dl(i)/rr+ &
               mu0*rr*rr*pd*(rcinv(i)-sinu(i)/rr)+ffp*(rcinv(i)+sinu(i)/rr) &
               -rr*bth*(mu0*rr*rr*pdd+ffpp))
          nu=fsi/(rr*rr*bth)
          nu1=-rcinv(i)-sinu(i)/rr-2.*psi2/psi1
          nu2=(sinu(i)/rr+2.*psi2/psi1)*(rcinv(i)+sinu(i)/rr)   &
               -3.*psi3/psi1+4.*psi2**2/psi1**2
          dnudpsi(i)=nu*ffp/fsi**2+nu*nu1/psi1
          t1(i)=nu/fsi+(nu*fsi/psi1**2)
          t2(i)=(nu*rr*rr/psi1**2)
          t3(i)=(nu/psi1)*(-rcinv(i)-sinu(i)/rr-(2./psi1)*(0.5*(rr*rcinv(i)+sinu(i))*bth))
          nup(k,i)=dnudpsi(i)*circumf(k)/(2.*pi)
          vpparr(i)=(nu1+2.*sinu(i)/rr)/(rr*bth**2)
!!$         vpp=vpp+dl*( &
!!$             1./(rpts(i)**2*bppts(i)**3)* &
!!$             (2.*uprimepts(i)*rpts(i)*bppts(i)+ &
!!$             rpts(i)**2*pprime_loc+ffprime_loc)+ &
!!$             1./(rpts(i-1)**2*bppts(i-1)**3)* &
!!$             (2.*uprimepts(i-1)*rpts(i-1)*bppts(i-1)+ &
!!$             rpts(i-1)**2*pprime_loc+ffprime_loc))
!         drhdpsi=1./psi_1
! dB^2/dpsi required for ballooning calculation....
          dbpdrho=rcinv(i)*bth-mu0*rr*pd-ffp/rr
          dbtdrho=(ffp*bth/fsi-fsi*sinu(i)/rr**2)
          dbsqdpsi(k,i)=2.*(bth*dbpdrho+fsi*dbtdrho/rr)/(rr*bth)
!
          ploty1(k,i)=nu*circumf(k)/(2.*pi)
          ploty2(k,i)=dnudpsi(i)*circumf(k)/(2.*pi)
          plotx(i)=(i-1.)*2.*pi/(npts-1.)
          if (k.eq.10) then
            nupplt(i)=dnudpsi(i)
            y1(i)=drdl
            y2(i)=dzdl
            y3(i)=d2rdl2
            y4(i)=d2zdl2
            y5(i)=rcinv(i)
            y6(i)=psi1
            chiplt(i)=(i-1.)*2.*pi/(npts-1.)
          end if
          d2nudpsi2(i)=(ffpp/fsi-ffp**2/fsi**3)*nu/fsi &
              +2.*ffp*nu*nu1/(psi1*fsi**2)+2.*nu*(nu2-nu1*psi2/psi1)/psi1**2
!  used to calculate  field integral round boundary
          bpsur(i)=sqrt(bsq)/bth
          bsqar(i)=bsq/bth
          bphiar(i)=bphi/bth
          rsqar(i)=rr*rr/bth
          rrar(i)=rr/bth
          riar(i)=1./(rr*bth)
          risar(i)=1./(rr*rr*bth)
          rar(i)=1./bth
          avblir(i)=1./((rr**4)*(bth**3))
          sfloc(i)=bphi/(rr*bth)
          sfloc(i)=nu
          bmod=sqrt(bsq)
          if (bmod.gt.bmax(k)) bmax(k)=bmod
          if (i.eq.1) then
            bmin(k)=bmod
          else
            if (bmod.lt.bmin(k)) bmin(k)=bmod
          end if
        end do
        do i=1,npts
          rr=rpts(k,i)
          zz=zpts(k,i)
          bphi=fsi/rr
          bth=bppts(k,i)
          bsq=bphi*bphi+bth*bth
          root=1.-sqrt(bsq)/bmax(k)
          erribm=-1.0e-10
          if (root.lt.erribm) then
           write(nw,*)'fatal error***root<0 in flxav, root=',root
           write(nw,*)' k=',k,' bsq=',bsq,' bmax=',bmax(k)
           write(nw,*)' bth=',bth,' bph=',i,bphi,' fsi=',fsi,' rr=',rr
           write(nw,*)' psi=',psi
           stop
          end if
          root=abs(root)
          colop(i)=sqrt(root)/bth
        end do
!  surface area
        if (k.eq.1) then
          call flxint(surfa,k,surar)
        end if
!  safety factor...
        call flxint(sfloc,k,ant)
        sfac(k)=ant/(2.*pi)
        call flxint(dnudpsi,k,ant)
        qp(k)=ant/(2.*pi)
        call flxint(t1,k,ant1)
        call flxint(t2,k,ant2)
        call flxint(t3,k,ant3)
        if (k.eq.100) write(6,*)' psi=',psiv(k)/umax,' T1=',ant1, &
        ' t2=',ant2,' T3=',ant3,' qp=',(ant1*ffp/fsi+ant2*mu0*pd+ant3)/(2.*pi)
        if (k.eq.100)   &
            write(6,*)' mu0*p-prime=',mu0*press(psiv(k),1),' f-prime=',fprof(psiv(k),1)/fprof(psiv(k),2),' f=',fprof(psiv(k),2)
        call flxint(vpparr,k,ant)
        vpp(k)=ant 
       !if (k.eq.5) write(6,*)' qp=',qp(k)
        call flxint(d2nudpsi2,k,ant)
        qpp(k)=ant/(2.*pi)
        call flxint(rar,k,rnor)
        call flxint(avblir,k,ant)
        avbli(k)=ant/rnor
        call flxint(bsqar,k,ant)
        bsqav(k)=ant/rnor
        call flxint(bphiar,k,ant)
        bphiav(k)=ant/rnor
        call flxint(rsqar,k,ant)
        rsqav(k)=ant/rnor
        call flxint(riar,k,ant)
        rinv(k)=ant/rnor
        call flxint(risar,k,ant)
        rsqinv(k)=ant/rnor
        call flxint(rrar,k,ant)
        rav(k)=ant/rnor
        rnorm(k)=rnor
        call flxint(colop,k,rint)
!  collisionality
        te=bk*tempe(psi,0)
        ti=bk*tempi(psi,1,0)
        tau=te/ti
        ne=dense(psi,0)
!  coulomb logarithm
        coolog=log(sqrt(ne*1.0d-6)*bk/te)
        coolog=24.-coolog
        bfac=bmax(k)*((bmax(k)+bmin(k))/(bmax(k)-bmin(k)))**2
        bfac=bfac*rint
!  electron--ion collision time
        epsfac=(eps0/eq)**2
        ratio=te/eq
        pro=0.
        do i=1,nimp+1
          zni=densi(psi,i,0)
          pro=pro+iz(i)*iz(i)*zni*eq
        end do
        rt=sqrt(me)*sqrt(te)
        botcol=pro*coolog
        topcol=rt*epsfac*ratio
        colte=(3.*2.*pi*sqrt(2.*pi))*topcol/botcol
        vthe=sqrt(2.*te/me)
        tnue(k)=bfac/(4.*colte*vthe)
        cnue(k)=sqrt(2.)*rnorm(k)*bmax(k)/(2.*pi*vthe*colte)
        if (te.eq.0.) then
!  Put collisionalities to arbitrary large number
          tnue(k)=1.0d10
          cnue(k)=1.0d10
        end if
!  ion--ion collision time
        epsfac=(eps0/eq)**2
        ratio=ti/eq
!  main species ion density
        zni=densi(psi,1,0)
        pro=zm*zm*zm*zm*zni*eq
        rt=sqrt(ti)
        botcol=pro*coolog
        topcol=rt*epsfac*ratio
        colti=(3.*4.*pi*sqrt(pi))*topcol/botcol
        vthi=sqrt(2.*ti)
        tnui(k)=bfac/(4.*colti*vthi)
        cnui(k)=sqrt(2.)*rnorm(k)*bmax(k)/(2.*pi*vthi*colti)
        if (ti.eq.0.) then
!  Put collisionalities to arbitrary large number
          tnui(k)=1.0d10
          cnui(k)=1.0d10
        end if
! integral of B/Btheta dl around the flux surface...
        call flxint(bpsur,k,ant)
        bdl(k)=ant
        bav(k)=ant/rnor
      end do
!  some axis averages are undefined, and we set to zero...
!      do k=2,ncon-1
!        write(6,*)' vpp=',vpp(k),' numerical=',(rnorm(k+1)-rnorm(k-1))/(psiv(k+1)-psiv(k-1))
!      end do
      fsi=fprof(0.0d0,2)
      avbli(ncon)=0.
      rnorm(ncon)=2.*pi*r0*r0*sfac(ncon)/fsi
      !  some axis values we know analytically....

      !Added by Bhavin
      bdl(ncon) = fsi/r0
      bav(ncon) = fsi/r0
      bsqav(ncon)=(fsi/r0)**2
      bphiav(ncon) = fsi/r0
      rsqav(ncon)=r0**2
      rinv(ncon)=1./r0
      rsqinv(ncon)=1./r0**2
      rav(ncon)=r0
!  calculate axis values by fitting on mid-plane...
      x1=rpts(ncon-1,1)
      x2=rpts(ncon-1,npts/2+1)
      x3=rpts(ncon-2,npts/2+1)
      q1=sfac(ncon-1)
      q2=sfac(ncon-1)
      q3=sfac(ncon-2)
      av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
      bv=(q1-q2)/(x1-x2)-av*(x1+x2)
      cv=q1-av*x1**2-bv*x1
      sfac(ncon)=av*r0**2+bv*r0+cv
      qp(ncon)=0.
      qpp(ncon)=0.
      q1=tnue(ncon-1)
      q2=tnue(ncon-1)
      q3=tnue(ncon-2)
      av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
      bv=(q1-q2)/(x1-x2)-av*(x1+x2)
      cv=q1-av*x1**2-bv*x1
      tnue(ncon)=av*r0**2+bv*r0+cv
      q1=cnue(ncon-1)
      q2=cnue(ncon-1)
      q3=cnue(ncon-2)
      av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
      bv=(q1-q2)/(x1-x2)-av*(x1+x2)
      cv=q1-av*x1**2-bv*x1
      cnue(ncon)=av*r0**2+bv*r0+cv
      q1=tnui(ncon-1)
      q2=tnui(ncon-1)
      q3=tnui(ncon-2)
      av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
      bv=(q1-q2)/(x1-x2)-av*(x1+x2)
      cv=q1-av*x1**2-bv*x1
      tnui(ncon)=av*r0**2+bv*r0+cv
      q1=cnui(ncon-1)
      q2=cnui(ncon-1)
      q3=cnui(ncon-2)
      av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
      bv=(q1-q2)/(x1-x2)-av*(x1+x2)
      cv=q1-av*x1**2-bv*x1
      cnui(ncon)=av*r0**2+bv*r0+cv
!   calculate fraction of trapped particles...
      do k=1,ncon
        psi=psiv(k)
        eps=epsv(k)
        fsi=fprof(psi,2)
        if (neo.ge.0) then
          if (k.eq.ncon) then
            ftrap(k)=0.
            ftrapd(k)=0.
          else
            binv=1./bmax(k)
            nst=100
            dst=binv/nst
            bigint=0.
            bigint2=0.
            do i=1,nst
              rla=(i-1)*dst+0.5*dst
              do l=1,npts
                rr=rpts(k,l)
                zz=zpts(k,l)
                bphi=fsi/rr
                bp=bppts(k,l)
                btot=sqrt(bp*bp+bphi*bphi)
                bla(l)=(sqrt(1.-rla*btot))/bp
                bla2(l)=btot/(bp*sqrt(1.-rla*btot))
              end do
              call flxint(bla,k,rlag)
              rlag=rlag/rnorm(k)
              bigint=bigint+rla*dst/rlag
              call flxint(bla2,k,rlag)
              rlag=rlag/rnorm(k)
              bigint2=bigint2+dst/rlag
            end do
!  fraction of trapped ptcles. on kth flux surface (k=1 => outside)
            ftrap(k)=1.-3.*bsqav(k)*bigint/4.
            ftrapd(k)=1.-3.*bsqav(k)*bigint2/2.
          end if
        end if
!  calculate fast ptcle bootstrap current
!        if (ifast.eq.1) call fastbs
!  calculate bootstrap current...
        if (k.eq.ncon) then
!!$          fc=1.
!!$          x=0.
!!$          pt=press(psi,0)
!!$          pd=press(psi,1)
!!$          ne=dense(psi,0)
!!$          te=tempe(psi,0)
!!$          pe=bk*ne*te
!!$          rj0=fsi*pe
!!$          bsq=bsqav(k)
!!$!  calculate correlation length
!!$          if (imat.eq.1) then
!!$            call conlen(k,fsi,conl,dotav)
!!$            if ((te.eq.0.).or.(ne.eq.0.).or.(ti.eq.0.)) then
!!$!  if the density or temperature=0 assume bootstrap current =0
!!$              bstrap=0.
!!$            else
!!$              call hirsig(psi,pt,pd,fc,fsi,conl,dotav,bsq,bstrap,eps,k)
!!$            end if
!!$          end if
          bstrap=0.
        else
           fc=1.-ftrap(k)
          x=ftrap(k)/fc
          pt=press(psi,0)
          pd=press(psi,1)
          ne=dense(psi,0)
          te=tempe(psi,0)
          pe=bk*ne*te
          rj0=fsi*pe
          bsq=bsqav(k)
!  calculate correlation length
          if (imat.eq.1) then
            call conlen(k,fsi,conl,dotav)
            if ((te.eq.0.).or.(ne.eq.0.).or.(ti.eq.0.)) then
!  if the density or temperature=0 assume bootstrap current =0
              bstrap=0.
            else
              call hirsig(psi,fc,fsi,conl,dotav,bsq,bstrap,eps,k)
            end if
          else
            if ((te.eq.0.).or.(ne.eq.0.).or.(ti.eq.0.)) then
!  if the density or temperature=0 assume bootstrap current =0
              bstrap=0.
            else
              call hirsh(psi,rj0,x,bsq,k,bstrap)
            end if
          end if
        end if
        bsj(k)=bstrap
      end do
!  neoclassical resistivity enhancement factor
      do k=1,ncon
        psi=psiv(k)
        eps=epsv(k)
        pt=press(psi,0)
        te=tempe(psi,0)
        ti=tempi(psi,1,0)
        tau=te/ti
        te=bk*te
        ne=dense(psi,0)
!  zeff
        zeff=zm
        if (imp.eq.1) then
          if (ne.gt.0.) then
            zeff=0.
            do l=1,nimp+1
              zni=densi(psi,l,0)
              zeff=zeff+(zni*iz(l)**2)/ne
            end do
!	  else
!	    write(nw,*)'error*** problem in flxav, ne=0'
!	    write(nw,*)'cannot evaluate zeff'
!	    stop
          end if
        end if
        coolog=log(sqrt(ne*1.0e-6)*bk/te)
        coolog=24.-coolog
!  electron--electron collision time
!	epsfac=(4.*pi*eps0/(eq))**2
        epsfac=(4.*pi*eps0/(zeff*eq))**2
        ratio=te/eq
        pro=ne*eq
        rt=sqrt(te)
        botcol=sqrt(2.*pi)*pro*coolog
        topcol=rt*epsfac*ratio
        colte=(3./4.)*topcol/botcol
        vthe=sqrt(2.*te)
!  safety factor
        qq=sfac(k)
!  electron collisionality (hirshman, hawryluk, birge)
        if (k.lt.ncon) then
          rnust=sqrt(2.)*r0*qq/(eps*colte*vthe*sqrt(eps))
          ft=ftrap(k)
        else
! set rnust is a dummy variable at the centre (ft=0)
          rnust=100.
          ft=0.
        end if
! If Te=0, conductivity -> so neoclassical enhancement is irrelevant
! -> rnust is a dummy variable which we set to be large
        if (te.eq.0.) rnust=1.0d9
        if (nco.eq.0) rnust=0.
!  neoclassical conductivity (hirshman, hawryluk, birge)
        call signeo(rnust,ft,zeff,sigfac)
        sighhb(k)=sigfac
      end do
      do k=2,ncon-1
        x1=psiv(k-1)
        x2=psiv(k)
        x3=psiv(k+1)
        q1=sfac(k-1)
        q2=sfac(k)
        q3=sfac(k+1)
        av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
        bv=(q1-q2)/(x1-x2)-av*(x1+x2)
        fq(k)=2.*av*x2+bv
        fqd(k)=2.*av
!!$        p1=press(x1,0)
!!$        p2=press(x2,0)
!!$        p3=press(x3,0)
!!$        av=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
!!$        bv=(p1-p2)/(x1-x2)-av*(x1+x2)
!!$        ppar(k)=2.*av*x2+bv
!!$        pppar(k)=2.*av
!!$        ppan(k)=1.03*press(x2,1)
!!$        pppan(k)=1.03*press(x2,2)
!        if (k.ge.90)   then
        !write(6,*)' k=',k,' psin=',psiv(k)/umax,' numerical shear=',fq(k),' analytic=',qp(k),' q=',sfac(k)
        !write(6,*)' k=',k,' shat=',psiv(k)*fq(k)/sfac(k),' s_hat=',psiv(k)*qp(k)/sfac(k),             &
        !' press=',press(psiv(k),0),' pp=',press(psiv(k),1)
!        end if
!        write(6,*)' q=',sfac(k),' psin=',psiv(k)/umax,' qp-norm=',qp(k)*umax
        do i=1,npts
!          x1=psiv(k-1)
!          x2=psiv(k)
!          x3=psiv(k+1)
          x1=sfac(k-1)
          x2=sfac(k)
          x3=sfac(k+1)
          q1=ploty1(k-1,i)
          q2=ploty1(k,i)
          q3=ploty1(k+1,i)
          av=((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/(x1-x3)
          bv=(q1-q2)/(x1-x2)-av*(x1+x2)
          ploty3(i)=(2.*av*x2+bv)*qp(k)
          ploty4(i)=ploty2(k,i)
        end do
!        call tstplt2(npts,plotx,ploty3,ploty4)
      end do
!!$      ppan(1)=ppan(2)
!!$      pppan(1)=pppan(2)
!!$      ppar(1)=ppar(2)
!!$      pppar(1)=pppar(2)
!!$      ppan(ncon)=ppan(ncon-1)
!!$      pppan(ncon)=pppan(ncon-1)
!!$      ppar(ncon)=ppar(ncon-1)
!!$      pppar(ncon)=pppar(ncon-1)
      fq(1)=0.
      fq(1)=(0.75*sfac(1)-sfac(2)+0.25*sfac(3))/  &
              (0.75*psiv(1)-psiv(2)+0.25*psiv(3))
      !write(6,*)' fqend=',fq(1)
      do k=2,ncon-2
        fqd(k)=(fq(k+1)-fq(k-1))/(psiv(k+1)-psiv(k-1))
      end do
      fqd(1)=(0.75*fq(1)-fq(2)+0.25*fq(3))/  &
              (0.75*psiv(1)-psiv(2)+0.25*psiv(3))
      fqd(ncon)=0.
      fqd(ncon-1)=0.
!!$       call tstplt(npts,chiplt,nupplt)
!!$       call tstplt(npts,chiplt,y1)
!!$       call tstplt(npts,chiplt,y2)
!!$       call tstplt(npts,chiplt,y3)
!!$       call tstplt(npts,chiplt,y4)
!!$       call tstplt(npts,chiplt,y5)
!!$       call tstplt(npts,chiplt,y6)
!       call tstplt2(ncon,psiv,ppar,ppan)
!       call tstplt2(ncon,psiv,pppar,pppan)
!!$       call tstplt(ncon,psiv,qpp)
!       call tstplt2(ncon,psiv,qpp,fqd)
!       call tstplt(ncon,psiv,sighhb)
!!$       call tstplt(ncon,psiv,rnorm)
!!$       call tstplt(ncon,psiv,rsqav)
!!$       call tstplt(ncon,psiv,rinv)
!!$       call tstplt(ncon,psiv,rsqinv)
!!$       call tstplt(ncon,psiv,rav)
!!$       call tstplt(ncon,psiv,ftrap)
!!$       call tstplt(ncon,psiv,bsj)
!!$       call tstplt(ncon,psiv,cnue)
!!$       call tstplt(ncon,psiv,cnui)
!!$       call tstplt(ncon,psiv,tnue)
!!$       call tstplt(ncon,psiv,tnui)

   end subroutine flxav
!
!*****************************************************************
!
      subroutine flxint(arr,k,ans)
!     ***************************
!
!  Calculates the flux surface integral of the variable held in
!  arr(i) calculated at flxr(k,i), flxz(k,i). k labels the flux
!  surface around which we are integrating
!
      use param
      implicit none
      double precision, intent(in) :: arr(npts)
      integer, intent(in) :: k
      double precision, intent(out) :: ans
      integer :: i
      double precision ::  dl
!
      ans=0.
      dl=circumf(k)/(npts-1.)
      do i=1,npts-1
        ans=ans+arr(i)*dl
      end do
   end subroutine flxint
!
!***************************************************************
 end module flux_average
