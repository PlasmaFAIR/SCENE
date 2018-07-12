      subroutine balloon(nf)
!     *********************
!
!  Performs a ballooning stability calculation for the nf flux surface
!  (nf=1 is the edge, nf=ncon is the magnetic axis
!  Note rpts and zpts are stored on equal arc length mesh, with npts
!  nup is calculated already in flxav

!
      use param
      use balpar
      implicit none
      integer nf,i,ip,im,k,j,ichi,jpass,kpass,ii
      double precision fprof,press,balfun,balfun2
      double precision fbal(2*npts*nturns,2),xipar(2*npts*nturns)
      double precision eq3(2*npts*nturns),eq4(2*npts*nturns)
      double precision fdbal(2*npts*nturns,2)
      double precision yp(npts),xp(npts)
      double precision eqq(npts),jac(npts),dt1dchi(npts),dt2dchi(npts)
      double precision dbsqdchi(npts)
      double precision pd,psi,bsq,fsi
      double precision chiv,bth,rr,dchi,dl,chi0,dchi0
      double precision chimax,chimin,bsqp,bsqm,t1p,t1m,t2p,t2m
      double precision hh,an,bn,cn,dn,betn,deln,cp,cm
      double precision lamold,lamnew,difold,difnew,err
      double precision x,y,yd,sumit,sumit1,s,alpha,hfun,hdfun
      double precision yplot(2*npts*nturns),y2(2*npts*nturns),nuarr(2*npts*nturns)
!
      jpass=0
      kpass=100
      psi=psiv(nf)
      pd=mu0*press(psi,1)
      fsi=fprof(psi,2)
!  Calculate and store equilibrium variables required....
      yp(1)=0.
      xp(1)=0.
      bsqp=(fprof(psiv(nf+1),2)/rpts(nf+1,1))**2+bppts(nf+1,1)**2
      bsqm=(fprof(psiv(nf-1),2)/rpts(nf-1,1))**2+bppts(nf-1,1)**2
      write(6,*)' calc=',dbsqdpsi(nf,1),' num=',       &
                (bsqp-bsqm)/(psiv(nf+1)-psiv(nf-1))

      dl=circumf(nf)/(npts-1)
      dchi=2.*pi/(npts-1)
      do i=1,npts
         bth=bppts(nf,i)
         jac(i)=dl/(bth*dchi)
      end do
      do i=1,npts
         ip=i+1
         if (ip.gt.npts) ip=2
         im=i-1
         if (im.lt.1) im=npts-1
         rr=rpts(nf,i)
         bth=bppts(nf,i)
         bsq=(fsi/rr)**2+bth**2
         if (i.gt.1) then
           xp(i)=(i-1)*dchi
           chiv=(i-1)*dchi
           yp(i)=yp(i-1)+0.5*(nup(nf,i-1)+nup(nf,i))*dchi-qp(nf)*dchi
!        write(6,*)' i=',i,' chi=',chiv,' yp=',yp(i),' nup=',nup(nf,i),' qp=',qp(nf)
!        write(6,*)' R=',rpts(nf,i),' Z=',zpts(nf,i),' bth=',bppts(nf,i)
         end if
         nuarr(i)=nup(nf,i)
         eqq(i)=(2.*pd+dbsqdpsi(nf,i))/bsq
         yplot(i)=-rcen*rcen*bppts(nf,1)*eqq(i)/2.
         bsqp=(fsi/rpts(nf,ip))**2+bppts(nf,ip)**2
         bsqm=(fsi/rpts(nf,im))**2+bppts(nf,im)**2
         dbsqdchi(i)=(bsqp-bsqm)/(2.*dchi)
         t1p=1./(jac(ip)*(rpts(nf,ip)*bppts(nf,ip))**2)
         t1m=1./(jac(im)*(rpts(nf,im)*bppts(nf,im))**2)
         dt1dchi(i)=(t1p-t1m)/(2.*dchi)
         t2p=(rpts(nf,ip)*bppts(nf,ip))**2/(jac(ip)*bsqp)
         t2m=(rpts(nf,im)*bppts(nf,im))**2/(jac(im)*bsqm)
         dt2dchi(i)=(t2p-t2m)/(2.*dchi)
      end do
!!$      call tstplt(npts,xp,yp,0.,0.)
!!$      call tstplt(npts,xp,yplot,0.,0.)
!      write(6,*)' yp1=',yp(1),' ypn=',yp(npts)
!      write(6,*)' bp1=',bppts(nf,1),' bpn=',bppts(nf,npts)
!      write(6,*)' R1=',rpts(nf,1),' Rn=',rpts(nf,npts)
      nchi=nturns*(npts-1)+1
      chimax=pi*nturns
      chimin=-pi*nturns
      dchi=(chimax-chimin)/(nchi-1)
      if (nchi0.gt.1) dchi0=2.*pi/(nchi0-1)
      do ichi=1,nchi0
        if (nchi0.eq.1) then
!  shift chi0 by pi as chi=0 is on inboard side...silly, but there you go!
          chi0=chi0val+pi
        else
          chi0=-pi+(ichi-1)*dchi0
        end if
        k=0
!      write(6,*)' *************psi=',psiv(nf),' ***************'
        do j=1,nturns+1
          do i=1,npts-1
            if (k.lt.nchi) then
            rr=rpts(nf,i)
            bth=bppts(nf,i)
            bsq=(fsi/rr)**2+bth**2
            k=k+1
            chi(k)=chimin+(k-1)*dchi
            eq1(k)=(1./(rr*bth)**2+((rr*bth)**2/bsq)*        &
                   (qp(nf)*(chi(k)-chi0)+yp(i))**2)/jac(i)
            eq2(k)=pd*(jac(i)*eqq(i)-(qp(nf)*(chi(k)-chi0)+yp(i))    &
                           *fsi*dbsqdchi(i)/bsq**2)
            eqd1(k)=2.*(qp(nf)*(chi(k)-chi0)+yp(i))*                    &
                    nup(nf,i)*(rr*bth)**2/(jac(i)*bsq)+dt1dchi(i)+     &
                         ((qp(nf)*(chi(k)-chi0)+yp(i))**2)*dt2dchi(i)
            end if
          end do
        end do
!  Do shooting to find eigenvalue, lambda....
!  j=1 is even solution, j=2 is odd
        fbal(nchi/2+1,1)=1.
        fdbal(nchi/2+1,1)=0.
        fbal(nchi/2+1,2)=0.
        fdbal(nchi/2+1,2)=1.
        hh=dchi/2.
        lam=lamges
        do
          jpass=jpass+1
          do j=1,2
            !  First shoot out from zero to max chi...
            do i=nchi/2+1,nchi-1
              x=chi(i)
              y=fbal(i,j)
              yd=fdbal(i,j)
              an=hh*balfun(x,i,y,yd)
              x=chi(i)+hh
              betn=hh*(yd+an/2.)
              y=fbal(i,j)+betn
              yd=fdbal(i,j)+an
              bn=hh*balfun(x,i,y,yd)
              yd=fdbal(i,j)+bn
              cn=hh*balfun(x,i,y,yd)
              deln=hh*(fdbal(i,j)+cn)
              y=fbal(i,j)+deln
              yd=fdbal(i,j)+2.*cn
              ii=i+1
              if (ii.ge.nchi) ii=nchi-1
              x=chi(i+1)
              dn=hh*balfun(x,ii,y,yd)
              fbal(i+1,j)=fbal(i,j)+dchi*(fdbal(i,j)+(an+bn+cn)/3.)
              fdbal(i+1,j)=fdbal(i,j)+(an+2.*bn+2.*cn+dn)/3.
!!$            write(61,*)' j=',j,' i+1=',i+1,' fdba=',fdbal(i+1,j)
            end do
          end do
          cp=-fbal(nchi,1)/fbal(nchi,2)
          do i=nchi/2+1,nchi
            fbal(i,1)=fbal(i,1)+cp*fbal(i,2)
            fdbal(i,1)=fdbal(i,1)+cp*fdbal(i,2)
          end do
          j=1
          !  and now back to min chi...
          do i=nchi/2+1,2,-1
            x=chi(i)
            y=fbal(i,j)
            yd=fdbal(i,j)
            an=-hh*balfun(x,i,y,yd)
            x=chi(i)-hh
            betn=-hh*(yd+an/2.)
            y=fbal(i,j)+betn
            yd=fdbal(i,j)+an
            bn=-hh*balfun(x,i,y,yd)
            yd=fdbal(i,j)+bn
            cn=-hh*balfun(x,i,y,yd)
            deln=-hh*(fdbal(i,j)+cn)
            y=fbal(i,j)+deln
            yd=fdbal(i,j)+2.*cn
            x=chi(i-1)
            ii=i-1
            if (ii.le.1) ii=2
            dn=-hh*balfun(x,i,y,yd)
            fbal(i-1,j)=fbal(i,j)-dchi*(fdbal(i,j)+(an+bn+cn)/3.)
            fdbal(i-1,j)=fdbal(i,j)+(an+2.*bn+2.*cn+dn)/3.
          end do
          difnew=fbal(1,1)
          !  Adjust lambda until cp=cm
          if (jpass.eq.1) then
            lamnew=0.999*lam
            lamold=lam
            difold=difnew
            lam=lamnew
            cycle
          else if (jpass.le.kpass) then
            lamnew=(difnew*lamold-difold*lam)/(difnew-difold)
            write(6,*)' lam=',lam,' difnew=',difnew
            err=abs((lam-lamnew)/(lam+lamnew))
            if (err.gt.1.0e-5) then
              lamold=lam
              lam=lamnew
              difold=difnew
              cycle
            end if
          else
            !  lambda=0 flags that could not find eigenvalue
            lam=0.
          end if
          exit
        end do

        write(92,35)chi0val,psiv(nf)/umax,lam
        write(6,*)' Flux surface=',nf,' psi=',psiv(nf)/umax,' lmbda=',lam
        lamges=lam
        do i=1,nchi
          fb(i)=fbal(i,1)
          write(61,*)' i=',i,' fdbal=',fdbal(i,1),fdbal(i,2)
          yplot(i)=fb(i)
          y2(i)=exp(-chi(i)**2/2.)
        end do
        sumit=0.
        sumit1=0.
        do k=1,nchi-1
          sumit=sumit+0.5*(eq2(k)*yplot(k)**2+eq2(k+1)*yplot(k+1)**2)
          sumit1=sumit1+0.5*(eq1(k)*fdbal(k,1)**2+eq1(k+1)*fdbal(k+1,1)**2)
        end do
! Now shoot to calculate xi_parallel, assuming it is odd in chi
        xipar(nchi/2+1)=0.
        k=nchi/2+1
        do i=nchi/2+1,nchi-1
          x=chi(i)
          if (x.lt.pi) then
            k=i
            xipar(k)=0.
            cycle
          end if
          an=dchi*balfun2(x,i,pd)
          x=chi(i)+dchi/2.
          bn=dchi*balfun2(x,i,pd)
          cn=bn
          x=chi(i+1)
          dn=dchi*balfun2(x,i,pd)
          xipar(i+1)=xipar(i)+(an+2.*bn+2.*cn+dn)/6.
        end do
        im=k
        do i=k+1,nchi
          im=im-1
          xipar(im)=-xipar(i)
        end do
        do i=im,nchi
          y2(i)=xipar(i)
        end do
        do i=1,im-1
          y2(i)=0.
        end do
!        call tstplt2(nchi,chi,yplot,y2)
!!$        call tstplt2(nchi,chi,eq2,eq4)
!        call tstplt(nchi,chi,eq3,0.,0.)
!        call tstplt(nchi,chi,eq4,0.,0.)
!        call tstplt(nchi,chi,eqd1,0.,0.)
!        call tstplt(nchi,chi,eq2,0.,0.)
!  end loop over chi0
      end do
 35   format(' chi0 (rad)=',f12.5,' flux=',e14.6,' lambda=',e14.6)
      return
  end subroutine balloon
!
!---------------------------------------------------------------------
!
      function balfun(x,k,y,yd)
!     -----------------------
!
!  Used for the solution to the ballooning equation
!
      use balpar
      implicit none
      integer k,pm
      double precision balfun,x,y,yd
      double precision rat,c1,c2,pi
      double precision alpha,s,hfun,hdfun
!
      pi=4.*atan(1.)
      pm=1
      if (x.lt.0.) pm=-1
      rat=(x-chi(k))/(chi(k+pm)-chi(k))
      if ((rat.lt.0.).or.(rat.gt.1.)) then
        write(6,*)' ERROR in balfun, rat=',rat
        write(6,*)' x=',x,' chik=',chi(k),' chik+1=',chi(k+1)
        stop
      end if
      c1=eqd1(k)/eq1(k)+rat*(eqd1(k+pm)/eq1(k+pm)-eqd1(k)/eq1(k))
      c2=eq2(k)/eq1(k)+rat*(eq2(k+pm)/eq1(k+pm)-eq2(k)/eq1(k))
      balfun=-c1*yd-lam*c2*y
!      balfun=(lam+x**2)*y
      return
  end function balfun
!---------------------------------------------------------------------
!
      function balfun2(x,k,pd)
!     -----------------------
!
!  Used to calculate xi_parallel from solution to ballooning equation
!  x= theta
!  k is the mesh point (x=x(k))
!  y(k) is the solution to the ballooning equation at x=x(k)
!
      use balpar
      implicit none
      integer k,pm
      double precision balfun2,x,pd,pi
      double precision rat,c1,c2
!
      pi=4.*atan(1.)
      pm=1
      if (x.lt.0.) pm=-1
      rat=(x-chi(k))/(chi(k+pm)-chi(k))
      if ((rat.lt.0.).or.(rat.gt.1.)) then
        write(6,*)' ERROR in balfun, rat=',rat
        write(6,*)' x=',x,' chik=',chi(k),' chik+1=',chi(k+1)
        stop
      end if
      c1=eq2(k)+rat*(eq2(k+pm)-eq2(k))
      c2=fb(k)+rat*(fb(k+pm)-fb(k))
      balfun2=c1*c2/pd
      return
  end function balfun2
