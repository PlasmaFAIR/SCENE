!*****************************************************************
!
      subroutine bdry(x,y,f,igo)
!     **********************
!
! Sets up the boundary of the flux surface
!
      use param
      implicit none
      double precision errg,dk,dth,kscl,th,rnor,znor,sg,s2g,ratr,rat
      double precision rv1,rv2,f1,f2,rup,rlo
      double precision x,y,f
      double precision xv(npts),yv(npts),gth(npts),ang(npts),rtemp(npts)
      double precision work(npts),ang2(npts)
      integer nit,i,j,ipit,jx,igo,iout,iv
!
      if (igo.eq.0) then
!  set up r(theta,psi)
        errg=1.0d-6
        nit=20
        jx=npts/4+1
        allocate(theta(npts),rth(npts,ncon))
        theta=0.; xv=0.; yv=0.; gth=0.; rth=0.; ang2=0.
!!$!  put X-points in using modified Bishop formula
        dk=kval/(ncon-1)
        if (kval.eq.0.) dk=1./(ncon-1.)
        dth=2.*pi/(npts-1.)
        ang(1)=-pi
        theta(1)=-pi
        do j=2,npts
          ang(j)=ang(1)+(j-1)*dth
!          write(91,*)' j=',j,' ang=',ang(j),' dth=',dth
          theta(j)=theta(1)+(j-1)*dth
        end do
        do i=ncon,2,-1
          if (kval.eq.0) then
! simple circle...
             kscl=dk*(i-1.)
             do j=1,npts
               gth(j)=1.
             end do
          else
            kscl=dk*(i-1.)
            gth(1)=sqrt(sqrt(4.+(kscl)**2)-2.)/(sqrt(1.+kscl)-1.)
            rv1=gth(1)*(sqrt(1.+kscl)-1.)
            th=pi
            f1=rv1**2*(4.+rv1**2-4.*rv1*sin(th)**2)-kscl**2
            gth(npts)=gth(1)
            do j=2,npts-1
              th=ang(j)
              rv1=gth(j-1)*(sqrt(1.+kscl)-1.)
              f1=rv1**2*(4.+rv1**2-4.*rv1*sin(th)**2)-kscl**2
              ipit=0
 10           ipit=ipit+1
              if (ipit.eq.1) then
                rv2=0.99*rv1
                f2=rv2**2*(4.+rv2**2-4.*rv2*sin(th)**2)-kscl**2
                rup=(f1*rv2-f2*rv1)/(f1-f2)
                rv1=rv2
                rv2=rup
                f1=f2
                goto 10
              else if (ipit.lt.nit) then
                f2=rv2**2*(4.+rv2**2-4.*rv2*sin(th)**2)-kscl**2
                rup=(f1*rv2-f2*rv1)/(f1-f2)
                rv1=rv2
                rv2=rup
                f1=f2
                if (abs((rv1-rv2)/(rv1+rv2)).gt.errg) goto 10
              else
                write(6,*)' Error*** could not converge r(theta)'
                stop
              end if
              gth(j)=rv2/(sqrt(1.+kscl)-1.)
            end do
          end if
          if (i.eq.ncon) then
            rnor=gth(1)
            znor=gth(jx)/rnor
          end if
          do j=1,npts
            gth(j)=gth(j)/rnor
           end do
          do j=1,npts
!            if (i.eq.2) write(91,*)' j h=',j,' ang=',ang(j)
            sg=sin(ang(j))
 !           if (i.eq.2) write(91,*)' j 4=',j,' ang=',ang(j),' sg=',sg
            s2g=sin(2.*ang(j))
            if (kval.gt.0.) then
              xv(j)=kscl*gth(j)*cos(ang(j)+tri*sg+quad*s2g)/kval
              yv(j)=kscl*gth(j)*elon*sg/(znor*kval)
!            if (i.eq.2) write(91,*)' j2=',j,' ang=',ang(j)
            else
              xv(j)=kscl*gth(j)*cos(ang(j)+tri*sg+quad*s2g)
              yv(j)=kscl*gth(j)*elon*sg/(znor)
            end if
!            if (i.eq.2) write(91,*)' j 3=',j,' ang=',ang(j)
!            if (i.eq.2) write(91,*)' j=',j,' gth before=',gth(j),' xv=',xv(j),' yv=',yv(j)
            gth(j)=sqrt(xv(j)**2+yv(j)**2)
!            if (i.eq.2) write(91,*)' j=',j,' gth=',gth(j),' ang=',ang(j),' sg=',sg
            ang2(j)=atan2(yv(j),xv(j))
          end do
          call spline1d(rtemp,theta,npts,gth,ang2,npts,work)
          do j=1,npts
            rth(j,i)=rtemp(j)
!            if (i.eq.2) write(91,*)' j=',j,' rth=',rth(j,i)
          end do
        end do
!  Need to scale rth to get aspect ratio correct
        rlo=rth(1,ncon)
        rup=rth(npts/2+1,ncon)
        do i=1,ncon
          do j=1,npts
            rth(j,i)=2.*tokeps*rcen*rth(j,i)/(rup+rlo)
!            write(91,*)' i=',i,' j=',j,' rth(j,i)=',rth(j,i)
          end do
        end do
        return
      end if
!  igo=1 evaluate f
      rup=sqrt((x-rcen)**2+y**2)
      th=atan2(y,x-rcen)
      dth=2.*pi/(npts-1.)
      j=int((th+pi)/dth)+1
      if (j.ge.npts) j=npts-1
      rat=(th-theta(j))/(theta(j+1)-theta(j))
      if ((rat.gt.1.).or.(rat.lt.0.)) then
        write(6,*)' Error in rat in bdry, rat=',rat
        write(6,*)' j=',j,' th=',th,' thetaj=',theta(j),' theta(j+1)=',theta(j+1)
        stop
      end if
      iv=0
      rv1=rth(j,ncon)+rat*(rth(j+1,ncon)-rth(j,ncon))
      rv2=rup
      iout=0
      if (rv1.lt.rv2) then
!  outside plasma
!         rv2=rv1-(rv2-rv1)
!         iout=1
        kscl=kval*rv2/rv1
        if (kval.eq.0.) kscl=rv2/rv1
      else
        if (rv2.lt.0.) then
          write(6,*)' rv2 leq 0***problem in bdry'
          write(6,*)' rup=',rup,' rv1=',rv1
          write(6,*)' x=',x,' rcen=',rcen,' y=',y
          stop
        end if
        do 20 i=1,ncon
          rv1=rth(j,i)+rat*(rth(j+1,i)-rth(j,i))
          if (rv1.gt.rv2) goto 20
          iv=iv+1
 20     continue
        if (iv.ge.ncon) iv=ncon-1
        if (iv.eq.0) iv=1
        dk=kval/(ncon-1.)
        if (kval.eq.0.) dk=1./(ncon-1.)
        rv1=rth(j,iv)+rat*(rth(j+1,iv)-rth(j,iv))
        ratr=(rv2-rv1)/(rth(j,iv+1)+rat*(rth(j+1,iv+1)-rth(j,iv+1))-rv1)
        if ((ratr.gt.1.).or.(ratr.lt.0.)) then
          write(6,*)' Error in ratr in bdry, rat=',ratr
          write(6,*)' iv=',iv,' rv1=',rv1,' rv2=',rv2
          stop
        end if
        kscl=dk*(iv+ratr-1.)
      end if
      if (kval.eq.0.) then
        if (iout.eq.1) kscl=2.-kscl
        f=1.-kscl
      else
        if (iout.eq.1) kscl=2.*kval-kscl
        f=kval-kscl
      end if
   end subroutine bdry

!
!*****************************************************************
