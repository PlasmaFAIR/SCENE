module equilibrium
  implicit none

contains

subroutine equil(niter)

  !----------------------------------------------------------------------
  !  Grad-Shafranov eqn solver
  !----------------------------------------------------------------------
  !
  use param
  use toms790, only : CSHEP2, CS2VAL
  implicit none
  integer i,j,k,n0,iax,jax,jm,nr1,nz1,ncon1,ij,ifail,npt,niter,igo
  double precision x,y,f,rm1,sf,dpsi
  double precision ur,ul,fd,fdd,rpeak,rl,upeak,upr
  double precision slst,r0lst,erres
  double precision px,py,pf,rat, Rmin
  double precision, dimension(:), allocatable:: psiold1,gst1
  double precision, dimension(:), allocatable:: r_old,z_old
  double precision, dimension(:,:), allocatable:: u1,grads
  double precision, dimension(:), allocatable:: xn,yn,fna,triang
  !
  integer, allocatable:: LCELL(:,:), LNEXT(:)
  double precision:: XMIN, YMIN, DX, DY, RMAX
  double precision, allocatable:: RW(:), A(:,:)
  integer:: NCC, NRR, NWW

  !
  !  If icont=-3 then we are iterating on ff' => do not need initialisation
  !  stuff
  if (icont.eq.-3) goto 30
  !  Initialize plasma edge routine...
  x=0
  y=0
  f=0.
  if (ibdry.eq.0) then
     igo=0
     call bdry(x,y,f,0)
  else
     igo=0
     call bdry2(x,y,f,0)

     if (igr .ge. 1 .and. ibdry .eq. 2) then
        !Display fourier of boundary
        call execute_command_line("ipython3 ~/SCENEv2/graphs/bdytest.py")

     end if
  end if
  ! ------------------------------------------------------------
  !           initialise ixout,idout arrays

  call initze

  write(6,*) z(1), z(10), z(nz), r(1), r(10), r(nr)
  !
  !  generate initial guess for U:
  do i=1,nr
     do j=1,nz
        if (ibdry.eq.0) then
           call bdry(r(i),z(j),f,1)
        else
           call bdry2(r(i),z(j),f,1)
        end if
        u(i,j)=f
     end do
  end do


  !  if ipass=0 uses input parameterisation of ff' (ie first run)
  !  if ipass=1 (set in ffdgen) uses mesh ff' calculated in ffdgen
  ipass=0
30 continue
  if (icont.eq.-3) ipass=1
  ! new start or continuation
  if (icont.lt.0) then
     !  continuation .....
     if (icont.eq.-1) then
        if (ipass.eq.0) then
           if (ipr.eq.0) write(nw,*)' continuation (input data)....'
           !  read in stored equilibrium
           read(7)title,scl,r0,u
           rewind 7

        end if
     else
        if (ipass.eq.0) then
           !  read in stored equilibrium
           if (ipr.eq.0) write(nw,*)' continuation (last equilibrium).....'
           read(12,*)ncon1,nr1,nz1
           read(12,*)scl,r0
           allocate(u1(nr1,nz1))
           allocate( r_old(nr1))
           allocate( z_old(nz1))
           allocate(psiold1(ncon1),gst1(ncon1))
           do i=1,nr1
              read(12,*)r_old(i)
           end do
           do j=1,nz1
              read(12,*)z_old(j)
           end do
           write(6,*) 'read in r and z successfully'
           do i=1,nr1
              do j=1,nz1
                 read(12,*)u1(i,j)
              end do
           end do
           write(6,*) 'read in u successfully'
           do i=1,ncon1
              read(12,*)psiold1(i),gst1(i)
           end do
           write(6,*) 'read in psi ff'' successfully'
           if (ncon1.ne.ncon) then

              !  Need to interpolate onto new mesh...make equi-distant in psi
              dpsi=1./(ncon-1.)
              do i=1,ncon
                 psiold(i)=1.-(i-1)*dpsi
                 if (i.eq.1) then
                    ij=1
                 else
                    ij=1
                    do 10 j=1,ncon1
                       if (psiold(i).gt.psiold1(j)) goto 10
                       ij=j
10                     continue
                    end if
                    if (ij.eq.ncon1) ij=ij-1
                    rat=(psiold(i)-psiold1(ij))/(psiold1(ij+1)-psiold1(ij))
                    gst(i)=gst1(ij)+rat*(gst1(ij+1)-gst1(ij))
                 end do
              else
              do i=1,ncon
                 psiold(i)=psiold1(i)
                 gst(i)=gst1(i)
                 write(6,*)  gst(i)
             end do
           end if
           !  And interpolate onto new R-Z mesh
           write(6,*)' nr1=',nr1,' nr=',nr,' nz1=',nz1,' nz=',nz
           if ((nr.eq.nr1).and.(nz.eq.nz1)) then
              write(6,*)' Same R-Z mesh as stored..'
              do i=1,nr
                 do j=1,nz
                    u(i,j)=u1(i,j)
                 end do
              end do
           else
              npt=nr1*nz1
              allocate( xn(npt),yn(npt),fna(npt),grads(2,npt),triang(7*npt) )
              grads=0.
              triang=0.
              k=0
              do  i=1,nr1
                 do  j=1,nz1
                    k=k+1
                    xn(k)=r_old(i)
                    yn(k)=z_old(j)
                    fna(k)=u1(i,j)
                 end do
              end do
              !              ifail=0
              ! Setup for interpolation
              NCC=MIN(17,npt-1)
              NWW=MIN(30,npt-1)
              NRR=ceiling( sqrt(dble(npt)/dble(3.0)) )
              allocate(LCELL(NRR,NRR))
              allocate(LNEXT(npt))
              allocate(RW(npt))
              allocate(A(9,npt))
              !              call e01saf(npt,xn,yn,fna,triang,grads,ifail)
              call CSHEP2(npt,xn,yn,fna,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
              if(ifail.NE.0) then
                 write(6,*) 'CSHEP2 ifail not zero!'
                 stop
              end if
              ! End of interpolation setup
              px=r0
              py=0.0d0
              ifail=1
              !                  call e01sbf(npt,xn,yn,fna,triang,grads,px,py,pf,ifail)
              pf=CS2VAL(px,py,npt,xn,yn,fna,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
              do i=1,nr
                 do j=1,nz
                    px=r(i)
                    py=z(j)
                    ifail=1
                    !                  call e01sbf(npt,xn,yn,fna,triang,grads,px,py,pf,ifail)
                    pf=CS2VAL(px,py,npt,xn,yn,fna,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
                    u(i,j)=pf
                    if (ixout(i,j).eq.1) then
                       if (u(i,j).lt.0.) u(i,j)=-u(i,j)
                    else
                       if (u(i,j).gt.0.) u(i,j)=-u(i,j)
                    end if
                 end do
              end do
              deallocate( xn,yn,fna,grads,triang )
              rewind 12

           end if
           ipass=1
        end if
     end if

  else
     if (ipr.eq.0) write(nw,*)' new start'
     !  new start .....
     ! initialise u to be value of bdry scaled.
     ! scale to correct order of magnitude
     rm1=rcen
     sf=1.0
50   if((rm1.ge.1.0).and.(rm1.lt.10.0))go to 70
     if(rm1.lt.1.0)go to 60
     rm1=rm1/10.0
     sf=sf*10.0
     goto 50
60   rm1=rm1*10.0
     sf=sf/10.0
     goto 50
70   continue
     do i=1,nr
        do j=1,nz
           if (ibdry.eq.0) then
              call bdry(r(i),z(j),f,1)
           else
              call bdry2(r(i),z(j),f,1)
           end if
           u(i,j)=f*sf
           if((ixout(i,j).eq.0).or.(ixout(i,j).eq.-3))u(i,j)=0.0
        end do
     end do
  end if

  !  calculate max. value of psi for the initial guess
  umax = -100.0
  do i=1,nr
     do j=1,nz
        if(u(i,j).ge.umax) then
           umax = u(i,j)
           jm=j
           if(u(i,j).ge.umax)iax=i
           if(u(i,j).ge.umax)jax=j
           if(u(i,j).ge.umax)umax=u(i,j)
        end if
     end do
  end do
  ! find rpeak by fitting quadratic
  ur=u(iax+1,jax)
  ul=u(iax-1,jax)
  fd=(ur-ul)/(2.0*dr)
  fdd=(ur-2*umax+ul)/(dr*dr)
  rpeak=r(iax) - fd/fdd
  r0=rpeak
  rl=rpeak - r(iax)
  upeak=umax + rl*fd + rl*rl/2.0 * fdd
  upr=umax
  ! umax labels max psi in plasma (i.e. value at magnetic axis)
  umax=upeak
  if (nouter.eq.0) then
     write(nw,*)'warning***no iteration done on initial eqbm!'
     write(nw,*)'if you want iteration reset nouter>0'
     goto 130
  end if
  do 120 n0=1,nouter
     call iter

     !  re-calculate maximum psi
     umax=-100.0
     do 111 j=1,nz
        do 110 i=1,nr
           if(u(i,j).ge.umax)iax=i
           if(u(i,j).ge.umax)jax=j
           if(u(i,j).ge.umax)umax=u(i,j)
110        continue
111        continue
           ! find rpeak by fitting quadratic
           ur=u(iax+1,jax)
           ul=u(iax-1,jax)
           fd=(ur-ul)/(2.0*dr)
           fdd=(ur-2*umax+ul)/(dr*dr)
           rpeak=r(iax) - fd/fdd
           r0=rpeak
           rl=rpeak - r(iax)
           upeak=umax + rl*fd + rl*rl/2.0 * fdd
           upr=umax
           ! umax labels max psi in plasma (i.e. value at magnetic axis)
           umax=upeak
           !  scale p'and ff' to obtain requested total current
           call constr(n0,0)
           if (n0.ne.1) then
              ! measure of convergence
              erres=abs((slst-scl)/(slst+scl))+abs((r0-r0lst)/(r0+r0lst))  &
                   +abs((cur-curtot)/(cur+curtot))
              ! test convergence of equilibrium

              if (mod(n0,10) .eq.0) write(6,135)errit,erres

              if (erres.lt.errit) goto 130
           end if
           ! store old values of s and r0 before iterating equilibrium
           slst=scl
           r0lst=r0
120        continue
130        continue
135        format('Required error=',e13.5,' achieved error=',e13.5)
           !
           !      if (ipass.eq.0) then
           if (niter.eq.1) then
              write(7)title,scl,r0,u
           end if
           rewind 7
           close(7)
           write(6,*)' is it writing data file...'
           if ((niter.gt.1).or.(npass.eq.1)) then
              write(6,*)' Writing eqbm data file'
              write(12,*)ncon,nr,nz
              write(12,*)scl,r0
              do i=1,nr
                 write(12,*)r(i)
              end do
              do j=1,nz
                 write(12,*)z(j)
              end do
              do i=1,nr
                 do j=1,nz
                    write(12,*)u(i,j)
                 end do
              end do
              do i=1,ncon
                 write(12,*)psiold(i),gst(i)
              end do
              rewind 12
              close(12)
              write(6,*)' Done writing data file'
           end if
           !  extrapolate outside plasma...uses Nag routines
           call extrap2
           ! calculate the value of poloidal field at centre of each mesh point
           call valbth


           do i=1,nr
              if (ixout(i,nsym).eq.0) cycle
              Rmin = r(i)
              exit
           end do
           do i=nr,1,-1
              if(ixout(i,nsym).eq.0) cycle
              Rmax = r(i)
              exit
           end do
           amin = (Rmax - Rmin)/2.

end subroutine equil
         !
         !*****************************************************************
!
      function bp(x,y)
!     ****************
!
!  This function interpolates a value at the point (x,y) from the
!  matrix bth. The interpolation is linear
!  (modified sykes routine)
      use param
      implicit none
      double precision x,y,xl,xr,xbase,yl,yr,ybase,bp,br,bz
      integer ii,i,jj,j
!
!  find the first x coord. not less than x
!
      i=0
      do 10 ii=1,nr-1
        if (x.lt.rcoord(ii)) goto 10
        i=ii
 10   continue
      if (i.eq.0) then
        write(nw,*)' ERROR***Asked for Btheta too close to box edge'
        write(nw,*)' x=',x,' rcoord(1)=',rcoord(1)
        stop
      end if
      xl=x-rcoord(i)
      xr=rcoord(i+1)-x
      xbase=xr+xl
!
!  find first y coord. not less than y
!
      do 20 jj=1,nz-1
        if(y.lt.zcoord(jj)) goto 20
        j=jj
 20   continue
      if (j.eq.0) then
        write(nw,*)' ERROR***Asked for Btheta too close to box edge'
        write(nw,*)' y=',y,' zcoord(1)=',zcoord(1)
        stop
      end if
      yl=y-zcoord(j)
      yr=zcoord(j+1)-y
      ybase=yr+yl
!
!  bilinear interpolation
!
!      bp=(yl*(xr*btheta(i,j+1)+xl*btheta(i+1,j+1)) &
!        +yr*(xr*btheta(i,j)+xl*btheta(i+1,j)))/(xbase*ybase)
      br=(yl*(xr*brcoord(i,j+1)+xl*brcoord(i+1,j+1)) &
        +yr*(xr*brcoord(i,j)+xl*brcoord(i+1,j)))/(xbase*ybase)
      bz=(yl*(xr*bzcoord(i,j+1)+xl*bzcoord(i+1,j+1)) &
        +yr*(xr*bzcoord(i,j)+xl*bzcoord(i+1,j)))/(xbase*ybase)
      bp=sqrt(br**2+bz**2)
!
   end function bp
!
!*****************************************************************
!
      subroutine bdry2(x,y,f,igo)
!     **********************
!
! Sets up the boundary of the flux surface
!
      use param
      implicit none
      double precision errg,dk,dth,kscl,th,rnor,znor,sg,s2g,rat
      double precision rv1,rv2,f1,f2,rup,rlo
      double precision x,y,f
      double precision xv(npts),yv(npts),gth(npts),ang(npts),rtemp(npts)
      double precision work(npts),ang2(npts)
      integer nit,i,j,ipit,jx,igo,iv,nh,k,izmid,nk,nextra,nreadpts
      integer izr2,ndmr,ndat
      double precision, dimension(:), allocatable:: rdim,zdim,rtmp,ztmp,thdim
      double precision intk,arg1,arg2,zval,rval,rofth
      character(len=30) string2
      character(len=2) string1
      !

!      write(6,*)' in bdry2, ibdry=',ibdry,' igo=',igo
      if (ibdry.eq.1) then
      if (igo.eq.0) then
!  set up r(theta,psi)
        errg=1.0d-6
        nit=20
        jx=npts/4+1
        allocate(theta(npts),rth(npts,ncon))
        theta=0.; xv=0.; yv=0.; gth=0.; rth=0.; ang2=0.
!!$!  put X-points in using modified Bishop formula
        dk=kval/(ncon-1)
        if (abs(kval).lt.1e-8) dk=1./(ncon-1.)
        dth=2.*pi/(npts-1.)
        ang(1)=-pi
        theta(1)=-pi
        do j=2,npts
          ang(j)=ang(1)+(j-1)*dth
!          write(91,*)' j=',j,' ang=',ang(j),' dth=',dth
          theta(j)=theta(1)+(j-1)*dth
        end do
        do i=ncon,2,-1
          if (abs(kval).lt.1e-8) then
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
      if ((rat.gt.1.).or.(rat.lt.-1.0d-10)) then
        write(6,*)' Error in rat in bdry2, rat=',rat
        write(6,*)' j=',j,' th=',th,' thetaj=',theta(j),' theta(j+1)=',theta(j+1)
        stop
      end if
      iv=0
      rv1=rth(j,ncon)+rat*(rth(j+1,ncon)-rth(j,ncon))
      rv2=rup
      kscl=kval*rup/rv1
      if (abs(kscl).lt.1e-8) kscl=rup/rv1
      if (abs(kval).lt.1e-8) then
!        if (iout.eq.1) kscl=2.-kscl
        f=1.-kscl
      else
!        if (iout.eq.1) kscl=2.*kval-kscl
        f=kval-kscl
      end if
!--------------------------------------------------------------
      else if (ibdry.eq.2) then
        if (igo.eq.0) then
          allocate(theta(npts))
          dth=2.*pi/(npts-1)
          do i=1,npts
            theta(i)=(i-1)*dth
          end do
!  Set up fourier components of plasma surface
          if (nfm.le.0) then
            write(6,*)' Must set no. Fourier modes >0, nfm=',nfm
            stop
          end if
          allocate(rfo(nfm),zfo(nfm))
! Read in free boundary: data stored anti-clockwise, starting on outboard mid-plane (theta=0). Last point is same as first (ie stored on 0</=theta</=2pi
          nh=39
          open(unit=nh,file='bdy.txt', &
                status='unknown',iostat=ios)
          if(ios.ne.0) then
            write(6,*) 'problem opening bdy.txt for boundary R,Z'
            stop
          endif
          read(nh,*)ndat
          allocate(rdim(ndat),zdim(ndat),rtmp(npts),ztmp(npts))
          do i=1,ndat
            read(nh,*)rdim(i),zdim(i)
          end do
          close(nh)
!  Find max and min R to define geometric axis
          r1=rdim(1)
          izmid=1
          r2=rdim(1)
          izr2=1
          do i=2,ndat
            if (r1.gt.rdim(i)) then
              r1=rdim(i)
              izmid=i
            end if
            if (r2.lt.rdim(i)) then
              r2=rdim(i)
              izr2=i
            end if
         end do
         print*, 'R1 and R2', r1, r2
          rcen=(r1+r2)/2.
          write(6,*)' rcen=',rcen,' izmid-',izmid,' izr2=',izr2
!          write(6,*)' zdim(1)=',zdim(1),' rdim(1)=',rdim(1)
          do i=1,ndat
            rdim(i)=rdim(i)-rcen
         end do

!---------------------------------------
          allocate (thdim(ndat),fm(nfm))
          do i=1,ndat-1
            thdim(i)=atan2(zdim(i),rdim(i))
            if (thdim(i).lt.0.) thdim(i)=thdim(i)+2.*pi
          end do
          thdim(ndat)=thdim(1)+2.*pi
!---------------------------------------
          do i=1,nfm
            k=i-1
            intk=0.
            do j=2,ndat
              dth=thdim(j)-thdim(j-1)
              arg1=sqrt(rdim(j-1)**2+zdim(j-1)**2)*cos(k*thdim(j-1))
              arg2=sqrt(rdim(j)**2+zdim(j)**2)*cos(k*thdim(j))
              intk=intk+0.5*(arg1+arg2)*dth
            end do
            if (k.eq.0) then
              fm(i)=intk/(2.*pi)
            else
              fm(i)=intk/pi
            end if
!            write(6,*)' i=',i,' fm=',fm(i)
          end do
          do j=1,npts
            th=theta(j)
            rofth=0.
            do i=1,nfm
              rofth=rofth+fm(i)*cos((i-1)*th)
            end do
            rtmp(j)=rofth*cos(th)
            ztmp(j)=rofth*sin(th)
          end do
! Evaluate the parameters related to plasma shape:
          r2=0.
          r1=0.
          do i=1,nfm
            r2=r2+fm(i)
            r1=r1+fm(i)*cos((i-1)*pi)*cos(pi)
          end do
          r2=rcen+r2
          r1=rcen+r1
          rcen=(r1+r2)/2.
          tokeps=(r2-r1)/(r1+r2)
          dth=(2.*pi)/(npts-1.)
          zs=0.
          rs=rcen
          do i=1,npts
            th=theta(i)
            zval=0.
            rval=rcen
            do j=1,nfm
              zval=zval+fm(j)*cos((j-1)*th)*sin(th)
              rval=rval+fm(j)*cos((j-1)*th)*cos(th)
            end do
            if (zs.lt.abs(zval)) then
              zs=abs(zval)
              rs=rval
            end if
          end do
          elon=2.*zs/(r2-r1)
          tri=(1.-rs/rcen)/tokeps


          !Replaces tstplt3
          write(6,*)' tokeps=',tokeps,' rcen=',rcen,' elon=',elon,' tri=',tri

          nh = nh+2

          open(unit=nh, file='bdytest.dat', status='unknown', iostat=ios)
          if (ios .ne. 0) then
             write(6,*) 'Failed to open bdytest in equil.f90'
             stop
          end if

          write(nh,*) npts

          do i=1,npts
             write(nh,*) theta(i), ztmp(i), rtmp(i)
          end do

          do i=1,ndat
             write(nh,*) thdim(i)
          end do

          close(nh)

          return
        else
!  Evaluate function
          th=atan2(y,x-rcen)
          rofth=0.
          do i=1,nfm
            rofth=rofth+fm(i)*cos((i-1)*th)
          end do
          f=1.-sqrt((x-rcen)**2+y**2)/rofth
          return
        end if
!--------------------------------------------------------------
      else
!        write(6,*)' in bdry Using Dmitrij's input file**********'
        if (igo.eq.0) then
          allocate(theta(npts))
!  Set up fourier components of plasma surface
          if (nfm.le.0) then
            write(6,*)' Must set no. Fourier modes >0, nfm=',nfm
            stop
          end if
          allocate(fm(nfm))
! Read in free boundary
          nh=39
          open(unit=nh,file='ctf07-bdry.dat', &
                status='unknown',iostat=ios)
          if(ios.ne.0) then
            write(6,*) 'problem opening ctf07-bdry.dat for boundary R,Z'
            stop
          endif
          do i=1,23
            read(nh,*)string1
          end do
      read(nh,*)string2,string2,string2,string2,string2,string2,rdim(1),zdim(1)
!          zdim(1)=-zdim(1)
          do i=1,232
            read(nh,*)string1
          end do
          ndmr=51
          do i=1,ndmr
            read(nh,*)zdim(i+1),rdim(i+1),rdim(2*ndmr+3-i)
          end do
          close(nh)
          nreadpts=2*ndmr+2
          zdim(ndmr+2)=-zdim(1)
          rdim(ndmr+2)=rdim(1)
          do i=1,ndmr
            zdim(ndmr+2+i)=zdim(ndmr-i+2)
          end do
!  Find max and min R to define geometric axis
          r1=rdim(1)
          izmid=1
          r2=rdim(1)
          do i=1,nreadpts
            if (r1.gt.rdim(i)) then
              r1=rdim(i)
              izmid=i
            end if
            if (r2.lt.rdim(i)) r2=rdim(i)
          end do
          rcen=(r1+r2)/2.
          do i=1,izmid
            rtmp(i)=rdim(i)
            ztmp(i)=zdim(i)
          end do
          do i=1,nreadpts-izmid+1
            rdim(i)=rdim(i+izmid-1)
            zdim(i)=zdim(i+izmid-1)
          end do
          j=0
          do i=nreadpts-izmid+2,nreadpts
            j=j+1
            rdim(i)=rtmp(j)
            zdim(i)=ztmp(j)
          end do
          do i=1,nreadpts
            rdim(i)=rdim(i)-rcen
            write(6,*)' i=',i,' R=',rdim(i),' Z=',zdim(i)
          end do
!  Load additional points:
!  Load additional points:
          k=1
          nextra=1
          rtmp(1)=rdim(1)
          ztmp(1)=zdim(1)
          do i=1,nreadpts-1
            do j=1,nextra
              k=k+1
              rtmp(k)=rdim(i)+j*(rdim(i+1)-rdim(i))/nextra
              ztmp(k)=zdim(i)+j*(zdim(i+1)-zdim(i))/nextra
            end do
          end do
          nk=k
!  Temporarily use theta for the input mesh
          do i=1,nk
            rdim(i)=rtmp(i)
            zdim(i)=ztmp(i)
            write(6,*)' k=',i,' R=',rdim(i),' Z=',zdim(i)
            theta(i)=atan2(zdim(i),rdim(i))
            write(6,*)' theta=',theta(i)
         end do

         nh=nh+2

         !replaces tstplt
          open(unit=nh, file='bdytest.dat', status='unknown', iostat=ios)
          if (ios .ne. 0) then
             write(6,*) 'Failed to open bdytest in equil.f90'
             stop
          end if

          write(nh,*) nk

          do i=1,npts
             write(nh,*) theta(i)
          end do

          close(nh)

!          if (igo.eq.0) stop
          do i=1,nfm
            k=i-1
            dth=theta(2)-theta(1)
            arg1=sqrt(rdim(1)**2+zdim(1)**2)*cos(k*theta(1))
            arg2=sqrt(rdim(2)**2+zdim(2)**2)*cos(k*theta(2))
            intk=0.5*(arg1+arg2)*dth
            do j=2,nk
              dth=theta(j)-theta(j-1)
              arg1=sqrt(rdim(j-1)**2+zdim(j-1)**2)*cos(k*theta(j-1))
              arg2=sqrt(rdim(j)**2+zdim(j)**2)*cos(k*theta(j))
              intk=intk+0.5*(arg1+arg2)*dth
            end do
            if (k.eq.0) then
              fm(i)=-intk/(2.*pi)
            else
              fm(i)=-intk/pi
            end if
          end do
          do j=1,nk
            th=theta(j)
            write(6,*)' j=',j,' th=',th
            rofth=0.
            do i=1,nfm
              rofth=rofth+fm(i)*cos((i-1)*th)
            end do
            rtmp(j)=rofth*cos(th)
            ztmp(j)=rofth*sin(th)
          end do
!          call tstplt2(nk,theta,rdim,rtmp)
!          call tstplt2(nk,theta,zdim,ztmp)
! Evaluate the parameters related to plasma shape:
          r2=0.
          r1=0.
          do i=1,nfm
            r2=r2+fm(i)
            r1=r1+fm(i)*cos((i-1)*pi)*cos(pi)
          end do
          r2=rcen+r2
          r1=rcen+r1
          rcen=(r1+r2)/2.
          tokeps=(r2-r1)/(r1+r2)
          dth=(2.*pi)/(npts-1.)
          zs=0.
          rs=rcen
          do i=1,npts
            theta(i)=-pi+(i-1)*dth
            th=theta(i)
            zval=0.
            rval=rcen
            do j=1,nfm
              zval=zval+fm(j)*cos((j-1)*th)*sin(th)
              rval=rval+fm(j)*cos((j-1)*th)*cos(th)
            end do
            if (zs.lt.abs(zval)) then
              zs=abs(zval)
              rs=rval
            end if
          end do
          elon=2.*zs/(r2-r1)
          tri=(1.-rs/rcen)/tokeps
          if (igo.eq.0) then
            do i=1,nfm
              write(6,*)' k=',i-1,' fm=',fm(i)
            end do
            write(6,*)' rcen=',rcen,' r1=',r1,' r2=',r2
            write(6,*)' zs=',zs,' rs=',rs
            write(6,*)' tokeps=',tokeps,' elon=',elon,' tri=',tri
            stop
          end if
        else
!  Evaluate function
          th=atan2(y,x-rcen)
          rofth=0.
          do i=1,nfm
            rofth=rofth+fm(i)*cos((i-1)*th)
          end do
          f=1.-sqrt((x-rcen)**2+y**2)/rofth
          return
        end if
     end if



!  ibdry not equal to 1 or 0=> use fourier description.
   end subroutine bdry2

!
!*****************************************************************
!
      subroutine bdry(x,y,f,igo)
!     **********************
!
! Sets up the boundary of the flux surface
!
      use param
      implicit none
      double precision rzero,rsb,rg,zdsq,zsbsq
      double precision al,ga,rb,zbsq
      double precision x,y,f
      integer igo
      logical debug
      save

      debug = .false.
!      data igo/0/
!
      if(igo.eq.1) go to 20
!  set up plasma boundary
      rg=0.5*(r1+r2)
      rzero=dsqrt(0.5*(r1*r1+r2*r2))
      rsb=rs/rzero
      al=((r2*r2+r1*r1)/(r2*r2-r1*r1))**2
      zdsq=sngl(zs*zs/(2.0*al*(1.0-rsb*rsb)))
      zsbsq=zs*zs/zdsq
      ga=(al*(rsb*rsb-1.0)**2 -1.0)/zsbsq
      ga=ga+rsb*rsb
      if (debug) write(6,*)' inside bdr and printing'
      if (ipr.eq.0) then
      write(nw,10) r1,r2,rs,zs
 10   format('d-shape boundary)',' r1,r2,rs,zs=',4e12.5)
      end if
!      write(6,*)' done writing'
!      igo=1
!      write(6,*)' returning from bdry'
      return
 20   rb=x/rzero
      zbsq=y*y/zdsq
      f=1.0-zbsq*(rb*rb-ga)-al*(rb*rb-1.0)**2
!
! special mod..
      if (x.lt.0.0) f=-1.0
  end subroutine bdry
!
!*****************************************************************
!
      subroutine constr(n0,ip)
!     ************************
!
! Scales scl to give correct total current
!
      use param
      implicit none
      double precision cval,ptot
      integer i,j,ip,n0
!
! evaluate total current
      curtot=0.0
      ptot=0.0
      !vv=0.
      do i=1,nr
    do 10 j=1,nz
      if(ixout(i,j).ne.1)goto 10
   cval=-rhs(u(i,j),r(i))/(r(i)*mu0)
!   read(6,*) test
   curtot=curtot+dr*dz*cval

!          psi=umax-u(i,j)
!          ptot=ptot+r(i)*press(psi,0)*dr*dz
!          vv=vv+r(i)*dr*dz
 10     continue
      end do
!      bv=mu0*rodi/(2.*pi*rcen)
!      betnow=(200.*mu0*ptot/(vv*bv**2))*tokeps*rcen*bv*1.0d6/cur
!!$!
!!$! rescale eigenvalue
!!$      if (betan.gt.0.) then
!!$!        if (abs((curtot-cur)/(curtot+cur)).lt.10.0*errit) then
!!$          write(6,*)' curtot=',curtot,' cur=',cur,' ratio=',curtot/cur
!!$          bpol=bpol*betan*curtot/(betnow*cur)
!!$!        end if
!!$      end if
      if(curtot.gt.(0.0))scl = scl*cur/curtot
      if (ip.eq.0) then
         if (ipr.eq.0) then
            if ( mod(n0,10) .eq. 0) then
      write(nw,30)n0,cur,curtot,scl,residl,r0
30    format(' it=',i3,' curr req=',e12.4,' cur got=',e12.4,     &
        ' eig s=',e12.4,' resid=',e12.4,' r0=',f8.4)
!      write(6,*)' bpol=',bpol,' betnow=',betnow
      end if
      !
      end if
      end if
      if (ipr.eq.0) then
      if(curtot.le.0.0)write(nw,40)
 40   format(' total curtot -ve..carry on,if remains -ve..',/,  &
       '  reduce chucking or swopping')
      end if
  end subroutine constr
!
!**********************************************************************
!
      subroutine extrap2
!     ******************
!
      use param
      use toms790, only : CSHEP2, CS2VAL
      implicit none
      double precision, dimension(:), allocatable::  x,y,f
      double precision, dimension(:,:), allocatable:: grads
      double precision px,py,pf
      integer, dimension(:), allocatable:: triang
      integer npt,k,i,j,ifail
      integer, allocatable:: LCELL(:,:), LNEXT(:)
      double precision:: XMIN, YMIN, DX, DY, RMAX
      double precision, allocatable:: RW(:), A(:,:)
      integer:: NCC, NRR, NWW
      logical :: debug

      debug = .false.
      
      if (debug) write(6,*) 'running extrap2'
!
!  Extrapolation routine to fill out the psi mesh on the R-Z grid
!  from the points calculated from the G-S solver
!
!
!  Load up points inside and just outside plasma
      npt=0
      do i=1,nr
        do j=1,nz
          if (ixout(i,j).ne.0) npt=npt+1
        end do
      end do
      allocate( x(npt),y(npt),f(npt),grads(2,npt),triang(7*npt) )
      k=0
      do  i=1,nr
        do  j=1,nz
          if (ixout(i,j).ne.0) then
            k=k+1
            x(k)=r(i)
            y(k)=z(j)
            f(k)=u(i,j)
          end if
        end do
      end do
!      ifail=0
!      write(6,*)' in e01saf'
      ! Setup for interpolation
      NCC=MIN(17,npt-1)
      NWW=MIN(30,npt-1)
      NRR=ceiling( sqrt(dble(npt)/dble(3.0)) )
      allocate(LCELL(NRR,NRR))
      allocate(LNEXT(npt))
      allocate(RW(npt))
      allocate(A(9,npt))
!      call e01saf(npt,x,y,f,triang,grads,ifail)
      call CSHEP2(npt,x,y,f,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
!      write(6,*)' out e01saf'
      do i=1,nr
        do j=1,nz
          if (ixout(i,j).eq.0) then
            px=r(i)
            py=z(j)
            ifail=1
!      write(6,*)' in e01saf'
!            call e01sbf(npt,x,y,f,triang,grads,px,py,pf,ifail)
!      write(6,*)' out e01saf'
            pf=CS2VAL(px,py,npt,x,y,f,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
            u(i,j)=pf
          end if
        end do
      end do
      deallocate( x,y,f,grads,triang )
  end subroutine extrap2
!
!!$!**********************************************************************
!
      subroutine extrap(iv,jv)
!     ************************
!
!
! special interpolation and extrapolation routine to get psi at
! -1 & -3 pts using a linear extrapolation technique
! revised 16/11/76   nearest 1 pts used to interpolate value of
! psi at boundary then next mesh-pt out used to extrapolate
! values to -1 pts. boundary then is never nearer than 1 mesh-pt
! from the pt used
!  sykes routine
!    n.b.   4/2/77   U=0.0 at bdry
!
      use param
      implicit none
      double precision fp,h,psip,fr,afr,psib,psiq,afp,hb
      integer itots,iv,jv,is,js,k,isum,ii,ivs,jvs,ivsp,jvsp
!
      if (ibdry.eq.0) then
        call bdry(r(iv),z(jv),fp,1)
      else
        call bdry2(r(iv),z(jv),fp,1)
      end if
     h=dr
     u(iv,jv)=0.0
     itots=0
     is=1
     js=1
!
!  p is the pt at which psi is required, r is the pt just inside
!  the plasma, and q is 1 mesh-pt further inside
!     move round all 8 neighbouring pts on subsequent its.
      do 10 k=1,8
        if(k.eq.2)js=0
        if(k.eq.3)js=-1
        if(k.eq.4)is=0
        if(k.eq.5)is=-1
        if(k.eq.6)js=0
        if(k.eq.7)js=1
        if(k.eq.8)is=0
        h=dr
        isum=abs(is)+abs(js)
        if (isum.eq.2) h=sqrt(2.)*dr
        ii=0
        ivs=iv+is
        if((ivs.gt.nr).or.(ivs.lt.1))go to 10
        jvs=jv+js
        if((jvs.gt.nz).or.(jvs.lt.1))go to 10
! if any n,ne,e,se,s,sw,w,nw pts are in plas then eval bdry f'n
! provided true bdry is not too close to any mesh pt then it is
! suitable for interp. etc take average of all such points.
        if(ixout(ivs,jvs).ne.1)go to 5
     ii=1
     ivsp=ivs+is
     jvsp=jvs+js
         psip=0.0
         if (ibdry.eq.0) then
           call bdry(r(ivs),z(jvs),fr,1)
         else
           call bdry2(r(ivs),z(jvs),fr,1)
         end if
         afr=abs(fr)
         psib=0.0
! interpolate to get hb (dist of boundary from r)
     afp=abs(fp)
     hb=afr*h/(afr+afp)
     psiq=u(ivsp,jvsp)
     psip=((hb-h)*psiq+2.0*h*psib)/(h+hb)
    5   continue
     u(iv,jv)=u(iv,jv)+ii*psip
         itots=itots+ii
 10   continue
!
      u(iv,jv)=u(iv,jv)/itots
   end subroutine extrap
!
!***************************************************************
!
    subroutine initze
!   ******************
!
!  sets up ixout and idout arrays which define the area occupied
! by the plasma (=1 for plasma, =-1,3 for plasma edge, =0 outside
! plasma.
!
      use param
      implicit none
      integer i,j,nrp1,nzp1
      integer im1,ip1,jm1,jp1,isum
      integer iprod
      double precision f,rb,zb
      integer sumx,sumd
!
      allocate( nix(nr*nz), niy(nr*nz) )
      nix=0; niy=0
      nrp1=nr-1
      nzp1=nz-1
! ------ first scan mesh for in/out points---------
      do i=1,nr
        do j=1,nz
          if (ibdry.eq.0) then
            call bdry(r(i),z(j),f,1)
          else
            call bdry2(r(i),z(j),f,1)
          end if
      if(f.gt.0.0)ixout(i,j)=1
      if(f.le.0.0)ixout(i,j)=0
        end do
      end do
!
! ------now replace just-outsiders with -1 or -3------
      do i=1,nr
        im1=i-1
        if(i.eq.1)im1=1
        ip1=i+1
        if (i.eq.nr) ip1=nr
    do 30 j=1,nz
      idout(i,j)=1
      if (ixout(i,j).eq.1) go to 30
          jm1=j-1
          if (j.eq.1) jm1=1
          jp1=j+1
          if (j.eq.nz) jp1=nz
          isum=ixout(i,jm1)+ixout(i,jp1)+ixout(ip1,j)+ixout(im1,j) &
             +ixout(im1,jm1)+ixout(ip1,jp1)+ixout(ip1,jm1)+ixout(im1,jp1)
          if(isum.eq.1)idout(i,j)=-3
          if(isum.gt.1)idout(i,j)=-1
          if(isum.eq.0)idout(i,j)=0
 30     continue
      end do
!
! ------now collect just-outsiders into n arrays---------
!        also work out sines,coses at these pts
      ntot=0
      do i=1,nr
        do 50 j=1,nz
          ixout(i,j)=idout(i,j)
          idout(i,j)=0
          if(ixout(i,j).ge.0)go to 50
          ntot=ntot+1
          nix(ntot)=i
          niy(ntot)=j
 50     continue
      end do
!
!
      call iout(' ixout  ')
!
! -----set up idout array of in/out cell-centres
!       (in if none of the 4 neighs is 0)
      do i=1,nrp1
    do 80 j=1,nzp1
      iprod=ixout(i,j)*ixout(i+1,j)*ixout(i,j+1)*ixout(i+1,j+1)
      idout(i,j)=0
      if(iprod.eq.0)go to 80
! definitely need j at this enclosed cell centre..but are we
! ...at edge?
      isum=ixout(i,j)+ixout(i+1,j)+ixout(i,j+1)+ixout(i+1,j+1)
      idout(i,j)=1
      if(isum.ne.4)idout(i,j)=-1
      rb = (r(i)+r(i+1))/2.0
      zb = (z(j)+z(j+1))/2.0
          if (ibdry.eq.0) then
            call bdry(rb,zb,f,1)
          else
            call bdry2(rb,zb,f,1)
          end if
      if(f.gt.(0.0))idout(i,j) = 1
 80     continue
      end do
      call iout(' idout  ')
!  check mesh is bounded by zeros
      sumx=0
      sumd=0
      do i=1,nr
    sumx=sumx+abs(ixout(i,1))
    sumd=sumd+abs(idout(i,1))
      end do
      do i=1,nr
    sumx=sumx+abs(ixout(i,nz))
    sumd=sumd+abs(idout(i,nz))
      end do
      do i=1,nz
    sumx=sumx+abs(ixout(1,i))
    sumd=sumd+abs(idout(1,i))
      end do
      do i=1,nz
    sumx=sumx+abs(ixout(nr,i))
    sumd=sumd+abs(idout(nr,i))
      end do
      if ((sumx.ne.0).or.(sumd.ne.0)) then
    write(nw,*)'fatal error***problem with mesh!'
    write(nw,*)' plasma boundary not surrounded by zeros'
    stop
      end if
!
   end subroutine initze
!
!**********************************************************************
!
      subroutine iout(word)
!     *********************
!
! writes out the ixout and idout arrays
!
      use param
      implicit none
      integer j,jy,i
      character(len=8) word
!
      if (ipr.eq.0) write(nw,10)word
  10  format(a)
      do 30 j=1,nz
    jy=nz+1-j
        if (word.eq.' ixout  ') then
          if (ipr.eq.0) write(nw,20)(ixout(i,jy),i=1,nr)
        else
          if (ipr.eq.0) write(nw,20)(idout(i,jy),i=1,nr)
        end if
 20     format(1x,64i2)
 30   continue
!
   end subroutine iout
!
!**********************************************************************
!
      subroutine iter
!     ***************
!
! attempt to iterate on total residuals ... fct problem
! uses nonlinear sor method to reduce residuals
!  modified sykes routine
!
      use param
      implicit none
      double precision drsq,rr,riph,rimh,ub,residb
      double precision zz,uc,rhsc,cc,others,ua,resida
      double precision uf,dist,del,resid
      double precision ur,ul,fd,fdd,rpeak,rl,upr,upeak
      integer ni,i,j,iax,jax,nn
!
!      dz=dr always
      if (abs(dr-dz) .gt. 1e-4) then
        write(nw,*)' fatal error***iter routine needs dr=dz'
        write(nw,*)' must set dr=dz'
        stop
      end if
      drsq=dr*dr
! --------------
       ni=0
      do ni=1,ninner
       do  i=1,nr
         rr=r(i)
         riph=rr+0.5*dr
         rimh=rr-0.5*dr
         do 10 j=1,nz
            if(ixout(i,j).ne.1)go to 10
           zz=z(j)
           uc=u(i,j)

           rhsc=rhs(uc,rr)
!  calculate lhs of g-s eqn.
!   coeff of u at the mesh point (i,j)
           cc=-2.0*(rr/rimh*rr/riph+1.0)
!   all terms involving other mesh points
           others=u(i,j+1)+u(i,j-1)+rr*(u(i-1,j)/rimh+u(i+1,j)/riph)
! note..still need to peturb uc even if zero
           ua=0.999*uc - 0.0001
           resida=(ua*cc+others)/drsq-rhsc
           ub=1.001*uc + 0.0001
           residb=(ub*cc+others)/drsq-rhsc
! find zero by reg falsi
           if(abs(residb-resida).lt.0.00001) go to 10
           uf=(residb*ua-resida*ub)/(residb-resida)
! dont move too far
! allow over-solution...
           dist=omega*(uf-uc)
           del=frac*uc
! ... but dont let u move too far
           if(dist.le.-del)dist=-del
           if(dist.gt.del)dist=del
           u(i,j)=uc+dist
 10      continue
       end do
! --------------------------
!
!  re-calculate maximum psi
       umax=-100.0
       iax=0
       jax=0
       do  j=1,nz
         do  i=1,nr
           if(u(i,j).ge.umax)iax=i
           if(u(i,j).ge.umax)jax=j
           if(u(i,j).ge.umax)umax=u(i,j)
         end do
       end do
! find rpeak by fitting quadratic
       ur=u(iax+1,jax)
       ul=u(iax-1,jax)
       fd=(ur-ul)/(2.0*dr)
       fdd=(ur-2*umax+ul)/(dr*dr)
       rpeak=r(iax) - fd/fdd
       r0=rpeak
       rl=rpeak - r(iax)
       upeak=umax + rl*fd + rl*rl/2.0 * fdd
       upr=umax
       umax=upeak
! extrapolate at edges to get u at -1 and -3 pts (just outsiders)
! needed to advance plas at next iteration
! extrap uses fact that u=0 at boundary
       do nn=1,ntot
         i=nix(nn)
         j=niy(nn)
         call extrap(i,j)
       end do
      end do
! --------------------
! calc euclidean norm of residual vector
      residl=0.0
      do i=1,nr
        rr=r(i)
        riph=rr+0.5*dr
        rimh=rr-0.5*dr
        do 40 j=1,nz
          if(ixout(i,j).ne.1)go to 40
          zz=z(j)
          uc=u(i,j)
          cc=-2.0*(rr/rimh*rr/riph+1.0)
          others=u(i,j+1)+u(i,j-1)+rr*(u(i-1,j)/rimh+u(i+1,j)/riph)
          resid=(uc*cc+others)/drsq-rhs(uc,rr)
          residl=residl+resid*resid
 40     continue
      end do
      residl=sqrt(residl)*dr*dz
  end subroutine iter
!
!**********************************************************************
!
      subroutine range(arr,lmax,alo,aup)
!     **********************************
!
      implicit none
      real arr(*)
      real alo,aup
      integer lmax,l

!
      alo=arr(1)
      aup=arr(1)
      do 10 l=2,lmax
    if (alo.gt.arr(l)) alo=arr(l)
    if (aup.lt.arr(l)) aup=arr(l)
 10   continue
   end subroutine range
!
!**********************************************************************
!
      function rhs(uu,rv)
!     *******************
!

      !  returns rhs of g-s equation
!
      use param
      use profiles_mod, only : fprof, press
      implicit none
      double precision pd,ffp,rjphi,rhs,psi,rv,uu
!
      ! Note -sign as u=psi_a-psi
      psi=umax-uu
      pd=-press(psi,1)
      ffp=-fprof(psi,1)

      rjphi=rv*rv*pd+ffp/mu0
      rhs=-mu0*rjphi

      return
  end function rhs
!
!**********************************************************************
!
      subroutine valbth
!     *****************
!
!  Calculates the value of poloidal magnetic field at centre
!  of each square in the mesh; rcoord and zcoord label the
!  r and z coordinates of these `central' points.
!  sykes routine (modified)
!
      use param
      implicit none
      integer nx,ny,je,inum,iss,js,ins,jn,k,inn,jnn,i,j
      double precision rc,dudz,dudr,brn,bzn,btot,bn,bnn,brnn,bznn
      double precision un, unn, utot
      double precision brtot,bztot
!
      if (icont.gt.-3) then
        allocate( rcoord(nr), zcoord(nz), btheta(nr,nz),   &
                brcoord(nr,nz),bzcoord(nr,nz),ucoord(nr,nz))
      end if
      nx = nr-1
      ny = nz-1
! load in new mesh
      do i=1,nr
        rcoord(i)=r(i)+0.5*dr
      end do
      do j=1,nz
        zcoord(j)=z(j)+0.5*dz
      end do
!  calculate btheta at new mesh points
      do i=1,nx
        rc=rcoord(i)
        do j = 1,ny
          dudz=0.5*(u(i+1,j+1)-u(i+1,j)+u(i,j+1)-u(i,j))/dz
          dudr=0.5*(u(i+1,j+1)-u(i,j+1)+u(i+1,j)-u(i,j))/dr

!         #Added by Bhavin 28.09.18 for RF calc by Simon freethy 
          ucoord(i,j)= (u(i+1,j+1)+u(i,j+1)+u(i+1,j)+u(i,j))/4.

          brcoord(i,j)=-dudz/rc
          bzcoord(i,j)=dudr/rc
          btheta(i,j)=sqrt(brcoord(i,j)**2+bzcoord(i,j)**2)
        end do
      end do
!   extrapolate to outside plasma
!   first extrapolate to -1 points (i.e. at boundary)
      je = 0
 50   continue
      je = je-1
      do i=1,nx
        do 70 j = 1,ny
          if (idout(i,j).ne.je) goto 70
          inum=0
          btot=0.
          brtot=0.
          bztot=0.
          iss=1
          js=1
          do 60 k=1,8
! scan mesh points around (i,j) to find ones further into plasma)
            if (k.eq.2) js=0
            if (k.eq.3) js=-1
            if (k.eq.4) iss=0
            if (k.eq.5) iss=-1
            if (k.eq.6) js=0
            if (k.eq.7) js=1
            if (k.eq.8) iss=0
            ins=i+iss
! if outside mesh cannot include
            if ((ins.gt.nx).or.(ins.lt.1)) goto 60
            jn=j+js
! if outside mesh cannot include
            if ((jn.gt.ny).or.(jn.lt.1)) goto 60
! if away from boundary cannot include
            if (idout(ins,jn).eq.0) goto 60
! if further out do not include
            if ((idout(ins,jn).le.je).and.(je.ne.0)) goto 60
            ! found mesh point further into plasma...
            un = ucoord(ins,jn)
            bn = btheta(ins,jn)
            brn = brcoord(ins,jn)
            bzn = bzcoord(ins,jn)
! continue further in in a straight line
            inn = ins + iss
            jnn = jn + js
            unn = ucoord(inn,jnn)
            bnn = btheta(inn,jnn)
            brnn = brcoord(inn,jnn)
            bznn = bzcoord(inn,jnn)
! linear extrapolation gives b (take average of all contributions)
            btot = btot + 2.0*bn - bnn
            brtot = brtot + 2.0*brn - brnn
            bztot = bztot + 2.0*bzn - bznn
            utot = utot + 2*un - unn
            inum = inum + 1
 60       continue
      if(inum.eq.0)goto 70
          ucoord(i,j) = utot/(1.0*inum)
          btheta(i,j) = btot/(1.0*inum)
          brcoord(i,j) = brtot/(1.0*inum)
          bzcoord(i,j) = bztot/(1.0*inum)
 70     continue
      end do
      if (je.eq.0) goto 80
!  now go back and extrapolate 0's next to boundary
        je = 1
        goto 50
 80   continue
!
!
 end subroutine valbth
!
!*****************************************************************
!
      subroutine xarea(arr,xint)
!     ********************************
!
!  calculates the array arr over the cross-sectional area of the plasma
!
      use param
      implicit none
      integer i,j
      double precision arr(nr,nz)
      double precision xint
!
      xint=0.0
      do i=1,nr
    do 10 j=1,nz
      if (ixout(i,j).le.0) goto 10
      xint=xint+arr(i,j)*dr*dz
 10     continue
      end do
   end subroutine xarea
!
!*****************************************************************
!
      subroutine valbth2(rg,zg,v,nrg,nzg)
!     ****************************************
!
!  Calculates the value of poloidal magnetic field at centre
!  of each square in the mesh; rcoord and zcoord label the
!  r and z coordinates of these `central' points.
!  sykes routine (modified)
!
      use param
      implicit none
      integer nx,ny,i,j
      double precision rc,dudz,dudr
      integer nrg,nzg
      double precision rg(nrg),zg(nzg),v(nrg,nzg)
      double precision drg,dzg,brc,bzc,bth,zc
!
      drg=rg(2)-rg(1)
      dzg=zg(2)-zg(1)
      nx = nrg-1
      ny = nzg-1
! load in new mesh
!      do j=1,nz
!        zcoord(j)=z(j)+0.5*dz
!      end do
!  calculate btheta at new mesh points
      do i=1,nx
        rc=rg(i)+0.5*drg
        do j = 1,ny
          zc=zg(j)+0.5*dzg
          dudz=0.5*(v(i+1,j+1)-v(i+1,j)+v(i,j+1)-v(i,j))/dzg
          dudr=0.5*(v(i+1,j+1)-v(i,j+1)+v(i+1,j)-v(i,j))/drg
          brc=dudz/(2.*pi*rc)
          bzc=-dudr/(2.*pi*rc)
          bth=sqrt(brc**2+bzc**2)
          if (v(i,j).le. 2.*pi*umax) then
            if (j.eq.nzg/2+1) then
            write(26,*)' R=',rg(i),' Z=',zg(j),' bth=',bth,' bthold=',bp(rc,zc)
            write(26,*)' R=',rg(i),' Z=',zg(j),' psi=',v(i,j)/(2.*pi)
            end if
          end if
        end do
      end do
      do i=1,nr
        write(26,*)' R=',r(i),' Z=',z(nsym),' psi=',umax-u(i,nsym)
      end do
!
!
 end subroutine valbth2
!
!*****************************************************************
end module equilibrium
