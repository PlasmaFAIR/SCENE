
      subroutine equil(niter)

!----------------------------------------------------------------------
!  Grad-Shafranov eqn solver
!----------------------------------------------------------------------
!
     use param
     implicit none
     integer i,j,k,n0,iax,jax,jm,nr1,nz1,ncon1,ij,ifail,npt,niter
     double precision x,y,f,rm1,sf,dpsi
     double precision ur,ul,fd,fdd,rpeak,rl,upeak,upr
     double precision slst,r0lst,erres,bp,z0,press
     double precision px,py,pf,rat
     double precision, dimension(:), allocatable:: psiold1,gst1
     double precision, dimension(:), allocatable:: r_old,z_old
     double precision, dimension(:,:), allocatable:: u1,grads
     double precision, dimension(:), allocatable:: xn,yn,fna,triang
     integer, allocatable:: LCELL(:,:), LNEXT(:)
     double precision:: XMIN, YMIN, DX, DY, RMAX
     double precision, allocatable:: RW(:), A(:,:)
     double precision, external:: CS2VAL
     integer:: NCC, NRR, NWW
!
!  If icont=-3 then we are iterating on ff' => do not need initialisation
!  stuff
      if (icont.eq.-3) goto 30
!  Initialize plasma edge routine...
      x=0.
      y=0.
      f=0.
      call bdry(x,y,f)
! ------------------------------------------------------------
!           initialise ixout,idout arrays
      call initze
!
!  generate initial guess for U:
      do i=1,nr
        do j=1,nz
	  call bdry(r(i),z(j),f)
	  u(i,j)=f
        end do
      end do
!  if ipass=0 uses input parameterisation of ff' (ie first run)
!  if ipass=1 (set in ffdgen) uses mesh ff' calculated in ffdgen
      ipass=0
 30   continue
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
            do i=1,nr1
              do j=1,nz1
                read(12,*)u1(i,j)
              end do
            end do
            do i=1,ncon1
              read(12,*)psiold1(i),gst1(i)
            end do
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
 10               continue
                end if
                if (ij.eq.ncon1) ij=ij-1
                rat=(psiold(i)-psiold1(ij))/(psiold1(ij+1)-psiold1(ij))
                gst(i)=gst1(ij)+rat*(gst1(ij+1)-gst1(ij))
              end do
            else
              do i=1,ncon
                psiold(i)=psiold1(i)
                gst(i)=gst1(i)
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
                  py=0.
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
 50     if((rm1.ge.1.0).and.(rm1.lt.10.0))go to 70
        if(rm1.lt.1.0)go to 60
        rm1=rm1/10.0
        sf=sf*10.0
        goto 50
 60     rm1=rm1*10.0
        sf=sf/10.0
        goto 50
 70     continue
        do i=1,nr
          do j=1,nz
            call bdry(r(i),z(j),f)
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
          end if
        end do
      end do
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
 110    continue
 111    continue
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
          write(6,135)errit,erres
	  if (erres.lt.errit) goto 130
	end if
! store old values of s and r0 before iterating equilibrium
	slst=scl
	r0lst=r0
 120  continue
!
!*****************************************************************
!
      subroutine bdry(x,y,f)
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
      save
      data igo/0/
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
      if (ipr.eq.0) then
      write(nw,10) r1,r2,rs,zs
 10   format('d-shape boundary)',' r1,r2,rs,zs=',4e12.5)
      end if
      igo=1
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
!!!$!**********************************************************************
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
      double precision fp,h,psip,fr,afr,psib,psiq,afp,aq,hb
      integer itots,iv,jv,is,js,k,isum,ii,ivs,jvs,ivsp,jvsp
!
      call bdry(r(iv),z(jv),fp)
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
         call bdry(r(ivs),z(jvs),fr)
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
      double precision uf,dist,del,resid,rhs
      double precision ur,ul,fd,fdd,rpeak,rl,upr,upeak
      integer ni,i,j,iax,jax,nn
!
!      dz=dr always
      if (dr.ne.dz) then
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
