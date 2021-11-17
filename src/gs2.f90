module gs2_output
  implicit none
contains
      subroutine gs2
!     --------------
!
!  Subroutine to calculate the various inputs required for the micro-
!  stability code GS2
!
!  HRW 22/05/03
!
!
      use equilibrium, only : bp, valbth2
      use param
      use profiles_mod, only : dense, densi, fprof, press, tempi, tempe, dpsidrho, shift
      use toms790, only : CSHEP2, CS2VAL
      use splines, only: spline, zspline
      implicit none
!
      integer :: ndsk,k,nspec,kk,ik
      double precision :: psi,zeff,zni19,ne19
      double precision :: dshafr, shafr, rmaj, rmin, epsil, beta_gs2
      double precision :: coolog, vss, ti, ni, ti_p, ni_p
      double precision :: coll, nref, tref, bcentr
      double precision :: zmag,rat,px,py,deltar,deltaz,pf,fsqedg
      double precision :: yp1,yp2
      double precision :: polvar(npts),aminor(ncon),fofpsi(ncon)
      double precision :: y2(ncon),fpint(ncon)


!      character(len=40), intent(in) :: filename
      double precision, dimension(:), allocatable :: ps,rg,zg
      double precision, dimension(:,:), allocatable :: mshvar,grads
      double precision, dimension(:), allocatable :: xn,yn,fna,triang
      double precision, dimension(:), allocatable :: temvar
      double precision, dimension(ncon) :: dpdrs,rhos

      integer :: nsurf,i,j,npt,ifail,nrg,nzg
!      type(vec2dl) :: sep
      integer, allocatable:: LCELL(:,:), LNEXT(:)
      double precision :: XMIN, YMIN, DX, DY, RMAX
      double precision, allocatable :: RW(:), A(:,:)
      integer :: NCC, NRR, NWW
!
      do k=2,ncon-1
         write(44,*)' psiN=',psiv(k)/umax
         write(44,*)' R        Z           Bp'
         do i=1,npts
           write(44,*)rpts(k,i),zpts(k,i),bp(rpts(k,i),zpts(k,i))
         end do
      end do
      ndsk=57
      open(unit=ndsk,file=runname(1:lrunname)//'_1.gs2', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_1.gs2'
         stop
      endif
!
      nsurf=128
      allocate(ps(nsurf),temvar(nsurf))
      do i=1,nsurf
        ps(i)=(dble(i-1)/dble(nsurf-1))*umax*2.*pi
      end do
!      write(ndsk,fmt='("GS2 input file",T30,"Produced by SCENE at:",a40)') date
      write(ndsk,fmt='("GS2 input file",T30,"Produced by SCENE at:",a40)')
   write(ndsk,fmt='(A80/,"GS2D Equilibrium Boundary description:",3(/A80))') &
       repeat('-',80),repeat('-',80),repeat('-',80),repeat('-',80)
   zmag=0.
   write(ndsk,fmt='(T2,"r0",T15,"a",T27,"rmag",T39,"zmag (m)"/1p,4e16.8)') rcen,tokeps*rcen,r0,zmag
   write(ndsk,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0(T)",T39,"ip(A)"/1p,4e16.8)')              &
   2.*pi*psiv(ncon),2.*pi*umax,mu0*rodi/(2.*pi*rcen),cur
   write(ndsk,fmt='("nfs"/I6)') nsurf
   write(ndsk,fmt='("Psi on 1d grid (for FS quantities) (T)")')
   write(ndsk,fmt='(1p,8e16.8)') ps
! take out 2*pi factor from ps to be consistent with psiv from SCENE
      do i=1,nsurf
        ps(i)=(dble(i-1)/dble(nsurf-1))*umax
      end do
   do k=1,ncon
     aminor(k)=0.5*(rpts(k,npts/2+1)-rpts(k,1))
!     write(6,*)' k=',k,' psi=',psiv(k),' aminor=',aminor(k)
   end do
!   call spline1d(temvar,ps,nsurf,aminor,psiv,ncon,work)
!  Linearly interpolate
   do k=1,nsurf
     ik=1
     do 10 kk=1,ncon
       if (ps(k).gt.psiv(ik)) goto 10
       ik=ik+1
 10  continue
     if (ik.eq.1) ik=2
     if (ik.gt.ncon) ik=ncon
     rat=(ps(k)-psiv(ik))/(psiv(ik-1)-psiv(ik))
     if ((rat.gt.1.).or.(rat.lt.0.)) then
       write(6,*)' Error in linear interpolation for aminor in gs2'
       write(6,*)' rat=',rat,' ps=',ps(k),' psiv=',psiv(ik),psiv(ik-1)
       write(6,*)' ik=',ik
       stop
     end if
     temvar(k)=aminor(ik)+rat*(aminor(ik-1)-aminor(ik))
   end do
   write(ndsk,fmt='("amin (m)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   !write(6,*)' written aminor'
!   call spline1d(temvar,ps,nsurf,sfac,psiv,ncon,work)
   do k=1,nsurf
     ik=1
     do 20 kk=1,ncon
       if (ps(k).gt.psiv(ik)) goto 20
       ik=ik+1
 20  continue
     if (ik.eq.1) ik=2
     if (ik.gt.ncon) ik=ncon
     rat=(ps(k)-psiv(ik))/(psiv(ik-1)-psiv(ik))
     if ((rat.gt.1.).or.(rat.lt.0.)) then
       write(6,*)' Error in linear interpolation for aminor in gs2'
       write(6,*)' rat=',rat,' ps=',ps(k),' psiv=',psiv(ik),psiv(ik-1)
     end if
     temvar(k)=sfac(ik)+rat*(sfac(ik-1)-sfac(ik))
   end do
   write(ndsk,fmt='("q")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   !write(6,*)' written q'
   do k=1,ncon
      psi=psiv(k)
      fofpsi(k)=fprof(psi,1)
   end do
!!$   yp1=fprof(psiv(1),4)
!!$   yp2=fprof(psiv(ncon),4)
   yp1=-2.d30
   yp2=-2.d30
   call spline(psiv,fofpsi,ncon,yp1,yp2,y2)
   psi=umax
   fsqedg=0.5*fprof(psi,2)**2
   call zspline(psiv,fofpsi,y2,ncon,fsqedg,1,fpint)
!   do k=1,ncon
!     fpint(k)=sqrt(2.*fpint(k))
!     write(6,*)' k=',k,' fpint=',fpint(k)
!   end do
!   do k=2,ncon-1
!     psi=psiv(k)
!     write(56,*)' psi=',psi,' fprof=',2.*fprof(psi,1),' fpdiff=',(fpint(k+1)-fpint(k-1))/(psiv(k+1)-psiv(k))
!   end do
   !write(6,*)' doing ssplies'
   do k=1,nsurf
!     write(6,*)' k=',k
     psi=ps(k)
!     call zsplint(psiv,fofpsi,y2,fpint,ncon,psi,yy,yp,fsq)
!     temvar(k)=sqrt(2.*fsq)
!    call spline1d(temvar,ps,nsurf,fpint,psiv,ncon,work)
!     temvar(k)=fprof(psi,2)
!     write(6,*)' ps=',ps(k)/umax,' f-prime=',fprof(psi,1)/fprof(psi,2)
     ik=1
     do 40 kk=1,ncon
       if (ps(k).gt.psiv(ik)) goto 40
       ik=ik+1
 40  continue
     if (ik.eq.1) ik=2
     if (ik.gt.ncon) ik=ncon
     rat=(ps(k)-psiv(ik))/(psiv(ik-1)-psiv(ik))
     if ((rat.gt.1.).or.(rat.lt.0.)) then
       write(6,*)' Error in linear interpolation for aminor in gs2'
       write(6,*)' rat=',rat,' ps=',ps(k),' psiv=',psiv(ik),psiv(ik-1)
     end if
     temvar(k)=fpint(ik)+rat*(fpint(ik-1)-fpint(ik))
     temvar(k)=sqrt(2.*temvar(k))
     write(56,*)' psi=',psi,' fprof=',fprof(psi,2),' temvar=',temvar(k)
!     write(56,*)' psi=',psi,' fprof=',fprof(psi,1),' ffprime=',yy
   end do
   write(ndsk,fmt='("f =r B_phi (Tm)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   !write(6,*)' written f'
   do k=1,nsurf
     psi=ps(k)
     temvar(k)=press(psi,0)
   end do
   write(ndsk,fmt='("p (Pa)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   !write(6,*)' written p'
   do k=1,nsurf
     psi=ps(k)
     temvar(k)=press(psi,1)/(2.*pi)
   end do
   write(ndsk,fmt='("dp/dpsi (Pa/Wb)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   !write(6,*)' written p-prime'
   write(ndsk,fmt='("No of points on LCFS=",I6)') npts
   k=0
   do j=npts/2+1,1,-1
     k=k+1
     polvar(k)=rpts(1,j)
   end do
   do j=npts-1,npts/2+1,-1
     k=k+1
     polvar(k)=rpts(1,j)
   end do
   write(ndsk,fmt='("r(j) (m) on LCFS")')
   write(ndsk,fmt='(1p,8e16.8)') polvar
   !write(6,*)' written rpts'
   k=0
   do j=npts/2+1,1,-1
     k=k+1
     polvar(k)=zpts(1,j)
   end do
   do j=npts-1,npts/2+1,-1
     k=k+1
     polvar(k)=zpts(1,j)
   end do
   write(ndsk,fmt='("z(j) (m) on LCFS")')
   write(ndsk,fmt='(1p,8e16.8)') polvar
   !write(6,*)' written zpts'
   nrg=nsurf
   nzg=nrg+1
   allocate(rg(nrg),zg(nzg),mshvar(nrg,nzg))
   write(ndsk,fmt='("NR",T14,"NZ"/2I6)') nrg, nzg
   write(ndsk,fmt='("rgrid (m)")')
   deltar=(r(nr)-r(1))/(nrg-1)
   do i=1,nrg
     rg(i)=r(1)+(i-1)*deltar
   end do
   !write(6,*)' writing rg'
   write(ndsk,fmt='(1p,8e16.8)') rg
   !WRITE(6,*)' WRITTEN RG'
   deltaz=(z(nz)-z(1))/(nzg-1)
   do j=1,nzg
     zg(j)=z(1)+(j-1)*deltaz
   end do
   write(ndsk,fmt='("zgrid (m)")')
   write(ndsk,fmt='(1p,8e16.8)') zg
   !write(6,*)' written zg'
!--------------------Interpolate-----------------------------------
      npt=nr*nz
      allocate( xn(npt),yn(npt),fna(npt),grads(2,npt),triang(7*npt) )
      grads=0.
      triang=0.
      k=0
      do i=1,nr
        do j=1,nz
          k=k+1
          xn(k)=r(i)
          yn(k)=z(j)
          fna(k)=2.*pi*(umax-u(i,j))
        end do
      end do
      ifail=0
      ! Setup for interpolation
      NCC=MIN(17,npt-1)
      NWW=MIN(30,npt-1)
      NRR=ceiling( sqrt(dble(npt)/dble(3.0)) )
      allocate(LCELL(NRR,NRR))
      allocate(LNEXT(npt))
      allocate(RW(npt))
      allocate(A(9,npt))
!      call e01saf(npt,xn,yn,fna,triang,grads,ifail)
      call CSHEP2(npt,xn,yn,fna,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
      px=r0
      py=0.
      ifail=1
!      call e01sbf(npt,xn,yn,fna,triang,grads,px,py,pf,ifail)
      pf=CS2VAL(px,py,npt,xn,yn,fna,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
      do i=1,nrg
        do j=1,nzg
          px=rg(i)
          py=zg(j)
          ifail=1
!          call e01sbf(npt,xn,yn,fna,triang,grads,px,py,pf,ifail)
          pf=CS2VAL(px,py,npt,xn,yn,fna,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)
          mshvar(i,j)=pf
        end do
      end do
      deallocate( xn,yn,fna,grads,triang )
      !write(6,*)' done interpolate'
      call valbth2(rg,zg,mshvar,nrg,nzg)
!---------------------------------------------------------------------
!   do i=1,nr
!     do j=1,nz
!       mshvar(i,j)=2.*pi*(umax-u(i,j))
!     end do
!   end do
   write(ndsk,fmt='("Psi on grid (Wb) : NB Bpol=(1/2pi r) grad(phi) x grad(psi), in (r, phi, z) rh system")')
   write(ndsk,fmt='(1p,8e16.8)') mshvar
   deallocate(rg,zg,mshvar)
   close(ndsk)
!write(6,*) 'parameters at psinorm=0.4'
!iloc=minloc(abs(ps-ps(1)-0.4d0*(ps(nsurf)-ps(1))))
!write(6,*) 'q=',q(iloc(1)),q(iloc(1)+1),q(iloc(1)-1)
!write(6,*) 'shat=',( q(iloc(1)+1)-q(iloc(1)-1) )/( ps(iloc(1)+1)-ps(iloc(1)-1) )*0.4d0*(psip-psmin)/q(iloc(1))
!write(6,*) 'dbeta/drho=',( pps(iloc(1)+1)-pps(iloc(1)-1) )/( ps(iloc(1)+1)-ps(iloc(1)-1) )*(psip-psmin)*2.0d0*mu0/b0**2
!
!
      ndsk=58
      open(unit=ndsk,file=runname(1:lrunname)//'_2.gs2', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_2.gs2'
         stop
      endif

      
      bcentr = mu0*rodi/(2.*pi*rcen)
      nref = dense(psiv(ncon),0)
      tref = tempe(psiv(ncon),0)

      write(ndsk,*) ' nref = ', nref, 'm^-3, tref = ',tref/1000.,'keV'
      write(ndsk,*) ' mref = ', mp, 'kg, lref = ',amin,'m'
      call dpsidrho(dpdrs,rhos)
      !
      do k=2,ncon-1
         psi=psiv(k)
         shafr = shift(k,0)
         dshafr = shift(k,1)
         rmaj = rcen + shafr
         rmin = maxval(rpts(k,:)) - rmaj

         epsil = rmin/rmaj
         beta_gs2 = 4.03d-3*nref*1.0d-19*                   &
              (tref/1000.)*(2.*pi*rcen/(mu0*rodi))**2
         
         write(ndsk,*)' rhoc:', psi/umax
         ne19=dense(psi,0)*1.0d-19
         write(ndsk,*)' beta:',beta_gs2
         write(ndsk,*)' B0:',(mu0*rodi)/(2.*pi*rcen)
         write(ndsk,*)' eps:', epsil
         write(ndsk,*)' pk:', 2*amin/(rmaj*sfac(k))
         write(ndsk,*)' epsl:', 2*amin/rmaj
         write(ndsk,*)' shift:',dshafr*umax/amin
         write(ndsk,*)' s_hat_input:', qp(k)*psi/sfac(k)
         !print*, 'Shear info: ', psi/umax,sfac(k),qp(k),rhos(k),dpdrs(k)
         !         write(ndsk,*)' beta_prime_input:',press(psi,1)*umax*beta_gs2/press(psi,0)
         write(ndsk,*)' beta_prime_input:',press(psi,1)*8.0d-7*pi*umax/bcentr**2
         zeff=zm
         if (imp.eq.1) then
            if (ne19.gt.0.) then
               zeff=0.
               do i=1,nimp+1
                  zni19=densi(psi,i,0)*1.0d-19
                  zeff=zeff+(zni19*iz(i)**2)/ne19
               end do
            end if
         end if
         write(ndsk,*)' zeff:',zeff

         !Include all species
         !  No. of species
         nspec=nimp+2
         !write(ndsk,*)' nspec:',nspec
         !do i=1,nspec-1

         !    !Temp and dens
         !    ti = tempi(psi,i,0)
         !    ni = densi(psi,i,0)

         !    ! Temp and dens gradient of psi
         !    !
         !    ! Already normalised to Lref = aminor
         !    ti_p = tempi(psi,i,1)*umax
         !    ni_p = densi(psi,i,1)*umax

         !    write(ndsk,*)' Ion species no.',i
         !    write(ndsk,*)' z:',iz(i)
         !    write(ndsk,*)' mass:',zmas(i)
         !    write(ndsk,*)' dens:',ni/nref
         !    write(ndsk,*)' temp:',ti/tref

         !    ! Needs temp and dens derivatives of minor radius
         !    write(ndsk,*)' tprim:',-ti_p/ti
         !    write(ndsk,*)' fprim:',-ni_p/ni

         !    coolog=log(sqrt(ni*1.0d-6)/ti)
         !    coolog=24.-coolog

         !    !From Wesson 2.15 (added by bhavin 21/03/18)
         !    coll = 6.6d17*zmas(i)**0.5 * (ti/1000)**1.5/(ni*iz(i)**4*coolog)
         !    vss = 1./coll
         !    !vss = sqrt(2.0)*pi*ni * iz(i)**4 * eq**4 * coolog &
         !    !	/ ( sqrt(zmas(i)*mp) * (ti*eq/bk)**1.5  * (4*pi*eps0)**2 )

         !    write(ndsk,*)' vnewk:',vss*amin/ sqrt(2*tref*bk/mp)

         ! end do


         ! Just Electron and ion
         nspec=2
         write(ndsk,*)' nspec:',2
     
         !Temp and dens - Enforce quasi-neutrality
         ti = tempi(psi,1,0)
         ni = dense(psi,0)

         ! Temp and dens gradient of psi
         !
         ! Already normalised to Lref = aminor
         ti_p = tempi(psi,1,1)*umax
         ni_p = densi(psi,1,1)*umax

         write(ndsk,*)' Ion species no.',1
         write(ndsk,*)' z:',iz(1)
         write(ndsk,*)' mass:',zmas(1)
         write(ndsk,*)' dens:',ni/nref
         write(ndsk,*)' temp:',ti/tref

         ! Needs temp and dens derivatives of minor radius
         write(ndsk,*)' tprim:',-ti_p/ti
         write(ndsk,*)' fprim:',-ni_p/ni

         coolog=log(sqrt(ni*1.0d-6)/ti)
         coolog=24.-coolog

         !From Wesson 2.15 (added by bhavin 21/03/18)
         coll = 6.6d17*zmas(1)**0.5 * (ti/1000)**1.5/(ni*iz(1)**4*coolog)
         vss = 1./coll
         !vss = sqrt(2.0)*pi*ni * iz(i)**4 * eq**4 * coolog &
         !	/ ( sqrt(zmas(i)*mp) * (ti*eq/bk)**1.5  * (4*pi*eps0)**2 )

         write(ndsk,*)' vnewk:',vss*amin/ sqrt(2*tref*bk/mp)

 
         ! Electrons
         write(ndsk,*)' Electron species'
         write(ndsk,*)' Z:',-1.
         write(ndsk,*)' mass:',me/mp
         write(ndsk,*)' dens:',dense(psi,0)/nref
         write(ndsk,*)' temp:',tempe(psi,0)/tref

         ! Make differential of rho, not psi
         write(ndsk,*)' tprim:',-tempe(psi,1)*umax/tempe(psi,0)
         write(ndsk,*)' fprim:',-dense(psi,1)*umax/dense(psi,0)

         coolog=log(sqrt(dense(psi,0)*1.0d-6)/tempe(psi,0))
         coolog=24.-coolog
         !vss = sqrt(2.0)*pi*dense(psi,0) * eq**4 * coolog &
         !	/ (sqrt(me) * (tempe(psi,0)/1000.)**1.5 * (4*pi*eps0)**2 )

         coll = 3.*(2.*pi)**1.5* eps0**2 * me**0.5 * (tempe(psi,0)*bk)**1.5 &
              / ( dense(psi,0) * eq**4 * coolog)
         vss = 1./coll

         write(ndsk,*)' vnewk:',vss*amin/ sqrt(2*tref*bk/mp)
         write(ndsk,*)'end'
      end do

      close(ndsk)
      !
      !end subroutine write_gs2

!!$      do k=ncon,1,-1
!!$        write(6,*)' psin=',psiv(k)/umax,' fprof=',fprof(psiv(k),1)/fprof(psiv(k),2),' p-prime=',press(psiv(k),1)
!!$      end do

      return
    end subroutine gs2


!!$      do k=1,nspec-1
!!$        write(ndsk,*)' Ion species no.',k
!!$        write(ndsk,*)' Z:',iz(k)
!!$        write(ndsk,*)' mass:',zmas(k)/zmas(1)
!!$        write(ndsk,*)' dens:',densi(psi,k,0)*1.0d-19/ne19
!!$        write(ndsk,*)' temp:',tempi(psi,k,0)/tempi(psi,1,0)
!!$        write(ndsk,*)' tprim:',(umax/tempi(psi,k,0))*tempi(psi,k,1)
!!$        write(ndsk,*)' fprim:',(umax/densi(psi,k,0))*densi(psi,k,1)
!!$        write(ndsk,*)' vnewk: waiting for formula'
!!$      end do
!!$      write(ndsk,*)' Electron species'
!!$      write(ndsk,*)' Z:',-1.
!!$      write(ndsk,*)' mass:',me/(zmas(1)*mp)
!!$      write(ndsk,*)' dens:',dense(psi,0)*1.0d-19/ne19
!!$      write(ndsk,*)' temp:',tempe(psi,0)/tempi(psi,1,0)
!!$      write(ndsk,*)' tprim:',(umax/tempe(psi,0))*tempe(psi,1)
!!$      write(ndsk,*)' fprim:',(umax/dense(psi,0))*dense(psi,1)
!!$      write(ndsk,*)' vnewk: waiting for formula'
!!$      close(ndsk)
end module gs2_output
