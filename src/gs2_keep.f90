      subroutine gs2
!     --------------
!
!  Subroutine to calculate the various inputs required for the micro-
!  stability code GS2
!
!  HRW 22/05/03
!
!
      use param
      implicit none
!
      integer ndsk,k,nspec,l,kk,ik
      double precision dense,psi,zeff,zni19,ne19,fprof,press,tempi,tempe
      double precision densi
      double precision zmag,rat
      double precision temvar(nr),work(2*ncon),polvar(npts),mshvar(nr,nz)
      double precision aminor(ncon)
      character(len=40) date
      integer, dimension(1):: iloc
!      character(len=40), intent(in) :: filename
      double precision, dimension(:), allocatable :: ps
      integer :: nsurf,i,j
!      type(vec2dl) :: sep
!

      ndsk=57
      open(unit=ndsk,file=runname(1:lrunname)//'_1.gs2', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'_1.gs2'
         stop
      endif
!
      nsurf=nr
      allocate(ps(nsurf))
      do i=1,nsurf
        ps(i)=(dble(i-1)/dble(nsurf-1))*umax*2.*pi
      end do
      write(ndsk,fmt='("GS2 input file",T30,"Produced by SCENE at:",a40)') date
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
     write(6,*)' k=',k,' psi=',psiv(k),' aminor=',aminor(k)
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
     rat=(ps(k)-psiv(ik))/(psiv(ik-1)-psiv(ik))
     if ((rat.gt.1.).or.(rat.lt.0.)) then
       write(6,*)' Error in linear interpolation for aminor in gs2'
       write(6,*)' rat=',rat,' ps=',ps(k),' psiv=',psiv(ik),psiv(ik-1)
     end if
     temvar(k)=aminor(ik)+rat*(aminor(ik-1)-aminor(ik))
   end do
   write(ndsk,fmt='("amin (m)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
!   call spline1d(temvar,ps,nsurf,sfac,psiv,ncon,work)
   do k=1,nsurf
     ik=1
     do 20 kk=1,ncon
       if (ps(k).gt.psiv(ik)) goto 20
       ik=ik+1
 20  continue
     if (ik.eq.1) ik=2
     rat=(ps(k)-psiv(ik))/(psiv(ik-1)-psiv(ik))
     if ((rat.gt.1.).or.(rat.lt.0.)) then
       write(6,*)' Error in linear interpolation for aminor in gs2'
       write(6,*)' rat=',rat,' ps=',ps(k),' psiv=',psiv(ik),psiv(ik-1)
     end if
     temvar(k)=sfac(ik)+rat*(sfac(ik-1)-sfac(ik))
   end do
   write(ndsk,fmt='("q")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   do k=1,nsurf
     psi=ps(k)/(2.*pi)
     temvar(k)=fprof(psi,2)
   end do
   write(ndsk,fmt='("f =r B_phi (Tm)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   do k=1,nsurf
     psi=ps(k)
     temvar(k)=press(psi,0)
   end do
   write(ndsk,fmt='("p (Pa)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   do k=1,nsurf
     psi=ps(k)
     temvar(k)=press(psi,1)/(2.*pi)
   end do
   write(ndsk,fmt='("dp/dpsi (Pa/Wb)")')
   write(ndsk,fmt='(1p,8e16.8)') temvar
   write(ndsk,fmt='("No of points on LCFS=",I6)') npts
   do j=1,npts
     polvar(j)=rpts(1,j)
   end do
   write(ndsk,fmt='("r(j) (m) on LCFS")')
   write(ndsk,fmt='(1p,8e16.8)') polvar
   do j=1,npts
     polvar(j)=zpts(1,j)
   end do
   write(ndsk,fmt='("z(j) (m) on LCFS")')
   write(ndsk,fmt='(1p,8e16.8)') polvar
   write(ndsk,fmt='("NR",T14,"NZ"/2I6)') NR, NZ
   write(ndsk,fmt='("rgrid (m)")')
   write(ndsk,fmt='(1p,8e16.8)') r
   write(ndsk,fmt='("zgrid (m)")')
   write(ndsk,fmt='(1p,8e16.8)') z
   do i=1,nr
     do j=1,nz
       mshvar(i,j)=2.*pi*(umax-u(i,j))
     end do
   end do
   write(ndsk,fmt='("Psi on grid (Wb) : NB Bpol=(1/2pi r) grad(phi) x grad(psi), in (r, phi, z) rh system")')
   write(ndsk,fmt='(1p,8e16.8)') mshvar
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
!
      k=ncon/2
      psi=psiv(k)
      write(ndsk,*)' normalised flux of chosen surface:',psi
      ne19=dense(psi,0)*1.0d-19
      write(ndsk,*)' beta:',4.03d-3*ne19*                   &
                   (tempi(psi,1,0)/1000.)*2.*pi*rcen/(mu0*rodi*1.0d6)
      zeff=zm
      if (imp.eq.1) then
        if (ne19.gt.0.) then
          zeff=0.
          do l=1,nimp+1
            zni19=densi(psi,l,0)*1.0d-19
            zeff=zeff+(zni19*iz(l)**2)/ne19
          end do
        end if
      end if
      write(ndsk,*)' zeff:',zeff
!  No. of species
      nspec=nimp+2
      write(ndsk,*)' nspec:',nspec
      do l=1,nspec-1
        write(ndsk,*)' Ion species no.',l
        write(ndsk,*)' Z:',iz(l)
        write(ndsk,*)' mass:',zmas(l)/zmas(1)
        write(ndsk,*)' dens:',densi(psi,l,0)*1.0d-19/ne19
        write(ndsk,*)' temp:',tempi(psi,l,0)/tempi(psi,1,0)
        write(ndsk,*)' tprim:',(umax/tempi(psi,l,0))*tempi(psi,l,1)
        write(ndsk,*)' fprim:',(umax/densi(psi,l,0))*densi(psi,l,1)
        write(ndsk,*)' vnewk: waiting for formula'
      end do
      write(ndsk,*)' Electron species'
      write(ndsk,*)' Z:',-1.
      write(ndsk,*)' mass:',me/(zmas(1)*mp)
      write(ndsk,*)' dens:',dense(psi,0)*1.0d-19/ne19
      write(ndsk,*)' temp:',tempe(psi,0)/tempi(psi,1,0)
      write(ndsk,*)' tprim:',(umax/tempe(psi,0))*tempe(psi,1)
      write(ndsk,*)' fprim:',(umax/dense(psi,0))*dense(psi,1)
      write(ndsk,*)' vnewk: waiting for formula'
      close(ndsk)
!
!end subroutine write_gs2



      return
 end subroutine gs2
