module nbi_mod
  implicit none
contains
subroutine nbicur()

  use param
  use profiles_mod, only : dense, tempe, tempi, densi, shift, elong, dpsidrho, dVdrho, rhotor
  implicit none

  double precision :: psi


  integer :: con, l, im,i,j, count, rr_index

  integer :: nbeams, nion, maxiter, inbfus, iflag

  real :: amb, zbeam, b0, volp
  double precision :: rr, nesum, tesum, neav,teav
  double precision :: bphi, bth, btot, fprof
  double precision, dimension(ncon) :: phi
  real :: zni, ne
  
  ! Input variables
  real, dimension(:), allocatable :: ebeam, pbeam, rtang, bwidth, bheigh, bgaussR, bgaussZ, bzpos

  integer, dimension(:), allocatable :: nbshape, nbptype

  real, dimension(:), allocatable :: aion, zion, ne20, tekev, tikev, zeff
  real, dimension(:,:), allocatable :: ni20, pwrfrac

  real, dimension(:), allocatable :: kappa, dkappa, shafr, dshafr

  real, dimension(:), allocatable :: vprime_norm, vprime, dvol, darea
  real, dimension(ncon) :: rhos

  double precision, dimension(:), allocatable :: dpdrs, rnormnb
  !Output variables

  real, dimension(:,:,:), allocatable :: hofr, srcfast, pitchangl

  real, dimension(:,:), allocatable :: shinethru, pnbbm, pnb

  real, dimension(:), allocatable ::  jnbTot, pnbe, pnbi, beamDens, beamVel, beamPress, beamFus
  real, dimension(:), allocatable :: pNBLoss, pNBAbsorb, pbfuse, pbfusi, snBeamDD, snBeamDT, nbcur
  real, dimension(:), allocatable :: etanb, gammanb, jnbfast, l31

  real :: pNBAbsorbTot, pNBLossTot, nbcurTot, etanbTot, beamBeta
  real :: beamFusTot, beamFusChTot, snDTTotal, snDDTotal

  integer :: iflagnb

  integer :: maxbeams, mxrho

  double precision, dimension(:), allocatable :: source
  maxbeams = 4
  mxrho = 51

    !J_nb = 0.
  !Nbi input param: no. of beams, mass of beams, charge, beam-fusion flag
  nbeams = 2
  amb = 2
  zbeam = 1
  inbfus = 1

  vol = 0.
  nesum =0.
  tesum = 0.
  count=0
  do i=1,nr

     rr = r(i)
     do j=1,nz
        if (ixout(i,j) .le. 0) cycle
        psi = umax - u(i,j)
        vol = vol+rr*dr*dz
        count = count +1
        nesum = nesum + dense(psi,0)
        tesum = tesum + tempe(psi,0)
     end do
  end do

  neav = nesum/count
  teav = tesum/count

  volp = sngl(vol*2*pi)

  !Iterations to find rho
  maxiter = 2

  open(unit=53,file='nbi.dat', iostat=ios)
  if (ios.ne.0) then
     write(nw,*) 'Error reading in nbi.dat file'
     stop
  end if
  
  read(53,*)
  read(53,*)
  
  read(53,*) nbeams
        

  allocate(ebeam(nbeams), pbeam(nbeams), rtang(nbeams), bwidth(nbeams), bheigh(nbeams), &
       nbptype(nbeams), bgaussR(nbeams), bgaussZ(nbeams), bzpos(nbeams), nbshape(nbeams) )

  allocate(pwrfrac(3, nbeams))


  read(53,*)
  read(53,*) ebeam

  read(53,*)
  read(53,*) pbeam


  read(53,*)
  read(53,*) rtang


  read(53,*)
  read(53,*) bzpos


  read(53,*)
  read(53,*) nbshape


  read(53,*)
  read(53,*) nbptype


  read(53,*)
  read(53,*) bwidth


  read(53,*)
  read(53,*) bheigh


  read(53,*)
  read(53,*) bgaussR


  read(53,*)
  read(53,*) bgaussZ


  read(53,*)
  read(53,*) pwrfrac(:,1)

  if (nbeams .eq. 2)  pwrfrac(:,2)=pwrfrac(:,1)
  close(53)

  paux=0.
  do l=1,nbeams
        paux = paux + pbeam(l)
  end do

  !Background plasma input quantities

  !Ions parameters
  if ( nimp+2 .gt. 6) then
     nion = 6
  else
     nion = nimp+1
  end if

  allocate(aion(nion), zion(nion), ne20(ncon), tekev(ncon), tikev(ncon), zeff(ncon), &
       dpdrs(ncon), rnormnb(ncon), dvol(ncon), darea(ncon), vprime(ncon), kappa(ncon), &
       dkappa(ncon), shafr(ncon), dshafr(ncon), vprime_norm(ncon), l31(ncon) )


  !aion = zmas
  !zion = iz
  !Toroidal magnetic field
  b0 = sngl(mu0*rodi/(2.*pi*r0))
  amin= (maxval(rpts(1,:))-minval(rpts(1,:)) )/2

  dpdrs=0.
  rnormnb=0.

  allocate( ni20(mxrho, nion))

  !Same order to NBEAM
  call dpsidrho(dpdrs, rnormnb)

  !Same order as NBEAM
  call dVdrho(dvol, darea, vprime)

  call rhotor(phi)
  !rho_tor = sqrt(Phi_tor)
  phi = sqrt(phi/phi(1))

  !call trapfact(l31)
  l31 = 0.

  do con=1,ncon
     
     psi = psiv(con)

     !Shafranov shift  and elongation + derivatives (all need to be reversed)
     !ncon-con+1 because array needs to start from axis
     shafr(ncon-con+1) = sngl(shift(con,0)/amin)

     dshafr(ncon-con+1) = sngl(shift(con,1)/amin * dpdrs(con))

     kappa(ncon-con+1) = sngl(elong(con,0))

     dkappa(ncon-con+1) = sngl(elong(con,1) * dpdrs(con))


     !Already in the correct order
     vprime_norm(ncon-con+1) = sngl(vprime(ncon-con+1) * dpdrs(con)/amin)


     !Electron and ion densities and temp (need to be reversed)
     ne = sngl(dense(psi,0))
     ne20(ncon-con+1) = ne*1.0e-20



     !Split ions into half D and half T
     ni20(ncon-con+1,1) = sngl(densi(psi,1,0)*1.0e-20)/2.
     aion(1) = 2.0
     zion(1) = 1.0

     ni20(ncon-con+1,nion) = sngl(densi(psi,1,0)*1.0e-20)/2.
     aion(nion) = 3.0
     zion(nion) = 1.0


     do im=2,nion-1

        !Half as D
        ni20(ncon-con+1,im) = sngl(densi(psi,im,0)*1.0e-20)
        aion(im) = real(zmas(im))
        zion(im) = iz(im)

     end do

     tekev(ncon-con+1) = sngl(tempe(psi,0)*1.0e-3)

     tikev(ncon-con+1) = sngl(tempi(psi,1,0)*1.0e-3)


     !Effective Z
     zeff(ncon-con+1)=sngl(zm)
     if (imp.eq.1) then
        if (ne.gt.0.) then
           zeff(ncon-con+1)=0.
           do l=1,nimp+1
              zni=real(densi(psi,l,0))
              zeff(ncon-con+1)=real(zeff(ncon-con+1)+(zni*iz(l)**2)/ne)
           end do

        end if
     end if


  end do

  !Need to use limit of V' at very small rho

  vprime_norm(2)=real(4.*pi*pi*rnormnb(ncon-1)*amin*kappa(2)*(rcen+amin*shafr(2)))

  !Set density to small value at edge otherwise NaNs
  ne20(ncon) = 0.01
  ni20(ncon,:) = 0.01

  rhos = sngl(rnormnb(ncon:1:-1))


  open(unit=55,file=runname(1:lrunname)//'.nbeams', &
       status='unknown', iostat=ios)
  if (ios.ne.0) then
     write(nw,*) 'Error opening nbeams file'
     stop
  end if

  write(55, *) nion, sngl(rcen), sngl(amin), b0, sngl(vol), ncon
  write(55, *) aion, zion

  do con=1,ncon
     write(55,*) rhos(con), vprime_norm(con), dvol(con), darea(con), kappa(con), dkappa(con), &
          shafr(con), dshafr(con), zeff(con), ne20(con), ni20(con, :), tekev(con), tikev(con), psiv(ncon-con+1)
  end do
  close(55)
  

  !Allocate arrays to outputs variables (need to do it for max no. of
  !beams/rhos as is allocated in nbeams
  allocate(hofr(mxrho, 3,maxbeams), srcfast(mxrho,3,maxbeams), pitchangl(mxrho, 3, maxbeams))

  allocate(shinethru(3,maxbeams))

  allocate(pNBLoss(maxbeams), pNBAbsorb(maxbeams), etanb(maxbeams), gammanb(maxbeams), nbcur(maxbeams))

  allocate(jnbTot(mxrho), pnbe(mxrho),pnbi(mxrho),beamDens(mxrho),beamVel(mxrho), &
       beamPress(mxrho), beamFus(mxrho), pbfuse(mxrho), pbfusi(mxrho),  &
       snBeamDD(mxrho), snBeamDT(mxrho), jnbfast(mxrho), pnbbm(mxrho,maxbeams), &
       pnb(ncon,nbeams)) 


  ! Call beam calc
  call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus, &
       rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ, &
       bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev, &
       tikev, zeff, sngl(rcen), sngl(amin), b0, sngl(vol), ncon, rhos, vprime_norm,  &
       dvol, darea, l31, srcfast,  &
       kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,  &
       pnbi, pnbbm,beamDens, beamVel, beamPress, beamFus, jnbfast, pbfuse, &
       pbfusi, snBeamDD, &
       snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot, &
       etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot, &
       beamFusChTot, snDTTotal, snDDTotal, pitchangl, iflagnb)

  print*, 'Completed NBI calc'
  ni20=0.0

  !Set total current to J_nb
  !Note need to reverse order of array again and use <B>/<B^2>
  ! from appropriate flxsur
  do i=1,ncon

     J_nb(ncon-i+1) = dble(jnbTot(i))*bav(ncon-i+1)/bsqav(ncon-i+1)

     !Account for poloidal field by multiplying Bt/B
     !Use poloidal field at OMP for each flux surface
     !Not great as some beam is off axis

     ! psi = psiv(ncon-i+1)
     ! rr = maxval(rpts(ncon-i+1,:))

     ! ! Calculate Bphi and Btheta
     ! rr_index = maxloc(rpts(ncon-i+1,:),1) 
     ! bphi = fprof(psi,2)/rr
     ! bth = bppts(ncon-i+1,rr_index)
     ! btot = sqrt(bth**2 + bphi**2)

     ! if (i .eq. 1) then
     !    bphi = fprof(0.0,2)/r0
     !    btot = bphi
     ! end if
     ! J_nb(ncon-i+1) = J_nb(ncon-i+1)*bphi/btot

  end do

  print*, 'flag', iflag
  print*, 'Total NB Current: ', nbcurTot
  print*, 'Total Power lost (MW): ',pNBLoss
  print*, 'Total Power absorbed (MW): ', pNBAbsorbTot
  !Shinethrough calc for all beams + components
  ! bmshine=0.
  ! do i=1,nbeams
  !    do j=1,3
  !       bmshine= bmshine +pwrfrac(j,i)*shinethru(j,i)
  !    end do
  ! end do

  !Sim power to e and i and reverse
  
  do i=1, nbeams
     ! Set Power to vol elements - reverse order (so rho 0->1)
     pnb(:,i) = dvol(ncon:1:-1)
     do j=1,ncon
        ! Times by beam power den (reverse order)
        pnb(j,i) = pnb(j,i) * pnbbm(ncon-j+1,i)

     end do
     
  end do
  
        
  allocate(source(ncon))
  call momsource(nbeams,dble(pnb),dble(ebeam),dble(rtang), source)
  nbmom = source
  
  
  bmshine = dble(sum(pNBLoss(1:nbeams)))
  bmfus = dble(beamFusTot)


   !Flx average quantities
   open(unit=48, file=runname(1:lrunname)//'_srcfast.dat', &
        status='unknown', iostat=ios)
   if (ios.ne.0) then
      write(6,*) 'problem opening ',runname(1:lrunname)//'_srcfast.dat'
      stop
   end if

   write(48,*) 'Fast source'

   do i = 1,ncon
      !Some along energy components
      write(48,28) srcfast(i,:,:nbeams)
   end do


28 format(8e13.6)
   close(48)
   
   !Flx average quantities
   open(unit=50, file=runname(1:lrunname)//'_flxav.dat', &
        status='unknown', iostat=ios)
   if (ios.ne.0) then
      write(6,*) 'problem opening ',runname(1:lrunname)//'_flxav.dat'
      stop
   end if

   write(50,*) 'Psi, NB cur (A/m^2), Pow to e (MW/m^3), Power to i (MW/m^3),  &
        &FI dens (/m^3), FI Vel (/m^2 s), flux sur vol (m^3), FI current (A/m^2)'
   do con=1,ncon
      write(50,30) psiv(ncon-con+1), jnbTot(con), pnbe(con), pnbi(con), beamDens(con), beamVel(con), dvol(con), &
           jnbfast(con), darea(con), pitchangl(con, 1, :2)

   end do
30 format(11e12.5)

   close(50)

end subroutine nbicur


subroutine momsource(nbeams,power,energy, rtan, source)

!Momentum from Beam energy (keV) and power (MW)
!Returns NBI momenum
  
  use param
  use profiles_mod, only : densi
  implicit none
  integer, intent(in) :: nbeams
  double precision, dimension(ncon,nbeams), intent(in) :: power  
  double precision, dimension(nbeams), intent(in) ::  rtan, energy
  double precision, dimension(ncon), intent(out) :: source
  double precision, dimension(ncon) :: bm_source,mass, vtor, angf
  double precision :: rmax

  integer :: i, con

  source = 0 
  
  do i=1,nbeams

     !Source from each beamline
     bm_source = power(:,i)*1e6*rtan(i)* &
          (energy(i)*1e3*eq/(2.*mp))**(-0.5)
     do con=1,ncon
        !Momentum source term for each flux surface
        rmax = maxval(rpts(con,:))
        if (con.eq.ncon) rmax = r0
        source(con) = source(con)+ bm_source(con)/rmax
        mass(con) = densi(psiv(con),1,0) *mp*zmai
     end do
  end do
  !Vtor = Mom*time/mass
  vtor = source  * taue /mass
  
  !Angular freq at R_outer
  angf = vtor/(maxval(rpts,2))

end subroutine momsource




subroutine trapfact(l31)
!     **********************************************
!
!  Trappec electron correct from Lin-Liu  & Hinton
!  Phys. Plasma 4, 4179 (1997)
!
      use param
      use profiles_mod, only : dense, densi, fprof
      implicit none      
      real, dimension(ncon), intent(out) :: l31
      double precision :: ne,zni,zeff,zb
      double precision :: fsi,dox
      double precision :: x,psi
      double precision :: fc
      integer :: l,k
!


      do k=1,ncon
         psi = psiv(k)
         ne=dense(psi,0)
         print*, ne
         fsi=fprof(psi,2)
         zeff=zm
         if (imp.eq.1) then
            if (ne.gt.0.) then
               zeff=0.
               do l=1,nimp+1
                  zni=densi(psi,l,0)
                  zeff=zeff+(zni*iz(l)**2)/ne
               end do
            else
               write(*,*) 'error*** problem in trapfact, ne=0'
               write(*,*) 'cannot evaluate zeff'
               stop
            end if
         end if
         zb=zeff
         print*, ftrap(k)
         fc=1.-ftrap(k)
         x=ftrap(k)/fc
         dox=1.414*zb+zb*zb+x*(0.754+2.657*zb+2.*zb*zb)+x*x*    &
              (0.348+1.243*zb+zb*zb)
         print*, dox
         !Flip the order for NBEAMS
         l31(ncon-k+1)=real(x*(0.754+2.21*zb+zb*zb+x*(0.348+1.243*zb+zb*zb))  &
              /dox)
      end do



    end subroutine trapfact

end module nbi_mod
