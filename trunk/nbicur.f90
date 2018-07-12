subroutine nbicur()

  use param
  implicit none

  double precision :: psi


  integer :: con, l, im,i,j, count

  integer :: nbeams, n, nion, maxiter, inbfus, iflag

  real :: amb, zbeam, b0, volp
  double precision :: rr, nesum, tesum, neav,teav

  double precision, dimension(ncon) :: phi
  real :: zni, ne

  ! input variables
  real, dimension(:), allocatable :: ebeam, pbeam, rtang, bwidth, bheigh, bgaussR, bgaussZ, bzpos

  integer, dimension(:), allocatable :: nbshape, nbptype

  real, dimension(:), allocatable :: aion, zion, ne20, tekev, tikev, zeff
  real, dimension(:,:), allocatable :: ni20, pwrfrac

  real, dimension(:), allocatable :: rnormnb, kappa, dkappa, shafr, dshafr

  real, dimension(:), allocatable :: dpdrs, vprime_norm, vprime, dvol, darea

  double precision :: dense, tempe, tempi, densi, shift, elong
  !Output variables

  real, dimension(:,:,:), allocatable :: hofr, jnbie

  real, dimension(:,:), allocatable :: shinethru, jnb

  real, dimension(:), allocatable ::  jnbTot, pnbe, pnbi, beamDens, beamVel, beamPress, beamFus, &
       pNBLoss, pNBAbsorb, pbfuse, pbfusi, snBeamDD, snBeamDT, nbcur, etanb, &
       gammanb, jnbfast, l31

  real :: pNBAbsorbTot, pNBLossTot, nbcurTot, etanbTot, beamBeta, beamFusTot, beamFusChTot, &
       snDTTotal, snDDTotal

  integer :: iflagnb

  integer :: maxbeams, mxrho

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
        !print*, dense(psi,0) , tempe(psi,0)
        tesum = tesum + tempe(psi,0)
     end do
  end do

  !print*, 'Tesum ', tesum, 'nesum', nesum, count
  neav = nesum/count
  teav = tesum/count

  volp = sngl(vol*2*pi)
  print*, 'Plasma volume', volp

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


  pwrfrac(:,2)=pwrfrac(:,1)
  close(53)

  paux=0.
  do l=1,nbeams
        paux = paux + pbeam(l)
  end do

  !print*, nbshape

  !Background plasma input quantities



  allocate(aion(nion), zion(nion), ne20(ncon), tekev(ncon), tikev(ncon), zeff(ncon), &
       dpdrs(ncon), rnormnb(ncon), dvol(ncon), darea(ncon), vprime(ncon), kappa(ncon), &
       dkappa(ncon), shafr(ncon), dshafr(ncon), vprime_norm(ncon), l31(ncon) )


  !Ions parameters
  if ( nimp+2 .gt. 6) then
     nion = 6
  else
     nion = nimp+1
  end if

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

  call trapfact(l31)

  do con=1,ncon

     psi = psiv(con)

     !Shafranov shift  and elongation + derivatives (all need to be reversed)
     !ncon-con+1 because array needs to start from axis
     shafr(ncon-con+1) = sngl(shift(con,0)/amin)

     dshafr(ncon-con+1) = sngl(shift(con,1)/amin * dpdrs(ncon-con+1))

     kappa(ncon-con+1) = sngl(elong(con,0))

     dkappa(ncon-con+1) = sngl(elong(con,1) * dpdrs(ncon-con+1))


     !Already in the correct order
     vprime_norm(ncon-con+1) = sngl(vprime(ncon-con+1) * dpdrs(ncon-con+1)/amin)


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
        aion(im) = zmas(im)
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
              zni=sngl(densi(psi,l,0))
              zeff(ncon-con+1)=sngl(zeff(ncon-con+1)+(zni*iz(l)**2)/ne)
           end do

        end if
     end if


  end do


  do i=1,ncon



  !   rr = maxval(rpts(ncon-i+1,:)) - r0
  !
  !   print*, rnormnb(i), vprime(i),  shafr(i), dshafr(i)
  end do


  !Need to use limit of V' at very small rho

  vprime_norm(2)=4.*pi*pi*rnormnb(2)*amin*kappa(2)*(rcen+amin*shafr(2))

  !Set density to small value at edge otherwise NaNs
  ne20(ncon) = 0.01
  ni20(ncon,:) = 0.01


  !Allocate arrays to outputs variables (need to do it for max no. of
  !beams/rhos as is allocated in nbeams
  allocate(hofr(mxrho, 3,maxbeams))

  allocate(shinethru(3,maxbeams))

  allocate(pNBLoss(maxbeams), pNBAbsorb(maxbeams), etanb(maxbeams), gammanb(maxbeams), nbcur(maxbeams))

  allocate(jnbTot(mxrho), pnbe(mxrho),pnbi(mxrho),beamDens(mxrho),beamVel(mxrho),beamPress(mxrho), &
       beamFus(mxrho), pbfuse(mxrho), pbfusi(mxrho), snBeamDD(mxrho), snBeamDT(mxrho), jnbfast(mxrho))


  ! Call beam calc
  call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus, &
       rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ, &
       bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev, &
       tikev, zeff, sngl(rcen), sngl(amin), b0, sngl(vol), ncon, rnormnb, vprime_norm,  &
       dvol, darea, l31,  &
       kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,  &
       pnbi, beamDens, beamVel, beamPress, beamFus, jnbfast, pbfuse, pbfusi, snBeamDD, &
       snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot, &
       etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot, &
       beamFusChTot, snDTTotal, snDDTotal, iflagnb)


  ni20=0.0

  !Set total current to J_nb
  !Note need to reverse order of array again and use <B>/<B^2>
  ! from appropriate flxsur
  do i=1,ncon
     !print*, i, beamVel(i), beamDens(i)
     J_nb(ncon-i+1) = dble(jnbTot(i))*bav(ncon-i+1)/bsqav(ncon-i+1)

  end do

  print*, 'flag', iflag
  print*, 'Total NB Current: ', nbcurTot

  !Shinethrough calc for all beams + components
  bmshine=0.
  do i=1,nbeams
     do j=1,3
        bmshine= bmshine +pwrfrac(j,i)*shinethru(j,i)
     end do
  end do

  bmshine = bmshine/nbeams
  bmfus = dble(beamFusTot)


   !Flx average quantities
   open(unit=50, file=runname(1:lrunname)//'_flxav.dat', &
        status='unknown', iostat=ios)
   if (ios.ne.0) then
      write(6,*) 'problem opening ',runname(1:lrunname)//'_flxav.dat'
      stop
   end if


   do con=1,ncon

      write(50,30) phi(ncon-con+1), jnbTot(con), pnbe(con), pnbi(con), beamDens(con), beamVel(con), dvol(con), jnbfast(con)

   end do
30 format(8e14.6)

   close(50)

  !print*, 'Charge beam fusion ', beamFusChTot
  !print*, 'Total beam target fusion ', beamFusTot

  !print*, 'Electron Average Temp ', teav, q' and density ', neav
end subroutine nbicur



subroutine trapfact(l31)
!     **********************************************
!
!  Trappec electron correct from Lin-Liu  & Hinton
!  Phys. Plasma 4, 4179 (1997)
!
      use param
      implicit none
      double precision dense,densi
      double precision ne,tau,zni,zeff,zb
      double precision pe,fsi,dox
      double precision x,rj0,psi,bsq,pt,pd,bstrap
      double precision rla, dst, btot,rr,zz,bphi,bp
      double precision bigint,fc, fprof
      integer i,l,k
!
      real, dimension(ncon) :: l31

      do k=1,ncon
         psi = psiv(k)
         ne=dense(psi,0)
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
               nw=10
               write(nw,*)'error*** problem in trapfact, ne=0'
               write(nw,*)'cannot evaluate zeff'
               stop
            end if
         end if
         zb=zeff

         fc=1.-ftrap(k)
         x=ftrap(k)/fc
         dox=1.414*zb+zb*zb+x*(0.754+2.657*zb+2.*zb*zb)+x*x*    &
              (0.348+1.243*zb+zb*zb)

         !Flip the order for NBEAMS
         l31(ncon-k+1)=real(x*(0.754+2.21*zb+zb*zb+x*(0.348+1.243*zb+zb*zb))  &
              /dox)
      end do



    end subroutine trapfact
