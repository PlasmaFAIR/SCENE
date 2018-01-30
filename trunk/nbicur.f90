subroutine nbicur()

  use param
  implicit none

  double precision :: psi


  integer :: con, l, im,i,j, count

  integer :: nbeams, n, nion, maxiter, inbfus, iflag

  real :: amb, zbeam, b0, volp
  double precision :: rr, nesum, tesum, neav,teav
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

  real, dimension(:), allocatable ::  jnbTot, pnbe, pnbi, beamDens, beamPress, beamFus, &
       pNBLoss, pNBAbsorb, pbfuse, pbfusi, snBeamDD, snBeamDT, nbcur, etanb, &
       gammanb

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
       nbptype(nbeams), bgaussR(nbeams), bgaussZ(nbeams), bzpos(nbeams))
  
  allocate(pwrfrac(3, nbeams))

	print*, nbeams
	
	read(53,*)
	read(53,*) ebeam

	read(53,*)
	read(53,*) pbeam


	read(53,*)
	read(53,*) rtang


	read(53,*)
	read(53,*) bzpos


	read(53,*)
	read(53,*) pwrfrac(:,1)


	
	pwrfrac(:,2)=pwrfrac(:,1)
	close(53)
	
	print*, pwrfrac
  
  nbshape = (/0,1/)
  bwidth = (/0.2,0.3/)
  bheigh = (/0.2,0.3/)
  nbptype = (/1,1/)

  bgaussR = (/0.1,0.1/)
  bgaussZ = (/0.1,0.1/)




  paux=0.
  do l=1,nbeams
     if (nbshape(l) .eq. 0) then
        paux = paux + pbeam(l)
     else
        paux = paux + 2*pbeam(l)
     end if
  end do
  
 

  !Background plasma input quantities



  allocate(aion(nion), zion(nion), ne20(ncon), tekev(ncon), tikev(ncon), zeff(ncon), &
       dpdrs(ncon), rnormnb(ncon), dvol(ncon), darea(ncon), vprime(ncon), kappa(ncon), &
       dkappa(ncon), shafr(ncon), dshafr(ncon), vprime_norm(ncon) )
  

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

  dpdrs=0.
  rnormnb=0.
  
  allocate( ni20(mxrho, nion))

	!Opposite order to NBEAM
  call dpsidrho(dpdrs, rnormnb)

  !Same order as NBEAM
  call dVdrho(dvol, darea, vprime)

  do con=1,ncon

     psi = psiv(con)
     

     !Shafranov shift  and elongation + derivatives
     !ncon-con+1 because array needs to start from axis
     shafr(ncon-con+1) = sngl(shift(con,0)/amin)

     dshafr(ncon-con+1) = sngl(shift(con,1)/amin) * dpdrs(con)

     kappa(ncon-con+1) = sngl(elong(con,0))

     dkappa(ncon-con+1) = sngl(elong(con,1)) * dpdrs(con)

     vprime_norm(ncon-con+1) = sngl(vprime(ncon-con+1) * dpdrs(con))

     !Electron and ion densities and temp
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
              zeff(ncon -con+1)=sngl(zeff(ncon-con+1))+(zni*iz(l)**2)/ne
           end do
          
        end if
     end if

     !print*,  ne20(ncon-con+1), ni20(ncon-con+1,1)
  end do

  !do i=1,ncon
	!		print*, vprime_norm(i), dvol(i), darea(i), psiv(i)/umax

  !
  !   rr = maxval(rpts(ncon-i+1,:)) - r0
  !   
  !   print*, rnormnb(i), vprime(i),  shafr(i), dshafr(i)
  !end do
  

  !Print*, 'calling nbeams'

!!!!! CHECK IF THIS IS VALID !!!!!!
  !Set density to small value at edge other NaNs
  ne20(ncon) = 0.01
  ni20(ncon,:) = 0.01
	
  
  !Allocate arrays to outputs variables (need to do it for max no. of
  !beams/rhos as is allocated in nbeams
  allocate(hofr(mxrho, 3,maxbeams))

  allocate(shinethru(3,maxbeams))

  allocate(pNBLoss(maxbeams), pNBAbsorb(maxbeams), etanb(maxbeams), gammanb(maxbeams), nbcur(maxbeams))

  allocate(jnbTot(mxrho), pnbe(mxrho),pnbi(mxrho),beamDens(mxrho),beamPress(mxrho), beamFus(mxrho), &
       pbfuse(mxrho), pbfusi(mxrho), snBeamDD(mxrho), snBeamDT(mxrho))

  ! Call beam calc
  call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus, &
       rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ, &
       bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev, &
       tikev, zeff, sngl(rcen), sngl(amin), b0, volp, ncon, rnormnb, vprime_norm, dvol, darea,  &
       kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,  &
       pnbi, beamDens, beamPress, beamFus, pbfuse, pbfusi, snBeamDD, &
       snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot, &
       etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot, &
       beamFusChTot, snDTTotal, snDDTotal, iflagnb)

  
  ni20=0.0

  !Set total current to J_nb
  !Note need to reverse order of array again and use <B>/<B^2>
  ! from appropriate flxsur
  do i=1,ncon

     J_nb(ncon-i+1) = dble(jnbTot(i))*bav(ncon-i+1)/bsqav(ncon-i+1)
	
		 print*, beamDens(ncon-i+1)
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

      psi = psiv(ncon-con+1)
           

      write(50,30) psi/umax, J_nb(ncon-con+1), pnbe(con), pnbi(con), beamDens(con)

   end do
30 format(8e14.6)

   close(50)

  print*, 'Charge beam fusion ', beamFusChTot
	print*, 'Total beam target fusion ', beamFusTot

  !print*, 'Electron Average Temp ', teav, q' and density ', neav
end subroutine nbicur

