c/    This program is a driver to test the NB module nbeams. It is not
c/    part of the nbeams package. Its purpose is to create a background
c/    plasma (geometry, MHD and plasma profiles) and then call the NB
c/    routines to calculate various neutral beam heating and current drive
c/    parameters.
c/
c/    The program reads the input file 'inbeams.dat' which contains a
c/    namelist with various input variables, calls the main routine of
c/    the NB module (calcBeams) and prints out some information to the
c/    file 'outbeams.dat'
c/
c/    Written by John Mandrekas, GIT, for the NTCC 08/29/00
c/
c/    Description of some of the input variables in the namelist
c/    (Only the local variables are included here. For a description
c/    of the variables in the calcBeams list see the extensive comments
c/    in the beginning of the calcbeams.f routine)
c/
c/    e0         - elongation at center
c/    ea         - elongation at edge
c/    shft0      - shift of magnetic axis (normalized to minor radius)
c/    teavg      - Density-weighted volume average electron temperature (keV)
c/    tea        - Electron temperature at the edge (separatrix) (keV)
c/    tiavg      - Density-weighted volume average ion temperature (keV)
c/    tia        - Ion temperature at the edge (separatrix) (keV)
c/    alfat      - Temperature profile exponent
c/    denav      - Average electron density (10^20 /m^3)
c/    edgavr     - edge-to-average ratio, ne(a)/<ne>
c/    alfan      - Density profile exponent
c/    fion       - Concentration of each ion species (ni/ne)

c/
      implicit none
      include 'nbparams.inc'

c/    Variables in the calcBeams variable list:
c/    ----------------------------------------
      integer iflag, ie, ib, nbeams, inbfus, n, nion, maxiter
      integer nbshape(maxBeams), nbptype(maxBeams)

      real amb, zbeam, r0, a, b0, volp, nbcurTot, etanbTot, beamBeta,
     .     pNBAbsorbTot, pNBLossTot, beamFusTot, beamFusChTot,
     .     snDDTotal, snDTTotal

      real ebeam(maxBeams), pbeam(maxBeams), rtang(maxBeams),
     .     bwidth(maxBeams), bheigh(maxBeams), bgaussR(maxBeams),
     .     bgaussZ(maxBeams), bzpos(maxBeams), pwrfrac(3, maxBeams),
     .     aion(maxIons), zion(maxIons), ne20(mxrho),
     .     ni20(mxrho,maxIons), tekev(mxrho), tikev(mxrho), zeff(mxrho),
     .     rnorm(mxrho), vprime(mxrho), dvol(mxrho), darea(mxrho),
     .     kappa(mxrho), dkappa(mxrho), shafr(mxrho), dshafr(mxrho),
     .     hofr(mxrho,3,maxBeams), shinethru(3,maxBeams), jnbTot(mxrho),
     .     pnbe(mxrho), pnbi(mxrho), beamDens(mxrho), beamPress(mxrho),
     .     pbfuse(mxrho), pbfusi(mxrho), snBeamDD(mxrho),beamFus(mxrho),
     .     snBeamDT(mxrho), nbcur(maxBeams), etanb(maxBeams),
     .     gammanb(maxBeams), pNBAbsorb(maxBeams), pNBLoss(maxBeams),
     .     l31(mxrho), srcfast(mxrho, 3, maxBeams), psiv(mxrho),
     .     pitchangl(mxrho, 3, maxBeams), beamVel(mxrho),
     .     jnbfast(mxrho), pnbbm(mxrho, maxBeams)

c/    Local variables:

      integer i, nbout, con, ncon, ios

      namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape,
     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
     .     bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0

      real, dimension(:), allocatable :: temp, tempion

      character(len=15) runname

      amb = 2
      zbeam = 1
      inbfus = 1
      maxiter = 2
      nbout = 20
      open (nbout, file = 'outbeams.dat', status = 'unknown')
      
      write(6,*) 'What is the name of the *.in input file?'
      write(6,*) '(without the .dat suffix):'
      read(5,*) runname

      open(unit=53,file='nbi.dat', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'Error reading in nbi.dat file'
         stop
      end if
  
      read(53,*)
      read(53,*)
      
      read(53,*) nbeams

      allocate(temp(nbeams))
      
      read(53,*)
      read(53,*) temp
      ebeam(:nbeams) = temp
      
      read(53,*)
      read(53,*) temp
      pbeam(:nbeams) = temp
      
      
      read(53,*)
      read(53,*) temp
      rtang(:nbeams) = temp
      
      
      read(53,*)
      read(53,*) temp
      bzpos(:nbeams) = temp
      

      read(53,*)
      read(53,*) temp
      nbshape(:nbeams) = int(temp)
      
      
      read(53,*)
      read(53,*) temp
      nbptype(:nbeams) = int(temp)
      
      
      read(53,*)
      read(53,*) temp
      bwidth(:nbeams) = temp


      read(53,*)
      read(53,*) temp
      bheigh(:nbeams) = temp
      
      
      read(53,*)
      read(53,*) temp
      bgaussR(:nbeams) = temp
      
      
      read(53,*)
      read(53,*) temp
      bgaussZ(:nbeams) = temp
      
      
      read(53,*)
      read(53,*) pwrfrac(:,1)
      
      if (nbeams .eq. 2)  pwrfrac(:,2)=pwrfrac(:,1)
      close(53)

      open(unit=55,file=trim(runname)//'.nbeams',status='unknown',
     . iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'Error opening nbeams file'
         stop
      end if
      
      read(55, *) nion, r0, a, b0, volp, ncon

      allocate(tempion(nion*2))
      
      read(55, *) tempion
      aion(:nion) = tempion(:nion)
      zion(:nion) = tempion(nion+1:)
      
      do con=1,ncon
         read(55,*) rnorm(con), vprime(con), dvol(con), darea(con),
     .        kappa(con), dkappa(con), shafr(con), dshafr(con),
     .        zeff(con), ne20(con),
     .        ni20(con, :nion), tekev(con), tikev(con), psiv(ncon-con+1)
      end do
      close(55)
  
c/    Call the beams package:

      call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus,
     .   rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ,
     .   bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev,
     .   tikev, zeff, r0, a, b0, volp, ncon, rnorm, vprime, dvol, darea,
     .   l31, srcfast,
     .   kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,
     .   pnbi, pnbbm, beamDens, beamVel, beamPress, beamFus, jnbfast, 
     .   pbfuse, pbfusi, snBeamDD,
     .   snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot,
     .   etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot,
     .   beamFusChTot, snDTTotal, snDDTotal, pitchangl, iflag)

c/    Write output quantities to file nbout:


c/    Global parameters
c/    -----------------
      write (nbout, 2000)
      write (nbout, 2100)
      write (nbout, 2200) pNBAbsorbTot, pNBLossTot, nbcurTot,
     .   etanbTot, beamBeta
      write (nbout, '(1x)')

      write (nbout, 2300)
      write (nbout, 2400)
      write (nbout, 2410) (pNBAbsorb(ib), ib = 1, nbeams)
      write (nbout, 2415) (pNBLoss(ib), ib = 1, nbeams)
      write (nbout, 2420) (nbcur(ib), ib = 1, nbeams)
      write (nbout, 2430) (etanb(ib), ib = 1, nbeams)
      write (nbout, 2440) (gammanb(ib), ib = 1, nbeams)
      write (nbout, '(1x)')
      write (nbout, 2450)
      write (nbout, 2460) (shinethru(1,ib), ib = 1, nbeams)
      write (nbout, 2470) (shinethru(2,ib), ib = 1, nbeams)
      write (nbout, 2480) (shinethru(3,ib), ib = 1, nbeams)
      write (nbout, '(1x)')

      if (inbfus.NE.0) then
         write (nbout, 2500)
         write (nbout, 2600) beamFusTot, beamFusChTot, snDTTotal,
     .      snDDTotal
      endif

      write (nbout, '(1x)')

c/    Write the deposition profile for each beamline and energy group:
c/    ---------------------------------------------------------------
      write (nbout, *) 'Neutral Beam Deposition Profiles'
      write (nbout, *) '--------------------------------'
      do ib = 1, nbeams
         write (nbout, 1000) ib
         do i = 1, ncon
            write (nbout,1100) rnorm(i), (hofr(i,ie,ib),ie=1,3)
         enddo
         write (nbout, '(1x)')
      enddo

c/    Write heating and current drive profile information:
c/    ----------------------------------------------------
      write (nbout, 3000)
      do i = 1, ncon
         write (nbout, 3100) rnorm(i), jnbTot(i), pnbe(i), pnbi(i),
     .   beamDens(i), beamPress(i), beamFus(i)
      enddo


 1000 format (1x, 'Beam ', i2/
     .   5x, 'rho', 6x, 'hofr_1', 7x, 'hofr_2', 7x, 'hofr_3')
 1100 format (1x, f8.5, 3(1x,e12.5))

 2000 format (1x, 'Global NB Heating and Current Drive Parameters',/
     .        1x, '----------------------------------------------')
 2100 format (1x, 'TOTALS:')
 2200 format (1x, 'Total Absorbed Power      = ', f9.4, ' MW'/
     .        1x, 'Total Lost Power          = ', f9.4, ' MW'/
     .        1x, 'Total NB Driven Current   = ', e12.5,' A'/
     .        1x, 'Total NBCD Efficiency     = ', e12.5,' A/W'/
     .        1x, 'Total Beam Beta           = ', e12.5)

 2300 format (1x, 'Information per Beamline and Energy Group:'/
     .        1x, '-----------------------------------------')
 2400 format (23x,'Beam 1', 10x, 'Beam 2', 10x, 'Beam 3'/
     .        23x,'------', 10x, '------', 10x, '------')
 2410 format (1x, 'Absorbed Power   ', 1x, e12.5, 2(4x, e12.5))
 2415 format (1x, 'Lost Power       ', 1x, e12.5, 2(4x, e12.5))
 2420 format (1x, 'NB driven current', 1x, e12.5, 2(4x, e12.5))
 2430 format (1x, 'NBCD efficiency  ', 1x, e12.5, 2(4x, e12.5))
 2440 format (1x, 'NBCD gamma       ', 1x, e12.5, 2(4x, e12.5))
 2450 format (1x, 'Beam Shinethrough')
 2460 format (1x, '   energy group 1', 1x, e12.5, 2(4x, e12.5))
 2470 format (1x, '   energy group 2', 1x, e12.5, 2(4x, e12.5))
 2480 format (1x, '   energy group 3', 1x, e12.5, 2(4x, e12.5))

 2500 format (1x, 'Beam-Target Parameters'/
     .        1x, '----------------------')
 2600 format (1x, 'Total Beam-Target Fusion Power   = ', e12.5, ' MW'/
     .        1x, 'Total Power to Charged Particles = ', e12.5, ' MW'/
     .        1x, 'Total DT Neutron Rate            = ', e12.5, ' n/s'/
     .        1x, 'Total DD Neutron Rate            = ', e12.5, ' n/s')

 3000 format (3x, 'rho', 5x, 'jnbtot', 6x, 'pNBe', 8x, 'pNBi', 8x,
     .   'nbfast', 6x, 'pressb', 6x, 'pfusb')
 3100 format (1x, f6.4, 6(1x, e11.4))

      close(nbout)
      
!     Flx average quantities
      open(unit=48, file=trim(runname)//'_srcfast.dat', 
     . status='unknown', iostat=ios)
      if (ios.ne.0) then
         write(6,*) 'problem opening ',
     .    runname//'_srcfast.dat'
         stop
      end if

      write(48,*) 'Fast source'

      do i = 1,ncon
!     Some along energy components
         write(48,28) srcfast(i,:,:nbeams)
      end do

 28   format(8e13.6)
      close(48)
   
!Flx average quantities
      open(unit=50, file=trim(runname)//'_flxav.dat',
     . status='unknown', iostat=ios)
      if (ios.ne.0) then
         write(6,*) 'problem opening ',runname//'_flxav.dat'
         stop
      end if

      write(50,*) "Psi, NB cur (A/m^2), Pow to e (MW/m^3), ",
     .     "Power to i (MW/m^3), FI dens (/m^3), FI Vel (/m^2 s),",
     .     "flux sur vol (m^3), FI current (A/m^2)"
      do con=1,ncon
         write(50,30) psiv(ncon-con+1), jnbTot(con), pnbe(con),
     .   pnbi(con), beamDens(con), beamVel(con), dvol(con), 
     .   jnbfast(con), darea(con), pitchangl(con, 1, :2)
         
      end do
 30   format(11e12.5)
      
      close(50)
      
      stop
      end
