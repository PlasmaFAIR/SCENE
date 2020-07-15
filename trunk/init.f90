  subroutine init
! ***************

!------------------------------------------------------------------
!  Initialises physical constants and the equilibrium parameters
!------------------------------------------------------------------
      use param
      implicit none

!
! initialise physical constants
      pi=4.0d0*atan(1.0d0)
      mu0=pi*4.0d-7
      bk=1.602d-19
      eps0=8.8542d-12
      mp=1.6726d-27
      me=9.1094d-31
      eq=1.6022d-19
! Define read and write units
      nread=11
      nw=6
! set default parameters
      call default
! read in desired parameters
      call parset
!!$! test density profile
!!$      umax=2.
!!$      dpsi=umax/(ncon-1.)
!!$      do i=1,ncon
!!$        psi=(i-1)*dpsi
!!$!        write(6,*)' psi=',psi,' ni=',densi(psi,1,0),' nid=',densi(psi,1,1)
!!$!        write(6,*)' psi=',psi,' ni2=',densi(psi,2,0),' nid2=',densi(psi,2,1)
!!$        if ((i.gt.1).and.(i.lt.ncon)) then
!!$           write(6,*)' psi=',psi,' nid=',densi(psi,2,1),            &
!!$               ' num=',(densi(psi+dpsi,2,0)-densi(psi-dpsi,2,0))/(2.*dpsi)
!!$           write(6,*)' nidd=',densi(psi,2,2),            &
!!$               ' num=',(densi(psi+dpsi,2,0)+densi(psi-dpsi,2,0)   &
!!$                        -2.*densi(psi,2,0))/(dpsi**2)
!!$           write(6,*)' psi=',psi,' Ted=',tempe(psi,1),            &
!!$               ' num=',(tempe(psi+dpsi,0)-tempe(psi-dpsi,0))/(2.*dpsi)
!!$           write(6,*)' Tedd=',tempe(psi,2),            &
!!$               ' num=',(tempe(psi+dpsi,0)+tempe(psi-dpsi,0)   &
!!$                        -2.*tempe(psi,0))/(dpsi**2)
!!$        end if
!!$      end do
!!$      write(6,*)' ane=',ane,' npo1=',npo1
end subroutine init

!*********************************************************************

  subroutine default
! ******************
!
!--------------------------------------------------------------------
!  Sets default parameter set
!--------------------------------------------------------------------
!
   use param
   use balpar
   implicit none

!  The number of points
!
      rcen=0.2
      ibdry=0
      nfm=3
      tokeps=0.9
      elon=1.5
      tri=0.3
      quad=0.
      kval=0.
      step=0.01
      ipr=0
      ipswtch=0
      ppow=2.
      p0=1.
      pa=0.
      pped=0.
      ape=2.1
      ppo1=3.0
      pedg=1.
      nipow=1.
      ni0=5.
      nia=0.
      naa=0.
      nbb=0.
      niped=0.
      niedg=1.
      tpoe=1.
      te0=1.
      tea=0.
      ifast=0
      teped=0.
      teedg=1.
      tpoi=1.
      ate=2.1
      tpo1=3.0
      ti0=1.
      tia=0.
      tiped=0.
      tiedg=1.
      fpow=0.1
      fpow1=1.
      fpow2=10.
      af0=1.
      af1=0.
      af2=0.
      fpow3=1.
      fpow4=10.
      af3=0.
      af4=0.
      powj=1.
      psic=0.
      bpol=1.
      npo1=3.
      ane=2.
      pfac=1.
      betan=-1.
      cur=100000.
      rodi=50000
      paux=0.
      zm=1.
      zmai=2.
      imat=0
      dil=1.
      imp=0
      icontour=0
      ibv=0
      itot=0
      fast=0
      neo=1
      nco=0
      ncon=15
      npts=200
      npass=10
      omega=1.6
      frac=0.1
      scl=1.0d8
      nouter=120
      ninner=4
      errit=1.0E-3
      errffd=0.01
      icont=1
!  Ballooning calculation variables
      ibal=0
      nbal=0
      nturns=10
      nchi=500
      nchi0=30
      chi0val=0.
      lamges=2.
      nbi = 0
      bmfus=0.

      !NB
      sig_r =0.08
      E_b = (/150,150/)
      sig_z = 0.08
      P_beam = (/10,40/)
      !I_0 = P_beam*1000./E_b
      P_frac = (/0.70,0.20,0.10/)
      R_t = (/0.8,1.1/)
      Z_beam = 0.0
      A_beam = 1
 end subroutine default
!
!**********************************************************************
!
      subroutine parset
!     *****************
!
!--------------------------------------------------------------------
!  Reads in user parameter set and does some compatibility checks
!---------------------------------------------------------------------
!
      use param
      use balpar
      implicit none
      character(len=4) word
      double precision val,rlo,zlo
      double precision press,psiax
      integer ival,i
      integer inr,inz
      integer ndatpt,nndatpt,ntdatpt
      integer nh,ndat
      double precision rmax,rmin,zmax,zmin,rval,zval
      double precision, dimension(:), allocatable:: psi_in,ne_in,te_in,work
      double precision, dimension(:), allocatable:: p_in,ffp_in
      double precision dpsi
      !
      logical :: debug

      debug = .false.
      
      write(6,*) 'What is the name of the *.in input file?'
      write(6,*) '(without the .dat suffix):'
      read(5,*) runname
      lrunname=index(runname,' ')-1
      open(unit=nread,file=runname(1:lrunname)//'.dat',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem opening file',runname(1:lrunname)//'.dat'
         stop
      endif
      read(nread,10)title
 10   FORMAT(A19)
! Read in equilibrium parameters
 60   read(nread,*)word,val
!  Avoid rounding error if VAL is an integer
      ival=int(val+0.1)
      if(val.LT.0.0) ival=int(val-0.1)
      if(word.eq.'rcen')rcen=val
      if(word.eq.'eps ')tokeps=val
      if(word.eq.'ibdr')ibdry=ival
      if(word.eq.'nfm ')nfm=ival
      if(word.eq.'elon')elon=val
      if(word.eq.'tri ')tri=val
      if(word.eq.'quad')quad=val
      if(word.eq.'kval')kval=val
      if(word.eq.'step')step=val
      if(word.eq.'ipr ')ipr=ival
      if(word.eq.'igr ')igr=ival
      if(word.eq.'fpow')fpow=val
      if(word.eq.'powj')powj=val
      if(word.eq.'fpo3')fpow3=val
      if(word.eq.'fpo4')fpow4=val
      if(word.eq.'af3 ')af3=val
      if(word.eq.'af4 ')af4=val
      if(word.eq.'fpo1')fpow1=val
      if(word.eq.'fpo2')fpow2=val
      if(word.eq.'af1 ')af1=val
      if(word.eq.'af0 ')af0=val
      if(word.eq.'af2 ')af2=val
      if(word.eq.'psic') psic=val
      if(word.eq.'ipsw')ipswtch=ival
      if(word.eq.'ppow')ppow=val
      if(word.eq.'p0  ')p0=val
      if(word.eq.'pa  ')pa=val
      if(word.eq.'pped')pped=val
      if(word.eq.'ppo1')ppo1=val
      if(word.eq.'ape ')ape=val
      if(word.eq.'pedg')pedg=val
      if(word.eq.'npow')nipow=val
      if(word.eq.'ni0 ')ni0=val
      if(word.eq.'nia ')nia=val
      if(word.eq.'naa ')naa=val
      if(word.eq.'nbb ')nbb=val
      if(word.eq.'nped')niped=val
      if(word.eq.'nedg')niedg=val
      if(word.eq.'ane ')ane=val
      if(word.eq.'npo1')npo1=val
      if(word.eq.'ate ')ate=val
      if(word.eq.'tpo1')tpo1=val
      if(word.eq.'tpoe')tpoe=val
      if(word.eq.'te0 ')te0=val
      if(word.eq.'tea ')tea=val
      if(word.eq.'tepd')teped=val
      if(word.eq.'tedg')teedg=val
      if(word.eq.'tpoi')tpoi=val
      if(word.eq.'ti0 ')ti0=val
      if(word.eq.'tia ')tia=val
      if(word.eq.'tipd')tiped=val
      if(word.eq.'tidg')tiedg=val
      if(word.eq.'cur ')cur=val
      if(word.eq.'paux')paux=val
      if(word.eq.'bpol')bpol=val
      if(word.eq.'pfac')pfac=val
      if(word.eq.'betn')betan=val
      if(word.eq.'dil ')dil=val
      if(word.eq.'rodi')rodi=val
      if(word.eq.'imat')imat=ival
      if(word.eq.'ifas')ifast=ival
      if(word.eq.'scl ')scl=val
      if(word.eq.'omeg')omega=val
      if(word.eq.'frac')frac=val
      if(word.eq.'zm  ')zm=val
      if(word.eq.'zmai')zmai=val
      if(word.eq.'ictr')icontour=ival
      if(word.eq.'imp ')imp=ival
      if(word.eq.'itot')itot=ival
      if(word.eq.'fast')fast=ival
      if(word.eq.'neo ')neo=ival
      if(word.eq.'nco ')nco=ival
      if(word.eq.'ncon')ncon=ival
      if(word.eq.'npts')npts=ival
      if(word.eq.'nout')nouter=ival
      if(word.eq.'ninn')ninner=ival
      if(word.eq.'npas')npass=ival
      if(word.eq.'erit')errit=val
      if(word.eq.'erfd')errffd=val
      if(word.eq.'icon')icont=ival
      if(word.eq.'ibal')ibal=ival
      if(word.eq.'nbal')nbal=ival
      if(word.eq.'ntrn')nturns=ival
      if(word.eq.'nchi')nchi=ival
      if(word.eq.'nci0')nchi0=ival
      if(word.eq.'chi0')chi0val=val
      if(word.eq.'lam ')lamges=val
      if(word.eq.'nbi ')nbi=ival
      if(word.eq.'fini')goto 70
      goto 60
 70   continue
      nturns=2*nturns
      if (chi0val.le.-pi) chi0val=chi0val+2.*pi
      if (chi0val.gt.pi) chi0val=chi0val-2.*pi
      inscl=scl
! Scale currents to be in Amps
      cur=cur*1.0d6
      rodi=rodi*1.0d6
! Set ibv ramp switch if neo=2
      if (neo.eq.2) then
        ibv=1
        neo=1
!!$        write(6,*)'Input R of peak bz-dot and R of bz-dot=0 (negative for constant BZ-dot)'
!!$        read(5,*)rpk,rm
      end if
! If ipswtch.eq.3 then read in grid of T and n data....Ti=Te assumed
      allocate(psi_m(ncon),ne_m(ncon),te_m(ncon),p_m(ncon))
      if (ipswtch.eq.3) then
        open(90,file='ntdat.dat')
        read(90,*)nndatpt,ntdatpt,nterp
        ndatpt=ntdatpt
        if (ndatpt.lt.nndatpt) ndatpt=nndatpt
        ndatpt=ndatpt+1
        deallocate(psi_m,ne_m,te_m)
        allocate(psi_in(nndatpt+1),ne_in(nndatpt+1),work(nndatpt+1))
        allocate(psi_m(nterp),ne_m(nterp),te_m(nterp),p_m(nterp))
! Equal spaced psi mesh
        dpsi=1./(nterp-1.)
        do i=1,nterp
          psi_m(i)=(i-1)*dpsi
        end do
        do i=2,nndatpt+1
          read(90,*)psi_in(i),ne_in(i)
        end do
        psi_in(1)=0.
        ne_in(1)=ne_in(2)
! Create fine-scale meshes, equal spaced in psi
        call spline1d(ne_m,psi_m,nterp,ne_in,psi_in,nndatpt+1,work)
        deallocate(psi_in,ne_in,work)
        allocate(psi_in(ntdatpt+1),te_in(ntdatpt+1),work(ntdatpt+1))
        do i=nterp,1,-1
          ne_m(i)=ne_m(i)/ne_m(1)
          write(6,*)' psi=',psi_m(i),' ne=',ne_m(i)
        end do
        do i=2,ntdatpt+1
          read(90,*)psi_in(i),te_in(i)
        end do
        psi_in(1)=0.
        te_in(1)=te0
! Create fine-scale meshes, equal spaced in psi
        call spline1d(te_m,psi_m,nterp,te_in,psi_in,ntdatpt+1,work)
        deallocate(psi_in,te_in)
      end if
!
!
! If ipswtch.eq.-1 then read in grid of p and ff' data
      if (ipswtch.eq.-1) then
        open(90,file='pdat.dat')
        read(90,*)nndatpt
        ndatpt=nndatpt
        nterp=ndatpt
        allocate(psi_in(ndatpt+1),p_in(ndatpt+1),ffp_in(ndatpt+1), &
                work(ndatpt+1))
        allocate(psi_m(nterp),p_m(nterp),ffp_m(nterp))
! Equal spaced psi mesh
        dpsi=1./(nterp-1.)
        do i=1,nterp
          psi_m(i)=(i-1)*dpsi
        end do
        do i=1,ndatpt
          read(90,*)psi_in(i),p_in(i),ffp_in(i)
        end do
! Create fine-scale meshes, equal spaced in psi
        call spline1d(p_m,psi_m,nterp,p_in,psi_in,ndatpt+1,work)
! Create fine-scale meshes, equal spaced in psi
        call spline1d(ffp_m,psi_m,nterp,ffp_in,psi_in,ndatpt+1,work)
        deallocate(psi_in,p_in,ffp_in)
      end if
!
!
! Put density parameters in SI
      ni0=ni0*1.0d19
      nia=nia*1.0d19
      niped=niped*1.0d19
      naa=naa*1.0d19
      nbb=nbb*1.0d19
! Check that number of poloidal points is odd...
      if (2*(npts/2).eq.npts) then
        write(6,*)'INPUT ERROR***npts must be odd'
        stop
      end if
! Quick check on consistency of edge parameters
      if ((imp.eq.0).and.(tea.lt.1.e-8).and.(tia.lt.1.e-8)) then
        write(6,*)'Profiles require non-zero edge temperature;'
        write(6,*)'set tea or tia non-zero (can be arbitrarily small)'
        stop
      end if
!!$      IF (ISO.EQ.1) THEN
!!$        tce=2.
!!$        IF (TCE.GE.1.) THEN
!!$          WRITE(6,*)'ERROR***TCE must be < 1...re-set in input file'
!!$          STOP
!!$        END IF
!!$      END IF
      if (imp.eq.1) then
!  Read in details of impurity ions
        read(nread,*)val
        val=val+0.1
        nimp=int(val)
        call allocimp
!  First element is main ion species
        iz(1)=int(zm)
        zmas(1)=zmai
        if (debug) write(6,*)' zmas1=',zmas(1),zmai
        ztpow(1)=tpoi
        zt0(1)=ti0
        zta(1)=tia
        ztped(1)=tiped
        ztedg(1)=tiedg
        znpow(1)=nipow
        zn0(1)=ni0
        zna(1)=nia
        znped(1)=niped
        znedg(1)=niedg
        if (nimp.gt.0) then
          do  i=1,nimp
            read(nread,*)iz(i+1),zmas(i+1)
            read(nread,*)ztpow(i+1),zt0(i+1),zta(i+1),ztped(i+1),ztedg(i+1)
            read(nread,*)znpow(i+1),zn0(i+1),zna(i+1),znped(i+1),znedg(i+1)
            zn0(i+1)=zn0(i+1)*1.0d19
            zna(i+1)=zna(i+1)*1.0d19
            znped(i+1)=znped(i+1)*1.0d19
          end do
        end if
      else
        nimp=0
        call allocimp
!  First element is main ion species
        iz(1)=int(zm)
        zmas(1)=zmai
        ztpow(1)=tpoi
        zt0(1)=ti0
        zta(1)=tia
        ztped(1)=tiped
        ztedg(1)=tiedg
        znpow(1)=nipow
        zn0(1)=ni0
        zna(1)=nia
        znped(1)=niped
        znedg(1)=niedg
      end if
!  If reading in mesh, calculate calculation boundary
      if (ibdry.eq.2) then
          nh=39
          open(unit=nh,file='bdy.txt', &
                status='unknown',iostat=ios)
          if(ios.ne.0) then
            write(6,*) 'problem opening bdy.txt for boundary R,Z'
            stop
          endif
          read(nh,*)ndat
          do i=1,ndat
            read(nh,*)rval,zval
            if (i.eq.1) then
              rmax=rval
              rmin=rval
              zmax=zval
              zmin=zval
              rs=rval
            else
              if (rmax.lt.rval) rmax=rval
              if (rmin.gt.rval) rmin=rval
              if (zmax.lt.zval) then
                zmax=zval
                rs=rval
              end if
              if (zmin.gt.zval) zmin=zval
            end if
          end do
          close(nh)
          if (abs(zmax+zmin).gt.1.0d-4) then
            write(6,*)'ERROR***input boundary not up-down symmetric'
            write(6,*)' zmax=',zmax,' zmin=',zmin
            stop
          end if
          rcen=(rmax+rmin)/2.
          tokeps=(rmax-rmin)/(rmax+rmin)
          elon=zmax/(rcen*tokeps)
          tri=asin((1.-rs/rcen)/tokeps)
      end if
!  magnetic axis guess...
      r0=rcen
!  vacuum part of f(psi)
      const=(mu0*rodi/(2.*pi))**2
!  calculate boundary parameters for sykes boundary...
      dr=step
      dz=step
      r1=rcen*(1.-tokeps)
      r2=rcen*(1.+tokeps)
      zs=rcen*tokeps*elon
      rs=rcen*(1.-tokeps*sin(tri))
      !rs=rcen*(1.-tokeps*tri)
!  Create mesh... 
      nr=int((r2-r1)/step+6)
      nz=int(2.*zs/step+6)
      inr=2*(nr/2)
! ensure nr and nz are odd
      if (inr.eq.nr) nr=nr+1
      inz=2*(nz/2)
      if (inz.eq.nz) nz=nz+1
! Allocate the mesh arrays...
      call allocmesh
! nsym labels z-coord of symmetry plane
      nsym=(nz+1)/2
! set up the box mesh
      rlo=rcen-((nr-1)/2)*dr
      zlo=-((nz-1)/2)*dz
      do i=1,nr
        r(i)=rlo+(i-1)*dr
      end do
      do  i=1,nz
        z(i)=zlo+(i-1)*dz
      end do
      if (abs(neo).eq.1) then
        if (itot.ne.0) then
          write(nw,*)'input warning***running cannot input total'
          write(nw,*)'current for neoclassical option'
          write(nw,*)'have reset itot=0 to specify external current'
          itot=0
        end if
      end if
!  Set the scaling parameters for temperature and pressure profiles
      if (tpoe.le.1.) then
        ten=0.
      else
        ten=(te0-tea-(teped-tea)*(exp(teedg)-exp(-teedg))/ &
                              (exp(teedg)+exp(-teedg))) &
           /(2**tpoe-1.-tpoe)
        if ((ipswtch.eq.6).or.(ipswtch.eq.7)) then
          ten=(te0-tea-(teped-tea)*(exp(teedg)-exp(-teedg))/ &
                              (exp(teedg)+exp(-teedg))) &
             /1.
          end if
      end if
      if (tpoi.le.1.) then
        tin=0.
      else
        tin=(ti0-tia-(tiped-tia)*(exp(tiedg)-exp(-tiedg))/ &
                              (exp(tiedg)+exp(-tiedg))) &
           /(2**tpoi-1.-tpoi)
      end if
      if (imp.eq.0) then
        if (ppow.gt.1.) then
        pn=(p0-pa-(pped-pa)*(exp(pedg)-exp(-pedg))/(exp(pedg)+exp(-pedg))) &
           /(2**ppow-1.-ppow)
        else
          pn=0.
        end if
      end if
      pscl=1.
      psiax=0.
      umax=1.
      pscl=abs(press(psiax,1)/(scl*bpol*pfac))
      if (debug) write(6,*)' pscl=',pscl

!  Allocate arrays for ff' iteration
      allocate( psiold(ncon),gst(ncon), J_nb(ncon), nbmom(ncon) )
!!$!  Run some checks on input parameters
!!$      if (tpoe.lt.1.) then
!!$        write(6,*)'INPUT ERROR****, must set tpoe>1'
!!$        stop
!!$      end if
!!$      if (tpoi.lt.1.) then
!!$        write(6,*)'INPUT ERROR****, must set tpoi>1'
!!$        stop
!!$      end if
!!$      if (ppow.lt.1.) then
!!$        write(6,*)'INPUT ERROR****, must set ppow>1'
!!$        stop
!!$      end if
!!$      if (nipow.lt.1.) then
!!$        write(6,*)'INPUT ERROR****, must set npow>1'
!!$        stop
!!$      end if
  end subroutine parset
!
!**********************************************************************
!
  subroutine allocimp
! *******************
!
! Allocates the impurity profiles arrays:
!
      use param
!
      allocate( ztpow(nimp+1),zt0(nimp+1),zta(nimp+1),ztped(nimp+1), &
                ztedg(nimp+1) )
      ztpow=0.; zt0=0.; zta=0.; ztped=0.; ztedg=0.
      allocate( znpow(nimp+1),zn0(nimp+1),zna(nimp+1),znped(nimp+1), &
                znedg(nimp+1))
      znpow=0.; zn0=0.; zna=0.; znped=0.; znedg=0.
      allocate( zmas(nimp+1),iz(nimp+1) )
      zmas=0.; iz=0
!
   end subroutine allocimp

!**********************************************************************
!
  subroutine allocmesh
! ********************
!
! Allocates the variable arrays on the R-Z eqbm mesh:
!
      use param
!
      allocate( r(nr), z(nz) )
      r=0.; z=0.
      allocate( u(nr,nz),ixout(nr,nz),idout(nr,nz) )
      u=0.; ixout=0; idout=0
!
   end subroutine allocmesh
