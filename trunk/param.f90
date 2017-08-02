module param
  implicit none
!  Job title
  character*50 runname
  integer lrunname,ios
  character*19 title
!
!physical constants...
  double precision pi,eps0,mu0,bk  ! pi, epsilon_0, mu_0 and Boltzman constant
  double precision mp,me,eq           ! proton and electron masses, 
                                      ! electron charge
! ibdry=0 uses original Sykes parameterisation of boundary
! 27/4/04 ibdry=1 option introduced to permit more general shapes to be
! modelled
   integer ibdry
! plasma shape parameters: Geometric centre, inverse aspec ratio, 
!                          elongation,triangularity, quadracity separatrix
  double precision rcen,tokeps, elon,tri,quad,kval, amin
!     kval=0 gives D, kval=1 gives separatrix (DND)
!
!  r(theta,psi) for boundary: R=r*cos(th), Z=r*sin(th) defined boundary on 
!  flux surface psi. There are nflux surfaces
     double precision, dimension(:), allocatable:: theta
     double precision, dimension(:,:), allocatable:: rth
!  Fourier components of flux surface (for ibdr=2)
     double precision, dimension(:), allocatable:: fm,rfo,zfo
     integer nfm
!
  double precision step     ! mesh step length
!
  integer ipr               ! controls printout to screen
  integer igr               ! controls graphics
!
! Profile parameters
     integer ipswtch         ! switch for profile forms...see profs.f90
!  pressure (scaled by mu0)
     double precision ppow   ! power of polynomial part
     double precision p0,pn     ! central pressure
     double precision ppeak     ! pressure peaking factor
     double precision sclcur    ! corrects psi for numerical errors in current
     double precision pa     ! edge pressure
     double precision pped   ! scales as pedestal pressure
     double precision pedg   ! pedestal gradient/width
     double precision pscl   ! scales p profile to get required beta-p
!  main ion density (10^19 m^-3)
     double precision nipow   ! power of polynomial part
     double precision ni0     ! central ion density
     double precision nia     ! edge ion density
     double precision niped   ! scales as pedestal density
     double precision niedg   !  pedestal gradient/width
     double precision naa,nbb   ! for quad ni+tanh (ipswtch=9)
!  electron temperature (ev)
     double precision tpoe   ! power of polynomial part
     double precision ate,tpo1  ! used to get flatter core T if ipswtch=6
     double precision ane,npo1  ! used to get flatter core n if ipswtch=6
     double precision te0,ten     ! central electron temperature
     double precision tea     ! edge electron temperature
     double precision teped   ! scales as pedestal temperature
     double precision teedg   ! pedestal gradient/width
!  ion temperature (eV)
     double precision tpoi   ! power of polynomial part
     double precision ti0,tin     ! central ion temperature
     double precision tia     ! edge ion temperature
     double precision tiped   ! scales as pedestal temperature
     double precision tiedg   ! pedestal gradient/width
!  ff-prime profile parameters...also used to parameterise J_parallel
     double precision fpow,fpow1,fpow2,fpow3,fpow4
     double precision af0,af1,af2,af3,af4
!  Leading exponent of current profile
     double precision powj 
!  if itot=0 and the user specifies the driven current profile, then 
!  a part of that current can be used to 'fill in' the bootstrap current  near
!  the axis; this driven current is stored in exph2. The extent in normalised 
!  psi over which the 'filling' occurs in specified in psic
     double precision psic
!
!  Some global parameters.....
     double precision bpol   ! measure of beta-poloidal
     double precision pfac   ! An extra boost to the pressure
     double precision betan ! Target value of betan (bpol re-scaled if betan>0)
     double precision cur    ! required plasma current
     double precision rodi   ! TF rod current
     double precision paux   ! auxiliary heating power
     double precision zm     ! charge on main ion species
     double precision zmai   ! mass of main ion species (cf p mass)
     double precision dil    ! ion dilution fraction: used only for fusion
                             ! power calculation when imp=0 is specified
     integer nbet       ! used to count number of beta iterations (if betan>0)
     
!

!
!  Switches for type of equilibrium
     integer imat         ! uses Hirshman formula if IMAT=0, or
                          ! Hirshman-Sigmar for IMAT=1
     integer imp          ! Set to 1 (0) to (not) include impurity ion species
     integer itot         ! If 1 derives equilibrium from ff-prime and 
                          ! p-prime profiles; If 0 derives equilibrium from
                          ! specified driven current profile
     integer neo          ! If set to 1, and ITOT set to 0, derives the 
                          ! equilibrium with neoclassical current profile
                          ! If set to 2, then ibv is set to 1, to indicate 
                          ! Bv ramp-induced toroidal electric field
     integer ibv
     double precision rpk,rm  ! R of peak Bz-dot and zero Bz-dot
     integer nco          ! Set to 1 (0), to (not) include collisions 
                          ! If set >1 then uses Shaing model for collisions
                          ! at tight aspect ratio (for bootstrap current)

!Fast particle parameters
     integer fast         ! Toggles fast particles off/on when fast=0/1.
     double precision fastb !Beta due to fast alpha (input as second impurity)

!  Equilibrium analysis...
     integer ncon        ! Required number of flux surfaces to analyse eqbm
     integer npts        ! Number of mesh points on each flux surface
     integer npass       ! Maximum number of equilibrium iterations in 
                         ! attempt to converge on required current profile
     integer icontour    ! When set to zero uses linear interpolation of 
                         ! psi(R,Z) grid: fine for R, Z and forst derivatives
                         ! Set to one to smooth the contours to give smooth 
                         ! second derivatives
                         ! Set to 2 to use Colin's contouring algorithm, which 
                         ! interpolates the psi mesh using a bi-cibic
                         ! spline
!
! Parameters used in Grad-Shafranov solver
     double precision errit ! error on convergence of Grad-Shafranov solver
     double precision errffd ! error on convergece of ff-prime (for itot=0)
     double precision omega ! control rate of convergence of G-S solver
     double precision frac ! control rate of convergence of G-S solver
     double precision inscl,scl ! used to get current correct
     integer ninner         ! control rate of convergence of G-S solver
     integer nouter  ! Maximum allowed number of iterations to solve G-S
     integer icont   ! If 1, uses a 'warm' start for the eqbm solver
!
!
! Read and write units:
     integer nread
     integer nw
!
!  Impurity ion information
     integer nimp      ! number of impurity ions
!  Impurity arrays are defined in an analagous way to n_e and Te profiles
     double precision, dimension(:), allocatable:: ztpow,zt0,zta,ztped,ztedg
     double precision, dimension(:), allocatable:: znpow,zn0,zna,znped,znedg
     double precision, dimension(:), allocatable:: zmas
     integer, dimension(:), allocatable:: iz

!----------------------------------------------------------------------
!
!  This ends the input variables.
!  We now list the equilibrium mesh variables
!------------------------------------------------------------------------

   double precision dr,dz    ! mesh size in R and Z directions (m)
   double precision r1,r2,rs,zs ! Inboard, outboard and 'corner' R-values
                                   ! Z-value of 'corner of D-shaped plasma (m)
   integer nr,nz         ! no. radial mesh points in R and Z directions
   integer nsym          ! labels the z-coordinate of the symmetry plane
   integer ipass         ! counts number of eqbm iterations
!
!  The arrays giving the major radius, and height of the eqbm mesh points
   double precision, dimension(:), allocatable:: r,z
!  The array containing the normalised flux calculated in the equbm
!  NOTE**** u=psi_a-psi is max on axis, and zero at the edge
  double precision, dimension(:,:), allocatable:: u  ! (psi-edge - psi)
!
! These arrays indicate points inside the plasma (when they are =1)
  integer, dimension(:,:), allocatable :: ixout  ! =1 if mesh point in plasma
  integer, dimension(:,:), allocatable :: idout  ! =1 if cell centre in plasma
  integer, dimension(:), allocatable :: nix,niy
  integer ntot

!  rcoord and zcoord label the R and Z coordinates of the centre of each
!  cell, and btheta is the poloidal magnetic field there
    double precision, dimension(:), allocatable:: rcoord,zcoord
    double precision, dimension(:,:), allocatable:: btheta,brcoord,bzcoord
!
!---------------------------------------------------------------------
!
!  Error tracking variables
!
!--------------------------------------------------------------------
    double precision residl,curtot
!
!---------------------------------------------------------------------
!
!  There now follows the list of equilibrium varables calculated in SCENE
!
!----------------------------------------------------------------------

   double precision r0     ! magnetic axis position
   double precision const  ! vacuum part of f(psi)
   double precision umax   ! u on axis, or psi at edge
   double precision surar  ! plasma surface area
!
!---------------------------------------------------------------------
!
!  And now we have the flux surface variables...
!
!  The poloidal flux (different definition to U...defined to be zero on axis)
! and increasing going outwards; eps = inverse aspect ratio of flux surfaces
!  First array element is plasma edge
  double precision, dimension(:), allocatable:: psiv,epsv
!  R and Z coordinates, and poloidal field
  double precision, dimension(:,:), allocatable:: rpts,zpts,bppts
!  flux surface circumference
   double precision, dimension(:), allocatable:: circumf 

! There follows a list of flux-surface averaged variables:
!  First we have q, q-prime and q-double-prime
  double precision, dimension(:), allocatable:: sfac,qp,qpp
!  Then we have   <B^2> and V''
  double precision, dimension(:), allocatable:: bsqav,vpp
!  And now: <R^2>, <R^{-1}>, <R^{-2}>, <R>, <1>
  double precision, dimension(:), allocatable:: rsqav,rinv,rsqinv,rav,rnorm
!  <R^{-4}*B_p^{-2}>  <B>
  double precision, dimension(:), allocatable:: avbli,bdl
!  collisionality variables
  double precision, dimension(:), allocatable:: cnue,tnue,cnui,tnui
!  trapped ptcle fraction, ftrap, and something used for fast ptcle b/strap
  double precision, dimension(:), allocatable:: ftrap,ftrapd
!  Neoclassical resistivity enhancement factor
  double precision, dimension(:), allocatable:: sighhb
!  bootstrap current <J.B>/sqrt<B^2>: thermal in bsj and alpha in ajdotb
  double precision, dimension(:), allocatable:: bsj,ajdotb
!
!  psiold is a normalised flux variable on which the current guess at
!  the ff' profile is evoluated (stored in gst)
 double precision, dimension(:), allocatable:: psiold,gst
!
!
!------------------------------------------------------------------------
!
!  Toroidal current profiles over the plasma cross-sectional area
!
!--------------------------------------------------------------------
!
!  Bootstrap current: thermal and alpha
  double precision, dimension(:,:), allocatable:: bsph,absph
!  Pfirsch-Schluter current
  double precision, dimension(:,:), allocatable:: psph
!  Diamagnetic current
  double precision, dimension(:,:), allocatable:: diph
!  Externally driven current (Ohmic + any auxiliary CD); exph2 is the 
!  part of the externally driven current which fills in the bootstrap
!  current on axis
  double precision, dimension(:,:), allocatable:: exph,exph2
!  Total current profile used in G-S eqn
  double precision, dimension(:,:), allocatable:: gradj
!  Spitzer current profile
  double precision, dimension(:,:), allocatable:: spit
  ! Neutral beam current
  double precision, dimension(:,:), allocatable:: nbph
!
!  And the total current to go with these....
 double precision totbs,totps,totdi,totex,totex2,totgs,spiti,alfbs
! Loop voltage (total, Spitzer and nc  without b/strap)
 double precision vloop,vspit,vnobs
!
!
!---------------------------------------------------------------
! 
!   Neutral Beam parameters
! Beam paramters: width in r and z, tangency radius and beam energy
 double precision:: sig_r, sig_z, R_t,E_b,Z_beam
 double precision :: A_beam, I_0

 
!
!---------------------------------------------------------------
!
!  Global equilibrium parameters
!
!---------------------------------------------------------------
!
! Stored energy, confinement time and volume averaged Zeff
  double precision conft,taue,zeffav
! beta variables
  double precision betai,betap,beta,betexp,betlim
! Volume averaged electron temperatrue and density
  double precision avt,avel
! Line averaged electron density (in 10**19 m**-3), and normalised to Greenwald
  double precision nebar,negw
! Plasma volume, cross-sectional area and internal inductance
  double precision vol,area,rli3,rli2
! Confinement time relative to scaling laws...
  double precision hipb98y1,hipb98y2
!
!----------------------------------------------------------
!
! Radiated power fraction for 3 ion species (1=main ion+ 2 impurities)
!
  double precision powtot(3),powsum
!
!----------------------------------------------------------
!
!----------------------------------------------------------
!
!  Quantities for burning device designs
!
!------------------------------------------------------------
!
!  Alpha heating power and He stored energy and confinement time
 double precision pfus,confa,tauh
!
!----------------------------------------------------------
!
!  Arrays for reading in density and temperature profiles from mesh
!  Equal-spaced in psi for speed
  double precision, dimension(:), allocatable:: psi_m,ne_m,te_m
! Pressure and ffp specified for ipswitch=-1
  double precision, dimension(:), allocatable:: p_m,ffp_m
!  The number of points
  integer nterp
!
! Rob wanted some data on total core current drive and edge current drive:
  double precision jcore,jtotal
!
!  electron viscosity coefficients
  double precision, dimension(:), allocatable:: mu1,mu2,mu3,nuee,nuei
!  fast alpha pressure
  double precision, dimension(:), allocatable:: fastp
!  fast particle fraction (vol average)
  double precision palpha
!
! If ifast=1 the bootstrap current includes the alpha-particle contribution
   integer ifast
!-----------------------------------------------------------
!
!  Store dnu/dpsi, which is useful for ballooning stability stuff
   double precision, dimension(:,:), allocatable:: nup,dbsqdpsi


 end module param


     
