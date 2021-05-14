module balpar

! Ballooning calculation variables

  implicit none
!
  integer nbal      ! set to 1 for single ballooning calculation >0 for all
                    ! flux surfaces
  integer ibal      ! The flux surface analysed for ballooning stability
                    ! when nbal=1
  integer  nturns    ! number of periods in 2*pi along field line for
                    ! ballooning calculation (nturns in positive and negative
                    ! directions)
  integer nchi0     ! number of chi0 values tested, evenly spaced between
                    ! -pi and pi (chi0=0 if nchi0=1)
  integer nchi      ! number of mesh points along field line for Runge-Kutte
                    ! calculation
  double precision chi0val  ! value of chi0 if nchi0=1 (0=outboard mid-pane)

     double precision, dimension(:), allocatable:: eq1,eq2,eqd1,chi,fb
     double precision, dimension(:), allocatable:: lambda
     double precision lam,lamges


 end module balpar
