module equil_data

  implicit none
!
  double precision, dimension(:,:), allocatable :: u  ! (psi-edge-psi)
  integer, dimension(:,:), allocatable :: ixout  ! =1 if mesh point in plasma
  integer, dimension(:,:), allocatable :: idout  ! =1 if cell centre in plasma
