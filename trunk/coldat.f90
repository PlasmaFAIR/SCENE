module coldat
  implicit none
!
!  Module contains data relating to Hirshman-Sigmar bootstrap model

!  Number of species
     integer nspec
!  mass, temperature, charge and density on each
     double precision, dimension(:), allocatable:: sm,ts,ch,dn
! collision times, ratio of thermal velocities and viscosity coeffs
     double precision, dimension(:,:), allocatable:: coltau,xab,vis
!  friction coefficients
     double precision, dimension(:,:), allocatable:: rl11,rl12,rl21,rl22

!
  end module coldat
