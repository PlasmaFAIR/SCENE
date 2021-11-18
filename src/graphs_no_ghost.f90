!> Stubs for the `flux_plot` module for build configurations without GHOST
module flux_plot
  implicit none
contains
  !> Setup GHOST grid file
  subroutine initialise_graphs
  end subroutine initialise_graphs

  !> Create all of the GHOST graphs
  subroutine make_all_graphs(debug)
    logical, optional, intent(in) :: debug
    if (present(debug)) continue
  end subroutine make_all_graphs
end module flux_plot
