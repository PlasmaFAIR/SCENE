module ghost_interface
  implicit none

  interface
    subroutine filnam(nmfile)
      character, intent(in) :: nmfile(*)
    end subroutine filnam

    subroutine filon
    end subroutine filon

    subroutine grend
    end subroutine grend

    subroutine map(xarea1, xarea2, yarea1, yarea2)
      real, intent(in) :: xarea1, xarea2, yarea1, yarea2
    end subroutine map
  end interface
  
end module ghost_interface
