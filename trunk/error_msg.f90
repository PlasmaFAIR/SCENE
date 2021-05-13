module scene_errors
  implicit none
contains
  subroutine error_msg(msg, flag)
    use param
    implicit none

    character(len=*), intent(in) :: msg
    integer, intent(in) :: flag
    integer :: unit

    open(newunit=unit, file='scene.error', status='unknown')

    write(unit, *) flag
    write(unit, *) msg

    close(unit)

  end subroutine error_msg
end module scene_errors
