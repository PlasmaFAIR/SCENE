subroutine error_msg(msg, flag)

  use param
  implicit none

  character(len=*), intent(in) :: msg
  integer, intent(in) :: flag
  integer :: unit

  unit = 10
  
  open(unit=unit, file='scene.error', status='unknown', iostat=ios)

  write(unit, *) flag
  write(unit, *) msg

  close(unit)

end subroutine error_msg
