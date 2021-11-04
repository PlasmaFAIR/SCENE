module nbi_mod
  implicit none
contains
  subroutine nbicur()
    error stop "Attempting to use NBI but SCENE has been compiled without NBEAMS support.&
         & Please reconfigure with -DSCENE_HAS_NBEAMS=on"
  end subroutine nbicur
end module nbi_mod
