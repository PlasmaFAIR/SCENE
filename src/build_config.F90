!> Compile-time parameters describing the built configuration
module build_config
  implicit none

#if SCENE_ENABLE_DOUBLE
  logical, parameter :: scene_enable_double = .true.
#else
  logical, parameter :: scene_enable_double = .false.
#endif

#if SCENE_USE_GHOST
  logical, parameter :: scene_use_ghost = .true.
#else
  logical, parameter :: scene_use_ghost = .false.
#endif

#if SCENE_USE_NBEAMS
  logical, parameter :: scene_use_nbeams = .true.
#else
  logical, parameter :: scene_use_nbeams = .false.
#endif

contains

  function formatted_build_config() result(output)
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
    character(len=:), allocatable :: output

    character(len=*), parameter :: nl = new_line('a')
    integer, parameter :: max_len = 2000

    allocate(character(len=max_len)::output)

    write(output, '(a, a)') "Built with ", compiler_version()
    write(output, '(a, a)') trim(output) // nl // "Build flags: ", compiler_options()
    write(output, '(a)') trim(output) // nl // nl // "Current build configuration: "
    write(output, '(a, l)') trim(output) // nl // "SCENE_ENABLE_DOUBLE: ", scene_enable_double
    write(output, '(a, l)') trim(output) // nl // "SCENE_USE_GHOST: ", scene_use_ghost
    write(output, '(a, l)') trim(output) // nl // "SCENE_USE_NBEAMS: ", scene_use_nbeams

    output = trim(output)
  end function formatted_build_config
end module build_config

