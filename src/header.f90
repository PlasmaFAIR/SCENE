!> SPDX-License-Identifier: BSD-3-Clause
!> author: Peter Hill
!> ---
!>
!> Contains functions for creating a standard header for files and output.
!>
!> Adapted from GS2, original author: Peter Hill
module header_mod
  use fuuid4, only : uuid_len
  implicit none

  !> A header for files and output.
  !>
  !> This is a type rather than a function so that we can use a
  !> consistent date/time and UUID across different files by reusing
  !> the same object.
  !>
  !> Call `header_type::to_string` to get the header as text.
  type :: header_type
    character(:), allocatable :: date_time
    character(len=uuid_len) :: run_uuid
  contains
    procedure :: to_string => header_to_string
  end type header_type

  interface header_type
    module procedure make_header
  end interface header_type

contains

  !> Return the UUID for this run
  !>
  !> Generated when first called, guaranteed not to change during the
  !> simulation.
  function simulation_run_uuid()
    use fuuid4, only : generate_uuid, uuid_len
    character(len=uuid_len) :: simulation_run_uuid
    character(len=uuid_len), save :: run_uuid
    logical, save :: first_run = .true.

    if (first_run) then
      run_uuid = generate_uuid()
      first_run = .false.
    end if

    simulation_run_uuid = run_uuid
  end function simulation_run_uuid

  !> Return the current date and time in ISO8601 format:
  !>     YYYY-MM-DDThh:mm:ss.ssssZhh:mm
  function date_iso8601()
    character(:), allocatable :: date_iso8601
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    call date_and_time(date, time, zone)

    date_iso8601 = date(1:4) // "-" // date(5:6) // "-" // date (7:8) &
         // "T" // time(1:2) // ":" // time(3:4) // ":" // time(5:10) &
         // "Z" // zone(1:3) // ":" // zone(4:5)
  end function date_iso8601

  !> Constructor for `header_type`.
  !>
  !> Stores the date using [[date_is08601]] and the UUID using
  !> [[simulation_run_uuid]]
  function make_header() result(this)
    type(header_type) :: this

    this%date_time = date_iso8601()
    this%run_uuid = simulation_run_uuid()
  end function make_header

  !> Return a multiline string with a standard header of the form:
  !>
  !>     Created by SCENE at 2021-02-02T15:01:26.370Z+00:00
  !>     Run UUID: 36310A48-6A73-4941-9366-410C5731027A
  !>     <optional file description>
  !>
  !> If \p comment_character is passed, it is prepended to the start
  !> of each line. Note that it is used as-is, including any whitespace.
  function header_to_string(this, comment_character, file_description) result(header)
    use git_version, only: get_git_version
    character(:), allocatable :: header

    class(header_type), intent(in) :: this
    !> Optional character(s) to start each line. Whitespace is not
    !> stripped or added. Default is empty string
    character(*), optional, intent(in) :: comment_character
    !> Optional file description. Treated as a single line, that is,
    !> new lines in this string will not begin with
    !> [[comment_character]]
    character(*), optional, intent(in) :: file_description

    ! Actual comment character to use, may be empty
    character(:), allocatable :: comment_char
    ! Note that description_line will contain the comment character if
    ! present
    character(:), allocatable :: description_line
    ! A literal new line `\n`, because Fortran
    character(*), parameter :: newline = new_line('a')

    if (present(comment_character)) then
      comment_char = comment_character
    else
      comment_char = ""
    end if

    if (present(file_description)) then
      ! If we add more lines to the header, we should add a newline here
      description_line = comment_char // file_description
    else
      description_line = ""
    end if

    header = comment_char // "Created by SCENE at " // this%date_time // newline &
         // comment_char // "SCENE version: " // get_git_version() // newline &
         // comment_char // "Run UUID: " // this%run_uuid // newline &
         // description_line

  end function header_to_string

end module header_mod
