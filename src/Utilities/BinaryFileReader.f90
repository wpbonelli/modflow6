module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream

  public :: BinaryFileHeaderType, BinaryFileReaderType

  type :: BinaryFileHeaderType
    integer(I4B) :: pos = 0
    integer(I4B) :: kper, kstp = 0
    real(DP) :: delt, pertim, totim = 0.0_DP
  contains
    procedure :: get_str
    procedure :: reset => reset_header
  end type BinaryFileHeaderType

  type, abstract :: BinaryFileIndexType
    integer(I4B) :: inunit
    logical(LGP) :: built = .false.
    integer(I4B) :: total = 0
    integer(I4B) :: current = 0
    class(BinaryFileHeaderType), allocatable :: headers(:)
  contains
    procedure(read_header_if), deferred :: read_header
    procedure :: reset
  end type BinaryFileIndexType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    type(BinaryFileHeaderType) :: headernext
    logical(LGP) :: endoffile
  contains
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
  end type BinaryFileReaderType

  abstract interface
    subroutine read_header_if(this, eof)
      import BinaryFileIndexType
      import I4B, LGP
      class(BinaryFileIndexType), intent(inout) :: this
      logical(LGP), intent(out) :: eof
    end subroutine read_header_if
  end interface

  abstract interface
    subroutine read_record_if(this, success, iout)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
    end subroutine read_record_if
  end interface
contains

  !> @brief Get a string representation of the header.
  function get_str(this) result(str)
    class(BinaryFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str

    write (str, '(*(G0))') &
      'Binary file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ')'
    str = trim(str)
  end function get_str

  subroutine reset(this)
    class(BinaryFileIndexType) :: this
    this%built = .false.
    this%total = 0
    this%current = 0
    if (allocated(this%headers)) deallocate(this%headers)
  end subroutine reset

  subroutine reset_header(this)
    class(BinaryFileHeaderType) :: this
    this%pos = 0
    this%kper = 0
    this%kstp = 0
    this%delt = 0.0_DP
    this%pertim = 0.0_DP
    this%totim = 0.0_DP
  end subroutine reset_header

  !> @brief Peek to see if another record is available.
  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat

    if (.not. this%endoffile) then
      read (this%inunit, iostat=iostat) this%headernext%kstp, this%headernext%kper
      if (iostat == 0) then
        call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
      else if (iostat < 0) then
        this%endoffile = .true.
      end if
    end if
  end subroutine peek_record

end module BinaryFileReaderModule
