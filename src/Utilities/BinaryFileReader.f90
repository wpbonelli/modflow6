module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream
  use BinaryFileHeaderModule, only: BinaryFileHeaderType, &
                                    BinaryFileHeaderWrapperType

  public :: BinaryFileReaderType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    type(BinaryFileHeaderType) :: header
    type(BinaryFileHeaderType) :: headernext
    integer(I4B) :: total = 0
    integer(I4B) :: current = 0
    logical(LGP) :: indexed = .false.
    logical(LGP) :: endoffile = .false.
    class(BinaryFileHeaderWrapperType), allocatable :: headers(:)
  contains
    procedure(build_index_if), deferred :: build_index
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_if(this, header, success, iout)
      import BinaryFileReaderType, BinaryFileHeaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      class(BinaryFileHeaderType), intent(out) :: header
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
    end subroutine read_record_if
  end interface

  abstract interface
    subroutine build_index_if(this, iout)
      import BinaryFileReaderType, I4B
      class(BinaryFileReaderType), intent(inout) :: this
      integer(I4B), intent(in), optional :: iout
    end subroutine build_index_if
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
