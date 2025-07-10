module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream

  public :: BinaryFileHeaderType, &
            BinaryFileHeaderWrapperType, &
            BinaryFileReaderType

  type :: BinaryFileHeaderType
    integer(I4B) :: pos
    integer(I4B) :: kper, kstp
    real(DP) :: delt, pertim, totim
  contains
    procedure :: get_str
  end type BinaryFileHeaderType

  type :: BinaryFileHeaderWrapperType
    class(BinaryFileHeaderType), allocatable :: header
  end type BinaryFileHeaderWrapperType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    logical(LGP) :: indexed = .false.
    logical(LGP) :: endoffile = .false.
    integer(I4B) :: current = 1
    type(BinaryFileHeaderWrapperType), allocatable :: headers(:)
    type(BinaryFileHeaderType) :: header
    type(BinaryFileHeaderType) :: headernext
  contains
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
    procedure :: build_index
    procedure :: this_header
    procedure :: next_header
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_if(this, success, iout, header_only)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
      logical(LGP), intent(in), optional :: header_only
    end subroutine read_record_if
  end interface
contains

  !> @brief Get a string representation of the header.
  function get_str(this) result(str)
    class(BinaryFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str
    ! local
    character(len=LENHUGELINE) :: temp

    write (temp, '(*(G0))') &
      'Binary file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ')'
    str = trim(temp)
  end function get_str

  !> @brief Peek to see if another record is available.
  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat

    if (this%indexed .and. allocated(this%headers)) then
      if (this%current == size(this%headers)) then
        this%endoffile = .true.
        return
      end if
    end if

    if (.not. this%endoffile) then
      read (this%inunit, iostat=iostat) this%headers(this%current + 1)%header%kstp, this%headers(this%current + 1)%header%kper
      if (iostat == 0) then
        this%current = this%current + 1
        call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
      else if (iostat < 0) then
        this%endoffile = .true.
      end if
    end if
  end subroutine peek_record

  subroutine build_index(this, iout)
    class(BinaryFileReaderType), intent(inout) :: this
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i
    logical(LGP) :: success

    if (this%indexed) return
    rewind (this%inunit)
    i = 0
    do
      call this%read_record(success, iout, header_only=.true.)
      i = i + 1
      if (this%endoffile) then
        if (success) exit
        call pstop(1, 'Error scanning record header')
      end if
    end do
    rewind (this%inunit)
    allocate (this%headers(i))
    i = 1
    do
      if (i > size(this%headers)) exit
      inquire(unit=this%inunit, pos=this%header%pos)
      call this%read_record(success, iout, header_only=.true.)
      if (.not. success) call pstop(1, 'Error reading record header')
      allocate (this%headers(i)%header, source=this%header)
      i = i + 1
    end do
    rewind (this%inunit)
    this%indexed = .true.
    this%endoffile = .false.
  end subroutine build_index

  function this_header(this) result(header)
    class(BinaryFileReaderType), intent(in) :: this
    type(BinaryFileHeaderType) :: header

    if (.not. this%indexed) call pstop(1, 'Index not built')
    if (this%current > size(this%headers)) call pstop(1, 'Current index out of bounds')
    header = this%headers(this%current)%header
  end function this_header

  function next_header(this) result(header)
    class(BinaryFileReaderType), intent(inout) :: this
    type(BinaryFileHeaderType) :: header

    if (.not. this%indexed) call pstop(1, 'Index not built')
    if (this%current >= size(this%headers)) call pstop(1, 'Current index out of bounds')
    this%current = this%current + 1
    if (this%current > size(this%headers)) then
      this%endoffile = .true.
      return
    end if
    header = this%headers(this%current)%header
  end function next_header

end module BinaryFileReaderModule
