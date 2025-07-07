module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream
  use BinaryFileHeaderModule, only: BinaryFileHeaderType

  public :: BinaryFileReaderType
  public :: read_record_if

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    integer(I4B) :: current = 1
    integer(I4B) :: total
    logical(LGP) :: built = .false.
    logical(LGP) :: endoffile = .false.
    type(BinaryFileHeaderType), allocatable :: entries(:)
  contains
    procedure :: build_index
    procedure :: destroy
    procedure :: peek_header
    procedure :: next_header
    procedure(read_record_if), deferred :: read_record
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_if(this, header, success, iout)
      import BinaryFileHeaderType
      import BinaryFileReaderType
      import I4B
      import LGP
      class(BinaryFileReaderType), intent(inout) :: this
      type(BinaryFileHeaderType), pointer, intent(inout) :: header
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
    end subroutine read_record_if
  end interface
contains

  subroutine build_index(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    type(BinaryFileHeaderType), pointer :: header
    logical(LGP) :: success
    integer(I4B) :: i, iostat, pos

    i = 0
    allocate (header)
    rewind (this%inunit)
    do
      call this%read_record(header, success)
      if (.not. success) call pstop(1, 'Error reading record')
      i = i + 1
    end do
    this%total = i
    allocate (this%entries(this%total))
    rewind (this%inunit)
    i = 1
    do
      if (i > this%total) exit
      inquire (this%inunit, pos=pos, iostat=iostat)
      call this%read_record(header, success)
      if (.not. success) call pstop(1, 'Error reading record')
      this%entries(i) = header
      this%entries(i)%pos = pos
      i = i + 1
    end do
    this%built = .true.
    this%current = 1
    deallocate (header)
  end subroutine build_index

  function peek_header(this) result(header)
    class(BinaryFileReaderType), intent(inout) :: this
    type(BinaryFileHeaderType) :: header

    header = this%entries(this%current)
  end function peek_header

  function next_header(this) result(header)
    class(BinaryFileReaderType), intent(inout) :: this
    type(BinaryFileHeaderType) :: header

    this%current = this%current + 1
    header = this%entries(this%current)
  end function next_header

  subroutine destroy(this)
    class(BinaryFileReaderType), intent(inout) :: this
    if (allocated(this%entries)) deallocate(this%entries)
  end subroutine destroy

end module BinaryFileReaderModule
