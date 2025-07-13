module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream
  use ConstantsModule, only: LENHUGELINE

  public :: BinaryFileHeaderType, BinaryFileReaderType

  type, abstract :: BinaryFileHeaderType
    integer(I4B) :: pos
    integer(I4B) :: kper, kstp
    real(DP) :: pertim, totim
  contains
    procedure :: get_str
    procedure(get_data_size_if), deferred :: get_data_size
    procedure(get_header_size_if), deferred :: get_header_size
  end type BinaryFileHeaderType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    integer(I4B) :: nrecords
    integer(I4B) :: current
    logical(LGP) :: indexed
    logical(LGP) :: endoffile
    integer(I4B), dimension(:), allocatable :: record_sizes
    class(BinaryFileHeaderType), allocatable :: header
    class(BinaryFileHeaderType), allocatable :: headernext

  contains
    procedure(read_header_if), deferred :: read_header
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
    procedure :: build_index
    procedure :: read_next
    procedure :: read_previous
  end type BinaryFileReaderType

  abstract interface
    function get_header_size_if(this) result(header_size)
      import BinaryFileHeaderType, I4B
      class(BinaryFileHeaderType), intent(in) :: this
      integer(I4B) :: header_size
    end function get_header_size_if
    function get_data_size_if(this) result(data_size)
      import BinaryFileHeaderType, I4B
      class(BinaryFileHeaderType), intent(in) :: this
      integer(I4B) :: data_size
    end function get_data_size_if
    subroutine read_header_if(this, success, iout, seek_next)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
      logical(LGP), intent(in), optional :: seek_next
    end subroutine read_header_if
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
        if (allocated(this%headernext)) deallocate (this%headernext)
      end if
    end if
  end subroutine peek_record

  !> @brief Build an index of data sizes in the file
  subroutine build_index(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    logical :: success

    this%indexed = .false.
    this%nrecords = 0

    ! first pass: count
    rewind (this%inunit)
    do
      ! call this%read_header(success, seek_next=.true.)
      call this%read_record(success)
      if (.not. success) exit
      this%nrecords = this%nrecords + 1
    end do
    if (allocated(this%record_sizes)) deallocate (this%record_sizes)
    allocate (this%record_sizes(this%nrecords))

    ! second pass: save
    rewind (this%inunit)
    do
      ! call this%read_header(success, seek_next=.true.)
      call this%read_record(success)
      if (.not. success) exit
      print *, "Indexed: ", this%header%get_str()
      this%record_sizes(this%nrecords) = this%header%get_header_size() + &
                                         this%header%get_data_size()
    end do

    this%indexed = .true.
    print *, "Index built with ", this%nrecords, " records."

    ! Reset to beginning
    rewind (this%inunit)
    this%current = 0
    this%endoffile = .false.
  end subroutine build_index

  !> @brief Move to the next record
  subroutine read_next(this, success)
    class(BinaryFileReaderType), intent(inout) :: this
    logical, intent(out) :: success

    call this%read_record(success)
    if (success) this%current = this%current + 1
  end subroutine read_next

  !> @brief Move to the previous record (requires index)
  subroutine read_previous(this, success)
    class(BinaryFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    ! local
    integer(I4B) :: iostat, record_size

    success = .false.
    if (.not. this%indexed) &
      call pstop(1, 'Index has not been built')
    if (this%current <= 1) return
    record_size = this%record_sizes(this%current - 1)
    call fseek_stream(this%inunit, -(record_size), 1, iostat)
    if (iostat /= 0) call pstop(1, 'Failed to seek backward')
    this%current = this%current - 1
    success = .true.
    call this%read_record(success)
    if (success) call fseek_stream(this%inunit, -(record_size), 1, iostat)
  end subroutine read_previous

end module BinaryFileReaderModule
