module HeadFileReaderModule

  use KindModule
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: LINELENGTH
  use BinaryFileHeaderModule, only: BinaryFileHeaderType
  use BinaryFileReaderModule, only: BinaryFileReaderType
  use InputOutputModule, only: fseek_stream

  implicit none

  private
  public :: HeadFileReaderType, HeadFileHeaderType

  type, extends(BinaryFileHeaderType) :: HeadFileHeaderType
    character(len=LINELENGTH) :: srcmodelname
    character(len=LINELENGTH) :: srcpackagename
    character(len=LINELENGTH) :: dstmodelname
    character(len=LINELENGTH) :: dstpackagename
  contains
    procedure :: get_str
  end type HeadFileHeaderType

  type, extends(BinaryFileReaderType) :: HeadFileReaderType
    integer(I4B) :: nlay
    real(DP), dimension(:), allocatable :: head
  contains
    procedure :: initialize
    procedure :: build_index
    procedure :: read_record
    procedure :: finalize
    procedure, private :: detect_nlay
  end type HeadFileReaderType

contains

  function get_str(this) result(str)
    class(HeadFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str

    write (str, '(*(G0))') &
      'Head file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ', text: ', trim(this%text), &
      ', srcmodelname: ', trim(this%srcmodelname), &
      ', srcpackagename: ', trim(this%srcpackagename), &
      ', dstmodelname: ', trim(this%dstmodelname), &
      ', dstpackagename: ', trim(this%dstpackagename), ')'
    str = trim(str)
  end function get_str

  subroutine detect_nlay(this, iout)
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: kstp_last, kper_last
    logical :: success
    !
    this%inunit = iu
    this%endoffile = .false.
    this%nlay = 0
    !
    ! -- Read the first head data record to set kstp_last, kstp_last
    call this%read_record(success)
    kstp_last = this%header%kstp
    kper_last = this%header%kper
    rewind (this%inunit)
    !
    ! -- Determine number of records within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading binary file to determine number of records per time step.'
    do
      call this%read_record(success, iout)
      if (.not. success) exit
      if (kstp_last /= this%header%kstp .or. kper_last /= this%header%kper) exit
      this%nlay = this%nlay + 1
    end do

    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nlay, &
      ' unique records in binary file.'
  end subroutine detect_nlay

  subroutine initialize(this, iu, iout)
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout

    this%inunit = iu
    this%endoffile = .false.
    this%nlay = 0
    call this%build_index(iout)
    call this%detect_nlay(iout)
  end subroutine initialize

  subroutine build_index(this, iout)
    class(HeadFileReaderType), intent(inout) :: this
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i
    logical(LGP) :: success
    type(HeadFileHeaderType) :: header

    if (this%indexed) return
    i = 1
    rewind (this%inunit)
    do
      call this%read_record(header, success, iout)
      if (this%endoffile) then
        if (success) exit
        call pstop(1, 'Error scanning record header')
      end if
      i = i + 1
    end do
    this%total = i
    allocate (this%headers(this%total))
    rewind (this%inunit)
    i = 1
    do
      if (i > this%total) exit
      call this%read_record(header, success, iout)
      if (.not. success) call pstop(1, 'Error reading record header')
      allocate (this%headers(i)%header, source=header)
      i = i + 1
    end do
    rewind (this%inunit)
    this%current = 1
    this%indexed = .true.
  end subroutine build_index

  subroutine read_record(this, header, success, iout)
    ! dummy
    class(HeadFileReaderType), intent(inout) :: this
    class(BinaryFileHeaderType), intent(out) :: header
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: iostat, iout_opt
    integer(I4B) :: ncol, nrow, ilay

    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    !
    this%header%kstp = 0
    this%header%kper = 0
    success = .true.
    this%headernext%kstp = 0
    this%headernext%kper = 0
    read (this%inunit, iostat=iostat) this%header%kstp, this%header%kper, &
      this%header%pertim, this%header%totim, this%text, ncol, nrow, ilay
    if (iostat /= 0) then
      success = .false.
      if (iostat < 0) this%endoffile = .true.
      return
    end if
    !
    ! -- allocate head to proper size
    if (.not. allocated(this%head)) then
      allocate (this%head(ncol * nrow))
    else
      if (size(this%head) /= ncol * nrow) then
        deallocate (this%head)
        allocate (this%head(ncol * nrow))
      end if
    end if

    ! read head array
    read (this%inunit) this%head

    ! set end of file flag
    call this%peek_record()
  end subroutine read_record

  subroutine finalize(this)
    class(HeadFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%head)) deallocate (this%head)
  end subroutine finalize

end module HeadFileReaderModule
