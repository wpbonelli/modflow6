module HeadFileReaderModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH
  use BinaryFileHeaderModule, only: BinaryFileHeaderType
  use BinaryFileReaderModule, only: BinaryFileReaderType

  implicit none

  private
  public :: HeadFileReaderType

  type, extends(BinaryFileReaderType) :: HeadFileReaderType
    character(len=16) :: text
    integer(I4B) :: nlay
    real(DP), dimension(:), allocatable :: head
  contains
    procedure :: initialize
    procedure :: finalize
  end type HeadFileReaderType

contains

  !< @brief initialize
  !<
  subroutine initialize(this, iu, iout)
    ! -- dummy
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: kstp_last, kper_last
    logical :: success
    type(BinaryFileHeaderType), pointer :: header
    !
    this%inunit = iu
    this%nlay = 0
    allocate(header)
    !
    ! -- Read the first head data record to set kstp_last, kstp_last
    call this%read_record(header, success)
    kstp_last = header%kstp
    kper_last = header%kper
    rewind (this%inunit)
    !
    ! -- Determine number of records within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading binary file to determine number of records per time step.'
    do
      call this%read_record(header, success, iout)
      if (.not. success) exit
      if (kstp_last /= header%kstp .or. kper_last /= header%kper) exit
      this%nlay = this%nlay + 1
    end do
    rewind (this%inunit)
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nlay, &
      ' unique records in binary file.'
  end subroutine initialize

  !< @brief read record
  !<
  subroutine read_record(this, header, success, iout)
    ! -- modules
    use InputOutputModule, only: fseek_stream
    ! -- dummy
    class(HeadFileReaderType), intent(inout) :: this
    type(BinaryFileHeaderType), pointer, intent(inout) :: header
    logical(LGP), intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! -- local
    integer(I4B) :: iostat, iout_opt
    integer(I4B) :: ncol, nrow, ilay
    !
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    !
    header%kstp = 0
    header%kper = 0
    success = .true.
    read (this%inunit, iostat=iostat) header%kstp, header%kper, header%pertim, &
      header%totim, this%text, ncol, nrow, ilay
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
    !
    ! -- read the head array
    read (this%inunit) this%head
    !
    if (this%current == this%total) this%endoffile = .true.
    this%current = this%current + 1
  end subroutine read_record

  !< @brief finalize
  !<
  subroutine finalize(this)
    class(HeadFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%head)) deallocate (this%head)
  end subroutine finalize

end module HeadFileReaderModule
