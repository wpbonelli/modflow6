module HeadFileReaderModule

  use KindModule
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: LINELENGTH
  use BinaryFileReaderModule, only: BinaryFileReaderType, &
                                    BinaryFileHeaderType, &
                                    BinaryFileIndexType

  implicit none

  private
  public :: HeadFileReaderType

  type, extends(BinaryFileHeaderType) :: HeadFileHeaderType
    character(len=16) :: text
    integer(I4B) :: ncol, nrow, ilay = 0
  contains
    procedure :: reset
  end type HeadFileHeaderType

  type, extends(BinaryFileIndexType) :: HeadFileIndexType
    type(HeadFileHeaderType) :: header
  contains
    procedure :: read_header
  end type HeadFileIndexType

  type, extends(BinaryFileReaderType) :: HeadFileReaderType
    type(HeadFileIndexType) :: index
    integer(I4B) :: nlay
    real(DP), dimension(:), allocatable :: head
  contains
    procedure :: initialize
    procedure :: read_record
    procedure :: finalize
  end type HeadFileReaderType

contains

  subroutine initialize(this, iu, iout)
    ! dummy
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    ! local
    integer(I4B) :: kstp_last, kper_last
    logical :: success
    
    this%inunit = iu
    this%index%inunit = iu
    this%endoffile = .false.
    this%nlay = 0
    
    ! Read the first head data record to set kstp_last, kstp_last
    call this%read_record(success)
    kstp_last = this%index%header%kstp
    kper_last = this%index%header%kper
    rewind (this%inunit)
    call this%index%reset() ! temporary
    
    ! Determine number of records within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading binary file to determine number of records per time step.'
    do
      call this%read_record(success, iout)
      if (.not. success) exit
      if (kstp_last /= this%index%header%kstp .or. kper_last /= this%index%header%kper) exit
      this%nlay = this%nlay + 1
    end do
    rewind (this%inunit)
    call this%index%reset() ! temporary
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nlay, &
      ' unique records in binary file.'
  end subroutine initialize

  subroutine read_header(this, eof)
    ! dummy
    class(HeadFileIndexType), intent(inout) :: this
    logical(LGP), intent(out) :: eof
    ! local
    integer(I4B) :: iostat
    
    eof = .false.
    inquire (this%inunit, pos=this%header%pos)
    read (this%inunit, iostat=iostat) this%header%kstp, this%header%kper, &
      this%header%pertim, this%header%totim, this%header%text, &
      this%header%ncol, this%header%nrow, this%header%ilay
    if (iostat /= 0) then
      if (iostat < 0) then
        eof = .true.
        return
      else
        print *, &
          'Error reading head file header ', this%current, ' at position ', this%header%pos
        call pstop(1, "error reading head file header")
      end if
    end if
    this%current = this%current + 1
  end subroutine read_header

  subroutine read_record(this, success, iout)
    ! dummy
    class(HeadFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: iout_opt
    
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    
    call this%index%read_header(success)
    if (.not. success) then 
      this%endoffile = .true.
      return     
    end if

    ! allocate head to proper size
    if (.not. allocated(this%head)) then
      allocate (this%head(this%index%header%ncol * &
                          this%index%header%nrow))
    else
      if (size(this%head) /= this%index%header%ncol * &
                             this%index%header%nrow) then
        deallocate (this%head)
        allocate (this%head(this%index%header%ncol * &
                            this%index%header%nrow))
      end if
    end if
    
    ! read the head array
    read (this%inunit) this%head
    
    call this%peek_record()
  end subroutine read_record

  subroutine finalize(this)
    class(HeadFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%head)) deallocate (this%head)
  end subroutine finalize

  subroutine reset(this)
    class(HeadFileHeaderType) :: this
    call this%BinaryFileHeaderType%reset()
    this%text = ''
    this%ncol = 0
    this%nrow = 0
    this%ilay = 0
  end subroutine reset

end module HeadFileReaderModule
