module BudgetFileReaderModule

  use KindModule
  use ErrorUtilModule, only: pstop
  use SimModule, only: store_error, store_error_unit
  use ConstantsModule, only: LINELENGTH
  use BinaryFileReaderModule, only: BinaryFileReaderType, &
                                    BinaryFileHeaderType, &
                                    BinaryFileIndexType

  implicit none

  private
  public :: BudgetFileReaderType

  type, extends(BinaryFileHeaderType) :: BudgetFileHeaderType
    character(len=16) :: budtxt
    integer(I4B) :: nval = 0
    integer(I4B) :: naux = 0
    integer(I4B) :: idum1 = 0
    integer(I4B) :: idum2 = 0
    integer(I4B) :: imeth = 0
  contains
    procedure :: reset
  end type BudgetFileHeaderType

  type, extends(BinaryFileIndexType) :: BudgetFileIndexType
    type(BudgetFileHeaderType) :: header
  contains
    procedure :: read_header
  end type BudgetFileIndexType

  type, extends(BinaryFileReaderType) :: BudgetFileReaderType
    type(BudgetFileIndexType) :: index
    logical :: hasimeth1flowja = .false.
    integer(I4B) :: nbudterms
    character(len=16), dimension(:), allocatable :: budtxtarray
    integer(I4B), dimension(:), allocatable :: imetharray
    character(len=16) :: srcmodelname
    character(len=16) :: srcpackagename
    integer(I4B) :: ndat
    integer(I4B), dimension(:), allocatable :: nauxarray
    character(len=16), dimension(:), allocatable :: auxtxt
    character(len=16), dimension(:, :), allocatable :: auxtxtarray
    integer(I4B) :: nlist
    real(DP), dimension(:), allocatable :: flowja
    integer(I4B), dimension(:), allocatable :: nodesrc
    integer(I4B), dimension(:), allocatable :: nodedst
    real(DP), dimension(:), allocatable :: flow
    real(DP), dimension(:, :), allocatable :: auxvar
    character(len=16) :: dstmodelname
    character(len=16) :: dstpackagename
    character(len=16), dimension(:), allocatable :: dstpackagenamearray
  contains
    procedure :: initialize
    procedure :: read_record
    procedure :: finalize
  end type BudgetFileReaderType

contains

  subroutine initialize(this, iu, iout, ncrbud)
    ! dummy
    class(BudgetFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(out) :: ncrbud
    ! local
    integer(I4B) :: ibudterm
    integer(I4B) :: kstp_last, kper_last
    integer(I4B) :: maxaux
    logical :: success
    
    this%inunit = iu
    this%index%inunit = iu
    this%endoffile = .false.
    this%nbudterms = 0
    ncrbud = 0
    maxaux = 0
    
    ! Determine number of budget terms within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading budget file to determine number of terms per time step.'
    
    ! Read through the first set of data for time step 1 and stress period 1
    do
      call this%read_record(success)
      if (.not. success) exit
      this%nbudterms = this%nbudterms + 1
      if (this%index%header%naux > maxaux) maxaux = &
        this%index%header%naux
      if (this%index%header%kstp /= this%headernext%kstp .or. &
          this%index%header%kper /= this%headernext%kper) &
        exit
    end do
    kstp_last = this%index%header%kstp
    kper_last = this%index%header%kper
    allocate (this%budtxtarray(this%nbudterms))
    allocate (this%imetharray(this%nbudterms))
    allocate (this%dstpackagenamearray(this%nbudterms))
    allocate (this%nauxarray(this%nbudterms))
    allocate (this%auxtxtarray(maxaux, this%nbudterms))
    this%auxtxtarray(:, :) = ''
    rewind (this%inunit)
    call this%index%reset() ! temporary
    
    ! Now read through again and store budget text names
    do ibudterm = 1, this%nbudterms
      call this%read_record(success, iout)
      if (.not. success) exit
      this%budtxtarray(ibudterm) = this%index%header%budtxt
      this%imetharray(ibudterm) = this%index%header%imeth
      this%dstpackagenamearray(ibudterm) = this%dstpackagename
      this%nauxarray(ibudterm) = this%index%header%naux
      if (this%index%header%naux > 0) then
        this%auxtxtarray(1:this%index%header%naux, ibudterm) = this%auxtxt(:)
      end if
      if (this%srcmodelname == this%dstmodelname) then
        if (allocated(this%nodesrc)) ncrbud = max(ncrbud, maxval(this%nodesrc))
      end if
    end do
    rewind (this%inunit)
    call this%index%reset() ! temporary
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nbudterms, &
      ' unique flow terms in budget file.'
  end subroutine initialize

  subroutine read_header(this, eof)
    ! dummy
    class(BudgetFileIndexType), intent(inout) :: this
    logical(LGP), intent(out) :: eof
    ! local
    integer(I4B) :: iostat
    
    eof = .false.
    inquire (this%inunit, pos=this%header%pos)
    read (this%inunit, iostat=iostat) this%header%kstp, this%header%kper, &
      this%header%budtxt, this%header%nval, &
      this%header%idum1, this%header%idum2
    if (iostat /= 0) then
      if (iostat < 0) then
        eof = .true.
        return
      else
        call pstop(1, "error reading budget file header")
      end if
    end if
    read (this%inunit) this%header%imeth, this%header%delt, &
      this%header%pertim, this%header%totim
  end subroutine read_header

  subroutine read_record(this, success, iout)
    ! dummy
    class(BudgetFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i, n, iout_opt
    character(len=LINELENGTH) :: errmsg
    
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if

    this%srcmodelname = ''
    this%srcpackagename = ''
    this%dstmodelname = ''
    this%dstpackagename = ''

    call this%index%read_header(success)
    if (.not. success) then 
      this%endoffile = .true.
      return     
    end if

    if (this%index%header%imeth == 1) then
      if (trim(adjustl(this%index%header%budtxt)) == 'FLOW-JA-FACE') then
        if (allocated(this%flowja)) deallocate (this%flowja)
        allocate (this%flowja(this%index%header%nval))
        read (this%inunit) this%flowja
        this%hasimeth1flowja = .true.
      else
        this%index%header%nval = &
          this%index%header%nval * &
          this%index%header%idum1 * &
          abs(this%index%header%idum2)
        if (allocated(this%flow)) deallocate (this%flow)
        allocate (this%flow(this%index%header%nval))
        if (allocated(this%nodesrc)) deallocate (this%nodesrc)
        allocate (this%nodesrc(this%index%header%nval))
        read (this%inunit) this%flow
        do i = 1, this%index%header%nval
          this%nodesrc(i) = i
        end do
      end if
    elseif (this%index%header%imeth == 6) then
      ! method code 6
      read (this%inunit) this%srcmodelname
      read (this%inunit) this%srcpackagename
      read (this%inunit) this%dstmodelname
      read (this%inunit) this%dstpackagename
      read (this%inunit) this%ndat
      this%index%header%naux = this%ndat - 1
      if (allocated(this%auxtxt)) deallocate (this%auxtxt)
      allocate (this%auxtxt(this%index%header%naux))
      read (this%inunit) this%auxtxt
      read (this%inunit) this%nlist
      if (allocated(this%nodesrc)) deallocate (this%nodesrc)
      allocate (this%nodesrc(this%nlist))
      if (allocated(this%nodedst)) deallocate (this%nodedst)
      allocate (this%nodedst(this%nlist))
      if (allocated(this%flow)) deallocate (this%flow)
      allocate (this%flow(this%nlist))
      if (allocated(this%auxvar)) deallocate (this%auxvar)
      allocate (this%auxvar(this%index%header%naux, this%nlist))
      read (this%inunit) (this%nodesrc(n), this%nodedst(n), this%flow(n), &
                          (this%auxvar(i, n), i=1, this%index%header%naux), n=1, this%nlist)
    else
      write (errmsg, '(a, a)') 'ERROR READING: ', &
        trim(this%index%header%budtxt)
      call store_error(errmsg)
      write (errmsg, '(a, i0)') 'INVALID METHOD CODE DETECTED: ', &
        this%index%header%imeth
      call store_error(errmsg)
      call store_error_unit(this%inunit)
    end if
    if (iout_opt > 0) then
      write (iout_opt, '(1pg15.6, a, 1x, a)') &
        this%index%header%totim, this%index%header%budtxt, &
        this%dstpackagename
    end if

    call this%peek_record()
  end subroutine read_record

  subroutine finalize(this)
    class(BudgetFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%auxtxt)) deallocate (this%auxtxt)
    if (allocated(this%flowja)) deallocate (this%flowja)
    if (allocated(this%nodesrc)) deallocate (this%nodesrc)
    if (allocated(this%nodedst)) deallocate (this%nodedst)
    if (allocated(this%flow)) deallocate (this%flow)
    if (allocated(this%auxvar)) deallocate (this%auxvar)
  end subroutine finalize

  subroutine reset(this)
    class(BudgetFileHeaderType) :: this
    call this%BinaryFileHeaderType%reset()
    this%budtxt = ''
    this%nval = 0
    this%naux = 0
    this%idum1 = 0
    this%idum2 = 0
  end subroutine reset

end module BudgetFileReaderModule
