module BudgetFileReaderModule

  use KindModule
  use ErrorUtilModule, only: pstop
  use SimModule, only: store_error, store_error_unit
  use ConstantsModule, only: LINELENGTH, LENBIGLINE
  use BinaryFileHeaderModule, only: BinaryFileHeaderType
  use BinaryFileReaderModule, only: BinaryFileReaderType
  use InputOutputModule, only: fseek_stream

  implicit none

  private
  public :: BudgetFileHeaderType, BudgetFileReaderType

  type, extends(BinaryFileHeaderType) :: BudgetFileHeaderType
    integer(I4B) :: naux
    integer(I4B) :: ndat
    integer(I4B) :: nval
    integer(I4B) :: nlist
    integer(I4B) :: idum1
    integer(I4B) :: idum2
    integer(I4B) :: imeth
    character(len=16) :: srcmodelname
    character(len=16) :: srcpackagename
    character(len=16) :: dstmodelname
    character(len=16) :: dstpackagename
    character(len=16), allocatable :: auxtxt(:)
  contains
    procedure :: get_str
  end type BudgetFileHeaderType

  type, extends(BinaryFileReaderType) :: BudgetFileReaderType
    logical :: hasimeth1flowja = .false.
    integer(I4B) :: nbudterms
    character(len=16), allocatable :: budtxtarray(:)
    integer(I4B), allocatable :: imetharray(:)
    integer(I4B), allocatable :: nauxarray(:)
    character(len=16), allocatable :: auxtxtarray(:, :)
    real(DP), allocatable :: flowja(:)
    integer(I4B), allocatable :: nodesrc(:)
    integer(I4B), allocatable :: nodedst(:)
    real(DP), allocatable :: flow(:)
    real(DP), allocatable :: auxvar(:, :)
    character(len=16) :: dstmodelname
    character(len=16) :: dstpackagename
    character(len=16), dimension(:), allocatable :: dstpackagenamearray
  contains
    procedure :: initialize
    procedure :: build_index
    procedure :: read_record
    procedure :: finalize
  end type BudgetFileReaderType

contains

  function get_str(this) result(str)
    class(BudgetFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str
    ! local
    character(len=LENBIGLINE) :: temp

    write (temp, '(*(G0))') &
      'Budget file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ', text: ', trim(this%text), &
      ', naux: ', this%naux, &
      ', ndat: ', this%ndat, &
      ', nval: ', this%nval, &
      ', nlist: ', this%nlist, &
      ', idum1: ', this%idum1, &
      ', idum2: ', this%idum2, &
      ', imeth: ', this%imeth, &
      ', srcmodelname: ', trim(this%srcmodelname), &
      ', srcpackagename: ', trim(this%srcpackagename), &
      ', dstmodelname: ', trim(this%dstmodelname), &
      ', dstpackagename: ', trim(this%dstpackagename), &
      ', auxtxt: ', this%auxtxt
    str = trim(temp)
  end function get_str

  subroutine initialize(this, iu, iout, ncrbud)
    ! -- dummy
    class(BudgetFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(out) :: ncrbud
    ! -- local
    integer(I4B) :: i, kstp_last, kper_last
    integer(I4B) :: maxaux

    this%inunit = iu
    this%nbudterms = 0
    ncrbud = 0
    maxaux = 0

    call this%build_index(iout=iout)

    if (iout > 0) &
      write (iout, '(a)') &
      'Reading budget file to determine number of terms per time step.'

    ! Look through the first set of headers for time step 1 and stress period 1
    do i = 1, this%total
      this%nbudterms = this%nbudterms + 1
      if (this%naux > maxaux) maxaux = this%naux
      if (this%header%kstp /= this%headernext%kstp .or. &
          this%header%kper /= this%headernext%kper) &
        exit
    end do
    kstp_last = this%header%kstp
    kper_last = this%header%kper
    allocate (this%budtxtarray(this%nbudterms))
    allocate (this%imetharray(this%nbudterms))
    allocate (this%dstpackagenamearray(this%nbudterms))
    allocate (this%nauxarray(this%nbudterms))
    allocate (this%auxtxtarray(maxaux, this%nbudterms))
    this%auxtxtarray(:, :) = ''

    ! Now look through again and store budget text names
    do i = 1, this%nbudterms
      select type (header => this%headers(i)%header)
      type is (BudgetFileHeaderType)
        this%budtxtarray(i) = header%text
        this%imetharray(i) = header%imeth
        this%dstpackagenamearray(i) = header%dstpackagename
        this%nauxarray(i) = header%naux
        if (header%naux > 0) then
          this%auxtxtarray(1:header%naux, i) = header%auxtxt(:)
        end if
        if (header%srcmodelname == header%dstmodelname) then
          if (allocated(this%nodesrc)) ncrbud = max(ncrbud, maxval(this%nodesrc))
        end if
      end select
    end do

    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nbudterms, &
      ' unique flow terms in budget file.'
  end subroutine initialize

  subroutine build_index(this, iout)
    class(BudgetFileReaderType), intent(inout) :: this
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i
    logical(LGP) :: success
    type(BudgetFileHeaderType) :: header

    if (this%indexed) return
    i = 0
    rewind (this%inunit)
    do
      call this%read_record(header, success, iout)
      i = i + 1
      if (this%endoffile) then
        if (success) exit
        call pstop(1, 'Error scanning record header')
      end if
    end do
    this%total = i
    allocate (this%headers(this%total))
    rewind (this%inunit)
    i = 1
    do
      if (i > this%total) exit
      call this%read_record(header, success, iout)
      print *, header%get_str()
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
    class(BudgetFileReaderType), intent(inout) :: this
    class(BinaryFileHeaderType), intent(out) :: header
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i, n, iostat, iout_opt
    character(len=LINELENGTH) :: errmsg

    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    !
    this%header%kstp = 0
    this%header%kper = 0
    this%budtxt = ''
    this%nval = 0
    this%naux = 0
    this%idum1 = 0
    this%idum2 = 0
    this%srcmodelname = ''
    this%srcpackagename = ''
    this%dstmodelname = ''
    this%dstpackagename = ''

    success = .true.
    this%headernext%kstp = 0
    this%headernext%kper = 0
    read (this%inunit, iostat=iostat) this%header%kstp, this%header%kper, &
      this%budtxt, this%nval, this%idum1, this%idum2
    if (iostat /= 0) then
      success = .false.
      if (iostat < 0) this%endoffile = .true.
      return
    end if
    read (this%inunit) this%imeth, this%header%delt, &
      this%header%pertim, this%header%totim
    if (this%imeth == 1) then
      if (trim(adjustl(this%budtxt)) == 'FLOW-JA-FACE') then
        if (allocated(this%flowja)) deallocate (this%flowja)
        allocate (this%flowja(this%nval))
        read (this%inunit) this%flowja
        this%hasimeth1flowja = .true.
      else
        this%nval = this%nval * this%idum1 * abs(this%idum2)
        if (allocated(this%flow)) deallocate (this%flow)
        allocate (this%flow(this%nval))
        if (allocated(this%nodesrc)) deallocate (this%nodesrc)
        allocate (this%nodesrc(this%nval))
        read (this%inunit) this%flow
        do i = 1, this%nval
          this%nodesrc(i) = i
        end do
      end if
    elseif (this%imeth == 6) then
      ! -- method code 6
      read (this%inunit) this%srcmodelname
      read (this%inunit) this%srcpackagename
      read (this%inunit) this%dstmodelname
      read (this%inunit) this%dstpackagename
      read (this%inunit) this%ndat
      this%naux = this%ndat - 1
      if (allocated(this%auxtxt)) deallocate (this%auxtxt)
      allocate (this%auxtxt(this%naux))
      read (this%inunit) this%auxtxt
      read (this%inunit) this%nlist
      if (allocated(this%nodesrc)) deallocate (this%nodesrc)
      allocate (this%nodesrc(this%nlist))
      if (allocated(this%nodedst)) deallocate (this%nodedst)
      allocate (this%nodedst(this%nlist))
      if (allocated(this%flow)) deallocate (this%flow)
      allocate (this%flow(this%nlist))
      if (allocated(this%auxvar)) deallocate (this%auxvar)
      allocate (this%auxvar(this%naux, this%nlist))
      read (this%inunit) (this%nodesrc(n), this%nodedst(n), this%flow(n), &
                          (this%auxvar(i, n), i=1, this%naux), n=1, this%nlist)
    else
      write (errmsg, '(a, a)') 'ERROR READING: ', trim(this%budtxt)
      call store_error(errmsg)
      write (errmsg, '(a, i0)') 'INVALID METHOD CODE DETECTED: ', this%imeth
      call store_error(errmsg)
      call store_error_unit(this%inunit)
    end if
    if (iout_opt > 0) then
      write (iout_opt, '(1pg15.6, a, 1x, a)') this%header%totim, this%budtxt, &
        this%dstpackagename
    end if
    !
    call this%peek_record()
  end subroutine read_record

  subroutine finalize(this)
    class(BudgetFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%flowja)) deallocate (this%flowja)
    if (allocated(this%nodesrc)) deallocate (this%nodesrc)
    if (allocated(this%nodedst)) deallocate (this%nodedst)
    if (allocated(this%flow)) deallocate (this%flow)
    if (allocated(this%auxvar)) deallocate (this%auxvar)
  end subroutine finalize

end module BudgetFileReaderModule
