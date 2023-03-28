module PrtOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMODELNAME, MNORMAL
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr
  use SimVariablesModule, only: errmsg, warnmsg
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryHelperModule, only: create_mem_path

  implicit none
  private
  public PrtOcType, oc_cr

  !> @ brief Output control for particle tracking models
  type, extends(OutputControlType) :: PrtOcType

    integer(I4B), pointer :: itrkout => null() !< binary output file
    integer(I4B), pointer :: itrkhdr => null() !< output header file
    integer(I4B), pointer :: itrkcsv => null() !< CSV output file
    integer(I4B), pointer :: itrkevent => null() !< track event option

  contains
    procedure :: oc_ar
    procedure :: oc_da => prt_oc_da
    procedure :: allocate_scalars => prt_oc_allocate_scalars
    procedure :: read_options => prt_oc_read_options

  end type PrtOcType

contains

  !> @ brief Create an output control object
  subroutine oc_cr(ocobj, name_model, inunit, iout)
    type(PrtOcType), pointer :: ocobj !< PrtOcType object
    character(len=*), intent(in) :: name_model !< name of the model
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output

    ! -- Create the object
    allocate (ocobj)

    ! -- Allocate scalars
    call ocobj%allocate_scalars(name_model)

    ! -- Save unit numbers
    ocobj%inunit = inunit
    ocobj%iout = iout

    ! -- Initialize block parser
    call ocobj%parser%Initialize(inunit, iout)
  end subroutine oc_cr

  subroutine prt_oc_allocate_scalars(this, name_model)
    class(PrtOcType) :: this
    character(len=*), intent(in) :: name_model !< name of model

    this%memoryPath = create_mem_path(name_model, 'OC')

    allocate (this%name_model)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%iperoc, 'IPEROC', this%memoryPath)
    call mem_allocate(this%iocrep, 'IOCREP', this%memoryPath)
    call mem_allocate(this%itrkout, 'ITRKOUT', this%memoryPath)
    call mem_allocate(this%itrkhdr, 'ITRKHDR', this%memoryPath)
    call mem_allocate(this%itrkcsv, 'ITRKCSV', this%memoryPath)
    call mem_allocate(this%itrkevent, 'ITRACKEVENT', this%memoryPath)

    this%name_model = name_model
    this%inunit = 0
    this%iout = 0
    this%ibudcsv = 0
    this%iperoc = 0
    this%iocrep = 0
    this%itrkout = 0
    this%itrkhdr = 0
    this%itrkcsv = 0
    this%itrkevent = -1
  end subroutine prt_oc_allocate_scalars

  !> @ brief Setup output control variables.
  subroutine oc_ar(this, mass, dis, dnodata)
    ! -- dummy
    class(PrtOcType) :: this !< PrtOcType object
    real(DP), dimension(:), pointer, contiguous, intent(in) :: mass !< particle mass
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization package
    real(DP), intent(in) :: dnodata !< no data value
    ! -- local
    integer(I4B) :: i, nocdobj, inodata
    type(OutputControlDataType), pointer :: ocdobjptr
    real(DP), dimension(:), pointer, contiguous :: nullvec => null()

    ! -- Initialize variables
    inodata = 0
    nocdobj = 2
    allocate (this%ocdobj(nocdobj))
    do i = 1, nocdobj
      call ocd_cr(ocdobjptr)
      select case (i)
      case (1)
        call ocdobjptr%init_dbl('BUDGET', nullvec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      case (2)
        call ocdobjptr%init_dbl('MASS', mass, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      end select
      this%ocdobj(i) = ocdobjptr
      deallocate (ocdobjptr)
    end do

    ! -- Read options or set defaults if this package not on
    if (this%inunit > 0) then
      call this%read_options()
    end if
  end subroutine oc_ar

  subroutine prt_oc_da(this)
    ! -- dummy
    class(PrtOcType) :: this
    ! -- local
    integer(I4B) :: i

    do i = 1, size(this%ocdobj)
      call this%ocdobj(i)%ocd_da()
    end do
    deallocate (this%ocdobj)

    deallocate (this%name_model)
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%iperoc)
    call mem_deallocate(this%iocrep)
    call mem_deallocate(this%itrkout)
    call mem_deallocate(this%itrkhdr)
    call mem_deallocate(this%itrkcsv)
    call mem_deallocate(this%itrkevent)
  end subroutine prt_oc_da

  subroutine prt_oc_read_options(this)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile, lowcase
    use ConstantsModule, only: LINELENGTH
    use TrackModule, only: TRACKHEADERS, TRACKTYPES
    use SimModule, only: store_error, store_error_unit
    use InputOutputModule, only: openfile, getunit
    ! -- dummy
    class(PrtOcType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: keyword2
    character(len=LINELENGTH) :: fname
    character(len=:), allocatable :: line
    integer(I4B) :: ierr
    integer(I4B) :: ipos
    logical :: isfound, found, endOfBlock
    type(OutputControlDataType), pointer :: ocdobjptr
    character(len=LINELENGTH) :: trkevent
    ! -- formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)

    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(/,1x,a,/)') 'PROCESSING OC OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        found = .false.
        select case (keyword)
        case ('BUDGETCSV')
          call this%parser%GetStringCaps(keyword2)
          if (keyword2 /= 'FILEOUT') then
            errmsg = "BUDGETCSV must be followed by FILEOUT and then budget &
              &csv file name.  Found '"//trim(keyword2)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          this%ibudcsv = GetUnit()
          call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                        filstat_opt='REPLACE')
          found = .true.
        case ('TRACK')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            ! parse filename
            call this%parser%GetString(fname)
            ! open binary track output file
            this%itrkout = getunit()
            call openfile(this%itrkout, this%iout, fname, 'DATA(BINARY)', &
                          form, access, filstat_opt='REPLACE', &
                          mode_opt=MNORMAL)
            write (this%iout, fmttrkbin) trim(adjustl(fname)), this%itrkout
            ! open and write ascii track header file
            this%itrkhdr = getunit()
            fname = trim(fname)//'.hdr'
            call openfile(this%itrkhdr, this%iout, fname, 'CSV', &
                          filstat_opt='REPLACE', mode_opt=MNORMAL)
            write (this%itrkhdr, '(a,/,a)') TRACKHEADERS, TRACKTYPES
          else
            call store_error('OPTIONAL TRACK KEYWORD MUST BE '// &
                             'FOLLOWED BY FILEOUT')
          end if
          found = .true.
        case ('TRACKCSV')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            ! parse filename
            call this%parser%GetString(fname)
            ! open CSV track output file and write headers
            this%itrkcsv = getunit()
            call openfile(this%itrkcsv, this%iout, fname, 'CSV', &
                          filstat_opt='REPLACE')
            write (this%iout, fmttrkcsv) trim(adjustl(fname)), this%itrkcsv
            write (this%itrkcsv, '(a)') TRACKHEADERS
          else
            call store_error('OPTIONAL TRACKCSV KEYWORD MUST BE &
              &FOLLOWED BY FILEOUT')
          end if
          found = .true.
        case ('TRACKEVENT')
          call this%parser%GetStringCaps(trkevent)
          select case (trkevent)
          case ('')
            this%itrkevent = -1
          case ('ALL')
            this%itrkevent = -1
          case ('RELEASE')
            this%itrkevent = 0
          case ('TRANSIT')
            this%itrkevent = 1
          case ('TIMESTEP')
            this%itrkevent = 2
          case ('TERMINATE')
            this%itrkevent = 3
          case ('WEAKSINK')
            this%itrkevent = 4
          case default
            write (errmsg, '(2a)') &
              'Looking for ALL, RELEASE, TRANSIT, TIMESTEP, &
              &TERMINATE, or WEAKSINK. Found: ', &
              trim(adjustl(trkevent))
            call store_error(errmsg, terminate=.TRUE.)
          end select
          found = .true.
        case default
          found = .false.
        end select

        if (.not. found) then
          do ipos = 1, size(this%ocdobj)
            ocdobjptr => this%ocdobj(ipos)
            if (keyword == trim(ocdobjptr%cname)) then
              found = .true.
              exit
            end if
          end do
          if (.not. found) then
            errmsg = "UNKNOWN OC OPTION '"//trim(keyword)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetRemainingLine(line)
          call ocdobjptr%set_option(line, this%parser%iuactive, this%iout)
        end if
      end do
      write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    end if
  end subroutine prt_oc_read_options

end module PrtOcModule
