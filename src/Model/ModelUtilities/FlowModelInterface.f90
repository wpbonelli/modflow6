module FlowModelInterfaceModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DONE, DZERO, DHALF, LINELENGTH, LENBUDTXT, &
                             LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error, store_error_unit
  use SimVariablesModule, only: errmsg
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use ListModule, only: ListType
  use BudgetFileReaderModule, only: BudgetFileReaderType
  use HeadFileReaderModule, only: HeadFileReaderType
  use GridFileReaderModule, only: GridFileReaderType
  use PackageBudgetModule, only: PackageBudgetType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr_bfr

  implicit none
  private
  public :: FlowModelInterfaceType

  type, extends(NumericalPackageType) :: FlowModelInterfaceType

    character(len=LENPACKAGENAME) :: text = '' !< text string for package
    logical, pointer :: flows_from_file => null() !< if .false., then flows come from GWF through GWF-Model exg
    type(ListType), pointer :: gwfbndlist => null() !< list of gwf stress packages
    integer(I4B), pointer :: iflowsupdated => null() !< flows were updated for this time step
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to this model ibound
    real(DP), dimension(:), pointer, contiguous :: gwfflowja => null() !< pointer to the GWF flowja array
    real(DP), dimension(:, :), pointer, contiguous :: gwfspdis => null() !< pointer to npf specific discharge array
    real(DP), dimension(:), pointer, contiguous :: gwfhead => null() !< pointer to the GWF head array
    real(DP), dimension(:), pointer, contiguous :: gwfsat => null() !< pointer to the GWF saturation array
    integer(I4B), dimension(:), pointer, contiguous :: ibdgwfsat0 => null() !< mark cells with saturation = 0 to exclude from dispersion
    integer(I4B), pointer :: idryinactive => null() !< mark cells with an additional flag to exclude from deactivation (gwe will simulate conduction through dry cells)
    real(DP), dimension(:), pointer, contiguous :: gwfstrgss => null() !< pointer to flow model QSTOSS
    real(DP), dimension(:), pointer, contiguous :: gwfstrgsy => null() !< pointer to flow model QSTOSY
    integer(I4B), pointer :: igwfstrgss => null() !< indicates if gwfstrgss is available
    integer(I4B), pointer :: igwfstrgsy => null() !< indicates if gwfstrgsy is available
    integer(I4B), pointer :: iubud => null() !< unit number GWF budget file
    integer(I4B), pointer :: iuhds => null() !< unit number GWF head file
    integer(I4B), pointer :: iumvr => null() !< unit number GWF mover budget file
    integer(I4B), pointer :: iugrb => null() !< unit number binary grid file
    integer(I4B), pointer :: nflowpack => null() !< number of GWF flow packages
    integer(I4B), dimension(:), pointer, contiguous :: igwfmvrterm => null() !< flag to indicate that gwf package is a mover term
    type(BudgetFileReaderType) :: bfr !< budget file reader
    type(HeadFileReaderType) :: hfr !< head file reader
    type(GridFileReaderType) :: gfr !< grid file reader
    type(PackageBudgetType), dimension(:), allocatable :: gwfpackages !< used to get flows between a package and gwf
    type(BudgetObjectType), pointer :: mvrbudobj => null() !< pointer to the mover budget object
    character(len=16), dimension(:), allocatable :: flowpacknamearray !< array of boundary package names (e.g. LAK-1, SFR-3, etc.)
    character(len=LENVARNAME) :: depvartype = ''

  contains

    procedure :: advance_bfr
    procedure :: advance_hfr
    procedure :: allocate_arrays
    procedure :: allocate_gwfpackages
    procedure :: allocate_scalars
    procedure :: deallocate_gwfpackages
    procedure :: finalize_bfr
    procedure :: finalize_hfr
    procedure :: finalize_gfr
    procedure :: fmi_ar
    procedure :: fmi_da
    procedure :: fmi_df
    procedure :: get_package_index
    procedure :: initialize_bfr
    procedure :: initialize_gwfterms_from_bfr
    procedure :: initialize_gwfterms_from_gwfbndlist
    procedure :: initialize_hfr
    procedure :: initialize_gfr
    procedure :: read_options
    procedure :: read_packagedata
  end type FlowModelInterfaceType

contains

  !> @brief Define the flow model interface
  !<
  subroutine fmi_df(this, dis, idryinactive)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: idryinactive
    ! -- formats
    character(len=*), parameter :: fmtfmi = &
      "(1x,/1x,'FMI -- FLOW MODEL INTERFACE, VERSION 2, 8/17/2023',            &
      &' INPUT READ FROM UNIT ', i0, //)"
    character(len=*), parameter :: fmtfmi0 = &
                    "(1x,/1x,'FMI -- FLOW MODEL INTERFACE,'&
                    &' VERSION 2, 8/17/2023')"
    !
    ! --print a message identifying the FMI package.
    if (this%iout > 0) then
      if (this%inunit /= 0) then
        write (this%iout, fmtfmi) this%inunit
      else
        write (this%iout, fmtfmi0)
        if (this%flows_from_file) then
          write (this%iout, '(a)') '  FLOWS ARE ASSUMED TO BE ZERO.'
        else
          write (this%iout, '(a)') '  FLOWS PROVIDED BY A GWF MODEL IN THIS &
            &SIMULATION'
        end if
      end if
    end if
    !
    ! -- Store pointers
    this%dis => dis
    !
    ! -- Read fmi options
    if (this%inunit /= 0) then
      call this%read_options()
    end if

    !
    ! -- Read packagedata options
    if (this%inunit /= 0 .and. this%flows_from_file) then
      call this%read_packagedata()
      call this%initialize_gwfterms_from_bfr()
    end if
    !
    ! -- If GWF-Model exchange is active, setup flow terms
    if (.not. this%flows_from_file) then
      call this%initialize_gwfterms_from_gwfbndlist()
    end if
    !
    ! -- Set flag that stops dry flows from being deactivated in a GWE
    !    transport model since conduction will still be simulated.
    !    0: GWE (skip deactivation step); 1: GWT (default: use existing code)
    this%idryinactive = idryinactive
  end subroutine fmi_df

  !> @brief Allocate the package
  !<
  subroutine fmi_ar(this, ibound)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    !
    ! -- store pointers to arguments that were passed in
    this%ibound => ibound
    !
    ! -- Allocate arrays
    call this%allocate_arrays(this%dis%nodes)
  end subroutine fmi_ar

  !> @brief Deallocate variables
  !<
  subroutine fmi_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(FlowModelInterfaceType) :: this

    ! -- todo: finalize hfr and bfr either here or in a finalize routine
    !
    ! -- deallocate any memory stored with gwfpackages
    call this%deallocate_gwfpackages()
    !
    ! -- deallocate fmi arrays
    deallocate (this%gwfpackages)
    deallocate (this%flowpacknamearray)
    call mem_deallocate(this%igwfmvrterm)
    call mem_deallocate(this%ibdgwfsat0)
    !
    if (this%flows_from_file) then
      call mem_deallocate(this%gwfstrgss)
      call mem_deallocate(this%gwfstrgsy)
    end if
    !
    ! -- special treatment, these could be from mem_checkin
    call mem_deallocate(this%gwfhead, 'GWFHEAD', this%memoryPath)
    call mem_deallocate(this%gwfsat, 'GWFSAT', this%memoryPath)
    call mem_deallocate(this%gwfspdis, 'GWFSPDIS', this%memoryPath)
    call mem_deallocate(this%gwfflowja, 'GWFFLOWJA', this%memoryPath)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%flows_from_file)
    call mem_deallocate(this%iflowsupdated)
    call mem_deallocate(this%igwfstrgss)
    call mem_deallocate(this%igwfstrgsy)
    call mem_deallocate(this%iubud)
    call mem_deallocate(this%iuhds)
    call mem_deallocate(this%iumvr)
    call mem_deallocate(this%iugrb)
    call mem_deallocate(this%nflowpack)
    call mem_deallocate(this%idryinactive)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine fmi_da

  !> @brief Allocate scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%flows_from_file, 'FLOWS_FROM_FILE', this%memoryPath)
    call mem_allocate(this%iflowsupdated, 'IFLOWSUPDATED', this%memoryPath)
    call mem_allocate(this%igwfstrgss, 'IGWFSTRGSS', this%memoryPath)
    call mem_allocate(this%igwfstrgsy, 'IGWFSTRGSY', this%memoryPath)
    call mem_allocate(this%iubud, 'IUBUD', this%memoryPath)
    call mem_allocate(this%iuhds, 'IUHDS', this%memoryPath)
    call mem_allocate(this%iumvr, 'IUMVR', this%memoryPath)
    call mem_allocate(this%iugrb, 'IUGRB', this%memoryPath)
    call mem_allocate(this%nflowpack, 'NFLOWPACK', this%memoryPath)
    call mem_allocate(this%idryinactive, "IDRYINACTIVE", this%memoryPath)
    !
    ! !
    ! -- Initialize
    this%flows_from_file = .true.
    this%iflowsupdated = 1
    this%igwfstrgss = 0
    this%igwfstrgsy = 0
    this%iubud = 0
    this%iuhds = 0
    this%iumvr = 0
    this%iugrb = 0
    this%nflowpack = 0
    this%idryinactive = 1
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine allocate_arrays(this, nodes)
    use MemoryManagerModule, only: mem_allocate
    !modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    integer(I4B) :: n
    !
    ! -- Allocate ibdgwfsat0, which is an indicator array marking cells with
    !    saturation greater than 0.0 with a value of 1
    call mem_allocate(this%ibdgwfsat0, nodes, 'IBDGWFSAT0', this%memoryPath)
    do n = 1, nodes
      this%ibdgwfsat0(n) = 1
    end do
    !
    ! -- Allocate differently depending on whether or not flows are
    !    being read from a file.
    if (this%flows_from_file) then
      call mem_allocate(this%gwfflowja, this%dis%con%nja, &
                        'GWFFLOWJA', this%memoryPath)
      call mem_allocate(this%gwfsat, nodes, 'GWFSAT', this%memoryPath)
      call mem_allocate(this%gwfhead, nodes, 'GWFHEAD', this%memoryPath)
      call mem_allocate(this%gwfspdis, 3, nodes, 'GWFSPDIS', this%memoryPath)
      do n = 1, nodes
        this%gwfsat(n) = DONE
        this%gwfhead(n) = DZERO
        this%gwfspdis(:, n) = DZERO
      end do
      do n = 1, size(this%gwfflowja)
        this%gwfflowja(n) = DZERO
      end do
      !
      ! -- allocate and initialize storage arrays
      if (this%igwfstrgss == 0) then
        call mem_allocate(this%gwfstrgss, 1, 'GWFSTRGSS', this%memoryPath)
      else
        call mem_allocate(this%gwfstrgss, nodes, 'GWFSTRGSS', this%memoryPath)
      end if
      if (this%igwfstrgsy == 0) then
        call mem_allocate(this%gwfstrgsy, 1, 'GWFSTRGSY', this%memoryPath)
      else
        call mem_allocate(this%gwfstrgsy, nodes, 'GWFSTRGSY', this%memoryPath)
      end if
      do n = 1, size(this%gwfstrgss)
        this%gwfstrgss(n) = DZERO
      end do
      do n = 1, size(this%gwfstrgsy)
        this%gwfstrgsy(n) = DZERO
      end do
      !
      ! -- If there is no fmi package, then there are no flows at all or a
      !    connected GWF model, so allocate gwfpackages to zero
      if (this%inunit == 0) call this%allocate_gwfpackages(this%nflowpack)
    end if
  end subroutine allocate_arrays

  !> @brief Read options from input file
  !<
  subroutine read_options(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, DEM6
    use InputOutputModule, only: getunit, openfile, urdaux
    use SimModule, only: store_error, store_error_unit
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING FMI OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SAVE_FLOWS')
          this%ipakcb = -1
        case default
          write (errmsg, '(a,3(1x,a))') &
            'UNKNOWN', trim(adjustl(this%text)), 'OPTION:', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF FMI OPTIONS'
    end if
  end subroutine read_options

  !> @brief Read packagedata block from input file
  !<
  subroutine read_packagedata(this)
    ! -- modules
    use OpenSpecModule, only: ACCESS, FORM
    use ConstantsModule, only: LINELENGTH, DEM6, LENPACKAGENAME
    use InputOutputModule, only: getunit, openfile, urdaux
    use SimModule, only: store_error, store_error_unit
    use DisModule, only: DisType
    use DisvModule, only: DisvType
    use DisuModule, only: DisuType
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword, fname
    integer(I4B) :: ierr
    integer(I4B) :: inunit
    integer(I4B) :: iapt
    logical :: isfound, endOfBlock
    logical :: blockrequired
    logical :: exist
    integer(I4B), allocatable :: idomain1d(:), idomain2d(:, :), idomain3d(:, :, :)
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and PRT Models do not have the same discretization for FMI&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  PRT Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
    character(len=*), parameter :: fmtidomerr = &
      "('GWF and PRT Models do not have the same discretization for FMI&
      & ',a,'.&
      &  GWF Model and PRT Model have different IDOMAIN arrays.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
    !
    ! -- initialize
    iapt = 0
    blockrequired = .true.
    !
    ! -- get packagedata block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              blockRequired=blockRequired, &
                              supportOpenClose=.true.)
    !
    ! -- parse packagedata block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING FMI PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('GWFBUDGET')
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEIN') then
            call store_error('GWFBUDGET KEYWORD MUST BE FOLLOWED BY '// &
                             '"FILEIN" then by filename.')
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          inunit = getunit()
          inquire (file=trim(fname), exist=exist)
          if (.not. exist) then
            call store_error('Could not find file '//trim(fname))
            call this%parser%StoreErrorUnit()
          end if
          call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                        ACCESS, 'UNKNOWN')
          this%iubud = inunit
          call this%initialize_bfr()
        case ('GWFHEAD')
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEIN') then
            call store_error('GWFHEAD KEYWORD MUST BE FOLLOWED BY '// &
                             '"FILEIN" then by filename.')
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          inquire (file=trim(fname), exist=exist)
          if (.not. exist) then
            call store_error('Could not find file '//trim(fname))
            call this%parser%StoreErrorUnit()
          end if
          inunit = getunit()
          call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                        ACCESS, 'UNKNOWN')
          this%iuhds = inunit
          call this%initialize_hfr()
        case ('GWFMOVER')
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEIN') then
            call store_error('GWFMOVER KEYWORD MUST BE FOLLOWED BY '// &
                             '"FILEIN" then by filename.')
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          inunit = getunit()
          call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                        ACCESS, 'UNKNOWN')
          this%iumvr = inunit
          call budgetobject_cr_bfr(this%mvrbudobj, 'MVT', this%iumvr, &
                                   this%iout)
          call this%mvrbudobj%fill_from_bfr(this%dis, this%iout)
        case ('GWFGRID')
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEIN') then
            call store_error('GWFGRID KEYWORD MUST BE FOLLOWED BY '// &
                             '"FILEIN" then by filename.')
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          inunit = getunit()
          call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                        ACCESS, 'UNKNOWN')
          this%iugrb = inunit
          call this%initialize_gfr()

          ! check node count
          if (this%dis%nodes /= this%gfr%get_nodes()) then
            write (errmsg, fmtdiserr) trim(this%text)
            call store_error(errmsg, terminate=.TRUE.)
          end if

          ! check idomain
          select case (this%gfr%grid_type)
          case ('DIS')
            call this%gfr%get_idomain(idomain3d)
            select type (dis => this%dis)
            type is (DisType)
              if (.not. all(dis%idomain == idomain3d)) then
                write (errmsg, fmtidomerr) trim(this%text)
                call store_error(errmsg, terminate=.TRUE.)
              end if
            end select
          case ('DISV')
            call this%gfr%get_idomain(idomain2d)
            select type (dis => this%dis)
            type is (DisvType)
              if (.not. all(dis%idomain == idomain2d)) then
                write (errmsg, fmtidomerr) trim(this%text)
                call store_error(errmsg, terminate=.TRUE.)
              end if
            end select
          case ('DISU')
            call this%gfr%get_idomain(idomain1d)
            select type (dis => this%dis)
            type is (DisuType)
              if (.not. all(dis%idomain == idomain1d)) then
                write (errmsg, fmtidomerr) trim(this%text)
                call store_error(errmsg, terminate=.TRUE.)
              end if
            end select
          end select
        case default
          write (errmsg, '(a,3(1x,a))') &
            'UNKNOWN', trim(adjustl(this%text)), 'PACKAGEDATA:', trim(keyword)
          call store_error(errmsg)
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF FMI PACKAGEDATA'
    end if
  end subroutine read_packagedata

  !> @brief Initialize the budget file reader
  subroutine initialize_bfr(this)
    class(FlowModelInterfaceType) :: this
    integer(I4B) :: ncrbud
    call this%bfr%initialize(this%iubud, this%iout, ncrbud)
    ! todo: need to run through the budget terms
    ! and do some checking
  end subroutine initialize_bfr

  !> @brief Advance the budget file reader
  !!
  !! Advance the budget file reader by reading the next chunk
  !! of information for the current time step and stress period.
  !<
  subroutine advance_bfr(this)
    ! -- modules
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    logical :: success
    integer(I4B) :: n
    integer(I4B) :: ipos
    integer(I4B) :: nu, nr
    integer(I4B) :: ip, i
    logical :: readnext
    ! -- format
    character(len=*), parameter :: fmtkstpkper = &
                 "(1x,/1x,'FMI READING BUDGET TERMS &
                 &FOR KSTP ', i0, ' KPER ', i0)"
    character(len=*), parameter :: fmtbudkstpkper = &
      "(1x,/1x, 'FMI SETTING BUDGET TERMS &
      &FOR KSTP ', i0, ' AND KPER ',        &
      &i0, ' TO BUDGET FILE TERMS FROM &
      &KSTP ', i0, ' AND KPER ', i0)"
    !
    ! -- If the latest record read from the budget file is from a stress
    ! -- period with only one time step, reuse that record (do not read a
    ! -- new record) if the running model is still in that same stress period,
    ! -- or if that record is the last one in the budget file.
    readnext = .true.
    if (kstp * kper > 1) then
      if (this%bfr%kstp == 1) then
        if (this%bfr%kpernext == kper + 1) then
          readnext = .false.
        else if (this%bfr%endoffile) then
          readnext = .false.
        end if
      else if (this%bfr%endoffile) then
        write (errmsg, '(4x,a)') 'REACHED END OF GWF BUDGET &
          &FILE BEFORE READING SUFFICIENT BUDGET INFORMATION FOR THIS &
          &GWT SIMULATION.'
        call store_error(errmsg)
        call store_error_unit(this%iubud)
      end if
    end if
    !
    ! -- Read the next record
    if (readnext) then
      !
      ! -- Write the current time step and stress period
      write (this%iout, fmtkstpkper) kstp, kper
      !
      ! -- loop through the budget terms for this stress period
      !    i is the counter for gwf flow packages
      ip = 1
      do n = 1, this%bfr%nbudterms
        call this%bfr%read_record(success, this%iout)
        if (.not. success) then
          write (errmsg, '(4x,a)') 'GWF BUDGET READ NOT SUCCESSFUL'
          call store_error(errmsg)
          call store_error_unit(this%iubud)
        end if
        !
        ! -- Ensure kper is same between model and budget file
        if (kper /= this%bfr%kper) then
          write (errmsg, '(4x,a)') 'PERIOD NUMBER IN BUDGET FILE &
            &DOES NOT MATCH PERIOD NUMBER IN TRANSPORT MODEL.  IF THERE &
            &IS MORE THAN ONE TIME STEP IN THE BUDGET FILE FOR A GIVEN &
            &STRESS PERIOD, BUDGET FILE TIME STEPS MUST MATCH GWT MODEL &
            &TIME STEPS ONE-FOR-ONE IN THAT STRESS PERIOD.'
          call store_error(errmsg)
          call store_error_unit(this%iubud)
        end if
        !
        ! -- if budget file kstp > 1, then kstp must match
        if (this%bfr%kstp > 1 .and. (kstp /= this%bfr%kstp)) then
          write (errmsg, '(4x,a)') 'TIME STEP NUMBER IN BUDGET FILE &
            &DOES NOT MATCH TIME STEP NUMBER IN TRANSPORT MODEL.  IF THERE &
            &IS MORE THAN ONE TIME STEP IN THE BUDGET FILE FOR A GIVEN STRESS &
            &PERIOD, BUDGET FILE TIME STEPS MUST MATCH GWT MODEL TIME STEPS &
            &ONE-FOR-ONE IN THAT STRESS PERIOD.'
          call store_error(errmsg)
          call store_error_unit(this%iubud)
        end if
        !
        ! -- parse based on the type of data, and compress all user node
        !    numbers into reduced node numbers
        select case (trim(adjustl(this%bfr%budtxt)))
        case ('FLOW-JA-FACE')
          !
          ! -- bfr%flowja contains only reduced connections so there is
          !    a one-to-one match with this%gwfflowja
          do ipos = 1, size(this%bfr%flowja)
            this%gwfflowja(ipos) = this%bfr%flowja(ipos)
          end do
        case ('DATA-SPDIS')
          do i = 1, this%bfr%nlist
            nu = this%bfr%nodesrc(i)
            nr = this%dis%get_nodenumber(nu, 0)
            if (nr <= 0) cycle
            this%gwfspdis(1, nr) = this%bfr%auxvar(1, i)
            this%gwfspdis(2, nr) = this%bfr%auxvar(2, i)
            this%gwfspdis(3, nr) = this%bfr%auxvar(3, i)
          end do
        case ('DATA-SAT')
          do i = 1, this%bfr%nlist
            nu = this%bfr%nodesrc(i)
            nr = this%dis%get_nodenumber(nu, 0)
            if (nr <= 0) cycle
            this%gwfsat(nr) = this%bfr%auxvar(1, i)
          end do
        case ('STO-SS')
          do nu = 1, this%dis%nodesuser
            nr = this%dis%get_nodenumber(nu, 0)
            if (nr <= 0) cycle
            this%gwfstrgss(nr) = this%bfr%flow(nu)
          end do
        case ('STO-SY')
          do nu = 1, this%dis%nodesuser
            nr = this%dis%get_nodenumber(nu, 0)
            if (nr <= 0) cycle
            this%gwfstrgsy(nr) = this%bfr%flow(nu)
          end do
        case default
          call this%gwfpackages(ip)%copy_values( &
            this%bfr%nlist, &
            this%bfr%nodesrc, &
            this%bfr%flow, &
            this%bfr%auxvar)
          do i = 1, this%gwfpackages(ip)%nbound
            nu = this%gwfpackages(ip)%nodelist(i)
            nr = this%dis%get_nodenumber(nu, 0)
            this%gwfpackages(ip)%nodelist(i) = nr
          end do
          ip = ip + 1
        end select
      end do
    else
      !
      ! -- write message to indicate that flows are being reused
      write (this%iout, fmtbudkstpkper) kstp, kper, this%bfr%kstp, this%bfr%kper
      !
      ! -- set the flag to indicate that flows were not updated
      this%iflowsupdated = 0
    end if
  end subroutine advance_bfr

  !> @brief Finalize the budget file reader
  subroutine finalize_bfr(this)
    class(FlowModelInterfaceType) :: this
    call this%bfr%finalize()
  end subroutine finalize_bfr

  !> @brief Initialize the head file reader
  subroutine initialize_hfr(this)
    class(FlowModelInterfaceType) :: this
    call this%hfr%initialize(this%iuhds, this%iout)
    ! todo: need to run through the head terms
    ! and do some checking
  end subroutine initialize_hfr

  !> @brief Advance the head file reader
  subroutine advance_hfr(this)
    ! modules
    use TdisModule, only: kstp, kper
    class(FlowModelInterfaceType) :: this
    integer(I4B) :: nu, nr, i, ilay
    integer(I4B) :: ncpl
    real(DP) :: val
    logical :: readnext
    logical :: success
    character(len=*), parameter :: fmtkstpkper = &
                         "(1x,/1x,'FMI READING HEAD FOR &
                         &KSTP ', i0, ' KPER ', i0)"
    character(len=*), parameter :: fmthdskstpkper = &
      "(1x,/1x, 'FMI SETTING HEAD FOR KSTP ', i0, ' AND KPER ',        &
      &i0, ' TO BINARY FILE HEADS FROM KSTP ', i0, ' AND KPER ', i0)"
    !
    ! -- If the latest record read from the head file is from a stress
    ! -- period with only one time step, reuse that record (do not read a
    ! -- new record) if the running model is still in that same stress period,
    ! -- or if that record is the last one in the head file.
    readnext = .true.
    if (kstp * kper > 1) then
      if (this%hfr%kstp == 1) then
        if (this%hfr%kpernext == kper + 1) then
          readnext = .false.
        else if (this%hfr%endoffile) then
          readnext = .false.
        end if
      else if (this%hfr%endoffile) then
        write (errmsg, '(4x,a)') 'REACHED END OF GWF HEAD &
          &FILE BEFORE READING SUFFICIENT HEAD INFORMATION FOR THIS &
          &GWT SIMULATION.'
        call store_error(errmsg)
        call store_error_unit(this%iuhds)
      end if
    end if
    !
    ! -- Read the next record
    if (readnext) then
      !
      ! -- write to list file that heads are being read
      write (this%iout, fmtkstpkper) kstp, kper
      !
      ! -- loop through the layered heads for this time step
      do ilay = 1, this%hfr%nlay
        !
        ! -- read next head chunk
        call this%hfr%read_record(success, this%iout)
        if (.not. success) then
          write (errmsg, '(4x,a)') 'GWF HEAD READ NOT SUCCESSFUL'
          call store_error(errmsg)
          call store_error_unit(this%iuhds)
        end if
        !
        ! -- Ensure kper is same between model and head file
        if (kper /= this%hfr%kper) then
          write (errmsg, '(4x,a)') 'PERIOD NUMBER IN HEAD FILE &
            &DOES NOT MATCH PERIOD NUMBER IN TRANSPORT MODEL.  IF THERE &
            &IS MORE THAN ONE TIME STEP IN THE HEAD FILE FOR A GIVEN STRESS &
            &PERIOD, HEAD FILE TIME STEPS MUST MATCH GWT MODEL TIME STEPS &
            &ONE-FOR-ONE IN THAT STRESS PERIOD.'
          call store_error(errmsg)
          call store_error_unit(this%iuhds)
        end if
        !
        ! -- if head file kstp > 1, then kstp must match
        if (this%hfr%kstp > 1 .and. (kstp /= this%hfr%kstp)) then
          write (errmsg, '(4x,a)') 'TIME STEP NUMBER IN HEAD FILE &
            &DOES NOT MATCH TIME STEP NUMBER IN TRANSPORT MODEL.  IF THERE &
            &IS MORE THAN ONE TIME STEP IN THE HEAD FILE FOR A GIVEN STRESS &
            &PERIOD, HEAD FILE TIME STEPS MUST MATCH GWT MODEL TIME STEPS &
            &ONE-FOR-ONE IN THAT STRESS PERIOD.'
          call store_error(errmsg)
          call store_error_unit(this%iuhds)
        end if
        !
        ! -- fill the head array for this layer and
        !    compress into reduced form
        ncpl = size(this%hfr%head)
        do i = 1, ncpl
          nu = (ilay - 1) * ncpl + i
          nr = this%dis%get_nodenumber(nu, 0)
          val = this%hfr%head(i)
          if (nr > 0) this%gwfhead(nr) = val
        end do
      end do
    else
      write (this%iout, fmthdskstpkper) kstp, kper, this%hfr%kstp, this%hfr%kper
    end if
  end subroutine advance_hfr

  !> @brief Finalize the head file reader
  subroutine finalize_hfr(this)
    class(FlowModelInterfaceType) :: this
    close (this%iuhds)
  end subroutine finalize_hfr

  !> @brief Initialize the grid file reader
  subroutine initialize_gfr(this)
    class(FlowModelInterfaceType) :: this
    call this%gfr%initialize(this%iugrb, this%iout)
  end subroutine initialize_gfr

  !> @brief Finalize the grid file reader
  subroutine finalize_gfr(this)
    class(FlowModelInterfaceType) :: this
    call this%gfr%finalize()
  end subroutine finalize_gfr

  !> @brief Initialize gwf terms from budget file
  !!
  !! initialize terms and figure out how many
  !! different terms and packages are contained within the file
  !<
  subroutine initialize_gwfterms_from_bfr(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use SimModule, only: store_error, store_error_unit, count_errors
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    integer(I4B) :: nflowpack
    integer(I4B) :: i, ip
    integer(I4B) :: naux
    logical :: found_flowja
    logical :: found_dataspdis
    logical :: found_datasat
    logical :: found_stoss
    logical :: found_stosy
    integer(I4B), dimension(:), allocatable :: imap
    !
    ! -- Calculate the number of gwf flow packages
    allocate (imap(this%bfr%nbudterms))
    imap(:) = 0
    nflowpack = 0
    found_flowja = .false.
    found_dataspdis = .false.
    found_datasat = .false.
    found_stoss = .false.
    found_stosy = .false.
    do i = 1, this%bfr%nbudterms
      select case (trim(adjustl(this%bfr%budtxtarray(i))))
      case ('FLOW-JA-FACE')
        found_flowja = .true.
      case ('DATA-SPDIS')
        found_dataspdis = .true.
      case ('DATA-SAT')
        found_datasat = .true.
      case ('STO-SS')
        found_stoss = .true.
        this%igwfstrgss = 1
      case ('STO-SY')
        found_stosy = .true.
        this%igwfstrgsy = 1
      case default
        nflowpack = nflowpack + 1
        imap(i) = 1
      end select
    end do
    !
    ! -- allocate gwfpackage arrays
    call this%allocate_gwfpackages(nflowpack)
    !
    ! -- Copy the package name and aux names from budget file reader
    !    to the gwfpackages derived-type variable
    ip = 1
    do i = 1, this%bfr%nbudterms
      if (imap(i) == 0) cycle
      call this%gwfpackages(ip)%set_name(this%bfr%dstpackagenamearray(i), &
                                         this%bfr%budtxtarray(i))
      naux = this%bfr%nauxarray(i)
      call this%gwfpackages(ip)%set_auxname(naux, this%bfr%auxtxtarray(1:naux, i))
      ip = ip + 1
    end do
    !
    ! -- Copy just the package names for the boundary packages into
    !    the flowpacknamearray
    ip = 1
    do i = 1, size(imap)
      if (imap(i) == 1) then
        this%flowpacknamearray(ip) = this%bfr%dstpackagenamearray(i)
        ip = ip + 1
      end if
    end do
    !
    ! -- Error if specific discharge, saturation or flowja not found
    if (.not. found_dataspdis) then
      write (errmsg, '(4x,a)') 'SPECIFIC DISCHARGE NOT FOUND IN &
                              &BUDGET FILE. SAVE_SPECIFIC_DISCHARGE AND &
                              &SAVE_FLOWS MUST BE ACTIVATED IN THE NPF PACKAGE.'
      call store_error(errmsg)
    end if
    if (.not. found_datasat) then
      write (errmsg, '(4x,a)') 'SATURATION NOT FOUND IN &
                              &BUDGET FILE. SAVE_SATURATION AND &
                              &SAVE_FLOWS MUST BE ACTIVATED IN THE NPF PACKAGE.'
      call store_error(errmsg)
    end if
    if (.not. found_flowja) then
      write (errmsg, '(4x,a)') 'FLOWJA NOT FOUND IN &
                              &BUDGET FILE. SAVE_FLOWS MUST &
                              &BE ACTIVATED IN THE NPF PACKAGE.'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine initialize_gwfterms_from_bfr

  !> @brief Initialize gwf terms from a GWF exchange
  subroutine initialize_gwfterms_from_gwfbndlist(this)
    ! -- modules
    use BndModule, only: BndType, GetBndFromList
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    ! -- local
    integer(I4B) :: ngwfpack
    integer(I4B) :: ngwfterms
    integer(I4B) :: ip
    integer(I4B) :: imover
    integer(I4B) :: ntomvr
    integer(I4B) :: iterm
    character(len=LENPACKAGENAME) :: budtxt
    class(BndType), pointer :: packobj => null()
    !
    ! -- determine size of gwf terms
    ngwfpack = this%gwfbndlist%Count()
    !
    ! -- Count number of to-mvr terms, but do not include advanced packages
    !    as those mover terms are not losses from the cell, but rather flows
    !    within the advanced package
    ntomvr = 0
    do ip = 1, ngwfpack
      packobj => GetBndFromList(this%gwfbndlist, ip)
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        ntomvr = ntomvr + 1
      end if
    end do
    !
    ! -- Allocate arrays in fmi of size ngwfterms, which is the number of
    !    packages plus the number of packages with mover terms.
    ngwfterms = ngwfpack + ntomvr
    call this%allocate_gwfpackages(ngwfterms)
    !
    ! -- Assign values in the fmi package
    iterm = 1
    do ip = 1, ngwfpack
      !
      ! -- set and store names
      packobj => GetBndFromList(this%gwfbndlist, ip)
      budtxt = adjustl(packobj%text)
      call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
      this%flowpacknamearray(iterm) = packobj%packName
      iterm = iterm + 1
      !
      ! -- if this package has a mover associated with it, then add another
      !    term that corresponds to the mover flows
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        budtxt = trim(adjustl(packobj%text))//'-TO-MVR'
        call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
        this%flowpacknamearray(iterm) = packobj%packName
        this%igwfmvrterm(iterm) = 1
        iterm = iterm + 1
      end if
    end do
  end subroutine initialize_gwfterms_from_gwfbndlist

  !> @brief Allocate budget packages
  !!
  !! gwfpackages is an array of PackageBudget objects.
  !! This routine allocates gwfpackages to the proper size and initializes some
  !! member variables.
  !<
  subroutine allocate_gwfpackages(this, ngwfterms)
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(FlowModelInterfaceType) :: this
    integer(I4B), intent(in) :: ngwfterms
    ! -- local
    integer(I4B) :: n
    character(len=LENMEMPATH) :: memPath
    !
    ! -- direct allocate
    allocate (this%gwfpackages(ngwfterms))
    allocate (this%flowpacknamearray(ngwfterms))
    !
    ! -- mem_allocate
    call mem_allocate(this%igwfmvrterm, ngwfterms, 'IGWFMVRTERM', this%memoryPath)
    !
    ! -- initialize
    this%nflowpack = ngwfterms
    do n = 1, this%nflowpack
      this%igwfmvrterm(n) = 0
      this%flowpacknamearray(n) = ''
      !
      ! -- Create a mempath for each individual flow package data set
      !    of the form, MODELNAME/FMI-FTn
      write (memPath, '(a, i0)') trim(this%memoryPath)//'-FT', n
      call this%gwfpackages(n)%initialize(memPath)
    end do
  end subroutine allocate_gwfpackages

  !> @brief Deallocate memory in the gwfpackages array
  subroutine deallocate_gwfpackages(this)
    class(FlowModelInterfaceType) :: this
    integer(I4B) :: n

    do n = 1, this%nflowpack
      call this%gwfpackages(n)%da()
    end do
  end subroutine deallocate_gwfpackages

  !> @brief Find the package index for the package with the given name
  subroutine get_package_index(this, name, idx)
    use BndModule, only: BndType, GetBndFromList
    class(FlowModelInterfaceType) :: this
    character(len=*), intent(in) :: name
    integer(I4B), intent(inout) :: idx
    ! -- local
    integer(I4B) :: ip
    !
    ! -- Look through all the packages and return the index with name
    idx = 0
    do ip = 1, size(this%flowpacknamearray)
      if (this%flowpacknamearray(ip) == name) then
        idx = ip
        exit
      end if
    end do
    if (idx == 0) then
      call store_error('Error in get_package_index.  Could not find '//name, &
                       terminate=.TRUE.)
    end if
  end subroutine get_package_index

end module FlowModelInterfaceModule
