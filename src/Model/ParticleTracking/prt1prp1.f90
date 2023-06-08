module PrtPrpModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DEM1, DONE, LENFTYPE, LINELENGTH, &
                             LENBOUNDNAME, LENPAKLOC, TABLEFT, TABCENTER, &
                             MNORMAL
  use BndModule, only: BndType
  use ObsModule, only: DefaultObsIdProcessor
  use TableModule, only: TableType, table_cr
  use TimeSeriesModule, only: TimeSeriesType
  use TimeSeriesRecordModule, only: TimeSeriesRecordType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  GetTimeSeriesLinkFromList
  use BlockParserModule, only: BlockParserType
  use PrtFmiModule, only: PrtFmiType
  use ParticleModule, only: ParticleListType
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_warning
  use SimVariablesModule, only: errmsg, warnmsg
  use ArrayHandlersModule, only: expandarray
  use GlobalDataModule
  use TrackDataModule, only: TrackDataType

  implicit none

  private
  public :: PrtPrpType
  public :: prp_create

  character(len=LENFTYPE) :: ftype = 'PRP'
  character(len=16) :: text = '             PRP'

  type, extends(BndType) :: PrtPrpType
    type(PrtFmiType), pointer :: fmi => null() !< flow model interface
    ! type(ParticleListType), dimension(:), pointer :: partlist       => null()  !< list of particle data
    type(ParticleListType), pointer :: partlist => null() !< list of particle data for the package
    integer(I4B), pointer :: nreleasepts => null() !< number of particle release points
    integer(I4B), pointer :: npart => null() !< number of particles in the particle data list
    integer(I4B), pointer :: npartmax => null() !< maximum number of particles in the particle data list
    real(DP), pointer :: stoptime => null() !< stop time for all particles
    real(DP), pointer :: stoptraveltime => null() !< stop travel time for all particles
    integer(I4B), pointer :: istopweaksink => null() !< weak sink option: 0 = do not stop, 1 = stop
    !< extend final steady state option: 0 = do not extend, 1 = extend
    ! integer(I4B), pointer :: iextendfinalss => null()
    integer(I4B), pointer :: istopzone => null() !< optional stop zone number; 0 = no stop zone
    integer(I4B), pointer :: ioutinactive => null() !< output for inactive particles: 0 = no output, 1 = output
    integer(I4B), pointer :: idrape => null() !< drape option: 0 = do not drape, 1 = drape to topmost active cell
    integer(I4B), dimension(:), pointer, contiguous :: noder => null() !< reduced node number of release point
    real(DP), dimension(:), pointer, contiguous :: x => null() !< x coordinate of particle release point
    real(DP), dimension(:), pointer, contiguous :: y => null() !< y coordinate of particle release point
    real(DP), dimension(:), pointer, contiguous :: z => null() !< z coordinate of particle release point
    ! real(DP), dimension(:), pointer, contiguous :: tbegin => null() !< begin time of particle release point
    ! real(DP), dimension(:), pointer, contiguous :: trepeat => null() !< repeat time interval of release point
    ! real(DP), dimension(:), pointer, contiguous :: tend => null() !< end time of particle release point
    !< stop time of particles released by particle release point
    ! kluge note: don't need this array if going with "global" stoptime value
    real(DP), dimension(:), pointer, contiguous :: tstop => null()
    character(len=LENBOUNDNAME), dimension(:), pointer, contiguous :: rptname &
                                                                      => null() !< release point name
    real(DP), dimension(:), pointer, contiguous :: massrls => null() !< mass released during time step
    ! real(DP), dimension(:), pointer, contiguous :: porosity => null() !< aquifer porosity
    ! real(DP), dimension(:), pointer, contiguous :: retfactor => null() !< retardation factor
    ! integer(I4B), dimension(:), pointer, contiguous :: izone => null() !< zone number
    integer(I4B), allocatable, dimension(:) :: kstp_list_rls !< allocatable time steps for releases in period
    integer(I4B), pointer :: ifreq_rls => null() !< release frequency (time steps) in period
    logical(LGP), pointer :: rls_first => null() !< flag for release on first time step in period
    logical(LGP), pointer :: rls_all => null() !< flag for release on all time steps in period
    logical(LGP), pointer :: rls_any => null() !< flag that indicates whether any release in period
    logical(LGP), pointer :: noperiodblocks => null() !< flag indicating if there are no period blocks in sim
    integer(I4B), pointer :: itrack1 => null() ! pointer to start of prp track data
    integer(I4B), pointer :: itrack2 => null() ! pointer to end of prp track data
    type(TrackDataType), pointer :: trackdata => null() ! pointer to model track data
    integer(I4B), pointer :: itrkout => null()
    integer(I4B), pointer :: itrkhdr => null()
    integer(I4B), pointer :: itrkcsv => null()

  contains

    procedure :: prp_allocate_arrays
    procedure :: prp_allocate_scalars
    procedure :: bnd_ar => prp_ar
    procedure :: bnd_ad => prp_ad
    procedure :: bnd_rp => prp_rp
    ! procedure :: bnd_cf => prp_cf
    ! procedure :: bnd_df => prp_df
    ! procedure :: bnd_fc => prp_fc
    procedure :: bnd_cq_simrate => prp_cq_simrate
    procedure :: bnd_da => prp_da
    procedure :: define_listlabel
    procedure :: prp_ot_trk
    procedure :: prp_set_pointers ! kluge?
    procedure :: bnd_options => prp_options
    procedure :: read_dimensions => prp_read_dimensions
    procedure :: prp_read_packagedata
    ! procedure :: read_data
    procedure :: sav_particles
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => prp_obs_supported
    procedure, public :: bnd_df_obs => prp_df_obs
    ! ! -- methods for time series
    ! procedure, public :: bnd_rp_ts => prp_rp_ts
  end type PrtPrpType

contains

  !> @brief Create a new particle release package
  subroutine prp_create(packobj, id, ibcnum, inunit, iout, namemodel, &
                        pakname, fmi)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    type(PrtFmiType), pointer :: fmi
    ! -- local
    type(PrtPrpType), pointer :: prpobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (prpobj)
    packobj => prpobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call prpobj%prp_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 4
    packobj%iscloc = 1
    !
    ! -- Store pointer to flow model interface
    prpobj%fmi => fmi
    !
    ! -- return
    return
  end subroutine prp_create

  !> @brief Deallocate memory
  !<
  subroutine prp_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(PrtPrpType) :: this
    !
    ! -- Deallocate parent package
    call this%BndType%bnd_da()
    !
    ! -- scalars
    call mem_deallocate(this%stoptime)
    call mem_deallocate(this%stoptraveltime)
    call mem_deallocate(this%istopweaksink)
    ! call mem_deallocate(this%iextendfinalss)
    call mem_deallocate(this%istopzone)
    call mem_deallocate(this%ioutinactive)
    call mem_deallocate(this%idrape)
    call mem_deallocate(this%nreleasepts)
    call mem_deallocate(this%ifreq_rls)
    call mem_deallocate(this%rls_first)
    call mem_deallocate(this%rls_all)
    call mem_deallocate(this%rls_any)
    call mem_deallocate(this%noperiodblocks)
    call mem_deallocate(this%npart)
    call mem_deallocate(this%npartmax)
    call mem_deallocate(this%itrkout)
    call mem_deallocate(this%itrkhdr)
    call mem_deallocate(this%itrkcsv)
    !
    ! -- arrays
    call mem_deallocate(this%noder)
    call mem_deallocate(this%x)
    call mem_deallocate(this%y)
    call mem_deallocate(this%z)
    ! call mem_deallocate(this%tbegin)
    ! call mem_deallocate(this%trepeat)
    ! call mem_deallocate(this%tend)
    call mem_deallocate(this%tstop)
    call mem_deallocate(this%rptname, 'RPTNAME', this%memoryPath)
    ! call mem_deallocate(this%porosity)
    ! call mem_deallocate(this%retfactor)
    ! call mem_deallocate(this%izone)
    deallocate (this%partlist%x) ! kluge note: use mem_deallocate for these arrays
    deallocate (this%partlist%y)
    deallocate (this%partlist%z)
    deallocate (this%partlist%iTrackingDomain)
    deallocate (this%partlist%iTrackingDomainBoundary)
    deallocate (this%partlist%trelease)
    deallocate (this%partlist%tstop)
    deallocate (this%partlist%ttrack)
    deallocate (this%partlist%istopweaksink)
    deallocate (this%partlist%istopzone)
    deallocate (this%partlist%istatus)
    deallocate (this%partlist%irpt)
    deallocate (this%partlist) ! kluge note: structure of arrays
    call mem_deallocate(this%massrls)
    !
    ! -- allocatable array (not pointer)
    if (allocated(this%kstp_list_rls)) deallocate (this%kstp_list_rls)
    !
    ! -- return
    return
  end subroutine prp_da

  !> @ brief Set pointers to model variables
  !<
  subroutine prp_set_pointers(this, ibound, itrack1, itrack2, trackdata)
    ! -- dummy variables
    class(PrtPrpType) :: this !< PrtPrpType object
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    integer(I4B), pointer :: itrack1
    integer(I4B), pointer :: itrack2
    type(TrackDataType), pointer :: trackdata
    !
    ! -- Set pointer to PRT model ibound
    this%ibound => ibound
    !
    ! -- Set pointers to track data
    this%itrack1 => itrack1
    this%itrack2 => itrack2
    this%trackdata => trackdata
    !
    ! -- return
    return
  end subroutine prp_set_pointers

  !> @brief Allocate arrays
  !<
  subroutine prp_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PrtPrpType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    integer(I4B) :: nps
    !
    ! ! -- call standard BndType allocate scalars
    ! call this%BndType%allocate_arrays()
    ! !
    ! -- Allocate
    call mem_allocate(this%noder, this%nreleasepts, 'NODER', this%memoryPath)
    call mem_allocate(this%x, this%nreleasepts, 'X', this%memoryPath)
    call mem_allocate(this%y, this%nreleasepts, 'Y', this%memoryPath)
    call mem_allocate(this%z, this%nreleasepts, 'Z', this%memoryPath)
    ! call mem_allocate(this%tbegin, this%nreleasepts, 'TBEGIN', this%memoryPath)
    ! call mem_allocate(this%trepeat, this%nreleasepts, 'TREPEAT', this%memoryPath)
    ! call mem_allocate(this%tend, this%nreleasepts, 'TEND', this%memoryPath)
    call mem_allocate(this%tstop, this%nreleasepts, 'TSTOP', this%memoryPath)
    call mem_allocate(this%rptname, LENBOUNDNAME, this%nreleasepts, &
                      'RPTNAME', this%memoryPath)
    ! call mem_allocate(this%porosity, this%dis%nodes, 'POROSITY',              &
    !                   this%memoryPath)
    ! call mem_allocate(this%retfactor, this%dis%nodes, 'RETFACTOR',            &
    !                   this%memoryPath)
    ! call mem_allocate(this%izone, this%dis%nodes, 'IZONE', this%memoryPath)
    allocate (this%partlist) ! kluge note: structure of arrays
    ! allocate(this%partlist%velmult(this%npartmax))
    allocate (this%partlist%x(this%npartmax)) ! kluge note: nprtmax is the initial max dimension
    allocate (this%partlist%y(this%npartmax)) ! kluge note: use mem_allocate for these arrays
    allocate (this%partlist%z(this%npartmax))
    ! kluge note: ditch crazy dims
    allocate (this%partlist%iTrackingDomain(this%npartmax, &
                                            levelMin:levelMax))
    ! kluge note: ditch crazy dims
    allocate (this%partlist%iTrackingDomainBoundary(this%npartmax, &
                                                    levelMin:levelMax))
    allocate (this%partlist%trelease(this%npartmax))
    allocate (this%partlist%tstop(this%npartmax))
    allocate (this%partlist%ttrack(this%npartmax))
    allocate (this%partlist%istopweaksink(this%npartmax))
    allocate (this%partlist%istopzone(this%npartmax))
    allocate (this%partlist%istatus(this%npartmax))
    allocate (this%partlist%irpt(this%npartmax))
    call mem_allocate(this%massrls, this%nreleasepts, 'MASSRLS', this%memoryPath)
    !
    ! -- Intialize
    do nps = 1, this%nreleasepts
      this%massrls(nps) = DZERO
    end do
    !
    ! -- The following array is allocatable (not a pointer) so it can be resized using
    if (allocated(this%kstp_list_rls)) deallocate (this%kstp_list_rls)
    allocate (this%kstp_list_rls(0))
    !
    ! -- Return
    return
  end subroutine prp_allocate_arrays

  !> @brief Allocate scalar members
  !<
  subroutine prp_allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PrtPrpType) :: this
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- Allocate scalars for this type
    call mem_allocate(this%stoptime, 'STOPTIME', this%memoryPath)
    call mem_allocate(this%stoptraveltime, 'STOPTRAVELTIME', this%memoryPath)
    call mem_allocate(this%istopweaksink, 'ISTOPWEAKSINK', this%memoryPath)
    ! call mem_allocate(this%iextendfinalss, 'IEXTENDFINALSS', this%memoryPath)
    call mem_allocate(this%istopzone, 'ISTOPZONE', this%memoryPath)
    call mem_allocate(this%ioutinactive, 'IOUTINACTIVE', this%memoryPath)
    call mem_allocate(this%idrape, 'IDRAPE', this%memoryPath)
    call mem_allocate(this%nreleasepts, 'NRELEASEPTS', this%memoryPath)
    call mem_allocate(this%ifreq_rls, 'IFREQ_RLS', this%memoryPath)
    call mem_allocate(this%rls_first, 'RLS_FIRST', this%memoryPath)
    call mem_allocate(this%rls_all, 'RLS_ALL', this%memoryPath)
    call mem_allocate(this%rls_any, 'RLS_ANY', this%memoryPath)
    call mem_allocate(this%noperiodblocks, 'NOPERIODBLOCKS', this%memoryPath)
    call mem_allocate(this%npart, 'NPART', this%memoryPath)
    call mem_allocate(this%npartmax, 'NPARTMAX', this%memoryPath)
    call mem_allocate(this%itrkout, 'ITRKOUT', this%memoryPath)
    call mem_allocate(this%itrkhdr, 'ITRKHDR', this%memoryPath)
    call mem_allocate(this%itrkcsv, 'ITRKCSV', this%memoryPath)
    !
    ! -- Set values
    this%stoptime = huge(1d0) ! kluge???
    this%stoptraveltime = huge(1d0) ! kluge???
    this%istopweaksink = 0
    ! this%iextendfinalss = 0
    this%istopzone = 0
    this%ioutinactive = 0
    this%idrape = 0
    this%nreleasepts = 0
    this%ifreq_rls = 0
    this%rls_first = .false.
    this%rls_all = .false.
    this%rls_any = .false.
    this%noperiodblocks = .false.
    this%npart = 0
    this%npartmax = 0
    this%itrkout = 0
    this%itrkhdr = 0
    this%itrkcsv = 0
    !
    ! -- return
    return
  end subroutine prp_allocate_scalars

  !> @ brief Allocate and read period data
  !<
  subroutine prp_ar(this)
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this !< PrtPrpType object
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- allocate and read observations
    call this%obs%obs_ar()
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
    !
    ! -- set boundname for each release point
    if (this%inamedbound /= 0) then
      do n = 1, this%nreleasepts
        this%boundname(n) = this%rptname(n)
      end do
    end if
    !
    ! -- copy noder into nodelist
    do n = 1, this%nreleasepts
      this%nodelist(n) = this%noder(n)
    end do
    !
    ! !
    ! ! -- setup pakmvrobj
    ! if (this%imover /= 0) then
    !   allocate(this%pakmvrobj)
    !   call this%pakmvrobj%ar(this%maxbound, this%maxbound, this%memoryPath)
    ! endif
    ! !
    ! -- return
    return
  end subroutine prp_ar

  !> @brief Advance a time step & release new particles according to period data
  !<
  subroutine prp_ad(this)
    ! -- modules
    use TdisModule, only: kstp, totimc
    ! -- dummy
    class(PrtPrpType) :: this
    ! -- local
    integer(I4B) :: i, n, ic
    integer(I4B) :: nps, np
    real(DP) :: trelease, tstop ! kluge?
    ! real(DP) :: top, bot, sat
    logical(LGP) :: isRelease
    !
    ! -- Reset particle mass released for time step
    do nps = 1, this%nreleasepts
      this%massrls(nps) = DZERO
    end do
    !
    ! -- Return early if there are no particles to release
    if (.not. this%rls_any) then
      return
    end if

    ! -- Check if there is to be a release at the start of this time step
    isRelease = .false.
    ! -- release all time steps
    if (this%rls_all) then
      isRelease = .true.
    else
      ! -- release on the first time step
      if (this%rls_first) then
        if (kstp == 1) isRelease = .true.
      end if
      ! -- release only after every this%ifreq_rls time steps elapse
      if (this%ifreq_rls > 0) then
        if ((kstp / this%ifreq_rls) * this%ifreq_rls == kstp) & ! kluge note: use modulo?
          isRelease = .true.
      end if
      ! -- ???
      n = size(this%kstp_list_rls)
      if (n > 0) then
        do i = 1, n
          ! kluge note: store advancing counter to avoid re-searching entire array each time?
          if (this%kstp_list_rls(i) == kstp) isRelease = .true.
        end do
      end if
    end if
    !
    ! -- Do the release, if there is one
    if (isRelease) then
      do nps = 1, this%nreleasepts
        ic = this%noder(nps) ! reduced node number (cell ID)
        ! -- If drape option activated, release particle in highest active
        ! -- cell vertically below release point. If no such active cell,
        ! -- do not release particle.
        if (this%idrape /= 0) then
          if (this%ibound(ic) == 0) then
            ! -- Search for highest active cell
            call this%dis%highest_active(ic, this%ibound)
            ! -- If returned cell is inactive, do not release particle
            if (this%ibound(ic) == 0) cycle ! kluge note: somehow record for the user that a particle was scheduled but not released?
          end if
        end if
        np = this%npart + 1 ! particle index
        this%npart = np ! ???
        trelease = totimc ! release time

        ! -- Set stopping time to earlier of times specified by STOPTIME and STOPTRAVELTIME
        if (this%stoptraveltime == huge(1d0)) then ! kluge huge?
          tstop = this%stoptime
        else
          tstop = trelease + this%stoptraveltime
          if (this%stoptime < tstop) tstop = this%stoptime
        end if

        this%partlist%x(np) = this%x(nps) ! kluge note: need check that specified coords are within the cell
        this%partlist%y(np) = this%y(nps)
        this%partlist%z(np) = this%z(nps)
        this%partlist%trelease(np) = trelease ! kluge
        this%partlist%tstop(np) = tstop
        this%partlist%ttrack(np) = trelease
        this%partlist%istopweaksink(np) = this%istopweaksink
        this%partlist%istopzone(np) = this%istopzone
        this%partlist%istatus(np) = 1
        this%partlist%irpt(np) = nps
        this%partlist%iTrackingDomain(np, 0) = 0 ! kluge???
        this%partlist%iTrackingDomainBoundary(np, 0) = 0 ! kluge???
        this%partlist%iTrackingDomain(np, 1) = 0 ! kluge???
        this%partlist%iTrackingDomainBoundary(np, 1) = 0
        this%partlist%iTrackingDomain(np, 2) = ic
        this%partlist%iTrackingDomainBoundary(np, 2) = 0
        this%partlist%iTrackingDomain(np, 3) = 0
        this%partlist%iTrackingDomainBoundary(np, 3) = 0

        ! -- Each particle currently assigned unit mass
        this%massrls(nps) = this%massrls(nps) + DONE
      end do
    end if
    !
    ! -- return
    return
  end subroutine prp_ad

  !> @ brief Read and prepare period data for particle input
  !<
  subroutine prp_rp(this)
    ! -- modules
    use TdisModule, only: kper, nper
    use InputOutputModule, only: urword
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this !< PrtPrpType object
    ! -- local variables
    integer(I4B) :: ierr
    integer(I4B) :: n
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    logical(LGP) :: rls_lsp
    character(len=LINELENGTH) :: keyword
    character(len=:), allocatable :: line
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
                      "('Looking for BEGIN PERIOD iper.  &
                      &Found ', a, ' instead.')"
    character(len=*), parameter :: fmt_steps = &
                                   "(6x,'TIME STEP(S) ',50(I0,' '))" ! kluge 50 (similar to STEPS in OC)?
    character(len=*), parameter :: fmt_freq = &
                                   "(6x,'EVERY ',I0,' TIME STEP(S)')"
    !
    ! -- Set ionper to the stress period number for which a new block of data
    !    will be read.
    if (this%inunit == 0) return
    !
    ! -- get stress period data
    if (this%ionper < kper) then
      !
      ! -- get period block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then
        !
        ! -- read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else
        !
        ! -- PERIOD block not found
        if (ierr < 0) then
          if (kper == 1) then
            ! -- End of file found; no period data for the simulation.
            this%noperiodblocks = .true.
          else
            ! -- End of file found; no more period data.
            this%ionper = nper + 1
          end if
        else
          ! -- Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
    end if
    !
    ! ! -- clear period data
    ! kluge note: releases occur only when explcitly specified for a
    ! period (consider a default exception for period 1)
    !
    ! if(allocated(this%kstp_list_rls)) deallocate(this%kstp_list_rls)
    ! allocate(this%kstp_list_rls(0))
    ! this%ifreq_rls = 0
    ! this%rls_first = .false.
    ! this%rls_all = .false.
    ! this%rls_any = .false.
    ! rls_lsp = .false.
    ! !
    ! -- if no period data for simulation, single release at beginning
    if (this%noperiodblocks) then
      if (kper == 1) then
        if (allocated(this%kstp_list_rls)) deallocate (this%kstp_list_rls)
        allocate (this%kstp_list_rls(0))
        this%ifreq_rls = 0
        this%rls_first = .true.
        this%rls_all = .false.
        this%rls_any = .true.
        rls_lsp = .false.
      else if (kper == 2) then
        this%rls_first = .false.
        this%rls_any = .false.
      end if
      ! -- else read data if ionper == kper
    else if (this%ionper == kper) then
      !
      ! -- clear period data
      if (allocated(this%kstp_list_rls)) deallocate (this%kstp_list_rls)
      allocate (this%kstp_list_rls(0))
      this%ifreq_rls = 0
      this%rls_first = .false.
      this%rls_all = .false.
      this%rls_any = .false.
      rls_lsp = .false.
      !
      ! -- loop to read records
      recordloop: do
        !
        ! -- Read the line
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- Process setting
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('ALL')
          this%rls_all = .true.
          this%rls_any = .true.
        case ('STEPS')
          call this%parser%GetRemainingLine(line)
          lloc = 1
          listsearch: do
            call urword(line, lloc, istart, istop, 2, ival, rval, -1, 0)
            if (ival > 0) then
              n = size(this%kstp_list_rls)
              call expandarray(this%kstp_list_rls)
              this%kstp_list_rls(n + 1) = ival
              cycle listsearch
            end if
            exit listsearch
          end do listsearch
          this%rls_any = .true.
        case ('FREQUENCY')
          ival = this%parser%GetInteger() ! kluge note: check for nonnegative
          this%ifreq_rls = ival
          this%rls_any = .true.
        case ('FIRST')
          this%rls_first = .true.
          this%rls_any = .true.
        case default
          write (errmsg, '(2a)') &
            'Looking for ALL, STEPS, FIRST, or FREQUENCY. Found: ', &
            trim(adjustl(keyword))
          call store_error(errmsg, terminate=.TRUE.)
        end select
        !
      end do recordloop
      !
      ! -- else repeat period settings
    else
      rls_lsp = .true.
    end if
    !
    if (.not. this%rls_any) then
      write (this%iout, "(1x,/1x,a)") 'NO PARTICLE RELEASES IN THIS STRESS '// &
        'PERIOD'
    else if (rls_lsp) then
      write (this%iout, "(1x,/1x,a)") 'REUSING PARTICLE RELEASE SETTINGS '// &
        'FROM LAST STRESS PERIOD'
    else if (this%rls_all) then
      write (this%iout, "(1x,/1x,a)") 'PARTICLE RELEASE SCHEDULED AT THE '// &
        'START OF ALL TIME STEPS IN STRESS PERIOD'
    else
      write (this%iout, "(1x,/1x,a)") 'PARTICLE RELEASE SCHEDULED AT THE '// &
        'START OF EACH TIME STEP THAT MATCHES ONE OR MORE OF THE FOLLOWING:'
      if (this%rls_first) write (this%iout, "(6x,a)") 'FIRST TIME STEP '// &
        '(START OF STRESS PERIOD)'
      if (this%ifreq_rls > 0) write (this%iout, fmt_freq) this%ifreq_rls
      n = size(this%kstp_list_rls)
      if (n > 0) write (this%iout, fmt_steps) this%kstp_list_rls
    end if
    !
    ! -- return
    return
  end subroutine prp_rp

  !> @ brief Calculate simrate.
  !!
  !!  Calculate the flow between package and the model and store in the
  !!  simvals variable.
  !!
  !<
  subroutine prp_cq_simrate(this, hnew, flowja, imover)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    class(PrtPrpType) :: this !< PrtPrpType object
    real(DP), dimension(:), intent(in) :: hnew ! kluge note: not needed but part of interface
    real(DP), dimension(:), intent(inout) :: flowja !< flow between package and model
    integer(I4B), intent(in) :: imover !< flag indicating if the mover package is active
    ! -- local variables
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: idiag
    real(DP) :: tled
    real(DP) :: rrate
    ! -- formats
    !
    ! -- Set reciprocal of time step length.
    tled = DONE / delt
    !
    ! -- If no boundaries, skip flow calculations.
    if (this%nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do i = 1, this%nbound
        node = this%nodelist(i)
        !
        ! -- If cell is no-flow or constant-head, then ignore it.
        rrate = DZERO
        ! kluge note: think about condition(s) under which to ignore cell
        if (node > 0) then
          idiag = this%dis%con%ia(node)
          ! kluge note: think about condition(s) under which to ignore cell
          ! if(this%ibound(node) > 0) then
          ! -- Calculate the flow rate into the cell.
          rrate = this%massrls(i) * tled
          ! end if
          flowja(idiag) = flowja(idiag) + rrate
        end if
        !
        ! -- Save simulated value to simvals array.
        this%simvals(i) = rrate
        !
      end do
    end if
    !
    ! -- return
    return
  end subroutine prp_cq_simrate

  !> @ brief Define list heading written to iout when PRINT_INPUT option is used
  !<
  subroutine define_listlabel(this) ! kluge note: update for PRT?
    class(PrtPrpType), intent(inout) :: this
    !
    ! -- create the header list label
    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'STRESS RATE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
    !
    ! -- return
    return
  end subroutine define_listlabel

  !> @brief Indicates whether observations are supported.
  !!
  !! Return true because PRP package supports observations.
  !! Overrides BndType%bnd_obs_supported().
  !<
  logical function prp_obs_supported(this) ! kluge note: want this???
    implicit none
    class(PrtPrpType) :: this
    prp_obs_supported = .true.
    return
  end function prp_obs_supported

  !> @brief Store observation type supported by PRP package.
  !!
  !! Overrides BndType%bnd_df_obs().
  !!
  !<
  subroutine prp_df_obs(this) ! kluge note: need this???
    implicit none
    ! -- dummy
    class(PrtPrpType) :: this
    ! -- local
    integer(I4B) :: indx
    call this%obs%StoreObsType('prp', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- return
    return
  end subroutine prp_df_obs

  !> @brief Set options specific to PrtPrpType (overrides BndType%bnd_options)
  !<
  subroutine prp_options(this, option, found)
    use OpenSpecModule, only: access, form
    use ConstantsModule, only: MAXCHARLEN, DZERO
    use InputOutputModule, only: urword, getunit, openfile
    use TrackDataModule, only: TRACKHEADERS, TRACKTYPES
    ! -- dummy
    class(PrtPrpType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    ! -- locals
    character(len=MAXCHARLEN) :: fname
    character(len=MAXCHARLEN) :: keyword
    ! -- formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    !
    ! ! -- reinitialize stoptime and stoptraveltime to huge    ! kluge?
    ! this%stoptime = huge(1d0)
    ! this%stoptraveltime = huge(1d0)
    ! !
    select case (option)
    case ('STOPTIME')
      this%stoptime = this%parser%GetDouble()
      found = .true.
    case ('STOPTRAVELTIME')
      this%stoptraveltime = this%parser%GetDouble()
      found = .true.
    case ('STOP_AT_WEAK_SINK')
      this%istopweaksink = 1
      found = .true.
      ! case ('EXTEND_FINAL_SS')
      !   this%iextendfinalss = 1
      !   found = .true.
      !   print *, "EXTEND_FINAL_SS option read in but not programmed yet"  ! kluge
      !   !!pause
      !   stop
    case ('ISTOPZONE')
      this%istopzone = this%parser%GetInteger()
      found = .true.
    case ('OUTPUT_FOR_INACTIVE')
      this%ioutinactive = 1
      found = .true.
      print *, "OUTPUT_FOR_INACTIVE option read in but not programmed yet" ! kluge
      ! pause
      stop
    case ('DRAPE')
      this%idrape = 1
      found = .true.
    case ('TRACK')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        ! parse filename
        call this%parser%GetString(fname)
        ! open binary output file
        this%itrkout = getunit()
        call openfile(this%itrkout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, filstat_opt='REPLACE', &
                      mode_opt=MNORMAL)
        write (this%iout, fmttrkbin) trim(adjustl(fname)), this%itrkout
        ! open and write ascii header spec file
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
        ! open CSV output file and write headers
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
    case default
      found = .false.
    end select
    !
    ! -- Return
    return
  end subroutine prp_options

  !> @brief Read the packagedata for this package
  !<
  subroutine prp_read_packagedata(this)
    ! use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(PrtPrpType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: cellid
    character(len=LENBOUNDNAME) :: bndName
    ! character(len=LENBOUNDNAME) :: bndNameTemp
    character(len=9) :: cno
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ival
    integer(I4B) :: n
    ! integer(I4B) :: j
    ! integer(I4B) :: ii
    ! integer(I4B) :: jj
    ! integer(I4B) :: ieqn
    ! integer(I4B) :: itmp
    integer(I4B) :: ierr
    ! integer(I4B) :: idx
    ! real(DP), pointer :: bndElem => null()
    ! -- local allocatable arrays
    character(len=LENBOUNDNAME), dimension(:), allocatable :: nametxt
    ! character(len=50), dimension(:, :), allocatable :: caux
    integer(I4B), dimension(:), allocatable :: nboundchk
    integer(I4B), dimension(:), allocatable :: noder
    real(DP), dimension(:), allocatable :: x
    real(DP), dimension(:), allocatable :: y
    real(DP), dimension(:), allocatable :: z
    ! real(DP), dimension(:), allocatable :: tbegin
    ! real(DP), dimension(:), allocatable :: trepeat
    ! real(DP), dimension(:), allocatable :: tend
    real(DP), dimension(:), allocatable :: tstop
    ! -- format
    character(len=*), parameter :: fmttend = &
      "('end time (', G0, ') must be greater than or equal to the              &
     &begin time (', G0, ').')"
    !
    ! -- allocate and initialize temporary variables
    allocate (noder(this%nreleasepts)) ! kluge ?
    allocate (x(this%nreleasepts))
    allocate (y(this%nreleasepts))
    allocate (z(this%nreleasepts))
    ! allocate(tbegin(this%nreleasepts))
    ! allocate(trepeat(this%nreleasepts))
    ! allocate(tend(this%nreleasepts))
    allocate (tstop(this%nreleasepts))
    allocate (nametxt(this%nreleasepts))
    !  if (this%naux > 0) then
    !    allocate(caux(this%naux, this%nreleasepts))    ! kluge !
    !  end if
    allocate (nboundchk(this%nreleasepts))
    !
    ! -- initialize temporary variables
    do n = 1, this%nreleasepts
      nboundchk(n) = 0
    end do
    !
    ! -- read particle release point data
    ! -- get particle release points block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportopenclose=.true.)
    !
    ! -- parse block if detected
    if (isfound) then
      ! write(this%iout,'(/1x,a)')                                                &
      !   'PROCESSING ' // trim(adjustl(this%text)) // ' PACKAGEDATA'
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%packName)) &
        //' PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        ival = this%parser%GetInteger()
        n = ival

        if (n < 1 .or. n > this%nreleasepts) then
          write (errmsg, '(a,1x,i0,a)') &
            'Release point number must be greater than 0 and less than ', &
            'or equal to', this%nreleasepts, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- node number
        call this%parser%GetCellid(this%dis%ndim, cellid)
        noder(n) = this%dis%noder_from_cellid(cellid, this%inunit, this%iout)
        !
        ! -- x, y, z coordinates
        x(n) = this%parser%GetDouble()
        y(n) = this%parser%GetDouble()
        z(n) = this%parser%GetDouble()
        !
        ! ! -- begin time for release point
        ! tbegin(n) = this%parser%GetDouble()
        ! !
        ! ! -- repeat time interval for release point
        ! trepeat(n) = this%parser%GetDouble()
        ! !
        ! ! -- end time for release point
        ! tend(n) = this%parser%GetDouble()
        ! !
        ! ! -- stop time for particles released by release point
        ! tstop(n) = this%parser%GetDouble()
        ! !
        ! ! -- get aux data
        ! do jj = 1, this%naux
        !   call this%parser%GetString(caux(jj, n))   ! kluge !
        ! end do
        !
        ! -- set default bndName
        write (cno, '(i9.9)') n
        bndName = 'PRP'//cno
        !
        ! ! -- read particle release point name
        ! if (this%inamedbound /= 0) then
        !   call this%parser%GetStringCaps(bndNameTemp)
        !   if (bndNameTemp /= '') then
        !     bndName = bndNameTemp        ! kluge !
        !   end if
        ! end if
        nametxt(n) = bndName
      end do

      ! write(this%iout,'(1x,a)')                                                  &
      !   'END OF ' // trim(adjustl(this%text)) // ' PACKAGEDATA'
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%packName))//' PACKAGEDATA'
      !
      ! -- check for duplicate or missing particle release points
      do n = 1, this%nreleasepts
        if (nboundchk(n) == 0) then
          write (errmsg, '(a,a,1x,i0,a)') 'No data specified for particle ', &
            'release point', n, '.'
          call store_error(errmsg)
        else if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'Data for particle release point', n, 'specified', nboundchk(n), &
            'times.'
          call store_error(errmsg)
        end if
      end do
    else
      call store_error('Required packagedata block not found.')
    end if
    !
    ! -- terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    !   ! -- allocate particle release point data          ! kluge !
    !   call this%prp_allocate_release point_arrays()
    !
    ! -- fill particle release point data with data stored in temporary local
    ! -- arrays
    do n = 1, this%nreleasepts
      this%noder(n) = noder(n)
      this%x(n) = x(n)
      this%y(n) = y(n)
      this%z(n) = z(n)
      ! this%tbegin(n) = tbegin(n)
      ! this%trepeat(n) = trepeat(n)
      ! this%tend(n) = tend(n)
      ! this%tstop(n) = tstop(n)
      this%rptname(n) = nametxt(n)
      !
      ! -- check for error condition
      ! if (this%tend(n) <= this%tbegin(n)) then
      !   write(cstr, fmttend) this%tend(n), this%tbegin(n)  ! kluge !
      !   write(*,'(A)') cstr    ! kluge
      !   !!pause
      !   stop
      ! end if
      !
      ! ! -- fill aux data
      ! do jj = 1, this%naux
      !   text = caux(jj, n)             ! kluge !
      !   ii = n
      !   bndElem => this%mauxvar(jj, ii)
      !   call read_value_or_time_series_adv(text, ii, jj, bndElem, this%packName,     &
      !                                      'AUX', this%tsManager, this%iprpak,   &
      !                                      this%auxname(jj))
      ! end do
    end do
    !
    ! -- deallocate local storage
    deallocate (noder)
    deallocate (x)
    deallocate (y)
    deallocate (z)
    ! deallocate(tbegin)
    ! deallocate(trepeat)
    ! deallocate(tend)
    deallocate (tstop)
    ! if (this%naux > 0) then     ! kluge !
    !   deallocate(caux)
    ! end if
    deallocate (nboundchk)
    !
    ! -- return
    return
  end subroutine prp_read_packagedata

  !> @brief Read package dimensions
  subroutine prp_read_dimensions(this)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(PrtPrpType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    !
    ! -- get dimension block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimension block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING PARTICLE INPUT DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NRELEASEPTS')
          this%nreleasepts = this%parser%GetInteger()
        case default
          write (errmsg, &
                 '(4x,a,a)') '****ERROR. UNKNOWN PARTICLE INPUT DIMENSION: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF PARTICLE INPUT DIMENSIONS'
    else
      call store_error('ERROR.  REQUIRED DIMENSIONS BLOCK NOT FOUND.')
    end if
    !
    ! -- set maxbound and nbound to nreleasepts
    this%maxbound = this%nreleasepts
    this%nbound = this%nreleasepts
    !
    ! -- set max number of particles
    this%npartmax = this%nreleasepts ! kluge note: does not account for repeating releases
    this%npartmax = this%npartmax * 10 ! kluge hardwire for some breathing space intially
    !
    ! -- allocate arrays for prp package
    call this%prp_allocate_arrays()
    !
    ! ! -- read porosity, retfactor, and izone data
    ! call this%read_data()
    ! !
    ! -- read packagedata
    call this%prp_read_packagedata()
    !
    ! -- return
    return
  end subroutine prp_read_dimensions

  !> @brief Save particle track data to an output file
  subroutine prp_ot_trk(this)
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this
    !
    ! -- write particle track data to binary output file
    if (this%itrkout /= 0) &
      call this%trackdata%save_track_data(this%itrkout, csv=.false., &
                                          itrack1=this%itrack1 + 1, &
                                          itrack2=this%itrack2)
    !
    ! -- write particle track data to CSV output file
    if (this%itrkcsv /= 0) &
      call this%trackdata%save_track_data(this%itrkcsv, csv=.true., &
                                          itrack1=this%itrack1 + 1, &
                                          itrack2=this%itrack2)
    !
    return
  end subroutine prp_ot_trk

  !> @brief Save particle information in binary format to icbcun
  !! todo: remove now that dedicated track output files are implemented
  !<
  subroutine sav_particles(this, icbcun)
    ! -- dummy
    class(PrtPrpType), intent(inout) :: this
    integer(I4B), intent(in) :: icbcun
    ! -- local
    character(len=16) :: text
    character(len=16), dimension(13) :: auxtxt
    real(DP), dimension(13) :: aux
    integer(I4B) :: nlist
    integer(I4B) :: itrack, irpt, icell
    integer(I4B) :: naux
    !
    if (icbcun /= 0) then
      !
      ! -- Write the header
      text = '      DATA-PRTCL'
      naux = 13
      auxtxt(:) = ['            kper', '            kstp', &
                   '            iprp', '            irpt', &
                   '           icell', '           izone', &
                   '         istatus', '         ireason', &
                   '        trelease', '               t', &
                   '               x', '               y', '               z']
      nlist = this%itrack2 - this%itrack1
      call this%dis%record_srcdst_list_header(text, &
                                              this%name_model, &
                                              this%packName, &
                                              this%name_model, &
                                              this%packName, &
                                              naux, &
                                              auxtxt, &
                                              icbcun, &
                                              nlist, &
                                              this%iout)
      !
      ! -- Write a zero for Q and particle data as aux variables
      do itrack = this%itrack1 + 1, this%itrack2
        irpt = this%trackdata%irpt(itrack)
        icell = this%trackdata%icell(itrack)
        aux(1) = this%trackdata%kper(itrack)
        aux(2) = this%trackdata%kstp(itrack)
        aux(3) = this%trackdata%iprp(itrack) ! todo: as above
        aux(4) = irpt
        aux(5) = icell
        aux(6) = this%trackdata%izone(itrack)
        aux(7) = this%trackdata%istatus(itrack)
        aux(8) = this%trackdata%ireason(itrack)
        aux(9) = this%trackdata%trelease(itrack)
        aux(10) = this%trackdata%t(itrack)
        aux(11) = this%trackdata%x(itrack)
        aux(12) = this%trackdata%y(itrack)
        aux(13) = this%trackdata%z(itrack)
        call this%dis%record_mf6_list_entry(icbcun, irpt, icell, DZERO, &
                                            naux, aux)
      end do
      !
    end if
    !
    ! -- return
    return
  end subroutine sav_particles

end module PrtPrpModule
