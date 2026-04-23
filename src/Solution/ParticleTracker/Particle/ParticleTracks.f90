!> @brief Particle track output module.
!!
!! Each particle's track consists of events reported as the particle is
!! advected through the model domain. Events are snapshots of particle
!! state, along with optional metadata, at a particular moment in time.
!!
!! Particles have no ID property. A particle can be uniquely identified
!! by the unique combination of its release attributes (model, package,
!! position, and time). This is possible because only one particle may
!! be released from a given point at a given time.
!!
!! This module consumes particle events and is responsible for writing
!! them to one or more track files, binary or CSV, and for logging the
!! events if requested. Each track file is associated with either a PRP
!! package or with the full PRT model (there may only be 1 such latter).
!!
!! Events are buffered via a pluggable strategy (memory or scratch file)
!! and flushed to disk only when a time step converges. See TrackBufferType,
!! MemoryTrackBufferType, and ScratchTrackBufferType.
!<
module ParticleTracksModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, DONE, DPIO180
  use ParticleModule, only: ParticleType, ACTIVE
  use ParticleEventModule, only: ParticleEventType
  use ReleaseEventModule, only: ReleaseEventType
  use TimeStepEventModule, only: TimeStepEventType
  use TerminationEventModule, only: TerminationEventType
  use WeakSinkEventModule, only: WeakSinkEventType
  use UserTimeEventModule, only: UserTimeEventType
  use CellExitEventModule, only: CellExitEventType
  use SubcellExitEventModule, only: SubcellExitEventType
  use DroppedEventModule, only: DroppedEventType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform

  implicit none
  public :: ParticleTrackFileType, &
            ParticleTracksType, &
            ParticleTrackEventSelectionType, &
            add_particle_event
  private :: save_record

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

  !> @brief Flat record of a particle track event.
  type :: TrackRecordType
    integer(I4B) :: kper, kstp, imdl, iprp, irpt, ilay, icu, izone
    integer(I4B) :: istatus, ireason
    real(DP) :: trelease, ttrack, x, y, z
    character(len=40) :: name
  end type TrackRecordType

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! Can be associated with a particle release point (PRP) package
  !! or with an entire model, and can be binary or text (CSV).
  !<
  type :: ParticleTrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type ParticleTrackFileType

  !> @brief Selection of particle events.
  type :: ParticleTrackEventSelectionType
    logical(LGP) :: release !< track release events
    logical(LGP) :: featexit !< track grid feature exits
    logical(LGP) :: timestep !< track timestep ends
    logical(LGP) :: terminate !< track termination events
    logical(LGP) :: weaksink !< track weak sink exit events
    logical(LGP) :: usertime !< track user-selected times
    logical(LGP) :: subfexit !< track subfeature exits
    logical(LGP) :: dropped !< track water table drops
  end type ParticleTrackEventSelectionType

  !> @brief Abstract track event buffer strategy.
  !!
  !! Buffers events between tracking runs so only the last successful
  !! run's events are written to disk. Implementations differ in where
  !! they stage the records: memory (MemoryTrackBufferType) or a
  !! temporary scratch file (ScratchTrackBufferType).
  !<
  type, abstract :: TrackBufferType
  contains
    procedure :: init => buffer_init_noop !< open/allocate resources (no-op default)
    procedure(buffer_append_if), deferred :: append !< buffer one record
    procedure(buffer_flush_if), deferred :: flush !< write buffered records, reset
    procedure(buffer_simple_if), deferred :: discard !< reset without writing
    procedure(buffer_simple_if), deferred :: destroy !< release resources
  end type TrackBufferType

  !> @brief In-memory track event buffer. Records are held in a
  !! dynamically growing array that doubles in capacity as needed.
  !<
  type, extends(TrackBufferType) :: MemoryTrackBufferType
    type(TrackRecordType), allocatable :: records(:) !< buffered records
    integer(I4B) :: nrecords = 0 !< number of records in buffer
  contains
    procedure :: append => memory_append
    procedure :: flush => memory_flush
    procedure :: discard => memory_discard
    procedure :: destroy => memory_destroy
  end type MemoryTrackBufferType

  !> @brief Scratch-file track event buffer. Records are written to
  !! an unformatted sequential scratch file that is rewound on discard
  !! and read back sequentially on flush. The OS page cache serves
  !! most I/O; physical disk I/O occurs only under memory pressure.
  !<
  type, extends(TrackBufferType) :: ScratchTrackBufferType
    integer(I4B) :: iun = 0 !< scratch file unit
    integer(I4B) :: nrecords = 0 !< records written since last discard/flush
  contains
    procedure :: init => scratch_init
    procedure :: append => scratch_append
    procedure :: flush => scratch_flush
    procedure :: discard => scratch_discard
    procedure :: destroy => scratch_destroy
  end type ScratchTrackBufferType

  !> @brief Particle track output manager. Handles printing as well as writing
  !! to files. One output unit can be configured for printing. Multiple files
  !! can be configured for writing, with each file optionally associated with
  !! a PRP package or with the full model. Events can be filtered by type, so
  !! that only certain event types are printed or written to files.
  !!
  !! Records are buffered via the configured strategy (memory or scratch file)
  !! and flushed to disk only when the time step converges. init_buffer must
  !! be called before any tracking occurs.
  !<
  type :: ParticleTracksType
    private
    integer(I4B), public :: iout = -1 !< log file unit
    integer(I4B), public :: ntrackfiles !< number of track files
    type(ParticleTrackFileType), public, allocatable :: files(:) !< track files
    type(ParticleTrackEventSelectionType), public :: selected !< event selection
    class(TrackBufferType), allocatable :: buffer !< buffer strategy
  contains
    procedure, public :: init_file
    procedure, public :: init_buffer
    procedure, public :: is_selected
    procedure, public :: select_events
    procedure, public :: buffer_event
    procedure, public :: flush_buffer
    procedure, public :: discard_buffer
    procedure, public :: destroy
    procedure :: expand_files
    procedure :: should_save
    procedure :: should_print
  end type ParticleTracksType

  abstract interface
    subroutine buffer_append_if(this, rec)
      import TrackBufferType, TrackRecordType
      class(TrackBufferType) :: this
      type(TrackRecordType), intent(in) :: rec
    end subroutine buffer_append_if

    subroutine buffer_flush_if(this, files, nfiles)
      import TrackBufferType, ParticleTrackFileType, I4B
      class(TrackBufferType) :: this
      type(ParticleTrackFileType), intent(in) :: files(:)
      integer(I4B), intent(in) :: nfiles
    end subroutine buffer_flush_if

    ! Shared interface for discard and destroy: both take only `this`.
    subroutine buffer_simple_if(this)
      import TrackBufferType
      class(TrackBufferType) :: this
    end subroutine buffer_simple_if
  end interface

contains

  !> @brief Initialize a binary or CSV output file.
  subroutine init_file(this, iun, csv, iprp)
    class(ParticleTracksType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    type(ParticleTrackFileType), pointer :: file

    if (.not. allocated(this%files)) then
      allocate (this%files(1))
    else
      call this%expand_files(increment=1)
    end if

    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp
    this%ntrackfiles = size(this%files)
    this%files(this%ntrackfiles) = file
  end subroutine init_file

  !> @brief Initialize the event buffer strategy.
  !! If use_scratch is .true., events are staged in a scratch file;
  !! otherwise they are held in memory. Must be called before tracking.
  subroutine init_buffer(this, use_scratch)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: use_scratch

    if (use_scratch) then
      allocate (ScratchTrackBufferType :: this%buffer)
    else
      allocate (MemoryTrackBufferType :: this%buffer)
    end if
    call this%buffer%init()
  end subroutine init_buffer

  !> @brief Destroy the particle track manager.
  subroutine destroy(this)
    class(ParticleTracksType) :: this
    call this%buffer%destroy()
    if (allocated(this%buffer)) deallocate (this%buffer)
    if (allocated(this%files)) deallocate (this%files)
  end subroutine destroy

  !> @brief Grow the array of track files.
  subroutine expand_files(this, increment)
    class(ParticleTracksType) :: this
    integer(I4B), optional, intent(in) :: increment
    integer(I4B) :: inclocal, isize, newsize
    type(ParticleTrackFileType), allocatable, dimension(:) :: temp

    if (present(increment)) then
      inclocal = increment
    else
      inclocal = 1
    end if

    if (allocated(this%files)) then
      isize = size(this%files)
      newsize = isize + inclocal
      allocate (temp(newsize))
      temp(1:isize) = this%files
      deallocate (this%files)
      call move_alloc(temp, this%files)
    else
      allocate (this%files(inclocal))
    end if
  end subroutine expand_files

  !> @brief Pick events to track.
  subroutine select_events(this, &
                           release, &
                           featexit, &
                           timestep, &
                           terminate, &
                           weaksink, &
                           usertime, &
                           subfexit, &
                           dropped)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: release
    logical(LGP), intent(in) :: featexit
    logical(LGP), intent(in) :: timestep
    logical(LGP), intent(in) :: terminate
    logical(LGP), intent(in) :: weaksink
    logical(LGP), intent(in) :: usertime
    logical(LGP), intent(in) :: subfexit
    logical(LGP), intent(in) :: dropped
    this%selected%release = release
    this%selected%featexit = featexit
    this%selected%timestep = timestep
    this%selected%terminate = terminate
    this%selected%weaksink = weaksink
    this%selected%usertime = usertime
    this%selected%subfexit = subfexit
    this%selected%dropped = dropped
  end subroutine select_events

  !> @brief Check if a given event code is selected for tracking.
  logical function is_selected(this, event) result(selected)
    class(ParticleTracksType), intent(inout) :: this
    class(ParticleEventType), intent(in) :: event

    select type (event)
    type is (ReleaseEventType)
      selected = this%selected%release
    type is (CellExitEventType)
      selected = this%selected%featexit
    type is (SubcellExitEventType)
      selected = this%selected%subfexit
    type is (TimeStepEventType)
      selected = this%selected%timestep
    type is (TerminationEventType)
      selected = this%selected%terminate
    type is (WeakSinkEventType)
      selected = this%selected%weaksink
    type is (UserTimeEventType)
      selected = this%selected%usertime
    type is (DroppedEventType)
      selected = this%selected%dropped
    class default
      call pstop(1, "unknown event type")
      selected = .false.
    end select
  end function is_selected

  !> @brief Buffer a particle event for deferred write.
  !! The buffer can be flushed to disk by flush_buffer,
  !! or discarded by discard_buffer (e.g., on a failed
  !! ATS retry or Picard iteration).
  subroutine buffer_event(this, particle, event)
    class(ParticleTracksType) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    type(TrackRecordType) :: rec

    rec%kper = event%kper; rec%kstp = event%kstp
    rec%imdl = event%imdl; rec%iprp = event%iprp
    rec%irpt = event%irpt; rec%ilay = event%ilay
    rec%icu = event%icu; rec%izone = event%izone
    rec%istatus = event%istatus; rec%ireason = event%get_code()
    rec%trelease = event%trelease; rec%ttrack = event%ttrack
    rec%x = event%x; rec%y = event%y; rec%z = event%z
    rec%name = trim(adjustl(particle%name))

    call this%buffer%append(rec)
  end subroutine buffer_event

  !> @brief Flush the event buffer to disk.
  !! Called from prt_ot after a time step converges.
  subroutine flush_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%flush(this%files, this%ntrackfiles)
  end subroutine flush_buffer

  !> @brief Discard buffered events without writing.
  !! Called from prt_ad at the start of each tracking attempt,
  !! so events from failed or superseded attempts are never written.
  subroutine discard_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%discard()
  end subroutine discard_buffer

  !> @brief Check whether a particle belongs in a given file.
  logical function should_save(this, particle, file) result(save)
    class(ParticleTracksType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    type(ParticleTrackFileType), intent(in) :: file
    save = (file%iun > 0 .and. &
            (file%iprp == -1 .or. file%iprp == particle%iprp))
  end function should_save

  !> @brief Is the output unit valid?
  logical function should_print(this)
    class(ParticleTracksType), intent(inout) :: this
    should_print = this%iout >= 0
  end function should_print

  subroutine buffer_init_noop(this)
    class(TrackBufferType) :: this
  end subroutine buffer_init_noop

  subroutine memory_append(this, rec)
    class(MemoryTrackBufferType) :: this
    type(TrackRecordType), intent(in) :: rec
    type(TrackRecordType), allocatable :: tmp(:)

    if (.not. allocated(this%records)) then
      allocate (this%records(64))
    else if (this%nrecords == size(this%records)) then
      allocate (tmp(size(this%records) * 2))
      tmp(1:this%nrecords) = this%records
      call move_alloc(tmp, this%records)
    end if
    this%nrecords = this%nrecords + 1
    this%records(this%nrecords) = rec
  end subroutine memory_append

  subroutine memory_flush(this, files, nfiles)
    class(MemoryTrackBufferType) :: this
    type(ParticleTrackFileType), intent(in) :: files(:)
    integer(I4B), intent(in) :: nfiles
    integer(I4B) :: n, i
    type(TrackRecordType) :: rec
    type(ParticleTrackFileType) :: file

    do n = 1, this%nrecords
      rec = this%records(n)
      do i = 1, nfiles
        file = files(i)
        if (file%iun > 0 .and. &
            (file%iprp == -1 .or. file%iprp == rec%iprp)) &
          call save_record(file%iun, rec, file%csv)
      end do
    end do
    this%nrecords = 0
  end subroutine memory_flush

  subroutine memory_discard(this)
    class(MemoryTrackBufferType) :: this
    this%nrecords = 0
  end subroutine memory_discard

  subroutine memory_destroy(this)
    class(MemoryTrackBufferType) :: this
    if (allocated(this%records)) deallocate (this%records)
  end subroutine memory_destroy

  subroutine scratch_init(this)
    class(ScratchTrackBufferType) :: this
    integer(I4B) :: istat
    open (newunit=this%iun, status='scratch', form='unformatted', &
          access='sequential', iostat=istat)
    if (istat /= 0) &
      call pstop(1, 'failed to open scratch track buffer file')
  end subroutine scratch_init

  subroutine scratch_append(this, rec)
    class(ScratchTrackBufferType) :: this
    type(TrackRecordType), intent(in) :: rec
    write (this%iun) rec
    this%nrecords = this%nrecords + 1
  end subroutine scratch_append

  subroutine scratch_flush(this, files, nfiles)
    class(ScratchTrackBufferType) :: this
    type(ParticleTrackFileType), intent(in) :: files(:)
    integer(I4B), intent(in) :: nfiles
    integer(I4B) :: n, i
    type(TrackRecordType) :: rec
    type(ParticleTrackFileType) :: file

    rewind (this%iun)
    do n = 1, this%nrecords
      read (this%iun) rec
      do i = 1, nfiles
        file = files(i)
        if (file%iun > 0 .and. &
            (file%iprp == -1 .or. file%iprp == rec%iprp)) &
          call save_record(file%iun, rec, file%csv)
      end do
    end do
    rewind (this%iun)
    this%nrecords = 0
  end subroutine scratch_flush

  subroutine scratch_discard(this)
    class(ScratchTrackBufferType) :: this
    rewind (this%iun)
    this%nrecords = 0
  end subroutine scratch_discard

  subroutine scratch_destroy(this)
    class(ScratchTrackBufferType) :: this
    if (this%iun > 0) then
      close (this%iun)
      this%iun = 0
    end if
    this%nrecords = 0
  end subroutine scratch_destroy

  !> @brief Save an event record to a binary or CSV file.
  subroutine save_record(iun, rec, csv)
    integer(I4B), intent(in) :: iun
    type(TrackRecordType), intent(in) :: rec
    logical(LGP), intent(in) :: csv

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, trim(rec%name)
    else
      write (iun) &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, rec%name
    end if
  end subroutine save_record

  !> @brief Add a particle event to be written to eligible
  !! files and printed to an output file unit if requested.
  !! This function should be subscribed as an event handler
  !! to particle event dispatchers. Events are buffered, to
  !! be written to output files upon successful completion
  !! of a time step when the framework OT hook is executed.
  function add_particle_event(context, particle, event) result(handled)
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: handled

    select type (context)
    type is (ParticleTracksType)
      if (context%should_print()) &
        call event%log(context%iout)
      if (context%is_selected(event)) &
        call context%buffer_event(particle, event)
      handled = .true.
    end select
  end function add_particle_event

end module ParticleTracksModule
