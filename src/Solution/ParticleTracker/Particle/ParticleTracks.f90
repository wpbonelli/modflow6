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
!! Events are buffered in memory or in a scratch file and are flushed to
!! disk only when a time step successfully finishes. All events are
!! buffered unconditionally; each file's own event selection and scope
!! filter are applied at flush time.
!<
module ParticleTracksModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType, ACTIVE
  use ParticleEventModule, only: ParticleEventType
  use ReleaseEventModule, only: ReleaseEventType
  use TimeStepEventModule, only: TimeStepEventType
  use TerminationEventModule, only: TerminationEventType
  use WeakSinkEventModule, only: WeakSinkEventType
  use UserTimeEventModule, only: UserTimeEventType
  use CellExitEventModule, only: CellExitEventType
  use SubcellExitEventModule, only: SubCellExitEventType
  use DroppedEventModule, only: DroppedEventType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform
  use ParticleTrackEventBufferModule
  use MemoryBufferModule, only: MemoryBufferType
  use ScratchFileBufferModule, only: ScratchFileBufferType

  implicit none
  public :: ParticleTracksType, &
            ParticleTrackEventSelectionType, &
            add_particle_event

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

  !> @brief Particle track output manager. Handles printing as well as writing
  !! to files. One output unit can be configured for printing. Multiple files
  !! can be configured for writing, with each file optionally associated with
  !! a PRP package or with the full model. Each file carries its own event
  !! selection. All events are buffered; per-file filtering happens at flush.
  !<
  type :: ParticleTracksType
    private
    integer(I4B), public :: iout = -1 !< log file unit
    integer(I4B), public :: ntrackfiles !< number of track files
    type(ParticleTrackFileType), public, allocatable :: files(:) !< track files
    class(ParticleTrackEventBufferType), allocatable :: buffer !< event buffer
  contains
    procedure, public :: init_file
    procedure, public :: init_buffer
    procedure, public :: buffer_event
    procedure, public :: flush_buffer
    procedure, public :: discard_buffer
    procedure, public :: destroy
    procedure :: expand_files
  end type ParticleTracksType

contains

  !> @brief Initialize a binary or CSV output file.
  subroutine init_file(this, iun, csv, iprp, selected)
    class(ParticleTracksType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    type(ParticleTrackEventSelectionType), intent(in), optional :: selected
    type(ParticleTrackFileType), pointer :: file

    call this%expand_files(increment=1)

    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp
    if (present(selected)) file%selected = selected
    this%ntrackfiles = size(this%files)
    this%files(this%ntrackfiles) = file
  end subroutine init_file

  !> @brief Initialize the event buffer strategy
  subroutine init_buffer(this, scratch)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: scratch

    if (scratch) then
      allocate (ScratchFileBufferType :: this%buffer)
    else
      allocate (MemoryBufferType :: this%buffer)
    end if
    call this%buffer%init()
  end subroutine init_buffer

  !> @brief Destroy the particle track manager.
  subroutine destroy(this)
    class(ParticleTracksType) :: this
    call this%buffer%destroy()
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

  !> @brief Buffer an event for deferred write.
  subroutine buffer_event(this, particle, event)
    class(ParticleTracksType) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    type(ParticleTrackRecordType) :: rec

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
  subroutine flush_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%flush(this%files)
  end subroutine flush_buffer

  !> @brief Discard buffered events without writing.
  subroutine discard_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%discard()
  end subroutine discard_buffer

  !> @brief Add a particle event to be written to eligible
  !! files and printed to an output file unit if requested.
  !! This function should be subscribed as an event handler
  !! to particle event dispatchers. Events are buffered to
  !! be written to output files upon successful completion
  !! of a time step when the framework OT hook is executed.
  !! All events are buffered; per-file filtering is deferred
  !! to flush time via each file's wants_record method.
  function add_particle_event(context, particle, event) result(handled)
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: handled

    select type (context)
    type is (ParticleTracksType)
      if (context%iout >= 0) &
        call event%log(context%iout)
      call context%buffer_event(particle, event)
      handled = .true.
    end select
  end function add_particle_event

end module ParticleTracksModule
