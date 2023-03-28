module TrackModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType

  implicit none

  private save_record
  public :: TrackControlType
  public :: TrackFileType

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! May be associated with a particle release point (PRP) package
  !! or with an entire model.
  !<
  type :: TrackFileType
    integer(I4B) :: iun = 0 !< file unit number
    logical(LGP) :: csv = .false. !< whether the file is binary or CSV
    integer(I4B) :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type TrackFileType

  !> @brief Manages particle track (i.e. pathline) files.
  !!
  !! Optionally filters events ("ireason" codes), e.g. release
  !! or cell-to-cell transitions. Events ("ireason" codes) are:
  !!
  !!   -1: ALL
  !!    0: RELEASE
  !!    1: TRANSIT
  !!    2: TIMESTEP
  !!    3: TERMINATE
  !!    4: WEAKSINK
  !!
  !! An arbitrary number of files can be managed, internal
  !! arrays are resized as needed.
  !<
  type :: TrackControlType
    private
    type(TrackFileType), public, allocatable :: trackfiles(:) !< output files
    integer(I4B), public :: ntrackfiles !< number of output files
    integer(I4B), public :: itrackevent !< track event selection
  contains
    procedure :: expand
    procedure, public :: init_track_file
    procedure, public :: save
    procedure, public :: set_track_event
  end type TrackControlType

  ! Track file headers
  character(len=*), parameter, public :: TRACKHEADERS = &
                'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,istatus,ireason,&
                &trelease,t,x,y,z,name'

  ! Track file dtypes
  character(len=*), parameter, public :: TRACKTYPES = &
                             '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
                             &<f8,<f8,<f8,<f8,<f8,|S40'

  ! Notes
  ! -----
  !
  ! Each particle's pathline consists of 1+ records reported as the particle
  ! is tracked over the model domain. Records are snapshots of the particle's
  ! state (e.g. tracking status, position) at a particular moment in time.
  !
  ! Particles have no ID property. Particles can be uniquely identified
  ! by composite key, i.e. combination of:
  !
  !   - imdl: originating model ID
  !   - iprp: originating PRP ID
  !   - irpt: particle release location ID
  !   - trelease: particle release time
  !
  ! Each record has an "ireason" property, which identifies the cause of
  ! the record. The user selects 1+ conditions or events for recording.
  ! Identical records (except "ireason") may be duplicated if multiple
  ! reporting conditions apply to particles at the same moment in time.
  ! Each "ireason" value corresponds to an OC "trackevent" option value:
  !
  !     0: particle released
  !     1: particle transitioned between cells
  !     2: current time step ended****
  !     3: particle terminated
  !     4: particle exited weak sink
  !
  ! Each record has an "istatus" property, which is the tracking status;
  ! e.g., awaiting release, active, terminated. A particle may terminate
  ! for several reasons. Status values greater than one imply termination.
  ! Particle status strictly increases over time, starting at zero:
  !
  !     0: pending release*
  !     1: active
  !     2: terminated at boundary face
  !     3: terminated in weak sink cell
  !     4: terminated in weak source cell**
  !     5: terminated in cell with no exit face
  !     6: terminated in cell with specified zone number
  !     7: terminated in inactive cell
  !     8: permanently unreleased***
  !     9: terminated for unknown reason*
  !
  ! PRT shares the same status enumeration as MODPATH 7. However, some
  ! don't apply to PRT; for instance, MODPATH 7 distinguishes forwards
  ! and backwards tracking, but status value 4 is not used by PRT.
  !
  ! --------------------------------------------------
  !   * is this necessary?
  !   ** unnecessary since PRT makes no distinction between forwards/backwards tracking
  !   *** e.g., released into an inactive cell, a stop zone cell, or a termination zone
  !   **** this may coincide with termination, in which case two events are reported

contains

  !> @brief Initialize a new track file
  subroutine init_track_file(this, iun, csv, iprp)
    ! -- dummy
    class(TrackControlType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    ! -- local
    type(TrackFileType), pointer :: file

    ! -- allocate or expand array
    if (.not. allocated(this%trackfiles)) then
      allocate (this%trackfiles(1))
    else
      call this%expand(increment=1)
    end if

    ! -- setup new file
    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp

    ! -- update array and counter
    this%ntrackfiles = size(this%trackfiles)
    this%trackfiles(this%ntrackfiles) = file

  end subroutine init_track_file

  !> @brief Expand the trackfile array, internal use only
  subroutine expand(this, increment)
    ! -- dummy
    class(TrackControlType) :: this
    integer(I4B), optional, intent(in) :: increment
    ! -- local
    integer(I4B) :: inclocal, isize, newsize
    type(TrackFileType), allocatable, dimension(:) :: temp

    ! -- initialize
    if (present(increment)) then
      inclocal = increment
    else
      inclocal = 1
    end if

    ! -- increase size of array by inclocal
    if (allocated(this%trackfiles)) then
      isize = size(this%trackfiles)
      newsize = isize + inclocal
      allocate (temp(newsize))
      temp(1:isize) = this%trackfiles
      deallocate (this%trackfiles)
      call move_alloc(temp, this%trackfiles)
    else
      allocate (this%trackfiles(inclocal))
    end if

  end subroutine expand

  !> @brief Save record to binary or CSV file, internal use only
  subroutine save_record(iun, particle, kper, kstp, reason, csv)
    ! -- dummy
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    logical(LGP), intent(in) :: csv
    ! -- local
    real(DP) :: x, y, z
    integer(I4B) :: status

    ! -- Get model (global) coordinates
    call particle%get_model_coords(x, y, z)

    ! -- Get status
    if (particle%istatus .lt. 0) then
      status = 1
    else
      status = particle%istatus
    end if

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        kper, &
        kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        reason, &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        trim(adjustl(particle%name))
    else
      write (iun) &
        kper, &
        kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        reason, &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        particle%name
    end if

  end subroutine

  !> @brief Save the particle's state to track output file(s).
  !!
  !! A record is saved to all enabled model-level files and to
  !! any PRP-level files with PRP index matching the particle's
  !! PRP index.
  !<
  subroutine save(this, particle, kper, kstp, reason, level)
    ! -- dummy
    class(TrackControlType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level
    ! -- local
    integer(I4B) :: i
    type(TrackFileType) :: file

    ! -- Only save if reporting is enabled for all events, or
    !    if the specified event matches the reporting reason.
    if (this%itrackevent > -1 .and. this%itrackevent /= reason) &
      return

    ! -- For now, only allow reporting from outside the tracking
    !    algorithm (e.g. release time), in which case level will
    !    not be provided, or if within the tracking solution, in
    !    subcells (level 3) only. This may change if the subcell
    !    ever needs to delegate to an even finer granularity, in
    !    which case the tracking solution would recurse 1+ calls
    !    deeper before advancing the particle and unwinding.
    if (present(level)) then
      if (level .ne. 3) return
    end if

    ! -- Save all enabled model-level files and any
    ! -- enabled and index-matching PRP-level files
    do i = 1, this%ntrackfiles
      file = this%trackfiles(i)
      if (file%iun > 0 .and. &
          (file%iprp == -1 .or. &
           file%iprp == particle%iprp)) &
        call save_record(file%iun, particle, &
                         kper, kstp, reason, csv=file%csv)
    end do

  end subroutine save

  !> @brief Set the event filter.
  !!
  !! Track event -1 indicates TRACKEVENT ALL and so on.
  !! If track event >= 0, only events of the given type
  !! will be recorded. Each non-negative tracking event
  !! number corresponds to an "ireason" code as appears
  !! in each row of output.
  !<
  subroutine set_track_event(this, itrackevent)
    class(TrackControlType) :: this
    integer(I4B), intent(in) :: itrackevent
    this%itrackevent = itrackevent
  end subroutine set_track_event

end module TrackModule
