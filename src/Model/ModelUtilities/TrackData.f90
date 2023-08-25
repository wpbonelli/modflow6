module TrackModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType

  implicit none

  private save_internal
  public :: TrackControlType
  public :: TrackFileType

  ! Types

  type :: TrackFileType
    ! Describes an output file containing all or some particle track data.
    ! May be associated with a particular particle release package (PRP).
    ! File is for entire model's pathlines if not associated with a PRP,
    ! otherwise tracks only particles in the particle release package.
    integer(I4B) :: iun = 0
    logical(LGP) :: csv = .false.
    integer(I4B) :: iprp = -1 ! prp 0 is reserved for exchange PRP
  end type TrackFileType

  type :: TrackControlType
    ! Manages file output for particle tracks (i.e. pathlines).
    ! Writes particle pathlines to file, optionally filtering by
    ! events of a given category, e.g. cell-to-cell transitions.
    ! Writes to 1 file of each kind, at each level, at a time. At most,
    ! this makes for 4. Kinds are: binary, csv. Levels are: model, prp.

    ! -1: ALL
    !  0: RELEASE
    !  1: TRANSIT
    !  2: TIMESTEP
    !  3: TERMINATE
    !  4: WEAKSINK
    integer(I4B), public :: itrackevent

    ! model-level output files
    type(TrackFileType), public :: trackfile = TrackFileType()
    type(TrackFileType), public :: trackcsvfile = TrackFileType(csv=.true.)
    !
    ! PRP-level output files
    type(TrackFileType), public :: prptrackfile = TrackFileType()
    type(TrackFileType), public :: prptrackcsvfile = TrackFileType(csv=.true.)

  contains
    procedure, public :: init_track_file
    procedure, public :: save_record
    procedure, public :: set_track_event
  end type TrackControlType

  ! Data model

  character(len=*), parameter, public :: TRACKHEADERS = &
                'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,istatus,ireason,&
                &trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKTYPES = &
                             '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
                             &<f8,<f8,<f8,<f8,<f8,|S40'

  ! Notes
  ! -----
  !
  ! Each particle's pathline consists of 1+ records reported as the particle
  ! is tracked over the model domain. Each record represents the particle's
  ! state (e.g. tracking status, position) at a particular moment in time.
  !
  ! Each record has an istatus property indicating the particle status, e.g.
  ! awaiting release, active, or terminated. A particle may terminate for
  ! several reasons. Any istatus value greater than 1 implies termination.
  ! The particle's lifecycle is a strictly increasing step function over
  ! istatus, starting at 0. For easier porting of models, PRT employs the
  ! same status enumeration as MODPATH 7 endpoint files. Some categories
  ! are not applicable to PRT, however, and will not appear in PRT output:
  ! for instance, MODPATH 7 distinguishes between forwards and backwards
  ! tracking, while PRT does not, so istatus value 4 is not used by PRT.
  !
  ! Each record has an ireason property identifying the reporting condition
  ! responsible for the record. The user selects 1+ reporting conditions.
  ! Identical records (excepting ireason) may be duplicated if multiple
  ! reporting conditions apply to particles at the same moment in time.
  ! Each ireason value corresponds to a PRT OC trackevent option value.
  !
  ! Particles have no ID property. Particles can be uniquely identified by
  ! composite key, i.e. combination of properties:
  !
  !   - imdl: originating model ID
  !   - iprp: originating PRP ID
  !   - irpt: particle release location ID
  !   - trelease: particle release time
  !
  ! Enumerations
  ! ------------
  !
  !   istatus: the particle's status at time of record
  !   -------
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
  !   ireason: the reason the record was reported
  !   -------
  !     0: particle released
  !     1: particle transitioned between cells
  !     2: current time step ended****
  !     3: particle terminated
  !     4: particle exited weak sink
  !     ...
  !
  !   * is this necessary?
  !   ** unnecessary since PRT makes no distinction between forwards/backwards tracking
  !   *** e.g., released into an inactive cell, a stop zone cell, or a termination zone
  !   **** this may coincide with termination, in which case two events are reported

contains

  !> @brief Initialize the appropriate track file
  subroutine init_track_file(this, iun, csv, iprp)
    ! -- dummy
    class(TrackControlType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    ! -- local
    logical(LGP) :: lcsv
    integer(I4B) :: lprp

    ! parse csv argument
    if (.not. present(csv)) then
      lcsv = .false.
    else
      lcsv = csv
    end if

    ! parse prp argument
    if (.not. present(iprp)) then
      lprp = -1
    else
      lprp = iprp
    end if

    ! setup appropriate track file
    if (lprp > -1) then
      if (lcsv) then
        this%prptrackcsvfile%iprp = lprp
        this%prptrackcsvfile%iun = iun
      else
        this%prptrackfile%iprp = lprp
        this%prptrackfile%iun = iun
      end if
    else
      if (lcsv) then
        this%trackcsvfile%iun = iun
      else
        this%trackfile%iun = iun
      end if
    end if
  end subroutine init_track_file

  subroutine save_internal(iun, particle, kper, kstp, reason, level, csv)
    ! -- dummy
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level
    logical(LGP), intent(in) :: csv
    ! -- local
    real(DP) :: xmodel, ymodel, zmodel
    integer(I4B) :: status

    ! -- Get model coordinates
    call particle%get_model_coords(xmodel, ymodel, zmodel)

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
        xmodel, &
        ymodel, &
        zmodel, &
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
        xmodel, &
        ymodel, &
        zmodel, &
        particle%name
    end if
  end subroutine

  !> @brief Save the particle's current state to output file(s).
  !!
  !! The record will be saved to both model-level files and
  !! any PRP-level files for which the particle's PRP index
  !! matches the track file's PRP index. This allows saving
  !! tracks for an entire model or for a particular package.
  subroutine save_record(this, particle, kper, kstp, reason, level)
    ! -- dummy
    class(TrackControlType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level

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
    if (present(level) .and. level .ne. 3) return

    ! save model-level file(s) if enabled
    if (this%trackfile%iun > 0) &
      call save_internal(this%trackfile%iun, particle, &
                         kper, kstp, reason, level, csv=.false.)
    if (this%trackcsvfile%iun > 0) &
      call save_internal(this%trackcsvfile%iun, particle, &
                         kper, kstp, reason, level, csv=.true.)

    ! save PRP-level file(s) if enabled
    if (particle%iprp == this%prptrackfile%iprp) then
      if (this%prptrackfile%iun > 0) &
        call save_internal(this%prptrackfile%iun, particle, &
                           kper, kstp, reason, level, csv=.false.)
      if (this%prptrackcsvfile%iun > 0) &
        call save_internal(this%prptrackcsvfile%iun, particle, &
                           kper, kstp, reason, level, csv=.true.)
    end if
  end subroutine save_record

  !> @brief Set the event to be reported.
  !!
  !! itrackevent -1 corresponds to TRACKEVENT ALL
  !! and so on. If itrackevent >= 0, the selected
  !! event type will be saved and others omitted.
  subroutine set_track_event(this, itrackevent)
    class(TrackControlType) :: this
    integer(I4B), intent(in) :: itrackevent
    this%itrackevent = itrackevent
  end subroutine set_track_event

end module TrackModule
