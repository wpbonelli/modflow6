module TrackDataModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType

  implicit none

  private save_internal
  public :: TrackDataType
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

  type :: TrackDataType
    ! Handles saving particle tracks (pathlines) to binary/CSV output.
    ! Writes to 1 file of each kind, at each level, at a time. At most,
    ! this makes for 4. Kinds are: binary, csv. Levels are: model, prp.

    ! model-level output files
    type(TrackFileType) :: trackfile = TrackFileType()
    type(TrackFileType) :: trackcsvfile = TrackFileType(csv=.true.)

    ! PRP-level output files
    type(TrackFileType) :: prptrackfile = TrackFileType()
    type(TrackFileType) :: prptrackcsvfile = TrackFileType(csv=.true.)
  contains
    procedure, public :: init_track_file
    procedure, public :: save_record
  end type TrackDataType

  ! Data model

  character(len=*), parameter, public :: TRACKHEADERS = &
                'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,istatus,ireason,&
                &trelease,t,x,y,z'

  character(len=*), parameter, public :: TRACKTYPES = &
                             '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
                             &<f8,<f8,<f8,<f8,<f8'

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
  !   * is this necessary?
  !   ** unnecessary since PRT makes no distinction between forwards/backwards tracking
  !   *** e.g., released into an inactive cell, a stop zone cell, or a termination zone
  !
  !   ireason: the reason the record was reported
  !   -------
  !     0: release
  !     1: cross spatial boundary (e.g. cell,  subcell)
  !     2: cross temporal boundary (e.g. time step end)
  !     3: exited weak sink
  !     ...

contains

  !> @brief Initialize the appropriate track file
  subroutine init_track_file(this, iun, csv, iprp)
    ! -- dummy
    class(TrackDataType) :: this
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
        zmodel
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
        zmodel
    end if

  end subroutine

  !> @brief Save the particle's current state to output file(s).
  !!
  !! The record will be saved to any model-level files on the
  !! TrackData instance and any PRP-level files for which the
  !! particle's PRP index matches the track file's PRP index.
  !! This allows for outputting particle tracks for an entire
  !! model, or independently for particular release packages.
  subroutine save_record(this, particle, kper, kstp, reason, level)
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level

    ! -- Only save record if particle is configured for all events or
    !    particle ievent matches provided reason for saving
    if (.not. (particle%ievent == -1 .or. particle%ievent == reason)) return

    ! -- If optional argument level is present and level isn't 3, return early
    !    (kluge note: adds after each subcell-level track)
    if (present(level) .and. level .ne. 3) return 

    ! save binary file(s) if enabled
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

end module TrackDataModule
