module TrackDataModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType
  use UtilMiscModule, only: transform_coords

  implicit none

  private
  public :: TrackDataType

  !> @brief Handles saving pathlines to binary and/or CSV output files.
  type :: TrackDataType
    integer(I4B), pointer :: ibinun => null()
    integer(I4B), pointer :: icsvun => null()
  contains
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

  !> @brief Save particle pathline record to a binary or CSV file.
  subroutine save_record(this, particle, kper, kstp, reason, level)
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: kper, kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level
    ! -- local
    logical(LGP) :: ladd
    real(DP) :: xmodel, ymodel, zmodel
    integer(I4B) :: status

    ! -- Determine whether to add track data
    if (.not. present(level)) then
      ! -- If optional argument level is not present, track data will be added
      ! -- by default
      ladd = .true.
    else
      ! If optional argument level is present, check criteria
      ladd = .false.
      if (level == 3) ladd = .true. ! kluge note: adds after each subcell-level track
    end if

    if (ladd) then

      ! -- Get model coordinates
      call particle%get_model_coords(xmodel, ymodel, zmodel)

      ! -- Get status
      if (particle%istatus .lt. 0) then
        status = 1
      else
        status = particle%istatus
      end if

      ! -- Save to CSV file
      if (this%icsvun /= 0) &
        write (this%icsvun, '(*(G0,:,","))') &
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

      ! -- Save to binary file
      if (this%ibinun /= 0) &
        write (this%ibinun) &
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
    return
  end subroutine save_record

end module TrackDataModule
