module TrackDataModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType
  use UtilMiscModule, only: transform_coords

  implicit none

  private
  public :: TrackDataType

  ! Notes
  ! -----
  !
  ! Each particle's track (or pathline) consists of 1+ records reported
  ! over the model domain while the particle is active.
  !
  ! Each record's istatus attribute indicates whether the particle is
  ! pending release, active, or terminated.
  !
  ! The record's ireason attribute indicates why it was reported. The user
  ! selects 1+ reporting conditions.
  !
  ! Identical particle states may be duplicated if multiple reporting
  ! conditions apply at a given time.
  !
  ! The particle's lifecycle is a strictly increasing step function over
  ! istatus, starting at 0.
  !
  ! There is no particle ID column. Particles can be uniquely identified by
  ! composite key, i.e. combination of column values:
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
  !   *is this necessary?
  !   **unnecessary since PRT makes no distinction between forwards/backwards tracking
  !   ***e.g., released into an inactive cell, a stop zone cell, or a termination zone
  !
  !   ireason: the reason the record was reported
  !   -------
  !     0: release
  !     1: cross spatial boundary (e.g. cell,  subcell)
  !     2: cross temporal boundary (e.g. time step end)
  !     3: exited weak sink
  !     ...

  character(len=*), parameter, public :: TRACKHEADERS = &
                'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,istatus,ireason,&
                &trelease,t,x,y,z'

  character(len=*), parameter, public :: TRACKTYPES = &
                             '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
                             &<f8,<f8,<f8,<f8,<f8'

  type :: TrackDataType
    integer(I4B), pointer :: ibinun => null()
    integer(I4B), pointer :: icsvun => null()
  contains
    procedure, public :: save_record
  end type TrackDataType

contains

  !> @brief Save pathline datum from a particle to a binary or CSV output file.
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
