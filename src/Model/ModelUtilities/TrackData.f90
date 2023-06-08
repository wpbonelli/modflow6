module TrackDataModule

  use KindModule, only: DP, I4B, LGP

  implicit none

  private
  public :: TrackDataType

  character(len=*), parameter, public :: TRACKHEADERS = &
                'kper,kstp,iprp,irpt,icell,izone,istatus,ireason,&
                &trelease,t,x,y,z'

  character(len=*), parameter, public :: TRACKTYPES = &
                             '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
                             &<f8,<f8,<f8,<f8,<f8'

  ! Structure of arrays to hold particle tracks. Arrays are in long format.
  ! Each particle's track across the simulation domain consists of 1+ rows.
  ! Particles can be uniquely identified by a combination of column values:
  !   - todo imdl: originating model ID
  !   - iprp: originating PRP ID
  !   - irpt: particle release location ID
  !   - trelease: particle release time (retrieved from partlist)
  type :: TrackDataType
    ! integer arrays
    integer(I4B), pointer :: nrows => null() ! total count of track data
    integer(I4B), dimension(:), pointer, contiguous :: kper ! stress period
    integer(I4B), dimension(:), pointer, contiguous :: kstp ! time step
    integer(I4B), dimension(:), pointer, contiguous :: irpt ! particle ID
    integer(I4B), dimension(:), pointer, contiguous :: iprp ! originating PRP ID
    integer(I4B), dimension(:), pointer, contiguous :: icell ! cell ID
    integer(I4B), dimension(:), pointer, contiguous :: izone ! todo zone number
    integer(I4B), dimension(:), pointer, contiguous :: istatus ! particle status
    integer(I4B), dimension(:), pointer, contiguous :: ireason ! reason for datum
    ! ireason can take values:
    !   0: release
    !   1: cross cell boundary (or subcell? or generic feature? worth distinguishing?)
    !   2: todo time series
    !   3: termination

    ! double arrays
    real(DP), dimension(:), pointer, contiguous :: trelease ! particle's release time
    real(DP), dimension(:), pointer, contiguous :: t ! current time
    real(DP), dimension(:), pointer, contiguous :: x ! current x coordinate
    real(DP), dimension(:), pointer, contiguous :: y ! current y coordinate
    real(DP), dimension(:), pointer, contiguous :: z ! current z coordinate

  contains
    procedure, public :: save_track_data
  end type TrackDataType

contains

  !> @brief Write track data to a binary or CSV output file.
  !!
  !! Arguments itrack1 and itrack2 may be provided to select a subset of
  !! track data to write to file. This can be used to write data for one
  !! or multiple contiguous PRPs, instead of all particles in the model.
  !!
  !<
  subroutine save_track_data(this, itrkun, csv, itrack1, itrack2)
    ! -- modules
    use ParticleModule, only: ParticleListType
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    integer(I4B), intent(in) :: itrkun
    logical(LGP), intent(in) :: csv
    integer(I4B), intent(in), optional :: itrack1, itrack2
    ! -- local
    integer(I4B) :: itrack, itrackmin, itrackmax
    integer(I4B) :: kper, kstp
    integer(I4B) :: iprp, irpt, icell, izone, istatus, ireason
    real(DP) :: trelease, t, x, y, z

    ! -- select subset of track data between itrack1 and itrack2
    !    if provided (can be used to select data for distinct PRPs)
    if (present(itrack1)) then
      itrackmin = itrack1
    else
      itrackmin = 1
    end if
    if (present(itrack2)) then
      itrackmax = itrack2
    else
      itrackmax = this%nrows
    end if

    ! -- write rows to file
    do itrack = itrackmin, itrackmax
      kper = this%kper(itrack)
      kstp = this%kstp(itrack)
      iprp = this%iprp(itrack)
      irpt = this%irpt(itrack)
      icell = this%icell(itrack)
      izone = this%izone(itrack)
      istatus = this%istatus(itrack)
      ireason = this%ireason(itrack)
      trelease = this%trelease(itrack)
      t = this%t(itrack)
      x = this%x(itrack)
      y = this%y(itrack)
      z = this%z(itrack)

      if (csv) then
        write (itrkun, '(*(G0,:,","))') &
          kper, kstp, iprp, irpt, icell, izone, istatus, ireason, &
          trelease, t, x, y, z
      else
        write (itrkun) &
          kper, kstp, iprp, irpt, icell, izone, istatus, ireason, &
          trelease, t, x, y, z
      end if
    end do

  end subroutine save_track_data

end module TrackDataModule
