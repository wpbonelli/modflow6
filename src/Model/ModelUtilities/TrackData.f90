module TrackDataModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ParticleModule, only: ParticleType
  use UtilMiscModule, only: transform_coords

  implicit none

  private
  public :: TrackDataType

  integer(I4B), parameter, public :: INITIAL_TRACK_SIZE = 1000

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
    character(len=:), pointer :: mempath => null() ! memory path of track data
    integer(I4B), pointer :: ntrack => null() ! total count of track data
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
    !   1: cross spatial boundary (cell? subcell? or generic feature? worth distinguishing?)
    !   2: cross temporal boundary (time step end)
    !   3: termination
    !   4: inactive

    ! double arrays
    real(DP), dimension(:), pointer, contiguous :: trelease ! particle's release time
    real(DP), dimension(:), pointer, contiguous :: t ! current time
    real(DP), dimension(:), pointer, contiguous :: x ! current x coordinate
    real(DP), dimension(:), pointer, contiguous :: y ! current y coordinate
    real(DP), dimension(:), pointer, contiguous :: z ! current z coordinate

  contains
    procedure, public :: allocate_arrays
    procedure, public :: deallocate_arrays
    procedure, public :: reallocate_arrays
    procedure, public :: add_track_data
    procedure, public :: reset_track_data
    procedure, public :: save_track_data
  end type TrackDataType

contains

  subroutine allocate_arrays(this, nt, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    integer(I4B), intent(in) :: nt
    character(len=*), intent(in) :: mempath
    !
    print *, 'allocating ', nt, ' slots for track data arrays'
    allocate (character(len=len(mempath)) :: this%mempath)
    ! call mem_allocate(this%mempath, size(mempath), 'TRACKMEMPATH', mempath)
    this%mempath = mempath ! kluge!!!
    ! call mem_setptr(mempath, 'TRACKMEMPATH', mempath)
    call mem_allocate(this%kper, nt, 'TRACKKPER', mempath)
    call mem_allocate(this%kstp, nt, 'TRACKKSTP', mempath)
    call mem_allocate(this%iprp, nt, 'TRACKIPRP', mempath)
    call mem_allocate(this%irpt, nt, 'TRACKIRPT', mempath)
    call mem_allocate(this%icell, nt, 'TRACKICELL', mempath)
    call mem_allocate(this%izone, nt, 'TRACKIZONE', mempath)
    call mem_allocate(this%istatus, nt, 'TRACKISTATUS', mempath)
    call mem_allocate(this%ireason, nt, 'TRACKIREASON', mempath)
    call mem_allocate(this%trelease, nt, 'TRACKTRELEASE', mempath)
    call mem_allocate(this%t, nt, 'TRACKT', mempath)
    call mem_allocate(this%x, nt, 'TRACKX', mempath)
    call mem_allocate(this%y, nt, 'TRACKY', mempath)
    call mem_allocate(this%z, nt, 'TRACKZ', mempath)
    !
    return
  end subroutine allocate_arrays

  subroutine deallocate_arrays(this, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    character(len=*), intent(in) :: mempath
    !
    deallocate (this%mempath) ! kluge!!!
    call mem_deallocate(this%kper)
    call mem_deallocate(this%kstp)
    call mem_deallocate(this%iprp)
    call mem_deallocate(this%irpt)
    call mem_deallocate(this%icell)
    call mem_deallocate(this%izone)
    call mem_deallocate(this%istatus)
    call mem_deallocate(this%ireason)
    call mem_deallocate(this%trelease)
    call mem_deallocate(this%t)
    call mem_deallocate(this%x)
    call mem_deallocate(this%y)
    call mem_deallocate(this%z)
    !
    return
  end subroutine deallocate_arrays

  subroutine reallocate_arrays(this, nt, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    integer(I4B), intent(in) :: nt
    character(len=*), intent(in) :: mempath
    !
    call mem_reallocate(this%kper, nt, 'TRACKKPER', mempath)
    call mem_reallocate(this%kstp, nt, 'TRACKKSTP', mempath)
    call mem_reallocate(this%iprp, nt, 'TRACKIPRP', mempath)
    call mem_reallocate(this%irpt, nt, 'TRACKIRPT', mempath)
    call mem_reallocate(this%icell, nt, 'TRACKICELL', mempath)
    call mem_reallocate(this%izone, nt, 'TRACKIZONE', mempath)
    call mem_reallocate(this%istatus, nt, 'TRACKISTATUS', mempath)
    call mem_reallocate(this%ireason, nt, 'TRACKIREASON', mempath)
    call mem_reallocate(this%trelease, nt, 'TRACKTRELEASE', mempath)
    call mem_reallocate(this%t, nt, 'TRACKT', mempath)
    call mem_reallocate(this%x, nt, 'TRACKX', mempath)
    call mem_reallocate(this%y, nt, 'TRACKY', mempath)
    call mem_reallocate(this%z, nt, 'TRACKZ', mempath)
    !
    return
  end subroutine reallocate_arrays

  !> @brief Add track data from a particle
  subroutine add_track_data(this, particle, reason, level)
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level
    ! -- local
    integer(I4B) :: ntrack, ntracksize
    integer(I4B) :: resizefactor, resizethresh
    real(DP) :: resizefraction
    logical(LGP) :: ladd
    real(DP) :: xmodel, ymodel, zmodel
    !
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
    !
    if (ladd) then
      !
      ! -- Expand track arrays by factor of 10 if at capacity.
      resizefactor = 10
      resizethresh = 100
      resizefraction = 0.01
      ntracksize = size(this%irpt)
      if ((ntracksize - this%ntrack) < 1) then
        ! print *, 'Expanding track arrays from ', ntracksize, &
        !   ' to ', ntracksize * resizefactor
        call this%reallocate_arrays(ntracksize * resizefactor, this%mempath)
      end if
      !
      ! -- Get model coordinates
      call particle%get_model_coords(xmodel, ymodel, zmodel)
      !
      ! -- Add track data
      ntrack = this%ntrack + 1
      this%ntrack = ntrack
      this%kper = kper
      this%kstp = kstp
      this%irpt(ntrack) = particle%irpt
      this%iprp = particle%iprp
      this%icell(ntrack) = particle%iTrackingDomain(2)
      this%izone = particle%izone
      this%istatus = particle%istatus
      this%ireason = reason
      this%trelease(ntrack) = particle%trelease
      this%t(ntrack) = particle%ttrack
      this%x(ntrack) = xmodel
      this%y(ntrack) = ymodel
      this%z(ntrack) = zmodel

    end if
    !
    return
  end subroutine add_track_data

  !> @brief Reset track data
  subroutine reset_track_data(this)
    ! -- dummy
    class(TrackDataType), intent(inout) :: this
    !
    this%ntrack = DZERO ! kluge note: zero out arrays, too, just for cleanliness?
    !
    return
  end subroutine reset_track_data

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
      itrackmax = this%ntrack
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
