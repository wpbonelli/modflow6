module test_prt_trackdata
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use ParticleModule, only: ParticleType, ParticleListType, &
                            create_particle, get_particle_id
  use TrackDataModule, only: TrackDataType
  implicit none
  private
  public :: collect_prt_TrackData
contains

  subroutine collect_prt_TrackData(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("reallocate", test_reallocate), &
                new_unittest("add_track_data", test_add_track_data) &
                ]
  end subroutine collect_prt_TrackData

  subroutine test_reallocate(error)
    ! -- modules
    use GlobalDataModule, only: levelMin, levelMax
    use ConstantsModule, only: LENMEMPATH, I4B

    ! -- dummy
    type(error_type), allocatable, intent(out) :: error

    ! -- locals
    type(ParticleType), pointer :: particle => null()
    type(TrackDataType), pointer :: trackdata => null()
    integer(I4B) :: nt1 = 10
    integer(I4B) :: nt2 = 20
    character(len=LENMEMPATH) :: mempath = "test_reallocate"

    ! allocate and initialize particle
    call create_particle(particle)
    particle%irpt = 0
    particle%iprp = 0
    particle%istopzone = 0
    particle%istopweaksink = 0
    particle%iTrackingDomain(levelMin:levelMax) = 0
    particle%iTrackingDomainBoundary(levelMin:levelMax) = 0
    particle%izone = 0
    particle%istatus = 0
    particle%x = 0
    particle%y = 0
    particle%z = 0
    particle%trelease = 0
    particle%tstop = 0
    particle%ttrack = 0
    particle%isTransformed = .false.
    particle%xOrigin = 0.0
    particle%yOrigin = 0.0
    particle%zOrigin = 0.0
    particle%sinrot = 0.0
    particle%cosrot = 0.0

    ! allocate track data
    allocate (trackdata)
    call trackdata%allocate_arrays(nt1, mempath)

    ! check initial array sizes
    call check(error, size(trackdata%kper) == nt1)
    call check(error, size(trackdata%kstp) == nt1)
    call check(error, size(trackdata%irpt) == nt1)
    call check(error, size(trackdata%iprp) == nt1)
    call check(error, size(trackdata%icell) == nt1)
    call check(error, size(trackdata%izone) == nt1)
    call check(error, size(trackdata%istatus) == nt1)
    call check(error, size(trackdata%ireason) == nt1)
    call check(error, size(trackdata%trelease) == nt1)
    call check(error, size(trackdata%t) == nt1)
    call check(error, size(trackdata%x) == nt1)
    call check(error, size(trackdata%y) == nt1)
    call check(error, size(trackdata%z) == nt1)
    if (allocated(error)) return

    ! resize particle arrays
    call trackdata%reallocate_arrays(nt2, mempath)

    ! check that arrays have been resized
    call check(error, size(trackdata%kper) == nt2)
    call check(error, size(trackdata%kstp) == nt2)
    call check(error, size(trackdata%irpt) == nt2)
    call check(error, size(trackdata%iprp) == nt2)
    call check(error, size(trackdata%icell) == nt2)
    call check(error, size(trackdata%izone) == nt2)
    call check(error, size(trackdata%istatus) == nt2)
    call check(error, size(trackdata%ireason) == nt2)
    call check(error, size(trackdata%trelease) == nt2)
    call check(error, size(trackdata%t) == nt2)
    call check(error, size(trackdata%x) == nt2)
    call check(error, size(trackdata%y) == nt2)
    call check(error, size(trackdata%z) == nt2)
    if (allocated(error)) return

  end subroutine test_reallocate

  subroutine test_add_track_data(error)
    ! -- modules
    use GlobalDataModule, only: levelMin, levelMax
    use ConstantsModule, only: LENMEMPATH, I4B

    ! -- dummy
    type(error_type), allocatable, intent(out) :: error

    ! -- locals
    type(ParticleType), pointer :: particle
    type(TrackDataType), pointer :: trackdata => null()
    integer(I4B) :: nt1 = 1
    integer(I4B) :: nt2
    integer(I4B) :: kper = 1, kstp = 1
    character(len=LENMEMPATH) :: mempath = "test_add_track_data"

    ! allocate and initialize particle
    print *, "allocating particle"
    call create_particle(particle)
    particle%irpt = 0
    particle%iprp = 0
    particle%istopzone = 0
    particle%istopweaksink = 0
    particle%iTrackingDomain(levelMin:levelMax) = 0
    particle%iTrackingDomainBoundary(levelMin:levelMax) = 0
    particle%izone = 0
    particle%istatus = 0
    particle%x = 0
    particle%y = 0
    particle%z = 0
    particle%trelease = 0
    particle%tstop = 0
    particle%ttrack = 0
    particle%isTransformed = .false.
    particle%xOrigin = 0.0
    particle%yOrigin = 0.0
    particle%zOrigin = 0.0
    particle%sinrot = 0.0
    particle%cosrot = 0.0

    ! allocate track data to initial size of 1
    print *, "allocating trackdata"
    allocate (trackdata)
    allocate (trackdata%ntrack)
    trackdata%ntrack = 0
    call trackdata%allocate_arrays(nt1, mempath)

    ! check initial array sizes
    print *, "checking initial trackdata size"
    call check(error, size(trackdata%kper) == nt1)
    call check(error, size(trackdata%kstp) == nt1)
    call check(error, size(trackdata%irpt) == nt1)
    call check(error, size(trackdata%iprp) == nt1)
    call check(error, size(trackdata%icell) == nt1)
    call check(error, size(trackdata%izone) == nt1)
    call check(error, size(trackdata%istatus) == nt1)
    call check(error, size(trackdata%ireason) == nt1)
    call check(error, size(trackdata%trelease) == nt1)
    call check(error, size(trackdata%t) == nt1)
    call check(error, size(trackdata%x) == nt1)
    call check(error, size(trackdata%y) == nt1)
    call check(error, size(trackdata%z) == nt1)
    if (allocated(error)) return

    ! add particle to track data
    print *, "adding first particle to track data"
    call trackdata%add_track_data(particle, kper=kper, &
                                  kstp=kstp, reason=0)

    ! check track data values
    print *, "checking initial track data values"
    call check(error, trackdata%kper(1) == 1)
    call check(error, trackdata%kstp(1) == 1)
    call check(error, trackdata%irpt(1) == 0)
    call check(error, trackdata%iprp(1) == 0)
    call check(error, trackdata%icell(1) == 0)
    call check(error, trackdata%izone(1) == 0)
    call check(error, trackdata%istatus(1) == 0)
    call check(error, trackdata%ireason(1) == 0)
    call check(error, trackdata%trelease(1) == 0)
    call check(error, trackdata%t(1) == 0)
    call check(error, trackdata%x(1) == 0)
    call check(error, trackdata%y(1) == 0)
    call check(error, trackdata%z(1) == 0)
    if (allocated(error)) return

    ! add another particle to track data
    print *, "adding another particle to track data"
    call trackdata%add_track_data(particle, kper=kper, &
                                  kstp=kstp, reason=0)

    ! check that arrays were automatically expanded by factor of 10
    print *, "checking arrays were expanded by factor of 10"
    nt2 = nt1 * 10
    call check(error, size(trackdata%kper) == nt2)
    call check(error, size(trackdata%kstp) == nt2)
    call check(error, size(trackdata%irpt) == nt2)
    call check(error, size(trackdata%iprp) == nt2)
    call check(error, size(trackdata%icell) == nt2)
    call check(error, size(trackdata%izone) == nt2)
    call check(error, size(trackdata%istatus) == nt2)
    call check(error, size(trackdata%ireason) == nt2)
    call check(error, size(trackdata%trelease) == nt2)
    call check(error, size(trackdata%t) == nt2)
    call check(error, size(trackdata%x) == nt2)
    call check(error, size(trackdata%y) == nt2)
    call check(error, size(trackdata%z) == nt2)
    if (allocated(error)) return

  end subroutine test_add_track_data

end module test_prt_trackdata
