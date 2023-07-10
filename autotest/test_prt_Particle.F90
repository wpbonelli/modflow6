module test_prt_particle
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use ParticleModule, only: ParticleType, ParticleListType, &
                            create_particle, get_particle_id
  implicit none
  private
  public :: collect_prt_Particle
contains

  subroutine collect_prt_Particle(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("dynamic_allocation", test_dynamic_allocation), &
                new_unittest("get_particle_id", test_get_particle_id) &
                ]
  end subroutine collect_prt_Particle

  subroutine test_get_particle_id(error)
    type(error_type), allocatable, intent(out) :: error
    type(ParticleType), pointer :: p

    call create_particle(p)
    p%iprp = 0
    p%irpt = 0
    p%trelease = 0.0
    print *, "Particle ID: ", get_particle_id(p)
    call check(error, get_particle_id(p) == "0-0-0.")

  end subroutine test_get_particle_id

  subroutine test_dynamic_allocation(error)
    ! -- modules
    use GlobalDataModule, only: levelMin, levelMax
    use ConstantsModule, only: LENMEMPATH, I4B

    ! -- dummy
    type(error_type), allocatable, intent(out) :: error

    ! -- locals
    type(ParticleListType), pointer :: partlist => null()
    integer(I4B) :: npartmax1 = 10
    integer(I4B) :: npartmax2 = 20
    integer(I4B) :: nlevels = levelMax - levelMin + 1
    character(len=LENMEMPATH) :: mempath = "test_dynamic_allocation"

    ! allocate particle arrays
    allocate (partlist)
    call partlist%allocate_arrays(npartmax1, levelMin, levelMax, mempath)

    ! check initial array sizes
    call check(error, size(partlist%iprp) == npartmax1)
    call check(error, size(partlist%irpt) == npartmax1)
    call check(error, size(partlist%istopweaksink) == npartmax1)
    call check(error, size(partlist%istopzone) == npartmax1)
    call check(error, size(partlist%iTrackingDomain, 1) == npartmax1)
    call check(error, size(partlist%iTrackingDomain, 2) == nlevels)
    call check(error, size(partlist%iTrackingDomainBoundary, 1) == npartmax1)
    call check(error, size(partlist%iTrackingDomainBoundary, 2) == nlevels)
    call check(error, size(partlist%icu) == npartmax1)
    call check(error, size(partlist%ilay) == npartmax1)
    call check(error, size(partlist%izone) == npartmax1)
    call check(error, size(partlist%istatus) == npartmax1)
    call check(error, size(partlist%x) == npartmax1)
    call check(error, size(partlist%y) == npartmax1)
    call check(error, size(partlist%z) == npartmax1)
    call check(error, size(partlist%trelease) == npartmax1)
    call check(error, size(partlist%tstop) == npartmax1)
    call check(error, size(partlist%ttrack) == npartmax1)
    if (allocated(error)) return

    ! resize particle arrays
    call partlist%reallocate_arrays(npartmax2, mempath)

    ! check that arrays have been resized
    call check(error, size(partlist%iprp) == npartmax2)
    call check(error, size(partlist%irpt) == npartmax2)
    call check(error, size(partlist%istopweaksink) == npartmax2)
    call check(error, size(partlist%istopzone) == npartmax2)
    call check(error, size(partlist%iTrackingDomain, 1) == npartmax2)
    call check(error, size(partlist%iTrackingDomain, 2) == nlevels)
    call check(error, size(partlist%iTrackingDomainBoundary, 1) == npartmax2)
    call check(error, size(partlist%iTrackingDomainBoundary, 2) == nlevels)
    call check(error, size(partlist%icu) == npartmax2)
    call check(error, size(partlist%ilay) == npartmax2)
    call check(error, size(partlist%istatus) == npartmax2)
    call check(error, size(partlist%x) == npartmax2)
    call check(error, size(partlist%y) == npartmax2)
    call check(error, size(partlist%z) == npartmax2)
    call check(error, size(partlist%trelease) == npartmax2)
    call check(error, size(partlist%tstop) == npartmax2)
    call check(error, size(partlist%ttrack) == npartmax2)
    if (allocated(error)) return

  end subroutine test_dynamic_allocation

end module test_prt_particle
