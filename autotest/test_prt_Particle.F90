module test_prt_particle
    use testdrive, only: error_type, unittest_type, new_unittest, check
    use ParticleModule, only: ParticleType, create_particle, get_particle_id
    implicit none
    private
    public :: collect_prt_Particle
    contains

    subroutine collect_prt_Particle(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("get_particle_id", test_get_particle_id)]
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
end module test_prt_particle
