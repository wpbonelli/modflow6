module ParticleModule

  use GlobalDataModule
  implicit none

  private
  public :: ParticleType
  public :: ParticleListType
  public :: create_particle

  ! -- Define the particle type (ParticleType)
  type ParticleType
    private
    integer, public :: ipart ! particle number (index)
    integer, public :: iprp ! particle release point number (index)
    ! integer, public :: igroup ! index of particle group to which the particle belongs
    ! integer, public :: imodel ! index of model to which the particle currently belongs
    ! double precision, public :: velmult ! velocity multiplier for the particle (determined by group)
    double precision, public :: x ! model x coord of particle
    double precision, public :: y ! model y coord of particle
    double precision, public :: z ! model z coord of particle
    double precision, public :: xlocal ! local x coord of particle
    double precision, public :: ylocal ! local y coord of particle
    double precision, public :: zlocal ! local z coord of particle
    integer, allocatable, public :: iTrackingDomain(:) ! array of indices for domains in the tracking domain hierarchy
    integer, allocatable, public :: iTrackingDomainBoundary(:) ! array of indices for tracking domain boundaries
    double precision, public :: trelease ! particle release time
    double precision, public :: tstop ! particle stop time
    double precision, public :: ttrack ! time to which particle has been tracked
    integer, public :: istopweaksink !< weak sink option: 0 = do not stop, 1 = stop
    integer, public :: istopzone !< stop zone number
    integer, public :: istatus !< particle status
    integer, public :: irpt !< release point number

  contains
    procedure, public :: destroy => destroy_particle ! destructor for the particle
  end type ParticleType

  ! -- Define the particle list type (ParticleListType)  ! kluge note: use separate module???
  type ParticleListType
    ! type(ParticleType), pointer :: particle
    ! double precision, allocatable, public :: velmult(:) ! velocity multiplier for the particle (determined by group)
    double precision, allocatable, public :: x(:) ! model x coord of particle
    double precision, allocatable, public :: y(:) ! model y coord of particle
    double precision, allocatable, public :: z(:) ! model z coord of particle
    ! double precision, allocatable, public :: xlocal(:) ! local x coord of particle
    ! double precision, allocatable, public :: ylocal(:) ! local y coord of particle
    ! double precision, allocatable, public :: zlocal(:) ! local z coord of particle
    integer, allocatable, public :: iTrackingDomain(:, :) ! array of indices for domains in the tracking domain hierarchy
    integer, allocatable, public :: iTrackingDomainBoundary(:, :) ! array of indices for tracking domain boundaries
    double precision, allocatable, public :: trelease(:) ! particle release time
    double precision, allocatable, public :: tstop(:) ! particle stop time
    double precision, allocatable, public :: ttrack(:) ! time to which particle has been tracked
    integer, allocatable, public :: istopweaksink(:) !< weak sink option: 0 = do not stop, 1 = stop
    integer, allocatable, public :: istopzone(:) !< stop zone number
    integer, allocatable, public :: istatus(:) !< particle status
    integer, allocatable, public :: irpt(:) !< release point number
  end type ParticleListType

contains

  !> @brief Create a new particle
  subroutine create_particle(particle)
    ! -- dummy
    type(ParticleType), pointer :: particle
    !
    allocate (particle)
    allocate (particle%iTrackingDomain(levelMin:levelMax))
    allocate (particle%iTrackingDomainBoundary(levelMin:levelMax))
    return
  end subroutine create_particle

  !> @brief Destructor for a particle
  subroutine destroy_particle(this)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    !
    deallocate (this%iTrackingDomain)
    deallocate (this%iTrackingDomainBoundary)
    return
  end subroutine destroy_particle

end module ParticleModule
