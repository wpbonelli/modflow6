module ParticleModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LENMEMPATH
  use GlobalDataModule
  use UtilMiscModule, only: transform_coords, modify_transf
  implicit none

  private
  public :: ParticleType
  public :: ParticleListType
  public :: create_particle
  public :: get_particle_id

  ! -- Define the particle type (ParticleType)
  type ParticleType
    private

    ! provenance/metadata
    ! integer, public :: imodel ! index of model to which the particle currently belongs
    integer(I4B), public :: irpt !< release point number
    integer(I4B), public :: iprp ! particle release point number (index)

    ! stop criteria
    integer(I4B), public :: istopweaksink !< weak sink option: 0 = do not stop, 1 = stop
    integer(I4B), public :: istopzone !< stop zone number

    ! tracking domain
    integer(I4B), allocatable, public :: iTrackingDomain(:) ! array of indices for domains in the tracking domain hierarchy
    integer(I4B), allocatable, public :: iTrackingDomainBoundary(:) ! array of indices for tracking domain boundaries

    ! track data
    integer(I4B), public :: izone !< current zone number
    integer(I4B), public :: istatus !< particle status
    real(DP), public :: x ! model x coord of particle
    real(DP), public :: y ! model y coord of particle
    real(DP), public :: z ! model z coord of particle
    real(DP), public :: trelease ! particle release time
    real(DP), public :: tstop ! particle stop time
    real(DP), public :: ttrack ! time to which particle has been tracked

    ! transformation
    logical(LGP), public :: isTransformed !< indicates whether coordinates have been transformed from model to local coordinates
    real(DP), public :: xOrigin !< x origin for coordinate transformation from model to local
    real(DP), public :: yOrigin !< y origin for coordinate transformation from model to local
    real(DP), public :: zOrigin !< z origin for coordinate transformation from model to local
    real(DP), public :: sinrot !< sine of rotation angle for coordinate transformation from model to local
    real(DP), public :: cosrot !< cosine of rotation angle for coordinate transformation from model to local

  contains
    procedure, public :: destroy => destroy_particle ! destructor for the particle
    procedure, public :: set_transf
    procedure, public :: reset_transf
    procedure, public :: transf_coords
    procedure, public :: get_model_coords
  end type ParticleType

  ! -- Define the particle list type (ParticleListType)  ! kluge note: use separate module???
  type ParticleListType

    ! provenance/metadata
    ! integer, public :: imodel ! index of model to which the particle currently belongs
    integer(I4B), allocatable, public :: irpt(:) !< release point number
    integer(I4B), allocatable, public :: iprp(:) ! particle release point number (index)

    ! stopping criteria
    integer(I4B), allocatable, public :: istopweaksink(:) !< weak sink option: 0 = do not stop, 1 = stop
    integer(I4B), allocatable, public :: istopzone(:) !< stop zone number

    ! tracking domain
    integer(I4B), allocatable, public :: iTrackingDomain(:, :) ! array of indices for domains in the tracking domain hierarchy
    integer(I4B), allocatable, public :: iTrackingDomainBoundary(:, :) ! array of indices for tracking domain boundaries

    ! track data
    integer(I4B), allocatable, public :: izone(:) !< current zone number
    integer(I4B), allocatable, public :: istatus(:) !< particle status
    real(DP), allocatable, public :: x(:) ! model x coord of particle
    real(DP), allocatable, public :: y(:) ! model y coord of particle
    real(DP), allocatable, public :: z(:) ! model z coord of particle
    real(DP), allocatable, public :: trelease(:) ! particle release time
    real(DP), allocatable, public :: tstop(:) ! particle stop time
    real(DP), allocatable, public :: ttrack(:) ! time to which particle has been tracked

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

  !> @brief Directly set transformation from model coordinates
  !! to particle's local coordinates
  !<
  subroutine set_transf(this, xOrigin, yOrigin, zOrigin, sinrot, cosrot)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    real(DP), intent(in) :: xOrigin
    real(DP), intent(in) :: yOrigin
    real(DP), intent(in) :: zOrigin
    real(DP), intent(in) :: sinrot
    real(DP), intent(in) :: cosrot
    ! -- local
    !
    this%xOrigin = xOrigin
    this%yOrigin = yOrigin
    this%zOrigin = zOrigin
    this%sinrot = sinrot
    this%cosrot = cosrot
    !
    this%isTransformed = .true.
    !
    return
  end subroutine set_transf

  subroutine transf_coords(this, xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                           sinrot_opt, cosrot_opt, invert_opt)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    real(DP), intent(in), optional :: xOrigin_opt
    real(DP), intent(in), optional :: yOrigin_opt
    real(DP), intent(in), optional :: zOrigin_opt
    real(DP), intent(in), optional :: sinrot_opt
    real(DP), intent(in), optional :: cosrot_opt
    logical(LGP), intent(in), optional :: invert_opt
    !
    ! -- Transform particle's coordinates
    call transform_coords(this%x, this%y, this%z, this%x, this%y, this%z, &
                          xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                          sinrot_opt, cosrot_opt, invert_opt)
    !
    ! -- Modify transformation from model coordinates to particle's new
    ! -- local coordinates by incorporating this latest transformation
    call modify_transf(this%xOrigin, this%yOrigin, this%zOrigin, &
                       this%sinrot, this%cosrot, &
                       xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                       sinrot_opt, cosrot_opt, invert_opt)
    ! -- Set isTransformed flag to true. Note that there is no check
    ! -- to see whether the modification brings the coordinates back
    ! -- to model coordinates (in which case the origin would be very
    ! -- close to zero and sinrot and cosrot would be very close to 0.
    ! -- and 1., respectively, allowing for roundoff error).
    this%isTransformed = .true.
    !
    return
  end subroutine transf_coords

  !> @brief Reset coordinate transformation to none
  !<
  subroutine reset_transf(this)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    ! -- local
    !
    this%xOrigin = DZERO
    this%yOrigin = DZERO
    this%zOrigin = DZERO
    this%sinrot = DZERO
    this%cosrot = DONE
    !
    this%isTransformed = .false.
    !
    return
  end subroutine reset_transf

  !> @brief Return the particle's model coordinates
  !<
  subroutine get_model_coords(this, xmodel, ymodel, zmodel)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    real(DP) :: xmodel
    real(DP) :: ymodel
    real(DP) :: zmodel
    ! -- local
    !
    if (this%isTransformed) then
      ! -- Transform back from local to model coordinates
      call transform_coords(this%x, this%y, this%z, xmodel, ymodel, zmodel, &
                            this%xOrigin, this%yOrigin, this%zOrigin, &
                            this%sinrot, this%cosrot, .true.)
    else
      ! -- Already in model coordinates
      xmodel = this%x
      ymodel = this%y
      zmodel = this%z
    end if
    !
    return
  end subroutine get_model_coords

  function get_particle_id(particle) result(id)
    ! -- dummy
    class(ParticleType), intent(in) :: particle
    character(len=LENMEMPATH) :: id
    !
    write (id, '(I0,"-",I0,"-",F0.0)') &
      particle%iprp, particle%irpt, particle%trelease
    !
    return
  end function get_particle_id

end module ParticleModule
