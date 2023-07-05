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
    procedure, public :: get_model_coords
    procedure, public :: reset_transf
    procedure, public :: set_transf
    procedure, public :: transf_coords
    procedure, public :: init_from_list
  end type ParticleType

  ! -- Define the particle list type (ParticleListType)  ! kluge note: use separate module???
  type ParticleListType

    ! provenance/metadata
    ! integer, public :: imodel ! index of model to which the particle currently belongs
    integer(I4B), dimension(:), pointer, contiguous :: irpt !< release point number
    integer(I4B), dimension(:), pointer, contiguous :: iprp ! particle release point number (index)

    ! stopping criteria
    integer(I4B), dimension(:), pointer, contiguous :: istopweaksink !< weak sink option: 0 = do not stop, 1 = stop
    integer(I4B), dimension(:), pointer, contiguous :: istopzone !< stop zone number

    ! tracking domain
    integer(I4B), dimension(:, :), allocatable :: iTrackingDomain ! array of indices for domains in the tracking domain hierarchy
    integer(I4B), dimension(:, :), allocatable :: iTrackingDomainBoundary ! array of indices for tracking domain boundaries

    ! track data
    integer(I4B), dimension(:), pointer, contiguous :: izone !< current zone number
    integer(I4B), dimension(:), pointer, contiguous :: istatus !< particle status
    real(DP), dimension(:), pointer, contiguous :: x ! model x coord of particle
    real(DP), dimension(:), pointer, contiguous :: y ! model y coord of particle
    real(DP), dimension(:), pointer, contiguous :: z ! model z coord of particle
    real(DP), dimension(:), pointer, contiguous :: trelease ! particle release time
    real(DP), dimension(:), pointer, contiguous :: tstop ! particle stop time
    real(DP), dimension(:), pointer, contiguous :: ttrack ! time to which particle has been tracked

  contains
    procedure, public :: allocate_arrays
    procedure, public :: deallocate_arrays
    procedure, public :: reallocate_arrays
    procedure, public :: update_from_particle
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

  !> @brief Allocate particle arrays
  subroutine allocate_arrays(this, np, lmin, lmax, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(ParticleListType), intent(inout) :: this
    integer(I4B), intent(in) :: np ! number of particles
    integer(I4B), intent(in) :: lmin ! minimum level in the tracking domain hierarchy
    integer(I4B), intent(in) :: lmax ! maximum level in the tracking domain hierarchy
    character(*), intent(in) :: mempath ! path to memory
    !
    call mem_allocate(this%irpt, np, 'PLIRPT', mempath)
    call mem_allocate(this%iprp, np, 'PLIPRP', mempath)
    ! -- kluge todo: update mem_allocate to allow custom range of indices?
    !    e.g. here we want to allocate 0-4 for trackdomain levels, not 1-5
    ! call mem_allocate(this%iTrackingDomain, np, lmax - lmin + 1, 'PLITD', mempath) ! kluge note: ditch crazy dims
    ! call mem_allocate(this%iTrackingDomainBoundary, np, &
    !                   lmax - lmin + 1, 'PLITDB', mempath) ! kluge note: ditch crazy dims
    allocate (this%iTrackingDomain(np, lmin:lmax))
    allocate (this%iTrackingDomainBoundary(np, lmin:lmax))
    call mem_allocate(this%izone, np, 'PLIZONE', mempath)
    call mem_allocate(this%istatus, np, 'PLISTATUS', mempath)
    call mem_allocate(this%x, np, 'PLX', mempath)
    call mem_allocate(this%y, np, 'PLY', mempath)
    call mem_allocate(this%z, np, 'PLZ', mempath)
    call mem_allocate(this%trelease, np, 'PLTRELEASE', mempath)
    call mem_allocate(this%tstop, np, 'PLTSTOP', mempath)
    call mem_allocate(this%ttrack, np, 'PLTTRACK', mempath)
    call mem_allocate(this%istopweaksink, np, 'PLISTOPWEAKSINK', mempath)
    call mem_allocate(this%istopzone, np, 'PLISTOPZONE', mempath)
    !
    return
  end subroutine allocate_arrays

  !> @brief Deallocate particle arrays
  subroutine deallocate_arrays(this, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(ParticleListType), intent(inout) :: this
    character(*), intent(in) :: mempath ! path to memory
    !
    call mem_deallocate(this%irpt, 'PLIRPT', mempath)
    call mem_deallocate(this%iprp, 'PLIPRP', mempath)
    ! call mem_deallocate(this%iTrackingDomain, 'PLITD', mempath)
    ! call mem_deallocate(this%iTrackingDomainBoundary, 'PLITDB', mempath)
    deallocate (this%iTrackingDomain)
    deallocate (this%iTrackingDomainBoundary)
    call mem_deallocate(this%izone, 'PLIZONE', mempath)
    call mem_deallocate(this%istatus, 'PLISTATUS', mempath)
    call mem_deallocate(this%x, 'PLX', mempath)
    call mem_deallocate(this%y, 'PLY', mempath)
    call mem_deallocate(this%z, 'PLZ', mempath)
    call mem_deallocate(this%trelease, 'PLTRELEASE', mempath)
    call mem_deallocate(this%tstop, 'PLTSTOP', mempath)
    call mem_deallocate(this%ttrack, 'PLTTRACK', mempath)
    call mem_deallocate(this%istopweaksink, 'PLISTOPWEAKSINK', mempath)
    call mem_deallocate(this%istopzone, 'PLISTOPZONE', mempath)
    !
    return
  end subroutine deallocate_arrays

  !> @brief Reallocate particle arrays
  subroutine reallocate_arrays(this, np, mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    use ArrayHandlersModule, only: ExpandArray2D
    ! -- dummy
    class(ParticleListType), intent(inout) :: this
    integer(I4B), intent(in) :: np ! number of particles
    character(*), intent(in) :: mempath ! path to memory
    !
    ! resize 1D arrays
    call mem_reallocate(this%irpt, np, 'PLIRPT', mempath)
    call mem_reallocate(this%iprp, np, 'PLIPRP', mempath)
    call mem_reallocate(this%izone, np, 'PLIZONE', mempath)
    call mem_reallocate(this%istatus, np, 'PLISTATUS', mempath)
    call mem_reallocate(this%x, np, 'PLX', mempath)
    call mem_reallocate(this%y, np, 'PLY', mempath)
    call mem_reallocate(this%z, np, 'PLZ', mempath)
    call mem_reallocate(this%trelease, np, 'PLTRELEASE', mempath)
    call mem_reallocate(this%tstop, np, 'PLTSTOP', mempath)
    call mem_reallocate(this%ttrack, np, 'PLTTRACK', mempath)
    call mem_reallocate(this%istopweaksink, np, 'PLISTOPWEAKSINK', mempath)
    call mem_reallocate(this%istopzone, np, 'PLISTOPZONE', mempath)
    !
    ! resize first dimension of 2D arrays
    call ExpandArray2D( &
      this%iTrackingDomain, &
      size(this%iTrackingDomain(:, 1) - np), &
      0)
    call ExpandArray2D( &
      this%iTrackingDomainBoundary, &
      size(this%iTrackingDomainBoundary(:, 1) - np), &
      0)
    !
    return
  end subroutine reallocate_arrays

  !> @brief Initialize particle from particle list
  subroutine init_from_list(this, partlist, im, iprp, irpt)
    ! -- dummy
    class(ParticleType), intent(inout) :: this
    type(ParticleListType), intent(in) :: partlist
    integer(I4B), intent(in) :: im ! model ID
    integer(I4B), intent(in) :: iprp ! particle release package ID
    integer(I4B), intent(in) :: irpt ! particle release point ID
    !
    this%iprp = iprp
    this%irpt = irpt ! kluge note: necessary to reset this here?
    this%istopweaksink = partlist%istopweaksink(irpt)
    this%istopzone = partlist%istopzone(irpt)
    this%iTrackingDomain(levelMin:levelMax) = &
      partlist%iTrackingDomain(irpt, levelMin:levelMax)
    this%iTrackingDomain(1) = im
    this%iTrackingDomainBoundary(levelMin:levelMax) = &
      partlist%iTrackingDomainBoundary(irpt, levelMin:levelMax)
    this%istatus = -1 ! need to reset this to -1 every timestep (for solution to proceed)
    this%x = partlist%x(irpt)
    this%y = partlist%y(irpt)
    this%z = partlist%z(irpt)
    this%trelease = partlist%trelease(irpt)
    this%tstop = partlist%tstop(irpt)
    this%ttrack = partlist%ttrack(irpt)
    !
    return
  end subroutine init_from_list

  !> @brief Update particle list from particle
  subroutine update_from_particle(this, particle, np)
    ! -- dummy
    class(ParticleListType), intent(inout) :: this
    type(ParticleType), intent(in) :: particle
    integer(I4B), intent(in) :: np
    !
    this%iTrackingDomain( &
      np, &
      levelMin:levelMax) = &
      particle%iTrackingDomain(levelMin:levelMax)
    this%iTrackingDomainBoundary( &
      np, &
      levelMin:levelMax) = &
      particle%iTrackingDomainBoundary(levelMin:levelMax)
    this%izone(np) = particle%izone
    this%istatus(np) = particle%istatus
    this%x(np) = particle%x
    this%y(np) = particle%y
    this%z(np) = particle%z
    this%ttrack(np) = particle%ttrack
    !
    return
  end subroutine update_from_particle

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
