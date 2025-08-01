module MethodSubcellModule
  use KindModule, only: DP, I4B
  use MethodCellModule, only: MethodCellType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType, &
                                 SubcellExitEventType

  private
  public :: MethodSubcellType

  type, abstract, extends(MethodCellType) :: MethodSubcellType
  contains
    procedure, public :: subcellexit
    procedure, public :: assess
  end type MethodSubcellType

contains

  !> @brief Particle exits a subcell.
  subroutine subcellexit(this, particle)
    class(MethodSubcellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event

    allocate (SubcellExitEventType :: event)
    select type (event)
    type is (SubcellExitEventType)
      event%exit_face = particle%iboundary(3)
      event%isc = particle%idomain(3)
    end select
    call this%events%dispatch(particle, event)
  end subroutine subcellexit

  subroutine assess(this, particle, cell_defn, tmax)
    ! dummy
    class(MethodSubcellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! noop
  end subroutine assess

end module MethodSubcellModule