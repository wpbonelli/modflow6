module MethodCellPassToBotModule

  use KindModule, only: DP, I4B
  use MethodModule
  use CellDefnModule
  use ParticleModule
  use TrackModule, only: TrackControlType
  implicit none

  private
  public :: MethodCellPassToBotType
  public :: create_methodCellPassToBot

  type, extends(MethodType) :: MethodCellPassToBotType
    private
    type(CellDefnType), pointer, public :: cellDefn => null() ! tracking domain for the method
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mCVP ! applies pass-to-bottom cell method
  end type MethodCellPassToBotType

contains

  !> @brief Create a new pass-to-bottom cell-method object
  subroutine create_methodCellPassToBot(methodCellPassToBot)
    ! -- dummy
    type(MethodCellPassToBotType), pointer :: methodCellPassToBot
    ! -- local
    !
    allocate (methodCellPassToBot)
    allocate (methodCellPassToBot%trackingDomainType)
    !
    ! -- This method delegates tracking to a submethod
    methodCellPassToBot%delegatesTracking = .FALSE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_cellDefn(methodCellPassToBot%cellDefn)
    methodCellPassToBot%trackingDomainType = "CellDefn" ! kluge???
    !
    return
    !
  end subroutine create_methodCellPassToBot

  !> @brief Destructor for a pass-to-bottom cell-method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    ! -- local
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  !> @brief Initialize a pass-to-bottom cell-method object
  subroutine init(this, particle, cellDefn, trackdata)
    ! -- dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(cellDefnType), pointer, intent(in) :: cellDefn
    type(TrackControlType), pointer :: trackdata
    !
    ! -- Set pointer to cell definition
    this%cellDefn => cellDefn
    !
    ! -- Set pointer to model track data
    this%trackdata => trackdata
    !
    return
    !
  end subroutine init

  !> @brief Apply pass-to-bottom method to a cell
  subroutine apply_mCVP(this, particle, tmax)
    ! -- modules
    use UtilMiscModule
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    !
    ! -- Update particle state, checking whether any reporting or
    ! -- termination conditions apply
    call this%update(particle, this%cellDefn)
    !
    ! -- Return early if particle is done advancing
    if (.not. particle%advancing) return
    !
    ! -- Pass particle vertically and instantaneously to cell bottom
    particle%z = this%cellDefn%bot
    particle%iTrackingDomainBoundary(2) = this%cellDefn%npolyverts + 2
    !
    ! -- Store track data
    call this%trackdata%save_record(particle, kper=kper, &
                                    kstp=kstp, reason=1)
    !
    return
    !
  end subroutine apply_mCVP

end module MethodCellPassToBotModule
