module MethodCellPassToBotModule

  use KindModule, only: DP, I4B
  use MethodModule
  use CellDefnModule
  use ParticleModule
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodCellPassToBotType
  public :: create_methodCellPassToBot

  ! -- Extend MethodType to the pass-to-bottom cell-method type
  ! -- (MethodCellPassToBot)
  type, extends(MethodType) :: MethodCellPassToBotType
    private
    type(CellDefnType), pointer, public :: cellDefn => null() ! tracking domain for the method
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mCVP ! applies pass-to-bottom cell method
!!    procedure, public :: pass => pass_mCVP                     ! passes the particle to the cell face
!!    procedure, public :: loadsub => loadsub_mCVP               ! loads the subcell method
!!    procedure, public :: load_subcell                          ! loads the lone subcell (subcell = cell)
  end type MethodCellPassToBotType

contains

  subroutine create_methodCellPassToBot(methodCellPassToBot)
! ******************************************************************************
! create_methodCellPassToBot -- Create a new pass-to-bottom cell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(MethodCellPassToBotType), pointer :: methodCellPassToBot
    ! -- local
! ------------------------------------------------------------------------------
    !
    allocate (methodCellPassToBot)
    allocate (methodCellPassToBot%trackingDomainType)
    !
    ! -- This method delegates tracking to a submethod
    methodCellPassToBot%delegatesTracking = .FALSE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_cellDefn(methodCellPassToBot%cellDefn)
!!    methodCellPassToBot%trackingDomain => methodCellPassToBot%cellDefn
    methodCellPassToBot%trackingDomainType = "CellDefn" ! kluge???
    !
    return
    !
  end subroutine create_methodCellPassToBot

  subroutine destroy(this)
! ******************************************************************************
! destroy -- Destructor for a pass-to-bottom cell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  subroutine init(this, particle, cellDefn, trackdata)
! ******************************************************************************
! init -- Initialize a pass-to-bottom cell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle ! kluge note: is particle needed???
    type(cellDefnType), pointer, intent(in) :: cellDefn
    type(TrackDataType), pointer :: trackdata
! ------------------------------------------------------------------------------
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
    particle%izone = this%cellDefn%izone
    !
    if (this%cellDefn%izone .ne. 0) then
      if (particle%istopzone .eq. this%cellDefn%izone) then
        ! -- Stop zone
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 6
      end if
    else if (this%cellDefn%inoexitface .ne. 0) then
      ! -- No exit face
      ! particle%iTrackingDomainBoundary(3) = 0
      particle%istatus = 5
    else if (particle%istopweaksink .ne. 0) then
      if (this%cellDefn%iweaksink .ne. 0) then
        ! -- Weak sink
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 3
      end if
    else
      !
      ! -- Pass particle vertically and instantaneously to cell bottom
      particle%z = this%cellDefn%bot
      particle%iTrackingDomainBoundary(2) = this%cellDefn%npolyverts + 2
      !
    end if
    !
    ! -- Store track data
    call this%trackdata%save_record(particle, kper=kper, &
                                    kstp=kstp, reason=1)
    !
    return
    !
  end subroutine apply_mCVP

end module MethodCellPassToBotModule
