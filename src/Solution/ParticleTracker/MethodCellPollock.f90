module MethodCellPollockModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE
  use MethodModule
  use MethodSubcellPoolModule
  use CellRectModule
  use SubcellRectModule
  use ParticleModule
  use TrackModule, only: TrackControlType
  implicit none

  private
  public :: MethodCellPollockType
  public :: create_methodCellPollock

  type, extends(MethodType) :: MethodCellPollockType
    private
    type(CellRectType), pointer, public :: cellRect => null() ! tracking domain for the method
    type(SubcellRectType), pointer :: subcellRect ! subcell object injected into subcell method
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mCP ! applies Pollock's cell method
    procedure, public :: pass => pass_mCP ! passes the particle to the cell face
    procedure, public :: loadsub => loadsub_mCP ! loads the subcell method
    procedure, public :: load_subcell ! loads the lone subcell (subcell = cell)
  end type MethodCellPollockType

contains

  !> @brief Create a new Pollock's cell-method object
  subroutine create_methodCellPollock(methodCellPollock)
    ! -- dummy
    type(MethodCellPollockType), pointer :: methodCellPollock
    ! -- local
    !
    allocate (methodCellPollock)
    !
    ! -- This method delegates tracking to a submethod
    methodCellPollock%delegatesTracking = .TRUE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_cellRect(methodCellPollock%cellRect)
    ! kluge note: "lazy" allocation here and in similar places
    methodCellPollock%trackingDomainType => methodCellPollock%cellRect%type
    !
    ! -- Create subdomain to be loaded and injected into the submethod
    call create_subcellRect(methodCellPollock%subcellRect)
    !
    return
    !
  end subroutine create_methodCellPollock

  !> @brief Destructor for a Pollock's cell-method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    ! -- local
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  !> @brief Initialize a Pollock's cell-method object
  subroutine init(this, particle, cellRect, trackctl)
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellRectType), pointer, intent(in) :: cellRect
    type(TrackControlType), pointer :: trackctl
    !
    ! -- Set pointer to cell definition
    this%cellRect => cellRect
    !
    ! -- Set pointer to model track output control
    this%trackctl => trackctl
    !
    return
    !
  end subroutine init

  !> @brief Load subcell method
  subroutine loadsub_mCP(this, particle, levelNext, submethod)
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    !
    ! -- Load rectangular subcell for injection into Pollock's subcell method
    call this%load_subcell(particle, levelNext, this%subcellRect)
    ! -- Select and initialize Pollock's subcell method and set subcell
    ! -- method pointer
    call methodSubcellPollock%init(this%subcellRect, this%trackctl)
    submethod => methodSubcellPollock
    !
    return
    !
  end subroutine loadsub_mCP

  !> @brief Having exited the lone subcell, pass the particle to the cell face
  !! In this case the lone subcell is the cell.
  subroutine pass_mCP(this, particle)
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! -- local
    integer :: exitFace, inface
    !
    exitFace = particle%iTrackingDomainBoundary(3)
    ! -- Map subcell exit face to cell face
    select case (exitFace) ! kluge note: exitFace uses Dave's iface convention
    case (0)
      inface = -1
    case (1)
      inface = 1
    case (2)
      inface = 3
    case (3)
      inface = 4
    case (4)
      inface = 2
    case (5)
      inface = 6 ! kluge note: inface=5 same as inface=1 due to wraparound
    case (6)
      inface = 7
    end select
    if (inface .eq. -1) then
      particle%iTrackingDomainBoundary(2) = 0
    else
      if ((inface .ge. 1) .and. (inface .le. 4)) then
        ! -- Account for local cell rotation
        inface = inface + this%cellRect%ipvOrigin - 1
        if (inface .gt. 4) inface = inface - 4
      end if
      particle%iTrackingDomainBoundary(2) = inface
    end if
    !
    return
    !
  end subroutine pass_mCP

  !> @brief Apply Pollock's method to a rectangular cell
  subroutine apply_mCP(this, particle, tmax)
    use UtilMiscModule
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! -- local
    double precision :: xOrigin, yOrigin, zOrigin, sinrot, cosrot
    !
    ! -- Update particle state, checking whether any reporting or
    ! -- termination conditions apply
    call this%update(particle, this%cellRect%cellDefn)
    !
    ! -- Return early if particle is done advancing
    if (.not. particle%advancing) return
    !
    ! -- If the particle is above the top of the cell (which is presumed to
    ! -- represent a water table above the cell bottom), pass the particle
    ! -- vertically and instantaneously to the cell top elevation and save
    ! -- the particle state to output file(s).
    if (particle%z > this%cellRect%cellDefn%top) then
      particle%z = this%cellRect%cellDefn%top
      call this%trackctl%save_record(particle, kper=kper, &
                                     kstp=kstp, reason=1) ! reason=1: cell transition
    end if
    !
    ! -- Transform particle location into local cell coordinates.
    ! -- Translated and rotated but not scaled relative to model.
    xOrigin = this%cellRect%xOrigin
    yOrigin = this%cellRect%yOrigin
    zOrigin = this%cellRect%zOrigin
    sinrot = this%cellRect%sinrot
    cosrot = this%cellRect%cosrot
    call particle%transf_coords(xOrigin, yOrigin, zOrigin, &
                                sinrot, cosrot, .false.)
    !
    ! -- Track across lone subcell
    call this%subtrack(particle, 2, tmax) ! kluge, hardwired to level 2
    !
    ! -- Transform particle location back to model coordinates
    call particle%transf_coords(xOrigin, yOrigin, zOrigin, &
                                sinrot, cosrot, .true.)
    !
    ! -- Reset transformation and eliminate accumulated roundoff error
    call particle%reset_transf()
    !
    return
    !
  end subroutine apply_mCP

  !> @brief Loads the lone rectangular subcell from the rectangular cell
  !! kluge note: is levelNext needed here and in similar "load" routines???
  subroutine load_subcell(this, particle, levelNext, subcellRect) !
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    type(SubcellRectType), pointer, intent(inout) :: subcellRect
    !
    ! -- Set subcell number to 1
    subcellRect%isubcell = 1
    particle%iTrackingDomain(levelNext) = 1 ! kluge note: is this the place to set this???
    !
    ! -- Subcell calculations will be done in local subcell coordinates
    subcellRect%dx = this%cellRect%dx
    subcellRect%dy = this%cellRect%dy
    subcellRect%dz = this%cellRect%dz
    subcellRect%sinrot = 0d0
    subcellRect%cosrot = 1d0 ! kluge note: rethink how/where to store subcell data???
    subcellRect%xOrigin = 0d0
    subcellRect%yOrigin = 0d0
    subcellRect%zOrigin = 0d0
    !
    ! -- Set subcell edge velocities
    subcellRect%vx1 = this%cellRect%vx1 ! kluge note: cell velocities now already account for retfactor and porosity
    subcellRect%vx2 = this%cellRect%vx2
    subcellRect%vy1 = this%cellRect%vy1
    subcellRect%vy2 = this%cellRect%vy2
    subcellRect%vz1 = this%cellRect%vz1
    subcellRect%vz2 = this%cellRect%vz2
    !
    return
    !
  end subroutine load_subcell

end module MethodCellPollockModule
