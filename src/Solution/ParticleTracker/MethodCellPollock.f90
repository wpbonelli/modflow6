module MethodCellPollockModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE
  use MethodModule
  use MethodSubcellPoolModule
  ! use CellModule
  use CellRectModule
  ! use CellDefnModule
  use SubcellRectModule
  use ParticleModule
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodCellPollockType
  public :: create_methodCellPollock

  ! -- Extend MethodType to the Pollock's cell-method type (MethodCellPollockType)
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
    ! methodCellPollock%trackingDomain => methodCellPollock%cellRect
    methodCellPollock%trackingDomainType => methodCellPollock%cellRect%type ! kluge note: "lazy" allocation here and in similar places
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
  subroutine init(this, particle, cellRect, trackdata)
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellRectType), pointer, intent(in) :: cellRect
    type(TrackDataType), pointer :: trackdata
    !
    ! -- Set pointer to cell definition
    this%cellRect => cellRect
    !
    ! -- Set pointer to model track data
    this%trackdata => trackdata
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
    call methodSubcellPollock%init(this%subcellRect)
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
    ! particle%iTrackingDomain(3) = 0            ! kluge for recursive  ! kluge note: set a "has_exited" attribute instead???
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
    ! if ((inface.ge.1).and.(inface.le.4)) then
    !   ! -- Account for local cell rotation
    !   inface = inface + this%cellRect%ipvOrigin - 1
    !   if (inface.gt.4) inface = inface - 4
    ! end if
    ! particle%iTrackingDomainBoundary(2) = inface
    if (inface .eq. -1) then
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      ! particle%iTrackingDomainBoundary(2) = 0
      ! particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))   ! kluge???
      particle%iTrackingDomainBoundary(2) = 0
    else
      if ((inface .ge. 1) .and. (inface .le. 4)) then
        ! -- Account for local cell rotation
        inface = inface + this%cellRect%ipvOrigin - 1
        if (inface .gt. 4) inface = inface - 4
      end if
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      particle%iTrackingDomainBoundary(2) = inface
      ! particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))   ! kluge???
    end if
    ! particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))   ! kluge???
    !
    return
    !
  end subroutine pass_mCP

  !> @brief Apply Pollock's method to a rectangular cell
  subroutine apply_mCP(this, particle, tmax)
    use UtilMiscModule
    ! -- dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! doubleprecision :: initialTime,maximumTime,t   ! kluge not in arg list yet
    ! -- local
    double precision :: xOrigin, yOrigin, zOrigin, sinrot, cosrot
    integer(I4B) :: ntrack
    !
    if (this%cellRect%cellDefn%izone .ne. 0) then
      if (particle%istopzone .eq. this%cellRect%cellDefn%izone) then
        ! -- Stop zone
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 6
        ! write(*,'(A,I,A,I)') "particle ", particle%ipart, " terminated in stop zone cell: ", particle%iTrackingDomain(2)  ! kluge
        ! return
      end if
    else if (this%cellRect%cellDefn%inoexitface .ne. 0) then
      ! -- No exit face
      ! particle%iTrackingDomainBoundary(3) = 0
      particle%istatus = 5
      ! write(*,'(A,I,A,I)') "particle ", particle%ipart, " terminated at cell w/ no exit face: ", particle%iTrackingDomain(2)  ! kluge
      ! return
    else if (particle%istopweaksink .ne. 0) then
      if (this%cellRect%cellDefn%iweaksink .ne. 0) then
        ! -- Weak sink
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 3
        ! write(*,'(A,I,A,I)')  "particle ", particle%ipart, " terminated at weak sink cell: ", particle%iTrackingDomain(2)  ! kluge
        ! return
      end if
    else
      !
      ! -- If the particle is above the top of the cell (which is presumed to
      ! -- represent a water table above the cell bottom), pass the particle
      ! -- vertically and instantaneously to the cell top elevation.
      if (particle%z > this%cellRect%cellDefn%top) then
        particle%z = this%cellRect%cellDefn%top
        ! -- Store track data
        ntrack = this%trackdata%ntrack + 1 ! kluge?
        this%trackdata%ntrack = ntrack
        this%trackdata%iptrack(ntrack) = particle%ipart
        this%trackdata%ictrack(ntrack) = particle%iTrackingDomain(2)
        this%trackdata%xtrack(ntrack) = particle%x
        this%trackdata%ytrack(ntrack) = particle%y
        this%trackdata%ztrack(ntrack) = particle%z
        this%trackdata%ttrack(ntrack) = particle%ttrack
      end if
      !
      ! -- Transform particle location into local cell coordinates
      xOrigin = this%cellRect%xOrigin
      yOrigin = this%cellRect%yOrigin
      zOrigin = this%cellRect%zOrigin
      sinrot = this%cellRect%sinrot
      cosrot = this%cellRect%cosrot
      call transform_coords(particle%x, particle%y, particle%z, &
                            xOrigin, yOrigin, zOrigin, sinrot, cosrot, .false., &
                            particle%xlocal, particle%ylocal, particle%zlocal)
      !
      ! -- Track across lone subcell
      call this%subtrack(particle, 2, tmax) ! kluge, hardwired to level 2
      !
      ! particle%iTrackingDomainBoundary(2) = 0
      ! !
      ! -- Transform particle location back to model coordinates
      call transform_coords(particle%xlocal, particle%ylocal, particle%zlocal, &
                            xOrigin, yOrigin, zOrigin, sinrot, cosrot, .true., &
                            particle%x, particle%y, particle%z)
      !
    end if
    !
    ! -- Store track data
    ntrack = this%trackdata%ntrack + 1 ! kluge?
    this%trackdata%ntrack = ntrack
    this%trackdata%iptrack(ntrack) = particle%ipart
    this%trackdata%ictrack(ntrack) = particle%iTrackingDomain(2)
    this%trackdata%xtrack(ntrack) = particle%x
    this%trackdata%ytrack(ntrack) = particle%y
    this%trackdata%ztrack(ntrack) = particle%z
    this%trackdata%ttrack(ntrack) = particle%ttrack
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
    ! -- local
    ! double precision :: velmult
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
    ! velmult = particle%velmult                      ! kluge note: apply velmult later, in subcell method???
    ! subcellRect%vx1 = velmult*this%cellRect%vx1     ! kluge note: assuming porosity=1. for now
    ! subcellRect%vx2 = velmult*this%cellRect%vx2
    ! subcellRect%vy1 = velmult*this%cellRect%vy1
    ! subcellRect%vy2 = velmult*this%cellRect%vy2
    ! subcellRect%vz1 = velmult*this%cellRect%vz1
    ! subcellRect%vz2 = velmult*this%cellRect%vz2
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
