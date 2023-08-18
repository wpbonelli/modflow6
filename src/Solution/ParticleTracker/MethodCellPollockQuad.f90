module MethodCellPollockQuadModule

  use KindModule, only: DP, I4B
  use MethodModule
  use MethodSubcellPoolModule
  use CellRectQuadModule
  use CellDefnModule
  use SubcellRectModule
  use ParticleModule
  use UtilMiscModule
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodCellPollockQuadType
  public :: create_methodCellPollockQuad

  type, extends(MethodType) :: MethodCellPollockQuadType
    private
    type(CellRectQuadType), pointer :: cellRectQuad => null() ! tracking domain for the method
    type(SubcellRectType), pointer :: subcellRect => null() ! subcell obj injected into method
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mCPQ ! applies Pollock's cell-quad method
    procedure, public :: pass => pass_mCPQ ! passes particle to next subcell or to cell face
    procedure, public :: loadsub => loadsub_mCPQ ! loads the subcell method
    procedure, public :: load_subcell ! loads the subcell
  end type MethodCellPollockQuadType

contains

  !> @brief Create a new Pollock's cell-quad-method object
  subroutine create_methodCellPollockQuad(methodCellPollockQuad)
    ! -- dummy
    type(MethodCellPollockQuadType), pointer :: methodCellPollockQuad
    !
    allocate (methodCellPollockQuad)
    !
    ! -- This method delegates tracking to a submethod
    methodCellPollockQuad%delegatesTracking = .TRUE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_cellRectQuad(methodCellPollockQuad%cellRectQuad)
    methodCellPollockQuad%trackingDomainType => &
      methodCellPollockQuad%cellRectQuad%type
    !
    ! -- Create subdomain to be loaded and injected into the submethod
    call create_subcellRect(methodCellPollockQuad%subcellRect)
    !
    return
    !
  end subroutine create_methodCellPollockQuad

  !> @brief Destructor for a Pollock's cell-quad-method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  !> @brief Initialize a Pollock's cell-quad-method object
  subroutine init(this, particle, cellRectQuad, trackdata)
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellRectQuadType), pointer, intent(in) :: cellRectQuad
    type(TrackDataType), pointer :: trackdata
    !
    ! -- Set pointer to cell definition
    this%cellRectQuad => cellRectQuad
    !
    ! -- Set pointer to model track data
    this%trackdata => trackdata
    !
    return
    !
  end subroutine init

  !> @brief Load subcell to inject into subcell method
  subroutine loadsub_mCPQ(this, particle, levelNext, submethod)
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    !
    ! -- Load subcell for injection into subcell method
    call this%load_subcell(particle, levelNext, this%subcellRect)
    ! -- Initialize subcell method and set subcell method pointer
    call methodSubcellPollock%init(this%subcellRect, this%trackdata)
    submethod => methodSubcellPollock
    !
    return
    !
  end subroutine loadsub_mCPQ

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mCPQ(this, particle)
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! -- local
    integer :: isc, exitFace, npolyverts, inface, infaceoff
    !
    exitFace = particle%iTrackingDomainBoundary(3)
    isc = particle%iTrackingDomain(3)
    npolyverts = this%cellRectQuad%cellDefn%npolyverts
    !
    select case (exitFace) ! kluge note: exitFace uses Dave's iface convention
    case (0)
      ! -- Subcell interior (cell interior)
      inface = -1
    case (1)
      select case (isc)
      case (1)
        ! -- W face, subcell 1 --> E face, subcell 4  (cell interior)
        particle%iTrackingDomain(3) = 4
        particle%iTrackingDomainBoundary(3) = 2
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (2)
        ! -- W face, subcell 2 --> E face, subcell 3 (cell interior)
        particle%iTrackingDomain(3) = 3
        particle%iTrackingDomainBoundary(3) = 2
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (3)
        ! -- W face, subcell 3 (cell face)
        inface = 1 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = 0
      case (4)
        ! -- W face, subcell 4 (cell face)
        inface = 2 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = -1
      end select
    case (2)
      select case (isc)
      case (1)
        ! -- E face, subcell 1 (cell face)
        inface = 3 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = 0
      case (2)
        ! -- E face, subcell 2 (cell face)
        inface = 4 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = -1
      case (3)
        ! -- E face, subcell 3 --> W face, subcell 2 (cell interior)
        particle%iTrackingDomain(3) = 2
        particle%iTrackingDomainBoundary(3) = 1
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (4)
        ! -- E face, subcell 4 --> W face subcell 1 (cell interior)
        particle%iTrackingDomain(3) = 1
        particle%iTrackingDomainBoundary(3) = 1
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      end select
    case (3)
      select case (isc)
      case (1)
        ! -- S face, subcell 1 --> N face, subcell 2 (cell interior)
        particle%iTrackingDomain(3) = 2
        particle%iTrackingDomainBoundary(3) = 4
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (2)
        ! -- S face, subcell 2 (cell face)
        inface = 4 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = 0
      case (3)
        ! -- S face, subcell 3 (cell face)
        inface = 1 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = -1
      case (4)
        ! -- S face, subcell 4 --> N face, subcell 3 (cell interior)
        particle%iTrackingDomain(3) = 3
        particle%iTrackingDomainBoundary(3) = 4
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      end select
    case (4)
      select case (isc)
      case (1)
        ! -- N face, subcell 1 (cell face)
        inface = 3 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = -1
      case (2)
        ! -- N face, subcell 2 --> S face, subcell 1 (cell interior)
        particle%iTrackingDomain(3) = 1
        particle%iTrackingDomainBoundary(3) = 3
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (3)
        ! -- N face, subcell 3 --> S face, subcell 4 (cell interior)
        particle%iTrackingDomain(3) = 4
        particle%iTrackingDomainBoundary(3) = 3
        inface = 0 ! kluge note: want Domain(2) unchanged; Boundary(2) = 0
      case (4)
        ! -- N face, subcell 4 (cell face)
        inface = 2 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
        infaceoff = 0
      end select
    case (5)
      ! -- Subcell bottom (cell bottom)
      inface = npolyverts + 2 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
    case (6)
      ! -- Subcell top (cell top)
      inface = npolyverts + 3 ! kluge note: want Domain(2) = -Domain(2); Boundary(2) = inface
    end select
    if (inface .eq. -1) then
      particle%iTrackingDomainBoundary(2) = 0
    else if (inface .eq. 0) then
      particle%iTrackingDomainBoundary(2) = 0
    else
      if ((inface .ge. 1) .and. (inface .le. 4)) then
        ! -- Account for local cell rotation
        inface = inface + this%cellRectQuad%irvOrigin - 1
        if (inface .gt. 4) inface = inface - 4
        inface = this%cellRectQuad%irectvert(inface) + infaceoff
        if (inface .lt. 1) inface = inface + npolyverts
      end if
      particle%iTrackingDomainBoundary(2) = inface
    end if
    !
    return
    !
  end subroutine pass_mCPQ

  !> @brief Apply Pollock's quad method to a rectangular-quad cell
  subroutine apply_mCPQ(this, particle, tmax)
    use UtilMiscModule
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! -- local
    double precision :: xOrigin, yOrigin, zOrigin, sinrot, cosrot
    !
    ! -- Update particle zone
    particle%izone = this%cellRectQuad%cellDefn%izone
    !
    if (this%cellRectQuad%cellDefn%izone .ne. 0) then
      if (particle%istopzone .eq. this%cellRectQuad%cellDefn%izone) then
        ! -- Stop zone
        particle%istatus = 6
      end if
    else if (this%cellRectQuad%cellDefn%inoexitface .ne. 0) then
      ! -- No exit face
      particle%istatus = 5
    else if (particle%istopweaksink .ne. 0) then
      if (this%cellRectQuad%cellDefn%iweaksink .ne. 0) then
        ! -- Weak sink
        particle%istatus = 3
      end if
    else
      !
      ! -- If the particle is above the top of the cell (which is presumed to
      ! -- represent a water table above the cell bottom), pass the particle
      ! -- vertically and instantaneously to the cell top elevation.
      if (particle%z > this%cellRectQuad%cellDefn%top) then
        particle%z = this%cellRectQuad%cellDefn%top
        ! -- Store track data
        call this%trackdata%save_record(particle, kper=kper, &
                                        kstp=kstp, reason=1)
      end if
      !
      ! -- Transform particle location into local cell coordinates.
      ! -- Translated and rotated but not scaled relative to model.
      xOrigin = this%cellRectQuad%xOrigin
      yOrigin = this%cellRectQuad%yOrigin
      zOrigin = this%cellRectQuad%zOrigin
      sinrot = this%cellRectQuad%sinrot
      cosrot = this%cellRectQuad%cosrot
      call particle%transf_coords(xOrigin, yOrigin, zOrigin, &
                                  sinrot, cosrot, .false.)
      !
      ! -- Track across subcells
      call this%subtrack(particle, 2, tmax) ! kluge, hardwired to level 2
      !
      ! -- Transform particle location back to model coordinates
      call particle%transf_coords(xOrigin, yOrigin, zOrigin, &
                                  sinrot, cosrot, .true.)
      !
      ! -- Reset transformation and eliminate accumulated roundoff error
      call particle%reset_transf()
      !
    end if
    !
    return
    !
  end subroutine apply_mCPQ

  !> @brief Load the rectangular subcell from the rectangular cell
  subroutine load_subcell(this, particle, levelNext, subcellRect)
    ! -- modules
    use ConstantsModule, only: DONE
    ! -- dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(SubcellRectType), intent(inout) :: subcellRect
    ! -- local
    double precision :: dx, dy, dz, areax, areay, areaz
    double precision :: dxprel, dyprel
    integer :: isc, npolyverts, m1, m2
    double precision :: qextl1, qextl2, qintl1, qintl2
    double precision :: factor, term
    !
    factor = DONE / this%cellRectQuad%cellDefn%retfactor
    factor = factor / this%cellRectQuad%cellDefn%porosity
    npolyverts = this%cellRectQuad%cellDefn%npolyverts
    !
    isc = particle%iTrackingDomain(3)
    ! -- Subcells 1, 2, 3, and 4 are Pollock's subcells A, B, C, and D,
    ! -- respectively
    !
    dx = this%cellRectQuad%dx
    dy = this%cellRectQuad%dy
    ! -- If not already known, determine subcell number
    if (isc .le. 0) then
      dxprel = particle%x / dx
      dyprel = particle%y / dy
      if (dxprel .lt. 5d-1) then
        if (dyprel .lt. 5d-1) then
          isc = 3
        else if (dyprel .gt. 5d-1) then
          isc = 4
        else
          ! kluge note: need to resolve this ambiguity based on flow direction
          print *, "particle initially on shared subcell edge"
          ! pause
          stop
        end if
      else if (dxprel .gt. 5d-1) then
        if (dyprel .lt. 5d-1) then
          isc = 2
        else if (dyprel .gt. 5d-1) then
          isc = 1
        else
          ! kluge note: need to resolve this ambiguity based on flow direction
          print *, "particle initially on shared subcell edge"
          ! pause
          stop
        end if
      else
        ! kluge note: need to resolve this ambiguity based on flow direction
        print *, "particle initially on shared subcell edge"
        ! pause
        stop
      end if
      subcellRect%isubcell = isc
      ! kluge note: as a matter of form, do we want to allow
      ! this subroutine to modify the particle???
      particle%iTrackingDomain(3) = isc
      ! kluge note: initial insubface is not currently being determined
    end if
    dx = 5d-1 * dx
    dy = 5d-1 * dy
    dz = this%cellRectQuad%cellDefn%top - &
         this%cellRectQuad%cellDefn%bot ! kluge note: need to account for partial saturation
    areax = dy * dz
    areay = dx * dz
    areaz = dx * dy
    qintl1 = this%cellRectQuad%qintl(isc)
    ! qintl list wraps around, so isc+1=5 is ok
    qintl2 = this%cellRectQuad%qintl(isc + 1)
    qextl1 = this%cellRectQuad%qextl1(isc)
    qextl2 = this%cellRectQuad%qextl2(isc)
    !
    subcellRect%dx = dx
    subcellRect%dy = dy
    subcellRect%dz = dz
    subcellRect%sinrot = 0d0
    subcellRect%cosrot = 1d0
    subcellRect%zOrigin = 0d0
    select case (isc)
    case (1)
      subcellRect%xOrigin = dx
      subcellRect%yOrigin = dy
      term = factor / areax
      subcellRect%vx1 = qintl1 * term
      subcellRect%vx2 = -qextl2 * term
      term = factor / areay
      subcellRect%vy1 = -qintl2 * term
      subcellRect%vy2 = -qextl1 * term
    case (2)
      subcellRect%xOrigin = dx
      subcellRect%yOrigin = 0d0
      term = factor / areax
      subcellRect%vx1 = -qintl2 * term
      subcellRect%vx2 = -qextl1 * term
      term = factor / areay
      subcellRect%vy1 = qextl2 * term
      subcellRect%vy2 = -qintl1 * term
    case (3)
      subcellRect%xOrigin = 0d0
      subcellRect%yOrigin = 0d0
      term = factor / areax
      subcellRect%vx1 = qextl2 * term
      subcellRect%vx2 = -qintl1 * term
      term = factor / areay
      subcellRect%vy1 = qextl1 * term
      subcellRect%vy2 = qintl2 * term
    case (4)
      subcellRect%xOrigin = 0d0
      subcellRect%yOrigin = dy
      term = factor / areax
      subcellRect%vx1 = qextl1 * term
      subcellRect%vx2 = qintl2 * term
      term = factor / areay
      subcellRect%vy1 = qintl1 * term
      subcellRect%vy2 = -qextl2 * term
    end select
    m1 = npolyverts + 2
    m2 = m1 + 1
    term = factor / areaz
    subcellRect%vz1 = 2.5d-1 * this%cellRectQuad%cellDefn%faceflow(m1) * term
    subcellRect%vz2 = -2.5d-1 * this%cellRectQuad%cellDefn%faceflow(m2) * term
    !
    return
    !
  end subroutine load_subcell

end module MethodCellPollockQuadModule
