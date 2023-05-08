module MethodCellTernaryModule

  use KindModule, only: DP, I4B
  use MethodModule
  use MethodSubcellPoolModule
  use CellPolyModule
  use CellDefnModule
  use SubcellTriModule
  use ParticleModule
  use Ternary ! kluge
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_methodCellTernary

  ! -- Extend MethodType to the ternary cell-method type (MethodCellTernaryType)
  type, extends(MethodType) :: MethodCellTernaryType
    private
    type(CellPolyType), pointer :: cellPoly => null() ! tracking domain for the method

    ! subcell object injected into subcell method
    ! kluge note: redundant to store these, kluge "99"
    type(subcellTriType), pointer :: subcellTri
    double precision :: x_vert(99), y_vert(99) ! cell vertex coordinates
    double precision :: xctr, yctr ! cell center coordinates
    double precision :: vx_vert_polygon(99), vy_vert_polygon(99) ! cell vertex velocities
    double precision :: vxctr, vyctr ! cell center velocities
    double precision :: ztop, zbot ! cell top and bottom elevations
    double precision :: dz ! cell thickness
    double precision :: vztop, vzbot ! cell top and bottom velocities
  contains
    procedure, public :: destroy => destroy ! destructor for the method
    procedure, public :: init => init ! initializes the method
    procedure, public :: apply => apply_mCT ! applies the ternary cell method
    procedure, public :: pass => pass_mCT ! passes particle to next subcell or to cell face
    procedure, public :: loadsub => loadsub_mCT ! loads the subcell method
    procedure, public :: load_subcell => load_subcell ! loads the subcell
  end type MethodCellTernaryType

contains

  !> @brief Create a new ternary cell-method object
  subroutine create_methodCellTernary(methodCellTernary)
    ! -- dummy
    type(MethodCellTernaryType), pointer :: methodCellTernary
    ! -- local
    !
    allocate (methodCellTernary)
    !
    ! -- This method delegates tracking to a submethod
    methodCellTernary%delegatesTracking = .TRUE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_cellPoly(methodCellTernary%cellPoly)
    ! methodCellTernary%trackingDomain => methodCellTernary%cellPoly
    methodCellTernary%trackingDomainType => methodCellTernary%cellPoly%type
    !
    ! -- Create subdomain to be loaded and injected into the submethod
    call create_subcellTri(methodCellTernary%subcellTri)
    !
    return
    !
  end subroutine create_methodCellTernary

  !> @brief Destructor for a ternary cell-method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! -- local
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  !> @brief Initialize a ternary cell-method object
  subroutine init(this, particle, cellPoly, trackdata)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellPolyType), pointer, intent(in) :: cellPoly
    type(TrackDataType), pointer :: trackdata
    !
    this%cellPoly => cellPoly
    !
    ! -- Set pointer to model track data
    this%trackdata => trackdata
    !
    return
    !
  end subroutine init

  !> @brief Load subcell to inject into subcell method
  subroutine loadsub_mCT(this, particle, levelNext, submethod)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    !
    ! -- Load subcell for injection into subcell method
    call this%load_subcell(particle, levelNext, this%subcellTri)
    ! -- Initialize subcell method and set subcell method pointer
    call methodSubcellTernary%init(this%subcellTri)
    submethod => methodSubcellTernary
    !
    return
    !
  end subroutine loadsub_mCT

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mCT(this, particle)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer :: isc, exitFace, inface, npolyverts
    !
    exitFace = particle%iTrackingDomainBoundary(3)
    isc = particle%iTrackingDomain(3)
    npolyverts = this%cellPoly%cellDefn%npolyverts
    !
    select case (exitFace)
    case (0)
      ! -- Subcell interior (cell interior)
      inface = -1
    case (1)
      ! -- Subcell face 1 (cell face)
      ! inface = npolyverts - isc        ! kluge note: is this general???
      inface = isc
      ! if (inface.eq.0) inface = 4
      if (inface .eq. 0) inface = npolyverts
    case (2)
      ! -- Subcell face --> next subcell in "cycle" (cell interior)
      isc = isc + 1
      if (isc .gt. npolyverts) isc = 1
      particle%iTrackingDomain(3) = isc
      particle%iTrackingDomainBoundary(3) = 3
      inface = 0
    case (3)
      ! -- Subcell face --> preceding subcell in "cycle" (cell interior)
      isc = isc - 1
      if (isc .lt. 1) isc = npolyverts
      particle%iTrackingDomain(3) = isc
      particle%iTrackingDomainBoundary(3) = 2
      inface = 0
    case (4)
      ! -- Subcell bottom (cell bottom)
      inface = npolyverts + 2
    case (5)
      ! -- Subcell top (cell top)
      inface = npolyverts + 3
    end select
    ! particle%iTrackingDomainBoundary(2) = inface
    ! if (inface.ne.0) particle%iTrackingDomain(3) = 0
    if (inface .eq. -1) then
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      ! particle%iTrackingDomainBoundary(2) = 0
      ! particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))   ! kluge???
      particle%iTrackingDomainBoundary(2) = 0
    else if (inface .eq. 0) then
      particle%iTrackingDomainBoundary(2) = 0
    else
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      particle%iTrackingDomainBoundary(2) = inface
      ! particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))   ! kluge???
    end if
    ! if (inface.ne.0) particle%iTrackingDomain(3) = -abs(particle%iTrackingDomain(3))
    !
    return
    !
  end subroutine pass_mCT

  !> @brief Apply the ternary method to a polygonal cell
  subroutine apply_mCT(this, particle, tmax)
    use ConstantsModule, only: DZERO, DONE, DHALF ! kluge???
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    ! double precision :: velmult    ! kluge???
    double precision :: retfactor
    integer :: npolyverts, iv, ivp1, ivm1
    double precision :: x0, y0, x1, y1, x2, y2, xsum, ysum
    double precision :: vxsum, vysum, flow0, flow1, v0x, v0y
    double precision :: d01x, d01y, d02x, d02y, det, area, term
    integer(I4B) :: ntrack
    !
    if (this%cellPoly%cellDefn%izone .ne. 0) then
      if (particle%istopzone .eq. this%cellPoly%cellDefn%izone) then
        ! -- Stop zone
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 6
        ! write(*,'(A,I,A,I)') "particle ", particle%ipart, &
        ! " terminated in stop zone cell: ", particle%iTrackingDomain(2)  ! kluge
        ! return
      end if
    else if (this%cellPoly%cellDefn%inoexitface .ne. 0) then
      ! -- No exit face
      ! particle%iTrackingDomainBoundary(3) = 0
      particle%istatus = 5
      ! write(*,'(A,I,A,I)') "particle ", particle%ipart, &
      ! " terminated at cell w/ no exit face: ", particle%iTrackingDomain(2)  ! kluge
      ! return
    else if (particle%istopweaksink .ne. 0) then
      if (this%cellPoly%cellDefn%iweaksink .ne. 0) then
        ! -- Weak sink
        ! particle%iTrackingDomainBoundary(3) = 0
        particle%istatus = 3
        ! write(*,'(A,I,A,I)')  "particle ", particle%ipart, &
        ! " terminated at weak sink cell: ", particle%iTrackingDomain(2)  ! kluge
        ! return
      end if
    else
      !
      ! -- If the particle is above the top of the cell (which is presumed to
      ! -- represent a water table above the cell bottom), pass the particle
      ! -- vertically and instantaneously to the cell top elevation.
      if (particle%z > this%cellPoly%cellDefn%top) then
        particle%z = this%cellPoly%cellDefn%top
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
      npolyverts = this%cellPoly%cellDefn%npolyverts
      !
      xsum = DZERO
      ysum = DZERO
      vxsum = DZERO
      vysum = DZERO
      area = DZERO
      this%ztop = this%cellPoly%cellDefn%top
      this%zbot = this%cellPoly%cellDefn%bot
      this%dz = this%ztop - this%zbot
      do iv = 1, npolyverts
        ivp1 = iv + 1
        if (ivp1 .gt. npolyverts) ivp1 = 1
        ivm1 = iv - 1
        if (ivm1 .lt. 1) ivm1 = npolyverts
        x0 = this%cellPoly%cellDefn%polyvert(iv)%x
        y0 = this%cellPoly%cellDefn%polyvert(iv)%y
        x2 = this%cellPoly%cellDefn%polyvert(ivp1)%x
        y2 = this%cellPoly%cellDefn%polyvert(ivp1)%y
        x1 = this%cellPoly%cellDefn%polyvert(ivm1)%x
        y1 = this%cellPoly%cellDefn%polyvert(ivm1)%y
        ! kluge note: assuming porosity=1. for now
        ! flow0 = this%cellPoly%cellDefn%faceflow(iv)/this%dz
        ! flow1 = this%cellPoly%cellDefn%faceflow(ivm1)/this%dz
        term = DONE / (this%cellPoly%cellDefn%porosity * this%dz)
        flow0 = this%cellPoly%cellDefn%faceflow(iv) * term
        flow1 = this%cellPoly%cellDefn%faceflow(ivm1) * term
        d01x = x1 - x0 ! kluge note: do this more efficiently, not recomputing things so much???
        d01y = y1 - y0
        d02x = x2 - x0
        d02y = y2 - y0
        ! kluge note: can det ever be zero, like maybe for a 180-deg vertex???
        ! oodet = DONE/(d01y*d02x - d02y*d01x)
        ! velmult = particle%velmult
        ! kluge note: "flow" is volumetric (face) flow rate per unit thickness, divided by porosity
        ! v0x = -velmult*oodet*(d02x*flow1 + d01x*flow0)
        ! v0y = -velmult*oodet*(d02y*flow1 + d01y*flow0)   !
        det = d01y * d02x - d02y * d01x
        retfactor = this%cellPoly%cellDefn%retfactor
        ! kluge note: can det ever be zero, like maybe for a 180-deg vertex???
        ! term = velfactor/det
        ! kluge note: can det ever be zero, like maybe for a 180-deg vertex???
        term = DONE / (retfactor * det)
        ! kluge note: "flow" here is volumetric flow rate (MODFLOW face flow)
        v0x = -term * (d02x * flow1 + d01x * flow0)
        ! per unit thickness, divided by porosity
        v0y = -term * (d02y * flow1 + d01y * flow0)
        this%vx_vert_polygon(iv) = v0x
        this%vy_vert_polygon(iv) = v0y
        xsum = xsum + x0
        ysum = ysum + y0
        vxsum = vxsum + v0x
        vysum = vysum + v0y
        this%x_vert(iv) = x0
        this%y_vert(iv) = y0
        area = area + x0 * y1 - x1 * y0
      end do
      ! kluge note: from get_cell2d_area
      !   a = 1/2 *[(x1*y2 + x2*y3 + x3*y4 + ... + xn*y1) -
      !             (x2*y1 + x3*y2 + x4*y3 + ... + x1*yn)]
      area = area * DHALF
      ! this%vzbot = velmult*this%cellPoly%cellDefn%faceflow(npolyverts+2)/area
      ! this%vztop = -velmult*this%cellPoly%cellDefn%faceflow(npolyverts+3)/area
      ! term = velfactor/(this%cellPoly%cellDefn%porosity*area)
      term = DONE / (retfactor * this%cellPoly%cellDefn%porosity * area)
      this%vzbot = this%cellPoly%cellDefn%faceflow(npolyverts + 2) * term
      this%vztop = -this%cellPoly%cellDefn%faceflow(npolyverts + 3) * term
      this%xctr = xsum / dble(npolyverts)
      this%yctr = ysum / dble(npolyverts)
      this%vxctr = vxsum / dble(npolyverts)
      this%vyctr = vysum / dble(npolyverts)
      !
      ! -- Track across subcells
      call this%subtrack(particle, 2, tmax) ! kluge, hardwired to level 2
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
  end subroutine apply_mCT

  !> @brief Loads the triangular subcell from the polygonal cell
  subroutine load_subcell(this, particle, levelNext, subcellTri)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(SubcellTriType), intent(inout) :: subcellTri
    ! -- local
    ! kluge note: maybe (in general) do tracking calc without velmult
    ! and divide exit time by velmult at the end???
    ! double precision :: velmult
    integer :: ic, isc, npolyverts
    integer :: iv0, iv1, ipv0, ipv1
    integer :: iv
    double precision :: x0, y0, x1, y1, x2, y2, x1rel, y1rel, x2rel, y2rel, xi, yi
    double precision :: di2, d02, d12, di1, d01, alphai, betai
    double precision :: betatol ! kluge
    !
    ic = this%cellPoly%cellDefn%icell
    subcellTri%icell = ic
    isc = particle%iTrackingDomain(3)
    npolyverts = this%cellPoly%cellDefn%npolyverts
    !
    ! -- Find subcell if not known       ! kluge note: from "find_init_triangle"
    if (isc .le. 0) then
      xi = particle%x
      yi = particle%y
      do iv = 1, npolyverts
        iv0 = iv
        iv1 = iv + 1
        if (iv1 .gt. npolyverts) iv1 = 1
        ! ipv0 = ivert_polygon(icell,iv0)
        ! ipv1 = ivert_polygon(icell,iv1)
        ipv0 = iv0 ! kluge???
        ipv1 = iv1
        x0 = this%x_vert(ipv0)
        y0 = this%y_vert(ipv0)
        x1 = this%x_vert(ipv1)
        y1 = this%y_vert(ipv1)
        x2 = this%xctr
        y2 = this%yctr
        x1rel = x1 - x0
        y1rel = y1 - y0
        x2rel = x2 - x0
        y2rel = y2 - y0
        di2 = xi * y2rel - yi * x2rel
        d02 = x0 * y2rel - y0 * x2rel
        d12 = x1rel * y2rel - y1rel * x2rel
        di1 = xi * y1rel - yi * x1rel
        d01 = x0 * y1rel - y0 * x1rel
        alphai = (di2 - d02) / d12
        betai = -(di1 - d01) / d12
        ! kluge note: can iTrackingDomainBoundary(2) be used to identify the subcell?
        betatol = -1e-7 ! kluge
        ! kluge note: think this handles points on triangle boundaries ok
        if ((alphai .ge. 0d0) .and. &
            (betai .ge. betatol) .and. &
            (alphai + betai .le. 1d0)) then
          isc = iv ! but maybe not!!!!!!!!!!!!
          exit ! kluge note: doesn't handle particle smack on cell center
        end if
      end do
      if (isc .le. 0) then
        write (*, '(A,I0,A,I0)') &
          "error -- initial triangle not found for particle ", &
          particle%ipart, " in cell ", ic ! kluge
        write (69, '(A,I0,A,I0)') &
          "error -- initial triangle not found for particle ", &
          particle%ipart, " in cell ", ic
        ! pause
        stop
      else
        ! subcellTri%isubcell = isc
        ! kluge note: as a matter of form, do we want to allow
        ! this subroutine to modify the particle???
        particle%iTrackingDomain(3) = isc
      end if
    end if
    subcellTri%isubcell = isc
    !
    ! -- Set coordinates and velocities at vertices of triangular subcell
    iv0 = isc
    iv1 = isc + 1
    if (iv1 .gt. npolyverts) iv1 = 1
    ! ipv0 = ivert_polygon(ic,iv0)
    ! ipv1 = ivert_polygon(ic,iv1)
    ipv0 = iv0 ! kluge???
    ipv1 = iv1
    subcellTri%x0 = this%x_vert(ipv0)
    subcellTri%y0 = this%y_vert(ipv0)
    subcellTri%x1 = this%x_vert(ipv1)
    subcellTri%y1 = this%y_vert(ipv1)
    subcellTri%x2 = this%xctr
    subcellTri%y2 = this%yctr
    subcellTri%v0x = this%vx_vert_polygon(iv0)
    subcellTri%v0y = this%vy_vert_polygon(iv0)
    subcellTri%v1x = this%vx_vert_polygon(iv1)
    subcellTri%v1y = this%vy_vert_polygon(iv1)
    subcellTri%v2x = this%vxctr
    subcellTri%v2y = this%vyctr
    subcellTri%ztop = this%ztop
    subcellTri%zbot = this%zbot
    subcellTri%dz = this%dz
    subcellTri%vzbot = this%vzbot
    subcellTri%vztop = this%vztop
    !
    return
    !
  end subroutine load_subcell

end module MethodCellTernaryModule
