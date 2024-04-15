module MethodCellTernaryModule

  use KindModule, only: DP, I4B
  use ErrorUtilModule, only: pstop
  use MethodModule
  use MethodSubcellPoolModule
  use CellPolyModule
  use CellDefnModule
  use SubcellTriModule, only: SubcellTriType, create_subcell_tri
  use ParticleModule
  use TrackModule, only: TrackFileControlType
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_method_cell_ternary

  ! input: number of vertices, vertex coordinates, and exterior edge velocities
  integer                                 :: nverts
  double precision, allocatable, dimension(:) :: xvert
  double precision, allocatable, dimension(:) :: yvert
  double precision, allocatable, dimension(:) :: vne
  ! intermediate: exterior edges
  double precision, allocatable, dimension(:) :: le
  double precision, allocatable, dimension(:) :: unex
  double precision, allocatable, dimension(:) :: uney
  ! intermediate: areas and divergence
  double precision               :: areacell
  double precision, allocatable, dimension(:) :: areasub
  double precision               :: divcell
  ! intermediate: interior edges
  double precision, allocatable, dimension(:) :: li
  double precision, allocatable, dimension(:) :: unix
  double precision, allocatable, dimension(:) :: uniy
  ! intermediate: interior edge midpoints
  double precision, allocatable, dimension(:) :: xmid
  double precision, allocatable, dimension(:) :: ymid
  double precision, allocatable, dimension(:) :: lm
  double precision, allocatable, dimension(:) :: umx
  double precision, allocatable, dimension(:) :: umy
  double precision, allocatable, dimension(:) :: kappax
  double precision, allocatable, dimension(:) :: kappay
  double precision, allocatable, dimension(:) :: vm0i
  double precision, allocatable, dimension(:) :: vm0e
  double precision, allocatable, dimension(:) :: vm1i
  double precision, allocatable, dimension(:) :: vm1e
  double precision, allocatable, dimension(:) :: vm0x
  double precision, allocatable, dimension(:) :: vm0y
  double precision, allocatable, dimension(:) :: vm1x
  double precision, allocatable, dimension(:) :: vm1y
  ! output: centroid coordinates and centroid and corner velocities
  double precision               :: xctrd
  double precision               :: yctrd
  double precision               :: vctrdx
  double precision               :: vctrdy
  double precision, allocatable, dimension(:) :: vc0x
  double precision, allocatable, dimension(:) :: vc0y
  double precision, allocatable, dimension(:) :: vc1x
  double precision, allocatable, dimension(:) :: vc1y
  ! intermediate: defined for convenience
  integer, allocatable, dimension(:)          :: iprev
  double precision, allocatable, dimension(:) :: xvertnext
  double precision, allocatable, dimension(:) :: yvertnext
  double precision, allocatable, dimension(:) :: unixnext
  double precision, allocatable, dimension(:) :: uniynext

  type, extends(MethodType) :: MethodCellTernaryType
    private
    real(DP), allocatable :: x_vert(:)
    real(DP), allocatable :: y_vert(:) !< cell vertex coordinates
    real(DP), allocatable :: v0x_vert_polygon(:)
    real(DP), allocatable :: v0y_vert_polygon(:)
    real(DP), allocatable :: v1x_vert_polygon(:)
    real(DP), allocatable :: v1y_vert_polygon(:) !< cell vertex velocities
    real(DP) :: xctr
    real(DP) :: yctr !< cell center coordinates
    real(DP) :: vxctr
    real(DP) :: vyctr !< cell center velocities
    real(DP) :: ztop
    real(DP) :: zbot !< cell top and bottom elevations
    real(DP) :: dz !< cell thickness
    real(DP) :: vztop
    real(DP) :: vzbot !< cell top and bottom velocities
    integer(I4B), public, pointer :: zeromethod
  contains
    procedure, public :: apply => apply_mct
    procedure, public :: destroy => destroy_mct
    procedure, public :: load => load_mct
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mct
    procedure         :: vertvelo_orig
    procedure         :: vertvelo
    !!!procedure         :: areapoly
    !!!procedure         :: hcsum
    !!!procedure         :: root_linear
  end type MethodCellTernaryType

contains

  !> @brief Create a tracking method
  subroutine create_method_cell_ternary(method)
    ! -- dummy
    type(MethodCellTernaryType), pointer :: method
    ! -- local
    type(CellPolyType), pointer :: cell
    type(SubcellTriType), pointer :: subcell

    allocate (method)
    allocate (method%zeromethod)
    call create_cell_poly(cell)
    method%cell => cell
    method%type => method%cell%type
    method%delegates = .true.
    call create_subcell_tri(subcell)
    method%subcell => subcell
    method%zeromethod = 0
  end subroutine create_method_cell_ternary

  !> @brief Destroy the tracking method
  subroutine destroy_mct(this)
    class(MethodCellTernaryType), intent(inout) :: this
    deallocate (this%type)
  end subroutine destroy_mct

  !> @brief Load subcell into tracking method
  subroutine load_mct(this, particle, next_level, submethod)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B), intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod

    select type (subcell => this%subcell)
    type is (SubcellTriType)
      call this%load_subcell(particle, subcell)
    end select
    call method_subcell_tern%init( &
      cell=this%cell, &
      subcell=this%subcell, &
      trackfilectl=this%trackfilectl, &
      tracktimes=this%tracktimes)
    submethod => method_subcell_tern
    method_subcell_tern%zeromethod = this%zeromethod
  end subroutine load_mct

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mct(this, particle)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer(I4B) :: isc
    integer(I4B) :: exitFace
    integer(I4B) :: inface
    integer(I4B) :: npolyverts

    exitFace = particle%iboundary(3)
    isc = particle%idomain(3)
    select type (cell => this%cell)
    type is (CellPolyType)
      npolyverts = cell%defn%npolyverts
    end select

    select case (exitFace)
    case (0)
      ! -- Subcell interior (cell interior)
      inface = -1
    case (1)
      ! -- Subcell face 1 (cell face)
      inface = isc
      if (inface .eq. 0) inface = npolyverts
    case (2)
      ! -- Subcell face --> next subcell in "cycle" (cell interior)
      isc = isc + 1
      if (isc .gt. npolyverts) isc = 1
      particle%idomain(3) = isc
      particle%iboundary(3) = 3
      inface = 0
    case (3)
      ! -- Subcell face --> preceding subcell in "cycle" (cell interior)
      isc = isc - 1
      if (isc .lt. 1) isc = npolyverts
      particle%idomain(3) = isc
      particle%iboundary(3) = 2
      inface = 0
    case (4)
      ! -- Subcell bottom (cell bottom)
      inface = npolyverts + 2
    case (5)
      ! -- Subcell top (cell top)
      inface = npolyverts + 3
    end select
    if (inface .eq. -1) then
      particle%iboundary(2) = 0
    else if (inface .eq. 0) then
      particle%iboundary(2) = 0
    else
      particle%iboundary(2) = inface
    end if
  end subroutine pass_mct

  !> @brief Apply the ternary method to a polygonal cell
  subroutine apply_mct(this, particle, tmax)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local

    ! -- Update particle state, checking whether any reporting or
    ! -- termination conditions apply
    call this%update(particle, this%cell%defn)

    ! -- Return early if particle is done advancing
    if (.not. particle%advancing) return

    ! -- If the particle is above the top of the cell (presumed water table)
    ! -- pass it vertically and instantaneously to the cell top and save the
    ! -- particle state to file
    if (particle%z > this%cell%defn%top) then
      particle%z = this%cell%defn%top
      call this%save(particle, reason=1) ! reason=1: cell transition
    end if
    
    ! -- Calculate vertex velocities
    if (particle%ivvorig > 0) then
      call this%vertvelo_orig()   ! kluge note: devoption for now
    else
      call this%vertvelo()
    end if

    ! -- Track across subcells
    call this%track(particle, 2, tmax) ! kluge, hardwired to level 2

  end subroutine apply_mct

  !> @brief Loads a triangular subcell from the polygonal cell
  subroutine load_subcell(this, particle, subcell)
    ! -- modules
    use ParticleModule, only: get_particle_id
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(SubcellTriType), intent(inout) :: subcell
    ! -- local
    integer(I4B) :: ic
    integer(I4B) :: isc
    integer(I4B) :: npolyverts
    integer(I4B) :: iv0
    integer(I4B) :: iv1
    integer(I4B) :: ipv0
    integer(I4B) :: ipv1
    integer(I4B) :: iv
    real(DP) :: x0
    real(DP) :: y0
    real(DP) :: x1
    real(DP) :: y1
    real(DP) :: x2
    real(DP) :: y2
    real(DP) :: x1rel
    real(DP) :: y1rel
    real(DP) :: x2rel
    real(DP) :: y2rel
    real(DP) :: xi
    real(DP) :: yi
    real(DP) :: di2
    real(DP) :: d02
    real(DP) :: d12
    real(DP) :: di1
    real(DP) :: d01
    real(DP) :: alphai
    real(DP) :: betai
    real(DP) :: betatol

    select type (cell => this%cell)
    type is (CellPolyType)
      ic = cell%defn%icell
      subcell%icell = ic
      isc = particle%idomain(3)
      npolyverts = cell%defn%npolyverts
      if (isc .le. 0) then
        xi = particle%x
        yi = particle%y
        do iv = 1, npolyverts
          iv0 = iv
          iv1 = iv + 1
          if (iv1 .gt. npolyverts) iv1 = 1
          ipv0 = iv0
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
          ! kluge note: can iboundary(2) be used to identify the subcell?
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
          print *, "error -- initial triangle not found for particle ", &
            get_particle_id(particle), " in cell ", ic
          call pstop(1)
        else
          ! subcellTri%isubcell = isc
          ! kluge note: as a matter of form, do we want to allow
          ! this subroutine to modify the particle???
          particle%idomain(3) = isc
        end if
      end if
      subcell%isubcell = isc

      ! -- Set coordinates and velocities at vertices of triangular subcell
      iv0 = isc
      iv1 = isc + 1
      if (iv1 .gt. npolyverts) iv1 = 1
      ipv0 = iv0
      ipv1 = iv1
      subcell%x0 = this%x_vert(ipv0)
      subcell%y0 = this%y_vert(ipv0)
      subcell%x1 = this%x_vert(ipv1)
      subcell%y1 = this%y_vert(ipv1)
      subcell%x2 = this%xctr
      subcell%y2 = this%yctr
      subcell%v0x = this%v0x_vert_polygon(iv0)
      subcell%v0y = this%v0y_vert_polygon(iv0)
      subcell%v1x = this%v1x_vert_polygon(iv0)   ! kluge note: the indices here actually refer to subcells, not vertices
      subcell%v1y = this%v1y_vert_polygon(iv0)
      subcell%v2x = this%vxctr
      subcell%v2y = this%vyctr
      subcell%ztop = this%ztop
      subcell%zbot = this%zbot
      subcell%dz = this%dz
      subcell%vzbot = this%vzbot
      subcell%vztop = this%vztop
    end select
  end subroutine load_subcell

  !> @brief Calculate vertex velocities the original way
  subroutine vertvelo_orig(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    integer(I4B) :: npolyverts
    integer(I4B) :: iv
    integer(I4B) :: ivp1
    integer(I4B) :: ivm1
    real(DP) :: retfactor
    real(DP) :: x0
    real(DP) :: y0
    real(DP) :: x1
    real(DP) :: y1
    real(DP) :: x2
    real(DP) :: y2
    real(DP) :: xsum
    real(DP) :: ysum
    real(DP) :: vxsum
    real(DP) :: vysum
    real(DP) :: flow0
    real(DP) :: flow1
    real(DP) :: v0x
    real(DP) :: v0y
    real(DP) :: d01x
    real(DP) :: d01y
    real(DP) :: d02x
    real(DP) :: d02y
    real(DP) :: det
    real(DP) :: area
    real(DP) :: term

    select type (cell => this%cell)
    type is (CellPolyType)

      npolyverts = cell%defn%npolyverts
      if (allocated(this%x_vert)) then
        deallocate (this%x_vert)
        deallocate (this%y_vert)
        deallocate (this%v0x_vert_polygon)
        deallocate (this%v0y_vert_polygon)
        deallocate (this%v1x_vert_polygon)
        deallocate (this%v1y_vert_polygon)
      end if
      allocate (this%x_vert(npolyverts))
      allocate (this%y_vert(npolyverts))
      allocate (this%v0x_vert_polygon(npolyverts))
      allocate (this%v0y_vert_polygon(npolyverts))
      allocate (this%v1x_vert_polygon(npolyverts))
      allocate (this%v1y_vert_polygon(npolyverts))

      xsum = DZERO
      ysum = DZERO
      vxsum = DZERO
      vysum = DZERO
      area = DZERO
      this%ztop = cell%defn%top
      this%zbot = cell%defn%bot
      this%dz = this%ztop - this%zbot
      do iv = 1, npolyverts
        ivp1 = iv + 1
        if (ivp1 .gt. npolyverts) ivp1 = 1
        ivm1 = iv - 1
        if (ivm1 .lt. 1) ivm1 = npolyverts
        x0 = cell%defn%polyvert(1, iv)
        y0 = cell%defn%polyvert(2, iv)
        x2 = cell%defn%polyvert(1, ivp1)
        y2 = cell%defn%polyvert(2, ivp1)
        x1 = cell%defn%polyvert(1, ivm1)
        y1 = cell%defn%polyvert(2, ivm1)
        term = DONE / (cell%defn%porosity * this%dz)
        flow0 = cell%defn%faceflow(iv) * term
        flow1 = cell%defn%faceflow(ivm1) * term
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
        retfactor = cell%defn%retfactor
        ! kluge note: can det ever be zero, like maybe for a 180-deg vertex???
        ! term = velfactor/det
        ! kluge note: can det ever be zero, like maybe for a 180-deg vertex???
        term = DONE / (retfactor * det)
        ! kluge note: "flow" here is volumetric flow rate (MODFLOW face flow)
        v0x = -term * (d02x * flow1 + d01x * flow0)
        ! per unit thickness, divided by porosity
        v0y = -term * (d02y * flow1 + d01y * flow0)
        this%v0x_vert_polygon(iv) = v0x
        this%v0y_vert_polygon(iv) = v0y
        this%v1x_vert_polygon(ivm1) = v0x   ! kluge note: the indices here actually refer to subcells, not vertices
        this%v1y_vert_polygon(ivm1) = v0y
        xsum = xsum + x0
        ysum = ysum + y0
        vxsum = vxsum + v0x
        vysum = vysum + v0y
        this%x_vert(iv) = x0
        this%y_vert(iv) = y0
        area = area + x0 * y1 - x1 * y0
      end do
      area = area * DHALF
      term = DONE / (retfactor * cell%defn%porosity * area)
      this%vzbot = cell%defn%faceflow(npolyverts + 2) * term
      this%vztop = -cell%defn%faceflow(npolyverts + 3) * term
      this%xctr = xsum / dble(npolyverts)
      this%yctr = ysum / dble(npolyverts)
      this%vxctr = vxsum / dble(npolyverts)
      this%vyctr = vysum / dble(npolyverts)

    end select
  end subroutine vertvelo_orig
  
  !> @brief Calculate vertex velocities
  subroutine vertvelo(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    integer(I4B) :: npolyverts
    integer(I4B) :: iv, ivm1
    real(DP) :: area
    real(DP) :: term

    integer                                  :: i
    double precision                         :: perturb
    double precision, allocatable, dimension(:)  :: xvals
    double precision, allocatable, dimension(:)  :: yvals
    double precision                         :: sixa
    double precision                         :: vm0i1
    double precision                         :: hcsumval
    double precision, allocatable, dimension(:) :: dx
    double precision, allocatable, dimension(:) :: dy
    double precision, allocatable, dimension(:) :: fact

    select type (cell => this%cell)
    type is (CellPolyType)

      npolyverts = cell%defn%npolyverts
      if (allocated(this%x_vert)) then
        deallocate (this%x_vert)
        deallocate (this%y_vert)
        deallocate (this%v0x_vert_polygon)
        deallocate (this%v0y_vert_polygon)
        deallocate (this%v1x_vert_polygon)
        deallocate (this%v1y_vert_polygon)
      end if
      allocate (this%x_vert(npolyverts))
      allocate (this%y_vert(npolyverts))
      allocate (this%v0x_vert_polygon(npolyverts))
      allocate (this%v0y_vert_polygon(npolyverts))
      allocate (this%v1x_vert_polygon(npolyverts))
      allocate (this%v1y_vert_polygon(npolyverts))
      
      nverts = npolyverts       ! kluge
      if (allocated(xvert)) then
        deallocate (xvert)
        deallocate (yvert)        
        deallocate (vne)
        deallocate (le)
        deallocate (unex)
        deallocate (uney)
        deallocate (areasub)
        deallocate (li)
        deallocate (unix)
        deallocate (uniy)
        deallocate (xmid)
        deallocate (ymid)
        deallocate (lm)
        deallocate (umx)
        deallocate (umy)
        deallocate (kappax)
        deallocate (kappay)
        deallocate (vm0i)
        deallocate (vm0e)
        deallocate (vm1i)
        deallocate (vm1e)
        deallocate (vm0x)
        deallocate (vm0y)
        deallocate (vm1x)
        deallocate (vm1y)
        deallocate (vc0x)
        deallocate (vc0y)
        deallocate (vc1x)
        deallocate (vc1y)
        deallocate (iprev)
        deallocate (xvertnext)
        deallocate (yvertnext)
        deallocate (unixnext)
        deallocate (uniynext)
      end if
      allocate (xvert(nverts))        ! x coordinates of cell vertices
      allocate (yvert(nverts))        ! y coordinates of cell vertices
      allocate (vne(npolyverts))          ! cell-edge normal velocities
      allocate (le(npolyverts))           ! lengths of exterior (cell) edges
      allocate (unex(npolyverts))         ! x components of unit normals to exterior edges
      allocate (uney(npolyverts))         ! y components of unit normals to exterior edges
      allocate (areasub(npolyverts))      ! subcell areas
      allocate (li(npolyverts))           ! lengths of interior edges ("spokes")
      allocate (unix(npolyverts))         ! x components of unit normals to interior edges
      allocate (uniy(npolyverts))         ! y components of unit normals to interior edges
      allocate (xmid(npolyverts))         ! x coordinates of midpoints
      allocate (ymid(npolyverts))         ! y coordinates of midpoints
      allocate (lm(npolyverts))           ! lengths of midpoint connectors
      allocate (umx(npolyverts))          ! x components of midpoint-connector (ccw) unit vectors
      allocate (umy(npolyverts))          ! y components of midpoint-connector (ccw) unit vectors
      allocate (kappax(npolyverts))       ! x components of kappa vectors
      allocate (kappay(npolyverts))       ! y components of kappa vectors
      ! Each subcell has two midpoint velocities (one on each interior edge): vm0 and vm1
      allocate (vm0i(npolyverts))         ! component of vm0 normal to the interior edge it's on
      allocate (vm0e(npolyverts))         ! component of vm0 in the direction normal to the corresponding exterior edge
      allocate (vm1i(npolyverts))         ! component of vm1 normal to the interior edge it's on
      allocate (vm1e(npolyverts))         ! component of vm1 in the direction normal to the corresponding exterior edge
      allocate (vm0x(npolyverts))         ! x component of vm0
      allocate (vm0y(npolyverts))         ! y component of vm0
      allocate (vm1x(npolyverts))         ! x component of vm1
      allocate (vm1y(npolyverts))         ! y component of vm1
      ! Each cell vertex (exterior corner) has two velocities: vc0 and vc1
      allocate (vc0x(npolyverts))         ! x component of vc0
      allocate (vc0y(npolyverts))         ! y component of vc0
      allocate (vc1x(npolyverts))         ! x component of vc1
      allocate (vc1y(npolyverts))         ! y component of vc1
      allocate (iprev(npolyverts))        ! vector of "previous" vextex indices defined for convenience
      allocate (xvertnext(npolyverts))    ! vector of "next" vertex x coordinates defined for convenience
      allocate (yvertnext(npolyverts))    ! vector of "next" vertex y coordinates defined for convenience
      allocate (unixnext(npolyverts))     ! vector of "next" interior unit-normal x coordinates defined for convenience
      allocate (uniynext(npolyverts))     ! vector of "next" interior unit-normal y coordinates defined for convenience
      
      do iv = 1, npolyverts
        this%x_vert(iv) = cell%defn%polyvert(1, iv)
        this%y_vert(iv) = cell%defn%polyvert(2, iv)
      end do

      xvert = this%x_vert
      yvert = this%y_vert

      this%ztop = cell%defn%top
      this%zbot = cell%defn%bot
      this%dz = this%ztop - this%zbot

      ! Shifted arrays for convenience
      do i = 1, nverts
          iprev(i) = i
      end do
      iprev = cshift(iprev, -1)
      xvertnext = cshift(xvert, 1)
      yvertnext = cshift(yvert, 1)
  
      ! Exterior edge unit normals (outward) and lengths
      allocate(dx(nverts), dy(nverts))
      dx = xvertnext - xvert
      dy = yvertnext - yvert
      le = dsqrt(dx * dx + dy * dy)
      unex = dy / le
      uney = -dx / le
      deallocate(dx, dy)
  
      ! Cell area
      areacell = areapoly(xvert, yvert)
      
      ! Cell centroid   ! kluge note: in general, this is NOT the average of the vertex coordinates
      sixa = areacell * 6.d0
      allocate(fact(nverts))
      fact = xvert * yvertnext - xvertnext * yvert
      xctrd = sum((xvert + xvertnext) * fact) / sixa
      yctrd = sum((yvert + yvertnext) * fact) / sixa
      deallocate(fact)
  
      ! Subcell areas
      allocate (xvals(3), yvals(3))
      do i = 1, nverts
          xvals(1) = xvert(i)
          xvals(2) = xvertnext(i)
          xvals(3) = xctrd
          yvals(1) = yvert(i)
          yvals(2) = yvertnext(i)
          yvals(3) = yctrd
          areasub(i) = areapoly(xvals, yvals)
      end do
      deallocate(xvals, yvals)
      
      ! Cell-edge normal velocities
      term = DONE / (cell%defn%porosity * cell%defn%retfactor)
      do i = 1, nverts
        area = le(i) * this%dz
        vne(i) = cell%defn%faceflow(i) * term / area
      end do
  
      ! Cell divergence (2D)
      divcell = sum(le * vne) / areacell
  
      ! Interior edge (ccw) unit normals and lengths
      dx = xvert - xctrd
      dy = yvert - yctrd
      li = dsqrt(dx * dx + dy * dy)
      unix = -dy / li
      uniy = dx / li
      ! Shifted arrays for convenience
      unixnext = cshift(unix, 1)
      uniynext = cshift(uniy, 1)
  
      ! Midpoints of interior edges
      xmid = 5.d-1 * (xvert + xctrd)
      ymid = 5.d-1 * (yvert + yctrd)
      
      ! Unit midpoint-connector (ccw) vectors and lengths
      dx = cshift(xmid, 1) - xmid
      dy = cshift(ymid, 1) - ymid
      lm = dsqrt(dx * dx + dy * dy)
      umx = dx / lm
      umy = dy / lm
  
      ! Kappa vectors (K tensor times unit midpoint-connector vectors)
      kappax = umx   ! kluge (isotropic K=1.)
      kappay = umy   ! kluge (isotropic K=1.)
  
      ! Use linearity to find vm0i[0] such that curl of the head gradient
      ! is zero
      vm0i1 = 0.d0
      perturb = 1.d-2
      vm0i1 = root_linear(hcsum, vm0i1, perturb)
      ! Evaluate head-cycle summation once more for its side effects
      hcsumval =  hcsum(vm0i1)
          
      ! Project linearly to get corner (vertex) velocities.
      vc0x = 2.d0 * vm0x - vctrdx
      vc0y = 2.d0 * vm0y - vctrdy
      vc1x = 2.d0 * vm1x - vctrdx
      vc1y = 2.d0 * vm1y - vctrdy
      
      do i = 1, nverts
        ivm1 = i - 1
        if (ivm1 .lt. 1) ivm1 = npolyverts
        this%v0x_vert_polygon(i) = vc0x(i)
        this%v0y_vert_polygon(i) = vc0y(i)
        this%v1x_vert_polygon(ivm1) = vc1x(ivm1)   ! kluge note: the indices here actually refer to subcells, not vertices
        this%v1y_vert_polygon(ivm1) = vc1y(ivm1)
      end do

      term = DONE / (cell%defn%retfactor * cell%defn%porosity * areacell)
      this%vzbot = cell%defn%faceflow(npolyverts + 2) * term
      this%vztop = -cell%defn%faceflow(npolyverts + 3) * term
      this%xctr = xctrd
      this%yctr = yctrd
      this%vxctr = vctrdx
      this%vyctr = vctrdy
 
    end select
  end subroutine vertvelo

function areapoly(xv, yv) result(area)
    ! dummy
    !!!class(MethodCellTernaryType), intent(inout) :: this
    double precision, dimension(:) :: xv
    double precision, dimension(:) :: yv
    ! result
    double precision               :: area

    area = 5.d-1 * sum(xv(:) * cshift(yv(:), 1) - cshift(xv(:), 1) * yv(:))

end function areapoly

double precision function hcsum(vm0i1)
    ! dummy
    !!!class(MethodCellTernaryType), intent(inout) :: this
    double precision vm0i1
    ! local
    double precision uprod(nverts), det(nverts), wt(nverts)
    double precision bi0x(nverts), be0x(nverts), bi0y(nverts), be0y(nverts)
    double precision bi1x(nverts), be1x(nverts), bi1y(nverts), be1y(nverts)
    double precision be01x(nverts), be01y(nverts)
    double precision emxx, emxy, emyx, emyy, rx, ry, emdet
    integer i, ip
    
    ! Set vm0i(1)
    vm0i(1) = vm0i1
    
    ! Get remaining vm0i's sequentially using divergence conditions
    do i = 2, nverts
        ip = iprev(i)
        vm0i(i) = (li(ip) * vm0i(ip) - le(ip) * vne(ip) + areasub(ip) * divcell) / li(i)
    end do

    ! Get vm1i's from vm0i's using continuity conditions
    vm1i = cshift(vm0i, 1)

    ! Get centroid velocity by setting up and solving 2x2 linear system
    uprod = unix * unex + uniy * uney
    det = 1.d0 - uprod * uprod
    bi0x = (unix - unex * uprod) / det
    be0x = (unex - unix * uprod) / det
    bi0y = (uniy - uney * uprod) / det
    be0y = (uney - uniy * uprod) / det
    uprod = unixnext * unex + uniynext * uney
    det = 1.d0 - uprod * uprod
    bi1x = (unixnext - unex * uprod) / det
    be1x = (unex - unixnext * uprod) / det
    bi1y = (uniynext - uney * uprod) / det
    be1y = (uney - uniynext * uprod) / det
    be01x = 5.d-1 * (be0x + be1x)
    be01y = 5.d-1 * (be0y + be1y)
    wt = 1.d0 / dble(nverts)      ! kluge (equal weights)
    emxx = 2.d0 - sum(wt * be01x * unex)
    emxy = -sum(wt * be01x * uney)
    emyx = -sum(wt * be01y * unex)
    emyy = 2.d0 - sum(wt * be01y * uney)
    rx = sum(wt * (bi0x * vm0i + bi1x * vm1i + be01x * vne))
    ry = sum(wt * (bi0y * vm0i + bi1y * vm1i + be01y * vne))
    emdet = emxx * emyy - emxy * emyx
    vctrdx = (emyy * rx - emxy * ry) / emdet
    vctrdy = (emxx * ry - emyx * rx) / emdet
    
    ! Get vm0e's using "known" conditions
    vm0e = 5.d-1 * (vne + unex * vctrdx + uney * vctrdy)
    
    ! Get vm1e's from uniformity along exterior edges
    vm1e = vm0e
    
    ! Transform vm0 and vm1 to (x, y) coordinates
    vm0x = bi0x * vm0i + be0x * vm0e
    vm0y = bi0y * vm0i + be0y * vm0e
    vm1x = bi1x * vm1i + be1x * vm0e
    vm1y = bi1y * vm1i + be1y * vm0e
    
    ! Calculate head-cycle summation (which is proportional to
    ! the curl of the head gradient)
    hcsum = sum(lm * (kappax * (vm0x + vm1x) + kappay * (vm0y + vm1y)))
    
    return

end function hcsum

function root_linear(func, x0, perturb) result(x)
    ! dummy
    !!!class(MethodCellTernaryType), intent(inout) :: this
    double precision func, x0, perturb
    external func   ! kluge
    ! result
    double precision x
    ! local
    double precision f0, f, jac

    ! Find the root of a linear function (func) that is complex
    ! to evaluate and would be difficult to invert explicitly
    
    ! Base evaluation of function
    f0 = func(x0)
    
    ! Calculate Jacobian (derivative)
    x = x0 + perturb
    f = func(x)
    jac = (f - f0) / perturb
    
    ! Solve by Newton-like linear projection
    x = x0 - f0 / jac
        
    return

end function root_linear

end module MethodCellTernaryModule

