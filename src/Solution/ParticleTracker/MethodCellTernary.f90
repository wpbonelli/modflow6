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
  use ArrayHandlersModule, only: ExpandArray
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_method_cell_ternary

  type, extends(MethodType) :: MethodCellTernaryType
    private
    integer(I4B), public, pointer :: zeromethod !< root-finding method to use
    integer(I4B) :: nverts !< number of vertices
    real(DP) :: xctr !< cell center x coord
    real(DP) :: yctr !< cell center y coord
    real(DP) :: vctrx !< cell center velocity x component
    real(DP) :: vctry !< cell center velocity y component
    real(DP) :: ztop !< cell top elevation
    real(DP) :: zbot !< cell bottom elevation
    real(DP) :: dz !< cell thickness
    real(DP) :: vztop !< cell top velocity
    real(DP) :: vzbot !< cell bottom velocity
    integer(I4B), allocatable, dimension(:) :: iprev !< shifted indices
    real(DP), allocatable, dimension(:) :: xvert !< cell vertex x coords
    real(DP), allocatable, dimension(:) :: yvert !< cell vertex y coords
    real(DP), allocatable, dimension(:) :: xvertnext !< shifted cell vertex x coords
    real(DP), allocatable, dimension(:) :: yvertnext !< shifted cell vertex y coords
    real(DP), allocatable, dimension(:) :: vne !< cell edge normal velocities
    real(DP), allocatable, dimension(:) :: vv0x
    real(DP), allocatable, dimension(:) :: vv0y
    real(DP), allocatable, dimension(:) :: vv1x
    real(DP), allocatable, dimension(:) :: vv1y !< cell vertex velocities
    ! below moved from vertvelo
    real(DP), allocatable, dimension(:) :: xvals
    real(DP), allocatable, dimension(:) :: yvals
    real(DP), allocatable, dimension(:) :: wk1
    real(DP), allocatable, dimension(:) :: wk2
    real(DP), allocatable, dimension(:) :: le !< lengths of exterior (cell) edges
    real(DP), allocatable, dimension(:) :: li !< lengths of interior edges ("spokes")
    real(DP), allocatable, dimension(:) :: unex !< x components of unit normals to exterior edges
    real(DP), allocatable, dimension(:) :: uney !< y components of unit normals to exterior edges
    real(DP), allocatable, dimension(:) :: areasub !< subcell areas
    real(DP), allocatable, dimension(:) :: unix !< x components of unit normals to interior edges
    real(DP), allocatable, dimension(:) :: uniy !< y components of unit normals to interior edges
    real(DP), allocatable, dimension(:) :: unixnext !< vector of "next" interior unit-normal x coords defined for convenience
    real(DP), allocatable, dimension(:) :: uniynext !< vector of "next" interior unit-normal y coords defined for convenience
    real(DP), allocatable, dimension(:) :: xmid !< x coords of midpoints
    real(DP), allocatable, dimension(:) :: ymid !< y coords of midpoints
    real(DP), allocatable, dimension(:) :: lm !< lengths of midpoint connectors
    real(DP), allocatable, dimension(:) :: umx !< x components of midpoint-connector (ccw) unit vectors
    real(DP), allocatable, dimension(:) :: umy !< y components of midpoint-connector (ccw) unit vectors
    real(DP), allocatable, dimension(:) :: kappax !< x components of kappa vectors
    real(DP), allocatable, dimension(:) :: kappay !< y components of kappa vectors
    real(DP), allocatable, dimension(:) :: vm0x !< x component of vm0
    real(DP), allocatable, dimension(:) :: vm0y !< y component of vm0
    real(DP), allocatable, dimension(:) :: vm1x !< x component of vm1
    real(DP), allocatable, dimension(:) :: vm1y !< y component of vm1
    ! below moved from hcsum
    real(DP), allocatable, dimension(:) :: vm0i
    real(DP), allocatable, dimension(:) :: vm0e
    real(DP), allocatable, dimension(:) :: vm1i
    real(DP), allocatable, dimension(:) :: vm1e
    real(DP), allocatable, dimension(:) :: uprod
    real(DP), allocatable, dimension(:) :: det
    real(DP), allocatable, dimension(:) :: wt
    real(DP), allocatable, dimension(:) :: bi0x
    real(DP), allocatable, dimension(:) :: be0x
    real(DP), allocatable, dimension(:) :: bi0y
    real(DP), allocatable, dimension(:) :: be0y
    real(DP), allocatable, dimension(:) :: bi1x
    real(DP), allocatable, dimension(:) :: be1x
    real(DP), allocatable, dimension(:) :: bi1y
    real(DP), allocatable, dimension(:) :: be1y
    real(DP), allocatable, dimension(:) :: be01x
    real(DP), allocatable, dimension(:) :: be01y
  contains
    procedure, public :: apply => apply_mct
    procedure, public :: destroy => destroy_mct
    procedure, public :: load => load_mct
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mct
    procedure :: vertvelo_orig
    procedure :: vertvelo
    procedure :: calc_thru_hcsum
  end type MethodCellTernaryType

contains

  !> @brief Create a tracking method
  subroutine create_method_cell_ternary(method)
    ! dummy
    type(MethodCellTernaryType), pointer :: method
    ! local
    type(CellPolyType), pointer :: cell
    type(SubcellTriType), pointer :: subcell

    allocate (method)
    allocate (method%zeromethod)
    allocate (method%xvals(3))
    allocate (method%yvals(3))
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
    ! dummy
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
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer(I4B) :: isc
    integer(I4B) :: exitFace
    integer(I4B) :: inface

    exitFace = particle%iboundary(3)
    isc = particle%idomain(3)  ! subcell id

    select case (exitFace)
    case (0)
      ! Subcell interior (cell interior)
      inface = -1
    case (1)
      ! Subcell face 1 (cell face)
      inface = isc
      if (inface .eq. 0) inface = this%nverts
    case (2)
      ! Subcell face --> next subcell in "cycle" (cell interior)
      isc = isc + 1
      if (isc .gt. this%nverts) isc = 1
      particle%idomain(3) = isc
      particle%iboundary(3) = 3
      inface = 0
    case (3)
      ! Subcell face --> preceding subcell in "cycle" (cell interior)
      isc = isc - 1
      if (isc .lt. 1) isc = this%nverts
      particle%idomain(3) = isc
      particle%iboundary(3) = 2
      inface = 0
    case (4)
      ! Subcell bottom (cell bottom)
      inface = this%nverts + 2
    case (5)
      ! Subcell top (cell top)
      inface = this%nverts + 3
    end select
    if (inface .eq. -1) then
      particle%iboundary(2) = 0
    else if (inface .eq. 0) then
      particle%iboundary(2) = 0
    else
      particle%iboundary(2) = inface
    end if

    print *, "passing from ", exitFace, " to ", inface
  end subroutine pass_mct

  !> @brief Apply the ternary method to a polygonal cell
  subroutine apply_mct(this, particle, tmax)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    integer(I4B) :: i
    integer(I4B) :: inc

    ! Update particle state, checking whether any reporting or
    ! termination conditions apply
    call this%update(particle, this%cell%defn)

    ! Return early if particle is done advancing
    if (.not. particle%advancing) return

    ! If the particle is above the top of the cell (presumed water table)
    ! pass it vertically and instantaneously to the cell top and save the
    ! particle state to file
    if (particle%z > this%cell%defn%top) then
      particle%z = this%cell%defn%top
      call this%save(particle, reason=1) ! reason=1: cell transition
    end if

    select type (cell => this%cell)
    type is (CellPolyType)
      ! Number of vertices
      this%nverts = cell%defn%npolyverts

      ! Expand arrays if needed
      if (size(this%xvert) < this%nverts) then
        inc = this%nverts - size(this%xvert)
        call ExpandArray(this%xvert, inc)
        call ExpandArray(this%yvert, inc)
        call ExpandArray(this%vne, inc)
        call ExpandArray(this%vv0x, inc)
        call ExpandArray(this%vv0y, inc)
        call ExpandArray(this%vv1x, inc)
        call ExpandArray(this%vv1y, inc)
        call ExpandArray(this%iprev, inc)
        call ExpandArray(this%xvertnext, inc)
        call ExpandArray(this%yvertnext, inc)
      end if

      ! Cell vertices
      do i = 1, this%nverts
        this%xvert(i) = cell%defn%polyvert(1, i)
        this%yvert(i) = cell%defn%polyvert(2, i)
      end do
      ! Top, bottom, and thickness
      this%ztop = cell%defn%top
      this%zbot = cell%defn%bot
      this%dz = this%ztop - this%zbot
      ! Shifted arrays
      do i = 1, this%nverts
        this%iprev(i) = i
      end do
      this%iprev = cshift(this%iprev, -1)
      this%xvertnext = cshift(this%xvert, 1)
      this%yvertnext = cshift(this%yvert, 1)
    end select

    ! Calculate vertex velocities
    if (particle%ivvorig > 0) then
      call this%vertvelo_orig() ! kluge note: devoption for now
    else
      print *, "calculating vertvelo"
      call this%vertvelo()
      print *, "calculated vertvelo"
    end if

    ! Track across subcells
    print *, "tracking over subcells"
    call this%track(particle, 2, tmax) ! kluge, hardwired to level 2
    print *, "tracked over subcells"

  end subroutine apply_mct

  !> @brief Loads a triangular subcell from the polygonal cell
  subroutine load_subcell(this, particle, subcell)
    ! modules
    use ParticleModule, only: get_particle_id
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(SubcellTriType), intent(inout) :: subcell
    ! local
    integer(I4B) :: ic
    integer(I4B) :: isc
    integer(I4B) :: iv0
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
      if (isc .le. 0) then
        xi = particle%x
        yi = particle%y
        do iv0 = 1, this%nverts
          x0 = this%xvert(iv0)
          y0 = this%yvert(iv0)
          x1 = this%xvertnext(iv0)
          y1 = this%yvertnext(iv0)
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
            isc = iv0 ! but maybe not!!!!!!!!!!!!
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

      ! Set coordinates and velocities at vertices of triangular subcell
      iv0 = isc
      subcell%x0 = this%xvert(iv0)
      subcell%y0 = this%yvert(iv0)
      subcell%x1 = this%xvertnext(iv0)
      subcell%y1 = this%yvertnext(iv0)
      subcell%x2 = this%xctr
      subcell%y2 = this%yctr
      subcell%v0x = this%vv0x(iv0)
      subcell%v0y = this%vv0y(iv0)
      subcell%v1x = this%vv1x(iv0) ! kluge note: the indices here actually refer to subcells, not vertices
      subcell%v1y = this%vv1y(iv0)
      subcell%v2x = this%vctrx
      subcell%v2y = this%vctry
      subcell%ztop = this%ztop
      subcell%zbot = this%zbot
      subcell%dz = this%dz
      subcell%vzbot = this%vzbot
      subcell%vztop = this%vztop
    end select
  end subroutine load_subcell

  !> @brief Calculate vertex velocities the original way    ! kluge note: temporary, for testing
  subroutine vertvelo_orig(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
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

      xsum = DZERO
      ysum = DZERO
      vxsum = DZERO
      vysum = DZERO
      area = DZERO
      this%ztop = cell%defn%top
      this%zbot = cell%defn%bot
      this%dz = this%ztop - this%zbot
      do iv = 1, this%nverts
        ivp1 = iv + 1
        if (ivp1 .gt. this%nverts) ivp1 = 1
        ivm1 = iv - 1
        if (ivm1 .lt. 1) ivm1 = this%nverts
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
        this%vv0x(iv) = v0x
        this%vv0y(iv) = v0y
        this%vv1x(ivm1) = v0x ! kluge note: the indices here actually refer to subcells, not vertices
        this%vv1y(ivm1) = v0y
        xsum = xsum + x0
        ysum = ysum + y0
        vxsum = vxsum + v0x
        vysum = vysum + v0y
        this%xvert(iv) = x0
        this%yvert(iv) = y0
        area = area + x0 * y1 - x1 * y0
      end do
      area = area * DHALF
      term = DONE / (retfactor * cell%defn%porosity * area)
      this%vzbot = cell%defn%faceflow(this%nverts + 2) * term
      this%vztop = -cell%defn%faceflow(this%nverts + 3) * term
      this%xctr = xsum / dble(this%nverts)
      this%yctr = ysum / dble(this%nverts)
      this%vctrx = vxsum / dble(this%nverts)
      this%vctry = vysum / dble(this%nverts)

    end select
  end subroutine vertvelo_orig

  !> @brief Calculate vertex velocities
  subroutine vertvelo(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    integer(I4B) :: iv
    integer(I4B) :: ivm1
    integer(I4B) :: i
    integer(I4B) :: inc
    real(DP) :: term
    real(DP) :: perturb
    real(DP) :: sixa
    real(DP) :: vm0i0
    real(DP) :: vm0ival
    real(DP) :: hcsum0
    real(DP) :: hcsum
    real(DP) :: jac
    real(DP) :: areacell
    real(DP) :: divcell

    select type (cell => this%cell)
    type is (CellPolyType)

      if (size(this%le) < this%nverts) then
        inc = this%nverts - size(this%le)
        call ExpandArray(this%le, inc)
        call ExpandArray(this%unex, inc)
        call ExpandArray(this%uney, inc)
        call ExpandArray(this%areasub, inc)
        call ExpandArray(this%li, inc)
        call ExpandArray(this%unix, inc)
        call ExpandArray(this%uniy, inc)
        call ExpandArray(this%xmid, inc)
        call ExpandArray(this%ymid, inc)
        call ExpandArray(this%lm, inc)
        call ExpandArray(this%umx, inc)
        call ExpandArray(this%umy, inc)
        call ExpandArray(this%kappax, inc)
        call ExpandArray(this%kappay, inc)
        call ExpandArray(this%vm0x, inc)
        call ExpandArray(this%vm0y, inc)
        call ExpandArray(this%vm1x, inc)
        call ExpandArray(this%vm1y, inc)
        call ExpandArray(this%unixnext, inc)
        call ExpandArray(this%uniynext, inc)
        call ExpandArray(this%wk1, inc)
        call ExpandArray(this%wk2, inc)
      end if

      ! Exterior edge unit normals (outward) and lengths
      this%wk1 = this%xvertnext - this%xvert
      this%wk2 = this%yvertnext - this%yvert
      this%le = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%unex = this%wk2 / this%le
      this%uney = -this%wk1 / this%le

      ! Cell area
      areacell = areapoly(this%xvert, this%yvert)

      ! Cell centroid   ! kluge note: in general, this is NOT the average of the vertex coordinates
      sixa = areacell * 6.d0
      this%wk1 = this%xvert * this%yvertnext - this%xvertnext * this%yvert
      this%xctr = sum((this%xvert + this%xvertnext) * this%wk1) / sixa
      this%yctr = sum((this%yvert + this%yvertnext) * this%wk1) / sixa

      ! Subcell areas
      do i = 1, this%nverts
        this%xvals(1) = this%xvert(i)
        this%xvals(2) = this%xvertnext(i)
        this%xvals(3) = this%xctr
        this%yvals(1) = this%yvert(i)
        this%yvals(2) = this%yvertnext(i)
        this%yvals(3) = this%yctr
        this%areasub(i) = areapoly(this%xvals, this%yvals)
      end do

      ! Cell-edge normal velocities
      term = DONE / (cell%defn%porosity * cell%defn%retfactor * this%dz)
      do i = 1, this%nverts
        this%vne(i) = cell%defn%faceflow(i) * term / this%le(i)
      end do

      ! Cell divergence (2D)
      divcell = sum(this%le * this%vne) / areacell

      ! Interior edge (ccw) unit normals and lengths
      this%wk1 = this%xvert - this%xctr
      this%wk2 = this%yvert - this%yctr
      this%li = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%unix = -this%wk2 / this%li
      this%uniy = this%wk1 / this%li
      ! Shifted arrays for convenience
      this%unixnext = cshift(this%unix, 1)
      this%uniynext = cshift(this%uniy, 1)

      ! Midpoints of interior edges
      this%xmid = 5.d-1 * (this%xvert + this%xctr)
      this%ymid = 5.d-1 * (this%yvert + this%yctr)

      ! Unit midpoint-connector (ccw) vectors and lengths
      this%wk1 = cshift(this%xmid, 1) - this%xmid
      this%wk2 = cshift(this%ymid, 1) - this%ymid
      this%lm = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%umx = this%wk1 / this%lm
      this%umy = this%wk2 / this%lm

      ! Kappa vectors (K tensor times unit midpoint-connector vectors)
      this%kappax = this%umx ! kluge (isotropic K=1.)
      this%kappay = this%umy ! kluge (isotropic K=1.)

      ! Use linearity to find vm0i[0] such that curl of the head gradient
      ! is zero
      perturb = 1.d-2 ! kluge?
      ! Calculations at base value
      vm0i0 = 0.d0
      call this%calc_thru_hcsum(vm0i0, divcell, hcsum0)
      ! Calculations at perturbed value
      vm0ival = vm0i0 + perturb
      call this%calc_thru_hcsum(vm0ival, divcell, hcsum)
      ! Calculations at root value
      jac = (hcsum - hcsum0) / perturb
      vm0ival = vm0i0 - hcsum0 / jac
      call this%calc_thru_hcsum(vm0ival, divcell, hcsum)

      ! Project linearly to get corner (vertex) velocities. Note that velocity
      ! vv1 is at the next vertex ccw from vv0, so vv0(i) and vv1(i) are the
      ! two vertex velocities used by triangular subcell i.
      this%vv0x = 2.d0 * this%vm0x - this%vctrx
      this%vv0y = 2.d0 * this%vm0y - this%vctry
      this%vv1x = 2.d0 * this%vm1x - this%vctrx
      this%vv1y = 2.d0 * this%vm1y - this%vctry

      ! Set top and bottom velocities
      term = DONE / (cell%defn%retfactor * cell%defn%porosity * areacell)
      this%vzbot = cell%defn%faceflow(this%nverts + 2) * term
      this%vztop = -cell%defn%faceflow(this%nverts + 3) * term

    end select
  end subroutine vertvelo

  subroutine calc_thru_hcsum(this, vm0ival, divcell, hcsum)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    real(DP) :: vm0ival
    real(DP) :: divcell
    real(DP) :: hcsum
    ! local
    integer(I4B) :: i
    integer(I4B) :: ip
    integer(I4B) :: inc
    real(DP) :: emxx
    real(DP) :: emxy
    real(DP) :: emyx
    real(DP) :: emyy
    real(DP) :: rx
    real(DP) :: ry
    real(DP) :: emdet

    if (size(this%vm0i) < this%nverts) then
      inc = this%nverts - size(this%vm0i)
      call ExpandArray(this%vm0i, inc)
      call ExpandArray(this%vm0e, inc)
      call ExpandArray(this%vm1i, inc)
      call ExpandArray(this%vm1e, inc)
      call ExpandArray(this%uprod, inc)
      call ExpandArray(this%det, inc)
      call ExpandArray(this%wt, inc)
      call ExpandArray(this%bi0x, inc)
      call ExpandArray(this%be0x, inc)
      call ExpandArray(this%bi0y, inc)
      call ExpandArray(this%be0y, inc)
      call ExpandArray(this%bi1x, inc)
      call ExpandArray(this%bi1y, inc)
      call ExpandArray(this%be1y, inc)
      call ExpandArray(this%be01x, inc)
      call ExpandArray(this%be01y, inc)
    end if

    ! Set vm0i(1)
    this%vm0i(1) = vm0ival

    ! Get remaining vm0i's sequentially using divergence conditions
    do i = 2, this%nverts
      ip = this%iprev(i)
      this%vm0i(i) = (this%li(ip) * this%vm0i(ip) - this%le(ip) * this%vne(ip) &
                      + this%areasub(ip) * divcell) / this%li(i)
    end do

    ! Get vm1i's from vm0i's using continuity conditions
    this%vm1i = cshift(this%vm0i, 1)

    ! Get centroid velocity by setting up and solving 2x2 linear system
    this%uprod = this%unix * this%unex + this%uniy * this%uney
    this%det = 1.d0 - this%uprod * this%uprod
    this%bi0x = (this%unix - this%unex * this%uprod) / this%det
    this%be0x = (this%unex - this%unix * this%uprod) / this%det
    this%bi0y = (this%uniy - this%uney * this%uprod) / this%det
    this%be0y = (this%uney - this%uniy * this%uprod) / this%det
    this%uprod = this%unixnext * this%unex + this%uniynext * this%uney
    this%det = 1.d0 - this%uprod * this%uprod
    this%bi1x = (this%unixnext - this%unex * this%uprod) / this%det
    this%be1x = (this%unex - this%unixnext * this%uprod) / this%det
    this%bi1y = (this%uniynext - this%uney * this%uprod) / this%det
    this%be1y = (this%uney - this%uniynext * this%uprod) / this%det
    this%be01x = 5.d-1 * (this%be0x + this%be1x)
    this%be01y = 5.d-1 * (this%be0y + this%be1y)
    this%wt = 1.d0 / dble(this%nverts) ! kluge (equal weights)
    emxx = 2.d0 - sum(this%wt * this%be01x * this%unex)
    emxy = -sum(this%wt * this%be01x * this%uney)
    emyx = -sum(this%wt * this%be01y * this%unex)
    emyy = 2.d0 - sum(this%wt * this%be01y * this%uney)
    rx = sum(this%wt * &
             (this%bi0x * this%vm0i + &
              this%bi1x * this%vm1i + &
              this%be01x * this%vne))
    ry = sum(this%wt * &
             (this%bi0y * this%vm0i + &
              this%bi1y * this%vm1i + &
              this%be01y * this%vne))
    emdet = emxx * emyy - emxy * emyx
    this%vctrx = (emyy * rx - emxy * ry) / emdet
    this%vctry = (emxx * ry - emyx * rx) / emdet

    ! Get vm0e's using "known" conditions
    this%vm0e = 5.d-1 * &
                (this%vne + this%unex * this%vctrx + this%uney * this%vctry)

    ! Get vm1e's from uniformity along exterior edges
    this%vm1e = this%vm0e

    ! Transform vm0 and vm1 to (x, y) coordinates
    this%vm0x = this%bi0x * this%vm0i + this%be0x * this%vm0e
    this%vm0y = this%bi0y * this%vm0i + this%be0y * this%vm0e
    this%vm1x = this%bi1x * this%vm1i + this%be1x * this%vm0e
    this%vm1y = this%bi1y * this%vm1i + this%be1y * this%vm0e

    ! Calculate head-cycle summation (which is proportional to
    ! the curl of the head gradient)
    hcsum = sum(this%lm * &
                (this%kappax * (this%vm0x + this%vm1x) + &
                 this%kappay * (this%vm0y + this%vm1y)))

  end subroutine calc_thru_hcsum

  function areapoly(xv, yv) result(area) ! kluge note: should this be packaged with other utilities?
    ! dummy
    double precision, dimension(:) :: xv
    double precision, dimension(:) :: yv
    ! result
    double precision :: area

    area = 5.d-1 * sum(xv(:) * cshift(yv(:), 1) - cshift(xv(:), 1) * yv(:))

  end function areapoly

end module MethodCellTernaryModule

