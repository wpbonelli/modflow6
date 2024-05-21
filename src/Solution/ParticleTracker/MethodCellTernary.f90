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
  use GeomUtilModule, only: area
  use ConstantsModule, only: DZERO, DONE, DTWO
  use ArrayHandlersModule, only: ExpandArray
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_method_cell_ternary

  type, extends(MethodType) :: MethodCellTernaryType
    private
    integer(I4B) :: nverts !< number of vertices
    real(DP), allocatable, dimension(:) :: xvert
    real(DP), allocatable, dimension(:) :: yvert !< cell vertex coordinates
    real(DP), allocatable, dimension(:) :: vne !< cell edge normal velocities
    real(DP), allocatable, dimension(:) :: vv0x
    real(DP), allocatable, dimension(:) :: vv0y
    real(DP), allocatable, dimension(:) :: vv1x
    real(DP), allocatable, dimension(:) :: vv1y !< cell vertex velocities
    real(DP) :: xctr
    real(DP) :: yctr !< cell center coordinates
    real(DP) :: vctrx
    real(DP) :: vctry !< cell center velocities
    real(DP) :: ztop
    real(DP) :: zbot !< cell top and bottom elevations
    real(DP) :: dz !< cell thickness
    real(DP) :: vztop
    real(DP) :: vzbot !< cell top and bottom velocities
    integer(I4B), allocatable, dimension(:) :: iprev !< array of shifted indices
    real(DP), allocatable, dimension(:) :: xvertnext
    real(DP), allocatable, dimension(:) :: yvertnext !< arrays of shifted cell vertex coordinates
    integer(I4B), public, pointer :: zeromethod
    ! vertvelo vars
    real(DP), allocatable, dimension(:) :: le
    real(DP), allocatable, dimension(:) :: unextx
    real(DP), allocatable, dimension(:) :: unexty
    real(DP), allocatable, dimension(:) :: areasub
    real(DP), allocatable, dimension(:) :: li
    real(DP), allocatable, dimension(:) :: unintx
    real(DP), allocatable, dimension(:) :: uninty
    real(DP), allocatable, dimension(:) :: xmid
    real(DP), allocatable, dimension(:) :: ymid
    real(DP), allocatable, dimension(:) :: lm
    real(DP), allocatable, dimension(:) :: umx
    real(DP), allocatable, dimension(:) :: umy
    real(DP), allocatable, dimension(:) :: kappax
    real(DP), allocatable, dimension(:) :: kappay
    real(DP), allocatable, dimension(:) :: vm0x
    real(DP), allocatable, dimension(:) :: vm0y
    real(DP), allocatable, dimension(:) :: vm1x
    real(DP), allocatable, dimension(:) :: vm1y
    real(DP), allocatable, dimension(:) :: unextxnext
    real(DP), allocatable, dimension(:) :: unextynext
    real(DP), allocatable, dimension(:) :: wk1
    real(DP), allocatable, dimension(:) :: wk2
    real(DP), allocatable, dimension(:) :: xvals
    real(DP), allocatable, dimension(:) :: yvals
    ! calc_thru_hsum vars
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
    procedure, public :: deallocate => destroy_mct
    procedure, public :: load => load_mct
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mct
    procedure :: vertvelo
    procedure :: calc_thru_hcsum
    procedure :: reallocate_vv
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
    ! vertvelo
    allocate (method%le(3)) ! lengths of exterior (cell) edges
    allocate (method%unextx(3)) ! x components of unit normals to exterior edges
    allocate (method%unexty(3)) ! y components of unit normals to exterior edges
    allocate (method%areasub(3)) ! subcell areas
    allocate (method%li(3)) ! lengths of interior edges ("spokes")
    allocate (method%unintx(3)) ! x components of unit normals to interior edges
    allocate (method%uninty(3)) ! y components of unit normals to interior edges
    allocate (method%xmid(3)) ! x coordinates of midpoints
    allocate (method%ymid(3)) ! y coordinates of midpoints
    allocate (method%lm(3)) ! lengths of midpoint connectors
    allocate (method%umx(3)) ! x components of midpoint-connector (cw) unit vectors
    allocate (method%umy(3)) ! y components of midpoint-connector (cw) unit vectors
    allocate (method%kappax(3)) ! x components of kappa vectors
    allocate (method%kappay(3)) ! y components of kappa vectors
    allocate (method%vm0x(3)) ! x component of vm0
    allocate (method%vm0y(3)) ! y component of vm0
    allocate (method%vm1x(3)) ! x component of vm1
    allocate (method%vm1y(3)) ! y component of vm1
    allocate (method%unextxnext(3)) ! vector of "next" interior unit-normal x coordinates defined for convenience
    allocate (method%unextynext(3)) ! vector of "next" interior unit-normal y coordinates defined for convenience
    allocate (method%wk1(3))
    allocate (method%wk2(3))
    allocate (method%xvals(3))
    allocate (method%yvals(3))
    ! calc_thru_hsum
    allocate (method%vm0i(3))
    allocate (method%vm0e(3))
    allocate (method%vm1i(3))
    allocate (method%vm1e(3))
    allocate (method%uprod(3))
    allocate (method%det(3))
    allocate (method%wt(3))
    allocate (method%bi0x(3))
    allocate (method%be0x(3))
    allocate (method%bi0y(3))
    allocate (method%be0y(3))
    allocate (method%bi1x(3))
    allocate (method%be1x(3))
    allocate (method%bi1y(3))
    allocate (method%be1y(3))
    allocate (method%be01x(3))
    allocate (method%be01y(3))
    call create_cell_poly(cell)
    method%cell => cell
    method%type => method%cell%type
    method%delegates = .true.
    call create_subcell_tri(subcell)
    method%subcell => subcell
  end subroutine create_method_cell_ternary

  !> @brief Destroy the tracking method
  subroutine destroy_mct(this)
    class(MethodCellTernaryType), intent(inout) :: this
    deallocate (this%type)
    ! vertvelo vars
    deallocate (this%le)
    deallocate (this%unextx)
    deallocate (this%unexty)
    deallocate (this%areasub)
    deallocate (this%li)
    deallocate (this%unintx)
    deallocate (this%uninty)
    deallocate (this%xmid)
    deallocate (this%ymid)
    deallocate (this%lm)
    deallocate (this%umx)
    deallocate (this%umy)
    deallocate (this%kappax)
    deallocate (this%kappay)
    deallocate (this%vm0x)
    deallocate (this%vm0y)
    deallocate (this%vm1x)
    deallocate (this%vm1y)
    deallocate (this%unextxnext)
    deallocate (this%unextynext)
    deallocate (this%wk1)
    deallocate (this%wk2)
    deallocate (this%xvals)
    deallocate (this%yvals)
    ! calc_thru_hsum vars
    deallocate (this%vm0i)
    deallocate (this%vm0e)
    deallocate (this%vm1i)
    deallocate (this%vm1e)
    deallocate (this%uprod)
    deallocate (this%det)
    deallocate (this%wt)
    deallocate (this%bi0x)
    deallocate (this%be0x)
    deallocate (this%bi0y)
    deallocate (this%be0y)
    deallocate (this%bi1x)
    deallocate (this%be1x)
    deallocate (this%bi1y)
    deallocate (this%be1y)
    deallocate (this%be01x)
    deallocate (this%be01y)
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
    isc = particle%idomain(3)

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
      ! (Re)allocate type-bound arrays
      if (allocated(this%xvert)) then
        deallocate (this%xvert)
        deallocate (this%yvert)
        deallocate (this%vne)
        deallocate (this%vv0x)
        deallocate (this%vv0y)
        deallocate (this%vv1x)
        deallocate (this%vv1y)
        deallocate (this%iprev)
        deallocate (this%xvertnext)
        deallocate (this%yvertnext)
      end if
      allocate (this%xvert(this%nverts))
      allocate (this%yvert(this%nverts))
      allocate (this%vne(this%nverts))
      allocate (this%vv0x(this%nverts))
      allocate (this%vv0y(this%nverts))
      allocate (this%vv1x(this%nverts))
      allocate (this%vv1y(this%nverts))
      allocate (this%iprev(this%nverts))
      allocate (this%xvertnext(this%nverts))
      allocate (this%yvertnext(this%nverts))
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
    call this%vertvelo()

    ! Track across subcells
    call this%track(particle, 2, tmax)

  end subroutine apply_mct

  !> @brief Loads a triangular subcell from the polygonal cell
  subroutine load_subcell(this, particle, subcell)
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
          ! assumes particle is in the cell, so no check needed for beta
          if ((alphai .ge. DZERO) .and. &
              (alphai + betai .le. DONE)) then
            isc = iv0
            exit
          end if
        end do
        if (isc .le. 0) then
          print *, "error -- initial triangle not found in cell ", ic, &
            " for particle at ", particle%x, particle%y, particle%z

          call pstop(1)
        else
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
      subcell%v1x = this%vv1x(iv0) ! the indices here actually refer to subcells, not vertices
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

  subroutine reallocate_vv(this)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    integer(I4B) :: inc

      inc = this%nverts
      call ExpandArray(this%le, inc) ! lengths of exterior (cell) edges
      call ExpandArray(this%unextx, inc) ! x components of unit normals to exterior edges
      call ExpandArray(this%unexty, inc) ! y components of unit normals to exterior edges
      call ExpandArray(this%areasub, inc) ! subcell areas
      call ExpandArray(this%li, inc) ! lengths of interior edges ("spokes")
      call ExpandArray(this%unintx, inc) ! x components of unit normals to interior edges
      call ExpandArray(this%uninty, inc) ! y components of unit normals to interior edges
      call ExpandArray(this%xmid, inc) ! x coordinates of midpoints
      call ExpandArray(this%ymid, inc) ! y coordinates of midpoints
      call ExpandArray(this%lm, inc) ! lengths of midpoint connectors
      call ExpandArray(this%umx, inc) ! x components of midpoint-connector (cw) unit vectors
      call ExpandArray(this%umy, inc) ! y components of midpoint-connector (cw) unit vectors
      call ExpandArray(this%kappax, inc) ! x components of kappa vectors
      call ExpandArray(this%kappay, inc) ! y components of kappa vectors
      call ExpandArray(this%vm0x, inc) ! x component of vm0
      call ExpandArray(this%vm0y, inc) ! y component of vm0
      call ExpandArray(this%vm1x, inc) ! x component of vm1
      call ExpandArray(this%vm1y, inc) ! y component of vm1
      call ExpandArray(this%unextxnext, inc) ! vector of "next" interior unit-normal x coordinates defined for convenience
      call ExpandArray(this%unextynext, inc) ! vector of "next" interior unit-normal y coordinates defined for convenience
      call ExpandArray(this%wk1, inc)
      call ExpandArray(this%wk2, inc)
      ! calc_thru_hsum
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
      call ExpandArray(this%be1x, inc)
      call ExpandArray(this%bi1y, inc)
      call ExpandArray(this%be1y, inc)
      call ExpandArray(this%be01x, inc)
      call ExpandArray(this%be01y, inc)
    end if
  end subroutine reallocate_vv

  !> @brief Calculate vertex velocities
  subroutine vertvelo(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    real(DP) :: term
    integer(I4B) :: i
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

      call this%reallocate_vv()     

      ! Exterior edge unit normals (outward) and lengths
      this%wk1 = this%xvertnext - this%xvert
      this%wk2 = this%yvertnext - this%yvert
      this%le = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%unextx = -this%wk2 / this%le
      this%unexty = this%wk1 / this%le

      ! Cell area
      areacell = area(this%xvert, this%yvert)

      ! Cell centroid (in general, this is NOT the average of the vertex coordinates)
      sixa = areacell * 6.d0
      this%wk1 = -(this%xvert * this%yvertnext - this%xvertnext * this%yvert)
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
        this%areasub(i) = area(this%xvals, this%yvals)
      end do

      ! Cell-edge normal velocities (outward)
      term = DONE / (cell%defn%porosity * cell%defn%retfactor * this%dz)
      do i = 1, this%nverts
        this%vne(i) = -cell%defn%faceflow(i) * term / this%le(i)
      end do

      ! Cell divergence (2D)
      divcell = sum(this%le * this%vne) / areacell

      ! Interior edge (cw) unit normals and lengths
      this%wk1 = this%xvert - this%xctr
      this%wk2 = this%yvert - this%yctr
      this%li = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%unintx = this%wk2 / this%li
      this%uninty = -this%wk1 / this%li
      ! Shifted arrays for convenience
      this%unextxnext = cshift(this%unintx, 1)
      this%unextynext = cshift(this%uninty, 1)

      ! Midpoints of interior edges
      this%xmid = 5.d-1 * (this%xvert + this%xctr)
      this%ymid = 5.d-1 * (this%yvert + this%yctr)

      ! Unit midpoint-connector (cw) vectors and lengths
      this%wk1 = cshift(this%xmid, 1) - this%xmid
      this%wk2 = cshift(this%ymid, 1) - this%ymid
      this%lm = dsqrt(this%wk1 * this%wk1 + this%wk2 * this%wk2)
      this%umx = this%wk1 / this%lm
      this%umy = this%wk2 / this%lm

      ! Kappa vectors (K tensor times unit midpoint-connector vectors) do not
      ! account for anisotropy, which is consistent with the way internal face
      ! flow calculations are done in MP7. The isotropic value of K does not
      ! matter in this case because it cancels out of the calculations, so
      ! K = 1 is assumed for simplicity.
      this%kappax = this%umx
      this%kappay = this%umy

      ! Use linearity to find vm0i[0] such that curl of the head gradient
      ! is zero
      perturb = 1.d-2
      ! Calculations at base value
      vm0i0 = 0.d0
      call this%calc_thru_hcsum(&
        vm0i0, &
        divcell, &
        areacell, &
        hcsum0)
      ! Calculations at perturbed value
      vm0ival = vm0i0 + perturb
      call this%calc_thru_hcsum(&
        vm0ival, &
        divcell, &
        areacell, &
        hcsum)
      ! Calculations at root value
      jac = (hcsum - hcsum0) / perturb
      vm0ival = vm0i0 - hcsum0 / jac
      call this%calc_thru_hcsum(&
        vm0ival, &
        divcell, &
        areacell, &
        hcsum)

      ! Project linearly to get corner (vertex) velocities. Note that velocity
      ! vv1 is at the next vertex cw from vv0, so vv0(i) and vv1(i) are the
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

  subroutine calc_thru_hcsum(this, vm0ival, divcell, areacell, hcsum)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    real(DP) :: vm0ival
    real(DP) :: divcell
    real(DP) :: hcsum
    real(DP) :: areacell
    ! local
    real(DP) :: emxx
    real(DP) :: emxy
    real(DP) :: emyx
    real(DP) :: emyy
    real(DP) :: rx
    real(DP) :: ry
    real(DP) :: emdet
    integer(I4B) :: i
    integer(I4B) :: ip

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
    this%uprod = this%unintx * this%unextx + this%uninty * this%unexty
    this%det = DONE - this%uprod * this%uprod
    this%bi0x = (this%unintx - this%unextx * this%uprod) / this%det
    this%be0x = (this%unextx - this%unintx * this%uprod) / this%det
    this%bi0y = (this%uninty - this%unexty * this%uprod) / this%det
    this%be0y = (this%unexty - this%uninty * this%uprod) / this%det
    this%uprod = this%unextxnext * this%unextx + this%unextynext * this%unexty
    this%det = DONE - this%uprod * this%uprod
    this%bi1x = (this%unextxnext - this%unextx * this%uprod) / this%det
    this%be1x = (this%unextx - this%unextxnext * this%uprod) / this%det
    this%bi1y = (this%unextynext - this%unexty * this%uprod) / this%det
    this%be1y = (this%unexty - this%unextynext * this%uprod) / this%det
    this%be01x = 5.d-1 * (this%be0x + this%be1x)
    this%be01y = 5.d-1 * (this%be0y + this%be1y)
    print *, "------", this%be01x
    this%wt = this%areasub / areacell
    emxx = DTWO - sum(this%wt * this%be01x * this%unextx)
    emxy = -sum(this%wt * this%be01x * this%unexty)
    emyx = -sum(this%wt * this%be01y * this%unextx)
    emyy = DTWO - sum(this%wt * this%be01y * this%unexty)
    rx = sum(this%wt * (this%bi0x * this%vm0i + this%bi1x * this%vm1i + this%be01x * this%vne))
    ry = sum(this%wt * (this%bi0y * this%vm0i + this%bi1y * this%vm1i + this%be01y * this%vne))
    emdet = emxx * emyy - emxy * emyx
    this%vctrx = (emyy * rx - emxy * ry) / emdet
    this%vctry = (emxx * ry - emyx * rx) / emdet

    ! Get vm0e's using "known" conditions
    this%vm0e = 5.d-1 * (this%vne + this%unextx * this%vctrx + this%unexty * this%vctry)

    ! Get vm1e's from uniformity along exterior edges
    this%vm1e = this%vm0e

    ! Transform vm0 and vm1 to (x, y) coordinates
    this%vm0x = this%bi0x * this%vm0i + this%be0x * this%vm0e
    this%vm0y = this%bi0y * this%vm0i + this%be0y * this%vm0e
    this%vm1x = this%bi1x * this%vm1i + this%be1x * this%vm0e
    this%vm1y = this%bi1y * this%vm1i + this%be1y * this%vm0e

    ! Calculate head-cycle summation (which is proportional to
    ! the curl of the head gradient)
    hcsum = sum(this%lm * (this%kappax * &
      (this%vm0x + this%vm1x) + &
      this%kappay * (this%vm0y + this%vm1y)))

  end subroutine calc_thru_hcsum

end module MethodCellTernaryModule

