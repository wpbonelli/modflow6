module MethodCellTernaryModule

  use KindModule, only: DP, I4B
  use ErrorUtilModule, only: pstop
  use MethodModule
  use MethodSubcellPoolModule
  use CellPolyModule
  use CellDefnModule
  use SubcellTriModule, only: SubcellTriType, create_subcell_tri
  use ParticleModule
  use TrackModule, only: TrackControlType
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_method_cell_ternary

  type, extends(MethodType) :: MethodCellTernaryType
    private
    double precision :: x_vert(99), y_vert(99) !< cell vertex coordinates (kluge 99, todo: allocate as needed)
    double precision :: xctr, yctr !< cell center coordinates
    double precision :: vx_vert_polygon(99), vy_vert_polygon(99) !< cell vertex velocities
    double precision :: vxctr, vyctr !< cell center velocities
    double precision :: ztop, zbot !< cell top and bottom elevations
    double precision :: dz !< cell thickness
    double precision :: vztop, vzbot !< cell top and bottom velocities
  contains
    procedure, public :: apply => apply_mct
    procedure, public :: destroy => destroy_mct
    procedure, public :: load => load_mct
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mct
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
  end subroutine destroy_mct

  !> @brief Load subcell into tracking method
  subroutine load_mct(this, particle, next_level, submethod)
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod

    select type (subcell => this%subcell)
    type is (SubcellTriType)
      call this%load_subcell(particle, subcell)
    end select
    call method_subcell_tern%init(subcell=this%subcell, trackctl=this%trackctl)
    submethod => method_subcell_tern
  end subroutine load_mct

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mct(this, particle)
    ! -- dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer :: isc, exitFace, inface, npolyverts

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
    use TdisModule, only: kper, kstp
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    double precision :: retfactor
    integer :: npolyverts, iv, ivp1, ivm1
    double precision :: x0, y0, x1, y1, x2, y2, xsum, ysum
    double precision :: vxsum, vysum, flow0, flow1, v0x, v0y
    double precision :: d01x, d01y, d02x, d02y, det, area, term

    select type (cell => this%cell)
    type is (CellPolyType)
      ! -- Update particle state, checking whether any reporting or
      ! -- termination conditions apply
      call this%update(particle, cell%defn)

      ! -- Return early if particle is done advancing
      if (.not. particle%advancing) return

      ! -- If the particle is above the top of the cell (which is presumed to
      ! -- represent a water table above the cell bottom), pass the particle
      ! -- vertically and instantaneously to the cell top elevation and save
      ! -- the particle state to output file(s).
      if (particle%z > cell%defn%top) then
        particle%z = cell%defn%top
        call this%trackctl%save(particle, kper=kper, &
                                kstp=kstp, reason=1) ! reason=1: cell transition
      end if

      npolyverts = cell%defn%npolyverts

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
      area = area * DHALF
      term = DONE / (retfactor * cell%defn%porosity * area)
      this%vzbot = cell%defn%faceflow(npolyverts + 2) * term
      this%vztop = -cell%defn%faceflow(npolyverts + 3) * term
      this%xctr = xsum / dble(npolyverts)
      this%yctr = ysum / dble(npolyverts)
      this%vxctr = vxsum / dble(npolyverts)
      this%vyctr = vysum / dble(npolyverts)

      ! -- Track across subcells
      call this%track(particle, 2, tmax) ! kluge, hardwired to level 2
    end select
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
    integer :: ic, isc, npolyverts
    integer :: iv0, iv1, ipv0, ipv1
    integer :: iv
    double precision :: x0, y0, x1, y1, x2, y2, x1rel, y1rel, x2rel, y2rel, xi, yi
    double precision :: di2, d02, d12, di1, d01, alphai, betai
    double precision :: betatol

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
      subcell%v0x = this%vx_vert_polygon(iv0)
      subcell%v0y = this%vy_vert_polygon(iv0)
      subcell%v1x = this%vx_vert_polygon(iv1)
      subcell%v1y = this%vy_vert_polygon(iv1)
      subcell%v2x = this%vxctr
      subcell%v2y = this%vyctr
      subcell%ztop = this%ztop
      subcell%zbot = this%zbot
      subcell%dz = this%dz
      subcell%vzbot = this%vzbot
      subcell%vztop = this%vztop
    end select
  end subroutine load_subcell

end module MethodCellTernaryModule
