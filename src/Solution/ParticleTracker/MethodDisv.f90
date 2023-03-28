module MethodDisvModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE
  use MethodModule
  use MethodCellPoolModule
  ! use CellModule
  use CellDefnModule
  use CellPolyModule
  use ParticleModule
  use PrtFmiModule, only: PrtFmiType
  use UtilMiscModule
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodDisvType
  public :: create_methodDisv

  ! -- Extend MethodType to the DISV-grid method type (MethodDisvType)
  type, extends(MethodType) :: MethodDisvType
    private
    type(PrtFmiType), pointer :: fmi => null() !< flow model interface
    ! type(CellDefnType), pointer :: cellDefn ! cellDefn object injected into cell method
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !< intercell flows
    type(CellPolyType), pointer :: cellPoly ! polygonal cell
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< pointer to aquifer porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor => null() !< pointer to retardation factor
    integer(I4B), dimension(:), pointer, contiguous :: izone => null() !< pointer to zone number
  contains
    procedure, public :: destroy ! destructor for the method      ! kluge note: must procedures like this be denoted as public (as and throughout)???
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mGDv ! applies the DISV-grid method
    procedure, public :: pass => pass_mGDv ! passes the particle to the next cell
    procedure, public :: loadsub => loadsub_mGDv ! loads the cell method
    procedure, public :: mapToNbrCell ! maps a location on the cell face to the shared face of a neighbor
    procedure, private :: get_npolyverts ! returns the number of polygon vertices for a cell in the grid
    procedure, private :: get_iatop ! returns index used to locate top elevation of cell in the grid
    procedure, private :: get_top ! returns top elevation based on index iatop
    procedure, public :: load_cellDefn ! loads cell definition from the grid
    procedure, public :: load_cellDefn_facenbr ! loads face neighbors to a cell object
    procedure, public :: load_cellDefn_polyverts ! loads polygon vertices to a cell object
    procedure, public :: load_cellDefn_flows ! loads flows to a cell object
    procedure, public :: load_cellDefn_ispv180 ! loads 180-degree vertex indicator to a cell object
    procedure, public :: load_cellDefn_basic ! loads basic components to a cell object from its grid
    procedure, public :: addBoundaryFlows_cellRect ! adds BoundaryFlows from the grid to the faceflow array of a rectangular cell
    procedure, public :: addBoundaryFlows_cellRectQuad ! adds BoundaryFlows from the grid to the faceflow array of a rectangular-quad cell
    procedure, public :: addBoundaryFlows_cellPoly ! adds BoundaryFlows from the grid to the faceflow array of a polygonal cell
  end type MethodDisvType

contains

  !> @brief Create a new DISV-grid method object
  subroutine create_methodDisv(methodDisv)
    ! -- dummy
    type(MethodDisvType), pointer :: methodDisv
    ! -- local
    !
    allocate (methodDisv)
    allocate (methodDisv%trackingDomainType)
    !
    ! -- This method delegates tracking to a submethod
    methodDisv%delegatesTracking = .TRUE.
    !
    methodDisv%trackingDomainType = "Disv" ! kluge???
    !
    ! call create_cellDefn(methodDisv%cellDefn)
    call create_cellPoly(methodDisv%cellPoly)
    !
    return
  end subroutine create_methodDisv

  !> @brief Destructor for a DISV-grid method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    ! -- local
    !
    deallocate (this%trackingDomainType)
    !
    return
  end subroutine destroy

  !> @brief Initialize a DISV-grid method object
  subroutine init(this, fmi, flowja, porosity, retfactor, izone, trackdata)
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(PrtFmiType), pointer :: fmi
    real(DP), dimension(:), pointer, contiguous :: flowja
    real(DP), dimension(:), pointer, contiguous :: porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor
    integer(I4B), dimension(:), pointer, contiguous :: izone
    type(TrackDataType), pointer :: trackdata
    ! -- local
    !
    this%fmi => fmi
    this%flowja => flowja
    this%porosity => porosity
    this%retfactor => retfactor
    this%izone => izone
    !
    ! -- Set pointer to model track data
    this%trackdata => trackdata
    !
    return
    !
  end subroutine init

  !> @brief Load the cell and the tracking method
  subroutine loadsub_mGDv(this, particle, levelNext, submethod)
    use CellRectModule
    use CellRectQuadModule
    ! use CellModule
    use CellUtilModule
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    ! -- local
    integer :: ic
    integer :: istatus
    type(CellRectType), pointer :: cellRect
    type(CellRectQuadType), pointer :: cellRectQuad
    ! class(CellType), pointer :: cell
    !
    ! -- Load cell
    ic = particle%iTrackingDomain(levelNext) ! kluge note: is cell number always known coming in?
    call this%load_cellDefn(ic, this%cellPoly%cellDefn)
    !
    ! -- Select and initialize cell method and set cell method pointer
    if (this%cellPoly%cellDefn%canBeCellRect) then
      call CellPolyToCellRect(this%cellPoly, cellRect, istatus)
      call methodCellPollock%init(particle, cellRect, this%trackdata)
      submethod => methodCellPollock
    else if (this%cellPoly%cellDefn%canBeCellRectQuad) then
      call CellPolyToCellRectQuad(this%cellPoly, cellRectQuad, istatus)
      call methodCellPollockQuad%init(particle, cellRectQuad, this%trackdata)
      submethod => methodCellPollockQuad
    else
      call methodCellTernary%init(particle, this%cellPoly, this%trackdata)
      submethod => methodCellTernary
    end if
    !
    return
    !
  end subroutine loadsub_mGDv

  !> @brief Pass a particle to the next cell, if there is one
  subroutine pass_mGDv(this, particle)
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! -- local
    integer :: inface, ipos, ic, inbr, idiag
    double precision :: z
    !
    inface = particle%iTrackingDomainBoundary(2)
    z = particle%z
    !
    inbr = this%cellPoly%cellDefn%facenbr(inface)
    if (inbr .eq. 0) then
      ! -- Exterior face; no neighbor to map to
      ! particle%iTrackingDomain(1) = 0
      ! particle%iTrackingDomain(2) = 0      ! kluge note: set a "has_exited" attribute instead???
      ! particle%iTrackingDomain(1) = -abs(particle%iTrackingDomain(1))   ! kluge???
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      particle%istatus = 2 ! kluge note: use -2 to allow check for transfer to another model???
      ! particle%iTrackingDomainBoundary(2) = -1
    else
      idiag = this%fmi%dis%con%ia(this%cellPoly%cellDefn%icell)
      ipos = idiag + inbr
      ic = this%fmi%dis%con%ja(ipos) ! kluge note: use PRT model's DIS instead of fmi's???
      particle%iTrackingDomain(2) = ic
      call this%mapToNbrCell(this%cellPoly%cellDefn, inface, z)
      particle%iTrackingDomainBoundary(2) = inface
      particle%iTrackingDomain(3:) = 0
      particle%iTrackingDomainBoundary(3:) = 0
      particle%z = z
      ! -- Update cell-cell flows of particle mass.
      !    Every particle is currently assigned unit mass.
      ! -- leaving old cell
      this%flowja(ipos) = this%flowja(ipos) - DONE
      ! -- entering new cell
      this%flowja(this%fmi%dis%con%isym(ipos)) &
        = this%flowja(this%fmi%dis%con%isym(ipos)) + DONE
    end if
    !
    return
    !
  end subroutine pass_mGDv

  !> @brief Map location on cell face to shared cell face of neighbor
  subroutine mapToNbrCell(this, cellDefnin, inface, z)
    ! dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefnin ! kluge???
    integer, intent(inout) :: inface
    double precision, intent(inout) :: z
    ! local
    integer :: icin, npolyvertsin
    integer :: ic, npolyverts, inbr, inbrnbr, j, m
    real(DP) :: zrel, topfrom, botfrom, top, bot, sat
    type(CellDefnType), pointer :: cellDefn ! kluge???
    !
    ! -- Map to shared cell face of neighbor
    inbr = cellDefnin%facenbr(inface)
    if (inbr .eq. 0) then ! kluge note: redundant check
      ! -- Exterior face; no neighbor to map to
      inface = -1 ! kluge???
    else
      ! -- Load definition for neighbor cell (neighbor with shared face)
      icin = cellDefnin%icell
      j = this%fmi%dis%con%ia(icin)
      ic = this%fmi%dis%con%ja(j + inbr)
      call create_cellDefn(cellDefn)
      call this%load_cellDefn(ic, cellDefn) ! kluge  ! kluge note: really only need to load facenbr and npolyverts for this purpose
      npolyvertsin = cellDefnin%npolyverts
      npolyverts = cellDefn%npolyverts
      if (inface .eq. npolyvertsin + 2) then
        ! -- Exits through bot, enters through top
        inface = npolyverts + 3
      else if (inface .eq. npolyvertsin + 3) then
        ! -- Exits through top, enters through bot
        inface = npolyverts + 2
      else
        ! -- Exits and enters through shared polygon face
        j = this%fmi%dis%con%ia(ic)
        do m = 1, npolyverts + 3 ! kluge note: use shared_edge in DisvGeom to find shared polygon face???
          inbrnbr = cellDefn%facenbr(m)
          if (this%fmi%dis%con%ja(j + inbrnbr) .eq. icin) then
            inface = m
            exit
          end if
        end do
        ! -- Map z between cells
        ! iatop = this%get_iatop(ic)
        ! top = this%get_top(iatop)
        topfrom = cellDefnin%top
        botfrom = cellDefnin%bot
        zrel = (z - botfrom) / (topfrom - botfrom)
        top = this%fmi%dis%top(ic) ! kluge note: use PRT model's DIS instead of fmi's???
        bot = this%fmi%dis%bot(ic)
        sat = this%fmi%gwfsat(ic)
        z = bot + zrel * sat * (top - bot)
      end if
      call cellDefn%destroy()
      deallocate (cellDefn)
    end if
    !
    return
    !
  end subroutine mapToNbrCell

  !> @brief Apply the DISV-grid method
  subroutine apply_mGDv(this, particle, tmax)
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! doubleprecision :: initialTime,maximumTime,t   ! kluge not in arg list yet
    ! -- local
    !
    ! -- Track across cells
    call this%subtrack(particle, 1, tmax) ! kluge, hardwired to level 1
    !
    return
    !
  end subroutine apply_mGDv

  !> @brief Return the number of polygon vertices for a cell in the grid
  function get_npolyverts(this, ic) result(npolyverts)
    use GwfDisvModule ! kluge???
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    integer, intent(in) :: ic
    ! -- local
    integer :: lay, icu, icu2d
    integer :: ncpl ! kluge???
    ! -- result
    integer :: npolyverts
    !
    select type (dis => this%fmi%dis) ! kluge??? ...
    type is (GwfDisvType)
      ncpl = dis%get_ncpl()
      icu = this%fmi%dis%get_nodeuser(ic)
      icu2d = icu - ((icu - 1) / ncpl) * ncpl ! kluge note: use MOD or MODULO???
      npolyverts = dis%iavert(icu2d + 1) - dis%iavert(icu2d) - 1
      if (npolyverts .le. 0) npolyverts = npolyverts + size(dis%javert) ! kluge???
    end select ! ... kluge???
    !
    return
  end function get_npolyverts

  !> @brief Get index used to locate top elevation of a cell in the grid
  function get_iatop(this, ic) result(iatop)
    use GwfDisvModule ! kluge???
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    integer, intent(in) :: ic
    ! -- result
    integer :: iatop
    ! -- local
    integer :: ncpl, icu
    !
    ! select type (dis => this%fmi%dis)                           ! kluge??? ...
    ! type is (GwfDisvType)
    !   ncpl = dis%get_ncpl()
    ! end select                                              ! ... kluge???
    ! if (ic.le.ncpl) then
    !   iatop = -ic
    ! else
    !   iatop = ic - ncpl
    ! end if
    ncpl = this%fmi%dis%get_ncpl()
    icu = this%fmi%dis%get_nodeuser(ic)
    if (icu .le. ncpl) then
      iatop = -1 ! kluge note: just returns -1 if in top layer and 1 otherwise
    else
      iatop = 1
    end if
    !
    return
  end function get_iatop

  !> @brief Get top elevation based on index iatop
  !! kluge note: not needed???
  function get_top(this, iatop) result(top) 
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    integer, intent(in) :: iatop
    ! -- result
    double precision :: top
    !
    if (iatop .lt. 0) then
      top = this%fmi%dis%top(-iatop)
    else
      top = this%fmi%dis%bot(iatop)
    end if
    !
    return
  end function get_top

  !> @brief Load cell definition from the grid
  subroutine load_cellDefn(this, ic, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    integer, intent(in) :: ic
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    !
    ! -- Load basic cell definition components
    call this%load_cellDefn_basic(ic, cellDefn)
    !
    ! -- Load polygon vertices
    call this%load_cellDefn_polyverts(cellDefn)
    !
    ! -- Load face neighbors
    call this%load_cellDefn_facenbr(cellDefn)
    !
    ! -- Load 180-degree indicator
    call this%load_cellDefn_ispv180(cellDefn)
    !
    ! -- Load flows
    ! -- (assumes face neighbors already loaded)
    call this%load_cellDefn_flows(cellDefn)
    !
    return
    !
  end subroutine load_cellDefn

  !> @brief Loads basic cell definition components from the grid
  subroutine load_cellDefn_basic(this, ic, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    integer, intent(in) :: ic
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    ! integer :: iatop,nnbrs
    integer :: iatop
    real(DP) :: top, bot, sat
    !
    cellDefn%icell = ic
    !
    ! ! -- Load nnbrs
    ! nnbrs = this%fmi%dis%con%ia(ic+1) - this%fmi%dis%con%ia(ic) - 1
    ! cellDefn%nnbrs = nnbrs          ! kluge note: is this actually needed anywhere???
    ! !
    ! -- Load npolyverts
    cellDefn%npolyverts = this%get_npolyverts(ic)
    !
    ! -- Load iatop, top, and bot
    iatop = this%get_iatop(ic)
    cellDefn%iatop = iatop
    ! cellDefn%top = this%get_top(iatop)
    top = this%fmi%dis%top(ic)
    bot = this%fmi%dis%bot(ic)
    sat = this%fmi%gwfsat(ic)
    top = bot + sat * (top - bot)
    cellDefn%top = top
    cellDefn%bot = bot
    !
    ! -- Load porosity, retfactor, and izone
    cellDefn%porosity = this%porosity(ic)
    cellDefn%retfactor = this%retfactor(ic)
    cellDefn%izone = this%izone(ic)
    !
    return
    !
  end subroutine load_cellDefn_basic

  !> @brief Loads polygon vertices to cell definition from the grid
  subroutine load_cellDefn_polyverts(this, cellDefn)
    use GwfDisvModule ! kluge???
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    integer :: ig, ic, npolyverts, iavert, iavertm1, m, j, icu, icu2d
    integer :: ncpl ! kluge???
    !
    ic = cellDefn%icell
    !
    if (cellDefn%npolyverts .eq. 0) cellDefn%npolyverts = this%get_npolyverts(ic)
    npolyverts = cellDefn%npolyverts
    !
    ! -- Load polygon vertices. Note that the polyverts array
    ! -- does not get reallocated if it is already allocated
    ! -- to a size greater than or equal to npolyverts+1.
    call allocate_as_needed(cellDefn%polyvert, npolyverts + 1)
    select type (dis => this%fmi%dis) ! kluge??? ...
    type is (GwfDisvType)
      ncpl = dis%get_ncpl()
      icu = this%fmi%dis%get_nodeuser(ic)
      icu2d = icu - ((icu - 1) / ncpl) * ncpl ! kluge note: use MOD or MODULO???
      iavert = dis%iavert(icu2d)
      iavertm1 = iavert - 1
      do m = 1, npolyverts
        j = dis%javert(iavertm1 + m)
        cellDefn%polyvert(m)%ivert = j
        cellDefn%polyvert(m)%x = dis%vertices(1, j)
        cellDefn%polyvert(m)%y = dis%vertices(2, j)
      end do
      ! -- List wraps around for convenience
      select type (vert => cellDefn%polyvert(npolyverts + 1))
      type is (VertexType)
        vert = cellDefn%polyvert(1)
      end select
    end select ! ... kluge???
    !
    return
    !
  end subroutine load_cellDefn_polyverts

  !> @brief Loads face neighbors to cell definition from the grid
  subroutine load_cellDefn_facenbr(this, cellDefn)
    use InputOutputModule ! kluge
    ! use DisvGeom             ! kluge
    use GwfDisvModule ! kluge
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    integer :: ic, npolyverts
    integer :: ic1, ic2, icu1, icu2, j1, j2, k1, k2, iloc, ipos
    integer :: istart1, istart2, istop1, istop2, iedgeface
    integer :: ncpl
    !
    ic = cellDefn%icell
    !
    if (cellDefn%npolyverts .eq. 0) cellDefn%npolyverts = this%get_npolyverts(ic)
    npolyverts = cellDefn%npolyverts
    !
    ! -- Load face neighbors. Note that the facenbr array
    ! -- does not get reallocated if it is already allocated
    ! -- to a size greater than or equal to npolyverts+3.
    call allocate_as_needed(cellDefn%facenbr, npolyverts + 3)
    select type (dis => this%fmi%dis) ! kluge type guard
    type is (GwfDisvType) ! kluge
      cellDefn%facenbr = 0
      ic1 = ic
      icu1 = dis%get_nodeuser(ic1)
      ncpl = dis%get_ncpl()
      call get_jk(icu1, ncpl, dis%nlay, j1, k1)
      istart1 = dis%iavert(j1)
      istop1 = dis%iavert(j1 + 1) - 1
      do iloc = 1, dis%con%ia(ic1 + 1) - dis%con%ia(ic1) - 1
        ipos = dis%con%ia(ic1) + iloc
        if (dis%con%mask(ipos) == 0) cycle ! kluge note: need mask here???
        ic2 = dis%con%ja(ipos)
        icu2 = dis%get_nodeuser(ic2)
        call get_jk(icu2, ncpl, dis%nlay, j2, k2)
        istart2 = dis%iavert(j2)
        istop2 = dis%iavert(j2 + 1) - 1
        call shared_edgeface(dis%javert(istart1:istop1), &
                             dis%javert(istart2:istop2), &
                             iedgeface)
        if (iedgeface /= 0) then
          ! -- Edge (polygon) face neighbor
          cellDefn%facenbr(iedgeface) = iloc
        else
          if (k2 > k1) then
            ! -- Bottom face neighbor
            cellDefn%facenbr(npolyverts + 2) = iloc
          else if (k2 < k1) then
            ! -- Top face neighbor
            cellDefn%facenbr(npolyverts + 3) = iloc
          else
            print *, "Error -- k2 should be <> k1, since no shared edge face"
          end if
        end if
      end do
      ! -- List of edge (polygon) faces wraps around
      cellDefn%facenbr(npolyverts + 1) = cellDefn%facenbr(1)
    end select
    !
    return
    !
  end subroutine load_cellDefn_facenbr

  !> @brief Find the edge face shared by two cells
  !!
  !! Find the shared edge face of cell1 shared by cell1 and cell2.
  !! isharedface will return with 0 if there is no shared edge
  !! face.  Proceed forward through ivlist1 and backward through
  !! ivlist2 as a clockwise face in cell1 must correspond to a
  !! counter clockwise face in cell2.
  !!
  !! kluge note: based on DisvGeom shared_edge
  !<
  subroutine shared_edgeface(ivlist1, ivlist2, iedgeface) 
    integer(I4B), dimension(:) :: ivlist1
    integer(I4B), dimension(:) :: ivlist2
    integer(I4B), intent(out) :: iedgeface
    integer(I4B) :: nv1
    integer(I4B) :: nv2
    integer(I4B) :: il1
    integer(I4B) :: il2
    logical :: found
    !
    found = .false.
    nv1 = size(ivlist1)
    nv2 = size(ivlist2)
    iedgeface = 0
    outerloop: do il1 = 1, nv1 - 1
      do il2 = nv2, 2, -1
        if (ivlist1(il1) == ivlist2(il2) .and. &
            ivlist1(il1 + 1) == ivlist2(il2 - 1)) then
          found = .true.
          iedgeface = il1
          exit outerloop
        end if
      end do
      if (found) exit
    end do outerloop
  end subroutine shared_edgeface

  !> @brief Load flows into the cell from the grid
  !!
  !! Load polygon face, top and bottom, and net distributed
  !! flows to cell definition from the grid.
  !!
  !<
  subroutine load_cellDefn_flows(this, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), intent(inout) :: cellDefn
    ! -- local
    integer :: ig, ic, npolyverts, m, n

    integer :: ioffset, nbf, m1, m2, mdiff ! kluge
    double precision :: qbf ! kluge
    !
    ic = cellDefn%icell
    !
    if (cellDefn%npolyverts .eq. 0) cellDefn%npolyverts = this%get_npolyverts(ic)
    npolyverts = cellDefn%npolyverts
    !
    ! -- Load face flows. Note that the faceflow array
    ! -- does not get reallocated if it is already allocated
    ! -- to a size greater than or equal to npolyverts+3.
    ! ! -- Also set noexitface flag.
    call allocate_as_needed(cellDefn%faceflow, npolyverts + 3)
    cellDefn%faceflow = 0d0
    ! cellDefn%inoexitface = 1
    ! -- As with polygon nbrs, polygon face flows wrap around for
    ! -- convenience at position npolyverts+1, and bot and top flows
    ! -- are tacked on the end of the list
    do m = 1, npolyverts + 3
      n = cellDefn%facenbr(m)
      if (n > 0) &
        cellDefn%faceflow(m) = this%fmi%gwfflowja(this%fmi%dis%con%ia(ic) + n)
        ! if (cellDefn%faceflow(m) < 0d0) cellDefn%inoexitface = 0
    end do
    call this%addBoundaryFlows_cellPoly(cellDefn)
    ! -- Set inoexitface flag
    cellDefn%inoexitface = 1
    do m = 1, npolyverts + 3 ! kluge note: can be streamlined with above code
      if (cellDefn%faceflow(m) < 0d0) cellDefn%inoexitface = 0
    end do
    !
    ! -- Add up net distributed flow
    cellDefn%distflow = this%fmi%SourceFlows(ic) + this%fmi%SinkFlows(ic) + &
                        this%fmi%StorageFlows(ic)
    !
    ! -- Set weak sink flag
    if (this%fmi%SinkFlows(ic) .ne. 0d0) then
      cellDefn%iweaksink = 1
    else
      cellDefn%iweaksink = 0
    end if
    !
    return
    !
  end subroutine load_cellDefn_flows

  !> @brief Load boundary flows from the grid into a rectangular cell
  subroutine addBoundaryFlows_cellRect(this, cellDefn) 
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), intent(inout) :: cellDefn
    ! -- local
    integer :: ig, ic, npolyverts, m, n
    integer :: ioffset
    !
    ic = cellDefn%icell
    npolyverts = cellDefn%npolyverts
    !
    ! kluge note - assignment of BoundaryFlows to faceflow below assumes vertex 1
    !  is at upper left of a rectangular cell and that BoundaryFlows still use old iface ordering
    ! ioffset = (ic - 1)*6
    ioffset = (ic - 1) * 10
    cellDefn%faceflow(1) = cellDefn%faceflow(1) + this%fmi%BoundaryFlows(ioffset + 4) ! kluge note: should these be additive (seems so)???
    cellDefn%faceflow(2) = cellDefn%faceflow(2) + this%fmi%BoundaryFlows(ioffset + 2)
    cellDefn%faceflow(3) = cellDefn%faceflow(3) + this%fmi%BoundaryFlows(ioffset + 3)
    cellDefn%faceflow(4) = cellDefn%faceflow(4) + this%fmi%BoundaryFlows(ioffset + 1)
    cellDefn%faceflow(5) = cellDefn%faceflow(1)
    ! cellDefn%faceflow(6) = cellDefn%faceflow(6) + this%fmi%BoundaryFlows(ioffset+5)
    ! cellDefn%faceflow(7) = cellDefn%faceflow(7) + this%fmi%BoundaryFlows(ioffset+6)
    cellDefn%faceflow(6) = cellDefn%faceflow(6) + this%fmi%BoundaryFlows(ioffset + 9)
    cellDefn%faceflow(7) = cellDefn%faceflow(7) + this%fmi%BoundaryFlows(ioffset + 10)
    !
    return
    !
  end subroutine addBoundaryFlows_cellRect

  !> @brief Load boundary flows from the grid into rectangular quadcell
  subroutine addBoundaryFlows_cellRectQuad(this, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), intent(inout) :: cellDefn
    ! -- local
    integer :: ig, ic, npolyverts, m, n, nn
    integer :: ioffset, nbf, m1, m2, mdiff
    double precision :: qbf
    integer :: irectvert(5) ! kluge
    !
    ic = cellDefn%icell
    npolyverts = cellDefn%npolyverts
    !
    ! kluge note - assignment of BoundaryFlows to faceflow below assumes vertex 1
    !  is at upper left of a rectangular cell and that BoundaryFlows still use old iface ordering
    ! ioffset = (ic - 1)*6
    ioffset = (ic - 1) * 10
    ! -- Polygon faces in positions 1 through npolyverts
    do n = 1, 4
      if (n .eq. 2) then
        nbf = 4
      else if (n .eq. 4) then
        nbf = 1
      else
        nbf = n
      end if
      qbf = this%fmi%BoundaryFlows(ioffset + nbf)
      ! call cellDefn%set_irectvert()    ! kluge
      nn = 0 ! kluge ...
      do m = 1, npolyverts
        if (.not. cellDefn%ispv180(m)) then
          nn = nn + 1
          irectvert(nn) = m
        end if
      end do
      irectvert(5) = irectvert(1) ! ... kluge
      m1 = irectvert(n)
      m2 = irectvert(n + 1)
      if (m2 .lt. m1) m2 = m2 + npolyverts
      mdiff = m2 - m1
      if (mdiff .eq. 1) then
        ! -- Assign BoundaryFlow to corresponding polygon face
        cellDefn%faceflow(m1) = cellDefn%faceflow(m1) + qbf
      else
        ! -- Split BoundaryFlow between two faces on quad-refined edge
        qbf = 5d-1 * qbf
        cellDefn%faceflow(m1) = cellDefn%faceflow(m1) + qbf
        cellDefn%faceflow(m1 + 1) = cellDefn%faceflow(m1 + 1) + qbf
      end if
    end do
    ! -- Wrap around to 1 in position npolyverts+1
    m = npolyverts + 1
    cellDefn%faceflow(m) = cellDefn%faceflow(1)
    ! -- Bottom in position npolyverts+2
    m = m + 1
    ! cellDefn%faceflow(m) = cellDefn%faceflow(m) + this%fmi%BoundaryFlows(ioffset+5)
    cellDefn%faceflow(m) = cellDefn%faceflow(m) + this%fmi%BoundaryFlows(ioffset + 9)
    ! -- Top in position npolyverts+3
    m = m + 1
    ! cellDefn%faceflow(m) = cellDefn%faceflow(m) + this%fmi%BoundaryFlows(ioffset+6)
    cellDefn%faceflow(m) = cellDefn%faceflow(m) + this%fmi%BoundaryFlows(ioffset + 10)
    !
    return
    !
  end subroutine addBoundaryFlows_cellRectQuad

  !> @brief Load boundary flows from the grid into a polygonal cell
  subroutine addBoundaryFlows_cellPoly(this, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), intent(inout) :: cellDefn
    ! -- local
    integer :: ic, npolyverts, ioffset, iv
    !
    ic = cellDefn%icell
    npolyverts = cellDefn%npolyverts
    !
    ! ioffset = (ic - 1)*6
    ioffset = (ic - 1) * 10 ! kluge note: hardwired for max 8 polygon faces plus top and bottom for now
    do iv = 1, npolyverts
      cellDefn%faceflow(iv) = cellDefn%faceflow(iv) + this%fmi%BoundaryFlows(ioffset+iv)  ! kluge note: should these be additive (seems so)???
    end do
    cellDefn%faceflow(npolyverts + 1) = cellDefn%faceflow(1)
    ! cellDefn%faceflow(npolyverts+2) = cellDefn%faceflow(npolyverts+2) + this%fmi%BoundaryFlows(ioffset+npolyverts+1)
    ! cellDefn%faceflow(npolyverts+3) = cellDefn%faceflow(npolyverts+3) + this%fmi%BoundaryFlows(ioffset+npolyverts+2)
    cellDefn%faceflow(npolyverts+2) = cellDefn%faceflow(npolyverts+2) + this%fmi%BoundaryFlows(ioffset+9)
    cellDefn%faceflow(npolyverts+3) = cellDefn%faceflow(npolyverts+3) + this%fmi%BoundaryFlows(ioffset+10)
    !
    return
    !
  end subroutine addBoundaryFlows_cellPoly

  !> @brief Load 180-degree vertex indicator array into cell.
  !!
  !! Loads 180-degree vertex indicator array to cell
  !! definition and sets flags that indicate how cell
  !! can be represented
  !!
  !<
  subroutine load_cellDefn_ispv180(this, cellDefn) ! kluge note: rename???
    implicit none
    ! -- dummy
    class(MethodDisvType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    ! integer :: ig,ic,npolyverts,iapnbr,iapv180
    integer :: npolyverts, m, j, m0, m1, m2
    integer :: ic
    integer :: num90, num180, numacute
    double precision :: x0, y0, x1, y1, x2, y2
    double precision :: epsang,epslen,s0x,s0y,s0mag,s2x,s2y,s2mag,sinang,dotprod
    logical last180
    !
    ic = cellDefn%icell
    !
    ! -- Set up polyverts if not done previously
    if (cellDefn%npolyverts .eq. 0) call this%load_cellDefn_polyverts(cellDefn)
    npolyverts = cellDefn%npolyverts
    !
    ! -- Load 180-degree indicator. Note that the ispv180 array
    ! -- does not get reallocated if it is already allocated
    ! -- to a size greater than or equal to npolyverts+1.
    ! -- Also, set flags that indicate how cell can be represented.
    call allocate_as_needed(cellDefn%ispv180, npolyverts + 1)
    cellDefn%ispv180(1:npolyverts + 1) = .false.
    cellDefn%canBeCellRect = .false.
    cellDefn%canBeCellRectQuad = .false.
    epsang = 1d-3 ! kluge hardwire, and using one value for all angles
    epslen = 1d-3 ! kluge hardwire
    num90 = 0
    num180 = 0
    numacute = 0
    last180 = .false.
    do m = 1, npolyverts ! kluge note: assumes non-self-intersecting polygon; no checks for self-intersection (e.g., star)
      m1 = m
      if (m1 .eq. 1) then
        m0 = npolyverts
        m2 = 2
      else if (m1 .eq. npolyverts) then
        m0 = npolyverts - 1
        m2 = 1
      else
        m0 = m1 - 1
        m2 = m1 + 1
      end if
      x0 = cellDefn%polyvert(m0)%x
      y0 = cellDefn%polyvert(m0)%y
      x1 = cellDefn%polyvert(m1)%x
      y1 = cellDefn%polyvert(m1)%y
      x2 = cellDefn%polyvert(m2)%x
      y2 = cellDefn%polyvert(m2)%y
      s0x = x0 - x1
      s0y = y0 - y1
      s0mag = dsqrt(s0x * s0x + s0y * s0y)
      s2x = x2 - x1
      s2y = y2 - y1
      s2mag = dsqrt(s2x * s2x + s2y * s2y)
      sinang = (s0x * s2y - s0y * s2x) / (s0mag * s2mag)
      if (dabs(sinang) .lt. epsang) then ! kluge note: is it better to check in terms of angle rather than sin{angle}???
        dotprod = s0x * s2x + s0y * s2y
        if (dotprod .gt. 0d0) then
          print *, "Cell ", ic, " has a zero angle" ! kluge
          print *, "      (tolerance epsang = ", epsang, ")"
          !!pause
          stop
        else
          if (last180) then
            print *, "Cell ", ic, " has consecutive 180-deg angles - not supported" ! kluge
            print *, "      (tolerance epsang = ", epsang, ")"
            !!pause
            stop
          else if (dabs((s2mag - s0mag) / max(s2mag, s0mag)) .gt. epslen) then
            print *, "Cell ", ic, " has a non-bisecting 180-deg vertex - not supported" ! kluge
            print *, "      (tolerance epslen = ", epslen, ")"
            !!pause
            stop
          end if
          num180 = num180 + 1 ! kluge note: want to evaluate 180-deg vertex using one criterion implemented in one place (procedure) to avoid potential disparities between multiple checks
          last180 = .true.
          cellDefn%ispv180(m) = .true.
        end if
      else if (sinang .gt. 0d0) then
        numacute = numacute + 1
        if (dabs(1d0 - sinang) .lt. epsang) num90 = num90 + 1
        last180 = .false.
      else
        print *, "Cell ", ic, " has an obtuse angle and so is nonconvex" ! kluge
        print *, "      (tolerance epsang = ", epsang, ")"
        !!pause
        stop
      end if
    end do
    if ((num90 .ne. 4) .and. (num180 .ne. 0)) then
      print *, "Cell ", ic, " is a non-rectangle with a 180-deg angle - not supported" ! kluge
      print *, "      (tolerance epsang = ", epsang, ")"
      ! pause
      stop
    end if
    ! -- List of 180-degree indicators wraps around for convenience
    cellDefn%ispv180(npolyverts + 1) = cellDefn%ispv180(1)
    !
    if (num90 .eq. 4) then
      if (num180 .eq. 0) then
        ! -- Can become CellRect
        cellDefn%canBeCellRect = .true.
      else
        ! -- Can become CellRectQuad
        cellDefn%canBeCellRectQuad = .true.
      end if
    end if
    !
    return
    !
  end subroutine load_cellDefn_ispv180

end module MethodDisvModule
