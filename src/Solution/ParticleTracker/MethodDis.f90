module MethodDisModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE
  use MethodModule
  use MethodCellPoolModule
  ! use CellModule
  use CellDefnModule
  use CellRectModule
  use ParticleModule
  use PrtFmiModule, only: PrtFmiType
  use UtilMiscModule
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodDisType
  public :: create_methodDis

  ! -- Extend MethodType to the DIS-grid method type (MethodDisType)
  type, extends(MethodType) :: MethodDisType
    private
    type(PrtFmiType), pointer :: fmi => null() !< flow model interface
    ! type(CellDefnType), pointer :: cellDefn ! cellDefn object injected into cell method
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !< intercell flows
    type(CellRectType), pointer :: cellRect ! rectangular cell
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< pointer to aquifer porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor => null() !< pointer to retardation factor
    integer(I4B), dimension(:), pointer, contiguous :: izone => null() !< pointer to zone number
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mGD ! applies the DIS-grid method
    procedure, public :: pass => pass_mGD ! passes the particle to the next cell
    procedure, public :: loadsub => loadsub_mGD ! loads the cell method
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
  end type MethodDisType

contains

  !> @brief Create a new DIS-grid method object
  subroutine create_methodDis(methodDis)
    ! -- dummy
    type(MethodDisType), pointer :: methodDis
    ! -- local
    !
    allocate (methodDis)
    allocate (methodDis%trackingDomainType)
    !
    ! -- This method delegates tracking to a submethod
    methodDis%delegatesTracking = .TRUE.
    !
    methodDis%trackingDomainType = "Dis" ! kluge???
    !
    call create_cellRect(methodDis%cellRect)
    !
    return
  end subroutine create_methodDis

  !> @brief Destructor for a DIS-grid method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    !
    deallocate (this%trackingDomainType)
    return
  end subroutine destroy

  !> @brief Initialize a DIS-grid method object
  subroutine init(this, fmi, flowja, porosity, retfactor, izone, trackdata)
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(PrtFmiType), pointer :: fmi
    real(DP), dimension(:), pointer, contiguous :: flowja
    real(DP), dimension(:), pointer, contiguous :: porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor
    integer(I4B), dimension(:), pointer, contiguous :: izone
    type(TrackDataType), pointer :: trackdata
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

  !> @brief Load the cell geometry and method (tracking strategy)
  subroutine loadsub_mGD(this, particle, levelNext, submethod)
    use ConstantsModule, only: DZERO, DONE
    use InputOutputModule
    use GwfDisModule, only: GwfDisType ! kluge???
    ! use CellRectModule ! kluge???
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    ! -- local
    integer :: ic, icu
    ! type(CellRectType), pointer :: cellRect
    integer :: irow, jcol, klay
    double precision :: areax, areay, areaz
    double precision :: dx, dy, dz
    double precision :: factor, term
    !
    ! -- Load cell
    ic = particle%iTrackingDomain(levelNext) ! kluge note: is cell number always known coming in?
    call this%load_cellDefn(ic, this%cellRect%cellDefn)
    select type (dis => this%fmi%dis) ! kluge???
    type is (GwfDisType)
      icu = dis%get_nodeuser(ic)
      call get_ijk(icu, dis%nrow, dis%ncol, dis%nlay, irow, jcol, klay)
      dx = dis%delr(jcol)
      dy = dis%delc(irow)
      dz = this%cellRect%cellDefn%top - this%cellRect%cellDefn%bot
      this%cellRect%dx = dx
      this%cellRect%dy = dy
      this%cellRect%dz = dz
      this%cellRect%sinrot = DZERO
      this%cellRect%cosrot = DONE
      this%cellRect%xOrigin = this%cellRect%cellDefn%polyvert(1)%x ! kluge note: could avoid using polyvert here
      this%cellRect%yOrigin = this%cellRect%cellDefn%polyvert(1)%y
      this%cellRect%zOrigin = this%cellRect%cellDefn%bot
      this%cellRect%ipvOrigin = 1
      areax = dx * dz
      areay = dy * dz
      areaz = dx * dy
      ! this%cellRect%vx1 = this%cellRect%cellDefn%faceflow(1)/areax ! kluge note: assuming porosity=1. and velmult=1. for now
      ! this%cellRect%vx2 = -this%cellRect%cellDefn%faceflow(3)/areax
      ! this%cellRect%vy1 = this%cellRect%cellDefn%faceflow(4)/areay
      ! this%cellRect%vy2 = -this%cellRect%cellDefn%faceflow(2)/areay
      ! this%cellRect%vz1 = this%cellRect%cellDefn%faceflow(6)/areaz
      ! this%cellRect%vz2 = -this%cellRect%cellDefn%faceflow(7)/areaz
      ! factor = this%cellRect%cellDefn%velfactor/this%cellRect%cellDefn%porosity
      factor = DONE / this%cellRect%cellDefn%retfactor
      factor = factor / this%cellRect%cellDefn%porosity
      term = factor / areax
      this%cellRect%vx1 = this%cellRect%cellDefn%faceflow(1) * term
      this%cellRect%vx2 = -this%cellRect%cellDefn%faceflow(3) * term
      term = factor / areay
      this%cellRect%vy1 = this%cellRect%cellDefn%faceflow(4) * term
      this%cellRect%vy2 = -this%cellRect%cellDefn%faceflow(2) * term
      term = factor / areaz
      this%cellRect%vz1 = this%cellRect%cellDefn%faceflow(6) * term
      this%cellRect%vz2 = -this%cellRect%cellDefn%faceflow(7) * term
    end select
    !
    ! -- Select and initialize cell method and set cell method pointer
    call methodCellPollock%init(particle, this%cellRect, this%trackdata)
    submethod => methodCellPollock
    !
    return
    !
  end subroutine loadsub_mGD

  !> @brief Pass a particle to the next cell, if there is one
  subroutine pass_mGD(this, particle)
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! -- local
    integer :: inface, ipos, ic, inbr, idiag
    real(DP) :: z, zrel, topfrom, botfrom, top, bot, sat
    !
    inface = particle%iTrackingDomainBoundary(2)
    z = particle%z
    !
    inbr = this%cellRect%cellDefn%facenbr(inface)
    if (inbr .eq. 0) then
      ! -- Exterior face; no neighbor to map to
      ! particle%iTrackingDomain(1) = 0
      ! particle%iTrackingDomain(2) = 0      ! kluge note: set a "has_exited" attribute instead???
      ! particle%iTrackingDomain(1) = -abs(particle%iTrackingDomain(1))   ! kluge???
      ! particle%iTrackingDomain(2) = -abs(particle%iTrackingDomain(2))   ! kluge???
      particle%istatus = 2 ! kluge note: use -2 to allow check for transfer to another model???
      ! particle%iTrackingDomainBoundary(2) = -1
    else
      idiag = this%fmi%dis%con%ia(this%cellRect%cellDefn%icell)
      ipos = idiag + inbr
      ic = this%fmi%dis%con%ja(ipos) ! kluge note: use PRT model's DIS instead of fmi's???
      particle%iTrackingDomain(2) = ic
      ! call this%mapToNbrCell(this%cellRect%cellDefn,inface,z)
      if (inface .eq. 1) then
        inface = 3
      else if (inface .eq. 2) then
        inface = 4
      else if (inface .eq. 3) then
        inface = 1
      else if (inface .eq. 4) then
        inface = 2
      else if (inface .eq. 6) then
        inface = 7
      else if (inface .eq. 7) then
        inface = 6
      end if
      particle%iTrackingDomainBoundary(2) = inface
      if (inface < 5) then
        ! -- Map z between cells
        ! iatop = this%get_iatop(ic)
        ! top = this%get_top(iatop)
        topfrom = this%cellRect%cellDefn%top
        botfrom = this%cellRect%cellDefn%bot
        zrel = (z - botfrom) / (topfrom - botfrom)
        top = this%fmi%dis%top(ic) ! kluge note: use PRT model's DIS instead of fmi's???
        bot = this%fmi%dis%bot(ic)
        sat = this%fmi%gwfsat(ic)
        z = bot + zrel * sat * (top - bot)
      end if
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
  end subroutine pass_mGD

  !> @brief Apply the method to a particle
  subroutine apply_mGD(this, particle, tmax)
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
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
  end subroutine apply_mGD

  !> @brief Return the number of polygon vertices for a cell in the grid
  function get_npolyverts(this, ic) result(npolyverts)
    use GwfDisModule ! kluge???
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    integer, intent(in) :: ic
    ! -- local
    ! -- result
    integer :: npolyverts
    !
    npolyverts = 4
    !
    return
  end function get_npolyverts

  !> @brief Get the index used to locate top elevation of a cell in the grid
  function get_iatop(this, ic) result(iatop)
    use GwfDisModule ! kluge???
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    integer, intent(in) :: ic
    ! -- result
    integer :: iatop
    ! -- local
    integer :: ncpl, icu
    !
    ! select type (dis => this%fmi%dis) ! kluge??? ...
    ! type is (GwfDisType)
    !   ncpl = dis%get_ncpl()
    ! end select ! ... kluge???
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

  !> @brief Returns a top elevation based on index iatop
  function get_top(this, iatop) result(top) ! kluge note: not needed???
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
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

  !> @brief Loads cell definition from the grid
  subroutine load_cellDefn(this, ic, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
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
    class(MethodDisType), intent(inout) :: this
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

  !> @brief Load polygon vertices into cell definition from the grid
  subroutine load_cellDefn_polyverts(this, cellDefn) ! kluge note: are polyverts even needed for MethodDis???
    use InputOutputModule ! kluge
    use GwfDisModule ! kluge???
    use ConstantsModule, only: DHALF
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    integer :: ic, icu, npolyverts
    integer :: irow, jcol, klay
    double precision :: cellx, celly, dxhalf, dyhalf
    !
    ic = cellDefn%icell
    !
    if (cellDefn%npolyverts .eq. 0) cellDefn%npolyverts = this%get_npolyverts(ic) ! kluge note: just set it to 4 ???
    npolyverts = cellDefn%npolyverts
    !
    ! -- Load polygon vertices. Note that the polyverts array
    ! -- does not get reallocated if it is already allocated
    ! -- to a size greater than or equal to npolyverts+1.
    call allocate_as_needed(cellDefn%polyvert, npolyverts + 1)
    select type (dis => this%fmi%dis) ! kluge??? ...
    type is (GwfDisType)
      icu = dis%get_nodeuser(ic)
      call get_ijk(icu, dis%nrow, dis%ncol, dis%nlay, irow, jcol, klay)
      cellx = dis%cellx(jcol)
      celly = dis%celly(irow)
      dxhalf = DHALF * dis%delr(jcol)
      dyhalf = DHALF * dis%delc(irow)
      ! -- SW vertex
      cellDefn%polyvert(1)%ivert = 0 ! kluge note: there's no vertex list to reference for DIS
      cellDefn%polyvert(1)%x = cellx - dxhalf
      cellDefn%polyvert(1)%y = celly - dyhalf
      ! -- NW vertex
      cellDefn%polyvert(2)%ivert = 0
      cellDefn%polyvert(2)%x = cellx - dxhalf
      cellDefn%polyvert(2)%y = celly + dyhalf
      ! -- NE vertex
      cellDefn%polyvert(3)%ivert = 0
      cellDefn%polyvert(3)%x = cellx + dxhalf
      cellDefn%polyvert(3)%y = celly + dyhalf
      ! -- SE vertex
      cellDefn%polyvert(4)%ivert = 0
      cellDefn%polyvert(4)%x = cellx + dxhalf
      cellDefn%polyvert(4)%y = celly - dyhalf
      ! -- List wraps around for convenience
      cellDefn%polyvert(5)%ivert = cellDefn%polyvert(1)%ivert
      cellDefn%polyvert(5)%x = cellDefn%polyvert(1)%x
      cellDefn%polyvert(5)%y = cellDefn%polyvert(1)%y
    end select ! ... kluge???
    !
    return
    !
  end subroutine load_cellDefn_polyverts

  !> @brief Loads face neighbors to cell definition from the grid
  subroutine load_cellDefn_facenbr(this, cellDefn)
    use InputOutputModule ! kluge
    ! use DisvGeom             ! kluge
    use GwfDisModule ! kluge
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    integer :: ic, npolyverts
    integer :: ic1, ic2, icu1, icu2, j1, j2, iloc, ipos
    integer :: irow1, irow2, jcol1, jcol2, klay1, klay2
    integer :: iedgeface
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
    type is (GwfDisType) ! kluge
      cellDefn%facenbr = 0
      ic1 = ic
      icu1 = dis%get_nodeuser(ic1)
      call get_ijk(icu1, dis%nrow, dis%ncol, dis%nlay, irow1, jcol1, klay1)
      ncpl = dis%get_ncpl()
      call get_jk(icu1, ncpl, dis%nlay, j1, klay1)
      do iloc = 1, dis%con%ia(ic1 + 1) - dis%con%ia(ic1) - 1
        ipos = dis%con%ia(ic1) + iloc
        if (dis%con%mask(ipos) == 0) cycle ! kluge note: need mask here???
        ic2 = dis%con%ja(ipos)
        icu2 = dis%get_nodeuser(ic2)
        call get_ijk(icu2, dis%nrow, dis%ncol, dis%nlay, irow2, jcol2, klay2)
        if (klay2 == klay1) then
          ! -- Edge (polygon) face neighbor
          if (irow2 > irow1) then
            ! Neighbor to the S
            iedgeface = 4 ! kluge note: make sure this numbering is consistent with numbering in cell method
          else if (jcol2 > jcol1) then
            ! Neighbor to the E
            iedgeface = 3
          else if (irow2 < irow1) then
            ! Neighbor to the N
            iedgeface = 2
          else
            ! Neighbor to the W
            iedgeface = 1
          end if
          cellDefn%facenbr(iedgeface) = iloc
        else if (klay2 > klay1) then
          ! -- Bottom face neighbor
          cellDefn%facenbr(npolyverts + 2) = iloc
        else
          ! -- Top face neighbor
          cellDefn%facenbr(npolyverts + 3) = iloc
        end if
      end do
      ! -- List of edge (polygon) faces wraps around
      cellDefn%facenbr(npolyverts + 1) = cellDefn%facenbr(1)
    end select
    !
    return
    !
  end subroutine load_cellDefn_facenbr

  !> @brief Find the shared edge face of two cells.
  !!
  !! Find the shared edge face of cell1 shared by cell1 and cell2.
  !! isharedface will return with 0 if there is no shared edge
  !! face.  Proceed forward through ivlist1 and backward through
  !! ivlist2 as a clockwise face in cell1 must correspond to a
  !! counter clockwise face in cell2
  !!
  !<
  subroutine shared_edgeface(ivlist1, ivlist2, iedgeface) ! kluge note: based on DisvGeom shared_edge
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

  !> @brief Load flows into the cell definition
  !!
  !! Loads polygon face, top and bottom, and net distributed
  !! flows from the grid into the cell definition.
  !!
  !>
  subroutine load_cellDefn_flows(this, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
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
    ! ! -- Also set inoexitface flag.
    call allocate_as_needed(cellDefn%faceflow, npolyverts + 3)
    cellDefn%faceflow = 0d0 ! kluge note: eventually use DZERO for 0d0 throughout
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
    ! -- Add BoundaryFlows to face flows
    call this%addBoundaryFlows_cellRect(cellDefn)
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

  !> @brief Add BoundaryFlows to the cell definition faceflow array
  subroutine addBoundaryFlows_cellRect(this, cellDefn)
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(CellDefnType), intent(inout) :: cellDefn
    ! -- local
    integer :: ig, ic, npolyverts, m, n
    integer :: ioffset
    !
    ic = cellDefn%icell
    npolyverts = cellDefn%npolyverts
    !
    ! kluge note - assignment of BoundaryFlows to faceflow below assumes vertex 1
    !  is at upper left of a rectangular cell
    ! ioffset = (ic - 1)*6
    ioffset = (ic - 1) * 10
    ! cellDefn%faceflow(1) = cellDefn%faceflow(1) + this%fmi%BoundaryFlows(ioffset+4)  ! kluge note: should these be additive (seems so)???
    ! cellDefn%faceflow(2) = cellDefn%faceflow(2) + this%fmi%BoundaryFlows(ioffset+2)
    ! cellDefn%faceflow(3) = cellDefn%faceflow(3) + this%fmi%BoundaryFlows(ioffset+3)
    ! cellDefn%faceflow(4) = cellDefn%faceflow(4) + this%fmi%BoundaryFlows(ioffset+1)
    cellDefn%faceflow(1) = cellDefn%faceflow(1) + this%fmi%BoundaryFlows(ioffset + 1) ! kluge note: should these be additive (seems so)???
    cellDefn%faceflow(2) = cellDefn%faceflow(2) + this%fmi%BoundaryFlows(ioffset + 2)
    cellDefn%faceflow(3) = cellDefn%faceflow(3) + this%fmi%BoundaryFlows(ioffset + 3)
    cellDefn%faceflow(4) = cellDefn%faceflow(4) + this%fmi%BoundaryFlows(ioffset + 4)
    cellDefn%faceflow(5) = cellDefn%faceflow(1)
    ! cellDefn%faceflow(6) = cellDefn%faceflow(6) + this%fmi%BoundaryFlows(ioffset+5)
    ! cellDefn%faceflow(7) = cellDefn%faceflow(7) + this%fmi%BoundaryFlows(ioffset+6)
    cellDefn%faceflow(6) = cellDefn%faceflow(6) + this%fmi%BoundaryFlows(ioffset + 9)
    cellDefn%faceflow(7) = cellDefn%faceflow(7) + this%fmi%BoundaryFlows(ioffset + 10)
    !
    return
    !
  end subroutine addBoundaryFlows_cellRect

  !> @brief Load vertex 180-degree indicator array into cell definition
  !!
  !! Loads 180-degree vertex indicator array to cell
  !! definition and sets flags that indicate how cell
  !! can be represented
  !!
  !<
  subroutine load_cellDefn_ispv180(this, cellDefn) ! kluge note: rename???
    implicit none
    ! -- dummy
    class(MethodDisType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: cellDefn
    ! -- local
    integer :: npolyverts
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
    cellDefn%canBeCellRect = .true.
    cellDefn%canBeCellRectQuad = .false.
    !
    return
    !
  end subroutine load_cellDefn_ispv180

end module MethodDisModule
