module CellRectQuadModule

  use CellDefnModule, only: CellDefnType
  implicit none

  private
  public :: CellRectQuadType
  public :: create_cellRectQuad

  type :: CellRectQuadType
    private
    character(len=40), pointer, public :: type ! character string that names the tracking domain type
    type(CellDefnType), pointer, public :: cellDefn => null() ! pointer to basic cell definition
    double precision, public :: dx ! dimension of cell in local x direction
    double precision, public :: dy ! dimension of cell in local y direction
    double precision, public :: dz ! dimension of cell in z direction
    double precision, public :: sinrot ! sine of rotation angle for local (x, y)
    double precision, public :: cosrot ! cosine of rotation angle for local (x, y)
    double precision, public :: xOrigin ! model x origin for local (x, y)
    double precision, public :: yOrigin ! model y origin for local (x, y)
    double precision, public :: zOrigin ! model z origin for local z
    integer, public :: irvOrigin ! origin rectangle vertex
    double precision, public :: qextl1(4), qextl2(4), qintl(5) ! external and internal subcell flows for the cell
    integer, allocatable, public :: irectvert(:) ! list of indices of the rectangle vertices
    integer, allocatable, public :: ipv4irv(:, :) ! list of the polygon vertex indices that correspond to the rectangle vertex indices
    double precision, allocatable, public :: rectflow(:, :) ! flow(s) for each rectangle face
  contains
    procedure, public :: destroy => destroy_cellRectQuad ! destructor for the cell
    procedure, public :: init => init_cellRectQuad ! initializes the rectangular-quad cell
    procedure :: load_irectvert ! loads list of indices of the rectangle vertices
    procedure :: get_irectvertSW ! gets index of southwest rectangle vertex
    procedure :: get_rectDimensionsRotation ! gets rectangular dimensions and rotation
    procedure, public :: get_rectflow ! returns a rectangle face flow
    procedure, public :: isRefinedFace ! returns whether a rectangle face is refined
  end type CellRectQuadType

contains

  !> @brief Create a new rectangular-quad cell
  subroutine create_cellRectQuad(cellRectQuad)
    ! -- dummy
    type(CellRectQuadType), pointer :: cellRectQuad
    !
    allocate (cellRectQuad)
    allocate (cellRectQuad%cellDefn)
    allocate (cellRectQuad%irectvert(5))
    allocate (cellRectQuad%ipv4irv(2, 4))
    allocate (cellRectQuad%rectflow(2, 4))
    allocate (cellRectQuad%type)
    cellRectQuad%type = 'CellRectQuad'
    !
    return
  end subroutine create_cellRectQuad

  !> @brief Destructor for a rectangular-quad cell
  subroutine destroy_cellRectQuad(this)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    !
    deallocate (this%cellDefn)
    deallocate (this%irectvert)
    deallocate (this%type)
    !
    return
  end subroutine destroy_cellRectQuad

  !> @brief Initialize a rectangular-quad cell
  subroutine init_cellRectQuad(this, cellDefn)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    type(CellDefnType), pointer :: cellDefn
    !
    ! -- Set pointer to cell definition
    this%cellDefn => cellDefn
    !
    ! -- Load the "rectangle vertices" for the cell
    call this%load_irectvert()
    !
    return
  end subroutine init_cellRectQuad

  !> @brief Load local polygon vertex indices
  !!
  !! Loads local polygon vertex indices of the four rectangle
  !! vertices of a rectangular-quad cell. Todo: rename?
  !!
  !<
  subroutine load_irectvert(this)
    implicit none
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    ! -- local
    integer :: npolyverts, n, m
    !
    npolyverts = this%cellDefn%get_npolyverts()
    !
    n = 0
    do m = 1, npolyverts
      if (.not. this%cellDefn%get_ispv180(m)) then
        n = n + 1
        this%irectvert(n) = m
        this%ipv4irv(1, n) = m
        this%rectflow(1, n) = this%cellDefn%get_faceflow(m)
        this%ipv4irv(2, n) = 0
        this%rectflow(2, n) = 0d0
      else
        if (n .ne. 0) then
          this%ipv4irv(2, n) = m
          this%rectflow(2, n) = this%cellDefn%get_faceflow(m)
        end if
      end if
    end do
    ! Wrap around for convenience
    this%irectvert(5) = this%irectvert(1)
    !
    return
  end subroutine load_irectvert

  !> @brief Get index of SW rectangle vertex
  !!
  !! Return the index (1, 2, 3, or 4) of the southwest
  !! rectangle vertex of a rectangular-quad cell
  !!
  !<
  function get_irectvertSW(this) result(irv1)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    integer :: irv1
    ! -- local
    integer :: irv, irv2, irv4, ipv1, ipv2, ipv4
    integer, dimension(4) :: irvnxt = (/2, 3, 4, 1/) ! kluge???
    double precision :: x1, y1, x2, y2, x4, y4
    !
    ! -- Find the "southwest" rectangle vertex by finding the vertex formed
    ! -- either by (1) a rectangle edge over which x decreases (going
    ! -- clockwise) followed by an edge over which x does not increase, or by
    ! -- (2) a rectangle edge over which y does not decrease (again going
    ! -- clockwise) followed by a rectangle edge over which y increases. In
    ! -- the end, ipv1 is the index (1, 2, 3, or 4) of the southwest
    ! -- rectangle vertex.
    do irv = 1, 4
      irv4 = irv
      irv1 = irvnxt(irv4)
      ipv4 = this%irectvert(irv4)
      ipv1 = this%irectvert(irv1)
      x4 = this%cellDefn%polyvert(ipv4)%x
      y4 = this%cellDefn%polyvert(ipv4)%y
      x1 = this%cellDefn%polyvert(ipv1)%x
      y1 = this%cellDefn%polyvert(ipv1)%y
      if (x1 .lt. x4) then
        irv2 = irvnxt(irv1)
        ipv2 = this%irectvert(irv2)
        x2 = this%cellDefn%polyvert(ipv2)%x
        if (x2 .le. x1) return
      else if (y1 .ge. y4) then
        irv2 = irvnxt(irv1)
        ipv2 = this%irectvert(irv2)
        y2 = this%cellDefn%polyvert(ipv2)%y
        if (y2 .gt. y1) return
      end if
    end do
    !
    return
    !
  end function get_irectvertSW

  !> @brief Get rectangular cell dimensions and rotation
  !!
  !! Compute rectangular dimensions and rotation of
  !! the cell using the specified rectangle vertex
  !! as the origin
  !!
  !<
  subroutine get_rectDimensionsRotation(this, irv1, xOrigin, yOrigin, zOrigin, &
                                        dx, dy, dz, sinrot, cosrot)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    integer :: irv1
    double precision :: xOrigin, yOrigin, zOrigin, dx, dy, dz, sinrot, cosrot
    ! -- local
    integer :: irv2, irv4, ipv1, ipv2, ipv4
    integer, dimension(4) :: irvnxt = (/2, 3, 4, 1/) ! kluge???
    double precision :: x1, y1, x2, y2, x4, y4, dx2, dy2, dx4, dy4
    !
    ! -- Get rectangle vertex neighbors irv2 and irv4
    irv2 = irvnxt(irv1)
    irv4 = irvnxt(irvnxt(irv2)) ! kluge
    !
    ! -- Get model coordinates at irv1, irv2, and irv4
    ipv1 = this%irectvert(irv1)
    x1 = this%cellDefn%polyvert(ipv1)%x
    y1 = this%cellDefn%polyvert(ipv1)%y
    ipv2 = this%irectvert(irv2)
    x2 = this%cellDefn%polyvert(ipv2)%x
    y2 = this%cellDefn%polyvert(ipv2)%y
    ipv4 = this%irectvert(irv4)
    x4 = this%cellDefn%polyvert(ipv4)%x
    y4 = this%cellDefn%polyvert(ipv4)%y
    !
    ! -- Compute rectangle dimensions
    xOrigin = x1
    yOrigin = y1
    zOrigin = this%cellDefn%bot
    dx2 = x2 - xOrigin
    dy2 = y2 - yOrigin
    dx4 = x4 - xOrigin
    dy4 = y4 - yOrigin
    dx = dsqrt(dx4 * dx4 + dy4 * dy4)
    dy = dsqrt(dx2 * dx2 + dy2 * dy2)
    dz = this%cellDefn%top - zOrigin ! kluge note: need to account for partial saturation
    !
    ! -- Compute sine and cosine of rotation angle (angle between "southern"
    ! -- rectangle side irv1-irv4 and the model x axis)
    sinrot = dy4 / dx
    cosrot = dx4 / dx
    !
    return
    !
  end subroutine get_rectDimensionsRotation

  !> @brief Return a rectangle face flow
  function get_rectflow(this, iq, irv) result(rectflow)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    integer :: iq, irv
    double precision :: rectflow
    !
    rectflow = this%rectflow(iq, irv)
    !
    return
  end function get_rectflow

  !> @brief Return whether a rectangle face is refined
  function isRefinedFace(this, irv)
    ! -- dummy
    class(CellRectQuadType), intent(inout) :: this
    integer :: irv
    logical :: isRefinedFace
    !
    if (this%ipv4irv(2, irv) .ne. 0) then
      isRefinedFace = .true.
    else
      isRefinedFace = .false.
    end if
    !
    return
  end function isRefinedFace

end module CellRectQuadModule
