module CellDefnModule

  implicit none

  private
  public :: CellDefnType
  public :: VertexType ! kluge public???
  public :: create_cellDefn

  ! -- Define the vertex type (VertexType)
  type VertexType
    integer :: ivert ! vertex number (index)
    double precision :: x ! x coord of vertex
    double precision :: y ! y coord of vertex
  end type VertexType

  ! -- Define the modflow cell definition type (CellDefnType)
  type CellDefnType
    private
    integer, public :: icell ! index of cell in source grid
    logical, public :: canBeCellRect ! flag to indicate whether cell can be represented as CellRectType
    logical, public :: canBeCellRectQuad ! flag to indicate whether cell can be represented as CellRectQuadType
    integer, public :: npolyverts ! number of vertices for cell polygon
    double precision, public :: porosity ! cell porosity
    double precision, public :: retfactor ! cell retardation factor
    integer, public :: izone ! cell zone number
    integer, public :: iweaksink ! weak sink indicator
    integer, public :: inoexitface ! no exit face indicator
    integer, public :: iatop ! kluge???               ! index of cell top in grid's top/bot arrays (<0 => top array)
    double precision, public :: top, bot ! top and bottom elevations of cell
    double precision, public :: sat ! cell saturation
    class(VertexType), allocatable, public :: polyvert(:) ! vertices for cell polygon
    logical, allocatable, public :: ispv180(:) ! indicator of 180-degree vertices (.true. = 180-degree angle at vertex)
    integer(kind=1), allocatable, public :: facenbr(:) ! neighbors that correspond to faces(/vertices)
    double precision, allocatable, public :: faceflow(:) ! flows that correspond to faces(/vertices)
    double precision, public :: distflow ! net distributed flow into cell
  contains
    procedure, public :: destroy => destroy_cellDefn ! destructor for the cell definition
    procedure, public :: init => init_cellDefn ! initializes the cell definition
    procedure, public :: get_npolyverts ! returns the number of polygon vertices
    procedure, public :: get_ispv180 ! returns 180-degree indicator for a vertex
    procedure, public :: get_botflow ! returns bottom flow
    procedure, public :: get_topflow ! returns top flow
    procedure, public :: get_distflow ! returns distributed flow
    procedure, public :: get_faceflow ! returns a face flow
  end type CellDefnType

contains

  !> @brief Create a new cell definition object
  subroutine create_cellDefn(cellDefn)
    ! -- dummy
    type(CellDefnType), pointer :: cellDefn
    !
    allocate (cellDefn)
    return
  end subroutine create_cellDefn

  !> @brief Destructor for a cell definition object
  subroutine destroy_cellDefn(this)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    !
    return
  end subroutine destroy_cellDefn

  !> @brief Initialize a cell definition object
  !! kluge note: needed???
  subroutine init_cellDefn(this)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    !
    return
  end subroutine init_cellDefn

  !> @brief Return the number of polygon vertices
  function get_npolyverts(this) result(npolyverts)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    integer :: npolyverts
    !
    npolyverts = this%npolyverts
    return
  end function get_npolyverts

  !> @brief Return 180-degree indicator for a vertex
  function get_ispv180(this, m) result(ispv180)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    integer :: m
    logical :: ispv180
    !
    ispv180 = this%ispv180(m)
    return
  end function get_ispv180

  !> @brief Return the bottom flow
  function get_botflow(this) result(botflow)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    double precision :: botflow
    !
    botflow = this%faceflow(this%npolyverts + 2)
    return
  end function get_botflow

  !> @brief Return the top flow
  function get_topflow(this) result(topflow)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    double precision :: topflow
    !
    topflow = this%faceflow(this%npolyverts + 3)
    return
  end function get_topflow

  !> @brief Return the distributed flow
  function get_distflow(this) result(distflow)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    double precision :: distflow
    !
    distflow = this%distflow
    return
  end function get_distflow

  !> @brief Return a face flow
  function get_faceflow(this, m) result(faceflow)
    ! -- dummy
    class(CellDefnType), intent(inout) :: this
    integer :: m
    double precision :: faceflow
    !
    faceflow = this%faceflow(m)
    return
  end function get_faceflow

end module CellDefnModule
