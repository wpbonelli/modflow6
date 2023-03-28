module CellRectModule

  ! use CellModule
  use CellDefnModule, only: CellDefnType
  implicit none

  private
  ! public :: CellType
  public :: CellRectType
  public :: create_cellRect

  type :: CellRectType
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
    integer, public :: ipvOrigin ! origin vertex
    double precision, public :: vx1 ! west-boundary local-x velocity
    double precision, public :: vx2 ! east-boundary local-x velocity
    double precision, public :: vy1 ! south-boundary local-y velocity
    double precision, public :: vy2 ! north-boundary local-y velocity
    double precision, public :: vz1 ! bottom-boundary z velocity
    double precision, public :: vz2 ! top-boundary z velocity
  contains
    procedure, public :: destroy => destroy_cellRect ! destructor for the cell
    procedure, public :: init_cellRect ! initializes the rectangular cell
  end type CellRectType

contains

  !> @brief Create a new rectangular cell
  subroutine create_cellRect(cellRect)
    ! -- dummy
    type(CellRectType), pointer :: cellRect
    !
    allocate (cellRect)
    allocate (cellRect%cellDefn)
    allocate (cellRect%type)
    cellRect%type = 'CellRect'
    return
  end subroutine create_cellRect

  !> @brief Destructor for a rectangular cell
  subroutine destroy_cellRect(this)
    ! -- dummy
    class(CellRectType), intent(inout) :: this
    !
    deallocate (this%cellDefn)
    deallocate (this%type)
    return
  end subroutine destroy_cellRect

  !> @brief Initialize a rectangular cell
  !! kluge note: needed???
  subroutine init_cellRect(this)
    ! -- dummy
    class(CellRectType), intent(inout) :: this
    !
    return
  end subroutine init_cellRect

end module CellRectModule
