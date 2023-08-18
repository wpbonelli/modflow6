module CellPolyModule

  use CellDefnModule, only: CellDefnType
  implicit none

  private
  public :: CellPolyType
  public :: create_cellPoly

  type :: CellPolyType
    private
    character(len=40), pointer, public :: type ! character string that names the tracking domain type
    type(CellDefnType), pointer, public :: cellDefn => null() ! pointer to basic cell definition
  contains
    procedure, public :: destroy => destroy_cellPoly ! destructor for the cell
    procedure, public :: init => init_cellPoly ! initializes the polygonal cell
  end type CellPolyType

contains

  !> @brief Create a new polygonal cell
  subroutine create_cellPoly(cellPoly)
    ! -- dummy
    type(CellPolyType), pointer :: cellPoly
    !
    allocate (cellPoly)
    allocate (cellPoly%cellDefn)
    allocate (cellPoly%type)
    cellPoly%type = 'CellPoly'
    !
    return
  end subroutine create_cellPoly

  !> @brief Destructor for a polygonal cell
  subroutine destroy_cellPoly(this)
    ! -- dummy
    class(CellPolyType), intent(inout) :: this
    !
    deallocate (this%cellDefn)
    deallocate (this%type)
    !
    return
  end subroutine destroy_cellPoly

  subroutine init_cellPoly(this)
    ! -- dummy
    class(CellPolyType), intent(inout) :: this
    !
    return
  end subroutine init_cellPoly

end module CellPolyModule
