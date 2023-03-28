module CellPolyModule

!!  use CellModule
  use CellDefnModule, only: CellDefnType
  implicit none

  private
  public :: CellPolyType
  public :: create_cellPoly

!!  ! -- Extend CellType to the polygonal cell type (CellPolyType)
!!  type, extends(CellType) :: CellPolyType
  ! -- Define the polygonal cell type (CellPolyType)
  type :: CellPolyType
    private
    character(len=40), pointer, public :: type ! character string that names the tracking domain type
    type(CellDefnType), pointer, public :: cellDefn => null() ! pointer to basic cell definition
  contains
    procedure, public :: destroy => destroy_cellPoly ! destructor for the cell
    procedure, public :: init => init_cellPoly ! initializes the polygonal cell
  end type CellPolyType

contains

  subroutine create_cellPoly(cellPoly)
! ******************************************************************************
! create_cellPoly -- Create a new polygonal cell
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(CellPolyType), pointer :: cellPoly
! ------------------------------------------------------------------------------
    !
    allocate (cellPoly)
    allocate (cellPoly%cellDefn) ! kluge note: use create_cellDefn???
    allocate (cellPoly%type)
    cellPoly%type = 'CellPoly'
    !
    return
  end subroutine create_cellPoly

  subroutine destroy_cellPoly(this)
! ******************************************************************************
! destroy_cellPoly -- Destructor for a polygonal cell
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(CellPolyType), intent(inout) :: this
! ------------------------------------------------------------------------------
    !
    deallocate (this%cellDefn)
    deallocate (this%type)
    !
    return
  end subroutine destroy_cellPoly

  subroutine init_cellPoly(this)
    ! ******************************************************************************
! init_cellPoly -- Initialize a polygonal cell
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(CellPolyType), intent(inout) :: this
! ------------------------------------------------------------------------------
    !
    ! kluge note: needed???
    !
    return
  end subroutine init_cellPoly

end module CellPolyModule
