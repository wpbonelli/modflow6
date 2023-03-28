module CellModule

  use CellDefnModule, only: CellDefnType
  implicit none
  private
  public :: CellType

  !> @brief A grid cell. Contains a cell-definition (composition over inheritance)
  type, abstract :: CellType
    character(len=40), pointer :: type ! tracking domain type
    type(CellDefnType), pointer :: defn => null() ! cell definition
  contains
    procedure(destroy), deferred :: destroy ! destructor for the cell
  end type CellType

  abstract interface
    subroutine destroy(this)
      import CellType
      class(CellType), intent(inout) :: this
    end subroutine
  end interface

end module CellModule
