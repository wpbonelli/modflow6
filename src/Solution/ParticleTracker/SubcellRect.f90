module SubcellRectModule

  implicit none

  private
  public :: SubcellRectType
  public :: create_subcellRect

  type SubcellRectType
    private
    character(len=40), pointer, public :: type ! character string that names the tracking domain type
    integer, public :: isubcell ! index of subcell in the cell
    integer, public :: icell ! index of cell in the source grid
    double precision, public :: sinrot ! sine of rotation angle for local (x, y)
    double precision, public :: cosrot ! cosine of rotation angle for local (x, y)
    double precision, public :: xOrigin ! cell x origin for local (x, y)
    double precision, public :: yOrigin ! cell y origin for local (x, y)
    double precision, public :: zOrigin ! cell z origin for local z
    double precision, public :: dx, dy, dz ! subcell dimensions
    double precision, public :: vx1, vx2, vy1, vy2, vz1, vz2 ! subcell face velocities
  contains
    procedure, public :: destroy => destroy_subcellRect ! destructor for the subcell
    procedure, public :: init => init_subcellRect ! initializes the rectangular subcell
  end type SubcellRectType

contains

  !> @brief Create a new rectangular subcell
  subroutine create_subcellRect(subcellRect)
    ! -- dummy
    type(SubcellRectType), pointer :: subcellRect
    !
    allocate (subcellRect)
    allocate (subcellRect%type)
    subcellRect%type = 'SubcellRect'
    return
  end subroutine create_subcellRect

  !> @brief Destructor for a rectangular subcell
  subroutine destroy_subcellRect(this)
    ! -- dummy
    class(SubcellRectType), intent(inout) :: this
    !
    deallocate (this%type)
    return
  end subroutine destroy_subcellRect

  !> @brief Initialize a rectangular subcell
  subroutine init_subcellRect(this)
    ! -- dummy
    class(SubcellRectType), intent(inout) :: this
    !
    return
  end subroutine init_subcellRect

end module SubcellRectModule
