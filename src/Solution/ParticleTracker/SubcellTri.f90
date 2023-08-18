module SubcellTriModule

  implicit none

  private
  public :: SubcellTriType
  public :: create_subcellTri

  type SubcellTriType
    private
    character(len=40), pointer, public :: type ! character string that names the tracking domain type
    integer, public :: isubcell ! index of subcell in the cell
    integer, public :: icell ! index of cell in the source grid
    double precision, public :: x0, y0, x1, y1, x2, y2 ! subcell corner coordinates
    double precision, public :: v0x, v0y, v1x, v1y, v2x, v2y ! subcell corner velocities
    double precision, public :: ztop, zbot ! subcell top and bottom elevations
    double precision, public :: dz ! subcell thickness
    double precision, public :: vztop, vzbot ! subcell top and bottom velocities
  contains
    procedure, public :: destroy => destroy_subcellTri ! destructor for the subcell
    procedure, public :: init => init_subcellTri ! initializes the triangular subcell
  end type SubcellTriType

contains

  !> @brief Create a new triangular subcell
  subroutine create_subcellTri(subcellTri)
    ! -- dummy
    type(SubcellTriType), pointer :: subcellTri
    !
    allocate (subcellTri)
    allocate (subcellTri%type)
    subcellTri%type = 'SubcellTri'
    !
    return
  end subroutine create_subcellTri

  !> @brief Destructor for a triangular subcell
  subroutine destroy_subcellTri(this)
    ! -- dummy
    class(SubcellTriType), intent(inout) :: this
    !
    deallocate (this%type)
    !
    return
  end subroutine destroy_subcellTri

  !> @brief Initialize a triangular subcell
  subroutine init_subcellTri(this)
    ! -- dummy
    class(SubcellTriType), intent(inout) :: this
    !
    return
  end subroutine init_subcellTri

end module SubcellTriModule
