module CellConnectionExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENEXCHANGENAME, LENMEMPATH
  use ListModule, only: ListType
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType

  implicit none

  private
  public :: CellConnectionExchangeType, AddCellConnectionExchangeToList, &
            GetCellConnectionExchangeFromList

  !> @brief CellConnectionExchangeType
  !!
  !! Base exchange type that stores the topological cell-to-cell connections
  !! between models. This type represents pure connectivity information without
  !! geometric data (connection lengths, areas, etc.). It stores:
  !!   - nexg: number of connected cell pairs
  !!   - nodem1, nodem2: node numbers in each model (reduced, 0-based)
  !!   - ihc: horizontal connection indicator (0=vertical, 1=horizontal)
  !!
  !! This is the foundation for more specialized exchange types that add
  !! geometric information (GeometricConnectionExchangeType) or numerical
  !! coupling methods (NumericalExchangeType).
  !<
  type, extends(BaseExchangeType) :: CellConnectionExchangeType
    integer(I4B), pointer :: nexg => null() !< number of connected cell pairs
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() !< node numbers in model 1 (reduced, 0-based)
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() !< node numbers in model 2 (reduced, 0-based)
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< horizontal connection indicator (0=vertical, 1=horizontal)
  contains
    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure :: allocate_scalars
    procedure :: allocate_arrays
  end type CellConnectionExchangeType

contains

  subroutine exg_df(this)
    class(CellConnectionExchangeType) :: this !< instance of cell connection exchange object
    ! override in child types
  end subroutine exg_df

  subroutine exg_ar(this)
    class(CellConnectionExchangeType) :: this !< instance of cell connection exchange object
    ! override in child types
  end subroutine exg_ar

  !> @brief Allocate scalar variables for this cell connection exchange
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(CellConnectionExchangeType) :: this !< instance of cell connection exchange object
    !
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    this%nexg = 0
  end subroutine allocate_scalars

  !> @brief Allocate array data for this cell connection exchange,
  !! using the number of connected nodes nexg
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(CellConnectionExchangeType) :: this !< instance of cell connection exchange object
    !
    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%ihc, this%nexg, 'IHC', this%memoryPath)
  end subroutine allocate_arrays

  !> @brief Deallocate memory associated with this cell connection exchange
  !<
  subroutine exg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(CellConnectionExchangeType) :: this !< instance of cell connection exchange object
    !
    call mem_deallocate(this%nodem1)
    call mem_deallocate(this%nodem2)
    call mem_deallocate(this%ihc)
    call mem_deallocate(this%nexg)
  end subroutine exg_da

  function CastAsCellConnectionExchangeClass(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(CellConnectionExchangeType), pointer :: res

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (CellConnectionExchangeType)
      res => obj
    end select
  end function CastAsCellConnectionExchangeClass

  subroutine AddCellConnectionExchangeToList(list, exchange)
    implicit none
    type(ListType), intent(inout) :: list
    class(CellConnectionExchangeType), pointer, intent(in) :: exchange
    class(*), pointer :: obj

    obj => exchange
    call list%Add(obj)
  end subroutine AddCellConnectionExchangeToList

  function GetCellConnectionExchangeFromList(list, idx) result(res)
    implicit none
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(CellConnectionExchangeType), pointer :: res
    class(*), pointer :: obj

    obj => list%GetItem(idx)
    res => CastAsCellConnectionExchangeClass(obj)
  end function GetCellConnectionExchangeFromList

end module CellConnectionExchangeModule
