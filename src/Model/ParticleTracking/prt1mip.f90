module PrtMipModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LINELENGTH
  use NumericalPackageModule, only: NumericalPackageType
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  !
  implicit none
  private
  public :: PrtMipType
  public :: mip_cr
  !
  type, extends(NumericalPackageType) :: PrtMipType
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< aquifer porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor => null() !< retardation factor
    integer(I4B), dimension(:), pointer, contiguous :: izone => null() !< zone number
  contains
    procedure :: mip_ar
    procedure :: mip_da
    procedure, private :: allocate_arrays
  end type PrtMipType

contains

  !> @brief Create a model input object
  subroutine mip_cr(mip, name_model, input_mempath, inunit, iout, dis)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    type(PrtMipType), pointer :: mip
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    ! -- locals
    logical(LGP) :: found_fname
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'NPF -- MODEL INPUT PACKAGE, VERSION 1, 08/08/2023', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! -- Create the object
    allocate (mip)
    !
    ! -- create name and memory path
    call mip%set_names(1, name_model, 'MIP', 'MIP')
    !
    ! -- Allocate scalars
    call mip%allocate_scalars()
    !
    ! -- Set variables
    mip%input_mempath = input_mempath
    mip%inunit = inunit
    mip%iout = iout
    !
    ! -- Set pointers
    mip%dis => dis
    !
    ! -- set name of input file
    call mem_set_value(mip%input_fname, 'INPUT_FNAME', mip%input_mempath, &
                       found_fname)
    !
    ! -- check if mip is enabled
    if (inunit > 0) then
      !
      ! -- Print a message identifying the model input package.
      write (iout, fmtheader) input_mempath
    end if
    !
    ! -- Return
    return
  end subroutine mip_cr

  !> @brief deallocate
  subroutine mip_da(this)
    ! -- modules
    use MemoryManagerExtModule, only: memorylist_remove
    use SimVariablesModule, only: idm_context
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(PrtMipType) :: this
    !
    ! -- Deallocate input memory
    call memorylist_remove(this%name_model, 'MIP', idm_context)
    !
    ! -- Deallocate parent package
    call this%NumericalPackageType%da()
    !
    ! -- scalars
    !
    ! -- arrays
    call mem_deallocate(this%porosity)
    call mem_deallocate(this%retfactor)
    call mem_deallocate(this%izone)
    !
    ! -- return
    return
  end subroutine mip_da

  !> @brief Allocate arrays
  subroutine allocate_arrays(this, nodes)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PrtMipType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    integer(I4B) :: i
    !
    ! -- Allocate
    call mem_allocate(this%porosity, nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%retfactor, nodes, 'RETFACTOR', this%memoryPath)
    call mem_allocate(this%izone, nodes, 'IZONE', this%memoryPath)
    !
    do i = 1, nodes
      this%porosity(i) = DZERO
      this%retfactor(i) = DONE
      this%izone(i) = 0
    end do
    !
    ! -- Return
    return
  end subroutine allocate_arrays

  !> @ brief Allocate and read model input
  subroutine mip_ar(this)
    ! -- modules
    use SimModule, only: store_error
    use MemoryManagerExtModule, only: mem_set_value
    use PrtMipInputModule, only: PrtMipParamFoundType
    ! -- dummy variables
    class(PrtMipType), intent(inout) :: this !< PrtMipType object
    ! -- local variables
    character(len=LINELENGTH) :: errmsg
    type(PrtMipParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map => null()
    !
    ! -- set map to convert user input data into reduced data
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! -- Allocate arrays
    call this%allocate_arrays(this%dis%nodes)
    !
    ! -- Source array inputs from IDM
    call mem_set_value(this%porosity, 'POROSITY', this%input_mempath, &
                       map, found%porosity)
    call mem_set_value(this%retfactor, 'RETFACTOR', this%input_mempath, &
                       map, found%retfactor)
    call mem_set_value(this%izone, 'IZONE', this%input_mempath, map, &
                       found%izone)
    !
    ! -- Ensure POROSITY was found
    if (.not. found%porosity) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: POROSITY not found'
      call store_error(errmsg)
    end if
    !
    ! -- return
    return
  end subroutine mip_ar

end module PrtMipModule
