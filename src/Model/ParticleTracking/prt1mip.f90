module PrtMipModule
  !
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE
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
    procedure, private :: read_options
    procedure :: read_data
  end type PrtMipType

contains

  subroutine mip_cr(mip, name_model, inunit, iout, dis)
! ******************************************************************************
! mip_cr -- Create a model input object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(PrtMipType), pointer :: mip
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
! ------------------------------------------------------------------------------
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
    mip%inunit = inunit
    mip%iout = iout
    !
    ! -- Set pointers
    mip%dis => dis
    !
    ! -- Initialize block parser
    call mip%parser%Initialize(mip%inunit, mip%iout)
    !
    ! -- Return
    return
  end subroutine mip_cr

  subroutine mip_da(this)
! ******************************************************************************
! mip_da -- deallocate
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(PrtMipType) :: this
! ------------------------------------------------------------------------------
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

  subroutine allocate_arrays(this, nodes)
! ******************************************************************************
! allocate_arrays
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PrtMipType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
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
  !!
  !!  Method to allocate and read model input for the MIP package.
  !!
  !<
  subroutine mip_ar(this)
    ! -- dummy variables
    class(PrtMipType), intent(inout) :: this !< PrtMipType object
    ! -- local variables
    !
    ! -- Print a message identifying the model input package.
    write (this%iout, 1) this%inunit
1   format(1x, /1x, 'MIP -- MODEL INPUT PACKAGE, VERSION X, X/XX/XXXX', & ! kluge update
           ' INPUT READ FROM UNIT ', i0)
    !
    ! -- Allocate arrays
    call this%allocate_arrays(this%dis%nodes)
    !
    ! -- Read options
    call this%read_options()
    !
    ! -- Read data
    call this%read_data()
    !
    ! -- return
    return
  end subroutine mip_ar

  subroutine read_options(this)
! ******************************************************************************
! read_options
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error
    ! -- dummy
    class(PrtMipType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
! ------------------------------------------------------------------------------
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING IC OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case default
          write (errmsg, '(4x,a,a)') 'Unknown IC option: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF IC OPTIONS'
    end if
    !
    ! -- Return
    return
  end subroutine read_options

  subroutine read_data(this)
! ******************************************************************************
! read_data
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error
    ! -- dummy
    class(PrtMipType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    character(len=:), allocatable :: line
    integer(I4B) :: istart, istop, lloc, ierr
    logical :: isfound, endOfBlock
! ------------------------------------------------------------------------------
    !
    ! -- get griddata block
    call this%parser%GetBlock('GRIDDATA', isfound, ierr)
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING GRIDDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        call this%parser%GetRemainingLine(line)
        lloc = 1
        select case (keyword)
        case ('POROSITY')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%porosity, &
                                        'POROSITY')
        case ('RETFACTOR')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%retfactor, &
                                        'RETFACTOR')
        case ('IZONE')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%izone, &
                                        'IZONE')
        case default
          write (errmsg, '(4x,a,a)') 'ERROR. UNKNOWN GRIDDATA TAG: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END PROCESSING GRIDDATA'
    else
      call store_error('ERROR.  REQUIRED GRIDDATA BLOCK NOT FOUND.')
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Return
    return
  end subroutine read_data

end module PrtMipModule
