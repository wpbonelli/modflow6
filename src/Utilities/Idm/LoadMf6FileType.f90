!> @brief This module contains the LoadMf6FileTypeModule
!!
!! This module contains the input data model routines for
!! loading the data from a MODFLOW 6 input file using the
!! block parser.
!!
!<
module LoadMf6FileTypeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENMEMPATH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use BlockParserModule, only: BlockParserType
  use ArrayReadersModule, only: ReadArray
  use InputOutputModule, only: parseline
  use InputDefinitionModule, only: InputParamDefinitionType
  use InputDefinitionSelectorModule, only: get_param_definition_type, &
                                           get_aggregate_definition_type
  use ModflowInputModule, only: ModflowInputType, getModflowInput
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use IdmLoggerModule, only: idm_log_var, idm_log_header, idm_log_close

  implicit none
  private
  public :: idm_load

  interface idm_load
    module procedure idm_load_from_blockparser
  end interface idm_load

contains

  !> @brief procedure to load a file
  !!
  !! Use parser to load information from an input file into the __INPUT__
  !! memory context location of the memory manager.
  !!
  !<
  subroutine idm_load_from_blockparser(parser, filetype, &
                                       component_type, subcomponent_type, &
                                       component_name, subcomponent_name, &
                                       subpackages, iout)
    use SimVariablesModule, only: idm_context
    type(BlockParserType), intent(inout) :: parser !< block parser
    character(len=*), intent(in) :: filetype !< file type to load, such as DIS6, DISV6, NPF6
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: component_name !< component name, such as MYGWFMODEL
    character(len=*), intent(in) :: subcomponent_name !< subcomponent name, such as MYWELLPACKAGE
    character(len=*), dimension(:), intent(in) :: subpackages !< array of subpackage types, such as ["TVK6", "OBS6"]
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B) :: iblock !< consecutive block number as defined in definition file
    type(ModflowInputType) :: mf6_input !< ModflowInputType
    character(len=LENMEMPATH) :: componentMemPath
    integer(I4B), dimension(:), contiguous, pointer :: mshape => null()
    !
    ! -- construct input object
    mf6_input = getModflowInput(filetype, component_type, &
                                subcomponent_type, component_name, &
                                subcomponent_name, subpackages)
    !
    ! -- model shape memory path
    componentMemPath = create_mem_path(component=mf6_input%component_name, &
                                       context=idm_context)
    !
    ! -- log lst file header
    call idm_log_header(mf6_input%component_name, &
                        mf6_input%subcomponent_name, iout)
    !
    ! -- process blocks
    do iblock = 1, size(mf6_input%p_block_dfns)
      call parse_block(parser, mf6_input, iblock, mshape, iout)
      !
      ! -- set model shape if discretion dimensions have been read
      if (mf6_input%p_block_dfns(iblock)%blockname == 'DIMENSIONS' .and. &
          filetype(1:3) == 'DIS') then
        call set_model_shape(mf6_input%file_type, componentMemPath, &
                             mf6_input%memoryPath, mshape)
      end if
    end do
    !
    ! -- close logging statement
    call idm_log_close(mf6_input%component_name, &
                       mf6_input%subcomponent_name, iout)
    !
    ! -- release allocated memory
    call mf6_input%destroy()
  end subroutine idm_load_from_blockparser

  !> @brief procedure to load a block
  !!
  !! Use parser to load information from a block into the __INPUT__
  !! memory context location of the memory manager.
  !!
  !<
  subroutine parse_block(parser, mf6_input, iblock, mshape, iout)
    use MemoryTypeModule, only: MemoryType
    use MemoryManagerModule, only: get_from_memorylist
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(ModflowInputType), intent(in) :: mf6_input !< ModflowInputType
    integer(I4B), intent(in) :: iblock !< consecutive block number as defined in definition file
    integer(I4B), dimension(:), contiguous, pointer, intent(inout) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    logical(LGP) :: isblockfound
    logical(LGP) :: endOfBlock
    logical(LGP) :: supportOpenClose
    integer(I4B) :: ierr
    logical(LGP) :: found
    type(MemoryType), pointer :: mt
    !
    ! -- disu vertices/cell2d blocks are contingent on NVERT dimension
    if (mf6_input%file_type == 'DISU6' .and. &
        (mf6_input%p_block_dfns(iblock)%blockname == 'VERTICES' .or. &
         mf6_input%p_block_dfns(iblock)%blockname == 'CELL2D')) then
      call get_from_memorylist('NVERT', mf6_input%memoryPath, mt, found, .false.)
      if (.not. found .or. mt%intsclr == 0) return
    end if
    !
    ! -- block open/close support
    supportOpenClose = (mf6_input%p_block_dfns(iblock)%blockname /= 'GRIDDATA')
    !
    ! -- parser search for block
    call parser%GetBlock(mf6_input%p_block_dfns(iblock)%blockname, isblockfound, &
                         ierr, supportOpenClose=supportOpenClose, &
                         blockRequired=mf6_input%p_block_dfns(iblock)%required)
    !
    ! -- process block
    if (isblockfound) then
      if (mf6_input%p_block_dfns(iblock)%aggregate) then
        !
        ! -- process block recarray type, set of variable 1d/2d types
        call parse_structarray_block(parser, mf6_input, iblock, mshape, iout)
      else
        do
          ! process each line in block
          call parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          !
          ! -- process line as tag(s)
          call parse_tag(parser, mf6_input, iblock, mshape, iout, .false.)
        end do
      end if
    end if

    return
  end subroutine parse_block

  !> @brief check subpackage
  !!
  !! Check and make sure that the subpackage is valid for
  !! this input file and load the filename of the subpackage
  !! into the memory manager.
  !!
  !<
  subroutine subpackage_check(parser, mf6_input, checktag, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(ModflowInputType), intent(in) :: mf6_input !< ModflowInputType
    character(len=LINELENGTH), intent(in) :: checktag !< subpackage string, such as TVK6
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=LINELENGTH) :: tag, fname_tag
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    integer(I4B) :: isubpkg

    do isubpkg = 1, size(mf6_input%subpackages)
      if (checktag == mf6_input%subpackages(isubpkg)) then
        fname_tag = trim(checktag)//'_FILENAME'
        call parser%GetStringCaps(tag)
        if (tag == 'FILEIN') then
          idt => get_param_definition_type(mf6_input%p_param_dfns, &
                                           mf6_input%component_type, &
                                           mf6_input%subcomponent_type, &
                                           fname_tag)
          call load_string_type(parser, idt, mf6_input%memoryPath, iout)
        else
          errmsg = 'Subpackage keyword must be followed by "FILEIN" '// &
                   'then by filename.'
          call store_error(errmsg)
        end if
      end if
    end do
  end subroutine subpackage_check

  !> @brief load an individual input record into memory
  !!
  !! Load an individual input record into the memory
  !! manager.  Allow for recursive calls in the case that multiple
  !! tags are on a single line.
  !!
  !<
  recursive subroutine parse_tag(parser, mf6_input, iblock, mshape, iout, &
                                 recursive_call)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(ModflowInputType), intent(in) :: mf6_input !< ModflowInputType
    integer(I4B), intent(in) :: iblock !< consecutive block number as defined in definition file
    integer(I4B), dimension(:), contiguous, pointer, intent(inout) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    logical(LGP), intent(in) :: recursive_call !< true if recursive call
    character(len=LINELENGTH) :: tag
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    !
    ! -- read tag name
    call parser%GetStringCaps(tag)
    if (recursive_call) then
      if (tag == '') then
        ! no data on line so return
        return
      end if
    end if
    !
    ! -- find keyword in input definition
    idt => get_param_definition_type(mf6_input%p_param_dfns, &
                                     mf6_input%component_type, &
                                     mf6_input%subcomponent_type, &
                                     tag)
    !
    ! -- allocate and load data type
    select case (idt%datatype)
    case ('KEYWORD')
      call load_keyword_type(parser, idt, mf6_input%memoryPath, iout)
      !
      ! -- load filename if subpackage tag
      call subpackage_check(parser, mf6_input, tag, iout)
      !
      ! -- set as dev option
      if (mf6_input%p_block_dfns(iblock)%blockname == 'OPTIONS' .and. &
          idt%tagname(1:4) == 'DEV_') then
        call parser%DevOpt()
      end if
    case ('STRING')
      call load_string_type(parser, idt, mf6_input%memoryPath, iout)
    case ('INTEGER')
      call load_integer_type(parser, idt, mf6_input%memoryPath, iout)
    case ('INTEGER1D')
      call load_integer1d_type(parser, idt, mf6_input%memoryPath, mshape, iout)
    case ('INTEGER3D')
      call load_integer3d_type(parser, idt, mf6_input%memoryPath, mshape, iout)
    case ('DOUBLE')
      call load_double_type(parser, idt, mf6_input%memoryPath, iout)
    case ('DOUBLE1D')
      call load_double1d_type(parser, idt, mf6_input%memoryPath, mshape, iout)
    case ('DOUBLE2D')
      call load_double2d_type(parser, idt, mf6_input%memoryPath, mshape, iout)
    case ('DOUBLE3D')
      call load_double3d_type(parser, idt, mf6_input%memoryPath, mshape, iout)
    case default
      write (errmsg, '(4x,a,a)') 'Failure reading data for tag: ', trim(tag)
      call store_error(errmsg)
      call parser%StoreErrorUnit()
    end select
    !
    ! -- continue line if in same record
    if (idt%in_record) then
      ! recursively call parse tag again to read rest of line
      call parse_tag(parser, mf6_input, iblock, mshape, iout, .true.)
    end if
    !
    ! --
    return
  end subroutine parse_tag

  !> @brief parse a structured array record into memory manager
  !!
  !! A structarray is similar to a numpy recarray.  It it used to
  !! load a list of data in which each column in the list may be a
  !! different type.  Each column in the list is stored as a 1d
  !! vector.
  !!
  !<
  subroutine parse_structarray_block(parser, mf6_input, iblock, mshape, iout)
    use StructArrayModule, only: StructArrayType, constructStructArray, &
                                 destructStructArray
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(ModflowInputType), intent(in) :: mf6_input !< ModflowInputType
    integer(I4B), intent(in) :: iblock !< consecutive block number as defined in definition file
    integer(I4B), dimension(:), contiguous, pointer, intent(inout) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    integer(I4B), pointer :: nrow
    integer(I4B) :: icol
    integer(I4B) :: ncol
    integer(I4B) :: nwords
    character(len=16), dimension(:), allocatable :: words
    type(StructArrayType), pointer :: struct_array
    character(len=:), allocatable :: parse_str
    !
    ! -- set input definition for this block
    idt => get_aggregate_definition_type(mf6_input%p_aggregate_dfns, &
                                         mf6_input%component_type, &
                                         mf6_input%subcomponent_type, &
                                         mf6_input%p_block_dfns(iblock)%blockname)
    !
    ! -- identify variable names, ignore first RECARRAY column
    parse_str = trim(idt%datatype)//' '
    call parseline(parse_str, nwords, words)
    ncol = nwords - 1
    !
    ! -- use shape to set the max num of rows
    call mem_setptr(nrow, idt%shape, mf6_input%memoryPath)
    !
    ! -- create a structured array
    struct_array => constructStructArray(ncol, nrow)
    do icol = 1, ncol
      !
      ! -- set pointer to input definition for this 1d vector
      idt => get_param_definition_type(mf6_input%p_param_dfns, &
                                       mf6_input%component_type, &
                                       mf6_input%subcomponent_type, &
                                       words(icol + 1))
      !
      ! -- allocate variable in memory manager
      call struct_array%mem_create_vector(icol, idt%datatype, idt%mf6varname, &
                                          mf6_input%memoryPath, idt%shape, &
                                          idt%preserve_case)
    end do
    !
    ! -- read the structured array
    call struct_array%read_from_parser(parser, iout)
    call parser%terminateblock()
    !
    ! -- destroy the structured array reader
    call destructStructArray(struct_array)
    !
    ! --
    return
  end subroutine parse_structarray_block

  !> @brief load type keyword
  !<
  subroutine load_keyword_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), pointer :: intvar
    call mem_allocate(intvar, idt%mf6varname, memoryPath)
    intvar = 1
    call idm_log_var(intvar, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_keyword_type

  !> @brief load type string
  !<
  subroutine load_string_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=LINELENGTH), pointer :: cstr
    integer(I4B) :: ilen
    ilen = LINELENGTH
    call mem_allocate(cstr, ilen, idt%mf6varname, memoryPath)
    call parser%GetString(cstr, (.not. idt%preserve_case))
    return
  end subroutine load_string_type

  !> @brief load type integer
  !<
  subroutine load_integer_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), pointer :: intvar
    call mem_allocate(intvar, idt%mf6varname, memoryPath)
    intvar = parser%GetInteger()
    call idm_log_var(intvar, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_integer_type

  !> @brief load type 1d integer
  !<
  subroutine load_integer1d_type(parser, idt, memoryPath, mshape, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B), pointer :: nsize1
    integer(I4B) :: nvals

    if (idt%shape == 'NODES') then
      nvals = product(mshape)
      call mem_allocate(int1d, nvals, idt%mf6varname, memoryPath)
    else
      call mem_setptr(nsize1, idt%shape, memoryPath)
      call mem_allocate(int1d, nsize1, idt%mf6varname, memoryPath)
    end if

    call read_grid_array(parser, mshape, idt%tagname, idt%layered, intarray=int1d)

    call idm_log_var(int1d, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_integer1d_type

  !> @brief load type 3d integer
  !<
  subroutine load_integer3d_type(parser, idt, memoryPath, mshape, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    integer(I4B) :: ndim
    integer(I4B) :: nsize1, nsize2, nsize3
    character(len=LINELENGTH) :: keyword

    ndim = size(mshape)

    ! set sizes
    if (ndim == 2) then
      nsize1 = mshape(2) ! NCPL
      nsize2 = 1
      nsize3 = mshape(1)
    elseif (ndim == 3) then
      nsize1 = mshape(3) ! NCOL
      nsize2 = mshape(2) ! NROW
      nsize3 = mshape(1) ! NLAY
    end if

    ! allocate the array using the memory manager
    call mem_allocate(int3d, nsize1, nsize2, nsize3, idt%mf6varname, memoryPath)

    ! fill the array from the file
    if (idt%blockname == 'GRIDDATA') then
      call parser%GetStringCaps(keyword)
      if (keyword == 'LAYERED') then
        ! read by layer
        call ReadArray(parser%iuactive, int3d(:, :, :), &
                       idt%mf6varname, ndim, nsize1, nsize2, &
                       nsize3, iout, 1, nsize3)
      else
        ! read full 3d array
        call ReadArray(parser%iuactive, int3d(:, :, :), idt%mf6varname, &
                       ndim, nsize1 * nsize2 * nsize3, iout)
      end if
    else
      ! read full 3d array
      call ReadArray(parser%iuactive, int3d(:, :, :), idt%mf6varname, &
                     ndim, nsize1 * nsize2 * nsize3, iout)
    end if

    call idm_log_var(int3d, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_integer3d_type

  !> @brief load type double
  !<
  subroutine load_double_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), pointer :: dblvar
    call mem_allocate(dblvar, idt%mf6varname, memoryPath)
    dblvar = parser%GetDouble()
    call idm_log_var(dblvar, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_double_type

  !> @brief load type 1d double
  !<
  subroutine load_double1d_type(parser, idt, memoryPath, mshape, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B), pointer :: nsize1
    integer(I4B) :: nvals

    if (idt%shape == 'NODES') then
      nvals = product(mshape)
      call mem_allocate(dbl1d, nvals, idt%mf6varname, memoryPath)
    else
      call mem_setptr(nsize1, idt%shape, memoryPath)
      call mem_allocate(dbl1d, nsize1, idt%mf6varname, memoryPath)
    end if

    call read_grid_array(parser, mshape, idt%tagname, idt%layered, dbl1d)
    call idm_log_var(dbl1d, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_double1d_type

  !> @brief load type 2d double
  !<
  subroutine load_double2d_type(parser, idt, memoryPath, mshape, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: ndim
    integer(I4B) :: nsize1, nsize2

    ndim = size(mshape)

    ! set sizes
    if (ndim == 2) then
      nsize1 = mshape(2) ! NCPL
      nsize2 = 1
    elseif (ndim == 3) then
      nsize1 = mshape(3) ! NCOL
      nsize2 = mshape(2) ! NROW
    end if

    ! allocate the array using the memory manager
    call mem_allocate(dbl2d, nsize1, nsize2, idt%mf6varname, memoryPath)

    ! fill the array from the file
    call ReadArray(parser%iuactive, dbl2d, idt%mf6varname, &
                   ndim, nsize1, nsize2, iout, 0)

    call idm_log_var(dbl2d, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_double2d_type

  !> @brief load type 3d double
  !<
  subroutine load_double3d_type(parser, idt, memoryPath, mshape, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    integer(I4B) :: ndim
    integer(I4B) :: nsize1, nsize2, nsize3
    character(len=LINELENGTH) :: keyword

    ndim = size(mshape)

    ! set sizes
    if (ndim == 2) then
      nsize1 = mshape(2) ! NCPL
      nsize2 = 1
      nsize3 = mshape(1)
    elseif (ndim == 3) then
      nsize1 = mshape(3) ! NCOL
      nsize2 = mshape(2) ! NROW
      nsize3 = mshape(1) ! NLAY
    end if

    ! allocate the array using the memory manager
    call mem_allocate(dbl3d, nsize1, nsize2, nsize3, idt%mf6varname, memoryPath)

    ! fill the array from the file
    if (idt%blockname == 'GRIDDATA') then
      call parser%GetStringCaps(keyword)
      if (keyword == 'LAYERED') then
        ! read by layer
        call ReadArray(parser%iuactive, dbl3d(:, :, :), &
                       idt%mf6varname, ndim, nsize1, nsize2, &
                       nsize3, iout, 1, nsize3)
      else
        ! read full 3d array
        call ReadArray(parser%iuactive, dbl3d(:, :, :), idt%mf6varname, &
                       ndim, nsize1 * nsize2 * nsize3, iout)
      end if
    else
      ! read full 3d array
      call ReadArray(parser%iuactive, dbl3d(:, :, :), idt%mf6varname, &
                     ndim, nsize1 * nsize2 * nsize3, iout)
    end if

    call idm_log_var(dbl3d, idt%mf6varname, memoryPath, iout)
    return
  end subroutine load_double3d_type

  !> @brief routine for setting the model shape
  !!
  !! The model shape must be set in the memory manager because
  !! individual packages need to know the shape of the arrays
  !! to read.
  !!
  !<
  subroutine set_model_shape(ftype, model_mempath, dis_mempath, model_shape)
    use MemoryTypeModule, only: MemoryType
    use MemoryManagerModule, only: get_from_memorylist
    character(len=*), intent(in) :: ftype
    character(len=*), intent(in) :: model_mempath
    character(len=*), intent(in) :: dis_mempath
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: model_shape
    integer(I4B), pointer :: ndim1
    integer(I4B), pointer :: ndim2
    integer(I4B), pointer :: ndim3

    select case (ftype)
    case ('DIS6')
      call mem_allocate(model_shape, 3, 'MODEL_SHAPE', model_mempath)
      call mem_setptr(ndim1, 'NLAY', dis_mempath)
      call mem_setptr(ndim2, 'NROW', dis_mempath)
      call mem_setptr(ndim3, 'NCOL', dis_mempath)
      model_shape = [ndim1, ndim2, ndim3]
    case ('DISV6')
      call mem_allocate(model_shape, 2, 'MODEL_SHAPE', model_mempath)
      call mem_setptr(ndim1, 'NLAY', dis_mempath)
      call mem_setptr(ndim2, 'NCPL', dis_mempath)
      model_shape = [ndim1, ndim2]
    case ('DISU6')
      call mem_allocate(model_shape, 1, 'MODEL_SHAPE', model_mempath)
      call mem_setptr(ndim1, 'NODES', dis_mempath)
      model_shape = [ndim1]
    end select

    return
  end subroutine set_model_shape

  !> @brief read an array that is the size of the model grid
  !<
  subroutine read_grid_array(parser, mshape, array_name, layered, dblarray, &
                             intarray)
    type(BlockParserType), intent(inout) :: parser !< block parser
    integer(I4B), dimension(:), intent(in) :: mshape !< model shape
    character(len=*), intent(in) :: array_name
    logical(LGP), intent(in) :: layered
    real(DP), dimension(:), optional, intent(inout) :: dblarray
    integer(I4B), dimension(:), optional, intent(inout) :: intarray
    integer(I4B) :: nvals
    integer(I4B) :: ndim
    integer(I4B) :: ndim1
    integer(I4B) :: ndim2
    integer(I4B) :: ndim3
    integer(I4B) :: k1
    integer(I4B) :: k2
    integer(I4B) :: iout !< unit number for output
    character(len=LINELENGTH) :: keyword

    ndim = size(mshape)
    if (present(dblarray)) then
      nvals = size(dblarray)
    end if
    if (present(intarray)) then
      nvals = size(intarray)
    end if
    iout = 0

    ! disu
    if (ndim == 1) then
      ndim1 = mshape(1) ! nodesuser
      ndim2 = 1 ! none
      ndim3 = 1 ! none
      k1 = 0
      k2 = 0

      ! disv
    else if (ndim == 2) then
      ndim1 = mshape(1) ! nlay
      ndim2 = 1 ! none
      ndim3 = mshape(2) ! ncpl
      k1 = 1
      k2 = ndim1

      ! dis
    else if (ndim == 3) then
      ndim1 = mshape(1) ! nlay
      ndim2 = mshape(2) ! nrow
      ndim3 = mshape(3) ! ncol
      k1 = 1
      k2 = ndim1
    end if

    call parser%GetStringCaps(keyword)
    if (keyword == 'LAYERED' .and. layered) then

      ! float array
      if (present(dblarray)) then
        call ReadArray(parser%iuactive, dblarray, &
                       array_name, ndim, ndim3, ndim2, &
                       ndim1, nvals, iout, k1, k2)
      end if

      ! integer array
      if (present(intarray)) then
        call ReadArray(parser%iuactive, intarray, &
                       array_name, ndim, ndim3, ndim2, &
                       ndim1, nvals, iout, k1, k2)
      end if

    else

      ! float array
      if (present(dblarray)) then
        call ReadArray(parser%iuactive, dblarray, array_name, &
                       ndim, nvals, iout, 0)
      end if

      ! integer array
      if (present(intarray)) then
        call ReadArray(parser%iuactive, intarray, array_name, &
                       ndim, nvals, iout, 0)
      end if

    end if

    return
  end subroutine read_grid_array

end module LoadMf6FileTypeModule