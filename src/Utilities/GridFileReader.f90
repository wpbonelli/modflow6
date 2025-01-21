module GridFileReaderModule

  use KindModule
  use SimModule, only: store_error, store_error_unit
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH
  use BaseDisModule, only: DisBaseType
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use DisuModule, only: DisuType
  use InputOutputModule, only: urword, upcase, openfile
  use Integer1dReaderModule, only: read_int1d
  use Integer2dReaderModule, only: read_int2d
  use Double1dReaderModule, only: read_dbl1d
  use Double2dReaderModule, only: read_dbl2d
  use HashTableModule, only: HashTableType, hash_table_cr, hash_table_da
  use DblHashTableModule, only: DblHashTableType, &
                                dbl_hash_table_cr, dbl_hash_table_da
  use ArrayHandlersModule, only: ExpandArray

  implicit none

  public :: GridFileReaderType, read_grb

  interface read_grb
    module procedure :: &
      read_int_0d, read_int_1d, read_int_2d, read_int_3d, &
      read_dbl_0d, read_dbl_1d, read_dbl_2d, read_dbl_3d
  end interface read_grb

  type :: GridFileReaderType
    private
    ! unit
    integer(I4B), public :: inunit
    ! header
    character(len=10), public :: grid_type
    integer(I4B), public :: version
    integer(I4B) :: ntxt
    integer(I4B) :: lentxt
    ! index
    type(HashTableType), pointer :: int_values !< maps variable name to integer scalar
    type(DblHashTableType), pointer :: dbl_values !< maps variable name to double scalar
    type(HashTableType), pointer :: arr_indices !< maps array variable name to index in file
    type(HashTableType), pointer :: arr_types !< maps array variable name to type of array (1=int, 2=double)
    type(HashTableType), pointer :: arr_shapes_index !< maps variable name to index in shapes
    integer(I4B), allocatable :: arr_shapes(:) !< flat array of array variable shapes
  contains
    procedure :: initialize
    procedure :: read_header
    procedure :: build_index
    procedure :: finalize
  end type GridFileReaderType

contains

  subroutine initialize(this, iu)
    class(GridFileReaderType) :: this
    integer(I4B), intent(in) :: iu

    this%inunit = iu

    call hash_table_cr(this%int_values)
    call dbl_hash_table_cr(this%dbl_values)
    call hash_table_cr(this%arr_indices)
    call hash_table_cr(this%arr_types)
    call hash_table_cr(this%arr_shapes_index)
    allocate (this%arr_shapes(0))

    call this%read_header()
    call this%build_index()
  end subroutine initialize

  subroutine read_header(this)
    ! dummy
    class(GridFileReaderType) :: this
    ! local
    character(len=50) :: line
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval

    ! grid type
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    if (line(istart:istop) /= 'GRID') then
      call store_error('Binary grid file must begin with "GRID". '//&
                       &'Found: '//line(istart:istop))
      call store_error_unit(this%inunit)
    end if
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    this%grid_type = line(istart:istop)
    call upcase(this%grid_type)

    ! version
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%version = ival

    ! ntxt
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%ntxt = ival

    ! lentxt
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%lentxt = ival

    print *, 'grid_type: ', trim(this%grid_type), &
      ' version: ', this%version, ' ntxt: ', this%ntxt, &
      ' lentxt: ', this%lentxt

  end subroutine read_header

  subroutine build_index(this)
    class(GridFileReaderType) :: this
    ! local
    character(len=50) :: line
    character(len=10) :: key, dtype
    integer(I4B) :: i, lloc, istart, istop, ndim, dim
    integer(I4B) :: shp_idx
    integer(I4B) :: ival
    real(DP) :: rval
    integer(I4B), allocatable :: shp(:)

    do i = 1, this%ntxt
      ! read line
      read (this%inunit) line

      ! parse variable name
      lloc = 1
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      key = line(istart:istop)
      call upcase(key)

      ! parse data type
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      dtype = line(istart:istop)
      call upcase(dtype)

      ! parse dimensions
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
      ndim = ival
      if (allocated(shp)) deallocate (shp)
      allocate (shp(ndim))

      ! parse shape
      do dim = 1, ndim
        call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
        shp(dim) = ival
      end do

      print *, 'key: ', trim(key), ' dtype: ', trim(dtype), &
        ' ndim: ', ndim, ' shp: ', shp

      ! if scalar, read and store it.
      ! if array, store shape/dimensions.
      if (ndim == 0) then
        call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
        if (dtype == "INTEGER") then
          call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
          call this%int_values%add(key, ival)
        else if (dtype == "DOUBLE") then
          call urword(line, lloc, istart, istop, 3, ival, rval, 0, 0)
          call this%dbl_values%add(key, rval)
        end if
        cycle
      else
        call this%arr_indices%add(key, i)
        if (dtype == "INTEGER") then
          call this%arr_types%add(key, 1)
        else if (dtype == "DOUBLE") then
          call this%arr_types%add(key, 2)
        end if
        shp_idx = size(this%arr_shapes)
        call ExpandArray(this%arr_shapes, increment=ndim)
        this%arr_shapes(shp_idx + 1:shp_idx + ndim) = shp
        call this%arr_shapes_index%add(key, shp_idx + 1)
      end if
    end do

    rewind (this%inunit)
  end subroutine build_index

  subroutine finalize(this)
    class(GridFileReaderType) :: this
    close (this%inunit)
    call hash_table_da(this%int_values)
    call dbl_hash_table_da(this%dbl_values)
    call hash_table_da(this%arr_indices)
    call hash_table_da(this%arr_types)
    call hash_table_da(this%arr_shapes_index)
    deallocate (this%arr_shapes)
  end subroutine finalize

  subroutine read_int_0d(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), intent(out) :: v
    v = this%int_values%get(key)
  end subroutine read_int_0d

  subroutine read_dbl_0d(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), intent(out) :: v
    v = this%dbl_values%get(key)
  end subroutine read_dbl_0d

  subroutine read_int_1d(this, key, v)
    ! dummy
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), allocatable, intent(out) :: v(:)
    ! local
    integer(I4B) :: i, shp

    ! check type
    if (this%arr_types%get(key) /= 1) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not an integer array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    shp = this%arr_shapes_index%get(key)
    allocate (v(shp))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    ! read the array
    read (this%inunit) v
  end subroutine read_int_1d

  subroutine read_dbl_1d(this, key, v)
    ! dummy
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), allocatable, intent(out) :: v(:)
    ! local
    integer(I4B) :: i, shp

    ! check type
    if (this%arr_types%get(key) /= 2) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not a double array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    shp = this%arr_shapes_index%get(key)
    allocate (v(shp))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    ! read the array
    read (this%inunit) v
  end subroutine read_dbl_1d

  subroutine read_int_2d(this, key, v)
    ! dummy
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), allocatable, intent(out) :: v(:, :)
    ! local
    integer(I4B) :: i
    integer(I4B), allocatable :: shp(:)

    ! check type
    if (this%arr_types%get(key) /= 1) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not an integer array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    allocate (shp(2))
    do i = this%arr_shapes_index%get(key), i + 1
      print *, 'setting dim ', i, ' shape to ', this%arr_shapes(i)
      shp(i) = this%arr_shapes(i)
    end do
    allocate (v(shp(1), shp(2)))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    print *, 'reading ', trim(key), ' of shape ', shp

    ! read the array
    read (this%inunit) v(:,:)

    deallocate (shp)
  end subroutine read_int_2d

  subroutine read_dbl_2d(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), allocatable, intent(out) :: v(:, :)
    ! local
    integer(I4B) :: i
    integer(I4B), allocatable :: shp(:)

    ! check type
    if (this%arr_types%get(key) /= 2) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not a double array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    allocate (shp(2))
    do i = this%arr_shapes_index%get(key), i + 1
      shp(i) = this%arr_shapes(i)
    end do
    allocate (v(shp(1), shp(2)))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    ! read the array
    read (this%inunit) v(:,:)

    deallocate (shp)
  end subroutine read_dbl_2d

  subroutine read_int_3d(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), allocatable, intent(out) :: v(:, :, :)
    ! local
    integer(I4B) :: i
    integer(I4B), allocatable :: shp(:)

    ! check type
    if (this%arr_types%get(key) /= 1) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not an integer array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    allocate (shp(3))
    do i = this%arr_shapes_index%get(key), i + 1
      shp(i) = this%arr_shapes(i)
    end do
    allocate (v(shp(1), shp(2), shp(3)))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    ! read the array
    read (this%inunit) v

    deallocate (shp)
  end subroutine read_int_3d

  subroutine read_dbl_3d(this, key, v)
    ! dummy
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), allocatable, intent(out) :: v(:, :, :)
    ! local
    integer(I4B) :: i
    integer(I4B), allocatable :: shp(:)

    ! check type
    if (this%arr_types%get(key) /= 2) then
      write (errmsg, '(a, a, a)') &
        'Variable ', trim(key), ' is not a double array.'
      call store_error(errmsg)
    end if

    ! get shape and allocate
    allocate (shp(3))
    do i = this%arr_shapes_index%get(key), i + 1
      shp(i) = this%arr_shapes(i)
    end do
    allocate (v(shp(1), shp(2), shp(3)))

    ! skip text lines
    do i = 1, this%ntxt + 4
      read (this%inunit)
    end do

    ! skip preceding arrays
    do i = 1, this%arr_indices%get(key) - 1
      read (this%inunit)
    end do

    ! read the array
    read (this%inunit) v

    deallocate (shp)
  end subroutine read_dbl_3d

end module GridFileReaderModule
