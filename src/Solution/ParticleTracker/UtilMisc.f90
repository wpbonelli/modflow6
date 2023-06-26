module UtilMiscModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use CellDefnModule, only: VertexType

  interface allocate_as_needed
    module procedure allocate_as_needed_int1, allocate_as_needed_dble, &
      allocate_as_needed_vert, allocate_as_needed_logical
  end interface

contains

  function FirstNonBlank(string) result(firstChar)
    implicit none
    character(len=*) string
    integer :: firstChar, length, n

    firstChar = 0
    length = len_trim(string)
    if (length .gt. 0) then
      do n = 1, length
        if (string(n:n) .ne. ' ') then
          firstChar = n
          return
        end if
      end do
    end if

  end function FirstNonBlank

  subroutine GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
    implicit none
    integer, intent(in) :: cellNumber, nlay, nrow, ncol
    integer, intent(inout) :: layer, row, column
    integer :: nrnc, r

    nrnc = nrow * ncol
    row = nrow
    column = ncol
    layer = cellNumber / nrnc
    r = mod(cellNumber, nrnc)
    if (r .gt. 0) then
      layer = layer + 1
      row = r / ncol
      r = mod(r, ncol)
      if (r .gt. 0) then
        row = row + 1
        column = r
      end if
    end if

  end subroutine GetLayerRowColumn

  subroutine ModExt(a, p, r, n)
    implicit none
    integer, intent(in) :: a, p
    integer, intent(inout) :: r, n

    r = 0
    n = 0
    if (a .eq. 0) return

    n = a / p
    r = mod(a, p)
    if (r .eq. 0) n = n + 1

  end subroutine ModExt

  !> @brief Find the minimum between StartValue and EndValue
  INTEGER FUNCTION FindMinimum(x, StartValue, EndValue)
    IMPLICIT NONE
    INTEGER, DIMENSION(1:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: StartValue, EndValue
    INTEGER :: Minimum
    INTEGER :: Location
    INTEGER :: i

    Minimum = x(StartValue) ! assume the first is the min
    Location = StartValue ! record its position
    DO i = StartValue + 1, EndValue ! start with next elements
      IF (x(i) < Minimum) THEN !   if x(i) less than the min?
        Minimum = x(i) !      Yes, a new minimum found
        Location = i !      record its position
      END IF
    END DO
    FindMinimum = Location ! return the position
  END FUNCTION FindMinimum

  !> @brief Swaps the values of the two arguments
  SUBROUTINE Swap(a, b)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: a, b
    INTEGER :: Temp

    Temp = a
    a = b
    b = Temp
  END SUBROUTINE Swap

  !> @brief Sort the array in ascending order
  SUBROUTINE Sort(x, Size)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Size
    INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
    INTEGER :: i
    INTEGER :: Location

    DO i = 1, Size - 1 ! except for the last
      Location = FindMinimum(x, i, Size) ! find min from this to last
      CALL Swap(x(i), x(Location)) ! swap this and the minimum
    END DO
  END SUBROUTINE Sort

  !> @brief Allocates or reallocates/resizes an integer(kind=1)
  !! array to meet or exceed a specified minimum dimension.
  !! kluge note: is there a better way than these new "allocate_as_needed" subroutines???
  subroutine allocate_as_needed_int1(i1array, mindim)
    implicit none
    ! -- dummy
    integer(kind=1), allocatable, intent(inout) :: i1array(:)
    integer, intent(in) :: mindim
    !
    if (.not. ALLOCATED(i1array)) then
      allocate (i1array(mindim))
    else if (SIZE(i1array) .lt. mindim) then
      deallocate (i1array)
      allocate (i1array(mindim))
    end if
    !
    return
  end subroutine allocate_as_needed_int1

  !> @brief Allocates or reallocates/resizes a double-precision
  !! array to meet or exceed a specified minimum dimension
  subroutine allocate_as_needed_dble(darray, mindim)
    implicit none
    ! -- dummy
    double precision, allocatable, intent(inout) :: darray(:)
    integer, intent(in) :: mindim
    !
    if (.not. ALLOCATED(darray)) then
      allocate (darray(mindim))
    else if (SIZE(darray) .lt. mindim) then
      deallocate (darray)
      allocate (darray(mindim))
    end if
    !
    return
  end subroutine allocate_as_needed_dble

  !> @brief Allocates or reallocates/resizes a vertex
  !! array to meet or exceed a specified minimum dimension
  subroutine allocate_as_needed_vert(varray, mindim)
    implicit none
    ! -- dummy
    class(VertexType), allocatable, intent(inout) :: varray(:)
    integer, intent(in) :: mindim
    !
    if (.not. ALLOCATED(varray)) then
      allocate (varray(mindim))
    else if (SIZE(varray) .lt. mindim) then
      deallocate (varray)
      allocate (varray(mindim))
    end if
    !
    return
  end subroutine allocate_as_needed_vert

  !> @brief Allocates or reallocates/resizes a logical
  !! array to meet or exceed a specified minimum dimension
  subroutine allocate_as_needed_logical(larray, mindim)
    implicit none
    ! -- dummy
    logical, allocatable, intent(inout) :: larray(:)
    integer, intent(in) :: mindim
    !
    if (.not. ALLOCATED(larray)) then
      allocate (larray(mindim))
    else if (SIZE(larray) .lt. mindim) then
      deallocate (larray)
      allocate (larray(mindim))
    end if
    !
    return
  end subroutine allocate_as_needed_logical

  !> @brief 3D translation and 2D rotation of coordinates
  !! with option to invert transformation
  subroutine transform_coords(xin, yin, zin, xout, yout, zout, &
                              xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                              sinrot_opt, cosrot_opt, invert_opt)
    implicit none
    ! -- dummy
    double precision :: xin, yin, zin
    double precision :: xout, yout, zout
    double precision, optional :: xOrigin_opt, yOrigin_opt, zOrigin_opt
    double precision, optional :: sinrot_opt, cosrot_opt
    logical, optional :: invert_opt
    ! -- local
    logical(LGP) :: isTranslation_add, isRotation_add, invert_add
    real(DP) :: xOrigin_add, yOrigin_add, zOrigin_add, sinrot_add, cosrot_add
    double precision :: x, y
    !
    ! -- Process option arguments and set defaults and flags
    call transf_opt_args_prep(xOrigin_add, yOrigin_add, zOrigin_add, &
                              sinrot_add, cosrot_add, invert_add, &
                              isTranslation_add, isRotation_add, &
                              xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                              sinrot_opt, cosrot_opt, invert_opt)
    !
    ! -- Apply transformation or its inverse
    if (.not. invert_add) then
      ! -- Apply transformation to coordinates
      if (isTranslation_add) then
        xout = xin - xOrigin_add
        yout = yin - yOrigin_add
        zout = zin - zOrigin_add
      else
        xout = xOrigin_add
        yout = yOrigin_add
        zout = zOrigin_add
      end if
      if (isRotation_add) then
        x = xout
        y = yout
        xout = x * cosrot_add + y * sinrot_add
        yout = -x * sinrot_add + y * cosrot_add
      end if
    else
      ! -- Apply inverse of transformation to coordinates
      if (isRotation_add) then
        x = xin * cosrot_add - yin * sinrot_add
        y = xin * sinrot_add + yin * cosrot_add
      else
        x = xin
        y = yin
      end if
      if (isTranslation_add) then
        xout = x + xOrigin_add
        yout = y + yOrigin_add
        zout = zin + zOrigin_add
      end if
    end if
    !
    return
    !
  end subroutine transform_coords

  !> @brief Modify transformation by applying an additional
  !! 3D translation and 2D rotation
  subroutine modify_transf(xOrigin, yOrigin, zOrigin, sinrot, cosrot, &
                           xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                           sinrot_opt, cosrot_opt, invert_opt)
    implicit none
    ! -- dummy
    double precision :: xOrigin, yOrigin, zOrigin
    double precision :: sinrot, cosrot
    double precision, optional :: xOrigin_opt, yOrigin_opt, zOrigin_opt
    double precision, optional :: sinrot_opt, cosrot_opt
    logical, optional :: invert_opt
    ! -- local
    logical(LGP) :: isTranslation_add, isRotation_add, invert_add
    real(DP) :: xOrigin_add, yOrigin_add, zOrigin_add, sinrot_add, cosrot_add
    real(DP) :: x0, y0, z0, s0, c0
    !
    ! -- Process option arguments and set defaults and flags
    call transf_opt_args_prep(xOrigin_add, yOrigin_add, zOrigin_add, &
                              sinrot_add, cosrot_add, invert_add, &
                              isTranslation_add, isRotation_add, &
                              xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                              sinrot_opt, cosrot_opt, invert_opt)
    !
    ! -- Copy existing transformation into working copy
    x0 = xOrigin
    y0 = yOrigin
    z0 = zOrigin
    s0 = sinrot
    c0 = cosrot
    !
    ! -- Modify transformation
    if (.not. invert_add) then
      !
      ! -- Apply additional transformation to existing transformation
      !
      if (isTranslation_add) then
        ! -- Calculate modified origin, XOrigin + R^T XOrigin_add, where
        ! -- XOrigin and XOrigin_add are the existing and additional origin
        ! -- vectors, respectively, and R^T is the transpose of the existing
        ! -- rotation matrix
        call transform_coords(xOrigin_add, yOrigin_add, zOrigin_add, &
                              xOrigin, yOrigin, zOrigin, &
                              x0, y0, z0, s0, c0, .true.)
      end if
      if (isRotation_add) then
        ! -- Calculate modified rotation matrix (represented by sinrot
        ! -- and cosrot) as R_add R, where R and R_add are the existing
        ! -- and additional rotation matrices, respectively
        sinrot = cosrot_add * s0 + sinrot_add * c0
        cosrot = cosrot_add * c0 - sinrot_add * s0
      end if
      !
    else
      !
      ! -- Apply inverse of additional transformation to existing transformation
      !
      ! -- Calculate modified origin, R^T (XOrigin + R_add XOrigin_add), where
      ! -- XOrigin and XOrigin_add are the existing and additional origin
      ! -- vectors, respectively, R^T is the transpose of the existing rotation
      ! -- matrix, and R_add is the additional rotation matrix
      if (isTranslation_add) then
        call transform_coords(-xOrigin_add, -yOrigin_add, zOrigin_add, &
                              x0, y0, z0, xOrigin, yOrigin, zOrigin, &
                              -sinrot_add, cosrot_add, .true.)
      end if
      xOrigin = c0 * x0 - s0 * y0
      yOrigin = s0 * x0 + c0 * y0
      zOrigin = z0
      if (isRotation_add) then
        ! -- Calculate modified rotation matrix (represented by sinrot
        ! -- and cosrot) as R_add^T R, where R and R_add^T are the existing
        ! -- rotation matirx and the transpose of the additional rotation
        ! -- matrix, respectively
        sinrot = cosrot_add * s0 - sinrot_add * c0
        cosrot = cosrot_add * c0 + sinrot_add * s0
      end if
      !
    end if
    !
    return
    !
  end subroutine modify_transf

  !> @brief Process option arguments and set defaults and flags
  !! for transformation operations
  subroutine transf_opt_args_prep(xOrigin_add, yOrigin_add, zOrigin_add, &
                                  sinrot_add, cosrot_add, invert_add, &
                                  isTranslation_add, isRotation_add, &
                                  xOrigin_opt, yOrigin_opt, zOrigin_opt, &
                                  sinrot_opt, cosrot_opt, invert_opt)
    implicit none
    ! -- dummy
    double precision :: xOrigin_add, yOrigin_add, zOrigin_add
    double precision :: sinrot_add, cosrot_add
    logical :: invert_add, isTranslation_add, isRotation_add
    double precision, optional :: xOrigin_opt, yOrigin_opt, zOrigin_opt
    double precision, optional :: sinrot_opt, cosrot_opt
    logical, optional :: invert_opt
    ! -- local
    !
    isTranslation_add = .false.
    xOrigin_add = DZERO
    if (present(xOrigin_opt)) then
      xOrigin_add = xOrigin_opt
      isTranslation_add = .true.
    end if
    yOrigin_add = DZERO
    if (present(yOrigin_opt)) then
      yOrigin_add = yOrigin_opt
      isTranslation_add = .true.
    end if
    zOrigin_add = DZERO
    if (present(zOrigin_opt)) then
      zOrigin_add = zOrigin_opt
      isTranslation_add = .true.
    end if
    isRotation_add = .false.
    sinrot_add = DZERO
    cosrot_add = DONE
    if (present(sinrot_opt)) then
      sinrot_add = sinrot_opt
      if (present(cosrot_opt)) then
        cosrot_add = cosrot_opt
      else
        ! -- If sinrot_opt is specified but cosrot_opt is not,
        ! -- default to corresponding non-negative cosrot_add
        cosrot_add = dsqrt(DONE - sinrot_add * sinrot_add)
      end if
      isRotation_add = .true.
    else if (present(cosrot_opt)) then
      cosrot_add = cosrot_opt
      ! -- cosrot_opt is specified but sinrot_opt is not, so
      ! -- default to corresponding non-negative sinrot_add
      sinrot_add = dsqrt(DONE - cosrot_add * cosrot_add)
      isRotation_add = .true.
    end if
    invert_add = .false.
    if (present(invert_opt)) invert_add = invert_opt
    !
    return
    !
  end subroutine transf_opt_args_prep

  !> @brief Adds a specified amount to an index and "wraps" the result if needed
  function add_wrap(istart, iadd, ilimit) result(iwrapped)
    implicit none
    ! -- dummy
    integer :: istart, iadd, ilimit, iwrapped
    !
    ! -- Add iadd to istart.  If the result exceeds ilimit, wrap around
    ! -- starting with 1. It is assumed (without checking) that iadd
    ! -- is positive and that no more than one wrap around is needed to ensure
    ! -- the final result is not greater than ilimit. For this to be true,
    ! -- it is sufficient to have istart no greater than ilimit and iadd
    ! -- no greater than ilimit.
    iwrapped = istart + iadd
    if (iwrapped .gt. ilimit) iwrapped = iwrapped - ilimit
    !
    return
    !
  end function add_wrap

  !> @brief Subtracts a specified amount from an index and "wraps" the
  !! result if necessary
  function subtr_wrap(istart, isubtr, ilimit) result(iwrapped)
    implicit none
    ! -- dummy
    integer :: istart, isubtr, ilimit, iwrapped
    !
    ! -- Subtract isubtr from istart.  If the result is less than 1, wrap
    ! -- around starting with ilimit. It is assumed (without checking) that
    ! -- isubtr is positive and that no more than one wrap around is needed
    ! -- to ensure the final result is not less than 1. For this to be true,
    ! -- it is sufficient to have istart no less than 1 and isubtr
    ! -- no greater than ilimit.
    iwrapped = istart - isubtr
    if (iwrapped .lt. 1) iwrapped = iwrapped + ilimit
    !
    return
    !
  end function subtr_wrap

  !> @brief Increments an index, "wrapping" the result if necessary
  function incr_wrap(istart, ilimit) result(iwrapped)
    implicit none
    ! -- dummy
    integer :: istart, ilimit, iwrapped
    !
    ! -- Increment istart. If istart is at ilimit, wrap around to 1. It is
    ! -- assumed (without checking) that istart is not greater than ilimit.
    if (istart .ne. ilimit) then
      iwrapped = istart + 1
    else
      iwrapped = 1
    end if
    !
    return
    !
  end function incr_wrap

  !> @brief Decrements an index, "wrapping" the result if necessary
  function decr_wrap(istart, ilimit) result(iwrapped)
    implicit none
    ! -- dummy
    integer :: istart, ilimit, iwrapped
    !
    ! -- Decrement istart. If istart is 1, wrap around to ilimit. It is
    ! -- assumed (without checking) that istart is not less than 1.
    if (istart .ne. 1) then
      iwrapped = istart - 1
    else
      iwrapped = ilimit
    end if
    !
    return
    !
  end function decr_wrap

end module UtilMiscModule
