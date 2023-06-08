module UtilMiscModule

  use CellDefnModule, only: VertexType

  interface allocate_as_needed
    module procedure allocate_as_needed_int1, allocate_as_needed_dble, &
      allocate_as_needed_vert, allocate_as_needed_logical
  end interface

contains

  function FirstNonBlank(string) result(firstChar)
    implicit none
    character(len=*), intent(in) :: string
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

!------------------------------------------------------
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

!------------------------------------------------------
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

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between StartValue and EndValue.
! --------------------------------------------------------------------

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

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

  SUBROUTINE Swap(a, b)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: a, b
    INTEGER :: Temp

    Temp = a
    a = b
    b = Temp
  END SUBROUTINE Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

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

  subroutine allocate_as_needed_int1(i1array, mindim)
! ******************************************************************************
! allocate_as_needed_int1 -- Allocates or reallocates/resizes an integer(kind=1)
!                            array to meet or exceed a specified minimum
!                            dimension
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    integer(kind=1), allocatable, intent(inout) :: i1array(:)
    integer, intent(in) :: mindim
! ------------------------------------------------------------------------------
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

  subroutine allocate_as_needed_dble(darray, mindim)
! ******************************************************************************
! allocate_as_needed_dble -- Allocates or reallocates/resizes a double-precision
!                            array to meet or exceed a specified minimum
!                            dimension
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    double precision, allocatable, intent(inout) :: darray(:)
    integer, intent(in) :: mindim
! ------------------------------------------------------------------------------
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

  subroutine allocate_as_needed_vert(varray, mindim)
! ******************************************************************************
! allocate_as_needed_vert -- Allocates or reallocates/resizes a vertex
!                            array to meet or exceed a specified minimum
!                            dimension
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    class(VertexType), allocatable, intent(inout) :: varray(:)
    integer, intent(in) :: mindim
! ------------------------------------------------------------------------------
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

  subroutine allocate_as_needed_logical(larray, mindim)
! ******************************************************************************
! allocate_as_needed_logical -- Allocates or reallocates/resizes a logical
!                               array to meet or exceed a specified minimum
!                               dimension
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    logical, allocatable, intent(inout) :: larray(:)
    integer, intent(in) :: mindim
! ------------------------------------------------------------------------------
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

  subroutine transform_coords(xin, yin, zin, xOrigin, yOrigin, zOrigin, &
                              sinrot, cosrot, invert, xout, yout, zout)
! ******************************************************************************
! transform_coords -- Translation and 2D rotation of coordinates,
!                     with option to invert transformation
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    double precision :: xin, yin, zin
    double precision :: xOrigin, yOrigin, zOrigin
    double precision :: sinrot, cosrot
    logical :: invert
    double precision :: xout, yout, zout
    ! -- local
    double precision :: x, y
! ------------------------------------------------------------------------------
    !
    if (.not. invert) then
      xout = xin - xOrigin
      yout = yin - yOrigin
      zout = zin - zOrigin
      if (sinrot .ne. 0d0) then
        x = xout
        y = yout
        xout = x * cosrot + y * sinrot
        yout = -x * sinrot + y * cosrot
      end if
    else
      if (sinrot .ne. 0d0) then
        x = xin * cosrot - yin * sinrot
        y = xin * sinrot + yin * cosrot
      else
        x = xin
        y = yin
      end if
      xout = x + xOrigin
      yout = y + yOrigin
      zout = zin + zOrigin
    end if
    !
    return
    !
  end subroutine transform_coords

  function add_wrap(istart, iadd, ilimit) result(iwrapped)
! ******************************************************************************
! add_wrap -- Adds a specified amount to an index and "wraps" the result
!             if necessary
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    integer :: istart, iadd, ilimit, iwrapped
! ------------------------------------------------------------------------------
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

  function subtr_wrap(istart, isubtr, ilimit) result(iwrapped)
! ******************************************************************************
! subtr_wrap -- Subtracts a specified amount from an index and "wraps" the
!               result if necessary
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    integer :: istart, isubtr, ilimit, iwrapped
! ------------------------------------------------------------------------------
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

  function incr_wrap(istart, ilimit) result(iwrapped)
! ******************************************************************************
! incr_wrap -- Increments an index, "wrapping" the result if necessary
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    integer :: istart, ilimit, iwrapped
! ------------------------------------------------------------------------------
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

  function decr_wrap(istart, ilimit) result(iwrapped)
! ******************************************************************************
! decr_wrap -- Decrements an index, "wrapping" the result if necessary
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    ! -- dummy
    integer :: istart, ilimit, iwrapped
! ------------------------------------------------------------------------------
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
