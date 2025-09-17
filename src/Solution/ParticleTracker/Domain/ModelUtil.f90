module ModelUtilModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use MethodModule, only: MethodType
  use MathUtilModule, only: is_close
  implicit none

  private
  public :: cell_is_dry
  public :: cell_is_sat

contains

  !> @brief Determine if a cell is dry
  logical function cell_is_dry(method, ic, tol)
    use MethodModule, only: MethodType
    ! dummy
    type(MethodType), intent(inout), pointer :: method
    integer(I4B), intent(in) :: ic
    real(DP), intent(in), optional :: tol
    ! local
    real(DP) :: tolerance

    if (present(tol)) then
      tolerance = tol
    else
      tolerance = DZERO
    end if

    cell_is_dry = is_close( &
      method%fmi%gwfsat(ic), DZERO, &
      atol=tolerance, &
      symmetric=.false.)

  end function cell_is_dry

  !> @brief Determine if a cell is saturated
  logical function cell_is_sat(method, ic, tol)
    use MethodModule, only: MethodType
    ! dummy
    type(MethodType), intent(inout), pointer :: method
    integer(I4B), intent(in) :: ic
    real(DP), intent(in), optional :: tol
    ! local
    real(DP) :: tolerance

    if (present(tol)) then
      tolerance = tol
    else
      tolerance = DZERO
    end if

    cell_is_sat = is_close( &
      method%fmi%gwfsat(ic), DONE, &
      atol=tolerance, &
      symmetric=.false.)

  end function cell_is_sat