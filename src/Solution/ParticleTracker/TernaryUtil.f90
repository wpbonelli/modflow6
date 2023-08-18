module TernaryUtil
  ! implicit none
  private
  public :: rotate
  public :: skew
  ! public :: num_vertices

contains
  !> @brief Rotate a vector
  subroutine rotate(rxx, rxy, ryx, ryy, xold, yold, xnew, ynew)
    implicit double precision(a - h, o - z)
    !
    xnew = rxx * xold + rxy * yold
    ynew = ryx * xold + ryy * yold
    !
    return
    !
  end subroutine

  !> @brief Skew a vector
  subroutine skew(iskew, sxx, sxy, syy, alp, bet)
    implicit double precision(a - h, o - z)
    !
    alpold = alp
    betold = bet
    if (iskew .eq. 1) then
      alp = sxx * alpold + sxy * betold
      bet = syy * betold
    else
      bet = betold / syy
      alp = (alpold - sxy * bet) / sxx
    end if
    !
    return
    !
  end subroutine

  !> @brief Count the number of vertices in a polygon
  ! integer function num_vertices(numvermx, ncpl, icell)
  !   use Ternary, only: ivert_polygon
  !   implicit double precision(a - h, o - z)
  !   !
  !   numver = numvermx
  !   do while (ivert_polygon(icell, numver) .eq. 0)
  !     numver = numver - 1
  !   end do
  !   !
  !   num_vertices = numver
  !   !
  !   return
  !   !
  ! end function

end module TernaryUtil
