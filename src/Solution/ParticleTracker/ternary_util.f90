!
subroutine rotate(rxx, rxy, ryx, ryy, xold, yold, xnew, ynew)
!
  implicit double precision(a - h, o - z)
  !
  ! -- Rotate a vector
  xnew = rxx * xold + rxy * yold
  ynew = ryx * xold + ryy * yold
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine skew(iskew, sxx, sxy, syy, alp, bet)
!
  implicit double precision(a - h, o - z)
  !
  ! -- Skew a vector
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
!
! ------------------------------------------------------------------------------
!
integer function num_vertices(numvermx, ncpl, icell)
!
  use ternarymod, only: ivert_polygon
  implicit double precision(a - h, o - z)
  !
  ! -- Count the number of vertices in a polygon
  !
  numver = numvermx
  do while (ivert_polygon(icell, numver) .eq. 0)
    numver = numver - 1
  end do
  !
  num_vertices = numver
  !
  return
  !
end function

