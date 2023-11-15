module MathUtilModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: MAXCHARLEN, LENHUGELINE, &
                             DZERO, DPREC, DSAME, &
                             LINELENGTH, LENHUGELINE, VSUMMARY

  implicit none
  private

  public :: is_same
  public :: is_close

contains

  !> @brief Function to determine if two reals are the same
  !!
  !! Function to evaluate if the difference between a and b are less than eps
  !! (i.e. a and b are the same).
  !!
  !<
  function is_same(a, b, eps) result(lvalue)
    ! -- return variable
    logical(LGP) :: lvalue !< boolean indicating if a and b are the same
    ! -- dummy variables
    real(DP), intent(in) :: a !< first number to evaluate
    real(DP), intent(in) :: b !< second number to evaluate
    real(DP), intent(in), optional :: eps !< optional maximum difference between a abd b (default=DSAME)
    ! -- local variables
    real(DP) :: epsloc
    real(DP) :: denom
    real(DP) :: rdiff
    !
    ! -- evaluate optioanl arguments
    if (present(eps)) then
      epsloc = eps
    else
      epsloc = DSAME
    end if
    lvalue = .FALSE.
    if (a == b) then
      lvalue = .TRUE.
    else
      if (abs(b) > abs(a)) then
        denom = b
      else
        denom = a
        if (abs(denom) == DZERO) then
          denom = DPREC
        end if
      end if
      rdiff = abs((a - b) / denom)
      if (rdiff <= epsloc) then
        lvalue = .TRUE.
      end if
    end if
    !
    ! -- return
    return
  end function is_same

  !> @brief Check if a value is approximately equal to another.
  !!
  !! By default this is asymmetric in a and b, with b taken to be
  !! the reference value, and the result computed by the formula:
  !! (abs(a - b) <= (atol + rtol * abs(b))). If symmetric is true
  !! then an alternate formula for symmetric equivalence is used:
  !! abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol). Default
  !! values for rtol and atol are 1e-05 and 0, respectively. Note
  !! that an atol of 0 is only suitable if a and b are not near 0.
  !<
  logical function is_close(a, b, rtol, atol, symmetric)
    ! modules
    use ConstantsModule, only: DSAME
    ! dummy
    real(DP), intent(in) :: a, b
    real(DP), intent(in), optional :: rtol, atol
    logical(LGP), intent(in), optional :: symmetric
    ! local
    real(DP) :: lrtol, latol
    logical(LGP) :: lsymmetric

    if (.not. present(rtol)) then
      lrtol = 1d-5
    else
      lrtol = rtol
    end if

    if (.not. present(atol)) then
      latol = DSAME
    else
      latol = atol
    end if

    if (.not. present(symmetric)) then
      lsymmetric = .false.
    else
      lsymmetric = symmetric
    end if

    if (lsymmetric) then
      ! symmetric, https://peps.python.org/pep-0485/
      is_close = abs(a - b) <= max(lrtol * max(abs(a), abs(b)), latol)
    else
      ! asymmetric, https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
      is_close = (abs(a - b) <= (latol + lrtol * abs(b)))
    end if
  end function is_close

end module MathUtilModule
