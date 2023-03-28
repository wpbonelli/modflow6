!
subroutine traverse_triangle(ntmax, nsave, diff, rdiff, isolv, tol, step, texit, &
    alpexit,betexit,itrifaceenter,itrifaceexit,rxx,rxy,ryx,ryy,sxx,sxy,syy,      &
                lbary, alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, vziodz, az)
!
  use ternarymod, only : xtrk,ytrk,ztrk,xtrktru,ytrktru,ztrktru,ca1,ca2,ca3,   &
   cb1, cb2, waa, wab, wba, wbb, alpp1, betp1, alppdiff, betpdiff, v0alp, v0bet, &
        v1alp, v1bet, v2alp, v2bet, icase, lenter, lsupout, tcpufindexit, nrptsoln
  implicit double precision(a - h, o - z)
  logical lface0, lbary
  intrinsic cpu_time
  !
  ! -- Compute elements of matrix W
  call get_w(alp0, bet0, alp1, bet1, alp2, bet2, lbary, waa, wab, wba, wbb)
  !
  ! -- Determine alpha and beta analytically as functions of time
  call solve_coefs(alpi, beti)
  if (.not. lsupout) then
    write (*, '(A,I2)') "   icase = ", icase
    write (69, '(A,I2)') "   icase = ", icase
  end if
  !
  ! -- Compute exit time (travel time to exit) and exit location
!!!  do irep=1,nrptsoln   ! kluge test
call find_exit_bary(isolv, 0, itrifaceenter, alp0, bet0, alp1, bet1, alpi, beti, &
                      rxx, rxy, ryx, ryy, tol, step, vziodz, az, &
                      texit0, alpexit0, betexit0)
call find_exit_bary(isolv, 1, itrifaceenter, alp1, bet1, alp2, bet2, alpi, beti, &
                      rxx, rxy, ryx, ryy, tol, step, vziodz, az, &
                      texit1, alpexit1, betexit1)
call find_exit_bary(isolv, 2, itrifaceenter, alp2, bet2, alp0, bet0, alpi, beti, &
                      rxx, rxy, ryx, ryy, tol, step, vziodz, az, &
                      texit2, alpexit2, betexit2)
  texit = min(texit0, texit1, texit2)
!!!  end do               ! kluge test
  ! -- Note that while the numbering of triangle faces is generally zero-based
  ! -- (0, 1, 2), itrifaceexit, which gets passed out, is one-based (1, 2, 3).
  if (texit .eq. texit0) then
    alpexit = alpexit0
    betexit = betexit0
    itrifaceexit = 1
  else if (texit .eq. texit1) then
    alpexit = alpexit1
    betexit = betexit1
    itrifaceexit = 2
  else if (texit .eq. texit2) then
    alpexit = alpexit2
    betexit = betexit2
    itrifaceexit = 3
  end if
  if (texit .eq. huge(1d0)) itrifaceexit = 0
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
logical function point_in_polygon(numvermx, nvert, ncpl, xi, yi, icell) ! kluge note: not used currently but will probably come in handy
  !
  use ternarymod, only: ivert_polygon, x_vert, y_vert
  implicit double precision(a - h, o - z)
  !
  ! -- Check whether point is in polygon
  point_in_polygon = .true.
  numver = num_vertices(numvermx, ncpl, icell)
  do iv = 1, numver
    iva = iv
    ivb = iv + 1
    if (ivb .gt. numver) ivb = 1
    ipva = ivert_polygon(icell, iva)
    ipvb = ivert_polygon(icell, ivb)
    xa = x_vert(ipva)
    ya = y_vert(ipva)
    xb = x_vert(ipvb)
    yb = y_vert(ipvb)
    term = (xb - xa) * (yi - ya) - (yb - ya) * (xi - xa)
    if (term .le. 0d0) then
      point_in_polygon = .false.
      exit
    end if
  end do
  !
  return
  !
end function
!
! ------------------------------------------------------------------------------
!
subroutine find_init_cell(numvermx,nvert,nlay,ncpl,xi,yi,zi,ilayeri,icelli)   ! kluge: move to "ternary_setup"???
!
  use ternarymod, only: ivert_polygon, x_vert, y_vert, z_bot
  implicit double precision(a - h, o - z)
  logical point_in_polygon
  external point_in_polygon
  !
  ! -- Find initial cell
  icelli = 0
  do icell = 1, ncpl
    if (point_in_polygon(numvermx, nvert, ncpl, xi, yi, icell)) then
      icelli = icell
      exit
    end if
  end do
  if (icelli .eq. 0) then
    write (*, '(A)') "error -- initial cell not found" ! kluge
    write (69, '(A)') "error -- initial cell not found"
    !!pause
    stop
  end if
  !
  ! -- Find initial layer
  ilayeri = 0
  if (zi .gt. z_bot(0, icelli)) then
    write (*, '(A)') "error -- initial z above top of model" ! kluge
    write (69, '(A)') "error -- initial z above top of model"
    !!pause
    stop
  end if
  do ilayer = 1, nlay
    if (zi .ge. z_bot(ilayer, icelli)) then
      ilayeri = ilayer
      exit
    end if
  end do
  if (ilayeri .eq. 0) then
    write (*, '(A)') "error -- initial z below bottom of model" ! kluge
    write (69, '(A)') "error -- initial z below bottom of model"
    !!pause
    stop
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine find_init_triangle(numvermx,nvert,ncpl,xi,yi,xctr,yctr,icell,itri)   ! kluge note: not currently used but maybe could be in load_subcell of MethodCellTernary???
!
  use ternarymod, only: ivert_polygon, x_vert, y_vert
  implicit double precision(a - h, o - z)
  !
  ! -- Find initial triangular subcell
  itri = 0
!!!  numver = num_vertices(numvermx,ncpl,icell)
  numver = 4 ! kluge
  do iv = 1, numver
    iv0 = iv
    iv1 = iv + 1
    if (iv1 .gt. numver) iv1 = 1
    ipv0 = ivert_polygon(icell, iv0)
    ipv1 = ivert_polygon(icell, iv1)
    x0 = x_vert(ipv0)
    y0 = y_vert(ipv0)
    x1 = x_vert(ipv1)
    y1 = y_vert(ipv1)
    x2 = xctr
    y2 = yctr
    x1rel = x1 - x0
    y1rel = y1 - y0
    x2rel = x2 - x0
    y2rel = y2 - y0
    di2 = xi * y2rel - yi * x2rel
    d02 = x0 * y2rel - y0 * x2rel
    d12 = x1rel * y2rel - y1rel * x2rel
    di1 = xi * y1rel - yi * x1rel
    d01 = x0 * y1rel - y0 * x1rel
    alphai = (di2 - d02) / d12
    betai = -(di1 - d01) / d12
if ((alphai .ge. 0d0) .and. (betai .ge. 0d0) .and. (alphai + betai .lt. 1d0)) then ! kluge note: think this handles points on triangle boundaries ok
      itri = iv ! but maybe not!!!!!!!!!!!!
      exit ! kluge note: doesn't handle particle smack on cell center
    end if
  end do
  if (itri .eq. 0) then
    write (*, '(A)') "error -- initial triangle not found" ! kluge
    write (69, '(A)') "error -- initial triangle not found"
    !!pause
    stop
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine canonical(x0, y0, x1, y1, x2, y2, v0x, v0y, v1x, v1y, v2x, v2y, &
  xi,yi,rxx,rxy,ryx,ryy,sxx,sxy,syy,lbary,alp0,bet0,alp1,bet1,alp2,bet2,alpi,beti)
!
  use ternarymod, only: v0alp, v0bet, v1alp, v1bet, v2alp, v2bet
  implicit double precision(a - h, o - z)
  logical lbary
  !
  ! -- Translate and rotate coordinates to "canonical" configuration
  x1diff = x1 - x0
  y1diff = y1 - y0
  x2diff = x2 - x0
  y2diff = y2 - y0
  baselen = dsqrt(x1diff * x1diff + y1diff * y1diff)
  oobaselen = 1d0 / baselen
  cosomega = x1diff * oobaselen
  sinomega = y1diff * oobaselen
  rxx = cosomega
  rxy = sinomega
  ryx = -sinomega
  ryy = cosomega
  alp0 = 0d0
  bet0 = 0d0
  alp1 = baselen
  bet1 = 0d0
  call rotate(rxx, rxy, ryx, ryy, x2diff, y2diff, alp2, bet2)
  call rotate(rxx, rxy, ryx, ryy, v0x, v0y, v0alp, v0bet)
  call rotate(rxx, rxy, ryx, ryy, v1x, v1y, v1alp, v1bet)
  call rotate(rxx, rxy, ryx, ryy, v2x, v2y, v2alp, v2bet)
  xidiff = xi - x0
  yidiff = yi - y0
  call rotate(rxx, rxy, ryx, ryy, xidiff, yidiff, alpi, beti)
  if (lbary) then
    sxx = 1d0 / alp1
    syy = 1d0 / bet2
    sxy = -alp2 * sxx * syy
    alp1 = 1d0
    alp2 = 0d0
    bet2 = 1d0
    call skew(1, sxx, sxy, syy, v0alp, v0bet)
    call skew(1, sxx, sxy, syy, v1alp, v1bet)
    call skew(1, sxx, sxy, syy, v2alp, v2bet)
    call skew(1, sxx, sxy, syy, alpi, beti)
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine get_w(alp0, bet0, alp1, bet1, alp2, bet2, lbary, waa, wab, wba, wbb)
!
  use ternarymod, only: v0alp, v0bet, v1alp, v1bet, v2alp, v2bet
  implicit double precision(a - h, o - z)
  logical lbary
  !
  ! -- Compute elements of W matrix
  ! -- Note: wab is the "alpha,beta" entry in matrix W
  !    and the alpha component of the w^(beta) vector
  v1alpdiff = v1alp - v0alp
!!!  v1betdiff = v1bet - v0bet
  v2alpdiff = v2alp - v0alp
  v2betdiff = v2bet - v0bet
  if (lbary) then
    waa = v1alpdiff
    wab = v2alpdiff
    wba = 0d0
    wbb = v2betdiff
  else
    ooalp1 = 1d0 / alp1
    oobet2 = 1d0 / bet2
    vterm = v1alpdiff * ooalp1
    waa = vterm
    wab = (v2alpdiff - alp2 * vterm) * oobet2
    wba = 0d0
    wbb = v2betdiff * oobet2
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine solve_coefs(alpi, beti)
!
  use ternarymod, only : ca1,ca2,ca3,cb1,cb2,waa,wab,wba,wbb,v0alp,v0bet,icase
  implicit double precision(a - h, o - z)
  !
  ! -- Compute analytical solution coefficients depending on case
  !
  zerotol = 1d-10 ! kluge
  !
  if (dabs(wbb) .gt. zerotol) then
    wratv = (wab / wbb) * v0bet
    acoef = v0alp - wratv
    bcoef = wratv + wab * beti
    afact = acoef / waa
    vfact = v0bet / wbb
    ! -- Coefs for beta do not depend on whether waa = 0 or not
    cb1 = -vfact ! const term in beta
    cb2 = vfact + beti ! coef for e(wbb*t) term in beta
    ! -- Coefs for alpha
    if (dabs(waa) .gt. zerotol) then
      ! -- Case waa <> 0, wbb <> 0
      if (dabs(wbb - waa) .gt. zerotol) then
        ! -- Subcase wbb <> waa
        bfact = bcoef / (wbb - waa)
        ca1 = -afact ! const term in alpha
        ca2 = alpi + afact - bfact ! coef for exp(waa*t) term in alpha
        ca3 = bfact ! coef for exp(wbb*t) term in alpha
        icase = 1
      else
        ! -- Subcase wbb = waa
        ca1 = -afact ! const term in alpha
        ca2 = alpi + afact ! coef for exp(waa*t) term in alpha
        ca3 = bcoef ! coef for t*exp(waa*t) term in alpha
        icase = -1
      end if
    else
      ! -- Case waa = 0, wbb <> 0
      bfact = bcoef / wbb
      ca1 = alpi - bfact ! const term in alpha
      ca2 = acoef ! coef for t term in alpha
      ca3 = bfact ! coef for exp(wbb*t) term in alpha
      icase = 2
    end if
  else
    ! -- Coefs for beta do not depend on whether waa = 0 or not
    cb1 = beti ! const term in beta
    cb2 = v0bet ! coef for t term in beta
    if (dabs(waa) .gt. zerotol) then
      ! -- Case waa <> 0, wbb = 0
      oowaa = 1d0 / waa
      vfact = (wab * oowaa) * v0bet
      ca1 = -oowaa * (v0alp + wab * beti + vfact) ! const term in alpha
      ca2 = -vfact ! coef for t term in alpha
      ca3 = alpi - ca1 ! coef for exp(waa*t) term in alpha
      icase = 3
    else
      ! -- Case waa = 0, wbb = 0
      ca1 = alpi ! const term in alpha
      ca2 = v0alp + wab * beti ! coef for t term in alpha
      ca3 = 5d-1 * wab * v0bet ! coef for t^2 term in alpha
      icase = 4
    end if
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine step_analytical(t, alp, bet)
!
  use ternarymod, only: ca1, ca2, ca3, cb1, cb2, waa, wab, wba, wbb, icase
  implicit double precision(a - h, o - z)
  !
  ! -- Step (evaluate) analytically depending on case
  !
  if (icase .eq. 1) then
    alp = ca1 + ca2 * dexp(waa * t) + ca3 * dexp(wbb * t)
    bet = cb1 + cb2 * dexp(wbb * t)
  else if (icase .eq. -1) then
    alp = ca1 + (ca2 + ca3 * t) * dexp(waa * t)
    bet = cb1 + cb2 * dexp(wbb * t)
  else if (icase .eq. 2) then
    alp = ca1 + ca2 * t + ca3 * dexp(wbb * t)
    bet = cb1 + cb2 * dexp(wbb * t)
  else if (icase .eq. 3) then
    alp = ca1 + ca2 * t + ca3 * dexp(waa * t)
    bet = cb1 + cb2 * t
  else if (icase .eq. 4) then
    alp = ca1 + (ca2 + ca3 * t) * t
    bet = cb1 + cb2 * t
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine step_euler(nt, step, vziodz, az, alpi, beti, t, alp, bet)
!
  use ternarymod, only: waa, wab, wba, wbb, v0alp, v0bet
  implicit double precision(a - h, o - z)
  !
  if (nt .eq. 0) then
    ! -- Initial location
    alp = alpi
    bet = beti
    t = 0d0
  else
    ! -- Step numerically
    valp = v0alp + waa * alp + wab * bet
    vbet = v0bet + wba * alp + wbb * bet
    if (step .lt. 0d0) then
      ! -- Compute time step based on abs value of step, interpreting the latter
      ! -- as a distance in canonical coordinates (alpha, beta, and scaled z)
      vz = vziodz * dexp(az * t)
      vmeasure = dsqrt(valp * valp + vbet * vbet + vz * vz)
      delt = -step / vmeasure
    else
      ! -- Set time step directly to step
      delt = step
    end if
    ikluge = 2 ! kluge
    if (ikluge .eq. 1) then
      t = t + delt
      alp = alp + valp * delt
      bet = bet + vbet * delt
    else
      rkn1 = valp
      rln1 = vbet
      thalf = t + 5d-1 * delt
      call step_analytical(thalf, alpproj, betproj)
      rkn2 = v0alp + waa * alpproj + wab * betproj
      rln2 = v0bet + wba * alpproj + wbb * betproj
      rkn3 = rkn2
      rln3 = rln2
      t = t + delt
      call step_analytical(t, alpproj, betproj)
      rkn4 = v0alp + waa * alpproj + wab * betproj
      rln4 = v0bet + wba * alpproj + wbb * betproj
      alp = alp + delt * (rkn1 + 2d0 * rkn2 + 2d0 * rkn3 + rkn4) / 6d0
      bet = bet + delt * (rln1 + 2d0 * rln2 + 2d0 * rln3 + rln4) / 6d0
    end if
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine find_exit_bary(isolv,itriface,itrifaceenter,alp1,bet1,alp2,bet2,alpi,beti,rxx,rxy,ryx,ryy,tol,step,vziodz,az,  &
                          texit, alpexit, betexit)
!
  use ternarymod, only : ca1,ca2,ca3,cb1,cb2,waa,wab,wba,wbb,v0alp,v0bet,v1alp,v1bet,v2alp,v2bet,  &
            alpp1, betp1, alpp2, betp2, alppdiff, betpdiff, icase, lenter, lsupout
  implicit double precision(a - h, o - z)
  logical lconstbet
  character convstrg * 11, facename(0:2) * 7, cfail * 60
  data(facename(itri), itri=0, 2)/"beta=0 ", "gamma=0", "alpha=0"/
  external fbary1, fbary2
  common / debug / ntdebug
  !
  ! -- Find the exit time and location in barycentric coordinates.
  ! -- Use iterative scheme or numerical integration indicated by isolv.
  !
  zerotol = 1d-10 ! kluge
!!!  texit = texitestim        ! kluge
  !
  if (itriface .eq. 0) then
    !
    ! -- Checking for exit on canonical face 0 (beta = 0)
    if (itrifaceenter .eq. 0) then
      ! -- Entrance face, so no exit. (Normal velocity is uniform along face 0,
      ! -- so it cannot be both an entrance and an exit.)
      texit = huge(1d0)
      cfail = "entered on face beta=0 (no exit)"
    else
      ! -- Not the entrance face, so check for outflow
      if (v0bet .ge. 0d0) then
        ! -- Inflow or no flow, so no exit
        texit = huge(1d0)
        cfail = "inflow or no flow on face beta=0 (no exit)"
      else
        ! -- Outflow, so check beta-velocity at the initial location,
        ! -- recognizing that it will never change sign along the
        ! -- trajectory (and will not be blocked from zero by an asymptote)
        vbeti = v0bet + wbb * beti
        if (vbeti .ge. 0d0) then
          ! -- Can't exit along beta = 0
          texit = huge(1d0)
          cfail = "cannot exit on face beta=0 (no exit)"
        else
          ! -- get alpt and check it
          call get_t_alpt(0d0, t, alpt)
          if ((alpt .ge. 0d0) .and. (alpt .le. 1d0)) then
            ! -- alpt within the edge, so exit found
            texit = t
            alpexit = alpt
            betexit = 0d0
            ntdebug = -111 ! kluge debug bludebug
          else
            ! -- alpt not within the edge, so not an exit
            texit = huge(1d0)
            cfail = "alpt not within face beta=0 (not an exit)"
          end if
        end if
      end if
    end if
    ! -- End canonical face 0 (beta = 0)
    !
  else
    !
    ! -- Checking for exit on canonical face 1 (gamma = 0.) or 2 (alpha = 0.)
    if (itriface .eq. 1) then
      ! -- Normal velocities (gamma components) at ends of canonical face 1
      v1n = -v1alp - v1bet
      v2n = -v2alp - v2bet
    else
      ! -- Normal velocities (alpha components) at ends of canonical face 2
      v1n = v0alp
      v2n = v2alp
    end if
    if ((v1n .ge. 0d0) .and. (v2n .ge. 0d0)) then
      ! -- No outflow at vn1 and vn2 corners; no outflow interval, so no exit.
      texit = huge(1d0)
      cfail = "no outflow interval (no exit)"
    else
      ! -- Find outflow interval
      call get_bet_outflow_bary(v1n, v2n, betoutlo, betouthi)
      ! -- Find trend of and limits on beta from beta{t} solution
      call get_bet_soln_limits(beti, betsollo, betsolhi, ibettrend)
      ! -- Look for exit
      if (ibettrend .eq. 0) then
        ! -- Beta is constant, so check if it's within the outflow interval;
        ! -- if not, no exit; if so, solve for t and alpha
        if ((beti .gt. betouthi) .or. (beti .lt. betoutlo)) then
          texit = huge(1d0)
          cfail = "beta=const not within outflow interval (not an exit)"
        else
          ! -- Check alpha-velocity at the initial location,
          ! -- recognizing that it will never change sign along the
          ! -- trajectory (and will not be blocked from zero by an asymptote)
          ! -- in this special case
          v0alpstar = v0alp + wab * beti
          valpi = v0alpstar + waa * alpi
          if ((itriface .eq. 1) .and. (valpi .le. 0d0)) then
            ! -- Can't exit along gamma = 0.
            texit = huge(1d0)
            cfail = "cannot exit on face gamma=0 (no exit)"
          else if ((itriface .eq. 2) .and. (valpi .ge. 0d0)) then
            ! -- Can't exit along alpha = 0.
            texit = huge(1d0)
            cfail = "cannot exit on face alpha=0 (no exit)"
          else
            ! -- get exit
            if (itriface .eq. 1) then
              alpexit = 1d0 - beti
            else
              alpexit = 0d0
            end if ! kluge note: seems like in this case (beta=const) this
            betexit = beti !   must be the ONLY exit; no need to check other edges??
            if (waa .ne. 0d0) then
              alplim = -v0alpstar / waa
              texit = dlog(alpexit - alplim / (alpi - alplim)) / waa
              ntdebug = -222 ! kluge debug bludebug
            else
              texit = (alpexit - alpi) / v0alpstar
              ntdebug = -333 ! kluge debug bludebug
            end if
          end if
        end if
        ! -- End constant-beta case
      else
        ! -- Beta varies along trajectory; combine outflow and soln limits on beta
        bethi = min(betouthi, betsolhi)
        betlo = max(betoutlo, betsollo)
        if (betlo .ge. bethi) then
          ! -- If bounds on bet leave no feasible interval, no exit
          texit = huge(1d0)
          cfail = "No feasible interval for beta"
        else
          ! -- Check sign of function value at beta bounds
          call get_t_alpt(bethi, thi, alphi)
          call get_t_alpt(betlo, tlo, alplo)
          if (itriface .eq. 1) then
            fax = 1d0 - betlo - alplo
            fbx = 1d0 - bethi - alphi
          else
            fax = alplo
            fbx = alphi
          end if
          if (fax * fbx .gt. 0d0) then
            ! -- Root not bracketed; no exit
            texit = huge(1d0)
            cfail = "Beta bounds betlo and bethi do not bracket a root"
          else
            if (isolv .eq. 1) then
              ! -- Use Brent's method with initial bounds on beta of betlo and bethi,
              ! -- assuming they bound the root
      call soln_brent(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail)
            else if (isolv .eq. 2) then
              ! -- Use Chandrupatla's method with initial bounds on beta of betlo and bethi,
              ! -- assuming they bound the root
      call soln_chand(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail)
            else if (isolv .eq. 3) then
              ! -- Use a test method with initial bounds on beta of betlo and bethi,
              ! -- assuming they bound the root
       call soln_test(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail)
            else if (isolv .eq. -1) then
              ! -- Use Euler integration to find exit
  call soln_euler(itriface, alpi, beti, step, vziodz, az, texit, alpexit, betexit)
            else
              write (*, '(A)') "Invalid isolv = ", isolv ! kluge
              write (69, '(A)') "Invalid isolv = ", isolv
              !!pause
              stop
            end if
          end if
        end if
        ! -- End variable-beta case
      end if
    end if
    ! -- End canonical face 1 (gamma = 0.) or 2 (alpha = 0.)
    !
  end if
  !
  if (texit .ne. huge(1d0)) then
    if (.not. lsupout) then
      write(*,'(A,A7,A,3(G14.5))') "   ", facename(itriface), ": exited at (t, alpha, beta) ", texit, alpexit, betexit
      write(69,'(A,A7,A,3(G14.5))') "   ", facename(itriface), ": exited at (t, alpha, beta) ", texit, alpexit, betexit
    end if
    if (texit .lt. 0d0) then
      write (*, '(A)') "texit is negative (unexpected)" ! kluge note: shouldn't get here
      write (69, '(A)') "texit is negative (unexpected)" ! kluge note: shouldn't get here
      stop
    end if
  else
    if (.not. lsupout) then
      write (*, '(A,A7,A,A)') "   ", facename(itriface), ": ", cfail
      write (69, '(A,A7,A,A)') "   ", facename(itriface), ": ", cfail
    end if
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
double precision function fbary1(bet)
!
  implicit double precision(a - h, o - z)
  external get_t_alpt
  !
  ! -- Return function for Brent's method applied to canonical face 1 (gamma = 0)
  !
  ! -- Evaluate gamma{t{beta}} = 1. - alpha{t{beta}} - beta
  call get_t_alpt(bet, t, alpt)
  fbary1 = 1d0 - alpt - bet
  !
  return
  !
end function
!
! ------------------------------------------------------------------------------
!
double precision function fbary2(bet)
!
  implicit double precision(a - h, o - z)
  external get_t_alpt
  !
  ! -- Return function for Brent's method applied to canonical face 2 (alpha = 0)
  !
  ! -- Evaluate alpha{t{beta}}
  call get_t_alpt(bet, t, alpt)
  fbary2 = alpt
  !
  return
  !
end function
!
! ------------------------------------------------------------------------------
!
subroutine get_t_alpt(bet, t, alp)
!
  use ternarymod, only: ca1, ca2, ca3, cb1, cb2, waa, wab, wba, wbb, icase
  implicit double precision(a - h, o - z)
  !
  ! -- Given beta evaluate t and alpha depending on case     ! kluge note: assumes cb2<>0, wbb<>0 as appropriate
  !
  zerotol = 1d-10 ! kluge
  !
  term = (bet - cb1) / cb2
  if (icase .eq. 1) then
    term = max(term, zerotol)
    t = dlog(term) / wbb
    alp = ca1 + ca2 * dexp(waa * t) + ca3 * dexp(wbb * t)
  else if (icase .eq. -1) then
    term = max(term, zerotol)
    t = dlog(term) / wbb
    alp = ca1 + (ca2 + ca3 * t) * dexp(waa * t)
  else if (icase .eq. 2) then
    term = max(term, zerotol)
    t = dlog(term) / wbb
    alp = ca1 + ca2 * t + ca3 * dexp(wbb * t)
  else if (icase .eq. 3) then
    t = term
    alp = ca1 + ca2 * t + ca3 * dexp(waa * t)
  else if (icase .eq. 4) then
    t = term
    alp = ca1 + (ca2 + ca3 * t) * t
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine get_bet_outflow_bary(vn1, vn2, betoutlo, betouthi)
!
  implicit double precision(a - h, o - z)
  !
  ! -- Find outflow interval
  !
  vndiff = vn2 - vn1
  !
  if (vn1 .lt. 0d0) then
    ! -- Outflow at vn1 corner
    betoutlo = 0d0
    if (vn2 .le. 0d0) then
      ! -- Outflow along entire edge (except possibly no-flow right at vn2 corner)
      betouthi = 1d0
    else
      ! -- Outflow along part of edge
      betouthi = -vn1 / vndiff
    end if
  else
    ! -- Outflow at vn2 corner
    betouthi = 1d0
    if (vn1 .le. 0d0) then
      ! -- Outflow along entire edge (except possibly no-flow right at vn1 corner)
      betoutlo = 0d0
    else
      ! -- Outflow along part of edge
      betoutlo = -vn1 / vndiff
    end if
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine get_bet_soln_limits(beti, betsollo, betsolhi, ibettrend)
!
  use ternarymod, only: v0bet, wbb
  implicit double precision(a - h, o - z)
  !
  ! -- Find trend of and limits on beta from beta{t} solution
  !
  if (wbb .gt. 0d0) then
    betlim = -v0bet / wbb
    if (beti .gt. betlim) then
      betsolhi = huge(1d0)
      betsollo = beti
      ibettrend = 1
    else if (beti .lt. betlim) then
      betsolhi = beti
      betsollo = -huge(1d0)
      ibettrend = -1
    else
      betsolhi = beti
      betsollo = beti
      ibettrend = 0
    end if
  else if (wbb .lt. 0d0) then
    betlim = -v0bet / wbb
    if (beti .gt. betlim) then
      betsolhi = beti
      betsollo = betlim
      ibettrend = -1
    else if (beti .lt. betlim) then
      betsolhi = betlim
      betsollo = beti
      ibettrend = 1
    else
      betsolhi = beti
      betsollo = beti
      ibettrend = 0
    end if
  else ! kluge note: use zerotol and elsewhere?
    if (v0bet .gt. 0d0) then
      betsolhi = huge(1d0)
      betsollo = beti
      ibettrend = 1
    else if (v0bet .lt. 0d0) then
      betsolhi = beti
      betsollo = -huge(1d0)
      ibettrend = -1
    else
      betsolhi = beti
      betsollo = beti
      ibettrend = 0
    end if
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine soln_brent(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail)
  !
  implicit double precision(a - h, o - z)
  character cfail * 60
  external fbary1, fbary2
  !
  ! -- Use Brent's method with initial bounds on beta of betlo and bethi,
  ! -- assuming they bracket the root
!!!  tol = 1d-7               ! kluge
  itmax = 50 ! kluge
  itact = itmax + 1 ! kluge
  blo = betlo
  bhi = bethi
  if (itriface .eq. 1) then
    betexit = zeroin(blo, bhi, fbary1, tol)
  else
    betexit = zeroin(blo, bhi, fbary2, tol)
  end if
  call get_t_alpt(betexit, texit, alpexit)
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine soln_chand(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail)
  !
  implicit double precision(a - h, o - z)
  character cfail * 60
  external fbary1, fbary2
  !
  ! -- Use Chandrupatla's method with initial bounds on beta of betlo and bethi,
  ! -- assuming they bracket the root
!!!  tol = 1d-7               ! kluge
  itmax = 50 ! kluge
  itact = itmax + 1 ! kluge
  blo = betlo
  bhi = bethi
  if (itriface .eq. 1) then
    betexit = zeroch(blo, bhi, fbary1, tol)
  else
    betexit = zeroch(blo, bhi, fbary2, tol)
  end if
  call get_t_alpt(betexit, texit, alpexit)
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine soln_test(itriface, betlo, bethi, tol, texit, alpexit, betexit, cfail) ! kluge test
  !
  implicit double precision(a - h, o - z)
  character cfail * 60
  external fbary1, fbary2
  !
  ! -- Use a test method with initial bounds on beta of betlo and bethi,
  ! -- assuming they bracket the root
!!!  tol = 1d-7               ! kluge
  itmax = 50 ! kluge
  itact = itmax + 1 ! kluge
  blo = betlo
  bhi = bethi
  if (itriface .eq. 1) then
    betexit = zerotest(blo, bhi, fbary1, tol)
  else
    betexit = zerotest(blo, bhi, fbary2, tol)
  end if
  call get_t_alpt(betexit, texit, alpexit)
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
subroutine soln_euler(itriface,alpi,beti,step,vziodz,az,texit,alpexit,betexit)
  !
  implicit double precision(a - h, o - z)
  common / debug / ntdebug ! kluge debug
  !
  ! -- Use Euler integration to find exit
  t = 0d0
  alp = alpi
  bet = beti
  if (itriface .eq. 1) gam = 1d0 - alpi - beti
  do nt = 1, 1000000000 ! kluge hardwired
    ! -- Save current time, alpha, and beta
    told = t
    alpold = alp
    betold = bet
    ! -- Step forward in time
!!!    t = dble(nt)*step
    call step_euler(nt, step, vziodz, az, alpi, beti, t, alp, bet)
!!!    if (nt.eq.0) then
!!!      znum = zi
!!!    else
!!!      vz = vzbot + az*(znum - zbot)
!!!      znum = znum + vz*delt           ! kluge note: can be smart about checking z
!!!    end if
    if (itriface .eq. 1) then
      ! -- If gamma has crossed zero, interpolate linearly
      ! -- to find crossing (exit) point
      gamold = gam
      gam = 1d0 - alp - bet
      if (gam .lt. 0d0) then
        wt = gamold / (gamold - gam)
        omwt = 1d0 - wt
        texit = omwt * told + wt * t
        alpexit = omwt * alpold + wt * alp
        betexit = omwt * betold + wt * bet
        exit
      end if
    else
      ! -- If alpha has crossed zero, interpolate linearly
      ! -- to find crossing (exit) point
      if (alp .lt. 0d0) then
        wt = alpold / (alpold - alp)
        omwt = 1d0 - wt
        texit = omwt * told + wt * t
        alpexit = omwt * alpold + wt * alp
        betexit = omwt * betold + wt * bet
        exit
      end if
    end if
    ! -- End time step loop
  end do
  ntdebug = nt ! kluge debug
  if (nt .gt. 1000000000) then ! kluge hardwired
    ! -- Exit not found after max number of time steps
    write (*, '(A)') "Didn't find exit in soln_euler" ! kluge note: shouldn't get here
    write (69, '(A)') "Didn't find exit in soln_euler" ! kluge note: shouldn't get here
    stop
  end if
  !
  return
  !
end subroutine
!
! ------------------------------------------------------------------------------
!
double precision function alpfun(t) ! kluge: not needed???
!
  use ternarymod, only : ca1,ca2,ca3,cb1,cb2,waa,wab,wba,wbb,icase,alpexitcopy
  implicit double precision(a - h, o - z)
  !
  ! -- Given t evaluate alpha depending on case     ! kluge note: assumes cb2<>0, wbb<>0 as appropriate
  !
  if (icase .eq. 1) then
    alpfun = ca1 + ca2 * dexp(waa * t) + ca3 * dexp(wbb * t)
  else if (icase .eq. -1) then
    alpfun = ca1 + (ca2 + ca3 * t) * dexp(waa * t)
  else if (icase .eq. 2) then
    alpfun = ca1 + ca2 * t + ca3 * dexp(wbb * t)
  else if (icase .eq. 3) then
    alpfun = ca1 + ca2 * t + ca3 * dexp(waa * t)
  else if (icase .eq. 4) then
    alpfun = ca1 + (ca2 + ca3 * t) * t
  end if
  alpfun = alpfun - alpexitcopy
  !
  return
  !
end function
!
! ------------------------------------------------------------------------------
!
double precision function zeroch(x0, x1, f, epsa)
!
  implicit double precision(a - h, o - z)
  !
  ! -- A zero of the function f{x} is computed in the interval (x0, x1)
  ! -- given tolerance epsa using Chandrupatla's method. FORTRAN code based
  ! -- generally on pseudocode in Scherer, POJ (2013) "Computational Physics:
  ! -- Simulation of Classical and Quantum Systems," 2nd ed., Springer, New York.
  !
!!  epsm = d1mach(4)
  epsm = epsilon(x0)
  !
  b = x0
  a = x1
  c = x1
  aminusb = a - b
  fb = f(b)
  fa = f(a)
  fc = f(c)
  t = 5d-1
  !
  do while (.true.)
!!    xt = a + t*(b - a)
    xt = a - t * aminusb
    ft = f(xt)
    if (sign(ft, fa) == ft) then
      c = a
      fc = fa
      a = xt
      fa = ft
    else
      c = b
      b = a
      a = xt
      fc = fb
      fb = fa
      fa = ft
    end if
    aminusb = a - b
    cminusb = c - b
    faminusfb = fa - fb
    fcminusfb = fc - fb
    xm = a
    fm = fa
    if (dabs(fb) < dabs(fa)) then
      xm = b
      fm = fb
    end if
    tol = 2d0 * epsm * dabs(xm) + epsa
!!    tl = tol/dabs(b - c)
    tl = tol / dabs(cminusb)
    if ((tl > 5d-1) .or. (fm == 0d0)) then
      zeroch = xm
      return
    end if
!!    xi = (a - b)/(c - b)
    xi = aminusb / cminusb
!!    phi = (fa - fb)/(fc - fb)
    phi = faminusfb / fcminusfb
    philo = 1d0 - dsqrt(1d0 - xi)
    phihi = dsqrt(xi)
    if ((phi > philo) .and. (phi < phihi)) then
!!!!      rab = fa/(fb - fa)
!!      rab = -fa/faminusfb
!!!!      rcb = fc/(fb - fc)
!!      rcb = -fc/fcminusfb
!!      rac = fa/(fc - fa)
!!!!      rbc = fb/(fc - fb)
!!      rbc = fb/fcminusfb
!!!!      t = rab*rcb + rac*rbc*(c - a)/(b - a)
!!      t = rab*rcb - rac*rbc*(c - a)/aminusb
      racb = fa / fcminusfb
      rcab = fc / faminusfb
      rbca = fb / (fc - fa)
      t = racb * (rcab - rbca * (c - a) / aminusb)
      if (t < tl) then
        t = tl
      else
        tlc = 1d0 - tl
        if (t > tlc) then
          t = tlc
        end if
      end if
    else
      t = 5d-1
    end if
!!    if (t < tl) t = tl
!!    if (t > 1d0 - tl) t = 1d0 - tl
  end do
  !
  return
  !
end function
!
! ------------------------------------------------------------------------------
!
double precision function zerotest(x0, x1, f, epsa) ! kluge test
!
  implicit double precision(a - h, o - z)
  logical retainedxa, retainedxb
  !
  ! -- A zero of the function f{x} is computed in the interval (x0, x1)
  ! -- given tolerance epsa using a test method.
  !
  epsm = epsilon(x0)
  !
  f0 = f(x0)
  if (f0 .eq. 0d0) then
    zerotest = x0
    return
  else if (f0 .lt. 0d0) then
    ya = x0
    yb = x1
    xa = f0
    xb = f(yb)
  else
    ya = x1
    yb = x0
    xa = f(ya)
    xb = f0
  end if
  ema = 1d0
  emb = 1d0
  retainedxa = .false.
  retainedxb = .false.
  !
  do while (.true.)
!!    yl = ya - xa*(yb - ya)/(xb - xa)
    yl = (ya * xb * emb - yb * xa * ema) / (xb * emb - xa * ema)
    tol = 4d0 * epsm * dabs(yl) + epsa
    if (dabs(yb - ya) .le. tol) then
      zerotest = yl
      return
    else
      xl = f(yl)
      if (xl .eq. 0d0) then
        zerotest = yl
        return
      else if (xl .gt. 0d0) then
        if (retainedxa) then
!!          ema = 1d0 - xl/xb
!!          if (ema <= 0d0) ema = 5d-1
          ema = 5d-1 ! kluge illinois
        else
          ema = 1d0
        end if
        emb = 1d0
        yb = yl
        xb = xl
        retainedxa = .true.
        retainedxb = .false.
      else
        if (retainedxb) then
!!          emb = 1d0 - xl/xa
!!          if (emb <= 0d0) emb = 5d-1
          emb = 5d-1 ! kluge illinois
        else
          emb = 1d0
        end if
        ema = 1d0
        ya = yl
        xa = xl
        retainedxa = .false.
        retainedxb = .true.
      end if
    end if
  end do
  !
  return
  !
end function
