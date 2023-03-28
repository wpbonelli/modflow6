module MethodSubcellTernaryModule
  use KindModule, only: DP, I4B
  use ErrorUtilModule, only: pstop
  use GeomUtilModule, only: skew
  use MethodModule, only: MethodType
  use CellModule, only: CellType
  use SubcellModule, only: SubcellType
  use SubcellTriModule, only: SubcellTriType, create_subcell_tri
  use ParticleModule, only: ParticleType, get_particle_id
  use TrackModule, only: TrackControlType
  use TernarySolveTrack, only: traverse_triangle, step_analytical, canonical
  use PrtFmiModule, only: PrtFmiType
  use BaseDisModule, only: DisBaseType
  use GwfDisvModule, only: GwfDisvType
  implicit none

  private
  public :: MethodSubcellTernaryType
  public :: create_method_subcell_ternary

  !> @brief Ternary triangular subcell tracking method
  type, extends(MethodType) :: MethodSubcellTernaryType
  contains
    procedure, public :: apply => apply_mst
    procedure, public :: destroy
    procedure, private :: track_subcell
  end type MethodSubcellTernaryType

contains

  !> @brief Create a new ternary subcell-method object
  subroutine create_method_subcell_ternary(method)
    ! -- dummy
    type(MethodSubcellTernaryType), pointer :: method
    ! -- local
    type(SubcellTriType), pointer :: subcell

    allocate (method)
    call create_subcell_tri(subcell)
    method%subcell => subcell
    method%type => method%subcell%type
    method%delegates = .false.
  end subroutine create_method_subcell_ternary

  !> @brief Destructor for a ternary subcell-method object
  subroutine destroy(this)
    class(MethodSubcellTernaryType), intent(inout) :: this
    deallocate (this%type)
  end subroutine destroy

  !> @brief Apply the ternary subcell method
  subroutine apply_mst(this, particle, tmax)
    class(MethodSubcellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax

    select type (subcell => this%subcell)
    type is (SubcellTriType)
      call this%track_subcell(subcell, particle, tmax)
    end select
  end subroutine apply_mst

  !> @brief Track a particle across a triangular subcell using the ternary method
  subroutine track_subcell(this, subcell, particle, tmax)
    ! modules
    use TdisModule, only: kper, kstp
    ! dummy
    class(MethodSubcellTernaryType), intent(inout) :: this
    class(SubcellTriType), intent(in) :: subcell
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    integer :: exitFace
    logical :: lbary ! kluge
    double precision :: x0, y0, x1, y1, x2, y2
    double precision :: v0x, v0y, v1x, v1y, v2x, v2y
    double precision :: xi, yi, zi, zirel, ztop, zbot, dz
    double precision :: rxx, rxy, ryx, ryy, sxx, sxy, syy
    double precision :: rot(2, 2), res(2), loc(2)
    double precision :: alp, bet, alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti
    double precision :: vzbot, vztop, vzi, vziodz, az, dtexitz
    double precision :: dt, t, dtexitxy, texit, x, y, z
    integer :: izstatus, itopbotexit
    integer :: ntmax, nsave, isolv, itrifaceenter, itrifaceexit
    double precision :: diff, rdiff, tol, step, dtexit, alpexit, betexit
    integer :: ntdebug ! kluge
    integer :: reason

    lbary = .true. ! kluge
    ntmax = 10000
    nsave = 1 ! needed???
    isolv = 1
    tol = 1d-7
    step = 1e-3 ! needed only for euler
    reason = -1

    ! -- Set some local variables for convenience
    xi = particle%x
    yi = particle%y
    zi = particle%z
    x0 = subcell%x0
    y0 = subcell%y0
    x1 = subcell%x1
    y1 = subcell%y1
    x2 = subcell%x2
    y2 = subcell%y2
    v0x = subcell%v0x
    v0y = subcell%v0y
    v1x = subcell%v1x
    v1y = subcell%v1y
    v2x = subcell%v2x
    v2y = subcell%v2y
    zbot = subcell%zbot
    ztop = subcell%ztop
    dz = subcell%dz
    vzbot = subcell%vzbot
    vztop = subcell%vztop

    ! -- Translate and rotate coordinates to "canonical" configuration
    call canonical(x0, y0, x1, y1, x2, y2, &
                   v0x, v0y, v1x, v1y, v2x, v2y, &
                   xi, yi, &
                   rxx, rxy, ryx, ryy, &
                   sxx, sxy, syy, &
                   alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, &
                   lbary)

    ! -- Do calculations related to analytical z solution, which can be done
    ! -- after traverse_triangle call if results not needed for adaptive time
    ! -- stepping during triangle (subcell) traversal
    ! kluge note: actually, can probably do z calculation just once for each cell
    zirel = (zi - zbot) / dz
    call calculate_dt(vzbot, vztop, dz, zirel, vzi, &
                      az, dtexitz, izstatus, &
                      itopbotexit)
    vziodz = vzi / dz

    ! -- Traverse triangular subcell
    ntdebug = -999 ! kluge debug bludebug
    itrifaceenter = particle%iboundary(3) - 1
    if (itrifaceenter .eq. -1) itrifaceenter = 999

    ! kluge note: can probably avoid calculating alpexit
    ! here in many cases and wait to calculate it later,
    ! once the final trajectory time is known
    call traverse_triangle(ntmax, nsave, diff, rdiff, &
                           isolv, tol, step, &
                           dtexitxy, alpexit, betexit, &
                           itrifaceenter, itrifaceexit, &
                           rxx, rxy, ryx, ryy, &
                           sxx, sxy, syy, &
                           alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, &
                           vziodz, az, lbary)

    ! -- Check for no exit face
    if ((itopbotexit .eq. 0) .and. (itrifaceexit .eq. 0)) then
      ! exitFace = 0
      ! particle%iboundary(3) = exitFace
      ! particle%istatus = 5
      ! return

      ! contact the developer situation (for now? always?)
      print *, "Subcell with no exit face: particle", get_particle_id(particle), &
        "cell", particle%idomain(2)
      call pstop(1)
    end if

    ! -- Determine (earliest) exit face and corresponding travel time to exit
    if (itopbotexit .eq. 0) then
      ! -- Exits through triangle face first
      exitFace = itrifaceexit
      dtexit = dtexitxy
    else if (itrifaceexit .eq. 0) then
      ! -- Exits through top/bottom first
      exitFace = 45
      dtexit = dtexitz
    else if (dtexitz .lt. dtexitxy) then
      ! -- Exits through top/bottom first
      exitFace = 45
      dtexit = dtexitz
    else
      ! -- Exits through triangle face first
      exitFace = itrifaceexit
      dtexit = dtexitxy
    end if
    if (exitFace .eq. 45) then
      if (itopbotexit .eq. -1) then
        exitFace = 4
      else
        exitFace = 5
      end if
    end if

    ! -- Compute exit time
    texit = particle%ttrack + dtexit

    if (texit .gt. tmax) then
      ! -- The computed exit time is greater than the maximum time, so set
      ! -- final time for particle trajectory equal to maximum time.
      t = tmax
      dt = t - particle%ttrack
      exitFace = 0
      particle%istatus = 1
      particle%advancing = .false.
      reason = 2 ! timestep end
    else
      ! -- The computed exit time is less than or equal to the maximum time,
      ! -- so set final time for particle trajectory equal to exit time.
      t = texit
      dt = dtexit
      reason = 1 ! cell transition
    end if

    ! -- Calculate final particle location
    ! -- kluge note: need to evaluate both alpha and beta here only
    ! -- for exitFace=0, otherwise just one or the other
    call step_analytical(dt, alp, bet)
    if (exitFace .eq. 1) then
      bet = 0d0
    else if (exitFace .eq. 2) then
      alp = 1d0 - bet
    else if (exitFace .eq. 3) then
      alp = 0d0
    end if
    loc = (/alp, bet/)
    if (lbary) loc = skew(loc, (/sxx, sxy, syy/), invert=.true.)
    rot = reshape((/rxx, rxy, ryx, ryy/), shape(rot))
    res = matmul(rot, loc) ! rotate vector
    x = res(1) + x0
    y = res(2) + y0
    if (exitFace .eq. 4) then
      z = zbot
    else if (exitFace .eq. 5) then
      z = ztop
    else
      if (izstatus .eq. 2) then ! kluge note: make this into a function
        ! -- vz uniformly zero
        z = zi
      else if (izstatus .eq. 1) then
        ! -- vz uniform, nonzero
        z = zi + vzi * dt
      else
        ! -- vz nonuniform
        z = zbot + (vzi * dexp(az * dt) - vzbot) / az
      end if
    end if

    ! -- Set final particle location in local (unscaled) subcell coordinates,
    ! -- final time for particle trajectory, and exit face
    particle%x = x
    particle%y = y
    particle%z = z
    particle%ttrack = t
    particle%iboundary(3) = exitFace

    ! -- Save particle track record
    if (reason > -1) &
      call this%trackctl%save(particle, kper=kper, &
                              kstp=kstp, reason=reason) ! reason=2: timestep
  end subroutine track_subcell

  !> @brief Do calculations related to analytical z solution
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !<
  subroutine calculate_dt(v1, v2, dx, xL, v, dvdx, &
                          dt, status, itopbotexit)
    doubleprecision, intent(in) :: v1, v2, dx, xL
    doubleprecision, intent(inout) :: v, dvdx, dt
    doubleprecision :: v2a, v1a, dv, dva, vv, vvv, zro, zrom, x, tol
    doubleprecision :: vr1, vr2, vr, v1v2
    integer :: status, itopbotexit
    logical :: noOutflow

    ! Initialize variables
    status = -1
    dt = 1.0d+20
    v2a = v2
    if (v2a .lt. 0d0) v2a = -v2a
    v1a = v1
    if (v1a .lt. 0d0) v1a = -v1a
    dv = v2 - v1
    dva = dv
    if (dva .lt. 0d0) dva = -dva

    ! Check for a uniform zero velocity in this direction.
    ! If so, set status = 2 and return (dt = 1.0d+20).
    tol = 1.0d-15
    if ((v2a .lt. tol) .and. (v1a .lt. tol)) then
      v = 0d0
      dvdx = 0d0
      status = 2
      itopbotexit = 0
      return
    end if

    ! Check for uniform non-zero velocity in this direction.
    ! If so, set compute dt using the constant velocity,
    ! set status = 1 and return.
    vv = v1a
    if (v2a .gt. vv) vv = v2a
    vvv = dva / vv
    if (vvv .lt. 1.0d-4) then
      zro = tol
      zrom = -zro
      v = v1
      x = xL * dx
      if (v1 .gt. zro) then
        dt = (dx - x) / v1
        itopbotexit = -2
      end if
      if (v1 .lt. zrom) then
        dt = -x / v1
        itopbotexit = -1
      end if
      dvdx = 0d0
      status = 1
      return
    end if

    ! Velocity has a linear variation.
    ! Compute velocity corresponding to particle position
    dvdx = dv / dx
    v = (1.0d0 - xL) * v1 + xL * v2

    ! If flow is into the cell from both sides there is no outflow.
    ! In that case, set status = 3 and return
    noOutflow = .true.
    if (v1 .lt. 0d0) noOutflow = .false.
    if (v2 .gt. 0d0) noOutflow = .false.
    if (noOutflow) then
      status = 3
      itopbotexit = 0
      return
    end if

    ! If there is a divide in the cell for this flow direction, check to see if the
    ! particle is located exactly on the divide. If it is, move it very slightly to
    ! get it off the divide. This avoids possible numerical problems related to
    ! stagnation points.
    if ((v1 .le. 0d0) .and. (v2 .ge. 0d0)) then
      if (abs(v) .le. 0d0) then
        v = 1.0d-20
        if (v2 .le. 0d0) v = -v
      end if
    end if

    ! If there is a flow divide, find out what side of the divide the particle
    ! is on and set the value of vr appropriately to reflect that location.
    vr1 = v1 / v
    vr2 = v2 / v
    vr = vr1
    itopbotexit = -1
    if (vr .le. 0d0) then
      vr = vr2
      itopbotexit = -2
    end if

    ! Check if velocity is in the same direction throughout cell (i.e. no flow divide).
    ! Check if product v1*v2 > 0 then the velocity is in the same direction throughout
    ! the cell (i.e. no flow divide). If so, set vr to reflect appropriate direction.
    v1v2 = v1 * v2
    if (v1v2 .gt. 0d0) then
      if (v .gt. 0d0) then
        vr = vr2
        itopbotexit = -2
      end if
      if (v .lt. 0d0) then
        vr = vr1
        itopbotexit = -1
      end if
    end if

    ! Compute travel time to exit face. Return with status = 0
    dt = log(vr) / dvdx
    status = 0
  end subroutine calculate_dt

end module MethodSubcellTernaryModule
