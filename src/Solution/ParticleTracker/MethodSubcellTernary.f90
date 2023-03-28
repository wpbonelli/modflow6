module MethodSubcellTernaryModule

  use KindModule, only: DP, I4B
  use MethodModule
  use SubcellTriModule
  use ParticleModule
  use ternarymod ! kluge
  implicit none

  private
  public :: MethodSubcellTernaryType
  public :: create_methodSubcellTernary

  ! -- Extend MethodType to the ternary subcell-method type (MethodSubcellTernaryType)
  type, extends(MethodType) :: MethodSubcellTernaryType
    private
    type(SubcellTriType), pointer, public :: subcellTri => null() ! tracking domain for the method
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mST ! applies the ternary subcell method (tracks particle)
  end type MethodSubcellTernaryType

contains

  subroutine create_methodSubcellTernary(methodSubcellTernary)
! ******************************************************************************
! create_methodSubcellTernary -- Create a new ternary subcell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(MethodSubcellTernaryType), pointer :: methodSubcellTernary
! ------------------------------------------------------------------------------
    !
    allocate (methodSubcellTernary)
    !
    ! -- This method does not delegate tracking to a submethod; it performs
    ! -- actual particle-tracking calculations
    methodSubcellTernary%delegatesTracking = .FALSE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_subcellTri(methodSubcellTernary%subcellTri)
!!    methodSubcellTernary%trackingDomain => methodSubcellTernary%subcellTri
   methodSubcellTernary%trackingDomainType => methodSubcellTernary%subcellTri%type
    !
    return
    !
  end subroutine create_methodSubcellTernary

  subroutine destroy(this)
! ******************************************************************************
! destroy -- Destructor for a ternary subcell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(MethodSubcellTernaryType), intent(inout) :: this
! ------------------------------------------------------------------------------
    !
    deallocate (this%trackingDomainType)
    !
    return
    !
  end subroutine destroy

  subroutine init(this, subcellTri)
! ******************************************************************************
! init -- Initialize a ternary subcell-method object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(MethodSubcellTernaryType), intent(inout) :: this
    type(SubcellTriType), pointer :: subcellTri
! ------------------------------------------------------------------------------
    !
    this%subcellTri => subcellTri
    !
    return
    !
  end subroutine init

  subroutine apply_mST(this, particle, tmax)
! ******************************************************************************
! apply_mST -- Apply the ternary method to a polygonal cell
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! dummy
    class(MethodSubcellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
!!!    doubleprecision :: initialTime,maximumTime,t   ! kluge not in arg list yet
    ! local
! ------------------------------------------------------------------------------
    !
!!!    initialTime = 0d0         ! kluge test
!!!    maximumTime = 9d99        ! kluge test
!!!    call track_sub(this%subcellTri,particle,initialTime,maximumTime,t)
    call track_sub(this%subcellTri, particle, tmax)
    !
    return
    !
  end subroutine apply_mST

!!!  subroutine track_sub(subcellTri,particle,initialTime,maximumTime,t)   ! kluge note: rename???
  subroutine track_sub(subcellTri, particle, tmax) ! kluge note: rename???
! ******************************************************************************
! track_sub -- Track a particle across a triangular subcell using the ternary
!              method
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! dummy
    class(SubcellTriType), intent(in) :: subcellTri
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
!!!    doubleprecision :: initialTime,maximumTime,t
    ! local
    integer :: exitFace
    logical :: lbary ! kluge
    double precision :: x0, y0, x1, y1, x2, y2, v0x, v0y, v1x, v1y, v2x, v2y
    double precision :: xi, yi, zi, zirel, ztop, zbot, dz
    double precision :: rxx, rxy, ryx, ryy, sxx, sxy, syy
    double precision :: alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti
    double precision :: vzbot, vztop, vzi, vziodz, az, dtexitz
    double precision :: alp, bet, dt, t, dtexitxy, texit, xc, yc, x, y, z
    integer :: izstatus, itopbotexit
    integer :: ntmax, nsave, isolv, itrifaceenter, itrifaceexit
    double precision :: diff,rdiff,tol,step,dtexit,alpexit,betexit,xcexit,ycexit
    double precision :: xexit, yexit, zexit
    integer :: ntdebug ! kluge
    ! ------------------------------------------------------------------------------
    !
    lbary = .true. ! kluge
    ntmax = 10000
    nsave = 1 ! needed???
!!    isolv = 1
    isolv = 1
!!    tol = 1d-7     ! adjustable
    tol = 1d-7
    step = 1e-3 ! needed only for euler
    !
    if (.not. allocated(ivert_polygon)) allocate (ivert_polygon(1, 4)) ! kluge, needed???
    ivert_polygon(1, 1) = 1
    ivert_polygon(1, 2) = 2
    ivert_polygon(1, 3) = 3
    ivert_polygon(1, 4) = 4
    !
    ! -- Set some local variables for convenience
    xi = particle%x
    yi = particle%y
    zi = particle%z
    x0 = subcellTri%x0
    y0 = subcellTri%y0
    x1 = subcellTri%x1
    y1 = subcellTri%y1
    x2 = subcellTri%x2
    y2 = subcellTri%y2
    v0x = subcellTri%v0x
    v0y = subcellTri%v0y
    v1x = subcellTri%v1x
    v1y = subcellTri%v1y
    v2x = subcellTri%v2x
    v2y = subcellTri%v2y
    zbot = subcellTri%zbot
    ztop = subcellTri%ztop
    dz = subcellTri%dz
    vzbot = subcellTri%vzbot
    vztop = subcellTri%vztop
    !
    ! -- Translate and rotate coordinates to "canonical" configuration
    call canonical(x0, y0, x1, y1, x2, y2, v0x, v0y, v1x, v1y, v2x, v2y, xi, yi, &
   rxx, rxy, ryx, ryy, sxx, sxy, syy, lbary, alp0, bet0, alp1, bet1, alp2, bet2, &
                   alpi, beti)
    !
    ! -- Do calculations related to analytical z solution, which can be done
    ! -- after traverse_triangle call if results not needed for adaptive time
    ! -- stepping during triangle (subcell) traversal   ! kluge note: actually, can probably do z calculation just once for each cell
    zirel = (zi - zbot) / dz
  call pr_CalculateDT_kluge(vzbot, vztop, dz, zirel, vzi, az, dtexitz, izstatus, &
                              itopbotexit)
    vziodz = vzi / dz
    !
    ! -- Traverse triangular subcell
    ntdebug = -999 ! kluge debug bludebug
    itrifaceenter = particle%iTrackingDomainBoundary(3) - 1
    if (itrifaceenter .eq. -1) itrifaceenter = 999
    call traverse_triangle(ntmax, nsave, diff, rdiff, isolv, tol, step, & ! kluge note: can probably avoid calculating alpexit
dtexitxy, alpexit, betexit, itrifaceenter, itrifaceexit, rxx, rxy, ryx, ryy, & !   here in many cases and wait to calculate it later,
sxx, sxy, syy, lbary, alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti, vziodz, az) !   once the final trajectory time is known
    !
    ! -- Check for no exit face
    if ((itopbotexit .eq. 0) .and. (itrifaceexit .eq. 0)) then
!!      exitFace = 0
!!      particle%iTrackingDomainBoundary(3) = exitFace
!!      particle%istatus = 5
!!      return
      print *, "======================================"
      print *, "Subcell with no exit face" ! kluge
      print *, "Particle ", particle%ipart
      print *, "Cell ", particle%iTrackingDomain(2)
      print *, "======================================"
      !!pause
      !!stop
    end if
    !
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
    !
    ! -- Compute exit time
    texit = particle%ttrack + dtexit
    !
    if (texit .gt. tmax) then
      ! -- The computed exit time is greater than the maximum time, so set
      ! -- final time for particle trajectory equal to maximum time.
      t = tmax
      dt = t - particle%ttrack
      exitFace = 0
      particle%istatus = 1
    else
      ! -- The computed exit time is less than or equal to the maximum time,
      ! -- so set final time for particle trajectory equal to exit time.
      t = texit
      dt = dtexit
    end if
    !
    ! -- Calculate final particle location
    call step_analytical(dt, alp, bet) ! kluge note: need to evaluate both alpha and beta here only for exitFace=0, otherwise just one or the other
    if (exitFace .eq. 1) then
      bet = 0d0
    else if (exitFace .eq. 2) then
      alp = 1d0 - bet
    else if (exitFace .eq. 3) then
      alp = 0d0
    end if
    xc = alp
    yc = bet
    if (lbary) call skew(-1, sxx, sxy, syy, xc, yc)
    call rotate(rxx, ryx, rxy, ryy, xc, yc, x, y)
    x = x + x0
    y = y + y0
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
    !
    ! -- Set final particle location in local (unscaled) subcell coordinates,
    ! -- final time for particle trajectory, and exit face
    particle%x = x
    particle%y = y
    particle%z = z
    particle%ttrack = t
    particle%iTrackingDomainBoundary(3) = exitFace
!!    write(*,*) itopbotexit, itrifaceexit, exitFace        ! kluge debug
!!    write(*,'(4G)') x, y, z, dt                           ! kluge debug
!!    write(69,*) itopbotexit, itrifaceexit, exitFace       ! kluge debug
!!    write(69,'(4G)') x, y, z, dt                          ! kluge debug
    !
    return
    !
  end subroutine track_sub

!-----------------------------------------------------------------------------
subroutine pr_CalculateDT_kluge(v1, v2, dx, xL, v, dvdx, dt, status, itopbotexit) ! kluge
    implicit none
    doubleprecision, intent(in) :: v1, v2, dx, xL
    doubleprecision, intent(inout) :: v, dvdx, dt
    doubleprecision :: v2a, v1a, dv, dva, vv, vvv, zro, zrom, x, tol
    doubleprecision :: vr1, vr2, vr, v1v2
    integer :: status, itopbotexit
    logical :: noOutflow
! ------------------------------------------------------------------------------
    !
    ! -- This subroutine consists partly or entirely of code written by
    ! -- David W. Pollock of the USGS for MODPATH 7. The authors of the present
    ! -- code are responsible for its appropriate application in this context
    ! -- and for any modifications or errors.
    !
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

    ! If there is a flow divide, this check finds out what side of the divide the particle
    ! is on and sets the value of vr appropriately to reflect that location.
    vr1 = v1 / v
    vr2 = v2 / v
    vr = vr1
    itopbotexit = -1
    if (vr .le. 0d0) then
      vr = vr2
      itopbotexit = -2
    end if

    ! Check to see if the velocity is in the same direction throughout the cell (i.e. no flow divide).
    ! Check to see if the product v1*v2 > 0 then the velocity is in the same direction throughout
    ! the cell (i.e. no flow divide). If so, set the value of vr to reflect the appropriate direction.
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

  end subroutine pr_CalculateDT_kluge

end module MethodSubcellTernaryModule
