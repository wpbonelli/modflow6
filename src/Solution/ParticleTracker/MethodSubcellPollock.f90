module MethodSubcellPollockModule
  use KindModule, only: DP, I4B
  use MethodModule
  use SubcellRectModule
  use ParticleModule
  use TrackDataModule, only: TrackDataType
  implicit none
  private
  public :: MethodSubcellPollockType
  public :: create_methodSubcellPollock
  public :: CalculateDT

  ! -- Extend MethodType to the Pollock's subcell-method type (MethodSubcellPollockType)
  type, extends(MethodType) :: MethodSubcellPollockType
    private
    type(SubcellRectType), pointer, public :: subcellRect => null() ! tracking domain for the method
    double precision, allocatable, public :: qextl1(:), qextl2(:), qintl(:) ! external and internal subcell flows for the cell
  contains
    procedure, public :: destroy ! destructor for the method
    procedure, public :: init ! initializes the method
    procedure, public :: apply => apply_mSP ! applies Pollock's subcell method (tracks particle)
  end type MethodSubcellPollockType

contains

  !> @brief Create a new Pollock's subcell-method object
  subroutine create_methodSubcellPollock(methodSubcellPollock)
    ! -- dummy
    type(MethodSubcellPollockType), pointer :: methodSubcellPollock
    !
    allocate (methodSubcellPollock)
    !
    ! -- This method does not delegate tracking to a submethod; it performs
    ! -- actual particle-tracking calculations
    methodSubcellPollock%delegatesTracking = .FALSE.
    !
    ! -- Create tracking domain for this method and set trackingDomain pointer
    call create_subcellRect(methodSubcellPollock%subcellRect)
    ! methodSubcellPollock%trackingDomain => methodSubcellPollock%subcellRect
    methodSubcellPollock%trackingDomainType => &
      methodSubcellPollock%subcellRect%type
    !
    return
  end subroutine create_methodSubcellPollock

  !> @brief Destructor for a Pollock's subcell-method object
  subroutine destroy(this)
    ! -- dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    !
    deallocate (this%trackingDomainType)
    return
  end subroutine destroy

  !> @brief Initialize a Pollock's subcell-method object
  subroutine init(this, subcellRect, trackdata)
    ! -- dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    type(SubcellRectType), pointer :: subcellRect
    type(TrackDataType), pointer :: trackdata
    !
    ! -- Set pointer to subcell definition
    this%subcellRect => subcellRect
    !
    ! -- Set pointer to particle track data
    this%trackdata => trackdata
    !
    return
  end subroutine init

  !> @brief Apply Pollock's method to a rectangular subcell
  subroutine apply_mSP(this, particle, tmax)
    use UtilMiscModule
    ! -- dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! doubleprecision :: initialTime,maximumTime,t   ! kluge not in arg list yet
    ! -- local
    double precision :: xOrigin, yOrigin, zOrigin, sinrot, cosrot
    !
    ! initialTime = 0d0         ! kluge test
    ! maximumTime = 9d99        ! kluge test
    ! !
    ! -- Transform particle location into local subcell coordinates
    xOrigin = this%subcellRect%xOrigin
    yOrigin = this%subcellRect%yOrigin
    zOrigin = this%subcellRect%zOrigin
    sinrot = this%subcellRect%sinrot
    cosrot = this%subcellRect%cosrot
    ! -- sinrot and cosrot should be 0. and 1., respectively, i.e., there
    ! -- is no rotation. Also no z translation; only x and y translations.
    call particle%transf_coords(xOrigin, yOrigin)
    !
    ! -- Track the particle across the subcell
    ! call track_sub(this%subcellRect,particle,initialTime,maximumTime,t)
    call track_sub(this%subcellRect, particle, tmax)
    !
    ! -- Transform particle location back to local cell coordinates
    call particle%transf_coords(xOrigin, yOrigin, invert_opt=.true.)
    !
    return
    !
  end subroutine apply_mSP

  !> @brief Track a particle across a rectangular subcell using Pollock's method
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !!
  !<
  subroutine track_sub(subcellRect, particle, tmax) ! kluge note: rename???
    ! dummy
    class(SubcellRectType), intent(in) :: subcellRect
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    doubleprecision :: vx1, vx2, vy1, vy2, vz1, vz2
    doubleprecision :: vx, dvxdx, vy, dvydy, vz, dvzdz
    doubleprecision :: dtexitx, dtexity, dtexitz, dtexit, texit, dt, t
    doubleprecision :: x, y, z
    integer :: statusVX, statusVY, statusVZ
    doubleprecision :: initialX, initialY, initialZ
    integer :: exitFace
    !
    ! -- Initial particle location in scaled subcell coordinates
    initialX = particle%x / subcellRect%dx
    initialY = particle%y / subcellRect%dy
    initialZ = particle%z / subcellRect%dz
    !
    ! -- Make local copies of face velocities for convenience
    vx1 = subcellRect%vx1
    vx2 = subcellRect%vx2
    vy1 = subcellRect%vy1
    vy2 = subcellRect%vy2
    vz1 = subcellRect%vz1
    vz2 = subcellRect%vz2
    !
    ! -- Compute time of travel to each possible exit face
    statusVX = CalculateDT(vx1, vx2, subcellRect%dx, initialX, vx, dvxdx, &
                           dtexitx)
    statusVY = CalculateDT(vy1, vy2, subcellRect%dy, initialY, vy, dvydy, &
                           dtexity)
    statusVZ = CalculateDT(vz1, vz2, subcellRect%dz, initialZ, vz, dvzdz, &
                           dtexitz)
    !
    ! -- Check for no exit face
    if ((statusVX .eq. 3) .and. (statusVY .eq. 3) .and. (statusVZ .eq. 3)) then
      ! exitFace = 0
      ! particle%iTrackingDomainBoundary(3) = exitFace
      ! particle%istatus = 5
      ! return

      ! kluge note: todo identify particle by "composite key" (not just irpt)
      print *, "======================================"
      print *, "Subcell with no exit face" ! kluge
      print *, "Particle ", particle%irpt
      print *, "Cell ", particle%iTrackingDomain(2)
      print *, "======================================"
      !!pause
      !!stop
    end if
    !
    ! -- Determine (earliest) exit face and corresponding travel time to exit
    exitFace = 0
    dtexit = 1.0d+30
    if ((statusVX .lt. 2) .or. (statusVY .lt. 2) .or. (statusVZ .lt. 2)) then
      ! -- Consider x-oriented faces
      dtexit = dtexitx
      if (vx .lt. 0d0) then
        exitFace = 1
      else if (vx .gt. 0) then
        exitFace = 2
      end if
      ! -- Consider y-oriented faces
      if (dtexity .lt. dtexit) then
        dtexit = dtexity
        if (vy .lt. 0d0) then
          exitFace = 3
        else if (vy .gt. 0d0) then
          exitFace = 4
        end if
      end if
      ! -- Consider z-oriented faces
      if (dtexitz .lt. dtexit) then
        dtexit = dtexitz
        if (vz .lt. 0d0) then
          exitFace = 5
        else if (vz .gt. 0d0) then
          exitFace = 6
        end if
      end if
    else
    end if
    !
    ! -- Compute exit time
    texit = particle%ttrack + dtexit
    !
    if (texit .gt. tmax) then
      ! -- The computed exit time is greater than the maximum time, so set
      ! -- final time for particle trajectory equal to maximum time and
      ! -- calculate particle location at that final time.
      t = tmax
      dt = t - particle%ttrack
      x = NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, subcellRect%dx, statusVX)
      y = NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, subcellRect%dy, statusVY)
      z = NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, subcellRect%dz, statusVZ)
      exitFace = 0
      particle%istatus = 1
    else
      ! -- The computed exit time is less than or equal to the maximum time,
      ! -- so set final time for particle trajectory equal to exit time and
      ! -- calculate exit location.
      t = texit
      dt = dtexit
      if ((exitFace .eq. 1) .or. (exitFace .eq. 2)) then
        x = 0d0
        y = NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, subcellRect%dy, statusVY)
        z = NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, subcellRect%dz, statusVZ)
        if (exitFace .eq. 2) x = 1.0d0
      else if ((exitFace .eq. 3) .or. (exitFace .eq. 4)) then
        x = NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, subcellRect%dx, statusVX)
        y = 0d0
        z = NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, subcellRect%dz, statusVZ)
        if (exitFace .eq. 4) y = 1.0d0
      else if ((exitFace .eq. 5) .or. (exitFace .eq. 6)) then
        x = NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, subcellRect%dx, statusVX)
        y = NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, subcellRect%dy, statusVY)
        z = 0d0
        if (exitFace .eq. 6) z = 1.0d0
      else
        print *, "something went wrong in track_sub" ! kluge
        !!pause                                                      ! kluge
        stop ! kluge
      end if
    end if
    !
    ! -- Set final particle location in local (unscaled) subcell coordinates,
    ! -- final time for particle trajectory, and exit face
    particle%x = x * subcellRect%dx
    particle%y = y * subcellRect%dy
    particle%z = z * subcellRect%dz
    particle%ttrack = t
    particle%iTrackingDomainBoundary(3) = exitFace
    !
    return
    !
  end subroutine track_sub

  !> @brief Calculate particle travel time to exit and exit status.
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !!
  !<
  function CalculateDT(v1, v2, dx, xL, v, dvdx, dt) result(status)
    ! dummy
    doubleprecision, intent(in) :: v1, v2, dx, xL
    doubleprecision, intent(inout) :: v, dvdx, dt
    ! result
    integer :: status
    ! local
    doubleprecision :: v2a, v1a, dv, dva, vv, vvv, zro, zrom, x, tol
    doubleprecision :: vr1, vr2, vr, v1v2
    logical :: noOutflow
    !
    ! -- Initialize variables.
    status = -1
    dt = 1.0d+20
    v2a = v2
    if (v2a .lt. 0d0) v2a = -v2a
    v1a = v1
    if (v1a .lt. 0d0) v1a = -v1a
    dv = v2 - v1
    dva = dv
    if (dva .lt. 0d0) dva = -dva
    !
    ! -- Check for a uniform zero velocity in this direction.
    ! -- If so, set status = 2 and return (dt = 1.0d+20).
    tol = 1.0d-15
    if ((v2a .lt. tol) .and. (v1a .lt. tol)) then
      v = 0d0
      dvdx = 0d0
      status = 2
      return
    end if
    !
    ! -- Check for uniform non-zero velocity in this direction.
    ! -- If so, set compute dt using the constant velocity,
    ! -- set status = 1 and return.
    vv = v1a
    if (v2a .gt. vv) vv = v2a
    vvv = dva / vv
    if (vvv .lt. 1.0d-4) then
      zro = tol
      zrom = -zro
      v = v1
      x = xL * dx
      if (v1 .gt. zro) dt = (dx - x) / v1
      if (v1 .lt. zrom) dt = -x / v1
      dvdx = 0d0
      status = 1
      return
    end if
    !
    ! -- Velocity has a linear variation.
    ! -- Compute velocity corresponding to particle position.
    dvdx = dv / dx
    v = (1.0d0 - xL) * v1 + xL * v2
    !
    ! -- If flow is into the cell from both sides there is no outflow.
    ! -- In that case, set status = 3 and return.
    noOutflow = .true.
    if (v1 .lt. 0d0) noOutflow = .false.
    if (v2 .gt. 0d0) noOutflow = .false.
    if (noOutflow) then
      status = 3
      return
    end if
    !
    ! -- If there is a divide in the cell for this flow direction, check to
    ! -- see if the particle is located exactly on the divide. If it is, move
    ! -- it very slightly to get it off the divide. This avoids possible
    ! -- numerical problems related to stagnation points.
    if ((v1 .le. 0d0) .and. (v2 .ge. 0d0)) then
      if (abs(v) .le. 0d0) then
        v = 1.0d-20
        if (v2 .le. 0d0) v = -v
      end if
    end if
    !
    ! -- If there is a flow divide, this check finds out what side of the
    ! -- divide the particle is on and sets the value of vr appropriately
    ! -- to reflect that location.
    vr1 = v1 / v
    vr2 = v2 / v
    vr = vr1
    if (vr .le. 0d0) then
      vr = vr2
    end if
    !
    ! -- If the product v1*v2 > 0, the velocity is in the same direction
    ! -- throughout the cell (i.e. no flow divide). If so, set the value
    ! -- of vr to reflect the appropriate direction.
    v1v2 = v1 * v2
    if (v1v2 .gt. 0d0) then
      if (v .gt. 0d0) vr = vr2
      if (v .lt. 0d0) vr = vr1
    end if
    !
    ! -- Check if vr is (very close to) zero.
    ! -- If so, set status = 2 and return (dt = 1.0d+20).
    if (dabs(vr) .lt. 1.0d-10) then
      v = 0d0
      dvdx = 0d0
      status = 2
      return
    end if
    !
    ! -- Compute travel time to exit face. Return with status = 0.
    dt = log(vr) / dvdx
    status = 0
    !
  end function CalculateDT

  !> @brief Update particle coordinates based on time increment
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !!
  !<
  function NewXYZ(v, dvdx, v1, v2, dt, x, dx, velocityProfileStatus) result(newX)
    ! dummy
    integer, intent(in) :: velocityProfileStatus
    doubleprecision, intent(in) :: v, dvdx, v1, v2, dt, x, dx
    ! result
    doubleprecision :: newX
    !
    newX = x
    select case (velocityProfileStatus)
    case (1)
      newX = newX + (v1 * dt / dx)
    case default
      if (v .ne. 0d0) then
        newX = newX + (v * (exp(dvdx * dt) - 1.0d0) / dvdx / dx)
      end if
    end select
    if (newX .lt. 0d0) newX = 0d0
    if (newX .gt. 1.0d0) newX = 1.0d0
    !
  end function NewXYZ

end module MethodSubcellPollockModule
