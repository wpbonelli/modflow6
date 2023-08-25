module MethodModule

  use KindModule, only: DP, I4B
  use GlobalDataModule
  use ParticleModule
  use CellDefnModule, only: CellDefnType
  use TrackModule, only: TrackControlType
  implicit none

  private
  public :: MethodType

  type, abstract :: MethodType
    ! private
    character(len=40), pointer, public :: trackingDomainType ! character string that names the tracking domain type
    logical, public :: delegatesTracking ! whether the method delegates tracking to an internal method
    type(TrackControlType), pointer :: trackctl ! the track output control object
  contains
    ! -- Implemented in all tracking methods
    procedure(apply), deferred :: apply ! applies the method
    ! -- Implemented (overridden) in tracking methods that delegate to submethods
    procedure :: pass ! passes a particle to the next tracking subdomain
    procedure :: loadsub ! loads tracking submethod
    ! -- Implemented in this base class
    procedure :: subtrack ! tracks the particle across subdomains
    procedure :: advance ! advances the particle
    procedure :: update ! update particle state, terminating if appropriate and reporting
  end type MethodType

  abstract interface

    !> @brief Must be implemented by all tracking methods
    subroutine apply(this, particle, tmax)
      import DP
      import MethodType
      import ParticleType
      class(MethodType), intent(inout) :: this
      type(ParticleType), pointer, intent(inout) :: particle
      real(DP), intent(in) :: tmax
    end subroutine apply

  end interface

contains

  !> @brief Track particle across subdomains
  recursive subroutine subtrack(this, particle, level, tmax)
    ! dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer :: level
    real(DP), intent(in) :: tmax
    ! local
    integer :: levelNext
    class(methodType), pointer :: submethod
    logical :: advancing
    !
    ! -- Set next level
    levelNext = level + 1
    !
    ! -- Track particle through next-level tracking domains (subdomains)
    ! -- where each domain is a grid cell (or subcell)
    advancing = .true.
    do while (advancing)
      ! -- Load subdomain tracking method (submethod)
      call this%loadsub(particle, levelNext, submethod)

      ! -- Delegate tracking to the submethod, which either does the tracking
      ! -- calculations or calls this subroutine (subtrack) recursively to
      ! -- continue delegating
      call submethod%apply(particle, tmax)

      ! -- Advance particle
      call advance(this, particle, levelNext, submethod, advancing)
    end do
    !
    return
    !
  end subroutine subtrack

  !> @brief Advance the particle
  subroutine advance(this, particle, levelNext, submethod, advancing)
    ! dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer :: levelNext
    class(MethodType), pointer :: submethod
    logical :: advancing
    !
    if (.not. particle%advancing) then
      particle%iTrackingDomainBoundary = 0
      advancing = .false.
    else
      call this%pass(particle)
      if (particle%iTrackingDomainBoundary(levelNext - 1) .ne. 0) &
        advancing = .false.
    end if
    !
    return
    !
  end subroutine advance

  !> @brief Load subdomain tracking method (submethod)
  subroutine loadsub(this, particle, levelNext, submethod)
    ! -- dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: levelNext
    class(MethodType), pointer, intent(inout) :: submethod
    !
    ! -- Implemented by all tracking methods that delegate tracking
    ! -- to a submethod
    !
    write (*, '(A)') "Type-bound procedure 'loadsub' is not implemented"
    write (*, '(A)') "in the 'Method' base class; it must be implemented by all"
    write (*, '(A)') "tracking methods that delegate tracking to a submethod"
    stop
    !
    return
    !
  end subroutine loadsub

  !> @brief Pass a particle to the next subcell
  subroutine pass(this, particle)
    ! -- dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    !
    ! -- Implemented by all tracking methods that delegate tracking
    ! -- to a submethod
    !
    write (*, '(A)') "Type-bound procedure 'pass' is not implemented in"
    write (*, '(A)') "the 'Method' base class; it must be implemented by all"
    write (*, '(A)') "tracking methods that delegate tracking to a submethod"
    stop
    !
    return
    !
  end subroutine pass

  !> @brief Update particle state and check termination conditions
  !!
  !! Update the particle's properties (e.g. advancing flag, zone number,
  !! status). If any termination conditions apply, the particle's status
  !! will be set to the appropriate termination value. If any reporting
  !! conditions apply, save particle state with the proper reason code.
  subroutine update(this, particle, celldefn)
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: celldefn
    particle%izone = celldefn%izone

    if (celldefn%izone .ne. 0) then
      if (particle%istopzone .eq. celldefn%izone) then
        particle%advancing = .false.
        particle%istatus = 6
        call this%trackctl%save_record(particle, kper=kper, &
                                       kstp=kstp, reason=3) ! reason=3: termination
      end if
    else if (celldefn%inoexitface .ne. 0) then
      particle%advancing = .false.
      particle%istatus = 5
      call this%trackctl%save_record(particle, kper=kper, &
                                     kstp=kstp, reason=3) ! reason=3: termination
    else if (celldefn%iweaksink .ne. 0) then
      if (particle%istopweaksink .ne. 0) then
        particle%advancing = .false.
        particle%istatus = 3
        call this%trackctl%save_record(particle, kper=kper, &
                                       kstp=kstp, reason=3) ! reason=3: termination
      else
        call this%trackctl%save_record(particle, kper=kper, &
                                       kstp=kstp, reason=4) ! reason=4: exited weak sink
      end if
    end if

  end subroutine update

end module MethodModule
