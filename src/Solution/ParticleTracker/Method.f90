module MethodModule

  use KindModule, only: DP, I4B
  use GlobalDataModule
  use ParticleModule ! kluge???
  use TrackDataModule, only: TrackDataType
  implicit none

  private
  public :: MethodType

  type, abstract :: MethodType
    ! private
    character(len=40), pointer, public :: trackingDomainType ! character string that names the tracking domain type
    logical, public :: delegatesTracking
    type(TrackDataType), pointer :: trackdata
  contains
    ! -- Implemented in all tracking methods
    procedure(apply), deferred :: apply ! applies the method
    ! -- Implemented (overridden) in tracking methods that delegate to submethods
    procedure :: pass ! passes a particle to the next tracking subdomain
    procedure :: loadsub ! loads tracking submethod
    ! -- Implemented in this base class
    procedure :: subtrack ! tracks the particle across subdomains    ! kluge, description, rename???
    procedure :: advance ! advances the particle                    ! kluge, description, rename???
  end type MethodType

  abstract interface

    !> @brief Must be implemented by all tracking methods
    subroutine apply(this, particle, tmax)
      import DP
      import MethodType
      import ParticleType
      class(MethodType), intent(inout) :: this
      type(ParticleType), pointer, intent(inout) :: particle ! kluge???
      real(DP), intent(in) :: tmax
      ! doubleprecision :: initialTime,maximumTime,t   ! kluge not in arg list yet
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
    ! doubleprecision :: initialTime,maximumTime,t ! kluge
    ! local
    integer :: levelNext
    class(methodType), pointer :: submethod
    logical :: isStillAdvancing
    character(len=40) :: formatString ! kluge???
    !
    ! write(formatString,'(A,I2,A)') "(1X,",(level+1)*4,"X,A,I)"                             ! kluge???
    ! write(*,formatString) TRIM(this%trackingDomain%type), particle%iTrackingDomain(level)  ! kluge???
    ! write(*,formatString) TRIM(this%trackingDomainType), particle%iTrackingDomain(level)  ! kluge???
    ! !
    ! -- Set next level
    levelNext = level + 1
    !
    ! -- Track particle through next-level tracking domains (subdomains)
    ! -- where each domain is a grid cell (or subcell)
    isStillAdvancing = .true.
    do while (isStillAdvancing)
      !
      ! -- Load subdomain tracking method (submethod)
      call this%loadsub(particle, levelNext, submethod)
      !
      ! -- Delegate tracking to the submethod, which either does the tracking
      ! -- calculations or calls this subroutine (subtrack) recursively to
      ! -- continue delegating
      call submethod%apply(particle, tmax)
      !
      ! -- Advance particle
      call advance(this, particle, levelNext, submethod, isStillAdvancing)
      !
    end do
    !
    return
    !
  end subroutine subtrack

  !> @brief Advance the particle
  subroutine advance(this, particle, levelNext, submethod, isStillAdvancing)
    ! dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer :: levelNext
    class(MethodType), pointer :: submethod
    logical :: isStillAdvancing
    ! local
    character(len=40) :: typeSubdomain ! kluge???
    integer :: ip, iSubdomain
    !
    ! ! -- Check whether the particle stopped advancing within the last tracking
    ! ! -- subdomain or has exited the tracking domain. If not, advance to the next
    ! ! -- subdomain.
    ! ip = particle%ipart
    ! typeSubdomain = submethod%trackingDomainType
    ! if (particle%iTrackingDomainBoundary(levelNext).eq.0) then
    !   particle%iTrackingDomainBoundary = 0        ! kluge???
    !   isStillAdvancing = .false.
    ! else
    !   call this%pass(particle)
    !   iSubdomain = particle%iTrackingDomain(levelNext)
    !   if (iSubdomain.lt.0) then
    !     isStillAdvancing = .false.
    !   end if
    ! end if
    ! -- Check whether the particle terminated within the last tracking
    ! -- subdomain. If not, advance to the next subdomain, if there is one.
    ip = particle%ipart
    ! if (particle%iTrackingDomainBoundary(levelNext).eq.0) then
    if (particle%istatus .ne. -1) then
      particle%iTrackingDomainBoundary = 0 ! kluge???
      isStillAdvancing = .false.
    else
      call this%pass(particle)
      if (particle%iTrackingDomainBoundary(levelNext - 1) .ne. 0) then
        isStillAdvancing = .false.
      end if
    end if
    ! ip = particle%ipart
    ! if (particle%istatus.eq.-1) then
    !   ! -- Particle not finished or terminated
    !   if (particle%iTrackingDomainBoundary(levelNext-1).eq.0) then  ! kluge note: pass in level instead of levelNext???
    !     ! -- Particle in domain interior, so pass to next subdomain
    !     call this%pass(particle)
    !     if (particle%istatus.ne.-1) then
    !       ! -- Particle finished or terminated, so no longer advancing
    !       ! -- in current domain
    !       isStillAdvancing = .false.
    !     else if (particle%iTrackingDomainBoundary(levelNext-1).ne.0) then
    !       ! -- Particle is at domain boundary, so no longer advancing
    !       ! -- in current domain
    !       isStillAdvancing = .false.
    !     end if
    !   else
    !     ! -- Particle is at domain boundary, so no longer advancing
    !     ! -- in current domain
    !     isStillAdvancing = .false.
    !   end if
    ! else
    !   ! -- Particle finished or terminated, so no longer advancing
    !   ! -- in current domain
    !   isStillAdvancing = .false.
    ! end if
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
    !!pause                                             ! kluge
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
    !!pause                                             ! kluge
    stop
    !
    return
    !
  end subroutine pass

end module MethodModule
