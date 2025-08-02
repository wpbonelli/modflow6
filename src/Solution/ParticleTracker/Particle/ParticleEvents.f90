module ParticleEventsModule
  use KindModule, only: DP, I4B, LGP
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, &
                                 ReleaseEventType, &
                                 CellExitEventType, &
                                 TimestepEventType, &
                                 TerminationEventType, &
                                 WeakSinkEventType, &
                                 UserTimeEventType
  implicit none

  private

  type, public, abstract :: ParticleEventConsumerType
  contains
    procedure(handle_event), deferred :: handle_event
  end type ParticleEventConsumerType

  type, public :: ParticleEventDispatcherType
    class(ParticleEventConsumerType), pointer :: consumer => null()
  contains
    procedure, public :: subscribe
    procedure, public :: unsubscribe
    procedure, public :: dispatch
    procedure :: destroy
  end type ParticleEventDispatcherType

  abstract interface
    subroutine handle_event(this, event)
      import ParticleEventConsumerType, ParticleType, ParticleEventType
      class(ParticleEventConsumerType), intent(inout) :: this
      class(ParticleEventType), pointer, intent(in) :: event
    end subroutine handle_event
  end interface

contains
  !> @brief Subscribe a consumer to the dispatcher.
  subroutine subscribe(this, consumer)
    class(ParticleEventDispatcherType), intent(inout) :: this
    class(ParticleEventConsumerType), target, intent(inout) :: consumer
    this%consumer => consumer
  end subroutine subscribe

  !> @brief Unsubscribe the consumer from the dispatcher.
  subroutine unsubscribe(this)
    class(ParticleEventDispatcherType), intent(inout) :: this
    if (associated(this%consumer)) then
      deallocate (this%consumer)
      this%consumer => null()
    end if
  end subroutine unsubscribe

  !> @brief Dispatch an event.
  subroutine dispatch(this, particle, event)
    use TdisModule, only: kper, kstp, totimc
    ! dummy
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(inout) :: event
    ! local
    integer(I4B) :: per, stp

    per = kper
    stp = kstp

    ! If tracking time falls exactly on a boundary between time steps,
    ! report the previous time step for this datum. This is to follow
    ! MP7's behavior, and because the particle will have been tracked
    ! up to this instant under the previous time step's conditions, so
    ! the time step we're about to start shouldn't get "credit" for it.
    if (particle%ttrack == totimc .and. (per > 1 .or. stp > 1)) then
      if (stp > 1) then
        stp = stp - 1
      else if (per > 1) then
        per = per - 1
        stp = 1
      end if
    end if

    event%particle_id%imdl = particle%imdl
    event%particle_id%iprp = particle%iprp
    event%particle_id%irpt = particle%irpt
    event%particle_id%trelease = particle%trelease
    event%particle_id%name = trim(adjustl(particle%name))
    event%tdis_coords%kper = per
    event%tdis_coords%kstp = stp
    event%grid_coords%ilay = particle%ilay
    event%grid_coords%icu = particle%icu
    event%grid_coords%izone = particle%izone
    event%time_coords%ttrack = particle%ttrack
    call particle%get_model_coords( &
      event%space_coords%x, &
      event%space_coords%y, &
      event%space_coords%z)
    event%istatus = particle%istatus
    event%ireason = event%get_code()

    call this%consumer%handle_event(event)
    deallocate (event)
  end subroutine dispatch

  !> @brief Destroy the dispatcher.
  subroutine destroy(this)
    class(ParticleEventDispatcherType), intent(inout) :: this
    if (associated(this%consumer)) deallocate (this%consumer)
  end subroutine destroy

end module ParticleEventsModule
