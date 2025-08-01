module MethodCellModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DONE, DZERO
  use MethodModule, only: MethodType
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, &
                                 CellExitEventType
  use CellDefnModule, only: CellDefnType
  use IteratorModule, only: IteratorType
  implicit none

  private
  public :: MethodCellType

  type, abstract, extends(MethodType) :: MethodCellType
  contains
    procedure, public :: cellexit
    procedure, public :: assess
  end type MethodCellType

contains

  !> @brief Particle exits a cell.
  subroutine cellexit(this, particle)
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event

    allocate (CellExitEventType :: event)
    select type (event)
    type is (CellExitEventType)
      event%exit_face = particle%iboundary(2)
      call check_cycle(particle, event)
      call store_event(particle, event)
    end select
    call this%events%dispatch(particle, event)
  end subroutine cellexit

  !> @brief Check for cycles in the particle's path.
  subroutine check_cycle(particle, event)
    class(ParticleType), intent(inout) :: particle !< particle
    class(CellExitEventType), pointer, intent(in) :: event !< particle event to check for
    ! local
    class(IteratorType), allocatable :: itr
    logical(LGP) :: found_cycle

    found_cycle = .false.
    itr = particle%history%Iterator()
    do while (itr%has_next())
      call itr%next()
      select type (prev => itr%value())
      class is (CellExitEventType)
        if (event%icu == prev%icu .and. &
            event%ilay == prev%ilay .and. &
            event%izone == prev%izone .and. &
            event%exit_face == prev%exit_face .and. &
            event%exit_face /= 0) then
          found_cycle = .true.
          exit
        end if
      end select
    end do

    if (found_cycle) then
      print *, "Cycle in particle path:"
      itr = particle%history%Iterator()
      do while (itr%has_next())
        call itr%next()
        select type (e => itr%value())
        class is (ParticleEventType)
          print *, e%get_str(time=.false.)
        end select
      end do
      call pstop(1, "Particle cycling detected")
    end if
  end subroutine check_cycle

  !> @brief Store the cell exit event in the particle's history.
  subroutine store_event(particle, event)
    class(ParticleType), intent(inout) :: particle
    class(CellExitEventType), pointer, intent(in) :: event
    ! local
    class(*), pointer :: p

    p => event
    call particle%history%Add(p)
    if (particle%history%Count() > 6) &
      call particle%history%remove_node_by_index(1, .true.)
  end subroutine store_event

  !> @brief Check reporting/terminating conditions before tracking
  !! the particle across the cell.
  !!
  !! Check a number of conditions determining whether to continue
  !! tracking the particle or terminate it, as well as whether to
  !! record any output data as per selected reporting conditions.
  !<
  subroutine assess(this, particle, cell_defn, tmax)
    ! modules
    use TdisModule, only: endofsimulation, totimc, totim
    use ParticleModule, only: TERM_WEAKSINK, TERM_NO_EXITS, &
                              TERM_STOPZONE, TERM_INACTIVE
    use ParticleEventModule, only: CELLEXIT, TERMINATE, &
                                   TIMESTEP, WEAKSINK, USERTIME
    ! dummy
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! local
    logical(LGP) :: dry_cell, dry_particle, no_exit_face, stop_zone, weak_sink
    integer(I4B) :: i
    real(DP) :: t, ttrackmax

    dry_cell = this%fmi%ibdgwfsat0(cell_defn%icell) == 0
    dry_particle = particle%z > cell_defn%top
    no_exit_face = cell_defn%inoexitface > 0
    stop_zone = cell_defn%izone > 0 .and. particle%istopzone == cell_defn%izone
    weak_sink = cell_defn%iweaksink > 0

    particle%izone = cell_defn%izone
    if (stop_zone) then
      call this%terminate(particle, status=TERM_STOPZONE)
      return
    end if

    if (no_exit_face .and. .not. dry_cell) then
      call this%terminate(particle, status=TERM_NO_EXITS)
      return
    end if

    if (weak_sink) then
      if (particle%istopweaksink > 0) then
        call this%terminate(particle, status=TERM_WEAKSINK)
        return
      else
        call this%weaksink(particle)
      end if
    end if

    if (dry_cell) then
      if (particle%idrymeth == 0) then
        ! drop to cell bottom. handled by pass
        ! to bottom method, nothing to do here
        no_exit_face = .false.
      else if (particle%idrymeth == 1) then
        ! stop
        call this%terminate(particle, status=TERM_INACTIVE)
        return
      else if (particle%idrymeth == 2) then
        ! stay
        particle%advancing = .false.
        no_exit_face = .false.

        ! we might report tracking times
        ! out of order here, but we want
        ! the particle termination event
        ! (if this is the last time step)
        ! to have the maximum tracking t,
        ! so we need to keep tabs on it.
        ttrackmax = totim

        ! update tracking time to time
        ! step end time and save record
        particle%ttrack = totim
        call this%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if

        ! terminate if last period/step
        if (endofsimulation) then
          particle%ttrack = ttrackmax
          call this%terminate(particle, status=TERM_NO_EXITS)
          return
        end if
      end if
    else if (dry_particle .and. this%name /= "passtobottom") then
      if (particle%idrymeth == 0) then
        ! drop to water table
        particle%z = cell_defn%top
        call this%cellexit(particle)
      else if (particle%idrymeth == 1) then
        ! stop
        call this%terminate(particle, status=TERM_INACTIVE)
        return
      else if (particle%idrymeth == 2) then
        ! stay
        particle%advancing = .false.
        no_exit_face = .false.

        ! we might report tracking times
        ! out of order here, but we want
        ! the particle termination event
        ! (if this is the last time step)
        ! to have the maximum tracking t,
        ! so we need to keep tabs on it.
        ttrackmax = totim

        ! update tracking time to time
        ! step end time and save record
        particle%ttrack = totim
        call this%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if
      end if
    end if

    if (no_exit_face) then
      particle%advancing = .false.
      particle%istatus = TERM_NO_EXITS
      call this%terminate(particle)
      return
    end if

  end subroutine assess

end module MethodCellModule
