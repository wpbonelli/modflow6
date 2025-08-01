module ParticleEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
  public :: CellExitEventType, TerminationEventType, ReleaseEventType
  public :: TimeStepEventType, WeakSinkEventType, UserTimeEventType
  public :: SubCellExitEventType
  public :: RELEASE, CELLEXIT, TIMESTEP, TERMINATE, WEAKSINK, USERTIME

  !> @brief Particle event enumeration.
  !!
  !! A number of events may occur to particles, each of which may (or may
  !! not) be of interest to the user. The user selects events to report.
  !<
  enum, bind(C)
    enumerator :: RELEASE = 0 !< particle was released
    enumerator :: CELLEXIT = 1 !< particle exited a cell
    enumerator :: TIMESTEP = 2 !< time step ended
    enumerator :: TERMINATE = 3 !< particle terminated
    enumerator :: WEAKSINK = 4 !< particle exited a weak sink
    enumerator :: USERTIME = 5 !< user-specified tracking time
    enumerator :: SUBCELLEXIT = 6 !< particle exited a subcell
  end enum

  !> @brief Base type for particle events.
  !!
  !! Events may be identical except for their type/code, reflecting the
  !! fact that several events of interest may occur at a given moment.
  type, abstract :: ParticleEventType
    type(ParticleType), pointer :: particle => null() ! particle causing the event
    integer(I4B) :: code = -1 ! event code
    integer(I4B) :: kper = 0, kstp = 0 ! period and step
    real(DP) :: time = 0.0_DP ! simulation time
    ! particle state snapshot for cycle detection
    integer(I4B) :: icu = 0 ! cell number at event
    integer(I4B) :: ilay = 0 ! layer at event
    integer(I4B) :: izone = 0 ! zone at event
    real(DP) :: x = 0.0_DP ! x coordinate at event
    real(DP) :: y = 0.0_DP ! y coordinate at event
    real(DP) :: z = 0.0_DP ! z coordinate at event
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_str
    procedure :: log
  end type ParticleEventType

  type, extends(ParticleEventType) :: CellExitEventType
    integer(I4B) :: exit_face = 0 ! exit face number
  end type CellExitEventType

  type, extends(ParticleEventType) :: TerminationEventType
  end type TerminationEventType

  type, extends(ParticleEventType) :: ReleaseEventType
  end type ReleaseEventType

  type, extends(ParticleEventType) :: TimeStepEventType
  end type TimestepEventType

  type, extends(ParticleEventType) :: WeakSinkEventType
  end type WeakSinkEventType

  type, extends(ParticleEventType) :: UserTimeEventType
  end type UserTimeEventType

  type, extends(ParticleEventType) :: SubCellExitEventType
    integer(I4B) :: exit_face = 0 ! exit face number
    integer(I4B) :: isc = 0 ! subcell number
  end type SubCellExitEventType

contains
  integer function get_code(this) result(code)
    class(ParticleEventType), intent(in) :: this

    select type (this)
    type is (ReleaseEventType); code = 0
    type is (CellExitEventType); code = 1
    type is (TimeStepEventType); code = 2
    type is (TerminationEventType); code = 3
    type is (WeakSinkEventType); code = 4
    type is (UserTimeEventType); code = 5
    type is (SubCellExitEventType); code = 6
    class default; call pstop(1, "unknown event type")
    end select
  end function get_code

  function get_verb(this) result(str)
    class(ParticleEventType), intent(in) :: this
    character(len=:), allocatable :: str

    select type (this)
    type is (ReleaseEventType); str = "released"
    type is (CellExitEventType); str = "exited cell"
    type is (SubCellExitEventType); str = "exited subcell"
    type is (TimeStepEventType); str = "completed timestep"
    type is (TerminationEventType); str = "terminated"
    type is (WeakSinkEventType); str = "exited weak sink"
    type is (UserTimeEventType); str = "user-specified tracking time"
    class default; call pstop(1, "unknown event type")
    end select
  end function get_verb

  function get_str(this, partid, cellid, space, time) result(str)
    class(ParticleEventType), intent(in) :: this
    logical(LGP), intent(in), optional :: partid, cellid, space, time
    character(len=:), allocatable :: str
    ! local
    character(LENHUGELINE) :: temp
    character(len=:), allocatable :: rptname_, &
                                     partid_, &
                                     cellid_, &
                                     space_, &
                                     time_, &
                                     status_
    logical(LGP) :: lpartid, lcellid, lspace, ltime

    if (present(partid)) then
      lpartid = partid
    else
      lpartid = .true.
    end if

    if (present(cellid)) then
      lcellid = cellid
    else
      lcellid = .true.
    end if

    if (present(space)) then
      lspace = space
    else
      lspace = .true.
    end if

    if (present(time)) then
      ltime = time
    else
      ltime = .true.
    end if

    rptname_ = trim(adjustl(this%particle%name))
    if (len_trim(rptname_) == 0) then
      rptname_ = ''
    else
      rptname_ = ' ('//rptname_//')'
    end if
    if (lpartid) then
      write (temp, '(*(G0))') &
        'from model ', this%particle%imdl, &
        ', package ', this%particle%iprp, &
        ', point ', this%particle%irpt, rptname_, &
        ', time ', this%particle%trelease
      partid_ = trim(adjustl(temp))
    else
      partid_ = ''
    end if
    write (temp, '(*(G0))') ' ', this%get_verb()
    if (lcellid) then
      write (temp, '(*(G0))') &
        ' in layer ', this%particle%ilay, &
        ', cell ', this%particle%icu, &
        ', zone ', this%particle%izone
      cellid_ = trim(adjustl(temp))
    else
      cellid_ = ''
    end if
    if (lspace) then
      write (temp, '(*(G0))') &
        ' at position x ', this%particle%x, &
        ', y ', this%particle%y, &
        ', z ', this%particle%z
      space_ = trim(adjustl(temp))
    else
      space_ = ''
    end if
    if (ltime) then
      write (temp, '(*(G0))') &
        ' at time t ', this%particle%ttrack, &
        ', period ', this%kper, &
        ', timestep ', this%kstp
      time_ = trim(adjustl(temp))
    else
      time_ = ''
    end if
    write (temp, '(*(G0))') ' with status ', this%particle%istatus
    status_ = trim(adjustl(temp))

    str = trim(adjustl('Particle '// &
                       partid_//' '// &
                       this%get_verb()//' '// &
                       cellid_//' '// &
                       space_//' '// &
                       time_//' '// &
                       status_))

  end function get_str

  subroutine log(this, iun, partid, cellid, space, time)
    class(ParticleEventType), intent(inout) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: partid, cellid, space, time

    if (iun < 0) return
    write (iun, '(*(G0))') get_str(this, partid, cellid, space, time)
  end subroutine log

end module ParticleEventModule
