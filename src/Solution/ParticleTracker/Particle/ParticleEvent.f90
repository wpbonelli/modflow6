module ParticleEventModule
  use ConstantsModule, only: DZERO
  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
  public :: CellExitEventType, TerminationEventType, ReleaseEventType
  public :: TimeStepEventType, WeakSinkEventType, UserTimeEventType
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
  end enum

  type :: ParticleIdType
    integer(I4B) :: imdl = 0 !< release model ID
    integer(I4B) :: iprp = 0 !< release package ID
    integer(I4B) :: irpt = 0 !< release point ID
    real(DP) :: trelease = DZERO !< release time
    character(len=:), allocatable :: name !< particle name
  end type ParticleIdType

  type :: GridCoordsType
    integer(I4B) :: ilay = 0 !< layer ID
    integer(I4B) :: icu = 0 !< cell ID
    integer(I4B) :: izone = 0 !< zone ID
  end type GridCoordsType

  type :: SpaceCoordsType
    real(DP) :: x = 0.0_DP !< x coordinate
    real(DP) :: y = 0.0_DP !< y coordinate
    real(DP) :: z = 0.0_DP !< z coordinate
  end type SpaceCoordsType

  type :: TimeCoordsType
    real(DP) :: ttrack = 0.0_DP !< time of the event
  end type TimeCoordsType

  type :: TdisCoordsType
    integer(I4B) :: kper = 0 !< stress period ID
    integer(I4B) :: kstp = 0 !< time step ID
  end type TdisCoordsType

  !> @brief Base type for particle events.
  !!
  !! Events may be identical except for their type/code, reflecting the
  !! fact that several events of interest may occur at a given moment.
  type, abstract :: ParticleEventType
    type(ParticleIdType) :: particle_id ! particle identity
    type(TdisCoordsType) :: tdis_coords ! stress period and time step
    type(GridCoordsType) :: grid_coords ! cell identity
    type(TimeCoordsType) :: time_coords ! time of the event
    type(SpaceCoordsType) :: space_coords ! particle location
    integer(I4B) :: istatus = -1 !< particle status
    integer(I4B) :: ireason = -1 ! event reason
  contains
    procedure :: get_code
    procedure :: get_verb
    ! procedure :: get_str
  end type ParticleEventType

  type, extends(ParticleEventType) :: CellExitEventType
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
    class default; call pstop(1, "unknown event type")
    end select
  end function get_code

  function get_verb(this) result(str)
    class(ParticleEventType), intent(in) :: this
    character(len=:), allocatable :: str

    select type (this)
    type is (ReleaseEventType); str = "released"
    type is (CellExitEventType); str = "exited cell"
    type is (TimeStepEventType); str = "completed timestep"
    type is (TerminationEventType); str = "terminated"
    type is (WeakSinkEventType); str = "exited weak sink"
    type is (UserTimeEventType); str = "user-specified tracking time"
    class default; call pstop(1, "unknown event type")
    end select
  end function get_verb
!
end module ParticleEventModule
