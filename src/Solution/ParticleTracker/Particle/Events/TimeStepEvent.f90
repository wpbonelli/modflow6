module TimeStepEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: TimeStepEventType

  type, extends(ParticleEventType) :: TimeStepEventType
  contains
    procedure :: get_code
    procedure :: get_verb
  end type TimeStepEventType

contains

  function get_code(this) result(code)
    class(TimeStepEventType), intent(in) :: this
    integer(I4B) :: code
    code = 2
  end function get_code

  function get_verb(this) result(verb)
    class(TimeStepEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'completed timestep'
  end function get_verb

end module TimeStepEventModule
