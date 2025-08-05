module TerminationEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: TerminationEventType

  type, extends(ParticleEventType) :: TerminationEventType
  contains
    procedure :: get_code
    procedure :: get_verb
  end type TerminationEventType

contains

  function get_code(this) result(code)
    class(TerminationEventType), intent(in) :: this
    integer(I4B) :: code
    code = 3
  end function get_code

  function get_verb(this) result(verb)
    class(TerminationEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'terminated'
  end function get_verb

end module TerminationEventModule
