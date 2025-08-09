module FeatExitEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, FEATEXIT
  implicit none

  private
  public :: FeatExitEventType

  type, extends(ParticleEventType) :: FeatExitEventType
  contains
    procedure :: get_code
    procedure :: get_verb
  end type FeatExitEventType

contains

  function get_code(this) result(code)
    class(FeatExitEventType), intent(in) :: this
    integer(I4B) :: code
    code = FEATEXIT
  end function get_code

  function get_verb(this) result(verb)
    class(FeatExitEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'exited grid feature'
  end function get_verb

end module FeatExitEventModule
