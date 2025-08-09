module ReleaseEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, RELEASE
  implicit none

  private
  public :: ReleaseEventType

  type, extends(ParticleEventType) :: ReleaseEventType
  contains
    procedure :: get_code
    procedure :: get_verb
  end type ReleaseEventType

contains

  function get_code(this) result(code)
    class(ReleaseEventType), intent(in) :: this
    integer(I4B) :: code
    code = RELEASE
  end function get_code

  function get_verb(this) result(verb)
    class(ReleaseEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'released'
  end function get_verb

end module ReleaseEventModule
