module ReleaseEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: ReleaseEventType

  type, extends(ParticleEventType) :: ReleaseEventType
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_str
  end type ReleaseEventType

contains

  function get_code(this) result(code)
    class(ReleaseEventType), intent(in) :: this
    integer(I4B) :: code
    code = 0
  end function get_code

  function get_verb(this) result(verb)
    class(ReleaseEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'released'
  end function get_verb

  function get_str(this) result(str)
    class(ReleaseEventType), intent(in) :: this
    character(len=:), allocatable :: str
    str = this%ParticleEventType%get_str()
  end function get_str

end module ReleaseEventModule
