module UserTimeEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: UserTimeEventType

  type, extends(ParticleEventType) :: UserTimeEventType
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_str
  end type UserTimeEventType

contains

  function get_code(this) result(code)
    class(UserTimeEventType), intent(in) :: this
    integer(I4B) :: code
    code = 5
  end function get_code

  function get_verb(this) result(verb)
    class(UserTimeEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'user-specified tracking time'
  end function get_verb

  function get_str(this) result(str)
    class(UserTimeEventType), intent(in) :: this
    character(len=:), allocatable :: str
    str = this%ParticleEventType%get_str()
  end function get_str

end module UserTimeEventModule
