module WeakSinkEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: WeakSinkEventType

  type, extends(ParticleEventType) :: WeakSinkEventType
  contains
    procedure :: get_code
    procedure :: get_verb
  end type WeakSinkEventType

contains

  function get_code(this) result(code)
    class(WeakSinkEventType), intent(in) :: this
    integer(I4B) :: code
    code = 4
  end function get_code

  function get_verb(this) result(verb)
    class(WeakSinkEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'exited weak sink'
  end function get_verb

end module WeakSinkEventModule
