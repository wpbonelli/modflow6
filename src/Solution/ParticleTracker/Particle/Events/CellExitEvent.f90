module CellExitEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: CellExitEventType

  type, extends(ParticleEventType) :: CellExitEventType
    integer(I4B) :: exit_face ! face through which the particle exited
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_str
  end type CellExitEventType

contains

  function get_code(this) result(code)
    class(CellExitEventType), intent(in) :: this
    integer(I4B) :: code
    code = 1
  end function get_code

  function get_verb(this) result(verb)
    class(CellExitEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'exited cell'
  end function get_verb

  function get_str(this) result(str)
    class(CellExitEventType), intent(in) :: this
    character(len=:), allocatable :: str
    character(len=LENHUGELINE) :: temp

    write (temp, '(*(G0))') &
      'Particle from model ', this%imdl, &
      ', package ', this%iprp, &
      ', point ', this%irpt, &
      ', time ', this%trelease, &
      ' '//this%get_verb()// &
      ' in layer ', this%ilay, &
      ', cell ', this%icu, &
      ', zone ', this%izone, &
      ' through face ', this%exit_face, &
      ' at x ', this%x, &
      ', y ', this%y, &
      ', z ', this%z, &
      ', time ', this%ttrack, &
      ', period ', this%kper, &
      ', timestep ', this%kstp, &
      ' with status ', this%istatus
    str = trim(adjustl(temp))
  end function get_str

end module CellExitEventModule
