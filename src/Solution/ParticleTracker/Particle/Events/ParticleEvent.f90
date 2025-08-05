module ParticleEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
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

  !> @brief Base type for particle events.
  !!
  !! Events may be identical except for their type/code, reflecting the
  !! fact that several events of interest may occur at a given moment.
  type, abstract :: ParticleEventType
    integer(I4B) :: imdl, iprp, irpt ! release model, package, and point
    real(DP) :: trelease = 0.0_DP ! release time
    integer(I4B) :: kper = 0, kstp = 0 ! period and step
    integer(I4B) :: ilay, icu, izone = 0
    real(DP) :: ttrack = 0.0_DP ! simulation time
    real(DP) :: x = 0.0_DP, y = 0.0_DP, z = 0.0_DP ! particle position
    integer(I4B) :: istatus = -1 ! status code
  contains
    procedure(get_code_if), deferred :: get_code
    procedure(get_verb_if), deferred :: get_verb
    procedure :: get_str
    procedure :: log
  end type ParticleEventType

  abstract interface
    function get_code_if(this) result(code)
      IMPORT I4B
      IMPORT ParticleEventType
      class(ParticleEventType), intent(in) :: this
      integer(I4B) :: code
    end function get_code_if

    function get_verb_if(this) result(verb)
      IMPORT ParticleEventType
      class(ParticleEventType), intent(in) :: this
      character(len=:), allocatable :: verb
    end function get_verb_if
  end interface

contains
  function get_str(this) result(str)
    class(ParticleEventType), intent(in) :: this
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
      ' at x ', this%x, &
      ', y ', this%y, &
      ', z ', this%z, &
      ', time ', this%ttrack, &
      ', period ', this%kper, &
      ', timestep ', this%kstp, &
      ' with status ', this%istatus
    str = trim(adjustl(temp))
  end function get_str

  subroutine log(this, iun)
    class(ParticleEventType), intent(inout) :: this
    integer(I4B), intent(in) :: iun

    if (iun >= 0) write (iun, '(*(G0))') this%get_str()
  end subroutine log

end module ParticleEventModule
