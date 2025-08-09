module ParticleEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
  public :: RELEASE
  public :: FEATEXIT
  public :: TIMESTEP
  public :: TERMINATE
  public :: WEAKSINK
  public :: USERTIME

  !> @brief Particle event enumeration.
  !!
  !! A number of events may occur to particles, each of which may (or may
  !! not) be of interest to the user. The user selects events to report.
  !<
  enum, bind(C)
    enumerator :: RELEASE = 0 !< particle was released
    enumerator :: FEATEXIT = 1 !< particle exited a grid feature
    enumerator :: TIMESTEP = 2 !< time step ended
    enumerator :: TERMINATE = 3 !< particle terminated
    enumerator :: WEAKSINK = 4 !< particle entered a weak sink
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
    procedure(get_int), deferred :: get_code
    procedure(get_str), deferred :: get_verb
    procedure :: get_text
    procedure :: log
  end type ParticleEventType

  abstract interface
    function get_int(this) result(ret)
      IMPORT I4B
      IMPORT ParticleEventType
      class(ParticleEventType), intent(in) :: this
      integer(I4B) :: ret
    end function get_int

    function get_str(this) result(ret)
      IMPORT ParticleEventType
      class(ParticleEventType), intent(in) :: this
      character(len=:), allocatable :: ret
    end function get_str
  end interface

contains

  function get_text(this) result(str)
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
  end function get_text

  subroutine log(this, iun)
    class(ParticleEventType), intent(inout) :: this
    integer(I4B), intent(in) :: iun

    if (iun >= 0) write (iun, '(*(G0))') this%get_text()
  end subroutine log

end module ParticleEventModule
