module MethodModelModule
  use KindModule, only: DP, I4B
  use MethodModule, only: MethodType, LEVEL_MODEL
  use CellModule, only: CellType
  use ParticleModule, only: ParticleType, TERM_BOUNDARY
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType

  private
  public :: MethodModelType
  public :: iface_to_iflowface
  public :: iflowface_to_iface

  type, abstract, extends(MethodType) :: MethodModelType
  contains
    procedure, public :: assess
    procedure, public :: get_level
    procedure, public :: check_internal_boundary
  end type MethodModelType

contains

  !> @brief Check particle reporting/termination status
  subroutine assess(this, particle, cell_defn, tmax)
    ! dummy
    class(MethodModelType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! noop
  end subroutine assess

  !> @brief Get the model method's level.
  function get_level(this) result(level)
    class(MethodModelType), intent(in) :: this
    integer(I4B) :: level
    level = LEVEL_MODEL
  end function get_level

  !> @brief Determine whether to stop at an internal assigned boundary face.
  !!
  !! This subroutine implements the logic for handling particles that reach
  !! internal boundary faces (faces assigned with IFLOWFACE). The behavior
  !! depends on the particle's INTERNAL_BOUNDARY_METHOD setting:
  !!
  !! - STOP: Terminate the particle at the face (default option)
  !! - IGNORE: Ignore boundary face assign, mimicking MODPATH 7
  !<
  subroutine check_internal_boundary(this, particle, cell)
    class(MethodModelType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(*), intent(inout) :: cell
    ! local
    integer(I4B) :: iflowface

    select type (cell)
    class is (CellType)
      iflowface = iface_to_iflowface(particle%iboundary(2), cell%defn%npolyverts)
      if (.not. this%fmi%is_boundary_face(cell%defn%icell, iflowface)) return
      select case (this%fmi%iffmeth)
      case (0) ! STOP
        call this%terminate(particle, status=TERM_BOUNDARY)
      case default
        ! no action needed for other methods, they are handled
        ! by differences in FMI's accounting of boundary flows
      end select
    end select
  end subroutine check_internal_boundary

  !> @brief Convert IFLOWFACE value to internal face index
  !! IFLOWFACE: 1-N = lateral faces, -1 = top, -2 = bottom
  !! Internal: 1-N = lateral faces, N+3 = top, N+2 = bottom
  function iflowface_to_iface(iflowface, npolyverts) result(iface)
    integer(I4B), intent(in) :: iflowface !< IFLOWFACE value
    integer(I4B), intent(in) :: npolyverts !< number of polygon vertices/lateral faces
    integer(I4B) :: iface

    if (iflowface == -1) then
      iface = npolyverts + 3 ! Top face
    else if (iflowface == -2) then
      iface = npolyverts + 2 ! Bottom face
    else if (iflowface >= 1 .and. iflowface <= npolyverts) then
      iface = iflowface ! Lateral faces (direct mapping)
    else
      iface = 0 ! Invalid IFLOWFACE
    end if

  end function iflowface_to_iface

  !> @brief Convert internal face index to IFLOWFACE value
  !! Internal: 1-N = lateral faces, N+2 = bottom, N+3 = top
  !! IFLOWFACE: 1-N = lateral faces, -2 = bottom, -1 = top
  function iface_to_iflowface(iface, npolyverts) result(iflowface)
    integer(I4B), intent(in) :: iface !< internal face index
    integer(I4B), intent(in) :: npolyverts !< number of polygon vertices/lateral faces
    integer(I4B) :: iflowface

    if (iface == npolyverts + 2) then
      iflowface = -2 ! Bottom face
    else if (iface == npolyverts + 3) then
      iflowface = -1 ! Top face
    else if (iface >= 1 .and. iface <= npolyverts) then
      iflowface = iface ! Lateral faces (direct mapping)
    else
      iflowface = 0 ! Invalid face index
    end if

  end function iface_to_iflowface

end module MethodModelModule
