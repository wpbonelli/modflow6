module PrtFmiModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, LENAUXNAME, LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error, store_warning
  use SimVariablesModule, only: errmsg, warnmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use BudgetObjectModule, only: BudgetObjectType
  use CellModule, only: MAX_POLY_CELLS
  use CellDefnModule, only: CellDefnType

  implicit none
  private
  public :: PrtFmiType
  public :: fmi_cr

  character(len=LENPACKAGENAME) :: text = '    PRTFMI'

  type, extends(FlowModelInterfaceType) :: PrtFmiType

    integer(I4B), pointer :: iffmeth => null() !< internal iflowface method: 0 = stop, 1 = ignore
    real(DP), allocatable, public :: SourceFlows(:) ! cell source flows
    real(DP), allocatable, public :: SinkFlows(:) ! cell sink flows
    real(DP), allocatable, public :: StorageFlows(:) ! cell storage flows
    real(DP), allocatable, public :: BoundaryFlows(:) ! cell boundary flows
    integer(I4B), allocatable, public :: BoundaryFaces(:) ! bitmask of cell iflowface, 0-30: lateral, 31 (-2): bottom, 32 (-1): top

  contains

    procedure :: fmi_ad
    procedure :: fmi_df => prtfmi_df
    procedure :: allocate_scalars => prtfmi_allocate_scalars
    procedure :: source_options => prtfmi_options
    procedure, private :: accumulate_flows
    procedure :: mark_boundary_face
    procedure :: is_boundary_face
    procedure :: validate_boundary_faces

  end type PrtFmiType

contains

  !> @brief Create a new PrtFmi object
  subroutine fmi_cr(fmiobj, name_model, input_mempath, inunit, iout)
    ! dummy
    type(PrtFmiType), pointer :: fmiobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(inout) :: inunit
    integer(I4B), intent(in) :: iout

    ! Create the object
    allocate (fmiobj)

    ! create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text

    ! Allocate scalars
    call fmiobj%allocate_scalars()

    ! Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout

    ! Assign dependent variable label
    fmiobj%depvartype = 'TRACKS          '

  end subroutine fmi_cr

  !> @brief Time step advance
  subroutine fmi_ad(this)
    ! modules
    use ConstantsModule, only: DHDRY
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: n
    character(len=15) :: nodestr
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE')"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a)"

    ! Set flag to indicated that flows are being updated.  For the case where
    ! flows may be reused (only when flows are read from a file) then set
    ! the flag to zero to indicated that flows were not updated
    this%iflowsupdated = 1

    ! If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) &
      call this%advance_bfr()

    ! If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) &
      call this%advance_hfr()

    ! If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) &
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)

    ! Accumulate flows
    call this%accumulate_flows()

    ! Validate IFLOWFACE assignments
    call this%validate_boundary_faces()

    ! if flow cell is dry, then set this%ibound = 0
    do n = 1, this%dis%nodes

      ! Calculate the ibound-like array that has 0 if saturation
      ! is zero and 1 otherwise
      if (this%gwfsat(n) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if

      ! Check if active model cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(n) == DHDRY) then
          ! cell should be made inactive
          this%ibound(n) = 0
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtdry) trim(nodestr)
        end if
      end if

      ! Convert dry model cell to active if flow has rewet
      if (this%ibound(n) == 0) then
        if (this%gwfhead(n) /= DHDRY) then
          ! cell is now wet
          this%ibound(n) = 1
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtrewet) trim(nodestr)
        end if
      end if
    end do
  end subroutine fmi_ad

  !> @brief Define the flow model interface
  subroutine prtfmi_df(this, dis, idryinactive)
    ! modules
    use SimModule, only: store_error
    ! dummy
    class(PrtFmiType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: idryinactive
    ! formats
    character(len=*), parameter :: fmtfmi = &
      "(1x,/1x,'FMI -- FLOW MODEL INTERFACE, VERSION 2, 8/17/2023',            &
      &' INPUT READ FROM MEMPATH: ', A, //)"
    character(len=*), parameter :: fmtfmi0 = &
                    "(1x,/1x,'FMI -- FLOW MODEL INTERFACE,'&
                    &' VERSION 2, 8/17/2023')"

    ! Print a message identifying the FMI package.
    if (this%iout > 0) then
      if (this%inunit /= 0) then
        write (this%iout, fmtfmi) this%input_mempath
      else
        write (this%iout, fmtfmi0)
        if (this%flows_from_file) then
          write (this%iout, '(a)') '  FLOWS ARE ASSUMED TO BE ZERO.'
        else
          write (this%iout, '(a)') '  FLOWS PROVIDED BY A GWF MODEL IN THIS &
            &SIMULATION'
        end if
      end if
    end if

    ! Store pointers
    this%dis => dis

    ! Read fmi options
    if (this%inunit /= 0) then
      call this%source_options()
    end if

    ! Read packagedata options
    if (this%inunit /= 0 .and. this%flows_from_file) then
      call this%source_packagedata()
      call this%initialize_gwfterms_from_bfr()
    end if

    ! If GWF-Model exchange is active, setup flow terms
    if (.not. this%flows_from_file) then
      call this%initialize_gwfterms_from_gwfbndlist()
    end if

    ! Set flag that stops dry flows from being deactivated
    ! TODO: consider if relevant to PRT
    this%idryinactive = idryinactive

    ! Allocate arrays
    allocate (this%StorageFlows(this%dis%nodes))
    allocate (this%SourceFlows(this%dis%nodes))
    allocate (this%SinkFlows(this%dis%nodes))
    allocate (this%BoundaryFlows(this%dis%nodes * MAX_POLY_CELLS))
    allocate (this%BoundaryFaces(this%dis%nodes))

  end subroutine prtfmi_df

  !> @brief Allocate scalars
  subroutine prtfmi_allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtFmiType) :: this

    call this%FlowModelInterfaceType%allocate_scalars()
    call mem_allocate(this%iffmeth, "IFFMETH", this%memoryPath)
    this%iffmeth = 0
  end subroutine prtfmi_allocate_scalars

  !> @brief Accumulate flows
  subroutine accumulate_flows(this)
    implicit none
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: j, i, ip, ib
    integer(I4B) :: ioffset, iflowface, iauxiflowface
    real(DP) :: qbnd
    character(len=LENAUXNAME) :: auxname
    integer(I4B) :: naux

    this%StorageFlows = DZERO
    if (this%igwfstrgss /= 0) &
      this%StorageFlows = this%StorageFlows + &
                          this%gwfstrgss
    if (this%igwfstrgsy /= 0) &
      this%StorageFlows = this%StorageFlows + &
                          this%gwfstrgsy
    this%SourceFlows = DZERO
    this%SinkFlows = DZERO
    this%BoundaryFlows = DZERO
    this%BoundaryFaces = 0
    do ip = 1, this%nflowpack
      iauxiflowface = 0
      naux = this%gwfpackages(ip)%naux
      if (naux > 0) then
        do j = 1, naux
          auxname = this%gwfpackages(ip)%auxname(j)
          if (trim(adjustl(auxname)) == "IFLOWFACE") then
            iauxiflowface = j
            exit
          end if
        end do
      end if
      do ib = 1, this%gwfpackages(ip)%nbound
        i = this%gwfpackages(ip)%nodelist(ib)
        if (i <= 0) cycle
        if (this%ibound(i) <= 0) cycle
        qbnd = this%gwfpackages(ip)%get_flow(ib)
        ! todo, after initial release: default iflowface values for different packages
        iflowface = 0
        if (iauxiflowface > 0) then
          iflowface = NINT(this%gwfpackages(ip)%auxvar(iauxiflowface, ib))
          if (this%iffmeth /= 1 .and. iflowface > 0) &
            call this%mark_boundary_face(i, iflowface)
          ! maps bot -2 -> MAX_POLY_CELLS - 1, top -1 -> MAX_POLY_CELLS
          if (iflowface < 0) iflowface = iflowface + MAX_POLY_CELLS + 1
        end if
        if (this%iffmeth /= 1 .and. iflowface > 0) then
          ioffset = (i - 1) * MAX_POLY_CELLS
          this%BoundaryFlows(ioffset + iflowface) = &
            this%BoundaryFlows(ioffset + iflowface) + qbnd
        else if (qbnd .gt. DZERO) then
          this%SourceFlows(i) = this%SourceFlows(i) + qbnd
        else if (qbnd .lt. DZERO) then
          this%SinkFlows(i) = this%SinkFlows(i) + qbnd
        end if
      end do
    end do
  end subroutine accumulate_flows

  !> @brief Mark a face as a boundary face. Maps IFLOWFACE to bit position.
  subroutine mark_boundary_face(this, n, iflowface)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: n, iflowface
    integer(I4B) :: bit_pos

    if (n <= 0 .or. n > size(this%BoundaryFaces)) return
    if (iflowface == 0) return
    if (iflowface > 0) then
      bit_pos = iflowface
    else
      bit_pos = 32 + iflowface ! -1: 31 (top), -2: 30 (bottom)
    end if
    if (bit_pos < 1 .or. bit_pos > 32) then
      print *, 'Invalid IFLOWFACE: ', iflowface
      call pstop(1)
    end if
    this%BoundaryFaces(n) = ibset(this%BoundaryFaces(n), bit_pos - 1)
  end subroutine mark_boundary_face

  !> @brief Check if a face is marked as boundary.
  function is_boundary_face(this, n, iflowface) result(is_boundary)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: n, iflowface
    ! local
    logical :: is_boundary
    integer(I4B) :: bit_pos

    is_boundary = .false.
    if (n <= 0 .or. n > size(this%BoundaryFaces)) return
    if (iflowface == 0) return
    if (iflowface > 0) then
      bit_pos = iflowface
    else
      bit_pos = 32 + iflowface ! -1: 31 (top), -2: 30 (bottom)
    end if
    if (bit_pos < 1 .or. bit_pos > 32) then
      print *, 'Invalid IFLOWFACE: ', iflowface
      call pstop(1)
    end if
    is_boundary = btest(this%BoundaryFaces(n), bit_pos - 1)
  end function is_boundary_face

  !> @brief Validate IFLOWFACE assignments.
  subroutine validate_boundary_faces(this)
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: n, iflowface, face_count
    integer(I4B), parameter :: max_faces = 30
    integer(I4B) :: assigned_faces(max_faces)
    ! formats
    character(len=*), parameter :: fmt_warn = &
      "('WARNING: Cell ', i0, ' &
      &  IFLOWFACE will be ignored: ', 20i3)"

    do n = 1, this%dis%nodes
      if (this%BoundaryFaces(n) == 0) cycle
      face_count = 0
      assigned_faces = 0
      ! positive (lateral)
      do iflowface = 1, max_faces
        if (btest(this%BoundaryFaces(n), iflowface - 1)) then
          face_count = face_count + 1
          if (face_count <= max_faces) assigned_faces(face_count) = iflowface
        end if
      end do
      ! negative (top/bottom)
      if (btest(this%BoundaryFaces(n), 30)) then ! -2 -> bit 30
        face_count = face_count + 1
        if (face_count <= max_faces) assigned_faces(face_count) = -2
      end if
      if (btest(this%BoundaryFaces(n), 31)) then ! -1 -> bit 31
        face_count = face_count + 1
        if (face_count <= max_faces) assigned_faces(face_count) = -1
      end if
      ! no validation per se just warning that IFLOWFACE
      ! is ignored if INTERNAL_BOUNDARY_METHOD is IGNORE
      if (this%iffmeth == 1 .and. face_count > 1) then
        write (warnmsg, fmt_warn) n, assigned_faces(1:min(face_count, 20))
        call store_warning(trim(adjustl(warnmsg)))
      end if
    end do
  end subroutine validate_boundary_faces

  !> @ brief Source input options for package
  subroutine prtfmi_options(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    class(PrtFmiType) :: this
    ! local
    logical(LGP) :: found_ipakcb, found_iffmeth
    character(len=LENVARNAME), dimension(2) :: iff_method = &
      &[character(len=LENVARNAME) :: 'STOP', 'IGNORE']
    ! formats
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO AND FLOW IMBALANCE CORRECTION ACTIVE.')"

    ! source base class options
    call this%FlowModelInterfaceType%source_options()

    ! source package input
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found_ipakcb)
    call mem_set_value(this%iffmeth, 'IFFMETH', this%input_mempath, &
                       iff_method, found_iffmeth)

    write (this%iout, '(1x,a)') 'PROCESSING FMI OPTIONS'

    if (found_ipakcb) then
      this%ipakcb = -1
      write (this%iout, fmtisvflow)
    end if

    if (found_iffmeth) then
      if (this%iffmeth == 0) then
        write (errmsg, '(a)') 'Unsupported internal boundary method. &
          &INTERNAL_BOUNDARY_METHOD must be "STOP" or "IGNORE"'
        call store_error(errmsg)
      else
        ! adjust for method zero indexing
        this%iffmeth = this%iffmeth - 1
      end if
    end if

    write (this%iout, '(1x,a)') 'END OF FMI OPTIONS'
  end subroutine prtfmi_options

end module PrtFmiModule
