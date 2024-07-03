!> @brief Print/save manager module.
module PrintSaveManagerModule

  use KindModule, only: DP, I4B, LGP
  use ArrayHandlersModule, only: expandarray
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use InputOutputModule, only: urword
  use TimeStepSelectModule, only: TimeStepSelectType

  implicit none
  private
  public :: PrintSaveManagerType
  public :: create_psm

  !> @brief Print/save manager type.
  !!
  !! Stores user settings as configured in an input file's period block
  !! and determines whether data should be printed to a list (log) file
  !! or saved to disk.
  !!
  !! The object should be initiated with the init() procedure.
  !!
  !! The rp() procedure will read a character string and configure the
  !! manager, where the character string may be of the following form:
  !!
  !!   PRINT ALL
  !!   PRINT STEPS 1 4 5 6
  !!   PRINT FIRST
  !!   PRINT LAST
  !!   PRINT FREQUENCY 4
  !!   SAVE ALL
  !!   SAVE STEPS 1 4 5 6
  !!   SAVE FIRST
  !!   SAVE LAST
  !!   SAVE FREQUENCY 4
  !!
  !! The should_print() and should_save() functions indicate whether
  !! to save or print during the current time step.
  !<
  type :: PrintSaveManagerType
    private
    type(TimeStepSelectType), pointer, public :: save_steps
    type(TimeStepSelectType), pointer, public :: print_steps
  contains
    procedure :: deallocate
    procedure :: init
    procedure :: rp
    procedure :: should_print
    procedure :: should_save
  end type PrintSaveManagerType

contains

  !> @brief Initialize or clear the print/save manager.
  function create_psm() result(psm)
    type(PrintSaveManagerType), pointer :: psm !< the print/save manager

    allocate (psm)
    allocate (psm%save_steps)
    allocate (psm%print_steps)
    call psm%save_steps%init()
    call psm%print_steps%init()
  end function create_psm

  subroutine init(this)
    class(PrintSaveManagerType) :: this
    if (associated(this%save_steps)) deallocate(this%save_steps)
    if (associated(this%print_steps)) deallocate(this%print_steps)
    allocate (this%save_steps)
    allocate (this%print_steps)
  end subroutine init

  subroutine deallocate(this)
    class(PrintSaveManagerType) :: this !< this instance
    deallocate (this%save_steps)
    deallocate (this%print_steps)
  end subroutine deallocate

  !> @ brief Read a line of input and prepare the manager.
  subroutine rp(this, linein, iout)
    ! -- dummy
    class(PrintSaveManagerType) :: this !< this instance
    character(len=*), intent(in) :: linein !< input line
    integer(I4B), intent(in) :: iout !< output file unit
    ! -- local
    character(len=len(linein)) :: line
    logical(LGP) lp, ls
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval
    
    line(:) = linein(:)
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    
    lp = .false.
    ls = .false.
    select case (line(istart:istop))
    case ('PRINT')
      lp = .true.
    case ('SAVE')
      ls = .true.
    case default
      write (errmsg, '(2a)') &
        'Looking for PRINT or SAVE. Found:', trim(adjustl(line))
      call store_error(errmsg, terminate=.TRUE.)
    end select

    if (lp) then
      call this%print_steps%read(line(istop + 2:))
      if (iout > 0) &
        call this%print_steps%log(iout, verb="PRINTED")
    else
      call this%save_steps%read(line(istop + 2:))
      if (iout > 0) &
        call this%save_steps%log(iout, verb="SAVED")
    end if
  end subroutine rp

  !> @ brief Determine if printing is enabled for this time step.
  logical function should_print(this, kstp, endofperiod)
    class(PrintSaveManagerType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< whether last step of stress period
    
    should_print = this%print_steps%is_selected(kstp, endofperiod)
  end function should_print

  !> @ brief Determine if saving is enabled for this time step.
  logical function should_save(this, kstp, endofperiod)
    class(PrintSaveManagerType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< whether last step of stress period
    
    should_save = this%save_steps%is_selected(kstp, endofperiod)
  end function should_save

end module PrintSaveManagerModule
