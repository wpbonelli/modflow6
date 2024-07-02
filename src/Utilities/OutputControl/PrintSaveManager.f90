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
  !! The time_to_print() and time_to_save() functions indicate whether
  !! to save or print during the current time step.
  !<
  type :: PrintSaveManagerType
    private
    type(TimeStepSelectType), pointer, public :: save_steps
    type(TimeStepSelectType), pointer, public :: print_steps
  contains
    procedure :: init
    procedure :: rp
    procedure :: should_print
    procedure :: should_save
  end type PrintSaveManagerType

contains

  !> @brief Initialize or clear the print/save manager.
  subroutine init(this)
    class(PrintSaveManagerType) :: this !< this instance
    call this%save_steps%init()
    call this%print_steps%init()
  end subroutine init

  !> @ brief Read a line of input and prepare the manager.
  subroutine rp(this, linein, iout)
    ! -- dummy
    class(PrintSaveManagerType) :: this !< this instance
    character(len=*), intent(in) :: linein !< input line
    integer(I4B), intent(in) :: iout !< output file unit
    ! -- local
    character(len=len(linein)) :: line
    logical lp, ls
    integer(I4B) :: n
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval
    ! -- formats
    character(len=*), parameter :: fmt_steps = &
      &"(6x,'THE FOLLOWING STEPS WILL BE ',A,': ',50(I0,' '))"
    character(len=*), parameter :: fmt_freq = &
      &"(6x,'THE FOLLOWING FREQUENCY WILL BE ',A,': ',I0)"
    
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
      call this%print_steps%read(line(istop:))
      if (iout > 0) then
        if (this%print_steps%all) then
          write (iout, "(6x,a)") 'ALL TIME STEPS WILL BE PRINTED'
        else if (this%print_steps%first) then
          write (iout, "(6x,a)") 'THE FIRST TIME STEP WILL BE PRINTED'
        else if (this%print_steps%last) then
          write (iout, "(6x,a)") 'THE LAST TIME STEP WILL BE PRINTED'
        else if (size(this%print_steps%steps) > 0) then
          write (iout, fmt_steps) 'PRINTED', this%print_steps%steps
        else if (this%print_steps%freq > 0) then
          write (iout, fmt_freq) 'PRINTED', this%print_steps%freq
        end if
      end if
    else
      call this%save_steps%read(line(istop:))
      if (iout > 0) then
        if (this%save_steps%all) then
          write (iout, "(6x,a)") 'ALL TIME STEPS WILL BE SAVED'
        else if (this%save_steps%first) then
          write (iout, "(6x,a)") 'THE FIRST TIME STEP WILL BE SAVED'
        else if (this%save_steps%last) then
          write (iout, "(6x,a)") 'THE LAST TIME STEP WILL BE SAVED'
        else if (size(this%save_steps%steps) > 0) then
          write (iout, fmt_steps) 'SAVED', this%save_steps%steps
        else if (this%save_steps%freq > 0) then
          write (iout, fmt_freq) 'SAVED', this%save_steps%freq
        end if
      end if
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
