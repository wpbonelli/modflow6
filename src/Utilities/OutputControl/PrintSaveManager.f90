!> @brief Print/save manager module.
module PrintSaveManagerModule

  use KindModule, only: DP, I4B, LGP
  use ArrayHandlersModule, only: expandarray
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use InputOutputModule, only: urword
  use PeriodBlockTimingModule, only: PeriodBlockTimingType

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
    type(PeriodBlockTimingType) :: save_timing
    type(PeriodBlockTimingType) :: print_timing
    logical(LGP) :: save_detected
    logical(LGP) :: print_detected
  contains
    procedure :: init
    procedure :: rp
    procedure :: kstp_to_print
    procedure :: kstp_to_save
  end type PrintSaveManagerType

contains

  !> @brief (Re-)initialize the manager.
  subroutine init(this)
    class(PrintSaveManagerType) :: this !< this instance
    
    call this%save_timing%init()
    call this%print_timing%init()
    this%save_detected = .false.
    this%print_detected = .false.
  end subroutine init

  !> @ brief Read a line of input and prepare the manager.
  subroutine rp(this, linein, iout)
    ! -- dummy
    class(PrintSaveManagerType) :: this !< psm object
    character(len=*), intent(in) :: linein !< character line of information
    integer(I4B), intent(in) :: iout !< unit number of output file
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
    
    this%save_detected = ls
    this%print_detected = lp

    if (lp) then
      call this%print_timing%read(line(istop:))
      if (iout > 0) then
        if (this%print_timing%all) then
          write (iout, "(6x,a)") 'ALL TIME STEPS WILL BE PRINTED'
        else if (this%print_timing%first) then
          write (iout, "(6x,a)") 'THE FIRST TIME STEP WILL BE PRINTED'
        else if (this%print_timing%last) then
          write (iout, "(6x,a)") 'THE LAST TIME STEP WILL BE PRINTED'
        else if (size(this%print_timing%kstp) > 0) then
          write (iout, fmt_steps) 'PRINTED', this%print_timing%kstp
        else if (this%print_timing%freq > 0) then
          write (iout, fmt_freq) 'PRINTED', this%print_timing%freq
        end if
      end if
    else
      call this%save_timing%read(line(istop:))
    end if

    ! -- set the steps to print or save
    case ('FREQUENCY')
      call urword(line, lloc, istart, istop, 2, ival, rval, -1, 0)
      if (lp) this%ifreq_print = ival
      if (ls) this%ifreq_save = ival
      if (iout > 0) then
        if (lp) write (iout, fmt_freq) 'PRINTED', this%ifreq_print
        if (ls) write (iout, fmt_freq) 'SAVED', this%ifreq_save
      end if
    case default
      write (errmsg, '(2a)') &
        'Looking for ALL, STEPS, FIRST, LAST, OR FREQUENCY. Found: ', &
        trim(adjustl(line))
      call store_error(errmsg, terminate=.TRUE.)
    end select
    !
    ! -- return
    return
  end subroutine rp

  !> @ brief Determine if it is time to print the data
  !!
  !!  Determine if data should be printed based on kstp and endofperiod
  !!
  !<
  logical function kstp_to_print(this, kstp, endofperiod)
    ! -- dummy
    class(PrintSaveManagerType) :: this !< psm object
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< flag indicating end of stress period
    ! -- local
    integer(I4B) :: i, n
    !
    kstp_to_print = .false.
    if (this%print_all) kstp_to_print = .true.
    if (kstp == 1 .and. this%print_first) kstp_to_print = .true.
    if (endofperiod .and. this%print_last) kstp_to_print = .true.
    if (this%ifreq_print > 0) then
      if (mod(kstp, this%ifreq_print) == 0) kstp_to_print = .true.
    end if
    n = size(this%kstp_list_print)
    if (n > 0) then
      do i = 1, n
        if (kstp == this%kstp_list_print(i)) then
          kstp_to_print = .true.
          exit
        end if
      end do
    end if
    !
    ! -- Return
    return
  end function kstp_to_print

  !> @ brief Determine if it is time to save the data
  !!
  !!  Determine if data should be saved based on kstp and endofperiod
  !!
  !<
  logical function kstp_to_save(this, kstp, endofperiod)
    ! -- dummy
    class(PrintSaveManagerType) :: this !< psm object
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< flag indicating end of stress period
    ! -- local
    integer(I4B) :: i, n
    !
    kstp_to_save = .false.
    if (this%save_all) kstp_to_save = .true.
    if (kstp == 1 .and. this%save_first) kstp_to_save = .true.
    if (endofperiod .and. this%save_last) kstp_to_save = .true.
    if (this%ifreq_save > 0) then
      if (mod(kstp, this%ifreq_save) == 0) kstp_to_save = .true.
    end if
    n = size(this%kstp_list_save)
    if (n > 0) then
      do i = 1, n
        if (kstp == this%kstp_list_save(i)) then
          kstp_to_save = .true.
          exit
        end if
      end do
    end if
    !
    ! -- Return
    return
  end function kstp_to_save

end module PrintSaveManagerModule
