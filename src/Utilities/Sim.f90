!> @brief Simulation utilities.
!!
!! This module contains simulation-scoped utilities for issuing messages
!! conditional on logging level, stopping the simulation, checking whether
!! the simulation has converged, and other functions. This module does not
!! depend on models, exchanges, solutions, or other simulation components.
!<
module SimModule

  use KindModule, only: DP, I4B, LGP
  use FileUtilModule, only: get_filename, close_all_files
  use ConstantsModule, only: MAXCHARLEN, LINELENGTH, LENHUGELINE &
                             DONE, &
                             IUSTART, IULAST, &
                             VSUMMARY, VALL, VDEBUG
  use MessageModule, only: Messages, write_message
  use SimVariablesModule, only: istdout, iout, isim_level, ireturnerr, &
                                iforcestop, warnmsg

  implicit none

  private
  public :: initial_message
  public :: final_message
  public :: sim_message
  public :: deprecation_warning
  public :: set_max_errors
  public :: store_error
  public :: store_error_unit
  public :: store_error_filename
  public :: store_note
  public :: store_warning
  public :: count_errors
  public :: count_notes
  public :: count_warnings
  public :: converge_reset
  public :: converge_check
  public :: ustop

  type(Messages) :: sim_errors
  type(Messages) :: sim_uniterrors
  type(Messages) :: sim_warnings
  type(Messages) :: sim_notes

contains

  !> @brief Return the number of stored errors.
  function count_errors() result(ncount)
    integer(I4B) :: ncount
    ncount = sim_errors%get_count()
  end function count_errors

  !> @brief Return the number of stored warnings.
  function count_warnings() result(ncount)
    integer(I4B) :: ncount
    ncount = sim_warnings%get_count()
  end function count_warnings

  !> @brief Return the number of stored notes.
  function count_notes() result(ncount)
    integer(I4B) :: ncount
    ncount = sim_notes%get_count()
  end function count_notes

  !> @brief Set the maximum number of errors to be stored in a simulation.
  subroutine set_max_errors(imax)
    integer(I4B), intent(in) :: imax !< maximum number of error messages that will be stored
    call sim_errors%set_max(imax)
  end subroutine set_max_errors

  !> @brief Store an error message for issuance at the end of the simulation.
  !!
  !! The terminate option controls whether to terminate the simulation.
  !! By default (unlike store_error_unit and store_error_filename), the
  !! simulation is not terminated.
  !<
  subroutine store_error(msg, terminate)
    ! -- dummy
    character(len=*), intent(in) :: msg !< error message
    logical, optional, intent(in) :: terminate !< whether to terminate
    ! -- local
    logical :: lterminate
    
    ! -- process optional variables
    if (present(terminate)) then
      lterminate = terminate
    else
      lterminate = .FALSE.
    end if
    
    ! -- store error
    call sim_errors%store(msg)
    
    ! -- terminate the simulation
    if (lterminate) call ustop()
  end subroutine store_error

  !> @brief Store an error associated with a file unit.
  !!
  !! The unit number is converted to filename for more readable errors.
  !! The terminate option controls whether to terminate the simulation.
  !! By default, the simulation is terminated.
  !<
  subroutine store_error_unit(iunit, terminate)
    ! -- dummy
    integer(I4B), intent(in) :: iunit !< open file unit number
    logical, optional, intent(in) :: terminate !< whether to terminate
    ! -- local
    logical :: lterminate
    character(len=LINELENGTH) :: fname
    character(len=LINELENGTH) :: errmsg
    
    ! -- process optional variables
    if (present(terminate)) then
      lterminate = terminate
    else
      lterminate = .TRUE.
    end if
    
    ! -- store error unit
    inquire (unit=iunit, name=fname)
    write (errmsg, '(3a)') &
      "Error occurred while reading file '", trim(adjustl(fname)), "'"
    call sim_uniterrors%store(errmsg)
    
    ! -- terminate the simulation
    if (lterminate) call ustop()
  end subroutine store_error_unit

  !> @brief Store an error associated with a filename.
  !!
  !! The terminate option controls whether to terminate the simulation.
  !! By default, the simulation is terminated.
  !<
  subroutine store_error_filename(filename, terminate)
    ! -- dummy
    character(len=*), intent(in) :: filename !< erroring file name
    logical, optional, intent(in) :: terminate !< whether to terminate
    ! -- local
    logical :: lterminate
    character(len=LINELENGTH) :: errmsg
    
    ! -- process optional variables
    if (present(terminate)) then
      lterminate = terminate
    else
      lterminate = .TRUE.
    end if
    
    ! -- store error unit
    write (errmsg, '(3a)') &
      "ERROR OCCURRED WHILE READING FILE '", trim(adjustl(filename)), "'"
    call sim_uniterrors%store(errmsg)
    
    ! -- terminate the simulation
    if (lterminate) call ustop()
  end subroutine store_error_filename

  !> @brief Store a warning for issuance at the end of the simulation.
  !! If the message matches an optional substring, it is not stored.
  subroutine store_warning(msg, substring)
    character(len=*), intent(in) :: msg !< warning message
    character(len=*), intent(in), optional :: substring !< prevent duplicates
    
    if (present(substring)) then
      call sim_warnings%store(msg, substring)
    else
      call sim_warnings%store(msg)
    end if
  end subroutine store_warning

  !> @brief Store a deprecation warning for issuance at the end of the simulation. 
  subroutine deprecation_warning(cblock, cvar, cver, endmsg, iunit)
    ! -- dummy variables
    character(len=*), intent(in) :: cblock !< block name
    character(len=*), intent(in) :: cvar !< variable name
    character(len=*), intent(in) :: cver !< version when variable was deprecated
    character(len=*), intent(in), optional :: endmsg !< optional user defined message to append
                                                     !! at the end of the deprecation warning
    integer(I4B), intent(in), optional :: iunit !< optional input file unit number with
                                                !! the deprecated variable
    ! -- local variables
    character(len=MAXCHARLEN) :: message
    character(len=LINELENGTH) :: fname
    
    ! -- build message
    write (message, '(a)') &
      trim(cblock)//" BLOCK VARIABLE '"//trim(cvar)//"'"
    if (present(iunit)) then
      call get_filename(iunit, fname)
      write (message, '(a,1x,3a)') &
        trim(message), "IN FILE '", trim(fname), "'"
    end if
    write (message, '(a)') &
      trim(message)//' WAS DEPRECATED IN VERSION '//trim(cver)//'.'
    if (present(endmsg)) then
      write (message, '(a,1x,2a)') trim(message), trim(endmsg), '.'
    end if
    
    ! -- store warning
    call sim_warnings%store(message)
  end subroutine deprecation_warning

  !> @brief Store a note for issuance at the end of the simulation.
  subroutine store_note(note)
    character(len=*), intent(in) :: note !< note
    call sim_notes%store(note)
  end subroutine store_note

  !> @brief Stop the simulation with an error, optionally issuing a message.
  subroutine ustop(stopmess, ioutlocal)
    character, optional, intent(in) :: stopmess * (*) !< optional message to print before stop
    integer(I4B), optional, intent(in) :: ioutlocal !< optional file to write final message to
    
    call print_final_message(stopmess, ioutlocal)
    call exit(ireturnerr)
  end subroutine ustop

  !> @brief Issue stored messages and optional final message, and close open files.
  subroutine print_final_message(stopmess, ioutlocal)
    ! -- dummy
    character, optional, intent(in) :: stopmess * (*) !< optional message to show before stop
    integer(I4B), optional, intent(in) :: ioutlocal !< optional file to write final message to
    ! -- local
    character(len=*), parameter :: fmt = '(1x,a)'
    character(len=*), parameter :: msg = 'Stopping due to error(s)'
    
    ! -- print stored messages
    call sim_notes%print_message('NOTES:', 'note(s)', &
                                 iunit=iout, level=VALL)
    call sim_warnings%print_message('WARNING REPORT:', 'warning(s)', &
                                    iunit=iout, level=VALL)
    call sim_errors%print_message('ERROR REPORT:', 'error(s)', iunit=iout)
    call sim_uniterrors%print_message('UNIT ERROR REPORT:', &
                                      'file unit error(s)', iunit=iout)
    
    ! -- write a stop message, if one is passed
    if (present(stopmess)) then
      if (stopmess .ne. ' ') then
        call sim_message(stopmess, fmt=fmt, iunit=iout)
        call sim_message(stopmess, fmt=fmt)
        if (present(ioutlocal)) then
          if (ioutlocal > 0 .and. ioutlocal /= iout) then
            write (ioutlocal, fmt) trim(stopmess)
            close (ioutlocal)
          end if
        end if
      end if
    end if
    
    ! -- flush buffered output
    flush (istdout)
    
    ! -- determine if an error condition has occurred
    if (sim_errors%get_count() > 0) then
      ireturnerr = 2
      if (present(ioutlocal)) then
        if (ioutlocal > 0 .and. ioutlocal /= iout) write (ioutlocal, fmt) msg
      end if
    end if
    
    ! -- close all open files
    call close_all_files()
    
  end subroutine print_final_message

  !> @brief Reset the simulation convergence flag.
  subroutine converge_reset()
    use SimVariablesModule, only: isimcnvg
    isimcnvg = 1
  end subroutine converge_reset

  !> @brief Check simulation convergence between time steps.
  !!
  !! By default, this routine indicates whether the simulation has converged.
  !! If CONTINUE is set in the simulation control file, the convergence flag
  !! and this routine's output will both indicate convergence, even if there
  !! is no convergence yet. The non-convergence counter is always incremented
  !! if the previous time step did not converge.
  !<
  subroutine converge_check(hasConverged)
    ! -- modules
    use SimVariablesModule, only: isimcnvg, numnoconverge, isimcontinue
    ! -- dummy
    logical, intent(inout) :: hasConverged !< whether simulation is considered converged
    ! -- format
    character(len=*), parameter :: fmtfail = &
      "(1x, 'Simulation convergence failure.', &
      &' Simulation will terminate after output and deallocation.')"
    
    ! -- initialize flag
    hasConverged = .true.
    
    ! -- increment convergence failure count if needed
    if (isimcnvg == 0) then
      numnoconverge = numnoconverge + 1
    end if
    
    ! -- continue if 'CONTINUE' specified in simulation control file
    if (isimcontinue == 1) then
      if (isimcnvg == 0) then
        isimcnvg = 1
      end if
    end if
    
    ! -- save simulation failure message
    if (isimcnvg == 0) then
      call sim_message('', fmt=fmtfail, iunit=iout)
      hasConverged = .false.
    end if
  end subroutine converge_check

  !> @brief Initialize message storage and issue initial simulation message.
  subroutine initial_message()
    ! -- modules
    use VersionModule, only: write_listfile_header
    use SimVariablesModule, only: simulation_mode
    
    ! -- initialize message storage
    call sim_errors%init()
    call sim_uniterrors%init()
    call sim_warnings%init()
    call sim_notes%init()
    
    ! -- show simulation banner
    call write_listfile_header(istdout, write_kind_info=.false., &
                               write_sys_command=.false.)
    
    ! -- indicate whether we are running in parallel
    if (simulation_mode == 'PARALLEL') &
      call sim_message('(MODFLOW runs in '//trim(simulation_mode)//' mode)', &
                       skipafter=1)
  end subroutine initial_message

  !> @brief Issue a final simulation message and terminate if necessary. 
  subroutine final_message()
    ! -- modules
    use SimVariablesModule, only: isimcnvg, numnoconverge, ireturnerr, &
                                  isimcontinue
    ! -- formats
    character(len=*), parameter :: fmtnocnvg = &
      &"(1x, 'Simulation convergence failure occurred ', i0, ' time(s).')"
    
    ! -- non-convergence warning if any timesteps failed to converge
    if (numnoconverge > 0) then
      write (warnmsg, fmtnocnvg) numnoconverge
      if (isimcontinue == 0) then
        call sim_errors%store(warnmsg)
      else
        call sim_warnings%store(warnmsg)
      end if
    end if
    
    ! -- write final message
    if (isimcnvg == 0) then
      call print_final_message('Premature termination of simulation.', iout)
    else
      call print_final_message('Normal termination of simulation.', iout)
    end if
    
    ! -- if convergence failure and CONTINUE isn't set, exit code 1
    if (isimcnvg == 0 .and. isimcontinue == 0) &
      ireturnerr = 1
    
    ! -- destroy messages
    call sim_errors%deallocate()
    call sim_uniterrors%deallocate()
    call sim_warnings%deallocate()
    call sim_notes%deallocate()
    
    ! -- return or halt
    if (iforcestop == 1) &
      call exit(ireturnerr)

  end subroutine final_message

  !> @brief Issue a message conditional on the simulation's logging level.
  subroutine sim_message(text, iunit, fmt, level, &
                         skipbefore, skipafter, advance)
    ! -- dummy
    character(len=*), intent(in) :: text !< message to write to iunit
    integer(I4B), intent(in), optional :: iunit !< optional file unit to write the message to (default=stdout)
    character(len=*), intent(in), optional :: fmt !< optional format to write the message (default='(a)')
    integer(I4B), intent(in), optional :: level !< optional level for the message (default=summary)
    integer(I4B), intent(in), optional :: skipbefore !< optional number of empty lines before message (default=0)
    integer(I4B), intent(in), optional :: skipafter !< optional number of empty lines after message (default=0)
    logical(LGP), intent(in), optional :: advance !< optional boolean indicating if advancing output (default is .TRUE.)
    ! -- local
    integer(I4B) :: iu
    character(len=*), parameter :: stdfmt = '(a)'

    if (present(inunit)) then
      iu = inunit
    else
      iu = istdout
    end if

    if (ilevel <= isim_level) &
      write_message(iu, text, stdfmt, skipbefore, skipafter, advance)
  end subroutine sim_message

end module SimModule
