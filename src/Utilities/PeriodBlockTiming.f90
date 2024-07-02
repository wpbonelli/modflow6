!> @brief Period block timing module.
module PeriodBlockTimingModule

  use KindModule, only: DP, I4B, LGP
  use ArrayHandlersModule, only: expandarray
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use InputOutputModule, only: urword

  implicit none
  private
  public :: PeriodBlockTimingType

  !> @brief Period block timing type.
  !!
  !! Represents a selection of time steps as configured in an input file's
  !! period block settings. The object should be initiated with the init()
  !! procedure. The read() procedure accepts a character string of form:
  !!
  !!   ALL
  !!   STEPS 1 4 5 6
  !!   FIRST
  !!   LAST
  !!   FREQUENCY 4
  !!
  !! The is_active(kstp) function indicates whether the given time step is
  !! selected. This function accepts an optional argument, indicating that
  !! the time step is the last in the stress period.
  !<
  type :: PeriodBlockTimingType
    logical(LGP) :: all
    logical(LGP) :: first
    logical(LGP) :: last
    integer(I4B) :: freq
    integer(I4B), allocatable :: kstp(:)
  contains
    procedure :: init
    procedure :: read
    procedure :: is_active
  end type PeriodBlockTimingType

contains

  !> @ brief Initialize the period block timing selection object.
  subroutine init(this)
    class(PeriodBlockTimingType) :: this !< this instance
    
    if (allocated(this%kstp)) deallocate (this%kstp)
    allocate (this%kstp(0))
    this%freq = 0
    this%first = .false.
    this%last = .false.
    this%all = .false.
  end subroutine init

  subroutine read(this, line)
    class(PeriodBlockTimingType) :: this !< this instance
    character(len=*), intent(in) :: line !< character line

    character(len=len(line)) :: l
    integer(I4B) :: n, lloc, istart, istop, ival
    real(DP) :: rval

    l(:) = line(:)
    lloc = 1
    
    call urword(l, lloc, istart, istop, 1, ival, rval, 0, 0)
    select case (l(istart:istop))
    case ('ALL')
      this%all = .true.
    case ('STEPS')
      listsearch: do
        call urword(l, lloc, istart, istop, 2, ival, rval, -1, 0)
        if (ival > 0) then
          n = size(this%kstp)
          call expandarray(this%kstp)
          this%kstp(n + 1) = ival
          cycle listsearch
        end if
        exit listsearch
      end do listsearch
    case ('FREQUENCY')
      call urword(l, lloc, istart, istop, 2, ival, rval, -1, 0)
      this%freq = ival
    case ('FIRST')
      this%first = .true.
    case ('LAST')
      this%last = .true.
    case default
      write (errmsg, '(2a)') &
        'Looking for ALL, STEPS, FIRST, LAST, OR FREQUENCY. Found: ', &
        trim(adjustl(line))
      call store_error(errmsg, terminate=.TRUE.)
    end select
  end subroutine read

  logical function is_active(this, kstp, endofperiod)
    class(PeriodBlockTimingType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in), optional :: endofperiod !< whether last step of stress period
    integer(I4B) :: i, n
    logical(LGP) :: lend
    
    if (present(endofperiod)) then
      lend = endofperiod
    else
      lend = .false.
    end if

    is_active = .false.
    if (this%all) is_active = .true.
    if (kstp == 1 .and. this%first) is_active = .true.
    if (lend .and. this%last) is_active = .true.
    if (this%freq > 0) then
      if (mod(kstp, this%freq) == 0) is_active = .true.
    end if
    n = size(this%kstp)
    if (n > 0) then
      do i = 1, n
        if (kstp == this%kstp(i)) then
          is_active = .true.
          exit
        end if
      end do
    end if
  end function is_active

end module PeriodBlockTimingModule