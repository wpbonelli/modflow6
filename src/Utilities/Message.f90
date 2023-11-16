!> @brief Generic utilities to store and issue messages to an output unit.
module MessageModule

  use KindModule, only: LGP, I4B, DP
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, DONE, &
                             VSUMMARY, LENHUGELINE
  use ArrayHandlersModule, only: ExpandArray

  implicit none

  public :: Messages
  public :: write_message
  public :: write_message_counter
  public :: write_message_centered

  !> @brief Container for related messages sharing a name and title.
  !!
  !! A maximum number of messages may be configured. Message storage
  !! arrays are dynamically reallocated up to the specified capacity.
  !<
  type :: Messages

    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: name
    character(len=MAXCHARLEN), allocatable, dimension(:) :: messages
    integer(I4B) :: num_messages !< number of messages stored
    integer(I4B) :: max_messages !< max number of messages to allow
    integer(I4B) :: max_exceeded !< excess requested beyond max allowed
    integer(I4B) :: expand_count !< amount to expand on dynamic reallocation

  contains

    procedure :: deallocate
    procedure :: get_count
    procedure :: init
    procedure :: set_max
    procedure :: store
    procedure :: write_all
    
  end type Messages

contains

  !> @brief Initialize message storage. Array allocation occurs on-the-fly.
  subroutine init(this, imax)
    class(Messages) :: this
    integer(I4B), intent(in), optional :: imax !< maximum message count
    this%num_messages = 0
    this%max_exceeded = 0
    this%expand_count = 100
    if (present(imax)) then
      this%max_messages = max(imax, 0)
    else
      this%max_messages = 1000
    end if
  end subroutine init

  !> @brief Return the number of messages currently stored.
  function get_count(this) result(nmessage)
    class(Messages) :: this
    integer(I4B) :: nmessage
    if (allocated(this%messages)) then
      nmessage = this%num_messages
    else
      nmessage = 0
    end if
  end function get_count

  !> @brief Set the maximum number of messages.
  subroutine set_max(this, imax)
    class(Messages) :: this
    integer(I4B), intent(in) :: imax
    this%max_messages = imax
  end subroutine set_max

  !> @ brief Deallocate message storage if needed.
  subroutine deallocate(this)
    class(Messages) :: this
    if (allocated(this%messages)) deallocate (this%messages)
  end subroutine deallocate

  !> @brief Store a message. Array allocation occurs on-the-fly.
  !! If the message matches a given substring, it is not stored.
  subroutine store(this, text, substring)
    ! -- dummy
    class(Messages) :: this
    character(len=*), intent(in) :: text
    character(len=*), intent(in), optional :: substring
    ! -- local
    integer(I4B) :: i

    ! -- if any message contains substring, don't duplicate
    if (present(substring)) then
      do i = 1, this%num_messages
        if (index(this%messages(i), substring) > 0) return
      end do
    end if

    ! -- if at capacity, increment excess counter and abort
    i = this%num_messages + 1 ! new message count
    if (i > this%max_messages) then
      this%max_exceeded = this%max_exceeded + 1
      return
    end if

    ! -- expand messages array if needed
    if (allocated(this%messages)) then
      if (this%num_messages == size(this%messages)) then
        call ExpandArray(this%messages, increment=this%expand_count)
        this%expand_count = int(this%expand_count * 1.1)
      end if
    end if

    ! store message and increment message counter
    this%messages(i) = text
    this%num_messages = i
    
  end subroutine store

  !> @brief Write all stored messages.
  subroutine write_all(this, title, name, iunit)
    ! -- dummy
    class(Messages) :: this !< MessageType object
    character(len=*), intent(in) :: title !< message title
    character(len=*), intent(in) :: name !< message name
    integer(I4B), intent(in), optional :: iunit !< optional file unit to save messages to
    ! -- local
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: cerr
    integer(I4B) :: iu
    integer(I4B) :: ilevel
    integer(I4B) :: i
    integer(I4B) :: isize
    integer(I4B) :: iwidth
    ! -- formats
    character(len=*), parameter :: stdfmt = "(/,A,/)"
    
    ! -- write the title and all message entries
    if (allocated(this%messages)) then
      isize = this%num_messages
      if (isize > 0) then
        
        ! -- calculate max width of prepended string for counter
        write (cerr, '(i0)') isize
        iwidth = len_trim(cerr) + 1
        
        ! -- write title for messages
        call write_message(iunit=iunit, text=title, fmt=stdfmt)
        
        ! -- write each message
        do i = 1, isize
          call write_message_counter(iunit=iunit, text=this%messages(i), icount=i, iwidth=iwidth)
        end do
        
        ! -- write the number of additional messages
        if (this%max_exceeded > 0) then
          write (errmsg, '(i0,3(1x,a))') &
            this%max_exceeded, 'additional', trim(name), &
            'detected but not printed.'
          call write_message(iunit=iunit, text=trim(errmsg), fmt='(/,1x,a)')
        end if
      end if
    end if
  end subroutine write_all

  ! standalone procedures

  !> @brief Configurable routine to write a message to an output unit.
  subroutine write_message(iunit, text, fmt, &
                           skipbefore, skipafter, advance)
    ! -- dummy
    integer(I4B), intent(in) :: iunit !< output unit to write the message to
    character(len=*), intent(in) :: text !< message to write
    character(len=*), intent(in), optional :: fmt !< format to write the message (default='(a)')
    integer(I4B), intent(in), optional :: skipbefore !< number of empty lines before message (default=0)
    integer(I4B), intent(in), optional :: skipafter !< number of empty lines after message (default=0)
    logical(LGP), intent(in), optional :: advance !< whether to use advancing output (default is .TRUE.)
    ! -- local
    character(len=3) :: cadvance
    integer(I4B) :: i
    integer(I4B) :: ilen
    integer(I4B) :: iu
    character(len=LENHUGELINE) :: simfmt
    character(len=*), parameter :: stdfmt = '(a)'
    character(len=*), parameter :: emptyfmt = '()'

    ! -- get message length
    ilen = len_trim(text)

    ! -- process optional arguments
    if (present(fmt)) then
      simfmt = fmt
    else
      if (ilen > 0) then
        simfmt = stdfmt
      else
        simfmt = emptyfmt
      end if
    end if
    if (present(advance)) then
      if (advance) then
        cadvance = 'YES'
      else
        cadvance = 'NO'
      end if
    else
      cadvance = 'YES'
    end if

    ! -- write empty line before message, if enabled
    if (present(skipbefore)) then
      do i = 1, skipbefore
        write (iu, *)
      end do
    end if

    ! -- write message if it isn't empty
    if (ilen > 0) then
      write (iu, trim(simfmt), advance=cadvance) text(1:ilen)
    else
      write (iu, trim(simfmt), advance=cadvance)
    end if

    ! -- write empty line after message, if enabled
    if (present(skipafter)) then
      do i = 1, skipafter
        write (iu, *)
      end do
    end if
  end subroutine write_message

  !> @brief Write a message, splitting across lines as needed.
  !!
  !! Subroutine that formats and writes a single message that
  !! may exceeed 78 characters in length. Messages longer than
  !! 78 characters are written across multiple lines. When a
  !! counter is passed in subsequent lines are indented.
  !<
  subroutine write_message_counter(iunit, text, icount, iwidth, &
                                  skipbefore, skipafter)
    ! -- dummy
    integer(I4B), intent(in) :: iunit !< the unit number to which the message is written
    character(len=*), intent(in) :: text !< message to be written
    integer(I4B), intent(in), optional :: icount !< counter to prepended to the message
    integer(I4B), intent(in), optional :: iwidth !< maximum width of the prepended counter
    integer(I4B), intent(in), optional :: skipbefore !< optional number of empty lines before message (default=0)
    integer(I4B), intent(in), optional :: skipafter !< optional number of empty lines after message (default=0)
    ! -- local
    integer(I4B), parameter :: len_line = 78
    character(len=LENHUGELINE) :: amessage
    character(len=len_line) :: line
    character(len=16) :: cfmt
    character(len=10) :: counter
    character(len=5) :: fmt_first
    character(len=20) :: fmt_cont
    logical(LGP) :: include_counter
    integer(I4B) :: isb
    integer(I4B) :: isa
    integer(I4B) :: jend
    integer(I4B) :: len_str1
    integer(I4B) :: len_str2
    integer(I4B) :: len_message
    integer(I4B) :: junit
    integer(I4B) :: ilevel
    integer(I4B) :: i
    integer(I4B) :: j

    ! -- abort if message is empty
    if (len_trim(text) < 1) return

    ! -- initialize local variables
    amessage = text
    counter = ''
    fmt_first = '(A)'
    fmt_cont = '(A)'
    len_str1 = 0
    len_str2 = len_line
    include_counter = .FALSE.
    junit = iunit
    j = 0

    ! -- process optional arguments
    if (present(skipbefore)) then
      isb = skipbefore
    else
      isb = 0
    end if
    if (present(skipafter)) then
      isa = skipafter
    else
      isa = 0
    end if

    ! -- create the counter to prepend to the start of the message,
    !    formats, and variables used to create strings
    if (present(iwidth) .and. present(icount)) then
      include_counter = .TRUE.

      ! -- write counter
      write (cfmt, '(A,I0,A)') '(1x,i', iwidth, ',".",1x)'
      write (counter, cfmt) icount

      ! -- calculate the length of the first and second string on a line
      len_str1 = len(trim(counter)) + 1
      len_str2 = len_line - len_str1

      ! -- write format for the continuation lines
      write (fmt_cont, '(a,i0,a)') &
        '(', len(trim(counter)) + 1, 'x,a)'
    end if

    ! -- calculate the length of the message
    len_message = len_trim(amessage)

    ! -- parse the message into multiple lines
5   continue
    jend = j + len_str2
    if (jend >= len_message) go to 100
    do i = jend, j + 1, -1
    if (amessage(i:i) .eq. ' ') then
    if (j == 0) then
    if (include_counter) then
      line = counter(1:len_str1)//amessage(j + 1:i)
    else
      line = amessage(j + 1:i)
    end if
    call write_message(iunit=junit, text=line, &
                     fmt=fmt_first, &
                     skipbefore=isb)
    else
    line = adjustl(amessage(j + 1:i))
    call write_message(iunit=junit, text=line, &
                     fmt=fmt_cont)
    end if
    j = i
    go to 5
    end if
    end do
    if (j == 0) then
    if (include_counter) then
      line = counter(1:len_str1)//amessage(j + 1:jend)
    else
      line = amessage(j + 1:jend)
    end if
    call write_message(iunit=junit, text=line, &
                     fmt=fmt_first, &
                     skipbefore=isb)
    else
    line = amessage(j + 1:jend)
    call write_message(iunit=junit, text=line, &
                     fmt=fmt_cont)
    end if
    j = jend
    go to 5

    ! -- last piece of amessage to write to a line
100 continue
    jend = len_message
    if (j == 0) then
    if (include_counter) then
      line = counter(1:len_str1)//amessage(j + 1:jend)
    else
      line = amessage(j + 1:jend)
    end if
    call write_message(iunit=junit, text=line, &
                     fmt=fmt_first, &
                     skipbefore=isb, skipafter=isa)
    else
    line = amessage(j + 1:jend)
    call write_message(iunit=junit, text=line, fmt=fmt_cont, &
                     skipafter=isa)
    end if
  end subroutine write_message_counter

  !> @brief Write horizontally centered text, left-padding as needed.
  subroutine write_message_centered(iunit, text, linelen)
    ! -- dummy
    integer(I4B), intent(in) :: iunit !< output unit to write text
    character(len=*), intent(in) :: text !< message to write to iunit
    integer(I4B), intent(in) :: linelen !< length of line to center text in
    ! -- local
    character(len=linelen) :: line
    character(len=linelen) :: blank
    integer(I4B) :: iu
    integer(I4B) :: len_message
    integer(I4B) :: jend
    integer(I4B) :: ipad
    integer(I4B) :: i
    integer(I4B) :: j

    ! -- initialize local variables
    blank = ''
    len_message = len_trim(adjustl(text))
    j = 0

    ! -- parse the amessage into multiple lines
    outer: do while (.true.)
      jend = j + linelen

      ! last line
      if (jend >= len_message) then
        jend = len_message
        line = text(j + 1:jend)
        ipad = ((linelen - len_trim(line)) / 2)
        call write_message(blank(1:ipad)//line, iunit=iunit)
        exit outer
      end if

      do i = jend, j + 1, -1
        if (text(i:i) .eq. ' ') then
          line = text(j + 1:i)
          ipad = ((linelen - len_trim(line)) / 2)
          call write_message(blank(1:ipad)//line, iunit=iunit)
          j = i
          cycle outer
        end if
      end do

      line = text(j + 1:jend)
      ipad = ((linelen - len_trim(line)) / 2)
      call write_message(blank(1:ipad)//line, iunit=iunit)
      j = jend
    end do outer

  end subroutine write_message_centered

end module MessageModule
