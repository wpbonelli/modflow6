module FileUtilModule

  use KindModule, only: DP, I4B, LGP
  use DefinedMacros, only: get_os
  use ConstantsModule, only: OSWIN, OSUNDEF, IUSTART, IULAST

  implicit none
  private

  ! members
  integer(I4B), public :: iunext = IUSTART !< next open file unit number to assign

  ! procedures
  public :: get_filename, get_fileunit, close_all_files

contains

  !> @brief Get the file name corresponding to an open file unit number.
  subroutine get_filename(iunit, fname)
    ! -- dummy variables
    integer(I4B), intent(in) :: iunit !< open file unit number
    character(len=*), intent(inout) :: fname !< file name
    ! -- local variables
    integer(I4B) :: ipos
    integer(I4B) :: ios
    integer(I4B) :: ilen

    ! -- get file name from unit number
    inquire (unit=iunit, name=fname)

    ! -- determine the operating system
    ios = get_os()

    ! -- extract filename from full path if needed,
    !    first checking forward slash (non-windows)
    if (ios /= OSWIN) &
      ipos = index(fname, '/', back=.TRUE.)
    ! -- check backslash on windows or undefined os
    !    if forward slash not found
    if ((ios == OSWIN .or. ios == OSUNDEF) .and. (ipos < 1)) &
      ipos = index(fname, '\', back=.TRUE.)

    ! -- strip path from file name
    if (ipos > 0) then
      ilen = len_trim(fname)
      write (fname, '(a)') fname(ipos + 1:ilen)//' '
    end if
  end subroutine get_filename

  !> @brief Retrieve the next unassigned unit number.
  function get_fileunit() result(iu)
    integer(I4B) :: iu
    integer(I4B) :: i
    logical :: opened

    do i = iunext, IULAST
      inquire (unit=i, opened=opened)
      if (.not. opened) exit
    end do
    iu = i
    iunext = iu + 1
  end function get_fileunit

  !> @brief Close all open files.
  subroutine close_all_files()
    ! -- local
    integer(I4B) :: i
    logical :: opened
    character(len=7) :: output_file

    ! -- close all open file units
    do i = iustart, iunext - 1
      ! -- determine if file is open, if not, skip it
      inquire (unit=i, opened=opened)
      if (.not. opened) cycle

      ! -- flush the file if it can be written to
      inquire (unit=i, write=output_file)
      if (trim(adjustl(output_file)) == 'YES') flush (i)

      ! -- close file
      close (i)
    end do
  end subroutine close_all_files

end module FileUtilModule
