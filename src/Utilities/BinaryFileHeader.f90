module BinaryFileHeaderModule

  use KindModule, only: I4B, LGP, DP
  use ConstantsModule, only: LINELENGTH
  implicit none

  public :: BinaryFileHeaderType, BinaryFileHeaderWrapperType

  type :: BinaryFileHeaderType
    integer(I4B) :: pos
    integer(I4B) :: kper, kstp
    real(DP) :: delt, pertim, totim
    character(len=16) :: text
  contains
    procedure :: get_str
  end type BinaryFileHeaderType

  type :: BinaryFileHeaderWrapperType
    class(BinaryFileHeaderType), allocatable :: header
  end type BinaryFileHeaderWrapperType

contains

  !> @brief Get a string representation of the header.
  function get_str(this) result(str)
    class(BinaryFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str

    write (str, '(*(G0))') &
      'Binary file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ', text: ', trim(this%text), &
      ')'
    str = trim(str)
  end function get_str

end module BinaryFileHeaderModule
