module TrackDataModule

  use KindModule, only: DP, I4B

  implicit none

  private
  public :: TrackDataType

  type :: TrackDataType
    integer(I4B), pointer :: ntrack => null() ! track data counter
    integer(I4B), dimension(:), pointer, contiguous :: iptrack ! kluge
    integer(I4B), dimension(:), pointer, contiguous :: ictrack ! kluge
    real(DP), dimension(:), pointer, contiguous :: xtrack ! kluge
    real(DP), dimension(:), pointer, contiguous :: ytrack ! kluge
    real(DP), dimension(:), pointer, contiguous :: ztrack ! kluge
    real(DP), dimension(:), pointer, contiguous :: ttrack ! kluge
  end type TrackDataType

end module TrackDataModule
