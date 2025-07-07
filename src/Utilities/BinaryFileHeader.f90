module BinaryFileHeaderModule
  use KindModule, only: I4B, DP, LGP
  implicit none

  type :: BinaryFileHeaderType
    integer(I4B) :: kstp
    integer(I4B) :: kper
    integer(I4B) :: kpstlast
    integer(I4B) :: kperlast
    integer(I4B) :: kstpnext
    integer(I4B) :: kpernext
    integer(I4B) :: pos
    real(DP) :: delt
    real(DP) :: pertim
    real(DP) :: totim
  end type BinaryFileHeaderType

end module BinaryFileHeaderModule