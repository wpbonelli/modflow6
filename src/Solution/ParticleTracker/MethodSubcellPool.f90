module MethodSubcellPoolModule ! kluge

  use MethodSubcellPollockModule
  use MethodSubcellTernaryModule
  implicit none

  private
  public :: create_methodSubcellPool ! create the subcell method pool object
  public :: destroy_methodSubcellPool ! destroy the subcell method pool object

type(MethodSubcellPollockType), pointer, public :: methodSubcellPollock => null() ! method for the method pool
type(MethodSubcellTernaryType), pointer, public :: methodSubcellTernary => null() ! method for the method pool

contains

  subroutine create_methodSubcellPool()
! ******************************************************************************
! create_methodSubcellPool -- Create the subcell method pool
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
! ------------------------------------------------------------------------------
    !
    ! -- Create subcell method pool
    call create_methodSubcellPollock(methodSubcellPollock)
    call create_methodSubcellTernary(methodSubcellTernary)
    !
    return
    !
  end subroutine create_methodSubcellPool

  subroutine destroy_methodSubcellPool()
! ******************************************************************************
! destroy_methodSubcellPool -- Destroy the subcell method pool
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
! ------------------------------------------------------------------------------
    !
    call methodSubcellPollock%destroy()
    deallocate (methodSubcellPollock)
    call methodSubcellTernary%destroy()
    deallocate (methodSubcellTernary)
    !
    return
    !
  end subroutine destroy_methodSubcellPool

end module MethodSubcellPoolModule
