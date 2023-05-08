module MethodSubcellPoolModule ! kluge

  use MethodSubcellPollockModule
  use MethodSubcellTernaryModule
  implicit none

  private
  public :: create_methodSubcellPool
  public :: destroy_methodSubcellPool

  type(MethodSubcellPollockType), pointer, public :: methodSubcellPollock &
                                                     => null()
  type(MethodSubcellTernaryType), pointer, public :: methodSubcellTernary &
                                                     => null()

contains

  !> @brief Create the subcell method pool
  subroutine create_methodSubcellPool()
    !
    call create_methodSubcellPollock(methodSubcellPollock)
    call create_methodSubcellTernary(methodSubcellTernary)
    !
    return
    !
  end subroutine create_methodSubcellPool

  !> @brief Destroy the subcell method pool
  subroutine destroy_methodSubcellPool()
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
