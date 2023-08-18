module MethodCellPoolModule

  use MethodCellPollockModule
  use MethodCellPollockQuadModule
  use MethodCellTernaryModule
  use MethodCellPassToBotModule
  implicit none

  private
  public :: create_methodCellPool
  public :: destroy_methodCellPool

  type(MethodCellPollockType), pointer &
    , public :: methodCellPollock => null() ! method for the method pool
  type(MethodCellPollockQuadType), pointer &
    , public :: methodCellPollockQuad => null() ! method for the method pool
  type(MethodCellTernaryType), pointer &
    , public :: methodCellTernary => null() ! method for the method pool
  type(MethodCellPassToBotType), pointer &
    , public :: methodCellPassToBot => null() ! method for the method pool

contains

  !> @brief Create the cell method pool
  subroutine create_methodCellPool()
    !
    ! -- Create cell method pool
    call create_methodCellPollock(methodCellPollock)
    call create_methodCellPollockQuad(methodCellPollockQuad)
    call create_methodCellTernary(methodCellTernary)
    call create_methodCellPassToBot(methodCellPassToBot)
    !
    return
    !
  end subroutine create_methodCellPool

  !> @brief Destroy the cell method pool
  subroutine destroy_methodCellPool()
    !
    call methodCellPollock%destroy()
    deallocate (methodCellPollock)
    call methodCellPollockQuad%destroy()
    deallocate (methodCellPollockQuad)
    call methodCellTernary%destroy()
    deallocate (methodCellTernary)
    call methodCellPassToBot%destroy()
    deallocate (methodCellPassToBot)
    !
    return
    !
  end subroutine destroy_methodCellPool

end module MethodCellPoolModule
