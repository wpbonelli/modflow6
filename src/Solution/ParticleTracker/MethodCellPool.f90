module MethodCellPoolModule ! kluge

  use MethodCellPollockModule
  use MethodCellPollockQuadModule
  use MethodCellTernaryModule
  implicit none

  private
  public :: create_methodCellPool ! create the cell method pool object
  public :: destroy_methodCellPool ! destroy the cell method pool object

  type(MethodCellPollockType), pointer, public :: methodCellPollock => null() ! method for the method pool
  type(MethodCellPollockQuadType), pointer, public :: methodCellPollockQuad   => null()   ! method for the method pool
  type(MethodCellTernaryType), pointer, public :: methodCellTernary => null() ! method for the method pool

contains

  subroutine create_methodCellPool()
! ******************************************************************************
! create_methodCellPool -- Create the cell method pool
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
! ------------------------------------------------------------------------------
    !
    ! -- Create cell method pool
    call create_methodCellPollock(methodCellPollock)
    call create_methodCellPollockQuad(methodCellPollockQuad)
    call create_methodCellTernary(methodCellTernary)
    !
    return
    !
  end subroutine create_methodCellPool

  subroutine destroy_methodCellPool()
! ******************************************************************************
! destroy_methodCellPool -- Destroy the cell method pool
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
! ------------------------------------------------------------------------------
    !
    call methodCellPollock%destroy()
    deallocate (methodCellPollock)
    call methodCellPollockQuad%destroy()
    deallocate (methodCellPollockQuad)
    call methodCellTernary%destroy()
    deallocate (methodCellTernary)
    !
    return
    !
  end subroutine destroy_methodCellPool

end module MethodCellPoolModule
