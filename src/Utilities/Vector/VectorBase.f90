module VectorBaseModule
  use KindModule, only: I4B, DP
  implicit none
  private

  type, public, abstract :: VectorBaseType
  contains
    procedure(create_if), deferred :: create
    procedure(destroy_if), deferred :: destroy
    procedure(get_array_if), deferred :: get_array

    procedure(zero_entries_if), deferred :: zero_entries
  end type VectorBaseType

  abstract interface
    subroutine create_if(this, n, name, mem_path)
      import VectorBaseType, I4B
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: n !< the nr. of elements in the (local) vector
      character(len=*) :: name !< the variable name (for access through memory manager)
      character(len=*) :: mem_path !< memory path for storing the underlying memory items
    end subroutine create_if
    subroutine destroy_if(this)
      import VectorBaseType
      class(VectorBaseType) :: this !< this vector
    end subroutine destroy_if
    function get_array_if(this) result(array)
      import VectorBaseType, DP
      class(VectorBaseType) :: this !< this vector
      real(DP), dimension(:), pointer, contiguous :: array
    end function get_array_if
    subroutine zero_entries_if(this)
      import VectorBaseType
      class(VectorBaseType) :: this
    end subroutine
  end interface

end module VectorBaseModule