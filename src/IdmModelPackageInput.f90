module ModelPackageInputModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LENFTYPE, LENPACKAGETYPE
  use SimModule, only: store_error, store_error_filename
  use gwfModule, only: GWF_NBASEPKG, GWF_NMULTIPKG, &
                         GWF_BASEPKG, GWF_MULTIPKG
  use gwtModule, only: GWT_NBASEPKG, GWT_NMULTIPKG, &
                         GWT_BASEPKG, GWT_MULTIPKG
  use prtModule, only: PRT_NBASEPKG, PRT_NMULTIPKG, &
                         PRT_BASEPKG, PRT_MULTIPKG

  implicit none
  private
  public :: supported_model_packages
  public :: multi_package_type

contains

  !> @brief set supported package types for model
  !!
  !! Allocate a list of package types supported
  !! by the model.  Base packages should be listed
  !! first as list determines load order.
  !!
  !<
  subroutine supported_model_packages(mtype, pkgtypes, numpkgs)
    ! -- modules
    ! -- dummy
    character(len=LENFTYPE), intent(in) :: mtype
    character(len=LENPACKAGETYPE), dimension(:), allocatable, &
      intent(inout) :: pkgtypes
    integer(I4B), intent(inout) :: numpkgs
    ! -- local
    !
    select case (mtype)
    case ('GWF6')
      numpkgs = GWF_NBASEPKG + GWF_NMULTIPKG
      allocate (pkgtypes(numpkgs))
      pkgtypes = [GWF_BASEPKG, GWF_MULTIPKG]
    case ('GWT6')
      numpkgs = GWT_NBASEPKG + GWT_NMULTIPKG
      allocate (pkgtypes(numpkgs))
      pkgtypes = [GWT_BASEPKG, GWT_MULTIPKG]
    case ('PRT6')
      numpkgs = PRT_NBASEPKG + PRT_NMULTIPKG
      allocate (pkgtypes(numpkgs))
      pkgtypes = [PRT_BASEPKG, PRT_MULTIPKG]
    case default
    end select
    !
    ! -- return
    return
  end subroutine supported_model_packages

  !> @brief Is the package multi-instance
  !<
  function multi_package_type(mtype_component, ptype_component, pkgtype) &
    result(multi_package)
    ! -- modules
    ! -- dummy
    character(len=LENFTYPE), intent(in) :: mtype_component
    character(len=LENFTYPE), intent(in) :: ptype_component
    character(len=LENFTYPE), intent(in) :: pkgtype
    ! -- return
    logical(LGP) :: multi_package
    ! -- local
    integer(I4B) :: n
    !
    multi_package = .false.
    !
    select case (mtype_component)
    case ('GWF')
      do n = 1, GWF_NMULTIPKG
        if (GWF_MULTIPKG(n) == pkgtype) then
          multi_package = .true.
          exit
        end if
      end do
    case ('GWT')
      do n = 1, GWT_NMULTIPKG
        if (GWT_MULTIPKG(n) == pkgtype) then
          multi_package = .true.
          exit
        end if
      end do
    case ('PRT')
      do n = 1, PRT_NMULTIPKG
        if (PRT_MULTIPKG(n) == pkgtype) then
          multi_package = .true.
          exit
        end if
      end do
    case default
    end select
    !
    ! -- return
    return
  end function multi_package_type

end module ModelPackageInputModule