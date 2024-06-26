! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmPrtDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use PrtNamInputModule
  use PrtDisInputModule
  use PrtDisvInputModule
  use PrtMipInputModule

  implicit none
  private
  public :: prt_param_definitions
  public :: prt_aggregate_definitions
  public :: prt_block_definitions
  public :: prt_idm_multi_package
  public :: prt_idm_subpackages
  public :: prt_idm_integrated

contains

  subroutine set_param_pointer(input_dfn, input_dfn_target)
    type(InputParamDefinitionType), dimension(:), pointer :: input_dfn
    type(InputParamDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_param_pointer

  subroutine set_block_pointer(input_dfn, input_dfn_target)
    type(InputBlockDefinitionType), dimension(:), pointer :: input_dfn
    type(InputBlockDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_block_pointer

  subroutine set_subpkg_pointer(subpkg_list, subpkg_list_target)
    character(len=16), dimension(:), pointer :: subpkg_list
    character(len=16), dimension(:), target :: subpkg_list_target
    subpkg_list => subpkg_list_target
  end subroutine set_subpkg_pointer

  function prt_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, prt_nam_param_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, prt_dis_param_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, prt_disv_param_definitions)
    case ('MIP')
      call set_param_pointer(input_definition, prt_mip_param_definitions)
    case default
    end select
    return
  end function prt_param_definitions

  function prt_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, prt_nam_aggregate_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, prt_dis_aggregate_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, prt_disv_aggregate_definitions)
    case ('MIP')
      call set_param_pointer(input_definition, prt_mip_aggregate_definitions)
    case default
    end select
    return
  end function prt_aggregate_definitions

  function prt_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, prt_nam_block_definitions)
    case ('DIS')
      call set_block_pointer(input_definition, prt_dis_block_definitions)
    case ('DISV')
      call set_block_pointer(input_definition, prt_disv_block_definitions)
    case ('MIP')
      call set_block_pointer(input_definition, prt_mip_block_definitions)
    case default
    end select
    return
  end function prt_block_definitions

  function prt_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = prt_nam_multi_package
    case ('DIS')
      multi_package = prt_dis_multi_package
    case ('DISV')
      multi_package = prt_disv_multi_package
    case ('MIP')
      multi_package = prt_mip_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="PRT"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function prt_idm_multi_package

  function prt_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, prt_nam_subpackages)
    case ('DIS')
      call set_subpkg_pointer(subpackages, prt_dis_subpackages)
    case ('DISV')
      call set_subpkg_pointer(subpackages, prt_disv_subpackages)
    case ('MIP')
      call set_subpkg_pointer(subpackages, prt_mip_subpackages)
    case default
    end select
    return
  end function prt_idm_subpackages

  function prt_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('DIS')
      integrated = .true.
    case ('DISV')
      integrated = .true.
    case ('MIP')
      integrated = .true.
    case default
    end select
    return
  end function prt_idm_integrated

end module IdmPrtDfnSelectorModule
