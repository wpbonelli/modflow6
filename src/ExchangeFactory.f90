module ExchangeFactoryModule
  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LINELENGTH
  use SimModule, only: store_error
  use SimVariablesModule, only: iout, model_names, model_loc_idx
  use CharacterStringModule, only: CharacterStringType
  use ArrayHandlersModule, only: ifind
  use GwfGwfExchangeModule, only: register_gwfgwf
  use GwfGwtExchangeModule, only: register_gwfgwt
  use GwfGweExchangeModule, only: register_gwfgwe
  use GwfPrtExchangeModule, only: register_gwfprt
  use GwtGwtExchangeModule, only: register_gwtgwt
  use GweGweExchangeModule, only: register_gwegwe
  use SwfGwfExchangeModule, only: register_swfgwf
  use VirtualGwfExchangeModule, only: register_virtual_gwfgwf
  use VirtualGwtExchangeModule, only: register_virtual_gwtgwt
  use VirtualGweExchangeModule, only: register_virtual_gwegwe
  use VirtualPrtExchangeModule, only: register_virtual_prtprt

  implicit none
  private
  public :: create_exchanges

contains

  subroutine create_exchanges(etypes, efiles, emnames_a, emnames_b, emempaths)
    ! -- dummy
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(in) :: etypes !< exg types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(in) :: efiles !< exg file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(in) :: emnames_a !< model a names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(in) :: emnames_b !< model b names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(in) :: emempaths
    ! -- local
    integer(I4B) :: exg_id, n
    integer(I4B) :: m1_id, m2_id
    logical(LGP) :: both_remote, both_local
    character(len=LINELENGTH) :: fname, name1, name2, exg_name
    character(len=LENMEMPATH) :: exg_mempath
    character(len=LINELENGTH) :: errmsg, exgtype
    ! -- formats
    character(len=*), parameter :: fmtmerr = "('Error in simulation control ', &
      &'file.  Could not find model: ', a)"

    exg_id = 0
    do n = 1, size(etypes)
      exgtype = etypes(n)
      fname = efiles(n)
      name1 = emnames_a(n)
      name2 = emnames_b(n)
      exg_mempath = emempaths(n)
      exg_id = exg_id + 1

      ! find model index in list
      m1_id = ifind(model_names, name1)
      if (m1_id < 0) then
        write (errmsg, fmtmerr) trim(name1)
        call store_error(errmsg, terminate=.true.)
      end if
      m2_id = ifind(model_names, name2)
      if (m2_id < 0) then
        write (errmsg, fmtmerr) trim(name2)
        call store_error(errmsg, terminate=.true.)
      end if

      ! both models on other process? then don't create it here...
      both_remote = (model_loc_idx(m1_id) == -1 .and. &
                     model_loc_idx(m2_id) == -1)
      both_local = (model_loc_idx(m1_id) > 0 .and. &
                    model_loc_idx(m2_id) > 0)
      if (.not. both_remote) write (iout, '(4x,a,a,i0,a,i0,a,i0)') &
        trim(exgtype), ' exchange ', exg_id, &
        ' will be created to connect model ', m1_id, &
        ' with model ', m2_id

      select case (exgtype)
      case ('GWF6-GWF6')
        write (exg_name, '(a,i0)') 'GWF-GWF_', exg_id
        if (.not. both_remote) &
          call register_gwfgwf( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
        call register_virtual_gwfgwf( &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id)
      case ('GWT6-GWT6')
        write (exg_name, '(a,i0)') 'GWT-GWT_', exg_id
        if (.not. both_remote) &
          call register_gwtgwt( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
        call register_virtual_gwtgwt( &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id)
      case ('GWE6-GWE6')
        write (exg_name, '(a,i0)') 'GWE-GWE_', exg_id
        if (.not. both_remote) &
          call register_gwegwe( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
        call register_virtual_gwegwe( &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id)
      case ('GWF6-GWT6')
        write (exg_name, '(a,i0)') 'GWF-GWT_', exg_id
        if (both_local) &
          call register_gwfgwt( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
      case ('GWF6-GWE6')
        write (exg_name, '(a,i0)') 'GWF-GWE_', exg_id
        if (both_local) &
          call register_gwfgwe( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
      case ('GWF6-PRT6')
        write (exg_name, '(a,i0)') 'GWF-PRT_', exg_id
        if (both_local) &
          call register_gwfprt( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
      case ('SWF6-GWF6')
        write (exg_name, '(a,i0)') 'SWF-GWF_', exg_id
        if (both_local) &
          call register_swfgwf( &
          fname, &
          exg_name, &
          exg_id, &
          m1_id, &
          m2_id, &
          exg_mempath)
      case default
        write (errmsg, '(a,a)') &
          'Unknown simulation exchange type: ', trim(exgtype)
        call store_error(errmsg, terminate=.true.)
      end select
    end do
  end subroutine create_exchanges

end module ExchangeFactoryModule
