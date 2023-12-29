module ModelsCreateModule
  use KindModule, only: I4B
  use ConstantsModule, only: LENMEMPATH, LINELENGTH, LENMODELNAME
  use ListsModule, only: basemodellist
  use SimModule, only: store_error, count_errors, &
                       store_error_filename, MaxErrors
  use SimVariablesModule, only: iout, model_names, model_ranks, model_loc_idx, &
                                idm_context, proc_id, nr_procs, simulation_mode
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_setptr, mem_allocate
  use CharacterStringModule, only: CharacterStringType

  implicit none
  private
  public :: models_create, create_load_mask

contains

  subroutine add_gwf_model(n, fname)
    ! -- modules
    use NumericalModelModule, only: NumericalModelType, &
                                    get_numerical_model_from_list
    use VirtualGwfModelModule, only: add_virtual_gwf_model
    use GwfModule, only: gwf_cr
    ! -- dummy
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: fname
    ! -- local
    class(NumericalModelType), pointer :: model

    model => null() ! can be null for remote models
    if (model_ranks(n) == proc_id) then
      write (iout, '(4x,2a,i0,a)') 'GWF6', ' model ', &
        n, ' will be created'
      call gwf_cr(fname, n, model_names(n))
      model => get_numerical_model_from_list(basemodellist, n)
      model_loc_idx(n) = n
    end if
    call add_virtual_gwf_model(n, model_names(n), model)

  end subroutine add_gwf_model

  subroutine add_gwt_model(n, fname)
    ! -- modules
    use NumericalModelModule, only: NumericalModelType, &
                                    get_numerical_model_from_list
    use VirtualGwtModelModule, only: add_virtual_gwt_model
    use GwtModule, only: gwt_cr
    ! -- dummy
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: fname
    ! -- local
    class(NumericalModelType), pointer :: model

    model => null() ! can be null for remote models
    if (model_ranks(n) == proc_id) then
      write (iout, '(4x,2a,i0,a)') 'GWT6', ' model ', &
        n, ' will be created'
      call gwt_cr(fname, n, model_names(n))
      model => get_numerical_model_from_list(basemodellist, n)
      model_loc_idx(n) = n
    end if
    call add_virtual_gwt_model(n, model_names(n), model)

  end subroutine add_gwt_model

  !> @brief Check that the model name is valid
  !<
  subroutine check_model_name(mtype, mname)
    ! -- dummy
    character(len=*), intent(in) :: mtype
    character(len=*), intent(inout) :: mname
    ! -- local
    integer :: ilen
    integer :: i
    character(len=LINELENGTH) :: errmsg
    logical :: terminate = .true.
    !
    ilen = len_trim(mname)
    if (ilen > LENMODELNAME) then
      write (errmsg, '(a,a)') 'Invalid model name: ', trim(mname)
      call store_error(errmsg)
      write (errmsg, '(a,i0,a,i0)') &
        'Name length of ', ilen, ' exceeds maximum length of ', &
        LENMODELNAME
      call store_error(errmsg, terminate)
    end if
    do i = 1, ilen
      if (mname(i:i) == ' ') then
        write (errmsg, '(a,a)') 'Invalid model name: ', trim(mname)
        call store_error(errmsg)
        write (errmsg, '(a)') &
          'Model name cannot have spaces within it.'
        call store_error(errmsg, terminate)
      end if
    end do

  end subroutine check_model_name

  !> @brief Create a load mask to determine which models
  !! should be loaded by idm on this process. This is in
  !! sync with models create. The mask array should be
  !! pre-allocated with size equal to the global number
  !! of models. It is returned as (1, 1, 0, 0, ... 0)
  !! with each entry being a load mask for the model
  !! at the corresponding location in the 'MNAME' array
  !< of the IDM.
  subroutine create_load_mask(mask_array)
    integer(I4B), dimension(:) :: mask_array
    ! local
    integer(I4B) :: i

    call create_load_balance(mask_array)
    do i = 1, size(mask_array)
      if (mask_array(i) == proc_id) then
        mask_array(i) = 1
      else
        mask_array(i) = 0
      end if
    end do

  end subroutine create_load_mask

  !> @brief Distribute the models over the available
  !! processes in a parallel run. Expects an array sized
  !< to the number of models in the global simulation
  subroutine create_load_balance(mranks)
    integer(I4B), dimension(:) :: mranks
    ! local
    integer(I4B) :: im, imm, ie, ip, cnt
    integer(I4B) :: nr_models, nr_gwf_models, nr_gwt_models
    integer(I4B) :: nr_exchanges
    integer(I4B) :: min_per_proc, nr_left
    integer(I4B) :: rank
    integer(I4B), dimension(:), allocatable :: nr_models_proc
    character(len=:), allocatable :: model_type_str
    character(len=LINELENGTH) :: errmsg
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: etypes !< exg types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_a !< model a names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_b !< model b names

    mranks = 0
    if (simulation_mode /= 'PARALLEL') return

    ! load IDM data
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)
    call mem_setptr(etypes, 'EXGTYPE', input_mempath)
    call mem_setptr(emnames_a, 'EXGMNAMEA', input_mempath)
    call mem_setptr(emnames_b, 'EXGMNAMEB', input_mempath)

    ! count flow models
    nr_models = size(mnames)
    nr_gwf_models = 0
    nr_gwt_models = 0
    do im = 1, nr_models
      if (mtypes(im) == 'GWF6') then
        nr_gwf_models = nr_gwf_models + 1
      else if (mtypes(im) == 'GWT6') then
        nr_gwt_models = nr_gwt_models + 1
      else
        model_type_str = mtypes(im)
        write (errmsg, *) 'Model type ', model_type_str, &
          ' not supported in parallel mode.'
        call store_error(errmsg, terminate=.true.)
      end if
    end do

    ! calculate nr of flow models for each rank
    allocate (nr_models_proc(nr_procs))
    min_per_proc = nr_gwf_models / nr_procs
    nr_left = nr_gwf_models - nr_procs * min_per_proc
    cnt = 1
    do ip = 1, nr_procs
      rank = ip - 1
      nr_models_proc(ip) = min_per_proc
      if (rank < nr_left) then
        nr_models_proc(ip) = nr_models_proc(ip) + 1
      end if
    end do

    ! assign ranks for flow models
    rank = 0
    do im = 1, nr_models
      if (mtypes(im) == 'GWF6') then
        if (nr_models_proc(rank + 1) == 0) then
          rank = rank + 1
        end if
        mranks(im) = rank
        nr_models_proc(rank + 1) = nr_models_proc(rank + 1) - 1
      end if
    end do

    ! match transport to flow
    nr_exchanges = size(etypes)

    do im = 1, nr_models
      if (.not. mtypes(im) == 'GWT6') cycle

      ! find match
      do ie = 1, nr_exchanges
        if (etypes(ie) == 'GWF6-GWT6' .and. mnames(im) == emnames_b(ie)) then
          ! this is the exchange, now find the flow model's rank
          rank = 0
          do imm = 1, nr_models
            if (mnames(imm) == emnames_a(ie)) then
              rank = mranks(imm)
              exit
            end if
          end do

          ! we have our rank, assign and go to next transport model
          mranks(im) = rank
          exit
        end if
      end do
    end do

    ! cleanup
    deallocate (nr_models_proc)

  end subroutine create_load_balance

  !> @brief Set the models to be used for the simulation
  !<
  subroutine models_create()
    ! -- locals
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mfnames !< model file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    character(len=LINELENGTH) :: errmsg, model_name, model_type, fname
    integer(I4B) :: im, n, nr_models_glob
    !
    ! -- set input memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    !
    ! -- set pointers to input context model attribute arrays
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mfnames, 'MFNAME', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)
    !
    ! -- allocate global arrays
    nr_models_glob = size(mnames)
    call mem_allocate(model_ranks, nr_models_glob, 'MRANKS', input_mempath)
    allocate (model_names(nr_models_glob))
    allocate (model_loc_idx(nr_models_glob))
    !
    ! -- assign models to cpu cores (in serial all to rank 0)
    call create_load_balance(model_ranks)
    !
    ! -- open model logging block
    write (iout, '(/1x,a)') 'READING SIMULATION MODELS'
    !
    ! -- create models
    im = 0
    do n = 1, size(mtypes)
      !
      ! -- attributes for this model
      model_type = mtypes(n)
      fname = mfnames(n)
      model_name = mnames(n)
      !
      call check_model_name(model_type, model_name)
      !
      ! increment global model id
      model_names(n) = model_name(1:LENMODELNAME)
      model_loc_idx(n) = -1
      !
      ! -- add a new (local or global) model
      select case (model_type)
      case ('GWF6')
        call add_gwf_model(n, fname)
        im = im + 1
      case ('GWT6')
        call add_gwt_model(n, fname)
        im = im + 1
      case default
        write (errmsg, '(a,a)') &
          'Unknown simulation model type: ', trim(model_type)
        call store_error(errmsg, terminate=.true.)
      end select
    end do
    !
    ! -- close model logging block
    write (iout, '(1x,a)') 'END OF SIMULATION MODELS'
    !
    ! -- sanity check
    if (simulation_mode == 'PARALLEL' .and. im == 0) then
      write (errmsg, '(a, i0)') &
        'No MODELS assigned to process ', proc_id
      call store_error(errmsg, terminate=.true.)
    end if

  end subroutine models_create

end module ModelsCreateModule
