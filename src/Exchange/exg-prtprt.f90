!> @brief PRT-PRT exchange module
module PrtPrtExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME, LINELENGTH, LENMEMPATH, &
                             LENMODELNAME, LENBOUNDNAME, LENAUXNAME
  use ListsModule, only: basemodellist, baseexchangelist
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use PrtModule, only: PrtModelType
  use ParticleModule, only: ParticleType, create_particle, &
                            create_particle_store, TERM_BOUNDARY
  use ParticleEventModule, only: ParticleEventType
  use MemoryHelperModule, only: create_mem_path
  use ExgPrtprtInputModule, only: ExgPrtprtParamFoundType
  use BndModule, only: BndType, GetBndFromList
  use CharacterStringModule
  use ExplicitModelModule, only: ExplicitModelType
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_deallocate
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use GeomUtilModule, only: get_ijk

  implicit none
  private

  public :: PrtPrtExchangeType
  public :: prtprt_cr, try_particle_transfer
  public :: GetPrtPrtExchangeFromList

  !> @brief Exchange that transfers particles between two PRT models.
  type, extends(BaseExchangeType) :: PrtPrtExchangeType
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file
    integer(I4B), pointer :: m1id => null() !< local index of model 1 in basemodellist
    integer(I4B), pointer :: m2id => null() !< local index of model 2 in basemodellist
    type(PrtModelType), pointer :: prtmodel1 => null() !< pointer to model 1
    type(PrtModelType), pointer :: prtmodel2 => null() !< pointer to model 2
    integer(I4B), pointer :: nexg => null() !< number of connected cell pairs
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    character(len=LENBOUNDNAME), dimension(:), &
      pointer, contiguous :: boundname => null() !< boundnames
    character(len=LENAUXNAME), dimension(:), &
      pointer, contiguous :: auxname => null() !< vector of auxname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< copy of vector auxname that can be stored in memory manager
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< array of auxiliary variable values
    integer(I4B), pointer :: iprpak => null() !< print input flag
    integer(I4B), pointer :: iprflow => null() !< print flag for cell by cell flows
    integer(I4B), pointer :: ipakcb => null() !< save flag for cell by cell flows
    integer(I4B), pointer :: inamedbound => null() !< flag to read boundnames
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() !< node numbers in model 1
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() !< node numbers in model 2
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< horizontal connection indicator array, size: nexg
    character(len=LENMODELNAME) :: gwfmodelname1 = '' !< name of gwf model for prt model 1
    character(len=LENMODELNAME) :: gwfmodelname2 = '' !< name of gwf model for prt model 2
    real(DP), dimension(:), pointer, contiguous :: gwfsimvals => null() !< gwf flows for exchange
    integer(I4B), pointer :: iflowface1 => null() !< index in auxvar for cell face number in model 1
    integer(I4B), pointer :: iflowface2 => null() !< index in auxvar for cell face number in model 2
    class(PrtPrtExchangeType), pointer :: self => null() !< pointer to self (for event subscriptions)
  contains
    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure :: exg_ad
    procedure :: connects_model
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_dimensions
    procedure, private :: source_data
    procedure, private :: noder
    procedure, private :: get_exchange_faces
  end type PrtPrtExchangeType

contains

  !> @brief Create a new PRT-PRT exchange and register it.
  subroutine prtprt_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! modules
    use SimVariablesModule, only: model_loc_idx
    ! dummy
    character(len=*), intent(in) :: filename !< exchange input file
    character(len=*), intent(in) :: name !< exchange name
    integer(I4B), intent(in) :: id !< exchange id
    integer(I4B), intent(in) :: m1_id !< global model 1 id
    integer(I4B), intent(in) :: m2_id !< global model 2 id
    character(len=*), intent(in) :: input_mempath !< IDM memory path for this exchange
    ! local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(PrtPrtExchangeType), pointer :: exchange => null()

    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    exchange%id = id
    exchange%name = name
    exchange%memoryPath = create_mem_path(exchange%name)
    exchange%input_mempath = input_mempath
    ! store a self-pointer for use in event subscriptions
    exchange%self => exchange

    ! allocate scalars and set defaults
    call exchange%allocate_scalars()
    exchange%filename = filename
    exchange%typename = 'PRT-PRT'

    ! NB: convert global model id to local index in basemodellist
    exchange%m1id = model_loc_idx(m1_id)
    exchange%m2id = model_loc_idx(m2_id)

    ! set model pointers
    call exchange%set_model_pointers()
  end subroutine prtprt_cr

  !> @brief Whether either of the models is connected by this exchange.
  logical function connects_model(this, model)
    class(PrtPrtExchangeType) :: this
    class(BaseModelType), pointer, intent(in) :: model

    select type (model)
    class is (PrtModelType)
      connects_model = (associated(this%prtmodel1, model) .or. &
                        associated(this%prtmodel2, model))
    end select
  end function connects_model

  !> @brief Define phase: store model pointers.
  subroutine exg_df(this)
    class(PrtPrtExchangeType) :: this
    ! local
    integer(I4B) :: iout
    class(*), pointer :: p

    iout = this%prtmodel1%iout
    call this%source_options(iout)
    call this%source_dimensions(iout)
    call this%allocate_arrays()
    call this%source_data(iout)
    !
    ! -- Validate that IFLOWFACE1 and IFLOWFACE2 are present
    if (this%iflowface1 == 0) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires IFLOWFACE1 as an auxiliary variable.'
      call store_error(errmsg)
    end if
    if (this%iflowface2 == 0) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires IFLOWFACE2 as an auxiliary variable.'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call store_error_filename(this%filename)
    end if
    !
    p => this%self ! self-pointer. can't point directly to `this`
    call this%prtmodel1%handlers%subscribe(try_particle_transfer, p)
    call this%prtmodel2%handlers%subscribe(try_particle_transfer, p)
  end subroutine exg_df

  !> @brief Allocate and read: load the connected-cell list from IDM.
  subroutine exg_ar(this)
    class(PrtPrtExchangeType) :: this
  end subroutine exg_ar

  !> @brief Deallocate memory.
  subroutine exg_da(this)
    class(PrtPrtExchangeType) :: this

    call mem_deallocate(this%m1id)
    call mem_deallocate(this%m2id)
    call mem_deallocate(this%nexg)
    call mem_deallocate(this%naux)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%iprflow)
    call mem_deallocate(this%ipakcb)
    call mem_deallocate(this%inamedbound)
    call mem_deallocate(this%iflowface1)
    call mem_deallocate(this%iflowface2)
    if (associated(this%nodem1)) then
      call mem_deallocate(this%nodem1, 'NODEM1', this%memoryPath)
    end if
    if (associated(this%nodem2)) then
      call mem_deallocate(this%nodem2, 'NODEM2', this%memoryPath)
    end if
    if (associated(this%ihc)) then
      call mem_deallocate(this%ihc, 'IHC', this%memoryPath)
    end if
    if (associated(this%auxvar)) then
      call mem_deallocate(this%auxvar, 'AUXVAR', this%memoryPath)
    end if
    if (associated(this%boundname)) then
      deallocate (this%boundname)
    end if
    if (associated(this%gwfsimvals)) then
      call mem_deallocate(this%gwfsimvals, 'GWFSIMVALS', this%memoryPath)
    end if
  end subroutine exg_da

  !> @brief Inject exchange boundary flows into FMI. Flows are
  !! consolidated, so that a coarse cell connected to multiple
  !! refined cells in another model will have as boundary flow
  !! the sum of flows associated with the cells connected to it.
  subroutine exg_ad(this)
    class(PrtPrtExchangeType) :: this
    integer(I4B) :: iexg, ic1, ic2, if1, if2
    real(DP) :: q

    do iexg = 1, this%nexg
      ic1 = this%nodem1(iexg)
      ic2 = this%nodem2(iexg)
      q = abs(this%gwfsimvals(iexg))
      call this%get_exchange_faces(iexg, if1, if2)
      call this%prtmodel1%fmi%add_boundary_flow(ic1, if1, -q)
      call this%prtmodel2%fmi%add_boundary_flow(ic2, if2, q)
    end do
  end subroutine

  subroutine get_exchange_faces(this, iexg, if1, if2)
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iexg !< exchange connection index
    integer(I4B), intent(out) :: if1, if2 !< face numbers
    !
    ! -- Get face numbers from auxiliary variables
    if1 = nint(this%auxvar(this%iflowface1, iexg))
    if2 = nint(this%auxvar(this%iflowface2, iexg))
  end subroutine

  subroutine set_model_pointers(this)
    class(PrtPrtExchangeType) :: this
    class(BaseModelType), pointer :: mb => null()

    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (PrtModelType)
      this%prtmodel1 => mb
    end select

    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      this%prtmodel2 => mb
    end select

    if (.not. associated(this%prtmodel1)) then
      write (errmsg, '(3a)') 'Problem with PRT-PRT exchange ', &
        trim(this%name), ': model 1 is not a PRT model.'
      call store_error(errmsg, terminate=.true.)
    end if
    if (.not. associated(this%prtmodel2)) then
      write (errmsg, '(3a)') 'Problem with PRT-PRT exchange ', &
        trim(this%name), ': model 2 is not a PRT model.'
      call store_error(errmsg, terminate=.true.)
    end if
  end subroutine set_model_pointers

  subroutine allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtPrtExchangeType) :: this

    allocate (this%filename)
    this%filename = ''

    call mem_allocate(this%m1id, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2id, 'M2ID', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    call mem_allocate(this%naux, 'NAUX', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%iprflow, 'IPRFLOW', this%memoryPath)
    call mem_allocate(this%ipakcb, 'IPAKCB', this%memoryPath)
    call mem_allocate(this%inamedbound, 'INAMEDBOUND', this%memoryPath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, &
                      'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, &
                      'AUXNAME_CST', this%memoryPath)
    call mem_allocate(this%iflowface1, 'IFLOWFACE1', this%memoryPath)
    call mem_allocate(this%iflowface2, 'IFLOWFACE2', this%memoryPath)

    this%m1id = 0
    this%m2id = 0
    this%nexg = 0
    this%naux = 0
    this%iprpak = 0
    this%iprflow = 0
    this%ipakcb = 0
    this%inamedbound = 0
    this%iflowface1 = 0
    this%iflowface2 = 0
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtPrtExchangeType) :: this

    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%ihc, this%nexg, 'IHC', this%memoryPath)
    ! auxname array is allocated while parsing
    call mem_allocate(this%auxvar, this%naux, this%nexg, &
                      'AUXVAR', this%memoryPath)
    ! allocate boundname
    if (this%inamedbound == 1) then
      allocate (this%boundname(this%nexg))
    else
      allocate (this%boundname(1))
    end if
    this%boundname(:) = ''
  end subroutine allocate_arrays

  !> @brief Read options from IDM.
  subroutine source_options(this, iout)
    use MemoryManagerExtModule, only: mem_set_value
    use ArrayHandlersModule, only: ifind
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! local
    type(ExgPrtprtParamFoundType) :: found
    logical(LGP) :: found_naux
    integer(I4B) :: n, ival

    call mem_set_value(this%gwfmodelname1, 'GWFMODELNAME1', this%input_mempath, &
                       found%gwfmodelname1)
    call mem_set_value(this%gwfmodelname2, 'GWFMODELNAME2', this%input_mempath, &
                       found%gwfmodelname2)
    call mem_set_value(this%naux, 'NAUX', this%input_mempath, found_naux)
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    call mem_set_value(this%iprpak, 'IPRPAK', this%input_mempath, found%iprpak)
    call mem_set_value(this%iprflow, 'IPRFLOW', this%input_mempath, found%iprflow)
    call mem_set_value(this%inamedbound, 'BOUNDNAMES', this%input_mempath, &
                       found%boundnames)

    ! reallocate aux arrays if aux variables provided
    if (found_naux .and. this%naux > 0) then
      call mem_reallocate(this%auxname, LENAUXNAME, this%naux, &
                          'AUXNAME', this%memoryPath)
      call mem_reallocate(this%auxname_cst, LENAUXNAME, this%naux, &
                          'AUXNAME_CST', this%memoryPath)
      call mem_set_value(this%auxname_cst, 'AUXILIARY', this%input_mempath, &
                         found%auxiliary)
      !
      do n = 1, this%naux
        this%auxname(n) = this%auxname_cst(n)
      end do
      !
      ! -- Look for IFLOWFACE1 and IFLOWFACE2 auxiliary variables
      ival = ifind(this%auxname, 'IFLOWFACE1')
      if (ival > 0) then
        this%iflowface1 = ival
      end if
      ival = ifind(this%auxname, 'IFLOWFACE2')
      if (ival > 0) then
        this%iflowface2 = ival
      end if
    end if

    if (.not. found%gwfmodelname1) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires that GWFMODELNAME1 be entered in the OPTIONS block.'
      call store_error(errmsg)
    end if
    if (.not. found%gwfmodelname2) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires that GWFMODELNAME2 be entered in the OPTIONS block.'
      call store_error(errmsg)
    end if

    if (found%ipakcb) then
      this%ipakcb = -1
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE SAVED TO BINARY BUDGET FILES.'
    end if

    if (found%iprpak) then
      write (iout, '(4x,a)') &
        'THE LIST OF EXCHANGES WILL BE PRINTED.'
    end if

    if (found%iprflow) then
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE PRINTED TO LIST FILES.'
    end if

    if (found%boundnames) then
      write (iout, '(4x,a)') 'EXCHANGE BOUNDARIES HAVE NAMES IN LAST COLUMN'
    end if
  end subroutine source_options

  !> @brief Read nexg from IDM.
  subroutine source_dimensions(this, iout)
    use MemoryManagerExtModule, only: mem_set_value
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    type(ExgPrtprtParamFoundType) :: found

    call mem_set_value(this%nexg, 'NEXG', this%input_mempath, found%nexg)

    if (.not. found%nexg) then
      write (errmsg, '(3a)') 'Required NEXG not found in exchange ', &
        trim(this%name), '.'
      call store_error(errmsg, terminate=.true.)
    end if
  end subroutine source_dimensions

  !> @brief Read cellidm1/cellidm2 from IDM and map to node numbers.
  subroutine source_data(this, iout)
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use GeomUtilModule, only: get_node
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm1 => null()
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm2 => null()
    integer(I4B) :: iexg, node1, node2, iaux
    integer(I4B), dimension(:), contiguous, pointer :: ihc
    real(DP), dimension(:, :), contiguous, pointer :: auxvar
    type(CharacterStringType), dimension(:), contiguous, pointer :: boundname

    call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
    call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)
    call mem_setptr(ihc, 'IHC', this%input_mempath)
    call mem_setptr(auxvar, 'AUXVAR', this%input_mempath)
    call mem_setptr(boundname, 'BOUNDNAME', this%input_mempath)

    do iexg = 1, this%nexg
      this%nodem1(iexg) = this%noder(this%prtmodel1, cellidm1(:, iexg), iout)
      this%nodem2(iexg) = this%noder(this%prtmodel2, cellidm2(:, iexg), iout)
      this%ihc(iexg) = ihc(iexg)
      do iaux = 1, this%naux
        this%auxvar(iaux, iexg) = auxvar(iaux, iexg)
      end do
      if (this%inamedbound == 1) then
        this%boundname(iexg) = boundname(iexg)
      end if
    end do
  end subroutine source_data

  !> @brief Returns reduced node number from user specified cell id.
  function noder(this, model, cellid, iout)
    ! modules
    use GeomUtilModule, only: get_node
    ! dummy
    class(PrtPrtExchangeType) :: this !< instance of exchange object
    type(PrtModelType), pointer, intent(in) :: model
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    integer(I4B) :: noder, node
    !
    if (model%dis%ndim == 1) then
      node = cellid(1)
    elseif (model%dis%ndim == 2) then
      node = get_node(cellid(1), 1, cellid(2), &
                      model%dis%mshape(1), 1, &
                      model%dis%mshape(2))
    else
      node = get_node(cellid(1), cellid(2), cellid(3), &
                      model%dis%mshape(1), &
                      model%dis%mshape(2), &
                      model%dis%mshape(3))
    end if
    noder = model%dis%get_nodenumber(node, 0)
  end function noder

  function try_particle_transfer(context, particle, event) result(transferred)
    use ModelExitEventModule, only: ModelExitEventType
    use PrtPrpModule, only: PrtPrpType, ExgPrtPrpType
    ! dummy
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: transferred
    ! local
    type(PrtPrtExchangeType), pointer :: exg
    class(BndType), pointer :: exgprp_obj
    type(ExgPrtPrpType), pointer :: exgprp
    integer(I4B) :: i, dest_cell
    logical(LGP) :: found

    transferred = .false.

    select type (context)
    type is (PrtPrtExchangeType)
      exg => context
    class default
      return
    end select

    select type (event)
    type is (ModelExitEventType)
      ! model1 -> model2
      found = .false.
      do i = 1, exg%nexg
        if (particle%icu == exg%nodem1(i)) then
          dest_cell = exg%nodem2(i)
          found = .true.
          exit
        end if
      end do
      if (found) then
        exgprp_obj => GetBndFromList( &
                      exg%prtmodel2%bndlist, &
                      exg%prtmodel2%bndlist%Count())
        select type (exgprp_obj)
        type is (ExgPrtPrpType)
          exgprp => exgprp_obj
          particle%itrdomain(1) = exg%prtmodel2%id
          particle%icu = dest_cell
          ! TODO: update particle coordinates
          exgprp%has_pending = .true.
          exgprp%nparticles = exgprp%nparticles + 1
          call exgprp%particles%resize(exgprp%nparticles, exgprp%memoryPath)
          call exgprp%particles%put(particle, exgprp%nparticles)
          transferred = .true.
          return
        end select
      end if

      ! model2 -> model1
      found = .false.
      do i = 1, exg%nexg
        if (particle%icu == exg%nodem2(i)) then
          dest_cell = exg%nodem1(i)
          found = .true.
          exit
        end if
      end do
      if (found) then
        exgprp_obj => GetBndFromList( &
                      exg%prtmodel1%bndlist, &
                      exg%prtmodel1%bndlist%Count())
        select type (exgprp_obj)
        type is (ExgPrtPrpType)
          exgprp => exgprp_obj
          particle%itrdomain(1) = exg%prtmodel1%id
          particle%icu = dest_cell
          ! TODO: update particle coordinates
          exgprp%has_pending = .true.
          exgprp%nparticles = exgprp%nparticles + 1
          call exgprp%particles%resize(exgprp%nparticles, exgprp%memoryPath)
          call exgprp%particles%put(particle, exgprp%nparticles)
          transferred = .true.
          return
        end select
      end if
    end select
  end function try_particle_transfer

  !> @brief Get PRT-PRT exchange from list
  function GetPrtPrtExchangeFromList(list, idx) result(res)
    use ListModule, only: ListType
    implicit none
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(PrtPrtExchangeType), pointer :: res
    class(*), pointer :: obj

    obj => list%GetItem(idx)
    res => CastAsPrtPrtExchange(obj)
  end function GetPrtPrtExchangeFromList

  !> @brief Cast object as PRT-PRT exchange
  function CastAsPrtPrtExchange(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(PrtPrtExchangeType), pointer :: res

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (PrtPrtExchangeType)
      res => obj
    end select
  end function CastAsPrtPrtExchange

end module PrtPrtExchangeModule
