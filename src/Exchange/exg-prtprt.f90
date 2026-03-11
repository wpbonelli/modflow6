!> @brief PRT-PRT Exchange
!!
!! Serial particle handoff between two PRT models sharing a domain
!! boundary.  Each model must correspond to a distinct GWF subdomain
!! that uses the same grid partitioning (1:1 cell alignment).
!!
!! During each time-step the ExplicitSolution drives an outer "handoff"
!! loop:
!!
!!   1. Every PRT model solves its own particles.  Particles that exit
!!      through a model-boundary face are placed in the model's outbox
!!      (status TERM_BOUNDARY) rather than being permanently terminated.
!!
!!   2. After all models have solved, each PrtPrtExchange calls
!!      do_transfer(), which moves outbox particles from each model into
!!      the inbox of its neighbour, updating the cell index using the
!!      connected-cell pair list.
!!
!!   3. The loop repeats while any model has inbox particles (ninbox > 0).
!!
!! Input file format
!! -----------------
!!   BEGIN EXCHANGEDATA
!!     <nodem1> <nodem2>
!!     ...
!!   END EXCHANGEDATA
!!
!! where nodem1 and nodem2 are 1-based user node numbers in the
!! respective models.
!<
module PrtPrtExchangeModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LINELENGTH, LENMEMPATH
  use ListsModule, only: basemodellist, baseexchangelist
  use SimModule, only: store_error, count_errors
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use PrtModule, only: PrtModelType
  use ParticleModule, only: ParticleType, create_particle, &
                            create_particle_store, TERM_BOUNDARY
  use MemoryHelperModule, only: create_mem_path
  use ExgPrtprtInputModule, only: ExgPrtprtParamFoundType

  implicit none
  private

  public :: PrtPrtExchangeType
  public :: prtprt_cr

  !> @brief Exchange that transfers particles between two PRT models.
  type, extends(BaseExchangeType) :: PrtPrtExchangeType

    integer(I4B), pointer :: m1id => null() !< local index of model 1 in basemodellist
    integer(I4B), pointer :: m2id => null() !< local index of model 2 in basemodellist
    type(PrtModelType), pointer :: prtmodel1 => null() !< pointer to model 1
    type(PrtModelType), pointer :: prtmodel2 => null() !< pointer to model 2
    integer(I4B), pointer :: nexg => null()  !< number of connected cell pairs
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() !< user node nos in model 1
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() !< user node nos in model 2

  contains

    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure :: do_transfer !< override base no-op with particle handoff
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: source_dimensions
    procedure, private :: source_data

  end type PrtPrtExchangeType

contains

  ! -----------------------------------------------------------------------
  !> @brief Create a new PRT-PRT exchange and register it.
  subroutine prtprt_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! -- modules
    use SimVariablesModule, only: model_loc_idx
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    character(len=*), intent(in) :: filename      !< exchange input file
    character(len=*), intent(in) :: name          !< exchange name
    integer(I4B), intent(in) :: id               !< exchange id
    integer(I4B), intent(in) :: m1_id            !< global model 1 id
    integer(I4B), intent(in) :: m2_id            !< global model 2 id
    character(len=*), intent(in) :: input_mempath !< IDM memory path for this exchange
    ! -- local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(PrtPrtExchangeType), pointer :: exchange => null()

    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    exchange%id = id
    exchange%name = name
    exchange%memoryPath = create_mem_path(exchange%name)
    exchange%input_mempath = input_mempath

    call exchange%allocate_scalars()

    ! NB: convert global model id to local index in basemodellist
    exchange%m1id = model_loc_idx(m1_id)
    exchange%m2id = model_loc_idx(m2_id)
  end subroutine prtprt_cr

  ! -----------------------------------------------------------------------
  !> @brief Define phase: store model pointers.
  subroutine exg_df(this)
    class(PrtPrtExchangeType) :: this
    call this%set_model_pointers()
  end subroutine exg_df

  ! -----------------------------------------------------------------------
  !> @brief Allocate and read: load the connected-cell list from IDM.
  subroutine exg_ar(this)
    class(PrtPrtExchangeType) :: this
    integer(I4B) :: iout

    iout = this%prtmodel1%iout
    call this%source_dimensions(iout)
    call this%source_data(iout)
  end subroutine exg_ar

  ! -----------------------------------------------------------------------
  !> @brief Transfer particles between models.
  !!
  !! Moves each particle in model 1's outbox into model 2's inbox (and
  !! vice-versa), updating the particle's cell index via the nodem1/nodem2
  !! mapping.  Outboxes are cleared after the transfer.
  !<
  subroutine do_transfer(this)
    class(PrtPrtExchangeType), intent(inout) :: this
    ! -- local
    type(ParticleType), pointer :: particle => null()
    integer(I4B) :: i, k, n1, n2
    integer(I4B) :: nout1, nout2
    integer(I4B) :: old_nin1, old_nin2
    character(len=LINELENGTH) :: mp1_inbox_path, mp2_inbox_path
    character(len=LINELENGTH) :: mp1_outbox_path, mp2_outbox_path

    mp1_inbox_path = create_mem_path(this%prtmodel1%name, 'PRTPRTIN')
    mp2_inbox_path = create_mem_path(this%prtmodel2%name, 'PRTPRTIN')
    mp1_outbox_path = create_mem_path(this%prtmodel1%name, 'PRTPRTOUT')
    mp2_outbox_path = create_mem_path(this%prtmodel2%name, 'PRTPRTOUT')

    nout1 = this%prtmodel1%noutbox
    nout2 = this%prtmodel2%noutbox

    ! Reuse a single particle for transfers
    call create_particle(particle)

    ! --- Transfer model1.outbox → model2.inbox ---
    if (nout1 > 0) then
      old_nin2 = this%prtmodel2%ninbox
      call this%prtmodel2%inbox%resize(old_nin2 + nout1, mp2_inbox_path)
      do i = 1, nout1
        ! Load particle, setting itrdomain(1) and imdl to model2
        call this%prtmodel1%outbox%get(particle, this%prtmodel2%id, &
                                       this%prtmodel1%outbox%iprp(i), i)
        ! Map the exit cell (nodem1(k)) to the entry cell (nodem2(k))
        do k = 1, this%nexg
          if (particle%icu == this%nodem1(k)) then
            particle%icu = this%nodem2(k)
            exit
          end if
        end do
        ! Mark as active so the receiving model will track it
        particle%istatus = 1 ! ACTIVE
        call this%prtmodel2%inbox%put(particle, old_nin2 + i)
      end do
      this%prtmodel2%ninbox = old_nin2 + nout1
    end if

    ! --- Transfer model2.outbox → model1.inbox ---
    if (nout2 > 0) then
      old_nin1 = this%prtmodel1%ninbox
      call this%prtmodel1%inbox%resize(old_nin1 + nout2, mp1_inbox_path)
      do i = 1, nout2
        call this%prtmodel2%outbox%get(particle, this%prtmodel1%id, &
                                       this%prtmodel2%outbox%iprp(i), i)
        do k = 1, this%nexg
          if (particle%icu == this%nodem2(k)) then
            particle%icu = this%nodem1(k)
            exit
          end if
        end do
        particle%istatus = 1 ! ACTIVE
        call this%prtmodel1%inbox%put(particle, old_nin1 + i)
      end do
      this%prtmodel1%ninbox = old_nin1 + nout2
    end if

    ! --- Clear both outboxes ---
    this%prtmodel1%noutbox = 0
    call this%prtmodel1%outbox%resize(0, mp1_outbox_path)
    this%prtmodel2%noutbox = 0
    call this%prtmodel2%outbox%resize(0, mp2_outbox_path)

    call particle%destroy()
    deallocate (particle)

    ! Map the exit cell (nodem1(k)) to the entry cell (nodem2(k))
    n1 = 0; n2 = 0 ! suppress unused variable warning
  end subroutine do_transfer

  ! -----------------------------------------------------------------------
  !> @brief Deallocate memory.
  subroutine exg_da(this)
    use MemoryManagerModule, only: mem_deallocate
    class(PrtPrtExchangeType) :: this

    call mem_deallocate(this%m1id)
    call mem_deallocate(this%m2id)
    call mem_deallocate(this%nexg)
    if (associated(this%nodem1)) then
      call mem_deallocate(this%nodem1, 'NODEM1', this%memoryPath)
    end if
    if (associated(this%nodem2)) then
      call mem_deallocate(this%nodem2, 'NODEM2', this%memoryPath)
    end if
  end subroutine exg_da

  ! -----------------------------------------------------------------------
  ! Private helpers

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
        trim(this%name), ': model 1 is not a PrtModelType.'
      call store_error(errmsg, terminate=.true.)
    end if
    if (.not. associated(this%prtmodel2)) then
      write (errmsg, '(3a)') 'Problem with PRT-PRT exchange ', &
        trim(this%name), ': model 2 is not a PrtModelType.'
      call store_error(errmsg, terminate=.true.)
    end if
  end subroutine set_model_pointers

  subroutine allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtPrtExchangeType) :: this

    call mem_allocate(this%m1id, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2id, 'M2ID', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)

    this%m1id = 0
    this%m2id = 0
    this%nexg = 0
  end subroutine allocate_scalars

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
    integer(I4B) :: iexg, node1, node2

    call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
    call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)

    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)

    do iexg = 1, this%nexg
      ! -- Convert multi-dim cell id to scalar user node number
      if (this%prtmodel1%dis%ndim == 1) then
        node1 = cellidm1(1, iexg)
      elseif (this%prtmodel1%dis%ndim == 2) then
        node1 = get_node(cellidm1(1, iexg), 1, cellidm1(2, iexg), &
                         this%prtmodel1%dis%mshape(1), 1, &
                         this%prtmodel1%dis%mshape(2))
      else
        node1 = get_node(cellidm1(1, iexg), cellidm1(2, iexg), &
                         cellidm1(3, iexg), &
                         this%prtmodel1%dis%mshape(1), &
                         this%prtmodel1%dis%mshape(2), &
                         this%prtmodel1%dis%mshape(3))
      end if
      this%nodem1(iexg) = this%prtmodel1%dis%get_nodenumber(node1, 0)

      if (this%prtmodel2%dis%ndim == 1) then
        node2 = cellidm2(1, iexg)
      elseif (this%prtmodel2%dis%ndim == 2) then
        node2 = get_node(cellidm2(1, iexg), 1, cellidm2(2, iexg), &
                         this%prtmodel2%dis%mshape(1), 1, &
                         this%prtmodel2%dis%mshape(2))
      else
        node2 = get_node(cellidm2(1, iexg), cellidm2(2, iexg), &
                         cellidm2(3, iexg), &
                         this%prtmodel2%dis%mshape(1), &
                         this%prtmodel2%dis%mshape(2), &
                         this%prtmodel2%dis%mshape(3))
      end if
      this%nodem2(iexg) = this%prtmodel2%dis%get_nodenumber(node2, 0)
    end do
  end subroutine source_data

end module PrtPrtExchangeModule
