module TrackingModelModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENBUDTXT, LENPAKLOC
  use ExplicitModelModule, only: ExplicitModelType
  use BaseDisModule, only: DisBaseType
  use SparseModule, only: sparsematrix
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType
  use ListModule, only: ListType
  use ParticleModule ! kluge???
  use MethodModule, only: MethodType
  use GlobalDataModule

  implicit none
  private
  public :: TrackingModelType, AddTrackingModelToList, &
            GetTrackingModelFromList

  type, extends(ExplicitModelType) :: TrackingModelType
    character(len=LINELENGTH), pointer :: filename => null() !input file name
    integer(I4B), pointer :: neq => null() !number of equations
    integer(I4B), pointer :: nja => null() !number of connections
    integer(I4B), pointer :: moffset => null() !offset of this model in the solution
    integer(I4B), pointer :: icnvg => null() !convergence flag
    ! integer(I4B), pointer :: npart => null() !number of particles in the model
    ! integer(I4B), pointer :: npartmax => null() !maximum number of particles in the model
    integer(I4B), dimension(:), pointer, contiguous :: ia => null() !csr row pointer
    integer(I4B), dimension(:), pointer, contiguous :: ja => null() !csr columns
    real(DP), dimension(:), pointer, contiguous :: x => null() !dependent variable (head, conc, etc)
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !right-hand side vector
    real(DP), dimension(:), pointer, contiguous :: cond => null() !conductance matrix
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() !pointer to position in solution matrix
    real(DP), dimension(:), pointer, contiguous :: xold => null() !dependent variable for previous timestep
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !intercell flows
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !ibound array
    !
    ! -- Derived types
    ! type(ParticleListType), dimension(:), pointer   :: partlist => null()        !list of particle data
    type(ParticleListType), pointer :: partlist => null() !list of particle data

  contains
    !
    ! -- Required for all models (override procedures defined in BaseModelType)
    procedure :: model_df
    procedure :: model_ar
    procedure :: model_fp
    procedure :: model_da
    !
    ! -- Methods specific to a tracking model
    ! procedure :: model_ac    ! kluge note: not needed???
    ! procedure :: model_mc    ! kluge note: not needed???
    procedure :: model_rp
    ! procedure :: model_ad
    procedure :: model_cf
    procedure :: model_fc
    procedure :: model_ptcchk
    procedure :: model_ptc
    procedure :: model_nr
    procedure :: model_cc
    procedure :: model_nur
    ! procedure :: model_cq
    ! procedure :: model_bd
    procedure :: model_bdcalc
    procedure :: model_bdsave
    procedure :: model_ot
    procedure :: model_bdentry
    ! procedure :: model_solve
    !
    ! -- Utility methods
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: set_moffset
    procedure :: set_idsoln
    procedure :: set_xptr
    procedure :: set_rhsptr
    procedure :: set_iboundptr
    procedure :: get_mrange
    procedure :: get_mcellid
    procedure :: get_mnodeu
    procedure :: get_iasym
    procedure :: get_method
  end type TrackingModelType

contains
  !
  ! -- Type-bound procedures for a tracking model
  !
  subroutine model_df(this)
    class(TrackingModelType) :: this
  end subroutine model_df

  ! subroutine model_ac(this, sparse)             ! kluge note: not needed???
  !   class(TrackingModelType) :: this
  !   type(sparsematrix), intent(inout) :: sparse
  ! end subroutine model_ac
  !
  ! subroutine model_mc(this, iasln, jasln)       ! kluge note: not needed???
  !   class(TrackingModelType) :: this
  !   integer(I4B), dimension(:), intent(in) :: iasln
  !   integer(I4B), dimension(:), intent(in) :: jasln
  ! end subroutine model_mc
  !
  subroutine model_ar(this)
    class(TrackingModelType) :: this
  end subroutine model_ar

  subroutine model_rp(this)
    class(TrackingModelType) :: this
  end subroutine model_rp

  ! subroutine model_ad(this)
  !   class(TrackingModelType) :: this
  ! end subroutine model_ad

  subroutine model_cf(this, kiter)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: kiter
  end subroutine model_cf

  subroutine model_fc(this, kiter, amatsln, njasln, inwtflag)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: njasln
    real(DP), dimension(njasln), intent(inout) :: amatsln
    integer(I4B), intent(in) :: inwtflag
  end subroutine model_fc

  subroutine model_ptcchk(this, iptc)
    class(TrackingModelType) :: this
    integer(I4B), intent(inout) :: iptc
    iptc = 0
  end subroutine model_ptcchk

  subroutine model_ptc(this, kiter, neqsln, njasln, &
                       ia, ja, x, rhs, amatsln, iptc, ptcf)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: neqsln
    integer(I4B), intent(in) :: njasln
    integer(I4B), dimension(neqsln + 1), intent(in) :: ia
    integer(I4B), dimension(njasln), intent(in) :: ja
    real(DP), dimension(neqsln), intent(in) :: x
    real(DP), dimension(neqsln), intent(in) :: rhs
    real(DP), dimension(njasln), intent(in) :: amatsln
    integer(I4B), intent(inout) :: iptc
    real(DP), intent(inout) :: ptcf
  end subroutine model_ptc

  subroutine model_nr(this, kiter, amatsln, njasln, inwtflag)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: njasln
    real(DP), dimension(njasln), intent(inout) :: amatsln
    integer(I4B), intent(in) :: inwtflag
  end subroutine model_nr

  subroutine model_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
  end subroutine model_cc

  subroutine model_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: neqmod
    real(DP), dimension(neqmod), intent(inout) :: x
    real(DP), dimension(neqmod), intent(in) :: xtemp
    real(DP), dimension(neqmod), intent(inout) :: dx
    integer(I4B), intent(inout) :: inewtonur
    real(DP), intent(inout) :: dxmax
    integer(I4B), intent(inout) :: locmax
  end subroutine model_nur

  ! subroutine model_cq(this, icnvg, isuppress_output)
  !   class(TrackingModelType) :: this
  !   integer(I4B),intent(in) :: icnvg
  !   integer(I4B), intent(in) :: isuppress_output
  ! end subroutine model_cq
  !
  ! subroutine model_bd(this, icnvg, isuppress_output)
  !   class(TrackingModelType) :: this
  !   integer(I4B),intent(in) :: icnvg
  !   integer(I4B), intent(in) :: isuppress_output
  ! end subroutine model_bd

  subroutine model_bdcalc(this, icnvg)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: icnvg
  end subroutine model_bdcalc

  subroutine model_bdsave(this, icnvg)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: icnvg
  end subroutine model_bdsave

  subroutine model_ot(this)
    class(TrackingModelType) :: this
  end subroutine model_ot

  subroutine model_bdentry(this, budterm, budtxt, rowlabel)
    class(TrackingModelType) :: this
    real(DP), dimension(:, :), intent(in) :: budterm
    character(len=LENBUDTXT), dimension(:), intent(in) :: budtxt
    character(len=*), intent(in) :: rowlabel
  end subroutine model_bdentry

  subroutine model_fp(this)
    class(TrackingModelType) :: this
  end subroutine model_fp

  subroutine model_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    class(TrackingModelType) :: this

    ! -- Scalars
    call mem_deallocate(this%neq)
    call mem_deallocate(this%nja)
    call mem_deallocate(this%icnvg)
    call mem_deallocate(this%moffset)
    deallocate (this%filename)
    ! call mem_deallocate(this%npart)
    ! call mem_deallocate(this%npartmax)
    !
    ! -- Arrays
    call mem_deallocate(this%xold)
    call mem_deallocate(this%flowja)
    call mem_deallocate(this%idxglo)
    ! deallocate(this%partlist)        ! kluge
    ! deallocate(this%partlist%velmult)
    ! deallocate(this%partlist%x)   ! kluge note: use mem_deallocate for these arrays
    ! deallocate(this%partlist%y)
    ! deallocate(this%partlist%z)
    ! deallocate(this%partlist%xlocal)
    ! deallocate(this%partlist%ylocal)
    ! deallocate(this%partlist%zlocal)
    ! deallocate(this%partlist%iTrackingDomain)
    ! deallocate(this%partlist%iTrackingDomainBoundary)
    ! deallocate(this%partlist%trelease)
    ! deallocate(this%partlist%tstop)
    ! deallocate(this%partlist%ttrack)
    ! deallocate(this%partlist%istopweaksink)
    ! deallocate(this%partlist%istopzone)
    ! deallocate(this%partlist%istatus)
    ! deallocate(this%partlist)      ! kluge note: structure of arrays - deallocate pointer elsewhere???
    call mem_deallocate(this%ibound)
    !
    ! -- nullify pointers
    ! call mem_deallocate(this%x, 'X', this%memoryPath)
    ! call mem_deallocate(this%rhs, 'RHS', this%memoryPath)
    if (associated(this%ibound)) &
      call mem_deallocate(this%ibound, 'IBOUND', this%memoryPath)
    !
    ! -- ExplicitModelType
    call this%ExplicitModelType%model_da()
    !
    !
    ! -- Return
    return
  end subroutine model_da

  ! subroutine model_solve(this)
  !   class(TrackingModelType), intent(inout) :: this
  ! end subroutine model_solve

  subroutine set_moffset(this, moffset)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: moffset
    this%moffset = moffset
  end subroutine set_moffset

  subroutine get_mrange(this, mstart, mend)
    class(TrackingModelType) :: this
    integer(I4B), intent(inout) :: mstart
    integer(I4B), intent(inout) :: mend
    mstart = this%moffset + 1
    mend = mstart + this%neq - 1
  end subroutine get_mrange

  subroutine set_idsoln(this, id)
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: id
    this%idsoln = id
  end subroutine set_idsoln

  subroutine allocate_scalars(this, modelname)
    use MemoryManagerModule, only: mem_allocate
    class(TrackingModelType) :: this
    character(len=*), intent(in) :: modelname
    !
    ! -- allocate basetype members
    call this%ExplicitModelType%allocate_scalars(modelname)
    !
    ! -- allocate members from this type
    call mem_allocate(this%neq, 'NEQ', this%memoryPath)
    call mem_allocate(this%nja, 'NJA', this%memoryPath)
    call mem_allocate(this%icnvg, 'ICNVG', this%memoryPath)
    call mem_allocate(this%moffset, 'MOFFSET', this%memoryPath)
    allocate (this%filename)
    ! call mem_allocate(this%npart, 'NPART', this%memoryPath)
    ! call mem_allocate(this%npartmax, 'NPARTMAX', this%memoryPath)
    !
    this%filename = ''
    this%neq = 0
    this%nja = 0
    this%icnvg = 0
    this%moffset = 0
    ! this%npart = 0
    ! this%npartmax = 0
    !
    ! -- return
    return
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
    use ConstantsModule, only: DZERO
    use MemoryManagerModule, only: mem_allocate
    class(TrackingModelType) :: this
    integer(I4B) :: i
    integer(I4B) :: np
    !
    call mem_allocate(this%xold, this%neq, 'XOLD', this%memoryPath)
    call mem_allocate(this%flowja, this%nja, 'FLOWJA', this%memoryPath)
    call mem_allocate(this%idxglo, this%nja, 'IDXGLO', this%memoryPath)
    ! call mem_allocate(this%partlist, this%npartmax, 'PARTLIST', this%memoryPath)   ! kluge note: update memory manager for this derived type???
    ! allocate(this%partlist(this%npartmax))                                          ! kluge
    ! do np=1,this%npartmax
    !   call create_particle(this%partlist(np)%particle)
    ! end do
    ! allocate(this%partlist)         ! kluge note: structure of arrays
    ! allocate(this%partlist%velmult(this%npartmax))
    ! allocate(this%partlist%x(this%npartmax))   ! kluge note: nprtmax is the initial max dimension
    ! allocate(this%partlist%y(this%npartmax))     ! kluge note: use mem_allocate for these arrays
    ! allocate(this%partlist%z(this%npartmax))
    ! allocate(this%partlist%xlocal(this%npartmax))
    ! allocate(this%partlist%ylocal(this%npartmax))
    ! allocate(this%partlist%zlocal(this%nreleasepts))
    ! allocate(this%partlist%iTrackingDomain(this%npartmax,levelMin:levelMax))         ! kluge note: ditch crazy dims
    ! allocate(this%partlist%iTrackingDomainBoundary(this%npartmax,levelMin:levelMax)) ! kluge note: ditch crazy dims
    ! allocate(this%partlist%trelease(this%npartmax))
    ! allocate(this%partlist%tstop(this%npartmax))
    ! allocate(this%partlist%ttrack(this%npartmax))
    ! allocate(this%partlist%istopweaksink(this%npartmax))
    ! allocate(this%partlist%istopzone(this%npartmax))
    ! allocate(this%partlist%istatus(this%npartmax))
    call mem_allocate(this%ibound, this%neq, 'IBOUND', this%memoryPath)
    !
    ! -- initialize
    do i = 1, size(this%flowja)
      this%flowja(i) = DZERO
    end do
    do i = 1, this%neq
      this%ibound(i) = 1 !default is active
    end do
    !
    ! -- return
    return
  end subroutine allocate_arrays

  subroutine set_xptr(this, xsln, varNameTgt, memPathTgt)
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(TrackingModelType) :: this
    real(DP), dimension(:), pointer, contiguous, intent(in) :: xsln
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    ! -- code
    this%x => xsln(this%moffset + 1:this%moffset + this%neq)
    call mem_checkin(this%x, 'X', this%memoryPath, varNameTgt, memPathTgt)
  end subroutine set_xptr

  subroutine set_rhsptr(this, rhssln, varNameTgt, memPathTgt)
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(TrackingModelType) :: this
    real(DP), dimension(:), pointer, contiguous, intent(in) :: rhssln
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    ! -- code
    this%rhs => rhssln(this%moffset + 1:this%moffset + this%neq)
    call mem_checkin(this%rhs, 'RHS', this%memoryPath, varNameTgt, memPathTgt)
  end subroutine set_rhsptr

  subroutine set_iboundptr(this, iboundsln, varNameTgt, memPathTgt) ! kluge note: not needed???
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(TrackingModelType) :: this
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: iboundsln
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    ! -- code
    this%ibound => iboundsln(this%moffset + 1:this%moffset + this%neq)
  call mem_checkin(this%ibound, 'IBOUND', this%memoryPath, varNameTgt, memPathTgt)
  end subroutine set_iboundptr

  subroutine get_mcellid(this, node, mcellid)
    use BndModule, only: BndType, GetBndFromList
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: node
    character(len=*), intent(inout) :: mcellid
    ! -- local
    character(len=20) :: cellid
    integer(I4B) :: ip, ipaknode, istart, istop
    class(BndType), pointer :: packobj

    if (node < 1) then
      cellid = ''
    else if (node <= this%dis%nodes) then
      call this%dis%noder_to_string(node, cellid)
    else
      cellid = '***ERROR***'
      ipaknode = node - this%dis%nodes
      istart = 1
      do ip = 1, this%bndlist%Count()
        packobj => GetBndFromList(this%bndlist, ip)
        if (packobj%npakeq == 0) cycle
        istop = istart + packobj%npakeq - 1
        if (istart <= ipaknode .and. ipaknode <= istop) then
          write (cellid, '(a, a, a, i0, a, i0, a)') '(', &
            trim(packobj%filtyp), '_', &
            packobj%ibcnum, '-', ipaknode - packobj%ioffset, ')'
          exit
        end if
        istart = istop + 1
      end do
    end if
    write (mcellid, '(i0, a, a, a, a)') this%id, '_', this%macronym, '-', &
      trim(adjustl(cellid))
    return
  end subroutine get_mcellid

  subroutine get_mnodeu(this, node, nodeu)
    use BndModule, only: BndType, GetBndFromList
    class(TrackingModelType) :: this
    integer(I4B), intent(in) :: node
    integer(I4B), intent(inout) :: nodeu
    ! -- local
    if (node <= this%dis%nodes) then
      nodeu = this%dis%get_nodeuser(node)
    else
      nodeu = -(node - this%dis%nodes)
    end if
    return
  end subroutine get_mnodeu

  function get_iasym(this) result(iasym)
    class(TrackingModelType) :: this
    integer(I4B) :: iasym
    iasym = 0
  end function get_iasym

  function get_method(this, particle) result(method)
    class(TrackingModelType) :: this
    type(ParticleType), pointer :: particle
    class(MethodType), pointer :: method
    method => null()
  end function get_method

  function CastAsTrackingModelClass(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(TrackingModelType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (TrackingModelType)
      res => obj
    end select
    return
  end function CastAsTrackingModelClass

  subroutine AddTrackingModelToList(list, model)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    class(TrackingModelType), pointer, intent(inout) :: model
    ! -- local
    class(*), pointer :: obj
    !
    obj => model
    call list%Add(obj)
    !
    return
  end subroutine AddTrackingModelToList

  function GetTrackingModelFromList(list, idx) result(res)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(TrackingModelType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsTrackingModelClass(obj)
    !
    return
  end function GetTrackingModelFromList

end module TrackingModelModule
