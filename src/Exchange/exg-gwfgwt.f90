module GwfGwtExchangeModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME, LINELENGTH
  use ListsModule, only: basemodellist, baseexchangelist, &
                         baseconnectionlist
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          get_smc_from_list
  use GwtGwtConnectionModule, only: GwtGwtConnectionType, CastAsGwtGwtConnection
  use GwfGwfConnectionModule, only: GwfGwfConnectionType, CastAsGwfGwfConnection
  use GwfGwfExchangeModule, only: GwfExchangeType, &
                                  GetGwfExchangeFromList
  use GwtGwtExchangeModule, only: GwtExchangeType
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use GwfModule, only: GwfModelType
  use GwtModule, only: GwtModelType
  use BndModule, only: BndType, GetBndFromList

  implicit none
  public :: GwfGwtExchangeType
  public :: gwfgwt_cr

  type, extends(BaseExchangeType) :: GwfGwtExchangeType

    integer(I4B), pointer :: m1_idx => null() !< index into the list of base exchanges for model 1
    integer(I4B), pointer :: m2_idx => null() !< index into the list of base exchanges for model 2
    character(len=LINELENGTH) :: filename !< the input file for the GWF-GWT exchange

  contains

    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: gwfbnd2gwtfmi
    procedure, private :: gwfconn2gwtconn
    procedure, private :: link_connections

  end type GwfGwtExchangeType

contains

  !> @brief Create a new GWF to GWT exchange object
  !<
  subroutine gwfgwt_cr(filename, id, m1_id, m2_id)
    ! -- modules
    use SimVariablesModule, only: model_loc_idx
    ! -- dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: m1_id
    integer(I4B), intent(in) :: m2_id
    ! -- local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(GwfGwtExchangeType), pointer :: exchange => null()
    character(len=20) :: cint
    !
    ! -- Create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)
    !
    ! -- Assign id and name
    exchange%id = id
    write (cint, '(i0)') id
    exchange%name = 'GWF-GWT_'//trim(adjustl(cint))
    exchange%memoryPath = exchange%name
    exchange%filename = filename
    !
    ! -- allocate scalars
    call exchange%allocate_scalars()
    !
    ! -- NB: convert from id to local model index in base model list
    exchange%m1_idx = model_loc_idx(m1_id)
    exchange%m2_idx = model_loc_idx(m2_id)
    !
    ! -- set model pointers
    call exchange%set_model_pointers()
  end subroutine gwfgwt_cr

  !> @brief Allocate and read
  !<
  subroutine set_model_pointers(this)
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    !
    ! -- set gwfmodel
    gwfmodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    gwtmodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Verify that gwf model is of the correct type
    if (.not. associated(gwfmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that gwt model is of the correct type
    if (.not. associated(gwtmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  Specified GWT Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Tell transport model fmi flows are not read from file
    gwtmodel%fmi%flows_from_file = .false.
    !
    ! -- Set a pointer to the GWF bndlist.  This will allow the transport model
    !    to look through the flow packages and establish a link to GWF flows
    gwtmodel%fmi%gwfbndlist => gwfmodel%bndlist
  end subroutine set_model_pointers

  !> @brief Define the GwfGwt Exchange object
  !<
  subroutine exg_df(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Check to make sure that flow is solved before transport and in a
    !    different IMS solution
    if (gwfmodel%idsoln >= gwtmodel%idsoln) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  The GWF model must be solved by a different IMS than the GWT model. &
        &Furthermore, the IMS specified for GWF must be listed in mfsim.nam &
        &before the IMS for GWT.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Set pointer to flowja
    gwtmodel%fmi%gwfflowja => gwfmodel%flowja
    call mem_checkin(gwtmodel%fmi%gwfflowja, &
                     'GWFFLOWJA', gwtmodel%fmi%memoryPath, &
                     'FLOWJA', gwfmodel%memoryPath)

    !
    ! -- Set the npf flag so that specific discharge is available for
    !    transport calculations if dispersion is active
    if (gwtmodel%indsp > 0) then
      gwfmodel%npf%icalcspdis = 1
    end if
  end subroutine exg_df

  !> @brief Allocate and read
  !<
  subroutine exg_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin, mem_allocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    logical :: subset_domain
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and GWT Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  GWT Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  Ensure discretization packages have the same shape.')"
    character(len=*), parameter :: fmtidomerr = &
      "('GWF and GWT Models do not have compatible discretizations for &
      &exchange ',a,'.&
      &  A cell that is inactive (IDOMAIN <= 0) in the GWF Model is active &
      &in the GWT Model.&
      &  The GWT active domain must be a subset of the GWF active domain.')"
    character(len=*), parameter :: fmtbuyvscerr = &
      "('GWF-GWT exchange ',a,': GWT''s active domain is a subset of &
      &GWF''s, but the GWF model has BUY or VSC active.  Coupling BUY/VSC &
      &to a GWT model with a smaller active domain is not yet supported.')"
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Check that the two models share the same underlying grid (same
    !    number of user nodes). The number of *active* (reduced) nodes may
    !    differ: the GWT active domain is allowed to be a subset of the
    !    GWF active domain (i.e. GWT's IDOMAIN may mark inactive cells that
    !    are active in GWF, but not the other way around).
    if (gwtmodel%dis%nodesuser /= gwfmodel%dis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), &
        gwfmodel%dis%nodesuser, &
        gwfmodel%dis%nodes, &
        gwtmodel%dis%nodesuser, &
        gwtmodel%dis%nodes
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- Check that GWT's active domain is a subset of GWF's, and if the
    !    two domains actually differ, build a node/connection map to
    !    translate between the two models' reduced numbering.
    call check_and_map_domains(this, gwfmodel, gwtmodel, fmtidomerr, &
                               subset_domain)
    !
    ! -- BUY/VSC read GWT's concentration/temperature using GWF's own node
    !    numbering; that coupling isn't mapped yet, so reject it rather
    !    than silently producing wrong buoyancy/viscosity effects.
    if (subset_domain .and. &
        (gwfmodel%inbuy > 0 .or. gwfmodel%invsc > 0)) then
      write (errmsg, fmtbuyvscerr) trim(this%name)
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- setup pointers to gwf variables allocated in gwf_ar
    if (subset_domain) then
      ! -- gwfflowja was already aliased to gwfmodel%flowja and checked in
      !    as 'GWFFLOWJA' in exg_df; that alias becomes the raw source, and
      !    gwfflowja itself becomes an owned, GWT-sized array, refreshed
      !    once per fmi_ad (see TspFmiType%translate_gwf_arrays)
      gwtmodel%fmi%gwfflowja_raw => gwtmodel%fmi%gwfflowja
      nullify (gwtmodel%fmi%gwfflowja)
      call mem_allocate(gwtmodel%fmi%gwfflowja, gwtmodel%dis%con%nja, &
                        'GWFFLOWJA_LOC', gwtmodel%fmi%memoryPath)
      !
      gwtmodel%fmi%gwfhead_raw => gwfmodel%x
      call mem_checkin(gwtmodel%fmi%gwfhead_raw, &
                       'GWFHEAD_RAW', gwtmodel%fmi%memoryPath, &
                       'X', gwfmodel%memoryPath)
      call mem_allocate(gwtmodel%fmi%gwfhead, gwtmodel%dis%nodes, &
                        'GWFHEAD', gwtmodel%fmi%memoryPath)
      !
      gwtmodel%fmi%gwfsat_raw => gwfmodel%npf%sat
      call mem_checkin(gwtmodel%fmi%gwfsat_raw, &
                       'GWFSAT_RAW', gwtmodel%fmi%memoryPath, &
                       'SAT', gwfmodel%npf%memoryPath)
      call mem_allocate(gwtmodel%fmi%gwfsat, gwtmodel%dis%nodes, &
                        'GWFSAT', gwtmodel%fmi%memoryPath)
      !
      gwtmodel%fmi%gwfspdis_raw => gwfmodel%npf%spdis
      call mem_checkin(gwtmodel%fmi%gwfspdis_raw, &
                       'GWFSPDIS_RAW', gwtmodel%fmi%memoryPath, &
                       'SPDIS', gwfmodel%npf%memoryPath)
      call mem_allocate(gwtmodel%fmi%gwfspdis, 3, gwtmodel%dis%nodes, &
                        'GWFSPDIS', gwtmodel%fmi%memoryPath)
    else
      gwtmodel%fmi%gwfhead => gwfmodel%x
      call mem_checkin(gwtmodel%fmi%gwfhead, &
                       'GWFHEAD', gwtmodel%fmi%memoryPath, &
                       'X', gwfmodel%memoryPath)
      gwtmodel%fmi%gwfsat => gwfmodel%npf%sat
      call mem_checkin(gwtmodel%fmi%gwfsat, &
                       'GWFSAT', gwtmodel%fmi%memoryPath, &
                       'SAT', gwfmodel%npf%memoryPath)
      gwtmodel%fmi%gwfspdis => gwfmodel%npf%spdis
      call mem_checkin(gwtmodel%fmi%gwfspdis, &
                       'GWFSPDIS', gwtmodel%fmi%memoryPath, &
                       'SPDIS', gwfmodel%npf%memoryPath)
    end if
    gwtmodel%fmi%igwfspdis = gwfmodel%npf%icalcspdis
    !
    ! -- setup pointers to the flow storage rates. GWF strg arrays are
    !    available after the gwf_ar routine is called.
    if (gwtmodel%inmst > 0) then
      if (gwfmodel%insto > 0) then
        gwtmodel%fmi%igwfstrgss = 1
        if (subset_domain) then
          gwtmodel%fmi%gwfstrgss_raw => gwfmodel%sto%strgss
          call mem_allocate(gwtmodel%fmi%gwfstrgss, gwtmodel%dis%nodes, &
                            'GWFSTRGSS', gwtmodel%fmi%memoryPath)
        else
          gwtmodel%fmi%gwfstrgss => gwfmodel%sto%strgss
        end if
        if (gwfmodel%sto%iusesy == 1) then
          gwtmodel%fmi%igwfstrgsy = 1
          if (subset_domain) then
            gwtmodel%fmi%gwfstrgsy_raw => gwfmodel%sto%strgsy
            call mem_allocate(gwtmodel%fmi%gwfstrgsy, gwtmodel%dis%nodes, &
                              'GWFSTRGSY', gwtmodel%fmi%memoryPath)
          else
            gwtmodel%fmi%gwfstrgsy => gwfmodel%sto%strgsy
          end if
        end if
      end if
    end if
    !
    ! -- Set a pointer to conc in buy
    if (gwfmodel%inbuy > 0) then
      call gwfmodel%buy%set_concentration_pointer(gwtmodel%name, gwtmodel%x, &
                                                  gwtmodel%ibound)
    end if
    !
    ! -- Set a pointer to conc (which could be a temperature) in vsc
    if (gwfmodel%invsc > 0) then
      call gwfmodel%vsc%set_concentration_pointer(gwtmodel%name, gwtmodel%x, &
                                                  gwtmodel%ibound)
    end if
    !
    ! -- transfer the boundary package information from gwf to gwt
    call this%gwfbnd2gwtfmi()
    !
    ! -- if mover package is active, then set a pointer to it's budget object
    if (gwfmodel%inmvr /= 0) then
      gwtmodel%fmi%mvrbudobj => gwfmodel%mvr%budobj
    end if
    !
    ! -- connect Connections
    call this%gwfconn2gwtconn(gwfmodel, gwtmodel)
  end subroutine exg_ar

  !> @brief Verify GWT's active domain is a subset of GWF's, and if the two
  !! domains actually differ, build the node and connection maps needed to
  !! translate between GWF's and GWT's reduced numbering.
  !<
  subroutine check_and_map_domains(this, gwfmodel, gwtmodel, fmtidomerr, &
                                   subset_domain)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    type(GwfModelType), pointer, intent(in) :: gwfmodel
    type(GwtModelType), pointer, intent(in) :: gwtmodel
    character(len=*), intent(in) :: fmtidomerr
    logical, intent(out) :: subset_domain
    ! -- local
    integer(I4B) :: nu, gn, pn, gm, pm, ipos, jpos
    !
    ! -- GWT's active domain must be a subset of GWF's: every user node
    !    active in GWT must also be active in GWF.
    subset_domain = (gwtmodel%dis%nodes /= gwfmodel%dis%nodes)
    do nu = 1, gwtmodel%dis%nodesuser
      if (gwtmodel%dis%get_nodenumber(nu, 0) /= 0 .and. &
          gwfmodel%dis%get_nodenumber(nu, 0) == 0) then
        write (errmsg, fmtidomerr) trim(this%name)
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end do
    !
    ! -- If the domains are identical, GWF's and GWT's reduced numbering
    !    coincide and no map is needed: leave fmi%gwf2loc/loc2gwf/loc2gwfja
    !    unassociated so downstream code takes the (existing) direct-index
    !    path.
    if (.not. subset_domain) return
    !
    ! -- Build the node maps
    call mem_allocate(gwtmodel%fmi%gwf2loc, gwfmodel%dis%nodes, &
                      'GWF2LOC', gwtmodel%fmi%memoryPath)
    call mem_allocate(gwtmodel%fmi%loc2gwf, gwtmodel%dis%nodes, &
                      'LOC2GWF', gwtmodel%fmi%memoryPath)
    do gn = 1, gwfmodel%dis%nodes
      nu = gwfmodel%dis%get_nodeuser(gn)
      gwtmodel%fmi%gwf2loc(gn) = gwtmodel%dis%get_nodenumber(nu, 0)
    end do
    do pn = 1, gwtmodel%dis%nodes
      nu = gwtmodel%dis%get_nodeuser(pn)
      gwtmodel%fmi%loc2gwf(pn) = gwfmodel%dis%get_nodenumber(nu, 0)
    end do
    !
    ! -- Build the connection map: for each GWT reduced connection position,
    !    find the corresponding position in GWF's ia/ja. This relies on the
    !    subset property just verified: every connection between two cells
    !    active in GWT is also a connection between the same two cells (in
    !    GWF's numbering) in GWF's larger active domain.
    call mem_allocate(gwtmodel%fmi%loc2gwfja, gwtmodel%dis%con%nja, &
                      'LOC2GWFJA', gwtmodel%fmi%memoryPath)
    do pn = 1, gwtmodel%dis%nodes
      gn = gwtmodel%fmi%loc2gwf(pn)
      ! -- diagonal position maps directly to GWF's diagonal position
      gwtmodel%fmi%loc2gwfja(gwtmodel%dis%con%ia(pn)) = gwfmodel%dis%con%ia(gn)
      do ipos = gwtmodel%dis%con%ia(pn) + 1, gwtmodel%dis%con%ia(pn + 1) - 1
        pm = gwtmodel%dis%con%ja(ipos)
        gm = gwtmodel%fmi%loc2gwf(pm)
        gwtmodel%fmi%loc2gwfja(ipos) = 0
        do jpos = gwfmodel%dis%con%ia(gn) + 1, gwfmodel%dis%con%ia(gn + 1) - 1
          if (gwfmodel%dis%con%ja(jpos) == gm) then
            gwtmodel%fmi%loc2gwfja(ipos) = jpos
            exit
          end if
        end do
        if (gwtmodel%fmi%loc2gwfja(ipos) == 0) then
          write (errmsg, '(a,a)') 'Programmer error: could not map a GWT &
            &connection onto the corresponding GWF connection for &
            &exchange ', trim(this%name)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end do
    end do
    !
    ! -- Warn about a real numerical caveat: unlike PRT (where "no exit
    !    face" at the edge of a subset domain is already a sound, existing
    !    concept), GWT's advection scheme assumes each active cell's water
    !    balance closes over its own connections. Where a GWT-active cell
    !    borders a cell that's active in GWF but excluded from GWT, GWF's
    !    real flux across that face is silently dropped rather than folded
    !    in as a boundary term, which can break that per-cell balance and
    !    produce local concentration artifacts (observed: unphysical
    !    overshoot at such a boundary cell in ad hoc testing). Warn rather
    !    than block, since a user excluding a genuinely low/no-flow region
    !    (the common ICBUND use case) may see no ill effect in practice.
    call warn_dropped_boundary_flux(this, gwfmodel, gwtmodel)
  end subroutine check_and_map_domains

  !> @brief Warn if any GWT-active cell borders a cell that's active in GWF
  !! but excluded from GWT -- see the caveat noted in check_and_map_domains.
  !<
  subroutine warn_dropped_boundary_flux(this, gwfmodel, gwtmodel)
    ! -- modules
    use SimVariablesModule, only: iout
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    type(GwfModelType), pointer, intent(in) :: gwfmodel
    type(GwtModelType), pointer, intent(in) :: gwtmodel
    ! -- local
    integer(I4B) :: gn, gm, ipos, ndropped
    character(len=*), parameter :: fmtwarn = &
      "(/1x,'WARNING: exchange ',a,': ',i0,' connection(s) between a GWT-&
      &active cell and a cell active in GWF but excluded from GWT were &
      &found. GWF flow across these connections is not represented in &
      &GWT''s transport equation, which can produce local mass-balance &
      &artifacts near the excluded region (concentrations outside the &
      &physically expected range). This is most likely to matter where &
      &GWF flow across the excluded boundary is not small.'/)"
    !
    ndropped = 0
    do gn = 1, gwfmodel%dis%nodes
      if (gwtmodel%fmi%gwf2loc(gn) <= 0) cycle
      do ipos = gwfmodel%dis%con%ia(gn) + 1, gwfmodel%dis%con%ia(gn + 1) - 1
        gm = gwfmodel%dis%con%ja(ipos)
        if (gwtmodel%fmi%gwf2loc(gm) <= 0) ndropped = ndropped + 1
      end do
    end do
    if (ndropped > 0) then
      write (iout, fmtwarn) trim(this%name), ndropped
    end if
  end subroutine warn_dropped_boundary_flux

  !> @brief Link GWT connections to GWF connections or exchanges
  !<
  subroutine gwfconn2gwtconn(this, gwfModel, gwtModel)
    ! -- modules
    use SimModule, only: store_error, store_error_filename, count_errors
    use SimVariablesModule, only: iout
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this !< this exchange
    type(GwfModelType), pointer :: gwfModel !< the flow model
    type(GwtModelType), pointer :: gwtModel !< the transport model
    ! -- local
    class(SpatialModelConnectionType), pointer :: conn => null()
    class(*), pointer :: objPtr => null()
    class(GwtGwtConnectionType), pointer :: gwtConn => null()
    class(GwfGwfConnectionType), pointer :: gwfConn => null()
    class(GwfExchangeType), pointer :: gwfExg => null()
    class(GwtExchangeType), pointer :: gwtExg => null()
    integer(I4B) :: ic1, ic2, iex
    integer(I4B) :: gwfConnIdx, gwfExIdx
    logical(LGP) :: areEqual
    !
    ! loop over all connections
    gwtloop: do ic1 = 1, baseconnectionlist%Count()
      !
      conn => get_smc_from_list(baseconnectionlist, ic1)
      if (.not. associated(conn%owner, gwtModel)) cycle gwtloop
      !
      ! start with a GWT conn.
      objPtr => conn
      gwtConn => CastAsGwtGwtConnection(objPtr)
      gwtExg => gwtConn%gwtExchange
      gwfConnIdx = -1
      gwfExIdx = -1
      !
      ! find matching GWF conn. in same list
      gwfloop: do ic2 = 1, baseconnectionlist%Count()
        conn => get_smc_from_list(baseconnectionlist, ic2)
        !
        if (associated(conn%owner, gwfModel)) then
          !
          objPtr => conn
          gwfConn => CastAsGwfGwfConnection(objPtr)
          gwfExg => gwfConn%gwfExchange
          !
          ! A model can have multiple exchanges, even connecting the same two
          ! models. We have a match if
          !  1. gwtgwt%model1 is connected to gwfgwf%model1
          !  2. gwtgwt%model2 is connected to gwfgwf%model2
          !  3. the list of connected nodes (nodem1, nodem2) is equivalent, such
          !     that it contains the same nodes, appearing in the same order in the
          !     exchange data block
          !
          if (gwfExg%v_model1%name /= gwtExg%gwfmodelname1) cycle
          if (gwfExg%v_model2%name /= gwtExg%gwfmodelname2) cycle
          !
          areEqual = (gwfExg%nexg == gwtExg%nexg)
          if (areEqual) then
            areEqual = all(gwfExg%nodem1 == gwtExg%nodem1)
            areEqual = areEqual .and. all(gwfExg%nodem2 == gwtExg%nodem2)
          end if
          if (areEqual) then
            ! same DIS, same exchange: link and go to next GWT conn.
            write (iout, '(/6a)') 'Linking exchange ', &
              trim(gwtExg%name), ' to ', trim(gwfExg%name), &
              ' (using interface model) for GWT model ', &
              trim(gwtModel%name)
            gwfConnIdx = ic2
            call this%link_connections(gwtConn, gwfConn)
            exit gwfloop
          end if
        end if
      end do gwfloop
      !
      ! fallback option: coupling to old gwfgwf exchange,
      ! the conditions are equal to what is used above
      ! (this will go obsolete at some point)
      if (gwfConnIdx == -1) then
        gwfloopexg: do iex = 1, baseexchangelist%Count()
          gwfExg => GetGwfExchangeFromList(baseexchangelist, iex)
          !
          if (.not. associated(gwfExg)) cycle gwfloopexg
          !
          if (associated(gwfExg%model1, gwfModel) .or. &
              associated(gwfExg%model2, gwfModel)) then
            !
            if (gwfExg%v_model1%name /= gwtExg%gwfmodelname1) cycle
            if (gwfExg%v_model2%name /= gwtExg%gwfmodelname2) cycle
            !
            areEqual = (gwfExg%nexg == gwtExg%nexg)
            !
            if (areEqual) then
              areEqual = all(gwfExg%nodem1 == gwtExg%nodem1)
              areEqual = areEqual .and. all(gwfExg%nodem2 == gwtExg%nodem2)
            end if
            if (areEqual) then
              ! link exchange to connection
              write (iout, '(/6a)') 'Linking exchange ', &
                trim(gwtExg%name), ' to ', trim(gwfExg%name), ' for GWT model ', &
                trim(gwtModel%name)
              gwfExIdx = iex
              if (gwtConn%owns_exchange) then
                gwtExg%gwfsimvals => gwfExg%simvals
                call mem_checkin(gwtExg%gwfsimvals, &
                                 'GWFSIMVALS', gwtExg%memoryPath, &
                                 'SIMVALS', gwfExg%memoryPath)
              end if
              !
              !cdl link up mvt to mvr
              if (gwfExg%inmvr > 0) then
                if (gwtConn%owns_exchange .and. &
                    gwtConn%gwtExchange%inmvt > 0) then
                  call gwtExg%mvt%set_pointer_mvrbudobj(gwfExg%mvr%budobj)
                end if
              end if
              !
              if (associated(gwfExg%model2, gwfModel)) gwtConn%exgflowSign = -1
              gwtConn%gwtInterfaceModel%fmi%flows_from_file = .false.
              !
              exit gwfloopexg
            end if
          end if
          !
        end do gwfloopexg
      end if
      !
      if (gwfConnIdx == -1 .and. gwfExIdx == -1) then
        ! none found, report
        write (errmsg, *) 'Cannot find GWF-GWF exchange when connecting'// &
          ' GWT model ', trim(gwtModel%name), ' with exchange ', &
          trim(gwtExg%name), ' to GWF model ', trim(gwfModel%name), &
          '. Note: GWF-GWF and GWT-GWT need identical exchange data '// &
          '(both in value and order) for the match to succeed.'
        call store_error(errmsg)
      end if
      !
    end do gwtloop
    !
    ! -- report errors
    if (count_errors() > 0) then
      call store_error_filename(this%filename)
    end if
  end subroutine gwfconn2gwtconn

  !> @brief Links a GWT connection to its GWF counterpart
  !<
  subroutine link_connections(this, gwtConn, gwfConn)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this !< this exchange
    class(GwtGwtConnectionType), pointer :: gwtConn !< GWT connection
    class(GwfGwfConnectionType), pointer :: gwfConn !< GWF connection
    !
    !gwtConn%exgflowja => gwfConn%exgflowja
    if (gwtConn%owns_exchange) then
      gwtConn%gwtExchange%gwfsimvals => gwfConn%gwfExchange%simvals
      call mem_checkin(gwtConn%gwtExchange%gwfsimvals, &
                       'GWFSIMVALS', gwtConn%gwtExchange%memoryPath, &
                       'SIMVALS', gwfConn%gwfExchange%memoryPath)
    end if
    !
    !cdl link up mvt to mvr
    if (gwfConn%gwfExchange%inmvr > 0) then
      if (gwtConn%owns_exchange .and. gwtConn%gwtExchange%inmvt > 0) then
        call gwtConn%gwtExchange%mvt%set_pointer_mvrbudobj( &
          gwfConn%gwfExchange%mvr%budobj)
      end if
    end if
    !
    if (associated(gwfConn%gwfExchange%model2, gwfConn%owner)) then
      gwtConn%exgflowSign = -1
    end if
    !
    ! fmi flows are not read from file
    gwtConn%gwtInterfaceModel%fmi%flows_from_file = .false.
    !
    ! set concentration pointer for buoyancy
    ! call gwfConn%gwfInterfaceModel%buy%set_concentration_pointer( &
    !   gwtConn%gwtModel%name, &
    !   gwtConn%conc, &
    !   gwtConn%icbound)
  end subroutine link_connections

  !> @brief Deallocate memory
  !<
  subroutine exg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    !
    call mem_deallocate(this%m1_idx)
    call mem_deallocate(this%m2_idx)
  end subroutine exg_da

  !> @brief Allocate package scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    !
    call mem_allocate(this%m1_idx, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2_idx, 'M2ID', this%memoryPath)
    this%m1_idx = 0
    this%m2_idx = 0
  end subroutine allocate_scalars

  !> @brief Call routines in FMI that will set pointers to the necessary flow
  !! data
  !<
  subroutine gwfbnd2gwtfmi(this)
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    integer(I4B) :: ngwfpack, ip, iterm, imover
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    class(BndType), pointer :: packobj => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Call routines in FMI that will set pointers to the necessary flow
    !    data (SIMVALS and SIMTOMVR) stored within each GWF flow package
    ngwfpack = gwfmodel%bndlist%Count()
    iterm = 1
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwfmodel%bndlist, ip)
      call gwtmodel%fmi%gwfpackages(iterm)%set_pointers( &
        'SIMVALS', &
        packobj%memoryPath, packobj%input_mempath)
      iterm = iterm + 1
      !
      ! -- If a mover is active for this package, then establish a separate
      !    pointer link for the mover flows stored in SIMTOMVR
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        call gwtmodel%fmi%gwfpackages(iterm)%set_pointers( &
          'SIMTOMVR', &
          packobj%memoryPath, packobj%input_mempath)
        iterm = iterm + 1
      end if
    end do
  end subroutine gwfbnd2gwtfmi

end module GwfGwtExchangeModule
