! Particle Tracking (PRT) Model

module PrtModule

  use KindModule, only: DP, I4B, LGP
  use InputOutputModule, only: ParseLine, upcase, lowcase
  use ConstantsModule, only: LENFTYPE, LENMEMPATH, DZERO, &
                             DONE, LENPAKLOC, LENBUDTXT, MNORMAL
  use VersionModule, only: write_listfile_header
  use TrackingModelModule, only: TrackingModelType
  use ExplicitModelModule, only: ExplicitModelType
  use BaseModelModule, only: BaseModelType
  use BndModule, only: BndType, AddBndToList, GetBndFromList
  use PrtPrpModule, only: PrtPrpType
  use PrtFmiModule, only: PrtFmiType
  use PrtMipModule, only: PrtMipType
  use PrtOcModule, only: PrtOcType
  use PrtObsModule, only: PrtObsType
  use BudgetModule, only: BudgetType
  use ListModule, only: ListType
  use MethodDisModule, only: MethodDisType
  use MethodDisvModule, only: MethodDisvType
  ! use MethodDisuModule, only: MethodDisuType
  use ParticleModule ! kluge
  use MethodModule
  use GlobalDataModule
  use TrackDataModule, only: TrackDataType
  use SimModule, only: count_errors, store_error, store_error_filename

  implicit none

  private
  public :: prt_cr
  public :: PrtModelType
  public :: CastAsPrtModel

  integer(I4B), parameter :: NBDITEMS = 1
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/'         STORAGE'/

  type, extends(TrackingModelType) :: PrtModelType

    type(PrtFmiType), pointer :: fmi => null() ! flow model interface
    type(PrtMipType), pointer :: mip => null() ! model input package
    type(PrtOcType), pointer :: oc => null() ! output control package
    type(PrtObsType), pointer :: obs => null() ! observation package
    type(BudgetType), pointer :: budget => null() ! budget object
    type(MethodDisType), pointer :: methodDis => null() ! method for dis grid        ! kluge?
    type(MethodDisvType), pointer :: methodDisv => null() ! method for disv grid
    ! type(MethodDisuType), pointer :: methodDisu => null() ! method for disu grid
    integer(I4B), pointer :: infmi => null() ! unit number FMI
    integer(I4B), pointer :: inmip => null() ! unit number MIP
    integer(I4B), pointer :: inmvt => null() ! unit number MVT
    integer(I4B), pointer :: inmst => null() ! unit number MST
    integer(I4B), pointer :: inadv => null() ! unit number ADV
    integer(I4B), pointer :: indsp => null() ! unit number DSP
    integer(I4B), pointer :: inssm => null() ! unit number SSM
    integer(I4B), pointer :: inoc => null() ! unit number OC
    integer(I4B), pointer :: inobs => null() ! unit number OBS
    integer(I4B), pointer :: nprp => null() ! number of PRP packages in the model
    real(DP), dimension(:), pointer, contiguous :: masssto => null() !< particle mass storage in cells, new value
    real(DP), dimension(:), pointer, contiguous :: massstoold => null() !< particle mass storage in cells, old value
    real(DP), dimension(:), pointer, contiguous :: ratesto => null() !< particle mass storage rate in cells
    type(TrackDataType), pointer :: trackdata

  contains

    procedure :: model_df => prt_df
    procedure :: model_ar => prt_ar
    procedure :: model_rp => prt_rp
    procedure :: model_ad => prt_ad
    procedure :: model_cf => prt_cf
    procedure :: model_fc => prt_fc
    procedure :: model_cc => prt_cc
    procedure :: model_cq => prt_cq
    procedure :: model_bd => prt_bd
    procedure :: model_ot => prt_ot
    procedure :: model_da => prt_da
    procedure :: model_solve => prt_solve
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: package_create
    procedure, private :: ftype_check
    procedure :: get_iasym => prt_get_iasym
    procedure, private :: prt_ot_flow
    procedure, private :: prt_ot_saveflow
    procedure, private :: prt_ot_printflow
    procedure, private :: prt_ot_dv
    procedure, private :: prt_ot_bdsummary
    procedure, private :: prt_ot_obs
    procedure, private :: prt_cq_sto
    procedure, private :: create_packages
    procedure, private :: create_bndpkgs
    procedure, private :: create_lstfile
    procedure, private :: log_namfile_options
    procedure :: get_method => prt_get_method

  end type PrtModelType

contains

  !> @brief Create a new particle tracking model object
  subroutine prt_cr(filename, id, modelname)
    ! -- modules
    use ListsModule, only: basemodellist
    use BaseModelModule, only: AddBaseModelToList
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CompilerVersion
    ! use MemoryManagerModule, only: mem_allocate, mem_set_print_option, mem_write_usage
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use GwfNamInputModule, only: GwfNamParamFoundType
    ! -- dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    character(len=*), intent(in) :: modelname
    ! -- local
    type(PrtModelType), pointer :: this
    class(BaseModelType), pointer :: model
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: lst_fname
    type(GwfNamParamFoundType) :: found
    !
    ! -- Allocate a new PRT Model (this)
    allocate (this)
    !
    ! -- Set this before any allocs in the memory manager can be done
    this%memoryPath = create_mem_path(modelname)
    !
    ! -- Allocate track data object
    allocate (this%trackdata)
    !
    ! -- Allocate scalars and add model to basemodellist
    call this%allocate_scalars(modelname)
    model => this
    call AddBaseModelToList(basemodellist, model)
    !
    ! -- Assign values
    this%filename = filename
    this%name = modelname
    this%macronym = 'PRT'
    this%id = id
    !
    ! -- set input model namfile memory path
    input_mempath = create_mem_path(modelname, 'NAM', idm_context)
    !
    ! -- copy option params from input context
    call mem_set_value(this%iprpak, 'PRINT_INPUT', input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', input_mempath, &
                       found%print_flows)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', input_mempath, found%save_flows)
    !
    ! -- create the list file
    call this%create_lstfile(lst_fname, filename, found%list)
    !
    ! -- activate save_flows if found
    if (found%save_flows) then
      this%ipakcb = -1
    end if
    !
    ! -- log set options
    if (this%iout > 0) then
      call this%log_namfile_options(found)
    end if
    !
    ! -- create model packages
    call this%create_packages()
    !
    ! -- return
    return
  end subroutine prt_cr

  !> @brief Define packages
  !
  ! (1) call df routines for each package
  ! (2) set variables and pointers
  !
  !<
  subroutine prt_df(this)
    ! -- modules
    use PrtPrpModule, only: PrtPrpType
    use ModelPackageInputsModule, only: NIUNIT_PRT
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Define packages and utility objects
    call this%dis%dis_df()
    call this%fmi%fmi_df(this%dis)
    call this%oc%oc_df()
    call this%budget%budget_df(NIUNIT_PRT, 'MASS', 'M')
    !
    ! -- Assign or point model members to dis members
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja
    !
    ! -- Define packages and assign iout for time series managers
    ! -- Set maximum number of particles
    ! this%npartmax = this%npart
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%neq, this%dis)
      packobj%TsManager%iout = this%iout
      packobj%TasManager%iout = this%iout
    end do
    !
    ! -- Allocate model arrays
    call this%allocate_arrays()
    !
    ! -- Store information needed for observations
    call this%obs%obs_df(this%iout, this%name, 'PRT', this%dis)
    !
    ! -- return
    return
  end subroutine prt_df

  !> @brief Allocate and read
  !
  ! (1) allocates and reads packages part of this model,
  ! (2) allocates memory for arrays part of this model object
  !
  !<
  subroutine prt_ar(this)
    ! -- modules
    use ConstantsModule, only: DHNOFLO
    use PrtPrpModule, only: PrtPrpType
    use PrtMipModule, only: PrtMipType
    ! -- dummy
    class(PrtModelType) :: this
    ! -- locals
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Allocate and read modules attached to model
    call this%fmi%fmi_ar(this%ibound)
    if (this%inmip > 0) call this%mip%mip_ar()
    !
    ! -- set up output control
    call this%oc%oc_ar(this%x, this%dis, DHNOFLO)
    call this%budget%set_ibudcsv(this%oc%ibudcsv)
    !
    ! -- Package input files now open, so allocate and read
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType)
        call packobj%prp_set_pointers(this%ibound, this%mip%izone)
      end select
      ! -- Read and allocate package
      call packobj%bnd_ar()
    end do
    !
    call this%methodDis%init(this%fmi, &
                             this%flowja, &
                             this%mip%porosity, &
                             this%mip%retfactor, &
                             this%mip%izone, &
                             this%trackdata)
    call this%methodDisv%init(this%fmi, &
                              this%flowja, &
                              this%mip%porosity, &
                              this%mip%retfactor, &
                              this%mip%izone, &
                              this%trackdata)
    !
    ! -- Set trackdata file units
    this%trackdata%ibinun = this%oc%itrkout
    this%trackdata%icsvun = this%oc%itrkcsv
    !
    ! -- return
    return
  end subroutine prt_ar

  !> @brief Read and prepare (calls package read and prepare routines)
  subroutine prt_rp(this)
    ! -- modules
    use TdisModule, only: readnewdata
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return
    !
    ! -- Read and prepare
    if (this%inoc > 0) call this%oc%oc_rp()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_rp()
      call packobj%bnd_rp_obs()
    end do
    !
    ! -- Return
    return
  end subroutine prt_rp

  !> @brief Time step advance (calls package advance subroutines)
  subroutine prt_ad(this)
    ! -- modules
    use SimVariablesModule, only: isimcheck, iFailedStepRetry
    ! -- dummy
    class(PrtModelType) :: this
    class(BndType), pointer :: packobj
    ! -- local
    integer(I4B) :: irestore
    integer(I4B) :: ip, n, i
    !
    ! -- Reset state variable
    irestore = 0
    if (iFailedStepRetry > 0) irestore = 1
    !
    ! -- Copy masssto into massstoold
    do n = 1, this%dis%nodes
      this%massstoold(n) = this%masssto(n)
    end do
    !
    ! -- Advance fmi
    call this%fmi%fmi_ad()
    !
    ! -- Advance
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ad()
      if (isimcheck > 0) then
        call packobj%bnd_ck()
      end if
    end do
    !
    ! -- Push simulated values to preceding time/subtime step
    call this%obs%obs_ad()
    !
    ! -- Initialize the flowja array.  Flowja is calculated each time,
    !    even if output is suppressed.  (Flowja represents flow of particle
    !    mass and is positive into a cell.  Currently, each particle is assigned
    !    unit mass.)  Flowja is updated continually as particles are tracked
    !    over the time step and at the end of the time step.  The diagonal
    !    position of the flowja array will contain the flow residual.
    do i = 1, this%nja
      this%flowja(i) = DZERO
    end do
    !
    ! -- return
    return
  end subroutine prt_ad

  !> @brief Calculate coefficients
  subroutine prt_cf(this, kiter)
    ! -- modules
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: kiter
    ! nothing to do
    return
  end subroutine prt_cf

  !> @brief Fill coefficients
  subroutine prt_fc(this, kiter, amatsln, njasln, inwtflag)
    ! -- modules
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: njasln
    real(DP), dimension(njasln), intent(inout) :: amatsln
    integer(I4B), intent(in) :: inwtflag
    ! nothing to do
    return
  end subroutine prt_fc

  !> @brief Final convergence check (calls package for cc routines)
  subroutine prt_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    ! nothing to do
    return
  end subroutine prt_cc

  !> @brief Calculate intercell flow (flowja)
  subroutine prt_cq(this, icnvg, isuppress_output)
    ! -- modules
    use SparseModule, only: csr_diagsum
    use TdisModule, only: delt
    use PrtPrpModule, only: PrtPrpType
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    real(DP) :: tled
    !
    ! -- Flowja is calculated each time, even if output is suppressed.
    !    Flowja represents flow of particle mass and is positive into a cell.
    !    Currently, each particle is assigned unit mass.
    !
    ! -- Reciprocal of time step size.
    tled = DONE / delt
    !
    ! -- Flowja was updated continually as particles were tracked over the
    !    time step.  At this point, flowja contains the net particle mass
    !    exchanged between cells during the time step.  To convert these to
    !    flow rates (particle mass per time), divide by the time step size.
    do i = 1, this%nja
      this%flowja(i) = this%flowja(i) * tled
    end do
    !
    ! -- Particle mass storage
    call this%prt_cq_sto()
    !
    ! -- Go through packages and call cq routines.  cf() routines are called
    !    first to regenerate non-linear terms to be consistent with the final
    !    conc solution.
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf(reset_mover=.false.) ! kluge note: cf call probably not needed
      call packobj%bnd_cq(this%x, this%flowja)
    end do
    !
    ! -- Finalize calculation of flowja by adding face flows to the diagonal.
    !    This results in the flow residual being stored in the diagonal
    !    position for each cell.
    call csr_diagsum(this%dis%con%ia, this%flowja)
    !
    ! -- Return
    return
  end subroutine prt_cq

  !> @brief Calculate particle mass storage
  subroutine prt_cq_sto(this)
    ! -- modules
    use TdisModule, only: delt
    use PrtPrpModule, only: PrtPrpType
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    integer(I4B) :: n
    integer(I4B) :: np
    integer(I4B) :: idiag
    integer(I4B) :: istatus
    real(DP) :: tled
    real(DP) :: rate
    !
    ! -- Reciprocal of time step size.
    tled = DONE / delt
    !
    ! -- Particle mass storage rate
    do n = 1, this%dis%nodes
      this%masssto(n) = DZERO
      this%ratesto(n) = DZERO
    end do
    do ip = 1, this%bndlist%Count() ! kluge note: could accumulate masssto on the fly in prt_solve instead
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType) ! kluge?
        do np = 1, packobj%npart
          istatus = packobj%partlist%istatus(np)
          if ((istatus > 0) .and. (istatus /= 8)) then ! kluge note: refine these conditions as necessary
            n = packobj%partlist%iTrackingDomain(np, 2)
            ! -- Each particle currently assigned unit mass
            this%masssto(n) = this%masssto(n) + DONE
          end if
        end do
      end select
    end do
    do n = 1, this%dis%nodes ! kluge note: set rate to zero and skip inactive nodes?
      rate = -(this%masssto(n) - this%massstoold(n)) * tled
      this%ratesto(n) = rate
      idiag = this%dis%con%ia(n)
      this%flowja(idiag) = this%flowja(idiag) + rate
    end do
    !
    ! -- Return
    return
  end subroutine prt_cq_sto

  !> @brief Calculate flows and budget
  !
  ! (1) Calculate intercell flows (flowja)
  ! (2) Calculate package contributions to model budget
  !
  !<
  subroutine prt_bd(this, icnvg, isuppress_output)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: rate_accumulator
    ! use ConstantsModule, only: DZERO
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- Save the solution convergence flag
    this%icnvg = icnvg
    !
    ! -- Budget routines (start by resetting).  Sole purpose of this section
    !    is to add in and outs to model budget.  All ins and out for a model
    !    should be added here to this%budget.  In a subsequent exchange call,
    !    exchange flows might also be added.
    call this%budget%reset()
    call rate_accumulator(this%ratesto, rin, rout)
    call this%budget%addentry(rin, rout, delt, budtxt(1), &
                              isuppress_output, '             PRT')
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd(this%budget)
    end do
    !
    ! -- Return
    return
  end subroutine prt_bd

  !> @brief Print and/or save model output
  subroutine prt_ot(this)
    ! -- modules
    use TdisModule, only: kstp, kper, tdis_ot, endofperiod
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    integer(I4B) :: idvsave
    integer(I4B) :: idvprint
    integer(I4B) :: icbcfl
    integer(I4B) :: icbcun
    integer(I4B) :: ibudfl
    integer(I4B) :: ipflag
    ! -- formats
    character(len=*), parameter :: fmtnocnvg = &
      "(1X,/9X,'****FAILED TO COMPLETE SOLUTION IN TIME STEP ', &
      &I0,' OF STRESS PERIOD ',I0,'****')"
    !
    ! -- Set write and print flags
    idvsave = 0
    idvprint = 0
    icbcfl = 0
    ibudfl = 0
    if (this%oc%oc_save('CONCENTRATION')) idvsave = 1
    if (this%oc%oc_print('CONCENTRATION')) idvprint = 1
    if (this%oc%oc_save('BUDGET')) icbcfl = 1
    if (this%oc%oc_print('BUDGET')) ibudfl = 1
    icbcun = this%oc%oc_save_unit('BUDGET')
    !
    ! -- Override ibudfl and idvprint flags for nonconvergence
    !    and end of period
    ibudfl = this%oc%set_print_flag('BUDGET', this%icnvg, endofperiod)
    idvprint = this%oc%set_print_flag('CONCENTRATION', this%icnvg, endofperiod)
    !
    ! -- Calculate and save observations
    call this%prt_ot_obs()
    !
    ! -- Save and print flows
    call this%prt_ot_flow(icbcfl, ibudfl, icbcun)
    !
    ! -- Save and print dependent variables
    call this%prt_ot_dv(idvsave, idvprint, ipflag)
    !
    ! -- Print budget summaries
    call this%prt_ot_bdsummary(ibudfl, ipflag)
    !
    ! -- Timing Output; if any dependendent variables or budgets
    !    are printed, then ipflag is set to 1.
    if (ipflag == 1) call tdis_ot(this%iout)
    !
    ! -- Write non-convergence message
    if (this%icnvg == 0) then
      write (this%iout, fmtnocnvg) kstp, kper
    end if
    !
    ! -- Return
    return
  end subroutine prt_ot

  !> @brief Calculate and save observations
  subroutine prt_ot_obs(this)
    class(PrtModelType) :: this
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! -- Calculate and save observations
    call this%obs%obs_bd()
    call this%obs%obs_ot()

    ! -- Calculate and save package obserations
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd_obs()
      call packobj%bnd_ot_obs()
    end do

  end subroutine prt_ot_obs

  !> @brief Save flows
  subroutine prt_ot_flow(this, icbcfl, ibudfl, icbcun)
    use PrtPrpModule, only: PrtPrpType
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! -- Save PRT flows
    call this%prt_ot_saveflow(this%nja, this%flowja, icbcfl, icbcun)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=0, icbcun=icbcun)
    end do

    ! -- Save advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=icbcfl, ibudfl=0)
    end do

    ! -- Print GWF flows
    ! no need to print flowja
    call this%prt_ot_printflow(ibudfl, this%flowja)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=ibudfl, icbcun=0)
    end do

    ! -- Print advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=0, ibudfl=ibudfl)
    end do

  end subroutine prt_ot_flow

  !> @brief Save intercell flows
  subroutine prt_ot_saveflow(this, nja, flowja, icbcfl, icbcun)
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: nja
    real(DP), dimension(nja), intent(in) :: flowja
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
    ! -- local
    integer(I4B) :: ibinun
    ! -- formats
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- Write the face flows if requested
    if (ibinun /= 0) then
      call this%dis%record_connection_array(flowja, ibinun, this%iout)
    end if
    !
    ! -- Return
    return
  end subroutine prt_ot_saveflow

  !> @brief Print intercell flows
  subroutine prt_ot_printflow(this, ibudfl, flowja)
    ! -- modules
    use TdisModule, only: kper, kstp
    use ConstantsModule, only: LENBIGLINE
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: ibudfl
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    character(len=LENBIGLINE) :: line
    character(len=30) :: tempstr
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm
    ! -- formats
    character(len=*), parameter :: fmtiprflow = &
                "(/,4x,'CALCULATED INTERCELL FLOW &
                &FOR PERIOD ', i0, ' STEP ', i0)"
    !
    ! -- Write flowja to list file if requested
    if (ibudfl /= 0 .and. this%iprflow > 0) then
      write (this%iout, fmtiprflow) kper, kstp
      do n = 1, this%dis%nodes
        line = ''
        call this%dis%noder_to_string(n, tempstr)
        line = trim(tempstr)//':'
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          m = this%dis%con%ja(ipos)
          call this%dis%noder_to_string(m, tempstr)
          line = trim(line)//' '//trim(tempstr)
          qnm = flowja(ipos)
          write (tempstr, '(1pg15.6)') qnm
          line = trim(line)//' '//trim(adjustl(tempstr))
        end do
        write (this%iout, '(a)') trim(line)
      end do
    end if
    !
    ! -- Return
    return
  end subroutine prt_ot_printflow

  !> @brief Print dependent variables
  subroutine prt_ot_dv(this, idvsave, idvprint, ipflag)
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    integer(I4B), intent(inout) :: ipflag
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! -- Print advanced package dependent variables
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_dv(idvsave, idvprint)
    end do

    ! -- save head and print head
    call this%oc%oc_ot(ipflag)

  end subroutine prt_ot_dv

  !> @brief Print budget summary
  subroutine prt_ot_bdsummary(this, ibudfl, ipflag)
    ! -- modules
    use TdisModule, only: kstp, kper, totim
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(inout) :: ipflag
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Package budget summary
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_bdsummary(kstp, kper, this%iout, ibudfl)
    end do
    !
    ! -- model budget summary
    if (ibudfl /= 0) then
      ipflag = 1
      ! -- model budget summary
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if

    ! -- Write to budget csv
    call this%budget%writecsv(totim)

  end subroutine prt_ot_bdsummary

  !> @brief Deallocate
  subroutine prt_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorylist_remove
    use SimVariablesModule, only: idm_context
    use MethodCellPoolModule, only: destroy_methodCellPool ! kluge
    use MethodSubcellPoolModule, only: destroy_methodSubcellPool ! kluge
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Deallocate idm memory
    call memorylist_remove(this%name, 'NAM', idm_context)
    call memorylist_remove(component=this%name, context=idm_context)
    !
    ! -- Internal flow packages deallocate
    call this%dis%dis_da()
    call this%fmi%fmi_da()
    call this%mip%mip_da()
    call this%budget%budget_da()
    call this%oc%oc_da()
    call this%obs%obs_da()
    !
    ! -- Internal package objects
    deallocate (this%dis)
    deallocate (this%fmi)
    deallocate (this%mip)
    deallocate (this%budget)
    deallocate (this%oc)
    deallocate (this%obs)
    !
    call this%methodDis%destroy
    call this%methodDisv%destroy
    ! call this%methodDisu%destroy
    !
    deallocate (this%methodDis)
    deallocate (this%methodDisv)
    ! deallocate(this%methodDisu)
    !
    ! -- Boundary packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_da()
      deallocate (packobj)
    end do
    !
    ! -- Scalars
    call mem_deallocate(this%infmi)
    call mem_deallocate(this%inmip)
    call mem_deallocate(this%inadv)
    call mem_deallocate(this%indsp)
    call mem_deallocate(this%inssm)
    call mem_deallocate(this%inmst)
    call mem_deallocate(this%inmvt)
    call mem_deallocate(this%inoc)
    call mem_deallocate(this%inobs)
    call mem_deallocate(this%nprp)
    call mem_deallocate(this%trackdata%ibinun)
    call mem_deallocate(this%trackdata%icsvun)
    !
    ! -- Arrays
    call mem_deallocate(this%masssto)
    call mem_deallocate(this%massstoold)
    call mem_deallocate(this%ratesto)
    !
    ! -- Track data object
    deallocate (this%trackdata)
    !
    ! -- TrackingModelType
    call this%TrackingModelType%model_da()
    !
    ! -- return
    return
  end subroutine prt_da

  !> @brief Return 1 if any package makes the matrix asymmetric, otherwise 0
  function prt_get_iasym(this) result(iasym)
    ! -- modules
    class(PrtModelType) :: this
    ! -- local
    integer(I4B) :: iasym
    !
    ! -- Start by setting iasym to zero
    iasym = 0
    !
    ! -- return
    return
  end function prt_get_iasym

  !> @brief Allocate memory for non-allocatable members
  subroutine allocate_scalars(this, modelname)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PrtModelType) :: this
    character(len=*), intent(in) :: modelname
    !
    ! -- allocate members from parent class
    call this%TrackingModelType%allocate_scalars(modelname)
    !
    ! -- allocate members that are part of model class
    call mem_allocate(this%infmi, 'INFMI', this%memoryPath)
    call mem_allocate(this%inmip, 'INMIP', this%memoryPath)
    call mem_allocate(this%inmvt, 'INMVT', this%memoryPath)
    call mem_allocate(this%inmst, 'INMST', this%memoryPath)
    call mem_allocate(this%inadv, 'INADV', this%memoryPath)
    call mem_allocate(this%indsp, 'INDSP', this%memoryPath)
    call mem_allocate(this%inssm, 'INSSM', this%memoryPath)
    call mem_allocate(this%inoc, 'INOC ', this%memoryPath)
    call mem_allocate(this%inobs, 'INOBS', this%memoryPath)
    call mem_allocate(this%nprp, 'NPRP', this%memoryPath) ! kluge?
    call mem_allocate(this%trackdata%ibinun, 'ITRKBIN', this%memoryPath)
    call mem_allocate(this%trackdata%icsvun, 'ITRKCSV', this%memoryPath)
    !
    this%infmi = 0
    this%inmip = 0
    this%inmvt = 0
    this%inmst = 0
    this%inadv = 0
    this%indsp = 0
    this%inssm = 0
    this%inoc = 0
    this%inobs = 0
    this%nprp = 0
    !
    ! -- return
    return
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  subroutine allocate_arrays(this)
    ! use ConstantsModule, only: DZERO
    use MemoryManagerModule, only: mem_allocate
    class(PrtModelType) :: this
    integer(I4B) :: n
    !
    ! -- Allocate arrays in TrackingModelType
    call this%TrackingModelType%allocate_arrays()
    !
    ! -- allocate and initialize arrays for mass storage
    call mem_allocate(this%masssto, this%dis%nodes, &
                      'MASSSTO', this%memoryPath)
    call mem_allocate(this%massstoold, this%dis%nodes, &
                      'MASSSTOOLD', this%memoryPath)
    call mem_allocate(this%ratesto, this%dis%nodes, &
                      'RATESTO', this%memoryPath)
    do n = 1, this%dis%nodes
      this%masssto(n) = DZERO
      this%massstoold(n) = DZERO
      this%ratesto(n) = DZERO
    end do
    !
    ! -- return
    return
  end subroutine allocate_arrays

  !> @brief Create boundary condition packages for this model
  subroutine package_create(this, filtyp, ipakid, ipaknum, pakname, inunit, &
                            iout)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use PrtPrpModule, only: prp_create
    use ApiModule, only: api_create
    ! -- dummy
    class(PrtModelType) :: this
    character(len=*), intent(in) :: filtyp
    character(len=LINELENGTH) :: errmsg
    integer(I4B), intent(in) :: ipakid
    integer(I4B), intent(in) :: ipaknum
    character(len=*), intent(in) :: pakname
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- local
    class(BndType), pointer :: packobj
    class(BndType), pointer :: packobj2
    integer(I4B) :: ip
    !
    ! -- This part creates the package object
    select case (filtyp)
    case ('PRP6')
      this%nprp = this%nprp + 1 ! kluge?
      call prp_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%fmi)
    case ('API6')
      call api_create(packobj, ipakid, ipaknum, inunit, iout, this%name, pakname)
    case default
      write (errmsg, *) 'Invalid package type: ', filtyp
      call store_error(errmsg, terminate=.TRUE.)
    end select
    !
    ! -- Packages is the bndlist that is associated with the parent model
    ! -- The following statement puts a pointer to this package in the ipakid
    ! -- position of packages.
    do ip = 1, this%bndlist%Count()
      packobj2 => GetBndFromList(this%bndlist, ip)
      if (packobj2%packName == pakname) then
        write (errmsg, '(a,a)') 'Cannot create package.  Package name  '// &
          'already exists: ', trim(pakname)
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end do
    call AddBndToList(this%bndlist, packobj)
    !
    ! -- return
    return
  end subroutine package_create

  !> @brief Check to make sure required input files have been specified
  subroutine ftype_check(this, indis)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: indis
    ! -- local
    character(len=LINELENGTH) :: errmsg
    !
    ! -- Check for DIS(u) and MIP. Stop if not present.
    if (indis == 0) then
      write (errmsg, '(1x,a)') &
        'Discretization (DIS6, DISV6, or DISU6) package not specified.'
      call store_error(errmsg)
    end if
    if (this%inmip == 0) then
      write (errmsg, '(1x,a)') &
        'Model input (MIP6) package not specified.'
      call store_error(errmsg)
    end if
    !
    if (count_errors() > 0) then
      write (errmsg, '(1x,a)') 'One or more required package(s) not specified.'
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end if
    !
    ! -- return
    return
  end subroutine ftype_check

  !> @brief Cast to PrtModelType
  function CastAsPrtModel(model) result(prtmodel)
    class(*), pointer :: model
    class(PrtModelType), pointer :: prtmodel
    !
    prtmodel => null()
    if (.not. associated(model)) return
    select type (model)
    class is (PrtModelType)
      prtmodel => model
    end select
  end function CastAsPrtModel

  !> @brief Solve the model
  subroutine prt_solve(this)
    ! -- modules
    use TdisModule, only: kper, kstp, totimc, totim, nper, nstp
    use PrtPrpModule, only: PrtPrpType
    ! -- dummy variables
    class(PrtModelType) :: this
    ! -- local variables
    integer(I4B) :: np, ip
    class(BndType), pointer :: packobj
    type(ParticleType), pointer :: particle
    class(MethodType), pointer :: method
    real(DP) :: tmax
    logical(LGP) :: limited
    integer(I4B) :: iprp

    call create_particle(particle)

    ! -- Loop over PRP packages
    iprp = 0
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType)

        iprp = iprp + 1

        ! -- Loop over particles in package
        do np = 1, packobj%npart

          ! -- If particle inactive, skip
          if (packobj%partlist%istatus(np) .ne. 1) cycle
          !
          ! -- Reset the particle's coordinate transformation
          call particle%reset_transf()

          ! -- Update particle properties from particle list
          call particle%load_from_list(packobj%partlist, this%id, iprp, np)

          ! -- Unless in last stress period and it has only one time step,
          ! -- limit max time to no later than end of time step
          tmax = particle%tstop
          limited = .true.
          if (kper == nper) then
            if (nstp(kper) == 1) then
              limited = .false.
            end if
          end if
          if (limited) then
            if (totim < particle%tstop) then
              tmax = totim
            end if
          end if

          ! -- Get the tracking method
          method => this%get_method(particle)

          ! -- If particle released during this time step, record its
          ! -- initial location in track data
          if (particle%trelease .ge. totimc) &
            call this%trackdata%save_record(particle, kper=kper, &
                                            kstp=kstp, reason=0)

          ! -- Apply the tracking method
          call method%apply(particle, tmax)

          ! -- Update particle in PRP package
          call packobj%partlist%update_from_particle(particle, np)
          !
        end do
      end select
    end do
    !
    call particle%destroy() ! kluge???
    deallocate (particle)
    !
    ! -- return
    return
  end subroutine prt_solve

  !> @brief Get the grid discretization method
  function prt_get_method(this, particle) result(method)
    !
    use GwfDisModule, only: GwfDisType
    use GwfDisvModule, only: GwfDisvType
    use GwfDisuModule, only: GwfDisuType
    !
    class(PrtModelType) :: this
    type(ParticleType), pointer :: particle ! kluge note: needed???
    class(MethodType), pointer :: method
    ! -- local
    !
    select type (dis => this%dis)
    type is (GwfDisType)
      method => this%methodDis
    type is (GwfDisvType)
      method => this%methodDisv
    type is (GwfDisuType)
      print *, "DISU grids not currently supported" ! kluge
      stop
    class default
      print *, "Grid type not recognized" ! kluge
      stop
    end select
    !
    return
    !
  end function prt_get_method

  !> @brief Source package info and begin to process
  !<
  subroutine create_bndpkgs(this, bndpkgs, pkgtypes, pkgnames, &
                            mempaths, inunits)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy
    class(PrtModelType) :: this
    integer(I4B), dimension(:), allocatable, intent(inout) :: bndpkgs
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgnames
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: mempaths
    integer(I4B), dimension(:), contiguous, &
      pointer, intent(inout) :: inunits
    ! -- local
    integer(I4B) :: ipakid, ipaknum
    character(len=LENFTYPE) :: pkgtype, bndptype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B) :: n

    if (allocated(bndpkgs)) then
      !
      ! -- create stress packages
      ipakid = 1
      bndptype = ''
      do n = 1, size(bndpkgs)
        !
        pkgtype = pkgtypes(bndpkgs(n))
        pkgname = pkgnames(bndpkgs(n))
        mempath = mempaths(bndpkgs(n))
        inunit => inunits(bndpkgs(n))
        !
        if (bndptype /= pkgtype) then
          ipaknum = 1
          bndptype = pkgtype
        end if
        !
        call this%package_create(pkgtype, ipakid, ipaknum, pkgname, inunit, &
                                 this%iout)
        ipakid = ipakid + 1
        ipaknum = ipaknum + 1
      end do
      !
      ! -- cleanup
      deallocate (bndpkgs)
    end if
    !
    ! -- return
    return
  end subroutine create_bndpkgs

  !> @brief Source package info and begin to process
  !<
  subroutine create_packages(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: expandarray
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    use GwfDisModule, only: dis_cr
    use GwfDisvModule, only: disv_cr
    use GwfDisuModule, only: disu_cr
    use BudgetModule, only: budget_cr
    use MethodDisModule, only: create_methodDis
    use MethodDisvModule, only: create_methodDisv
    ! use MethodDisuModule, only: create_methodDisu
    use MethodCellPoolModule, only: create_methodCellPool
    use MethodSubcellPoolModule, only: create_methodSubcellPool
    use PrtMipModule, only: mip_cr
    use PrtFmiModule, only: fmi_cr
    use PrtOcModule, only: oc_cr
    use PrtObsModule, only: prt_obs_cr
    ! -- dummy
    class(PrtModelType) :: this
    ! -- local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    integer(I4B), dimension(:), contiguous, &
      pointer :: inunits => null()
    character(len=LENMEMPATH) :: model_mempath
    character(len=LENFTYPE) :: pkgtype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B), dimension(:), allocatable :: bndpkgs
    integer(I4B) :: n
    integer(I4B) :: indis = 0 ! DIS enabled flag
    character(len=LENMEMPATH) :: mempathmip = ''
    !
    ! -- set input memory paths, input/model and input/model/namfile
    model_mempath = create_mem_path(component=this%name, context=idm_context)
    !
    ! -- set pointers to model path package info
    call mem_setptr(pkgtypes, 'PKGTYPES', model_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', model_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', model_mempath)
    call mem_setptr(inunits, 'INUNITS', model_mempath)
    !
    do n = 1, size(pkgtypes)
      !
      ! attributes for this input package
      pkgtype = pkgtypes(n)
      pkgname = pkgnames(n)
      mempath = mempaths(n)
      inunit => inunits(n)
      !
      ! -- create dis package first as it is a prerequisite for other packages
      select case (pkgtype)
      case ('DIS6')
        indis = 1
        call dis_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DISV6')
        indis = 1
        call disv_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DISU6')
        indis = 1
        call disu_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('MIP6')
        this%inmip = 1
        mempathmip = mempath
      case ('FMI6')
        this%infmi = inunit
      case ('OC6')
        this%inoc = inunit
      case ('OBS6')
        this%inobs = inunit
      case ('PRP6')
        call expandarray(bndpkgs)
        bndpkgs(size(bndpkgs)) = n
      case default
        ! TODO
      end select
    end do
    !
    ! -- Create utility objects
    call budget_cr(this%budget, this%name)
    call create_methodDis(this%methodDis)
    call create_methodDisv(this%methodDisv)
    ! call create_methodDisu(this%methodDisu)
    call create_methodCellPool() ! kluge
    call create_methodSubcellPool() ! kluge
    !
    ! -- Create packages that are tied directly to model
    call mip_cr(this%mip, this%name, mempathmip, this%inmip, this%iout, this%dis)
    call fmi_cr(this%fmi, this%name, this%infmi, this%iout)
    call oc_cr(this%oc, this%name, this%inoc, this%iout)
    call prt_obs_cr(this%obs, this%inobs)
    !
    ! -- Check to make sure that required ftype's have been specified
    call this%ftype_check(indis)
    !
    call this%create_bndpkgs(bndpkgs, pkgtypes, pkgnames, mempaths, inunits)
    !

  end subroutine create_packages

  subroutine create_lstfile(this, lst_fname, model_fname, defined)
    ! -- modules
    use KindModule, only: LGP
    use InputOutputModule, only: openfile, getunit
    ! -- dummy
    class(PrtModelType) :: this
    character(len=*), intent(inout) :: lst_fname
    character(len=*), intent(in) :: model_fname
    logical(LGP), intent(in) :: defined
    ! -- local
    integer(I4B) :: i, istart, istop
    !
    ! -- set list file name if not provided
    if (.not. defined) then
      !
      ! -- initialize
      lst_fname = ' '
      istart = 0
      istop = len_trim(model_fname)
      !
      ! -- identify '.' character position from back of string
      do i = istop, 1, -1
        if (model_fname(i:i) == '.') then
          istart = i
          exit
        end if
      end do
      !
      ! -- if not found start from string end
      if (istart == 0) istart = istop + 1
      !
      ! -- set list file name
      lst_fname = model_fname(1:istart)
      istop = istart + 3
      lst_fname(istart:istop) = '.lst'
    end if
    !
    ! -- create the list file
    this%iout = getunit()
    call openfile(this%iout, 0, lst_fname, 'LIST', filstat_opt='REPLACE')
    !
    ! -- write list file header
    call write_listfile_header(this%iout, 'PARTICLE TRACKING MODEL (PRT)')
    !
    ! -- return
    return
  end subroutine create_lstfile

  !> @brief Write model namfile options to list file
  !<
  subroutine log_namfile_options(this, found)
    use GwfNamInputModule, only: GwfNamParamFoundType
    class(PrtModelType) :: this
    type(GwfNamParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'NAMEFILE OPTIONS:'

    if (found%newton) then
      write (this%iout, '(4x,a)') &
        'NEWTON-RAPHSON method enabled for the model.'
      if (found%under_relaxation) then
        write (this%iout, '(4x,a,a)') &
          'NEWTON-RAPHSON UNDER-RELAXATION based on the bottom ', &
          'elevation of the model will be applied to the model.'
      end if
    end if

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'STRESS PACKAGE INPUT WILL BE PRINTED '// &
        'FOR ALL MODEL STRESS PACKAGES'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'PACKAGE FLOWS WILL BE PRINTED '// &
        'FOR ALL MODEL PACKAGES'
    end if

    if (found%save_flows) then
      write (this%iout, '(4x,a)') &
        'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END NAMEFILE OPTIONS:'
  end subroutine log_namfile_options

end module PrtModule
