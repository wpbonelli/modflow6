!> @brief Particle event buffering strategies.

module ParticleTrackEventBufferModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ParticleEventModule, only: RELEASE, FEATEXIT, SUBFEXIT, TIMESTEP, &
                                 TERMINATE, WEAKSINK, USERTIME, DROPPED

  implicit none

  !> @brief Selection of particle events to write to a track file.
  type :: ParticleTrackEventSelectionType
    logical(LGP) :: release = .false. !< track release events
    logical(LGP) :: featexit = .false. !< track grid feature exits
    logical(LGP) :: subfexit = .false. !< track sub-grid-scale feature exits
    logical(LGP) :: timestep = .false. !< track timestep ends
    logical(LGP) :: terminate = .false. !< track termination events
    logical(LGP) :: weaksink = .false. !< track weak sink exit events
    logical(LGP) :: usertime = .false. !< track user-selected times
    logical(LGP) :: dropped = .false. !< track water table drops
  end type ParticleTrackEventSelectionType

  !> @brief Flat record of a particle track event.
  type :: ParticleTrackRecordType
    integer(I4B) :: kper, kstp, imdl, iprp, irpt, ilay, icu, izone
    integer(I4B) :: istatus, ireason
    real(DP) :: trelease, ttrack, x, y, z
    character(len=40) :: name
  end type ParticleTrackRecordType

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! Can be associated with a particle release point (PRP) package
  !! or with an entire model, and can be binary or text (CSV). An
  !! event selection controls which event types are written. Scope
  !! is controlled by iprp: -1 writes all particles; any other value
  !! writes only particles originating from that PRP package index.
  !<
  type :: ParticleTrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
    type(ParticleTrackEventSelectionType), public :: selected !< event selection
  contains
    procedure, public :: wants_record
  end type ParticleTrackFileType

  !> @brief Event buffering strategy
  type, abstract :: ParticleTrackEventBufferType
    integer(I4B) :: nrecords = 0 !< number of records stored
  contains
    procedure :: init => buffer_init !< open/allocate resources
    procedure(buffer_append), deferred :: append !< buffer one record
    procedure(buffer_flush), deferred :: flush !< write buffered records, reset
    procedure(buffer_simple), deferred :: discard !< reset without writing
    procedure(buffer_simple), deferred :: destroy !< release resources
  end type ParticleTrackEventBufferType

  abstract interface
    subroutine buffer_append(this, rec)
      import ParticleTrackEventBufferType, ParticleTrackRecordType
      class(ParticleTrackEventBufferType) :: this
      type(ParticleTrackRecordType), intent(in) :: rec
    end subroutine buffer_append

    subroutine buffer_flush(this, files)
      import ParticleTrackEventBufferType, ParticleTrackFileType
      class(ParticleTrackEventBufferType) :: this
      type(ParticleTrackFileType), intent(in) :: files(:)
    end subroutine buffer_flush

    ! Shared interface for discard and destroy: both take only `this`.
    subroutine buffer_simple(this)
      import ParticleTrackEventBufferType
      class(ParticleTrackEventBufferType) :: this
    end subroutine buffer_simple
  end interface

contains

  subroutine buffer_init(this)
    class(ParticleTrackEventBufferType) :: this
  end subroutine buffer_init

  !> @brief Whether this file wants a given record.
  !!
  !! A record is wanted if both scope and event selection match:
  !! scope matches when iprp is -1 (model-level) or equals the
  !! record's iprp; event selection matches when the corresponding
  !! flag in the file's selected struct is .true.
  logical function wants_record(this, rec)
    class(ParticleTrackFileType), intent(in) :: this
    type(ParticleTrackRecordType), intent(in) :: rec
    logical :: event_ok

    if (this%iun <= 0) then
      wants_record = .false.
      return
    end if

    if (this%iprp /= -1 .and. this%iprp /= rec%iprp) then
      wants_record = .false.
      return
    end if

    select case (rec%ireason)
    case (RELEASE)
      event_ok = this%selected%release
    case (FEATEXIT)
      event_ok = this%selected%featexit
    case (SUBFEXIT)
      event_ok = this%selected%subfexit
    case (TIMESTEP)
      event_ok = this%selected%timestep
    case (TERMINATE)
      event_ok = this%selected%terminate
    case (WEAKSINK)
      event_ok = this%selected%weaksink
    case (USERTIME)
      event_ok = this%selected%usertime
    case (DROPPED)
      event_ok = this%selected%dropped
    case default
      event_ok = .false.
    end select

    wants_record = event_ok
  end function wants_record

  !> @brief Save an event record to a binary or CSV file.
  subroutine save_record(iun, rec, csv)
    integer(I4B), intent(in) :: iun
    type(ParticleTrackRecordType), intent(in) :: rec
    logical(LGP), intent(in) :: csv

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, trim(rec%name)
    else
      write (iun) &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, rec%name
    end if
  end subroutine save_record

end module ParticleTrackEventBufferModule
