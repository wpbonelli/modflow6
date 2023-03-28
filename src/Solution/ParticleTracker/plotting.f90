module plottingmodule

  double precision, public :: xscale, yscale, zscale
  integer, public :: ntrk
  double precision, public :: xtrk(0:1000000), ytrk(0:1000000), ztrk(0:1000000)
  logical :: plot1stgridonly, isgridplotted(10, 10), doneplottinggrids

contains

  subroutine x3dstart()
!
    implicit double precision(a - h, o - z)
    !
    ! -- Write X3D header info
    write (160, '(A/2A/4A)') '<?xml version="1.0" encoding="UTF-8"?>', &
      '<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.0', &
      '//EN" "http://www.web3d.org/specifications/x3d-3.0.dtd">', &
      "<X3D profile='Interchange' version='3.0' ", &
      "xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' ", &
      "xsd:noNamespaceSchemaLocation=", &
      "'http://www.web3d.org/specifications/x3d-3.0.xsd'>"
    write (160, '(2X,A/4X,A/2X,A)') "<head>", &
      "<meta content='example.x3d' name='title'/>", "</head>"
    write (160, '(2X,A)') "<Scene>"
    !
    return
    !
  end subroutine

! ------------------------------------------------------------------------------

  subroutine vptbck(xmin, xmax, ymin, ymax, zmin, zmax)
    !
    implicit double precision(a - h, o - z)
    !
    ! -- Write X3D viewpoint and background
    xmid = 5d-1 * (xmin + xmax)
    ymid = 5d-1 * (ymin + ymax)
    zmid = 5d-1 * (zmin + zmax)
    xview = xmid
    yview = ymid
    xspan = xmax - xmin
    yspan = ymax - ymin
    zview = zmax + 2d0 * max(xspan, yspan)
    ovecx = 0d0
    ovecy = 0d0
    ovecz = 1d0
    oangl = 0d0
    write (160, '(4X,A/6X,A,3(2X,G),A/6X,A,4(2X,G),A/6X,A,3(2X,G),A)') &
      "<Viewpoint description='centered_above'", &
      "position='", xview, yview, zview, "'", &
      "orientation='", ovecx, ovecy, ovecz, oangl, "'", &
      "centerOfRotation='", xmid, ymid, zmid, "'/>"
    write (160, '(4X,A)') "<Background skyColor='1. 1. 1.'/>"
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine drwgrid(dis)
!
    use CellDefnModule
    use UtilMiscModule
    use BaseDisModule, only: DisBaseType
    use GwfDisvModule
    use GwfDisModule
    use InputOutputModule ! kluge
    use ConstantsModule, only: DHALF
    implicit double precision(a - h, o - z) ! kluge
    class(DisBaseType), pointer :: dis
    type(CellDefnType), pointer :: cellDefn
    integer :: ic, npolyverts, icu, icu2d
    integer :: ncpl ! kluge???
    integer :: irow, jcol, klay
    double precision :: cellx, celly, dxhalf, dyhalf
    !
    ! -- Draw all cells in grid.
    !
    red = 0.
    grn = 0.
    blu = 0.
    trn = 0.9
    !
    call create_cellDefn(cellDefn)
    !
    ncells = dis%nodes
    ncpl = dis%get_ncpl()
    !
    do ic = 1, ncells
      !
      icu = dis%get_nodeuser(ic)
      !
      select type (dis)
      type is (GwfDisvType)
        !
        icu2d = icu - ((icu - 1) / ncpl) * ncpl ! kluge note: use MOD or MODULO???
        npolyverts = dis%iavert(icu2d + 1) - dis%iavert(icu2d) - 1
        if (npolyverts .le. 0) npolyverts = npolyverts + size(dis%javert) ! kluge???
        cellDefn%npolyverts = npolyverts
        ! -- Load polygon vertices. Note that the polyverts array
        ! -- does not get reallocated if it is already allocated
        ! -- to a size greater than or equal to npolyverts+1.
        call allocate_as_needed(cellDefn%polyvert, npolyverts + 1)
        iavert = dis%iavert(icu2d)
        iavertm1 = iavert - 1
        do m = 1, npolyverts
          j = dis%javert(iavertm1 + m)
          cellDefn%polyvert(m)%ivert = j ! kluge note: is ivert ever used???
          cellDefn%polyvert(m)%x = dis%vertices(1, j) + dis%xorigin
          cellDefn%polyvert(m)%y = dis%vertices(2, j) + dis%yorigin ! kluge note: dis%angrot is not taken into account here
        end do
        ! -- List wraps around for convenience
        cellDefn%polyvert(npolyverts + 1) = cellDefn%polyvert(1)
        !
      type is (GwfDisType)
        !
        cellDefn%npolyverts = 4
        npolyverts = cellDefn%npolyverts
        ! -- Load polygon vertices. Note that the polyverts array
        ! -- does not get reallocated if it is already allocated
        ! -- to a size greater than or equal to npolyverts+1.
        call allocate_as_needed(cellDefn%polyvert, npolyverts + 1)
        call get_ijk(icu, dis%nrow, dis%ncol, dis%nlay, irow, jcol, klay)
        cellx = dis%cellx(jcol) + dis%xorigin
        celly = dis%celly(irow) + dis%yorigin ! kluge note: dis%angrot is not taken into account here
        dxhalf = DHALF * dis%delr(jcol)
        dyhalf = DHALF * dis%delc(irow)
        ! -- SW vertex
        cellDefn%polyvert(1)%ivert = 0 ! kluge note: there's no vertex list to reference for DIS
        cellDefn%polyvert(1)%x = cellx - dxhalf
        cellDefn%polyvert(1)%y = celly - dyhalf
        ! -- NW vertex
        cellDefn%polyvert(2)%ivert = 0
        cellDefn%polyvert(2)%x = cellx - dxhalf
        cellDefn%polyvert(2)%y = celly + dyhalf
        ! -- NE vertex
        cellDefn%polyvert(3)%ivert = 0
        cellDefn%polyvert(3)%x = cellx + dxhalf
        cellDefn%polyvert(3)%y = celly + dyhalf
        ! -- SE vertex
        cellDefn%polyvert(4)%ivert = 0
        cellDefn%polyvert(4)%x = cellx + dxhalf
        cellDefn%polyvert(4)%y = celly - dyhalf
        ! -- List wraps around for convenience
        cellDefn%polyvert(5)%ivert = cellDefn%polyvert(1)%ivert
        cellDefn%polyvert(5)%x = cellDefn%polyvert(1)%x
        cellDefn%polyvert(5)%y = cellDefn%polyvert(1)%y
        !
      end select
      !
      if (icu .le. ncpl) then
        cellDefn%iatop = -1
      else
        cellDefn%iatop = 1
      end if
      cellDefn%top = dis%top(ic)
      cellDefn%bot = dis%bot(ic)
      !
      call drwcell(cellDefn, red, grn, blu, trn)
      !
    end do
    !
    call cellDefn%destroy()
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine drwcell(cellDefn, red, grn, blu, trn)
!
    use CellDefnModule
    implicit double precision(a - h, o - z)
    type(CellDefnType), pointer :: cellDefn
    type(VertexType), pointer :: verts(:)
    !
    ! -- Draw top and/or bottom polygon of cell.
    ! -- (If drawing top, do that first.)
    !
    npolyverts = cellDefn%npolyverts
    verts => cellDefn%polyvert
    iatop = cellDefn%iatop
    if (iatop .lt. 0) then ! kluge layered
      top = cellDefn%top
      call drwpoly(verts, npolyverts, top, red, grn, blu, trn)
    end if
    bot = cellDefn%bot
    call drwpoly(verts, npolyverts, bot, red, grn, blu, trn)
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine drwpoly(verts, numver, z, red, grn, blu, trn)
!
    use CellDefnModule, only: VertexType
    implicit double precision(a - h, o - z)
    type(VertexType) :: verts(numver)
    !
    ! -- Write polygon as indexed line set
    call shapestart(red, grn, blu, trn)
    write (160, '(8X,A)') "<IndexedLineSet coordIndex='"
    write (160, '(10x,999I)') (m, m=0, numver), -1
    ! -- Point coordinates for indexed line set
    write (160, '(10X,A/10X,A)') "'>", "<Coordinate point='"
    do m = 1, numver
      write (160, '(10x,3G)') verts(m)%x * xscale, verts(m)%y * yscale, z * zscale
    end do
    m = 1
    write (160, '(10x,3G)') verts(m)%x * xscale, verts(m)%y * yscale, z * zscale
    write (160, '(10X,A)') "'/>"
    write (160, '(8X,A)') "</IndexedLineSet>"
    call shapeend()
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine drwtrk(xorigin, yorigin, red, grn, blu)
    !
    implicit double precision(a - h, o - z)
    !
    ! -- Draw a particle track
    !
    ! -- Transparency
    trn = 0.
    !
    ! -- Write track as indexed line set
!!  write(160,'(3(A,I,", "),A)') "<!-- PARTICLE ",ip," CELL ",icell," GRID ",igrid," -->"
    write (160, '(6X,A)') "<Shape>"
    ! -- Appearance
    write (160, '(8X,A/10X,A,A,3(F6.3,1X),A,A,F6.3,A,                  &
  &    /10X,A,A,A/8X,A)') &
      "<Appearance>", "<Material diffuseColor='0 0 0' ", &
      "emissiveColor='", red, grn, blu, "' ", &
      "transparency='", trn, "'/>", &
      "<LineProperties applied='true' linetype='1' ", &
      "linewidthScaleFactor='2' metadata='X3DMetadataObject' >", &
      "</LineProperties>", "</Appearance>"
    ! -- Indexed line set
    write (160, '(8X,A)') "<IndexedLineSet coordIndex='"
    write (160, '(10X,999I)') (iv, iv=0, ntrk), -1
    ! -- Point coordinates for indexed line set
    write (160, '(10X,A/10X,A)') "'>", "<Coordinate point='"
    do iv = 0, ntrk
      x = (xorigin + xtrk(iv)) * xscale
      y = (yorigin + ytrk(iv)) * yscale ! kluge note: angrot is not taken into account here
      z = ztrk(iv) * zscale
      write (160, '(10X,3G)') x, y, z
    end do
    write (160, '(10X,A)') "'/>"
    write (160, '(8X,A/6X,A)') "</IndexedLineSet>", "</Shape>"
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine shapestart(red, grn, blu, trn)
!
    implicit double precision(a - h, o - z)
    !
    ! -- Start shape
    write (160, '(6X,A)') "<Shape>"
    ! -- Appearance
    write (160, '(8X,A/10X,A,A,3(F6.3,1X),A,A,F6.3,A/8X,A)') &
      "<Appearance>", "<Material diffuseColor='0 0 0' ", &
      "emissiveColor='", red, grn, blu, "' ", &
      "transparency='", trn, "'/>", "</Appearance>"
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine shapeend()
!
    implicit double precision(a - h, o - z)
    !
    ! -- End shape
    write (160, '(6X,A)') "</Shape>"
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine x3dend()
    !
    implicit double precision(a - h, o - z)
    !
    ! -- Close remaining X3D tags
    write (160, '(2X,A/A)') "</Scene>", "</X3D>"
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine plot_open()
    !
    implicit double precision(a - h, o - z)
    logical :: od
    !
    ! -- Open and begin x3d plot file
    inquire (file="./plot.x3d", OPENED=od)
    if (.not. od) then
      open (unit=160, file="./plot.in", status="old")
      read (160, *) xmin, xmax, xscale
      read (160, *) ymin, ymax, yscale
      read (160, *) zmin, zmax, zscale
      read (160, *) plot1stgridonly
      close (160)
      open (unit=160, file="./plot.x3d", status="unknown")
      call x3dstart()
!!    call vptbck(0.d0,700.d0,0.d0,700.d0,-100.d0,0.d0) ! test006
!!    call vptbck(0.d0,2.d0,0.d0,1.d0,0.d0,1.d0) ! test201_gwtbuy-henryCHD
      call vptbck(xmin, xmax, ymin, ymax, zmin, zmax)
      isgridplotted = .false.
      doneplottinggrids = .false.
    end if
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine plot_close()
    !
    implicit double precision(a - h, o - z)
    logical :: od
    !
    ! -- End and close x3d plot file
    inquire (file="plot.x3d", OPENED=od)
    if (od) then
      call x3dend()
      close (160)
    end if
    !
    return
    !
  end subroutine

!-------------------------------------------------------------------------------

  subroutine plot_track(xorigin, yorigin, kper)
    !
    implicit double precision(a - h, o - z)
    !
    if (2 * (kper / 2) .eq. kper) then
      red = 0d0
      grn = 4d-1
      blu = 2d-1
    else
      red = 0d0
      grn = 3d-1
      blu = 8d-1
    end if
    call drwtrk(xorigin, yorigin, red, grn, blu)
    !
    return
    !
  end subroutine

  subroutine plot_init_track(x, y, z)
    !
    implicit double precision(a - h, o - z)
    !
    ntrk = 0
    xtrk(ntrk) = x
    ytrk(ntrk) = y
    ztrk(ntrk) = z
    !
    return
    !
  end subroutine

  subroutine plot_update_track(levelNext, x, y, z)
    !
    implicit double precision(a - h, o - z)
    !
    if (levelNext .eq. 2) then
      ntrk = ntrk + 1
      if (ntrk .gt. size(xtrk) - 1) then
        print *, "ntrk exceeds hardwired size of ", size(xtrk) - 1 ! kluge
!!      pause
        stop
      end if
      xtrk(ntrk) = x
      ytrk(ntrk) = y
      ztrk(ntrk) = z
    end if
    !
    return
    !
  end subroutine

end module plottingmodule
