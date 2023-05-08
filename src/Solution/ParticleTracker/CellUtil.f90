module CellUtilModule

  implicit none

  private
  public :: CellPolyToCellRect
  public :: CellPolyToCellRectQuad
  public :: MetricForPointInCellPolygon

contains

  !> @brief Convert CellPoly representation to CellRect if possible
  subroutine CellPolyToCellRect(cellPoly, cellRect, istatus)
    use ConstantsModule, only: DONE
    use CellRectModule, only: CellRectType, create_cellRect
    use CellPolyModule, only: CellPolyType
    use CellDefnModule, only: CellDefnType ! kluge???
    ! -- dummy
    type(CellPolyType), pointer :: cellPoly
    type(CellRectType), pointer :: cellRect
    integer :: istatus
    ! -- local
    type(CellDefnType), pointer :: cellDefn
    integer :: ipv, ipv1, ipv2, ipv3, ipv4
    integer, dimension(4) :: ipvnxt = (/2, 3, 4, 1/)
    double precision :: x1, y1, x2, y2, x4, y4
    double precision :: dx2, dy2, dx4, dy4, areax, areay, areaz
    double precision :: xOrigin, yOrigin, zOrigin, dx, dy, dz, sinrot, cosrot
    double precision :: factor, term
    !
    call create_cellRect(cellRect)
    cellDefn => cellPoly%cellDefn
    ! -- kluge note: no check whether conversion is possible; assumes it is
    !
    ! -- Translate and rotate the rectangular cell into local coordinates
    ! -- with x varying from 0 to dx and y varying from 0 to dy. Choose the
    ! -- "south-west" vertex to be the local origin so that the rotation
    ! -- angle is zero if the cell already aligns with the model x and y
    ! -- coordinates. The "southwest" vertex is found by finding the vertex
    ! -- formed either by (1) an edge over which x decreases (going
    ! -- clockwise) followed by an edge over which x does not increase, or
    ! -- by (2) an edge over which y does not decrease (again going
    ! -- clockwise) followed by an edge over which y increases. In the end,
    ! -- ipv1 is the local vertex number (within the cell, taking a value
    ! -- of 1, 2, 3, or 4) of the southwest vertex, and ipv2, ipv3, and
    ! -- ipv4 are the local vertex numbers of the remaining three vertices
    ! -- going clockwise.
    do ipv = 1, 4
      ipv4 = ipv
      ipv1 = ipvnxt(ipv4)
      x4 = cellDefn%polyvert(ipv4)%x
      y4 = cellDefn%polyvert(ipv4)%y
      x1 = cellDefn%polyvert(ipv1)%x
      y1 = cellDefn%polyvert(ipv1)%y
      if (x1 .lt. x4) then
        ipv2 = ipvnxt(ipv1)
        x2 = cellDefn%polyvert(ipv2)%x
        if (x2 .le. x1) then
          y2 = cellDefn%polyvert(ipv2)%y
          exit
        end if
      else if (y1 .ge. y4) then
        ipv2 = ipvnxt(ipv1)
        y2 = cellDefn%polyvert(ipv2)%y
        if (y2 .gt. y1) then
          x2 = cellDefn%polyvert(ipv2)%x
          exit
        end if
      end if
    end do
    ipv3 = ipvnxt(ipv2)
    !
    ! -- Compute upper bounds on the local coordinates (the rectangular
    ! -- dimensions of the cell) and the sine and cosine of the rotation
    ! -- angle, and store local origin information
    xOrigin = x1
    yOrigin = y1
    zOrigin = cellDefn%bot
    dx2 = x2 - xOrigin
    dy2 = y2 - yOrigin
    dx4 = x4 - xOrigin
    dy4 = y4 - yOrigin
    dx = dsqrt(dx4 * dx4 + dy4 * dy4)
    dy = dsqrt(dx2 * dx2 + dy2 * dy2)
    dz = cellDefn%top - zOrigin ! kluge note: need to account for partial saturation
    sinrot = dy4 / dx
    cosrot = dx4 / dx
    cellRect%cellDefn = cellPoly%cellDefn ! kluge???
    cellRect%dx = dx
    cellRect%dy = dy
    cellRect%dz = dz
    cellRect%sinrot = sinrot
    cellRect%cosrot = cosrot
    cellRect%xOrigin = xOrigin
    cellRect%yOrigin = yOrigin
    cellRect%zOrigin = zOrigin
    cellRect%ipvOrigin = ipv1
    !
    ! -- Compute (unscaled) cell edge velocities from face flows
    areax = dx * dz
    areay = dy * dz
    areaz = dx * dy
    ! cellRect%vx1 = cellDefn%faceflow(ipv1)/areax     ! kluge note: assuming porosity=1. and velmult=1. for now
    ! cellRect%vx2 = -cellDefn%faceflow(ipv3)/areax
    ! cellRect%vy1 = cellDefn%faceflow(ipv4)/areay
    ! cellRect%vy2 = -cellDefn%faceflow(ipv2)/areay
    ! cellRect%vz1 = cellDefn%faceflow(6)/areaz
    ! cellRect%vz2 = -cellDefn%faceflow(7)/areaz
    ! factor = cellDefn%velfactor/cellDefn%porosity
    factor = DONE / (cellDefn%retfactor * cellDefn%porosity)
    term = factor / areax
    cellRect%vx1 = cellDefn%faceflow(ipv1) * term
    cellRect%vx2 = -cellDefn%faceflow(ipv3) * term
    term = factor / areay
    cellRect%vy1 = cellDefn%faceflow(ipv4) * term
    cellRect%vy2 = -cellDefn%faceflow(ipv2) * term
    term = factor / areaz
    cellRect%vz1 = cellDefn%faceflow(6) * term
    cellRect%vz2 = -cellDefn%faceflow(7) * term
    !
    istatus = 1 ! kluge note: needed???
    return
    !
  end subroutine CellPolyToCellRect

  !> @brief Convert CellPoly representation to CellRectQuad if possible
  subroutine CellPolyToCellRectQuad(cellPoly, cellRectQuad, istatus)
    use CellRectQuadModule, only: CellRectQuadType, create_cellRectQuad
    use CellPolyModule, only: CellPolyType
    use UtilMiscModule, only: add_wrap, incr_wrap
    ! -- dummy
    type(CellPolyType), pointer :: cellPoly
    type(CellRectQuadType), pointer :: cellRectQuad
    integer :: istatus
    ! -- local
    integer :: i, irv, isc
    double precision :: qhalf, qdisttopbot, q1, q2, q4
    !
    call create_cellRectQuad(cellRectQuad)
    call cellRectQuad%init(cellPoly%cellDefn)
    ! kluge note: no check whether conversion is possible; assumes it is
    ! -- Translate and rotate the rect-quad cell into local coordinates with
    ! -- x varying from 0 to dx and y varying from 0 to dy. Choose the "south-
    ! -- west" rectangle vertex to be the local origin so that the rotation
    ! -- angle is zero if the cell already aligns with the model x and y
    ! -- coordinates.
    cellRectQuad%irvOrigin = cellRectQuad%get_irectvertSW() ! kluge note: no need to pass all that stuff in call below -- set internally in CellRectQuad
    call cellRectQuad%get_rectDimensionsRotation( &
      cellRectQuad%irvOrigin, cellRectQuad%xOrigin, &
      cellRectQuad%yOrigin, cellRectQuad%zOrigin, &
      cellRectQuad%dx, cellRectQuad%dy, &
      cellRectQuad%dz, cellRectQuad%sinrot, &
      cellRectQuad%cosrot)
    !
    ! -- Set the external and internal face flows used for subcells
    do i = 0, 3
      irv = add_wrap(cellRectQuad%irvOrigin, i, 4)
      isc = add_wrap(i, 3, 4)
      if (.not. cellRectQuad%isRefinedFace(irv)) then
        qhalf = 5d-1 * cellRectQuad%get_rectflow(1, irv)
        cellRectQuad%qextl2(isc) = qhalf
        isc = incr_wrap(isc, 4)
        cellRectQuad%qextl1(isc) = qhalf
      else
        cellRectQuad%qextl2(isc) = cellRectQuad%get_rectflow(1, irv)
        isc = incr_wrap(isc, 4)
        cellRectQuad%qextl1(isc) = cellRectQuad%get_rectflow(2, irv)
      end if
    end do
    qdisttopbot = 2.5d-1 * (cellRectQuad%cellDefn%get_distflow() &
                            + cellRectQuad%cellDefn%get_botflow() &
                            + cellRectQuad%cellDefn%get_topflow())
    q1 = qdisttopbot + cellRectQuad%qextl1(1) + cellRectQuad%qextl2(1)
    q2 = qdisttopbot + cellRectQuad%qextl1(2) + cellRectQuad%qextl2(2)
    q4 = qdisttopbot + cellRectQuad%qextl1(4) + cellRectQuad%qextl2(4)
    cellRectQuad%qintl(1) = -5d-1 * (q1 + 5d-1 * (q2 - q4))
    cellRectQuad%qintl(2) = cellRectQuad%qintl(1) + q1
    cellRectQuad%qintl(3) = cellRectQuad%qintl(2) + q2
    cellRectQuad%qintl(4) = cellRectQuad%qintl(1) - q4
    cellRectQuad%qintl(5) = cellRectQuad%qintl(1)
    !
    istatus = 1 ! kluge note: needed???
    return
    !
  end subroutine CellPolyToCellRectQuad

  !> @brief Get a metric indicating if a point is inside a cell polygon
  !!
  !! Return a metric for whether or not a 2D point
  !! is within a cell's polygon, including the cell
  !! boundary (negative value, no; otherwise, yes)
  !!
  !! kluge note: probably won't be needed in the long run
  !<
  function MetricForPointInCellPolygon(xpt, ypt, cellDefn) result(value)
    ! -- modules
    use CellDefnModule, only: CellDefnType
    ! -- dummy
    double precision :: xpt, ypt
    type(CellDefnType) :: cellDefn
    ! -- result
    double precision :: value
    ! -- local
    integer :: npolyverts, iv, iva, ivb
    double precision :: xa, ya, xb, yb, valueiv
    !
    npolyverts = cellDefn%npolyverts
    value = huge(1d0) ! kluge
    do iv = 1, npolyverts
      iva = iv
      ivb = iv + 1
      if (ivb .gt. npolyverts) ivb = 1
      xa = cellDefn%polyvert(iva)%x
      ya = cellDefn%polyvert(iva)%y
      xb = cellDefn%polyvert(ivb)%x
      yb = cellDefn%polyvert(ivb)%y
      valueiv = -(xb - xa) * (ypt - ya) + (yb - ya) * (xpt - xa)
      value = min(value, valueiv)
    end do
    !
    return
    !
  end function MetricForPointInCellPolygon

end module CellUtilModule
