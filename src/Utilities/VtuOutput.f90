!> @brief VTK (VTU) output of results on polyhedral cells
!!
!! Used with the FULL_3D option of the DISU Package.  The cell polyhedra are
!! read from a geometry file (extension .v3d, written by VORO2MESH) and, at
!! each time step for which the dependent variable is saved, an ASCII VTU
!! file is written together with a PVD file that collects the time series.
!!
!! Geometry file format (free format):
!!   NCELLS n
!!   NPOINTS p
!!   x y z                (p lines)
!!   nf                   (for each cell, n times)
!!   nv i1 ... inv        (for each of the nf faces, 1-based point numbers)
!<
module VtuOutputModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error

  implicit none
  private
  public :: VtuOutputType

  type :: VtuOutputType
    character(len=LINELENGTH) :: root = '' !< root name of the output files
    integer(I4B) :: ncells = 0
    integer(I4B) :: npoints = 0
    integer(I4B) :: nstep = 0
    integer(I4B) :: nfacestream = 0
    real(DP), allocatable :: pts(:, :) !< (3,npoints)
    integer(I4B), allocatable :: nfaces(:) !< faces of each cell
    integer(I4B), allocatable :: facestream(:) !< VTK face stream of all cells
    integer(I4B), allocatable :: faceoff(:) !< end of the face stream of each cell
    integer(I4B), allocatable :: conn(:) !< unique points of each cell
    integer(I4B), allocatable :: connoff(:) !< end of connectivity of each cell
    real(DP), allocatable :: times(:)
    character(len=LINELENGTH), allocatable :: files(:)
    real(DP), allocatable :: vtimes(:) !< time series of the cell-center vector files
    character(len=LINELENGTH), allocatable :: vfiles(:)
  contains
    procedure :: init => vtu_init
    procedure :: write_step => vtu_write_step
    procedure :: da => vtu_da
  end type VtuOutputType

contains

  subroutine vtu_init(this, geofile, root)
    class(VtuOutputType) :: this
    character(len=*), intent(in) :: geofile
    character(len=*), intent(in) :: root
    integer(I4B) :: iu, ios, i, j, k, nf, nv, ip, ic, ifs, icn, np
    character(len=16) :: key
    integer(I4B), allocatable :: tmp(:), mark(:)
    !
    this%root = root
    open (newunit=iu, file=trim(geofile), status='old', action='read', &
          iostat=ios)
    if (ios /= 0) then
      write (errmsg, '(a)') 'Could not open the VTU geometry file "'// &
        trim(geofile)//'".'
      call store_error(errmsg, terminate=.true.)
    end if
    read (iu, *, iostat=ios) key, this%ncells
    read (iu, *, iostat=ios) key, this%npoints
    if (ios /= 0 .or. this%ncells < 1 .or. this%npoints < 4) then
      write (errmsg, '(a)') 'Invalid header in the VTU geometry file "'// &
        trim(geofile)//'". Expected NCELLS and NPOINTS.'
      call store_error(errmsg, terminate=.true.)
    end if
    allocate (this%pts(3, this%npoints))
    do i = 1, this%npoints
      read (iu, *, iostat=ios) this%pts(1:3, i)
      if (ios /= 0) call bad(geofile)
    end do
    allocate (this%nfaces(this%ncells), this%faceoff(this%ncells))
    allocate (this%connoff(this%ncells))
    allocate (mark(this%npoints))
    mark = 0
    ! first pass storage: grow arrays with a simple doubling strategy
    allocate (this%facestream(max(1024, 8 * this%ncells)))
    allocate (this%conn(max(1024, 8 * this%ncells)))
    ifs = 0
    icn = 0
    do ic = 1, this%ncells
      read (iu, *, iostat=ios) nf
      if (ios /= 0) call bad(geofile)
      this%nfaces(ic) = nf
      call put_fs(this%facestream, ifs, nf)
      np = 0
      do j = 1, nf
        read (iu, *, iostat=ios) nv
        if (ios /= 0) call bad(geofile)
        backspace (iu)
        if (allocated(tmp)) deallocate (tmp)
        allocate (tmp(nv + 1))
        read (iu, *, iostat=ios) tmp
        if (ios /= 0) call bad(geofile)
        call put_fs(this%facestream, ifs, nv)
        do k = 1, nv
          ip = tmp(k + 1)
          if (ip < 1 .or. ip > this%npoints) call bad(geofile)
          call put_fs(this%facestream, ifs, ip - 1)
          if (mark(ip) /= ic) then
            mark(ip) = ic
            call put_fs(this%conn, icn, ip - 1)
          end if
        end do
      end do
      this%faceoff(ic) = ifs
      this%connoff(ic) = icn
    end do
    close (iu)
    this%nfacestream = ifs
    call shrink(this%facestream, ifs)
    call shrink(this%conn, icn)
    allocate (this%times(0), this%files(0))
  end subroutine vtu_init

  subroutine bad(geofile)
    character(len=*), intent(in) :: geofile
    write (errmsg, '(a)') 'Error reading the VTU geometry file "'// &
      trim(geofile)//'".'
    call store_error(errmsg, terminate=.true.)
  end subroutine bad

  subroutine put_fs(a, n, v)
    integer(I4B), allocatable, intent(inout) :: a(:)
    integer(I4B), intent(inout) :: n
    integer(I4B), intent(in) :: v
    integer(I4B), allocatable :: b(:)
    if (n == size(a)) then
      allocate (b(2 * size(a)))
      b(1:n) = a(1:n)
      call move_alloc(b, a)
    end if
    n = n + 1
    a(n) = v
  end subroutine put_fs

  subroutine shrink(a, n)
    integer(I4B), allocatable, intent(inout) :: a(:)
    integer(I4B), intent(in) :: n
    integer(I4B), allocatable :: b(:)
    allocate (b(n))
    b = a(1:n)
    call move_alloc(b, a)
  end subroutine shrink

  !> @brief Write one VTU file with the scalar (and optional vector) results
  !!
  !! nodereduced maps the user cell number to the reduced number (<= 0 for
  !! inactive cells, which receive the value 1.0e30)
  !<
  subroutine vtu_write_step(this, totim, nodereduced, sname, svals, vname, &
                            vvals, vol, xyz)
    class(VtuOutputType) :: this
    real(DP), intent(in) :: totim
    integer(I4B), dimension(:), intent(in) :: nodereduced !< size 0: no reduction
    character(len=*), intent(in) :: sname
    real(DP), dimension(:), intent(in) :: svals
    character(len=*), intent(in), optional :: vname
    real(DP), dimension(:, :), intent(in), optional :: vvals
    real(DP), dimension(:), intent(in), optional :: vol
    real(DP), dimension(:, :), intent(in), optional :: xyz !< cell centers (3,reduced nodes)
    integer(I4B) :: iu, i, nr, ic, ios
    character(len=LINELENGTH) :: fname
    real(DP), parameter :: inactive = 1.0e30_DP
    !
    this%nstep = this%nstep + 1
    write (fname, '(a,a,i6.6,a)') trim(this%root), '_', this%nstep, '.vtu'
    open (newunit=iu, file=trim(fname), status='replace', action='write', &
          iostat=ios)
    if (ios /= 0) then
      write (errmsg, '(a)') 'Could not open the VTU output file "'// &
        trim(fname)//'".'
      call store_error(errmsg, terminate=.true.)
    end if
    write (iu, '(a)') '<?xml version="1.0"?>'
    write (iu, '(a)') '<VTKFile type="UnstructuredGrid" version="1.0" '// &
      'byte_order="LittleEndian">'
    write (iu, '(a)') '<UnstructuredGrid>'
    write (iu, '(a,i0,a,i0,a)') '<Piece NumberOfPoints="', this%npoints, &
      '" NumberOfCells="', this%ncells, '">'
    write (iu, '(a)') '<Points>'
    write (iu, '(a)') '<DataArray type="Float64" NumberOfComponents="3" '// &
      'format="ascii">'
    do i = 1, this%npoints
      write (iu, '(3(es24.16e3,1x))') this%pts(1:3, i)
    end do
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '</Points>'
    write (iu, '(a)') '<Cells>'
    write (iu, '(a)') '<DataArray type="Int64" Name="connectivity" '// &
      'format="ascii">'
    call write_ints(iu, this%conn)
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '<DataArray type="Int64" Name="offsets" format="ascii">'
    call write_ints(iu, this%connoff)
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
    call write_ints(iu, [(42, ic=1, this%ncells)])
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '<DataArray type="Int64" Name="faces" format="ascii">'
    call write_ints(iu, this%facestream)
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '<DataArray type="Int64" Name="faceoffsets" '// &
      'format="ascii">'
    call write_ints(iu, this%faceoff)
    write (iu, '(a)') '</DataArray>'
    write (iu, '(a)') '</Cells>'
    write (iu, '(a)') '<CellData>'
    ! -- scalar
    write (iu, '(a)') '<DataArray type="Float64" Name="'//trim(sname)// &
      '" format="ascii">'
    do ic = 1, this%ncells
      nr = reduced(ic)
      if (nr > 0) then
        write (iu, '(es24.16e3)') svals(nr)
      else
        write (iu, '(es24.16e3)') inactive
      end if
    end do
    write (iu, '(a)') '</DataArray>'
    ! -- vector
    if (present(vname) .and. present(vvals)) then
      write (iu, '(a)') '<DataArray type="Float64" Name="'//trim(vname)// &
        '" NumberOfComponents="3" format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) then
          write (iu, '(3(es24.16e3,1x))') vvals(1:3, nr)
        else
          write (iu, '(3(es24.16e3,1x))') 0.0_DP, 0.0_DP, 0.0_DP
        end if
      end do
      write (iu, '(a)') '</DataArray>'
    end if
    if (present(vol)) then
      write (iu, '(a)') '<DataArray type="Float64" Name="CELL_VOLUME" '// &
        'format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) then
          write (iu, '(es24.16e3)') vol(nr)
        else
          write (iu, '(es24.16e3)') 0.0_DP
        end if
      end do
      write (iu, '(a)') '</DataArray>'
    end if
    write (iu, '(a)') '</CellData>'
    write (iu, '(a)') '</Piece>'
    write (iu, '(a)') '</UnstructuredGrid>'
    write (iu, '(a)') '</VTKFile>'
    close (iu)
    ! -- time series
    call append_series(this, totim, fname)
    call write_pvd(this%root, this%times, this%files)
    ! -- vectors at the cell centers (point cloud, for glyphs in ParaView)
    if (present(vname) .and. present(vvals) .and. present(xyz)) then
      call write_vectors()
    end if
  contains
    subroutine write_vectors()
      character(len=LINELENGTH) :: vfname
      integer(I4B) :: iv, nact, k
      real(DP), allocatable :: t(:)
      character(len=LINELENGTH), allocatable :: f(:)
      !
      write (vfname, '(a,a,i6.6,a)') trim(this%root), '_vectors_', this%nstep, &
        '.vtu'
      open (newunit=iv, file=trim(vfname), status='replace', action='write', &
            iostat=ios)
      if (ios /= 0) return
      nact = 0
      do ic = 1, this%ncells
        if (reduced(ic) > 0) nact = nact + 1
      end do
      write (iv, '(a)') '<?xml version="1.0"?>'
      write (iv, '(a)') '<VTKFile type="UnstructuredGrid" version="1.0" '// &
        'byte_order="LittleEndian">'
      write (iv, '(a)') '<UnstructuredGrid>'
      write (iv, '(a,i0,a,i0,a)') '<Piece NumberOfPoints="', nact, &
        '" NumberOfCells="', nact, '">'
      write (iv, '(a)') '<Points>'
      write (iv, '(a)') '<DataArray type="Float64" NumberOfComponents="3" '// &
        'format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) write (iv, '(3(es24.16e3,1x))') xyz(1:3, nr)
      end do
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '</Points>'
      write (iv, '(a)') '<Cells>'
      write (iv, '(a)') '<DataArray type="Int64" Name="connectivity" '// &
        'format="ascii">'
      call write_ints(iv, [(k - 1, k=1, nact)])
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '<DataArray type="Int64" Name="offsets" format="ascii">'
      call write_ints(iv, [(k, k=1, nact)])
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
      call write_ints(iv, [(1, k=1, nact)])
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '</Cells>'
      write (iv, '(a)') '<PointData Vectors="'//trim(vname)//'" Scalars="'// &
        trim(sname)//'">'
      write (iv, '(a)') '<DataArray type="Float64" Name="'//trim(vname)// &
        '" NumberOfComponents="3" format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) write (iv, '(3(es24.16e3,1x))') vvals(1:3, nr)
      end do
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '<DataArray type="Float64" Name="'//trim(vname)// &
        '_MAGNITUDE" format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) write (iv, '(es24.16e3)') sqrt(sum(vvals(1:3, nr)**2))
      end do
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '<DataArray type="Float64" Name="'//trim(sname)// &
        '" format="ascii">'
      do ic = 1, this%ncells
        nr = reduced(ic)
        if (nr > 0) write (iv, '(es24.16e3)') svals(nr)
      end do
      write (iv, '(a)') '</DataArray>'
      write (iv, '(a)') '</PointData>'
      write (iv, '(a)') '</Piece>'
      write (iv, '(a)') '</UnstructuredGrid>'
      write (iv, '(a)') '</VTKFile>'
      close (iv)
      ! -- time series of the vector files
      k = 0
      if (allocated(this%vtimes)) k = size(this%vtimes)
      allocate (t(k + 1), f(k + 1))
      if (k > 0) then
        t(1:k) = this%vtimes
        f(1:k) = this%vfiles
      end if
      t(k + 1) = totim
      f(k + 1) = vfname
      call move_alloc(t, this%vtimes)
      call move_alloc(f, this%vfiles)
      call write_pvd(trim(this%root)//'_vectors', this%vtimes, this%vfiles)
    end subroutine write_vectors

    integer(I4B) function reduced(icell)
      integer(I4B), intent(in) :: icell
      if (size(nodereduced) == 0) then
        reduced = icell
      else
        reduced = nodereduced(icell)
      end if
    end function reduced
  end subroutine vtu_write_step

  subroutine append_series(this, totim, fname)
    class(VtuOutputType) :: this
    real(DP), intent(in) :: totim
    character(len=*), intent(in) :: fname
    real(DP), allocatable :: t(:)
    character(len=LINELENGTH), allocatable :: f(:)
    integer(I4B) :: n, k
    n = size(this%times)
    allocate (t(n + 1), f(n + 1))
    t(1:n) = this%times
    do k = 1, n
      f(k) = this%files(k)
    end do
    t(n + 1) = totim
    f(n + 1) = fname
    call move_alloc(t, this%times)
    call move_alloc(f, this%files)
  end subroutine append_series

  subroutine write_pvd(root, times, files)
    character(len=*), intent(in) :: root
    real(DP), dimension(:), intent(in) :: times
    character(len=LINELENGTH), dimension(:), intent(in) :: files
    integer(I4B) :: iu, k, ios, p
    character(len=LINELENGTH) :: base
    open (newunit=iu, file=trim(root)//'.pvd', status='replace', &
          action='write', iostat=ios)
    if (ios /= 0) return
    write (iu, '(a)') '<?xml version="1.0"?>'
    write (iu, '(a)') '<VTKFile type="Collection" version="0.1" '// &
      'byte_order="LittleEndian">'
    write (iu, '(a)') '<Collection>'
    do k = 1, size(times)
      base = files(k)
      p = max(index(base, '/', back=.true.), index(base, '\', back=.true.))
      write (iu, '(a,es24.16e3,a)') '<DataSet timestep="', times(k), &
        '" group="" part="0" file="'//trim(base(p + 1:))//'"/>'
    end do
    write (iu, '(a)') '</Collection>'
    write (iu, '(a)') '</VTKFile>'
    close (iu)
  end subroutine write_pvd

  subroutine write_ints(iu, a)
    integer(I4B), intent(in) :: iu
    integer(I4B), dimension(:), intent(in) :: a
    integer(I4B) :: i, j
    do i = 1, size(a), 12
      j = min(i + 11, size(a))
      write (iu, '(12(i0,1x))') a(i:j)
    end do
  end subroutine write_ints

  subroutine vtu_da(this)
    class(VtuOutputType) :: this
    if (allocated(this%pts)) deallocate (this%pts)
    if (allocated(this%nfaces)) deallocate (this%nfaces)
    if (allocated(this%facestream)) deallocate (this%facestream)
    if (allocated(this%faceoff)) deallocate (this%faceoff)
    if (allocated(this%conn)) deallocate (this%conn)
    if (allocated(this%connoff)) deallocate (this%connoff)
    if (allocated(this%times)) deallocate (this%times)
    if (allocated(this%files)) deallocate (this%files)
    if (allocated(this%vtimes)) deallocate (this%vtimes)
    if (allocated(this%vfiles)) deallocate (this%vfiles)
  end subroutine vtu_da

end module VtuOutputModule
