!
module Ternary
  !
  ! -- Particle track arrays
  double precision, allocatable :: xtrk(:), ytrk(:), ztrk(:)
  ! -- "True" particle track arrays for comparison
  double precision, allocatable :: xtrktru(:), ytrktru(:), ztrktru(:)
  ! -- Vertex-coordinate arrays
  double precision, allocatable :: x_vert(:), y_vert(:)
  ! -- Node-coordinate arrays (for display purposes only)
  double precision, allocatable :: x_node(:), y_node(:)
  ! -- Cell-bottom elevation array (layer 0 = top)
  double precision, allocatable :: z_bot(:, :)
  ! -- Initial particle coordinate arrays
  double precision, allocatable :: xip(:), yip(:), zip(:)
  ! -- Final particle coordinate arrays
  double precision, allocatable :: xfp(:), yfp(:), zfp(:)
  double precision, allocatable :: xfptru(:), yfptru(:), zfptru(:)
  ! -- Travel time and path length arrays for particles
  double precision, allocatable :: traveltime(:), &
    traveltimetru(:), &
    pathlength(:), &
    pathlengthtru(:)
  ! -- Array of flows across polygon edges of each cell and across cell bottoms (layer 0 = top)
  double precision, allocatable :: flow_polygon(:, :, :), flux_bot(:, :)
  ! -- Arrays for components (x and y) of velocity at the vertices of a polygon
  double precision, allocatable :: vx_vert_polygon(:), vy_vert_polygon(:)
  ! -- Arrays that hold the two cell numbers associated with each interface
  integer, allocatable :: icell1_inter(:), icell2_inter(:)
  ! -- Arrays that hold the two vertex numbers associated with each interface
  integer, allocatable :: ivert1_inter(:), ivert2_inter(:)
  ! -- Arrays that hold the vertex and interface numbers for the polygon in each cell
  integer, allocatable :: ivert_polygon(:, :), iinter_polygon(:, :)
  ! -- "IFACE" arrays that hold cell numbers and face numbers for CHD conditions
  integer, allocatable :: ifacechd_layer(:), ifacechd_cell(:), ifacechd_face(:)
  ! -- Arrays to save entrance information
  integer, allocatable :: ilayersaved(:), &
                          icellsaved(:), &
                          itrisaved(:), &
                          numversaved(:)
  !
  ! -- Analytical solution coefficients
  double precision ca1, ca2, ca3, cb1, cb2
  ! -- Elements of the "velocity matrix," W
  double precision waa, wab, wba, wbb
  ! -- Coordinates of adjacent polygon-edge vertices and coordinate offsets,
  ! -- used in exit-point calculation
  double precision alpp1, betp1, alpp2, betp2, alppdiff, betpdiff
  ! -- "Canonical" velocity components at corners of triangular subcell
  double precision v0alp, v0bet, v1alp, v1bet, v2alp, v2bet
  ! -- Cumulative CPU time for finding exits
  double precision tcpufindexit
  ! -- Case index for analytical solution
  integer icase
  ! -- Number of solution repeats for timing
  integer nrptsoln
  ! -- Flag that indicates that the entrance face is being processed
  logical lenter
  ! -- Flag for suppression of output during timing
  logical :: lsupout = .true. ! kluge
  ! -- Copies of beti and alpexit     ! kluge
  double precision beticopy, alpexitcopy
  !
end module
