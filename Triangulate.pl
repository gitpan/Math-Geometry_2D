#!/ap/tools/bin/perl

#         This program is an implementation of a fast polygon
# triangulation algorithm based on the paper "A simple and fast
# incremental randomized algorithm for computing trapezoidal
# decompositions and for triangulating polygons" by Raimund Seidel.
#
#         The algorithm handles simple polygons with holes. The input is
# specified as contours. The outermost contour is anti-clockwise, while
# all the inner contours must be clockwise. No point should be repeated
# in the input. A sample input file 'data_1' is provided.
#
#         The output is a reference to a list of triangles. Each triangle
# is ar reference to an array fo three points, each point is a reference
# to an array holdign the x and y coordinates of the point.
# The number of output triangles produced for a polygon with n points is,
#         (n - 2) + 2*(#holes)
#
#         The program is a translation to perl of the C program written by
# Narkhede A. and Manocha D., Fast polygon triangulation algorithm based
# on Seidel's Algorithm, UNC-CH, 1994.
# Note that in this perl version, there are no statically allocated arrays
# so the only limit is the amount of (virtual) memory available.
#
# See also:
#
#   R. Seidel
#     "A simple and Fast Randomized Algorithm for Computing Trapezoidal
#      Decompositions and for Triangulating Polygons"
#     "Computational Geometry Theory & Applications"
#      Number = 1, Year 1991, Volume 1, Pages 51-64
#
#   J. O'Rourke
#     "Computational Geometry in {C}"
#      Cambridge University Press  - 1994
#
# Input specified as a contour with the restrictions mentioned above:
#   - first polygon is the outer shape and must be anti-clockwise.
#   - next polygons are inner shapels (holes) must be clockwise.
#   - Inner and outer shapes must be simple .
#
# Every contour is specified by giving all its points in order. No
# point shoud be repeated. i.e. if the outer contour is a square,
# only the four distinct endpoints shopudl be specified in order.
#
# Returns a reference to an array holding the triangles.
#

use strict;
use Carp;
use POSIX;

my $C_EPS     = 1e-10; # tolerance value: Used for making
                       # all decisions about collinearity or
                       # left/right of segment. Decrease
                       # this value if the input points are
                       # spaced very close together

my $INFINITY = 1<<29;

my $TRUE  = 1;
my $FALSE = 0;

my $T_X    = 1;
my $T_Y    = 2;
my $T_SINK = 3;

my $ST_VALID   = 1;
my $ST_INVALID = 2;

my $FIRSTPT = 1;
my $LASTPT  = 2;

my $S_LEFT  = 1;
my $S_RIGHT = 2;

my $SP_SIMPLE_LRUP =  1; # for splitting trapezoids
my $SP_SIMPLE_LRDN =  2;
my $SP_2UP_2DN     =  3;
my $SP_2UP_LEFT    =  4;
my $SP_2UP_RIGHT   =  5;
my $SP_2DN_LEFT    =  6;
my $SP_2DN_RIGHT   =  7;
my $SP_NOSPLIT     = -1;

my $TRI_LHS = 1;
my $TRI_RHS = 2;
my $TR_FROM_UP = 1;    # for traverse-direction
my $TR_FROM_DN = 2;

my $choose_idx = 1;
my @permute;
my $q_idx;
my $tr_idx;
my @qs;       # Query structure
my @tr;       # Trapezoid structure
my @seg;      # Segment table

my @mchain; # Table to hold all the monotone
            # polygons . Each monotone polygon
            # is a circularly linked list
my @vert;   # chain init. information. This
            # is used to decide which
            # monotone polygon to split if
            # there are several other
            # polygons touching at the same
            # vertex
my @mon;    # contains position of any vertex in
            # the monotone chain for the polygon
my @visited;
my @op;     # contains the resulting list of triangles
            # and their vertex number
my ($chain_idx, $op_idx, $mon_idx);

sub TriangulatePolygon {

  $choose_idx = 1;
  @seg = ();
  @mchain = ();
  @vert = ();
  @mon = ();
  @visited = ();
  @op = ();

  my ($polygonrefs) = @_;
  my @polygons = @{$polygonrefs};

  my $ccount = 0;
  my $i = 1;
  while ($ccount < @polygons) {
    my @vertexarray = @{$polygons[$ccount]};
    my $npoints     = @vertexarray;
    my $first = $i;
    my $last  = $first + $npoints - 1;
    for (my $j = 0; $j < $npoints; $j++, $i++) {
      my @vertex = @{$vertexarray[$j]};
      $seg[$i]{v0}{x} = $vertex[0];
      $seg[$i]{v0}{y} = $vertex[1];
      if ($i == $last) {
        $seg[$i]{next} = $first;
        $seg[$i]{prev} = $i-1;
        my %tmp = %{$seg[$i]{v0}};
        $seg[$i-1]{v1} = \%tmp;
      } elsif ($i == $first) {
        $seg[$i]{next} = $i+1;
        $seg[$i]{prev} = $last;
        my %tmp = %{$seg[$i]{v0}};
        $seg[$last]{v1} = \%tmp;
      } else {
        $seg[$i]{prev} = $i-1;
        $seg[$i]{next} = $i+1;
        my %tmp = %{$seg[$i]{v0}};
        $seg[$i-1]{v1} = \%tmp;
      }
      $seg[$i]{is_inserted} = $FALSE;
    }
    $ccount++;
  }

  my $n = $i-1;

  _generate_random_ordering($n);
  _construct_trapezoids($n);
  my $nmonpoly = _monotonate_trapezoids($n);
  my $ntriangles = _triangulate_monotone_polygons($n, $nmonpoly);
  # now get the coordinates for all the triangles
  my @result;
  for (my $i = 0; $i < $ntriangles; $i++) {
    my @vertices = @{$op[$i]};
    my $triangle = [[$seg[$vertices[0]]{v0}{x},$seg[$vertices[0]]{v0}{y}],
                    [$seg[$vertices[1]]{v0}{x},$seg[$vertices[1]]{v0}{y}],
                    [$seg[$vertices[2]]{v0}{x},$seg[$vertices[2]]{v0}{y}]];
    push @result,$triangle;
  }
  return [@result];;
}

# Generate a random permutation of the segments 1..n
sub _generate_random_ordering {
  @permute = ();
  my ($n) = @_;
  my @input;
  for (my $i = 1 ; $i <= $n ; $i++) {
    $input[$i] = $i;
  }
  my $i = 1;
  for (my $i = 1 ; $i <= $n ; $i++) {
    my $m = int rand($#input) + 1;
    $permute[$i] = $input[$m];
    splice @input,$m,1;
  }
}

# Return the next segment in the generated random ordering of all the
# segments in S
sub _choose_segment {
  return $permute[$choose_idx++];
}

# Return a new node to be added into the query tree
sub _newnode {
  return $q_idx++;
}

# Return a free trapezoid
sub _newtrap {
  $tr[$tr_idx]{lseg} = -1;
  $tr[$tr_idx]{rseg} = -1;
  $tr[$tr_idx]{state} = $ST_VALID;
  return $tr_idx++;
}

# Floating point number comparison
sub _fp_equal {
  my ($X, $Y, $POINTS) = @_;
  my ($tX, $tY);
  $tX = sprintf("%.${POINTS}g", $X);
  $tY = sprintf("%.${POINTS}g", $Y);
  return $tX eq $tY;
}

# Return the maximum of the two points
sub _max {
  my ($v0_ref, $v1_ref) = @_;
  my %v0   = %{$v0_ref};
  my %v1   = %{$v1_ref};
  if ($v0{y} > $v1{y} + $C_EPS) {
    return \%v0;
  } elsif (_fp_equal($v0{y}, $v1{y}, $precision)) {
    if ($v0{x} > $v1{x} + $C_EPS) {
      return \%v0;
    } else {
      return \%v1;
    }
  } else {
    return \%v1;
  }
}

# Return the minimum of the two points
sub _min {
  my ($v0_ref, $v1_ref) = @_;
  my %v0   = %{$v0_ref};
  my %v1   = %{$v1_ref};
  if ($v0{y} < $v1{y} - $C_EPS) {
    return \%v0;
  } elsif (_fp_equal($v0{y}, $v1{y}, $precision)) {
    if ($v0{x} < $v1{x}) {
      return \%v0;
    } else {
      return \%v1;
    }
  } else {
    return \%v1;
  }
}

sub _greater_than {
  my ($v0_ref, $v1_ref) = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  if ($v0{y} > $v1{y} + $C_EPS) {
    return 1;
  } elsif ($v0{y} < $v1{y} - $C_EPS) {
    return 0;
  } else {
    return ($v0{x} > $v1{x});
  }
}

sub _equal_to {
  my ($v0_ref, $v1_ref) = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  return ( _fp_equal($v0{y}, $v1{y}, $precision) &&
           _fp_equal($v0{x}, $v1{x}, $precision) );
}

sub _greater_than_equal_to {
  my ($v0_ref, $v1_ref) = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  if ($v0{y} > $v1{y} + $C_EPS) {
    return 1;
  } elsif ($v0{y} < $v1{y} - $C_EPS) {
    return 0;
  } else {
    return ($v0{x} >= $v1{x});
  }
}

sub _less_than {
  my ($v0_ref, $v1_ref) = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  if ($v0{y} < $v1{y} - $C_EPS) {
    return 1;
  } elsif ($v0{y} > $v1{y} + $C_EPS) {
    return 0;
  } else {
    return ($v0{x} < $v1{x});
  }
}

# Initilialise the query structure (Q) and the trapezoid table (T)
# when the first segment is added to start the trapezoidation. The
# query-tree starts out with 4 trapezoids, one S-node and 2 Y-nodes
#
#                4
#   -----------------------------------
#                 \
#       1          \        2
#                   \
#   -----------------------------------
#                3
#

sub _init_query_structure {
  my ($segnum) = @_;

  my ($i1,$i2,$i3,$i4,$i5,$i6,$i7,$root);
  my ($t1,$t2,$t3,$t4);

  @qs = ();
  @tr = ();

  $q_idx  =  $tr_idx = 1;

  $i1 = _newnode();
  $qs[$i1]{nodetype} = $T_Y;

  my %tmpmax = %{_max($seg[$segnum]{v0}, $seg[$segnum]{v1})}; # root
  $qs[$i1]{yval} = {x => $tmpmax{x} , y => $tmpmax{y}};
  $root = $i1;

  $qs[$i1]{right} = $i2 = _newnode();
  $qs[$i2]{nodetype} = $T_SINK;
  $qs[$i2]{parent} = $i1;

  $qs[$i1]{left} = $i3 = _newnode();
  $qs[$i3]{nodetype} = $T_Y;
  my %tmpmin = %{_min($seg[$segnum]{v0}, $seg[$segnum]{v1})}; # root
  $qs[$i3]{yval} = {x => $tmpmin{x} , y => $tmpmin{y}};
  $qs[$i3]{parent} = $i1;

  $qs[$i3]{left} = $i4 = _newnode();
  $qs[$i4]{nodetype} = $T_SINK;
  $qs[$i4]{parent} = $i3;

  $qs[$i3]{right} = $i5 = _newnode();
  $qs[$i5]{nodetype} = $T_X;
  $qs[$i5]{segnum} = $segnum;
  $qs[$i5]{parent} = $i3;

  $qs[$i5]{left} = $i6 = _newnode();
  $qs[$i6]{nodetype} = $T_SINK;
  $qs[$i6]{parent} = $i5;

  $qs[$i5]{right} = $i7 = _newnode();
  $qs[$i7]{nodetype} = $T_SINK;
  $qs[$i7]{parent} = $i5;

  $t1 = _newtrap();    # middle left
  $t2 = _newtrap();    # middle right
  $t3 = _newtrap();    # bottom-most
  $t4 = _newtrap();    # topmost

  $tr[$t1]{hi} = {x => $qs[$i1]{yval}{x} , y => $qs[$i1]{yval}{y}};
  $tr[$t2]{hi} = {x => $qs[$i1]{yval}{x} , y => $qs[$i1]{yval}{y}};
  $tr[$t4]{lo} = {x => $qs[$i1]{yval}{x} , y => $qs[$i1]{yval}{y}};
  $tr[$t1]{lo} = {x => $qs[$i3]{yval}{x} , y => $qs[$i3]{yval}{y}};
  $tr[$t2]{lo} = {x => $qs[$i3]{yval}{x} , y => $qs[$i3]{yval}{y}};
  $tr[$t3]{hi} = {x => $qs[$i3]{yval}{x} , y => $qs[$i3]{yval}{y}};
  $tr[$t4]{hi} = {x =>      $INFINITY , y =>      $INFINITY};
  $tr[$t3]{lo} = {x => -1 * $INFINITY , y => -1 * $INFINITY};
  $tr[$t1]{rseg} = $tr[$t2]{lseg} = $segnum;
  $tr[$t1]{u0} = $tr[$t2]{u0} = $t4;
  $tr[$t1]{d0} = $tr[$t2]{d0} = $t3;
  $tr[$t4]{d0} = $tr[$t3]{u0} = $t1;
  $tr[$t4]{d1} = $tr[$t3]{u1} = $t2;

  $tr[$t1]{sink} = $i6;
  $tr[$t2]{sink} = $i7;
  $tr[$t3]{sink} = $i4;
  $tr[$t4]{sink} = $i2;

  $tr[$t1]{state} = $tr[$t2]{state} = $ST_VALID;
  $tr[$t3]{state} = $tr[$t4]{state} = $ST_VALID;

  $qs[$i2]{trnum} = $t4;
  $qs[$i4]{trnum} = $t3;
  $qs[$i6]{trnum} = $t1;
  $qs[$i7]{trnum} = $t2;

  $seg[$segnum]{is_inserted} = $TRUE;
  return $root;
}

# Update the roots stored for each of the endpoints of the segment.
# This is done to speed up the location-query for the endpoint when
# the segment is inserted into the trapezoidation subsequently
#
sub _find_new_roots {
  my ($segnum) = @_;

  return if ($seg[$segnum]{is_inserted});

  $seg[$segnum]{root0} = _locate_endpoint($seg[$segnum]{v0}, $seg[$segnum]{v1}, $seg[$segnum]{root0});
  $seg[$segnum]{root0} = $tr[$seg[$segnum]{root0}]{sink};

  $seg[$segnum]{root1} = _locate_endpoint($seg[$segnum]{v1}, $seg[$segnum]{v0}, $seg[$segnum]{root1});
  $seg[$segnum]{root1} = $tr[$seg[$segnum]{root1}]{sink};
}

# Main routine to perform trapezoidation
sub _construct_trapezoids {
  my ($nseg) = @_; #

  # Add the first segment and get the query structure and trapezoid
  # list initialised

  my $root = _init_query_structure(_choose_segment());

  for (my $i = 1 ; $i <= $nseg; $i++) {
    $seg[$i]{root0} = $seg[$i]{root1} = $root;
  }
  for (my $h = 1; $h <= _math_logstar_n($nseg); $h++) {
    for (my $i = _math_N($nseg, $h -1) + 1; $i <= _math_N($nseg, $h); $i++) {
      _add_segment(_choose_segment());
    }
    # Find a new root for each of the segment endpoints
    for (my $i = 1; $i <= $nseg; $i++) {
      _find_new_roots($i);
    }
  }
  for (my $i = _math_N($nseg, _math_logstar_n($nseg)) + 1; $i <= $nseg; $i++) {
    _add_segment(_choose_segment());
  }
}

# Add in the new segment into the trapezoidation and update Q and T
# structures. First locate the two endpoints of the segment in the
# Q-structure. Then start from the topmost trapezoid and go down to
# the  lower trapezoid dividing all the trapezoids in between .
#

sub _add_segment {
  my ($segnum) = @_;

  my ($tu, $tl, $sk, $tfirst, $tlast, $tnext);
  my ($tfirstr, $tlastr, $tfirstl, $tlastl);
  my ($i1, $i2, $t, $t1, $t2, $tn);
  my $tritop = 0;
  my $tribot = 0;
  my $is_swapped = 0;
  my $tmptriseg;
  my %s = %{$seg[$segnum]};

  if (_greater_than($s{v1}, $s{v0})) { # Get higher vertex in v0
    my %tmp;
    %tmp   = %{$s{v0}};
    $s{v0} = {x => $s{v1}{x} , y => $s{v1}{y}};
    $s{v1} = {x =>   $tmp{x} , y =>   $tmp{y}};
    my $tmp   = $s{root0};
    $s{root0} = $s{root1};
    $s{root1} = $tmp;
    $is_swapped = 1;
  }

  if (($is_swapped) ? !_inserted($segnum, $LASTPT) :
       !_inserted($segnum, $FIRSTPT)) { # insert v0 in the tree
    my $tmp_d;

    $tu = _locate_endpoint($s{v0}, $s{v1}, $s{root0});
    $tl = _newtrap();          # tl is the new lower trapezoid
    $tr[$tl]{state} = $ST_VALID;
    my %tmp = %{$tr[$tu]};
    my %tmphi = %{$tmp{hi}};
    my %tmplo = %{$tmp{lo}};
    $tr[$tl] = \%tmp;
    $tr[$tl]{hi} = {x => $tmphi{x} , y => $tmphi{y}};
    $tr[$tl]{lo} = {x => $tmplo{x} , y => $tmplo{y}};
    $tr[$tu]{lo} = {x => $s{v0}{x} , y => $s{v0}{y}};
    $tr[$tl]{hi} = {x => $s{v0}{x} , y => $s{v0}{y}};
    $tr[$tu]{d0} = $tl;
    $tr[$tu]{d1} = 0;
    $tr[$tl]{u0} = $tu;
    $tr[$tl]{u1} = 0;

    if ((($tmp_d = $tr[$tl]{d0}) > 0) && ($tr[$tmp_d]{u0} == $tu)) {
      $tr[$tmp_d]{u0} = $tl;
    }
    if ((($tmp_d = $tr[$tl]{d0}) > 0) && ($tr[$tmp_d]{u1} == $tu)) {
      $tr[$tmp_d]{u1} = $tl;
    }

    if ((($tmp_d = $tr[$tl]{d1}) > 0) && ($tr[$tmp_d]{u0} == $tu)) {
      $tr[$tmp_d]{u0} = $tl;
    }
    if ((($tmp_d = $tr[$tl]{d1}) > 0) && ($tr[$tmp_d]{u1} == $tu)) {
      $tr[$tmp_d]{u1} = $tl;
    }

    # Now update the query structure and obtain the sinks for the
    # two trapezoids

    $i1 = _newnode();          # Upper trapezoid sink
    $i2 = _newnode();          # Lower trapezoid sink
    $sk = $tr[$tu]{sink};

    $qs[$sk]{nodetype} = $T_Y;
    $qs[$sk]{yval}     = {x => $s{v0}{x} , y=> $s{v0}{y}};
    $qs[$sk]{segnum}   = $segnum;  # not really reqd ... maybe later
    $qs[$sk]{left}     = $i2;
    $qs[$sk]{right}    = $i1;

    $qs[$i1]{nodetype} = $T_SINK;
    $qs[$i1]{trnum}    = $tu;
    $qs[$i1]{parent}   = $sk;

    $qs[$i2]{nodetype} = $T_SINK;
    $qs[$i2]{trnum}    = $tl;
    $qs[$i2]{parent}   = $sk;

    $tr[$tu]{sink} = $i1;
    $tr[$tl]{sink} = $i2;
    $tfirst = $tl;
  } else {  # v0 already present
            # Get the topmost intersecting trapezoid
    $tfirst = _locate_endpoint($s{v0}, $s{v1}, $s{root0});
    $tritop = 1;
  }


  if (($is_swapped) ? !_inserted($segnum, $FIRSTPT) :
       !_inserted($segnum, $LASTPT)) { # insert v1 in the tree
    my $tmp_d;

    $tu = _locate_endpoint($s{v1}, $s{v0}, $s{root1});
    $tl = _newtrap();         # tl is the new lower trapezoid
    $tr[$tl]{state} = $ST_VALID;
    my %tmp = %{$tr[$tu]};
    my %tmphi = %{$tmp{hi}};
    my %tmplo = %{$tmp{lo}};
    $tr[$tl] = \%tmp;
    $tr[$tl]{hi} = {x => $tmphi{x} , y => $tmphi{y}};
    $tr[$tl]{lo} = {x => $tmplo{x} , y => $tmplo{y}};
    $tr[$tu]{lo} = {x => $s{v1}{x} , y => $s{v1}{y}};
    $tr[$tl]{hi} = {x => $s{v1}{x} , y => $s{v1}{y}};
    $tr[$tu]{d0} = $tl;
    $tr[$tu]{d1} = 0;
    $tr[$tl]{u0} = $tu;
    $tr[$tl]{u1} = 0;

    if ((($tmp_d = $tr[$tl]{d0}) > 0) && ($tr[$tmp_d]{u0} == $tu)) {
      $tr[$tmp_d]{u0} = $tl;
    }
    if ((($tmp_d = $tr[$tl]{d0}) > 0) && ($tr[$tmp_d]{u1} == $tu)) {
      $tr[$tmp_d]{u1} = $tl;
    }

    if ((($tmp_d = $tr[$tl]{d1}) > 0) && ($tr[$tmp_d]{u0} == $tu)) {
      $tr[$tmp_d]{u0} = $tl;
    }
    if ((($tmp_d = $tr[$tl]{d1}) > 0) && ($tr[$tmp_d]{u1} == $tu)) {
      $tr[$tmp_d]{u1} = $tl;
    }

    # Now update the query structure and obtain the sinks for the
    # two trapezoids

    $i1 = _newnode();          # Upper trapezoid sink
    $i2 = _newnode();          # Lower trapezoid sink
    $sk = $tr[$tu]{sink};

    $qs[$sk]{nodetype} = $T_Y;
    $qs[$sk]{yval}     = {x => $s{v1}{x} , y => $s{v1}{y}};
    $qs[$sk]{segnum}   = $segnum;   # not really reqd ... maybe later
    $qs[$sk]{left}     = $i2;
    $qs[$sk]{right}    = $i1;

    $qs[$i1]{nodetype} = $T_SINK;
    $qs[$i1]{trnum}    = $tu;
    $qs[$i1]{parent}   = $sk;

    $qs[$i2]{nodetype} = $T_SINK;
    $qs[$i2]{trnum}    = $tl;
    $qs[$i2]{parent}   = $sk;

    $tr[$tu]{sink} = $i1;
    $tr[$tl]{sink} = $i2;
    $tlast = $tu;
  } else {  # v1 already present
            # Get the lowermost intersecting trapezoid
    $tlast = _locate_endpoint($s{v1}, $s{v0}, $s{root1});
    $tribot = 1;
  }

  # Thread the segment into the query tree creating a new X-node
  # First, split all the trapezoids which are intersected by s into
  # two

  $t = $tfirst;               # topmost trapezoid

  while (($t > 0) &&
         _greater_than_equal_to($tr[$t]{lo}, $tr[$tlast]{lo})) {
                              # traverse from top to bot
    my ($t_sav, $tn_sav);
    $sk = $tr[$t]{sink};
    $i1 = _newnode();          # left trapezoid sink
    $i2 = _newnode();          # right trapezoid sink

    $qs[$sk]{nodetype} = $T_X;
    $qs[$sk]{segnum}   = $segnum;
    $qs[$sk]{left}     = $i1;
    $qs[$sk]{right}    = $i2;

    $qs[$i1]{nodetype} = $T_SINK;   # left trapezoid (use existing one)
    $qs[$i1]{trnum}    = $t;
    $qs[$i1]{parent}   = $sk;

    $qs[$i2]{nodetype} = $T_SINK;   # right trapezoid (allocate new)
    $qs[$i2]{trnum}    = $tn = _newtrap();
    $tr[$tn]{state}    = $ST_VALID;
    $qs[$i2]{parent}   = $sk;

    if ($t == $tfirst) {
      $tfirstr = $tn;
    }
    if (_equal_to($tr[$t]{lo}, $tr[$tlast]{lo})) {
      $tlastr = $tn;
    }

    my %tmp = %{$tr[$t]};
    my %tmphi = %{$tmp{hi}};
    my %tmplo = %{$tmp{lo}};
    $tr[$tn] = \%tmp;
    $tr[$tn]{hi} = {x => $tmphi{x} , y => $tmphi{y}};
    $tr[$tn]{lo} = {x => $tmplo{x} , y => $tmplo{y}};
    $tr[$t]{sink} = $i1;
    $tr[$tn]{sink} = $i2;
    $t_sav  = $t;
    $tn_sav = $tn;

    # error

    if (($tr[$t]{d0} <= 0) && ($tr[$t]{d1} <= 0)) {  # case cannot arise
      print "add_segment: error\n";

    # only one trapezoid below. partition t into two and make the
    # two resulting trapezoids t and tn as the upper neighbours of
    # the sole lower trapezoid

    } elsif (($tr[$t]{d0} > 0) && ($tr[$t]{d1} <= 0)) { # Only one trapezoid below
      if (($tr[$t]{u0} > 0) && ($tr[$t]{u1} > 0)) {     # continuation of a chain from abv.
        if ($tr[$t]{usave} > 0) {                 # three upper neighbours
          if ($tr[$t]{uside} == $S_LEFT) {
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = -1;
            $tr[$tn]{u1} = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0}  = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
            $tr[$tr[$tn]{u1}]{d0} = $tn;
          } else {                                # intersects in the right
            $tr[$tn]{u1} = -1;
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = $tr[$t]{u0};
            $tr[$t]{u0}  = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0} = $t;
            $tr[$tr[$t]{u1}]{d0} = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
          }

          $tr[$t]{usave} = $tr[$tn]{usave} = 0;
        } else {                                  # No usave.... simple case
          $tr[$tn]{u0} = $tr[$t]{u1};
          $tr[$t]{u1}  = $tr[$tn]{u1} = -1;
          $tr[$tr[$tn]{u0}]{d0} = $tn;
        }
      } else {                              # fresh seg. or upward cusp
        my $tmp_u = $tr[$t]{u0};
        my ($td0, $td1);
        if ((($td0 = $tr[$tmp_u]{d0}) > 0) &&
            (($td1 = $tr[$tmp_u]{d1}) > 0)) {  # upward cusp
          if (($tr[$td0]{rseg} > 0) &&
              !_is_left_of($tr[$td0]{rseg}, $s{v1})) {
            $tr[$t]{u0} = $tr[$t]{u1} = $tr[$tn]{u1} = -1;
            $tr[$tr[$tn]{u0}]{d1} = $tn;
          } else {   # cusp going leftwards
            $tr[$tn]{u0} = $tr[$tn]{u1} = $tr[$t]{u1} = -1;
            $tr[$tr[$t]{u0}]{d0} = $t;
          }
        } else {     # fresh segment
          $tr[$tr[$t]{u0}]{d0} = $t;
          $tr[$tr[$t]{u0}]{d1} = $tn;
        }
      }

      if (_fp_equal($tr[$t]{lo}{y}, $tr[$tlast]{lo}{y}, $precision) &&
          _fp_equal($tr[$t]{lo}{x}, $tr[$tlast]{lo}{x}, $precision) && $tribot) {
        # bottom forms a triangle

        if ($is_swapped) {
          $tmptriseg = $seg[$segnum]{prev};
        } else {
          $tmptriseg = $seg[$segnum]{next};
        }

        if (($tmptriseg > 0) && _is_left_of($tmptriseg, $s{v0})) { # L-R downward cusp
          $tr[$tr[$t]{d0}]{u0} = $t;
          $tr[$tn]{d0} = $tr[$tn]{d1} = -1;
        } else { # R-L downward cusp
          $tr[$tr[$tn]{d0}]{u1} = $tn;
          $tr[$t]{d0} = $tr[$t]{d1} = -1;
        }
      } else {
        if (($tr[$tr[$t]{d0}]{u0} > 0) && ($tr[$tr[$t]{d0}]{u1} > 0)) {
          if ($tr[$tr[$t]{d0}]{u0} == $t) {  # passes thru LHS
            $tr[$tr[$t]{d0}]{usave} = $tr[$tr[$t]{d0}]{u1};
            $tr[$tr[$t]{d0}]{uside} = $S_LEFT;
          } else {
            $tr[$tr[$t]{d0}]{usave} = $tr[$tr[$t]{d0}]{u0};
            $tr[$tr[$t]{d0}]{uside} = $S_RIGHT;
          }
        }
        $tr[$tr[$t]{d0}]{u0} = $t;
        $tr[$tr[$t]{d0}]{u1} = $tn;
      }

      $t = $tr[$t]{d0};

    } elsif (($tr[$t]{d0} <= 0) && ($tr[$t]{d1} > 0)) {  # Only one trapezoid below
      if (($tr[$t]{u0} > 0) && ($tr[$t]{u1} > 0)) {      # continuation of a chain from abv.
        if ($tr[$t]{usave} > 0) {     # three upper neighbours
          if ($tr[$t]{uside} == $S_LEFT) {
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = -1;
            $tr[$tn]{u1} = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0}  = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
            $tr[$tr[$tn]{u1}]{d0} = $tn;
          } else {  # intersects in the right
            $tr[$tn]{u1} = -1;
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = $tr[$t]{u0};
            $tr[$t]{u0}  = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0}  = $t;
            $tr[$tr[$t]{u1}]{d0}  = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
          }

          $tr[$t]{usave} = $tr[$tn]{usave} = 0;

        } else {  # No usave.... simple case
          $tr[$tn]{u0} = $tr[$t]{u1};
          $tr[$t]{u1} = $tr[$tn]{u1} = -1;
          $tr[$tr[$tn]{u0}]{d0} = $tn;
        }
      } else {  # fresh seg. or upward cusp
        my $tmp_u = $tr[$t]{u0};
        my ($td0,$td1);
        if ((($td0 = $tr[$tmp_u]{d0}) > 0) &&
            (($td1 = $tr[$tmp_u]{d1}) > 0)) {    # upward cusp
          if (($tr[$td0]{rseg} > 0) &&
              !_is_left_of($tr[$td0]{rseg}, $s{v1})) {
              $tr[$t]{u0} = $tr[$t]{u1} = $tr[$tn]{u1} = -1;
              $tr[$tr[$tn]{u0}]{d1} = $tn;
          } else {
            $tr[$tn]{u0} = $tr[$tn]{u1} = $tr[$t]{u1} = -1;
            $tr[$tr[$t]{u0}]{d0} = $t;
          }
        } else {  # fresh segment
          $tr[$tr[$t]{u0}]{d0} = $t;
          $tr[$tr[$t]{u0}]{d1} = $tn;
        }
      }

      if (_fp_equaL($tr[$t]{lo}{y}, $tr[$tlast]{lo}{y}, $precision) &&
          _fp_equal($tr[$t]{lo}{x}, $tr[$tlast]{lo}{x}, $precision) && $tribot) {
        # bottom forms a triangle
        my $tmpseg;

        if ($is_swapped) {
          $tmptriseg = $seg[$segnum]{prev};
        } else {
          $tmptriseg = $seg[$segnum]{next};
        }

        if (($tmpseg > 0) && _is_left_of($tmpseg, $s{v0})) {
          # L-R downward cusp
          $tr[$tr[$t]{d1}]{u0} = $t;
          $tr[$tn]{d0} = $tr[$tn]{d1} = -1;
        } else {
          # R-L downward cusp
          $tr[$tr[$tn]{d1}]{u1} = $tn;
          $tr[$t]{d0} = $tr[$t]{d1} = -1;
        }
      } else {
        if (($tr[$tr[$t]{d1}]{u0} > 0) && ($tr[$tr[$t]{d1}]{u1} > 0)) {
          if ($tr[$tr[$t]{d1}]{u0} == $t) { # passes thru LHS
            $tr[$tr[$t]{d1}]{usave} = $tr[$tr[$t]{d1}]{u1};
            $tr[$tr[$t]{d1}]{uside} = $S_LEFT;
          } else {
            $tr[$tr[$t]{d1}]{usave} = $tr[$tr[$t]{d1}]{u0};
            $tr[$tr[$t]{d1}]{uside} = $S_RIGHT;
          }
        }
        $tr[$tr[$t]{d1}]{u0} = $t;
        $tr[$tr[$t]{d1}]{u1} = $tn;
      }

      $t = $tr[$t]{d1};

    # two trapezoids below. Find out which one is intersected by
    # this segment and proceed down that one

    } else {
      my $tmpseg = $tr[$tr[$t]{d0}]{rseg};
      my ($y0,$yt);
      my %tmppt;
      my ($tnext, $i_d0, $i_d1);

      $i_d0 = $i_d1 = $FALSE;
      if (_fp_equal($tr[$t]{lo}{y}, $s{v0}{y}, $precision)) {
        if ($tr[$t]{lo}{x} > $s{v0}{x}) {
          $i_d0 = $TRUE;
        } else {
          $i_d1 = $TRUE;
        }
      } else {
        $tmppt{y} = $y0 = $tr[$t]{lo}{y};
        $yt       = ($y0 - $s{v0}{y})/($s{v1}{y} - $s{v0}{y});
        $tmppt{x} = $s{v0}{x} + $yt * ($s{v1}{x} - $s{v0}{x});

        if (_less_than(\%tmppt, $tr[$t]{lo})) {
          $i_d0 = $TRUE;
        } else {
          $i_d1 = $TRUE;
        }
      }

      # check continuity from the top so that the lower-neighbour
      # values are properly filled for the upper trapezoid

      if (($tr[$t]{u0} > 0) && ($tr[$t]{u1} > 0)) {  # continuation of a chain from abv.
        if ($tr[$t]{usave} > 0) {  # three upper neighbours
          if ($tr[$t]{uside} == $S_LEFT) {
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = -1;
            $tr[$tn]{u1} = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0}  = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
            $tr[$tr[$tn]{u1}]{d0} = $tn;
          } else {                    # intersects in the right
            $tr[$tn]{u1} = -1;
            $tr[$tn]{u0} = $tr[$t]{u1};
            $tr[$t]{u1}  = $tr[$t]{u0};
            $tr[$t]{u0}  = $tr[$t]{usave};

            $tr[$tr[$t]{u0}]{d0}  = $t;
            $tr[$tr[$t]{u1}]{d0}  = $t;
            $tr[$tr[$tn]{u0}]{d0} = $tn;
          }

          $tr[$t]{usave} = $tr[$tn]{usave} = 0;
        } else {                      # No usave.... simple case
          $tr[$tn]{u0} = $tr[$t]{u1};
          $tr[$tn]{u1} = -1;
          $tr[$t]{u1}  = -1;
          $tr[$tr[$tn]{u0}]{d0} = $tn;
        }
      } else {                        # fresh seg. or upward cusp
        my $tmp_u = $tr[$t]{u0};
        my ($td0, $td1);
        if ((($td0 = $tr[$tmp_u]{d0}) > 0) &&
            (($td1 = $tr[$tmp_u]{d1}) > 0)) {  # upward cusp
          if (($tr[$td0]{rseg} > 0) &&
              !_is_left_of($tr[$td0]{rseg}, $s{v1})) {
            $tr[$t]{u0} = $tr[$t]{u1} = $tr[$tn]{u1} = -1;
            $tr[$tr[$tn]{u0}]{d1} = $tn;
          } else {
            $tr[$tn]{u0} = $tr[$tn]{u1} = $tr[$t]{u1} = -1;
            $tr[$tr[$t]{u0}]{d0} = $t;
          }
        } else {                               # fresh segment
          $tr[$tr[$t]{u0}]{d0} = $t;
          $tr[$tr[$t]{u0}]{d1} = $tn;
        }
      }

      if (_fp_equal($tr[$t]{lo}{y}, $tr[$tlast]{lo}{y}, $precision) &&
          _fp_equal($tr[$t]{lo}{x}, $tr[$tlast]{lo}{x}, $precision) && $tribot) {
        # this case arises only at the lowest trapezoid.. i.e.
        # tlast, if the lower endpoint of the segment is
        # already inserted in the structure

        $tr[$tr[$t]{d0}]{u0} = $t;
        $tr[$tr[$t]{d0}]{u1} = -1;
        $tr[$tr[$t]{d1}]{u0} = $tn;
        $tr[$tr[$t]{d1}]{u1} = -1;

        $tr[$tn]{d0} = $tr[$t]{d1};
        $tr[$t]{d1} = $tr[$tn]{d1} = -1;

        $tnext = $tr[$t]{d1};
      } elsif ($i_d0) {               # intersecting d0
        $tr[$tr[$t]{d0}]{u0} = $t;
        $tr[$tr[$t]{d0}]{u1} = $tn;
        $tr[$tr[$t]{d1}]{u0} = $tn;
        $tr[$tr[$t]{d1}]{u1} = -1;

        # new code to determine the bottom neighbours of the
        # newly partitioned trapezoid

        $tr[$t]{d1} = -1;

        $tnext = $tr[$t]{d0};
      } else {                        # intersecting d1
        $tr[$tr[$t]{d0}]{u0} = $t;
        $tr[$tr[$t]{d0}]{u1} = -1;
        $tr[$tr[$t]{d1}]{u0} = $t;
        $tr[$tr[$t]{d1}]{u1} = $tn;

        # new code to determine the bottom neighbours of the
        # newly partitioned trapezoid

        $tr[$tn]{d0} = $tr[$t]{d1};
        $tr[$tn]{d1} = -1;

        $tnext = $tr[$t]{d1};
      }

      $t = $tnext;
    }

    $tr[$t_sav]{rseg} = $tr[$tn_sav]{lseg}  = $segnum;
  } # end-while

  # Now combine those trapezoids which share common segments. We can
  # use the pointers to the parent to connect these together. This
  # works only because all these new trapezoids have been formed
  # due to splitting by the segment, and hence have only one parent

  $tfirstl = $tfirst;
  $tlastl  = $tlast;
  merge_trapezoids($segnum, $tfirstl, $tlastl, $S_LEFT);
  merge_trapezoids($segnum, $tfirstr, $tlastr, $S_RIGHT);

  $seg[$segnum]{is_inserted} = $TRUE;
}

# Returns true if the corresponding endpoint of the given segment is
# already inserted into the segment tree. Use the simple test of
# whether the segment which shares this endpoint is already inserted

sub _inserted {
  my ($segnum, $whichpt) = @_;
  if ($whichpt == $FIRSTPT) {
    return $seg[$seg[$segnum]{prev}]{is_inserted};
  } else {
    return $seg[$seg[$segnum]{next}]{is_inserted};
  }
}

# This is query routine which determines which trapezoid does the
# point v lie in. The return value is the trapezoid number.
#

sub _locate_endpoint {
  my ($v_ref, $vo_ref, $r) = @_;
  my %v    = %{$v_ref};
  my %vo   = %{$vo_ref};
  my %rptr = %{$qs[$r]};

  SWITCH: {
    ($rptr{nodetype} == $T_SINK) && do {
      return $rptr{trnum};
    };
    ($rptr{nodetype} == $T_Y) && do {
      if (_greater_than(\%v, $rptr{yval})) { # above
        return _locate_endpoint(\%v, \%vo, $rptr{right});
      } elsif (_equal_to(\%v, $rptr{yval})) { # the point is already
                                              # inserted.
          if (_greater_than(\%vo, $rptr{yval})) {          # above
            return _locate_endpoint(\%v, \%vo, $rptr{right});
          } else {
            return _locate_endpoint(\%v, \%vo, $rptr{left}); # below
          }
      } else {
        return _locate_endpoint(\%v, \%vo, $rptr{left});     # below
      }
    };
    ($rptr{nodetype} == $T_X) && do {
      if (_equal_to(\%v, $seg[$rptr{segnum}]{v0}) ||
          _equal_to(\%v, $seg[$rptr{segnum}]{v1})) {
        if (_fp_equal($v{y}, $vo{y}, $precision)) { # horizontal segment
          if ($vo{x} < $v{x}) {
            return _locate_endpoint(\%v, \%vo, $rptr{left});  # left
          } else {
            return _locate_endpoint(\%v, \%vo, $rptr{right}); # right
          }
        } elsif (_is_left_of($rptr{segnum}, \%vo)) {
            return _locate_endpoint(\%v, \%vo, $rptr{left});  # left
        } else {
            return _locate_endpoint(\%v, \%vo, $rptr{right}); # right
        }
      } elsif (_is_left_of($rptr{segnum}, \%v)) {
        return _locate_endpoint(\%v, \%vo, $rptr{left});  # left
      } else {
        return _locate_endpoint(\%v, \%vo, $rptr{right}); # right
      }
    };
    # default
    croak("Haggu !!!!!");
  }
}

# Thread in the segment into the existing trapezoidation. The
# limiting trapezoids are given by tfirst and tlast (which are the
# trapezoids containing the two endpoints of the segment. Merges all
# possible trapezoids which flank this segment and have been recently
# divided because of its insertion
#

sub merge_trapezoids {
  my ($segnum, $tfirst, $tlast, $side) = @_;
  my ($t, $tnext, $cond);
  my $ptnext;

  # First merge polys on the LHS
  $t = $tfirst;
  # while (($t > 0) && _greater_than_equal_to($tr[$t]{lo}, $tr[$tlast]{lo})) {
  while ($t > 0) {
    last if (! _greater_than_equal_to($tr[$t]{lo}, $tr[$tlast]{lo}));
    if ($side == $S_LEFT) {
      $cond = (((($tnext = $tr[$t]{d0}) > 0) && ($tr[$tnext]{rseg} == $segnum)) ||
               ((($tnext = $tr[$t]{d1}) > 0) && ($tr[$tnext]{rseg} == $segnum)));
    } else {
      $cond = (((($tnext = $tr[$t]{d0}) > 0) && ($tr[$tnext]{lseg} == $segnum)) ||
               ((($tnext = $tr[$t]{d1}) > 0) && ($tr[$tnext]{lseg} == $segnum)));
    }
    if ($cond) {
      if (($tr[$t]{lseg} == $tr[$tnext]{lseg}) &&
          ($tr[$t]{rseg} == $tr[$tnext]{rseg})) { # good neighbours
                                                  # merge them
        # Use the upper node as the new node i.e. t
        $ptnext = $qs[$tr[$tnext]{sink}]{parent};
        if ($qs[$ptnext]{left} == $tr[$tnext]{sink}) {
          $qs[$ptnext]{left} = $tr[$t]{sink};
        } else {
          $qs[$ptnext]{right} = $tr[$t]{sink};     # redirect parent
        }
        # Change the upper neighbours of the lower trapezoids
        if (($tr[$t]{d0} = $tr[$tnext]{d0}) > 0) {
          if ($tr[$tr[$t]{d0}]{u0} == $tnext) {
            $tr[$tr[$t]{d0}]{u0} = $t;
          } elsif ($tr[$tr[$t]{d0}]{u1} == $tnext) {
            $tr[$tr[$t]{d0}]{u1} = $t;
          }
        }
        if (($tr[$t]{d1} = $tr[$tnext]{d1}) > 0) {
          if ($tr[$tr[$t]{d1}]{u0} == $tnext) {
            $tr[$tr[$t]{d1}]{u0} = $t;
          } elsif ($tr[$tr[$t]{d1}]{u1} == $tnext) {
            $tr[$tr[$t]{d1}]{u1} = $t;
          }
        }
        $tr[$t]{lo} = {x => $tr[$tnext]{lo}{x} , y=> $tr[$tnext]{lo}{y}};
        $tr[$tnext]{state} = 2; # invalidate the lower
                                # trapezium
      } else {            #* not good neighbours
        $t = $tnext;
      }
    } else {              #* do not satisfy the outer if
        $t = $tnext;
    }
  } # end-while
}

# Retun TRUE if the vertex v is to the left of line segment no.
# segnum. Takes care of the degenerate cases when both the vertices
# have the same y--cood, etc.
#

sub _is_left_of {
  my ($segnum, $v_ref) = @_;
  my %s = %{$seg[$segnum]};
  my $area;
  my %v = %{$v_ref};

  if (_greater_than($s{v1}, $s{v0})) { # seg. going upwards
    if (_fp_equal($s{v1}{y}, $v{y}, $precision)) {
      if ($v{x} < $s{v1}{x}) {
        $area = 1;
      } else {
        $area = -1;
      }
    } elsif (_fp_equal($s{v0}{y}, $v{y}, $precision)) {
      if ($v{x} < $s{v0}{x}) {
        $area = 1;
      } else{
        $area = -1;
      }
    } else {
      $area = _Cross($s{v0}, $s{v1}, \%v);
    }
  } else {                        # v0 > v1
    if (_fp_equal($s{v1}{y}, $v{y}, $precision)) {
      if ($v{x} < $s{v1}{x}) {
        $area = 1;
      } else {
        $area = -1;
      }
    } elsif (_fp_equal($s{v0}{y}, $v{y}, $precision)) {
      if ($v{x} < $s{v0}{x}) {
        $area = 1;
      } else {
        $area = -1;
      }
    } else {
      $area = _Cross($s{v1}, $s{v0}, \%v);
    }
  }
  if ($area > 0) {
    return $TRUE;
  } else {
    return $FALSE;
  };
}

sub _Cross {
  my ($v0_ref, $v1_ref, $v2_ref) = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  my %v2 = %{$v2_ref};
  return ( ($v1{x} - $v0{x}) * ($v2{y} - $v0{y}) -
           ($v1{y} - $v0{y}) * ($v2{x} - $v0{x}) );
}

# Get log*n for given n
sub _math_logstar_n {
  my ($n) = @_;
  my $i = 0;
  for ($i = 0 ; $n >= 1 ; $i++) {
    $n = log($n)/log(2);  # log2
  }
  return ($i - 1);
}

sub _math_N {
  my ($n,$h) = @_;
  my $v = $n;
  for (my $i = 0 ; $i < $h; $i++) {
    $v = log($v)/log(2);  # log2
  }
  return (ceil($n/$v));
}

# This function returns TRUE or FALSE depending upon whether the
# vertex is inside the polygon or not. The polygon must already have
# been triangulated before this routine is called.
# This routine will always detect all the points belonging to the
# set (polygon-area - polygon-boundary). The return value for points
# on the boundary is not consistent!!!
#

sub is_point_inside_polygon {
  my @vertex = @_;
  my %v;
  my ($trnum, $rseg);

  %v = {x => $vertex[0] , y => $vertex[1]};

  $trnum = _locate_endpoint(&v, &v, 1);
  my %t = %{$tr[$trnum]};

  if ($t{state} == $ST_INVALID) {
    return $FALSE;
  }

  if (($t{lseg} <= 0) || ($t{rseg} <= 0)) {
    return $FALSE;
  }
  $rseg = $t{rseg};
  return _greater_than_equal_to($seg[$rseg]{v1}, $seg[$rseg]{v0});
}

sub _Cross_Sine {
  my ($v0_ref, $v1_ref)  = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  return ($v0{x} * $v1{y} - $v1{x} * $v0{y});
}

sub _Length {
  my ($v0_ref)  = @_;
  my %v0 = %{$v0_ref};
  return (sqrt($v0{x} * $v0{x} + $v0{y} * $v0{y}));
}

sub _Dot {
  my ($v0_ref, $v1_ref)  = @_;
  my %v0 = %{$v0_ref};
  my %v1 = %{$v1_ref};
  return ($v0{x} * $v1{x} + $v0{y} * $v1{y})
}

# Function returns TRUE if the trapezoid lies inside the polygon
sub inside_polygon {
  my ($t_ref) = @_;
  my %t = %{$t_ref};
  my $rseg = $t{rseg};
  if ($t{state} == $ST_INVALID) {
    return 0;
  }
  if (($t{lseg} <= 0) || ($t{rseg} <= 0)) {
    return 0;
  }
  if ((($t{u0} <= 0) && ($t{u1} <= 0)) ||
      (($t{d0} <= 0) && ($t{d1} <= 0)))  { # triangle
    return (_greater_than($seg[$rseg]{v1}, $seg[$rseg]{v0}));
  }
  return 0;
}

# return a new mon structure from the table
sub _newmon {
  return ++$mon_idx;
}

# return a new chain element from the table
sub _new_chain_element {
  return ++$chain_idx;
}

sub _get_angle {
  my ($vp0_ref, $vpnext_ref, $vp1_ref) = @_;
  my %vp0    = %{$vp0_ref};
  my %vpnext = %{$vpnext_ref};
  my %vp1    = %{$vp1_ref};

  my ($v0, $v1);

  $v0 = {x => $vpnext{x} - $vp0{x} , y => $vpnext{y} - $vp0{y}};
  $v1 = {x => $vp1{x}    - $vp0{x} , y => $vp1{y}    - $vp0{y}};

  if (_Cross_Sine($v0, $v1) >= 0) { # sine is positive
    return _Dot($v0, $v1)/_Length($v0)/_Length($v1);
  } else {
    return (-1 * _Dot($v0, $v1)/_Length($v0)/_Length($v1) - 2);
  }
}

# (v0, v1) is the new diagonal to be added to the polygon. Find which
# chain to use and return the positions of v0 and v1 in p and q
sub _get_vertex_positions {
  my ($v0, $v1) = @_;

  my (%vp0, %vp1);
  my ($angle, $temp);
  my ($tp, $tq);

  %vp0 = %{$vert[$v0]};
  %vp1 = %{$vert[$v1]};

  # p is identified as follows. Scan from (v0, v1) rightwards till
  # you hit the first segment starting from v0. That chain is the
  # chain of our interest

  $angle = -4.0;
  for (my $i = 0; $i < 4; $i++) {
    if ($vp0{vnext}[$i] <= 0) {
      next;
    }
    if (($temp = _get_angle($vp0{pt}, $vert[$vp0{vnext}[$i]]{pt}, $vp1{pt})) > $angle) {
      $angle = $temp;
      $tp = $i;
    }
  }

  # $ip_ref = \$tp;

  # Do similar actions for q

  $angle = -4.0;
  for (my $i = 0; $i < 4; $i++)
    {
      if ($vp1{vnext}[$i] <= 0) {
        next;
      }
      if (($temp = _get_angle($vp1{pt}, $vert[$vp1{vnext}[$i]]{pt}, $vp0{pt})) > $angle) {
        $angle = $temp;
        $tq = $i;
      }
    }

  # $iq_ref = \$tq;

  return ($tp,$tq);

}

# v0 and v1 are specified in anti-clockwise order with respect to
# the current monotone polygon mcur. Split the current polygon into
# two polygons using the diagonal (v0, v1)
#
sub _make_new_monotone_poly {
  my ($mcur, $v0, $v1) = @_;

  my ($p, $q, $ip, $iq);
  my $mnew = _newmon;
  my ($i, $j, $nf0, $nf1);

  my %vp0 = %{$vert[$v0]};
  my %vp1 = %{$vert[$v1]};

  ($ip,$iq) = _get_vertex_positions($v0, $v1);

  $p = $vp0{vpos}[$ip];
  $q = $vp1{vpos}[$iq];

  # At this stage, we have got the positions of v0 and v1 in the
  # desired chain. Now modify the linked lists

  $i = _new_chain_element;        # for the new list
  $j = _new_chain_element;

  $mchain[$i]{vnum} = $v0;
  $mchain[$j]{vnum} = $v1;

  $mchain[$i]{next} = $mchain[$p]{next};
  $mchain[$mchain[$p]{next}]{prev} = $i;
  $mchain[$i]{prev} = $j;
  $mchain[$j]{next} = $i;
  $mchain[$j]{prev} = $mchain[$q]{prev};
  $mchain[$mchain[$q]{prev}]{next} = $j;

  $mchain[$p]{next} = $q;
  $mchain[$q]{prev} = $p;

  $nf0 = $vp0{nextfree};
  $nf1 = $vp1{nextfree};

  $vert[$v0]{vnext}[$ip] = $v1;

  $vert[$v0]{vpos}[$nf0] = $i;
  $vert[$v0]{vnext}[$nf0] = $mchain[$mchain[$i]{next}]{vnum};
  $vert[$v1]{vpos}[$nf1] = $j;
  $vert[$v1]{vnext}[$nf1] = $v0;

  $vert[$v0]{nextfree}++;
  $vert[$v1]{nextfree}++;

  $mon[$mcur] = $p;
  $mon[$mnew] = $i;
  return $mnew;
}

# Main routine to get monotone polygons from the trapezoidation of
# the polygon.
#

sub _monotonate_trapezoids {
  my ($n) = @_;

  my $tr_start;

  # First locate a trapezoid which lies inside the polygon
  # and which is triangular
  my $i;
  for ($i = 1; $i < $#tr; $i++) {
    if (inside_polygon($tr[$i])) {
      last;
    }
  }
  $tr_start = $i;

  # Initialise the mon data-structure and start spanning all the
  # trapezoids within the polygon

  for (my $i = 1; $i <= $n; $i++) {
    $mchain[$i]{prev} = $seg[$i]{prev};
    $mchain[$i]{next} = $seg[$i]{next};
    $mchain[$i]{vnum} = $i;
    $vert[$i]{pt} = {x => $seg[$i]{v0}{x} , y => $seg[$i]{v0}{y}};
    $vert[$i]{vnext}[0] = $seg[$i]{next}; # next vertex
    $vert[$i]{vpos}[0] = $i;              # locn. of next vertex
    $vert[$i]{nextfree} = 1;
  }

  $chain_idx = $n;
  $mon_idx = 0;
  $mon[0] = 1;                       # position of any vertex in the first chain

  # traverse the polygon
  if ($tr[$tr_start]{u0} > 0) {
    _traverse_polygon(0, $tr_start, $tr[$tr_start]{u0}, $TR_FROM_UP);
  } elsif ($tr[$tr_start]{d0} > 0) {
    _traverse_polygon(0, $tr_start, $tr[$tr_start]{d0}, $TR_FROM_DN);
  }

  # return the number of polygons created
  return _newmon;
}

# recursively visit all the trapezoids
sub _traverse_polygon {
  my ($mcur, $trnum, $from, $dir) = @_;

  if (!$trnum) {  # patch dvdp
    return 0;
  }
  my %t = %{$tr[$trnum]};
  my ($howsplit, $mnew);
  my ($v0, $v1, $v0next, $v1next);
  my ($retval, $tmp);
  my $do_switch = $FALSE;

  if (($trnum <= 0) || $visited[$trnum]) {
    return 0;
  }

  $visited[$trnum] = $TRUE;

  # We have much more information available here.
  # rseg: goes upwards
  # lseg: goes downwards

  # Initially assume that dir = TR_FROM_DN (from the left)
  # Switch v0 and v1 if necessary afterwards

  # special cases for triangles with cusps at the opposite ends.
  # take care of this first
  if (($t{u0} <= 0) && ($t{u1} <= 0)) {
    if (($t{d0} > 0) && ($t{d1} > 0)) { # downward opening triangle
      $v0 = $tr[$t{d1}]{lseg};
      $v1 = $t{lseg};
      if ($from == $t{d1}) {
        $do_switch = $TRUE;
        $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
        _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
        _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
      } else {
        $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
        _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
        _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
      }
    } else {
      $retval = $SP_NOSPLIT;        # Just traverse all neighbours
      _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
      _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
      _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
      _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
    }
  } elsif (($t{d0} <= 0) && ($t{d1} <= 0)) {
    if (($t{u0} > 0) && ($t{u1} > 0)) { # upward opening triangle
      $v0 = $t{rseg};
      $v1 = $tr[$t{u0}]{rseg};
      if ($from == $t{u1}) {
        $do_switch = $TRUE;
        $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
        _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
      } else {
        $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
        _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
      }
    } else {
      $retval = $SP_NOSPLIT;        # Just traverse all neighbours
      _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
      _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
      _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
      _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
    }
  } elsif (($t{u0} > 0) && ($t{u1} > 0)) {
    if (($t{d0} > 0) && ($t{d1} > 0)) { # downward + upward cusps
      $v0 = $tr[$t{d1}]{lseg};
      $v1 = $tr[$t{u0}]{rseg};
      $retval = $SP_2UP_2DN;
      if ((($dir == $TR_FROM_DN) && ($t{d1} == $from)) ||
          (($dir == $TR_FROM_UP) && ($t{u1} == $from))) {
        $do_switch = $TRUE;
        $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
        _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
        _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
      } else {
        $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
        _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
        _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
      }
    } else {                      #* only downward cusp
      if (_equal_to($t{lo}, $seg[$t{lseg}]{v1})) {
        $v0 = $tr[$t{u0}]{rseg};
        $v1 = $seg[$t{lseg}]{next};

        $retval = $SP_2UP_LEFT;
        if (($dir == $TR_FROM_UP) && ($t{u0} == $from)) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
        }
      } else {
        $v0 = $t{rseg};
        $v1 = $tr[$t{u0}]{rseg};
        $retval = $SP_2UP_RIGHT;
        if (($dir == $TR_FROM_UP) && ($t{u1} == $from)) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
        }
      }
    }
  } elsif (($t{u0} > 0) || ($t{u1} > 0)) { # no downward cusp
    if (($t{d0} > 0) && ($t{d1} > 0)) { # only upward cusp
      if (_equal_to($t{hi}, $seg[$t{lseg}]{v0})) {
        $v0 = $tr[$t{d1}]{lseg};
        $v1 = $t{lseg};
        $retval = $SP_2DN_LEFT;
        if (!(($dir == $TR_FROM_DN) && ($t{d0} == $from))) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
        }
      } else {
        $v0 = $tr[$t{d1}]{lseg};
        $v1 = $seg[$t{rseg}]{next};

        $retval = $SP_2DN_RIGHT;
        if (($dir == $TR_FROM_DN) && ($t{d1} == $from)) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
        }
      }
    } else { # no cusp
      if (_equal_to($t{hi}, $seg[$t{lseg}]{v0}) &&
          _equal_to($t{lo}, $seg[$t{rseg}]{v0})) {
        $v0 = $t{rseg};
        $v1 = $t{lseg};
        $retval = $SP_SIMPLE_LRDN;
        if ($dir == $TR_FROM_UP) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
        }
      } elsif (_equal_to($t{hi}, $seg[$t{rseg}]{v1}) &&
               _equal_to($t{lo}, $seg[$t{lseg}]{v1})) {
        $v0 = $seg[$t{rseg}]{next};
        $v1 = $seg[$t{lseg}]{next};

        $retval = $SP_SIMPLE_LRUP;
        if ($dir == $TR_FROM_UP) {
          $do_switch = $TRUE;
          $mnew = _make_new_monotone_poly($mcur, $v1, $v0);
          _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{d0}, $trnum, $TR_FROM_UP);
        } else {
          $mnew = _make_new_monotone_poly($mcur, $v0, $v1);
          _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
          _traverse_polygon($mnew, $t{u0}, $trnum, $TR_FROM_DN);
          _traverse_polygon($mnew, $t{u1}, $trnum, $TR_FROM_DN);
        }
      } else { # no split possible
        $retval = $SP_NOSPLIT;
        _traverse_polygon($mcur, $t{u0}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mcur, $t{d0}, $trnum, $TR_FROM_UP);
        _traverse_polygon($mcur, $t{u1}, $trnum, $TR_FROM_DN);
        _traverse_polygon($mcur, $t{d1}, $trnum, $TR_FROM_UP);
      }
    }
  }

  return $retval;
}

# For each monotone polygon, find the ymax and ymin (to determine the
# two y-monotone chains) and pass on this monotone polygon for greedy
# triangulation.
# Take care not to triangulate duplicate monotone polygons

sub _triangulate_monotone_polygons {
  my ($nvert, $nmonpoly) = @_;

  my ($ymax, $ymin);
  my ($p, $vfirst, $posmax, $posmin, $v);
  my ($vcount, $processed);

  $op_idx = 0;
  for (my $i = 0; $i < $nmonpoly; $i++) {
    $vcount = 1;
    $processed = $FALSE;
    $vfirst = $mchain[$mon[$i]]{vnum};
    $ymax = {x => $vert[$vfirst]{pt}{x} , y => $vert[$vfirst]{pt}{y}};
    $ymin = {x => $vert[$vfirst]{pt}{x} , y => $vert[$vfirst]{pt}{y}};
    $posmax = $posmin = $mon[$i];
    $mchain[$mon[$i]]{marked} = $TRUE;
    $p = $mchain[$mon[$i]]{next};
    while (($v = $mchain[$p]{vnum}) != $vfirst) {
      if ($mchain[$p]{marked}) {
        $processed = $TRUE;
        last;                # break from while
      } else {
        $mchain[$p]{marked} = $TRUE;
      }

      if (_greater_than($vert[$v]{pt}, $ymax)) {
        $ymax = {x => $vert[$v]{pt}{x} , y => $vert[$v]{pt}{y}};
        $posmax = $p;
      }
      if (_less_than($vert[$v]{pt}, $ymin)) {
        $ymin = {x => $vert[$v]{pt}{x} , y => $vert[$v]{pt}{y}};
        $posmin = $p;
      }
      $p = $mchain[$p]{next};
      $vcount++;
    }

    if ($processed) {              # Go to next polygon
      next;
    }

    if ($vcount == 3) {            # already a triangle
      $op[$op_idx][0] = $mchain[$p]{vnum};
      $op[$op_idx][1] = $mchain[$mchain[$p]{next}]{vnum};
      $op[$op_idx][2] = $mchain[$mchain[$p]{prev}]{vnum};
      $op_idx++;
    } else {                      # triangulate the polygon
      $v = $mchain[$mchain[$posmax]{next}]{vnum};
      if (_equal_to($vert[$v]{pt}, $ymin)) {  # LHS is a single line
        _triangulate_single_polygon($nvert, $posmax, $TRI_LHS);
      } else {
        _triangulate_single_polygon($nvert, $posmax, $TRI_RHS);
      }
    }
  }

  return $op_idx;
}

# A greedy corner-cutting algorithm to triangulate a y-monotone
# polygon in O(n) time.
# Joseph O-Rourke, Computational Geometry in C.
#
sub _triangulate_single_polygon {
  my ($nvert, $posmax, $side) = @_;

  my $v;
  my @rc;
  my $ri = 0;        # reflex chain
  my ($endv, $tmp, $vpos);

  if ($side == $TRI_RHS) {   # RHS segment is a single segment
    $rc[0] = $mchain[$posmax]{vnum};
    $tmp   = $mchain[$posmax]{next};
    $rc[1] = $mchain[$tmp]{vnum};
    $ri = 1;

    $vpos = $mchain[$tmp]{next};
    $v = $mchain[$vpos]{vnum};

    if (($endv = $mchain[$mchain[$posmax]{prev}]{vnum}) == 0) {
      $endv = $nvert;
    }
  } else {                              # LHS is a single segment
    $tmp = $mchain[$posmax]{next};
    $rc[0] = $mchain[$tmp]{vnum};
    $tmp = $mchain[$tmp]{next};
    $rc[1] = $mchain[$tmp]{vnum};
    $ri = 1;

    $vpos = $mchain[$tmp]{next};
    $v = $mchain[$vpos]{vnum};

    $endv = $mchain[$posmax]{vnum};
  }

  while (($v != $endv) || ($ri > 1)) {
    if ($ri > 0) {              # reflex chain is non-empty
      if (_Cross($vert[$v]{pt}, $vert[$rc[$ri - 1]]{pt}, $vert[$rc[$ri]]{pt}) > 0) {
        # convex corner: cut if off
        $op[$op_idx][0] = $rc[$ri - 1];
        $op[$op_idx][1] = $rc[$ri];
        $op[$op_idx][2] = $v;
        $op_idx++;
        $ri--;
      } else {     # non-convex
                   # add v to the chain
        $ri++;
        $rc[$ri] = $v;
        $vpos = $mchain[$vpos]{next};
        $v = $mchain[$vpos]{vnum};
      }
    } else {       # reflex-chain empty: add v to the
                   # reflex chain and advance it
      $rc[++$ri] = $v;
      $vpos = $mchain[$vpos]{next};
      $v = $mchain[$vpos]{vnum};
    }
  } # end-while

  # reached the bottom vertex. Add in the triangle formed
  $op[$op_idx][0] = $rc[$ri - 1];
  $op[$op_idx][1] = $rc[$ri];
  $op[$op_idx][2] = $v;
  $op_idx++;
  $ri--;

}

1;
