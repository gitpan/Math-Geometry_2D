#
# Geometry_2D.pm, version 0.99 July 2002
#
# Copyright (c) 2002 Danny Van de Pol - Alcatel Telecom Belgium
# danny.vandepol@alcatel.be
#
# Free usage under the same Perl Licence condition.
#

package Math::Geometry_2D;

use vars qw($VERSION $precision);
$VERSION   = '1.00';
$precision = 7;

require Exporter;
@ISA='Exporter';
@EXPORT = qw/SegmentLength Determinant DotProduct CrossProduct
             TriangleArea Colinear
             SegmentIntersection LineIntersection
             Perpendicular PerpendicularFoot
             DistanceToLine DistanceToSegment
             Gpc2Polygons GpcClip
             CircleToPoly ArcToPoly
            /;

use strict;
use GPC;
use Carp;

require "Triangulate.pl";

=pod

=head1 NAME

Geometry_2D - Module with 2D geometry functions

=head1 SYNOPSIS

 use Math::Geometry_2D;
 $polygon = Math::Geometry_2D=>new; creates a new polygon object;
 $contour = Math::Geometry_2D=>new; creates a new contour object;

=head4 Formats

A point is a refenece to an array holding the x and y coordinates of the point.

 $point = [$x_coord,$y_coord];

A polygon is a reference to an (ordered) array of points.  The first point is the
begin and end point of the polygon. The points can be given in any direction
(clockwise or counter clockwise).

 $points = [[$x1,$y1],[$x2,$y2], ... ];
 $polygon->points($points);                    # assign points to polygon object
 $points1 = [[$x1,$y1],[$x2,$y2], ... ];
 $points2 = [[ax1,by1],[ax2,by2], ... ];
 $contour->polygons([$points1,$points2, ...]); # assign polgyons to contour object

A contour is a reference to an array of polygons.  By convention, the first polygon
is the outer shape, all other polygons represent holes in the outer shape.  The outer
shape must enclose all holes !
Using this convention, the points can be given in any direction, however, keep
in mind that some functions (e.g. triangulation) require that the outer polygons
are entered in counter clockwise order and the inner polygons (holes) in clock
wise order.  The points, polygons, add_polygons methods will automatically set the
right order of points.

 $contour = [$poly1,$poly2], ... ];

=head1 METHODS

The available methods are:

=head4 $polygon->points(arg);

 Returns the polygon points if no argument is entered
 If the argument is a refence to a points array, sets the points for a polygon object

=head4 $contour->polygons(arg);

 Returns the contour polygons if no argument is entered
 If the argument is a refence to a polygons array, sets the polygons for a contour object

=head4 $contour->num_polygons;

 Returns the total number of polygons in the contour.

=head4 $contour->add_polygons(arg);

 Adds a list of polygons to a contour object (if the contour object doesn't have any
 polygons yet, the very first polygon reference from the list is used as the outer
 shape.  Returns the total number of polygons in the contour.

=head4 $contour->get_polygons(arg_1,arg_2, ... );

 Returns a list of polygons where each element of the list corresponds to the polygon
 at index arg_x - starting at 0.  If the index arg_x is out of range, the corresponding
 value in the result list wil be undefined.  If no argument is entered, a full list of
 all polygons will be returned. Please note that this method returns a list rather
 then a reference.

=head4 $polygon->cleanup;

 Remove colinear points from the polygon/contour.

=head4 $polygon->isconvex;

 Returns true if the polygon/contour is convex (a contour is considered to be convex if
 the outer shape is convex)

=head4 $polygon->issimple;

 Returns true if the polygon/contour is simple (a contour is considered to be simple if
 all it's polygons are simple)

=head4 $polygon->perimeter;

 Returns the perimeter of the polygon/contour (the perimeter of a contour is the perimeter
 of the outer shape)

=head4 $polygon->area;

 Returns the signed area of the polygon/contour (positive if the points are in counter
 clockwise order) (the area of a contour is the area of the outer shape minus the tota
 area the holes)

=head4 $polygon->centroid;

 Returns the centroid of the polygon/contour

=head4 $polygon->isinside($point);

 Returns true if point is inside the polygon/contour (a point is inside a contour if
 it is inside the outer polygon and not inside a hole)

=head4 $polygon->rotate($angle,$center);

 Returns polygon/contour rotated $angle (in radians) around $center

=head4 $polygon->move($dx,$dy);

 Returns polygon/contour moved $dx in x direction and $dy in y direction

=head4 $polygon->mirrorx($center);

 Returns polygon/contour mirrored in x direction
 with (vertical) axis of reflection through point $center

=head4 $polygon->mirrory($center);

 Returns polygon/contour mirrored in y direction
 with (horizontal) axis of reflection through point $center

=head4 $polygon->mirror($axos);

 Returns polygon mirrored/contour along axis $axis (= array with 2 points defining
 axis of reflection))

=head4 $polygon->scale($csale,$center);

 Returns polygon/contour scaled by a factor $scale, center of scaling is $scale

=head4 $polygon->bbox;

 Returns the polygon's/contour's bounding box

=head4 $polygon->minrectangle;

 Returns the polygon's/contour's minimal (area) enclosing rectangle

=head4 $polygon->convexhull;

 Returns a polygon representing the convex hull of the polygon/contour

=head4 $polygon->convexhull2;

 Returns a polygon representing the convex hull of an arbitrary set of points
 (works also on a contour, however a contour is a set of polygons and polygons
  are ordered sets of points so the method above will be faster)

=head4 $polygon->triangulate;

 Triangulates a polygon/contour based on Raimund Seidel's algorithm:
 'A simple and fast incremental randomized algorithm for computing trapezoidal
 decompositions and for triangulating polygons'
 Returns a reference to a list of triangles

=head4 $polygon->convert2gpc;

 Converts a polygon/contour to a gpc structure and returns the resulting gpc structure

=head1 EXPORTS

=head4 SegmentLength[$p1,$p2];

 Returns the length of the segment (vector) p1p2

=head4 Determinant(x1,y1,x2,y2);

 Returns | x1 y1 | which is x1*y2 - y1*x2
         | x2 y2 |

=head4 DotProduct($p1,$p2,$p3,$p4);

 Returns the vector dot product of vectors p1p2 and p3p4
 or the dot product of p1p2 and p2p3 if $p4 is ommited from the argument list

=head4 CrossProduct($p1,$p2,$p3);

 Returns the vector cross product of vectors p1p2 and p1p3

=head4 TriangleArea($p1,$p2,$p3);

 Returns the signed area of the triangle p1p2p3

=head4 Colinear($p1,$p2,$p3);

 Returns true if p1,p2 and p3 are colinear

=head4 SegmentIntersection($p1,$p2,$p3,$p4);

 Returns false if segments don't intersect
 Returns the intersection point of segments p1p2 and p3p4

=head4 LineIntersection($p1,$p2,$p3,$p4);

 Returns false if lines don't intersect (parallel lines)
 Returns the intersection point of lines p1p2 and p3p4

=head4 Perpendicular($p1,$p2,$p3,$p4);

 Returns true if lines (segments) p1p2 and p3p4 are perpendicular

=head4 PerpendicularFoot($p1,$p2,$p3);

 Returns the perpendicular foot of p3 on line p1p2

=head4 DistanceToLine($p1,$p2,$p3);

 Returns the perpendicular dostance of p3 to line p1p2

=head4 DistanceToSegment($p1,$p2,$p3);

 Returns the perpendicular distance of p3 to segment p1p2

=head4 Gpc2Polygons($gpc_contour);

 Comverts a gpc contour structure to an array of contours and returns the array

=head4 GpcClip($operation,$gpc_contour_1,$gpc_contour_2);

 $operation is DIFFERENCE, INTERSECTION, XOR or UNION
 $gpc_polygon_1 is the source polygon
 $gpc_polygon_2 is the clip polygon
 Returns a gpc polygon structure which is the result of the gpc clipping operation

=head4 CircleToPoly($i,$p1,$p2,$p3);

 Converts the circle through points p1p2p3 to a polygon with i segments

=head4 CircleToPoly($i,$center,$p1);

 Converts the circle with center through points p1 to a polygon with i segments

=head4 CircleToPoly($i,$center,$radius);

 Converts the circle with center and radius to a polygon with i segments

=head4 ArcleToPoly($i,$p1,$p2,$p3);

 Converts the arc with begin point p1, intermediate point p2 and end point p3
 to a (non-closed !) polygon with i segments

=head4 ArcleToPoly($i,$center,$p1,$p2,$direction);

 Converts the arc with center, begin point p1 and end point p2 to a
 (non-closed !) polygon with i segments.  If direction is 0, the arc
 is traversed counter clockwise from p1 to p2, clockwise if direction is 1

=cut


require 5.005;

my $delta = 10 ** (-$precision);

################################################################################
#
# calculate length of a line segment
#
# args : reference to array with 2 points defining line segment
#
sub SegmentLength {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 2) {
    carp("Need 2 points for a segment length calculation");
    return;
  }
  my @a = @{$points[0]};
  my @b = @{$points[1]};
  my $length = sqrt(DotProduct([$points[0],$points[1],$points[0],$points[1]]));
  return $length;
}
################################################################################
#  
#  The determinant for the matrix  | x1 y1 |
#                                  | x2 y2 |
#
# args : x1,y1,x2,y2
#
sub Determinant {
  my ($x1,$y1,$x2,$y2) = @_;
  return ($x1*$y2 - $x2*$y1);
}
################################################################################
#
# vector dot product
# calculates dotproduct vectors p1p2 and p3p4
# The dot product of a and b  is written as a.b and is
# defined by a.b = |a|*|b|*cos q 
#
# args : reference to an array with 4 points p1,p2,p3,p4 defining 2 vectors
#        a = vector p1p2 and b = vector p3p4
#        or
#        reference to an array with 3 points p1,p2,p3 defining 2 vectors
#        a = vector p1p2 and b = vector p1p3
#
sub DotProduct {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  my (@p1,@p2,@p3,@p4);
  if (@points == 4) {
    @p1 = @{$points[0]};
    @p2 = @{$points[1]};
    @p3 = @{$points[2]};
    @p4 = @{$points[3]};
  } elsif (@points == 3) {
    @p1 = @{$points[0]};
    @p2 = @{$points[1]};
    @p3 = @{$points[0]};
    @p4 = @{$points[2]};
  } else {
    carp("Need 3 or 4 points for a dot product");
    return;
  }
  return ($p2[0]-$p1[0])*($p4[0]-$p3[0]) + ($p2[1]-$p1[1])*($p4[1]-$p3[1]);
}
################################################################################
#
# returns vector cross product of vectors p1p2 and p1p3
# using Cramer's rule
#
# args : reference to an array with 3 points p1,p2 and p3
#
sub CrossProduct {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 3) {
    carp("Need 3 points for a cross product");
    return;
  }
  my @p1 = @{$points[0]};
  my @p2 = @{$points[1]};
  my @p3 = @{$points[2]};
  my $det_p2p3 = &Determinant($p2[0], $p2[1], $p3[0], $p3[1]);
  my $det_p1p3 = &Determinant($p1[0], $p1[1], $p3[0], $p3[1]);
  my $det_p1p2 = &Determinant($p1[0], $p1[1], $p2[0], $p2[1]);
  return ($det_p2p3-$det_p1p3+$det_p1p2);
}
################################################################################
#
#  The Cramer's Rule for area of a triangle is
#                                  | x1 y1 1 |
#                        A = 1/2 * | x2 y2 1 |
#                                  | x3 y3 1 |
# Which is 'half of the cross product of vectors ab and ac.
# The cross product of the vectors ab and ac is a vector perpendicular to the
# plane defined by ab and bc with a magnitude equal to the area of the
# parallelogram defined by a, b, c and ab + bc (vector sum)
# Don't forget that:  (ab x ac) = - (ac x ab)  (x = cross product)
# Which just means that if you reverse the vectors in the cross product,
# the resulting vector points in the opposite direction
# The direction of the resulting vector can be found with the "right hand rule"
# This can be used to determine the order of points a, b and c:
# clockwise or counter clockwise
#
# args : reference to an array with 3 points p1.p2,p3
#
sub TriangleArea {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 3) {  # need 3 points for a triangle ...
    carp("A triangle should have exactly 3 points");
    return;
  }
  return CrossProduct($pointsref)/2;
}
################################################################################
# 
# Check if 3 points are colinear
# Points are colinear if triangle area is 0
# Triangle area is crossproduct/2 so we can check the crossproduct instead
#
# args : reference to an array with 3 points p1.p2,p3
#
sub Colinear {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 3) {
    carp("Colinear only checks colinearity for 3 points");
    return;
  }
  # check the area of the triangle to find
  return (abs(CrossProduct($pointsref)) < $delta);
}
################################################################################
#
# calculate intersection point of 2 line segments
# returns false if segments don't intersect
# The theory:
#
#  Parametric representation of a line
#    if p1 (x1,y1) and p2 (x2,y2) are 2 points on a line and
#       P1 is the vector from (0,0) to (x1,y1)
#       P2 is the vector from (0,0) to (x2,y2)
#    then the parametric representation of the line is P = P1 + k (P2 - P1)
#    where k is an arbitrary scalar constant.
#    for a point on the line segement (p1,p2)  value of k is between 0 and 1
#
#  for the 2 line segements we get
#      Pa = P1 + k (P2 - P1)
#      Pb = P3 + l (P4 - P3)
#
#  For the intersection point Pa = Pb so we get the following equations
#      x1 + k (x2 - x1) = x3 + l (x4 - x3)
#      y1 + k (y2 - y1) = y3 + l (y4 - y3)
#  Which using Cramer's Rule results in
#          (x4 - x3)(y1 - y3) - (y4 - x3)(x1 - x3)
#      k = ---------------------------------------
#          (y4 - y3)(x2 - x1) - (x4 - x3)(y2 - y1)
#   and
#          (x2 - x1)(y1 - y3) - (y2 - y1)(x1 - x3)
#      l = ---------------------------------------
#          (y4 - y3)(x2 - x1) - (x4 - x3)(y2 - y1)
#
#  Note that the denominators are equal.  If the denominator is 9,
#  the lines are parallel.  Intersection is detected by checking if
#  both k and l are between 0 and 1.
#
#  The intersection point p5 (x5,y5) is:
#     x5 = x1 + k (x2 - x1)
#     y5 = y1 + k (y2 - y1)
#
# 'Touching' segments are considered as not intersecting
#
# args : reference to an array with 4 points p1,p2,p3,p4
#
sub SegmentIntersection {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 4) {
    carp("SegmentIntersection needs 4 points");
    return;
  }
  my @p1 = @{$points[0]}; # p1,p2 = segment 1
  my @p2 = @{$points[1]};
  my @p3 = @{$points[2]}; # p3,p4 = segment 2
  my @p4 = @{$points[3]};
  my @p5;
  my $n1 = Determinant(($p3[0]-$p1[0]),($p3[0]-$p4[0]),($p3[1]-$p1[1]),($p3[1]-$p4[1]));
  my $n2 = Determinant(($p2[0]-$p1[0]),($p3[0]-$p1[0]),($p2[1]-$p1[1]),($p3[1]-$p1[1]));
  my $d  = Determinant(($p2[0]-$p1[0]),($p3[0]-$p4[0]),($p2[1]-$p1[1]),($p3[1]-$p4[1]));
  if ($d == 0) {
    return 0; # parallel
  }
  if (!(($n1/$d < 1) && ($n2/$d < 1) &&
        ($n1/$d > 0) && ($n2/$d > 0))) {
    return 0;
  }
  $p5[0] = $p1[0] + $n1/$d * ($p2[0] - $p1[0]);
  $p5[1] = $p1[1] + $n1/$d * ($p2[1] - $p1[1]);
  return \@p5; # intersection point
}
################################################################################
#
# Intersection point of 2 lines - (almost) identical as for Segments
# each line is defined by 2 points
# 
# args : reference to an array with 4 points p1,p2,p3,p4
#
sub LineIntersection {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points < 4) {
    carp("LineIntersection needs 4 points");
    return;
  }
  my @p1 = @{$points[0]}; # p1,p2 = line 1
  my @p2 = @{$points[1]};
  my @p3 = @{$points[2]}; # p3,p4 = line 2
  my @p4 = @{$points[3]};
  my @p5;
  my $n1 = Determinant(($p3[0]-$p1[0]),($p3[0]-$p4[0]),($p3[1]-$p1[1]),($p3[1]-$p4[1]));
  my $n2 = Determinant(($p2[0]-$p1[0]),($p3[0]-$p1[0]),($p2[1]-$p1[1]),($p3[1]-$p1[1]));
  my $d  = Determinant(($p2[0]-$p1[0]),($p3[0]-$p4[0]),($p2[1]-$p1[1]),($p3[1]-$p4[1]));
  if ($d == 0) {
    return 0; # parallel
  }
  $p5[0] = $p1[0] + $n1/$d * ($p2[0] - $p1[0]);
  $p5[1] = $p1[1] + $n1/$d * ($p2[1] - $p1[1]);
  return \@p5; # intersection point
}
################################################################################
#
# returns true if 2 lines (segments) are perpendicular
# Lines are perpendicular if dot product is 0
# 
# args : reference to an array with 4 points p1,p2,p3,p4
#        p1p2 = line 1
#        p3p4 = line 2
#
sub Perpendicular {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 4) {
    carp("Perpendicular needs 4 points defining 2 lines or segments");
    return;
  }
  return (abs(DotProduct([$points[0],$points[1],$points[2],$points[3]])) < $delta);
}
################################################################################
#
# Calculates the 'perpendicular foot' of a point on a line
#
# args: reference to array with 3 points p1,p2,p3
#       p1p2 = line
#       p3   = point for which perpendicular foot is to be calculated
#
sub PerpendicularFoot {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points != 3) {
    carp("PerpendicularFoot needs 3 points defining a line and a point");
    return;
  }
  my @p1 = @{$points[0]}; # p1,p2 = line
  my @p2 = @{$points[1]};
  my @p3 = @{$points[2]}; # p3 point
  # vector penpenidular to line
  my @v;
  $v[0] =     $p2[1] - $p1[1];  # y2-y1
  $v[1] =  - ($p2[0] - $p1[0]); # -(x2-x1);
  # p4 = v + p3 is a second point of the line perpendicular to p1p2 going through p3
  my @p4;
  $p4[0] =  $p3[0] + $v[0];
  $p4[1] =  $p3[1] + $v[1];
  return LineIntersection([\@p1,\@p2,\@p3,\@p4]);
}
################################################################################
#
# Calculate distance from point p to line segment p1p2
#
# args: reference to array with 3 points: p1,p2,p3
#       p1p2 = segment
#       p3   = point for which distance is to be calculated
# returns distance from p3 to line segment p1p2
#         which is the smallest value from:
#            distance p3p1
#            distance p3p2
#            perpendicular distance from p3 to line p1p2
#
sub DistanceToSegment {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points < 3) {
    carp("DistanceToSegment needs 3 points defining a segment and a point");
    return;
  }
  # the perpendicular distance is the height of the parallelogram defined
  # by the 3 points devided by the base
  # Note the this is a signed value so it can be used to check at which
  # side the point is located
  # we use dot products to find out where point is located1G/dotpro
  my $d1 = DotProduct([$points[0],$points[1],$points[0],$points[2]]);
  my $d2 = DotProduct([$points[0],$points[1],$points[0],$points[1]]);
  my $dp = CrossProduct([$points[2],$points[0],$points[1]]) / sqrt $d2;
  if ($d1 <= 0) {
    return SegmentLength([$points[2],$points[0]]);
  } elsif ($d2 <= $d1) {
    return SegmentLength([$points[2],$points[1]]);
  } else {
    return $dp;
  }
}
################################################################################
#
# Calculate distance from point p to line p1p2
#
# args: reference to array with 3 points: p1,p2,p3
#       p1p2 = line
#       p3   = point for which distance is to be calculated
# returns 2 numbers
#   - perpendicular distance from p3 to line p1p2
#   - distance from p3 to line segment p1p2
#     which is the smallest value from:
#            distance p3p1
#            distance p3p2
#
sub DistanceToLine {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points < 3) {
    carp("DistanceToLine needs 3 points defining a line and a point");
    return;
  }
  # the perpendicular distance is the height of the parallelogram defined
  # by the 3 points devided by the base
  # Note the this is a signed value so it can be used to check at which
  # side the point is located
  # we use dot products to find out where point is located1G/dotpro
  my $d  = DotProduct([$points[0],$points[1],$points[0],$points[1]]);
  my $dp = CrossProduct([$points[2],$points[0],$points[1]]) / sqrt $d;
  return $dp;
}
################################################################################
#
# Initializer
#
sub new {
  my $invocant = shift;
  my $class = ref($invocant) || $invocant;
  my $self = { @_ };
  bless($self,$class);
  return $self;
}
################################################################################
#
# args: reference to polygon object
#
sub points {
  my Math::Geometry_2D $self = shift;
  if (@_) {
    if ($self->get_polygons) {
      carp("Object is a contour - can't add points");
      return;
    } else {
      # delete existing info
      $self->{points} = ();
      my $pointsref = shift;
      # normalize (a single polygon has only an outer shape make
      # -> the points order clockwise)
      if (PolygonArea($pointsref) > 0) {
        $self->{points} = $pointsref;
      } else {
        $self->{points} = [reverse @{$pointsref}];
      }
    }
  }
  return $self->{points};
}
################################################################################
#
# args: reference to polygon object
#
sub polygons {
  my Math::Geometry_2D $self = shift;
  if (@_) {
    if ($self->points) {
      carp("Object is a polygon - can't add polygons");
      return;
    } else {
      # delete existing info
      $self->{polygons} = ();
      my $polygons = shift;
      my @polygonrefs = @{$polygons};
      $self->add_polygons(@polygonrefs);
    }
  }
  return $self->{polygons};
}
################################################################################
#
# args: none
# returns the number of polygons in the contour
#
sub num_polygons {
  my Math::Geometry_2D $self = shift;
  my $polygons = $self->{polygons};
  return 0 if (! $polygons);
  return scalar @{$polygons};
}
################################################################################
#
# args: list of references to polygons
# returns the number of polygons in the contour
#
sub add_polygons {
  my Math::Geometry_2D $self = shift;
  return if (! @_); # nothing to add
  # can't add polygons to a polygon object
  if ($self->points) {
    carp("Object is a polygon - can't add polygons");
    return;
  }
  # first polygon is outer polygon
  if (! $self->num_polygons) {
    my $outer = shift;
    # counter clockwise for outer polygon
    if (PolygonArea($outer) < 0) {
      push @{$self->{polygons}}, [reverse @{$outer}];
    } else {
      push @{$self->{polygons}}, $outer;
    }
  }
  # inner polygon(s)
  while (@_) {
    # clockwise for inner polygon
    my $inner = shift;
    if (PolygonArea($inner) > 0) {
      push @{$self->{polygons}}, [reverse @{$inner}];
    } else {
      push @{$self->{polygons}}, $inner;
    }
  }
  return scalar @{$self->{polygons}};
}
################################################################################
#
# args: list of indices
# returns list of polygons indicated by indices
#         (list value at position n is undefined if the index at position
#          n is out of range)
#         list of all polygons indicated by indices
#
sub get_polygons {
  my Math::Geometry_2D $self = shift;
  my @result;
  my $polygons = $self->{polygons};
  return if (! $polygons);
  my $i = 0;
  if (@_) {
    while (@_) {
      my $index = int shift;
      if ($index >= 0 && $index < num_polygons($self)) {
        $result[$i] = ${$polygons}[$index];
      } else {
        $result[$i] = undef;
      }
      $i++;
    }
    return @result;
  } else {
    return @{$polygons};
  }
}
################################################################################
# cleanup polygon = remove colinear points
#
# args: reference to polygon or contour object
#
sub cleanup {
  my ($self) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {    # polygon object
    my @points = @$pointsref;
    for (my $i=0 ; $i< @points && @points > 2 ;$i++) {
      if (Colinear([$points[$i-2],$points[$i-1],$points[$i]])) {
        splice @points,$i-1,1;
        $i--;
      }
    }
    # replace polygon points
    $self->points([@points]);
    return [@points];
  } else {             # contour object
    my @polygonrefs = $self->get_polygons;
    for (my $j = 0; $j < @polygonrefs; $j++) {
      $pointsref = $polygonrefs[$j];
      my @points = @$pointsref;
      for (my $i=0 ; $i< @points && @points > 2 ;$i++) {
        if (Colinear([$points[$i-2],$points[$i-1],$points[$i]])) {
          splice @points,$i-1,1;
          $i--;
        }
      }
      $polygonrefs[$j] = [@points];
    }
    $self->polygons([@polygonrefs]);
    return [@polygonrefs];
  }
}
################################################################################
#
# Ah - more vector algebra
# We consider every set of 3 subsequent points p1,p2,p3 on the polygon and calculate
# the vector product of the vectors  p1p2 and p1p3.  All these products should
# have the same sign.  If the sign changes, the polygon is not convex
#
# make sure to remove colinear points first before calling perimeter
# (I prefer not to include the call to cleanup)
#
# args: reference to polygon or contour object
#       (for a contour we only check the outer shape)
#
sub isconvex {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  return 1 if (@points < 5); # every poly with a less then 5 points is convex
  my $prev = 0;
  for (my $i = 0 ; $i < @points ; $i++) {
    my $tmp = CrossProduct([$points[$i-2],$points[$i-1],$points[$i]]);
    # check if sign is different from pervious one(s)
    if ( ($prev < 0 && $tmp > 0) ||
         ($prev > 0 && $tmp < 0) ) {
      return 0;
    }
    $prev = $tmp;
  }
  return 1;
}
################################################################################
#
# Brute force attack:
# just check intersection for every segment versus every other segment
# so for a polygon with n ponts this will take n**2 intersection calculations
# I added a few simple improvements: to boost speed:
#   - don't check adjacant segments
#   - don't check against 'previous' segments (if we checked segment x versus y,
#     we don't need to check y versus x anymore)
# Results in (n-2)*(n-1)/2 - 1 checks  which is close to n**2/2 for large n
#
# make sure to remove colinear points first before calling perimeter
# (I prefer not to include the call to cleanup)
#
# args: reference to polygon or contour object
#       (a contour is considered to be simple if all it's shapes are simple)
#
sub IsSimplePolygon {
  my ($pointsref) = @_;
  my @points = @$pointsref;
  return 1 if (@points < 4); # triangles are simple polygons ...
  for (my $i = 0 ; $i < @points-2 ; $i++) {
    # check versus all next non-adjacant edges
    for (my $j = $i+2 ; $j < @points ; $j++) {
      # don't check first versus last segment (adjacant)
      next if ($i == 0 && $j == @points-1);
      if (SegmentIntersection([$points[$i-1],$points[$i],$points[$j-1],$points[$j]])) {
        return 0;
      }
    }
  }
  return 1;
}
################################################################################
#
# Check if polyogn or contour is simple
sub issimple {
  my ($self) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    return IsSimplePolygon($pointsref);
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      return 0 if (! IsSimplePolygon($_));
    }
    return 1;
  }
}
################################################################################
# makes only sense for simple polygons
# make sure to remove colinear points first before calling perimeter
# (I prefer not to include the call to colinear)
#
# args: reference to polygon or contour object
# returns the perimeter of the polygon or the perimeter of the outer shape of
# the contour
#
sub perimeter {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  my $perimeter = 0;
  if ($pointsref) {
    my @points = @$pointsref;
    if (@points < 3) { # no perimeter for lines and points
      carp("Can't calculate perimeter: polygon should have at least 3 points");
      return;
    }
    for (my $index=0;$index < @points; $index++) {
      $perimeter += SegmentLength([$points[$index-1],$points[$index]]);
    }
  }
  return $perimeter;
}
################################################################################
# makes only sense for simple polygons
# make sure to remove colinear points first before calling area
# returns a signed value, can be used to find out whether
# the order of points is clockwise or counter clockwise
# (I prefer not to include the call to colinear)
#
# args: reference to an array of points
#
sub PolygonArea {
  my $pointsref = $_[0];
  my @points = @$pointsref;
  if (@points < 3) { # no area for lines and points
    carp("Can't calculate area: polygon should have at least 3 points");
    return;
  }
  my $area = 0;
  while(@points >= 3){
    $area+=TriangleArea([$points[0],$points[1],$points[2]]);
    splice @points,1,1;
  }
  return $area;
}
################################################################################
# Calculates the area of a polygon or a contour
# Makes only sense for simple polygons
# Returns a signed value so it can be used to find out whether
# the order of points in a polygon is clockwise or counter
# clockwise.
#
# args: reference to polygon or contour object
#
sub area {
  my ($self) = @_;
  my $pointsref = $self->points;
  my $area = 0;
  if ($pointsref) {
    $area = PolygonArea($pointsref);
  } else {
    my @polygonrefs = $self->get_polygons;
    foreach (@polygonrefs) {
      $area += PolygonArea($_);
    }
  }
  return $area;
}
################################################################################
#
# calculate the centroid of a polygon or contour
# (a.k.a. the center of mass a.k.a. the center of gravity)
#
# The centroid is calculated as the weighted sum of the centroids
# of a partition of the polygon into triangles. The centroid of a
# triangle is simply the average of its three vertices, i.e., it
# has coordinates (x1 + x2 + x3)/3 and (y1 + y2 + y3)/3. 
# In fact, the triangulation need not be a partition, but rather
# can use positively and negatively oriented triangles (with positive
# and negative areas), as is used when computing the area of a polygon
#
# makes only sense for simple polygons
# make sure to remove colinear points first before calling centroid
# (I prefer not to include the call to cleanup)
#
# args: reference to polygon object
#
sub centroid {
  my ($self) = @_;
  my $trianglesref = $self->triangulate;
  my @triangles = @{$trianglesref};

  if (! @triangles) { # no result from triangulation
    carp("Nothing to calculate centroid for");
    return;
  }

  my @c;
  my $total_area;
  # triangulate
  foreach my $triangleref (@triangles) {
    my @triangle = @{$triangleref};
    my $area = TriangleArea([$triangle[0],$triangle[1],$triangle[2]]);
    # weighted centroid = area * centroid = area * sum / 3
    # we postpone division by 3 till we divide by total area to
    # minimize number of calculations
    $c[0] += ($triangle[0][0]+$triangle[1][0]+$triangle[2][0]) * $area;
    $c[1] += ($triangle[0][1]+$triangle[1][1]+$triangle[2][1]) * $area;
    $total_area += $area;
  }
  $c[0] = $c[0]/($total_area*3);
  $c[1] = $c[1]/($total_area*3);
  return \@c;
}
################################################################################
#
# The winding number method has been cused here.  Seems to
# be the most accurate one and, if well written, it matches
# the performance of the crossing number method.
# The winding number method counts the number of times a polygon
# winds around the point.  If the result is 0, the points is outside
# the polygon.
#
# args: reference to polygon object
#       reference to a point
#
sub IsInsidePolygon {
  my ($pointsref,$pointref) = @_;
  my @points = @$pointsref;
  if (@points < 3) { # polygon should at least have 3 points ...
    carp("Can't run inpolygon: polygon should have at least 3 points");
    return;
  }
  if (! $pointref) {
    carp("Can't run inpolygon: no point entered");
    return;
  }
  my @point = @$pointref;
  my $wn;  # thw winding number counter
  for (my $i = 0 ; $i < @points ; $i++) {
    if ($points[$i-1][1] <= $point[1]) { # start y <= P.y
      if ($points[$i][1] > $point[1]) {  # // an upward crossing
        if (CrossProduct([$points[$i-1],$points[$i],$pointref]) > 0) {
          # point left of edge
          $wn++;                         # have a valid up intersect
        }
      }
    } else {                             # start y > P.y (no test needed)
      if ($points[$i][1] <= $point[1]) { # a downward crossing
        if (CrossProduct([$points[$i-1],$points[$i],$pointref]) < 0) {
          # point right of edge
          $wn--;                         # have a valid down intersect
        }
      }
    }
  }
  return $wn;
}
################################################################################
#
# Check if polygon inside polygon or contour
# (for a contour, a point is inside when it's within the outer shape and
#  not within one of the inner shapes (holes) )
sub isinside {
  my ($self,$pointref) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    return IsInsidePolygon($pointsref,$pointref);
  } else {
    my @polygonrefs = $self->get_polygons;
    return 0 if (! IsInsidePolygon($polygonrefs[0],$pointref));
    my @result;
    for (my $i = 1; $i <@polygonrefs; $i++) {
      return 0 if (IsInsidePolygon($polygonrefs[$i],$pointref));
    }
    return 1;
  }
}
################################################################################
#
# a counter clockwise rotation over an angle a is given by the formula
#
#  / x2 \      /  cos(a)  -sin(a) \  / x1 \
#  |    |   =  |                  |  |    |
#  \ y2 /      \  sin(a)   cos(a) /  \ y1 /
#
# args: reference to polygon object
#       angle (in radians)
#       reference to center point (use origin if no center point entered)
#
sub RotatePolygon {
  my ($pointsref,$angle,$center) = @_;
  my $xc = 0;
  my $yc = 0;
  if ($center) {
    my @point = @$center;
    $xc = $point[0];
    $yc = $point[1];
  }
  if ($pointsref) {
    my @points = @$pointsref;
    for (my $i = 0 ; $i < @points ; $i++) {
      my $x = $xc + cos($angle)*($points[$i][0] - $xc) - sin($angle)*($points[$i][1] - $yc);
      my $y = $yc + sin($angle)*($points[$i][0] - $xc) + cos($angle)*($points[$i][1] - $yc);
      $points[$i][0] = $x;
      $points[$i][1] = $y;
    }
    return [@points];
  }
}
################################################################################
#
# rotate jpolygon or contour
#
sub rotate {
  my ($self,$angle,$center) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(RotatePolygon($pointsref,$angle,$center));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      push @{$self->{polygons}}, RotatePolygon($_,$angle,$center);
    }
  }
}
################################################################################
#
# move a polygon over a distance in x and y direction
#
# args: reference to polygon object
#       X offset
#       y offset
#
sub MovePolygon {
  my ($pointsref,$dx,$dy) = @_;
  if ($pointsref) {
    my @points = @$pointsref;
    for (my $i = 0 ; $i < @points ; $i++) {
      $points[$i][0] = $points[$i][0] + $dx;
      $points[$i][1] = $points[$i][1] + $dy;
    }
    return [@points];
  }
}
################################################################################
#
# Move polygon or contour
#
sub move {
  my ($self,$dx,$dy) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(MovePolygon($pointsref,$dx,$dy));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      $self->add_polygon(MovePolygon($_,$dx,$dy));
    }
  }
}
################################################################################
#
# mirror in x direction - vertical axis through point referenced by $center
# if no center entered, use y axis
#
# args: reference to polygon object
#       reference to center
#
sub MirrorXPolygon {
  my ($pointsref,$center) = @_;
  my @points = @$pointsref;
  if (@points == 0) { # nothing to mirror
    carp("Nothing to mirror ...");
    return;
  }
  my $xc = 0;
  my $yc = 0;
  if ($center) {
    my @point = @$center;
    $xc = $point[0];
    $yc = $point[1];
  }
  for (my $i = 0 ; $i < @points ; $i++) {
    $points[$i][0] = 2*$xc - $points[$i][0];
  }
  return [@points];
}
################################################################################
#
# mirror polygon or contour in x direction
#    (vertical axis through point referenced by $center)
sub mirrorx {
  my ($self,$dx,$dy) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(MirrorXPolygon($pointsref,$dx,$dy));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      $self->add_polygons(MirrorXPolygon($_,$dx,$dy));
    }
  }
}
################################################################################
#
# mirror in y direction - horizontal axis through point referenced by $center
# if no center entered, use x axis
#
# args: reference to polygon object
#       reference to center
#
sub MirrorYPolygon {
  my ($pointsref,$center) = @_;
  my @points = @$pointsref;
  if (@points == 0) { # nothing to mirror
    carp("Nothing to mirror ...");
    return;
  }
  my $xc = 0;
  my $yc = 0;
  if ($center) {
    my @point = @$center;
    $xc = $point[0];
    $yc = $point[1];
  }
  for (my $i = 0 ; $i < @points ; $i++) {
    $points[$i][1] = 2*$yc - $points[$i][1];
  }
  return [@points];
}
################################################################################
#
# mirror polygon or contour in x direction
#    (vertical axis through point referenced by $center)
sub mirrory {
  my ($self,$dx,$dy) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(MirrorYPolygon($pointsref,$dx,$dy));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      $self->add_polygons(MirrorYPolygon($_,$dx,$dy));
    }
  }
}
################################################################################
#
# mirror around axis determined by 2 points (p1p2)
#
# args: reference to polygon object
#       reference to array with to points defining reflection axis
#
sub MirrorPolygon {
  my ($pointsref,$axisref) = @_;
  my @points = @$pointsref;
  my @axis   = @$axisref;
  if (@axis != 2) { # need 2 points defining axis
    carp("Can't mirror: 2 points need to define axis");
    return;
  }
  my $p1ref = $axis[0];
  my $p2ref = $axis[1];
  my @p1 = @$p1ref;
  my @p2 = @$p2ref;
  if (@points == 0) { # nothing to mirror
    carp("Nothing to mirror ...");
    return;
  }
  for (my $i = 0 ; $i < @points ; $i++) {
    my $footref = PerpendicularFoot([\@p1,\@p2,$points[$i]]);
    my @foot = @$footref;
    $points[$i][0] = $foot[0] - ($points[$i][0] - $foot[0]);
    $points[$i][1] = $foot[1] - ($points[$i][1] - $foot[1]);
  }
  return [@points];
}
################################################################################
#
# mirror polygon or contour around axis determined by 2 points (p1p2)
#
sub mirror {
  my ($self,$axisref) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(MirrorPolygon($pointsref,$axisref));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      $self->add_polygons(MirrorPolygon($_,$axisref));
    }
  }
}
################################################################################
#
# scale polygon from center
# I would choose the centroid ...
#
# args: reference to polygon object
#       scale factor
#       reference to center point
#
sub ScalePolygon {
  my ($pointsref,$scale,$center) = @_;
  my @points = @$pointsref;
  if (@points == 0) { # nothing to scale
    carp("Nothing to scale ...");
    return;
  }
  my $xc = 0;
  my $yc = 0;
  if ($center) {
    my @point = @$center;
    $xc = $point[0];
    $yc = $point[1];
  }
  # subtract center, scale and add center again
  for (my $i = 0 ; $i < @points ; $i++) {
    $points[$i][0] = $scale * ($points[$i][0] - $xc) + $xc;
    $points[$i][1] = $scale * ($points[$i][1] - $yc) + $yc;
  }
  return [@points];
}
################################################################################
#
# scale polygon from center
# I would choose the centroid ...
#
sub scale {
  my ($self,$scale,$center) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    $self->points(ScalePolygon($pointsref,$scale,$center));
  } else {
    my @polygonrefs = $self->get_polygons;
    $self->{polygons} = ();
    my @result;
    foreach (@polygonrefs) {
      $self->add_polygons(ScalePolygon($_,$scale,$center));
    }
  }
}
################################################################################
#
# The "bounding box" of a set of points is the box with horizontal
# and vertical edges that contains all points
#
# args: reference to array of points or a contour
# returns reference to array of 4 points representing bounding box
#
sub bbox {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  if (@points < 3) { # polygon should at least have 3 points ...
    carp("Can't determine bbox: polygon should have at least 3 points");
    return;
  }
  my $min_x = $points[0][0];
  my $min_y = $points[0][1];
  my $max_x = $points[0][0];
  my $max_y = $points[0][1];
  for (my $i = 1 ; $i < @points ; $i++) {
     $min_x = $points[$i][0] if ($points[$i][0] < $min_x);
     $min_y = $points[$i][1] if ($points[$i][1] < $min_y);
     $max_x = $points[$i][0] if ($points[$i][0] > $max_x);
     $max_y = $points[$i][1] if ($points[$i][1] > $max_y);
  }
  return [[$min_x,$min_y],
          [$min_x,$max_y],
          [$max_x,$max_y],
          [$max_x,$min_y]];
}
################################################################################
#
# The "minimal enclosing rectangle" of a set of points is the box with minimal area
# that contains all points.
# We'll use the rotating calipers method here which works only on convex polygons
# so before calling minbbox, create the convex hull first for the set of points
# (taking into account whether or not the set of points represents a polygon).
#
# args: reference to array of points representing a convex polygon
# returns reference to array of 4 points representing minimal bounding rectangle
#
sub minrectangle {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  if (@points < 3) { # polygon should at least have 3 points ...
    carp("Can't determine minrectangle: polygon should have at least 3 points");
    return;
  }
  my $d;
  # scan all segments and for each segment, calculate the area of the bounding
  # box that has one side coinciding with the segment
  my $min_area = 0;
  my @indices;
  for (my $i = 0 ; $i < @points ; $i++) {
    # for each segment, find the point (vertex) at the largest perpendicular distance
    # the opposite side of the current rectangle runs through this point
    my $mj;       # index of point at maximum distance
    my $maxj = 0; # maximum distance (squared)
    # Get coefficients of the implicit line equation ax + by +c = 0
    # Do NOT normalize since scaling by a constant
    # is irrelevant for just comparing distances.
    my $a = $points[$i-1][1] - $points[$i][1];
    my $b = $points[$i][0] - $points[$i-1][0];
    my $c = $points[$i-1][0] * $points[$i][1] - $points[$i][0] * $points[$i-1][1];
    # loop through point array testing for max distance to current segment
    for (my $j = -1 ; $j < @points-1 ; $j++) {
      next if ($j == $i || $j == $i-1); # exclude points of current segment
      # just use dist squared (sqrt not needed for comparison)
      # since the polygon is convex, all points are at the same side
      # so we don't need to take the absolute value for dist
      my $dist = $a * $points[$j][0] + $b * $points[$j][1] + $c;
      if ($dist > $maxj) {    # this point is further
          $mj   = $j;         # so have a new maximum
          $maxj = $dist;
      }
    }
    # the line -bx+ay+c=0 is perpendicular to ax+by+c=0
    # now find index of extreme points corresponding to perpendicular line
    # initialize to first point (note that points of current segment could
    # be one or even both of the extreme points)
    my $mk = 0;
    my $ml = 0;
    my $mink = -$b * $points[0][0] + $a * $points[0][1] + $c;
    my $maxl = -$b * $points[0][0] + $a * $points[0][1] + $c;
    for (my $j = 1 ; $j < @points ; $j++) {
      # use signed dist to get extreme points
      my $dist = -$b * $points[$j][0] + $a * $points[$j][1] + $c;
      if ($dist < $mink) {    # this point is further
          $mk   = $j;         # so have a new maximum
          $mink = $dist;
      }
      if ($dist > $maxl) {    # this point is further
          $ml   = $j;         # so have a new maximum
          $maxl = $dist;
      }
    }
    # now $maxj/sqrt(a**2+b**2) is the height of the current rectangle
    # and (|$mink| + |$maxl|)/sqrt(a**2+b**2) is the width
    # since area is width*height we can waste the costly sqrt function
    my $area = abs($maxj * ($mink-$maxl)) / ($a**2 +$b**2);
    if ($area < $min_area || ! $min_area) {
      $min_area = $area;
      @indices = ($i,$mj,$mk,$ml);
    }
  }
  my ($i,$j,$k,$l) = @indices;
  # Finally, get the corners of the minimum enclosing rectangle
  my $p1 = PerpendicularFoot([$points[$i-1],$points[$i],$points[$k]]);
  my $p2 = PerpendicularFoot([$points[$i-1],$points[$i],$points[$l]]);
  # now we calculate the second point on the line parallel to
  # the segment i going through the vertex j
  my $p  = [$points[$j][0]+$points[$i-1][0]-$points[$i][0],
            $points[$j][1]+$points[$i-1][1]-$points[$i][1]];
  my $p3 = PerpendicularFoot([$points[$j],$p,$points[$l]]);
  my $p4 = PerpendicularFoot([$points[$j],$p,$points[$k]]);
  return [$p1,$p2,$p3,$p4];
}
################################################################################
#
# triangulate polygon or contour
#
# args: polygon or contour object
# returns a reference to an array triangles
#
sub triangulate {
  my ($self) = @_;
  my $pointsref = $self->points;
  if ($pointsref) {
    return TriangulatePolygon([$pointsref]);
  } else {
    my $polygonrefs = $self->polygons;
    if ($polygonrefs) {
      return TriangulatePolygon($polygonrefs);
    }
  }
}
################################################################################
#
# convexhull using the Melkman algorithm 
# (the set of input points represent a polygon and are thus ordered
#
# args: reference to ordered array of points representing a polygon
#       or contour (for a contour, we calculate the hull for the
#       outer shape)
# returns a reference to an array of the convex hull vertices
#
sub convexhull {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  return ([@points]) if (@points < 5);            # need at least 5 points
  # initialize a deque D[] from bottom to top so that the
  # 1st tree vertices of V[] are a counterclockwise triangle
  my @result;
  my $bot = @points-2;
  my $top = $bot+3;           # initial bottom and top deque indices
  $result[$bot] = $points[2]; # 3rd vertex is at both bot and top
  $result[$top] = $points[2]; # 3rd vertex is at both bot and top
  if (CrossProduct([$points[0], $points[1], $points[2]]) > 0) {
    $result[$bot+1] = $points[0];
    $result[$bot+2] = $points[1];       # ccw vertices are: 2,0,1,2
  } else {
    $result[$bot+1] = $points[1];
    $result[$bot+2] = $points[0];       # ccw vertices are: 2,1,0,2
  }

  # compute the hull on the deque D[]
  for (my $i=3; $i < @points; $i++) {   # process the rest of vertices
    # test if next vertex is inside the deque hull
    if ((CrossProduct([$result[$bot], $result[$bot+1], $points[$i]]) > 0) &&
      (CrossProduct([$result[$top-1], $result[$top], $points[$i]]) > 0) ) {
        last;         # skip an interior vertex
    }

    # incrementally add an exterior vertex to the deque hull
    # get the rightmost tangent at the deque bot
    while (CrossProduct([$result[$bot], $result[$bot+1], $points[$i]]) <= 0) {
      ++$bot;                      # remove bot of deque
      }
    $result[--$bot] = $points[$i]; # insert $points[i] at bot of deque

    # get the leftmost tangent at the deque top
    while (CrossProduct([$result[$top-1], $result[$top], $points[$i]]) <= 0) {
      --$top;                      # pop top of deque
      }
    $result[++$top] = $points[$i]; #/ push $points[i] onto top of deque
  }

  # transcribe deque D[] to the output hull array H[]
  my @returnval;
  for (my $h = 0; $h <= ($top-$bot-1); $h++) {
    $returnval[$h] = $result[$bot + $h];
  }

  return [(@returnval)];
}
################################################################################
#
# convexhull using Andrew's monotone chain 2D convex hull algorithm
# returns a reference to an array of the convex hull vertices
#
# args: reference to array of points (doesn't really need to be a polygon)
#       (also works for a contour - however, since a contour should consist
#       of polygons - which are ordered sets of points - the algorithm
#       above will be faster)
# returns a reference to an array of the convex hull vertices
#
sub convexhull2 {
  my ($self) = @_;
  my $pointsref = $self->points;
  if (! $pointsref) {
    $pointsref = ($self->get_polygons(0))[0];
    return if (! $pointsref); # empty object
  }
  my @points = @$pointsref;
  return ([@points]) if (@points < 5);            # need at least 5 points
  # first, sort the points by increasing x and y-coordinates
  @points = sort ByXY (@points);
  # Get the indices of points with min x-coord and min|max y-coord
  my @hull;
  my $bot = 0;
  my $top = -1;
  my $minmin = 0;
  my $minmax;
  my $xmin = $points[0][0];
  for (my $i = 1 ; $i < @points ; $i++) {
    if ($points[$i][0] != $xmin) {
      $minmax = $i - 1;
      last
    }
  }
  if ($minmax == @points-1) {      # degenerate case: all x-coords == xmin
    $hull[++$top] = $points[$minmin];
    if ($points[$minmax][1] != $points[$minmin][1]) { # a nontrivial segment
      $hull[$==$top] = $points[$minmax];
      return [@points];
    }
  }

  # Get the indices of points with max x-coord and min|max y-coord
  my $maxmin = 0;
  my $maxmax = @points - 1;
  my $xmax = $points[@points-1][0];
  for (my $i = @points - 2 ; $i >= 0 ; $i--) {
    if ($points[$i][0] != $xmax) {
      $maxmin = $i + 1;
      last;
    }
  }

  # Compute the lower hull on the stack @lower
  $hull[++$top] = $points[$minmin];    # push minmin point onto stack
  my $i = $minmax;
  while (++$i <= $maxmin) {
    # the lower line joins points[minmin] with points[maxmin]
    if (CrossProduct([$points[$minmin],$points[$maxmin],$points[$i]]) >= 0 && $i < $maxmin) {
      next;  # ignore points[i] above or on the lower line
    }
    while ($top > 0) {           # there are at least 2 points on the stack
      # test if points[i] is left of the line at the stack top
      if (CrossProduct([$hull[$top-1], $hull[$top], $points[$i]]) > 0) {
        last;                    # points[i] is a new hull vertex
      } else {
        $top--;
      }
    }
    $hull[++$top] = $points[$i]; # push points[i] onto stack
  }

  # Next, compute the upper hull on the stack H above the bottom hull
  if ($maxmax != $maxmin) {       # if distinct xmax points
    push @hull,$points[$maxmax];  # push maxmax point onto stack
  }
  $bot = $top;
  $i = $maxmin;
  while (--$i >= $minmax) {
    # the upper line joins points[maxmax] with points[minmax]
    if (CrossProduct([$points[$maxmax],$points[$minmax],$points[$i]]) >= 0 && $i > $minmax) {
      next;                        # ignore points[i] below or on the upper line
    }
    while ($top > $bot) {          # at least 2 points on the upper stack
      # test if points[i] is left of the line at the stack top
      if (CrossProduct([$hull[$top-1],$hull[$top],$points[$i]]) > 0) {
        last;                      # points[i] is a new hull vertex
      } else {
        $top--;
      }
    }
    $hull[++$top] = $points[$i];   # push points[i] onto stack
  }
  if ($minmax == $minmin) {
    shift @hull;                   # remove joining endpoint from stack
  }
  return [@hull];
}
################################################################################
#
# Sorting function to surt points first by X coordinate, then by Y coordinate
#
sub ByXY {
  my @p1 = @$a;
  my @p2 = @$b;
  my $result = $p1[0] <=> $p2[0];
  if ($result){
    return $result;
  } else {
    return $p1[1] <=> $p2[1];
  }
}
################################################################################
#
# convert polygon/contour to gpc contour
#
sub convert2gpc {
  my ($self,$dx,$dy) = @_;
  my @polygons;
  my $pointsref = $self->points;
  if ($pointsref) {
    push @polygons,$pointsref;
  } else {
    @polygons = $self->get_polygons;
  }
  foreach (@polygons) {
    my @points = @{$_};
    if (@points < 3) { # need at least 3 points
      carp("Can't convert to gpc structure: polygon should have at least 3 points");
      return;
    }
  }
  my $contour = GPC::new_gpc_polygon();
  GPC::gpc_polygon_num_contours_set($contour,scalar(@polygons));
  # array for hole pointers
  my $hole_array = GPC::int_array(scalar(@polygons));
  GPC::gpc_polygon_hole_set($contour,$hole_array);
  my $vlist = GPC::new_gpc_vertex_list();
  for (my $i = 0; $i < @polygons; $i++) {
    if ($i == 0) {
      GPC::int_set($hole_array,$i,0);
    } else {
      GPC::int_set($hole_array,$i,1);
    }
    my @points = @{$polygons[$i]};
    my @gpc_vertexlist;
    foreach my $vertex (@points) {
      my $v = GPC::new_gpc_vertex();
      GPC::gpc_vertex_x_set($v,$$vertex[0]);
      GPC::gpc_vertex_y_set($v,$$vertex[1]);
      push @gpc_vertexlist,$v;
    }
    my $va = create_gpc_vertex_array(@gpc_vertexlist);
    my $vl = GPC::new_gpc_vertex_list();
    GPC::gpc_vertex_list_vertex_set($vl,$va);
    GPC::gpc_vertex_list_num_vertices_set($vl,scalar(@points));
    GPC::gpc_vertex_list_set($vlist,$i,$vl);
  }
  GPC::gpc_polygon_contour_set($contour,$vlist);
  return $contour;
}
################################################################################
#
# convert gpc object to a set of contours
# A gpc contour object can consist of multiple outer shapes each having holes,
#
sub Gpc2Polygons {
  my ($gpc) = @_;
  my @result; # array with contours
  my @inner;  # array holding the inner polygons
  my @outer;  # array holding the outer polygons
  my $num_contours = GPC::gpc_polygon_num_contours_get($gpc);
  my $contour      = GPC::gpc_polygon_contour_get($gpc);
  my $hole_array   = GPC::gpc_polygon_hole_get($gpc);
  # for each shape of the gpc object
  for (my $i = 0 ; $i < $num_contours ; $i++) {
    my @polygon;
    # get the hole flag
    my $hole = GPC::int_get($hole_array,$i);
    # get the vertices
    my $vl = GPC::gpc_vertex_list_get($contour,$i);
    my $num_vertices = GPC::gpc_vertex_list_num_vertices_get($vl);
    my $va = GPC::gpc_vertex_list_vertex_get($vl);
    for (my $j = 0 ; $j < $num_vertices ; $j++) {
      my $v = GPC::gpc_vertex_get($va,$j);
      my $x = GPC::gpc_vertex_x_get($v);
      my $y = GPC::gpc_vertex_y_get($v);
      push @polygon,[$x,$y];
    }
    # create lists of inner and outer shapes
    if ($hole) {
      push @inner,[@polygon];
    } else {
      push @outer,[@polygon];
    }
  }
  # shortcut: if there is only one outer shape, we're done
  if (@outer == 1) {
    my $obj = Math::Geometry_2D->new;
    $obj->add_polygons(@outer,@inner);
    push @result,$obj;
  } else {
    foreach (@outer) {
      # create contour for each outer shape
      my $obj = Math::Geometry_2D->new;
      $obj->polygons([$_]);
      push @result,$obj;
      # if an inner shape has at least one point inside this
      # outer shape, it belongs to this outer shape (so all
      # points are inside it)
      my $i = 0;
      while ($i < @inner) {
        my @polygon = @{$inner[$i]};
        if ($obj->isinside($polygon[0])) {
          $obj->add_polygons($inner[$i]);
          splice @inner,$i,1;
        }
        $i++;
      }
    }
  }
  return @result;
}
################################################################################
#
# gpc polygon clipping operatins
#
sub GpcClip {
  my ($op,$gpc_poly_1,$gpc_poly_2) = @_;
  my $result = GPC::new_gpc_polygon();
  SWITCH: {
    ($op eq "DIFFERENCE") && do {
      GPC::gpc_polygon_clip(0,$gpc_poly_1,$gpc_poly_2,$result);
      return $result;
    };
    ($op eq "INTERSECTION") && do {
      GPC::gpc_polygon_clip(1,$gpc_poly_1,$gpc_poly_2,$result);
      return $result;
    };
    ($op eq "XOR") && do {
      GPC::gpc_polygon_clip(2,$gpc_poly_1,$gpc_poly_2,$result);
      return $result;
    };
    ($op eq "UNION") && do {
      GPC::gpc_polygon_clip(3,$gpc_poly_1,$gpc_poly_2,$result);
      return $result;
    };
    return;
  }
}
###############################################################################
#
# create gpc vertex array pointer
#
sub create_gpc_vertex_array {
  my $len = scalar(@_);
  my $va = GPC::gpc_vertex_array($len);
  for (my $i=0; $i<$len; $i++) {
    my $val = shift;
    GPC::gpc_vertex_set($va,$i,$val);
  }
  return $va;
}
################################################################################
#
my $pi = atan2(1,1) * 4;
#
################################################################################
#
# convert a circle to a polygon
# arguments: first argument is the number of segments,
#            the other arguments are:
#  p1,p2,p3       : 3 points on the circle
# or
#  center,p1      : center and a point on the circle
# or
#  center,radius  : the center and the radius of the circle
#
sub CircleToPoly {
  my @args = @_;
  my @result;
  my ($segments,$p1,$p2,$p3,$center,$radius);
  if (@args == 4) {      # 3 points
    ($segments,$p1,$p2,$p3) = @args;
    $center = CalcCenter($p1,$p2,$p3);
    $radius = SegmentLength([$p1,$center]);
  } elsif (@args == 3) {
    if (ref $args[2]) {  # center + 1 point
      ($segments,$center,$p1) = @args;
      $radius = SegmentLength([$p1,$center]);
    } else {             # center + radius
      ($segments,$center,$radius) = @args;
    }
  } else {
    return;
  }
  my $angle = ($pi * 2) / $segments;
  for (my $i = 0 ; $i < $segments ; $i++) {
    push @result, [${$center}[0] + $radius * cos($angle * $i),
                   ${$center}[1] + $radius * sin($angle * $i)]
  }
  return [@result];
}
################################################################################
#
# convert an arc to a polygon
# arguments: first argument is the number of segments,
#            the other arguments are:
#  p1,p2,p3          : startpoint, intermediate point, endpoint
# or
#  $center,p1,p2,$dir : center, startpoint, endpoint,  direction
#                       direction 0 counter clockwise
#                                 1 clockwise
# Note: the return value is a set of points, NOT a closed polygon !!!
#
sub ArcToPoly {
  my @args = @_;
  my @result;
  my ($segments,$p1,$p2,$p3,$center,$direction);
  my ($radius,$angle);
  my ($start_angle, $end_angle);
  if (@args == 4) {      # 3 points
    ($segments,$p1,$p2,$p3) = @args;
    $center = CalcCenter($p1,$p2,$p3);
    $radius = SegmentLength([$p1,$center]);
    # calculate start and end angles
    $start_angle  = CalcAngle($center,$p1);
    my $mid_angle = CalcAngle($center,$p2);
    $end_angle    = CalcAngle($center,$p3);
    if ( (($mid_angle   < $start_angle) && ($start_angle < $end_angle)) ||
         (($start_angle < $end_angle)   && ($end_angle   < $mid_angle)) ||
         (($end_angle   < $mid_angle)   && ($mid_angle   < $start_angle)) ) {
      $direction = 1;
    }
    $angle = $end_angle - $start_angle;
  } elsif (@args == 5) {  # center, begin, end, direction
    ($segments,$center,$p1,$p3,$direction) = @args;
    $radius = SegmentLength([$p1,$center]);
    # calculate start and end angles
    $start_angle = CalcAngle($center,$p1);
    $end_angle   = CalcAngle($center,$p3);
    $angle = $end_angle - $start_angle;
  } else {
    return;
  }

  if ($direction) {  # clockwise
    if ($angle > 0) {
      $angle = $angle - ($pi * 2);
    }
  } else {
    if ($angle < 0) {
      $angle = $angle + ($pi * 2);
    }
  }
  $angle = $angle / $segments;

  push @result,$p1; # start point
  for (my $i = 1 ; $i < $segments ; $i++) {
    push @result, [${$center}[0] + $radius * cos($start_angle + $angle * $i),
                   ${$center}[1] + $radius * sin($start_angle + $angle * $i)]
  }
  push @result,$p3; # end point
  return [@result];
}
################################################################################
#
# Calculate the center of a circle going through 3 points
#
sub CalcCenter {
  my ($p1_ref, $p2_ref, $p3_ref) = @_;
  my ($x1,$y1) = @{$p1_ref};
  my ($x2,$y2) = @{$p2_ref};
  my ($x3,$y3) = @{$p3_ref};
  # calculate midpoints of line segments p1p2 p2p3
  my $u1 = ($x1 + $x2)/2;
  my $v1 = ($y1 + $y2)/2;
  my $u2 = ($x2 + $x3)/2;
  my $v2 = ($y2 + $y3)/2;
  # linear equations y = a + bx
  my ($a1,$a2);
  my ($b1,$b2);
  # intersect (center) coordinates
  my ($xi,$yi);
  # slope of perpendicular = -1/slope
  if ($y1 != $y2) {
    $b1 = - ($x1 - $x2)/($y1 - $y2);
    $a1 = $v1 - $b1 * $u1;
  } else {
    $xi = $u1;
  }
  if ($y2 != $y3) {
    $b2 = - ($x2 - $x3)/($y2 - $y3);
    $a2 = $v2 - $b2 * $u2;
  } else {
    $xi = $u2;
  }
  # parallel lines (colinear is also parallel)
  return if ($b1 == $b2 || (!$b1 && !$b2));
  $xi = - ($a1 - $a2)/($b1 - $b2) if (!$xi);
  $yi = $a1 + $b1 * $xi if ($b1);
  $yi = $a2 + $b2 * $xi if ($b1);
  return [($xi,$yi)];
}
################################################################################
#
# calculate angel of vector p1p2
#
sub CalcAngle {
  my ($p1_ref,$p2_ref) = @_;
  my ($x1,$y1) = @{$p1_ref};
  my ($x2,$y2) = @{$p2_ref};
  return   0   if ($y1 == $y2 && $x1 == $x2);
  return   0   if ($y1 == $y2 && $x1 < $x2);
  return $pi   if ($y1 == $y2 && $x1 > $x2);
  return $pi/2 if ($x1 == $x2 && $y1 < $y2);
  return ($pi *3)/2 if ($x1 == $x2 && $y1 > $y2);
  my $angle = atan2($y2-$y1,$x2-$x1);
  return $angle;
}
################################################################################
1;
