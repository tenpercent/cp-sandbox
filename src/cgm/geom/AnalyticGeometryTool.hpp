//-------------------------------------------------------------------------
// COPYRIGHT 1998  CATERPILLAR INC.  ALL RIGHTS RESERVED
//
// Filename      : AnalyticGeometryTool.hpp
//
// Purpose       : This class performs calculations on analytic geometry
//                 (points, lines, arcs, planes, polygons).  Capabilities 
//                 include vector and point math, matrix operations, 
//                 measurements, intersections and comparison/containment checks.
//
// Special Notes : As most of these functions were taken from Caterpillar's
//                 Cmath library with minimal modification for this class,
//                 the vectors, planes, etc. are represented in non-CUBIT
//                 format.  These could be converted with some effort.
// 
//                 Other Notes:
//                 ------------
//                 1. Points and vectors are represented as double[3].
//
//                 2. Lines (unless denoted as 'segments') are defined as 
//                    a point and a vector (point-direction form).  The parametric
//                    equations of a line are:
//
//                         x = xo + t*a
//                         y = yo + t*a
//                         z = zo + t*a
//
//                 3. Planes are also defined as a point and a vector
//                    (point-normal form).  The equation of a plane is:
//
//                         a*(x-xo) + b*(y-yo) + c*(z-zo) = 0
//                      
//                     For example, the equation of a plane passing through the 
//                     point (3,-1,7) and perpendicular to the vector n = (4,2,-5) 
//                     is: 4(x-3) + 2(y+1) -5(z-7) = 0
//
//                     This can be reduced to:
//
//                         A*x + B*y + C*z + D = 0
//  
//                     For the example, the equation is: 4x + 2y -5z + 25 = 0
//
//                 4. Most functions which can accept a 3x3 matrix have
//                    an analogous 4x4 function.  This is for convenience
//                    only (the extra rows/columns are ignored within the function).
//                 
//
// Creator       : Steve Storm
//
// Creation Date : 10/16/98
//-------------------------------------------------------------------------

#ifndef ANALYTIC_GEOMETRY_TOOL_HPP
#define ANALYTIC_GEOMETRY_TOOL_HPP

#include <math.h>
#include "CubitDefines.h"
#include "CubitGeomConfigure.h"

#include <iostream>

class CubitBox;
class CubitPlane;
class CubitVector;
template <class X> class DLIList;

  // Multiplication constants
  const double AGT_E =             2.7182818284590452354;	// e
  const double AGT_LOG_2E =        1.4426950408889634074;	// log2(e)
  const double AGT_LOG_10E =       0.43429448190325182765; // log10(e)
  const double AGT_LN_2 =          0.69314718055994530942; // ln(2)
  const double AGT_LN_10 =         2.30258509299404568402; // ln(10)
  const double AGT_PI  =           3.14159265358979323846; // PI
  const double AGT_PI_DIV_2 =      1.57079632679489661923; // PI/2
  const double AGT_PI_DIV_4 =      0.78539816339744830962; // PI/4
  const double AGT_1_DIV_PI =      0.31830988618379067154; // 1/PI
  const double AGT_2_DIV_PI =      0.63661977236758134308; // 2/PI
  const double AGT_2_DIV_SQRT_PI = 1.12837916709551257390; // 2/(sqrt(PI))
  const double AGT_SQRT_2 =        1.41421356237309504880; // sqrt(2)
  const double AGT_SQRT_1_DIV_2 =  0.70710678118654752440; // sqrt(1/2)
  const double AGT_PI_DIV_180 =    0.017453292519943295;   // PI/180
  const double AGT_180_DIV_PI =    57.295779513082323;     // 180/PI 

  const double DEGtoRAD = AGT_PI_DIV_180;
  const double RADtoDEG = AGT_180_DIV_PI;

  typedef struct {
     double Vec1[3];	/* First vector that defines plane of arc */
     double Vec2[3];	/* Second vector that defines plane of arc */
     double Origin[3];	/* Center of arc (on plane) */
     double StartAngle;	/* Starting angle (in radians) of arc */
     double EndAngle;	/* End angle (in radians) of arc */
     double Radius;	/* Radius of arc */
  } AGT_3D_Arc;

  // Return values */
  enum AgtEquality {
     AGT_EQUAL       = 0, // define for comparisons
     AGT_LESSTHAN    = 1, // define for comparisons
     AGT_GREATERTHAN = 2  // define for comparisons
  };

  enum AgtSide {
     AGT_ON_PLANE =  0,  // define for which side of plane
     AGT_POS_SIDE =  1,  // define for which side of plane
     AGT_NEG_SIDE =  2,  // define for which side of plane
     AGT_INT_PLANE = 3,  // define for intersects plane
     AGT_CROSS = 4       // define for line crossing plane
  };

  enum AgtOrientation {
    AGT_OPP_DIR = 0,  // define for vector direction comparision
    AGT_SAME_DIR = 1  // define for vector direction comparision
  };

  enum AgtDistanceMethod {
    AGT_FRACTION = 0,  // define for distance along a curve
    AGT_ABSOLUTE = 1   // define for distance along a curve
  };

// Macros
#define number_sign(number) (number<0 ? -1 : 1)
#define copy_vec(vec1,vec2) copy_pnt(vec1,vec2)
#define set_vec(pnt,x,y,z) set_pnt(pnt,x,y,z)

///////////////////////////////////////////////////////////////////////////////
//                               MAGIC SOFTWARE                              //
// See note later in file regarding MAGIC softare.                           //
///////////////////////////////////////////////////////////////////////////////
// Minimal 2D bounding box - parameters
const int maxPartition = 32;
const int invInterp = 32;
const double angleMin = -0.25*AGT_PI;
const double angleMax = +0.25*AGT_PI;
const double angleRange = 0.5*AGT_PI;

typedef struct
{
    double x, y;
}
Point2;

typedef struct
{
    double x, y, z;
}
Point3;

typedef struct
{
    Point2 center;
    Point2 axis[2];
    double extent[2];
}
OBBox2;

typedef struct
{
    Point3 center;
    Point3 axis[3];
    double extent[3];
}
OBBox3;

// threshold on determinant to conclude lines/rays/segments are parallel
const double par_tolerance = 1e-06;
#ifndef INFINITY
const double INFINITY = CUBIT_DBL_MAX;
#endif
const Point3 PT_INFINITY = { INFINITY, INFINITY, INFINITY };
#define DIST(x) (Sqrt(Abs(x)))

// Line in 2D
typedef struct
{
    // line is Dot(N,X) = c, N not necessarily unit length
    Point2 N;
    double c;
}
AgtLine;

// Triangle in 2D
typedef struct
{
    Point2 v[3];
}
Triangle;

// Linked list of triangle
typedef struct AgtTriList
{
    Triangle* tri;
    AgtTriList* next;
}
AgtTriList;

//---------------------------------------------------------------------------
// Lines in 3D
//---------------------------------------------------------------------------
typedef struct
{
    // Line is L(t) = b+t*m for any real-valued t
    // Ray has constraint t >= 0, b is the origin of the ray
    // Line segment has constraint 0 <= t <= 1, b and b+m are end points

    Point3 b, m;
}
Line3;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Circles in 3D
//---------------------------------------------------------------------------
typedef struct
{
    // Plane containing circle is Dot(N,X-C) = 0 where X is any point in the
    // plane.  Vectors U, V, and N form an orthonormal right-handed set
    // (matrix [U V N] is orthonormal and has determinant 1).  Circle within
    // the plane is parameterized by X = C + R*(cos(A)*U + sin(A)*V) where
    // A is an angle in [0,2*pi).

    Point3 U, V, N;  // coordinate system of plane containing circle
    Point3 C;  // center
    double R;  // radius > 0
}
Circle3;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Triangles in 3D
//---------------------------------------------------------------------------
typedef struct
{
    // Triangle points are tri(s,t) = b+s*e0+t*e1 where 0 <= s <= 1,
    // 0 <= t <= 1, and 0 <= s+t <= 1.

    // If the vertices of the triangle are v0, v1, and v2, then b = v0,
    // e0 = v1 - v0, and e1 = v2 - v0.

    Point3 b, e0, e1;
}
Triangle3;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Parallelograms in 3D
//---------------------------------------------------------------------------
typedef struct
{
    // Rectoid points are rect(s,t) = b+s*e0+t*e1 where 0 <= s <= 1
    // and 0 <= t <= 1.  Could have called this a parallelogram, but
    // I do not like typing long words...

    Point3 b, e0, e1;
}
Rectoid3;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Planes in 3D
//---------------------------------------------------------------------------
typedef struct
{
    // plane is Dot(N,X) = d
    Point3 N;
    double d;
}
Plane3;
//---------------------------------------------------------------------------

// END OF MAGIC SOFTWARE
///////////////////////////////////////////////////////////////////////////////

class CUBIT_GEOM_EXPORT AnalyticGeometryTool
{
public:

   ~AnalyticGeometryTool();
   static AnalyticGeometryTool* instance();

   //*********************************************************
   // Double numbers
   //*********************************************************
   CubitBoolean dbl_eq( double val_1, double val_2 );
   // Check to see if double numbers are equal within epsilon.
   // Values are equal if fabs(val_1 - val_2) < epsilon.
   // Epsilon is determined by next two functions.
   double get_epsilon();
   double set_epsilon( double new_epsilon );
   // Get/set epsilon used for double number comparisons.
   // Default value is 1e-6.
#ifdef BOYD14
   void swap_dbl(double& dbl1,double& dbl2);
   // Swap two double number memory locations.
   double min_array( const double *array, int len, int &pos );
   double max_array( const double *array, int len, int &pos );
   // Find minimum or maximum of an array of double numbers.
   // Position in array is returned as well.
#endif
   void round_near_val( double& val );
   // Round value to either -1, 0, or 1 if within epsilon to
   // to one of these.  Epsilon determined by get_epsilon/set_epsilon
   // functions.

   //**************************************************************************
   // Matrices & Transforms
   //**************************************************************************
   // Note: For these functions the matrix format is defined as follows:
   //    
   //    Consider the transformation from [x,y,z] to [x',y',z']:    
   //                             _          _
   //                            | x1 y1 z1 0 |
   //    [x',y',z',1] = [x,y,z,1]| x2 y2 z2 0 |
   //                            | x3 y3 z3 0 |
   //			                     | ox oy oz 1 |
   //			                      -          -
   void transform_pnt( double m[4][4], double pin[3], double pout[3] );
   // This functions applies the transformation matrix to the specified
   // point and returns the new point coordinates through the call list.
   // Point in and point out can be the same memory address or different.
#ifdef BOYD14
   void transform_pnt_arr( double m[4][4], double pin_arr[][3], 
                           size_t num_pnts, double pout_arr[][3] );
   // Transform a point array.  Points can be transformed in place.
#endif
   void transform_vec( double m3[3][3], double vin[3], double vout[3] );
   void transform_vec( double m4[4][4], double vin[3], double vout[3] );
   // This routine applies a transformation matrix to a specified vector
   // and returns the new vector orientation.  Vector in and vector out can
   // be the same memory address or different.
#ifdef BOYD14
   double transform_dist( double m3[3][3], double distance, double v[3] );
   double transform_dist( double m4[4][4], double distance, double v[3] );
   // Transforms a distance with given transform matrix
   // The routine can accept the vector direction (in non-transformed system) 
   // of the distance.  This is probably a very unusual case - meaning that,
   // compared to the old system, the new system is not square (a different
   // scale factor can apply to each component direction - ie., the distance
   // 6 units covers in the x-direction may not equal the distance 6 units
   // covers in the y-direction - a rectangular system!).
#endif
   void transform_line( double rot_mtx[4][4], double origin[3], double axis[3] );
   void transform_line( double rot_mtx[4][4], CubitVector &origin, CubitVector &axis );
   // Transforms a line using the given transformation mtx.  The mtx is typically
   // built using add_origin_to_rotation_mtx.
   void copy_mtx( double from[3][3],double to[3][3] );
   void copy_mtx( double from[4][4], double to[4][4] );
   void copy_mtx( double from[4][4], double to[3][3] );
   void copy_mtx(double from[3][3], double to[4][4], double *origin = NULL );
   // This routine simply copies one matrix to another.  If a NULL is passed in
   // for the from matrix, then the identity matrix is copied into the out matrix.
   // Last function allows you to specify 1st 3 doubles of row 4 (origin) of the 
   // 4x4 matrix.
   void create_rotation_mtx( double theta, double v[3], double mtx3x3[3][3] );
   void create_rotation_mtx( double theta, double v[3], double mtx4x4[4][4] );
   // This routine determines the tranformation matrix given the theta and
   // the vector to cross through.
   void add_origin_to_rotation_mtx( double rot_mtx[4][4], double origin[3] );
   // Adds origin to rotation matrix, so you can rotate points about a line.
   // Line is defined as original vector used in create_rotation_mtx and
   // the origin.
   void identity_mtx( double mtx3x3[3][3] );
   void identity_mtx( double mtx4x4[4][4] );
   // Simply sets the given matrix to the identity matrix.
   void mtx_to_angs( double mtx3x3[3][3], double &ax, double &ay, double &az );
   void mtx_to_angs( double mtx4x4[4][4], double &ax, double &ay, double &az );
   // Gets rotation angles to rotate one system to another - returned rotation 
   // angles are about global system.
   // To use the result from this function, align the unoriented object's origin 
   // with the global origin, then apply the rotation angles returned from this 
   // routine to the object in the order of ax,ay,az about the global origin.  
   // The object will be oriented such that its xyz axes will point in the same
   // direction as the object whose transformation matrix was given.  The object
   // can then be translated.
   void rotate_system( double mtx[3][3], double sys_in[3][3], 
                       double sys_out[3][3] );
   void rotate_system( double mtx[4][4], double sys_in[4][4], 
                       double sys_out[4][4] );
   // This routine rotates a 3x3 system given a transformation matrix.  The 
   // matrix can be rotated in place by sending same variable in & out, or a 
   // new matrix can be created.  In the 4x4 case, the translation portion
   // of the matrix is unchanged.
#ifdef BOYD14
   double mtx_to_ratio( double m3[3][3] );
   double mtx_to_ratio( double m4[4][4] );
   // This function returns the ratio of a length in the second system to a 
   // length in the first system.  Multiply the length in system 1 by this 
   // ratio to get the length in system 2.
#endif
   double det_mtx( double m[3][3] );
   // Find determinant of matrix.
   void mult_mtx( double a[3][3],double b[3][3], double d[3][3] );
   void mult_mtx( double a[4][4], double b[4][4], double d[4][4] );
   // Multiply matrices together.  If any input is NULL, the identity
   // matrix is used in its place.
   CubitStatus inv_mtx_adj( double mtx[3][3], double inv_mtx[3][3] );
   //This routine inverts a 3x3 matrix using the adjoint method.  If NULL 
   // is sent in for first matrix, the second matrix is set to the identity 
   // matrix.  If the same memory address is passed in for both incoming and 
   // outgoing matrices, the matrix is inverted in place.  Returns CUBIT_FAILURE
   // if no inverse exists.
   CubitStatus inv_trans_mtx( double transf[4][4], double inv_transf[4][4] );
   // This routine inverts a 4x4 transformation matrix.  If NULL is sent in 
   // for first matrix, the second matrix is set to the identity matrix.  If 
   // the same memory address is passed in for both incoming and outgoing matrices, 
   // the matrix is inverted in place.  Uses LU decomposition.
   void vecs_to_mtx( double xvec[3], double yvec[3], double zvec[3], 
                     double matrix[3][3] );
   void vecs_to_mtx( double xvec[3], double yvec[3], double zvec[3], 
                     double origin[3], double matrix[4][4] );
   // Creates a transformation matrix given three vectors (and origin
   // for 4x4 case).
#ifdef BOYD14
   void mtx_to_vecs( double matrix[3][3], double xvec[3], double yvec[3], 
                     double zvec[3] );
   void mtx_to_vecs( double matrix[4][4], double xvec[3], double yvec[3], 
                     double zvec[3], double origin[3] );
   // Extract vectors from a matrix.  NULL can be input for any of the
   // arguments and the routine will not extract that vector.  NULL can
   // also be input for the matrix, in which case the function assumes
   // the identity matrix.
#endif

   void print_mtx( double mtx[3][3] );
   void print_mtx( double mtx[4][4] );
   // Prints matrix values, for debugging.

   //**************************************************************************
   // 3D Points
   //**************************************************************************
   void copy_pnt( double pnt_in[3], double pnt_out[3] );
   // Copy one double[3] point to another double[3] point (uses memcpy).
   // If first point in NULL then second point is set to 0,0,0.
   void copy_pnt( double pnt_in[3], CubitVector &cubit_vec );
   void copy_pnt( CubitVector &cubit_vec, double pnt_out[3] );
   // For going back and forth from CubitVectors
   void set_pnt( double pnt[3], double x, double y, double z );
   // Sets the value of pnt to x,y,z (inline function).
#ifdef BOYD14
   void swap_pnt(double pnt1[3],double pnt2[3]);
   // Swaps two double[3] points with each other.
#endif
   CubitBoolean pnt_eq( double pnt1[3],double pnt2[3] );
   // Compares two points to determine if they are equivalent.  The
   // equivalence is based on epsilon (get_epsilon, set_epsilon).
   // If the distance between them is less than epsilon, the points
   // are considered equal.
#ifdef BOYD14
   CubitBoolean pnt_arr_eq( double **pnt_arr1, long num1, double **pnt_arr2,
                           long num2, long start_pos1 = 0, long start_pos2 = 0,
                           long num_to_check = -1 );
   // Compares two point arrays or portions of two point arrays to 
   // determine if they are equivalent within epsilon (get_epsilon,
   // set_epsilon).  Checking will stop when either of the array ends
   // is encountered or num_to_check is exceeded.  Note that he function
   // can return CUBIT_TRUE even if the arrays are of a different length.
#endif
   CubitStatus mirror_pnt( double pnt[3], double pln_orig[3], 
                           double pln_norm[3], double m_pnt[3]);
   // Function to mirror a point about a plane.  The mirror point
   // is the same distance from the plane as the original, but on 
   // the opposite side of the plane.  Function returns CUBIT_FAILURE
   // if the point is on the plane - in which case the point is
   // just copied.
#ifdef BOYD14
   double ** mirror_pnt_arr( double **pnt_arr, long num_pnts, 
                             double pln_orig[3], double pln_norm[3],
                             double **mirror_arr );
   // Mirrors an array of points about a plane.  If mirror_arr is
   // non-NULL, function assumes it was allocated to correct lenght.
   // If NULL, function allocates the array.
#endif
   CubitStatus next_pnt( double str_pnt[3], double vec_dir[3], double len,
                         double new_pnt[3]);
   // Given start pnt, vec dir and length find next point in space.  The
   // vector direction is first unitized then the following formula is used:
   //  new_point[] = start_point[] + (length * unit_vector[])
   // new_pnt can be the same as input point (overwrites old location).
#ifdef BOYD14
   CubitStatus pnt_between2pnts( double pnt1[3], double pnt2[3], 
                                 AgtDistanceMethod dist_meth, double dist,
                                 double out_pnt[3] );
   // This routine creates a point between two other points, at a specified
   // distance between the two given points.  The distance can be either
   // absolute (meaning an actual distance, ie. in mm), or a normalized distance
   // (ie., 50% or 0.5 between the other two points).  The reference point for
   // the distance is always the first point. 
   // - The normalized distance does not need to be between 0.0 & 1.0.  For example,
   // a 1.25 normalized distance would result in a point 25% beyond the second 
   // point (not actually "between" the two points).   
   // - The distances can also be negative.  For example, a -1.25 normalized
   // distance would result in a point 25% before the first point.  A -50.0
   // absolute distance would extend 50.0 from the first point in the opposite
   // direction of the second point.
#endif

   //***************************************************************************
   // 3D Vectors
   //***************************************************************************
   CubitStatus unit_vec( double vin[3], double vout[3] );
   // Finds unit vector of input vector (vout can equal vin - converts in place).
   double dot_vec( double uval[3], double vval[3] );
   // Dots two vectors.  Property of dot product is:
   //   angle between vectors acute if dot product > 0
   //   angle between vectors obtuse if dot product < 0
   //   angle between vectors 90 deg if dot product = 0
   void cross_vec( double uval[3], double vval[3], double cross[3] );
   // Finds cross product of two vectors:
   //    cross[0] = u[1] * v[2] - u[2] * v[1];
   //    cross[1] = u[2] * v[0] - u[0] * v[2];
   //    cross[2] = u[0] * v[1] - u[1] * v[0];
   void cross_vec_unit( double uval[3], double vval[3], double cross[3] );
   // Finds unit vector that is cross product of two vectors.
   void orth_vecs( double uvect[3], double vvect[3], double wvect[3] );
   // Finds 2 arbitrary orthoganal vectors to the first.
   double mag_vec( double vec[3] );
   // Finds the magnitude of a vector:
   //    magnitude = sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
   CubitStatus get_vec ( double str_pnt[3], double stp_pnt[3],
                         double vector_out[3] );
   // Finds a vector from str_pnt to stp_pnt.
   // Returns CUBIT_FAILURE if points are coincident.
   CubitStatus get_vec_unit( double str_pnt[3], double stp_pnt[3],
                             double uv_out[3] );
   // Finds a unit vector from str_pnt to stp_pnt.
   // Returns CUBIT_FAILURE if points are coincident.
   void mult_vecxconst( double constant, double vec[3], double vec_out[3] );
   // Multiples a vector by a constant (vec_out can equal vec).
#ifdef BOYD14
   void add_vec_vec( double vec_1[3], double vec_2[3], double vec_out[3]);
   // Adds two vectors (vec_out can equal vec_1 or vec_2)
   void sub_vec_vec( double vec_1[3], double vec_2[3], double vec_out[3] );
   // Subtracts two vectors (vec_out can equal vec_1 or vec_2).
#endif
#ifdef BOYD14
   void rotate_vec( double v1[3], double v2[3], double angle, double v_out[3] );
   // Rotate a vector about another vector a given angle.  v_out can equal
   // v1 (v_out replaces v1) or be a new vector.
#endif
   void reverse_vec( double vin[3],double vout[3] );
   // Multiples -1.0 by a vector (vout can equal vin).
   double angle_vec_vec( double v1[3],double v2[3] );
   // Finds angle in radians between two vectors.  Angle will always be
   // between 0.0 and PI.  For accuracy, the routine uses the cosine for large
   // angles and the sine for small angles.
   
   //***************************************************************************
   // Distances
   //***************************************************************************
   double dist_pnt_pnt( double pnt1[3], double pnt2[3] );
   // Determines distance between two points in space.
#ifdef BOYD14
   double dist_ln_ln( double ln_1_orig[3], double ln_1_vec[3],
                      double ln_2_orig[3], double ln_2_vec[3],
                      unsigned short *num_ints = NULL, double int_pnt1[3] = NULL,
                      double int_pnt2[3] = NULL );
   // Determines shortest distance between two infinite lines.  Optional
   // output also includes intersection points.  The num_ints will be zero
   // if the lines are parallel, one if they physically cross and two if
   // they are non-parallel and non-crossing.
   double dist_pnt_pln( double pnt[3], double pln_orig[3], double pln_norm[3], 
                        AgtSide *side = NULL, double pln_int[3] = NULL );
   // This routine determines the shortest distance between a point and a
   // plane.  This is the perpendicular distance from the point to the plane.
   // Optionally you can find which side of the plane the point is on and 
   // the plane intersection with the point.
#endif
   double dist_pln_pln( double pln_1_orig[3], double pln_1_norm[3], 
                        double pln_2_orig[3], double pln_2_norm[3],
                        AgtSide *side = NULL, AgtOrientation *orien = NULL,
                        unsigned short *status = NULL );
   // This routine determines the shortest distance between two planes.  This 
   // is the perpendicular distance between the two planes.  Note that if the
   // two planes are not parallel, the returned distance is zero.  The routine
   // can also be used to determine which side of the first plane the second
   // plane is on, and the orientation of the two planes to each other.
#ifdef BOYD14
   double dist_ln_pln( double ln_orig[3], double ln_vec[3], double pln_orig[3],
                       double pln_norm[3], AgtSide *side = NULL, 
                       unsigned short *status = NULL );
   // This routine determines the shortest distance between a line and a plane.  
   // This is the perpendicular distance between the line and the plane.  Note 
   // that if the line is not parallel to the plane, the returned distance is 
   // zero.  You can also find which side of the plane the line is on and
   // the status (= 0 if line & plane not parallel).
   double dist_pnt_ln( double pnt[3], double ln_orig[3], double ln_vec[3],
                       double int_pnt[3] = NULL );
   // This routine determines the shortest distance between a point and a line.  
   // This is the perpendicular distance between the point and the line.  You
   // can optionally get the intersection point of the point and line.
#endif
   
   //***************************************************************************
   // Intersections
   //***************************************************************************
   CubitStatus int_ln_pln( double ln_orig[3], double ln_vec[3], double pln_orig[3],
                           double pln_norm[3], double int_pnt[3] );
   // This routine calculates the intersection point of a line and a plane.
   // Returns CUBIT_FAILURE if no intersection exists (line is parallel to the
   // plane).
   int int_ln_ln( double p1[3], double v1[3], double p2[3], double v2[3], 
                  double int_pnt1[3], double int_pnt2[3] );
   // This function finds the intersection of two lines.  If the lines cross,
   // it finds the one unique point.  If the lines do not cross (but are non-
   // parallel), it finds the nearest intersection points (a line drawn from 
   // each of these intersections would be perpendicular to both lines).
   // Returns number of intersections found:
   //   0 intersections found, lines are parallel
   //   1 intersection found, lines physically intersect
	//   2 intersections found, lines do not intersect but nearest 
	//     intersections found
   int int_pnt_ln( double pnt[3], double ln_orig[3],  double ln_vec[3],
                            double int_pnt[3] );
   // This routine finds the perpendicular intersection of a point & line.  This
   // point will lie on the line.
   // Returns 0 if point is not on the line, otherwise 1.
   int int_pnt_pln( double pnt[3], double pln_orig[3], double pln_norm[3], 
                             double pln_int[3] );
   // This routine determines the perpendicular intersection of a point and a 
   // plane.  This is the perpendicular intersection point on the plane.
   // Returns 0 if point is not on the plane, otherwise 1.  
   CubitStatus int_pln_pln( double pln_1_orig[3], double pln_1_norm[3],
                            double pln_2_orig[3], double pln_2_norm[3],
                            double ln_orig[3], double ln_vec[3] );
   // This routine finds intersection of two planes.  This intersection results
   // in a line.  Returns CUBIT_FAILURE if planes are parallel.
#ifdef BOYD14
   CubitStatus int_pln_pln_pln( double pln_1_orig[3], double pln_1_norm[3],
                                double pln_2_orig[3], double pln_2_norm[3],
                                double pln_3_orig[3], double pln_3_norm[3],
                                double int_pnt[3] );
   // This routine finds the intersection of three orthoganal planes.  This 
   // intersection is a point.  Returns CUBIT_FAILURE if no point intersection
   // exists.
#endif
   
   //**************************************************************************
   // Comparison/Containment Tests
   //**************************************************************************
   CubitBoolean is_vec_par( double vec_1[3], double vec_2[3] );
   // This routine checks to see if two vectors are parallel.  Two vectors
   // are considered parallel if they are pointing in the same direction or
   // opposite directions.
   CubitBoolean is_vec_perp( double vec_1[3],double vec_2[3]);
   // This routine checks to see if two vectors are perpendicular.  Two vectors
   // are considered perpendicular if they are PI/2 radians to each other.
   CubitBoolean is_vecs_same_dir( double vec_1[3], double vec_2[3] );
   // This routine checks to see if two vectors point in the same direction.
   // The vector magnitudes do not have to be equal for a successful return.
   CubitBoolean is_pnt_on_ln( double pnt[3], double ln_orig[3], double ln_vec[3] );
   // This routined determines if a specified point is on a given infinite line
   // defined by a point and a vector.
   CubitBoolean is_pnt_on_ln_seg( double pnt1[3], double end1[3], double end2[3] );
   // This routine determines if a specified point is on a given line segment
   // defined by two endpoints.  The point is on the line only if it lies on
   // the line between or on the two endpoints.
   CubitBoolean is_pnt_on_pln( double pnt[3], double pln_orig[3], double pln_norm[3] );
   // This routined determines if a specified point is on a given plane
   // which is defined by an origin and three vectors.
   CubitBoolean is_ln_on_pln( double ln_orig[3], double ln_vec[3],
                              double pln_orig[3], double pln_norm[3] );
   // This routine determines if a specified line (defined by a point and
   // a vector) is in a given plane (defined by the two vectors in the plane 
   // and the normal to that plane).
#ifdef BOYD14
   CubitBoolean do_rays_converge( double ray_1_orig[3], double ray_1_vec[3],
                                  double ray_2_orig[3], double ray_2_vec[3] );
   // This routine checks to see if two rays converge or diverge towards each
   // other.  A ray is a line that is not infinite (consists of an origin and
   // a vector direction).  Rays are defined to converge if they intersect
   // within a common plane to each ray (ie., if the rays do not lie in a 
   // common plane, one of the rays is temporarily translated into the common
   // plane for checking).  If the rays are parallel, no one common plane exists,
   // so they are checked as-is.  See the following examples for a clearer 
   // picture of how this routine works.  Assume that the rays in the first row
   // are in planes that are parallel to this document but not in planes that are
   // necessarily coincident.  Assume that the rays in the second row are in the
   // same plane.  Note that rays may converge even if they are not in the same 
   // plane.
   //    
   //    Examples:
   //    
   //      *---->          <----*           <----*                      
   //              ^               ^             ^            ^
   //              |               |             |            |      
   //              |               |             |            |
   //              *               *             *            *---->
   //	    
   //      converge        diverge           converge        converge   
   //      
   //      
   //      *--->  *--->    *--->  <---*    *--->            *--->    
   //                                            <---*         *--->
   //        converge        converge        diverge         diverge 
   CubitBoolean is_ln_on_ln( double ln_orig1[3], double ln_vec1[3],
                             double ln_orig2[3], double ln_vec2[3] );
   // This routine checks to see if two lines are colinear.  Two lines are 
   // considered colinear if they are parallel and the origin of one line is on 
   // the other line.
#endif
   CubitBoolean is_pln_on_pln( double pln_orig1[3], double pln_norm1[3],
                               double pln_orig2[3], double pln_norm2[3] );
   // This routine checks to see if two infinite planes are coincident.

   //**************************************************************************
   // Arcs/Circles
   //**************************************************************************
   void setup_arc( double vec_1[3], double vec_2[3], double origin[3],
                  double start_angle, double end_angle, // Angles in radians
                  double radius, AGT_3D_Arc &arc );
   void setup_arc( CubitVector& vec_1, CubitVector& vec_2, 
                   CubitVector origin, double start_angle, // Angles in radians
                   double end_angle, double radius,
                   AGT_3D_Arc &arc );
   // Functions to populate an arc structure.  The arc is defined with two
   // orthogonal vectors in a plane, the arc origin, the beginning angle
   // in radians to rotate from the start vector (towards second vector), the 
   // ending angle in radians to rotate from the start vector, and the radius 
   // of the arc.  The arc is parameterized from 0 to 1.
   void get_arc_xyz( AGT_3D_Arc &arc, double param, double pnt[3] );
   void get_arc_xyz( AGT_3D_Arc &arc, double param, CubitVector& pnt );
   // Given a parameter value from 0 to 1 on the arc, return the xyz location.
   // Arc is assumed to be parameterized from 0 to 1.

   int get_num_circle_tess_pnts( double radius, double len_tolerance = 1e-4 );
   // Get number of tessellation points on the circle required to 
   // approximate the length of the circle within len_tolerance.  Can
   // be used to find the number of tessellations points to display a 
   // circle - smaller circles will have fewer tessellation points, larger 
   // circles will have more tessellation points with the same tolerance.
   // Minimum number of tessellations points returned is 8.

   //**************************************************************************
   // Miscellaneous
   //**************************************************************************
   void get_pln_orig_norm( double A, double B, double C, double D, 
                           double pln_orig[3], double pln_norm[3] = NULL );
   // Finds an origin-normal format plane from reduced form
   //   A*x + B*y + C*z + D = 0
   void get_box_corners( double box_min[3],double box_max[3],double c[8][3] );
   // Find 8 corners of a box given minimum and maximum points.  Box is
   // defined starting from left-bottom-front clockwise (at front of box), 
   // then to the rear in same order.
   CubitStatus min_pln_box_int_corners( const CubitPlane& plane, 
                                        const CubitBox& box,
                                        int extension_type, double extension,
                                        CubitVector& p1, CubitVector& p2,
                                        CubitVector& p3, CubitVector& p4,
                                        CubitBoolean silent = CUBIT_FALSE );
   CubitStatus min_pln_box_int_corners( CubitVector& vec1,
                                        CubitVector& vec2,
                                        CubitVector& vec3,
                                        CubitVector& box_min, CubitVector& box_max,
                                        int extension_type, double extension,
                                        CubitVector& p1, CubitVector& p2, 
                                        CubitVector& p3, CubitVector& p4,
                                        CubitBoolean silent = CUBIT_FALSE );
   CubitStatus min_pln_box_int_corners( double plane_norm[3], double plane_coeff,
                                       double box_min[3], double box_max[3],
                                       int extension_type, double extension,
                                       double p1[3], double p2[3],   // OUT
                                       double p3[3], double p4[3],   // OUT
                                       CubitBoolean silent = CUBIT_FALSE); 
   // Finds the 4 corner points of the input infinite plane that just 
   // intersects the input box (defined by bottom-left and top-right corners, 
   // axis aligned with the cubit coordinate system) - plane should have
   // minimal area.  The resultant plane's corner points can be extended out 
   // by a percentage of diagonal or absolute (making it bigger than the minimal
   // area plane).  extension_type - 0=none, 1=percentage, 2=absolute
   // If silent is CUBIT_TRUE, no error messages are given.

   CubitStatus get_tight_bounding_box( DLIList<CubitVector*> &point_list,                                       
                                       CubitVector& p1, CubitVector& p2,
                                       CubitVector& p3, CubitVector& p4,
                                       CubitVector& p5, CubitVector& p6,
                                       CubitVector& p7, CubitVector& p8);
   CubitStatus get_tight_bounding_box( DLIList<CubitVector*> &point_list,                                              
                                       CubitVector &center,
                                       CubitVector axes[3],
                                       CubitVector &extension );
   // Finds the minimum size bounding box that fits around the points.  Box
   // will not necessarily be oriented in xyz directions.

   CubitStatus min_cyl_box_int( double radius, CubitVector& axis, CubitVector& center,
                                CubitBox& box,
                                int extension_type, double extension,
                                CubitVector &start, CubitVector &end,
                                int num_tess_pnts = 2048 );
   CubitStatus min_cyl_box_int( double radius, double axis[3], double center[3], 
                                double box_min[3], double box_max[3], 
                                int extension_type, double extension, 
                                double start[3], double end[3],
                                int num_tess_pnts = 2048 );
   // Finds the start and end of a cylinder that just intersects the input
   // box (defined by bottom-left and top-right corners, axis aligned with 
   // the cubit coordinate system).  The resultant cylinder's start and end
   // points can be extended out by a percentage of cylinder's length or
   // absolute (making it longer in both directions).
   //   extension_type - 0=none, 1=percentage, 2=absolute
   // The num_tess_pnts is the number of line points along the circle the 
   // function projects to the box to calculate the minimum intersection
   // (more points could give a more accurate result for non-axis aligned
   // cylinders).

   double MinTriangleTriangle (Triangle3& tri0, Triangle3& tri1, double& s, double& t, double& u, double& v);
   // Finds minimum distance between two triangles in 3D (MAGIC)
   
   double AreaIntersection (Triangle &tri1, Triangle &tri2 );
   // Finds area of intersection between two triangles (MAGIC)

   double Angle( Triangle3& tri1, Triangle3& tri2 );
   // Finds angle between normals of the two triangles (radians)

   double MinPointTriangle (const Point3& p, const Triangle3& tri, double& s, double& t);
   // Finds minimum distance beween a triange and point in 3D (MAGIC)

#ifdef BOYD14
   void Center( Triangle3& tri, double center[3] );
   // Finds center position of given triangle
#endif

   void Normal( Triangle3& tri, double normal[3] );
   // Finds normal vector of given triangle

   double ProjectedOverlap( Triangle3& tri1, Triangle3& tri2 );
   // Projects tri2 to the plane of tri1 and finds the overlap area

protected:
   AnalyticGeometryTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.
   
private:
   
   static AnalyticGeometryTool* instance_;

   double agtEpsilon;  // default = 1e-6


///////////////////////////////////////////////////////////////////////////////
//                               MAGIC SOFTWARE
//This code was obtained from Dave Eberly, at www.magic-software.com.  It
//has been modified only to be placed into AnalyticGeometryTool.  This was
//done for convenience only.  An alternate solution would be to create a 
//separate library for these functions.
//
//Steve Storm, storm@CAT.com, 3-27-99
///////////////////////////////////////////////////////////////////////////////
//MAGIC is an acronym for My Alternate Graphics and Image Code. The
//initial code base originated in 1991 as an attempt to answer questions that
//arise in my favorite news group, comp.graphics.algorithms, and has been
//growing ever since. Magic is intended to provide free source code for solving
//problems that commonly arise in computer graphics and image analysis. While
//the code at this web site is free, additional conditions are: 
//
// *  You may distribute the original source code to others at no charge. You
//    got it for free, so do not charge others for it. 
// *  You may modify the original source code and distribute it to others at no
//    charge. The modified code must be documented to indicate that it is not
//    part of the original package. I do not want folks to blame me for bugs
//    introduced by modifications. I do accept blame for bugs that are my
//    doing and will do my best to fix them. 
// *  You may use this code for non-commercial purposes. You may also
//    incorporate this code into commercial packages. However, you may not
//    sell any of your source code which contains my original and/or modified
//    source code (see items 1 and 2 in this list of conditions). In such a case,
//    you would need to factor out my code and freely distribute it. Send me
//    email if you use it. I am always interested in knowing where and how
//    my code is being used. 
// *  The original code comes with absolutely no warranty and no guarantee is
//    made that the code is bug-free. Caveat emptor. 
//
// Dave Eberly, eberly@magic-software.com, www.magic-software.com
///////////////////////////////////////////////////////////////////////////////
//FILE: minbox2.h
   double Area (int N, Point2* pt, double angle);
   void MinimalBoxForAngle (int N, Point2* pt, double angle, OBBox2& box);
   OBBox2 MinimalBox2 (int N, Point2* pt);
   // Functions to find minimal area rectangle that surrounds a set of points.

#if 1
//FILE: minbox3.h
   void MatrixToAngleAxis( double** R, double& angle, double axis[3] );
   void AngleAxisToMatrix( double angle, double axis[3], double R[3][3] );
   double Volume( int N, Point3* pt, double angle[3] );
   void MinimalBoxForAngles( int N, Point3* pt, double angle[3], OBBox3& box );
   void GetInterval( double A[3], double D[3], double& tmin, double& tmax );
   void Combine( double result[3], double A[3], double t, double D[3] );
   double MinimizeOnInterval( int N, Point3* pt, double A[3], double D[3] );
   double MinimizeOnLattice( int N, Point3* pt, double A[3], int layers,
                             double thickness );
   void InitialGuess( int N, Point3* pt, double angle[3] );
   OBBox3 MinimalBox3( int N, Point3* pt );
   // Functions to find minimal volume box that surrounds a set of points.
#endif

   double Abs (double x);
   double ACos (double x);
   double ATan2 (double y, double x);
   double Cos (double x);
   double Sign (double x);
   double Sin (double x);
   double Sqrt (double x);
   double UnitRandom ();
   double Dot (const Point2& p, const Point2& q);
#ifdef BOYD14
   Point2 Perp (const Point2& p);
#endif
   double Length (const Point2& p);
   CubitBoolean Unitize (Point2& p );
   double Dot (const Point3& p, const Point3& q);
   Point3 Cross (const Point3& p, const Point3& q);
   double Length (const Point3& p);
   CubitBoolean Unitize (Point3& p );
   // Supporting code for triangle calculations

#ifdef BOYD14
// FILE: lin3lin3.cpp
   double MinLineLine (const Line3& line0, const Line3& line1, double& s, double& t);
   double MinLineRay (const Line3& line, const Line3& ray, double& s, double& t);
   double MinLineLineSegment (const Line3& line, const Line3& seg, double& s, double& t);
   double MinRayRay (const Line3& ray0, const Line3& ray1, double& s, double& t);
   double MinRayLineSegment (const Line3& ray, const Line3& seg, double& s, double& t);
#endif
   double MinLineSegmentLineSegment (const Line3& seg0, const Line3& seg1, double& s, double& t);

   //get intersection between bounding box and plane
   int get_plane_bbox_intersections( double box_min[3],
                                     double box_max[3],
                                     double pln_orig[3],
                                     double pln_norm[3],
                                     double int_array[6][3] );
// FILE: pt3tri3.cpp


// FILE: lin3tri3.cpp
   double MinLineSegmentTriangle (const Line3& seg, const Triangle3& tri, double& r, double& s, double& t);
   // Code to support distance between triangles

//FILE: triasect.cpp
   AgtLine EdgeToLine (Point2* v0, Point2* v1);
   void TriangleLines (Triangle* T, AgtLine line[3]);
   AgtTriList* SplitAndDecompose (Triangle* T, AgtLine* line);
   AgtTriList* Intersection (Triangle* T0, Triangle* T1);
   double AreaTriangle (Triangle* T);
   double Area (AgtTriList* list);
   double Area (Triangle &tri1, Triangle &tri2 );
   // Code to support area of overlap of two triangles
};

#if 1
// FILE: eigen.h
class mgcEigen
{
public:
        mgcEigen (int _size);
        ~mgcEigen ();

        mgcEigen& Matrix (float** inmat);

#ifdef BOYD14
        // solve eigensystem
        void EigenStuff2 ();  // uses TriDiagonal2
        void EigenStuff3 ();  // uses TriDiagonal3
        void EigenStuff4 ();  // uses TriDiagonal4
        void EigenStuffN ();  // uses TriDiagonalN
        void EigenStuff  ();  // uses switch statement
#endif

#ifdef BOYD14
        // solve eigensystem, use decreasing sort on eigenvalues
        void DecrSortEigenStuff2 ();
        void DecrSortEigenStuff3 ();
        void DecrSortEigenStuff4 ();
        void DecrSortEigenStuffN ();
        void DecrSortEigenStuff  ();
#endif

#ifdef BOYD14
        // solve eigensystem, use increasing sort on eigenvalues
        void IncrSortEigenStuff2 ();
        void IncrSortEigenStuff3 ();
        void IncrSortEigenStuff4 ();
        void IncrSortEigenStuffN ();
        void IncrSortEigenStuff  ();
#endif

private:
        int size;
        float** mat;
        float* diag;
        float* subd;

        // Householder reduction to tridiagonal form
        void Tridiagonal2 (float** pmat, float* pdiag, float* psubd);
        void Tridiagonal3 (float** pmat, float* pdiag, float* psubd);
        void Tridiagonal4 (float** pmat, float* pdiag, float* psubd);
        void TridiagonalN (int n, float** mat, float* diag, float* subd);

        // QL algorithm with implicit shifting, applies to tridiagonal matrices
        void QLAlgorithm (int n, float* pdiag, float* psubd, float** pmat);

        // sort eigenvalues from largest to smallest
        void DecreasingSort (int n, float* eigval, float** eigvec);

        // sort eigenvalues from smallest to largest
        void IncreasingSort (int n, float* eigval, float** eigvec);

// error handling
public:
        static int verbose1;
        static unsigned error;
        static void Report (std::ostream& ostr);
private:
        static const unsigned invalid_size;
        static const unsigned allocation_failed;
        static const unsigned ql_exceeded;
        static const char* message[3];
        static int Number (unsigned single_error);
        static void Report (unsigned single_error);
};

// FILE: eigen.h
class mgcEigenD
{
public:
        mgcEigenD (int _size);
        ~mgcEigenD ();

    // set the matrix for eigensolving
        double& Matrix (int row, int col) { return mat[row][col]; }
        mgcEigenD& Matrix (double** inmat);

        // get the results of eigensolving
        double Eigenvector (int row, int col) { return mat[row][col]; }
        const double** Eigenvector () { return (const double**) mat; }

#ifdef BOYD14
        // solve eigensystem
        void EigenStuff2 ();  // uses TriDiagonal2
#endif
        void EigenStuff3 ();  // uses TriDiagonal3
#ifdef BOYD14
        void EigenStuff4 ();  // uses TriDiagonal4
        void EigenStuffN ();  // uses TriDiagonalN
        void EigenStuff  ();  // uses switch statement
#endif

#ifdef BOYD14
        // solve eigensystem, use decreasing sort on eigenvalues
        void DecrSortEigenStuff2 ();
        void DecrSortEigenStuff3 ();
        void DecrSortEigenStuff4 ();
        void DecrSortEigenStuffN ();
        void DecrSortEigenStuff  ();
#endif

#ifdef BOYD14
        // solve eigensystem, use increasing sort on eigenvalues
        void IncrSortEigenStuff2 ();
        void IncrSortEigenStuff3 ();
        void IncrSortEigenStuff4 ();
        void IncrSortEigenStuffN ();
        void IncrSortEigenStuff  ();
#endif

private:
        int size;
        double** mat;
        double* diag;
        double* subd;

        // Householder reduction to tridiagonal form
        void Tridiagonal2 (double** mat, double* diag, double* subd);
        void Tridiagonal3 (double** mat, double* diag, double* subd);
        void Tridiagonal4 (double** mat, double* diag, double* subd);
        void TridiagonalN (int n, double** mat, double* diag, double* subd);

        // QL algorithm with implicit shifting, applies to tridiagonal matrices
        void QLAlgorithm (int n, double* diag, double* subd, double** mat);

        // sort eigenvalues from largest to smallest
        void DecreasingSort (int n, double* eigval, double** eigvec);

        // sort eigenvalues from smallest to largest
        void IncreasingSort (int n, double* eigval, double** eigvec);

// error handling
public:
        static int verbose1;
        static unsigned error;
        static void Report (std::ostream& ostr);
private:
        static const unsigned invalid_size;
        static const unsigned allocation_failed;
        static const unsigned ql_exceeded;
        static const char* message[3];
        static int Number (unsigned single_error);
        static void Report (unsigned single_error);
};
#endif

///// MAGIC SOFTWARE - INLINE FUNCTIONS
//---------------------------------------------------------------------------
// Wrapped math functions
//---------------------------------------------------------------------------
inline double AnalyticGeometryTool::Abs (double x)
{
    return double(fabs(x));
}

inline double AnalyticGeometryTool::ATan2 (double y, double x)
{
    return double(atan2(y,x));
}

inline double AnalyticGeometryTool::Sqrt (double x)
{
    return double(sqrt(x));
}

inline double AnalyticGeometryTool::UnitRandom ()
{
    return double(rand())/double(RAND_MAX);
}

//---------------------------------------------------------------------------
// Points in 2D
//---------------------------------------------------------------------------

inline double AnalyticGeometryTool::Dot (const Point3& p, const Point3& q)
{
    return double(p.x*q.x + p.y*q.y + p.z*q.z);
}


inline double AnalyticGeometryTool::set_epsilon( double new_epsilon )
{ 
   double old_epsilon = agtEpsilon;
   agtEpsilon = new_epsilon;
   return old_epsilon;
}

inline CubitBoolean AnalyticGeometryTool::dbl_eq(double val_1,double val_2)
{ 
   CubitBoolean result;
   if (fabs(val_1 - val_2) < agtEpsilon)
      result = CUBIT_TRUE;
   else
      result = CUBIT_FALSE;
   return result;
}

#endif

