//- Class: Cubit2DPoint
//-
//- Description: This file defines the Cubit2DPoint class which is a
//- standard two-dimensional point in space. 
//-
//- Owner: Steve Storm

#ifndef CUBIT2DPOINT_HPP
#define CUBIT2DPOINT_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

template <class X> class DLIList;

class CUBIT_UTIL_EXPORT Cubit2DPoint
{
public:
  
    //- Heading: Constructors and Destructor
  Cubit2DPoint();  //- Default constructor.
  
  Cubit2DPoint( const double x, const double y );
    //- Constructor: create point from two components

  Cubit2DPoint( const double xy[2] );
    //- Constructor: create point from array

  Cubit2DPoint( const Cubit2DPoint& copy_from );  //- Copy Constructor
  Cubit2DPoint( const Cubit2DPoint* copy_from );  //- Copy Constructor
  
    //- Heading: Set and Inquire Functions
  void set( const double x, const double y );
    //- Change point components to {x}, {y}
  
  void set( const double xy[2] );
    //- Change point components to xy[0], xy[1]
  
  void set( const Cubit2DPoint& to_copy );
    //- Same as operator=(const Cubit2DPoint&)
  
  double x() const; //- Return x component of point
  double y() const; //- Return y component of point
  
  void get_xy( double &x, double &y ); //- Get x, y components
  void get_xy( double xy[3] ); //- Get xy array
  
  void x( const double x ); //- Set x component of point
  void y( const double y ); //- Set y component of point

#ifdef BOYD15
  double distance_between(const Cubit2DPoint& test_point);
    //- Calculate the distance between this point and another.
#endif

  void print_me();
    //- Prints out the coordinates of this point.

#ifdef BOYD15
  CubitBoolean within_tolerance(const Cubit2DPoint &pnt2,
                                double tol = 1e-10) const;
    //- Compare two points to see if they are spatially equal
    //- within the tolerance given.
#endif

  void min_max( const Cubit2DPoint &pnt2,
                double &xmin, double &xmax,
                double &ymin, double &ymax ) const;
    //- Get min & max coordinates between this point & the
    //- passed-in point.  Analogous to a 2D bounding box.

  void update_min_max( Cubit2DPoint &min, Cubit2DPoint &max ) const;
    //- Include this points coordinates in the given min-max bounding box.

  CubitBoolean is_on_line_segment( const Cubit2DPoint &end1,
                                   const Cubit2DPoint &end2,
                                   double tol = 1e-10 ) const;
    //- Determine if this point is on the line segment defined by end1 
    //- and end2.  Uses a triangle method for speed.  Note that in this 
    //- implementation the tolerance is not how close the point is to 
    //- the line - rather, it is how small a triangle area needs to be
    //- in order to consider its three points co-linear.

  Cubit2DPoint &operator=(const Cubit2DPoint &from);

private:
  
  double xVal;  //- x component of point.
  double yVal;  //- y component of point.
};

inline Cubit2DPoint::Cubit2DPoint(const Cubit2DPoint& copy_from)
    : xVal(copy_from.xVal), yVal(copy_from.yVal)
{}

inline Cubit2DPoint::Cubit2DPoint(const Cubit2DPoint* copy_from)
    : xVal(copy_from->xVal), yVal(copy_from->yVal)
{}

inline Cubit2DPoint::Cubit2DPoint()
    : xVal(0.0), yVal(0.0)
{}

inline Cubit2DPoint::Cubit2DPoint( const double x,
                                   const double y )
    : xVal(x), yVal(y)
{}

inline double Cubit2DPoint::x() const
{ return xVal; }
inline double Cubit2DPoint::y() const
{ return yVal; }

inline void Cubit2DPoint::x( const double x )
{ xVal = x; }
inline void Cubit2DPoint::y( const double y )
{ yVal = y; }

inline void Cubit2DPoint::set( const double x,
                               const double y )
{
  xVal = x;
  yVal = y;
}

inline void Cubit2DPoint::set( const double xy[2] )
{
  xVal = xy[0];
  yVal = xy[1];
}

inline Cubit2DPoint& Cubit2DPoint::operator=(const Cubit2DPoint &from)  
{
  xVal = from.xVal; 
  yVal = from.yVal;
  return *this;
}

inline void Cubit2DPoint::set(const Cubit2DPoint& to_copy)
{
  *this = to_copy;
}

#endif


