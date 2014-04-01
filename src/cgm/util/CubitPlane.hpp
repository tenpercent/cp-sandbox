//- Class: CubitPlane
//-
//- Description: This file defines the CubitPlane class which is a
//- three-dimensional planar surface defined by the equation  
//- {Ax + By + Cz + D = 0}. Plane normal is normalized
//- Only basic functionality required by other Cubit classes is 
//- currently supported.
//-
//- Owner: Greg Sjaardema
//- Checked by: Tony Edwards 8/18/94
//- Version: $Id: 

#ifndef CUBITPLANE_HPP
#define CUBITPLANE_HPP

template <class X> class DLIList;
#include "CubitVector.hpp"
#include "CubitPlaneStruct.h"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitPlane
{
public:
  
    //- Heading: Constructors and Destructor
  CubitPlane(); 
    //- Default constructor.
  
  CubitPlane(const double A, const double B,
             const double C, const double D);
    //- Constructor: create plane from three components of normal and 
    //- coefficient
  
  CubitPlane(const CubitVector &Normal, const double D);
    //- Constructor: create plane from plane normal and coefficient.
  
  CubitPlane(const CubitVector &Normal, const CubitVector &point);
    //- Constructor: create plane from plane normal that passes through point.
  
  CubitPlane(DLIList<CubitVector*> &positions);
    //- Constructor: create plane closest to the set of points in {positions} 
    //- using Newell's Method.
  
  CubitPlane(const CubitPlane& copy_from);  //- Copy Constructor
  
  CubitPlane(const CubitPlaneStruct& from);  
  
  int mk_plane_with_points(const CubitVector& vector1,
                           const CubitVector& vector2,
                           const CubitVector& vector3);
    //- Create a plane given three points represented as 
    //- CubitVectors.
    //- Return a CUBIT_FAILURE if the points are collinear.
  
  const CubitVector& normal() const;
    //- Return Plane normal (normalized)

  void normal(const CubitVector &temp_normal);
    //- set the normal for this plane
  
  double coefficient() const;
    //- Return the coefficient
  
  void coefficient(const double temp_coeff);
    //- set the coefficient

  void set(const CubitVector &Normal, const CubitVector &point);
    //- redefine plane using normal and point on plane
  
    //- Heading: Other Functions
  CubitVector point_on_plane() const;
    //- Returns a random point on the plane.
  
#ifdef BOYD15
  double report_plane_error(DLIList<CubitVector*> &points);
    //- Reports the average squared distance of the points
    //- to the plane. Deviation is calculated as  {sum(diff^2)/points.size()}
    //- of {points}
#endif
  
  double distance(const CubitVector &vector) const;
    //- Calculates the distance from {vector} to plane. If the point lies
    //- behind the plane (the opposite direction than the normal points),
    //- the returned distance will be negative.
  
  CubitVector intersect(const CubitVector &base,
                        const CubitVector &direction) const;
    //- Calculate intersection of line from {base} in direction {direction}
    //- and plane {this}.
    //- Returns coordinates of intersection in a CubitVector.
  
  int intersect(const CubitPlane &other_plane,
                CubitVector &origin, CubitVector &vector) const;
    //- Calculate the intersection of {this} with {other_plane}
    //- Returns a point on the intersection line in {origin} and the
    //- direction of the intersection line in {vector}.
    //- Returns CUBIT_FALSE if the planes are coplanar
    
  CubitVector project( const CubitVector& point ) const;
    //- Project a point onto plane
  
  CubitPlane& operator=(const CubitPlane &plane);
    //- assignment
  
  CubitPlane &operator=(const CubitPlaneStruct &from);

  operator CubitPlaneStruct() 
    {
      CubitPlaneStruct to;
      to.normal_ = normal_;
      to.d_ = d_;
      return to;
    }

private:
  
  CubitVector normal_;  //- Normal to plane.
  double d_;            //- Coefficient
  
};

inline const CubitVector& CubitPlane::normal() const { return normal_; }

inline void CubitPlane::normal(const CubitVector &temp_normal)
{normal_ = temp_normal;}

inline double CubitPlane::coefficient() const { return d_; }
inline void CubitPlane::coefficient(const double temp_coeff)
{d_ = temp_coeff; }

inline CubitPlane::CubitPlane(const CubitPlaneStruct &from)  
{
  normal_ = from.normal_;
  d_ = from.d_;
}

#endif

