//- Class: CubitPlane
//- Description: This file defines the CubitPlane class.
//- Owner: Greg Sjaardema
//- Checked by:


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CubitPlane.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"

CubitPlane::CubitPlane(const double A, const double B,
                       const double C, const double D) : normal_(A, B, C)
{
  double length = normal_.length();
  d_ = D / length;
  normal_ /= length;
}

CubitPlane::CubitPlane() : normal_(0.0, 0.0, 0.0)
{
  d_ = 0.0;
}

CubitPlane::CubitPlane(const CubitVector &Normal, const double D) :
  normal_(Normal)
{
  double length = normal_.length();
  d_ = D / length;
  normal_ /= length;
}
	
CubitPlane::CubitPlane(const CubitVector &Normal, const CubitVector &point) :
  normal_(Normal)
{
// plane_D = -(plane_normal.x()*X + plane_normal.y()*Y + plane_normal.z()*Z)
  normal_.normalize();
  d_ = -1.0 * (normal_.x()*point.x() + normal_.y()*point.y() +
	       normal_.z()*point.z());
}

void CubitPlane::set(const CubitVector &Normal, const CubitVector &point)
{
  normal_ = Normal;
  normal_.normalize();
  d_ = -1.0 * (normal_.x()*point.x() + normal_.y()*point.y() +
	       normal_.z()*point.z());
}

CubitPlane::CubitPlane(DLIList<CubitVector*> &positions) 
{
// Newells method for determining approximate plane equation.
// Plane equation is of the form:
// 0 = plane_D +plane_normal.x()*X + plane_normal.y()*Y + plane_normal.z()*Z

  CubitVector vector_diff;
  CubitVector vector_sum;
  CubitVector ref_point = CubitVector (0.0, 0.0, 0.0);
  CubitVector tmp_vector;
  normal_ = CubitVector(0.0, 0.0, 0.0);
  CubitVector *vector1, *vector2;
  
  for (int i=positions.size(); i > 0; i--)
  {
    vector1 = positions.get_and_step();
    vector2 = positions.get();
    vector_diff = (*vector2) - (*vector1);
    ref_point += (*vector1);
    vector_sum = (*vector1) + (*vector2);
    
    tmp_vector.set(vector_diff.y() * vector_sum.z(),
                   vector_diff.z() * vector_sum.x(),
                   vector_diff.x() * vector_sum.y());
    normal_ += tmp_vector;
  }
  double divisor = positions.size() * normal_.length();
  d_ = (ref_point % normal_) / divisor;
  normal_.normalize();
  normal_ *= -1.0;
}

CubitPlane::CubitPlane(const CubitPlane& copy_from)
{
  normal_ = copy_from.normal_;
  d_ = copy_from.d_;
}

int CubitPlane::mk_plane_with_points( const CubitVector& V0, 
                                      const CubitVector& V1,
                                      const CubitVector& V2 )
{
    // First check to make sure that the points are not collinear
  
    // Vector going from Vertex 0 to Vertex 1	    
  CubitVector vector01 = V1 - V0;
  vector01.normalize();
  
    // Vector going from Vertex 0 to Vertex 2	  
  CubitVector vector02 = V2 - V0;
  vector02.normalize();
  
    // If the 3 points are collinear, then the cross product of these two
    // vectors will yield a null vector (one whose length is zero).
  normal_ = vector01 * vector02;
  
  if (normal_.length() < CUBIT_RESABS)
  {
    PRINT_ERROR("Points are collinear.\n"
                "       Cannot create a CubitPlane object.\n");
    d_ = 0.0;
    return CUBIT_FAILURE;
  }
  
  normal_.normalize();
  
  double D0 =   -(normal_ % V0);
  double D1 =   -(normal_ % V1);
  double D2 =   -(normal_ % V2);
  d_ = (D0 + D1 + D2) / 3.0;
  
  return CUBIT_SUCCESS;
}

CubitVector CubitPlane::point_on_plane() const
{
  if (normal_.x() != 0)
    return CubitVector(-d_ / normal_.x(), 0, 0);
  else if (normal_.y() != 0)
    return CubitVector(0, -d_ / normal_.y(), 0);
  else if (normal_.z() != 0)
    return CubitVector(0, 0, -d_ / normal_.z());
    // If A B and C are all zero, the plane is invalid,
    // Just return <0,0,0>
  return CubitVector(0,0,0);
}

//-  Plane assignment operator.
CubitPlane& CubitPlane::operator=(const CubitPlane &plane)
{
  if (this != &plane)
  {
    normal_ = plane.normal_;
    d_ = plane.d_;
  }
  return *this;
}

double CubitPlane::distance(const CubitVector &vector) const
{
  return normal_ % vector + d_;
}

CubitVector CubitPlane::intersect(const CubitVector &base,
                                  const CubitVector &direction) const
{
  double t = -(normal_ % base + d_) / (normal_ % direction);
  return (base + direction * t);
}

int CubitPlane::intersect(const CubitPlane &plane_2,
                          CubitVector &origin, CubitVector &vector) const
{
    // Code from Graphics Gems III, Georgiades
    // Calculate the line of intersection between two planes. 
    // Initialize the unit direction vector of the line of intersection in
    // xdir.
    // Pick the point on the line of intersection on the coordinate plane most
    // normal to xdir.
    // Return TRUE if successful, FALSE otherwise (indicating that the planes
    // don't intersect). The order in which the planes are given determines the
    // choice of direction of xdir.
    //
    // int GetXLine(vect4 *pl1, vect4 *plane_2, vect3 *vector, vect3 *xpt)
  double invdet;  // inverse of 2x2 matrix determinant
  vector = normal() * plane_2.normal();
  CubitVector dir2(vector.x()*vector.x(),
                   vector.y()*vector.y(),
                   vector.z()*vector.z());
  CubitVector plane1n = normal();
  CubitVector plane2n = plane_2.normal();
  
  if (dir2.z() > dir2.y() && dir2.z() > dir2.x() && dir2.z() > CUBIT_RESABS)
  {
      // then get a point on the XY plane
    invdet = 1.0 / vector.z();
    
      //solve < pl1.x * origin.x + pl1.y * origin.y = - pl1.w >
      //      < plane2n.x * origin.x + plane2n.y * origin.y = - plane2n.w >
    origin.set(plane1n.y() * plane_2.coefficient() -
               plane2n.y() * coefficient(),
               plane2n.x() * coefficient() -
               plane1n.x() * plane_2.coefficient(),
               0.0);
  }
  else if (dir2.y() > dir2.x() && dir2.y() > CUBIT_RESABS)
  {
      // then get a point on the XZ plane
    invdet = -1.0 / vector.y();
    
      // solve < pl1.x * origin.x + pl1.z * origin.z = -pl1.w >
      //       < plane2n.x * origin.x + plane2n.z * origin.z = -plane2n.w >
    origin.set(plane1n.z() * plane_2.coefficient() -
               plane2n.z() * coefficient(),
               0.0,
               plane2n.x() * coefficient() -
               plane1n.x() * plane_2.coefficient());
  }
  else if (dir2.x() > CUBIT_RESABS)
  {
      // then get a point on the YZ plane
    invdet = 1.0 / vector.x();
    
      // solve < pl1.y * origin.y + pl1.z * origin.z = - pl1.w >
      //             < plane2n.y * origin.y + plane2n.z * origin.z = - plane2n.w >
    origin.set(0.0,
               plane1n.z() * plane_2.coefficient() -
               plane2n.z() * coefficient(),
               plane2n.y() * coefficient() -
               plane1n.y() * plane_2.coefficient());
  }
  else // vector is zero, then no point of intersection exists
    return CUBIT_FALSE;
  
  origin *= invdet;
  invdet = 1.0 / (float)sqrt(dir2.x() + dir2.y() + dir2.z());
  vector *= invdet;
  return CUBIT_TRUE;
}

CubitVector CubitPlane::project( const CubitVector& point ) const
{
  return point - distance(point) * normal();
}
