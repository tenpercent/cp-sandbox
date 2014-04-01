//- Class: CubitBox
//- Description: This file defines the CubitBox class.
//- Owner: Greg Sjaardema
//- Checked by:

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "CubitDefines.h"

CubitBox::~CubitBox()
{}

CubitBox::CubitBox() :
    minimum_(0.0, 0.0, 0.0),
    maximum_(0.0, 0.0, 0.0)
{}

CubitBox::CubitBox(const CubitVector &min, const CubitVector &max)
{
  minimum_.set (CUBIT_MIN (min.x(), max.x()),
                CUBIT_MIN (min.y(), max.y()),
                CUBIT_MIN (min.z(), max.z()));
  maximum_.set (CUBIT_MAX (min.x(), max.x()),
                CUBIT_MAX (min.y(), max.y()),
                CUBIT_MAX (min.z(), max.z()));
}

CubitBox::CubitBox(const CubitVector &min_max) :
    minimum_(min_max),
    maximum_(min_max)
{}

CubitBox::CubitBox(const CubitBox& copy_from):
    minimum_(copy_from.minimum_),
    maximum_(copy_from.maximum_)
{}

void CubitBox::reset(const CubitVector& vector)
{
  minimum_ = vector;
  maximum_ = vector;
}

void CubitBox::reset(const CubitVector &min, const CubitVector &max)
{
  minimum_.set (CUBIT_MIN (min.x(), max.x()),
		CUBIT_MIN (min.y(), max.y()),
		CUBIT_MIN (min.z(), max.z()));
  maximum_.set (CUBIT_MAX (min.x(), max.x()),
		CUBIT_MAX (min.y(), max.y()),
		CUBIT_MAX (min.z(), max.z()));
}

void CubitBox::reset(const CubitBox &box)
{
  minimum_ = box.minimum_;
  maximum_ = box.maximum_;
}

CubitVector CubitBox::minimum() const
{
  return CubitVector(minimum_);
}

CubitVector CubitBox::maximum() const
{
  return CubitVector(maximum_);
}

double CubitBox::max_x() const
{
  return maximum_.x();
}

double CubitBox::max_y() const
{
  return maximum_.y();
}

double CubitBox::max_z() const
{
  return maximum_.z();
}

double CubitBox::min_x() const
{
  return minimum_.x();
}

double CubitBox::min_y() const
{
  return minimum_.y();
}

double CubitBox::min_z() const
{
  return minimum_.z();
}

CubitVector CubitBox::center() const
{
  return CubitVector(minimum_ + maximum_) / 2.0;
}

CubitVector CubitBox::diagonal() const
{
  return CubitVector(maximum_ - minimum_);
}
  
double CubitBox::x_range() const
{ return (maximum_.x() - minimum_.x()); }

double CubitBox::y_range() const
{ return (maximum_.y() - minimum_.y()); }

double CubitBox::z_range() const
{ return (maximum_.z() - minimum_.z()); }

CubitBox& CubitBox::operator=(const CubitBox &box)
{
  if (this != &box)
    {
      minimum_ = box.minimum_;
      maximum_ = box.maximum_;
    }
  return *this;
}

bool CubitBox::overlap( double tolerance , const CubitBox& box ) const
{
  //     | - note the '!'.  This condition checks that the boxes
  //     |   do NOT intersect and negates the result. 
  return ! ( (box.minimum_.x() - maximum_.x() > tolerance) ||
             (box.minimum_.y() - maximum_.y() > tolerance) ||
             (box.minimum_.z() - maximum_.z() > tolerance) ||
             (minimum_.x() - box.maximum_.x() > tolerance) ||
             (minimum_.y() - box.maximum_.y() > tolerance) ||
             (minimum_.z() - box.maximum_.z() > tolerance) );
}
             
  
  

// Union of this and box
CubitBox& CubitBox::operator|=(const CubitBox& box)
{
  minimum_.x(CUBIT_MIN(minimum_.x(), box.minimum_.x()));
  minimum_.y(CUBIT_MIN(minimum_.y(), box.minimum_.y()));
  minimum_.z(CUBIT_MIN(minimum_.z(), box.minimum_.z()));

  maximum_.x(CUBIT_MAX(maximum_.x(), box.maximum_.x()));
  maximum_.y(CUBIT_MAX(maximum_.y(), box.maximum_.y()));
  maximum_.z(CUBIT_MAX(maximum_.z(), box.maximum_.z()));
  return *this;
}

// Union of this and vector
CubitBox& CubitBox::operator|=(const CubitVector& vector)
{
  minimum_.x(CUBIT_MIN(minimum_.x(), vector.x()));
  minimum_.y(CUBIT_MIN(minimum_.y(), vector.y()));
  minimum_.z(CUBIT_MIN(minimum_.z(), vector.z()));

  maximum_.x(CUBIT_MAX(maximum_.x(), vector.x()));
  maximum_.y(CUBIT_MAX(maximum_.y(), vector.y()));
  maximum_.z(CUBIT_MAX(maximum_.z(), vector.z()));
  return *this;
}

// Intersection of this and box
inline CubitBox& CubitBox::operator&=(const CubitBox& box)
{
  minimum_.x(CUBIT_MAX(minimum_.x(), box.minimum_.x()));
  minimum_.y(CUBIT_MAX(minimum_.y(), box.minimum_.y()));
  minimum_.z(CUBIT_MAX(minimum_.z(), box.minimum_.z()));

  maximum_.x(CUBIT_MIN(maximum_.x(), box.maximum_.x()));
  maximum_.y(CUBIT_MIN(maximum_.y(), box.maximum_.y()));
  maximum_.z(CUBIT_MIN(maximum_.z(), box.maximum_.z()));

  if (minimum_.x() > maximum_.x() ||
      minimum_.y() > maximum_.y() ||
      minimum_.z() > maximum_.z())
  {
    minimum_.set(0.0, 0.0, 0.0);
    maximum_.set(0.0, 0.0, 0.0);
  }
  return *this;
}

CubitBox& CubitBox::operator*=(double scale)
{
  CubitVector center_vec = center();
  *this -= center_vec;
  minimum_ *= scale;
  maximum_ *= scale;
  *this += center_vec;
  return *this;
}

CubitBox& CubitBox::operator/=(double scale)
{
  assert(scale != 0.0);
  *this *= 1/scale;
  return *this;
}

CubitBox& CubitBox::operator+=(const CubitVector& offset)
{
  minimum_ += offset;
  maximum_ += offset;
  return *this;
}

CubitBox& CubitBox::operator-=(const CubitVector& offset)
{
  minimum_ -= offset;
  maximum_ -= offset;
  return *this;
}

int CubitBox::operator<(const CubitBox& box) const
{
  return (box.minimum_.x() < minimum_.x() &&
          box.minimum_.y() < minimum_.y() &&
          box.minimum_.z() < minimum_.z() &&
          box.maximum_.x() > maximum_.x() &&
          box.maximum_.y() > maximum_.y() &&
          box.maximum_.z() > maximum_.z());
}

int CubitBox::operator<=(const CubitBox& box) const
{
  return (box.minimum_.x() <= minimum_.x() &&
          box.minimum_.y() <= minimum_.y() &&
          box.minimum_.z() <= minimum_.z() &&
          box.maximum_.x() >= maximum_.x() &&
          box.maximum_.y() >= maximum_.y() &&
          box.maximum_.z() >= maximum_.z());
}

int CubitBox::operator>(const CubitBox& box) const
{
  return (box.minimum_.x() > minimum_.x() &&
          box.minimum_.y() > minimum_.y() &&
          box.minimum_.z() > minimum_.z() &&
          box.maximum_.x() < maximum_.x() &&
          box.maximum_.y() < maximum_.y() &&
          box.maximum_.z() < maximum_.z());
}

int CubitBox::operator>=(const CubitBox& box) const
{
  return (box.minimum_.x() >= minimum_.x() &&
          box.minimum_.y() >= minimum_.y() &&
          box.minimum_.z() >= minimum_.z() &&
          box.maximum_.x() <= maximum_.x() &&
          box.maximum_.y() <= maximum_.y() &&
          box.maximum_.z() <= maximum_.z());
}

int CubitBox::operator>(const CubitVector& vect) const
{
  return (vect.x() > minimum_.x() && vect.x() < maximum_.x() &&
          vect.y() > minimum_.y() && vect.y() < maximum_.y() &&
          vect.z() > minimum_.z() && vect.z() < maximum_.z() );
}

int CubitBox::operator>=(const CubitVector& vect) const
{
  return (vect.x() >= minimum_.x() && 
          vect.x() <= maximum_.x() &&
          vect.y() >= minimum_.y() && 
          vect.y() <= maximum_.y() &&
          vect.z() >= minimum_.z() && 
          vect.z() <= maximum_.z() );
}

CubitBox operator|(const CubitBox& lhs, const CubitBox& rhs)
{
  return CubitBox(lhs) |= rhs;
}

CubitBox operator&(const CubitBox& lhs, const CubitBox& rhs)
{
  return CubitBox(lhs) &= rhs;
}

CubitBox operator*(double lhs, const CubitBox& rhs)
{
  return CubitBox(rhs) *= lhs;
}

CubitBox operator*(const CubitBox& lhs, double rhs)
{
  return CubitBox(lhs) *= rhs;
}

CubitBox operator/(const CubitBox& lhs, double rhs)
{
  return CubitBox(lhs) /= rhs;
}

CubitBox operator+(const CubitBox& lhs, const CubitVector& rhs)
{
  return CubitBox(lhs) += rhs;
}

CubitBox operator-(const CubitBox& lhs, const CubitVector& rhs)
{
  return CubitBox(lhs) -= rhs;
}

double CubitBox::distance_squared( const CubitVector& point ) const
{
  return (point - closest_point(point)).length_squared();
}



//-------------------------------------------------------------------------
// Purpose       : Find the closest point on the CubitBox
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/98
//-------------------------------------------------------------------------
CubitVector CubitBox::closest_point( const CubitVector& point ) const
{
  CubitVector result;
  
  if( point.x() < minimum_.x() )
    result.x( minimum_.x() );
  else if( point.x() > maximum_.x() )
    result.x( maximum_.x() );
  else
    result.x( point.x() );
  
  if( point.y() < minimum_.y() )
    result.y( minimum_.y() );
  else if( point.y() > maximum_.y() )
    result.y( maximum_.y() );
  else
    result.y( point.y() );
  
  if( point.z() < minimum_.z() )
    result.z( minimum_.z() );
  else if( point.z() > maximum_.z() )
    result.z( maximum_.z() );
  else
    result.z( point.z() );
  
  return result;
}

CubitBox::CubitBox(const double min[3], const double max[3])
{
  minimum_.set (CUBIT_MIN (min[0], max[0]),
                CUBIT_MIN (min[1], max[1]),
                CUBIT_MIN (min[2], max[2]));
  maximum_.set (CUBIT_MAX (min[0], max[0]),
                CUBIT_MAX (min[1], max[1]),
                CUBIT_MAX (min[2], max[2]));
}

void CubitBox::reset(const double min[3], const double max[3])
{
  minimum_.set (CUBIT_MIN (min[0], max[0]),
                CUBIT_MIN (min[1], max[1]),
                CUBIT_MIN (min[2], max[2]));
  maximum_.set (CUBIT_MAX (min[0], max[0]),
                CUBIT_MAX (min[1], max[1]),
                CUBIT_MAX (min[2], max[2]));
}

void CubitBox::get_corners( CubitVector vectors[8] ) const
{
  vectors[0] = minimum_;
  vectors[1] = CubitVector (maximum_.x(), minimum_.y(), minimum_.z());
  vectors[2] = CubitVector (maximum_.x(), maximum_.y(), minimum_.z());
  vectors[3] = CubitVector (minimum_.x(), maximum_.y(), minimum_.z());
  
  vectors[4] = CubitVector (minimum_.x(), minimum_.y(), maximum_.z());
  vectors[5] = CubitVector (maximum_.x(), minimum_.y(), maximum_.z());
  vectors[6] = maximum_;
  vectors[7] = CubitVector (minimum_.x(), maximum_.y(), maximum_.z());
}

int CubitBox::operator&&(const CubitBox& box) const 
{
    // Return false if there is no overlap
    // along at least one of the 3 axes.
  if (minimum_.x() > box.maximum_.x() ||
      maximum_.x() < box.minimum_.x() ||
      minimum_.y() > box.maximum_.y() ||
      maximum_.y() < box.minimum_.y() ||
      minimum_.z() > box.maximum_.z() ||
      maximum_.z() < box.minimum_.z() )
    return CUBIT_FALSE;
  
    // If you didn't return false...
  return CUBIT_TRUE;
}

int CubitBox::operator||(const CubitBox& box) const 
{
    // Return false if there is no overlap
    // along at least one of the 3 axes.
  if (minimum_.x() >= box.maximum_.x() ||
      maximum_.x() <= box.minimum_.x() ||
      minimum_.y() >= box.maximum_.y() ||
      maximum_.y() <= box.minimum_.y() ||
      minimum_.z() >= box.maximum_.z() ||
      maximum_.z() <= box.minimum_.z() )
    return CUBIT_FALSE;
  
    // If you didn't return false...
  return CUBIT_TRUE;
}

int CubitBox::operator<=(const CubitVector& vect) const
{
  return (vect.x() <= minimum_.x() ||
          vect.x() >= maximum_.x() ||
          vect.y() <= minimum_.y() ||
          vect.y() >= maximum_.y() ||
          vect.z() <= minimum_.z() ||
          vect.z() >= maximum_.z() );
}
