//- Class: CubitBox
//-
//- Description: This file defines the CubitBox class which represents
//- an axis-aligned rectangular box which can be uses as a bounding box.
//-
//- Owner: Greg Sjaardema
//- Checked by: 
//- Version: $Id: 

#ifndef CUBITBOX_HPP
#define CUBITBOX_HPP

#include "CubitVector.hpp"
#include "CubitBoxStruct.h"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitBox
{
public:
    //- Heading: Constructors and Destructor
  CubitBox(); 
    //- Default constructor.
  
  CubitBox(const CubitVector &min, const CubitVector &max);
    //- Constructor: create box from two CubitVectors
  
  CubitBox(const double min[3], const double max[3] );
    //- Constructor: create box from two coordinates
  
  CubitBox(const CubitVector &min_max);
    //- Constructor: create box from one CubitVector
  
  CubitBox(const CubitBox& copy_from);  //- Copy Constructor
  
  CubitBox(const CubitBoxStruct& from);
  
  ~CubitBox(); 
    //- destructor
    
  CubitBox& bounding_box();
  
  
  void reset(const CubitVector &vector);
  void reset(const CubitVector &min, const CubitVector &max);
  void reset(const CubitBox &box);
  void reset(const double min[3], const double max[3]);
    //- reset ranges

  double max_x() const;
  double max_y() const;
  double max_z() const;

  double min_x() const;
  double min_y() const;
  double min_z() const;

  CubitVector minimum()  const;
  CubitVector maximum()  const;
  CubitVector center()   const;
  CubitVector diagonal() const;
    //- Return Box minimum/maximum/center
  
  void get_corners ( CubitVector corners[8] ) const;
    //- Fills 'corners' with the corners of this box.
    //- The order is:
    //-   0) minimum()
    //-   1-3) Front face (Constant minimum z-plane), normal out of box
    //-        using right hand rule, including corner[0].
    //-   4-7) Same as 0-3, but offset to back face
    //-        (constant maximum z-plane).  Normal of these last 4 points
    //-        is into box relative to back plane (same direction as
    //-        normal w/ first 4 points).  Maximum ends up at index 6.
  
  double x_range() const;
  double y_range() const;
  double z_range() const;
    //- x, y, and z range of the box (max - min)

  bool overlap( double tolerance, const CubitBox& other_box ) const;
    //- Check if boxes are within passed tolerance of each other.
    //- If tolerance is 0, use && or || operator.
  
    //- Heading: Operators

    // Operators that modify {this}
  CubitBox& operator=(const CubitBox &box);
  CubitBox& operator|=(const CubitBox& box);
  CubitBox& operator|=(const CubitVector& vector);
  CubitBox& operator&=(const CubitBox& box);
  CubitBox& operator*=(double scale);
  CubitBox& operator/=(double scale);
  CubitBox& operator+=(const CubitVector& offset);
  CubitBox& operator-=(const CubitVector& offset);
    //- {=}  - Assignment
    //- {|=} - Union of {this} and {box}
    //- {&=} - Intersection (overlap) of {this} and {box}
    //- {*=} - Scale {this} about box center
    //- {/=} - Scale {this} about box center
    //- {+=} - Move {this} by {offset}  CubitVector
    //- {-=} - Move {this} by {-offset} CubitVector
  
    // Operators that check for containment
  int operator< (const CubitBox& box) const;
  int operator<=(const CubitBox& box) const;
  int operator> (const CubitBox& box) const;
  int operator>=(const CubitBox& box) const;
  int operator> (const CubitVector& vect) const;
  int operator>=(const CubitVector& vect) const;
  int operator<=(const CubitVector& vect) const;
  int operator&&(const CubitBox& box) const;
  int operator||(const CubitBox& box) const;
    //- {<}  - Is {this} completely surrounded by {box}?
    //- {>}  - Does {this} completely surround {box}?
    //- {<=,>=} - As above, but inner box may touch
    //-           boundary of outer box.
    //- {>}  - Is {vect} contained within {this}, but not on boundary?
    //- {>=} - Is {vect} contained within or on the boundary of {this}?
    //- {<=} - Is {vect} outside or on boundary of {this}?
    //- {&&} - Do {this} and {box} intersect?  Just butting against each
    //-        other also counts as an intersection. See {||}.
    //- {||} - Do {this} and {box} intersect?  Just butting against each
    //-        other does NOT count as an intersection.  See {&&}.
  
  CubitBox &operator=(const CubitBoxStruct &from);

  operator CubitBoxStruct() 
    {
      CubitBoxStruct to;
      to.minimum_ = minimum_;
      to.maximum_ = maximum_;
      return to;
    }

    // Operators that return a modification of {this}.
    // {this} itself is not modified.
  friend CubitBox operator|(const CubitBox& lhs, const CubitBox& rhs);
  friend CubitBox operator|(const CubitBox& lhs, const CubitVector& rhs);
  friend CubitBox operator&(const CubitBox& lhs, const CubitBox& rhs);
  friend CubitBox operator*(const CubitBox& lhs, double rhs);
  friend CubitBox operator*(double rhs, const CubitBox& lhs);
  friend CubitBox operator/(const CubitBox& lhs, double rhs);
  friend CubitBox operator+(const CubitBox& lhs, const CubitVector& rhs);
  friend CubitBox operator-(const CubitBox& lhs, const CubitVector& rhs);
  
#ifdef BOYD15
  double distance( const CubitVector& position ) const;
    //R double
    //R- The shortest distance from the passed position to the
    //R- box, if the point is outside the box.  Zero if the
    //R- passed position is within the box or on its boundary.
    //I position
    //I- A position from which to evaluate the shortest distance
    //I- to the box.
#endif
	
	double distance_squared( const CubitVector& position ) const;
  
  CubitVector closest_point( const CubitVector& position ) const;
    //R CubitVector
    //R- The closest point on the box to the passed position.
    //R- The passed position will be returned if it is within
    //R- the box.
    //I- A position from which to evaluate the closest point
    //I- on the box.
  
#ifdef BOYD15
  CubitBoolean within_2_dim( CubitBox& box ) const;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I box
    //I- The box to test this relative to.
    //- This method tests if this box is within the passed box
    //- in at least two dimensions (i.e. within the projection of
    //- box on the xz-plane, the xy-plane, or the yz-plane.)
#endif
  
private:
  
  CubitVector minimum_; //- X, Y, and Z position of minimum corner
  CubitVector maximum_; //- X, Y, and Z position of maximum corner
};

inline CubitBox& CubitBox::operator=(const CubitBoxStruct &from)  
{
  minimum_ = from.minimum_;
  maximum_ = from.maximum_;
  return *this;
}

inline CubitBox::CubitBox(const CubitBoxStruct &from)  
{
  minimum_ = from.minimum_;
  maximum_ = from.maximum_;
}
#endif

