//-------------------------------------------------------------------------
// Filename      : Point.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef POINT_HPP
#define POINT_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "GeometryEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT Point : public GeometryEntity
{
public :
  Point() ;
    //- The default constructor
  
  virtual ~Point() ;
    //- The destructor
  
  virtual const type_info& topology_entity_type_info() const;
  
  virtual CubitVector coordinates() const = 0;  // pure virtual
    //R CubitVector
    //R- Contains the coordinate values {x y z} of this Point
    //- Returns the spatial coordinates of this Point.
  
  virtual double measure() 
    {return 1.0;}
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
  
  virtual GeometryType geometry_type()
    {return UNDEFINED_POINT_TYPE;}
    //R GeometryType
    //R- the enumerated type of the geometric represnetation.

  virtual CubitBox bounding_box() const;
  
protected: 
  
private:
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

