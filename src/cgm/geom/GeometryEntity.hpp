//-------------------------------------------------------------------------
// Filename      : GeometryEntity.hpp
//
// Purpose       : This file contains the declarations of the base class 
//                 GeometryEntity.  This is the root of the class hierarchy
//                 the encapsulates all geometry information attached to
//                 topological entities.  The topological entities are
//                 independent of their underlying geometric representation.
//
// Special Notes : 
//
// Creator       : Xuechen Liu 
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef GEOMETRY_ENTITY_HPP
#define GEOMETRY_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "TopologyBridge.hpp"
#include "GeometryDefines.h"
#include "CubitBox.hpp"
#include "CubitObservable.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity ;
template <class X> class DLIList;
class GeometryQueryEngine ;
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN ENUM DEFINITIONS       **********
// ********** END ENUM DEFINITIONS         **********

class CUBIT_GEOM_EXPORT GeometryEntity : public TopologyBridge
{
public :
  
  GeometryEntity() : myId(0) {}
    //- The default constructor
  
  virtual ~GeometryEntity();
    //- The destructor
    //- The destructor is made pure virtual to prevent instantiation
    //- of this class.
  
  virtual CubitBox bounding_box() const = 0;
  
  virtual double measure() = 0;
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
  
  virtual GeometryType geometry_type() = 0;
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
    //- Returns the type of geometric representation, e.g.
    //- STRAIGHT_CURVE_TYPE, SPHERE_SURFACE_TYPE, etc

  int get_saved_id() const { return myId; }
  void set_saved_id( int value ) { myId = value; }
  
protected:
private:

  int myId;
  
};

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
#endif

