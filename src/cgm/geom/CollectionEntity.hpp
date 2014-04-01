//-------------------------------------------------------------------------
// Filename      : CollectionEntity.hpp
//
// Purpose       : This class is the base class for all the other solid 
//                 model entities in the Geometry subsystem. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

#ifndef COLLECTION_ENTITY_HPP
#define COLLECTION_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include "CubitString.hpp"
#include "CubitGeomConfigure.h"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class GeometryQueryEngine;
class RefCollection;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT CollectionEntity
{
public:
  
  CollectionEntity();
    // Default constructor.
  
  virtual ~CollectionEntity() ; 
    // Destructor.
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const = 0;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual const type_info& entity_type_info() const
    { return typeid(CollectionEntity); }
    //- return the type for this CollectionEntity
  
  virtual CubitString name() const = 0;
  virtual void name(CubitString) = 0;
    // Collection Entities aren't RefEntities, but they may have names.
    // Note: This is just here to keep the initial Pro/E
    //       implementation happy.

  CubitStatus set_owner_attribute(RefCollection* owner)
    {
      if (myOwner)
        return CUBIT_FAILURE;
      myOwner = owner;
      return CUBIT_SUCCESS;
    }
  RefCollection* get_owner_attribute() const
    { return myOwner; }
  
protected:
  RefCollection* myOwner;
private:
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

