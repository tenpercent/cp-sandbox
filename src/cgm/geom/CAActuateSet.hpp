//-------------------------------------------------------------------------
// Filename      : CAActuateSet.hpp
//
// Purpose       : Maintain the list of entities for which attributes are
//                 being actuated such that any entities destroyed 
//                 during actuation (e.g. merging) get removed from the
//                 list.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------

#ifndef CA_ACTUATE_SET_HPP
#define CA_ACTUATE_SET_HPP

#include "DLIList.hpp"
#include "CubitObserver.hpp"
#include "DagType.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;

class CUBIT_GEOM_EXPORT CAActuateSet : public CubitObserver
{
  public:
  
    CAActuateSet( DLIList<RefEntity*>& actuate_list );
      //- Constructor
      //- Passed list of entities such that attributes will
      //- be actuated for the passed entities and all child
      //- entities.
    
    ~CAActuateSet();
      //- Destructor
    
    void set_current_dimension( int dimension );
      //- Populate the "current list" with all entities
      //- of the specified dimension (4 = body, 3 = volume, etc.)
      //- "current list" is populated with all entities of
      //- the specified dimension from the actuate_list passed
      //- to the constuctor, and children of entities in the
      //- actuate_list with the passed dimension.
    
    inline RefEntity* remove_next();
      //- Pop one entity from the "current list".  Returns NULL
      //- when there are no more entities.
    
    virtual CubitStatus notify_observer( CubitObservable* observable,
                                  const CubitEvent& observer_event,
                                  CubitBoolean from_observable );
      //- Update for destroyed entities.
      
  private:
  
    void append_to_current( DLIList<RefEntity*>& query_results );
      //- Append passed entities to currentList
      
    static DagType get_type_id( int dimension );
  
    DLIList<RefEntity*> typeList[5];
      //- Lists containing the entities passed to the constructor,
      //- grouped and indexed by dimension (body = 4).
    
    DLIList<RefEntity*> currentList;
      //- The list of entities of the same dimension D, including
      //- the entities from the typeList of dimension D and children
      //- of dimension D for all entities in the typeLists above 
      //- with dimension > D.
    
    int currentDimension;
      //- The dimension of entities in currentList.
};
 
//-------------------------------------------------------------------------
// Purpose       : Pop one entity from the current list.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
inline RefEntity* CAActuateSet::remove_next()
{
  return currentList.size() ? currentList.pop() : 0;
}

#endif

