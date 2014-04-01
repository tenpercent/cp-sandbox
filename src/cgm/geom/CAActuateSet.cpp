//-------------------------------------------------------------------------
// Filename      : CAActuateSet.cpp
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

#include "CAActuateSet.hpp"
#include "RefEntity.hpp"
#include "ModelEntity.hpp"
#include "ModelQueryEngine.hpp"
#include "Body.hpp"
#include "BasicTopologyEntity.hpp"

#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "RefVolume.hpp"
#include "RefGroup.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
CAActuateSet::CAActuateSet( DLIList<RefEntity*>& actuate_list )
  : currentDimension(-1)
{
    // register this as a static observer so that we can
    // remove entities from the lists as they are destroyed
  register_static_observer( this );
  
    // put all entities in the actuate_list into the typeList
    // for the appropriate dimension of entity.
  for( int i = actuate_list.size(); i--; )
  {
    RefEntity* entity_ptr = actuate_list.get_and_step();
    int dimension = entity_ptr->dimension();
    if( dimension < 0 ) // body
      dimension = 4;
    typeList[dimension].append( entity_ptr );
  }
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
CAActuateSet::~CAActuateSet()
{
    // remove from static observer list
  unregister_static_observer( this );
}

//-------------------------------------------------------------------------
// Purpose       : Populate currentList with entities of the specified
//                 dimension.  This includes children of any entities 
//                 of a higher dimension that are in the lists managed
//                 by this object.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
void CAActuateSet::set_current_dimension( int dimension )
{
  int i;
  DLIList<ModelEntity*> query_source, query_target;
  DLIList<RefEntity*> temp_list;
  
    // Clean out current list before adding new entities
  currentList.clean_out();
  
    // Get the target type to query for.
  DagType type = get_type_id( dimension );

    // Get children of higher-order entities
  for( i = 4; i > dimension; i-- )
  {
    query_source.clean_out();
    query_target.clean_out();
    
    CAST_LIST( typeList[i], query_source, ModelEntity );
    ModelQueryEngine::instance()
      ->query_model( query_source, type, query_target );
    
    temp_list.clean_out();
    CAST_LIST( query_target, temp_list, RefEntity );
    
    append_to_current( temp_list );
  }
  
    // Add lcoal entities of current dimension
  append_to_current( typeList[dimension] );

    // Save current dimension
  currentDimension = dimension;
}

//-------------------------------------------------------------------------
// Purpose       : Add entities to currentList
//
// Special Notes : Make sure no duplicates are added
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
void CAActuateSet::append_to_current( DLIList<RefEntity*>& list )
{
  int i;
  
    // Set marks on all new entities
  for( i = list.size(); i--; )
    list.get_and_step()->marked(1);
  
    // Clear marks on entities in current list, 
    // including those also in the new list.
  for( i = currentList.size(); i--; )
    currentList.get_and_step()->marked(0);
  
    // Any entities in the new list that are still
    // marked are not already in the current list.
    // Add them.
  for( i = list.size(); i--; )
  {
    RefEntity* entity_ptr = list.get_and_step();
    if( entity_ptr->marked() )
    {
      currentList.append( entity_ptr );
      entity_ptr->marked(0);
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Remove deleted entities from lists
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/02
//-------------------------------------------------------------------------
CubitStatus CAActuateSet::notify_observer( CubitObservable* observable,
                               const CubitEvent& observer_event,
                               CubitBoolean )
{
  RefEntity* entity_ptr;
  int dimension;
  
    // only care about entities that are destroyed
  if( observer_event.get_event_type() != MODEL_ENTITY_DESTRUCTED &&
      observer_event.get_event_type() != ENTITY_DESTRUCTED )
    return CUBIT_SUCCESS;


    // is it a body?
  if( (entity_ptr = dynamic_cast<Body*>(observable) ) != NULL )
    dimension = 4;

    // is it some other topology entity
  else if( (entity_ptr = dynamic_cast<BasicTopologyEntity*>(observable) )
           != NULL )
    dimension = entity_ptr->dimension();

    // otherwise we don't care about it
  else
    return CUBIT_FAILURE;
  

    // if it exists in the type list, remove it
  if( typeList[dimension].move_to( entity_ptr ) )
    typeList[dimension].extract();
  
    // if it exists in the current list, remove it.
  if( dimension == currentDimension && 
      currentList.move_to( entity_ptr ) )
    currentList.extract();
    
  return CUBIT_SUCCESS;
}

DagType CAActuateSet::get_type_id( int dimension )
{
  switch( dimension )
  {
    case 4:  return DagType::body_type();     
    case 3:  return DagType::ref_volume_type();
    case 2:  return DagType::ref_face_type();  
    case 1:  return DagType::ref_edge_type();  
    case 0:  return DagType::ref_vertex_type();
    default: assert(0); 
             return DagType::invalid_type();
  }
}
