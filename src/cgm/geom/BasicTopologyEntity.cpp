///-------------------------------------------------------------------------
// Filename      : BasicTopologyEntity.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 BasicTopologyEntity. 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/15/96 
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "SenseEntity.hpp"
#include "BasicTopologyEntity.hpp"
#include "GroupingEntity.hpp"
#include "GeometryEntity.hpp"
#include "DLIList.hpp"
#include "ModelQueryEngine.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
BasicTopologyEntity::~BasicTopologyEntity()
{
  while (firstSenseEntity && remove_sense_entity(firstSenseEntity));
  while (firstGroupingEntity && remove_grouping_entity(firstGroupingEntity));
  assert(!firstSenseEntity && !firstGroupingEntity);
}


//-------------------------------------------------------------------------
// Purpose       : This function gets the list of GroupingEntities of this 
//                 BasicTopologyEntity.
//
// Special Notes : In the DAG, the GroupingEntities associated with the 
//                 current BasicTopologyEntity are linked from one GroupingEntity
//                 to the next in order to maintain the order in
//                 which the GroupingEntities were added to the BTE.
//
//                 Complete reimplementation - jk, July 2003
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::get_grouping_entity_list(
    DLIList<GroupingEntity*>& list) const
{
  for (GroupingEntity* ptr = firstGroupingEntity; ptr; ptr = ptr->next())
  {
    assert(ptr->get_basic_topology_entity_ptr() == this);
    list.append(ptr);
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get the parent SenseEntities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/08/99
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::get_sense_entity_list(
    DLIList<SenseEntity*>& list) const
{
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
  {
    assert(ptr->get_basic_topology_entity_ptr() == this);
    list.append(ptr);
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Attach a child grouping entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::add_grouping_entity( GroupingEntity* gpe_ptr) 
{
  if (gpe_ptr->dag_type().parent() != dag_type())
    return CUBIT_FAILURE;
  
  if (gpe_ptr->get_basic_topology_entity_ptr())
    return CUBIT_FAILURE;
    
  assert(!gpe_ptr->next());
  
  if (firstGroupingEntity)
  {
    if (!gpe_ptr->insert_after(lastGroupingEntity))
      return CUBIT_FAILURE;
  }
  else
  {
    firstGroupingEntity = gpe_ptr;
  }  
    
  lastGroupingEntity = gpe_ptr;
  gpe_ptr->set_basic_topology_entity_ptr(this);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Disconnect child grouping entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::remove_grouping_entity( GroupingEntity* gpe_ptr) 
{
  if (gpe_ptr->get_basic_topology_entity_ptr() != this)
    return CUBIT_FAILURE;
  
  GroupingEntity* next_gpe_ptr = gpe_ptr->next();
  GroupingEntity* prev_gpe_ptr = gpe_ptr->previous();
  
  if (!gpe_ptr->remove_from_list())
    return CUBIT_FAILURE;
  
  if (firstGroupingEntity == gpe_ptr)
    firstGroupingEntity = next_gpe_ptr;
  if (lastGroupingEntity == gpe_ptr)
    lastGroupingEntity = prev_gpe_ptr;

  gpe_ptr->set_basic_topology_entity_ptr(0);
  return CUBIT_SUCCESS;
}
  
//-------------------------------------------------------------------------
// Purpose       : Update list of child grouping entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/12/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::set_grouping_entity_list( 
                                        DLIList<GroupingEntity*>& list,
                                        DLIList<GroupingEntity*>& removed )
{
  int i;
  
    // Remove all?
  if (list.size() == 0)
  {
    get_grouping_entity_list( removed );
    return disconnect_all_children();
  }
  
    // Check for error conditions before modifying anything
  list.reset();
  for (i = list.size(); i--; )
  {
    GroupingEntity* grouping_entity = list.get_and_step();
    
      // Check to make sure that we are getting the correct type of 
      // GroupingEntity.
    if ( dag_type() != grouping_entity->dag_type().parent() )
       return CUBIT_FAILURE ;
   
      // Check that the grouping entity is not already in some other
      // basic topology entity
    if ( grouping_entity->get_basic_topology_entity_ptr() &&
         grouping_entity->get_basic_topology_entity_ptr() != this )
      return CUBIT_FAILURE;
  }
  
    // Special case for first entity in list.
  list.reset();
  GroupingEntity* new_first = list.get_and_step();
    // No sense entities currently attached...
  if (!firstGroupingEntity)
  {
    firstGroupingEntity = lastGroupingEntity = new_first;
    new_first->set_basic_topology_entity_ptr(this);
  }
    // Already attached, but not first in list...
  else if( firstGroupingEntity != new_first )
  {
    if (!new_first->get_basic_topology_entity_ptr())
      new_first->set_basic_topology_entity_ptr(this);
    else
    {
      if (lastGroupingEntity == new_first)
        lastGroupingEntity = new_first->previous();
      new_first->remove_from_list();
    }
      
    new_first->insert_before(firstGroupingEntity);
    firstGroupingEntity = new_first;
  }
  
    // Now loop through remaining sense entities.
  GroupingEntity* prev = new_first;
  for (i = list.size() - 1; i--; )
  {
    GroupingEntity* curr = list.get_and_step();

      // If next grouping entity in input list is not
      // next grouping entity in this BTE...
    if (prev->next() != curr)
    {
      if (!curr->get_basic_topology_entity_ptr())
        curr->set_basic_topology_entity_ptr(this);
      else
      {
        if (lastGroupingEntity == curr)
          lastGroupingEntity = curr->previous();
        curr->remove_from_list();
      }
      curr->insert_after(prev);
    }
    
      // update lastSenseEntity if necessary...
    if (lastGroupingEntity == prev)
      lastGroupingEntity = curr;
      
      // iterate
    prev = curr;
  }
  
    // Disconnect any grouping entities in this BTE
    // that were not in in the input list (they should now
    // be at the end of the list of grouping entities in this)
    // and pass them back in the 'removed' list.
  CubitStatus result = CUBIT_SUCCESS;
  while (prev != lastGroupingEntity)
  {
    removed.append(prev->next());
    if (!remove_grouping_entity(prev->next()))
    {
      assert(0);
      result = CUBIT_FAILURE;
      prev = prev->next();
    }
  }
  
  return result;
}
  



//-------------------------------------------------------------------------
// Purpose       : Attach a parent sense entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::add_sense_entity( SenseEntity* se_ptr ) 
{
  if (se_ptr->dag_type() != dag_type().parent())
    return CUBIT_FAILURE;
  
  if (se_ptr->get_basic_topology_entity_ptr())
    return CUBIT_FAILURE;
    
  assert(!se_ptr->next_on_bte());
  if (firstSenseEntity)
  {
    assert(!lastSenseEntity->next_on_bte());
    lastSenseEntity->set_bte_next(se_ptr);
  }
  else
  {
    firstSenseEntity = se_ptr;
  }
  
  lastSenseEntity = se_ptr;
  se_ptr->set_basic_topology_entity_ptr(this);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Disconnect parent sense entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::remove_sense_entity( SenseEntity* se_ptr) 
{
  if (se_ptr->get_basic_topology_entity_ptr() != this)
    { assert(0); return CUBIT_FAILURE; }
  
  if (!firstSenseEntity)
    return CUBIT_FAILURE;
  
  if (firstSenseEntity == se_ptr)
  {
    if (lastSenseEntity == se_ptr)
    {
      firstSenseEntity = lastSenseEntity = 0;
    }
    else
    {
      firstSenseEntity = se_ptr->next_on_bte();
    }
  }
  else
  {
    SenseEntity* prev = firstSenseEntity;
    while (prev->next_on_bte() != se_ptr)
    {
      prev = prev->next_on_bte();
      if (!prev)
        return CUBIT_FAILURE;
    }
    
    prev->set_bte_next( se_ptr->next_on_bte() );
    if (lastSenseEntity == se_ptr)
      lastSenseEntity = prev;
  }
  
  se_ptr->set_bte_next(0);
  se_ptr->set_basic_topology_entity_ptr(0);
  return CUBIT_SUCCESS;
}
   
//-------------------------------------------------------------------------
// Purpose       : Find connecting sense entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
SenseEntity* BasicTopologyEntity::find_sense_entity(GroupingEntity* gpe) const
{
  SenseEntity *ptr, *result = 0;
  for (ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
  {
    if (ptr->get_grouping_entity_ptr() == gpe)
    {
      if (result)
        return 0;
      result = ptr;
    }
  }
  
  return result;
}
SenseEntity* BasicTopologyEntity::find_sense_entity(BasicTopologyEntity* bte) const
{
  SenseEntity *ptr, *result = 0;
  for (ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
  {
    GroupingEntity* gpe = ptr->get_grouping_entity_ptr();
    BasicTopologyEntity* tmp_bte = gpe ? gpe->get_basic_topology_entity_ptr() : 0;
    if (tmp_bte == bte)
    {
      if (result)
        return 0;
      result = ptr;
    }
  }
  
  return result;
}

CubitStatus BasicTopologyEntity::get_sense_entities( 
                                DLIList<SenseEntity*>& result,
                                GroupingEntity* in_this )
{
  int input_size = result.size();
  
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
    if (!in_this || ptr->get_grouping_entity_ptr() == in_this)
      result.append(ptr);
  
  return result.size() > input_size ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus BasicTopologyEntity::get_sense_entities( 
                                DLIList<SenseEntity*>& result,
                                BasicTopologyEntity* in_this )
{
  int input_size = result.size();
  
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
    if (!in_this || ptr->get_parent_basic_topology_entity_ptr() == in_this)
      result.append(ptr);
  
  return result.size() > input_size ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Check if this entity is nonmanifold in the parent
//                 sense entity.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitBoolean BasicTopologyEntity::is_nonmanifold(GroupingEntity* gpe)
{
  int count = 0;
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
    if (ptr->get_grouping_entity_ptr() == gpe)
      count++;
  
  assert(count);
  return count == 1 ? CUBIT_FALSE : CUBIT_TRUE;
}

GeometryType BasicTopologyEntity::geometry_type() const
{
     //- returns type of underlying geometry representation
     //- (see GeometryEntity.hpp for list of types)
   return get_geometry_entity_ptr()->geometry_type();
}

//-------------------------------------------------------------------------
// Purpose       : These functions return the GeometryEntity of this 
//                 BasicTopologyEntity.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------
GeometryEntity* BasicTopologyEntity::get_geometry_entity_ptr() const
{
  TopologyBridge* bridge = bridge_manager()->topology_bridge();
  return dynamic_cast<GeometryEntity*>(bridge);
}


CubitBox BasicTopologyEntity::bounding_box()
{
   return get_geometry_entity_ptr()->bounding_box() ;
}

//-------------------------------------------------------------------------
// Purpose       : Set the GeometryEntity pointer of this BTE
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/24/96
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::set_geometry_entity_ptr(
  GeometryEntity* GE_ptr) 
{
//   if (dag_type().dimension() != GE_ptr.dimension())
//   {
//     PRINT_ERROR("Internal Error: %s:%d: Mismatched BTE/GeometryEntity.\n",
//                  __FILE__,__LINE__);
//     return CUBIT_FAILURE;
//   }
   return TopologyEntity::set_topology_bridge(GE_ptr);
}

double BasicTopologyEntity::measure()
{ return get_geometry_entity_ptr()->measure(); }

//-------------------------------------------------------------------------
// Purpose       : get parent sense entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int BasicTopologyEntity::get_parents( DLIList<ModelEntity*>* list ) const
{
  if (!firstSenseEntity)
    return 0;
  
  int count = 0;
  SenseEntity* ptr;
  
  for (ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
  {
    count++;
    if (list)
      list->append(ptr);
  }

  return count;
}



//-------------------------------------------------------------------------
// Purpose       : get child grouping entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int BasicTopologyEntity::get_children( DLIList<ModelEntity*>* list ) const
{
  if (!firstGroupingEntity)
    return 0;
  
  int count = 0;
  GroupingEntity* ptr;
  
  for (ptr = firstGroupingEntity; ptr; ptr = ptr->next())
  {
    count++;
    if (list)
      list->append(ptr);
  }

  return count;
}

//-------------------------------------------------------------------------
// Purpose       : Functions to support ModelQueryEngine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/24/03
//-------------------------------------------------------------------------
CubitBoolean BasicTopologyEntity::query_append_parents( DLIList<ModelEntity*>& list )
{
  if (!firstSenseEntity)
    return CUBIT_FALSE;
  
  
  CubitBoolean found_some = CUBIT_FALSE;
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next_on_bte())
    if (!mqe->encountered(ptr))
    {
      list.append(ptr);
      found_some = CUBIT_TRUE;
    }
  
  return found_some;  
}
CubitBoolean BasicTopologyEntity::query_append_children( DLIList<ModelEntity*>& list )
{
  if (!firstGroupingEntity)
    return 0;
  
  CubitBoolean found_some = CUBIT_FALSE;
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  GroupingEntity* ptr;
  
  for (ptr = firstGroupingEntity; ptr; ptr = ptr->next())
  {
    if (!mqe->encountered(ptr))
    {
      found_some = CUBIT_TRUE;
      list.append(ptr);
    }
  }

  return found_some;  
}

//-------------------------------------------------------------------------
// Purpose       : Generic function to remove child grouping entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::remove_child_link(ModelEntity* entity_ptr)
{
  GroupingEntity* gpe_ptr = dynamic_cast<GroupingEntity*>(entity_ptr);
  if (!gpe_ptr)
    return CUBIT_FAILURE;
  
  return remove_grouping_entity(gpe_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Remove from parent SenseEntitys
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::disconnect_all_parents( DLIList<ModelEntity*>* list )
{
  while (firstSenseEntity)
  {
    SenseEntity* se_ptr = firstSenseEntity;
    
    if (!remove_sense_entity(se_ptr))
      return CUBIT_FAILURE;

    if (list)
      list->append(se_ptr);
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Disconnect child GroupingEntities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus BasicTopologyEntity::disconnect_all_children( DLIList<ModelEntity*>* list )
{
  while (firstGroupingEntity)
  {
    GroupingEntity* gpe_ptr = firstGroupingEntity;
    
    if (!remove_grouping_entity(gpe_ptr))
      return CUBIT_FAILURE;

    if (list)
      list->append(gpe_ptr);
  }
  return CUBIT_SUCCESS;
}



// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

