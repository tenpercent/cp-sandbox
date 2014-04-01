//-------------------------------------------------------------------------
// Filename      : GroupingEntity.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 GroupingEntity. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96 
//
// Owner         :  Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "GroupingEntity.hpp"
#include "SenseEntity.hpp"
#include "BasicTopologyEntity.hpp"

#include "DLIList.hpp"
#include "ModelQueryEngine.hpp"

// ********** END CUBIT INCLUDES           **********


// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 09/03/96
//-------------------------------------------------------------------------
GroupingEntity::~GroupingEntity()
{
  if (myParent)
    myParent->remove_grouping_entity(this);
  while (firstSenseEntity && remove_sense_entity(firstSenseEntity));
    
  assert (!myParent && !nextInParent && !firstSenseEntity);
}

//-------------------------------------------------------------------------
// Purpose       : This function returns a list of SenseEntity
//                 pointers associated with this grouping entity.
//
// Special Notes : Complete reimplementation - j.k. July 2003
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/28/96
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::get_sense_entity_list( DLIList<SenseEntity*>& list)
{
  if (!firstSenseEntity)
    return CUBIT_SUCCESS;
  
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next())
  {
    assert(ptr->get_grouping_entity_ptr() == this);
    list.append(ptr);
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function adds a sense entity to the list of 
//                 sense entities of a grouping entity.
//
// Special Notes : In the DAG, the sense entities associated with the 
//                 current grouping entity are linked from one sense
//                 entity to another in order to maintain the order in
//                 which the sense entities were added to the grouping
//                 entity.
//
//                 Complete reimplementation - j.k. July 2003
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::add_sense_entity(SenseEntity *sense_entity,
                                             SenseEntity *after_this) 
{
    // Check to make sure that we are getting the correct type of 
    // SenseEntity.
  if ( dag_type() != sense_entity->dag_type().parent() )
     return CUBIT_FAILURE ;
   
    // Check that the sense entity is not already in some other
    // grouping entity
  if ( sense_entity->get_grouping_entity_ptr() )
    return CUBIT_FAILURE;
  
    // prev and next ptrs should be NULL if sense entity is not
    // in a grouping entity
  assert (!sense_entity->next() && !sense_entity->previous());
  
  if (after_this)
  {
    if (after_this->get_grouping_entity_ptr() != this )
      return CUBIT_FAILURE;
  
    if (!sense_entity->gpe_insert_after(after_this))
      return CUBIT_FAILURE;
      
    if (after_this == lastSenseEntity)
      lastSenseEntity = sense_entity;
  }
  else if (lastSenseEntity)
  {
    if (!sense_entity->gpe_insert_after(lastSenseEntity))
      return CUBIT_FAILURE;
    lastSenseEntity = sense_entity;
  }
  else
  {
    firstSenseEntity = lastSenseEntity = sense_entity;
  }
  
  sense_entity->set_grouping_entity_ptr(this);
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Remove a child sense entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::remove_sense_entity( SenseEntity* sense_entity_ptr )
{
  if (sense_entity_ptr->get_grouping_entity_ptr() != this)
    return CUBIT_FAILURE;
  
  if (firstSenseEntity == sense_entity_ptr)
    firstSenseEntity = sense_entity_ptr->next();
  if (lastSenseEntity == sense_entity_ptr)
    lastSenseEntity = sense_entity_ptr->previous();
  
  sense_entity_ptr->gpe_remove();
  sense_entity_ptr->set_grouping_entity_ptr(NULL);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Change/update/re-order child SenseEntity list
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/03/03
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::set_sense_entity_list( 
                                        DLIList<SenseEntity*>& list,
                                        DLIList<SenseEntity*>& removed )
{
  int i;
  
    // Remove all?
  if (list.size() == 0)
  {
    get_sense_entity_list( removed );
    return disconnect_all_children();
  }
  
    // Check for error conditions before modifying anything
  list.reset();
  for (i = list.size(); i--; )
  {
    SenseEntity* sense_entity = list.get_and_step();
    
      // Check to make sure that we are getting the correct type of 
      // SenseEntity.
    if ( dag_type() != sense_entity->dag_type().parent() )
       return CUBIT_FAILURE ;
   
      // Check that the sense entity is not already in some other
      // grouping entity
    if ( sense_entity->get_grouping_entity_ptr() &&
         sense_entity->get_grouping_entity_ptr() != this )
      return CUBIT_FAILURE;
  }
  
    // Special case for first entity in list.
  list.reset();
  SenseEntity* new_first = list.get_and_step();
    // No sense entities currently attached...
  if (!firstSenseEntity)
  {
    firstSenseEntity = lastSenseEntity = new_first;
    new_first->set_grouping_entity_ptr(this);
  }
    // Already attached, but not first in list...
  else if( firstSenseEntity != new_first )
  {
    if (!new_first->get_grouping_entity_ptr())
      new_first->set_grouping_entity_ptr(this);
    else
    {
      if (lastSenseEntity == new_first)
        lastSenseEntity = new_first->previous();
      new_first->gpe_remove();
    }
      
    new_first->gpe_insert_before(firstSenseEntity);
    firstSenseEntity = new_first;
  }
  
    // Now loop through remaining sense entities.
  SenseEntity* prev = new_first;
  for (i = list.size() - 1; i--; )
  {
    SenseEntity* curr = list.get_and_step();

      // If next sense entity in input list is not
      // next sense entity in this GroupingEntity...
    if (prev->next() != curr)
    {
      if (!curr->get_grouping_entity_ptr())
        curr->set_grouping_entity_ptr(this);
      else
      {
        if (lastSenseEntity == curr)
          lastSenseEntity = curr->previous();
        curr->gpe_remove();
      }
      curr->gpe_insert_after(prev);
    }
    
      // update lastSenseEntity if necessary...
    if (lastSenseEntity == prev)
      lastSenseEntity = curr;
      
      // iterate
    prev = curr;
  }
  
    // Disconnect any sense entities in this GroupingEntity
    // that were not in in the input list (they should now
    // be at the end of the list of sense entities in this)
    // and pass them back in the 'removed' list.
  CubitStatus result = CUBIT_SUCCESS;
  while (prev != lastSenseEntity)
  {
    removed.append(prev->next());
    if (!remove_sense_entity(prev->next()))
    {
      assert(0);
      result = CUBIT_FAILURE;
      prev = prev->next();
    }
  }
  
  return result;
}
  

//-------------------------------------------------------------------------
// Purpose       : Invert
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
void GroupingEntity::reverse_direction()
{
  SenseEntity* ptr;
  if (!firstSenseEntity)
    return;
  
    // For each child sense entity
  for (ptr = firstSenseEntity; ptr; ptr = ptr->previous())
  {
      // change linked list pointers
    ptr->swap_gpe_list_ptrs();
      // change sense
    ptr->reverse_sense();
  }
  
    // Change first pointer the old last entity
    // (Preserves order as returnd by old DLIList rep and 
    //  makes this work to reverse a chain.)
  ptr = firstSenseEntity;
  firstSenseEntity = lastSenseEntity;
  lastSenseEntity = ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Get parent basic topology entity.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int GroupingEntity::get_parents( DLIList<ModelEntity*>* list ) const
{
  if (!myParent)
    return 0;
  
  if (list)
    list->append(myParent);
  
  return 1;
}

//-------------------------------------------------------------------------
// Purpose       : Get child sense entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int GroupingEntity::get_children( DLIList<ModelEntity*>* list ) const
{
  if (!firstSenseEntity)
    return 0;
  
  int count = 0;
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next() )
  {
    assert(ptr->get_grouping_entity_ptr() == this);
    if(list)
      list->append(ptr);
    count++;
  }
  
  return count;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : Remove a child sense entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::remove_child_link( ModelEntity* child_ptr )
{
  SenseEntity* se = dynamic_cast<SenseEntity*>(child_ptr);
  if (!se)
    return CUBIT_FAILURE;
  
  return remove_sense_entity(se);
}

//-------------------------------------------------------------------------
// Purpose       : Remove from all parent BasicTopologyEntities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::disconnect_all_parents( DLIList<ModelEntity*>* list )
{
  if (!myParent)
    return CUBIT_SUCCESS;
  if (list) 
    list->append(myParent);
  return myParent->remove_grouping_entity(this);
}

//-------------------------------------------------------------------------
// Purpose       : Remove all child SenseEntitys
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus GroupingEntity::disconnect_all_children( DLIList<ModelEntity*>* list)
{
  while (firstSenseEntity)
  {
    if (list)
      list->append(firstSenseEntity);
    if (!remove_sense_entity(firstSenseEntity))
      return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
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
CubitBoolean GroupingEntity::query_append_parents( DLIList<ModelEntity*>& list )
{
  if (myParent && !ModelQueryEngine::instance()->encountered(myParent))
  {
    list.append(myParent);
    return CUBIT_TRUE;
  }
  
  return CUBIT_FALSE;
}
CubitBoolean GroupingEntity::query_append_children( DLIList<ModelEntity*>& list )
{
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  CubitBoolean found_some = CUBIT_FALSE;
  
  if (!firstSenseEntity)
    return CUBIT_FALSE;
  
  for (SenseEntity* ptr = firstSenseEntity; ptr; ptr = ptr->next())
  {
    if (!mqe->encountered(ptr))
    {
      list.append(ptr);
      found_some = CUBIT_TRUE;
    }
  }
    
  return found_some;
}
