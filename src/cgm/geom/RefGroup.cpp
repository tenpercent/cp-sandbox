//-- File: RefGroup.cc

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "RefGroup.hpp"
#include "CubitObservable.hpp"
#include "MergeEvent.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"  
#include "CastTo.hpp"
#include "GeometryQueryTool.hpp"
#include "RefEntityFactory.hpp"

RefGroup::RefGroup( const char* name, int id)
{
  recursionMark = 0;

  if( id == 0 )
    entityId = RefEntityFactory::instance()->next_ref_group_id();
  else 
    entityId = id;

  // assign default names
  assign_default_name();
  
  if (name != NULL)
    entity_name( name );
}

RefGroup::RefGroup (DLIList<RefEntity*>& entity_list)
{
  recursionMark = 0;
  //entityList = entity_list;
  int i;
  for ( i = entity_list.size(); i > 0; i-- )
  {
    RefEntity *ent_ptr = entity_list.get_and_step();
    if ( !entityList.move_to(ent_ptr) )
      entityList.append(ent_ptr);
  }  
  entityList.reset();
  for (i=entityList.size(); i > 0; i--)
    register_observable(entityList.get_and_step());

  entityId = RefEntityFactory::instance()->next_ref_group_id();

    // Notify Model about the creation of this object
  CubitObserver::notify_static_observers(this, MODEL_ENTITY_CONSTRUCTED) ;

  // assign default names
  assign_default_name();

}

RefGroup::RefGroup(int /*proe_type*/)
{
   // This constructor was created just to get around not doing a notify,
   // so we can use it for Pro/E parts and assemblies.
  recursionMark = 0;
}

RefGroup::~RefGroup ()
{
   remove_all_ref_entities();
}

CubitStatus RefGroup::add_ref_entity(RefEntity *ref_entity)
{
  if (!entityList.move_to(ref_entity)) 
  {
    entityList.append(ref_entity);
    register_observable(ref_entity);
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

CubitStatus RefGroup::add_ref_entity(DLIList<RefEntity*>& entity_list)
{
  RefEntity* entity = NULL;
  for (int i=entity_list.size(); i > 0; i--) 
  {
    entity = entity_list.get_and_step();
    if ( entityList.move_to(entity) == CUBIT_FALSE )
    {
       entityList.append(entity);
       register_observable(entity);
    }
  }

  return CUBIT_SUCCESS ;
}

CubitStatus RefGroup::remove_ref_entity(RefEntity *entity,
                                        const CubitBoolean from_observable)
{
  
  if (entityList.remove(entity) != NULL) {
    unregister_observable(entity, from_observable);
    return CUBIT_SUCCESS;
  }
  else return CUBIT_FAILURE;
}

CubitStatus RefGroup::remove_ref_entity(DLIList<RefEntity*> &entity_list,
                                        const CubitBoolean from_observable)
{
  
  int i;
  RefEntity *entity;
  for (i = entity_list.size(); i > 0; i--) {
    entity = entity_list.get_and_step();
    if (entityList.remove(entity) != NULL) {
      unregister_observable(entity, from_observable);
    }
  }
    
  return CUBIT_SUCCESS;
}

int RefGroup::remove_all_ref_entities()
{
  int num_entities = entityList.size();
  for (int i=entityList.size(); i > 0; i--) {
    unregister_observable(entityList.get_and_step());
  }
  entityList.clean_out();
  return num_entities;
}

void RefGroup::get_child_ref_entities(DLIList<RefEntity*>& entity_list)
{
  entityList.reset();
  entity_list.merge_unique(entityList, CUBIT_TRUE );
}

void RefGroup::get_child_entities(DLIList<CubitEntity*>& cub_entity_list)
{
  entityList.reset();
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( entityList, temp_list );
  cub_entity_list.merge_unique(temp_list, CUBIT_TRUE);
//  cub_entity_list.merge_unique(entityList, CUBIT_TRUE);
}

void RefGroup::expand_group( DLIList<RefEntity*> & entity_list )
{
    //This function will get the list of child entities for this group and 
    //if there are any groups in its child list, it will call this fucntion
    //again on those groups.  The end result should be that the returned
    //list does not contain groups... just basic_topology entities
    //(the other ref_entities) 
  DLIList<RefEntity*> group_list, mixed_list;
  RefEntity *ref_entity_ptr;
  RefGroup *group;
  mixed_list.append( this );

  // don't swap the upper and lower bounds on ii.
  for ( int ii = 0; ii < mixed_list.size(); ii++ )
  {
    //Note, the list gets added to, so don't do a simple get_and_step
    mixed_list.reset();
    mixed_list.step(ii);
    ref_entity_ptr = mixed_list.get();
    
    group = CAST_TO( ref_entity_ptr, RefGroup );
    if (group)
      {
	// avoid infinite recursion of groups A and B containing each other
        if (!group_list.move_to(group)) {
	  group->get_child_ref_entities(mixed_list);
	  group_list.append(group);
        }
      }
    // In the future if we have other RefEntities that are grouping
    // entities then we should have another statement here...
    else
      entity_list.append( ref_entity_ptr );
  }

  return;
}

void RefGroup::get_parent_ref_entities(DLIList<RefEntity*>&)
{
  //- appends all ref entities that own this to entity_list.
  //- Goes up just one dimension.
  // There is nothing above us. Do nothing
}

int RefGroup::maximum_dimension()
{
  recursionMark = 1;
  // This routine returns the maximum dimension of its owned subentities.
  // The 'only' kludge is that if one of its subentities is a RefGroup, it
  // must call 'maximum_dimension' on that entity instead of dimension().
  
  RefEntity *entity;
  int test_dim = 0;
  int dimension = 0;
  
  for (int i = entityList.size(); i > 0; i--)  
  {
    entity = entityList.get_and_step();
    RefGroup* group;
    if ( (group = CAST_TO(entity,RefGroup) ) != NULL )
    {
      if (group->recursionMark == 0)
	test_dim = group->maximum_dimension();
    } 
    else 
    {
      test_dim = entity->dimension();
    }
    dimension = CUBIT_MAX(dimension, test_dim);
  }
  recursionMark = 0;
  return dimension;
}

CubitBox RefGroup::bounding_box()
{
  recursionMark = 1;
 
  CubitBox super_box = CubitBox();
  int super_box_defined = CUBIT_FALSE;
  
  for (int i = entityList.size(); i > 0; i--) {
      
    RefEntity *entity = entityList.get_and_step();
    RefGroup *group = CAST_TO(entity, RefGroup);
    if (!group || group->recursionMark == 0) {
      CubitBox entity_box = entity->bounding_box();
      
      // "Concatenate" this box with the super_box, creating a bounding
      // box that bounds the entities (from the list), processed so far.
      if (super_box_defined)
	super_box |= entity_box;
      else {
	super_box = entity_box;
	super_box_defined = CUBIT_TRUE;
      }
    }
  }
  recursionMark = 0;

  return super_box;
}

void RefGroup::get_sub_entities(DLIList<RefEntity*> &entity_list)
{  
  recursionMark = 1;

  //- appends all ref entities owned by this entity on entity_list
  //- and recurses all the way down to dimension 0
  DLIList<RefEntity*> local_entity_list;
  get_child_ref_entities(local_entity_list);

  // *need* to merge now if a group
  // else more efficient to merge later
  entity_list.merge_unique(local_entity_list, CUBIT_TRUE);
DLIList<RefEntity*> temp_list, temp_list2;


  for (int i=local_entity_list.size(); i > 0; i--) {
    // take some care to avoid infinite recursion
    RefEntity *child = local_entity_list.get_and_step();
    RefGroup *group = CAST_TO(child, RefGroup);
    if (!group || group->recursionMark == 0) {
	temp_list2.clean_out();
        child->get_all_child_ref_entities(temp_list2);  
 	temp_list.merge_unique(temp_list2);
     }
  }

  entity_list.merge_unique(temp_list);
  recursionMark = 0;
}

void RefGroup::is_mergeable(AutoMergeStatus val)
{
  recursionMark = 1;
  autoMergeStatus = val;
  DLIList<RefEntity*> children;
  get_child_ref_entities( children );
  for ( int i = children.size(); i > 0; i-- ) {
    RefEntity *child = children.get_and_step();
    RefGroup *group = CAST_TO(child, RefGroup);
    if (!group || group->recursionMark == 0)
      child->is_mergeable(val);
  }
  recursionMark = 0;
}

CubitVector RefGroup::center_point()
{ return bounding_box().center(); }


int RefGroup::subtract(RefGroup *group_to_subtract, RefGroup *target_group)
{
  // Get the list of RefEntities associated with each group
  DLIList<RefEntity*> final_list;
  DLIList<RefEntity*> entity_list_2;

  get_child_ref_entities(final_list);
  group_to_subtract->get_child_ref_entities(entity_list_2);

  // At this point, we have three groups, all non-null. Note that
  // some or all of these may be the same group, so don't destroy
  // the target group until all information is generated.
  
  // 1. Final = group 2
  // 2. Remove from final all items that are in both final and group 1
  for (int i=entity_list_2.size(); i > 0; i--) {
    RefEntity *entity = entity_list_2.get_and_step();
    final_list.remove(entity);
  }

  target_group->remove_all_ref_entities();
  target_group->add_ref_entity(final_list);

  return CUBIT_SUCCESS;
}

int RefGroup::intersect(RefGroup *other_group, RefGroup *target_group)
{
  // Get the list of RefEntities associated with each group
  DLIList<RefEntity*> entity_list_1;
  DLIList<RefEntity*> entity_list_2;

  get_child_ref_entities(entity_list_1);
  other_group->get_child_ref_entities(entity_list_2);

  DLIList<RefEntity*> final_list;

  // At this point, we have three groups, all non-null. Note that
  // some or all of these may be the same group, so don't destroy
  // the target group until all information is generated.
  
  // Final = all items in both group 1 and group 2.
  for (int i=entity_list_2.size(); i > 0; i--) {
    RefEntity *entity = entity_list_2.get_and_step();
    if (entity_list_1.move_to(entity)) {
      final_list.append_unique(entity);
    }
  }

  target_group->remove_all_ref_entities();
  target_group->add_ref_entity(final_list);

  return CUBIT_SUCCESS;
}
  
int RefGroup::unite(RefGroup *other_group, RefGroup *target_group)
{
  // Get the list of RefEntities associated with each group
  DLIList<RefEntity*> final_list;
  DLIList<RefEntity*> entity_list_2;

  get_child_ref_entities(final_list);
  other_group->get_child_ref_entities(entity_list_2);

  // At this point, we have three groups, all non-null. Note that
  // some or all of these may be the same group, so don't destroy
  // the target group until all information is generated.
  
  // 1. Final = group 1
  // 2. Add all items that are in group 2, but not already in group 1
  for (int i=entity_list_2.size(); i > 0; i--) {
    RefEntity *entity = entity_list_2.get_and_step();
    final_list.append_unique(entity);
  }
  
  target_group->remove_all_ref_entities();
  target_group->add_ref_entity(final_list);

  return CUBIT_SUCCESS;
}

int RefGroup::validate()
{
  // NOTE: RefGroup::validate() should not call RefEntity::validate()
  //       directly since the contained entities will make that
  //       call in their respective validate() functions.

  recursionMark = 1;
  int error = 0;
  for (int i = entityList.size(); i > 0; i--) {
    RefEntity *entity = entityList.get_and_step();
    RefGroup *group = CAST_TO(entity, RefGroup);
    if (!group || group->recursionMark == 0) 
      error += entity->validate();
  }
  recursionMark = 0;
  return error;
}

CubitStatus RefGroup::delete_group(RefGroup *group_ptr, CubitBoolean propagate)
{
   // This function will delete the corresponding group from the Model
   // If propagate is CUBIT_TRUE, the contained groups are deleted also.

   group_ptr->notify_all_observers( MODEL_ENTITY_DESTRUCTED );

   if( propagate )
   {
      DLIList<RefGroup*> contained_groups;
      get_contained_groups ( group_ptr, contained_groups );
      for( int i=0; i<contained_groups.size(); i++ )
      {
         delete_group( contained_groups.get_and_step(), false );
      }
   }
   
   group_ptr->remove_all_ref_entities();
   delete group_ptr;

   return CUBIT_SUCCESS;
}

void RefGroup::delete_all_groups()
{
   DLIList<RefGroup*> ref_groups;
   GeometryQueryTool::instance()->ref_groups(ref_groups);
   int num_groups = ref_groups.size();
   for( int i=0; i<num_groups; i++ )
      delete_group( ref_groups.get_and_step() );

}

void RefGroup::get_contained_groups (RefGroup *group_ptr, DLIList<RefGroup*> &contained_groups)
{
   DLIList<RefEntity*> group_list, mixed_list;
   RefEntity *ref_entity_ptr;
   RefGroup *group;
   
   mixed_list.append( group_ptr );
   
   // Using expand_group algorithm to get sub_groups of group_ptr
   // don't swap the upper and lower bounds on ii.
   for ( int ii = 0; ii < mixed_list.size(); ii++ )
   {
      // Note, the list gets added to, so don't do a simple get_and_step
      mixed_list.reset();
      mixed_list.step(ii);
      ref_entity_ptr = mixed_list.get();
      
      group = CAST_TO( ref_entity_ptr, RefGroup ); //which it should be the first time around
      if (group)
      {
         // avoid infinite recursion of groups A and B containing each other
         if (!group_list.move_to(group)) {
            group->get_child_ref_entities(mixed_list);
            group_list.append(group);
         }  
      }
   }
   CAST_LIST(group_list, contained_groups, RefGroup);
}

void RefGroup::get_groups_within( CubitEntity* cubit_entity_ptr, 
                                  DLIList<RefGroup*> &groups_within,
                                  const CubitBoolean recursive)
{
  RefEntity* ref_entity_ptr = CAST_TO( cubit_entity_ptr, RefEntity );
  if( ref_entity_ptr != NULL )
  {
    get_groups_within( ref_entity_ptr, groups_within, recursive);
    return;
  }

  return;
}

void RefGroup::get_groups_within( RefEntity* ref_entity_ptr, 
                                  DLIList<RefGroup*> &groups_within,
                                  const CubitBoolean recursive)
{
  RefGroup* ref_group_ptr;
  
  // Get the observer's of this entity
  DLIList <CubitObserver*> observer_list;
  ref_entity_ptr->get_observer_list( observer_list );
  
  // Now, loop through the observer list, adding them
  for (int i=0; i<observer_list.size(); i++)
  {
    ref_group_ptr = CAST_TO (observer_list.get(), RefGroup);
    
    if( ref_group_ptr != NULL )
    {
      // Add this observer
      groups_within.append_unique( ref_group_ptr );
      
      // Add all the groups the observer is within
      if (recursive) get_groups_within( ref_group_ptr, groups_within );
    }
    observer_list.step();
  }
}

// Here we need to move through all the model's group's, since observer's
// aren't recursive (i.e., if group 2 contains group 3 which contains
// hex 1, then hex 1's observer is only group 3.  Group 2 actually
// contains hex 1 however, but through group 3).
void RefGroup::get_groups_within( RefGroup* ref_group_ptr, 
                                  DLIList<RefGroup*> &groups_within,
                                  const CubitBoolean recursive)
{
  // Get the model's groups
  DLIList<RefGroup*> model_group_list;
  GeometryQueryTool::instance()->ref_groups(model_group_list);
  int i, j;
  RefGroup* model_group;
  RefGroup* contained_group;
  
  // Go through all the groups, finding which ones contain
  // the entity, at any level below (checks groups within
  // groups).
  for( i=0; i<model_group_list.size(); i++ )
  {
    model_group = model_group_list.get_and_step();
    
    // See if it exists directly in model_group
    DLIList<RefEntity*> contained_ref_entities;
    model_group->get_child_ref_entities( contained_ref_entities );
    if( contained_ref_entities.move_to( (RefEntity*)ref_group_ptr ) )
      groups_within.append_unique( model_group );

    if (!recursive) continue;
    
    // Now check all of model group's contained groups
    DLIList<RefGroup*> contained_group_list;
    model_group->get_contained_groups( model_group, contained_group_list );
    contained_group_list.remove( model_group );
    for( j=0; j<contained_group_list.size(); j++ )
    {
      contained_group = contained_group_list.get_and_step();
      
      contained_ref_entities.clean_out();
      contained_group->get_child_ref_entities( contained_ref_entities );
      
      if( contained_ref_entities.move_to( (RefEntity*)ref_group_ptr ) )
        groups_within.append_unique( model_group );
    }
  }
}

CubitStatus RefGroup::notify_observer(CubitObservable *observable,
                                      const CubitEvent &observer_event,
                                      CubitBoolean from_observable)
{
  RefEntity *entity = CAST_TO(observable, RefEntity);
  if (entity != NULL) {
    int event = observer_event.get_event_type();
    if (event == MODEL_ENTITY_CONSTRUCTED ||
        event == ENTITY_CONSTRUCTED) {
      add_ref_entity(entity);
      return CUBIT_SUCCESS;
    }
    else if (event == MODEL_ENTITY_DESTRUCTED ||
             event == ENTITY_DESTRUCTED) {
      remove_ref_entity(entity, from_observable);
      return CUBIT_SUCCESS;
    }
    else if ( event == ENTITIES_MERGED )
    {
        // Out with the old...
      remove_ref_entity(entity, from_observable);

        // ...in with the new.
      MergeEvent *merge_event = (MergeEvent *) &observer_event;
      RefEntity *kept_entity = 
        CAST_TO( merge_event->get_kept_entity(), RefEntity );
      add_ref_entity(kept_entity);
    }
  }
  
  else {
  
    // else it's an error, unrecognized entity
    PRINT_ERROR("Notify called with an observable that's unrecognized.\n");
    assert(0);
  }
  return CUBIT_FAILURE;
}

/*
void RefGroup::draw (int color)
{
   recursionMark = 1;
   if (is_visible())
   {
      for (int i = entityList.size(); i > 0; i--)
      {
         RefEntity *entity = entityList.get_and_step();
         RefGroup *group = CAST_TO(entity, RefGroup);
         if (!group || group->recursionMark == 0)
            entity->draw(color);
      }
   }
   recursionMark = 0;
}
*/

