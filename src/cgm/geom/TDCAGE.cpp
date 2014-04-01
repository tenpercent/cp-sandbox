#include "TDCAGE.hpp"
#include "RefGroup.hpp"
#include "RefEntity.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"

int TDCAGE::is_cage(const ToolData* td)
{
  return (CAST_TO(const_cast<ToolData*>(td), TDCAGE) != NULL);
}


int TDCAGE::td_sequence_number(const RefGroup *group)
{
    // go through the group list and get the sequence number;
    // return -1 if the group isn't in the list

  groupList.reset();
  sequenceList.reset();
  for (int i = groupList.size(); i > 0; i--) {
    if (group == groupList.get_and_step()) return sequenceList.get();

    else sequenceList.step();
  }
  
  return -1;
}

void TDCAGE::insert_entity(RefEntity *entity, const int seq_num,
                           RefGroup *into_group)
{
    //- insert the entity into the group following the sequence number
    //- given
  TDCAGE *td_cage;
  CubitBoolean inserted = CUBIT_FALSE;
  DLIList<RefEntity*> copy_group_list;
  copy_group_list = into_group->entityList;
  
  into_group->entityList.last();
  for (int i = into_group->entityList.size(); i > 0; i--) {
    RefEntity *into_entity = into_group->entityList.get();
    td_cage = (TDCAGE *) into_entity->get_TD(&TDCAGE::is_cage);
    if (td_cage && td_cage->td_sequence_number(into_group) < seq_num) {
	  if ( !copy_group_list.move_to(entity) )
        into_group->entityList.insert(entity);
      inserted = CUBIT_TRUE;
      break;
    }
    else into_group->entityList.back();
  }
  
    // if we got here and didn't insert, either the list is empty
    // or it isn't and our seq_num is less than all items; either
    // way, we insert at beginning
  if (inserted == CUBIT_FALSE){
	  if ( !copy_group_list.move_to(entity) )
	    into_group->entityList.insert_first(entity);
  }

    // now register the observable entity in the group
  into_group->register_observable(entity);

    // now, make sure the inserted entity has a td_cage and this group
    // and seq number are registered
  td_cage = (TDCAGE *) entity->get_TD(&TDCAGE::is_cage);
  if (!td_cage) {
    td_cage = new TDCAGE(-1);
    entity->add_TD(td_cage);
  }
  
  td_cage->insert_group(into_group, seq_num);
}

int TDCAGE::group_sequence_number(RefGroup *group, const RefEntity *entity)
{
    //- return the sequence number of the given entity in the given group

    // built for speed; don't accumulate index, use loop counter for that
  group->entityList.reset();
  for (int i = group->entityList.size(); i > 0; i--) {
    if (group->entityList.get_and_step() == entity)
        // index in reverse direction, so subtract to get real index
      return group->entityList.size() - i;
  }
  
  return -1;
}

void TDCAGE::initialize_group_sequence_list(RefEntity *entity) 
{
    //- get the groups owning entity and build the sequence lists
  DLIList<RefGroup*> ref_groups;
  RefGroup::get_groups_within(entity, ref_groups, CUBIT_FALSE);

  groupList.clean_out();
  sequenceList.clean_out();
  
  for (int i = ref_groups.size(); i > 0; i--) {
    RefGroup *ref_group = ref_groups.get_and_step();
    groupList.append(ref_group);
    sequenceList.append(group_sequence_number(ref_group, entity));
  }
}

    
  
