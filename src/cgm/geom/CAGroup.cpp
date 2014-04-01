//- Class:          CAGroup
//- Owner:          Greg Nielson
//- Description:    Cubit Attribute for groups.
//- Checked By:
//- Version:

#include "CAGroup.hpp"
#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "RefGroup.hpp"
#include "CubitObserver.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "TDCAGE.hpp"
#include "GeometryQueryTool.hpp"

#include <stdlib.h>
#include <time.h>

// initialize this CA's static members
CubitBoolean CAGroup::initialize_rand = CUBIT_TRUE;

CubitAttrib* CAGroup_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAGroup *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAGroup(entity);
  }
  else
  {
    new_attrib = new CAGroup(entity, p_csa);
  }

  return new_attrib;
}

CAGroup::CAGroup(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  initialize();
}

CAGroup::CAGroup(RefEntity* new_attrib_owner,
                 CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{

  initialize();
  
  DLIList<int*> *i_list = csa_ptr->int_data_list();
  DLIList<CubitString*> *cs_list = csa_ptr->string_data_list();

    // first, the ints
  i_list->reset();
  
  int num_groups = *(i_list->get_and_step());
  
    // groupID
  int i;
  for (i = num_groups; i > 0; i--)
    groupID.append(*(i_list->get_and_step()));

      // uniqueID
  for (i = num_groups; i > 0; i--)
    uniqueID.append(*(i_list->get_and_step()));

      // sequenceNumbers
  for (i = num_groups; i > 0; i--)
    sequenceNumbers.append(*(i_list->get_and_step()));

    // numOwningGroups
  for (i = num_groups; i > 0; i--)
    numOwningGroups.append(*(i_list->get_and_step()));

    // total number of owning groups
  int total_owning_groups = *(i_list->get_and_step());
  
    // owningGroupID
  for (i = total_owning_groups; i > 0; i--)
    owningGroupID.append(*(i_list->get_and_step()));
    
    // owningUniqueID
  for (i = total_owning_groups; i > 0; i--)
    owningUniqueID.append(*(i_list->get_and_step()));

    // owningSequenceNumbers
  for (i = total_owning_groups; i > 0; i--)
    owningSequenceNumbers.append(*(i_list->get_and_step()));

    // ancestor groups added after first implementation of CAGroup,
    // therefore not all attributes may have ancestor groups
    // setting total_ancestor_groups to zero first will short-circuit
    // loops below if there aren't any ancestor groups, so no need for
    // an 'if' statement

    // total number of ancestor groups
  int total_ancestor_groups = 0;

  if (i_list->size() > (4*num_groups + 3*total_owning_groups + 2))
    total_ancestor_groups = *(i_list->get_and_step());
  
    // ancestorGroupID
  for (i = total_ancestor_groups; i > 0; i--)
    ancestorGroupID.append(*(i_list->get_and_step()));
    
    // ancestorUniqueID
  for (i = total_ancestor_groups; i > 0; i--)
    ancestorUniqueID.append(*(i_list->get_and_step()));
    
    // ancestorOwnedGroupUid
  for (i = total_ancestor_groups; i > 0; i--)
    ancestorOwnedGroupUid.append(*(i_list->get_and_step()));
    
     // ancestorSequenceNumbers
  for (i = total_ancestor_groups; i > 0; i--)
    ancestorSequenceNumbers.append(*(i_list->get_and_step()));
    
    // now, doubles (none)

    // now, strings
  cs_list->reset();

    // attribute internal name (just pop the list)
  cs_list->step();

    // groupNames
  for (i = num_groups; i > 0; i--)
    groupNames.append(new CubitString(*(cs_list->get_and_step())));

    // owningGroupNames
  for (i = total_owning_groups; i > 0; i--)
    owningGroupNames.append(new CubitString(*(cs_list->get_and_step())));

    // ancestorGroupName
  for (i = total_ancestor_groups; i > 0; i--)
    ancestorGroupName.append(new CubitString(*(cs_list->get_and_step())));

    // ok, we're done
}

CAGroup::~CAGroup()
{
  int i;
  for ( i = groupNames.size(); i > 0; i--)
    delete groupNames.get_and_step();
  
  for (i = owningGroupNames.size(); i > 0; i--)
    delete owningGroupNames.get_and_step();
}

void CAGroup::initialize()
{
  if(initialize_rand == CUBIT_TRUE)
  {
#if defined NT || defined _NT || defined CUBIT_LINUX
    srand((unsigned)time(NULL));
#else
    srand48(time(NULL));
#endif
    initialize_rand = CUBIT_FALSE;
  }
}



CubitStatus CAGroup::actuate()
{
  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;
  
    // need to reset all the lists to keep them in sync
  groupID.reset();
  uniqueID.reset();
  groupNames.reset();
  sequenceNumbers.reset();
  owningGroupID.reset();
  owningUniqueID.reset();
  owningGroupNames.reset();
  owningSequenceNumbers.reset();
  numOwningGroups.reset();
  
    // go through all the groups on this CA
  int i;
  for (i = groupID.size(); i > 0; i--) {
      // pop the data for this group off the lists
    int group_id = groupID.get_and_step();
    int unique_id = uniqueID.get_and_step();
    CubitString *group_name = groupNames.get_and_step();
    int seq_num = sequenceNumbers.get_and_step();

    RefGroup *ref_group =
        assign_group(attribOwnerEntity, group_id, unique_id, 
                     group_name, seq_num);
    
    
      // check for groups owning this group
    int owning_groups = numOwningGroups.get_and_step();
  
    for (int j = owning_groups; j > 0; j--)
    {
      int owning_group_id = owningGroupID.get_and_step();
      int owning_unique_id = owningUniqueID.get_and_step();
      CubitString *owning_group_name = owningGroupNames.get_and_step();
      seq_num = owningSequenceNumbers.get_and_step();

      assign_group(ref_group, owning_group_id, owning_unique_id,
                   owning_group_name, seq_num);
    }

  } // loop over all groups in this CAGroup

    // now do ancestors
  ancestorGroupID.reset();
  ancestorUniqueID.reset();
  ancestorOwnedGroupUid.reset();
  ancestorGroupName.reset();
  ancestorSequenceNumbers.reset();
  
  for (i = ancestorGroupID.size(); i > 0; i--) {
    /*RefGroup *ancestor_group = */
      assign_ancestor_group(ancestorGroupID.get_and_step(), 
                            ancestorUniqueID.get_and_step(),
                            ancestorGroupName.get_and_step(),
                            ancestorOwnedGroupUid.get_and_step(),
                            ancestorSequenceNumbers.get_and_step());
  }
  
  deleteAttrib = CUBIT_TRUE;
  hasActuated = CUBIT_TRUE;
  
  return CUBIT_SUCCESS;
}

RefGroup *CAGroup::assign_group(RefEntity *owned_entity,
                                const int group_id, const int unique_id,
                                CubitString *group_name,
                                const int seq_num)
{
  RefGroup* parent_group = NULL;
  CubitBoolean group_id_exists = CUBIT_FALSE;

    // search for the group corresponding to this id and unique id
  RefGroup *ref_group = GeometryQueryTool::instance()->get_last_ref_group();
  for(int i = GeometryQueryTool::instance()->num_ref_groups(); i > 0 && !parent_group; i--)
  {
    ref_group = GeometryQueryTool::instance()->get_next_ref_group();
    if(ref_group->id() == group_id)
    {
      group_id_exists = CUBIT_TRUE;
    }

    ToolData *td_temp = ref_group->get_TD(&TDCAGE::is_cage);
    TDCAGE *td_cagroup = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);
    if (td_cagroup != NULL && td_cagroup->unique_id() == unique_id)
      parent_group = ref_group;
  }
  
  if(!parent_group)
  {
      // else make a new group, and assign id and unique id
      // also assign group name
    parent_group = RefEntityFactory::instance()->construct_RefGroup();
    if(!group_id_exists) {
      parent_group->set_id(0);
      parent_group->set_id(group_id);
    }
    else
      PRINT_INFO("Creating group %d to hold attribute group %d\n",
                 parent_group->id(), group_id);

      // put a td on this new group with the right unique_id
    TDCAGE *td_cagroup = new TDCAGE(unique_id);
    parent_group->add_TD(td_cagroup);

      // add the attribOwnerEntity to the group and name the group
    parent_group->entity_name(*group_name);
  }

    // add the entity to the group with the proper sequence number
  TDCAGE::insert_entity(owned_entity, seq_num, parent_group);
  
  return parent_group;
}
  
RefGroup *CAGroup::assign_ancestor_group(const int ancestor_id,
                                         const int ancestor_uid,
                                         const CubitString *ancestor_name,
                                         const int owned_group_uid,
                                         const int seq_num)
{
    // for each call of this function, we:
    // - get an ancestor group with ancestor_id, ancestor_uid and
    //   ancestor_name (make one if it doesn't exist)
    // - add group with owned_group_uid to that ancestor group
    //   (owned_group_uid should exist, error if not)
    // search for the group corresponding to this id and unique id
  RefGroup* ancestor_group = NULL;
  RefGroup *owned_group = NULL;
  CubitBoolean ancestor_id_exists = CUBIT_FALSE;

  RefGroup *ref_group = GeometryQueryTool::instance()->get_last_ref_group();
  
  for(int i = GeometryQueryTool::instance()->num_ref_groups(); i > 0 && 
        (!ancestor_group || !owned_group); i--)
  {
    ref_group = GeometryQueryTool::instance()->get_next_ref_group();
  
    if(ref_group->id() == ancestor_id)
    {
      ancestor_id_exists = CUBIT_TRUE;
    }
    
    ToolData *td_temp = ref_group->get_TD(&TDCAGE::is_cage);
    TDCAGE *td_cagroup = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);
    if (td_cagroup != NULL && td_cagroup->unique_id() == ancestor_uid) {
      ancestor_group = ref_group;
    }
        
    if (td_cagroup != NULL && td_cagroup->unique_id() == owned_group_uid) {
      owned_group = ref_group;
    }
  }

  assert(owned_group != 0);
  
  if (!ancestor_group)
  {
      // make a new group, and assign id and unique id
      // also assign group name
    ancestor_group = RefEntityFactory::instance()->construct_RefGroup();
;
    if(!ancestor_id_exists) {
      ancestor_group->set_id(0);
      ancestor_group->set_id(ancestor_id);
    }
    else
      PRINT_INFO("Creating group %d to hold attribute group %d\n",
                 ancestor_group->id(), ancestor_id);

      // put a td on this new group with the right unique_id
    TDCAGE *td_cagroup = new TDCAGE(ancestor_uid);
    ancestor_group->add_TD(td_cagroup);

      // add the owned group to the group and name the group
    ancestor_group->entity_name(*ancestor_name);
  }

    // add the entity to the group with the proper sequence number
  TDCAGE::insert_entity(owned_group, seq_num, ancestor_group);

  return ancestor_group;
}
  
CubitStatus CAGroup::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // get the groups containing attribOwnerEntity
  RefGroup* ref_group;
  DLIList<RefGroup*> ref_group_list;
  int i;

  RefGroup::get_groups_within(attribOwnerEntity, ref_group_list, CUBIT_FALSE);

  if( ref_group_list.size() == 0)
  {
    delete_attrib(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }

    // else, this entity is owned by groups
  RefGroup* parent_ref_group;

    // get a td_cage onto the attribOwnerEntity, to keep sequence numbers
    // for the owning groups
  TDCAGE *td_entity = (TDCAGE *) attribOwnerEntity->get_TD(&TDCAGE::is_cage);
  if (!td_entity) {
    td_entity = new TDCAGE(-1);
    attribOwnerEntity->add_TD(td_entity);
  }

    // now, initialize that tdcage
  td_entity->initialize_group_sequence_list(attribOwnerEntity);
  
    // ok, now write the data for the groups to this attribute
  for(i = ref_group_list.size(); i > 0; i--)
  {
      // get the refgroup which gets assigned to this CAGroup
    ref_group = ref_group_list.get_and_step();

      // First, make sure there's a TDCAGE on the RefGroup
    ToolData *td_temp = ref_group->get_TD(&TDCAGE::is_cage);
    TDCAGE *td_cagroup = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);
    if (td_cagroup == NULL) {
#if defined NT || defined _NT || defined CUBIT_LINUX
      td_cagroup = new TDCAGE(rand());
#else
      td_cagroup = new TDCAGE((int)lrand48());
#endif
      ref_group->add_TD(td_cagroup);
      td_cagroup->initialize_group_sequence_list(ref_group);
    }
    
      // append to this CAGroup the id and the unique id of the group;
      // also append the group name
    groupID.append(ref_group->id());
    uniqueID.append(td_cagroup->unique_id());
    groupNames.append(new CubitString(ref_group->entity_name()));

      // get and append the sequence number of the attribOwnerEntity
      // in this group
    int seq_number = td_entity->td_sequence_number(ref_group);
    assert(seq_number != -1);
    sequenceNumbers.append(seq_number);

      // check this group for containing (parent) groups
    DLIList<RefGroup*> parent_ref_group_list;
    RefGroup::get_groups_within(ref_group, parent_ref_group_list, CUBIT_FALSE);

      // append the number of parent groups to the right list
    numOwningGroups.append(parent_ref_group_list.size());
    
      // for each parent group, do essentially the same thing, adding the
      // data to the owningGroup lists
    for(int j = parent_ref_group_list.size(); j>0; j--)
    {
      parent_ref_group = parent_ref_group_list.get_and_step();
      td_temp = parent_ref_group->get_TD(&TDCAGE::is_cage);
      TDCAGE *td_parent = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);
      if (td_parent == NULL) {
#if defined NT || defined _NT || defined CUBIT_LINUX
        td_parent = new TDCAGE(rand());
#else
        td_parent = new TDCAGE((int)lrand48());
#endif
        parent_ref_group->add_TD(td_parent);
        td_parent->initialize_group_sequence_list(parent_ref_group);
      }

        // append to this CAGroup the id and the unique id of the group;
        // also append the group name
      owningGroupID.append(parent_ref_group->id());
      owningUniqueID.append(td_parent->unique_id());
      owningGroupNames.append(new CubitString(parent_ref_group->entity_name()));

        // get and append the sequence number of the group in the parent
        // group
      seq_number = td_cagroup->td_sequence_number(parent_ref_group);
      assert(seq_number != -1);
      owningSequenceNumbers.append(seq_number);

        // finally, build a list of distant ancestors, in case there are
        // groups more than twice removed from any real entities; make
        // it a recursive function, so that it goes all the way up the chain
        // of ancestors
      build_ancestor_list(parent_ref_group);

    } // loop over parent groups
  } // loop over groups containing attribOwnerEntity

  return CUBIT_SUCCESS;
}

CubitStatus CAGroup::reset()
{
  groupID.clean_out();
    //- group ids containing attribOwnerEntity
  
  uniqueID.clean_out();
    //- unique ids of groups containing attribOwnerEntity

  int i;
  for (i = groupNames.size(); i > 0; i--)
    delete groupNames.get_and_step();
  groupNames.clean_out();
    //- names of groups containing attribOwnerEntity

  sequenceNumbers.clean_out();
    //- sequence numbers of this entity in the groups

  numOwningGroups.clean_out();
    //- for each group in groupID, number of groups owning those groups

  owningGroupID.clean_out();
    //- group ids containing groups containing attribOwnerEntity
  
  owningUniqueID.clean_out();
    //- unique ids of groups containing groups containing attribOwnerEntity

  for (i = owningGroupNames.size(); i > 0; i--)
    delete owningGroupNames.get_and_step();
  owningGroupNames.clean_out();
    //- names of groups containing groups containing attribOwnerEntity

  owningSequenceNumbers.clean_out();
    //- sequence numbers of groups in owning groups
  
    //- for each ancestor (a group which owns only other groups, with those
    //- those groups owning only other groups), we store the group id, uid,
    //- name, and the uid of the owned group to which this is an ancestor
  ancestorGroupID.clean_out();
  ancestorUniqueID.clean_out();
  for (i = ancestorGroupName.size(); i > 0; i--)
    delete ancestorGroupName.get_and_step();
  ancestorGroupName.clean_out();
  ancestorOwnedGroupUid.clean_out();
  ancestorSequenceNumbers.clean_out();
  
  return CUBIT_SUCCESS;
}

void CAGroup::build_ancestor_list(RefGroup *parent_ref_group)
{
  DLIList<RefGroup*> ancestor_ref_group_list;
  RefGroup::get_groups_within(parent_ref_group, ancestor_ref_group_list, CUBIT_FALSE);

    // now, recursively work on ancestor list, adding owning group id,
    // uid and name, and owned group uid, to lists
  ToolData *td_temp = parent_ref_group->get_TD(&TDCAGE::is_cage);
  TDCAGE *td_parent = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);

  for (int j = ancestor_ref_group_list.size(); j > 0; j--) {
    RefGroup *ancestor = ancestor_ref_group_list.get_and_step();
    td_temp = ancestor->get_TD(&TDCAGE::is_cage);
    TDCAGE *td_cagroup = (td_temp ? CAST_TO(td_temp, TDCAGE) : NULL);
    if (td_cagroup == NULL) {
#if defined NT || defined _NT || defined CUBIT_LINUX
      td_cagroup = new TDCAGE(rand());
#else
      td_cagroup = new TDCAGE((int)lrand48());
#endif
      ancestor->add_TD(td_cagroup);
      td_cagroup->initialize_group_sequence_list(ancestor);
    }

    ancestorGroupID.append(ancestor->id());
    ancestorUniqueID.append(td_cagroup->unique_id());
    ancestorGroupName.append(new CubitString(ancestor->entity_name()));
    ancestorOwnedGroupUid.append(td_parent->unique_id());
    
      // get and append the sequence number of the group in the parent
      // group
    int seq_number = td_parent->td_sequence_number(ancestor);
    assert(seq_number != -1);
    ancestorSequenceNumbers.append(seq_number);

    build_ancestor_list(ancestor);
  }
}
      
CubitSimpleAttrib* CAGroup::cubit_simple_attrib()
{
  DLIList<CubitString*> cs_list;
  DLIList<double> d_list;
  DLIList<int> i_list;

    // first, the ints
    // groupID
  groupID.reset();
  i_list.append(groupID.size());
  int i;
  for (i = groupID.size(); i > 0; i--)
    i_list.append(groupID.get_and_step());
    
    // uniqueID
  uniqueID.reset();
  for (i = uniqueID.size(); i > 0; i--)
    i_list.append(uniqueID.get_and_step());
    
    // sequenceNumbers
  sequenceNumbers.reset();
  for (i = sequenceNumbers.size(); i > 0; i--)
    i_list.append(sequenceNumbers.get_and_step());
    
    // numOwningGroups
  numOwningGroups.reset();
  for (i = numOwningGroups.size(); i > 0; i--)
    i_list.append(numOwningGroups.get_and_step());

    // size of owningGroupID
  i_list.append(owningGroupID.size());
  
    // owningGroupID
  owningGroupID.reset();
  for (i = owningGroupID.size(); i > 0; i--)
    i_list.append(owningGroupID.get_and_step());
    
    // owningUniqueID
  owningUniqueID.reset();
  for (i = owningUniqueID.size(); i > 0; i--)
    i_list.append(owningUniqueID.get_and_step());
    
    // owningSequenceNumbers
  owningSequenceNumbers.reset();
  for (i = owningSequenceNumbers.size(); i > 0; i--)
    i_list.append(owningSequenceNumbers.get_and_step());
    
    // size of ancestorGroupID
  i_list.append(ancestorGroupID.size());
  
    // ancestorGroupID
  ancestorGroupID.reset();
  for (i = ancestorGroupID.size(); i > 0; i--)
    i_list.append(ancestorGroupID.get_and_step());

    // ancestorUniqueID
  ancestorUniqueID.reset();
  for (i = ancestorUniqueID.size(); i > 0; i--)
    i_list.append(ancestorUniqueID.get_and_step());

    // ancestorOwnedGroupUid
  ancestorOwnedGroupUid.reset();
  for (i = ancestorOwnedGroupUid.size(); i > 0; i--)
    i_list.append(ancestorOwnedGroupUid.get_and_step());

    // ancestorSequenceNumbers
  ancestorSequenceNumbers.reset();
  for (i = ancestorSequenceNumbers.size(); i > 0; i--)
    i_list.append(ancestorSequenceNumbers.get_and_step());

    // now, doubles (none)

    // now, strings
    // attribute internal name
  cs_list.append(new CubitString(att_internal_name()));

    // groupNames
  groupNames.reset();
  for (i = groupID.size(); i > 0; i--)
    cs_list.append(new CubitString(*(groupNames.get_and_step())));

    // owningGroupNames
  owningGroupNames.reset();
  for (i = owningGroupNames.size(); i > 0; i--)
    cs_list.append(new CubitString(*(owningGroupNames.get_and_step())));

    // ancestorGroupName
  ancestorGroupName.reset();
  for (i = ancestorGroupName.size(); i > 0; i--)
    cs_list.append(new CubitString(*(ancestorGroupName.get_and_step())));

  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list,
                                                          &d_list,
                                                          &i_list);

  for (i = cs_list.size(); i--; ) delete cs_list.get_and_step();
  
  return csattrib_ptr;
}

void CAGroup::has_written(CubitBoolean has_written)
{
    //- overloaded has_written function, resets td_cage on owner
  if (has_written == CUBIT_TRUE && hasWritten == CUBIT_FALSE)
      // reset the td_cage on the owner
    attribOwnerEntity->delete_TD(&TDCAGE::is_cage);
    
  hasWritten = has_written;
}
