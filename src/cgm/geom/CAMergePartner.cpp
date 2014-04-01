//- Class:          CAMergePartner
//- Owner:          Greg Nielson
//- Description:    Cubit Attribute for merge partners.
//- Checked By:
//- Version:

#include "CAMergePartner.hpp"
#include "DLIList.hpp"
#include "SDLCAMergePartnerList.hpp"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "TDCompare.hpp"
#include "TDUniqueId.hpp"
#include "MergeTool.hpp"
#include "CastTo.hpp"
#include "CubitAttribUser.hpp"
#include "GeometryQueryTool.hpp"
#include "RefVolume.hpp"
#include "TopologyBridge.hpp"

#include <stdlib.h>
#include <time.h>

// factory function
CubitAttrib* CAMergePartner_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAMergePartner *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAMergePartner(entity);
  }
  else
  {
    new_attrib = new CAMergePartner(entity, p_csa);
  }

  return new_attrib;
}

// initialize this CA's static members
CAMergePartner::CAMergePartner (RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  initialize();
  
}

CAMergePartner::CAMergePartner(RefEntity* new_attrib_owner,
                               CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  initialize();
  
  assert(csa_ptr->int_data_list() != NULL);
  int * i_temp = csa_ptr->int_data_list()->get_and_step();
  if(i_temp != NULL)
    mergeID = *i_temp;
  if( csa_ptr->int_data_list()->size() > 1 )
  {
    i_temp = csa_ptr->int_data_list()->get_and_step();
    switch( *i_temp ) {
      case -1: bridge_sense_ = CUBIT_REVERSED; break;
      case  1: bridge_sense_ = CUBIT_FORWARD;  break;
      case  0: bridge_sense_ = CUBIT_UNKNOWN; break;
      default: bridge_sense_ = CUBIT_UNKNOWN; assert(0);
    }
  }
}

CAMergePartner::~CAMergePartner() 
{
    //
  mergeID = 0;
}


void CAMergePartner::initialize()
{
  mergeID = -1;
  bridge_sense_ = CUBIT_UNKNOWN;
  isSurvivor = 0;
}


CubitStatus CAMergePartner::actuate()
{
  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;
  
  BasicTopologyEntity * bte_ptr = CAST_TO(attribOwnerEntity,BasicTopologyEntity); 
  if (bte_ptr == NULL || bte_ptr->deactivated() == CUBIT_TRUE)
    return CUBIT_FAILURE;
  
  DLIList<RefEntity*> merge_list;
  merge_prepare(merge_list);
  merge_list.append_unique(attribOwnerEntity);

  CubitStatus result = CUBIT_SUCCESS;

  if (merge_list.size() > 1)
    result = (NULL == MergeTool::instance()->force_merge(merge_list))
            ? CUBIT_FAILURE : CUBIT_SUCCESS;

  if( result )
  {
    hasActuated = CUBIT_TRUE;
  }        

  return result;
}

void CAMergePartner::merge_prepare(DLIList<RefEntity*> &merge_list)
{
  DLIList<CubitAttrib*> my_ca_list;
  CAMergePartner *my_camp_ptr;
  RefEntity* re_ptr;

    // get all the merge partner attributes that are on my owner
  attribOwnerEntity->find_cubit_attrib_type(CA_MERGE_PARTNER, my_ca_list);
  merge_list.clean_out();
  DLIList<ToolDataUser*> td_list, temp_td_list;
  int i;
  for (i = my_ca_list.size(); i > 0; i--)
  {
    my_camp_ptr = CAST_TO(my_ca_list.get(),CAMergePartner);
    my_ca_list.step();
    td_list.clean_out();
      // get all the objects with this unique id (which is also the merge id)
    TDUniqueId::find_td_unique_id(my_camp_ptr->merge_id(), temp_td_list);
    td_list += temp_td_list;
  }
  
    // now put those entities into the merge_list
  for (i = td_list.size(); i > 0; i--) {
    re_ptr = CAST_TO(td_list.get(), RefEntity);
    if (re_ptr) merge_list.append(re_ptr);
    td_list.step();
  }
  
    // Now get bridge sense for each entity in list.
    // Add this entity to list, too.
  merge_list.append( attribOwnerEntity );
  for( i = merge_list.size(); i--; )
  {
    RefEntity* ent = merge_list.get_and_step();
    TopologyEntity* te = dynamic_cast<TopologyEntity*>(ent);
    if( te->bridge_manager()->number_of_bridges() != 1 )
      continue;
    
    my_ca_list.clean_out();
    ent->find_cubit_attrib_type(CA_MERGE_PARTNER, my_ca_list);
    assert( my_ca_list.size() < 2);
    if( !my_ca_list.size() )
      continue;
    
    my_camp_ptr = dynamic_cast<CAMergePartner*>(my_ca_list.get());
    if( my_camp_ptr->bridge_sense() == CUBIT_UNKNOWN )
      continue;
  }  
  merge_list.pop(); // take attribOwnerEntity back off list
    
  return;
}
  
CubitStatus CAMergePartner::actuate_list(DLIList<RefEntity*> entity_list)
{

    // given a list of ref entities (usually all entities of a given type),
    // actuate the camp's on those entities
  RefEntity *ref_ent, *keeper;
  DLIList<CubitAttrib*> ca_list;
  DLIList<TopologyBridge*> bridge_list(entity_list.size());
  SDLCAMergePartnerList sorted_camp_list;
  int i;
  for(i = entity_list.size(); i > 0; i--)
  {
    ref_ent = entity_list.get_and_step();
    ca_list.clean_out();
    ref_ent->find_cubit_attrib_type(CA_MERGE_PARTNER, ca_list);
    assert(ca_list.size() < 2); // There should only be one
                                //  merge partner per entity
    if(ca_list.size() > 0)
    {
      CAMergePartner* attrib = dynamic_cast<CAMergePartner*>(ca_list.get());
      sorted_camp_list.append( attrib );
      TopologyEntity* te = dynamic_cast<TopologyEntity*>(ref_ent);
      TopologyBridge* bridge = te->bridge_manager()->topology_bridge();
      bridge_list.append(bridge);
    }
  }
  sorted_camp_list.sort();
  sorted_camp_list.reset();

  if (DEBUG_FLAG(90)) {
    for (i = sorted_camp_list.size(); i > 0; i--) {
      CAMergePartner *camp_ptr = sorted_camp_list.get_and_step();
      ref_ent = camp_ptr->attrib_owner();
      PRINT_DEBUG_90("%s %d, unique id = %d\n", ref_ent->class_name(),
                     ref_ent->id(), camp_ptr->merge_id());
    }
  }
  
    // now go through all the camp's for the entity list; peel off
    // camp's with the same id, and merge the associated entities together
  while (sorted_camp_list.size() > 0) 
  {
    DLIList<RefEntity*> refent_list;
    DLIList<CubitAttrib*> camp_list;
    keeper = NULL;

      // get the next list of entities with the same camp id
    sorted_camp_list.last();
    CAMergePartner *camp_ptr = sorted_camp_list.remove();
    sorted_camp_list.back();
    camp_list.append(camp_ptr);
    int current_id = camp_ptr->merge_id();
    while (sorted_camp_list.size() > 0 &&
           sorted_camp_list.get()->merge_id() == current_id) {
      camp_list.append(sorted_camp_list.remove());
      sorted_camp_list.back();
    }
    
    if (camp_list.size() == 1) continue;
    
    CubitBoolean has_actuated = camp_list.get()->has_actuated();

      // check the has actuated flag; if one is set, they all should be;
      // also, compile list of ref entities while we're at it
    for (current_id = camp_list.size(); current_id > 0; current_id--) {
      ref_ent = camp_list.get()->attrib_owner();
      refent_list.append(ref_ent);
      if (!keeper || ref_ent->id() < keeper->id()) keeper = ref_ent;
      assert(camp_list.get()->has_actuated() == has_actuated);
      camp_list.step();
    }

      // if they have already actuated, go on to next ones
    if (has_actuated == CUBIT_TRUE) continue;

      // otherwise merge
    if(refent_list.size() > 1)
      MergeTool::instance()->force_merge(refent_list);

        // remove the cubit attribute from the surviving parent
    keeper->remove_cubit_attrib(CA_MERGE_PARTNER);

  } // loop over existing camp's
  
  return CUBIT_SUCCESS;
}

CubitStatus CAMergePartner::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
    // set the updated flag
  hasUpdated = CUBIT_TRUE;
  BasicTopologyEntity* bte_ptr = CAST_TO(attribOwnerEntity, BasicTopologyEntity);

  if( (bte_ptr == NULL) || (dynamic_cast<RefVolume*>(bte_ptr) != NULL))
  {
    delete_attrib(CUBIT_TRUE);
  }
  else if( (bte_ptr->bridge_manager()->number_of_bridges() == 1) &&
           (bte_ptr->bridge_manager()->topology_bridge()->bridge_sense() != CUBIT_REVERSED) )
  {
    delete_attrib(CUBIT_TRUE);
  }
  else
  {
      // get the merge id from the TDUniqueId for the owner entity
    mergeID = TDUniqueId::get_unique_id(attribOwnerEntity);
    bridge_sense_ = CUBIT_UNKNOWN;
    isSurvivor = 0;
  }

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CAMergePartner::cubit_simple_attrib()
{
  DLIList<CubitString*> cs_list;
  DLIList<double> d_list;
  DLIList<int> i_list;

  i_list.append(mergeID);
  cs_list.append(new CubitString(att_internal_name()));

  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list,
                                                          &d_list,
                                                          &i_list);
  int i;
  for ( i = cs_list.size(); i--;) delete cs_list.get_and_step();
  
  return csattrib_ptr;
}

void CAMergePartner::set_survivor( CubitSimpleAttrib* csa, int is_survivor )
{
  //get the list we want to modify from the CSA
  DLIList<int*>* data = csa->int_data_list();
  
  //change or append?
  if( data->size() >= 4 )
  {
    data->reset();
    data->step();
    data->step();
    *(data->next()) = is_survivor;
  }
  else 
  {
    assert( data->size() == 3 );
    data->append( new int(is_survivor) );
    data->reset();
  }
}

void CAMergePartner::set_bridge_sense( CubitSimpleAttrib* csa, CubitSense sense )
{
  //encode/decode sense as:
  // CUBIT_FORWARD :  1
  // CUBIT_REVERSE : -1
  // CUBIT_UNKNOWN :  0
  int i = 0;
  switch( sense ) { 
    case CUBIT_FORWARD  : i =  1; break;
    case CUBIT_REVERSED : i = -1; break;
    case CUBIT_UNKNOWN  : i =  0; break;
    default: assert(0);
  }
  
  //get the list we want to modify from the CSA
  DLIList<int*>* data = csa->int_data_list();
  
  //change or append?
  if( data->size() >= 2 )
  {
    data->reset();
    *(data->next()) = i;
  }
  else 
  {
    assert( data->size() == 1 );
    data->append( new int(i) );
    data->reset();
  }
}

CubitBoolean CAMergePartner::is_survivor( CubitSimpleAttrib* csa )
{
  DLIList<int*>* data = csa->int_data_list();
  if( data && (data->size() >= 4) )
  {
    data->reset();
    data->step();
    data->step();
    int i = *(data->next());
    if( i == 1 )
      return true;
  }
  return false;
}


CubitSense CAMergePartner::get_bridge_sense( CubitSimpleAttrib* csa )
{
  DLIList<int*>* data = csa->int_data_list();
  if( data && (data->size() >= 2) )
  {
    data->reset();
    int i = *(data->next());
    switch(i) {
      case  1: return CUBIT_FORWARD;
      case -1: return CUBIT_REVERSED;
      default: assert(0);  // if -DNDEBUG, fall through to unknown
      case  0: return CUBIT_UNKNOWN;
    }
  }

  return CUBIT_UNKNOWN;
}

void CAMergePartner::set_saved_id( CubitSimpleAttrib* csa, int id )
{
   //get the list we want to modify from the CSA
  DLIList<int*>* data = csa->int_data_list();

    // ID goes after bridge sense, so save bridge sense first
  if (data->size() == 1)
    set_bridge_sense( csa, CUBIT_UNKNOWN );
  
    // change?
  if (data->size() > 2 )
  {
    data->reset();
    *data->next(2) = id;
  }
    // set?
  else
  {
    assert(data->size() == 2);
    data->append( new int(id) );
  }
}

int CAMergePartner::get_saved_id( CubitSimpleAttrib* csa )
{
   //get the list we want to modify from the CSA
  DLIList<int*>* data = csa->int_data_list();

  if (data->size() < 3) 
    return 0;
  
  data->reset();
  return *data->next(2);
}
    
  

void CAMergePartner::print() 
{
  PRINT_INFO("Attribute MERGE_PARTNER, %s %d: mergeId = %d, sense = %d.\n",
             attribOwnerEntity->class_name(), attribOwnerEntity->id(),
             mergeID, bridge_sense_);
}
