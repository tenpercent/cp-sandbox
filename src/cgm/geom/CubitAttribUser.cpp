//- class:          CubitAttribUser
//- Owner:          Greg Nielson
//- Description:    implementation of CubitAttribUser class.
//- Checked By:
//- Version: $Id:


// #include "CubitSimpleAttrib.hpp"
#include "CubitAttribUser.hpp"
#include "CastTo.hpp"
#include "Body.hpp"
#include "GeometryEntity.hpp"
#include "BasicTopologyEntity.hpp"
#include "BodySM.hpp"
#include "CAGroup.hpp"
#include "CAMergePartner.hpp"
//#include "CADeferredAttrib.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
// #include "DLIList.hpp"
// #include "CubitString.hpp"
#include "GeometryModifyTool.hpp"
#include "BridgeManager.hpp"


CubitAttribUser::CubitAttribUser(CubitAttrib* cubit_attrib_ptr)
{
  headAttrib = cubit_attrib_ptr;
}

CubitAttribUser::~CubitAttribUser()
{
  //delete all ToolData's chained off this user.
  CubitAttrib *ca_ptr = headAttrib;
  
  while ( ca_ptr) {
    CubitAttrib *next = ca_ptr->next_attrib();
    delete ca_ptr;
    ca_ptr = next;
  }

  headAttrib = NULL;
}

CubitAttrib* CubitAttribUser::get_cubit_attrib (int attrib_type,
                                                CubitBoolean create_if_missing)
{
  CubitAttrib* cubit_attrib_ptr = NULL;
  RefEntity* entity = NULL;
  DLIList<CubitAttrib*> attrib_list;
  find_cubit_attrib_type (attrib_type, attrib_list);
  if (attrib_list.size() > 0)
    cubit_attrib_ptr =  attrib_list.get();
  else if ( create_if_missing == CUBIT_TRUE )
  {
    entity = CAST_TO(this, RefEntity);
    cubit_attrib_ptr = CGMApp::instance()->attrib_manager()->create_cubit_attrib(attrib_type, entity, NULL);
  }

  return cubit_attrib_ptr;
}

CubitStatus CubitAttribUser::add_cubit_attrib (CubitAttrib* cubit_attrib_ptr)
{
  if (cubit_attrib_ptr == NULL)
    return CUBIT_FAILURE;

  CubitStatus cubit_attribute_status = CUBIT_FAILURE;
  CubitAttrib *temp_ptr;
  
  if (headAttrib == NULL)
  {
    headAttrib = cubit_attrib_ptr;
    cubit_attribute_status = CUBIT_SUCCESS;
  }
  else
  {
    for(temp_ptr = headAttrib;
        temp_ptr->next_attrib() != NULL;
        temp_ptr = temp_ptr->next_attrib());
    temp_ptr->set_next_attrib(cubit_attrib_ptr);
    cubit_attribute_status = CUBIT_SUCCESS;
  }
  return cubit_attribute_status;
}

CubitStatus CubitAttribUser::put_simple_attrib 
  (CubitSimpleAttrib* new_csattrib_ptr, CubitBoolean append_it)
{
  Body* Body_ptr;
  Body_ptr = CAST_TO(this, Body);
  BasicTopologyEntity* BTE_ptr;
  BTE_ptr = CAST_TO(this, BasicTopologyEntity);
  if ( Body_ptr != NULL )
  {
      // Get the OSME pointer
    BodySM* OSME_ptr = Body_ptr->get_body_sm_ptr();
    
      // check for duplicates
    if (DEBUG_FLAG(94))
    {
      DLIList<CubitSimpleAttrib*> cs_list;
      OSME_ptr->get_simple_attribute(cs_list);
      for (int i = cs_list.size(); i > 0; i--)
      {
        CubitSimpleAttrib *cs_attrib = cs_list.get_and_step();
        if (CubitSimpleAttrib::equivalent(cs_attrib, new_csattrib_ptr))
          PRINT_INFO("Trying to add equivalent attribute of type %s on Body %d.\n",
                     new_csattrib_ptr->character_type().c_str(), Body_ptr->id());
        else if (cs_attrib->character_type() == new_csattrib_ptr->character_type())
          PRINT_INFO("Trying to add attribute of same type %s on Body %d.\n",
                     new_csattrib_ptr->character_type().c_str(), Body_ptr->id());
      }
    }
    
//________  Change Code by DZ of Cat,  3/11/99 12:25:30 PM  ________
      // Attach this name to it
    if ( append_it )
       append_simple_attribute(OSME_ptr, new_csattrib_ptr);
//________  Change End by DZ of Cat,  3/11/99 12:25:30 PM  ________
  }
    // Deal with BasicTopologyEntities
  else if ( BTE_ptr != NULL )
  {
      // Get the GeometryEntity pointer
    GeometryEntity* GE_ptr =
      BTE_ptr->get_geometry_entity_ptr();
    
      // check for duplicates
    if (DEBUG_FLAG(94))
    {
      DLIList<CubitSimpleAttrib*> cs_list;
      GE_ptr->get_simple_attribute(cs_list);
      for (int i = cs_list.size(); i > 0; i--)
      {
        CubitSimpleAttrib *cs_attrib = cs_list.get_and_step();
        if (CubitSimpleAttrib::equivalent(cs_attrib, new_csattrib_ptr))
          PRINT_INFO("Trying to add equivalent attribute of type %s on %s %d.\n",
                     new_csattrib_ptr->character_type().c_str(), 
                     BTE_ptr->class_name(), BTE_ptr->id());
        else if (cs_attrib->character_type() == new_csattrib_ptr->character_type())
          PRINT_INFO("Trying to add attribute of same type %s on %s %d.\n",
                     new_csattrib_ptr->character_type().c_str(), 
                     BTE_ptr->class_name(), BTE_ptr->id());
      }
    }

//_________  Add Code by DZ of Cat,  3/11/99 11:10:29 AM  _________
      // Attach this name to it
    if ( append_it )
       append_simple_attribute(GE_ptr, new_csattrib_ptr);
//_________  Code End by DZ of Cat,  3/11/99 11:10:29 AM  _________
   
  }
    // Deal with all other RefEntities.  As these do not have any
    // underlying entities (such as solid model entities) associated
    // with them, this block is a do-nothing block.
  else
  {
  }
  
  if (DEBUG_FLAG(90))
  {
    RefEntity *ref_entity = CAST_TO(this, RefEntity);
    PRINT_DEBUG_90( "Putting simple attribute of type %s on"
                " %s %d.\n", new_csattrib_ptr->character_type().c_str(),
                ref_entity->class_name(), ref_entity->id());
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus CubitAttribUser::clear_simple_attribs()
{

  CubitAttrib *cubit_attrib_ptr = NULL;
  CubitAttrib *next_attrib_ptr = NULL;
  for(cubit_attrib_ptr = headAttrib;
      cubit_attrib_ptr != NULL;)
  {
    //ignore Assembly and Name attributes
    next_attrib_ptr = cubit_attrib_ptr->next_attrib();
    if( cubit_attrib_ptr->int_attrib_type() != CA_ENTITY_NAME && 
        cubit_attrib_ptr->int_attrib_type() != CA_ASSEMBLY_DATA )
    {
      remove_cubit_attrib( cubit_attrib_ptr );
      delete cubit_attrib_ptr;
    }
    cubit_attrib_ptr = next_attrib_ptr;
  }
 
  if (DEBUG_FLAG(94))
  {
    PRINT_DEBUG_94("CubitAttribUser::clear_simple_attribs()\n");
  }
  TopologyEntity* te_ptr = dynamic_cast<TopologyEntity*>(this);
  if( !te_ptr )
    return CUBIT_FAILURE;
  
  remove_all_simple_attribute(te_ptr->bridge_manager()->topology_bridge());
  write_cubit_attrib_by_type(CA_ASSEMBLY_DATA);
  write_cubit_attrib_by_type( CA_ENTITY_NAME );
  set_written_flag(CUBIT_FALSE);
  return CUBIT_SUCCESS; 
}

CubitStatus CubitAttribUser::write_specific_cubit_attrib(CubitAttrib* cubit_attrib_ptr)
{
  CubitStatus result = CUBIT_SUCCESS;
  
    // don't write if this attrib is marked as deleted or it has already been written
    // Also, only write if the write flag is on for this attribute type
  if ((cubit_attrib_ptr->delete_attrib() != CUBIT_TRUE) &&
//      (cubit_attrib_ptr->has_written() != CUBIT_TRUE) &&
      (CGMApp::instance()->attrib_manager()->auto_write_flag(cubit_attrib_ptr->int_attrib_type()) == CUBIT_TRUE)) 
  {

    if (DEBUG_FLAG(90)) {
      RefEntity *this_entity = CAST_TO(this, RefEntity);
      PRINT_DEBUG_90( "Writing attribute for %s %d, type %s\n",
                  this_entity->class_name(), this_entity->id(),
                  cubit_attrib_ptr->att_internal_name());
    }

      // set has_written flag here, before actually writing - a trick to
      // allow some attributes to detect they're being written, and do any
      // final setup they need to do before the actual write takes place (used
      // in CAMeshContainer, for example)
    cubit_attrib_ptr->has_written(CUBIT_TRUE);

    CubitSimpleAttrib *csa_ptr = cubit_attrib_ptr->cubit_simple_attrib();
    assert(csa_ptr != 0);
    
    result = put_simple_attrib(csa_ptr);

      // if the write wasn't successful, reset the write flag
    if (result != CUBIT_SUCCESS) cubit_attrib_ptr->has_written(CUBIT_FALSE);

    delete csa_ptr;
  }
  
  return result;
}
  
CubitStatus CubitAttribUser::write_cubit_attrib_by_type(int attrib_type)
{
  DLIList<CubitAttrib*> attrib_list;
  find_cubit_attrib_type(attrib_type, attrib_list);
  CubitStatus write_status = write_cubit_attrib_list(attrib_list);
  return write_status;
}
  
CubitStatus CubitAttribUser::write_cubit_attrib_list(DLIList<CubitAttrib*> 
                                                       attrib_list)
{
  CubitStatus write_status = CUBIT_SUCCESS;
  attrib_list.reset();
  for(int i = attrib_list.size(); i > 0; i--)
  {
    CubitAttrib* cubit_attrib_ptr = attrib_list.get_and_step();
    if (write_specific_cubit_attrib(cubit_attrib_ptr) == CUBIT_FAILURE)
      write_status = CUBIT_FAILURE;
  }
  return write_status;
}
  
CubitStatus CubitAttribUser::write_cubit_attribs()
{
    // write the attributes
  CubitStatus write_status = CUBIT_SUCCESS;
  DLIList<CubitAttrib*> attrib_list;
  CubitAttrib *attrib;
  get_cubit_attrib_list(attrib_list);
  attrib_list.reset();
  int i;
  for(i = attrib_list.size(); i != 0; i--)
  {
    attrib = attrib_list.get_and_step();
    if (write_specific_cubit_attrib(attrib) == CUBIT_FAILURE)
      write_status = CUBIT_FAILURE;
  }
  
  return write_status;
}
 
void CubitAttribUser::split_owner(DLIList<CubitSimpleAttrib*> &csa_list)
{
  
    //- if owner is to be split, get simple attribs for new entity

    // first auto create attributes
  /*CubitStatus success = */
  CGMApp::instance()->attrib_manager()->auto_update_attribs(CAST_TO(this, RefEntity));

    // call split for each attribute, getting back new csa's
  DLIList<CubitAttrib*> ca_list;
  get_cubit_attrib_list(ca_list);
  int i;
  for ( i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get_and_step();
    CubitSimpleAttrib *csa_ptr = ca_ptr->split_owner();
    if (csa_ptr != NULL) csa_list.append(csa_ptr);
  }

    // now, check delete flag for each ca_ptr, and delete if necessary
  for (i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get_and_step();
    if (ca_ptr->delete_attrib()) {
      remove_cubit_attrib(ca_ptr);
      delete ca_ptr;
    }
  }

    // ok, we're done
}

void CubitAttribUser::merge_owner(RefEntity *deletable_entity)
{
    //- if owner is to be merged, combine attribs from deletable_entity with this

    // if no deletable entity, keep the entity name and delete the other
    // attribs on this entity
  int i;
  if (deletable_entity == NULL) {
    DLIList<CubitAttrib*> ca_list;
    get_cubit_attrib_list(ca_list);
    for (i = ca_list.size(); i > 0; i--) {
      CubitAttrib *ca_ptr = ca_list.get_and_step();
      if (ca_ptr->int_attrib_type() != CA_ENTITY_NAME &&
          ca_ptr->int_attrib_type() != CA_ASSEMBLY_DATA)
      {
        remove_cubit_attrib(ca_ptr);
        delete ca_ptr;
      }
    }

    return;
  }
    
    // first auto create attributes
  /*CubitStatus success = */
  CGMApp::instance()->attrib_manager()->auto_update_attribs(CAST_TO(this, RefEntity));

    // copy any attr's on deletable_entity but not on this
  auto_create_for_merge(deletable_entity);

    // call merge for each attribute, passing other as argument
  DLIList<CubitAttrib*> ca_list, deletable_ca_list;
  get_cubit_attrib_list(ca_list);
  deletable_entity->get_cubit_attrib_list(deletable_ca_list);
  assert(ca_list.size() >= deletable_ca_list.size());
  for (i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get_and_step();
    deletable_ca_list.reset();
    CubitAttrib *deletable_ca_ptr = NULL;
    if (deletable_ca_list.size() > 0) {
      deletable_ca_ptr = deletable_ca_list.get();
        // get the corresponding deletable attribute, then extract it from
        //   the list
      while (ca_ptr->int_attrib_type() != deletable_ca_ptr->int_attrib_type() &&
             !deletable_ca_list.is_at_end()) {
        deletable_ca_list.step();
        deletable_ca_ptr = deletable_ca_list.get();
      }
      if (ca_ptr->int_attrib_type() == deletable_ca_ptr->int_attrib_type())
        deletable_ca_list.extract();
      else deletable_ca_ptr = NULL;
    }
    ca_ptr->merge_owner(deletable_ca_ptr);
  }
  
    // now, check delete flag for each ca_ptr, and delete if necessary
  for (i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get_and_step();
    if (ca_ptr->delete_attrib()) {
      remove_cubit_attrib(ca_ptr);
      delete ca_ptr;
    }
  }

  // ok, we're done
}

void CubitAttribUser::transf_owner(const CubitVector &matrow1,
                                   const CubitVector &matrow2,
                                   const CubitVector &matrow3,
                                   const CubitVector &translate_vec,
                                   const double scale_factor)
{
  
    //- called if owner is to be transformed; simply passes information
    //- to attribs on this entity; does *not* autocreate
  
  DLIList<CubitAttrib*> ca_list;
  get_cubit_attrib_list(ca_list);
  for (int i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get_and_step();
    ca_ptr->transf_owner(matrow1, matrow2, matrow3,
                         translate_vec, scale_factor);
  }

    // ok, we're done
}

void CubitAttribUser::auto_create_for_merge(RefEntity *deletable_entity)
{
    // copy any attr's on deletable_entity but not on this
  DLIList<CubitAttrib*> deletable_ca_list;
  deletable_entity->get_cubit_attrib_list(deletable_ca_list);
  if (deletable_ca_list.size() == 0) return;
  
  DLIList<CubitAttrib*> new_list;
  RefEntity *entity = CAST_TO(this, RefEntity);
  for (int i = deletable_ca_list.size(); i > 0; i--) {
    new_list.clean_out();
    CubitAttrib *deletable_ca_ptr = deletable_ca_list.get_and_step();
    find_cubit_attrib_type(deletable_ca_ptr->int_attrib_type(), new_list);
    if (new_list.size() == 0)
    {
      CGMApp::instance()->attrib_manager()->create_cubit_attrib(
                                                deletable_ca_ptr->int_attrib_type(),
                                                entity, NULL);
    }
  }
}

CubitStatus CubitAttribUser::actuate_cubit_attrib (CubitAttrib* cubit_attrib_ptr)
{
  CubitStatus return_value = CUBIT_FAILURE;
  if (cubit_attrib_ptr != NULL )
  {
    if (cubit_attrib_ptr->has_actuated() == CUBIT_TRUE ||
        cubit_attrib_ptr->delete_attrib() == CUBIT_TRUE) {
      return CUBIT_FAILURE;
    }
    
    return_value = cubit_attrib_ptr->actuate();
    if (cubit_attrib_ptr->delete_attrib() == CUBIT_TRUE)
    {
      return_value = remove_cubit_attrib(cubit_attrib_ptr);
      delete cubit_attrib_ptr;
    }
  }
  return return_value;
}

CubitStatus CubitAttribUser::actuate_cubit_attrib (int attrib_type)
{
  DLIList<CubitAttrib*> attrib_list;
  find_cubit_attrib_type(attrib_type, attrib_list);
  CubitStatus actuate_status = actuate_cubit_attrib(attrib_list);
  return actuate_status;
}

CubitStatus CubitAttribUser::actuate_cubit_attrib(DLIList<CubitAttrib*>
                                                      attrib_list)
{
  CubitStatus actuate_status = CUBIT_SUCCESS;
  attrib_list.reset();
  for(int i = 0; i < attrib_list.size(); i++)
  {
    CubitAttrib* cubit_attrib_ptr = attrib_list.get_and_step();
    if (actuate_cubit_attrib(cubit_attrib_ptr) == CUBIT_FAILURE)
      actuate_status = CUBIT_FAILURE;
  }
  return actuate_status;
} 

CubitStatus CubitAttribUser::actuate_cubit_attrib(DLIList<RefEntity*> refent_list,
                                                  int attrib_type)
                                                  
{
  CubitStatus actuate_status = CUBIT_SUCCESS;
  DLIList<CubitAttrib*> attrib_list;
  RefEntity *ref_ent;
  for(int i = refent_list.size(); i > 0; i--)
  {
    ref_ent = refent_list.get_and_step();
    ModelEntity* me_ptr = dynamic_cast<ModelEntity*>(ref_ent);
    if( me_ptr && me_ptr->deactivated() )
      continue;
    ref_ent->find_cubit_attrib_type(attrib_type, attrib_list);
    if (attrib_list.size() > 0)
      return ((attrib_list.get())->actuate_list(refent_list));
  }
  return actuate_status;
}

//CubitStatus CubitAttribUser::actuate_cubit_attrib ()
//{
//  CubitStatus actuate_status = CUBIT_SUCCESS;
//  DLIList<CubitAttrib*> attrib_list;
//  CubitAttrib *attrib;
//  get_cubit_attrib_list(attrib_list);
//  attrib_list.reset();
//
//  for(int i = attrib_list.size(); i != 0; i--)
//  {
//    attrib = attrib_list.get_and_step();
//      if (attrib->actuate() == CUBIT_FAILURE)
//      {
//        actuate_status = CUBIT_FAILURE;
//      }
//      if (attrib->delete_attrib() == CUBIT_TRUE)
//      {
//        remove_cubit_attrib(attrib);
//        delete attrib;
//      }
//  }
//  return actuate_status;
//}

CubitStatus CubitAttribUser::auto_actuate_cubit_attrib (CubitBoolean from_constructor,
                                                        CubitBoolean after_geom_changes)
{
  CubitStatus actuate_status = CUBIT_SUCCESS;
  CubitAttrib *attrib;
  DLIList<CubitAttrib*> attrib_list;
  get_cubit_attrib_list(attrib_list);
  attrib_list.reset();
  int i;
  for( i = attrib_list.size(); i != 0; i--)
  {
    attrib = attrib_list.get_and_step();
      // check first for deletable attribute; this attribute shouldn't really be here,
      // but until we figure out why it is...
    if (attrib->delete_attrib() == CUBIT_TRUE) {
      PRINT_DEBUG_90("Trying to auto actuate deletable attribute - this is bad...\n");
    }
    else if (attrib->auto_actuate_flag() == CUBIT_TRUE &&
        !attrib->has_actuated() &&
        (!from_constructor || attrib->actuate_in_constructor()) &&
        (after_geom_changes || !attrib->actuate_after_geom_changes()))
    {

      PRINT_DEBUG_90("Actuating attribute type %s for %s %d\n",
                     attrib->att_internal_name(), attrib->attrib_owner()->class_name(),
                     attrib->attrib_owner()->id());
      
      if (attrib->actuate() == CUBIT_FAILURE)
      {
        actuate_status = CUBIT_FAILURE;
      }

        // need to check again for delete flag, since it might have been set
        // in actuate function
      if (attrib->delete_attrib() == CUBIT_TRUE)
      {
        remove_cubit_attrib(attrib);
        delete attrib;
      }
    }
  }

//  if (CADeferredAttrib::cleanup_cadas(from_constructor, after_geom_changes) == CUBIT_FAILURE)
//    actuate_status = CUBIT_FAILURE;

  return actuate_status;
}

//CubitStatus CubitAttribUser::update_cubit_attrib (CubitAttrib* cubit_attrib_ptr)
//{
//  if (cubit_attrib_ptr != NULL)
//      return cubit_attrib_ptr->update();
//  return CUBIT_FAILURE;
//}

CubitStatus CubitAttribUser::update_cubit_attrib (int attrib_type)
{
  DLIList<CubitAttrib*> attrib_list;
  CubitAttrib* new_attrib = NULL;
  find_cubit_attrib_type(attrib_type, attrib_list);
  assert(attrib_list.size() <= 1);
  
  if(attrib_list.size() == 0)
    new_attrib = get_cubit_attrib(attrib_type);
  else
    new_attrib = attrib_list.get();

  assert(new_attrib != 0);
  CubitStatus update_status = new_attrib->update();
  
  return update_status;
}

//CubitStatus CubitAttribUser::update_cubit_attrib(DLIList<CubitAttrib*>
//                                                      attrib_list)
//{
//  CubitStatus update_status = CUBIT_SUCCESS;
//  attrib_list.reset();
//  for(int i = attrib_list.size(); i > 0; i--)
//  {
//    CubitAttrib* cubit_attrib_ptr = attrib_list.get_and_step();
//    if (cubit_attrib_ptr->update() == CUBIT_FAILURE)
//      update_status = CUBIT_FAILURE;
//  }
//  return update_status;
//} 

//CubitStatus CubitAttribUser::update_cubit_attrib()
//{
//  return CUBIT_FAILURE;
//} 

CubitStatus CubitAttribUser::auto_update_cubit_attrib ()
{
    // for this cau, automatically create and update ca's

    // first, create ca's for any attribute type which has its auto
    // update flag set and which isn't present yet on this entity
  /*CubitStatus success = */
  CGMApp::instance()->attrib_manager()->auto_update_attribs(CAST_TO(this, RefEntity));

    // now, update all attributes present on this entity
  CubitStatus update_status = CUBIT_SUCCESS;
  DLIList<CubitAttrib*> attrib_list;
  CubitAttrib *attrib;
  get_cubit_attrib_list(attrib_list);
  attrib_list.reset();
  int i;
  for( i = attrib_list.size(); i != 0; i--)
  {
    attrib = attrib_list.get_and_step();
    if (!attrib->has_updated()) {
        // if this attribute has written already, reset the information in it
        // so it gets a "clean" update (otherwise information can be added to
        // lists on the attrib more than once)
      if (CUBIT_TRUE == attrib->has_written())
         attrib->reset();
    
        // reset the delete flag here, so we don't need to do it in every attribute
        // (it'll get set back to delete if the update isn't successful)
      attrib->delete_attrib(CUBIT_FALSE);
      PRINT_DEBUG_90("Updating attribute type %s for %s %d, delete = ",
                     attrib->att_internal_name(), attrib->attrib_owner()->class_name(),
                     attrib->attrib_owner()->id());
      update_status = attrib->update();

      PRINT_DEBUG_90("%s\n",
                     (attrib->delete_attrib() == CUBIT_FALSE ? "NO" : "YES"));

        // if update was successful, reset the written flag (assumes that all
        // updates are done before writing starts)
        // (don't reset if it's the entity_name attrib, since this attrib is
        // written automatically every time it's updated)
      if (update_status == CUBIT_SUCCESS && attrib->int_attrib_type() != CA_ENTITY_NAME)
        attrib->has_written(CUBIT_FALSE);
    }
  }
/*
  RefEntity* re_ptr;
  re_ptr = CAST_TO(this, RefEntity);
  DLIList<RefEntity*> children;
  re_ptr->get_child_ref_entities( children );
  for ( i = children.size(); i > 0; i-- )
  {
    children.get()->auto_update_cubit_attrib();
    children.step();
  }
*/
  return update_status;
}

void CubitAttribUser::auto_reset_cubit_attrib (DLIList<RefEntity*> ref_ents)
{
    // set the update flag back off for all attribs on these entities and their children
  DLIList<RefEntity*> children, temp_list;
  int i;
  for (i = ref_ents.size(); i > 0; i--) {
    temp_list.clean_out();
    ref_ents.get_and_step()->get_all_child_ref_entities(temp_list);
    children += temp_list;
  }
  ref_ents.merge_unique(children);
  
    // ok, have a unique'd list if entities; now reset on each of them
  ref_ents.reset();
  for (i = ref_ents.size(); i > 0; i--) {
    ref_ents.get()->set_written_flag(CUBIT_FALSE);
    ref_ents.get_and_step()->set_updated_flag(CUBIT_FALSE);
  }
}

//void CubitAttribUser::auto_reset_cubit_attrib ()
//{
//    // set the update flag back off for all attribs on these entities and their children
//  DLIList<CubitAttrib*> attrib_list;
//  CubitAttrib *attrib;
//  get_cubit_attrib_list(attrib_list);
//  attrib_list.reset();
//  int i;
//  for( i = attrib_list.size(); i != 0; i--)
//  {
//    attrib = attrib_list.get_and_step();
//    if (attrib->has_updated() && attrib->has_written())
//      attrib->has_updated(CUBIT_FALSE);
//  }
//}

CubitStatus CubitAttribUser::auto_update_cubit_attrib (DLIList<RefEntity*> &entity_list,
                                                       CubitBoolean write_too)
{
    //- for entity_list, auto create, update and write ca's

    // now, reset the update and write flags for these entities
    // and their children
  auto_reset_cubit_attrib(entity_list);

    // need to update all then write all in separate loops,
    // to prevent duplication of attributes
  int i;
  for ( i = entity_list.size(); i > 0; i--) 
     entity_list.get_and_step()->auto_update_cubit_attrib();

  if (write_too) {
    for (i = entity_list.size(); i > 0; i--) 
      entity_list.get_and_step()->write_cubit_attribs();
  }
  
  
  return CUBIT_SUCCESS;
}

CubitStatus CubitAttribUser::clear_all_simple_attrib( DLIList<RefEntity*>& entity_list )
{
  CubitStatus result = CUBIT_SUCCESS;
  for( int i = entity_list.size(); i--; )
    if( entity_list.get_and_step()->clear_simple_attribs() != CUBIT_SUCCESS )
      result = CUBIT_FAILURE;
  return result;
}


void CubitAttribUser::find_cubit_attrib_type (int type,
                                              DLIList<CubitAttrib*>& attrib_list) const
{
  for(CubitAttrib* cubit_attrib_ptr = headAttrib;
      cubit_attrib_ptr != NULL;
      cubit_attrib_ptr = cubit_attrib_ptr->next_attrib())
  {
    if (type == cubit_attrib_ptr->int_attrib_type())
      attrib_list.append_unique(cubit_attrib_ptr);
  }
}

CubitStatus CubitAttribUser::remove_cubit_attrib(DLIList<CubitAttrib*>
                                                      attrib_list)
{
  CubitStatus remove_status = CUBIT_SUCCESS;
  attrib_list.reset();
  for(int i = attrib_list.size(); i != 0; i--)
  {
    CubitAttrib* cubit_attrib_ptr = attrib_list.get_and_step();
    if (remove_cubit_attrib(cubit_attrib_ptr) == CUBIT_FAILURE)
      remove_status = CUBIT_FAILURE;
    delete cubit_attrib_ptr;
  }
  return remove_status;
} 

CubitStatus CubitAttribUser::remove_cubit_attrib(int attrib_type)
{
  DLIList<CubitAttrib*> attrib_list;
  find_cubit_attrib_type(attrib_type, attrib_list);
  CubitStatus remove_status = remove_cubit_attrib(attrib_list);
  return remove_status;
} 

CubitStatus CubitAttribUser::remove_cubit_attrib (CubitAttrib*
                                                  cubit_attrib_ptr)
{
  if (cubit_attrib_ptr == NULL)
    return CUBIT_FAILURE;
  CubitStatus remove_geom = remove_attrib_geometry_entity(cubit_attrib_ptr);
  //CubitStatus remove_cubit = CUBIT_FAILURE;
  CubitBoolean once = CUBIT_FALSE;
  if (headAttrib == cubit_attrib_ptr) {
    headAttrib = cubit_attrib_ptr->next_attrib();
    once = CUBIT_TRUE;
  }
  
  CubitAttrib* temp_cubit_attrib_ptr = headAttrib;
  while (temp_cubit_attrib_ptr != NULL) {
    if (cubit_attrib_ptr == temp_cubit_attrib_ptr->next_attrib()) {
      if (once) 
        PRINT_DEBUG_90("Removing attribute more than once.\n");
      temp_cubit_attrib_ptr->set_next_attrib(cubit_attrib_ptr->
                                             next_attrib());
      once = CUBIT_TRUE;
    }
      
    temp_cubit_attrib_ptr = temp_cubit_attrib_ptr->next_attrib();
  }

  return remove_geom;
}

//CubitStatus CubitAttribUser::remove_cubit_attrib ()
//{
//  CubitAttrib *cubit_attrib_ptr = NULL;
//  for(cubit_attrib_ptr = headAttrib;
//      cubit_attrib_ptr != NULL;)
//  {
//    headAttrib = cubit_attrib_ptr->next_attrib();
//    delete cubit_attrib_ptr;
//    cubit_attrib_ptr = headAttrib;
//  }
//  remove_attrib_geometry_entity();
//  return CUBIT_SUCCESS;
//}
  
  

CubitStatus CubitAttribUser::auto_read_cubit_attrib()
{
    // auto read all simple attributes on this entity & create CA's for them;
    // checks global and CA-specific auto read flag, and puts CSA's back on
    // entity if auto read for that type is off
  
  CubitStatus cubit_assign_status = CUBIT_FAILURE;
  
    // Get the GeometryEntity of this RefEntity (it could be the OSME
    // if it is a Body) and get its name, if it exists
  
  DLIList<CubitSimpleAttrib*> csattrib_list;

    // Deal with a Body entity
  Body* Body_ptr = CAST_TO(this, Body);
  BasicTopologyEntity* BTE_ptr = CAST_TO(this, BasicTopologyEntity);
  if ( Body_ptr != NULL )
  {
    BodySM* OSME_ptr = Body_ptr->get_body_sm_ptr();
    OSME_ptr->get_simple_attribute(csattrib_list);
    remove_all_simple_attribute(OSME_ptr);
  }
  else if ( BTE_ptr != NULL )
  {
    GeometryEntity* GE_ptr = BTE_ptr->get_geometry_entity_ptr();
    GE_ptr->get_simple_attribute(csattrib_list);
    remove_all_simple_attribute(GE_ptr);
  }
  else return CUBIT_SUCCESS;
  
  csattrib_list.reset();
  for(int i = csattrib_list.size(); i != 0; i--)
  {
    CubitSimpleAttrib* cubit_simple_attrib_ptr = csattrib_list.get_and_step();
    int csa_type = CGMApp::instance()->attrib_manager()->attrib_type(cubit_simple_attrib_ptr);
    CubitAttrib *new_attrib = NULL;
    if (CGMApp::instance()->attrib_manager()->auto_read_flag(csa_type))
        // auto create this CA
      new_attrib = CGMApp::instance()->attrib_manager()->create_cubit_attrib(csa_type,
                                                    CAST_TO(this, RefEntity), cubit_simple_attrib_ptr);
    
    if (new_attrib != NULL && new_attrib->delete_attrib() == CUBIT_FALSE) 
      cubit_assign_status = CUBIT_SUCCESS;
    else if (new_attrib != NULL && new_attrib->delete_attrib() == CUBIT_TRUE) {
      remove_cubit_attrib(new_attrib);
      delete new_attrib;
    }
    else
      put_simple_attrib(cubit_simple_attrib_ptr);

    delete cubit_simple_attrib_ptr;
  }
  
  return cubit_assign_status;
}

CubitStatus CubitAttribUser::read_cubit_attrib(int attrib_type)
{

  CubitStatus read_status;
  
    // get all simple attrib's
  DLIList<CubitSimpleAttrib*> csattrib_list;

    // Deal with a Body entity
  Body* Body_ptr = CAST_TO(this, Body);
  BasicTopologyEntity* BTE_ptr = CAST_TO(this, BasicTopologyEntity);
  if ( Body_ptr != NULL )
  {
    BodySM* OSME_ptr = Body_ptr->get_body_sm_ptr();
    read_status = OSME_ptr->get_simple_attribute(csattrib_list);
    remove_all_simple_attribute(OSME_ptr);
  }
  else if ( BTE_ptr != NULL )
  {
    GeometryEntity* GE_ptr = BTE_ptr->get_geometry_entity_ptr();
    read_status = GE_ptr->get_simple_attribute(csattrib_list);
    remove_all_simple_attribute(GE_ptr);
  }
  else return CUBIT_SUCCESS;
  
  csattrib_list.reset();
  for(int i = csattrib_list.size(); i != 0; i--)
  {
    CubitSimpleAttrib* cubit_simple_attrib_ptr = csattrib_list.get_and_step();
    int csa_type = CGMApp::instance()->attrib_manager()->attrib_type(cubit_simple_attrib_ptr);
    if (attrib_type == CA_ALL_ATTRIBUTES || csa_type == attrib_type) {
        // create this CA
//       CubitAttrib *new_attrib =
      CGMApp::instance()->attrib_manager()->create_cubit_attrib(csa_type,
                                                                CAST_TO(this, RefEntity), cubit_simple_attrib_ptr);
    }
    else {
        // otherwise we don't want to read this one; since the get_simple_attribute
        // took the attribute off, we'll have to put this one back on
      put_simple_attrib(cubit_simple_attrib_ptr);
    }
    delete cubit_simple_attrib_ptr;
  }
      
  return read_status;
}
  
//CubitStatus CubitAttribUser::read_cubit_attrib(CubitBoolean read_children)
//{
//  CubitStatus read_status = CUBIT_SUCCESS;
//  if (read_cubit_attrib(CA_ALL_ATTRIBUTES) == CUBIT_FAILURE)
//      read_status = CUBIT_FAILURE;
//
//  if (read_children) {
//      // now, call read for all children
//    DLIList<RefEntity*> children;
//    RefEntity* re_ptr;
//    re_ptr = CAST_TO(this, RefEntity);
//    re_ptr->get_all_child_ref_entities( children );
//    CubitStatus temp_status;
//    for ( int i = children.size(); i > 0; i-- )
//    {
//      temp_status = children.get()->read_cubit_attrib(CA_ALL_ATTRIBUTES);
//      if (temp_status == CUBIT_FAILURE) read_status = CUBIT_FAILURE;
//      children.step();
//    }
//  }
//  
//  return read_status;
//}
 
void CubitAttribUser::get_cubit_attrib_list (DLIList<CubitAttrib*>& attrib_list)
{
  for(CubitAttrib* cubit_attrib_ptr = headAttrib;
      cubit_attrib_ptr != NULL;
      cubit_attrib_ptr = cubit_attrib_ptr->next_attrib()) {
    assert(NULL != cubit_attrib_ptr);
    attrib_list.append_unique(cubit_attrib_ptr);
  }
}

int CubitAttribUser::num_cubit_attrib()
{
  int number = 0;
  for(CubitAttrib* cubit_attrib_ptr = headAttrib;
      cubit_attrib_ptr != NULL;
      cubit_attrib_ptr = cubit_attrib_ptr->next_attrib())
    number++;

  return number;
}

CubitStatus CubitAttribUser::remove_attrib_geometry_entity (CubitAttrib*
                                                            cubit_attrib_ptr)
{
  CubitStatus removed = CUBIT_FAILURE;
  CubitSimpleAttrib *csattrib_ptr;
  
  if (cubit_attrib_ptr != NULL &&
      (csattrib_ptr = cubit_attrib_ptr->cubit_simple_attrib()) != NULL)
  {

    if (DEBUG_FLAG(90)) {
      RefEntity *ref_entity = CAST_TO(this, RefEntity);
      PRINT_DEBUG_90( "Removing simple attribute of type %s on"
                  " %s %d.\n", csattrib_ptr->character_type().c_str(),
                  ref_entity->class_name(), ref_entity->id());
    }
    
    TopologyEntity* topo_ptr = CAST_TO(this, TopologyEntity);
    if (topo_ptr && topo_ptr->bridge_manager()->topology_bridge())
    {
        remove_simple_attribute(topo_ptr->bridge_manager()->topology_bridge(), csattrib_ptr);
      removed = CUBIT_SUCCESS;
    }

    delete csattrib_ptr;
  }
  return removed;
}

//CubitStatus CubitAttribUser::remove_attrib_geometry_entity ()
//{
//  CubitStatus removed = CUBIT_FAILURE;
//  Body* Body_ptr = CAST_TO(this, Body);
//    
//  BasicTopologyEntity* BTE_ptr = CAST_TO(this, BasicTopologyEntity);
//      //Deal with Bodies
//  if ( Body_ptr != NULL )
//  {
//      // Get the OSME pointer
//    BodySM* OSME_ptr = Body_ptr->get_body_sm_ptr();
//    remove_all_simple_attribute(OSME_ptr);
//    removed = CUBIT_SUCCESS;
//  }
//  
//    // Deal with BasicTopologyEntities
//  else if ( BTE_ptr != NULL )
//  {
//      // Get the GeometryEntity pointer
//    GeometryEntity* GE_ptr = 
//      BTE_ptr->get_geometry_entity_ptr();
//    remove_all_simple_attribute(GE_ptr);
//    removed = CUBIT_SUCCESS;
//  }
//
//    // All other Entities
//  else
//  {
//  }
//  return removed;
//}

void CubitAttribUser::set_written_flag(CubitBoolean flag)
{
  DLIList<CubitAttrib*> ca_list;
  get_cubit_attrib_list(ca_list);
  for (int i = ca_list.size(); i > 0; i--) {
    if (flag == CUBIT_FALSE && ca_list.get()->has_written() == CUBIT_TRUE)
      ca_list.get()->reset();
    ca_list.get_and_step()->has_written(flag);
  }
}

void CubitAttribUser::set_updated_flag(CubitBoolean flag)
{
  DLIList<CubitAttrib*> ca_list;
  get_cubit_attrib_list(ca_list);
  for (int i = ca_list.size(); i > 0; i--)
    ca_list.get_and_step()->has_updated(flag);
}

//void CubitAttribUser::print_attribs() 
//{
//  PRINT_INFO("Attributes on %s %d:\n", 
//             CAST_TO(this, RefEntity)->entity_name().c_str(),
//             CAST_TO(this, RefEntity)->id());
//  DLIList<CubitAttrib*> attrib_list;
//  get_cubit_attrib_list(attrib_list);
//  int i;
//  for (i = attrib_list.size(); i > 0; i--)
//    attrib_list.get_and_step()->print();
//  
//  if (attrib_list.size() == 0) PRINT_INFO("(none)\n");
//}

void CubitAttribUser::append_simple_attribute(TopologyBridge *bridge, CubitSimpleAttrib* attrib_ptr)
{
    // for merged objects, put attribute on other entities
  if (CubitSimpleAttrib::get_push_attribs() == CUBIT_TRUE ||
    CGMApp::instance()->attrib_manager()->auto_update_flag(CA_MERGE_PARTNER) == CUBIT_TRUE) {
  
    DLIList<TopologyBridge*> tb_list;
    bridge->bridge_manager()->get_bridge_list(tb_list);

      // if this entity is merged, it will have > 1 entity in the bridge
      // list
    tb_list.reset();
    assert(tb_list.size() == 0 || tb_list.get() == bridge);

      // Special handling of MergePartner attribute.  
      // Need to store the bridge sense in the attribute, which
      // is potentially different for different TopologyBridges.
    int type = CGMApp::instance()->attrib_manager()->attrib_type_from_internal_name(attrib_ptr->character_type().c_str());
    if ( type == CA_MERGE_PARTNER )
    {
        // now adjust the attrib on each bridge to hold the 
        // relative sense for that bridge.
      for (int i = tb_list.size(); i > 0; i--) {
        TopologyBridge *temp_tb = tb_list.get_and_step();

        CAMergePartner::set_bridge_sense( attrib_ptr, temp_tb->bridge_sense() );

        GeometryEntity *geom_ptr = dynamic_cast<GeometryEntity*>(temp_tb);
        if (geom_ptr ) 
        {
          //if we're copying a merged entity, saved id should be zero
          if( GeometryModifyTool::instance()->get_copy_entity() ) 
            CAMergePartner::set_saved_id( attrib_ptr, 0 );
          else
            CAMergePartner::set_saved_id( attrib_ptr, geom_ptr->get_saved_id() );

          //First bridge should be marked as "survivor"
          if( i == tb_list.size() )
          {
            CAMergePartner::set_survivor( attrib_ptr, 1 );
          }
          else
            CAMergePartner::set_survivor( attrib_ptr, 0 );
        }
        append_attrib_internal(temp_tb, attrib_ptr);
      }
    }
      // For anything other than CAMergePartner, just append the
      // unmodified attribute to each bridge.
    else
    {
      for ( int i = tb_list.size(); i--; )
        append_attrib_internal( tb_list.get_and_step(), attrib_ptr);
    }
  }
  else
  {
      // Append this name to the primary bridge
    append_attrib_internal(bridge, attrib_ptr);  
  }
}

void CubitAttribUser::append_attrib_internal(TopologyBridge *bridge, CubitSimpleAttrib* attrib_ptr)
{
  DLIList<CubitSimpleAttrib*> others(1);
  
    // Check for duplicates
  if ( attrib_ptr->character_type() != "DEFERRED_ATTRIB" )
  {
    bridge->get_simple_attribute( attrib_ptr->character_type().c_str(), others );
    while ( others.size() )
    {
      CubitSimpleAttrib* dup_attrib = others.pop();
      bridge->remove_simple_attribute_virt( dup_attrib );
      delete dup_attrib;
    }
  }
  else
  {
    attrib_ptr->string_data_list()->reset();
    CubitString real_name = *attrib_ptr->string_data_list()->next();
    bridge->get_simple_attribute("DEFERRED_ATTRIB", others);
    while ( others.size() )
    {
      CubitSimpleAttrib* dup_attrib = others.pop();
      dup_attrib->string_data_list()->reset();
      if ( *dup_attrib->string_data_list()->next() == real_name )
        bridge->remove_simple_attribute_virt(dup_attrib);
      delete dup_attrib;
    }
  }
  
    // append attribute
  bridge->append_simple_attribute_virt(attrib_ptr);  
}

void CubitAttribUser::remove_all_simple_attribute(TopologyBridge* bridge)
{
    // remove this name from the primary object
  bridge->remove_all_simple_attribute_virt();

    // for merged objects, remove all attributes from other entities
  if (CubitSimpleAttrib::get_push_attribs() == CUBIT_TRUE ||
      CGMApp::instance()->attrib_manager()->auto_update_flag(CA_MERGE_PARTNER) == CUBIT_TRUE) {
    DLIList<TopologyBridge*> tb_list;
    bridge->bridge_manager()->get_bridge_list(tb_list);
    tb_list.reset();
    assert(tb_list.size() == 0 || tb_list.get() == bridge);
    for (int i = tb_list.size(); i > 0; i--) {
      TopologyBridge *temp_tb = tb_list.get_and_step();
      if (temp_tb != bridge) temp_tb->remove_all_simple_attribute_virt();
    }
  }
}

void CubitAttribUser::remove_simple_attribute(TopologyBridge* bridge, CubitSimpleAttrib* attrib_ptr)
{
    // remove this name from the primary object
  bridge->remove_simple_attribute_virt(attrib_ptr);
  
    // for merged objects, remove attribute from other entities
  if (CubitSimpleAttrib::get_push_attribs() == CUBIT_TRUE ||
      CGMApp::instance()->attrib_manager()->auto_update_flag(CA_MERGE_PARTNER) == CUBIT_TRUE) {
    DLIList<TopologyBridge*> tb_list;
    BridgeManager *bm_ptr = bridge->bridge_manager();
    if (bm_ptr != NULL) bm_ptr->get_bridge_list(tb_list);
    tb_list.reset();
    assert(tb_list.size() == 0 || tb_list.get() == bridge);
    for (int i = tb_list.size(); i > 0; i--) {
      TopologyBridge *temp_tb = tb_list.get_and_step();
      if (temp_tb != bridge) temp_tb->remove_simple_attribute_virt(attrib_ptr);
    }
  }
}

