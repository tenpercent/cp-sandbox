//- Class:          CAEntityName
//- Owner:          Greg Nielson
//- Description:    Cubit Attribute for entity names.
//- Checked By:
//- Version:

#include "CAEntityName.hpp"
#include "BasicTopologyEntity.hpp"
#include "Body.hpp"
#include "RefEntityName.hpp"
#include "CastTo.hpp"

CubitAttrib* CAEntityName_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAEntityName *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAEntityName(entity);
  }
  else
  {
    new_attrib = new CAEntityName(entity, p_csa);
  }

  return new_attrib;
}

CAEntityName::CAEntityName(RefEntity* new_attrib_owner,
                           CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  PRINT_DEBUG_94( "Creating ENTITY_NAME attribute from CSA for %s %d\n",
              (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
              (attribOwnerEntity ? attribOwnerEntity->id() : 0));
  
  DLIList<CubitString*> *cs_list = csa_ptr->string_data_list();
  cs_list->reset();

    // step over the attribute type
  cs_list->step();

    // now read name / option pairs
  for (int i = cs_list->size()-1; i > 0; i--)
  {
    CubitString *cs = cs_list->get_and_step();
    if (*cs == CubitString(""))
      PRINT_WARNING("Empty name attribute for %s %d.\n",
                    attribOwnerEntity->class_name(),
                    attribOwnerEntity->id());
    else
      entityNames.append(new CubitString(*cs));
  }

  if (entityNames.size() == 0)
    deleteAttrib = CUBIT_TRUE;
}

CAEntityName::CAEntityName(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  PRINT_DEBUG_94( "Creating ENTITY_NAME attribute for %s %d\n",
              (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
              (attribOwnerEntity ? attribOwnerEntity->id() : 0));
}

CAEntityName::~CAEntityName()
{
  PRINT_DEBUG_94("Deleting ENTITY_NAME attribute\n");
  for (int i = entityNames.size(); i > 0; i--) {
    delete entityNames.get_and_step();
  }
  
  entityNames.clean_out();
}

CubitStatus CAEntityName::actuate()
{
  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;

  PRINT_DEBUG_94( "Actuating ENTITY_NAME attribute for %s %d\n",
              attribOwnerEntity->class_name(), attribOwnerEntity->id());

  entityNames.reset();

  CubitBoolean update_attribs = (CubitBoolean) 
      (RefEntityName::instance()->get_generate_default_names()!=1);
  
  RefEntityName::instance()->add_refentity_name(attribOwnerEntity,entityNames,
                                                update_attribs);

  hasActuated = CUBIT_TRUE;
  return CUBIT_SUCCESS;
}

CubitStatus CAEntityName::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
  PRINT_DEBUG_94( "Updating ENTITY_NAME attribute for %s %d\n",
              attribOwnerEntity->class_name(), attribOwnerEntity->id());

    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // first, remove this attrib in its old form from the geometry entity
  if (hasWritten == CUBIT_TRUE) {
    attribOwnerEntity->remove_attrib_geometry_entity(this);
    hasWritten = CUBIT_FALSE;
  }
  
  DLIList<CubitString*> names;
  int num_names = RefEntityName::instance()->
      get_refentity_name(attribOwnerEntity,names);
  if( num_names == 0)
  {
    delete_attrib(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      // get the names; since RefEntityName passes back pointers to
      // its own strings, create new ones here
    int i;

      // first, delete all the old names on the list
    for (i = entityNames.size(); i > 0; i--)
      delete entityNames.get_and_step();
    entityNames.clean_out();
    
    for (i = names.size(); i > 0; i--) {
      entityNames.append(new CubitString(*names.get_and_step()));
    }
    
      // reset the delete flag if it was set before
    delete_attrib(CUBIT_FALSE);
    
      // now, write to geometry entity
    attribOwnerEntity->write_specific_cubit_attrib(this);
  }
  return CUBIT_SUCCESS;
}

CubitStatus CAEntityName::reset()
{
  PRINT_DEBUG_94("CAEntityName::reset()\n");
  
    //- reset function, cleans out name list
  int i;
  for ( i = entityNames.size(); i > 0; i--)
    delete entityNames.get_and_step();

  entityNames.clean_out();

    // need to reset hasUpdated flag too, so next update will do something
  hasUpdated = CUBIT_FALSE;
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib *CAEntityName::split_owner()
{
    // if this entity is to be split, pass back a simple attribute with
    // duplicate name data to be put on new entity
  PRINT_DEBUG_94("CAEntityName::split_owner()\n");
  update();
  return cubit_simple_attrib();
}

void CAEntityName::merge_owner(CubitAttrib *deletable_attrib)
{
    // if this entity is to be merged, copy names over from deletable entity
  CAEntityName *caen_ptr = CAST_TO(deletable_attrib, CAEntityName);
  if (caen_ptr)
  {
    DLIList<CubitString*> &other_names = caen_ptr->entityNames;
    other_names.reset();
    for (int i = other_names.size(); i--; )
    {
      entityNames.append(new CubitString(*(other_names.get_and_step())));
    }
  }
}

CubitSimpleAttrib* CAEntityName::cubit_simple_attrib()
{
  PRINT_DEBUG_94("CAEntityName::cubit_simple_attrib()\n");
  
  DLIList<CubitString*> cs_list;

    // pack the string list:
    // character type of this CA
  cs_list.append(new CubitString(att_internal_name()));

    // name, option pairs
  int i;
  for ( i = entityNames.size(); i > 0; i--) {
    cs_list.append(new CubitString(*entityNames.get_and_step()));
  }

  CubitSimpleAttrib* csattrib_ptr =
      new CubitSimpleAttrib(&cs_list, NULL, NULL);

    // since CubitSimpleAttrib constructor creates new names from
    // the list passed in, we should delete the strings created in
    // this function
  for (i = cs_list.size(); i > 0; i--)
    delete cs_list.get_and_step();

  return csattrib_ptr;
}


void CAEntityName::print()
{
    // print info on this attribute
  entityNames.reset();
  
  PRINT_INFO("CAEntityName: owner = %s %d; names: ",
             attribOwnerEntity->class_name(), attribOwnerEntity->id());
  for (int i = entityNames.size(); i > 0; i--)
    PRINT_INFO("%s ", entityNames.get_and_step()->c_str());

  PRINT_INFO("\n");
}
