//- Class:          CAUniqueId
//- Owner:          Tim Tautges
//- Description:    Cubit Attribute for unique ids
//- Checked By:
//- Version:

#include "CAUniqueId.hpp"
#include "CubitSimpleAttrib.hpp"
#include "RefEntity.hpp"
#include "TDUniqueId.hpp"

DLIList<CAUniqueId *> CAUniqueId::allCAUniqueIds;
bool CAUniqueId::autoUniqueId = false;

CubitAttrib* CAUniqueId_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAUniqueId *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAUniqueId(entity);
  }
  else
  {
    new_attrib = new CAUniqueId(entity, p_csa);
  }

  return new_attrib;
}

CAUniqueId::~CAUniqueId()
{
  if (allCAUniqueIds.move_to(this))
    allCAUniqueIds.extract();
}

CAUniqueId::CAUniqueId(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  uniqueId = -1;
  allCAUniqueIds.append(this);
}

CAUniqueId::CAUniqueId(RefEntity* new_attrib_owner,
                               CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  uniqueId = *csa_ptr->int_data_list()->get();
  allCAUniqueIds.append(this);
}

CubitStatus CAUniqueId::actuate()
{

  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;

    // create a TDUniqueId for the entity, if it doesn't already
    // exist
  TDUniqueId *uid = (TDUniqueId *) attrib_owner()->get_TD(&TDUniqueId::is_unique_id);

  if (uid != NULL) {
      // check to make sure it's the same unique id
    if (uid->unique_id() != uniqueId) {
      PRINT_ERROR("Different unique id found for %s %d.\n",
                  attrib_owner()->class_name(), attrib_owner()->id());
      return CUBIT_FAILURE;
    }
  }
  else {
      // else make a new one
    uid = new TDUniqueId(attrib_owner(), uniqueId);
  }
  
  delete_attrib(CUBIT_TRUE);
  hasActuated = CUBIT_TRUE;
  
  return CUBIT_SUCCESS;
}

CubitStatus CAUniqueId::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // if the owner has a unique id, save it, otherwise delete this one
  TDUniqueId *td_uid = (TDUniqueId *) attrib_owner()->get_TD(&TDUniqueId::is_unique_id);

  if (NULL == td_uid && autoUniqueId) 
    td_uid = new TDUniqueId(attrib_owner());

  if (td_uid == NULL) {
    delete_attrib(CUBIT_TRUE);
  }

  else {
    uniqueId = td_uid->unique_id();
    if (delete_attrib() == CUBIT_TRUE) delete_attrib(CUBIT_FALSE);
  }
  
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CAUniqueId::cubit_simple_attrib()
{
  return new CubitSimpleAttrib(att_internal_name(), "", "",
                               uniqueId, 0.0);
}

CubitStatus CAUniqueId::actuate_all()
{
    //- actuate all the CAUI's on the list, then empty the list
  for (int i = allCAUniqueIds.size(); i > 0; i--) {
    CAUniqueId *cauid = allCAUniqueIds.get();
    if (cauid->actuate() == CUBIT_SUCCESS) allCAUniqueIds.extract();
    else allCAUniqueIds.step();
  }
  
  return CUBIT_SUCCESS;
}

void CAUniqueId::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAUniqueId: owner = %s %d: uid=%d\n",
             attribOwnerEntity->class_name(), attribOwnerEntity->id(),
             uniqueId);
}
