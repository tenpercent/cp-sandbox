
#include "CubitAttribManager.hpp"
#include "CubitAttribUser.hpp"
#include "CastTo.hpp"
#include "assert.h"
#include "CADefines.hpp"

#include "RefEntity.hpp"

CubitAttribManager::CubitAttribManager() 
    : silentFlag(false)
{
}

CubitAttribManager::~CubitAttribManager()
{
/*
  std::map<int, CARegisterData*>::iterator iter = mTypeToCAData.begin();

  while (mTypeToCAData.end() != iter)
  {
    delete iter->second;
    iter++;
  }
*/
  int i;
  assert(mTypeNames.size() == mInternalNames.size());


  for (i=mTypeNames.size(); i>0; i--)
  {
    delete [] mTypeNames.get_and_step();
    delete [] mInternalNames.get_and_step();
  }

}

CubitStatus
CubitAttribManager::register_attrib_type( const int att_type,
                                         const char* att_type_name,
                                      const char* att_internal_name,
                                      CACreateFunction p_create_function,
                                      CubitBoolean auto_actuate_flag,
                                      CubitBoolean auto_update_flag,
                                      CubitBoolean auto_write_flag,
                                      CubitBoolean auto_read_flag,
                                      CubitBoolean actuate_in_constructor,
                                      CubitBoolean actuate_after_geom_changes)
{
  assert(CA_UNDEFINED != att_type);
  assert(CA_ALL_ATTRIBUTES != att_type);

/*
  unsigned int flags = 0;
  if (auto_actuate_flag) flags |= CA_AUTOACTUATE;
  if (auto_update_flag) flags |= CA_AUTOUPDATE;
  if (auto_write_flag) flags |= CA_AUTOWRITE;
  if (auto_read_flag) flags |= CA_AUTOREAD;
  if (actuate_in_constructor) flags |= CA_ACTUATEINCONSTRUCTOR;
  if (actuate_after_geom_changes) flags |= CA_ACTUATEAFTERGEOMCHANGES;
  CARegisterData *p_entry = new CARegisterData(att_type_name, att_internal_name,
                                               p_create_function, flags);
  
  std::pair<std::map<int, CARegisterData*>::iterator, bool> result_pair;
  result_pair = mTypeToCAData.insert( std::pair<int, CARegisterData*>(att_type, p_entry) );

  if (!result_pair.second) // failed to insert because the type was already registered
  {
    assert(false);
    delete p_entry;
    return CUBIT_FAILURE;
  }
*/

  int index = mTypes.where_is_item(att_type);
  if (-1 != index) // type already registered
  {
    if (!silentFlag) PRINT_ERROR("Attribute type %d already registered\n.", att_type);
    return CUBIT_FAILURE;
  }

  int size = mTypes.size();
  assert (mTypeNames.size() == size);
  assert (mCreatorFunctions.size() == size);
  assert (mInternalNames.size() == size);
  assert (mAutoActuateFlags.size() == size);
  assert (mAutoUpdateFlags.size() == size);
  assert (mAutoWriteFlags.size() == size);
  assert (mAutoReadFlags.size() == size);
  assert (mActuateInConstructor.size() == size);
  assert (mActuateAfterGeomChanges.size() == size);
//  assert (mAttribImported.size() == size);

  mTypes.append(att_type);

  int namesize = strlen(att_type_name);
  char *stemp = new char[namesize+1];
  strcpy(stemp, att_type_name);
  mTypeNames.append(stemp);

  namesize = strlen(att_internal_name);
  stemp = new char[namesize+1];
  strcpy(stemp, att_internal_name);
  mInternalNames.append(stemp);

  mCreatorFunctions.append(p_create_function);

  mAutoActuateFlags.append(auto_actuate_flag);
  mAutoUpdateFlags.append(auto_update_flag);
  mAutoWriteFlags.append(auto_write_flag);
  mAutoReadFlags.append(auto_read_flag);
  mActuateInConstructor.append(actuate_in_constructor);
  mActuateAfterGeomChanges.append(actuate_after_geom_changes);
//  mAttribImported.append(CUBIT_FALSE); // TODO - get rid of mAttribImported

  return CUBIT_SUCCESS;
}

CubitAttrib*
CubitAttribManager::create_cubit_attrib(const int attrib_type,
                                        RefEntity *attrib_user,
                                        CubitSimpleAttrib *p_csa)
{
  CubitAttrib* new_attrib = NULL;
/*
  CARegisterData *p_entry = NULL;
  CACreateFunction p_creator;
  std::map<int, CARegisterData*>::iterator iter;

  iter = mTypeToCAData.find(attrib_type);
  if (mTypeToCAData.end() == iter)
  {
    assert(false);
    return NULL;
  }

  p_entry = iter->second;
  p_creator = p_entry->mCreateFunction;
  new_attrib = (*p_creator)(attrib_user);
  return new_attrib;
*/

  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) // type not registered
  {
    if (!silentFlag) PRINT_ERROR("Attribute type %d can't be created because it hasn't been registered.\n", 
                attrib_type);
    return NULL;
  }


  CACreateFunction p_creator;
  mCreatorFunctions.reset();
  p_creator = mCreatorFunctions.next(index);

  new_attrib = (*p_creator)(attrib_user, p_csa);
  return new_attrib;

}

/*
CubitStatus CubitAttribManager::actuate_list(DLIList<RefEntity*> entity_list)
{
  int i, j;
  RefEntity * ref_ent;
  for(i = entity_list.size(); i > 0; i--)
  {
    ref_ent = entity_list.get_and_step();

    mTypes.reset();
    for (j = mTypes.size(); j>0; j--)
    {
      ref_ent->actuate_cubit_attrib(mTypes.get_and_step());
    }
  }
  return CUBIT_SUCCESS;
}
*/
CubitStatus CubitAttribManager::auto_update_attribs(RefEntity *cau)
{
    //- create attribs whose auto update flag is set

  CubitStatus status = CUBIT_SUCCESS;
  DLIList<CubitAttrib*> attrib_list;

/*
  std::map<int, CARegisterData*>::iterator iter = mTypeToCAData.begin();

  while (mTypeToCAData.end() != iter)
  {
    CubitBoolean auto_update = ((iter->second)->flags) & CA_AUTOUPDATE;
    int attrib_type = iter->first;

    if (auto_update)
    {
      attrib_list.clean_out();
      cau->find_cubit_attrib_type(attrib_type, attrib_list);
      if (attrib_list.size() == 0) {
        create_cubit_attrib(attrib_type, cau);
        if (status == CUBIT_FAILURE) break;
      }
    }
    iter++;
  }
*/
  mAutoUpdateFlags.reset();
  mTypes.reset();
  assert(mTypes.size() == mAutoUpdateFlags.size());
  int index;
  for (index=mAutoUpdateFlags.size(); index>0; index--) {
      // check the auto update flag first, if not set we can go
    if (mAutoUpdateFlags.get() == CUBIT_TRUE)
    {

        // else we must create a CA of type if there's not one already there
      attrib_list.clean_out();
      cau->find_cubit_attrib_type(mTypes.get(), attrib_list);
      if (attrib_list.size() == 0) {
        create_cubit_attrib(mTypes.get(), cau, NULL);
        if (status == CUBIT_FAILURE) break;
      }
    }
    mAutoUpdateFlags.step();
    mTypes.step();
  }

  return status;
}

void CubitAttribManager::set_all_auto_actuate_flags(CubitBoolean value)
{

//  std::map<int, CARegisterData*>::iterator iter = mTypeToCAData.begin();

//  while (mTypeToCAData.end() != iter)
//  {

  mAutoActuateFlags.reset();
  for (int i = mAutoActuateFlags.size(); i>0; i--)
  {
    mAutoActuateFlags.change_to(value);
    mAutoActuateFlags.step();
  }
}

CubitStatus CubitAttribManager::set_auto_actuate_flag(int attrib_type, CubitBoolean value)
{
  assert(CA_UNDEFINED != attrib_type);
  assert(CA_ALL_ATTRIBUTES != attrib_type);

  int index = mTypes.where_is_item(attrib_type);

  if (-1 == index)
    return CUBIT_FAILURE;

  mAutoActuateFlags.reset();
  mAutoActuateFlags.step(index);
  mAutoActuateFlags.change_to(value);
  return CUBIT_SUCCESS;
}

void CubitAttribManager::set_all_auto_update_flags(CubitBoolean value)
{
  mAutoUpdateFlags.reset();
  for (int i = mAutoUpdateFlags.size(); i>0; i--)
  {
    mAutoUpdateFlags.change_to(value);
    mAutoUpdateFlags.step();
  }
}

CubitStatus CubitAttribManager::set_auto_update_flag(int attrib_type, CubitBoolean value)
{
  assert(CA_UNDEFINED != attrib_type);
  assert(CA_ALL_ATTRIBUTES != attrib_type);

  int index = mTypes.where_is_item(attrib_type);

  if (-1 == index)
    return CUBIT_FAILURE;

  mAutoUpdateFlags.reset();
  mAutoUpdateFlags.step(index);
  mAutoUpdateFlags.change_to(value);
  return CUBIT_SUCCESS;
}

void CubitAttribManager::set_all_auto_write_flags(CubitBoolean value)
{
  mAutoWriteFlags.reset();
  for (int i = mAutoWriteFlags.size(); i>0; i--)
  {
    mAutoWriteFlags.change_to(value);
    mAutoWriteFlags.step();
  }
}

CubitStatus CubitAttribManager::set_auto_write_flag(int attrib_type, CubitBoolean value)
{
  assert(CA_UNDEFINED != attrib_type);
  assert(CA_ALL_ATTRIBUTES != attrib_type);

  int index = mTypes.where_is_item(attrib_type);

  if (-1 == index)
    return CUBIT_FAILURE;

  mAutoWriteFlags.reset();
  mAutoWriteFlags.step(index);
  mAutoWriteFlags.change_to(value);
  return CUBIT_SUCCESS;
}

void CubitAttribManager::set_all_auto_read_flags(CubitBoolean value)
{
  mAutoReadFlags.reset();
  for (int i = mAutoReadFlags.size(); i>0; i--)
  {
    mAutoReadFlags.change_to(value);
    mAutoReadFlags.step();
  }
}

CubitStatus CubitAttribManager::set_auto_read_flag(int attrib_type, CubitBoolean value)
{
  assert(CA_UNDEFINED != attrib_type);
  assert(CA_ALL_ATTRIBUTES != attrib_type);

  int index = mTypes.where_is_item(attrib_type);

  if (-1 == index)
    return CUBIT_FAILURE;

  mAutoReadFlags.reset();
  mAutoReadFlags.step(index);
  mAutoReadFlags.change_to(value);
  return CUBIT_SUCCESS;
}

void CubitAttribManager::auto_flag(int flag) 
{
  if (flag == -1) {
    flag = auto_flag();
    if (flag == -1) {
      if (!silentFlag) PRINT_ERROR("Can't change attribute flag with toggle, "
                  "some are already set.\n");
      return;
    }
  }

  CubitBoolean set_flag = (flag == 1 ? CUBIT_TRUE : CUBIT_FALSE);
  
  mAutoUpdateFlags.reset();
  mAutoActuateFlags.reset();
  assert (mAutoUpdateFlags.size() == mAutoActuateFlags.size());
  for (int i = mAutoUpdateFlags.size(); i>0; i--)
  {
    mAutoUpdateFlags.change_to(set_flag);
    mAutoUpdateFlags.step();
    mAutoActuateFlags.change_to(set_flag);
    mAutoActuateFlags.step();
  }

  if (!set_flag) {
      // make sure entity_name flag isn't set false here

    int index = mTypes.where_is_item(CA_ENTITY_NAME);
    assert(-1 != index);

    mAutoUpdateFlags.reset();
    mAutoUpdateFlags.step(index);
    mAutoUpdateFlags.change_to(CUBIT_TRUE);

    mAutoActuateFlags.reset();
    mAutoActuateFlags.step(index);
    mAutoActuateFlags.change_to(CUBIT_TRUE);
  }
}
  
int CubitAttribManager::auto_flag() 
{
  int sum = 0;
  CubitBoolean b_temp;
  mAutoUpdateFlags.reset();
  mAutoActuateFlags.reset();
  assert (mAutoUpdateFlags.size() == mAutoActuateFlags.size());
  for (int i = mAutoUpdateFlags.size(); i>0; i--)
  {
    b_temp = mAutoUpdateFlags.get_and_step();
    if (b_temp) sum++;
    b_temp = mAutoActuateFlags.get_and_step();
    if (b_temp) sum++;
  }
    
  if (2*mAutoUpdateFlags.size() == sum) return 0;
  else if (sum == 0) return 1;
  else return -1;
}

int CubitAttribManager::attrib_type(const char* name)
{
  int i;
  mTypeNames.reset();
  assert (mTypes.size() == mTypeNames.size());
  for (i=0; i<mTypeNames.size(); i++)
  {
    if (!strcmp(name, mTypeNames.get_and_step()))
    {
      mTypes.reset();
      mTypes.step(i);
      return mTypes.get();
    }
  }

  return CA_UNDEFINED;
}

int CubitAttribManager::attrib_type(CubitSimpleAttrib *csa_ptr)
{
  CubitString char_type = csa_ptr->character_type();
  return attrib_type_from_internal_name(char_type.c_str());
}

int CubitAttribManager::attrib_type_from_internal_name(const char* name)
{
  int i;
  mInternalNames.reset();
  assert (mTypes.size() == mInternalNames.size());
  for (i=0; i<mInternalNames.size(); i++)
  {
    if (!strcmp(name, mInternalNames.get_and_step()))
    {
      mTypes.reset();
      mTypes.step(i);
      return mTypes.get();
    }
  }

  return CA_UNDEFINED;
}

CubitBoolean CubitAttribManager::auto_actuate_flag(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mAutoActuateFlags.reset();
  mAutoActuateFlags.step(index);
  return mAutoActuateFlags.get();
}

CubitBoolean CubitAttribManager::auto_update_flag(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mAutoUpdateFlags.reset();
  mAutoUpdateFlags.step(index);
  return mAutoUpdateFlags.get();
}

CubitBoolean CubitAttribManager::auto_write_flag(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mAutoWriteFlags.reset();
  mAutoWriteFlags.step(index);
  return mAutoWriteFlags.get();
}

CubitBoolean CubitAttribManager::auto_read_flag(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mAutoReadFlags.reset();
  mAutoReadFlags.step(index);
  return mAutoReadFlags.get();
}

//- return the internal name of this CA given the enumerated attribute type
const char * CubitAttribManager::att_internal_name(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return NULL;
  }

  mInternalNames.reset();
  mInternalNames.step(index);
  return mInternalNames.get();
} 

//- return the name of this CA given the enumerated attribute type
const char * CubitAttribManager::att_name(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return NULL;
  }

  mTypeNames.reset();
  mTypeNames.step(index);
  return mTypeNames.get();
} 

CubitBoolean CubitAttribManager::actuate_in_constructor(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mActuateInConstructor.reset();
  mActuateInConstructor.step(index);
  return mActuateInConstructor.get();
}

CubitBoolean CubitAttribManager::actuate_after_geom_changes(int attrib_type)
{
  int index = mTypes.where_is_item(attrib_type);
  if (-1 == index) {
    if (!silentFlag) PRINT_ERROR("Attribute type %d not recognized.\n", attrib_type);
    return false;
  }

  mActuateAfterGeomChanges.reset();
  mActuateAfterGeomChanges.step(index);
  return mActuateAfterGeomChanges.get();
}

void CubitAttribManager::get_registered_types(DLIList<int> &types)
{
  types = mTypes;
  types.reset();
}


