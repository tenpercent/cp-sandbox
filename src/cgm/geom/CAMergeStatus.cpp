//-------------------------------------------------------------------------
// Filename      : CAMergeStatus.cpp
//
// Purpose       : Save RefEntity::autoMergeStatus
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/06/02
//-------------------------------------------------------------------------

#include "RefEntity.hpp"
#include "CAMergeStatus.hpp"
#include "CADefines.hpp"

CubitAttrib* CAMergeStatus_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAMergeStatus *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAMergeStatus(entity);
  }
  else
  {
    new_attrib = new CAMergeStatus(entity, p_csa);
  }

  return new_attrib;
}

CAMergeStatus::CAMergeStatus( RefEntity* owner, CubitSimpleAttrib* csa )
  : CubitAttrib( owner )
{ 
  assert( csa && csa->int_data_list()->size() == 1 );
  int i = *(csa->int_data_list()->get());
  assert( i == 0 || i == 1 || i == 2 );
  status = (AutoMergeStatus)i;
}

CAMergeStatus::CAMergeStatus( RefEntity* owner )
  : CubitAttrib( owner ), status(AUTO_MERGE_AUTO)
{ }

CAMergeStatus::~CAMergeStatus() 
{ }

CubitStatus CAMergeStatus::actuate() 
{
  if( hasActuated )
    return CUBIT_SUCCESS;
  
  if( !attribOwnerEntity )
    return CUBIT_FAILURE;
    
  DLIList<CubitAttrib*> att_list;
  attribOwnerEntity->find_cubit_attrib_type(CA_MERGE_PARTNER, att_list);
  if( att_list.size() )
    return CUBIT_FAILURE;
  
  assert( status == AUTO_MERGE_ON || status == AUTO_MERGE_OFF );  
  attribOwnerEntity->is_mergeable(status);
  
  deleteAttrib = CUBIT_FALSE;
  hasActuated = CUBIT_TRUE;

  return CUBIT_SUCCESS;
}

CubitStatus CAMergeStatus::update()
{
  if( hasUpdated )
    return CUBIT_SUCCESS;
    
  if( !attribOwnerEntity )
    return CUBIT_FAILURE;
  
  status = attribOwnerEntity->merge_status();
  if( status == AUTO_MERGE_AUTO )
    delete_attrib(CUBIT_TRUE);
  else
    hasUpdated = CUBIT_TRUE;
  
  return CUBIT_SUCCESS;
}

CubitStatus CAMergeStatus::reset()
{
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CAMergeStatus::cubit_simple_attrib()
{
  if( deleteAttrib )
    return 0;
  
  assert( status != AUTO_MERGE_AUTO );
  CubitSimpleAttrib* result = 0;

  CubitString name( att_internal_name() );
  DLIList<CubitString*> string_list(1);
  string_list.append( &name );
  
  int int_data = (int)status;
  DLIList<int> int_list;
  int_list.append( int_data );
  
  result = new CubitSimpleAttrib( &string_list, 0, &int_list );
  return result;
}

int CAMergeStatus::int_attrib_type()
{
  return CA_MERGE_STATUS;
}

const type_info& CAMergeStatus::entity_type_info() const
{
  return typeid(CAMergeStatus);
}

