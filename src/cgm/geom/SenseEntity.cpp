//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : SenseEntity.C 
//
// Purpose       : This file contains the implementation of the class 
//                 SenseEntity. 
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

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"
#include "SenseEntity.hpp"
#include "GroupingEntity.hpp"
#include "BasicTopologyEntity.hpp"
#include "DLIList.hpp"
#include "ModelQueryEngine.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
SenseEntity::~SenseEntity()
{
  if (myParent)
    myParent->remove_sense_entity(this);
  if (myChild)
    myChild->remove_sense_entity(this);
  assert (!myParent && !myChild);
}

//-------------------------------------------------------------------------
// Purpose       : This function is used to attach a BasicTopologyEntity as a
//                 child of the current SenseEntity in the DAG. 
//
// Special Notes : Complete reimplementation - jk, July 2003
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/28/96
//-------------------------------------------------------------------------
CubitStatus SenseEntity::attach_basic_topology_entity(
    BasicTopologyEntity* basic_topology_entity_ptr)
{
  return basic_topology_entity_ptr->add_sense_entity( this );
}



//-------------------------------------------------------------------------
// Purpose       : Change the BTE attached to this SenseEntity
//
// Special Notes : derived from attach_basic_topology_entity()
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/01
//-------------------------------------------------------------------------
CubitStatus SenseEntity::switch_basic_topology_entity( BasicTopologyEntity* new_bte )
{
  if (new_bte->dag_type().parent() != dag_type())
    return CUBIT_FAILURE;
  
  if (!myChild->remove_sense_entity(this))
    return CUBIT_FAILURE;
    
  if (!new_bte->add_sense_entity(this))
    return CUBIT_FAILURE; 
    
  return CUBIT_SUCCESS;
}



 
//-------------------------------------------------------------------------
// Purpose       : Reverse
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/27/01
//-------------------------------------------------------------------------
void SenseEntity::reverse_sense()
{ 
  switch (mySense)
  {
    case CUBIT_FORWARD : mySense = CUBIT_REVERSED; break;
    case CUBIT_REVERSED: mySense = CUBIT_FORWARD ; break;
    default            : mySense = CUBIT_UNKNOWN ; break;
  }
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : get parent grouping entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int SenseEntity::get_parents( DLIList<ModelEntity*>* list ) const
{
  if (!myParent)
    return 0;
  
  if (list)
    list->append(myParent);
  
  return 1;
}

//-------------------------------------------------------------------------
// Purpose       : get child basic topology entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
int SenseEntity::get_children( DLIList<ModelEntity*>* list ) const
{
  if (!myChild)
    return 0;
  
  if (list)
    list->append(myChild);
  
  return 1;
}

//-------------------------------------------------------------------------
// Purpose       : Remove from child basic topology entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus SenseEntity::remove_child_link( ModelEntity* entity_ptr )
{
  if (entity_ptr != myChild)
    return CUBIT_FAILURE;
  
  BasicTopologyEntity* bte_ptr = static_cast<BasicTopologyEntity*>(entity_ptr);
  return bte_ptr->remove_sense_entity(this);
}

//-------------------------------------------------------------------------
// Purpose       : Get parent BasicTopologyEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
BasicTopologyEntity* SenseEntity::get_parent_basic_topology_entity_ptr()
{
  return myParent ? myParent->get_basic_topology_entity_ptr() : 0;
}

//-------------------------------------------------------------------------
// Purpose       : Remove from parent GroupingEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus SenseEntity::disconnect_all_parents( DLIList<ModelEntity*>* list )
{
  if (!myParent)
    return CUBIT_SUCCESS;
  
  if (list)
    list->append(myParent);
  
  return myParent->remove_sense_entity(this);
}

//-------------------------------------------------------------------------
// Purpose       : Disconnect form child BasicTopologyEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CubitStatus SenseEntity::disconnect_all_children( DLIList<ModelEntity*>* list )
{
  if (!myChild)
    return CUBIT_SUCCESS;
  
  if (list)
    list->append(myChild);
  
  return myChild->remove_sense_entity(this);
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
CubitBoolean SenseEntity::query_append_parents( DLIList<ModelEntity*>& list )
{
  if (myParent && !ModelQueryEngine::instance()->encountered(myParent))
  {
    list.append(myParent);
    return CUBIT_TRUE;
  }
  
  return CUBIT_FALSE;
}
CubitBoolean SenseEntity::query_append_children( DLIList<ModelEntity*>& list )
{
  if (myChild && !ModelQueryEngine::instance()->encountered(myChild))
  {
    list.append(myChild);
    return CUBIT_TRUE;
  }
  
  return CUBIT_FALSE;
}


// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

