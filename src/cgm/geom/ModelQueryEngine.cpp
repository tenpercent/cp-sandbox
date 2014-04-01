//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//
// Filename      : ModelQueryEngine.C 
//
// Purpose       : This file contains the implementation of the class 
//                 ModelQueryEngine.
//
// Special Notes : 
//
// Creator       : Xuechen Liu 
//
// Creation Date : 06/08/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********

#include <stdio.h>

// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitMessage.hpp"
#include "ModelQueryEngine.hpp"
#include "ModelEntity.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "RefVertex.hpp"
#include "CoVertex.hpp"
#include "Chain.hpp"
#include "RefEdge.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "RefFace.hpp"
#include "CoFace.hpp"
#include "Shell.hpp"
#include "Body.hpp"
#include "CoVolume.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********

ModelQueryEngine* ModelQueryEngine::instance_ = NULL;

// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Controls the access and creation of the sole instance
//                 of this class.
//
// Special Notes :
//
// Creator       : Xuechen Liu 
//
// Creation Date : 06/08/96
//-------------------------------------------------------------------------
ModelQueryEngine* ModelQueryEngine::instance()
{
  if (instance_ == NULL)
  {
    instance_ = new ModelQueryEngine;
  }
  return instance_;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor of the ModelQueryEngine class. 
//
// Special Notes :
//
// Creator       : Xuechen Liu 
//
// Creation Date : 06/08/96
//-------------------------------------------------------------------------
ModelQueryEngine::~ModelQueryEngine()
{
   // Set static instance_ to zero to indicated that we are dead.
   instance_ = NULL;
}


// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : Constructor of the ModelQueryEngine class. 
//
// Special Notes :
//
// Creator       : Xuechen Liu 
//
// Creation Date : 06/08/96
//-------------------------------------------------------------------------
ModelQueryEngine::ModelQueryEngine() : queryCallStackDepth(0)
{
}


//-------------------------------------------------------------------------
// Purpose       : Query model
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_model_and_append( 
                                          ModelEntity& source_object,
                                          DagType target_type,
                                          DLIList<ModelEntity*>& result_set )
{
  BeginQuery lock;
  
  DagType source_type = source_object.dag_type();

  if (!source_type.is_valid() || !target_type.is_valid())
    return CUBIT_FAILURE;

  else if (source_type < target_type)
    return query_append_parents( source_object, target_type, result_set );

  else if (source_type > target_type)
    return query_append_children( source_object, target_type, result_set );  

  else // same type
  {
    assert(source_type == target_type);
    if (!encountered(&source_object))
      result_set.append(&source_object);
    return CUBIT_SUCCESS;
  }
}


//-------------------------------------------------------------------------
// Purpose       : Query model
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_model( 
                                          ModelEntity& source_object,
                                          DagType target_type,
                                          DLIList<ModelEntity*>& result_set )
{
  result_set.clean_out();
  return query_model_and_append( source_object, target_type, result_set );
}


//-------------------------------------------------------------------------
// Purpose       : Query model
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_model( 
                                    DLIList<ModelEntity*>& source_set,
                                    DagType target_type,
                                    DLIList<ModelEntity*>& result_set )
{
  result_set.clean_out();
  return query_model_and_append( source_set, target_type, result_set );
}


//-------------------------------------------------------------------------
// Purpose       : Query model
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_model_and_append( 
                                    DLIList<ModelEntity*>& source_set,
                                    DagType target_type,
                                    DLIList<ModelEntity*>& result_set )
{
  BeginQuery lock;
  CubitStatus result = CUBIT_SUCCESS;
  
  source_set.reset();
  for (int i = source_set.size(); i--; )
  {
    ModelEntity& source = *source_set.get_and_step();
    if (!query_model_and_append( source, target_type, result_set ))
      result = CUBIT_FAILURE;
  }
  
  return result;
}

// ********** END PROTECTED FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Query downwards
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_append_children ( 
                                        ModelEntity& source_object,
                                        DagType target_type,
                                        DLIList<ModelEntity*>& result_set )
{
  BeginQuery lock;
  
  DagType current_type = source_object.dag_type();
  assert(current_type.is_valid() && target_type.is_valid());
  assert(current_type > target_type);
  
  intermediateNodeSets[0].clean_out();
  intermediateNodeSets[0].append(&source_object);
  int current_index = 0;
  
  while (current_type > target_type)
  {
    DLIList<ModelEntity*>& current_set = intermediateNodeSets[ current_index];
    DLIList<ModelEntity*>&    next_set = intermediateNodeSets[!current_index];
    
    next_set.clean_out();
    current_set.reset();
    for (int i = current_set.size(); i--; )
    {
      ModelEntity* current_ptr = current_set.get_and_step();
      current_ptr->query_append_children(next_set);
    }

    current_index = !current_index;
    current_type--;
  }
  
  result_set += intermediateNodeSets[current_index];
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Query upwards
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus ModelQueryEngine::query_append_parents( 
                                        ModelEntity& source_object,
                                        DagType target_type,
                                        DLIList<ModelEntity*>& result_set )
{
  BeginQuery lock;
  
  DagType current_type = source_object.dag_type();
  assert(current_type.is_valid() && target_type.is_valid());
  assert(current_type < target_type);
  
  intermediateNodeSets[0].clean_out();
  intermediateNodeSets[0].append(&source_object);
  int current_index = 0;
  
  while (current_type < target_type)
  {
    DLIList<ModelEntity*>& current_set = intermediateNodeSets[ current_index];
    DLIList<ModelEntity*>&    next_set = intermediateNodeSets[!current_index];
    
    next_set.clean_out();
    current_set.reset();
    for (int i = current_set.size(); i--; )
    {
      ModelEntity* current_ptr = current_set.get_and_step();
      current_ptr->query_append_parents(next_set);
    }

    current_index = !current_index;
    current_type++;
  }
  
  result_set += intermediateNodeSets[current_index];
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Mark node encountered
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/24/03
//-------------------------------------------------------------------------
bool ModelQueryEngine::encountered( ModelEntity* node_ptr )
{
  if (node_ptr->encountered_)
    return true;
  
  encounteredSet.append(node_ptr);
  node_ptr->encountered_ = true;
  return false;
}

void ModelQueryEngine::inc_query_call_stack()
{
  queryCallStackDepth++;
}

void ModelQueryEngine::dec_query_call_stack()
{
  assert(queryCallStackDepth > 0);
  queryCallStackDepth--;
  if (queryCallStackDepth == 0)
    while (encounteredSet.size())
      encounteredSet.pop()->encountered_ = CUBIT_FALSE;
}

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********
 
// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
 

