//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : DAG.cc
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********

// ********** BEGIN OPEN INVENTOR INCLUDES    **********
// ********** END OPEN INVENTOR INCLUDES      **********

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "DAG.hpp"
#include "CastTo.hpp"
#include "ModelEntity.hpp"
#include "CubitObserver.hpp"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN EXTERN DECLARATIONS       **********
// ********** END EXTERN DECLARATIONS         **********

// ********** BEGIN STATIC DECLARATIONS       **********

DAG* DAG::instance_ = NULL ;

// ********** END STATIC DECLARATIONS         **********

// ********** BEGIN PUBLIC FUNCTIONS          **********

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the only instance of the class.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------
DAG* DAG::instance()
{
  if (instance_ == NULL ) 
  {
    instance_ = new DAG;
  }
  return instance_ ;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------
DAG::~DAG()
{
  instance_ = NULL ;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a pointer to a deactivated DAG 
//                 node and adds it to the list of deactivated DAG nodes.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------

void DAG::add_deactivated_DAG_node(ModelEntity* deactivatedDAGNodePtr) 
{
  deactivatedDAGNodeList_.append(deactivatedDAGNodePtr) ;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a pointer to a deactivated DAG 
//                 node and removes it to the list of deactivated DAG nodes.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------

void DAG::remove_deactivated_DAG_node(ModelEntity* deactivatedDAGNodePtr) 
{
  if (deactivatedDAGNodeList_.move_to(deactivatedDAGNodePtr))
    deactivatedDAGNodeList_.extract();
}

//-------------------------------------------------------------------------
// Purpose       : This function deletes all the deactivated DAG nodes
//                 and cleans the list of deactivated DAG nodes.
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------
void DAG::cleanout_deactivated_DAG_nodes() 
{
  ModelEntity* tempDAGNode = NULL;
  while ( deactivatedDAGNodeList_.size() > 0)
  {
    tempDAGNode = NULL;
    deactivatedDAGNodeList_.last();
    tempDAGNode = deactivatedDAGNodeList_.get();
    if( tempDAGNode != NULL )
    {
      delete tempDAGNode;
    }
  }
   
  deactivatedDAGNodeList_.clean_out() ;
}

//-------------------------------------------------------------------------
// Purpose       :  This function is used to notify the DAG class that an
//                  event has occurred.  The input object that initiated the
//                  the event is a DAG Node.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/03/96
//-------------------------------------------------------------------------
void DAG::notify(ModelEntity* DAGNodePtr, CubitEventType event)
{
    // The input DAG Node was destructed.
  if ( event == DAG_NODE_DESTRUCTED )
  {
      // If this input pointer exists in the deactivated DAG Node list,
      // remove it
    deactivatedDAGNodeList_.remove(DAGNodePtr);
  }
}


//-------------------------------------------------------------------------
// Purpose	 : Find the children of a node at a specified level in the 
//		   DAG, relative to the node.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 09/16/97
//-------------------------------------------------------------------------
CubitStatus DAG::get_children_at_level( ModelEntity* parent, int level,
                                        DLIList<ModelEntity*>& result_set )
{
  DLIList<ModelEntity*> parent_list;  // nodes at current level
	
    // Start with a parent_list with one node, the passed parent
  parent_list.append( parent );
	
    // Call with list
  return get_children_at_level( parent_list, level, result_set );
}

//-------------------------------------------------------------------------
// Purpose	 : Find the parents of a node at a specified level in the 
//		   DAG, relative to the node.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 09/16/97
//-------------------------------------------------------------------------
CubitStatus DAG::get_parents_at_level( ModelEntity* child, int level,
                                       DLIList<ModelEntity*>& result_set )
{
  DLIList<ModelEntity*> child_list;  // nodes at current level
	
    // Start with a parent_list with one node, the passed parent
  child_list.append( child );
	
    // Call with list
  return get_parents_at_level( child_list, level, result_set );
}


//-------------------------------------------------------------------------
// Purpose	 : Find the children of a list of nodes at a specified 
//		   level in the DAG, relative to the node.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 09/16/97
//-------------------------------------------------------------------------
CubitStatus DAG::get_children_at_level( DLIList<ModelEntity*>& parent_list, 
                                        int level,
                                        DLIList<ModelEntity*>& result_set )
{
  assert( level > 0 );
  if( parent_list.size( ) < 1 ) return CUBIT_FAILURE;
	
  DLIList<ModelEntity*> child_list;   // children of all parents
  DLIList<ModelEntity*> current_node_children;
	
    //Get a list of all the children one level down from 
    //the nodes in parent_list.
  for( int i = 0; i < parent_list.size( ); i++ )
  {
    ModelEntity* current = parent_list.get_and_step( );
		
    current_node_children.clean_out( );
    current->get_children( &current_node_children );
		
    for( int j = 0; j < current_node_children.size( ); j++ )
    {
      child_list.append_unique(
        current_node_children.get_and_step( ) );
    }
  }
	
    //If we found no children, return failure.
  if( child_list.size( ) < 1 ) return CUBIT_FAILURE;
	
    //If level 1, then return the list of children.
  if( level == 1 )
  {
    result_set = child_list;
    return CUBIT_SUCCESS;
  }
    //Otherwise get next level down by recursion.
  else
  {
    return get_children_at_level( child_list, level - 1, 
                                  result_set );
  }	
}

//-------------------------------------------------------------------------
// Purpose	 : Find the parents of a list of nodes at a specified level
//		   in the DAG, relative to the node.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 09/16/97
//-------------------------------------------------------------------------
CubitStatus DAG::get_parents_at_level( DLIList<ModelEntity*>& child_list, 
                                       int level,
                                       DLIList<ModelEntity*>& result_set )
{
  assert( level > 0 );
  if( child_list.size( ) < 1 ) return CUBIT_FAILURE;
	
  DLIList<ModelEntity*> parent_list;   // children of all parents
  DLIList<ModelEntity*> current_node_parents;
	
    //Get a list of all the parents one level up from 
    //the nodes in child_list.
  for( int i = 0; i < child_list.size( ); i++ )
  {
    ModelEntity* current = child_list.get_and_step( );
		
    current_node_parents.clean_out( );
    current->get_parents( &current_node_parents );
		
    for( int j = 0; j < current_node_parents.size( ); j++ )
    {
      parent_list.append_unique(
        current_node_parents.get_and_step( ) );
    }
  }
	
    //If we found no parents, return failure.
  if( parent_list.size( ) < 1 ) return CUBIT_FAILURE;
	
    //If level 1, then return the list of parents.
  if( level == 1 )
  {
    result_set = parent_list;
    return CUBIT_SUCCESS;
  }
    //Otherwise get next level down by recursion.
  else
  {
    return get_parents_at_level( parent_list, level - 1, 
                                 result_set );
  }	
}




// ********** END PUBLIC FUNCTIONS            **********

// ********** BEGIN PROTECTED FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The default constructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------
DAG::DAG() : deactivatedDAGNodeList_() 
{
}

// ********** END PROTECTED FUNCTIONS         **********

// ********** BEGIN PRIVATE FUNCTIONS         **********
// ********** END PRIVATE FUNCTIONS           **********

// ********** BEGIN HELPER CLASSES            **********
// ********** END HELPER CLASSES              **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN STATIC FUNCTIONS          **********
// ********** END STATIC FUNCTIONS            **********

