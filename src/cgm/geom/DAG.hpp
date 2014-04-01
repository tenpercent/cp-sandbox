//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : DAG.hpp
//
// Purpose       : The DAG class performs general operations related
//                 to DAGs and ModelEntitys. Its current functionalities:
//                   1)Maintain a list of deactivated ModelEntitys
//		 Note: The "nodes" of the DAG are type ModelEntity.
//
// Special Notes : The class is a singleton.
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

#ifndef DAG_HPP
#define DAG_HPP

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********

// ********** BEGIN OPEN INVENTOR INCLUDES    **********
// ********** END OPEN INVENTOR INCLUDES      **********

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "CubitDefines.h"
#include "CubitEventDefines.h"
#include "CubitGeomConfigure.h"

#include "DLIList.hpp"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS      **********
class ModelEntity;
// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********
// ********** END ENUM DEFINITIONS            **********

class CUBIT_GEOM_EXPORT DAG
{

// ********** BEGIN FRIEND CLASS DECLARATIONS **********
// ********** END FRIEND CLASS DECLARATIONS   **********

public:

  static DAG* instance(); 
    //- Controlled access and creation of the sole instance of this class.
    //- Singleton pattern

  ~DAG() ;
    //- Destructor

  void add_deactivated_DAG_node(ModelEntity* deactivatedDAGNodePtr) ;
    //R void
    //I deactivatedDAGNodePtr
    //I- A pointer to a deactivated DAG node (ModelEntity)
    //- This function takes a pointer to a deactivated DAG node and
    //- adds it to the list of deactivated DAG nodes.

  void remove_deactivated_DAG_node(ModelEntity* deactivatedDAGNodePtr) ;
    //R void
    //I deactivatedDAGNodePtr
    //I- A pointer to a deactivated DAG node
    //- This function takes a pointer to a deactivated DAG node and
    //- removes it from the list of deactivated DAG nodes.

  void cleanout_deactivated_DAG_nodes() ;
    //R void
    //- This function deletes all the deactivated DAG nodes and cleans
    //- the list of deactivated DAG nodes. 
    //- N.B. deletion of objects connected to the DAGNodes are not
    //-      done at base class. The classes derived from DAGNode
    //-      must do that.
      
#ifdef BOYD14
  void remove_and_cleanout_deactivated_DAG_nodes() ;
    //R void
    //- This function removes the deactivated DAG nodes from the DAG
    //- then calls the function cleanout_deactivated_DAG_nodes.
#endif
      
  void notify(ModelEntity* DAGNodePtr, CubitEventType event);
    //R void
    //I DAGNodePtr
    //I- The input DAG Node pointer to the DAG Node that initiated
    //I- the event
    //I event
    //I- The type of event that occurred -- enum in CubitDefines
    //- This function is used to notify the DAG class that an
    //- event has occurred.  The input object that initiated the
    //- the event is a DAG Node.

// Methods added by Jason Kraftcheck on September 16, 1996.

  CubitStatus get_children_at_level( ModelEntity* parent, int level, 
                                     DLIList<ModelEntity*>& result_set );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I parent
    //I- A pointer to the DAG node to get the children of.
    //I level
    //I- Levels down to traverse DAG.
    //O result_set
    //O- The list of DAGNodes found.
    //- These methods return the children of a DAG node at the
    //- specified level down from the passed node.
    //-
    //- Note: These methods are implemented here rather than in
    //-       DAGNode because a breadth-first search will be
    //-	much more efficient than a depth-first search in
    //-	DAGNode.

  CubitStatus get_parents_at_level( ModelEntity* child, int level, 
                                    DLIList<ModelEntity*>& result_set );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I child
    //I- A pointer to the DAG node to get the parents of.
    //I level
    //I- Levels up to traverse DAG.
    //O result_set
    //O- The list of DAGNodes found.
    //- These methods return the parents of a DAG node at the
    //- specified level up from the passed node.
    //-
    //- Note: These methods are implemented here rather than in
    //-       DAGNode because a breadth-first search will be
    //-	much more efficient than a depth-first search in
    //-	DAGNode.
      
  CubitStatus get_children_at_level( DLIList<ModelEntity*>& parents, int level, 
                                     DLIList<ModelEntity*>& result_set );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I parents
    //I- A list of DAG nodes to get the children of.
    //I level
    //I- Levels down to traverse DAG.
    //O result_set
    //O- The list of DAGNodes found.
    //- These methods return the children of a DAG node at the
    //- specified level down from the passed node.
    //-
    //- Note: It is assumed that the parents passed are all at the
    //-	same level in the DAG.  The behavior of this method
    //-	is undefined if this it not the case.

  CubitStatus get_parents_at_level( DLIList<ModelEntity*>& children, int level, 
                                    DLIList<ModelEntity*>& result_set );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I child
    //I- A list of DAG nodes to get the parents of.
    //I level
    //I- Levels up to traverse DAG.
    //O result_set
    //O- The list of DAGNodes found.
    //- These methods return the parents of a DAG node at the
    //- specified level up from the passed node.
    //-
    //- Note: It is assumed that the children passed are all at the
    //-	same level in the DAG.  The behavior of this method
    //-	is undefined if this it not the case.
      
// End methods added by Jason Kraftcheck on September 16, 1996.

protected: 

  DAG() ;
    //- The default constructor

private:

  static DAG* instance_ ;
    //- Pointer to the only instance of the class

  DLIList<ModelEntity*> deactivatedDAGNodeList_ ;
    //- List of deactivated snodes
} ;


// ********** BEGIN INLINE FUNCTIONS          **********
// ********** END INLINE FUNCTIONS            **********

// ********** BEGIN FRIEND FUNCTIONS          **********
// ********** END FRIEND FUNCTIONS            **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END HELPER CLASS DECLARATIONS   **********

#endif

