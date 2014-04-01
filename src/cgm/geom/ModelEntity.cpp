//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : ModelEntity.cc
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/15/96
//
// Owner         : Malcolm J. Panthaki
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

#include "ModelEntity.hpp"
#include "ModelQueryEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "CastTo.hpp"
#include "CubitObservable.hpp"
#include "CubitObserver.hpp"
#include "DAG.hpp"
#include "RefVolume.hpp"
#include "Loop.hpp"
#include "RefEdge.hpp"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN EXTERN DECLARATIONS       **********
// ********** END EXTERN DECLARATIONS         **********

// ********** BEGIN STATIC DECLARATIONS       **********
// ********** END STATIC DECLARATIONS         **********

// ********** BEGIN PUBLIC FUNCTIONS          **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor
//
// Special Notes : 
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/15/96
//-------------------------------------------------------------------------

ModelEntity::ModelEntity() : 
    deactivatedStatus_(CUBIT_FALSE),
    encountered_(CUBIT_FALSE)
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/15/96
//-------------------------------------------------------------------------

ModelEntity::~ModelEntity()
{
     // Make sure that the ModelEntity class does not have a pointer to this
     // ModelEntity in its list.  If it does, then remove the instances
     // of the pointer
   DAG::instance()->notify(this, DAG_NODE_DESTRUCTED);
  
}


//-------------------------------------------------------------------------
// Purpose       : Determine whether this object can be deleted from
//                 memory or not.
//
// Special Notes : This function will recursively remove all the
//                 removable child entities from the DAG and call their
//                 destructors. Destructors of other objects, like
//                 GeometryEntity, need to be called from the destructors 
//                 of different ModEnts. Virtual destructors put to use :) :)
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/05/96
//-------------------------------------------------------------------------

CubitStatus ModelEntity::remove_from_DAG(CubitBoolean recurse_flag)
{
     // This counter will be used in the recursion to test whether
     // the call is from outside or from the function itself. When
     // the call comes from outside, the counter should always be
     // zero.

  CubitBoolean this_recurse = recurse_flag;
  if (recurse_flag == CUBIT_FALSE) recurse_flag = CUBIT_TRUE;
  
  DLIList<ModelEntity*> childModEntList;
   
     // Check to see if there are no parents of this object. 
   if ( get_parents() == 0 )
   {
      if (this_recurse == CUBIT_FALSE)
      {
         // Since we are not recursing, this is a top-level entity.
         // Notify the static observers that a top-level entity is being
         // destructed.  This must be done before children are disconnected.
       CubitObservable *top_level = CAST_TO(this, CubitObservable);
       CubitObserver::notify_static_observers(top_level, TOP_LEVEL_ENTITY_DESTRUCTED);
      }

        // Go through all the children and remove their link to
        // the current object.
      
      ModelEntity* tempModEntPtr = NULL ;
      ModelEntity* childModEntPtr = NULL ;
      
      childModEntList.clean_out();
      disconnect_all_children(&childModEntList);
      
        // The following while conditional may not work...it depends on
        // what is_at_end does when you step over the end...CHECK THIS
      int i;
      for( i = 0 ; i < childModEntList.size() ; i++ )
      {
           // Get the next ModelEnti in the child list and make sure its
           // pointer to its parent is removed.
         tempModEntPtr = childModEntList.get_and_step();
          
           // Try remove() on the child ModEnt. If it comes back with
           // a success, then delete it.
         childModEntPtr = tempModEntPtr;
         
         if ( childModEntPtr->remove_from_DAG(recurse_flag) == CUBIT_SUCCESS )
         {
              // Now deactivate the child ModEnt
            childModEntPtr->deactivated(CUBIT_TRUE) ;

              // remove it from observables, just before we go to delete it
            CubitObservable *observable = CAST_TO(childModEntPtr, CubitObservable);
            if (observable) 
            {
              if( !observable->notify_observers( MODEL_ENTITY_DESTRUCTED ) )
                 return CUBIT_FAILURE;
            }
         }
      }
      
      
        // If this is the top of the recursion, then clean out all the deactivated
        // entities.
      if (this_recurse == CUBIT_FALSE)
      {
         this->deactivated(CUBIT_TRUE) ;

         // remove it from observables, just before we go to delete it
         CubitObservable *observable = CAST_TO(childModEntPtr, CubitObservable);
         if (observable) 
         {
           if( !observable->notify_observers( MODEL_ENTITY_DESTRUCTED ) )
              return CUBIT_FAILURE;

         }
         GeometryQueryTool::instance()->cleanout_deactivated_geometry() ;
      }
      
      return CUBIT_SUCCESS ;
   }
   else
   {
      return CUBIT_FAILURE ;
   }
}

void ModelEntity::disconnect_from_DAG()
{
    // disconnects this entity from any others
    // to which it is connected in the DAG; does not delete the DAGNode

  disconnect_all_children();
  disconnect_all_parents();
}

//-------------------------------------------------------------------------
// Purpose       : Get and set functions for the deactivated flag.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------

void ModelEntity::deactivated(CubitBoolean flag)
{
   if (deactivatedStatus_ != flag)
   {
     deactivatedStatus_ = flag ;
     if (flag == CUBIT_TRUE)
     {
        DAG::instance()->add_deactivated_DAG_node(this) ;
     }
     else
     {
        DAG::instance()->remove_deactivated_DAG_node(this) ;
     }
   }
}

CubitBoolean ModelEntity::deactivated() const
{
   return (CubitBoolean)deactivatedStatus_ ;
}




// ********** END PUBLIC FUNCTIONS            **********
// ********** BEGIN PROTECTED FUNCTIONS       **********
// ********** END PROTECTED FUNCTIONS         **********
// ********** BEGIN PRIVATE FUNCTIONS         **********
















// ********** END PRIVATE FUNCTIONS           **********

// ********** BEGIN HELPER CLASSES            **********
// ********** END HELPER CLASSES              **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN STATIC FUNCTIONS          **********
// ********** END STATIC FUNCTIONS            **********

