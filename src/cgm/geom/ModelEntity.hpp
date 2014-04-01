//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : ModelEntity.hpp
//
// Purpose       : This class serves as the pure abstract class of
//                 all the Model-related classes. Each ModelEntity Object
//		   can be considered a Node in the DAG.
//
// Special Notes : This is a pure virtual class to prevent instantiation.
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/15/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef MODEL_ENTITY_HPP
#define MODEL_ENTITY_HPP

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
#include "DagType.hpp"
#include "CubitGeomConfigure.h"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS      **********

template <class X> class DLIList;

// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********
// ********** END ENUM DEFINITIONS            **********

class CUBIT_GEOM_EXPORT ModelEntity
{

// ********** BEGIN FRIEND CLASS DECLARATIONS **********

   friend class ModelQueryEngine ;
   friend class GroupingEntity;
// ********** END FRIEND CLASS DECLARATIONS   **********

 public:

   ModelEntity();
      //- Default constructor

   virtual ~ModelEntity();

   virtual int get_children( DLIList<ModelEntity*>* list = 0 ) const = 0;
   virtual int get_parents ( DLIList<ModelEntity*>* list = 0 ) const = 0;
      //- Returns number of parents/children.  
      //- If a list is passed, pointers to parents/children will
      //- be appended to the passed list.
      
   virtual DagType dag_type() const = 0;
      
   void deactivated(CubitBoolean flag) ;
   CubitBoolean deactivated() const ;
      //- Get and set funtions for the deactivation flag.
	
   CubitStatus remove_from_DAG(CubitBoolean recursive_flag = CUBIT_FALSE);
   void disconnect_from_DAG();
   virtual CubitStatus remove_child_link( ModelEntity* child_ptr ) = 0;

   virtual CubitStatus 
    disconnect_all_children( DLIList<ModelEntity*>* children = 0 ) = 0;
   virtual CubitStatus
    disconnect_all_parents( DLIList<ModelEntity*>* parents = 0 ) = 0;

protected:

   virtual CubitBoolean query_append_parents( DLIList<ModelEntity*>& list ) = 0;
   virtual CubitBoolean query_append_children(DLIList<ModelEntity*>& list ) = 0;
   
   cBit deactivatedStatus_ : 1;
     //- The deactivated flag
   
   cBit encountered_ : 1;
    // Indicates if it has already been found in 
    // in a DAG query.  



private:

   ModelEntity const& operator=(ModelEntity const& aModEnt) ;
   ModelEntity(const ModelEntity& aModEnt);
    // These are declared here to prevent the compiler
    // from generating default ones.  These have no implementation.

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

