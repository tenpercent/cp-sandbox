//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : CoFace.cc
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Jihong Ma
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

#include "CoFace.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "GroupingEntity.hpp"
#include "CastTo.hpp"
#include "GeometryDefines.h"
#include "Shell.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********


// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
CoFace::CoFace() 
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
CoFace::~CoFace() 
{
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to a face and the
//                 sense of this CoFace. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
CoFace::CoFace(RefFace* facePtr, CubitSense sense) 
{
   attach_basic_topology_entity(facePtr) ;
   set_sense(sense) ;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns a pointer to the RefFace which
//                 the current CoFace is associated with.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
RefFace* CoFace::get_ref_face_ptr()  
{
   // Call the generic function defined in the SenseEntity class to
   // do the real work
   BasicTopologyEntity* BTEPtr = get_basic_topology_entity_ptr();

   // Cast the returned pointer to RefFace and return it
   return CAST_TO( BTEPtr, RefFace );   
}

//-------------------------------------------------------------------------
// Purpose       : Get the parent Shell.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/09/99
//-------------------------------------------------------------------------
Shell* CoFace::get_shell_ptr()  
{
   GroupingEntity* gpe_ptr = get_grouping_entity_ptr();

   // Cast the returned pointer to RefFace and return it
   return CAST_TO( gpe_ptr, Shell );   
}

//-------------------------------------------------------------------------
// Purpose       : This function returns a pointer to the RefVolume which
//                 the current CoFace is associated with.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 03/14/97
//-------------------------------------------------------------------------
RefVolume* CoFace::get_ref_volume()
{
   BasicTopologyEntity* bte_ptr = get_parent_basic_topology_entity_ptr();

   // Cast the returned pointer to RefFace and return it
   return dynamic_cast<RefVolume*>(bte_ptr);
}


//-------------------------------------------------------------------------
// Purpose       : This function is called after a child of a ModelEntity 
//                 is switched. The sense of a CoFace may change if one 
//                 of its RefFace changes. This function takes care of 
//                 that. If the sense of the RefFaces that were switched 
//                 is same, nothing is done. If the RefFaces are of 
//                 opposite sense, the sense of this object is switched, 
//                 i.e. if it was FORWARD, it is made REVERSE, and vice
//                 versa.
//
// Special Notes : This function assumes the RefFaces are equivalent...
//
// Creator       : David White
//
// Creation Date : 4/3/97
//-------------------------------------------------------------------------

void CoFace::switch_child_notify(ModelEntity const* newChild, 
                                 ModelEntity const* oldChild ) 
{
   // Make sure the entities being switched are RefFaces. If not,
   // get out of this function.
  ModelEntity* tmp_new_child = const_cast<ModelEntity*>(newChild);
  ModelEntity* tmp_old_child = const_cast<ModelEntity*>(oldChild);
   RefFace* new_child_ref_face = CAST_TO(tmp_new_child, RefFace) ;
   RefFace* old_child_ref_face = CAST_TO(tmp_old_child, RefFace) ;

   if ( ( new_child_ref_face == NULL ) || ( old_child_ref_face == NULL ) )
   {
      return ;
   }

   CubitSense sense =
       old_child_ref_face->compare_alignment( new_child_ref_face );
   
     // If the sense of the old RefFace relative to the new RefFace is 
     // same, nothing needs to be done. However, if the relative sense
     // is reversed, switch the sense of the CoFace.
   if ( sense == CUBIT_REVERSED ) 
   {
     if ( get_sense() == CUBIT_FORWARD )
     {
       set_sense(CUBIT_REVERSED) ;
     }
     else
     {
       set_sense(CUBIT_FORWARD) ;
     }
   }
}


// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

