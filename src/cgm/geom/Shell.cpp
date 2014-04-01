//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Shell.C
//
// Purpose       : This file contains the implementation of the class
//                  Shell. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Xuechen Liu
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
#include "Shell.hpp"
#include "CoFace.hpp"
#include "CastTo.hpp"
#include "DLIList.hpp"
#include "RefFace.hpp"
#include "CubitBox.hpp"
#include "RefVolume.hpp"
#include "ShellSM.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
Shell::Shell()
{
}

//-------------------------------------------------------------------------
// Purpose       : A constructor woth a pointer to an other solid
//                 model entity.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
Shell::Shell(ShellSM* OSMEPtr)
{
   bridge_manager()->add_bridge(OSMEPtr) ;
}



RefVolume* Shell::get_ref_volume_ptr()
{
	return CAST_TO( get_basic_topology_entity_ptr(), RefVolume );
}

//-------------------------------------------------------------------------
// Purpose       : Find the bounding box of the Shell
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/27/99
//-------------------------------------------------------------------------
CubitBox Shell::bounding_box()
{
	CubitBox result;
	DLIList<RefFace*> face_list;
	ref_faces( face_list );
	if( face_list.size() > 0 )
		result = face_list.get_and_step()->bounding_box();
	for( int i = face_list.size(); i > 1; i-- )
		result |= face_list.get_and_step()->bounding_box();
	return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get ShellSM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/23/03
//-------------------------------------------------------------------------
ShellSM* Shell::get_shell_sm_ptr() const
{
  return dynamic_cast<ShellSM*>(bridge_manager()->topology_bridge());
}

//-------------------------------------------------------------------------
// Purpose       : Check if shell is a sheet.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/15/04
//-------------------------------------------------------------------------
CubitBoolean Shell::is_sheet()
{
  DLIList<RefFace*> faces;
  ref_faces(faces);
  while (faces.size())
    if ( ! faces.pop()->is_nonmanifold(this) )
      return CUBIT_FALSE;
  return CUBIT_TRUE;
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

