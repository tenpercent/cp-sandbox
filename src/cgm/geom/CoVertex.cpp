//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : CoVertex.cc
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

#include "CoVertex.hpp"
#include "RefVertex.hpp"
#include "DLIList.hpp"
#include "Chain.hpp"
#include "RefEdge.hpp"

#include "CastTo.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********


// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The default constructor
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------

CoVertex::CoVertex() 
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------

CoVertex::~CoVertex() 
{
}

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
CoVertex::CoVertex(RefVertex* vertexPtr) 
{
   attach_basic_topology_entity(vertexPtr) ;
   set_sense(CUBIT_FORWARD) ;
}

//-------------------------------------------------------------------------
// Purpose       : Get the RefVertex associated with this CoVertex
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
RefVertex* CoVertex::get_ref_vertex_ptr()  
{
   // Call the generic function defined in the SenseEntity class to
   // do the real work
   BasicTopologyEntity* BTEPtr = get_basic_topology_entity_ptr();

   // Cast the returned pointer to RefVertex and return it
   return CAST_TO( BTEPtr, RefVertex );   
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
