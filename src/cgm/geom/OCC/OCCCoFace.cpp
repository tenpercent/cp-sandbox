//-------------------------------------------------------------------------
// Filename      : OCCCoFace.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 04/28/08
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********
// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CastTo.hpp"
#include "OCCCoFace.hpp"
#include "OCCLoop.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitUtil.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********


OCCCoFace::OCCCoFace( OCCSurface *surf_ptr, 
	  	      OCCShell *shell_ptr, CubitSense sense )
	 	    : mySurface(surf_ptr), 
 	  	      myShell(shell_ptr),faceSense(sense)
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Jane HU  
//
// Creation Date : 04/28/08
//-------------------------------------------------------------------------
OCCCoFace::~OCCCoFace()
{
}

GeometryQueryEngine*
OCCCoFace::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}


void OCCCoFace::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  parents.append(shell());
}

void OCCCoFace::get_children_virt( DLIList<TopologyBridge*>& children )
{
  children.append(surface());
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

