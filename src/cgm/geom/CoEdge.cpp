//-------------------------------------------------------------------------
// Filename      : CoEdge.cpp
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

// ********** BEGIN CUBIT INCLUDES         **********
#include "CoEdge.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "Loop.hpp"
#include "DLIList.hpp"
#include "GeometryDefines.h"
#include "GeometryQueryTool.hpp"
#include "CastTo.hpp"
#include "CoEdgeSM.hpp"
#include "Curve.hpp"

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
CoEdge::CoEdge() 
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
CoEdge::~CoEdge() 
{
  remove_from_observers();
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to a edge and the
//                 sense of this CoEdge. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
CoEdge::CoEdge(RefEdge* edgePtr, CubitSense sense) 
{
   attach_basic_topology_entity(edgePtr) ;
   set_sense(sense) ;
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CoEdge::CoEdge( CoEdgeSM* osme_ptr )
{
  set_co_edge_sm_ptr(osme_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Get the RefEdge associated with this CoEdge.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
RefEdge* CoEdge::get_ref_edge_ptr()  
{
   // Call the generic function defined in the SenseEntity class to
   // do the real work
   BasicTopologyEntity* BTEPtr = get_basic_topology_entity_ptr();

   // Cast the returned pointer to RefEdge and return it
   return CAST_TO( BTEPtr, RefEdge );   
}

//-------------------------------------------------------------------------
// Purpose       : Get the parent Loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/23/99
//-------------------------------------------------------------------------
Loop* CoEdge::get_loop_ptr()
{
	return CAST_TO( get_grouping_entity_ptr(), Loop );
}


//-------------------------------------------------------------------------
// Purpose       : Get the ref_face associated with this CoEdge.
//
// Special Notes : Will assert if more than one ref_face but will
//                 return NULL if no RefFaces are associated with it.
//
// Creator       : David White
//
// Creation Date : 3/14/97
//-------------------------------------------------------------------------
RefFace* CoEdge::get_ref_face()
{
  return dynamic_cast<RefFace*>(get_parent_basic_topology_entity_ptr());
}
  


//-------------------------------------------------------------------------
// Purpose       : This function is called after a child of a ModelEntity 
//                 is switched. The sense of a CoEdge may change if one 
//                 of its RefEdges changes. This function takes care of 
//                 that. If the sense of the RefEdges that were switched 
//                 is same, nothing is done. If the RefEdges are of 
//                 opposite sense, the sense of this object is switched, 
//                 i.e. if it was FORWARD, it is made REVERSE, and vice
//                 versa. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------

void CoEdge::switch_child_notify(ModelEntity const* newChild, 
                                 ModelEntity const* oldChild) 
{
    // Make sure the entities being switched are RefEdges. If not,
    // get out of this function.

  ModelEntity* tmp_new_child = const_cast<ModelEntity*>(newChild);
  ModelEntity* tmp_old_child = const_cast<ModelEntity*>(oldChild);
  RefEdge * newChildRefEdge = CAST_TO(tmp_new_child, RefEdge) ;
  RefEdge * oldChildRefEdge = CAST_TO(tmp_old_child, RefEdge) ;
  
  if ( ( newChildRefEdge == NULL ) || ( oldChildRefEdge == NULL ) )
  {
    return ;
  }

    // If the children are RefEdges, get the sense of the old RefEdge 
    // relative to the new one.

  CubitSense sense;
  CubitBoolean spatially_equal;
    //We really don't care at this point if the two are spatially equal.  
    // The user could have done a force merge and in which case they are not.
    // So just ignore that spatial comparison.

//When this function is called, the merge has already taken place, so force
//  merge actually makes more sense to use here...
  CubitBoolean force_merge = CUBIT_TRUE;

  /*CubitStatus stat = */
  newChildRefEdge->relative_sense( 
     oldChildRefEdge, GeometryQueryTool::get_geometry_factor(), &sense,
     spatially_equal, force_merge );
  
    // If the sense of the old RefEdge relative to the new RefEdge is 
    // same, nothing needs to be done. However, if the relative sense
    // is reversed, switch the sense of the CoEdge.
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

//-------------------------------------------------------------------------
// Purpose       : Get CoEdgeSM pointer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/23/03
//-------------------------------------------------------------------------
CoEdgeSM* CoEdge::get_co_edge_sm_ptr() const
{
  return dynamic_cast<CoEdgeSM*>(bridge_manager()->topology_bridge());
}


CubitStatus CoEdge::set_co_edge_sm_ptr( CoEdgeSM* ptr )
{
  if (!bridge_manager()->topology_bridge())
    return bridge_manager()->add_bridge(ptr);
  else if(bridge_manager()->topology_bridge() != ptr)
    return CUBIT_FAILURE;
  else
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Check if child RefEdges are equal and have correct
//                 relative sense.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/01/04
//-------------------------------------------------------------------------
CubitBoolean CoEdge::about_spatially_equal( CoEdge* other_coedge,
                                            CubitSense relative_sense,
                                            double tolerance_factor,
                                            CubitBoolean notify_refEntity )
{
  RefEdge* this_ref_edge = get_ref_edge_ptr();
  RefEdge* other_ref_edge = other_coedge->get_ref_edge_ptr();
  CubitSense edge_sense;
  
  if (!this_ref_edge->about_spatially_equal( other_ref_edge,
                                             tolerance_factor,
                                             &edge_sense,
                                             notify_refEntity ))
    return CUBIT_FALSE;
  
  if (this_ref_edge->get_curve_ptr()->geometry_type() == POINT_CURVE_TYPE ||
      other_ref_edge->get_curve_ptr()->geometry_type() == POINT_CURVE_TYPE)
    return CUBIT_TRUE;
  
  if (edge_sense == CUBIT_UNKNOWN)
  {
    PRINT_WARNING("Failed to determine relative sense of curves.\n");
    return CUBIT_TRUE;
  }
  
  bool coedges_reversed = get_sense() != other_coedge->get_sense();
  bool want_reversed = edge_sense != relative_sense;
  if (coedges_reversed == want_reversed)
    return CUBIT_TRUE;
  
 // if (notify_refEntity)
 // {
 //   this_ref_edge->remove_compare_data();
 //   other_ref_edge->remove_compare_data();
 // }
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Get start/end vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/04
//-------------------------------------------------------------------------
RefVertex* CoEdge::start_vertex() 
{
  return get_sense() == CUBIT_FORWARD ? 
         get_ref_edge_ptr()->start_vertex() :
         get_ref_edge_ptr()->end_vertex();
}
RefVertex* CoEdge::end_vertex() 
{
  return get_sense() == CUBIT_REVERSED ? 
         get_ref_edge_ptr()->start_vertex() :
         get_ref_edge_ptr()->end_vertex();
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

