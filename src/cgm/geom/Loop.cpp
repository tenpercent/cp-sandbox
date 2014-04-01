//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Loop.C
//
// Purpose       : This file contains the implementation of the class
//                  Loop. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
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

#include "Loop.hpp"

#include "CoEdge.hpp"
#include "SenseEntity.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"
#include "CubitUtil.hpp"

#include "GeometryQueryEngine.hpp"

#include "DLIList.hpp"
#include "CastTo.hpp"

#include "Curve.hpp"
#include "LoopSM.hpp"

#include "ModelEntity.hpp"

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
Loop::Loop()
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
Loop::Loop(LoopSM* OSMEPtr)
{
   set_topology_bridge(OSMEPtr) ;
}

//-------------------------------------------------------------------------
// Purpose       : Gets the angle metric for this Loop.
//                 
// Special Notes : The actual computation is done by the underlying geometric
//                 modeling engine as this is a geometric computation.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 1/10/97
//-------------------------------------------------------------------------
CubitStatus Loop::get_angle_metric(double& angle_metric)
{
    // If the OSME knows its angle metric (some modelers cache
    // some data along these lines), return that value
  LoopSM* loop_sm = get_loop_sm_ptr();
  if (loop_sm && loop_sm->get_angle_metric(angle_metric))
    return CUBIT_SUCCESS;
  
    // First, get the list of coedges
  DLIList<CoEdge*> coedges;
  ordered_co_edges(coedges);
  
    // Now build a polygon approximation of the curves in the loop.
    // The polygon is a straight-line approximation to and is
    // topologically equivalent to the loop. At any given point,
    // the polygon may be far from the actual loop, however.
    //  - samitch
  
    // If there are no coedges, this is an empty loop
  if (coedges.size() == 0)
  {
    angle_metric = 0;
    return CUBIT_FAILURE;
  }
  
    // Loop through each coedge
  int i, j;
  DLIList<CubitVector*> polygon_points;
  DLIList<CubitVector*> interior_points;
  CoEdge* cur_coedge = NULL;
  coedges.reset();
  for (i = coedges.size(); i--; )
  {
      // Get the first point on this curve
    cur_coedge = coedges.get_and_step();
    RefEdge* cur_refedge = cur_coedge->get_ref_edge_ptr();
    polygon_points.append(new CubitVector
                          (cur_coedge->get_sense() == CUBIT_FORWARD ?
                           cur_refedge->start_vertex()->coordinates() :
                           cur_refedge->end_vertex()->coordinates()));
    
      // Get the interior points for approximation
    CubitSense return_sense;
    interior_points.clean_out();
    cur_refedge->get_interior_extrema(interior_points, return_sense);
      // Now put the points into the polygon.
      // We don't need to re-allocate any CubitVectors because we are just
      // copying pointers to dynamically allocated CubitVectors.
    if (cur_coedge->get_sense() == return_sense)
    {
      interior_points.reset();
      for (j = interior_points.size(); j--; )
        polygon_points.append(interior_points.get_and_step());
    }
    else
    {
      interior_points.last();
      for (j = interior_points.size(); j--; )
        polygon_points.append(interior_points.get_and_back());
    }
  }
  
    // Now that we have all of the points, compute and sum up
    // the internal angles on the polygon approximation.
  double angle, angle_sum = 0;
  RefFace* surface = cur_coedge->get_ref_face();
  CubitVector *point[3], t[2], normal;
  point[0] = polygon_points.get_and_step();
  point[1] = polygon_points.get_and_step();
  t[0] = *point[1] - *point[0];
  for (i = polygon_points.size(); i--; )
  {
      // Determine proper internal surface angle at point[1]
    point[2] = polygon_points.get_and_step();
    normal = surface->normal_at(*point[1]);
    t[1] = *point[1] - *point[2] ;
    angle = normal.vector_angle(t[1], t[0]);
    
      // Add up the total
    angle_sum += angle;
    
      // Iterate
    point[1] = point[2];
    t[0] = -t[1];
  }
  angle_metric = angle_sum / CUBIT_PI - polygon_points.size();
  
    // Clean up dynamically allocated vectors
  for (i = polygon_points.size(); i>0; i--) 
    delete polygon_points.get_and_step();
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : gets the ref edges in order with respect to this Loop.
//                 
// Special Notes : appends the edges to the list, with respect to coedges.
//
// Creator       : David White
//
// Creation Date : 03/25/97
//-------------------------------------------------------------------------
CubitStatus Loop::ordered_ref_edges(DLIList<RefEdge*>& ordered_edge_list )
{
  CubitStatus status = CUBIT_SUCCESS;

  DLIList<SenseEntity*> sense_entity_list;

  status = this->get_sense_entity_list(sense_entity_list);
  if ( status == CUBIT_FAILURE )
  {
    PRINT_ERROR("In Loop::ordered_ref_edges\n");
    PRINT_ERROR("       Problem getting the CoEdges of this Loop.\n");
    return CUBIT_FAILURE;
  }
  
    //Get the ref_edges associated with each co_edge.
  for ( int ii = sense_entity_list.size(); ii > 0; ii-- )
  {
    SenseEntity* se_ptr = sense_entity_list.get_and_step();
    BasicTopologyEntity* bte_ptr = se_ptr->get_basic_topology_entity_ptr();
    ordered_edge_list.append( dynamic_cast<RefEdge*>(bte_ptr) );
  }

  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------------
// Purpose       : gets the ref edges in order with respect to this Loop.
//                 
// Special Notes : appends the edges to the list, with respect to coedges.
//
// Creator       : David White
//
// Creation Date : 03/25/97
//-------------------------------------------------------------------------
CubitStatus Loop::ordered_co_edges(DLIList<CoEdge*>& ordered_coedge_list )
{
  CubitStatus status = CUBIT_SUCCESS;
  
  DLIList<SenseEntity*> sense_entity_list;
  
  status = this->get_sense_entity_list(sense_entity_list);
  if ( status == CUBIT_FAILURE )
  {
    PRINT_ERROR("In Loop::ordered_co_edges\n"
                "       Problem getting the CoEdges of this Loop.\n");
    return CUBIT_FAILURE;
  }
  
    // Cast the SenseEntity list to a CoEdge list
  DLIList<CoEdge*> co_edge_list(sense_entity_list.size());
  CAST_LIST( sense_entity_list, co_edge_list , CoEdge);
  
  ordered_coedge_list += co_edge_list;
  return CUBIT_SUCCESS;
}



RefFace* Loop::get_ref_face_ptr()
{
	return CAST_TO( get_basic_topology_entity_ptr(), RefFace );
}

//-------------------------------------------------------------------------
// Purpose       : Get LoopSM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/23/03
//-------------------------------------------------------------------------
LoopSM* Loop::get_loop_sm_ptr() const
{
  return dynamic_cast<LoopSM*>(bridge_manager()->topology_bridge());
}

//-------------------------------------------------------------------------
// Purpose       : Check if loops are spatially equal.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/01/04
//-------------------------------------------------------------------------
CubitBoolean Loop::about_spatially_equal( DLIList<CoEdge*>& other_coedges,
                                          CubitSense relative_sense,
                                          double tolerance_factor,
                                          CubitBoolean notify_refEntity )
{
  DLIList<CoEdge*> this_coedges(other_coedges.size());
  
    // Loops must have same number of coedges to match.
  this->ordered_co_edges( this_coedges );
  if (this_coedges.size() != other_coedges.size())
    return CUBIT_FALSE;
  
    // Want to compare coedges in order, so make sure we have
    // them in the correct order.
  if (relative_sense == CUBIT_REVERSED)
    this_coedges.reverse();
  
    // Try to match all coedges.  Begin with the first coedge
    // in this loop.  For each coedge in the other loop that 
    // it matches, check if all the other coedges match in the
    // correct order.
  int other_loop_index = 0;
  this_coedges.reset();
  other_coedges.reset();
  CoEdge* this_coedge = this_coedges.get_and_step();
  for (int i = other_coedges.size(); i--; )
  {
      // Loop until we find a matching CoEdge
    CoEdge* other_coedge = other_coedges.get_and_step();
    if (!this_coedge->about_spatially_equal( other_coedge,
                                             relative_sense,
                                             tolerance_factor,
                                             notify_refEntity ))
      continue;
    
      // Found a matching coedge.  Now try to match all the
      // others in the correct order.
    bool match = true;
    other_loop_index = other_coedges.get_index();
    for (int j = other_coedges.size() - 1; j-- && match; )
    {
      this_coedge = this_coedges.get_and_step();
      other_coedge = other_coedges.get_and_step();
      match = this_coedge->about_spatially_equal( other_coedge,
                                                  relative_sense,
                                                  tolerance_factor,
                                                  notify_refEntity );
    }
    
      // Matched all coedges, in order.  Done.
    if (match)
      return CUBIT_TRUE;
    
     // Try again, as perhaps the first coedge of this loop
     // also matches some other one in the second loop and
     // if we start with that one, the remaining coedges will
     // also match.
    this_coedges.reset();
    this_coedge = this_coedges.get_and_step();
    other_coedges.reset();
    other_coedges.step( other_loop_index );
  }
  
    // If here, loops didn't match.
  return CUBIT_FALSE;
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


