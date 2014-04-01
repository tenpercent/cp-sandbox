//-------------------------------------------------------------------------
// Filename      : Loop.hpp
//
// Purpose       : Represents a loop of edges of a face of a model.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef LOOP_HPP
#define LOOP_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitBox.hpp"
#include "GroupingEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class CoEdge;
class RefFace;
class RefEdge;
class RefVertex;
class CubitVector;
class LoopSM;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT Loop : public GroupingEntity
{
public :
  
  Loop() ;
    //- The default constructor.
  
  Loop(LoopSM* OSMEPtr) ;
    //- The constructor with a pointer to an other solid model entity.

  DagType dag_type() const { return DagType::loop_type(); }
  
  CubitBoolean is_external() ;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Returns CUBIT_TRUE if the Loop is an external Loop and CUBIT_FALSE
    //- otherwise.
  
  CubitStatus get_angle_metric(double& angle_metric);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O angle_metric
    //O- Angle metric value
    //- Gets the angle metric for this Loop. The value of the angle  
    //- metric returned is
    //- ( sum(interior_angles_at_vertices)/pi - n ).
    //- For planar faces this metric is 2.0 for internal
    //- Loops and -2.0 for
    //- external Loops.
  
  CubitStatus ordered_ref_edges( DLIList<RefEdge*>& ordered_edge_list );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O ordered_edge_list
    //O- Ordered list of ref edges with respect to this loop.
    //- This function works through the co_edges to get the correct
    //- ref_edges with respect to this loop.

  CubitStatus ordered_co_edges (DLIList<CoEdge*>& ordered_coedge_list);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O ordered_edge_list
    //O- Ordered list of co edges with respect to this loop.
    //- This function is differnt than the simple co_edges function.
    //- You must use this function if you want the co_edges to respect
    //- the loop in their list.
		
	RefFace* get_ref_face_ptr( );
	//R RefFace*
	//R- A pointer to the parent RefFace.
	
	CoEdge* get_co_edge_ptr( RefEdge* ref_edge_ptr );
	//R CoEdge* 
	//R- A pointer to a child CoEdge
	//I ref_edge_ptr
	//I- A pointer to a child RefEdge
	//- This method returns the child CoEdge connecting the
	//- passed RefEdge to this Loop.  If the RefEdge occurs
	//- more than once in the Loop, only the first CoEdge
	//- encountered will be returned.
	
      
  CubitBox bounding_box();
    //R CubitBox
    //R- The bounding box of this loop.
    //- This method returns the bounding box of a Loop
    //- as the union of the bounding boxes of its child
    //- RefEdges.

	CubitStatus tangents( RefVertex* vertex_ptr, 
		                    CubitVector& first_vector,
	                      CubitVector& second_vector ) ;
		//R CubitStatus
		//I veretx_ptr
		//I- A pointer to a vertex in this loop.
		//O first_vector
		//O- For the first edge on the loop which contains the passed
		//O- vertex, this is (an approximation of) the tangent of that
		//O- edge at the location of the vertex, pointing outwards from
		//O- the vertex.
		//O second_vector
		//O- For the second edge on the loop which contains the passed
		//O- vertex, this is (an approximation of) the tangent of that
		//O- edge at the location of the vertex, pointing outwards from
		//O- the vertex.
		//- Get the tangent vectors of the loop at the location of the
		//- passed vertex, pointing outwards from the vertex.
		//- 
		//- NOTE:  Both returned vectors point outward from the vertex.
		//-        To have both tangents pointing in the direction of
		//-        the loop, reverse first_vector.
	
  CubitBoolean validate() ;
		//R CubitBoolean
		//R- CUBIT_TRUE if the Loop is topologically valid.
		//R- CUBIT_FALSE otherwise.
		//- Validate the Loop.
		//- Returns CubitFalse if:
		//-  * The Loop has no CoEdges.
		//-  * Any CoEdge has no child RefEdge.
		//-  * Any CoEdge has an undefined sense (sense is CUBIT_UNKNOWN).
		//-  * The Loop is not topologically closed.
		//-  * The sense of a CoEdge does not correspond to
		//-    the order of the vertices on its child RefEdge.
		//-  * The Loop has two child RefEdges which form a topological
		//-    figure-eight (there is only one child vertex of the Loop.)

  LoopSM* get_loop_sm_ptr() const;
  
    CubitBoolean about_spatially_equal( DLIList<CoEdge*>& other_loop_coedges,
                                        CubitSense relative_sense,
                                        double tolerance_factor = 1.0,
                                        CubitBoolean notify_refEntity 
                                           = CUBIT_FALSE );
                                        

   protected: 

   private:
    Loop( const Loop& );
    void operator=( const Loop&);
} ;

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

