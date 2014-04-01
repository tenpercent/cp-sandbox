//-------------------------------------------------------------------------
// Filename      : GeometryUtil.hpp
//
// Purpose       : This class groups simple geometric utility functions.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/31/98
//-------------------------------------------------------------------------

#ifndef GEOMETRY_UTILITY_HPP
#define GEOMETRY_UTILITY_HPP

#include "CubitDefines.h"
#include "CubitGeomConfigure.h"


class Surface;
class Curve;
class Point;
class Loop;
class CubitVector;
class CoEdge;
class RefEdge;
class RefFace;
class RefVolume;
class RefVertex;
class Shell;
class Body;
template <class X> class DLIList;
class CoFace;
class RefEdge;
class TopologyEntity;
class TopologyEntity;
class RefEntity;
class TopologyBridge;

class CUBIT_GEOM_EXPORT GeometryUtil
{
public:
  
  double default_angle_tol;     //CUBIT_PI      == NONE
  double default_midpoint_tol;  //CUBIT_DBL_MAX == NONE
  
  static GeometryUtil* instance()
    {
      if( instance_ == NULL ) instance_ = new GeometryUtil;
      return instance_;
    }
  
  CubitBoolean is_position_within_loop( const CubitVector& position,
                                        Loop* loop_ptr,
                                        CubitVector* closest_on_loop = NULL,
                                        CubitVector* normal = NULL );
    //R CubitBoolean    
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I position
    //I- The position to test.  This position will be moved to the
    //I- surface if it does not lie on the surface.
    //I loop_ptr
    //I- A pointer to the loop to test.
    //O closest_on_loop
    //O- The closest point on the loop.
    //I normal
    //I- A pointer to a vector containing the normal on the surface at
    //I- the passed position.  If a normal is passed, the function assumes
    //I- that the point is on the surface.  
    //I- If no normal is specified, the function will query the loop_ptr 
    //I- for its parent face, moved the passed position to that face, and
    //I- evaluate the normal at that position.
    //- This method tests if a point lies within a loop.
    //-
    //- NOTE: This method will return the expected value only when
    //-       the loop is an outside loop!  If the loop represents
    //-       a hole, CUBIT_TRUE will be returned when the point
    //-       is OUTSIDE the loop. (or on the bounded RefFace)
  
  
  CubitStatus closest_point_trimmed( RefFace* face_ptr, 
                                     const CubitVector& from_pt,
                                     CubitVector& result_pt,
                                     CubitPointContainment* cont = 0 );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I face_ptr
    //I- The surface to evaluate on.
    //I from_pt
    //I- The position to evaluate from.
    //O result_pt
    //O- The closest position on the surface.
    //- This is an implementation of closest_point_timmed which uses
    //- the modeling engine's closest_point method, but uses
    //- loop angle metric to determine if a point lies in a loop.

    double surface_cpu_time;
    double other_cpu_time;
    //Stats kept by closest_point_trimmed.
  
  CubitBoolean is_point_on_face( RefFace* ref_face_ptr,
                                 const CubitVector& point,
                                 CubitBoolean check_on_surface = CUBIT_TRUE );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I ref_face_ptr
    //I- A pointer to the RefFace to test the point on
    //I point
    //I- The point to test
    //I check_on_surface
    //I- Check if the point is within CUBIT_RESABS of the surface.
    //- Test if the closest point on the passed RefFace to the passed
    //- point is within the boundary of the RefFace, and if 
    //- check_on_surface is true, check that the passed point is within
    //- CUBIT_RESABS of the surface.

  CubitStatus make_linearized_curve( Surface* surface_ptr,
                                     const CubitVector& start_pt,
                                     const CubitVector& end_pt,
                                     DLIList<CubitVector*>& segment_points,
                                     double arc_angle_tol = CUBIT_DBL_MAX,
                                     double midpoint_dist_tol = CUBIT_DBL_MAX,
                                     const CubitVector* mid_point = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I surface_ptr
    //I- A surface the curve is to lie on.
    //I start_pt
    //I- The start point of the curve.
    //I end_pt
    //I- The end point of the curve.
    //O segment_points
    //O- The curve segments.
    //I arc_angle_tol
    //I- The maximum angle between surface normal vectors at the start and
    //I- end of a segment in a plane tangent to the segment.
    //I midpoint_dist_tol
    //I- The maximum distance bewteen the midpoint of a segment
    //I- and the surface.  
    //- This method creates a linear segment
    //- beween the start and end points, and then moves that line
    //- onto the surface by splitting it, and moveing the split 
    //- points onto the surface.  Each segment is bisected and the split
    //- point is moved to the surface, until the tolerances are met.
  
  double loop_area( Loop* loop_ptr );
    //R double
    //R- A (very) rough approximation of the area of a loop.
    //- Do a very rough approximation of the area of a loop.

  CoEdge* closest_loop_coedge( Loop* loop_ptr, 
                               const CubitVector& from_pt,
                               CoEdge*& other_coedge,
                               CubitVector* closest_pt = NULL );
  CoEdge* closest_face_coedge( RefFace* face_ptr, 
                               const CubitVector& from_pt,
                               CoEdge*& other_coedge,
                               CubitVector* closest_pt = NULL );
    //R CoEdge*
    //R- The closest coedge
    //I loop_ptr
    //I- The loop to check the coedges of.
    //I from_pt
    //I- The position to check from
    //O other_coedge
    //O- If the closest point was at a vertex, the other coedge
    //O- at that vertex.  This coedge follows the returned coedge
    //O- in the loop.
    //O closest_pt
    //O- The closest point on the loop.
    //- This method finds the closest point on a loop, and returnes the
    //- coedge it occured on, and optionally the position on that coedge.
    
  CubitStatus closest_shell_coface( Shell* shell_ptr,
                                    const CubitVector& from_pt,
                                    DLIList<CoFace*>& result_set,
                                    CubitVector* closest_pt = NULL );
    //R CubitStatus
    //I shell_ptr
    //I- A pointer to the shell to evaluate.
    //I from_pt
    //I- The position to evaluate.
    //O result_set
    //O- The closest coface to the passed position, or several CoFaces
    //O- if they are equidistant from the position.
    //O closest_pt
    //O- The closest point on the returned coface.
    //- This method finds the closest point on a shell, and returnes the
    //- coface it occured on, and optionally the position on that coface.
    
  CubitBoolean is_position_within_shell( const CubitVector& position,
                                         Shell* shell_ptr,
                                         CubitVector* closest_on_shell = NULL );
    //R CubitBoolean    
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I position
    //I- The position to test.  
    //I shell_ptr
    //I- A pointer to the shell to test.
    //O closest_on_shell
    //O- The closest point on the shell.
    //- This method tests if a point lies within a shell.
    //-
    //- NOTE: This method will return the expected value only when
    //-       the shell is an outside loop!  If the shell represents
    //-       a void, CUBIT_TRUE will be returned when the point
    //-       is OUTSIDE the shell. (or within the bounded RefVolume)
  
  CubitBoolean is_position_in_volume( const CubitVector& position,
                                      RefVolume* ref_volume_ptr,
                                      CubitVector* closest_on_bounary = NULL );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I ref_volume_ptr
    //I- A pointer to the RefVolume to test the point on
    //I position
    //I- The point to test
    //- Test if the passed position is within the bounded RefVolume.
    
  CubitBoolean is_shell_a_void( Shell* shell_ptr );
    //R CubitBoolean
    //R- CUBIT_TRUE  if the passed shell is a void, or
    //R- CUBIT_FALSE if the passed shell is a bounded region.
    //I shell_ptr
    //I- A pointer to the Shell to test
    //- Test if a Shell is a void or an enclosed region.

  RefEdge* find_connected_ref_edge_by_coord( CubitVector& coords, 
                                             DLIList<RefEdge*>& ref_edge_list,
                                             int& start_flg );
    // R RefEdge - pointer to refedge containing coords as start or end location
    // R- NULL if no refedge found
    // I coords 
    // I- xyz position to search for
    // I ref_edge_list
    // I- list of refedges to search through
    // O start_flg
    // O- A flag indicating whether coords was on start or end of refedge
    // O-  1 = on start
    // O-  0 = on end

  void find_connected_ref_edges_by_coord( CubitVector& coords, 
                                          DLIList<RefEdge*>& ref_edge_list,
                                          DLIList<RefEdge*>& connected_edges_list );
    //- From a list of ref_edges, find all that connect by coordinate.

  CubitStatus form_ref_edge_loop_by_coord( DLIList<RefEdge*>& ref_edge_list );
    //- Sorts a list of RefEdges such that they form a loop.  The RefEdges
    //- need not be connected by vertices - the function uses coordinates
    //- to find connected edges.  The order is preserved if they already form 
    //- a loop.  Otherwise, the order is determined from the first->second 
    //- ref_edge, or if these are not connected, from the start->end of first 
    //- ref_edge.  Function returns CUBIT_FAILURE if no loop can be found.

  CubitStatus form_ref_edge_chain_by_coord( DLIList<RefEdge*>& ref_edge_list );
    //- Sorts a list of RefEdges such that they form a chain.  The order is
    //- preserved if they already form a chain.  Otherwise, the order is
    //- determined from the first->second ref_edge, or if these are not
    //- connected, from the end of the first ref_edge whose coordinates are
    //- not shared by other ref_edges in the list.
    //- Function returns CUBIT_FAILURE if no single chain can be found.

  CubitStatus check_valid_chain_by_coord( DLIList<RefEdge*>& ref_edge_list, 
                                          DLIList<RefVertex*>& problem_vertices );
  CubitStatus check_valid_loop_by_coord( DLIList<RefEdge*>& ref_edge_list,
                                         DLIList<RefVertex*>& problem_vertices );
    //- Functions to ensure that a valid chain or loop was found.  Checks
    //- number of connected edges to each vertex

  // - Same functions as above, but connections by vertex rather than coord
  RefEdge* find_connected_ref_edge( RefVertex* ref_vertex_ptr, 
                                    DLIList<RefEdge*>& ref_edge_list,
                                    int& start_flg );
  //- From a list of ref_edges, finds one that connects by vertex.
  //- List is searched sequentially.  Also returns whether the
  //- start or end of the found curve matches the coordinates.
  //- Returns NULL if no matching RefEdge found.

  void find_connected_ref_edges( RefVertex* ref_vertex_ptr, 
                                 DLIList<RefEdge*>& ref_edge_list,
                                 DLIList<RefEdge*>& connected_edges_list );
  //- From a list of ref_edges, find all that connect by vertex.
  
  CubitStatus form_ref_edge_chain( DLIList<RefEdge*>& ref_edge_list );
  //- Sorts a list of RefEdges such that they form a chain.  The order is
  //- preserved if they already form a chain.  Otherwise, the order is
  //- determined from the first->second ref_edge, or if these are not
  //- connected, from the end of the first ref_edge whose coordinates are
  //- not shared by other ref_edges in the list.
  //- Function returns CUBIT_FAILURE if no single chain can be found.

  CubitStatus form_ref_edge_loop( DLIList<RefEdge*>& ref_edge_list );
  //- Sorts a list of RefEdges such that they form a loop.  The RefEdges
  //- need not be connected by vertices - the function uses coordinates
  //- to find connected edges.  The order is preserved if they already form 
  //- a loop.  Otherwise, the order is determined from the first->second 
  //- ref_edge, or if these are not connected, from the start->end of first 
  //- ref_edge.  Function returns CUBIT_FAILURE if no loop can be found.

  CubitStatus check_valid_chain( DLIList<RefEdge*>& ref_edge_list, 
                                 DLIList<RefVertex*>& problem_vertices );
  CubitStatus check_valid_loop( DLIList<RefEdge*>& ref_edge_list,
                                DLIList<RefVertex*>& problem_vertices );
  //- Functions to ensure that a valid chain or loop was found.  Checks
  //- number of connected edges to each vertex
    
  CubitBoolean valid_edge( RefEdge* edge_ptr, 
                           CubitBoolean print_error = CUBIT_TRUE );
    //- Check if the vertices are in the correct order on the edge wrt
    //- to the specified start and end parameters of the edge.
    //- Check if the tangent vectors at the end vertices are in the
    //- correct direction.
  CubitBoolean valid_loop_coedges( Loop* loop_ptr,
                                   CubitBoolean print_error = CUBIT_TRUE );
    //- Check that the loop is closed and that the sense of each coedge
    //- is logical wrt the next coedge.
  CubitBoolean valid_shell_cofaces( Shell* shell_ptr,
                                    CubitBoolean print_error = CUBIT_TRUE );
    //- Check the relative sense of each coface wrt to its ajacent cofaces.
  CubitBoolean valid_topology( TopologyEntity* topo_ptr,
                               CubitBoolean print_error = CUBIT_TRUE,
                               DLIList<TopologyEntity*>* invalid_list = NULL );
    //- Call the above validation methods for the passed TopologyEntity 
    //- and all of its children.       
    
  CubitBoolean valid_sm_topology( DLIList<RefEntity*>& entity_list,
                                  CubitBoolean print_error = CUBIT_TRUE );
    //- Verify that the VGI DAG matches the solid modeler's topology
    //- graph.  Added for verifying correct unmerge.     
  CubitBoolean valid_sm_topology( Body*     body_ptr, CubitBoolean print );
    //- called by valid_sm_topology(DLIList<RefEntity*>&)
  CubitBoolean valid_sm_topology( RefVolume* vol_ptr, CubitBoolean print );
    //- called by valid_sm_topology(DLIList<RefEntity*>&)
  CubitBoolean valid_sm_topology( Shell*   shell_ptr, CubitBoolean print );
    //- called by valid_sm_topology(RefVolume*)
  CubitBoolean valid_sm_topology( RefFace*  face_ptr, CubitBoolean print );
    //- called by valid_sm_topology(DLIList<RefEntity*>&)
  CubitBoolean valid_sm_topology( Loop*     loop_ptr, CubitBoolean print );
    //- called by valid_sm_topology(RefFace*)
  CubitBoolean valid_sm_topology( CoEdge* coedge_ptr, CubitBoolean print );
    //- called by valid_sm_topology(Loop*)
  CubitBoolean valid_sm_topology( RefEdge*  edge_ptr, CubitBoolean print );
    //- called by valid_sm_topology(DLIList<RefEntity*>&)
  CubitBoolean valid_sm_topology( RefVertex* vtx_ptr, CubitBoolean print );
    //- called by valid_sm_topology(DLIList<RefEntity*>&)
                          
  static void list_SM_topology( TopologyBridge* bridge, int depth );
  
    
protected:

  static void list_SM_topology( TopologyBridge* bridge, int depth, 
                                int indent, int counts[8] );
  static int print_topo_bridge( TopologyBridge* bridge, int indent );

  CubitBoolean inside_of_curve( const CubitVector& curve_tangent,
                                const CubitVector& curve_position,
                                const CubitVector& surf_position,
                                const CubitVector& surf_normal );
  
  CubitStatus recursive_make_curve( Surface* surface_ptr,
                                    const CubitVector& start_pt,
                                    const CubitVector& end_pt,
                                    DLIList<CubitVector*>& segment_points,
                                    double arc_angle_tol,
                                    double midpoint_dist_tol );
    //R CubitStatus
    //I surface_ptr
    //I- The surface the curve should lie on.
    //I start_pt, end_pt
    //I- The end points of the curve.
    //O segment_points
    //O- The list of positions defining the linear polyline result.
    //I cos_arc_angle_tol
    //I- The arc angle tolerance.
    //I midpoint_dist_tol_sqr
    //I- The midpoint distance tolerance.
    //- This is the recursive method used by 
    //- GeometryUtil::make_linearized_curve() to create
    //- a segmented curve on a surface.
    //-
    //- This method creates a line between the start and end points,
    //- and moves that line onto the surface by recursively bisecting 
    //- the line and moving the midpoint to the surface. 
  
  GeometryUtil();
    //- constructor

private:
  
  static GeometryUtil* instance_;
  
  int linearized_curve_debug_count_;
  
};

#endif

