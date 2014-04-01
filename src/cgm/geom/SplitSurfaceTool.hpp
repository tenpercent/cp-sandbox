//-------------------------------------------------------------------------
// Filename      : SplitSurfaceTool.hpp
//
// Purpose       :
//
// Special Notes : 
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------

#ifndef SPLITSURFACETOOL_HPP
#define SPLITSURFACETOOL_HPP

#include "CubitDefines.h"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
class Cubit2DPoint;
class RefVertex;
class RefEdge;
class RefFace;
class TDSplitSurface;
//class Curve;
class GeometryModifyEngine;

template <class X> class DLIList;

class CUBIT_GEOM_EXPORT SplitSurfaceTool
{

  public:

    SplitSurfaceTool();
    ~SplitSurfaceTool(){}

    CubitStatus preview( RefFace *ref_face_ptr,
                         DLIList<CubitVector*> &locations,
                         DLIList<DLIList<CubitVector*>*> &vec_lists,
                         CubitBoolean create_ref_edges_flg );
    //- Preview function for simple surface splitting.  Temporary Curves are
    //- created on the RefFace and displayed.  The curves are created from 
    //- the input vec_lists (straight lines, arcs or splines).  The input 
    //- locations are the original input locations from the user, which are 
    //- displayed as well as the Curves.

    CubitStatus calculate_split_curves( RefFace *ref_face_ptr,
                                        DLIList<CubitVector*> &locations,
                                        DLIList<DLIList<CubitVector*>*> &vec_lists,
                                        DLIList<Curve*>& curve_list );
    //- Calculate the split curves based on the surface and user parameters.
    //- The curves are created from the input vec_lists (straight lines, 
    //- arcs or splines). The input locations are the original input locations 
    //- from the user.

    CubitStatus calculate_split_curves( DLIList<RefFace*> &ref_face_list, 
                                        int num_segs, double fraction, 
                                        double distance, RefEdge *from_curve_ptr,
                                        DLIList<RefVertex*> &corner_vertex_list,
                                        DLIList<RefVertex*> &through_vertex_list,
                                        RefEdge *curve_dir_ptr,
                                        CubitBoolean preview_flg,
                                        CubitBoolean create_ref_edges_flg,
                                        CubitBoolean just_curves_flg,
                                        DLIList<DLIList<Curve*>*> &curve_lists_list ); 
    //- Calculates the curves to split a chain of surfaces.
    //- ref_face_list - chain of surfaces to split.  Can be given in any order
    //-                 as long as they are connected.  Can be a continuous
    //-                 loop of surfaces.  If a single surface is given, by
    //-                 default it will be split along the narrowest aspect
    //-                 ratio of the surface (i.e., along a fillet).  Periodic
    //-                 surfaces are not handled.
    //- num_segs - the number of segments to create (must be >= 2 );
    //- fraction - the fraction along the surfaces to split, not valid if
    //-            num_segs > 2.  This value is not used if a distance is
    //-            specified instead.  
    //- distance - if 2 segments, allow the split to be at a user specified
    //-            distance across the surface.  Specify -1.0 to use the
    //-            fraction instead.
    //- from_curve_ptr - (OPTIONAL) if user specified a fraction or distance,
    //-                  orient from this curve.  If not specified, a default
    //-                  is chosen automatically.
    //- corner_vertex_list - (OPTIONAL) the user can specify the corners of the,
    //-                      chain, typically if the detected corners are incorrect.
    //-                      The split direction is from corners 0-1 to 
    //-                      corners 2-3 (for single surfaces, if not overridden
    //-                      by curve_dir_ptr).
    //- through_vertex_list - (OPTIONAL) user specifies forced vertices for the split
    //-                       to run through (on curves).  Not valid with
    //-                       more than 2 segments.
    //- curve_dir_ptr - (OPTIONAL) for single surfaces, the direction of split can be
    //-                 specified by picking a curve on the surface.
    //- preview_flg - if CUBIT_TRUE, just draw the curves that will be used to split 
    //-               instead of actually splitting.              
    //- create_ref_edges_flg - valid only if preview_flg=CUBIT_TRUE.  If CUBIT_TRUE,
    //-                    create RefEdges *instead* of splitting.
    //- Note this function also utilizes two additional settings which are 
    //- supplied through the set_tolerance and set_parametric_flg functions.    

    CubitStatus split_surface( RefFace *ref_face_ptr,
                               DLIList<CubitVector*> &locations,
                               DLIList<DLIList<CubitVector*>*> &vec_lists );
    //- Split function for simple surface splitting.  Temporary Curves are
    //- created on the RefFace and used to split the RefFace.  The curves are 
    //- created from the input vec_lists (straight lines, arcs or splines).
    //- The input locations are the original input locations from the user,
    //- which are displayed for user reference.

    CubitStatus split_surface( RefFace *ref_face_ptr, DLIList<Curve*> &curve_list);
    //- Split function for simple surface splitting using a list of predefined curves
    
    CubitStatus split_surfaces( DLIList<RefFace*> &ref_face_list,
                                int num_segs,
                                double fraction,
                                double distance,
                                RefEdge *from_curve_ptr,
                                DLIList<RefVertex*> &corner_vertex_list,
                                DLIList<RefVertex*> &through_vertex_list,
                                RefEdge *curve_dir_ptr = NULL,
                                CubitBoolean preview_flg = CUBIT_FALSE,
                                CubitBoolean create_ref_edges_flg = CUBIT_FALSE );
    //- Splits the chain of surfaces.
    //- ref_face_list - chain of surfaces to split.  Can be given in any order
    //-                 as long as they are connected.  Can be a continuous
    //-                 loop of surfaces.  If a single surface is given, by
    //-                 default it will be split along the narrowest aspect
    //-                 ratio of the surface (i.e., along a fillet).  Periodic
    //-                 surfaces are not handled.
    //- num_segs - the number of segments to create (must be >= 2 );
    //- fraction - the fraction along the surfaces to split, not valid if
    //-            num_segs > 2.  This value is not used if a distance is
    //-            specified instead.  
    //- distance - if 2 segments, allow the split to be at a user specified
    //-            distance across the surface.  Specify -1.0 to use the
    //-            fraction instead.
    //- from_curve_ptr - (OPTIONAL) if user specified a fraction or distance,
    //-                  orient from this curve.  If not specified, a default
    //-                  is chosen automatically.
    //- corner_vertex_list - (OPTIONAL) the user can specify the corners of the,
    //-                      chain, typically if the detected corners are incorrect.
    //-                      The split direction is from corners 0-1 to 
    //-                      corners 2-3 (for single surfaces, if not overridden
    //-                      by curve_dir_ptr).
    //- through_vertex_list - (OPTIONAL) user specifies forced vertices for the split
    //-                       to run through (on curves).  Not valid with
    //-                       more than 2 segments.
    //- curve_dir_ptr - (OPTIONAL) for single surfaces, the direction of split can be
    //-                 specified by picking a curve on the surface.
    //- preview_flg - if CUBIT_TRUE, just draw the curves that will be used to split 
    //-               instead of actually splitting.              
    //- create_ref_edges_flg - valid only if preview_flg=CUBIT_TRUE.  If CUBIT_TRUE,
    //-                    create RefEdges *instead* of splitting.
    //- Note this function also utilizes two additional settings which are 
    //- supplied through the set_tolerance and set_parametric_flg functions.    

    CubitStatus split_surfaces( DLIList<RefFace*> &ref_face_list,
                                int num_segs,
                                double fraction,
                                double distance,
                                RefEdge *from_curve_ptr,
                                DLIList<RefVertex*> &corner_vertex_list,
                                DLIList<RefVertex*> &through_vertex_list,
                                RefEdge *curve_dir_ptr,
                                DLIList<DLIList<Curve*>*> &curve_lists_list );
    //- Same as above except this method simply returns a list of lists of 
    //- Curves (one list per input RefFace).  This method is intended for
    //- use by a programmer (rather than a user from GeometryCommands).
    //- Note it is the calling functions responsibility to free the memory
    //- allocated in curve_lists_list, as well as the Curves (note this can
    //- be accomplished with the function free_curves_lists).

    CubitStatus split_surfaces( DLIList<RefFace*> &ref_face_list,
                                int num_segs,
                                double fraction,
                                double distance,
                                RefEdge *from_curve_ptr,
                                DLIList<RefVertex*> &corner_vertex_list,
                                DLIList<RefVertex*> &through_vertex_list,
                                RefEdge *curve_dir_ptr,
                                DLIList<Curve*> &curve_list );
    //- Same as above except this method simply returns a simple list of Curves.

    void free_curves_lists( DLIList<DLIList<Curve*>*> &curve_lists_list,
                            CubitBoolean free_curves_flg = CUBIT_TRUE );
    //- Free the curves and lists memory.  If free_curves_flg is CUBIT_TRUE,
    //- free the Curves as well (however, the function checks if a RefEdge is 
    //- attached to the Curve - if a RefEdge is attached the Curve is not freed).

    static void set_tolerance( double tol );
    //- Set the tolerance (function returns previous tolerance)
    static double get_tolerance();
    //- Get the tolerance.  In layman's terms, the tolerance controls how 
    //- closely the split actually follows the center of the surface. In the
    //- code, it determines how accurately to facet the bounding curves (the 
    //- facets are used to interpolate points on the surface using a mapping 
    //- concept).  It is also used in an iterative procedure in fitting splines
    //- to the line segments generated by the mapping (just fitting a spline 
    //- through the interpolated points was found to give very bad results - 
    //- instead we check the deviation of the spline from the line segments and
    //- add additional points until the resultant spline is within tolerance to
    //- the line segments.

    static void set_auto_detect_triangles_flg( CubitBoolean auto_detect_flg );
    //- Set the auto detect triangles flag
    static CubitBoolean get_auto_detect_triangles_flg();
    //- Get the auto detect triangles flag

    static void set_side_angle_threshold( double angle );
    //- Set the side angle threshold
    static double get_side_angle_threshold();
    //- Get the side angle threshold.  If there is a corner angle within this
    // threshold to 180 and another corner angle less than the point threshold, a
    // triangle is created.

    static void set_point_angle_threshold( double tol );
    //- Set the point angle threshold
    static double get_point_angle_threshold();
    //- Get the point angle threshold.  If there is a corner angle less than
    // this value and another corner angle within the side threshold to 180, a
    // triangle is created.

    static void set_parametric_flg( CubitBoolean parametric_flg );
    //- Set the parametric flag
    static CubitBoolean get_parametric_flg();
    //- Get the parametric flag.  If the parametric_flg is CUBIT_TRUE, find
    //- spline locations in the parametric space of the surface, otherwise
    //- initially find the spline locations in 3D space then project back to
    //- the surface.  Typically, the parametric space gives better results
    //- (on curvy surfaces, it will result in a spline that is much closer to
    //- the middle of the surface).  However, sometimes the parameter space
    //- is bad or the mapping algorigm gets confused (frequently on conic
    //- surfaces, resulting in points revolved 180 deg away), so the 3D method
    //- gives a better result.  Note the default is thus 3D.

    static void initialize_settings();
    //- Initialize settings for this class

    CubitStatus draw_preview( DLIList<Curve*> &curve_list, int color=CUBIT_BLUE );
    //- preview the curves

  private:

    CubitStatus split_surfaces( DLIList<RefFace*> &ref_face_list,
                                int num_segs,
                                double fraction,
                                double distance,
                                RefEdge *from_curve_ptr,
                                DLIList<RefVertex*> &corner_vertex_list,
                                DLIList<RefVertex*> &through_vertex_list,
                                RefEdge *curve_dir_ptr,
                                CubitBoolean preview_flg,
                                CubitBoolean create_ref_edges_flg,
                                CubitBoolean just_curves_flg,
                                DLIList<DLIList<Curve*>*> &curve_lists_list );
    //- Main workhorse function for automatic splitting.  See above for 
    //- description of inputs.  Note if just_curves_flg==CUBIT_TRUE, the 
    //- preview_flg and create_ref_edges_flg are ignored.  See above for
    //- description of the inputs.

    CubitStatus check_valid_faces();
    //- Check each face for multiple loops, hardlines and against looping on
    //- themselves (if more than one face).

    CubitStatus order_face_list();
    //- Order the refFaceChain from one end of surface chain to the other.
    //- This function also detects the isLoop condition (a continuous loop
    //- of surfaces).

    CubitStatus get_neighbors( RefFace *seed_ref_face, 
                               DLIList<RefFace*> &input_ref_faces,
                               DLIList<RefFace*> &neighbor_ref_faces );
    //- Get all the neighbors of the input seed RefFace.  Note second
    //- list is intended to be copied in.

    CubitStatus get_outer_loops();
    //- Get all of the outer loops - outerCurveLoop, outerCoEdgeLoop,
    //- outerVertexList.  Make sure they start at the start surface
    //- in the chain.

    CubitStatus get_outer_coedge_loop();
    //- Given the ordered refFaceChain, finds the outerCoEdgeLoop.  

    CubitStatus find_loop_start( CoEdge *&start_co_edge_ptr );
    //- Finds the start of the loop when constructing the outerCoEdgeLoop

    CubitStatus get_outer_vertices( DLIList<RefVertex*> &ref_vertex_list );
    //- Fill the given vertex list using the outerCoEdgeLoop.  This is
    //- simply an ordered list of the vertices in the loop.

    CubitStatus get_outer_curves( DLIList<RefEdge*> &ref_edge_list  );
    //- Fill the given refedge list using the outerCoEdgeLoop.  This is 
    //- simply an ordered list of refedges in the loop.
    
    CubitStatus pick_4_corners();
    //- For a chain of surfaces.  Uses angle angle criteria to pick the 
    //- corner vertices so that a mapping algorithm can be used to find 
    //- the split locations on the surfaces.  Two vertices must be on
    //- the start surface, and two on the end.

    CubitStatus pick_expanded_corners( RefFace *ref_face_ptr,
                                       DLIList<CoEdge*> &co_edge_list,
                                       int &offset, double epsilon, 
                                       int &A_i, int &B_i );
    //- Used by pick_4_corners.  Finds the two best corners on the given
    //- RefFace (which will be either the start or the end surface of the
    //- chain).  If there are multiple vertices with the same angle, it
    //- will find the pair that are the farthest from the center.

    CubitStatus pick_4_corners_simple();
    //- For a single surface.  Just finds the 4 corners closest to PI/2.

    CubitStatus update_corners_for_triangle();
    //- Determines whether any corners are close to 180 degrees (above
    //- angleTolerance).  On a single surface, the corner that is 
    //- closest to 180 and above angleTolerance is removed.  The next
    //- angle farthest from PI/2 is chosen as the corner of the triangle.
    //- For a chain, both ends of the chain are checked.

    CubitStatus remove_corner( int corner_to_remove, int collapsed_corner,
                               CubitBoolean set_collapsed_first );
    //- Removes a corner making a triangle.  If set_collapsed_first is
    //- CUBIT_TRUE, then the collapsed side is set to the first corner,
    //- meaning the split will originate there.

    CubitStatus fill_corners( int best_corner_1, int best_corner_2,
                              int best_corner_3, int best_corner_4 );
    //- Fill the cornerCoEdge array given the position of the corners

    CubitStatus fill_side_intervals( int best_corner_1, int best_corner_2,
                                     int best_corner_3, int best_corner_4 );
    //- Fill sideInterval array given teh position of the corners

    CubitStatus order_selected_corners( DLIList<RefVertex*> &outer_vertex_list );
    //- If the user specified corners, this function will reorient the lists
    //- so that the desired split direction will occur.  It also does error
    //- checking on the users selected corners.

    int number_in_list( DLIList<RefVertex*> &corner_vertex_list, 
                        RefVertex *ref_vertex_ptr );
    //- Count the number of times the given RefVertex exists in the list

    CubitBoolean is_in_outer_loop( RefVertex *ref_vertex_ptr );
    //- Determines whether the given vertex exists in the outer loop

    CubitStatus compute_angles( double *&angles, double &turn_angle_sum );
    //- This finds an array of angles around the outer loop for finding
    //- corners.  Note the angles array is allocated - the calling function
    //- must free it.  turn_angle_sum is currently not used.

    void order_corners( int &corner_1, int &corner_2,
				                int &corner_3, int &corner_4 );
    void order_corners( int &corner_1, int &corner_2,
				                int &corner_3 );
    //- Given the best corners in any order, these functions just sort them
    //- from lowest to highest.

    CubitStatus adjust_for_split_direction( RefEdge *curve_dir_ptr,
                                            RefEdge *from_curve_ptr,
                                            DLIList<RefVertex*> &corner_vertex_list );
    //- Adjust the lists so that the split direction is proper.  See function
    //- for logic details (i.e., always split along chains of surfaces, attempt
    //- to split along highest aspect ratio for single surfaces, etc.).

    double get_side_length( int side );
    //- Determines the length of a given side of the logical rectangle.

    double compute_next_angle( DLIList<CoEdge*> &co_edge_list );
    //- Finds the next angle when walking around outerCoEdgeLoop finding 
    //- angles for corner picking.  It steps the list.

    CubitStatus reorient_loop( int start_offset );
    //- Reorients the global loop (outerCoEdgeLoop) so that it starts at the beginning.

    CubitStatus shuffle_corners_forward();
    //- This shuffles the outer loop, side interval and corner array forward
    //- by sideInterval[0] (ie., useful when changing the split direction for
    //- a single surface).

    CubitStatus shuffle_zero_side_forward();
    //- This shuffles the zero interval side (present on triangles) forward.
    //- Side interval and corner array affected.  This function is necessary
    //- if user specifies a 'direction' or 'from' curve when splitting a triangle.

    CubitStatus populate_curve_loops();
    //- Populates the CoEdge lists in the surface tooldatas.  The CoEdge lists
    //- are organized by side of the logical rectangle on each surface in
    //- refFaceChain.

    CubitBoolean is_edge_in_list( RefEdge *ref_edge_ptr, DLIList<CoEdge*> &co_edge_list );
    //- Determines if the given RefEdge is the child of any of the CoEdges in co_edge_list

    CubitStatus check_face_loop( RefFace *ref_face_ptr, CubitBoolean &is_loop );
    //- Checks a single surface to see if it is a loop.

    CubitStatus check_for_loop( RefFace *ref_face_ptr1, RefFace *ref_face_ptr2, 
                                CubitBoolean &is_loop, 
                                CoEdge *&start_co_edge_ptr );
    //- Checks a two-surface chain to see if it is a loop.  If it is a loop, it can
    //- return the starting coedge as well (warning - if not a loop, the starting 
    //- coedge is not set).

    CubitStatus ordered_co_edges( RefFace *ref_face_ptr, DLIList<CoEdge*> &co_edge_list );
    //- Get an ordered list of coedges from the RefFace.  The RefFace must only
    //- contain one loop.

    CoEdge *get_next_co_edge( CoEdge *prev_co_edge_ptr, DLIList<CoEdge*> &co_edge_list );
    //- Given the previous CoEdge, get the next CoEdge using the CoEdges in 
    //- the given list.

    RefVertex *start_vertex( CoEdge *co_edge_ptr );
    //- Get the starting vertex of the given CoEdge.

    RefVertex *end_vertex( CoEdge *co_edge_ptr );
    //- Get the ending vertex of the given CoEdge.

    CoEdge *get_complimentary_co_edge( CoEdge *co_edge_ptr, 
                                       DLIList<CoEdge*> co_edge_list );
    //- Find the CoEdge from co_edge_list that shares a common RefEdge with
    //- the input co_edge_ptr.  Returns NULL if no complimentary CoEdge is
    //- found.  Note the input co_edge_list is copied so that it's list 
    //- position is not changed.

    CubitStatus get_attached_coedges_start( CoEdge *co_edge_ptr,
                                            DLIList<CoEdge*> &co_edge_list,
                                            DLIList<CoEdge*> &attached_co_edge_list );
    //- Get coedges from co_edge_list whose start vertices are attached to the start 
    //- vertex of co_edge_ptr

    CubitStatus remove_other_corner_coedges( CoEdge *co_edge_1, CoEdge *co_edge_2, 
                                             DLIList<CoEdge*> &all_co_edge_list,
                                             DLIList<CoEdge*> &att_co_edge_list );
    //- For isLoop situation, removes CoEdges from the "attached_co_edge_list" 
    //- that are connected to the other shared corner.  The inputs are the 
    //- previous and next CoEdges surrounding this (not the other) corner.

    CubitStatus add_tooldata( DLIList<CoEdge*> &co_edge_list, int vertex_type );
    //- Add a tooldata to each of the CoEdges in co_edge_list containing the 
    //- vertex_type at the start of the CoEdge.

    CubitBoolean is_vertex_on_side( RefVertex *ref_vertex_ptr, int side );
    //- Determines if the given vertex lies on the given side.  Note the
    //- vertex is checked spatially, not literally as if in the vertices
    //- in the outer loop.

    CubitBoolean is_vertex_in_surface( RefVertex *ref_vertex_ptr, RefFace *ref_face );
    //- Determines if the given vertex is in the given RefFace.

    CubitBoolean is_curve_in_outer_loop( RefEdge *ref_edge_ptr );
    //- Determines if the given curve lies in the outer loop of curves.

    CubitBoolean is_curve_on_side( RefEdge *ref_edge_ptr, int side );
    //- Determines if the given curve lies on the given side.

    CubitBoolean is_chain_two_triangles();
    //- Determines wheither the chain of surfaces is actually two connected
    //- triangles (see diagram in function).

    CubitBoolean is_triangle();
    //- Determines if a single surface chain is a triangle

    CubitStatus check_through_vertices( char *type );
    //- Determines if the through vertices are on valid sides of the surface.
    //- Call this function after adjusting for a Curve direction.
    //- This function only checks - it gives an error if through vertices are 
    //- on invalid sides.  Type is either "direction" or "from", for the error
    //- message ("direction" if the direction exists and was checked previously, 
    //- "from" if only a from curve was given).

    CoEdge* prev_co_edge( CoEdge *co_edge_ptr );
    //- From the outerCoEdgeLoop, get the CoEdge previous to the one given.

    CubitStatus position_co_edge_list( int i, DLIList<CoEdge*> &co_edge_list );
    //- Position the co_edge_list to the start of the loop

    CubitStatus get_a_coedges( DLIList<CoEdge*> &co_edge_list, 
                               DLIList<CoEdge*> &a_coedges );
    CubitStatus get_b_coedges( DLIList<CoEdge*> &co_edge_list, 
                               DLIList<CoEdge*> &b_coedges );
    CubitStatus get_c_coedges( DLIList<CoEdge*> &co_edge_list, 
                               DLIList<CoEdge*> &c_coedges );
    CubitStatus get_d_coedges( DLIList<CoEdge*> &co_edge_list, 
                               int num_so_far, DLIList<CoEdge*> &d_coedges );
    //- These functions find the coedges on each side of the logical rectangle,
    //- using the vertex classifications

    CubitStatus find_spline_curves( RefFace *ref_face_ptr, int num_segs, 
                                    double distance, 
                                    DLIList<Curve*> *curve_list_ptr,
                                    double tolerance,
                                    CubitBoolean parametric_flg = CUBIT_FALSE,
                                    CubitBoolean preview_flg = CUBIT_FALSE,
                                    CubitBoolean create_ref_edges_flg = CUBIT_FALSE );
    //- Uses a mapping algorithm to calculate the interior points to pass
    //- spline(s) through and creates the Curve(s) for the surface.  By default, 
    //- it will do the mapping in the 3D space of the surface.  It can be
    //- overridden to use the parametric space if a surface parametric space
    //- exists.  If preview_flg is TRUE, then the interior spline points are
    //- drawn for the user to query.  This setting, along with the create_ref_edges_flg
    //- also helps to determine whether to project the curves to the surface.
    //- In ACIS, projecting is very slow, so we skip it when previewing; however
    //- we use it if the user is ultimately creating RefEdges (note the RefEdges
    //- will not be created within this function).

    CubitStatus fill_boundary_coords( TDSplitSurface *tdss, int nr, int nc, 
                                      CubitVector ***coords );
    //- Fills the boundary array for the mapping algorithm in 3D

    CubitStatus draw_boundary_coords( TDSplitSurface *tdss );
    //- Debug function - just draws the coords for each side, if debug 100 is on

    CubitStatus smooth_interior_coords( RefFace *ref_face_ptr,
                                        TDSplitSurface *tdss,
                                        double tolerance,
                                        int nr, int nc, double distance,
                                        CubitVector ***coords );
    //- Smooths the interior coords in 3D so that the split location is correct.
    //- Respect distance if not -1.0.

    CubitStatus fill_interior_coords( TDSplitSurface *tdss, int nr, int nc,
                                      CubitVector ***coords );
    //- Fills the interior coordinates using a mapping algorithm in 3D

    CubitVector *make_interior_coord( CubitVector ***coords, int nr, int nc,
                                      double ada, double tse, int j, int i );
    //- Calculates an interior coordinate using the mapping algorithm in 3D

    CubitStatus fill_boundary_coords( TDSplitSurface *tdss, int nr, int nc, 
                                      Cubit2DPoint ***coords );
    //- Fills the boundary array for the mapping algorithm in 2D

    Cubit2DPoint *get_uv_point( RefFace *ref_face_ptr, CubitVector *vec_ptr );
    //- Given an input 3D point gets a UV point on the given surface

    CubitStatus fill_interior_coords( TDSplitSurface *tdss, int nr, int nc,
                                      Cubit2DPoint ***coords );
    //- Fills the interior coordinates using a mapping algorithm in 2D

    CubitStatus get_tse_array( TDSplitSurface *tdss, int tse_ints,
                               DLIList<double> &tse_array );
    //- For the mapping algorithm, finds an array of tse which are the
    //- average fractional distances across the AC direction of the surface

    Cubit2DPoint *make_interior_coord( Cubit2DPoint ***coords, int nr, int nc,
                                       double ada, double tse, int j, int i );
    //- Calculates an interior coordinate using the mapping algorithm in 2D

    CubitStatus draw_points( DLIList<CubitVector*> &pnt_list, int color, 
                             int flush = CUBIT_TRUE );
    CubitStatus draw_point( CubitVector &pnt, int color, int flush = CUBIT_FALSE );
    //- Draw points

    CubitStatus draw_preview( Curve *curve_ptr, CubitBoolean flush=CUBIT_TRUE,
      int color=CUBIT_BLUE );
    //- Draw the curves

    CubitStatus create_ref_edges( DLIList<Curve*> &curve_list );
    //- Create RefEdges from the given Curves

    int number_coedges( RefEdge *ref_edge_ptr, RefFace *ref_face_ptr );
    //- Finds the number of coedges for the given edge on the given surface

    void delete_surf_tooldatas( DLIList<RefFace*> &ref_face_list );
    //- Deletes tooldatas from the given ref_face_list

    void delete_coedge_tooldatas( DLIList<CoEdge*> &co_edge_list );
    //- Deletes tooldatas from the given co_edge_list

    void delete_vertex_tooldatas( DLIList<RefVertex*> &ref_vertex_list );
    //- Deletes tooldatas from the given ref_vertex_list

    void list_sides_debug();
    //- Lists the sides for debugging purposes

    Curve *create_curve( DLIList<CubitVector*> &vec_list,
                         Surface *surf_ptr,                       
                         double tolerance = GEOMETRY_RESABS,
                         CubitBoolean iterate = CUBIT_FALSE,
                         CubitBoolean draw_pnts = CUBIT_FALSE,
                         CubitBoolean project_curve = CUBIT_TRUE );
    //- Creates a curve given the input points.  The function will attempt
    //- to use a straight line or arc, if possible.  In order to create
    //- a straight line or arc, all of the input points must lie within
    //- 0.5*tolerance of a line or arc (0.5 to minimize the total deviation,
    //- exception - if tolerance==GEOMETRY_RESABS do not take 0.5), and the
    //- line or arc must lie within resabs of the surface (so that imprinting
    //- can occur).  If iterate is TRUE, if a line or arc cannot be created,
    //- the function will iterate until a spline is created that lies within 
    //- tolerance of the line segments input in vec_list (i.e., the resultant
    //- spline approximates a series of line segments).  If draw_pnts is TRUE,
    //- the function will draw the resultant spline points.  If project_curve
    //- is CUBIT_TRUE, project the curve to the surface (this option was 
    //- added to speed up previewing - for that we skip the project).

    CubitBoolean check_points_straight( Surface *surf_ptr,
                                        DLIList<CubitVector*> &point_list,
                                        double tolerance );
    //- Check to see if the points approximate a straight line, within
    //- tolerance.  Also, if surface is not a plane, all points are projected
    //- to the resultant line and checked to make sure they lie within resabs
    //- (GEOMETRY_RESABS or a user specified tolerance) to the surf_ptr so
    //- that an imprint will be successful.  If only two points are input, 
    //- additional locations along the line are checked.

    Curve* check_points_arc( Surface *surf_ptr,
                             GeometryModifyEngine *gme, 
                             DLIList<CubitVector*> &point_list,
                             double resabs,
                             double tolerance );
    //- Check to see if the input points approximate an arc.  If there are
    //- only two input points, try to create an arc if the surface is a
    //- cone, sphere or torus.

    Curve* check_if_arc( GeometryModifyEngine *gme,
                         Surface *surf_ptr,
                         DLIList<CubitVector*> &point_list,
                         CubitVector &start_pnt,
                         CubitVector &mid_pnt,
                         CubitVector &end_pnt,
                         double resabs,
                         double tolerance,
                         double &deviation );
    //- Supplemental function for check_points_arc.  This function can be
    //- called multiple times from check_points_arc as we try to iterate
    //- to a solution.

    CubitStatus get_arc_mid_pnt( Surface *surf_ptr,
                                 CubitVector &start_pnt,
                                 CubitVector &end_pnt,
                                 CubitVector &mid_pnt,
                                 double tolerance );
    //- Given the start and end points, and the surface, this function finds
    //- a third point on the arc.

    Curve* create_arc_two( GeometryModifyEngine *gme,
                           Surface *surf_ptr,
                           CubitVector &start_pnt,
                           CubitVector &end_pnt,
                           double resabs,
                           double tolerance );
    //- Given the surface and two points, attempt to find an arc that lies
    //- on the surface.

    Curve* create_arc_three( GeometryModifyEngine *gme,
                             Surface *surf_ptr,
                             CubitVector &start_pnt,
                             CubitVector &mid_pnt,
                             CubitVector &end_pnt,
                             double resabs );
    //- Given the surface and three points, attempt to find an arc that lies
    //- on the surface.

    CubitStatus reflect_arc_pnt( CubitVector &vec1,
                                 CubitVector &vec2,
                                 CubitVector &vec3,
                                 CubitVector &pnt_to_reflect,
                                 CubitVector &out_pnt );
    //- Given three points and another input point, create an arc and
    //- reflect the input point about the arc center.

    CubitBoolean is_point_on_surface( Surface *surf_ptr, CubitVector &pnt, 
                                      double resabs );
    //- Determine if the input point is within the given resabs of the surface.

    int count_surfaces_in_owning_body( RefFace *ref_face_ptr, Body *&body_ptr );
    //- Returns: Number of surfaces in the owning body of the RefFace
    //-          -1 - RefFace is free (no owning body)
    //-          -2 - RefFace is owned by more than one body
    //- Also returns the owning body via argument list

    int count_surfaces_in_body( Body *body_ptr );
    //- Return number of surfaces in the body

    int count_curves_in_body( Body *body_ptr );
    //- Return number of curves in the body

    static CubitBoolean parametricFlg; // Split parametrically or not (3D)
    static double splitTolerance; // Length tolerance for tessellations, etc.
    static CubitBoolean autoDetectTriangles; // Detect triangles automatically?
    static double sideAngleThreshold;  // Closest corner within this threshold
                                       // to 180 is removed (if pointAngleThreshold
                                       // criteria also met)
    static double pointAngleThreshold; // Corner with angle below this becomes
                                       // the triangle point (if sideAngleThreshold
                                       // criteria also met)

    DLIList<RefFace*> refFaceChain;
    DLIList<RefVertex*> throughVertexList;
    DLIList<CoEdge*> outerCoEdgeLoop;
    CoEdge *cornerCoEdge[4]; // Corner at start of given CoEdge
    CubitBoolean isLoop; // A continuous loop of surfaces?
    int sideInterval[4];
};

inline void
SplitSurfaceTool::set_tolerance( double tol )
{ splitTolerance = tol; }

inline double
SplitSurfaceTool::get_tolerance()
{ return splitTolerance; }

inline void
SplitSurfaceTool::set_auto_detect_triangles_flg( CubitBoolean flg )
{ autoDetectTriangles = flg; }

inline CubitBoolean
SplitSurfaceTool::get_auto_detect_triangles_flg()
{ return autoDetectTriangles; }

inline void
SplitSurfaceTool::set_side_angle_threshold( double tol )
{ sideAngleThreshold = tol; }

inline double
SplitSurfaceTool::get_side_angle_threshold()
{ return sideAngleThreshold; }

inline void
SplitSurfaceTool::set_point_angle_threshold( double tol )
{ pointAngleThreshold = tol; }

inline double
SplitSurfaceTool::get_point_angle_threshold()
{ return pointAngleThreshold; }

inline void
SplitSurfaceTool::set_parametric_flg( CubitBoolean parametric_flg )
{ parametricFlg = parametric_flg; }

inline CubitBoolean
SplitSurfaceTool::get_parametric_flg()
{ return parametricFlg; }

#endif

