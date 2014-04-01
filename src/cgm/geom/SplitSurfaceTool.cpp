//-------------------------------------------------------------------------
// Filename      : SplitSurfaceTool.cpp
//
// Purpose       : Split a single or chain of surfaces (ie., split a fillet 
//                 down the middle so that a mesh sweep can occur).  Used by
//                 the "split surface" commands.
//
//   Split Surface <id> Across [Pair] Location <options multiple locs>
//       [Preview [Create]]
//   Split Surface <id> Across Location <multiple locs> Onto Curve <id>
//       [Preview [Create]]
//
//   Split Surface <id_list> [Corner Vertex <id_list>] [Direction Curve <id>]
//       [Segment <val> | Fraction|Distance <val> [From Curve <id>]]
//       [Through Vertex <id_list>] [Parametric <on|OFF>] [Tolerance <value>]
//       [Preview [Create]]
//
// Special Notes : 
//
// Creator       : Steve Storm
//
// Creation Date : 10/06/2002
//-------------------------------------------------------------------------

#include "SplitSurfaceTool.hpp"
#include "AnalyticGeometryTool.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "Point.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "Loop.hpp"
#include "CubitMessage.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitUtil.hpp"
#include "DLIList.hpp"
#include "TDSplitSurface.hpp"
#include "GfxDebug.hpp"
#include "GfxPreview.hpp"
#include "Cubit2DPoint.hpp"
#include "GMem.hpp"
#include "SettingHandler.hpp"

#ifdef CUBIT_GUI
#ifndef NO_USAGE_TRACKING
#include "GUIInterface.h"
// #include "stdafx.h"
// #include "UsageTrackingTool.h"
#endif
#endif

CubitBoolean SplitSurfaceTool::parametricFlg = CUBIT_FALSE;
double SplitSurfaceTool::splitTolerance = 1.0;
CubitBoolean SplitSurfaceTool::autoDetectTriangles = CUBIT_TRUE;
double SplitSurfaceTool::sideAngleThreshold = 27.0; // From 180
double SplitSurfaceTool::pointAngleThreshold = 45.0; // Below is a point

SplitSurfaceTool::SplitSurfaceTool()
{
  int i;
  for(i=0; i<4; i++)
  {
    sideInterval[i] = 0;
    cornerCoEdge[i] = NULL;
  }
  isLoop = CUBIT_FALSE;
}

//Initialize all settings in this class
void
SplitSurfaceTool::initialize_settings()
{
  SettingHandler::instance()->add_setting("Split Surface Tolerance", 
                                          SplitSurfaceTool::set_tolerance, 
					                                SplitSurfaceTool::get_tolerance);

  SettingHandler::instance()->add_setting("Split Surface Parametric", 
                                          SplitSurfaceTool::set_parametric_flg, 
					                                SplitSurfaceTool::get_parametric_flg);

  SettingHandler::instance()->add_setting("Split Surface Auto Detect Triangles", 
                                          SplitSurfaceTool::set_auto_detect_triangles_flg, 
					                                SplitSurfaceTool::get_auto_detect_triangles_flg);

  SettingHandler::instance()->add_setting("Split Surface Side Angle Threshold",
                                          SplitSurfaceTool::set_side_angle_threshold, 
					                                SplitSurfaceTool::get_side_angle_threshold);

  SettingHandler::instance()->add_setting("Split Surface Point Angle Threshold", 
                                          SplitSurfaceTool::set_point_angle_threshold, 
					                                SplitSurfaceTool::get_point_angle_threshold);
}

CubitStatus                                
SplitSurfaceTool::preview( RefFace *ref_face_ptr,
                           DLIList<CubitVector*> &locations,
                           DLIList<DLIList<CubitVector*>*> &vec_lists,
                           CubitBoolean create_ref_edges_flg )
{
  // Create curves from the input vec_lists (locations are just the original
  // locations the user specified - these need to be drawn)
  int i;
  Curve *curve_ptr;
  DLIList<CubitVector*> *vec_list_ptr;
  vec_lists.reset();

  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();

  for( i=vec_lists.size(); i--; )
  {
    vec_list_ptr = vec_lists.get_and_step();

    vec_list_ptr->reset();
    if( vec_list_ptr->size() < 2 )
    {
      PRINT_ERROR( "Unable to create a curve from less than two locations.\n" );
      continue;
    }

    curve_ptr = create_curve( *vec_list_ptr, surf_ptr );

    if( curve_ptr )
    {
      if( create_ref_edges_flg == CUBIT_TRUE )
      {
        GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
      }
      else
      {
        draw_preview( curve_ptr, CUBIT_FALSE );
         curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
      }
    }
  }
  
  // Draw locations too
  draw_points( locations, CUBIT_BLUE );
  
  return CUBIT_SUCCESS;
}

// Create the curves used by the splitting algorithm.
// Note: calling code is responsible for cleaning up the
//       created curves when it is done with them.
CubitStatus                                
SplitSurfaceTool::calculate_split_curves( RefFace *ref_face_ptr,
                                 DLIList<CubitVector*> &locations,
                                 DLIList<DLIList<CubitVector*>*> &vec_lists,
                                 DLIList<Curve*>& curve_list )
{
  // Create curves from the input vec_lists (locations are just the original
  // locations the user specified - these need to be drawn).  NOTE: drawing
  // the locations was removed 6/16/05 per Sandia request.
  int i;
  Curve *curve_ptr;
  DLIList<CubitVector*> *vec_list_ptr;
  vec_lists.reset();

  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();

  for( i=vec_lists.size(); i--; )
  {
    vec_list_ptr = vec_lists.get_and_step();

    vec_list_ptr->reset();
    if( vec_list_ptr->size() < 2 )
    {
      PRINT_ERROR( "Unable to create a curve from less than two locations.\n" );
      continue;
    }

    curve_ptr = create_curve( *vec_list_ptr, surf_ptr );

    if( curve_ptr )
    {
      curve_list.append( curve_ptr );
      if( DEBUG_FLAG( 100 ) )
        draw_preview( curve_ptr, CUBIT_TRUE, CUBIT_RED );
    }
  }

  if( !curve_list.size() )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitStatus                                
SplitSurfaceTool::split_surface( RefFace *ref_face_ptr, DLIList<Curve*> &curve_list)
{
  // Count number of surfaces in owning body prior to split - this is used to
  // determine if split is actually successful.
  Body *body_ptr;
  int num_surfaces_prior = count_surfaces_in_owning_body( ref_face_ptr, body_ptr );
  int num_curves_prior = count_curves_in_body( body_ptr );

  if( num_surfaces_prior == -1 )
  {
    PRINT_ERROR( "Cannot split a surface that is not part of a volume\n" );
    return CUBIT_FAILURE;
  }
  else if( num_surfaces_prior == -2 )
  {
    PRINT_ERROR( "Cannot split a surface that is contained by multiple volumes\n" );
    return CUBIT_FAILURE;
  }

  // Perform the split on the real (e.g. ACIS) geometry. 
  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();

  DLIList<Surface*> surface_list;
  surface_list.append( surf_ptr );
  DLIList<DLIList<Curve*>*> curve_lists_list;
  curve_lists_list.append( &curve_list );

  Body *new_body_ptr;

  if( GeometryModifyTool::instance()->imprint( surface_list,
    curve_lists_list, new_body_ptr ) == CUBIT_FAILURE )
  {
    while( curve_list.size() ) 
      delete curve_list.pop();
    return CUBIT_FAILURE;
  }

  while( curve_list.size() ) 
    delete curve_list.pop();
  
  // TODO: Need to find a better place to do the drawing.  Should not be 
  // using debug graphics to do previews or highlights from within CGM. 
  // This could create unwanted dependencies between CGM and the graphics. 

  // Draw locations too
  // draw_points( locations, CUBIT_BLUE );
  
  int num_surfaces_after = count_surfaces_in_body( body_ptr );
  
  if( num_surfaces_after > num_surfaces_prior )
    return CUBIT_SUCCESS;
  else
  {
    int num_curves_after = count_curves_in_body( body_ptr );
    if( num_curves_after > num_curves_prior )
      return CUBIT_SUCCESS;
  }
  
  PRINT_ERROR( "Split failed - surface %d was not split\n", ref_face_ptr->id() );
  return CUBIT_FAILURE;
}

CubitStatus                                
SplitSurfaceTool::split_surface( RefFace *ref_face_ptr,
                                 DLIList<CubitVector*> &locations,
                                 DLIList<DLIList<CubitVector*>*> &vec_lists )
{
  // Count number of surfaces in owning body prior to split - this is used to
  // determine if split is actually successful.
  Body *body_ptr;
  int num_surfaces_prior = count_surfaces_in_owning_body( ref_face_ptr, body_ptr );
  int num_curves_prior = count_curves_in_body( body_ptr );
  if( num_surfaces_prior == -1 )
  {
    PRINT_ERROR( "Cannot split a surface that is not part of a volume\n" );
    return CUBIT_FAILURE;
  }
  else if( num_surfaces_prior == -2 )
  {
    PRINT_ERROR( "Cannot split a surface that is contained by multiple volumes\n" );
    return CUBIT_FAILURE;
  }

  // Find the splitting curves
  DLIList<Curve*> curve_list;
  CubitStatus err  = calculate_split_curves( ref_face_ptr, locations, 
                                             vec_lists, curve_list );
  if( err == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Perform the split on the real (e.g. ACIS) geometry. 
  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();

  DLIList<Surface*> surface_list;
  surface_list.append( surf_ptr );
  DLIList<DLIList<Curve*>*> curve_lists_list;
  curve_lists_list.append( &curve_list );

  Body *new_body_ptr;

  if( GeometryModifyTool::instance()->imprint( surface_list,
    curve_lists_list, new_body_ptr ) == CUBIT_FAILURE )
  {
    while( curve_list.size() ) 
      delete curve_list.pop();
    return CUBIT_FAILURE;
  }

  while( curve_list.size() ) 
    delete curve_list.pop();

  // TODO: Need to find a better place to do the drawing.  Should not be 
  // using debug graphics to do previews or highlights from within CGM. 
  // This could create unwanted dependencies between CGM and the graphics. 

  // Draw locations too
  //draw_points( locations, CUBIT_BLUE );

  int num_surfaces_after = count_surfaces_in_body( body_ptr );
  
  if( num_surfaces_after > num_surfaces_prior )
    return CUBIT_SUCCESS;
  else
  {
    int num_curves_after = count_curves_in_body( body_ptr );
    if( num_curves_after > num_curves_prior )
      return CUBIT_SUCCESS;
  }
  
  PRINT_ERROR( "Split failed - surface %d was not split\n", ref_face_ptr->id() );
  return CUBIT_FAILURE;
}

CubitStatus
SplitSurfaceTool::draw_points( DLIList<CubitVector*> &pnt_list, int color, 
                               int flush )
{
  int i;
  for( i=pnt_list.size(); i--; )
  {
    CubitVector *pnt_ptr = pnt_list.get_and_step();
    draw_point( *pnt_ptr, color );
  }

  if( flush )
  {
    GfxPreview::flush();
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::draw_point( CubitVector &pnt, int color, int flush )
{
  GfxPreview::draw_point( pnt, color );
  if( flush )
  {
    GfxPreview::flush();
  }
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::split_surfaces( DLIList<RefFace*> &ref_face_list,
                                 int num_segs,
                                 double fraction,
                                 double distance,
                                 RefEdge *from_curve_ptr,
                                 DLIList<RefVertex*> &corner_vertex_list,
                                 DLIList<RefVertex*> &through_vertex_list,
                                 RefEdge *curve_dir_ptr,
                                 CubitBoolean preview_flg,
                                 CubitBoolean create_ref_edges_flg )
{
  CubitStatus status;

  CubitBoolean just_curves_flg = CUBIT_FALSE;
  DLIList<DLIList<Curve*>*> curve_lists_list;

  // Call the primary function
  status = split_surfaces( ref_face_list, num_segs, fraction, distance, 
                           from_curve_ptr, corner_vertex_list, 
                           through_vertex_list, curve_dir_ptr,
                           preview_flg, create_ref_edges_flg,
                           just_curves_flg, curve_lists_list );

  free_curves_lists( curve_lists_list );

#ifdef CUBIT_GUI
#ifndef NO_USAGE_TRACKING
  if(!preview_flg)
      GUIInterface::instance()->log_split_surface_op();//UsageTrackingTool::instance()->LogSplitSurfaceOp();
#endif
#endif // CUBIT_GUI

  return status;
}

CubitStatus
SplitSurfaceTool::split_surfaces( DLIList<RefFace*> &ref_face_list,
                                  int num_segs,
                                  double fraction,
                                  double distance,
                                  RefEdge *from_curve_ptr,
                                  DLIList<RefVertex*> &corner_vertex_list,
                                  DLIList<RefVertex*> &through_vertex_list,
                                  RefEdge *curve_dir_ptr,
                                  DLIList<DLIList<Curve*>*> &curve_lists_list )
{
  CubitStatus status;
  CubitBoolean preview_flg = CUBIT_FALSE;
  CubitBoolean create_ref_edges_flg = CUBIT_FALSE;
  CubitBoolean just_curves_flg = CUBIT_TRUE;

  // Call the primary function
  status = split_surfaces( ref_face_list, num_segs, fraction, distance, 
                           from_curve_ptr, corner_vertex_list, 
                           through_vertex_list, curve_dir_ptr,
                           preview_flg, create_ref_edges_flg,
                           just_curves_flg, curve_lists_list );

  return status;
}

CubitStatus
SplitSurfaceTool::split_surfaces( DLIList<RefFace*> &ref_face_list,
                                  int num_segs,
                                  double fraction,
                                  double distance,
                                  RefEdge *from_curve_ptr,
                                  DLIList<RefVertex*> &corner_vertex_list,
                                  DLIList<RefVertex*> &through_vertex_list,
                                  RefEdge *curve_dir_ptr,
                                  DLIList<Curve*> &curve_list )
{
  CubitStatus status;
  CubitBoolean preview_flg = CUBIT_FALSE;
  CubitBoolean create_ref_edges_flg = CUBIT_FALSE;
  CubitBoolean just_curves_flg = CUBIT_TRUE;
  DLIList<DLIList<Curve*>*> curve_lists_list;

  // Call the primary function
  status = split_surfaces( ref_face_list, num_segs, fraction, distance, 
                           from_curve_ptr, corner_vertex_list, 
                           through_vertex_list, curve_dir_ptr,
                           preview_flg, create_ref_edges_flg,
                           just_curves_flg, curve_lists_list );

  // Convert the curves to the simple curve_list
  int i, j;
  curve_lists_list.reset();
  for( i=curve_lists_list.size(); i--; )
  {
    DLIList<Curve*> *curve_list_ptr = curve_lists_list.get_and_step();

    curve_list_ptr->reset();
    for( j=curve_list_ptr->size(); j--; )
    {
      Curve *curve_ptr = curve_list_ptr->get_and_step();
      curve_list.append( curve_ptr );
    }
  }

  free_curves_lists( curve_lists_list, CUBIT_FALSE );

  return status;
}

CubitStatus
SplitSurfaceTool::calculate_split_curves( DLIList<RefFace*> &ref_face_list, 
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
                                          DLIList<DLIList<Curve*>*> &curve_lists_list )
{
  int i;
  refFaceChain = ref_face_list;
  throughVertexList = through_vertex_list;

  // Get tolerance and parametric_flg from member variables
  double tolerance = splitTolerance;
  CubitBoolean parametric_flg = parametricFlg;

  // Check for valid fraction
  if( fraction<0.0 || fraction>1.0 )
  {
    PRINT_ERROR( "Fraction must be between 0.0 and 1.0 - aborting.\n" );
    return CUBIT_FAILURE;
  }

  // Check for valid number of segments
  if( num_segs < 2 )
  {
    PRINT_ERROR( "Number of specified segments must be >= 2\n" );
    return CUBIT_FAILURE;
  }

  // Not valid to specify through vertices with multiple segments
  if( num_segs > 2 && throughVertexList.size() )
  {
    PRINT_ERROR( "Through vertices specified - not valid if number of segments > 2\n" );
    return CUBIT_FAILURE;
  }

  // Check number of through vertices
  if( throughVertexList.size() > refFaceChain.size()+1 )
  {
    PRINT_ERROR( "Too many 'through' vertices specified.\n"
      "       Can only be equal to one more than the number of surfaces to split.\n" );
    return CUBIT_FAILURE;
  }

  // If a distance was specified, set the parametric flag to FALSE
  if( distance != -1.0 )
    parametric_flg = CUBIT_FALSE;

  // Check the individual surfaces for errors (multiple loops, hardlines)
  if( check_valid_faces() == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Order the face list from one end of the chain to the other.  This function
  // also checks for the isLoop condition.
  if( order_face_list() == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  if( DEBUG_FLAG(154) )
  {
    DLIList<CubitEntity*> cubit_faces;
    CAST_LIST(refFaceChain, cubit_faces, CubitEntity);
    CubitUtil::list_entity_ids( "\nOrdered surface list: ", 
      cubit_faces, 80, "\n", CUBIT_FALSE );
  }

  // Get all the outer loops, in the proper order.  This also makes sure 
  // they start at the beginning of the first surface in the chain.
  if( get_outer_loops() == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Make sure direction and from curves are valid, if specified
  if( curve_dir_ptr && !is_curve_in_outer_loop( curve_dir_ptr ) )
  {
    PRINT_ERROR( "Specified direction curve %d not found on side of logical rectangle\n",
      curve_dir_ptr->id() );
    return CUBIT_FAILURE;
  }
  if( from_curve_ptr && !is_curve_in_outer_loop( from_curve_ptr ) )
  {
    PRINT_ERROR( "Specified from curve %d not found on side of logical rectangle\n",
      from_curve_ptr->id() );
    return CUBIT_FAILURE;
  }

  // Find the corners
  if( corner_vertex_list.size() )
  {
    if( corner_vertex_list.size() != 4 )
    {
      PRINT_ERROR( "You must select exactly 4 corner vertices\n" );
      return CUBIT_FAILURE;
    }

    if( isLoop )
    {
      PRINT_WARNING( "Ignoring specified corners since continuous loop situation\n" );
    }
    else
    {
      // Now order these vertices properly w/respect to the outerCoEdgeLoop
      if( order_selected_corners( corner_vertex_list ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;
    }
  }
  else
  {
    // Pick the corners based on angle criteria
    if( refFaceChain.size() > 1 )
    {
      if( pick_4_corners() == CUBIT_FAILURE )
        return CUBIT_FAILURE;
    }
    else
    {
      // Single surface case
      if( pick_4_corners_simple() == CUBIT_FAILURE )
        return CUBIT_FAILURE;
    }
  }

  // Adjust the surface to the split direction.  This function also adjusts
  // and checks for errors if a throughVertexList was specified.
  if( adjust_for_split_direction( curve_dir_ptr, from_curve_ptr, corner_vertex_list )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // If previewing, draw the corners in RED
  if( preview_flg == CUBIT_TRUE )
  {
    CubitVector temp_vec = start_vertex(cornerCoEdge[0])->coordinates();
  
    draw_point( temp_vec, CUBIT_RED );
    temp_vec = start_vertex(cornerCoEdge[1])->coordinates();
    draw_point( temp_vec, CUBIT_RED );
    temp_vec = start_vertex(cornerCoEdge[2])->coordinates();
    draw_point( temp_vec, CUBIT_RED );
    temp_vec = start_vertex(cornerCoEdge[3])->coordinates();
    draw_point( temp_vec, CUBIT_RED, CUBIT_TRUE );
  }

  if( DEBUG_FLAG(154) )
  {
    DLIList<RefEdge*> curve_list;
    get_outer_curves( curve_list );
    DLIList<CubitEntity*> cubit_edges;
    CAST_LIST(curve_list, cubit_edges, CubitEntity);
    CubitUtil::list_entity_ids( "Ordered curves outer loop: ", 
      cubit_edges, 80, "\n", CUBIT_FALSE );
    
    PRINT_INFO( "Best corner 1 = Vertex %d\n", start_vertex(cornerCoEdge[0])->id() );
    PRINT_INFO( "Best corner 2 = Vertex %d\n", start_vertex(cornerCoEdge[1])->id() );
    PRINT_INFO( "Best corner 3 = Vertex %d\n", start_vertex(cornerCoEdge[2])->id() );
    PRINT_INFO( "Best corner 4 = Vertex %d\n", start_vertex(cornerCoEdge[3])->id() );

    PRINT_INFO( "Side interval #1 = %d\n", sideInterval[0] );
    PRINT_INFO( "Side interval #2 = %d\n", sideInterval[1] );
    PRINT_INFO( "Side interval #3 = %d\n", sideInterval[2] );
    PRINT_INFO( "Side interval #4 = %d\n", sideInterval[3] );
  }

  // Add a tooldata to each surface for later use
  RefFace *ref_face_ptr;
  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    ref_face_ptr->add_TD( new TDSplitSurface( ref_face_ptr ) );
  }

  // Add curve loops to each tooldata.  These curve loops will form the
  // boundary of a mapped mesh algorithm that will find the interior
  // points for the split locations.
  if( populate_curve_loops() == CUBIT_FAILURE )
  {
    delete_surf_tooldatas( refFaceChain );
    return CUBIT_FAILURE;
  }

  // Now build up the coordinates and get the curve(s) for each surface
  // curve_lists_list will contain a list of curve lists for each surface

  TDSplitSurface *tdss;
  refFaceChain.reset();
  for(i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    tdss = (TDSplitSurface *)ref_face_ptr->
      get_TD(&TDSplitSurface::is_split_surface);

    // Setup matched curve tessellations on all curves (the tessellations
    // are "matched" from side to side - the points on each side are
    // at identical parameter values along the side).  This is not a 
    // perfect solution but seems to work okay for the typical geometry
    // types this tool is used on.
    if( tdss->tessellate_sides( tolerance, fraction, distance, num_segs,
        throughVertexList ) == CUBIT_FAILURE )
    {
      delete_vertex_tooldatas( through_vertex_list );
      delete_surf_tooldatas( refFaceChain );
      return CUBIT_FAILURE;
    }

    // curve_list_ptr will contain a list of curves for the individual surface
    DLIList<Curve*> *curve_list_ptr = new DLIList<Curve*>;
    curve_lists_list.append( curve_list_ptr );

    // Use a mapping concept to get the spline coordinates of the curve to
    // split with.  The function find_spline_curves will populate curve_list_ptr
    // with potentially multiple curves (typically there will only be one curve 
    // though).
    if( find_spline_curves( ref_face_ptr, num_segs, distance, curve_list_ptr,
                            tolerance, parametric_flg, preview_flg,
                            create_ref_edges_flg ) == CUBIT_FAILURE )
    {
      delete_vertex_tooldatas( throughVertexList );
      delete_surf_tooldatas( refFaceChain );
      return CUBIT_FAILURE;
    }

    if( just_curves_flg == CUBIT_TRUE )
      continue;

    if( preview_flg == CUBIT_TRUE )
    {
      if( create_ref_edges_flg == CUBIT_FALSE )
        draw_preview( *curve_list_ptr );
      else
      {
        create_ref_edges( *curve_list_ptr );

        // This just draws each curve as we go, same as preview does
        GfxPreview::flush();
      }
    }
  }

  // Determine if all of the 'through' vertices were used - if not give a warning
  if( throughVertexList.size() )
  {
    DLIList<RefVertex*> non_used_through_vertex_list;
    RefVertex *ref_vertex_ptr;
    throughVertexList.reset();
    for( i=throughVertexList.size(); i--; )
    {
      ref_vertex_ptr = throughVertexList.get_and_step();
      tdss = (TDSplitSurface *)ref_vertex_ptr->
        get_TD(&TDSplitSurface::is_split_surface);
      if( !tdss )
        non_used_through_vertex_list.append( ref_vertex_ptr );
    }
    
    if( non_used_through_vertex_list.size() )
    {
      DLIList<CubitEntity*> cubit_verts;
      CAST_LIST(non_used_through_vertex_list, cubit_verts, CubitEntity);
      CubitUtil::list_entity_ids( "WARNING - unused 'through' vertices: ", 
        cubit_verts, 80, "\n", CUBIT_FALSE );
      PRINT_INFO("         They were not found on the split path.\n" );
    }
  }

  // Remove tooldatas
  delete_vertex_tooldatas( throughVertexList );
  delete_surf_tooldatas( refFaceChain );

  return CUBIT_SUCCESS;
}

//- Begin Private Functions
CubitStatus
SplitSurfaceTool::split_surfaces( DLIList<RefFace*> &ref_face_list, 
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
                                  DLIList<DLIList<Curve*>*> &curve_lists_list )
{

  // Find the splitting curves
  if( calculate_split_curves( ref_face_list, num_segs, fraction, distance,
                              from_curve_ptr, corner_vertex_list, through_vertex_list,
                              curve_dir_ptr, preview_flg, create_ref_edges_flg,
                              just_curves_flg, curve_lists_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;                         
  
  // Do the splitting
  if( preview_flg==CUBIT_FALSE && just_curves_flg==CUBIT_FALSE )
  {
    int i;
    DLIList<Surface*> surface_list;
    refFaceChain.reset();
    for( i=refFaceChain.size(); i--; )
      surface_list.append( refFaceChain.get_and_step()->get_surface_ptr() );

    Body *new_body_ptr;
    if( GeometryModifyTool::instance()->imprint( surface_list,
      curve_lists_list, new_body_ptr ) == CUBIT_FAILURE )
    {
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

void
SplitSurfaceTool::free_curves_lists( DLIList<DLIList<Curve*>*> &curve_lists_list,
                                     CubitBoolean free_curves_flg)
{
  while( curve_lists_list.size() )
  {
    DLIList<Curve*> *curve_list_ptr = curve_lists_list.pop();
    
    if( free_curves_flg == CUBIT_TRUE )
    {
      while( curve_list_ptr->size() )
      {
        Curve *curve_ptr = curve_list_ptr->pop();
        
        // If there is no RefEdge attached to this Curve, delete it
        RefEdge* ref_edge_ptr = dynamic_cast<RefEdge*>(curve_ptr->topology_entity());
        if( !ref_edge_ptr )
           curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
      }
    }

    delete curve_list_ptr;
  }
}

void
SplitSurfaceTool::delete_surf_tooldatas( DLIList<RefFace*> &ref_face_list )
{
  int i;
  for( i=ref_face_list.size(); i--; )
  {
    RefFace *ref_face_ptr = ref_face_list.get_and_step();
    ref_face_ptr->delete_TD( &TDSplitSurface::is_split_surface );
  }
}

void
SplitSurfaceTool::delete_coedge_tooldatas( DLIList<CoEdge*> &co_edge_list )
{
  int i;
  for( i=co_edge_list.size(); i--; )
  {
    CoEdge *co_edge_ptr = co_edge_list.get_and_step();
    co_edge_ptr->delete_TD( &TDSplitSurface::is_split_surface );
  }
}

void
SplitSurfaceTool::delete_vertex_tooldatas( DLIList<RefVertex*> &ref_vertex_list )
{
  int i;
  for( i=ref_vertex_list.size(); i--; )
  {
    RefVertex *ref_vertex_ptr = ref_vertex_list.get_and_step();
    ref_vertex_ptr->delete_TD( &TDSplitSurface::is_split_surface );
  }
}

CubitStatus
SplitSurfaceTool::check_valid_faces()
{
  // Make sure each face doesn't have multiple loops or hardlines
  int i;
  RefFace *ref_face_ptr;
  CubitBoolean is_loop;
  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    if( check_face_loop( ref_face_ptr, is_loop ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    if( refFaceChain.size()>1 && is_loop==CUBIT_TRUE )
    {
      PRINT_ERROR( "Surface %d loops back on itself - cannot split chain.\n",
        ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::order_face_list()
{
  int i;
  RefFace *ref_face_ptr = NULL;

  if( refFaceChain.size() > 1 )
  {
    // First make sure surfaces are at least connected
    refFaceChain.reset();
    for( i=refFaceChain.size(); i--; )
    {
      ref_face_ptr = refFaceChain.get_and_step();

      DLIList<RefFace*> neighbor_ref_faces;
      get_neighbors( ref_face_ptr, refFaceChain, neighbor_ref_faces );
      if( neighbor_ref_faces.size() == 0 )
      {
        PRINT_ERROR( "You must select a continuous chain of surfaces.\n"
          "       Surface %d is not attached to the other surfaces.\n",
          ref_face_ptr->id() );
        return CUBIT_FAILURE;
      }
    } 
  }

  if( refFaceChain.size() < 3 )
  {
    if( refFaceChain.size() == 1 )
    {
      // Check for isLoop situation.  Assume we have a loop if some of
      // the RefEdges have 2 CoEdges that don't immediately loop back
      // on themselves (this is a hardline).
      ref_face_ptr = refFaceChain.get();
      if( check_face_loop( ref_face_ptr, isLoop ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;

      return CUBIT_SUCCESS;
    }

    if( refFaceChain.size() == 2 )
    {
      // Check for isLoop situation.  Assume we have a loop if two separate
      // chains of curves from the first surface are shared by the second
      // surface.
      refFaceChain.reset();
      RefFace *ref_face_ptr1 = refFaceChain.get_and_step();
      RefFace *ref_face_ptr2 = refFaceChain.get();
      CoEdge *start_co_edge_ptr = NULL; // Dummy
      if( check_for_loop( ref_face_ptr1, ref_face_ptr2, isLoop, start_co_edge_ptr )
        == CUBIT_FAILURE )
        return CUBIT_FAILURE;
    }
    return CUBIT_SUCCESS;
  }

  // Get the face list going from one end to the other
  DLIList<RefFace*> ordered_face_list;

  // Find an end - would have only one attached face
  int found = 0;
  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();

    DLIList<RefFace*> neighbor_ref_faces;
    get_neighbors( ref_face_ptr, refFaceChain, neighbor_ref_faces );
    if( neighbor_ref_faces.size() == 1 )
    {
      found = 1;
      break;
    }
    else if( neighbor_ref_faces.size() == 0 )
    {
      PRINT_ERROR( "You must select a continuous chain of surfaces.\n"
        "       Surface %d is not attached to the other surfaces.\n",
        ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }

  // Check for continuous loop of surfaces
  if( !found )
  {
    // This is a continuous loop of surfaces - just use the first one the 
    // user picked as starting surface of the patch.
    isLoop = CUBIT_TRUE;
    refFaceChain.reset();
    ref_face_ptr = refFaceChain.get();
  }

  // Walk across the surfaces
  ordered_face_list.append( ref_face_ptr );
  DLIList<RefFace*> remaining_face_list = refFaceChain;
  remaining_face_list.remove( ref_face_ptr );
  remaining_face_list.reset();
  int num_faces = refFaceChain.size()-1;

  for( i=num_faces; i--; )
  {
    DLIList<RefFace*> neighbor_ref_faces;
    get_neighbors( ref_face_ptr, remaining_face_list, neighbor_ref_faces );
    if( (!isLoop && neighbor_ref_faces.size() != 1) ||
        (i==num_faces-1 && isLoop && neighbor_ref_faces.size() != 2) ||
        (i!=num_faces-1 && isLoop && neighbor_ref_faces.size() != 1) )
    {
      PRINT_ERROR( "Selected surfaces do not appear to form logical rectangle\n" );
      return CUBIT_FAILURE;
    }
    neighbor_ref_faces.reset();
    ref_face_ptr = neighbor_ref_faces.get();
    ordered_face_list.append( ref_face_ptr );
    remaining_face_list.remove( ref_face_ptr );
  }

  refFaceChain.clean_out();
  ordered_face_list.reset();
  refFaceChain = ordered_face_list;

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::check_face_loop( RefFace *ref_face_ptr, CubitBoolean &is_loop )
{
  is_loop = CUBIT_FALSE;
  
  // Get list of coedges for this surface
  DLIList<CoEdge*> co_edge_list;
  if( ordered_co_edges( ref_face_ptr, co_edge_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  int i;
  co_edge_list.reset();
  CoEdge *co_edge_ptr;
  DLIList<RefEdge*> ref_edge_list;
  RefEdge *ref_edge_ptr;
  int prev_i = -99;
  RefEdge *prev_ref_edge_ptr = NULL;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();

    if( ref_edge_list.is_in_list( ref_edge_ptr ) )
    {
      // This may be a loop.
      is_loop = CUBIT_TRUE;

      // Make sure it is not a hardline.  Check for situation where loop
      // turns back on itself.

      if( ref_edge_ptr==prev_ref_edge_ptr && prev_i==i+1 )
      {
        PRINT_ERROR( "Detected hardline in Surface %d - cannot split.\n", 
          ref_face_ptr->id() );
        return CUBIT_FAILURE;
      }
    }

    prev_i = i;
    prev_ref_edge_ptr = ref_edge_ptr;

    ref_edge_list.append( ref_edge_ptr );
  }

  return CUBIT_SUCCESS;
}

CubitStatus 
SplitSurfaceTool::check_for_loop( RefFace *ref_face_ptr1, RefFace *ref_face_ptr2, 
                                  CubitBoolean &is_loop, CoEdge *&start_co_edge_ptr )
{
  // This algorithm is not perfect, but should cover most cases.
  RefEdge *ref_edge_ptr;

  is_loop = CUBIT_FALSE;
  
  // Get list of coedges for the first surface
  DLIList<CoEdge*> co_edge_list1;
  if( ordered_co_edges( ref_face_ptr1, co_edge_list1 ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Get list of refedges for the second surface
  DLIList<RefEdge*> ref_edge_list2;
  DLIList<Loop*> loop_list;
  ref_face_ptr2->ordered_loops( loop_list );
  if( loop_list.size() > 1 )
  {
    PRINT_ERROR( "Only surfaces with a single loop are allowed.\n"
      "       Surface %d has %d loops.\n", ref_face_ptr2->id(), 
      loop_list.size() );
    return CUBIT_FAILURE;
  }
  loop_list.get()->ordered_ref_edges( ref_edge_list2 );

  // Walk around the first surface

  // First, position the list at a starting coedge - one that is shared with
  // the other surface.
  int i, j;
  co_edge_list1.reset();
  for( i=0; i<co_edge_list1.size(); i++ )
  {
    ref_edge_ptr = co_edge_list1.get()->get_ref_edge_ptr();
    
    if( ref_edge_list2.is_in_list( ref_edge_ptr ) )
    {
      if( i==0 ) // Backup until we are at the start
      {
        for( j=co_edge_list1.size(); j--; )
        {
          co_edge_list1.back();
          ref_edge_ptr = co_edge_list1.get()->get_ref_edge_ptr();
          if( ref_edge_list2.is_in_list( ref_edge_ptr ) )
            continue;
          co_edge_list1.step();
        }
      }
      break;
    }
    co_edge_list1.step();
  }
  
  // Now the list is positioned at the start of one of the shared chains of 
  // curves.
  start_co_edge_ptr = co_edge_list1.get();

  // Walk along the first surface
  int found_break = 0;
  for( i=co_edge_list1.size(); i--; )
  {
    ref_edge_ptr = co_edge_list1.get_and_step()->get_ref_edge_ptr();

    if( found_break && ref_edge_list2.is_in_list( ref_edge_ptr ) )
    {
      is_loop = CUBIT_TRUE;
      return CUBIT_SUCCESS;
    }

    if( !ref_edge_list2.is_in_list( ref_edge_ptr ) )
      found_break = 1;
  }

  return CUBIT_SUCCESS;
}

CubitStatus 
SplitSurfaceTool::get_neighbors( RefFace *seed_ref_face,
                                 DLIList<RefFace*> &input_ref_faces,
                                 DLIList<RefFace*> &neighbor_ref_faces )
{
  int i, j;

  // Copy the input_ref_faces_copy list so as to not change it
  DLIList<RefFace*> input_ref_faces_copy(input_ref_faces.size());
  input_ref_faces_copy = input_ref_faces;

  // Get the edges
  DLIList<RefEdge*> ref_edge_list;
  seed_ref_face->ref_edges( ref_edge_list );

  // Get all ref_faces attached to these edges
  RefEdge *ref_edge_ptr;
  RefFace *ref_face_ptr;
  for( i=ref_edge_list.size(); i--; )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();

    DLIList<RefFace*> attached_ref_faces;
    ref_edge_ptr->ref_faces( attached_ref_faces );

    for( j=attached_ref_faces.size(); j--; )
    {
      ref_face_ptr = attached_ref_faces.get_and_step();

      if( ref_face_ptr == seed_ref_face )
        continue;

      // Don't consider ref_faces that aren't in the input list
      if( input_ref_faces_copy.is_in_list( ref_face_ptr ) )
      {
        neighbor_ref_faces.append_unique( ref_face_ptr );
      }
    }
  }

  // Respect the incoming order of ref_faces
  if( neighbor_ref_faces.size() == 2 )
  {
    input_ref_faces_copy.reset();
    neighbor_ref_faces.reset();
    RefFace *first_neighbor = neighbor_ref_faces.get_and_step();
    RefFace *second_neighbor = neighbor_ref_faces.get_and_step();
    for( i=input_ref_faces_copy.size(); i--; )
    {
      ref_face_ptr = input_ref_faces_copy.get_and_step();
      if( ref_face_ptr == first_neighbor )
        break;
      else if( ref_face_ptr == second_neighbor )
      {
        // Reverse the order
        neighbor_ref_faces.reverse();
        break;
      }
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_outer_loops()
{
  // Get the outer curve loop
  if( get_outer_coedge_loop() == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // If we found a loop, we already setup the rest of the datastructures
  if( isLoop )
    return CUBIT_SUCCESS;

  // Make sure the loop has at least 3 vertices, to form at least a 
  // degenerate rectangle
  if( outerCoEdgeLoop.size() < 3 )
  {
    PRINT_ERROR( "Selected surfaces do not appear to form logical rectangle\n" );
    return CUBIT_FAILURE;
  }

  // Orient the loops so they start at the first CoEdge of the 
  // starting RefFace in the chain
  if( refFaceChain.size() > 1 )
  {
    refFaceChain.reset();
    RefFace *first_ref_face = refFaceChain.get();

    // Move forward in the vertex list until we are on a surface
    // other than the first one
    int i;
    CoEdge *co_edge_ptr = NULL;
    outerCoEdgeLoop.reset();
    for( i=outerCoEdgeLoop.size(); i--; )
    {
      co_edge_ptr = outerCoEdgeLoop.get_and_step();
      if( co_edge_ptr->get_ref_face() != first_ref_face )
        break;
    }

    // Now move forward in the list until we are on the first surface
    for( i=outerCoEdgeLoop.size(); i--; )
    {
      co_edge_ptr = outerCoEdgeLoop.get();
      if( co_edge_ptr->get_ref_face() == first_ref_face )
        break;
      outerCoEdgeLoop.step();
    }

    PRINT_DEBUG_154( "First vertex on first surface = %d\n", start_vertex( co_edge_ptr )->id() );
    PRINT_DEBUG_154( "Index of vertex list = %d\n", outerCoEdgeLoop.get_index() );

    // Reorient the lists appropriately
    reorient_loop( outerCoEdgeLoop.get_index() );

    if( DEBUG_FLAG(154) )
    {
      DLIList<RefVertex*> ref_vertex_list;
      get_outer_vertices( ref_vertex_list );
      DLIList<CubitEntity*> cubit_verts;
      CAST_LIST(ref_vertex_list, cubit_verts, CubitEntity);
      CubitUtil::list_entity_ids( "New vertex list: ", 
        cubit_verts, 80, "\n", CUBIT_FALSE );
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_outer_coedge_loop()
{

  // Walk around the loop

  // Start by gathering a list of all coedges
  int i;
  DLIList<CoEdge*> co_edge_list;

  RefFace *ref_face_ptr;
  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    if( ordered_co_edges( ref_face_ptr, co_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

//  PRINT_INFO( "Entire coedge loop, starting at beginning -\n" );
//  PRINT_INFO( "-------------------------------------------\n" );
//  co_edge_list.reset();
//  for( i=co_edge_list.size(); i--; )
//  {
//    CoEdge *co_edge_ptr = co_edge_list.get_and_step();
//    PRINT_INFO( " Curve %d on Surface %d\n",
//      co_edge_ptr->get_ref_edge_ptr()->id(), 
//      co_edge_ptr->get_ref_face()->id() );
//  }

  CoEdge *start_co_edge_ptr;
  if( find_loop_start( start_co_edge_ptr ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to find loop start\n" ); // TODO: better message
    return CUBIT_FAILURE;
  }

//  PRINT_INFO( "Loop start = Curve %d on Surface %d\n", 
//  start_co_edge_ptr->get_ref_edge_ptr()->id(),
//  start_co_edge_ptr->get_ref_face()->id() );

  // Remove CoEdges on curves shared between surfaces (except for the starting
  // chain for isLoop situation).  For isLoop, lets create a list of those to
  // be added back in later (keep_co_edge_list).
  DLIList<CoEdge*> keep_co_edge_list;
  CoEdge *shared_co_edge_ptr;
  CoEdge *co_edge_ptr;  
  if( isLoop )
  {
    // Method - traverse co_edge_list, adding coedges to keep_co_edge_list 
    // until we don't find a coedge that shares a RefEdge with another coedge.  
    // We use the term "complimentary" CoEdge to refer to the other CoEdge that
    // shares the common RefEdge.
    //    
    // To illustrate, in the diagram below, assume S1 is the first surface and
    // S2 is the last surface in the surface loop we are splitting (note that
    // S1 can equal S2).  C1, C2 and C3 make up the "chain" of curvs, or "seam"
    // that is between S1 and S2.  Note that C1 has 2 CoEdges, one owned by S1 
    // and one owned by S2.  So, starting at C1 on S1, we step through the list
    // until we encounter a CoEdge that does not have a "complimentary" CoEdge.  
    // Just keep track of the CoEdges as we go.
    //
    //     ------------+------------
    //    .            |            .
    //    .            |C1          .
    //    .            |            .
    //    .        ^   +   |        .
    //    .        |   |   |        .
    //    .   S2   |   |C2 |   S1   .
    //    .        |   |   v        .
    //    .            +            .
    //    .            |            .
    //    .            |C3          .
    //    .            |            .
    //     ------------+------------      

    // For this code, it is important that all the ordered coedges from the 
    // starting surface are at the beginning of the co_edge_list (otherwise, we 
    // will traverse partway up the first chain and then back).  So we rebuild 
    // co_edge_list - putting all the ordered coedges from the first surface at
    // the beginning of co_edge_list. 
    DLIList<CoEdge*> temp_co_edge_list = co_edge_list;
    temp_co_edge_list.move_to( start_co_edge_ptr );
    co_edge_list.clean_out();

    RefFace *start_ref_face_ptr = start_co_edge_ptr->get_ref_face();

    for( i=temp_co_edge_list.size(); i--; )
    {
      co_edge_ptr = temp_co_edge_list.get();
      if( co_edge_ptr->get_ref_face() == start_ref_face_ptr )
      {
        co_edge_list.append( co_edge_ptr );
        temp_co_edge_list.change_to( NULL );
      }
      temp_co_edge_list.step();
    }
    temp_co_edge_list.remove_all_with_value( NULL );

    co_edge_list += temp_co_edge_list;

//    PRINT_INFO( "Entire coedge loop -\n" );
//    PRINT_INFO( "--------------------\n" );
//    co_edge_list.reset();
//    for( i=co_edge_list.size(); i--; )
//    {
//      CoEdge *co_edge_ptr = co_edge_list.get_and_step();
//      PRINT_INFO( " Curve %d on Surface %d\n",
//        co_edge_ptr->get_ref_edge_ptr()->id(), 
//        co_edge_ptr->get_ref_face()->id() );
//    }

//    PRINT_INFO( "Coedges to keep - \n" );
    co_edge_list.reset();
    co_edge_ptr = co_edge_list.get();
    RefFace *ref_ref_face_ptr = co_edge_ptr->get_ref_face();
    RefFace *shared_ref_face_ptr;
    shared_co_edge_ptr = get_complimentary_co_edge( co_edge_ptr, co_edge_list );
    RefFace *ref_shared_ref_face_ptr = shared_co_edge_ptr->get_ref_face();
    for( i=co_edge_list.size(); i--; )
    {
      // Keep track of the starting and ending chain which we will lose in the
      // code below that removes all coedges that share a refedge. To do this, 
      // add coedges until we don't find a shared one, or the reference surfaces
      // change (ref_ref_face_ptr || ref_shared_ref_face_ptr).  Why look at
      // the reference surfaces?  Because there could be a triangle as in the
      // diagram below.
      // 
      //     ------------+-------------
      //    .            |             .
      //    .            |C1           .
      //    .            |             .
      //    .        ^   +   |    S1   .
      //    .        |   |   |         .
      //    .   S2   |   |C2 |       / .
      //    .        |   |   v     /   .
      //    .            +       /     .
      //    .            |     /       .
      //    .            |C3 /     S3  .
      //    .            | /           .
      //     ------------+-------------      
      co_edge_ptr = co_edge_list.get_and_step();
      ref_face_ptr = co_edge_ptr->get_ref_face();
      shared_co_edge_ptr = get_complimentary_co_edge( co_edge_ptr, co_edge_list );
      if( !shared_co_edge_ptr )
        break;
      shared_ref_face_ptr = shared_co_edge_ptr->get_ref_face();
      if( ref_face_ptr==ref_ref_face_ptr && 
          shared_ref_face_ptr==ref_shared_ref_face_ptr  )
      {
        keep_co_edge_list.append( co_edge_ptr );
//        PRINT_INFO( " Keeping Curve %d on Surface %d\n", 
//          co_edge_ptr->get_ref_edge_ptr()->id(), 
//          co_edge_ptr->get_ref_face()->id() );
        keep_co_edge_list.append( shared_co_edge_ptr );
//        PRINT_INFO( " Keeping Curve %d on Surface %d\n", 
//          shared_co_edge_ptr->get_ref_edge_ptr()->id(), 
//          shared_co_edge_ptr->get_ref_face()->id() );
      }
      else
        break;
    }
  }

  // Remove all coedges that share a RefEdge
  for( i=co_edge_list.size(); i--; )
  {
    // Add coedges until we don't find a shared one - this adds just the
    // first chain.
    co_edge_ptr = co_edge_list.get();
    if( co_edge_ptr == NULL )
    {
      co_edge_list.step();
      continue;
    }
    shared_co_edge_ptr = get_complimentary_co_edge( co_edge_ptr, co_edge_list );
    if( shared_co_edge_ptr )
    {
      // This is not very efficient, but okay since we have relatively short
      // lists.
      co_edge_list.move_to( shared_co_edge_ptr );
//      PRINT_INFO( " Removing Curve %d on Surface %d\n",
//        shared_co_edge_ptr->get_ref_edge_ptr()->id(), 
//          shared_co_edge_ptr->get_ref_face()->id() );

      co_edge_list.change_to( NULL );
      co_edge_list.move_to( co_edge_ptr );
//      PRINT_INFO( " Removing Curve %d on Surface %d\n",
//        co_edge_ptr->get_ref_edge_ptr()->id(), 
//          co_edge_ptr->get_ref_face()->id() );

      co_edge_list.change_to( NULL );
    }
    co_edge_list.step();
  }

  co_edge_list.remove_all_with_value( NULL );

  // Add back in the removed coedges if isLoop
  if( isLoop )
  {
    // Put them at the beginning for a slight efficiency gain
    DLIList<CoEdge*> temp_co_edge_list = co_edge_list;
    co_edge_list.clean_out();
    co_edge_list += keep_co_edge_list;
    co_edge_list += temp_co_edge_list;
  }

  // Walk around the loop until we encounter the start coedge again.
  // Note loop is positioned at starting coedge.
  CoEdge *prev_co_edge_ptr = start_co_edge_ptr;
  outerCoEdgeLoop.append( start_co_edge_ptr );
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = get_next_co_edge( prev_co_edge_ptr, co_edge_list );

    if( co_edge_ptr == NULL )
    {
      PRINT_ERROR( "Selected surfaces do not appear to form logical rectangle\n" );
      return CUBIT_FAILURE;
    }

    if( co_edge_ptr == start_co_edge_ptr ) // We're done
      break;

    outerCoEdgeLoop.append( co_edge_ptr );
    prev_co_edge_ptr = co_edge_ptr;
    co_edge_list.remove( co_edge_ptr );
  }

  // Set the corners and such if isLoop
  if( isLoop )
  {
    // Start of first curve is first corner
    int best_corner_1, best_corner_2, best_corner_3=-1, best_corner_4=-1;

    best_corner_1 = 0;
    best_corner_2 = keep_co_edge_list.size()/2;
    outerCoEdgeLoop.reset();
    for( i=0; i<outerCoEdgeLoop.size()+1; i++ )
    {
      co_edge_ptr = outerCoEdgeLoop.get_and_step();
      if( i == best_corner_1 )
      {
        cornerCoEdge[0] = co_edge_ptr;
        continue;
      }

      if( i == best_corner_2 )
      {
        cornerCoEdge[1] = co_edge_ptr;
        continue;
      }

      if( start_vertex(co_edge_ptr) == start_vertex(cornerCoEdge[1]) )
      {
        best_corner_3 = i;
        cornerCoEdge[2] = co_edge_ptr;
        continue;
      }

      if( start_vertex(co_edge_ptr) == start_vertex(cornerCoEdge[0]) )
      {
        if( i==outerCoEdgeLoop.size() )
          best_corner_4 = i-1;
        else
          best_corner_4 = i;
        cornerCoEdge[3] = co_edge_ptr;
        break;
      }
    }

    if( best_corner_3==-1 || best_corner_4==-1 )
    {
      PRINT_ERROR( "unable to find corners\n" );
      return CUBIT_FAILURE;
    }

    fill_side_intervals( best_corner_1, best_corner_2, best_corner_3, best_corner_4 );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::find_loop_start( CoEdge *&start_co_edge_ptr )
{
  // For non-isLoop situation, start of loop will be first non-shared curve
  // on surface loop
  
  int i, j;
  RefFace *ref_face_ptr;
  if( !isLoop )
  {
    // First, get the outer curve loop of the surfaces
    refFaceChain.reset();
    ref_face_ptr = refFaceChain.get_and_step();

    DLIList<CoEdge*> co_edge_list;
    if( ordered_co_edges( ref_face_ptr, co_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    // If we have one surface just return
    co_edge_list.reset();
    if( refFaceChain.size() == 1 )
    {
      start_co_edge_ptr = co_edge_list.get();
      return CUBIT_SUCCESS;
    }

    // Get the list of coedges in the next surface in the chain
    ref_face_ptr = refFaceChain.get();
    DLIList<CoEdge*> next_co_edge_list;
    if( ordered_co_edges( ref_face_ptr, next_co_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    // Find a coedge on first surface that shares a curve with the next surface
    CoEdge *co_edge_ptr;
    RefEdge *ref_edge_ptr;
    co_edge_list.reset();
    int found = 0;
    for( i=co_edge_list.size(); i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      if( is_edge_in_list( ref_edge_ptr, next_co_edge_list ) )
      {
        found = 1;
        break;
      }
    }

    // This should never happen
    if( !found )
    {
      PRINT_ERROR( "Unable to find curves shared between surfaces\n" );
      return CUBIT_FAILURE;
    }

    // Loop until we find a curve not shared with the other surface - this
    // is the start
    for( i=co_edge_list.size(); i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      if( !is_edge_in_list( ref_edge_ptr, next_co_edge_list ) )
      {
        start_co_edge_ptr = co_edge_ptr;
        return CUBIT_SUCCESS;
      }
    }

    // This should never happen
    PRINT_ERROR( "Unable to find start of loop - please report this\n" );
    return CUBIT_FAILURE;
  }

  // isLoop situation
  if( refFaceChain.size() == 1 )
  {
    // Find first edge that has two coedges on this surface - this is the
    // start of the list
    refFaceChain.reset();
    ref_face_ptr = refFaceChain.get_and_step();

    DLIList<CoEdge*> co_edge_list;
    if( ordered_co_edges( ref_face_ptr, co_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    co_edge_list.reset();
    for( i=0; i<co_edge_list.size(); i++ )
    {
      RefEdge *ref_edge_ptr = co_edge_list.get()->get_ref_edge_ptr();
      if( number_coedges( ref_edge_ptr, ref_face_ptr ) == 2 )
      {
        if( i==0 ) // Backup until we are at the start
        {
          for( j=co_edge_list.size(); j--; )
          {
            co_edge_list.back();
            ref_edge_ptr = co_edge_list.get()->get_ref_edge_ptr();
            if( number_coedges( ref_edge_ptr, ref_face_ptr ) == 1 )
              break;
            co_edge_list.step();
          }
        }
        break;
      }
      co_edge_list.step();
    }

    start_co_edge_ptr = co_edge_list.get();
    return CUBIT_SUCCESS;
  }

  if( refFaceChain.size() == 2 )
  {
    refFaceChain.reset();
    RefFace *ref_face_ptr1 = refFaceChain.get_and_step();
    RefFace *ref_face_ptr2 = refFaceChain.get();
    CubitBoolean is_loop; // Dummy
    return check_for_loop( ref_face_ptr1, ref_face_ptr2, is_loop, start_co_edge_ptr );
  }

  // Case refFaceChain.size() > 2

  // Use first and last surfaces.  Also, adjust order of refFaceChain to
  // account for user's selection of start and end vertex, if selected.
  refFaceChain.reset();
  RefFace *start_surf = refFaceChain.get_and_back();
  RefFace *end_surf = refFaceChain.get();
  
  // Get list of coedges for the first surface
  DLIList<CoEdge*> start_co_edge_list;
  if( ordered_co_edges( start_surf, start_co_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

  // Get list of refedges for the last surface
  DLIList<RefEdge*> end_ref_edge_list;
  DLIList<Loop*> loop_list;
  end_surf->ordered_loops( loop_list );
  loop_list.get()->ordered_ref_edges( end_ref_edge_list );

  // Walk around the first surface
  RefEdge *ref_edge_ptr;

  // First, position the list at a starting coedge - one that is shared with
  // the other surface.
  start_co_edge_list.reset();
  for( i=0; i<start_co_edge_list.size(); i++ )
  {
    ref_edge_ptr = start_co_edge_list.get()->get_ref_edge_ptr();
    
    if( end_ref_edge_list.is_in_list( ref_edge_ptr ) )
    {
      if( i==0 ) // Backup until we are at the start
      {
        for( j=start_co_edge_list.size(); j--; )
        {
          start_co_edge_list.back();
          ref_edge_ptr = start_co_edge_list.get()->get_ref_edge_ptr();
          if( end_ref_edge_list.is_in_list( ref_edge_ptr ) )
            continue;
          start_co_edge_list.step();
        }
      }
      break;
    }
    start_co_edge_list.step();
  }
  
  // Now the list is positioned at the start of one of the shared chains of 
  // curves.
  start_co_edge_ptr = start_co_edge_list.get();

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_outer_vertices( DLIList<RefVertex*> &ref_vertex_list )
{
  int i;

  outerCoEdgeLoop.reset();
  for( i=outerCoEdgeLoop.size(); i--; )
    ref_vertex_list.append( start_vertex( outerCoEdgeLoop.get_and_step() ) );

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_outer_curves( DLIList<RefEdge*> &ref_edge_list )
{
  int i;
  
  outerCoEdgeLoop.reset();
  for( i=outerCoEdgeLoop.size(); i--; )
    ref_edge_list.append( outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr() );
  
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::populate_curve_loops()
{
  int i;
  RefFace *ref_face_ptr;
  TDSplitSurface *tdss;

  if( refFaceChain.size() == 1 )
  {
    ref_face_ptr = refFaceChain.get();
    tdss = (TDSplitSurface *)ref_face_ptr->get_TD(&TDSplitSurface::is_split_surface);
    tdss->add_coedges( outerCoEdgeLoop, sideInterval );
    return CUBIT_SUCCESS;
  }
  
  // Steps
  // 1. Go through each CoEdge on outside, adding TD to each CoEdge
  //    (with data for the start vertex of the CoEdge), marking them
  //    per the code in the diagram below. 
  //
  //    6    5   4
  //    +----+---+ 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //   7+----+---+3 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //   7+-+---+--+3 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //    +---+----+ 
  //    0   1    2 
  //
  //    
  // 2. Traverse through the curves in each surface, assigning 
  //    curves to proper side of surface utilizing vertex classification
  //    from above.
  
  // First gather a list of all coedges
  DLIList<CoEdge*> all_co_edge_list;
  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    ordered_co_edges( ref_face_ptr, all_co_edge_list );
  }
  
  refFaceChain.reset();
  RefFace *start_surf = refFaceChain.get();
  refFaceChain.last();
  RefFace *end_surf = refFaceChain.get();
  
  // STEP 1
  CoEdge *co_edge_ptr;
  RefVertex *start_vertex_ptr;
  int type = 0;
  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.back();
  CoEdge *prev_co_edge_ptr = outerCoEdgeLoop.get_and_step();
  for( i=0; i<outerCoEdgeLoop.size(); i++ )
  {
    co_edge_ptr = outerCoEdgeLoop.get_and_step();
    
    // Get all CoEdges whose start vertices are attached to the start of this
    // CoEdge.
    DLIList<CoEdge*> att_co_edge_list;
    get_attached_coedges_start( co_edge_ptr, all_co_edge_list, att_co_edge_list );
    // This should never happen, but check anyway
    if( att_co_edge_list.size() == 0 )
    {
      PRINT_ERROR( "Unable to find any curves attached to the start of Curve %d\n",
        co_edge_ptr->get_ref_edge_ptr()->id() );
      delete_coedge_tooldatas( all_co_edge_list );
      return CUBIT_FAILURE;
    }
    
    // Remove CoEdges from outerCoEdgeLoop
    att_co_edge_list -= outerCoEdgeLoop;
    
    // Add back in the current co_edge_ptr
    att_co_edge_list.append( co_edge_ptr );
    
    // Get starting vertex on this coedge
    start_vertex_ptr = start_vertex( co_edge_ptr );
    
    // Special case - need to remove possible wrong CoEdges.  If isLoop remove
    // all of the other CoEdges that are on the outerCoEdgeLoop (particularly
    // necessary as we traverse up and down the "seam", and at the corners).
    // If a corner, we need to deal with a possible triangle there by keeping
    // the proper CoEdge.
    if( isLoop )
    {
      // Find corner number from i
      int corner = -1;
      if( i == 0 )
        corner = 0;
      else if( i == sideInterval[0] )
        corner = 1;
      else if( i == sideInterval[0] + sideInterval[1] )
        corner = 2;
      else if( i== sideInterval[0] + sideInterval[1] + sideInterval[2] )
        corner = 3;
      
      if( corner != -1 )
      {        
        // Remove coedges that belong to the other corner shared by this corner
        if( remove_other_corner_coedges( prev_co_edge_ptr, co_edge_ptr,
          all_co_edge_list, att_co_edge_list ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
      }
    }
    
    // DEBUG
    //      if( att_co_edge_list.size() > 1 )
    //      {
    //        PRINT_INFO( "Att_co_edge_list.size() = %d\n", att_co_edge_list.size() );
    //        att_co_edge_list.reset();
    //        for( int k=att_co_edge_list.size(); k--; )
    //          PRINT_INFO( "  Adding TD to Curve = %d, Surface = %d\n", att_co_edge_list.get()->
    //          get_ref_edge_ptr()->id(), att_co_edge_list.get_and_step()->get_ref_face()->id() );
    //      }
    //      else
    //        PRINT_INFO( " Adding TD to Curve %d on Surface %d\n",
    //        att_co_edge_list.get()->get_ref_edge_ptr()->id(),
    //        att_co_edge_list.get()->get_ref_face()->id() );
    
    switch (type)
    {
    case 0:
      if( add_tooldata( att_co_edge_list, 0 ) == CUBIT_FAILURE )
      {
        delete_coedge_tooldatas( all_co_edge_list );
        return CUBIT_FAILURE;
      }
      if( cornerCoEdge[0] != cornerCoEdge[1] )
        type = 1;
      else
        type = 3;
      break;
    case 1:
      if( co_edge_ptr != cornerCoEdge[1] )
      {
        if( add_tooldata( att_co_edge_list, 1 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
      }
      else
      {
        if( add_tooldata( att_co_edge_list, 2 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
        if( cornerCoEdge[1] != cornerCoEdge[2] )
          type = 3;
        else
          type = 5;
      }
      break;
    case 3:
      if( co_edge_ptr != cornerCoEdge[2] )
      {
        if( add_tooldata( att_co_edge_list, 3 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
      }
      else
      {
        if( add_tooldata( att_co_edge_list, 4 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
        if( cornerCoEdge[2] != cornerCoEdge[3] )
          type = 5;
        else
          type = 7;
      }
      break;
    case 5:
      if( co_edge_ptr != cornerCoEdge[3] )
      {
        if( add_tooldata( att_co_edge_list, 5 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
      }
      else
      {
        if( add_tooldata( att_co_edge_list, 6 ) == CUBIT_FAILURE )
        {
          delete_coedge_tooldatas( all_co_edge_list );
          return CUBIT_FAILURE;
        }
        type = 7;
      }
      break;
    case 7:
      if( add_tooldata( att_co_edge_list, 7 ) == CUBIT_FAILURE )
      {
        delete_coedge_tooldatas( all_co_edge_list );
        return CUBIT_FAILURE;
      }
      break;
    }
    
    prev_co_edge_ptr = co_edge_ptr;
    
  }
    
  // Debug
  if( DEBUG_FLAG(154) )
  {
    outerCoEdgeLoop.reset();
    for( i=outerCoEdgeLoop.size(); i--; )
    {
      co_edge_ptr = outerCoEdgeLoop.get_and_step();
      TDSplitSurface *tdssv = (TDSplitSurface *)co_edge_ptr->
        get_TD(&TDSplitSurface::is_split_surface);
      
      PRINT_INFO( "Vertex %d - type %d\n", start_vertex( co_edge_ptr )->id(),
        tdssv->get_type() );
    }
  }
  
  // STEP 2
  refFaceChain.reset();
  for( i=0; i<refFaceChain.size(); i++ )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    tdss = (TDSplitSurface *)ref_face_ptr->
      get_TD(&TDSplitSurface::is_split_surface);
    
    // Get loop of coedges on this surface
    DLIList<CoEdge*> co_edge_list;
    ordered_co_edges( ref_face_ptr, co_edge_list );
    
    // Now, get to the start of the loop.  This will either be on vertex 0 or
    // the last vertex 7 (see diagram above).
    position_co_edge_list( i, co_edge_list );
    
    PRINT_DEBUG_154( "Surface %d: start of loop = curve %d\n", 
      ref_face_ptr->id(), co_edge_list.get()->get_ref_edge_ptr()->id() );
    
    // Now that we have the loop positioned properly, we can add the proper
    // curves (classified by side) to the TD on the surface.      
    DLIList<CoEdge*> side_coedges;
    int num_so_far = 0;
    // Check for special case - side A collapsed - triangle tip at start of
    // split
    if( ref_face_ptr == start_surf && cornerCoEdge[0] == cornerCoEdge[1] )
      ;
    else
      get_a_coedges( co_edge_list, side_coedges );
    num_so_far = side_coedges.size();
    
    // Side A could be collapsed (in case of triangle tip at start of split)
    if( side_coedges.size() )
      tdss->add_a_coedges( side_coedges );
    else
      tdss->add_a_coedges( side_coedges, start_vertex( cornerCoEdge[0] ) );
    
    side_coedges.clean_out();
    get_b_coedges( co_edge_list, side_coedges );
    num_so_far += side_coedges.size();
    // Side B could be collapsed (if we are on a triangle)
    if( side_coedges.size() )
      tdss->add_b_coedges( side_coedges );
    else
      tdss->add_b_coedges( side_coedges, start_vertex( co_edge_list.get() ) );
    
    side_coedges.clean_out();
    // Check for special case - side C collapsed - triangle tip at end of
    // split
    if( ref_face_ptr == end_surf && cornerCoEdge[2] == cornerCoEdge[3] )
      ;
    else
      get_c_coedges( co_edge_list, side_coedges );
    num_so_far += side_coedges.size();
    
    // Side C could be collapsed (in case of triangle tip at end of split)
    if( side_coedges.size() )
      tdss->add_c_coedges( side_coedges );
    else
      tdss->add_c_coedges( side_coedges, start_vertex( cornerCoEdge[2] ) );
    
    side_coedges.clean_out();
    if( get_d_coedges( co_edge_list, num_so_far, 
      side_coedges ) == CUBIT_FAILURE )
    {
      // Remove tooldatas from coedges
      delete_coedge_tooldatas( all_co_edge_list );
      return CUBIT_FAILURE;
    }
    
    // Side D could be collapsed (if we are on a triangle)
    if( side_coedges.size() )
      tdss->add_d_coedges( side_coedges );
    else
      tdss->add_d_coedges( side_coedges, start_vertex( co_edge_list.get() ) );
  }
  
  // Remove tooldatas from coedges
  delete_coedge_tooldatas( all_co_edge_list );

  return CUBIT_SUCCESS;  
}

CubitStatus
SplitSurfaceTool::ordered_co_edges( RefFace *ref_face_ptr, 
                                    DLIList<CoEdge*> &co_edge_list )
{
  DLIList<Loop*> loop_list;
  ref_face_ptr->ordered_loops( loop_list );

  if( loop_list.size() > 1 )
  {
    PRINT_ERROR( "Only surfaces with a single loop are allowed.\n"
      "       Surface %d has %d loops.\n", ref_face_ptr->id(), 
      loop_list.size() );
    return CUBIT_FAILURE;
  }

  if( loop_list.size() == 0 )
  {
    PRINT_ERROR( "Surface %d does not contain any loops - cannot split\n",
      ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  loop_list.get()->ordered_co_edges( co_edge_list );

  return CUBIT_SUCCESS;
}

CoEdge *
SplitSurfaceTool::get_next_co_edge( CoEdge *prev_co_edge_ptr, 
                                    DLIList<CoEdge*> &co_edge_list )
{
  RefVertex *end_vertex_ptr = end_vertex( prev_co_edge_ptr );

  // Find all potential connected coedges
  int i;
  CoEdge *co_edge_ptr;
  DLIList<CoEdge*> connected_co_edge_list;
  co_edge_list.reset();
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    if( end_vertex_ptr == start_vertex( co_edge_ptr ) )
      connected_co_edge_list.append( co_edge_ptr );
  }

  // We should always have one or two connected coedges - if not error
  if( connected_co_edge_list.size() == 0 )
  {
    PRINT_ERROR( "Didn't find connected coedge after Curve %d on Surface %d\n", 
      prev_co_edge_ptr->get_ref_edge_ptr()->id(), 
      prev_co_edge_ptr->get_ref_face()->id() );
    return NULL;
  }
  else if( connected_co_edge_list.size() == 1 )
  {
    return connected_co_edge_list.get();
  }
  else if( connected_co_edge_list.size() > 2 )
  {
    PRINT_ERROR( "Found %d coedges after Curve %d on Surface %d\n", 
      connected_co_edge_list.size(),
      prev_co_edge_ptr->get_ref_edge_ptr()->id(),
      prev_co_edge_ptr->get_ref_face()->id() );
    return NULL;
  }
  else if( connected_co_edge_list.size() == 2 )
  {
    // Avoid looping back on ourselves
    connected_co_edge_list.reset();
    if( prev_co_edge_ptr->get_ref_edge_ptr() == 
      connected_co_edge_list.get()->get_ref_edge_ptr() )
      return connected_co_edge_list.step_and_get();
    else
      return connected_co_edge_list.get();
  }

  return NULL;  
}

CubitBoolean
SplitSurfaceTool::is_edge_in_list( RefEdge *ref_edge_ptr, DLIList<CoEdge*> &co_edge_list )
{
  int i;
  CoEdge *co_edge_ptr;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    if( co_edge_ptr->get_ref_edge_ptr() == ref_edge_ptr )
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

RefVertex *
SplitSurfaceTool::start_vertex( CoEdge *co_edge_ptr )
{
  if( co_edge_ptr == NULL )
    return NULL;

  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  
  if ( co_edge_ptr->get_sense() == CUBIT_REVERSED )
    return ref_edge_ptr->end_vertex();
  else
    return ref_edge_ptr->start_vertex();
}

RefVertex *
SplitSurfaceTool::end_vertex( CoEdge *co_edge_ptr )
{
  if( co_edge_ptr == NULL )
    return NULL;

  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  
  if ( co_edge_ptr->get_sense() == CUBIT_REVERSED )
    return ref_edge_ptr->start_vertex();
  else
    return ref_edge_ptr->end_vertex();
}

CoEdge *
SplitSurfaceTool::get_complimentary_co_edge( CoEdge *co_edge_ptr, 
                                             DLIList<CoEdge*> co_edge_list )
{
  // Note: co_edge_list must be copied so that the calling code's list
  //       position is not changed.
  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  CoEdge *shared_co_edge_ptr;
  int i;
  for( i=co_edge_list.size(); i--; )
  {
    shared_co_edge_ptr = co_edge_list.get_and_step();
    if( shared_co_edge_ptr == NULL )
      continue;
    if( shared_co_edge_ptr == co_edge_ptr )
      continue;
    if( shared_co_edge_ptr->get_ref_edge_ptr() == ref_edge_ptr )
      return shared_co_edge_ptr;
  }
  return NULL;
}

CubitStatus
SplitSurfaceTool::get_attached_coedges_start( CoEdge *co_edge_ptr,
                                              DLIList<CoEdge*> &co_edge_list,
                                              DLIList<CoEdge*> &attached_co_edge_list )
{
  int i;
  CoEdge *attached_co_edge_ptr;

  RefVertex *start_vertex_ptr;
  start_vertex_ptr = start_vertex( co_edge_ptr );

  RefVertex *attached_start_vertex_ptr;

  co_edge_list.reset();
  for( i=co_edge_list.size(); i--; )
  {
    attached_co_edge_ptr = co_edge_list.get_and_step();

    attached_start_vertex_ptr = start_vertex( attached_co_edge_ptr );

    if( attached_start_vertex_ptr == start_vertex_ptr )
      attached_co_edge_list.append( attached_co_edge_ptr );
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus 
SplitSurfaceTool::remove_other_corner_coedges( CoEdge *co_edge_1,
                                               CoEdge *co_edge_2,
                                               DLIList<CoEdge*> &all_co_edge_list,
                                               DLIList<CoEdge*> &att_co_edge_list )
{
  // If only one attached coedge there is nothing to do
  if( att_co_edge_list.size() == 1 )
  {
    // It should be co_edge_2
    if( att_co_edge_list.get() != co_edge_2 )
    {
      PRINT_ERROR( "Internal error in split surface function - please report this.\n" );
      PRINT_DEBUG_154( "       Function remove_other_coedges.\n" );
      PRINT_DEBUG_154( "       Using curves %d and %d on surface %d, got curve %d\n",
        co_edge_1->get_ref_edge_ptr()->id(), co_edge_2->get_ref_edge_ptr()->id(),
        co_edge_1->get_ref_face()->id(), att_co_edge_list.get()->get_ref_edge_ptr()->id() );
      return CUBIT_FAILURE;
    }
    return CUBIT_SUCCESS;
  }

  // If we got this far, there are CoEdges we must examine to determine if they
  // belong with this corner or the other corner that shares this vertex.

  DLIList<CoEdge*> keep_co_edge_list;

  // Always keep co_edge_2 (the tooldata represents the start vertex on this coedge)
  keep_co_edge_list.append( co_edge_2 );
  att_co_edge_list.remove( co_edge_2 );

  // Traverse from co_edge_2
  DLIList<CoEdge*> keep_co_edge_list2;
  int done = 0;
  CoEdge *comp_co_edge_ptr = co_edge_2;
  CoEdge *co_edge_ptr;
  while( !done )
  {
    RefFace *ref_face_ptr = comp_co_edge_ptr->get_ref_face();
    DLIList<CoEdge*> co_edge_list;
    ordered_co_edges( ref_face_ptr, co_edge_list );
    co_edge_list.move_to( comp_co_edge_ptr );
    co_edge_list.back();
    co_edge_ptr = co_edge_list.get();
    if( co_edge_ptr == co_edge_1 )
      break;

    comp_co_edge_ptr = get_complimentary_co_edge( co_edge_ptr, all_co_edge_list );

    keep_co_edge_list2.append( comp_co_edge_ptr );

    // Avoid infinite loop, in case of something we didn't anticipate
    if( keep_co_edge_list2.size() > att_co_edge_list.size() )
    {
      PRINT_ERROR( "Internal error in split surface function - please report this.\n" );
      PRINT_DEBUG_154( "       Function remove_other_coedges #1.\n" );
      return CUBIT_FAILURE;
    }
  }

  int i;
  for( i=keep_co_edge_list2.size(); i--; )
  {
    co_edge_ptr = keep_co_edge_list2.get_and_step();
    if( !att_co_edge_list.is_in_list( co_edge_ptr ) )
    {
      PRINT_ERROR( "Internal error in split surface function - please report this.\n" );
      PRINT_DEBUG_154( "       Function remove_other_coedges #2.\n" );
      return CUBIT_FAILURE;
    }
  }

  keep_co_edge_list += keep_co_edge_list2;

  att_co_edge_list.clean_out();

  att_co_edge_list = keep_co_edge_list;

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::add_tooldata( DLIList<CoEdge*> &co_edge_list, int vertex_type )
{
  int i;
  CoEdge *co_edge_ptr;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();

    // Something's wrong if we have are trying to put a tooldata on a CoEdge that
    // already has one.
    if( (TDSplitSurface *)co_edge_ptr->get_TD(&TDSplitSurface::is_split_surface) )
    {
      PRINT_ERROR( "Internal error in split surface function - please report this.\n" );
      PRINT_DEBUG_154( "       Function add_tooldata.\n" );
      PRINT_DEBUG_154( "       Curve = %d on surface %d\n", co_edge_ptr->get_ref_edge_ptr()->id(),
        co_edge_ptr->get_ref_face()->id() );
      return CUBIT_FAILURE;
    }

    co_edge_ptr->add_TD( new TDSplitSurface( vertex_type ) );
  }
  return CUBIT_SUCCESS;
}

CubitStatus 
SplitSurfaceTool::pick_4_corners()
{
  // Special case - for the loop situation the corners are already picked
  if( isLoop )
    return CUBIT_SUCCESS;

  // Find 4 reasonable corners preferably close to PI/2.  We find at least 
  // 2 corners from the starting surface and 2 from the ending surface.
  
  // Note that at this point, the vertex loop starts at the beginning
  // of the first surface in the chain.

  // Check for trivial case
  if( outerCoEdgeLoop.size() == 4 )
  {
    outerCoEdgeLoop.reset();
    cornerCoEdge[0] = outerCoEdgeLoop.get_and_step();
    cornerCoEdge[1] = outerCoEdgeLoop.get_and_step();
    cornerCoEdge[2] = outerCoEdgeLoop.get_and_step();
    cornerCoEdge[3] = outerCoEdgeLoop.get();
    
    sideInterval[0] = 1;
    sideInterval[1] = 1;
    sideInterval[2] = 1;
    sideInterval[3] = 1;

    return CUBIT_SUCCESS;
  }

  // Initialize variables
  refFaceChain.reset();
  RefFace *start_surf = refFaceChain.get();
  refFaceChain.last();
  RefFace *end_surf = refFaceChain.get();
  outerCoEdgeLoop.reset();
  CoEdge *co_edge_ptr;
  int A_i=-1, B_i=-1, D_i=-1, E_i=-1;

  // Walk along the loop, to find the best corners for the first and last
  // surfaces.  The leading edge (B) tries to push forward, but the rear (A) 
  // tries to hang back.  Move the lead forward, then check to see if the rear
  // should move into it's old spot.

  // Note that the outerCoEdgeLoop starts at the first vertex of the starting
  // surface (shown after A in this diagram).
  // 
  // Final result for this patch results in corners where the 0's are shown.

  //
  //          E     D
  //          0-----+-------+------0 
  //          |     |       |      |
  //        F |     |       |      | 
  //   +------+     |       |      +------+
  //   |            |       |             |   ^
  //   |    End     |       |    Start    |   | Loop
  //   |  Surface   |       |   Surface   |   | Direction
  //   |            |       |             |   |
  //   +------+     |       |      +------+
  //          |     |       |      | C
  //          |     |       |      |
  //          0-----+-------+------0 B
  //                        A
  //

  double epsilon = .001745; // Angles considered comparable within about .1 degrees.
  int offset = 0;

  pick_expanded_corners( start_surf, outerCoEdgeLoop, offset, epsilon, A_i, B_i );
  PRINT_DEBUG_154( "First pass picked %d and %d on surf %d\n",
    A_i, B_i, start_surf->id() );

  // Advance the vertex and coedge lists until we are on the end_surf

  // Special case - if vertices are shared between start and end surface, we
  // want the opportunity to use the start vertex on the end surface, even
  // though we may have already picked it. We want to consider it again so
  // that we can find the best corners on the end surface.  Since we will 
  // always find at least two corners on a surface, if we did not consider
  // all vertices on the end surface, we might find a garbage vertex.  If
  // we do pick up a duplicate, we will clean it out and find the next best
  // vertex later.  If we did not do this, in the figure below, the first
  // pass through pick_expanded corners would find C1 and C2.  The next
  // pass would find C3 and C4.  We don't want C4, but we would pick it 
  // just because pick_expanded_corners must pick exactly two corners.  
  // Instead, lets pick C2 and C3 the next time around and pickup the lower
  // left corner as a cleanup step later.
  //                
  //               C4
  //   C1 +--------+--------------------+ C3
  //      |         \                   |
  //      |           \                 |
  //      |             \      end      |
  //      |               \             |
  //      |                 \           |
  //      |    start          \         |
  //      |                     \       |
  //      |                       \     |
  //      |                         \   |
  //      |                           \ |
  //      +--------+--------------------+ C2
  
  // Backup (to consider last vertex again, as explained above)
  outerCoEdgeLoop.back();
  offset--;
  
  int i;
  for( i=0; i<outerCoEdgeLoop.size(); i++ )
  {
    co_edge_ptr = outerCoEdgeLoop.get();
    if( is_vertex_in_surface( start_vertex(co_edge_ptr), end_surf ) )
      break;
    
    offset++;
    outerCoEdgeLoop.step();
  }

  pick_expanded_corners( end_surf, outerCoEdgeLoop, offset, epsilon, D_i, E_i );

  PRINT_DEBUG_154( "Second pass picked %d and %d on surf %d\n",
    D_i, E_i, end_surf->id() );

  if( E_i == outerCoEdgeLoop.size() )
    E_i = 0;

  // Error checking
  if( A_i==B_i )
    B_i = -1;
  if( A_i==D_i )
    D_i = -1;
  if( A_i==E_i )
    E_i = -1;
  if( B_i==D_i )
    D_i = -1;
  if( B_i==E_i )
    E_i = -1;
  if( D_i==E_i )
    E_i = -1;
  
  if( A_i==-1 || B_i==-1 || D_i==-1 || E_i==-1 )
  {
    // Most likely error case is that we only found three corners.  Here
    // is an example where this can happen.  Notice that the upper left
    // corner did not get selected.  This is because when processing the
    // left surface, the best corners were found that were "expanded out"
    // the farthest - thus we missed the "inside" corner.
    // 
    //    
    //    +-----------------------------+ C
    //    |                            /|
    //    |                          /  |
    //    |                        /    |
    //    |                      /      |
    //    |                    /        |
    //    |                  /          |
    //    |                /            |
    //    |              /              |
    //    |            /                |
    //    |          /                  |
    //  C +--------+--------------------+ C

    //
    // In this case, lets just cop out and find the next best corner(s)
    int num = 0;
    if( A_i==-1 )
      num++;
    if( B_i==-1 )
      num++;
    if( D_i==-1 )
      num++;
    if( E_i==-1 )
      num++;

    if( num>2 )
    {
      PRINT_ERROR( "Unable to pick corners - please specify manually.\n"
        "       Type 'help split surface' for proper syntax\n" );
      return CUBIT_FAILURE;
    }

    int best_corner_1 = -1, best_corner_2 = -1;

    double diff;
    double pi2 = CUBIT_PI/2.0;
    double min = CUBIT_DBL_MAX;
    outerCoEdgeLoop.reset();
    for( i=0; i<outerCoEdgeLoop.size(); i++ )
    {
      co_edge_ptr = outerCoEdgeLoop.get();

      if( i==A_i || i==B_i || i==D_i || i==E_i ||
          (!is_vertex_in_surface( start_vertex(co_edge_ptr), start_surf ) &&
           !is_vertex_in_surface( start_vertex(co_edge_ptr), end_surf ) ) )
      {
        outerCoEdgeLoop.step();
        continue;
      }

      diff = fabs( pi2 - compute_next_angle( outerCoEdgeLoop ) );

      if( diff < min )
      {
        best_corner_1 = i;
        min = diff;
      }
    }

    if( best_corner_1 == -1 )
    {
      PRINT_ERROR( "Unable to pick corners - please specify manually.\n"
        "       Type 'help split surface' for proper syntax\n" );
      return CUBIT_FAILURE;
    }

    if( num>1 )
    {
      min = CUBIT_DBL_MAX;
      outerCoEdgeLoop.reset();
      for( i=0; i<outerCoEdgeLoop.size(); i++ )
      {
        co_edge_ptr = outerCoEdgeLoop.get();
        
        if( i==A_i || i==B_i || i==D_i || i==E_i ||
          i==best_corner_1 ||
          (!is_vertex_in_surface( start_vertex(co_edge_ptr), start_surf ) &&
           !is_vertex_in_surface( start_vertex(co_edge_ptr), end_surf ) ) )
        {
          outerCoEdgeLoop.step();
          continue;
        }
        
        diff = fabs( pi2 - compute_next_angle( outerCoEdgeLoop ) );
        
        if( diff < min )
        {
          best_corner_2 = i;
          min = diff;
        }
      }
      if( best_corner_2 == -1 )
      {
        PRINT_ERROR( "Unable to pick corners - please specify manually.\n"
          "       Type 'help split surface' for proper syntax\n" );
        return CUBIT_FAILURE;
      }
    }

    int best_corner[6];

    // Find the four corners we have
    i = 0;
    if( A_i != -1 )
      best_corner[i++] = A_i;
    if( B_i != -1 )
      best_corner[i++] = B_i;
    if( D_i != -1 )
      best_corner[i++] = D_i;
    if( E_i != -1 )
      best_corner[i++] = E_i;
    if( best_corner_1 != -1 )
      best_corner[i++] = best_corner_1;
    if( best_corner_2 != -1 )
      best_corner[i++] = best_corner_2;

    if( i>4 )
    {
      PRINT_ERROR( "Internal error in split surface function - please report\n" );
      return CUBIT_FAILURE;
    }

    order_corners( best_corner[0], best_corner[1], best_corner[2], best_corner[3] );

    A_i = best_corner[0];
    B_i = best_corner[1];
    D_i = best_corner[2];
    E_i = best_corner[3];
  }

  // Fill cornerCoEdge array
  fill_corners( A_i, B_i, D_i, E_i );

  // Make sure we found at least two vertices on start surf and two vertices
  // on end surf
  if( !is_vertex_in_surface( start_vertex(cornerCoEdge[0]), start_surf ) ||
      !is_vertex_in_surface( start_vertex(cornerCoEdge[1]), start_surf ) ||
      !is_vertex_in_surface( start_vertex(cornerCoEdge[2]), end_surf ) ||
      !is_vertex_in_surface( start_vertex(cornerCoEdge[3]), end_surf ) )
  {
    PRINT_ERROR( "Unable to pick corners - please specify manually.\n"
      "       Type 'help split surface' for proper syntax\n" );
    return CUBIT_FAILURE;
  }

  // Fill sideInterval array
  fill_side_intervals( A_i, B_i, D_i, E_i );

  // Get loops in proper position (from start)
  reorient_loop( A_i );

  if( autoDetectTriangles == CUBIT_TRUE )
  {
    // Check to see if any corners are close to 180 and should be removed
    // to create a triangle or triangles
    if( update_corners_for_triangle() == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::pick_expanded_corners( RefFace *ref_face_ptr,
                                         DLIList<CoEdge*> &co_edge_list,
                                         int &offset, double epsilon, int &A_i, int &B_i )
{
  CoEdge *co_edge_ptr;
  double pi2 = CUBIT_PI/2.0; // Target angle (find vertices closest to this)
  double A = -1.0, B = -1.0, C;
  int old_B_i;
  double old_B;

  A_i = -1;
  B_i = -1;

  int initial_offset = offset;

  int i;
  for( i=0; i<co_edge_list.size(); i++ )
  {
    co_edge_ptr = co_edge_list.get();

    // If no longer in the surface we are done
    if( !is_vertex_in_surface( start_vertex(co_edge_ptr), ref_face_ptr ) )
    {
      A_i += initial_offset;
      B_i += initial_offset;
      return CUBIT_SUCCESS;
    }
    
    // Calculate difference of angle from PI/2.  Looking for smaller C as 
    // compared to previous A & B.  This also steps the lists.
    C = fabs( pi2 - compute_next_angle( co_edge_list ) );
    
    if( A == -1.0 ) // Initialization of A
    {
      A_i = i;
      A = C;
    }
    else if( B == -1.0 ) // Initialization of B
    {
      B_i = i;
      B = C;
    }
    // Check if B should move forward (it wants to!)
    else if( fabs(C - B) < epsilon || C < B ) // C == B || C < B
    {                                         // C <= B (effectively)
      // Push B forward to C's spot.
      old_B = B;
      old_B_i = B_i;
      B = C;
      B_i = i;
      
      // Check if A should move into B's old spot (it doesn't want to!)
      if( !(fabs(old_B - A)<epsilon) && old_B < A ) // old_B != A && old_B < A
      {                                             // old_B < A (effectively)
        // Push A forward to B's old spot
        A = old_B;
        A_i = old_B_i;
      }
    }
    // Also look for smaller C as compared to A
    else if( !(fabs(C - A)<epsilon) && C < A ) // C != A && C < A
    {                                          // C < A (effectively)
      // Push B forward to C's spot.
      old_B = B;
      old_B_i = B_i;
      B = C;
      B_i = i;
      
      // Push A forward to B's old spot
      A = old_B;
      A_i = old_B_i;
    }
    offset++;
  }

  A_i += initial_offset;
  B_i += initial_offset;
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::pick_4_corners_simple()
{
  // Special case - for the loop situation the corners are already picked
  if( isLoop )
    return CUBIT_SUCCESS;

  // We simply find the 4 corners closest to PI/2 - used
  // when we only have one surface.  Results sometimes 
  // aren't ideal....but user can pick corners in these
  // situations.

  // Store best quadruple
  int best_corner_1, best_corner_2, best_corner_3, best_corner_4;

  const int number_vertices = outerCoEdgeLoop.size();
  if( number_vertices == 3 )
  {
    // Triangle - special case (default to split to first tip of triangle)
    best_corner_1 = 0;
    best_corner_2 = 0;
    best_corner_3 = 1;
    best_corner_4 = 2;
  }
  else
  {
    const int num_angles = number_vertices;  
    
    // Compute angles
    double *angles;
    double turn_angle_sum = 0.;
    compute_angles( angles, turn_angle_sum );
    
    // For now just find the 4 vertices that are closest to PI/2
    int i;
    double diff;
    double pi2 = CUBIT_PI/2.0;
    double min = CUBIT_DBL_MAX;
    for( i=0; i<num_angles; i++ )
    {
      diff = fabs( pi2 - angles[i] );
      if( diff < min )
      {
        best_corner_1 = i;
        min = diff;
      }
    }
    
    min = CUBIT_DBL_MAX;
    for( i=0; i<num_angles; i++ )
    {
      diff = fabs( pi2 - angles[i] );
      if( diff < min && i != best_corner_1 )
      {
        best_corner_2 = i;
        min = diff;
      }
    }
    
    min = CUBIT_DBL_MAX;
    for( i=0; i<num_angles; i++ )
    {
      diff = fabs( pi2 - angles[i] );
      if( diff < min && i != best_corner_1 && i != best_corner_2 )
      {
        best_corner_3 = i;
        min = diff;
      }
    }
    
    min = CUBIT_DBL_MAX;
    for( i=0; i<num_angles; i++ )
    {
      diff = fabs( pi2 - angles[i] );
      if( diff < min && i != best_corner_1 && i != best_corner_2 && i != best_corner_3 )
      {
        best_corner_4 = i;
        min = diff;
      }
    }

    delete [] angles;
    
    order_corners( best_corner_1, best_corner_2, best_corner_3, best_corner_4 );
  }

  // Fill cornerCoEdge and sideInterval arrays
  fill_corners( best_corner_1, best_corner_2, best_corner_3, best_corner_4 );
  fill_side_intervals( best_corner_1, best_corner_2, best_corner_3, best_corner_4 );

  // Get loops in proper position (from start)
  reorient_loop( best_corner_1 );

  if( autoDetectTriangles == CUBIT_TRUE )
  {
    // Check to see if any corners are close to 180 and should be removed
    // to create a triangle
    if( update_corners_for_triangle() == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::update_corners_for_triangle()
{
  // Get four angles, in degrees
  double angle[4];
  outerCoEdgeLoop.reset();
  angle[0] = 180.0/CUBIT_PI * compute_next_angle( outerCoEdgeLoop );
  outerCoEdgeLoop.step( sideInterval[0]-1 );
  angle[1] = 180.0/CUBIT_PI * compute_next_angle( outerCoEdgeLoop );
  outerCoEdgeLoop.step( sideInterval[1]-1 );
  angle[2] = 180.0/CUBIT_PI * compute_next_angle( outerCoEdgeLoop );
  outerCoEdgeLoop.step( sideInterval[2]-1 );
  angle[3] = 180.0/CUBIT_PI * compute_next_angle( outerCoEdgeLoop );

  if( refFaceChain.size() > 1 )
  {
    // If a chain, do for both sides.  Check end surface first (otherwise
    // arrays could be previously modified when checking the other side)
    int corner_to_remove = -1;
    int collapsed_corner = -1;

    if( fabs(180.0-angle[2])<=sideAngleThreshold && angle[3] <= pointAngleThreshold )
    {
      corner_to_remove = 2;
      collapsed_corner = 3;
    }
    else if( fabs(180.0-angle[3])<=sideAngleThreshold && angle[2] <= pointAngleThreshold )
    {
      corner_to_remove = 3;
      collapsed_corner = 2;
    }

    if( corner_to_remove != -1 )
      remove_corner( corner_to_remove, collapsed_corner, CUBIT_FALSE );

    // Now check start surface
    corner_to_remove = -1;
    collapsed_corner = -1;

    if( fabs(180.0-angle[0])<=sideAngleThreshold && angle[1] <= pointAngleThreshold )
    {
      corner_to_remove = 0;
      collapsed_corner = 1;
    }
    else if( fabs(180.0-angle[1])<=sideAngleThreshold && angle[0] <= pointAngleThreshold )
    {
      corner_to_remove = 1;
      collapsed_corner = 0;
    }

    if( corner_to_remove != -1 )
      remove_corner( corner_to_remove, collapsed_corner, CUBIT_TRUE );
  }
  else
  {
    // Determine angle within sideAngleThreshold and closest to 180
    double diff = CUBIT_DBL_MAX;
    int best_side_corner = -1;
    int i;
    for( i=0; i<4; i++ )
    {
      if( fabs(180.0-angle[i])<=sideAngleThreshold && fabs(180.0-angle[i])<diff )
      {
        best_side_corner = i;
        diff = fabs(180.0-angle[i]);
      }
    }

    if( best_side_corner == -1 )
    {
      //PRINT_INFO( "No angle within criteria above sideAngleThreshold\n" );
      return CUBIT_SUCCESS;
    }

    //PRINT_INFO( "Best side corner = %d, angle = %f\n", best_side_corner, 
    //  angle[best_side_corner] );

    // Determine angle most below pointAngleThreshold
    diff = CUBIT_DBL_MAX;
    int best_point_corner = -1;
    for( i=0; i<4; i++ )
    {
      if( angle[i] <= pointAngleThreshold && (pointAngleThreshold-angle[i]) < diff )
      {
        best_point_corner = i;
        diff = pointAngleThreshold - angle[i];
      }
    }

    if( best_point_corner == -1 )
    {
      //PRINT_INFO( "No angle within criteria below pointAngleThreshold\n" );
      return CUBIT_SUCCESS;
    }

    //PRINT_INFO( "Best point corner = %d, angle = %f\n", best_point_corner, 
    //  angle[best_point_corner] );

    if( best_side_corner != -1 && best_point_corner != -1 )
      remove_corner( best_side_corner, best_point_corner, CUBIT_TRUE );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::remove_corner( int corner_to_remove, int collapsed_corner,
                                 CubitBoolean set_collapsed_first )
{
  // Need to adjust:
  //  outerCoEdgeLoop
  //  cornerCoEdge
  //  sideInterval

  // Shuffle the corners so that corner to remove is 2nd corner
  switch( corner_to_remove )
  {
  case 0:
    shuffle_corners_forward();
    shuffle_corners_forward();
    shuffle_corners_forward();
    collapsed_corner += 1;

    // Ensure collapsed_corner does not wrap the array (note we could have used
    // the remainder function instead, as follows:
    // int remainder = 0;
    // collapsed_corner = (remainder=collapsed_corner%4) ? remainder : collapsed_corner;
    if( collapsed_corner == 4 )
      collapsed_corner = 0;    

    break;
  case 1:
    break;
  case 2:
    shuffle_corners_forward();
    collapsed_corner -= 1;
    break;
  case 3:
    shuffle_corners_forward();
    shuffle_corners_forward();
    collapsed_corner -= 2;
    break;
  }

  if( collapsed_corner < 0 )
    collapsed_corner += 4;

  switch( collapsed_corner )
  {
  case 0:
    sideInterval[1] = sideInterval[0] + sideInterval[1];
    sideInterval[0] = 0;
    break;
  case 1:
    PRINT_ERROR( "Algorithm error in auto-detecting triangle - please report.\n" );
    return CUBIT_FAILURE;
  case 2:
    sideInterval[0] = sideInterval[0] + sideInterval[1];
    sideInterval[1] = 0;
    break;
  case 3:
    sideInterval[0] = sideInterval[0] + sideInterval[1];
    sideInterval[1] = sideInterval[2];
    sideInterval[2] = 0;
  }

  // Now redo cornerCoEdge array
  outerCoEdgeLoop.reset();
  cornerCoEdge[0] = outerCoEdgeLoop.get();
  outerCoEdgeLoop.step( sideInterval[0] );
  cornerCoEdge[1] = outerCoEdgeLoop.get();
  outerCoEdgeLoop.step( sideInterval[1] );
  cornerCoEdge[2] = outerCoEdgeLoop.get();
  outerCoEdgeLoop.step( sideInterval[2] );
  cornerCoEdge[3] = outerCoEdgeLoop.get();

  if( set_collapsed_first == CUBIT_TRUE )
  {
    // Adjust the arrays so that collapsed_corner is at the front
    switch( collapsed_corner )
    {
    case 2:
      shuffle_corners_forward();
      break;
    case 3:
      shuffle_corners_forward();
      shuffle_corners_forward();
      shuffle_corners_forward();
    }
  }
  else
  {
    // Return to original position
    switch( corner_to_remove )
    {
    case 0:
      shuffle_corners_forward();
      break;
    case 1:
      break;
    case 2:
      shuffle_corners_forward();
      shuffle_corners_forward();
      shuffle_corners_forward();
      break;
    case 3:
      shuffle_corners_forward();
      shuffle_corners_forward();
      break;
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::fill_corners( int best_corner_1, int best_corner_2,
                                int best_corner_3, int best_corner_4 )
{
  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.step(best_corner_1);
  cornerCoEdge[0] = outerCoEdgeLoop.get();

  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.step(best_corner_2);
  cornerCoEdge[1] = outerCoEdgeLoop.get();

  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.step(best_corner_3);
  cornerCoEdge[2] = outerCoEdgeLoop.get();

  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.step(best_corner_4);
  cornerCoEdge[3] = outerCoEdgeLoop.get();

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::fill_side_intervals( int best_corner_1, int best_corner_2,
                                       int best_corner_3, int best_corner_4 )
{
  sideInterval[0] = best_corner_2 - best_corner_1;
  sideInterval[1] = best_corner_3 - best_corner_2;
  sideInterval[2] = best_corner_4 - best_corner_3;
  sideInterval[3] = best_corner_1 + outerCoEdgeLoop.size() - best_corner_4;

  return CUBIT_SUCCESS;
}

CubitStatus 
SplitSurfaceTool::compute_angles( double *&angles, double &turn_angle_sum )
{
  const int num_angles = outerCoEdgeLoop.size();
  angles = new double[ num_angles ];
  turn_angle_sum = 0.;
  outerCoEdgeLoop.reset();
  int i;
  for(i=0; i<num_angles; i++) 
  {
    angles[i] = compute_next_angle( outerCoEdgeLoop );
    //turn_angle_sum += CUBIT_PI - angles[i]; // Not currently used
  }
  return CUBIT_SUCCESS;
}

double
SplitSurfaceTool::compute_next_angle( DLIList<CoEdge*> &co_edge_list )
{  
  CoEdge *co_edge_1 = co_edge_list.prev();
  CoEdge *co_edge_2 = co_edge_list.get_and_step();

  // Coordinates of common point
  CubitVector vertex_point = start_vertex( co_edge_2 )->coordinates();
  
  // Find normal to the the face at the common vertex
  RefFace *ref_face1 = co_edge_1->get_ref_face();
  CubitVector normal1 = ref_face1->normal_at(vertex_point, NULL);

  RefFace *ref_face2 = co_edge_2->get_ref_face();
  CubitVector normal2 = ref_face2->normal_at(vertex_point, NULL);

  CubitVector normal = (normal1 + normal2)/2; // Average

  // Find directed tangents to determine interior angle
  // Use sense of edge with respect to this face's loop.
  CubitVector tangent_1, tangent_2;
  co_edge_1->get_ref_edge_ptr()->tangent( vertex_point, tangent_1 );
  co_edge_2->get_ref_edge_ptr()->tangent( vertex_point, tangent_2 );
  
  if ( co_edge_1->get_sense() == CUBIT_REVERSED )
    tangent_1 = -tangent_1;
  if ( co_edge_2->get_sense() == CUBIT_REVERSED )
    tangent_2 = -tangent_2;
  
  //  At this point we have the tangents going in the correct loop
  //  sense.
  //  Now get tangent pointing away from the center for the correct
  //  angle
  tangent_1 = -tangent_1;

  // Return angle from given tangents and normal to face
  double angle = normal.vector_angle( tangent_2, tangent_1 );

  PRINT_DEBUG_154( " Angle at vertex %d = %f\n", start_vertex( co_edge_2 )->id(), 
                   angle*180.0/CUBIT_PI );
  return angle;
}

void
SplitSurfaceTool::order_corners( int &corner_1, int &corner_2,
				                         int &corner_3, int &corner_4 )
{
  // sort the corners by indices
  int swap_temp;
#define SWAP( a, b ) swap_temp = (a); (a) = (b); (b) = swap_temp
  if ( corner_1 > corner_2 ) {
    SWAP( corner_1, corner_2 );
  }
  if ( corner_2 > corner_3 ) {
    SWAP( corner_2, corner_3 );
  }
  if ( corner_3 > corner_4 ) {
    SWAP( corner_3, corner_4 );
  }
  // 4 is now set
  order_corners( corner_1, corner_2, corner_3 );
#undef SWAP
}

void
SplitSurfaceTool::order_corners( int &corner_1, int &corner_2, 
                                 int &corner_3 )
{
  // sort the corners by indices
  int swap_temp;
#define SWAP( a, b ) swap_temp = (a); (a) = (b); (b) = swap_temp
  if ( corner_1 > corner_2 ) {
    SWAP( corner_1, corner_2 );
  }
  if ( corner_2 > corner_3 ) {
    SWAP( corner_2, corner_3 );
  }
  // 3 is now set
  if ( corner_1 > corner_2 ) {
    SWAP( corner_1, corner_2 );
  }
  // 2 is now set
  // 1 is set as its the only one left
#undef SWAP
}

CubitStatus 
SplitSurfaceTool::order_selected_corners( DLIList<RefVertex*> &corner_vertex_list )
{
  int i;

  RefVertex *ref_vertex_ptr;
  CoEdge *co_edge_ptr;

  // Verify that these corners exist in the outerCoEdgeLoop
  for( i = corner_vertex_list.size(); i--; )
  {
    ref_vertex_ptr = corner_vertex_list.get_and_step();

    if( is_in_outer_loop( ref_vertex_ptr ) == CUBIT_FALSE )
    {
      PRINT_ERROR( "Selected corner, Vertex %d, not found in outer loop of patch.\n",
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }
  }

  // We need to start at the first corner the user selected - shuffle the
  // lists to do this.
  corner_vertex_list.reset();
  RefVertex *start_corner = corner_vertex_list.get();

  int offset = 0;
  outerCoEdgeLoop.reset();
  for( i=outerCoEdgeLoop.size(); i--; )
  {
    ref_vertex_ptr = start_vertex( outerCoEdgeLoop.get_and_step() );
    if( ref_vertex_ptr == start_corner )
      break;
    offset++;
  }

  // Now shuffle the lists to get them into the proper order
  if( offset )
  {
    outerCoEdgeLoop.reset();
    DLIList<CoEdge*> outer_co_edge_loop = outerCoEdgeLoop;
    
    outerCoEdgeLoop.clean_out();

    outer_co_edge_loop.reset();
    
    // Set the loop at start offset
    outer_co_edge_loop.step( offset );
    for( i=outer_co_edge_loop.size(); i--; )
    {
      outerCoEdgeLoop.append( outer_co_edge_loop.get_and_step() );
    }
  }

  // Cruise through the outer coedge list to get the proper corner order 
  int corner = 0;
  int num = 0;
  int side_interval = 0;
  outerCoEdgeLoop.reset();

  co_edge_ptr = outerCoEdgeLoop.get_and_step();

  // Special case for first corner (always at first corner to start)
  num = number_in_list( corner_vertex_list, start_vertex(co_edge_ptr) );
  if( num > 2 )
  {
    PRINT_ERROR( "Vertex %d was selected for %d of the corners - max is 2\n",
      start_vertex(co_edge_ptr)->id(), num );
    return CUBIT_FAILURE;
  }

  cornerCoEdge[corner++] = co_edge_ptr;
  if( num == 2 )
  {
    sideInterval[corner-1] = 0;
    cornerCoEdge[corner++] = co_edge_ptr;
  }

  for( i=1; i<outerCoEdgeLoop.size(); i++ )
  {
    co_edge_ptr = outerCoEdgeLoop.get_and_step();
    side_interval++;

    num = number_in_list( corner_vertex_list, start_vertex(co_edge_ptr) );
    if( num > 2 )
    {
      PRINT_ERROR( "Vertex %d was selected for %d of the corners - max is 2\n",
        start_vertex(co_edge_ptr)->id(), num );
      return CUBIT_FAILURE;
    }

    if( num > 0 )
    {
      sideInterval[corner-1] = side_interval;
      side_interval = 0;
      cornerCoEdge[corner++] = co_edge_ptr;

      if( num == 2 )
      {
        sideInterval[corner-1] = 0;
        cornerCoEdge[corner++] = co_edge_ptr;
      }
    }
  }

  sideInterval[3] = outerCoEdgeLoop.size() - sideInterval[0] - sideInterval[1] - sideInterval[2];

  return CUBIT_SUCCESS;
}

int
SplitSurfaceTool::number_in_list( DLIList<RefVertex*> &corner_vertex_list, 
                                  RefVertex *ref_vertex_ptr )
{
  int num = 0;
  int i;
  corner_vertex_list.reset();
  RefVertex *check_ptr;
  for( i=4; i--; )
  {
    check_ptr = corner_vertex_list.get_and_step();
    if( check_ptr == ref_vertex_ptr )
      num++;
  }

  return num;
}

CubitBoolean
SplitSurfaceTool::is_in_outer_loop( RefVertex *ref_vertex_ptr )
{
  int i;
  outerCoEdgeLoop.reset();

  for( i=outerCoEdgeLoop.size(); i--; )
  {
    if( ref_vertex_ptr == start_vertex( outerCoEdgeLoop.get_and_step() ) )
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitStatus
SplitSurfaceTool::adjust_for_split_direction( RefEdge *curve_dir_ptr,
                                              RefEdge *from_curve_ptr,
                                              DLIList<RefVertex*> &corner_vertex_list )
{
  // Adjust the loops to the split direction. The split direction is
  // determined as follows:
  // 1) For multiple surfaces - always across the chain.
  //
  //          split direction
  //   c2        <-------->            c1
  //   +---------+-----------+-----------+
  //   |         |           |           |
  //   |_  _  _  |_  _  _  _ | _  _  _  _|
  //   |         |  split    |           |side 0
  //   |         |           |           |
  //   +---------+-----------+-----------+
  //   c3            side 3              c0
  //
  // 2) If a "from" curve is specfied, put the "from" curve on 
  //    side 3, so that the fraction and distance goes from the
  //    right side.


  // For single surfaces, the logic is a little more tricky:
  // 1) From optional "direction" curve given by user (parallel to split) 
  // 2) From optional "from" curve given by user (parallel to split)
  // 3) From optional through vertex list given by user.
  // 4) From optional corner vertices given by user, if given in a 
  //    sensible order (ie., not criss-crossed). Split direction
  //    will be perpendicular to the side from first to second
  //    selected corner (ie., perpendicular to c0-c1 in the 
  //    diagram below).
  // 5) From aspect ratio of patch - split along "skinny" direction.
 
  //          split direction
  //   c2        <-------->            c1
  //   +-------------------------------+
  //   |                               |
  //   |_  _  _  _  _  _  _  _  _  _  _|
  //   |           split               |
  //   |                               |
  //   +-------------------------------+
  //   c3                              c0
  //
  //   In all cases check for incompatible inputs (i.e., conflicting 
  //   "direction" and "from" curves, conflicting "through" vertices
  //   on single surfaces, etc.)

  corner_vertex_list.reset();

  if( refFaceChain.size() > 1 && corner_vertex_list.size() )
  {
    refFaceChain.reset();
    RefFace *start_surf = refFaceChain.get();
    refFaceChain.last();
    RefFace *end_surf = refFaceChain.get();

    // We might have to shuffle the selected corners to orient them to split 
    // across the surfaces.  To be valid, we must have the two starting 
    // corners in the start surface and the two ending corners in the end 
    // surface.  Shuffle the corners until this condition is met.  If
    // we shuffle three times and the condition is still not met, we have
    // been given invalid vertices.
    //
    // Here is a case we must be able to handle, with corners specified
    // as shown.  We must shuffle three times.

    //  C3                        C2
    //   +-----+------------------+
    //   |\     \                 |
    //   |  \    \                |
    //   |    \    \     end      |
    //   |      \   \             |
    //   |        \   \           |
    //   |          \   \         |
    //   |            \   \       |
    //   |              \  \      |
    //   |                \  \    |
    //   |     start        \ \   |
    //   |                    \ \ |
    //   |                      \\|
    //   +------------------------+
    //  C0                        C1

    if( !is_vertex_in_surface( start_vertex(cornerCoEdge[0]), start_surf ) ||
        !is_vertex_in_surface( start_vertex(cornerCoEdge[1]), start_surf ) ||
        !is_vertex_in_surface( start_vertex(cornerCoEdge[2]), end_surf ) ||
        !is_vertex_in_surface( start_vertex(cornerCoEdge[3]), end_surf ) )
    {
      shuffle_corners_forward();
      if( !is_vertex_in_surface( start_vertex(cornerCoEdge[0]), start_surf ) ||
          !is_vertex_in_surface( start_vertex(cornerCoEdge[1]), start_surf ) ||
          !is_vertex_in_surface( start_vertex(cornerCoEdge[2]), end_surf ) ||
          !is_vertex_in_surface( start_vertex(cornerCoEdge[3]), end_surf ) )
      {
        shuffle_corners_forward();
        if( !is_vertex_in_surface( start_vertex(cornerCoEdge[0]), start_surf ) ||
            !is_vertex_in_surface( start_vertex(cornerCoEdge[1]), start_surf ) ||
            !is_vertex_in_surface( start_vertex(cornerCoEdge[2]), end_surf ) ||
            !is_vertex_in_surface( start_vertex(cornerCoEdge[3]), end_surf ) )
        {
          shuffle_corners_forward();
          if( !is_vertex_in_surface( start_vertex(cornerCoEdge[0]), start_surf ) ||
              !is_vertex_in_surface( start_vertex(cornerCoEdge[1]), start_surf ) ||
              !is_vertex_in_surface( start_vertex(cornerCoEdge[2]), end_surf ) ||
              !is_vertex_in_surface( start_vertex(cornerCoEdge[3]), end_surf ) )
          {
            PRINT_ERROR( "Selected vertices are invalid.  Two vertices must be on\n"
              "       start surface and two on end surface of chain.\n" );
            return CUBIT_FAILURE;
          }
        }
      }
    }
  }
  if( refFaceChain.size() > 1 )
  {
    if( curve_dir_ptr )
    {
      // Check for rare case with two triangles - that case allows a switch in
      // direction
      if( is_chain_two_triangles() )
      {
        if( is_curve_on_side( curve_dir_ptr, 0 ) ||
            is_curve_on_side( curve_dir_ptr, 2 ) )
        {
          shuffle_corners_forward();
        }

        if( throughVertexList.size() > 2 )
        {
          PRINT_ERROR( "For two triangle case, you can only specify 2 'through' vertices.\n" );
          return CUBIT_FAILURE;
        }

        // Make sure through vertices (if any) are on side 0 or 2
        if( check_through_vertices( "'direction'" ) == CUBIT_FAILURE )
          return CUBIT_FAILURE;
      }
      else
        PRINT_WARNING( "Curve direction ignored - not necessary for a chain of surfaces.\n" );
    }
      
    // Flip sides if user selected a from curve
    if( from_curve_ptr )
    {
      // Check for rare case with two triangles - that case allows a switch in
      // direction
      if( is_chain_two_triangles() )
      {
        if( is_curve_on_side( from_curve_ptr, 0 ) ||
            is_curve_on_side( from_curve_ptr, 2 ) )
        {
          if( curve_dir_ptr )
          {
            PRINT_ERROR( "From curve must be parallel to direction of split.\n" );
            return CUBIT_FAILURE;
          }
          else
            shuffle_corners_forward();
        }

        // Next two checks are redundant if a direction was specified, but that 
        // doesn't hurt anything.
        if( throughVertexList.size() > 2 )
        {
          PRINT_ERROR( "For two triangle case, you can only specify 2 'through' vertices.\n" );
          return CUBIT_FAILURE;
        }

        // Make sure through vertices (if any) are on side 0 or 2
        if( check_through_vertices( "'from'" ) == CUBIT_FAILURE )
          return CUBIT_FAILURE;
      }
      else
      {
        // Make sure this curve is on side 3
        if( is_curve_on_side( from_curve_ptr, 0 ) ||
            is_curve_on_side( from_curve_ptr, 2 ) )
        {
          PRINT_ERROR( "From curve must be parallel to direction of split.\n" );
          return CUBIT_FAILURE;
        }
        if( is_curve_on_side( from_curve_ptr, 1 ) )
        {
          shuffle_corners_forward();
          shuffle_corners_forward();
        }
      }   
    }
    return CUBIT_SUCCESS;
  }
  else
    // SINGLE surface
  {
    // Special case for triangle
    if( is_triangle() && (curve_dir_ptr || from_curve_ptr ) )
    {
      // For a triangle, curve_dir_ptr and from_curve_ptr mean the same
      // thing.  Do a check.
      if( curve_dir_ptr && from_curve_ptr && (from_curve_ptr != curve_dir_ptr) )
      {
        PRINT_ERROR( "For triangle case, specifying both a 'direction' and 'from'\n"
          "       curve is redundant - these have the same meaning.\n" );
        return CUBIT_FAILURE;
      }

      if( from_curve_ptr && !curve_dir_ptr )
        curve_dir_ptr = from_curve_ptr;

      // Desired state: 0 interval side opposite the side that curve_dir_ptr
      // is on.  Then ensure that side past selected curve is side 0.  A
      // tricky bit of logic!        

      // Start out by setting zero side to side 0
      if( sideInterval[1] == 0 )
      {
        shuffle_corners_forward();          
      }
      else if (sideInterval[2] == 0 )
      {
        shuffle_corners_forward();
        shuffle_corners_forward();
      }
      else if( sideInterval[3] == 0 )
      {
        shuffle_corners_forward();
        shuffle_corners_forward();
        shuffle_corners_forward();
      }

      // Move zero side and ensure that side past selected curve is side 0.
      if( is_curve_on_side( curve_dir_ptr, 1 ) )
      {
        // Zero side should be 3
        shuffle_zero_side_forward();
        shuffle_zero_side_forward();

        shuffle_corners_forward();
      }
      else if( is_curve_on_side( curve_dir_ptr, 2 ) )
      {
        shuffle_corners_forward();
        shuffle_corners_forward();
        shuffle_corners_forward();
      }
      else if( is_curve_on_side( curve_dir_ptr, 3 ) )
      {
        shuffle_zero_side_forward();
      }

      // Check through vertices to make sure they are valid
      if( check_through_vertices( "'direction'" ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;

      return CUBIT_SUCCESS;
    }
    else if( curve_dir_ptr )
    {
      // Curve should be on side 1 or 3 - otherwise adjust the loops forward
      if( is_curve_on_side( curve_dir_ptr, 0 ) )
        shuffle_corners_forward();
      else if( is_curve_on_side( curve_dir_ptr, 2 ) )
        shuffle_corners_forward();

      // Now check if there is a from curve too, along with curve direction
      if( from_curve_ptr )
      {
        if( is_curve_on_side( from_curve_ptr, 0 ) ||
          is_curve_on_side( from_curve_ptr, 2 ) )
        {
          PRINT_ERROR( "'Direction' and 'From' curves must be on same or opposite\n"
            "       sides of logical rectangle\n" );
          return CUBIT_FAILURE;
        }
        // Must be on side 3
        if( is_curve_on_side( from_curve_ptr, 1 ) )
        {
          shuffle_corners_forward();
          shuffle_corners_forward();
        }
      }

      // Check through vertices to make sure they are valid
      if( check_through_vertices( "'direction'" ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;

      return CUBIT_SUCCESS;
    }

    // Single surface WITHOUT a direction given
    if( from_curve_ptr )
    {
      // Curve should be on side 3 - otherwise adjust the loops forward
      if( is_curve_on_side( from_curve_ptr, 2 ) )
      {
        shuffle_corners_forward();
        shuffle_corners_forward();
        shuffle_corners_forward();
      }
      else if( is_curve_on_side( from_curve_ptr, 1 ) )
      {
        shuffle_corners_forward();
        shuffle_corners_forward();
      }
      else if( is_curve_on_side( from_curve_ptr, 0 ) )
        shuffle_corners_forward();

      // Check through vertices to make sure they are valid
      if( check_through_vertices( "'from'" ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;

      return CUBIT_SUCCESS;
    }

    // Single surface WITHOUT a 'direction' or 'from' curve given - if through
    // vertices let them determine the split direction.
    if( throughVertexList.size() > 2 )
    {
      PRINT_ERROR( "For a single surface, you can only specify 2 'through' vertices.\n" );
      return CUBIT_FAILURE;
    }

    RefVertex *start_vertex_ptr = NULL;
    RefVertex *end_vertex_ptr = NULL;
    throughVertexList.reset();
    if( throughVertexList.size() > 0 )
      start_vertex_ptr = throughVertexList.get_and_step();
    if( throughVertexList.size() > 1 )
      end_vertex_ptr = throughVertexList.get();

    // First check for start and end vertices
    if( start_vertex_ptr || end_vertex_ptr )
    {
      if( start_vertex_ptr )
      {
        // Check if on side A
        if( !is_vertex_on_side( start_vertex_ptr, 0 ) )
        {
          // Shuffle forward and check again
          shuffle_corners_forward();
          if( !is_vertex_on_side( start_vertex_ptr, 0 ) )
          {
            // Shuffle forward and check again
            shuffle_corners_forward();
            if( !is_vertex_on_side( start_vertex_ptr, 0 ) )
            {
              // Shuffle forward and check again
              shuffle_corners_forward();
              if( !is_vertex_on_side( start_vertex_ptr, 0 ) )
              {
                PRINT_ERROR( "Specified 'through' vertex %d is not valid\n",
                  start_vertex_ptr->id() );
                return CUBIT_FAILURE;
              }
            }
          }
        }
      }
      if( end_vertex_ptr )
      {
        // Check if on side C
        if( !is_vertex_on_side( end_vertex_ptr, 2 ) )
        {
          // Shuffle forward and check again
          shuffle_corners_forward();
          if( !is_vertex_on_side( end_vertex_ptr, 2 ) )
          {
            // Shuffle forward and check again
            shuffle_corners_forward();
            if( !is_vertex_on_side( end_vertex_ptr, 2 ) )
            {
              // Shuffle forward and check again
              shuffle_corners_forward();
              if( !is_vertex_on_side( end_vertex_ptr, 2 ) )
              {
                PRINT_ERROR( "Specified 'through' vertex %d is not valid\n",
                  end_vertex_ptr->id() );
                return CUBIT_FAILURE;
              }
            }
          }
        }
      }
      if( start_vertex_ptr && end_vertex_ptr )
      {
        if( !is_vertex_on_side( start_vertex_ptr, 0 ) ||
            !is_vertex_on_side( end_vertex_ptr, 2 ) )
        {
          PRINT_ERROR( "Specified 'through' vertices invalid - both not on split path.\n" );
          return CUBIT_FAILURE;
        }
      }

      return CUBIT_SUCCESS;

    }
    
    // If second corner is 1 or 3, use the direction implied.  Otherwise,
    // user gave corners in a criss-crossed manner and we can't gain any
    // implied direction from them.
    else if( corner_vertex_list.size() && 
             corner_vertex_list.step_and_get() == start_vertex( cornerCoEdge[1]) )
      return CUBIT_SUCCESS;
    else if( corner_vertex_list.size() && 
             corner_vertex_list.get() == start_vertex( cornerCoEdge[3]) )
      ;
    else
    {
      // We use aspect ratio to determine split direction.  Split along
      // narrowest direction.

      // Skip if we have a triangle
      if( is_triangle() )
        return CUBIT_SUCCESS;

      double len0 = get_side_length( 0 );
      double len1 = get_side_length( 1 );
      double len2 = get_side_length( 2 );
      double len3 = get_side_length( 3 );

      double ratio = (len0+len2) / (len1+len3);
      PRINT_DEBUG_154( "Ratio = %f\n", ratio );

      if( ratio < 1.0001 )
        return CUBIT_SUCCESS;
      else
        PRINT_DEBUG_154( "Ratio deems we adjust forward\n" );
    }
    
    // If we got this far, we need to adjust the split direction - move
    // loops forward by sideInterval[0] (start at corner 1)
    shuffle_corners_forward();
  }

  return CUBIT_SUCCESS;
}

CubitBoolean
SplitSurfaceTool::is_chain_two_triangles()
{
  if( refFaceChain.size() != 2 )
    return CUBIT_FALSE;

  //  C3-0-1-2                 C2-3-0-1
  //   +------------------------+
  //   |\                       |
  //   |  \                     |
  //   |    \          end      |
  //   |      \                 |
  //   |        \               |
  //   |          \             |
  //   |            \           |
  //   |              \         |
  //   |                \       |
  //   |     start        \     |
  //   |                    \   |
  //   |                      \ |
  //   +------------------------+
  //  C0-1-2-3                 C1-2-3-0

  // 4 possibilities (see diagram above - corners 0-1-2-3 or 1-2-3-0 etc
  // going ccw from bottom left)

  CoEdge *c0 = cornerCoEdge[0];
  CoEdge *c1 = cornerCoEdge[1];
  CoEdge *c2 = cornerCoEdge[2];
  CoEdge *c3 = cornerCoEdge[3];

  if( (prev_co_edge(c0)->get_ref_face() == c0->get_ref_face() &&
      prev_co_edge(c1)->get_ref_face() != c1->get_ref_face() &&
      prev_co_edge(c2)->get_ref_face() == c2->get_ref_face() &&
      prev_co_edge(c3)->get_ref_face() != c3->get_ref_face() ) ||

      (prev_co_edge(c1)->get_ref_face() == c1->get_ref_face() &&
      prev_co_edge(c2)->get_ref_face() != c2->get_ref_face() &&
      prev_co_edge(c3)->get_ref_face() == c3->get_ref_face() &&
      prev_co_edge(c0)->get_ref_face() != c0->get_ref_face() ) )
  {
    // 2-triangle case
    return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

CubitBoolean
SplitSurfaceTool::is_triangle()
{
  if( sideInterval[0]==0 || sideInterval[1]==0 ||
    sideInterval[2]==0 || sideInterval[3]==0 )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

CubitStatus
SplitSurfaceTool::check_through_vertices( char *type )
{
  if( throughVertexList.size() )
  {
    RefVertex *start_vertex_ptr = NULL;
    RefVertex *end_vertex_ptr = NULL;
    throughVertexList.reset();
    start_vertex_ptr = throughVertexList.get_and_step();
    if( throughVertexList.size() > 1 )
      end_vertex_ptr = throughVertexList.get();
    
    if( start_vertex_ptr )
    {
      if( is_vertex_on_side( start_vertex_ptr, 0 ) || 
          is_vertex_on_side( start_vertex_ptr, 2 ) )
        ;
      else
      {
        if( throughVertexList.size() == 1 )
          PRINT_ERROR( "Through vertice(s) not compatible with %s curve given\n",
          type );
        else
          PRINT_ERROR( "Through vertices not compatible with %s curve given\n",
          type);
        return CUBIT_FAILURE;
      }
    }
    if( end_vertex_ptr )
    {
      if( is_vertex_on_side( end_vertex_ptr, 0 ) || 
          is_vertex_on_side( end_vertex_ptr, 2 ) )
        ;
      else
      {
        PRINT_ERROR( "Through vertices not compatible with %s curve given\n" );
        return CUBIT_FAILURE;
      }
    }
  }
  return CUBIT_SUCCESS;  
}

CoEdge *
SplitSurfaceTool::prev_co_edge( CoEdge *co_edge_ptr )
{
  outerCoEdgeLoop.move_to( co_edge_ptr );
  return outerCoEdgeLoop.prev();
}

double
SplitSurfaceTool::get_side_length( int side )
{
  int i;
  RefEdge *ref_edge_ptr;
  double length = 0.0;
  outerCoEdgeLoop.reset();
  switch( side )
  {
  case 0:
    for( i=sideInterval[0]; i--; )
    {
      ref_edge_ptr = outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr();
      length += ref_edge_ptr->measure();
    }
    break;

  case 1:
    outerCoEdgeLoop.step( sideInterval[0] );
    for( i=sideInterval[1]; i--; )
    {
      ref_edge_ptr = outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr();
      length += ref_edge_ptr->measure();
    }
    break;

  case 2:

    outerCoEdgeLoop.step( sideInterval[0]+sideInterval[1] );
    for( i=sideInterval[2]; i--; )
    {
      ref_edge_ptr = outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr();
      length += ref_edge_ptr->measure();
    }
    break;

  case 3:

    outerCoEdgeLoop.step( sideInterval[0]+sideInterval[1]+sideInterval[2] );
    for( i=sideInterval[3]; i--; )
    {
      ref_edge_ptr = outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr();
      length += ref_edge_ptr->measure();
    }
    break;
  }

  return length;
}

CubitStatus
SplitSurfaceTool::reorient_loop( int start_offset )
{
  if( start_offset == 0 )
    return CUBIT_SUCCESS;

  DLIList<CoEdge*> outer_co_edge_loop( outerCoEdgeLoop.size() );

  outerCoEdgeLoop.reset();
  outerCoEdgeLoop.step( start_offset );

  int i;
  for( i=outerCoEdgeLoop.size(); i--; )
    outer_co_edge_loop.append( outerCoEdgeLoop.get_and_step() );

  outerCoEdgeLoop.clean_out();

  outer_co_edge_loop.reset();

  outerCoEdgeLoop = outer_co_edge_loop;

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::shuffle_corners_forward()
{
  // Simply move loops forward by sideInterval[0] (start at corner 1)
  reorient_loop( sideInterval[0] );
  
  // Shuffle sideInterval and cornerCoEdge too
  int side_interval_temp = sideInterval[0];
  sideInterval[0] = sideInterval[1];
  sideInterval[1] = sideInterval[2];
  sideInterval[2] = sideInterval[3];
  sideInterval[3] = side_interval_temp;

  CoEdge *corner_coedge_temp = cornerCoEdge[0];
  cornerCoEdge[0] = cornerCoEdge[1];
  cornerCoEdge[1] = cornerCoEdge[2];
  cornerCoEdge[2] = cornerCoEdge[3];
  cornerCoEdge[3] = corner_coedge_temp;
  
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::shuffle_zero_side_forward()
{
  if( sideInterval[0] == 0 )
  {
    cornerCoEdge[1] = cornerCoEdge[2];
    sideInterval[0] = sideInterval[1];
    sideInterval[1] = 0;
  }
  else if( sideInterval[1] == 0 )
  {
    cornerCoEdge[2] = cornerCoEdge[3];
    sideInterval[1] = sideInterval[2];
    sideInterval[2] = 0;
  }
  else if( sideInterval[2] == 0 )
  {
    cornerCoEdge[3] = cornerCoEdge[0];
    sideInterval[2] = sideInterval[3];
    sideInterval[3] = 0;
  }
  else if( sideInterval[3] == 0 )
  {
    cornerCoEdge[0] = cornerCoEdge[1];
    sideInterval[3] = sideInterval[0];
    sideInterval[0] = 0;
  }
  return CUBIT_SUCCESS;
}

CubitBoolean
SplitSurfaceTool::is_vertex_on_side( RefVertex *ref_vertex_ptr, int side )
{
  // Note this will return CUBIT_TRUE if the input vertex is on one of the
  // specified corners.

  // Position outerCoEdge loop to start of corner (note switch statement
  // falls through - case 3 will step 3 times).
  outerCoEdgeLoop.reset();
  switch( side )
  {
  case 3:
    outerCoEdgeLoop.step( sideInterval[2] );
  case 2:
    outerCoEdgeLoop.step( sideInterval[1] );
  case 1:
    outerCoEdgeLoop.step( sideInterval[0] );
  }

  CoEdge *co_edge_ptr = outerCoEdgeLoop.get();
  CubitVector ref_coords = ref_vertex_ptr->coordinates();

  if( sideInterval[side] == 0 )
  {
    // Compare coordinates
    RefVertex *start_vertex_ptr = start_vertex( co_edge_ptr );
    if( start_vertex_ptr == ref_vertex_ptr )
      return CUBIT_TRUE;
    if( ref_coords.about_equal( start_vertex_ptr->coordinates() ) )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }

  // Check if it is on each coedge
  CubitPointContainment pnt_containment;
  int i;
  for( i=sideInterval[side]; i--; )
  {
    co_edge_ptr = outerCoEdgeLoop.get_and_step();
    pnt_containment = co_edge_ptr->get_ref_edge_ptr()->
      point_containment( ref_coords );
    if( pnt_containment == CUBIT_PNT_ON )
      return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

CubitBoolean
SplitSurfaceTool::is_vertex_in_surface( RefVertex *ref_vertex_ptr, 
                                        RefFace *ref_face_ptr )
{
  // Base on finding RefFaces since there will typically be fewer RefFaces
  // to search than RefVertices
  DLIList<RefFace*> ref_face_list;
  ref_vertex_ptr->ref_faces( ref_face_list );
  return ref_face_list.is_in_list( ref_face_ptr );
}

CubitBoolean
SplitSurfaceTool::is_curve_in_outer_loop( RefEdge *ref_edge_ptr )
{
  outerCoEdgeLoop.reset();
  int i;
  for( i=outerCoEdgeLoop.size(); i--; )
  {
    if( ref_edge_ptr == outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr() )
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitBoolean
SplitSurfaceTool::is_curve_on_side( RefEdge *ref_edge_ptr, int side )
{
  if( sideInterval[side] == 0 )
    return CUBIT_FALSE;

  // Position outerCoEdge loop to start of corner (note switch statement
  // falls through - case 3 will step 3 times).
  outerCoEdgeLoop.reset();
  switch( side )
  {
  case 3:
    outerCoEdgeLoop.step( sideInterval[2] );
  case 2:
    outerCoEdgeLoop.step( sideInterval[1] );
  case 1:
    outerCoEdgeLoop.step( sideInterval[0] );
  }

  // Compare against each coedge on the given side
  int i;
  for( i=sideInterval[side]; i--; )
  {
    if( ref_edge_ptr == outerCoEdgeLoop.get_and_step()->get_ref_edge_ptr() )
      return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

CubitStatus
SplitSurfaceTool::position_co_edge_list( int i, DLIList<CoEdge*> &co_edge_list )
{
  //    6    5   4
  //    +----+---+ 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //   7+----+---+3 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //   7+-+---+--+3 
  //    |        | 
  //    |        | 
  //   7+        +3 
  //    |        | 
  //    +---+----+ 
  //    0   1    2 
  // Now, get to the start of the loop.  This will either be on vertex 0 or
  // the last vertex 7 (see diagram above).

  CoEdge *co_edge_ptr;
  CoEdge *next_co_edge_ptr;
  co_edge_list.reset();
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    next_co_edge_ptr = co_edge_list.get_and_back();
    
    // Get the tooldata from the start vertex
    TDSplitSurface *tdssv_start = (TDSplitSurface *)co_edge_ptr->
      get_TD(&TDSplitSurface::is_split_surface);
    
    if( tdssv_start )
    {
      if( tdssv_start->get_type() == 0 )
        break; // This is the start
      else if (tdssv_start->get_type() == 7)
      {
        // Check end of curve - if not a 7 or 0 we are done
        TDSplitSurface *tdssv_end = (TDSplitSurface *)next_co_edge_ptr->
          get_TD(&TDSplitSurface::is_split_surface);
        if( !tdssv_end )
          break;
        if( tdssv_end->get_type() != 7 && tdssv_end->get_type() != 0 )
          break;
      }
      else if (tdssv_start->get_type() == 6)
      {
        // Check end of curve - if not a 7 or 0 we are done
        //  (it will typically be NULL or type 3 - see top surface
        //   in diagram below)
        //
        //    6   5    4    
        //    +---+----+
        //    | \      |
        //    |   \    |
        //    |     +  |
        //    |       \|
        //   7+-+---+--+3 
        //    |        | 
        //    |        | 
        //   7+        +3 
        //    |        | 
        //    +---+----+ 
        //    0   1    2 
        
        TDSplitSurface *tdssv_end = (TDSplitSurface *)next_co_edge_ptr->
          get_TD(&TDSplitSurface::is_split_surface);
        if( !tdssv_end )
          break;
        if( tdssv_end->get_type() != 7 && tdssv_end->get_type() != 0 )
          break;
      }
    }
    co_edge_list.step();
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_a_coedges( DLIList<CoEdge*> &co_edge_list, 
                                 DLIList<CoEdge*> &a_coedges )
{
  // Get all curves until end vertex type = 2 or 3 or 4 (3 or 4 in case of triangle)
  int i;
  CoEdge *co_edge_ptr;
  CoEdge *next_co_edge_ptr;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    a_coedges.append( co_edge_ptr );

    next_co_edge_ptr = co_edge_list.get();
    
    // Get the tooldata from the next coedge (need value on end vertex of 
    // first coedge).
    TDSplitSurface *tdss = (TDSplitSurface *)next_co_edge_ptr->
      get_TD(&TDSplitSurface::is_split_surface);

    if( !tdss )
      continue;

    if( tdss->get_type() == 2 || tdss->get_type() == 3 
        || tdss->get_type() == 4 )
      break;
  }
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_b_coedges( DLIList<CoEdge*> &co_edge_list, 
                                 DLIList<CoEdge*> &b_coedges )
{
  // Keep getting curves as long as end vertex type = 3 or 4 (stop at 4)
  int i;
  CoEdge *co_edge_ptr;
  CoEdge *next_co_edge_ptr;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    next_co_edge_ptr = co_edge_list.get_and_back();
    
    // Get the tooldata from the end vertex
    TDSplitSurface *tdss = (TDSplitSurface *)next_co_edge_ptr->
      get_TD(&TDSplitSurface::is_split_surface);

    if( !tdss )
      return CUBIT_SUCCESS;

    if( tdss->get_type() == 4 )
    {
      b_coedges.append( co_edge_ptr );
      co_edge_list.step();
      return CUBIT_SUCCESS;
    }

    if( tdss->get_type() == 3 )
    {
      b_coedges.append( co_edge_ptr );
      co_edge_list.step();
      continue;
    }

    return CUBIT_SUCCESS;
  }
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_c_coedges( DLIList<CoEdge*> &co_edge_list, 
                                 DLIList<CoEdge*> &c_coedges )
{
  // Keep getting curves until end vertex type = 6 or 7 or 0 (for triangle)
  int i;
  CoEdge *co_edge_ptr;
  CoEdge *next_co_edge_ptr;
  for( i=co_edge_list.size(); i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    next_co_edge_ptr = co_edge_list.get_and_back();
    
    // Get the tooldata from the end vertex
    TDSplitSurface *tdss = (TDSplitSurface *)next_co_edge_ptr->
      get_TD(&TDSplitSurface::is_split_surface);

    if( !tdss )
    {
      c_coedges.append( co_edge_ptr );
      co_edge_list.step();
      continue;
    }

    if( tdss->get_type() == 6 || tdss->get_type() == 7 ||
        tdss->get_type() == 0 )
    {
      c_coedges.append( co_edge_ptr );
      co_edge_list.step();
      return CUBIT_SUCCESS;
    }
    
    c_coedges.append( co_edge_ptr );
    co_edge_list.step();

  }
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_d_coedges( DLIList<CoEdge*> &co_edge_list, 
                                int num_so_far, DLIList<CoEdge*> &d_coedges )
{
  // Get remaining curves
  if( co_edge_list.size() == num_so_far )
    return CUBIT_SUCCESS;

  if( num_so_far > co_edge_list.size() )
  {
    PRINT_ERROR( "Unexpected error in algorithm; aborting.\n" );
    PRINT_DEBUG_154( "       Surface = %d, num_so_far = %d\n", 
      co_edge_list.get()->get_ref_face()->id(), num_so_far );
    return CUBIT_FAILURE;
  }

  int i;
  CoEdge *co_edge_ptr;
  for( i=co_edge_list.size()-num_so_far; i--; )
  {
    co_edge_ptr = co_edge_list.get_and_step();
    d_coedges.append( co_edge_ptr );
  }

  return CUBIT_SUCCESS;
}

void
SplitSurfaceTool::list_sides_debug()
{
  if( !DEBUG_FLAG(154) )
    return;

  int i, j;
  RefFace *ref_face_ptr;
  TDSplitSurface *tdss;

  refFaceChain.reset();
  for( i=refFaceChain.size(); i--; )
  {
    ref_face_ptr = refFaceChain.get_and_step();
    tdss = (TDSplitSurface *)ref_face_ptr->
      get_TD(&TDSplitSurface::is_split_surface);
    
    PRINT_INFO( "Surface %d:\n", ref_face_ptr->id() );
    
    DLIList<CoEdge*> *side_co_edge_list;
    side_co_edge_list = tdss->get_a_coedges();
    DLIList<RefEdge*> side_ref_edge_list;
    side_co_edge_list->reset();
    for( j=side_co_edge_list->size(); j--; )
      side_ref_edge_list.append( side_co_edge_list->get_and_step()->get_ref_edge_ptr() );
    
    DLIList<CubitEntity*> cubit_edges;
    CAST_LIST(side_ref_edge_list, cubit_edges, CubitEntity);
    CubitUtil::list_entity_ids( "Side 0: ", 
      cubit_edges, 80, "\n", CUBIT_FALSE );
    
    side_co_edge_list = tdss->get_b_coedges();
    side_ref_edge_list.clean_out();
    side_co_edge_list->reset();
    for( j=side_co_edge_list->size(); j--; )
      side_ref_edge_list.append( side_co_edge_list->get_and_step()->get_ref_edge_ptr() );
    
    CAST_LIST(side_ref_edge_list, cubit_edges, CubitEntity);
    CubitUtil::list_entity_ids( "Side 1: ", 
      cubit_edges, 80, "\n", CUBIT_FALSE );
    
    side_co_edge_list = tdss->get_c_coedges();
    side_ref_edge_list.clean_out();
    side_co_edge_list->reset();
    for( j=side_co_edge_list->size(); j--; )
      side_ref_edge_list.append( side_co_edge_list->get_and_step()->get_ref_edge_ptr() );
    
    CAST_LIST(side_ref_edge_list, cubit_edges, CubitEntity);
    CubitUtil::list_entity_ids( "Side 2: ", 
      cubit_edges, 80, "\n", CUBIT_FALSE );
    
    side_co_edge_list = tdss->get_d_coedges();
    side_ref_edge_list.clean_out();
    side_co_edge_list->reset();
    for( j=side_co_edge_list->size(); j--; )
      side_ref_edge_list.append( side_co_edge_list->get_and_step()->get_ref_edge_ptr() );
    
    CAST_LIST(side_ref_edge_list, cubit_edges, CubitEntity);
    CubitUtil::list_entity_ids( "Side 3: ", 
      cubit_edges, 80, "\n", CUBIT_FALSE );
  }
}

CubitStatus
SplitSurfaceTool::find_spline_curves( RefFace *ref_face_ptr, int num_segs, 
                                      double distance, 
                                      DLIList<Curve*> *curve_list_ptr,
                                      double tolerance,
                                      CubitBoolean parametric_flg,
                                      CubitBoolean preview_flg,
                                      CubitBoolean create_ref_edges_flg )
{
  int i, j;

  TDSplitSurface *tdss = (TDSplitSurface *)ref_face_ptr->
    get_TD(&TDSplitSurface::is_split_surface);

  TopologyBridge* bridge = 0;
//   GeometryModifyEngine* gme = 
  GeometryModifyTool::instance()->get_engine( ref_face_ptr, &bridge );
  Surface* surf_ptr = dynamic_cast<Surface*>(bridge);

  // Interpolate to find the spline points
  //
  //                          sideC
  //                        a[nr-1][i]
  //          _________________________________________
  //    (0,1)|                                         |(1,1)
  //         |                                         |
  //         |                                         |
  //         |                                         |
  //         |                                         |
  //         |                                         |
  // a[j][0] |                                         |
  //         |                                         | a[j][nc-1] 
  //  sideD  |                                         |
  //         | ada (row)                               | sideB
  //         | ^                                       |
  //         | |                                       |
  //         | | (tse,ada)                             |
  //         | |                                       |
  //         | +----->tse (column)                     |
  //    (0,0)|_________________________________________|(1,0)
  //                           a[0][i]    
  //                            sideA
  
  // Split direction is vertical in the above diagram.
  
  if( ref_face_ptr->is_parametric() && parametric_flg==CUBIT_TRUE )
  {
    PRINT_DEBUG_154( "Using 2D mapping to find interior points\n" );
    // Do the mapping in parametric space - we have found this is
    // more accurate (especially since we don't have the benefit of
    // any elemental smoothing when we are done).  However, we have
    // also found that occasionally, especially on conic surfaces,
    // the mapping will create points in the wrong space of the 
    // surface (i.e., on the "other side".).  For this reason we
    // default to 3D mapping.
    // Allocate a matrix of CubitVector pointers
    int nr = tdss->coord_list_size_b(); // number of rows
    int nc = num_segs+1;                // number of columns
    Cubit2DPoint ***coords;
    coords = new Cubit2DPoint **[nr];
    
    for( j=0; j<nr; j++ )
      coords[j] = new Cubit2DPoint *[nc];
    
    // Initialize      
    for( i=0; i<nr; i++ )
    {
      for( j=0; j<nc; j++ )
        coords[i][j] = NULL;
    }
    
    // Fill the boundary coordinates into coords
    fill_boundary_coords( tdss, nr, nc, coords );
    
    // Fill the interior coords
    fill_interior_coords( tdss, nr, nc, coords );
    
    // Generate the 3D vectors
    CubitVector vec;

    for( j=1; j<nc-1; j++ )
    {
      DLIList<CubitVector*> spline_points;
      for( i=0; i<nr; i++ )
      {
        vec = ref_face_ptr->position_from_u_v( coords[i][j]->x(), coords[i][j]->y() );
        spline_points.append( new CubitVector( vec ) );
      }

      // Create the curve from the vectors
      CubitBoolean project_curve = CUBIT_TRUE;
      if( preview_flg == CUBIT_TRUE && create_ref_edges_flg == CUBIT_FALSE )
        project_curve = CUBIT_FALSE;
      Curve *curve_ptr = create_curve( spline_points, surf_ptr, tolerance, 
                                       CUBIT_TRUE, preview_flg, project_curve );      

      // Free spline points
      while( spline_points.size() ) delete spline_points.pop();
      
      if( curve_ptr == NULL )
      {
        // Just fail if a curve can't be created
        PRINT_ERROR( "Unable to create curve on Surface %d\n", ref_face_ptr->id() );
        while( curve_list_ptr->size() ) delete curve_list_ptr->pop();
        return CUBIT_FAILURE;
      }
      else
      {
        curve_list_ptr->append( curve_ptr );
      }
    }
    
    // Free memory - all of the Cubit2DPoints were allocated      
    for( i=0; i<nr; i++ )
    {
      for( j=0; j<nc; j++ )
        delete coords[i][j];
    }
    
    // Free matrix memory
    for( i=0; i<nr; i++ )
      delete coords[i];
    delete coords;
  }
  else
  {
    // Need to do the mapping in 3D space and project back to surface
    PRINT_DEBUG_154( "Using 3D mapping to find interior points\n" );
    
    // Allocate a matrix of CubitVector pointers
    int nr = tdss->coord_list_size_b(); // number of rows
    int nc = num_segs+1;                // number of columns
    CubitVector ***coords;
    coords = new CubitVector **[nr];
    
    for( i=0; i<nr; i++ )
      coords[i] = new CubitVector *[nc];
    
    // Initialize      
    for( i=0; i<nr; i++ )
    {
      for( j=0; j<nc; j++ )
        coords[i][j] = NULL;
    }
    
    // Fill the boundary coordinates into coords
    fill_boundary_coords( tdss, nr, nc, coords );
    
    // Fill the interior coords
    fill_interior_coords( tdss, nr, nc, coords );

    // Smooth the 3D points (note - don't check result - the function only gives
    // a warning if it fails).
    smooth_interior_coords( ref_face_ptr, tdss, tolerance, nr, nc, distance, coords );
    
    // Create split curves
    for( j=1; j<nc-1; j++ )
    {
      DLIList<CubitVector*> spline_points;
      for( i=0; i<nr; i++ )
      {
        ref_face_ptr->move_to_surface( *(coords[i][j]) );
        spline_points.append( new CubitVector( *(coords[i][j]) ) );
      }

      // Create the curve from the vectors
      CubitBoolean project_curve = CUBIT_TRUE;
      if( preview_flg == CUBIT_TRUE && create_ref_edges_flg == CUBIT_FALSE )
        project_curve = CUBIT_FALSE;
      Curve *curve_ptr = create_curve( spline_points, surf_ptr, tolerance,
        CUBIT_TRUE, preview_flg, project_curve ); 

      // Free spline points
      while( spline_points.size() ) delete spline_points.pop();
      
      if( curve_ptr == NULL )
      {
        // Just fail if a curve can't be created
        PRINT_ERROR( "Unable to create curve on Surface %d\n", ref_face_ptr->id() );
        while( curve_list_ptr->size() ) delete curve_list_ptr->pop();
        return CUBIT_FAILURE;
      }
      else
        curve_list_ptr->append( curve_ptr );
    }
    
    // Free interior points (the boundary locations were not allocated for 
    // the matrix)     
    for( i=1; i<nr-1; i++ )
    {
      for( j=1; j<nc-1; j++ )
        delete coords[i][j];
    }
    
    // Free matrix memory
    for( i=0; i<nr; i++ )
      delete coords[i];
    delete coords;
  }

  return CUBIT_SUCCESS; 
}

CubitStatus
SplitSurfaceTool::fill_boundary_coords( TDSplitSurface *tdss, int nr, 
                                        int nc, CubitVector ***coords )
{
  //                                sideC
  //                            coords[nr-1][c]
  //               _________________________________________
  //         (0,1)|                                         |(1,1)
  //              |                                         |
  //              |                                         |
  //              |                                         |
  //              |                                         |
  //              |                                         |
  // coords[r][0] |                                         |
  //              |                                         | coords[r][nc-1] 
  //       sideD  |                                         |
  //              | ada (row)                               | sideB
  //              | ^                                       |
  //              | |                                       |
  //              | | (tse,ada)                             |
  //              | |                                       |
  //              | +----->tse (column)                     |
  //         (0,0)|_________________________________________|(1,0)
  //                            coords[0][c]    
  //                               sideA
  
  // Split direction is vertical in the above diagram.

  // Debug (with DEBUG_FLAG_154)
  draw_boundary_coords( tdss );

  int r, c;

  // Note: sideB and sideD contain the corner coords (sideA and sideC just 
  //       contain the interior coords)
  tdss->coord_list_reset_a();
  tdss->coord_list_reset_b();
  tdss->coord_list_last_c(); // Will need to be traversed backwards
  tdss->coord_list_last_d(); // Will need to be traversed backwards

  for( r=0; r<nr; r++ )
  {
    // Rows go up B and D (note these contain the corner coords)
    coords[r][0] = tdss->coord_list_get_and_back_d();
    coords[r][nc-1] = tdss->coord_list_get_and_step_b();
  }

  for( c=1; c<nc-1; c++ )
  {
    // Columns go up A and C (note these DO NOT contain the corner coords)
    coords[0][c] = tdss->coord_list_get_and_step_a();
    coords[nr-1][c] = tdss->coord_list_get_and_back_c();
  }

  return CUBIT_SUCCESS;
}

// Debug function
CubitStatus
SplitSurfaceTool::draw_boundary_coords( TDSplitSurface *tdss )
{
  if( !DEBUG_FLAG(154) )
    return CUBIT_SUCCESS;

  int i;
  CubitVector *vec_ptr;
  CoEdge *co_edge_ptr;
  DLIList<CoEdge*> *co_edge_list_ptr;
  DLIList<RefEdge*> curve_list;
  DLIList<CubitEntity*> cubit_curves;

  // Note: sideB and sideD contain the corner coords (sideA and sideC just 
  //       contain the interior coords)

  PRINT_INFO( "Surface %d, Side A, %d coords, YELLOW\n", 
    tdss->ref_face_ptr()->id(), tdss->coord_list_size_a() );
  
  co_edge_list_ptr = tdss->get_a_coedges();
  co_edge_list_ptr->reset();
  for( i=co_edge_list_ptr->size(); i--; )
  {
    co_edge_ptr = co_edge_list_ptr->get_and_step();
    curve_list.append( co_edge_ptr->get_ref_edge_ptr() );
  }
  CAST_LIST(curve_list, cubit_curves, CubitEntity);
  CubitUtil::list_entity_ids( " Curves: ", 
    cubit_curves, 80, "\n", CUBIT_FALSE );

  tdss->coord_list_reset_a();
  for( i=tdss->coord_list_size_a(); i--; )
  {
    vec_ptr = tdss->coord_list_get_and_step_a();
    PRINT_DEBUG_100( " create vertex %f %f %f color yellow\n", vec_ptr->x(),
      vec_ptr->y(), vec_ptr->z() );
    GfxDebug::draw_point( *vec_ptr, CUBIT_YELLOW );
  }

  PRINT_INFO( "Surface %d, Side B, %d coords, MAGENTA\n", 
    tdss->ref_face_ptr()->id(), tdss->coord_list_size_b() );

  co_edge_list_ptr = tdss->get_b_coedges();
  co_edge_list_ptr->reset();
  curve_list.clean_out();
  for( i=co_edge_list_ptr->size(); i--; )
  {
    co_edge_ptr = co_edge_list_ptr->get_and_step();
    curve_list.append( co_edge_ptr->get_ref_edge_ptr() );
  }
  cubit_curves.clean_out();
  CAST_LIST(curve_list, cubit_curves, CubitEntity);
  CubitUtil::list_entity_ids( " Curves: ", 
    cubit_curves, 80, "\n", CUBIT_FALSE );

  tdss->coord_list_reset_b();
  for( i=tdss->coord_list_size_b(); i--; )
  {
    vec_ptr = tdss->coord_list_get_and_step_b();
    PRINT_DEBUG_100( " create vertex %f %f %f color magenta\n", vec_ptr->x(),
      vec_ptr->y(), vec_ptr->z() );
    GfxDebug::draw_point( *vec_ptr, CUBIT_MAGENTA );
  }

  PRINT_INFO( "Surface %d, Side C, %d coords, GREEN\n", 
    tdss->ref_face_ptr()->id(), tdss->coord_list_size_c() );

  co_edge_list_ptr = tdss->get_c_coedges();
  co_edge_list_ptr->reset();
  curve_list.clean_out();
  for( i=co_edge_list_ptr->size(); i--; )
  {
    co_edge_ptr = co_edge_list_ptr->get_and_step();
    curve_list.append( co_edge_ptr->get_ref_edge_ptr() );
  }
  cubit_curves.clean_out();
  CAST_LIST(curve_list, cubit_curves, CubitEntity);
  CubitUtil::list_entity_ids( " Curves: ", 
    cubit_curves, 80, "\n", CUBIT_FALSE );

  tdss->coord_list_reset_c();
  for( i=tdss->coord_list_size_c(); i--; )
  {
    vec_ptr = tdss->coord_list_get_and_step_c();
    PRINT_DEBUG_100( " create vertex %f %f %f color green\n", vec_ptr->x(),
      vec_ptr->y(), vec_ptr->z() );
    GfxDebug::draw_point( *vec_ptr, CUBIT_GREEN );
  }

  PRINT_INFO( "Surface %d, Side D, %d coords, RED\n", 
    tdss->ref_face_ptr()->id(), tdss->coord_list_size_d() );

  co_edge_list_ptr = tdss->get_d_coedges();
  co_edge_list_ptr->reset();
  curve_list.clean_out();
  for( i=co_edge_list_ptr->size(); i--; )
  {
    co_edge_ptr = co_edge_list_ptr->get_and_step();
    curve_list.append( co_edge_ptr->get_ref_edge_ptr() );
  }
  cubit_curves.clean_out();
  CAST_LIST(curve_list, cubit_curves, CubitEntity);
  CubitUtil::list_entity_ids( " Curves: ", 
    cubit_curves, 80, "\n", CUBIT_FALSE );

  tdss->coord_list_reset_d();
  for( i=tdss->coord_list_size_d(); i--; )
  {
    vec_ptr = tdss->coord_list_get_and_step_d();
    PRINT_DEBUG_100( " create vertex %f %f %f color red\n", vec_ptr->x(),
      vec_ptr->y(), vec_ptr->z() );
    GfxDebug::draw_point( *vec_ptr, CUBIT_RED );
  }
  
  GfxDebug::flush();

  return CUBIT_SUCCESS;
}

// Cloned from MapToolSupport, except MapToolSupport does not accurately
// calculate ada and tse
CubitStatus
SplitSurfaceTool::fill_interior_coords( TDSplitSurface *tdss,
                                        int nr, int nc, CubitVector ***coords )
{
  int r, c;

  int ada_ints = nr-1;
  int tse_ints = nc-1;

  double ada, tse;

  DLIList<double> tse_array;
  if( get_tse_array( tdss, tse_ints, tse_array ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  
  // Create all the new coords in the interior of the array
  if( tdss )
  {
    // Initialize list position

    // Must use side that is not collapsed
    if( tdss->is_b_collapsed() && tdss->is_d_collapsed() )
      ; // Nothing to do
    else if( tdss->is_b_collapsed() == CUBIT_FALSE )
    {
      tdss->param_list_reset_b();
      ada = tdss->param_list_get_and_step_b();
    }
    else
    {
      // Must go backwards on this list
      tdss->param_list_last_d();
      ada = tdss->param_list_get_and_back_d();
    }
  }
  for( r=1; r<ada_ints; r++ )
  {
    if( tdss )
    {
      if( tdss->is_b_collapsed() && tdss->is_d_collapsed() )
        ada = (double)r/(double)ada_ints;
      else if( tdss->is_b_collapsed() )
        ada = tdss->param_list_get_and_back_d()/tdss->length_d();
      else
        ada = tdss->param_list_get_and_step_b()/tdss->length_b();
    }
    else
    {
      ada = (double)r/(double)ada_ints;
    }

    tse_array.reset();
    for( c=1; c<tse_ints; c++ )
    {
      tse = tse_array.get_and_step();
      
      // Create a new coord at the mapped location on the face
      coords[r][c] = make_interior_coord( coords, nr, nc, ada, tse, r, c );
    }  
  }   
  
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::smooth_interior_coords( RefFace *ref_face_ptr,
                                          TDSplitSurface *tdss,
                                          double tolerance,
                                          int nr, int nc, double distance,
                                          CubitVector ***coords )
{
  int r, c; // Row, Column
  CubitVector *vec_ptr;

  // Preserve original coordinates and restore to those if there is an error 
  // (just give warning).
  DLIList<CubitVector*> backup_coords( nr*nc );
  for( r=0; r<nr; r++ )
  {
    for( c=0; c<nc; c++ )
      backup_coords.append( new CubitVector( *coords[r][c] ) );
  }

  // Problem: after doing a 3D map, when the interior coordinates are projected
  // back to the surface they may not lie at the correct percentage across the
  // surface (thus a fillet may not be split exactly in the middle).

  // Solution: create spline curves on the surface (perpendicular to the split
  // direction) at each coordinate.  Move the coordinate to the proper fraction
  // distance along the spline curve.  The trick is to find the proper fraction
  // distance - we use a mapping concept to determine the fraction.  The real 
  // need for the "fraction map" is that point A and point B in the diagram
  // below can be at different percentages along sideA and sideC, so the 
  // internal fractions need to be "smoothed out".  An alternative to this
  // "spline" technique might be to do something like real Winslow smoothing,
  // but that would require a lot more work to implement...

  // Setup fraction "map" (actual surface split is vertical in diagram below) -
  // note the splines perpendicular to the split direction).
  //
  //                           sideC            A
  //          _________________________+________+____________
  //    (0,1)|                                               |(1,1)
  //         |                                               |
  //         |                                               |
  //         +............spline......+........+.............+ 
  //         |                                               |
  //         |                                               |
  //         |                                               |
  //         |                                               | 
  //  sideD  +............spline.....+........+..............+ 
  //         |                                               |
  //         | y (row)                                       | sideB (tessellations)
  //         | ^                                             |
  //         | |                                             |
  //         | | (x,y)                                       |
  //         | |                                             |
  //         | +----->x (column)              B              |
  //    (0,0)|_____________________+________ +______________ |(1,0)
  //
  //                            sideA (split locations)

  Cubit2DPoint ***frac_coords;
  frac_coords = new Cubit2DPoint **[nr];
  
  for( r=0; r<nr; r++ )
    frac_coords[r] = new Cubit2DPoint *[nc];
  
  // Initialize      
  for( r=0; r<nr; r++ )
  {
    for( c=0; c<nc; c++ )
      frac_coords[r][c] = NULL;
  }
    
  // Fill the boundary "fractions" into frac_coords

  frac_coords[0][0] = new Cubit2DPoint( 0.0, 0.0 );
  frac_coords[0][nc-1] = new Cubit2DPoint( 1.0, 0.0 );
  frac_coords[nr-1][0] = new Cubit2DPoint( 0.0, 1.0 );
  frac_coords[nr-1][nc-1] = new Cubit2DPoint( 1.0, 1.0 );

  // Note: sideB and sideD contain the corner coords (sideA and sideC just 
  //       contain the interior coord)
  tdss->param_list_reset_a();
  tdss->param_list_reset_b();
  tdss->param_list_get_and_step_b(); // Move off the corner
  tdss->param_list_last_c(); // Will need to be traversed backwards
  tdss->param_list_last_d(); // Will need to be traversed backwards
  tdss->param_list_get_and_back_d(); // Move off the corner

  double frac_a, frac_b, frac_c, frac_d;
  
  for( r=1; r<nr-1; r++ )
  {
    // Rows go up B and D (note these contain the corner coords)
    if( tdss->is_b_collapsed() && tdss->is_d_collapsed() )
    {
      frac_b = double(r)/double(nr-1);
      frac_d = frac_b;
    }
    else if( tdss->is_d_collapsed() )
    {
      // Use B for both sides
      frac_b = tdss->param_list_get_and_step_b()/tdss->length_b();
      frac_d = frac_b;
    }
    else if( tdss->is_b_collapsed() )
    {
      // Use D for both sides
      frac_d = (tdss->length_d()-tdss->param_list_get_and_back_d())/tdss->length_d();
      frac_b = frac_d;
    }
    else
    {
      frac_b = tdss->param_list_get_and_step_b()/tdss->length_b();
      frac_d = (tdss->length_d()-tdss->param_list_get_and_back_d())/tdss->length_d();
    }

    frac_coords[r][0] = new Cubit2DPoint( 0.0, frac_d );
    frac_coords[r][nc-1] = new Cubit2DPoint( 1.0, frac_b );

//    PRINT_INFO( "frac_coords[%d][0] = %f, %f\n", r, frac_coords[r][0]->x(),
//    frac_coords[r][0]->y() );
//    PRINT_INFO( "frac_coords[%d][%d] = %f, %f\n", r, nc-1, frac_coords[r][nc-1]->x(),
//    frac_coords[r][nc-1]->y() );
  }

  for( c=1; c<nc-1; c++ )
  {
    // Columns go up A and C (note these DO NOT contain the corner coords)
    if( tdss->is_a_collapsed() && tdss->is_c_collapsed() )
    {
      frac_a = double(c)/double(nc-1);
      frac_c = frac_a;
    }
    else if( tdss->is_a_collapsed() )
    {
      // Use C for both sides
      frac_c = (tdss->length_c()-tdss->param_list_get_and_back_c())/tdss->length_c();
      frac_a = frac_c;
    }
    else if( tdss->is_c_collapsed() )
    {
      // Use A for both sides
      frac_a = tdss->param_list_get_and_step_a()/tdss->length_a();
      frac_c = frac_a;
    }
    else
    {
      frac_a = tdss->param_list_get_and_step_a()/tdss->length_a();
      frac_c = (tdss->length_c()-tdss->param_list_get_and_back_c())/tdss->length_c();
    }

    frac_coords[0][c] = new Cubit2DPoint( frac_a, 0.0 );
    frac_coords[nr-1][c] = new Cubit2DPoint( frac_c, 1.0 );

//    PRINT_INFO( "frac_coords[0][%d] = %f, %f\n", c, frac_coords[0][c]->x(), 
//    frac_coords[0][c]->y() );
//    PRINT_INFO( "frac_coords[%d][%d] = %f, %f\n", nr-1, c, frac_coords[nr-1][c]->x(),
//    frac_coords[nr-1][c]->y() );
  }

  // Fill interior coordinates
  fill_interior_coords( tdss, nr, nc, frac_coords );

  // Now create the splines

  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme = 
    GeometryModifyTool::instance()->get_engine( ref_face_ptr, &bridge );
//   Surface* surf_ptr = dynamic_cast<Surface*>(bridge);

  for( r=1; r<nr-1; r++ )
  {
    // Move interior points to surface
    for( c=1; c<nc-1; c++ )
      ref_face_ptr->move_to_surface( *(coords[r][c]) );

    // Make a spline across the surface
    DLIList<CubitVector*> spline_points;
    for( c=0; c<nc; c++ )
      spline_points.append( coords[r][c] );

    spline_points.reset();
    CubitVector *start_pnt = spline_points.get_and_back();
    CubitVector *end_pnt = spline_points.get();

    Point *start_Point = gme->make_Point( *start_pnt );
    Point *end_Point = gme->make_Point( *end_pnt );
    if( start_Point == NULL || end_Point == NULL )
    {
      PRINT_WARNING( "Unable to adjust split positions - split may be inaccurate\n" );

      // Restore ALL original points & move to surface
      backup_coords.reset();
      for( r=0; r<nr; r++ )
      {
        for( c=0; c<nc; c++ )
        {
          vec_ptr = backup_coords.get_and_step();
          ref_face_ptr->move_to_surface( *vec_ptr );
          coords[r][c]->set( *vec_ptr );
          delete vec_ptr;
        }
      }

      // Free memory
      if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
      if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
      for( r=0; r<nr; r++ )
      {
        for( c=0; c<nc; c++ )
          delete frac_coords[r][c];
      }
      for( r=0; r<nr; r++ )
        delete frac_coords[r];
      delete frac_coords;
      
      return CUBIT_FAILURE;
    }

    // Here, note it is not necessary to create the curve on the surface
    Curve *curve_ptr = gme->make_Curve( SPLINE_CURVE_TYPE, start_Point, end_Point, 
                                        spline_points );

    // Debug
    if( DEBUG_FLAG( 100 ) )
      draw_preview( curve_ptr, CUBIT_TRUE, CUBIT_RED );

    if( curve_ptr == NULL )
    {
      PRINT_WARNING( "Unable to adjust split positions - split may be inaccurate\n" );

      // Restore ALL original points & move to surface
      backup_coords.reset();
      for( r=0; r<nr; r++ )
      {
        for( c=0; c<nc; c++ )
        {
          vec_ptr = backup_coords.get_and_step();
          ref_face_ptr->move_to_surface( *vec_ptr );
          coords[r][c]->set( *vec_ptr );
          delete vec_ptr;
        }
      }

      // Free memory
      for( r=0; r<nr; r++ )
      {
        for( c=0; c<nc; c++ )
          delete frac_coords[r][c];
      }
      for( r=0; r<nr; r++ )
        delete frac_coords[r];
      delete frac_coords;
      
      return CUBIT_FAILURE;
    }

    // Space interior points along the Curve per the original fractions
    spline_points.reset();
    spline_points.step();
    double curve_length = curve_ptr->measure();
    for ( c=1; c<nc-1; c++ )
    {
      vec_ptr = spline_points.get_and_step();
      
      // Find total length of curve and multiply by the fraction
      double mv_distance;
      if( distance == -1.0 )
        mv_distance = curve_length*(frac_coords[r][c]->x());
      else
        mv_distance = distance;
      
      CubitVector curve_position;
      if( curve_ptr->point_from_arc_length ( *start_pnt, mv_distance, curve_position )
        == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to adjust split positions - split may be inaccurate\n" );

        // Restore ALL original points & move to surface
        backup_coords.reset();
        for( r=0; r<nr; r++ )
        {
          for( c=0; c<nc; c++ )
          {
            vec_ptr = backup_coords.get_and_step();
            ref_face_ptr->move_to_surface( *vec_ptr );
            coords[r][c]->set( *vec_ptr );
            delete vec_ptr;
          }
        }

        // Free memory
        for( r=0; r<nr; r++ )
        {
          for( c=0; c<nc; c++ )
            delete frac_coords[r][c];
        }
        
        // Free matrix memory
        for( r=0; r<nr; r++ )
          delete frac_coords[r];
        delete frac_coords;

        // Free spline
        curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
        return CUBIT_FAILURE;
      }

      vec_ptr->set( curve_position );
    }

    curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    //GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
  }

  // Free memory - all of the frac_coords were allocated      
  for( r=0; r<nr; r++ )
  {
    for( c=0; c<nc; c++ )
      delete frac_coords[r][c];
  }
  
  // Free matrix memory
  for( r=0; r<nr; r++ )
    delete frac_coords[r];
  delete frac_coords;

  // Free backup coords
  while( backup_coords.size() )
    delete backup_coords.pop();

  return CUBIT_SUCCESS;
}

CubitVector*
SplitSurfaceTool::make_interior_coord( CubitVector ***coords, 
                                       int nr, int nc,
                                       double ada, double tse,
                                       int r, int c )
{
  double temp_x, temp_y, temp_z;
  CubitVector *vec;
  int ada_ints = nr-1;
  int tse_ints = nc-1;

  temp_x =   (1.0 - ada)*coords[0][c]->x()
          + ada*coords[ada_ints][c]->x()
          + (1.0 - tse)*coords[r][0]->x()
          + tse*coords[r][tse_ints]->x()
          - (1.0 - tse)*(1.0 - ada)*coords[0][0]->x()
          - (1.0 - tse)*ada*coords[ada_ints][0]->x()
          - tse*(1.0 - ada)*coords[0][tse_ints]->x()
          - tse*ada*coords[ada_ints][tse_ints]->x();

  temp_y =   (1.0 - ada)*coords[0][c]->y() 
          + ada*coords[ada_ints][c]->y()
          + (1.0 - tse)*coords[r][0]->y() 
          + tse*coords[r][tse_ints]->y()
          - (1.0 - tse)*(1.0 - ada)*coords[0][0]->y()
          - (1.0 - tse)*ada*coords[ada_ints][0]->y()
          - tse*(1.0 - ada)*coords[0][tse_ints]->y()
          - tse*ada*coords[ada_ints][tse_ints]->y();

  temp_z =   (1.0 - ada)*coords[0][c]->z() 
          + ada*coords[ada_ints][c]->z()
          + (1.0 - tse)*coords[r][0]->z() 
          + tse*coords[r][tse_ints]->z()
          - (1.0 - tse)*(1.0 - ada)*coords[0][0]->z()
          - (1.0 - tse)*ada*coords[ada_ints][0]->z()
          - tse*(1.0 - ada)*coords[0][tse_ints]->z()
          - tse*ada*coords[ada_ints][tse_ints]->z();

  vec = new CubitVector( temp_x, temp_y, temp_z );

  return vec;
}

CubitStatus
SplitSurfaceTool::fill_boundary_coords( TDSplitSurface *tdss, int nr, 
                                        int nc, Cubit2DPoint ***coords )
{
  //                                sideC
  //                            coords[nr-1][c]
  //               _________________________________________
  //         (0,1)|                                         |(1,1)
  //              |                                         |
  //              |                                         |
  //              |                                         |
  //              |                                         |
  //              |                                         |
  // coords[r][0] |                                         |
  //              |                                         | coords[r][nc-1] 
  //       sideD  |                                         |
  //              | ada (row)                               | sideB
  //              | ^                                       |
  //              | |                                       |
  //              | | (tse,ada)                             |
  //              | |                                       |
  //              | +----->tse (column)                     |
  //         (0,0)|_________________________________________|(1,0)
  //                            coords[0][c]    
  //                               sideA
  
  // Split direction is vertical in the above diagram.

  // Debug (with DEBUG_FLAG_154)
  draw_boundary_coords( tdss );

  int r, c;
  RefFace *ref_face_ptr = tdss->ref_face_ptr();

  // Note: sideB and sideD contain the corner coords (sideA and sideC just 
  //       contain the interior coords)
  tdss->coord_list_reset_a();
  tdss->coord_list_reset_b();
  tdss->coord_list_last_c(); // Will need to be traversed backwards
  tdss->coord_list_last_d(); // Will need to be traversed backwards

  for( r=0; r<nr; r++ )
  {
    // Rows go up B and D (note these contain the corner coords)
    coords[r][0] = get_uv_point( ref_face_ptr, tdss->coord_list_get_and_back_d() );
    coords[r][nc-1] = get_uv_point( ref_face_ptr, tdss->coord_list_get_and_step_b() );
  }

  for( c=1; c<nc-1; c++ )
  {
    // Columns go up A and C (note these DO NOT contain the corner coords)
    coords[0][c] = get_uv_point( ref_face_ptr, tdss->coord_list_get_and_step_a() );
    coords[nr-1][c] = get_uv_point( ref_face_ptr, tdss->coord_list_get_and_back_c() );
  }

  return CUBIT_SUCCESS;
}

Cubit2DPoint *
SplitSurfaceTool::get_uv_point( RefFace *ref_face_ptr, CubitVector *vec_ptr )
{
  double u, v;
  ref_face_ptr->u_v_from_position( *vec_ptr, u, v );
  return( new Cubit2DPoint( u, v ) );
}

// Cloned from MapToolSupport, except MapToolSupport does not accurately
// calculate ada and tse
CubitStatus
SplitSurfaceTool::fill_interior_coords( TDSplitSurface *tdss,
                                        int nr, int nc, Cubit2DPoint ***coords )
{
  int r, c;

  int ada_ints = nr-1;
  int tse_ints = nc-1;

  double ada, tse;

  DLIList<double> tse_array;
  if( get_tse_array( tdss, tse_ints, tse_array ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  
  // Create all the new coords in the interior of the array
  if( tdss )
  {
    // Initialize list position

    // Must use side that is not collapsed
    if( tdss->is_b_collapsed() && tdss->is_d_collapsed() )
      ; // Nothing to do
    else if( tdss->is_b_collapsed() == CUBIT_FALSE )
    {
      tdss->param_list_reset_b();
      ada = tdss->param_list_get_and_step_b();
    }
    else
    {
      // Must go backwards on this list
      tdss->param_list_last_d();
      ada = tdss->param_list_get_and_back_d();
    }
  }
  for( r=1; r<ada_ints; r++ )
  {
    if( tdss )
    {
      if( tdss->is_b_collapsed() && tdss->is_d_collapsed() )
        ada = (double)r/(double)ada_ints;
      else if( tdss->is_b_collapsed() )
        ada = tdss->param_list_get_and_back_d()/tdss->length_d();
      else
        ada = tdss->param_list_get_and_step_b()/tdss->length_b();
    }
    else
    {
      ada = (double)r/(double)ada_ints;
    }

    tse_array.reset();
    for( c=1; c<tse_ints; c++ )
    {
      tse = tse_array.get_and_step();
      
      // Create a new coord at the mapped location on the face
      coords[r][c] = make_interior_coord( coords, nr, nc, ada, tse, r, c );
    }  
  }   
  
  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::get_tse_array( TDSplitSurface *tdss, int tse_ints,
                                 DLIList<double> &tse_array )
{
  int c;
  double tse;

  tdss->param_list_reset_a();
  tdss->param_list_last_c();

  for( c = 1; c < tse_ints; c++ )
  {
    if( tdss )
    {
      if( tdss->is_a_collapsed() && tdss->is_c_collapsed() )
      {
        tse = (double)c/(double)tse_ints;
      }
      else if( tdss->is_a_collapsed() )
      {
        // Just use C (go backwards on side C)
        tse =  tdss->param_list_get_and_back_c()/tdss->length_c();
      }
      else if( tdss->is_c_collapsed() )
      {
        // Just use A
        tse = tdss->param_list_get_and_step_a()/tdss->length_a();
      }
      else
      {
        // Take average fractional location across the surface
        tse = (tdss->param_list_get_and_step_a()/tdss->length_a() + 1.0 - 
          tdss->param_list_get_and_back_c()/tdss->length_c())/2.0;
      }      
    }
    else
    {
      tse = (double)c/(double)tse_ints;
    }
    
    tse_array.append( tse );

  }

  return CUBIT_SUCCESS;
}

Cubit2DPoint*
SplitSurfaceTool::make_interior_coord( Cubit2DPoint ***coords, 
                                       int nr, int nc,
                                       double ada, double tse,
                                       int r, int c )
{
  double temp_x, temp_y;
  Cubit2DPoint *pnt;
  int ada_ints = nr-1;
  int tse_ints = nc-1;

  temp_x =   (1.0 - ada)*coords[0][c]->x()
          + ada*coords[ada_ints][c]->x()
          + (1.0 - tse)*coords[r][0]->x()
          + tse*coords[r][tse_ints]->x()
          - (1.0 - tse)*(1.0 - ada)*coords[0][0]->x()
          - (1.0 - tse)*ada*coords[ada_ints][0]->x()
          - tse*(1.0 - ada)*coords[0][tse_ints]->x()
          - tse*ada*coords[ada_ints][tse_ints]->x();

  temp_y =   (1.0 - ada)*coords[0][c]->y() 
          + ada*coords[ada_ints][c]->y()
          + (1.0 - tse)*coords[r][0]->y() 
          + tse*coords[r][tse_ints]->y()
          - (1.0 - tse)*(1.0 - ada)*coords[0][0]->y()
          - (1.0 - tse)*ada*coords[ada_ints][0]->y()
          - tse*(1.0 - ada)*coords[0][tse_ints]->y()
          - tse*ada*coords[ada_ints][tse_ints]->y();

  pnt = new Cubit2DPoint( temp_x, temp_y );

  return pnt;
}

CubitStatus
SplitSurfaceTool::draw_preview( DLIList<Curve*> &curve_list, int color )
{
  int i;
  Curve *curve_ptr;
  curve_list.reset();  

  // clear any old previews first
  GfxPreview::clear();

  for( i=curve_list.size(); i--; )
  {
    curve_ptr = curve_list.get_and_step();
    draw_preview( curve_ptr, CUBIT_FALSE, color );
  }

  GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::draw_preview( Curve *curve_ptr, CubitBoolean flush,
                                int color )
{
  int num_points;
  CubitStatus result;
  GMem g_mem;

  // clear any old previews first
  GfxPreview::clear();
  
  // get the graphics 
  result = curve_ptr->get_geometry_query_engine()->
    get_graphics( curve_ptr, num_points, &g_mem );
  
  if (result==CUBIT_FAILURE || num_points == 0)
  {
    PRINT_WARNING("Unable to preview a curve\n" );
  }
  
  // Draw the polyline
  GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, color );
  if( flush )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus
SplitSurfaceTool::create_ref_edges( DLIList<Curve*> &curve_list )
{
  int i;
  Curve *curve_ptr;
  curve_list.reset();  
  for( i=curve_list.size(); i--; )
  {
    curve_ptr = curve_list.get_and_step();
    GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
  }
  return CUBIT_SUCCESS;
}

int
SplitSurfaceTool::number_coedges( RefEdge *ref_edge_ptr, RefFace *ref_face_ptr )
{
  DLIList<CoEdge*> co_edge_list;
  ref_edge_ptr->get_co_edges( co_edge_list, ref_face_ptr );
  return co_edge_list.size();
}

Curve *
SplitSurfaceTool::create_curve( DLIList<CubitVector*> &vec_list,
                                Surface *surf_ptr,
                                double tolerance,
                                CubitBoolean iterate,
                                CubitBoolean draw_pnts,
                                CubitBoolean project_curve )
{
  Curve *curve_ptr = NULL;

  if( vec_list.size() < 2 )
  {
    PRINT_ERROR( "Unable to create a curve from less than two locations.\n" );
    return NULL;
  }

  GeometryModifyEngine* gme = 
    GeometryModifyTool::instance()->get_engine( surf_ptr );

  if (NULL == gme)
  {
    PRINT_ERROR("No geometry modify engine available.  Unable to create split curve.");
    return NULL;
  }

  GeometryQueryEngine* gqe = surf_ptr->get_geometry_query_engine();
  // Note the user can set the following value through
  //  set geometry accuracy <val>
  double resabs = gqe->get_sme_resabs_tolerance();

  // Find start and end locations
  vec_list.reset();
  CubitVector *start_pnt = vec_list.get_and_back();
  CubitVector *end_pnt = vec_list.get();
  
  // Check to see if a simple straight line should be created.  Use 1/2 the
  // tolerance because points can deviate on either side of the straight line
  // (unless tolerance is resabs - then just use it).
  if( check_points_straight( surf_ptr, vec_list, 
                             tolerance==resabs ? resabs : 0.5*tolerance  ) )
  {
    // Create a straight line on the surface
    Point *start_Point = gme->make_Point( *start_pnt );
    Point *end_Point = gme->make_Point( *end_pnt );
    if( start_Point == NULL || end_Point == NULL )
    {
      if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
      if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
      return NULL;
    }
    curve_ptr = gme->make_Curve( start_Point, end_Point, surf_ptr );
    if( draw_pnts == CUBIT_TRUE )
    {
      // Drawing here is a bit of a hack but it saves storing the point
      // list
      draw_point( *start_pnt, CUBIT_BLUE );
      draw_point( *end_pnt, CUBIT_BLUE );
    }

    // Free memory
   // start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point );
  //  end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point );

    return curve_ptr;
  }

  // Check to see if an arc should be created.  Use 1/2 the tolerance because 
  // points can deviate on either side of the arc (unless tolerance is resabs -
  // then just use resabs).
  if( (curve_ptr = check_points_arc( surf_ptr, gme, vec_list, resabs, 
                                    tolerance==resabs ? resabs : 0.5*tolerance ) ) != NULL )
  {
    if( draw_pnts == CUBIT_TRUE )
    {
      // Drawing here is a bit of a hack but it saves storing the point
      // lists
      draw_point( *start_pnt, CUBIT_BLUE );
      draw_point( *end_pnt, CUBIT_BLUE );
    }
    return curve_ptr;
  }

  // Ok - its not a straight line or arc, so create a spline

  // Copy the vectors since we are going to add to the list
  DLIList<CubitVector*> spline_points;
  int i;
  vec_list.reset();
  for( i=vec_list.size(); i--; )
    spline_points.append( new CubitVector(*vec_list.get_and_step()) );

  int times_through_loop = 0;
  int done = 0;
  while( !done )
  {
    times_through_loop++;

    // Create a spline.  Do not create on the surface since that can be
    // computationally expensive and we are going to iterate.  Later we
    // will project the final result to the surface.  Note: the preceeding
    // statement is not true - the gme->make_Curve function would not 
    // actually create the curve on the surface - it only projects the 
    // points to the surface then makes a free spline - this is puzzling.
    Point *start_Point = gme->make_Point( *start_pnt );
    Point *end_Point = gme->make_Point( *end_pnt );
    if( start_Point == NULL || end_Point == NULL )
    {
      while( spline_points.size() ) delete spline_points.pop();
      if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
      if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
      return NULL;
    }
    curve_ptr = gme->make_Curve( SPLINE_CURVE_TYPE, start_Point, end_Point, 
                                 spline_points );

    if( curve_ptr == NULL )
    {
      while( spline_points.size() ) delete spline_points.pop();
      start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
      end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
      return NULL;
    }

    if( iterate == CUBIT_FALSE )
      break;
    
    // Check to see if spline is close enough
    CubitVector *start_vec_ptr, *end_vec_ptr;
    spline_points.reset();
    CubitVector mid_vec;
    double dist;
    int insert_num = 0;
    for( i=spline_points.size()-1; i--; )
    {
      start_vec_ptr = spline_points.get_and_step();
      end_vec_ptr = spline_points.get();
      mid_vec = (*start_vec_ptr + *end_vec_ptr)/2.0;
      
      // Check distance from mid_vec to Curve
      CubitVector mid_vec_crv;
      curve_ptr->closest_point( mid_vec, mid_vec_crv );
      
      double x2 = mid_vec.x()-mid_vec_crv.x();
      x2 = x2*x2;
      double y2 = mid_vec.y()-mid_vec_crv.y();
      y2 = y2*y2;
      double z2 = mid_vec.z()-mid_vec_crv.z();
      z2 = z2*z2;
      
      dist = sqrt( x2+y2+z2 );
      
      if( dist > tolerance )
      {
        // Insert the midpoint location into the list
        spline_points.back();
        spline_points.insert( new CubitVector( mid_vec ) );
        insert_num++;
      }
      
    }

    if( insert_num == 0 )
      done++;

    if( times_through_loop > 100 )
    {
      PRINT_WARNING( "Unable to closely approximate spline\n" );
      done++;
    }

    if( !done )
    {
       curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    }
    
    if( insert_num > 0 )
      PRINT_DEBUG_154( "Inserted %d points when creating spline\n", insert_num );
  }

  if( surf_ptr->geometry_type() != PLANE_SURFACE_TYPE && 
      project_curve == CUBIT_TRUE )
  {
    // Project the spline to the surface (if the surface is not planar)
    DLIList<Surface*> surf_list;
    surf_list.append( surf_ptr );
    DLIList<Curve*> curve_list;
    curve_list.append( curve_ptr );
    DLIList<Curve*> curve_list_new;

    if( gme->project_edges( surf_list, curve_list, curve_list_new ) == CUBIT_FAILURE )
    {
      PRINT_WARNING( "Unable to project curve to surface - split may fail\n" );
    }
    else
    {
        curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
        curve_ptr = curve_list_new.get();
    }
  }

  if( curve_ptr && draw_pnts == CUBIT_TRUE )
  {
    // Drawing here is a bit of a hack but it saves storing the point
    // lists
    draw_points( spline_points, CUBIT_BLUE, CUBIT_FALSE );
  }

  // Free memory
  while( spline_points.size() ) delete spline_points.pop();

  return curve_ptr;
}

CubitBoolean
SplitSurfaceTool::check_points_straight( Surface *surf_ptr,
                                         DLIList<CubitVector*> &point_list,
                                         double tolerance )
{
  GeometryType geom_type = surf_ptr->geometry_type();

  // Get start and end point
  point_list.reset();
  CubitVector start_pnt( *point_list.get_and_back() );
  CubitVector end_pnt( *point_list.get() );

  // If 2 or 3 points, check points along a straight line to make sure they
  // lie on the surface.  If we only have 2 points and we meet this criteria
  // we are done.
  if( point_list.size() < 4 )
  {
    if( point_list.size()==2 && geom_type == PLANE_SURFACE_TYPE )
      return CUBIT_TRUE;

    // Check several points on the straight line - see if on surface

    CubitVector check_pnt;
    check_pnt.set( start_pnt.x() + .25 * (end_pnt.x() - start_pnt.x()),
                   start_pnt.y() + .25 * (end_pnt.y() - start_pnt.y()),
                   start_pnt.z() + .25 * (end_pnt.z() - start_pnt.z()) );
    if( is_point_on_surface( surf_ptr, check_pnt, GEOMETRY_RESABS ) == CUBIT_FALSE ) 
      return CUBIT_FALSE;

    check_pnt.set( start_pnt.x() + .5 * (end_pnt.x() - start_pnt.x()),
                   start_pnt.y() + .5 * (end_pnt.y() - start_pnt.y()),
                   start_pnt.z() + .5 * (end_pnt.z() - start_pnt.z()) );
    if( is_point_on_surface( surf_ptr, check_pnt, GEOMETRY_RESABS ) == CUBIT_FALSE )
      return CUBIT_FALSE;

    check_pnt.set( start_pnt.x() + .75 * (end_pnt.x() - start_pnt.x()),
                   start_pnt.y() + .75 * (end_pnt.y() - start_pnt.y()),
                   start_pnt.z() + .75 * (end_pnt.z() - start_pnt.z()) );
    if( is_point_on_surface( surf_ptr, check_pnt, GEOMETRY_RESABS ) == CUBIT_FALSE )
      return CUBIT_FALSE;

    if( point_list.size() == 2 )
      return CUBIT_TRUE;
  }

  // Get vector along the line
  CubitVector ln_vec = end_pnt - start_pnt;
  ln_vec.normalize();
  
  int i;
  point_list.reset();
  point_list.step();
  for( i=point_list.size()-2; i--; )
  {
    CubitVector *pnt = point_list.get_and_step();
    
    // Get vector from point to line origin
    CubitVector vec = *pnt - start_pnt;
    
    // Calculate distance from point to origin
    double distpnt_orig = vec.length();
    
    // Take dot product of vec & unit vector to get distance to intersection 
    // along the line
    double dot = fabs( vec % ln_vec );
    
    // Calculate distance from point to line w/sqrt (use fabs to account for
    // tiny roundoff error that could result in a negative number)
    double distance = sqrt(fabs(distpnt_orig*distpnt_orig - dot*dot));
    
    if( distance > tolerance )
      return CUBIT_FALSE;
    
    if( geom_type != PLANE_SURFACE_TYPE )
    {      
      // Find the actual intersection point and determine if it is within
      // GEOMETRY_RESABS of the surface
      
      // Travel along the line to get the intersection point
      CubitVector int_pnt;
      start_pnt.next_point( ln_vec, dot, int_pnt );
      
      // Find the closest point on the surface
      CubitVector closest_loc_on_surf;
      if( surf_ptr->closest_point( int_pnt, &closest_loc_on_surf ) == CUBIT_FAILURE )
      {
        return CUBIT_FALSE;
      }
      
      // Find the distance between these two points
      double dist = int_pnt.distance_between( closest_loc_on_surf );
      if( dist > GEOMETRY_RESABS )
          return CUBIT_FALSE;
      }
    }

  return CUBIT_TRUE;
}

Curve *
SplitSurfaceTool:: check_points_arc( Surface *surf_ptr,
                                     GeometryModifyEngine *gme,
                                     DLIList<CubitVector*> &point_list,
                                     double resabs,
                                     double tolerance )
{
  if( point_list.size() < 2 )
    return 0;

  point_list.reset();
  CubitVector start_pnt( *point_list.get_and_back() );
  CubitVector end_pnt( *point_list.get() );

  if( point_list.size() == 2 )
    return create_arc_two( gme, surf_ptr, start_pnt, end_pnt, resabs,
                           tolerance );

  // Get a mid point location
  CubitVector mid_pnt;
  int mid = point_list.size()/2;
  point_list.reset();
  point_list.step( mid );
  mid_pnt = *point_list.get();

  // Check if 3 points are colinear
  double a = (start_pnt-mid_pnt).length();
  double b = (mid_pnt-end_pnt).length();
  double c = (end_pnt-start_pnt).length();
  
  double denom = sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
  if( denom < GEOMETRY_RESABS )
    return 0;

  // Special case if 3 points
  if( point_list.size() == 3 )
    return create_arc_three( gme, surf_ptr, start_pnt, mid_pnt, end_pnt, resabs );

  // Get a potentially more accurate arc mid_pnt
  get_arc_mid_pnt( surf_ptr, start_pnt, end_pnt, mid_pnt, tolerance );

  // Create an actual arc and do some comparisons
  double deviation;
  Curve *curve_ptr = check_if_arc( gme, surf_ptr, point_list, start_pnt,
                                   mid_pnt, end_pnt, resabs, tolerance, deviation );

  if( curve_ptr )
    return curve_ptr;

  // Try harder to find an arc - use each point as the center point
  int i;
  point_list.reset();
  point_list.step();
  CubitVector *mid_pnt_ptr;
  for( i=point_list.size()-2; i--; )
  {
    mid_pnt_ptr = point_list.get_and_step();

    curve_ptr = check_if_arc( gme, surf_ptr, point_list, start_pnt,
                              *mid_pnt_ptr, end_pnt, resabs, tolerance,
                              deviation );

    if( curve_ptr )
      return curve_ptr;

    if( deviation > tolerance )
    {
      return 0;
    }
  }

  return 0;
}

CubitStatus
SplitSurfaceTool::get_arc_mid_pnt( Surface *surf_ptr,
                                   CubitVector &start_pnt,
                                   CubitVector &end_pnt,
                                   CubitVector &mid_pnt,
                                   double tolerance )
{
  // Possible Surface types.
  //  CONE_SURFACE_TYPE
  //  PLANE_SURFACE_TYPE
  //  SPHERE_SURFACE_TYPE
  //  SPLINE_SURFACE_TYPE
  //  TORUS_SURFACE_TYPE
  //  BEST_FIT_SURFACE_TYPE
  //  FACET_SURFACE_TYPE
  //  UNDEFINED_SURFACE_TYPE
  GeometryType geom_type = surf_ptr->geometry_type();

  // If surface is a cone, sphere or torus, we will *try* to get a more accurate
  // mid point location.  This can often result in a better arc.
  if( (geom_type == CONE_SURFACE_TYPE || 
      geom_type == SPHERE_SURFACE_TYPE ||
      geom_type == TORUS_SURFACE_TYPE) &&
      surf_ptr->is_parametric() )
  {
    // Determine if arc traverses approximately along the uv space of the surface
    // Using the surface uv space to create the arc should be the most accurate
    // method of creating the arc.

    double start_u, start_v;
    if( surf_ptr->u_v_from_position ( start_pnt, start_u, start_v )
      == CUBIT_FAILURE )
      return CUBIT_FAILURE;
    
    double end_u, end_v;
    if( surf_ptr->u_v_from_position ( end_pnt, end_u, end_v )
      == CUBIT_FAILURE )
      return CUBIT_FAILURE;
    
    // Use average.  Note this could go to the opposite side of the surface
    // - I'm not sure how to handle this other than checking the resultant 
    // point and flipping to the other side if it is off.
    double mid_u, mid_v;
    mid_u = (start_u+end_u)/2.0;
    mid_v = (start_v+end_v)/2.0;

    CubitVector new_mid_pnt = surf_ptr->position_from_u_v( mid_u, mid_v );

    CubitPointContainment pnt_on_flg;
    pnt_on_flg = surf_ptr->point_containment( mid_u, mid_v );
    if( pnt_on_flg == CUBIT_PNT_OUTSIDE ||
        pnt_on_flg == CUBIT_PNT_UNKNOWN )
    {
      CubitVector mid_pnt_ref;

      if( reflect_arc_pnt( start_pnt, new_mid_pnt, end_pnt, new_mid_pnt, mid_pnt_ref ) ==
        CUBIT_SUCCESS )
      {
        // Make sure the point is on the surface
        CubitVector mid_pnt_ref_on;
        surf_ptr->closest_point( mid_pnt_ref, &mid_pnt_ref_on );
        if( mid_pnt_ref.distance_between( mid_pnt_ref_on ) < tolerance )
        {
          pnt_on_flg = surf_ptr->point_containment( mid_pnt_ref );
          if( pnt_on_flg == CUBIT_PNT_INSIDE ||
              pnt_on_flg == CUBIT_PNT_BOUNDARY )
          {
            mid_pnt = mid_pnt_ref_on;
            return CUBIT_SUCCESS;
          }
        }
      }
    }
    else
    {
      mid_pnt = new_mid_pnt;
      return CUBIT_SUCCESS;
    }
  }

  return CUBIT_FAILURE;
}

Curve*
SplitSurfaceTool::create_arc_two( GeometryModifyEngine *gme,
                                  Surface *surf_ptr,
                                  CubitVector &start_pnt,
                                  CubitVector &end_pnt,
                                  double resabs,
                                  double tolerance )
{
  CubitVector mid_pnt;
  if( get_arc_mid_pnt( surf_ptr, start_pnt, end_pnt, mid_pnt, tolerance )
    == CUBIT_FAILURE )
    return 0;

  // Check if 3 points are colinear
  double a = (start_pnt-mid_pnt).length();
  double b = (mid_pnt-end_pnt).length();
  double c = (end_pnt-start_pnt).length();
  
  double denom = sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
  if( denom < GEOMETRY_RESABS )
    return 0;

  Point *start_Point = gme->make_Point( start_pnt );
  Point *end_Point = gme->make_Point( end_pnt );
  Point *mid_Point = gme->make_Point( mid_pnt );
  if( start_Point == NULL || mid_Point == NULL || end_Point == NULL )
  {
    if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
    if( mid_Point ) mid_Point->get_geometry_query_engine()->delete_solid_model_entities(mid_Point);
    if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
    return 0;
  }

  // Create a 3-pt arc with the GeometryModifyEngine
  Curve *curve_ptr = gme->create_arc_three( start_Point, 
                                            mid_Point,
                                            end_Point );

  if( !curve_ptr )
    return 0;

  // Check 1/4 and 3/4 location - make sure on surface
  CubitVector check_loc;
  if( curve_ptr->position_from_fraction( .25, check_loc ) == CUBIT_FAILURE )
  {
    curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  // Now make sure within GEOMETRY_RESABS to the surface
  if( is_point_on_surface( surf_ptr, check_loc, resabs ) == CUBIT_FALSE )
  {
    curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  if( curve_ptr->position_from_fraction( .75, check_loc ) == CUBIT_FAILURE )
  {
    curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  // Now make sure within GEOMETRY_RESABS to the surface
  if( is_point_on_surface( surf_ptr, check_loc, resabs ) == CUBIT_FALSE )
  {
    curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  return curve_ptr;
}

Curve*
SplitSurfaceTool::create_arc_three( GeometryModifyEngine *gme,
                                    Surface *surf_ptr,
                                    CubitVector &start_pnt,
                                    CubitVector &mid_pnt,
                                    CubitVector &end_pnt,
                                    double resabs )
{
  Point *start_Point = gme->make_Point( start_pnt );
  Point *end_Point = gme->make_Point( end_pnt );
  Point *mid_Point = gme->make_Point( mid_pnt );
  if( start_Point == NULL || mid_Point == NULL || end_Point == NULL )
  {
    if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
    if( mid_Point ) mid_Point->get_geometry_query_engine()->delete_solid_model_entities(mid_Point);
    if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
    return 0;
  }

  // Create a 3-pt arc with the GeometryModifyEngine
  Curve *curve_ptr = gme->create_arc_three( start_Point, 
                                            mid_Point,
                                            end_Point );
  if( !curve_ptr )
    return 0;

  // Check 1/4 and 3/4 location - make sure on surface
  CubitVector check_loc;
  if( curve_ptr->position_from_fraction( .25, check_loc ) == CUBIT_FAILURE )
  {
     curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  // Now make sure within GEOMETRY_RESABS to the surface
  if( is_point_on_surface( surf_ptr, check_loc, resabs ) == CUBIT_FALSE )
  {
     curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  if( curve_ptr->position_from_fraction( .75, check_loc ) == CUBIT_FAILURE )
  {
     curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  // Now make sure within GEOMETRY_RESABS to the surface
  if( is_point_on_surface( surf_ptr, check_loc, resabs ) == CUBIT_FALSE )
  {
     curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
    return 0;
  }

  return curve_ptr;
}

Curve *
SplitSurfaceTool::check_if_arc( GeometryModifyEngine *gme,
                                Surface *surf_ptr,
                                DLIList<CubitVector*> &point_list,
                                CubitVector &start_pnt,
                                CubitVector &mid_pnt,
                                CubitVector &end_pnt,
                                double resabs,
                                double tolerance,
                                double &deviation )
{
  // Check if 3 points are colinear
  double a = (start_pnt-mid_pnt).length();
  double b = (mid_pnt-end_pnt).length();
  double c = (end_pnt-start_pnt).length();
  
  double denom = sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
  if( denom < GEOMETRY_RESABS )
  {
    deviation = 99999.0;
    return 0;
  }
  
  // Create an actual arc and do some comparisons

  Point *start_Point = gme->make_Point( start_pnt );
  Point *end_Point = gme->make_Point( end_pnt );
  Point *mid_Point = gme->make_Point( mid_pnt );
  if( start_Point == NULL || mid_Point == NULL || end_Point == NULL )
  {
    if( start_Point ) start_Point->get_geometry_query_engine()->delete_solid_model_entities(start_Point);
    if( mid_Point ) mid_Point->get_geometry_query_engine()->delete_solid_model_entities(mid_Point);
    if( end_Point ) end_Point->get_geometry_query_engine()->delete_solid_model_entities(end_Point);
    return 0;
  }

  // Create a 3-pt arc with the GeometryModifyEngine
  Curve *curve_ptr = gme->create_arc_three( start_Point, 
                                            mid_Point,
                                            end_Point );
  if( curve_ptr == NULL )
    return 0;

  // Now check all points to see if they are within tolerance to the arc
  // At the same time, check to make sure that each point on the arc is
  // within GEOMETRY_RESABS of the surface (otherwise an imprint will likely
  // not occur).
  int i;
  CubitVector closest_loc;
  point_list.reset();
  double max_deviation = -1.0;
  for( i=point_list.size(); i--; )
  {
    CubitVector *vec_ptr = point_list.get_and_step();
    if( curve_ptr->closest_point( *vec_ptr, closest_loc ) == CUBIT_FAILURE )
    {
       curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
      return 0;
    }

    deviation = vec_ptr->distance_between( closest_loc );
    if( deviation > max_deviation )
      max_deviation = deviation;

    if( deviation > tolerance )
    {
       curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
      return 0;
    }

    // Now make sure within GEOMETRY_RESABS (actually resbas) to the surface
    if( is_point_on_surface( surf_ptr, closest_loc, resabs ) == CUBIT_FALSE )
    {
       curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr );
      return 0;
    }
  }

  deviation = max_deviation;
  return curve_ptr;
}

CubitStatus
SplitSurfaceTool::reflect_arc_pnt( CubitVector &vec1,
                                   CubitVector &vec2,
                                   CubitVector &vec3,
                                   CubitVector &pnt_to_reflect,
                                   CubitVector &out_pnt )

{
  double a = (vec1-vec2).length();
  double b = (vec2-vec3).length();
  double c = (vec3-vec1).length();
  
  double denom = sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
  if(denom < GEOMETRY_RESABS )
    return CUBIT_FAILURE;

  AnalyticGeometryTool *agt = AnalyticGeometryTool::instance();
  
  // Find the circle radius
  double rad = (a*b*c)/denom;
  
  // Find a coordinate system in the plane of the arc
  // Start by getting normal to arc
  CubitVector z = (vec1-vec2)*(vec2-vec3);
  z.normalize();
  CubitVector x;
  CubitVector y;
  z.orthogonal_vectors( x,y );

  // Now we are going to transform the points to a local coordinate
  // system in the plane of the arc - thats where we can work to
  // determine the equation of the circle (in 2D).
  double mtx_global_to_local[4][4];
  double mtx_local_to_global[4][4];

  double xt[3], yt[3], zt[3], origint[3];
  agt->copy_pnt( x, xt );
  agt->copy_pnt( y, yt );
  agt->copy_pnt( z, zt );
  agt->copy_pnt( vec1, origint );

  agt->vecs_to_mtx( xt, yt, zt, origint, mtx_local_to_global );
  agt->inv_trans_mtx( mtx_local_to_global, mtx_global_to_local );

  double vec1_loc[3], vec2_loc[3], vec3_loc[3];
  agt->copy_pnt( vec1, vec1_loc );
  agt->copy_pnt( vec2, vec2_loc );
  agt->copy_pnt( vec3, vec3_loc );

  agt->transform_pnt( mtx_global_to_local, vec1_loc, vec1_loc );
  agt->transform_pnt( mtx_global_to_local, vec2_loc, vec2_loc );
  agt->transform_pnt( mtx_global_to_local, vec3_loc, vec3_loc );

//  PRINT_INFO( "vec1_loc = %f, %f, %f\n", vec1_loc[0], vec1_loc[1], vec1_loc[2] );
//  PRINT_INFO( "vec2_loc = %f, %f, %f\n", vec2_loc[0], vec2_loc[1], vec2_loc[2] );
//  PRINT_INFO( "vec3_loc = %f, %f, %f\n", vec3_loc[0], vec3_loc[1], vec3_loc[2] );

  double x_1 = vec1_loc[0];
  double y_1 = vec1_loc[1];
  double x_2 = vec2_loc[0];
  double y_2 = vec2_loc[1];
  double x_3 = vec3_loc[0];
  double y_3 = vec3_loc[1];
  
  double N_1_mat[2][2];
  N_1_mat[0][0] = (x_2*x_2 + y_2*y_2-(x_1*x_1 + y_1*y_1));
  N_1_mat[0][1] = y_2-y_1;
  N_1_mat[1][0] = (x_3*x_3 + y_3*y_3-(x_1*x_1 + y_1*y_1));
  N_1_mat[1][1] = y_3-y_1;
  double N_2_mat[2][2];
  N_2_mat[0][0] = x_2-x_1;
  N_2_mat[0][1] = (x_2*x_2 + y_2*y_2-(x_1*x_1 + y_1*y_1));
  N_2_mat[1][0] = x_3-x_1;
  N_2_mat[1][1] = (x_3*x_3 + y_3*y_3-(x_1*x_1 + y_1*y_1));
  double D_mat[2][2];
  D_mat[0][0] = x_2-x_1;
  D_mat[0][1] = y_2-y_1;
  D_mat[1][0] = x_3-x_1;
  D_mat[1][1] = y_3-y_1;
  
  double N_1 = N_1_mat[0][0]*N_1_mat[1][1] - N_1_mat[0][1]*N_1_mat[1][0];
  double N_2 = N_2_mat[0][0]*N_2_mat[1][1] - N_2_mat[0][1]*N_2_mat[1][0];
  double D = D_mat[0][0]*D_mat[1][1] - D_mat[0][1]*D_mat[1][0];

  double center_loc[3];
  center_loc[0] = N_1/(2.0*D);
  center_loc[1] = N_2/(2.0*D);
  center_loc[2] = 0.0;

//  double centert[3];
//  agt->transform_pnt( mtx_local_to_global, center_loc, centert );
//  
//  CubitVector center( centert );
//  PRINT_INFO( "Arc Center = %f, %f, %f\n", center.x(), center.y(), center.z() );
//
//  GfxDebug::draw_point( center, CUBIT_BLUE );

  // Now that we have the center, it is simple matter to find the reflected
  // point from pnt_to_reflect by traversing 2*rad along a vector from
  // pnt_to_reflect to center.
  double pnt_to_reflectt[3];
  agt->copy_pnt( pnt_to_reflect, pnt_to_reflectt );
  double pnt_to_reflect_loc[3];
  agt->transform_pnt( mtx_global_to_local, pnt_to_reflectt, pnt_to_reflect_loc );
  double trav_vec[3];
  if( agt->get_vec( pnt_to_reflect_loc, center_loc, trav_vec ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  agt->unit_vec( trav_vec, trav_vec );

  double ref_pnt_loc[3];
  agt->next_pnt( pnt_to_reflect_loc, trav_vec, 2.0*rad, ref_pnt_loc );

  // Now transform this point to the global csys
  double ref_pntt[3];
  agt->transform_pnt( mtx_local_to_global, ref_pnt_loc, ref_pntt );

  out_pnt.set( ref_pntt[0], ref_pntt[1], ref_pntt[2] );

  return CUBIT_SUCCESS;
}

CubitBoolean
SplitSurfaceTool::is_point_on_surface( Surface *surf_ptr, CubitVector &pnt, 
                                       double resabs )
{
  // Make sure within GEOMETRY_RESABS to the surface
  CubitVector closest_loc_on_surf;
  if( surf_ptr->closest_point( pnt, &closest_loc_on_surf ) == CUBIT_FAILURE )
    return CUBIT_FALSE;

  double dist = pnt.distance_between( closest_loc_on_surf );
  if( dist > resabs ) // Was GEOMETRY_RESABS
    return CUBIT_FALSE;
  
  return CUBIT_TRUE;
}

int
SplitSurfaceTool::count_surfaces_in_owning_body( RefFace *ref_face_ptr, 
                                                 Body *&body_ptr )
{
  // Returns: Number of surfaces in the owning body of the RefFace
  //          -1 - RefFace is free (no owning body)
  //          -2 - RefFace is owned by more than one body
  body_ptr = 0;
  DLIList<Body*> body_list;

  ref_face_ptr->bodies( body_list );
  if( body_list.size() == 0 )
    return -1;
  else if( body_list.size() > 1 )
    return -2;
  
  body_ptr = body_list.get();
  DLIList<RefFace*> ref_face_list;
  body_ptr->ref_faces( ref_face_list );
  return ref_face_list.size();
}

int
SplitSurfaceTool::count_surfaces_in_body( Body *body_ptr )
{
  // Return number of surfaces in the body
  DLIList<RefFace*> ref_face_list;
  body_ptr->ref_faces( ref_face_list );
  return ref_face_list.size();
}

int
SplitSurfaceTool::count_curves_in_body( Body *body_ptr )
{
  // Return number of edges in the body
  DLIList<RefEdge*> ref_edge_list;
  body_ptr->ref_edges( ref_edge_list );
  return ref_edge_list.size();
}

