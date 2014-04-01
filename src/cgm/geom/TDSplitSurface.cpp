//-Class: TDSplitSurface.cpp

#include "TDSplitSurface.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "GMem.hpp"

TDSplitSurface::TDSplitSurface( int vertex_type )
{
  vertexType = vertex_type;
  sideA = NULL;
  sideB = NULL;
  sideC = NULL;
  sideD = NULL;
}

TDSplitSurface::TDSplitSurface( RefFace *ref_face_ptr )
{
  refFacePtr = ref_face_ptr;

  sideA = NULL;
  sideB = NULL;
  sideC = NULL;
  sideD = NULL;
}

TDSplitSurface::~TDSplitSurface()
{
  if( sideA )
    delete sideA;
  if( sideB )
    delete sideB;
  if( sideC )
    delete sideC;
  if( sideD )
    delete sideD;
}

CubitStatus
TDSplitSurface::add_coedges( DLIList<CoEdge*> &co_edge_list,
                             int side_interval[] )
{
  int i;
  co_edge_list.reset();
  CoEdge *co_edge_ptr;
  DLIList<CoEdge*> chain;
  RefVertex *start_vertex_ptr;

  // SIDE A
  if( side_interval[0] == 0 )
  {
    co_edge_ptr = co_edge_list.get();
    start_vertex_ptr = start_vertex( co_edge_ptr );
    CubitVector tmp(start_vertex_ptr->coordinates());
    sideA = new SSSide( refFacePtr, chain, &tmp );
  }
  else
  {
    for( i=side_interval[0]; i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      chain.append( co_edge_ptr );
    }
    sideA = new SSSide( refFacePtr, chain );
  }
  chain.clean_out();

  // SIDE B
  if( side_interval[1] == 0 )
  {
    co_edge_ptr = co_edge_list.get();
    start_vertex_ptr = start_vertex( co_edge_ptr );
    CubitVector tmp(start_vertex_ptr->coordinates());
    sideB = new SSSide( refFacePtr, chain, &tmp );
  }
  else
  {
    for( i=side_interval[1]; i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      chain.append( co_edge_ptr );
    }
    sideB = new SSSide( refFacePtr, chain );
  }
  chain.clean_out();

  // SIDE C
  if( side_interval[2] == 0 )
  {
    co_edge_ptr = co_edge_list.get();
    start_vertex_ptr = start_vertex( co_edge_ptr );
    CubitVector tmp(start_vertex_ptr->coordinates());
    sideC = new SSSide( refFacePtr, chain, &tmp );
  }
  else
  {
    for( i=side_interval[2]; i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      chain.append( co_edge_ptr );
    }
    sideC = new SSSide( refFacePtr, chain );
  }
  chain.clean_out();

  // SIDE D
  if( side_interval[3] == 0 )
  {
    co_edge_ptr = co_edge_list.get();
    start_vertex_ptr = start_vertex( co_edge_ptr );
    CubitVector tmp(start_vertex_ptr->coordinates());
    sideD = new SSSide( refFacePtr, chain, &tmp );
  }
  else
  {
    for( i=side_interval[3]; i--; )
    {
      co_edge_ptr = co_edge_list.get_and_step();
      chain.append( co_edge_ptr );
    }
    sideD = new SSSide( refFacePtr, chain );
  }

  return CUBIT_SUCCESS;
}

CubitStatus 
TDSplitSurface::add_a_coedges( DLIList<CoEdge*> &a_coedges,
                               RefVertex *start_vertex_ptr )
{
  if( start_vertex_ptr == NULL )
    sideA = new SSSide( refFacePtr, a_coedges );
  else
  {
    // Collapsed side (for triangle)
    CubitVector tmp(start_vertex_ptr->coordinates()); 
    sideA = new SSSide( refFacePtr, a_coedges, &tmp );
  }
  return CUBIT_SUCCESS;
}

CubitStatus 
TDSplitSurface::add_b_coedges( DLIList<CoEdge*> &b_coedges,
                               RefVertex *start_vertex_ptr )
{
  if( start_vertex_ptr == NULL )
    sideB = new SSSide( refFacePtr, b_coedges );
  else
  {
    // Collapsed side (for triangle)
    CubitVector tmp(start_vertex_ptr->coordinates()); 
    sideB = new SSSide( refFacePtr, b_coedges, &tmp );
  }
  return CUBIT_SUCCESS;
}

CubitStatus 
TDSplitSurface::add_c_coedges( DLIList<CoEdge*> &c_coedges,
                               RefVertex *start_vertex_ptr )
{
  if( start_vertex_ptr == NULL )
    sideC = new SSSide( refFacePtr, c_coedges );
  else
  {
    // Collapsed side (for triangle)
    CubitVector tmp(start_vertex_ptr->coordinates()); 
    sideC = new SSSide( refFacePtr, c_coedges, &tmp );
  }
  return CUBIT_SUCCESS;
}

CubitStatus 
TDSplitSurface::add_d_coedges( DLIList<CoEdge*> &d_coedges,
                               RefVertex *start_vertex_ptr )
{
  if( start_vertex_ptr == NULL )
    sideD = new SSSide( refFacePtr, d_coedges );
  else
  {
    // Collapsed side (for triangle)
    CubitVector tmp(start_vertex_ptr->coordinates()); 
    sideD = new SSSide( refFacePtr, d_coedges, &tmp );
  }
  return CUBIT_SUCCESS;
}

DLIList<CoEdge*> *
TDSplitSurface::get_a_coedges()
{
  return sideA->co_edges();
}

DLIList<CoEdge*> *
TDSplitSurface::get_b_coedges()
{
  return sideB->co_edges();
}

DLIList<CoEdge*> *
TDSplitSurface::get_c_coedges()
{
  return sideC->co_edges();
}

DLIList<CoEdge*> *
TDSplitSurface::get_d_coedges()
{
  return sideD->co_edges();
}

CubitStatus
TDSplitSurface::tessellate_sides( double tol, double fraction, double distance,
                                  int num_segs, 
                                  DLIList<RefVertex*> &through_vertex_list )
{
  // Sides B and D will use the graphics tessellation to build the 
  // param lists.
  if( sideB->build_param_list_from_facets( tol ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  if( sideD->build_param_list_from_facets( tol ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Syncronize the lists so that we have the same number and evenly
  // spaced tessellations on B and D sides
  if( sideB->syncronize_lists( sideD, tol ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to interpolate split location.\n" );
    return CUBIT_FAILURE;
  }

  // Sides A and C typically only need to retrieve the 50% location,
  // but may have a location different than 50% or multiple locations.
  // Note if num_segs>2 the fraction is ignored.

  // Populate lists for side A
  if( sideA->build_param_list( fraction, distance, num_segs, through_vertex_list ) 
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Populate lists for side C
  if( distance != -1.0 )
    distance = sideC->length()-distance;

  if( sideC->build_param_list( 1.0-fraction, distance, num_segs, 
                               through_vertex_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

RefVertex *
TDSplitSurface::start_vertex( CoEdge *co_edge_ptr )
{
  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  
  if ( co_edge_ptr->get_sense() == CUBIT_REVERSED )
    return ref_edge_ptr->end_vertex();
  else
    return ref_edge_ptr->start_vertex();
}

//================================================================================
// Description: A trivial class to hold two parameter values so that they 
//              can be stored in a DLIList.  They are the min and max parameter
//              space along the composite curve. (Split Surface Param)
// Author     : Steve Storm
// Date       : 2/3/2004
//================================================================================
SSParam::SSParam( double min, double max )
{
  uMin = min;
  uMax = max;
}

SSParam::~SSParam()
{
}

//================================================================================
// Description: This class (Split Surface Side) holds a chain of curves on one 
//              side of the surface.  It is needed to handle queries using a 
//              composite curve concept.
// Author     : Steve Storm
// Date       : 2/3/2004
//================================================================================
SSSide::SSSide( RefFace *ref_face_ptr, DLIList<CoEdge*> &co_edges,
                const CubitVector *collapsed_loc_ptr )
{
  refFacePtr = ref_face_ptr;
  isCollapsed = CUBIT_FALSE;

  coEdgeChain = co_edges;
  double param_low = 0.0;
  paramHigh = 0.0;
  int i;
  CoEdge *co_edge_ptr;
  for( i=coEdgeChain.size(); i--; )
  {
    co_edge_ptr = coEdgeChain.get_and_step();
    paramHigh += co_edge_ptr->get_ref_edge_ptr()->measure();
    coEdgeParamList.append( new SSParam( param_low, paramHigh) );
    param_low = paramHigh;
  }

  if( !coEdgeChain.size() )
  {
    // Need to add a single point to hold the location of this side
    // ie., the side is collapsed (point of a triangle)
    isCollapsed = CUBIT_TRUE;

    assert( collapsed_loc_ptr != NULL );

    paramHigh = 0.0;
    coordList.append( new CubitVector( *collapsed_loc_ptr ) );
    paramList.append( 0.0 );
  }
}

SSSide::~SSSide()
{
  // Free memory
  while( coEdgeParamList.size() ) 
    delete coEdgeParamList.pop();
  while( coordList.size() ) 
    delete coordList.pop();
}

CubitStatus
SSSide::position_from_u( double u_value, CubitVector &output_position)
{
  if( isCollapsed )
  {
    output_position = *coordList.get();
    return CUBIT_SUCCESS;
  }

  // Determine which CoEdge the given u is on
  int i;
  coEdgeChain.reset();
  coEdgeParamList.reset();
  CoEdge *co_edge_ptr = NULL;
  SSParam *param_ptr = NULL;
  double coedge_param_max;
  for( i=coEdgeChain.size(); i--; )
  {
    co_edge_ptr = coEdgeChain.get_and_step();
    param_ptr = coEdgeParamList.get_and_step();
    coedge_param_max = param_ptr->umax();
    if( u_value <= coedge_param_max )
      break;
  }

  // We have found the correct coedge.  Get it's RefEdge. 
  // Subtract it's start parameter.  
  // We now have the distance along the curve to traverse.
  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  double i_dist; // individual distance
  i_dist = u_value - param_ptr->umin();

  double ui_min, ui_max;
  ref_edge_ptr->get_param_range( ui_min, ui_max );

  if( co_edge_ptr->get_sense() == CUBIT_REVERSED ) 
    i_dist = ref_edge_ptr->measure() - i_dist;

  double ui = ref_edge_ptr->u_from_arc_length( ui_min, i_dist );

  return ref_edge_ptr->position_from_u( ui, output_position );
}

CubitStatus
SSSide::u_from_position( const CubitVector &input_position, double &u )
{
  if( isCollapsed )
  {
    u = 0.0;
    return CUBIT_SUCCESS;
  }

  // First check each curve to determine which one the input position is on
  int i;
  coEdgeChain.reset();
  coEdgeParamList.reset();
  CoEdge *co_edge_ptr = NULL;
  SSParam *param_ptr = NULL;
  int found = 0;
  CubitPointContainment pnt_containment;
  for( i=coEdgeChain.size(); i--; )
  {
    co_edge_ptr = coEdgeChain.get_and_step();
    param_ptr = coEdgeParamList.get_and_step();
    pnt_containment = co_edge_ptr->get_ref_edge_ptr()->point_containment( input_position );
    if( pnt_containment == CUBIT_PNT_ON )
    {
      found = 1;
      break;
    }
  }

  if( !found )
  {
    PRINT_ERROR( "Position %f, %f, %f not found on any curve.\n",
      input_position.x(), input_position.y(), input_position.z() );
    return CUBIT_FAILURE;
  }

  // Now we know which CoEdge it is on.  Get the curve. Get the parameter of 
  // the individual curve.
  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  double ui = ref_edge_ptr->u_from_position( input_position );
  
  // Get the parameter of the composite curve

  // Add the distance along the curve to its start param (umin)
  double sum = param_ptr->umin();

  double ui_min, ui_max;
  ref_edge_ptr->get_param_range( ui_min, ui_max );

  CubitSense sense = co_edge_ptr->get_sense();
  double root_param = (sense == CUBIT_FORWARD) ? ui_min : ui_max;
  double lfu = ( root_param < ui ) ?
    ref_edge_ptr->length_from_u( root_param, ui ) :
    ref_edge_ptr->length_from_u( ui, root_param );

  u = sum + fabs( lfu );

  return CUBIT_SUCCESS;
}

// Same as above function but curve is known
CubitStatus
SSSide::u_from_position( const CubitVector &input_position, 
                         CoEdge *co_edge_ptr, SSParam *param,
                         double &u )
{
  if( isCollapsed )
  {
    u = 0.0;
    return CUBIT_SUCCESS;
  }

  RefEdge *ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
  double ui = ref_edge_ptr->u_from_position( input_position );

  // Get the parameter of the composite curve
  double sum = param->umin();

  double ui_min, ui_max;
  ref_edge_ptr->get_param_range( ui_min, ui_max );

  CubitSense sense = co_edge_ptr->get_sense();
  double root_param = (sense == CUBIT_FORWARD) ? ui_min : ui_max;
  double lfu = ( root_param < ui ) ?
    ref_edge_ptr->length_from_u( root_param, ui ) :
    ref_edge_ptr->length_from_u( ui, root_param );

  u = sum + fabs( lfu );

  return CUBIT_SUCCESS;
}

CubitBoolean
SSSide::is_vertex_on( RefVertex *ref_vertex_ptr )
{
  if( isCollapsed )
  {
    // Compare coordinates
    CubitVector *side_coord = coordList.get();
    
    if( side_coord->about_equal( ref_vertex_ptr->coordinates() ) )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }

  // Check if it is on each coedge
  CubitVector ref_coords = ref_vertex_ptr->coordinates();
  CubitPointContainment pnt_containment;
  int i;
  CoEdge *co_edge_ptr;
  coEdgeChain.reset();
  for( i=coEdgeChain.size(); i--; )
  {
    co_edge_ptr = coEdgeChain.get_and_step();
    pnt_containment = co_edge_ptr->get_ref_edge_ptr()->
      point_containment( ref_coords );
    if( pnt_containment == CUBIT_PNT_ON )
      return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

CubitStatus
SSSide::build_param_list_from_facets( double tolerance )
{
  if( isCollapsed )
    return CUBIT_SUCCESS;

  int i, j, num_pnts;
  CoEdge *co_edge_ptr;
  SSParam *param_ptr;
  double param;
  CubitVector vec, *vec_ptr;

  coEdgeChain.reset();
  coEdgeParamList.reset();
  co_edge_ptr = coEdgeChain.get();
  param_ptr = coEdgeParamList.get();

  // Get the first coordinate
  if( co_edge_ptr->get_sense() == CUBIT_FORWARD )
    vec = co_edge_ptr->get_ref_edge_ptr()->start_coordinates();
  else
    vec = co_edge_ptr->get_ref_edge_ptr()->end_coordinates();

  // Add to list
  if( u_from_position( vec, co_edge_ptr, param_ptr, param ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  paramList.append( param );

  for( i=coEdgeChain.size(); i--; )
  {
    DLIList<CubitVector*> temp_vec_list;
    co_edge_ptr = coEdgeChain.get_and_step();
    param_ptr = coEdgeParamList.get_and_step();
    GMem *g_mem = new GMem;
    
    co_edge_ptr->get_ref_edge_ptr()->get_graphics( *g_mem, tolerance );
    num_pnts = g_mem->pointListCount;
    
    GPoint* point_list = g_mem->point_list();
    
    for( j=0; j<num_pnts; j++ )
    {
      vec_ptr = new CubitVector(
        point_list[j].x, point_list[j].y, point_list[j].z );
      temp_vec_list.append( vec_ptr );
    }

    delete g_mem;

    if( co_edge_ptr->get_sense() == CUBIT_REVERSED )
      temp_vec_list.reverse();

    // Calculate corresponding parameter values
    temp_vec_list.reset();
    // Skip the first point since already in the list
    vec_ptr = temp_vec_list.get_and_step();
    for( j=1; j<num_pnts; j++ )
    {
      vec_ptr = temp_vec_list.get_and_step();
      if( u_from_position( *vec_ptr, co_edge_ptr, param_ptr, param )
        == CUBIT_FAILURE )
        return CUBIT_FAILURE;
      paramList.append( param );
    }

    while( temp_vec_list.size() )
      delete temp_vec_list.pop();
  }

  // Extremes of paramList should be 0.0 and paramHigh - force them to
  // these values exactly.  This avoids some slight roundoff errors that
  // can cause an extra split point close to the start or end of a surface
  // that won't even be cleaned up in syncronize_lists.  This causes 
  // unexpected results at the start or end of the split.
  paramList.reset();
  double first_param = paramList.get_and_back();
  double last_param = paramList.get();
  paramList.reset();
  if( first_param > last_param )
  {
    paramList.change_to( paramHigh );
    paramList.back();
    paramList.change_to( 0.0 );
  }
  else
  {
    paramList.change_to( 0.0 );
    paramList.back();
    paramList.change_to( paramHigh );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SSSide::build_param_list( double fraction, double distance, int num_segs,
                          DLIList<RefVertex*> &through_vertex_list )
{
  int i;
  double frac, param;
  // Add paramList.size to account for potential collapsed edge that already
  // has an item in the list
  for( i=1+paramList.size(); i<num_segs; i++ )
  {
    // Example locations (we need to fill middle only): 0 1 2 3 4
    if( num_segs > 2 )
      frac = (double)i/num_segs;
    else if( through_vertex_list.size() )
    {
      int j;
      int found = 0;
      through_vertex_list.reset();
      for( j=through_vertex_list.size(); j--; )
      {
        RefVertex *ref_vertex_ptr = through_vertex_list.get_and_step();
        if( is_vertex_on( ref_vertex_ptr ) )
        {
          coordList.append( new CubitVector( ref_vertex_ptr->coordinates() ) );
          if( u_from_position( ref_vertex_ptr->coordinates(), param ) == CUBIT_FAILURE )
            return CUBIT_FAILURE;
          paramList.append( param );
          found = 1;

          // Keep track of vertices that were used - we can give a warning when
          // done if some weren't used.
          TDSplitSurface *tdss = (TDSplitSurface *)ref_vertex_ptr->
            get_TD(&TDSplitSurface::is_split_surface);
          if( !tdss )
            ref_vertex_ptr->add_TD( new TDSplitSurface( 0 ) );
          break;
        }
      }
      if( found )
        continue;
      else
        frac = fraction;
    }
    else
      frac = fraction;

    if( distance != -1.0 )
    {
      if( distance > paramHigh )
      {
        PRINT_ERROR( "Surface %d is not wide enough to support split distance of %f\n",
          refFacePtr->id(), distance );
        return CUBIT_FAILURE;
      }
      param = distance;
    }
    else
      param = frac*paramHigh;

    paramList.append( param );

    CubitVector vec;
    if( position_from_u( param, vec ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
    
    coordList.append( new CubitVector( vec ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
SSSide::syncronize_lists( SSSide *other_side, double param_tol )
{
  // This will syncronize paramList in "this" and other_side, also
  // update coordList correspondingly in both.

  // Note: sideB and sideD (which we are operating on here) will contain the
  // corner coords (sideA and sideC just contain the interior coords)

  // This should never happen but check anyway
  if( paramHigh == 0.0 && other_side->paramHigh == 0.0 )
    return CUBIT_FAILURE;

  paramList.reset();
  other_side->paramList.reset();

  DLIList<double> other_param_list = other_side->paramList;
  int i;
  double param;

  // Reverse other_param_list (which is on side D) so that it is going in the
  // same direction as this paramList (side B)
  other_param_list.reset();
  for( i=other_param_list.size(); i--; )
  {
    param = other_param_list.get();
    other_param_list.change_to( other_side->paramHigh - param );
    other_param_list.step();
  }
  
  other_param_list.reverse(); // May make sort faster

  // Parameter lists must be "scaled" to match
  double this_factor = 1.0;
  double other_factor = 1.0;
  double scale_factor = 1.0;
  
  // Only scale if neither side is a zero length side and sides are not equal
  // in length
  if( !isCollapsed && !other_side->is_collapsed() &&
      other_side->paramHigh != paramHigh )
  {
    // Use longer side as "baseline" - it may have more "features"
    // so we will be sure to capture those
    if( paramHigh > other_side->paramHigh )
    {
      other_param_list.reset();

      other_factor = paramHigh/other_side->paramHigh;
      scale_factor = other_factor;
      for( i=other_param_list.size(); i--; )
      {
        param = other_param_list.get();
        other_param_list.change_to( param*other_factor );
        other_param_list.step();
      }
    }
    else
    {
      paramList.reset();

      this_factor = other_side->paramHigh/paramHigh;
      scale_factor = this_factor;
      for( i=paramList.size(); i--; )
      {
        param = paramList.get();
        paramList.change_to( param*this_factor );
        paramList.step();
      }
    }
  }

  DLIList<double> combined_param_list = paramList;
  combined_param_list.merge_unique( other_param_list );

  // This should never happen, but check anyway
  if( combined_param_list.size() < 2 )
    return CUBIT_FAILURE;

  // Sort the list from low to high
  combined_param_list.sort();

  // Remove near values
  param_tol = param_tol*scale_factor; // Scale tolerance to length
  double prev_param;
  combined_param_list.reset();
  prev_param = combined_param_list.get_and_step();
  for( i=combined_param_list.size()-2; i--; )
  {
    param = combined_param_list.get();
    if( fabs(param - prev_param) < param_tol )
    {
      combined_param_list.change_to( -1.0 );
    }
    else
    {
      prev_param = param;
    }
    combined_param_list.step();
  }
  combined_param_list.remove_all_with_value(-1.0);

  // We need to check if there are near coincident points at the end of the list
  combined_param_list.last();
  prev_param = combined_param_list.get_and_back();
  param = combined_param_list.get();
  // Use 10.0 * GEOMETRY_RESABS because Granite cannot handle points
  // any closer together when making a spline
  if( fabs( prev_param-param ) < 10.0*GEOMETRY_RESABS*scale_factor )
  {
    // Remove the second to last value
    combined_param_list.last();
    combined_param_list.back();
    combined_param_list.remove();

    if( combined_param_list.size() < 2 )
      return CUBIT_FAILURE;
  }

  // Also check if there are near coincident points anywhere along the 
  // other (shorter) side (very doubtful, but check anyway)
  if( scale_factor > 1.0 )
  {
    param_tol = 10.0*GEOMETRY_RESABS*scale_factor;
    double prev_param;
    combined_param_list.reset();
    prev_param = combined_param_list.get_and_step();
    for( i=combined_param_list.size()-2; i--; )
    {
      param = combined_param_list.get();
      if( fabs(param - prev_param) < param_tol )
      {
        combined_param_list.change_to( -1.0 );
      }
      else
      {
        prev_param = param;
      }
      combined_param_list.step();
    }
    combined_param_list.remove_all_with_value(-1.0);
  }
  
  // Update this paramList
  if( isCollapsed )
  {
    // This is a collapsed triangle side, so just stuff it with the proper
    // number of params and coordinates
    paramList.clean_out();
    for( i=combined_param_list.size(); i--; )
      paramList.append( 0.0 );

    // It will already have one coordinate - copy it to the others
    CubitVector *vec;
    vec = coordList.get();
    for( i=combined_param_list.size()-1; i--; )
      coordList.append( new CubitVector( *vec ) );
  }
  else
  {
    paramList.clean_out();
    paramList = combined_param_list;

    // Generate this coordList.
    CubitVector vec;
    paramList.reset();
    for( i=paramList.size(); i--; )
    {
      param = paramList.get();
      param = param/this_factor; // Normalize it
      paramList.change_to( param ); // Set proper value in paramList
      position_from_u( param, vec );
      coordList.append( new CubitVector( vec ) );
      paramList.step();
    }
  }

  // Do the same for the other side
  if( other_side->is_collapsed() )
  {
    // This is a collapsed triangle side, so just stuff it with the proper
    // number of params and coordinates
    other_side->paramList.clean_out();
    for( i=combined_param_list.size(); i--; )
      other_side->paramList.append( 0.0 );

    // It will already have one coordinate - copy it to the others
    CubitVector *vec;
    vec = other_side->coordList.get();
    for( i=combined_param_list.size()-1; i--; )
      other_side->coordList.append( new CubitVector( *vec ) );
  }
  else
  {
    combined_param_list.reset();
    for( i=combined_param_list.size(); i--; )
    {
      param = combined_param_list.get();
      combined_param_list.change_to( other_side->paramHigh-param/other_factor );
      combined_param_list.step();
    }
    combined_param_list.reverse();

    CubitVector vec;
    other_side->paramList.clean_out();
    other_side->paramList = combined_param_list;
    combined_param_list.reset();
    for( i=combined_param_list.size(); i--; )
    {
      param = combined_param_list.get_and_step();
      other_side->position_from_u( param, vec );
      other_side->coordList.append( new CubitVector( vec ) );
    }
  }

  return CUBIT_SUCCESS;
}
