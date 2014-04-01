//-------------------------------------------------------------------------
// Filename      : GeometryUtil.cc
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/03/98
//-------------------------------------------------------------------------
#define COMPSURF_KEEP_CPT_STATS

#include "GeometryUtil.hpp"
#include "GeometryQueryEngine.hpp"

#include "CubitUtil.hpp"
#include "CubitString.hpp"
#include "CubitVector.hpp"

#include "BodySM.hpp"
#include "Lump.hpp"
#include "ShellSM.hpp"
#include "Surface.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"
#include "Curve.hpp"
#include "Point.hpp"

#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "Body.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "CoFace.hpp"
#include "CoEdge.hpp"
#include "Chain.hpp"
#include "CoVertex.hpp"

#include "DLIList.hpp"

#include "CpuTimer.hpp"
#include "GfxDebug.hpp"

GeometryUtil* GeometryUtil::instance_ = NULL;



//-------------------------------------------------------------------------
// Purpose       : Test if a point lies within a loop.
//
// Special Notes : Uses loop angle metric.
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 08/04/98
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::is_position_within_loop( 
                                              const CubitVector& position,
                                              Loop* loop_ptr,
                                              CubitVector* closest_on_loop,
                                              CubitVector* passed_normal )
{
  CubitVector surf_point, normal, pt_on_curve;
  if( passed_normal )
  {
    surf_point = position;
    normal = *passed_normal;
  }
  else
  {
    loop_ptr->get_ref_face_ptr()->get_surface_ptr()->
      closest_point( position, &surf_point, &normal );
  }
  
  CoEdge* other_coedge = 0;
  CoEdge* closest_coedge = 
    closest_loop_coedge( loop_ptr, position, other_coedge, &pt_on_curve );
  //If we got back NULL, either the loop has no CoEdges, or all Curves
  //occur multiple times.
  if( !closest_coedge ) return CUBIT_TRUE;
  
  if( closest_on_loop ) *closest_on_loop = pt_on_curve;
  
  CubitVector coe_normal, coe_cross, tangent1, tangent2, junk;

  if( ! other_coedge )
  {
    Curve* curve_ptr = closest_coedge->get_ref_edge_ptr()->get_curve_ptr();
    double u = curve_ptr->u_from_position( pt_on_curve );
    
    /**** Special case: closest point on loop at G1 discontinuity ****/
    if( curve_ptr->G1_discontinuous( u, &tangent1, &tangent2 ) )
    {
      if( closest_coedge->get_sense() == CUBIT_REVERSED )
      {
        tangent1 *= -1.0;
        tangent2 *= -1.0;
      }
      loop_ptr->get_ref_face_ptr()->get_surface_ptr()->
        closest_point( pt_on_curve, NULL, &coe_normal );
      coe_cross = tangent1 * tangent2;
      double sum  = (coe_cross + coe_normal).length_squared();
      double diff = (coe_cross - coe_normal).length_squared();
    
      CubitBoolean inside1 = 
        inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
      CubitBoolean inside2 = 
        inside_of_curve( tangent2, pt_on_curve, surf_point, normal );
        
      if( (sum > diff) || (diff < CUBIT_DBL_MIN) ) 
        //discontinuity is at a convexity
      {
        //the point must be inside of both sub-edges
        return (inside1 && inside2) ? CUBIT_TRUE : CUBIT_FALSE;
      }
      else //discontinuity is at a concavity
      {
        //the point must be inside of one of the sub-edges
        return (inside1 || inside2) ? CUBIT_TRUE : CUBIT_FALSE;
      }
    }

    else /**** This is the normal, non-special case part ****/
    { 
      curve_ptr->closest_point( pt_on_curve, junk, &tangent1 );
      if( closest_coedge->get_sense() == CUBIT_REVERSED )
        tangent1 *= -1.0;
        
      return inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
    }
  }
    /**** Special case: closest point on loop at vertex ****/
  else
  {
    Curve* curve1_ptr = closest_coedge->get_ref_edge_ptr()->get_curve_ptr();
    Curve* curve2_ptr =   other_coedge->get_ref_edge_ptr()->get_curve_ptr();
    curve1_ptr->closest_point( pt_on_curve, junk, &tangent1 );
    curve2_ptr->closest_point( pt_on_curve, junk, &tangent2 );

    if( closest_coedge->get_sense() == CUBIT_REVERSED )
      tangent1 *= -1.0;
    if(   other_coedge->get_sense() == CUBIT_REVERSED )
      tangent2 *= -1.0;

    loop_ptr->get_ref_face_ptr()->get_surface_ptr()->
      closest_point( pt_on_curve, NULL, &coe_normal );
    
    coe_cross = tangent1 * tangent2;
    double sum  = (coe_cross + coe_normal).length_squared();
    double diff = (coe_cross - coe_normal).length_squared();
    
    CubitBoolean inside1 = 
      inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
    CubitBoolean inside2 = 
      inside_of_curve( tangent2, pt_on_curve, surf_point, normal );
        
    if( (sum > diff) || (diff < CUBIT_DBL_MIN) ) 
      //the common vertex is at a convexity
    {
      //the point must be inside of both coedges
      return (inside1 && inside2) ? CUBIT_TRUE : CUBIT_FALSE;
    }
    else //the common vertex is at a concavity
    {
      //the point must be inside of one of the coedges
      return (inside1 || inside2) ? CUBIT_TRUE : CUBIT_FALSE;
    }
  }
}



                            
//-------------------------------------------------------------------------
// Purpose       : Attempt to calculate the area of a loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/19/99
//-------------------------------------------------------------------------
/*
double GeometryUtil::loop_area( Loop* loop_ptr )
{
  int i;
  
  DLIList<CoEdge*> coedges;
  loop_ptr->co_edges( coedges );
  
  RefFace* ref_face_ptr = loop_ptr->get_ref_face_ptr();
  
  
    // Get list of segment points around loop
  DLIList<CubitVector*> point_list, interior_edge_points;
  coedges.reset();
  for( i = coedges.size(); i > 0; i-- )
  {
    CoEdge* coedge_ptr = coedges.get_and_step();
    RefEdge* edge_ptr = coedge_ptr->get_ref_edge_ptr();
    CubitSense list_sense;
    interior_edge_points.clean_out();
    edge_ptr->get_interior_extrema( interior_edge_points, list_sense );
    if( coedge_ptr->get_sense() == CUBIT_REVERSED )
    {
      list_sense = CubitUtil::opposite_sense( list_sense );
      point_list.append( new CubitVector( 
        edge_ptr->end_vertex()->coordinates() ) );
    }
    else point_list.append( new CubitVector( 
      edge_ptr->start_vertex()->coordinates() ) );
      
    if( list_sense == CUBIT_REVERSED )
      interior_edge_points.reverse();
    
    point_list += interior_edge_points;
  }
  
    // Pick any point for which the refface has a valid normal
  const double STEP = 1e2 * GEOMETRY_RESABS;
  const double SMALL = CUBIT_RESABS * CUBIT_RESABS;
  CubitVector normal, *point, *prev, *current;
  for( i = point_list.size(); i--; )
  {
    point = point_list.get_and_step();
    normal = ref_face_ptr->normal_at( *point );
    if( normal.length_squared() > SMALL )
      break;
  }
  
    // If we did not find a valid normal at any point (this should
    // almost never happen), try stepping a small distance along 
    // each segment
  for( i = point_list.size(); (normal.length_squared() <= SMALL) && (i > 0); i-- )
  {
    prev = point;
    point = point_list.get_and_step();
    CubitVector direction = *prev - *point;
    double length = direction.length();
    if( length < CUBIT_RESABS )
      continue;
    double step_fraction = STEP / length;
    if( step_fraction >= 1.0 )
      continue;
    
    CubitVector position = *point + step_fraction * direction;
    normal = ref_face_ptr->normal_at( position );
    if( normal.length_squared() > SMALL )
    {
      point = new CubitVector( position );
      point_list.back();
      point_list.insert( point );
      break;
    }
  }
  
    // RefFace does not have a valid normal anywhere??
  if( normal.length_squared() <= SMALL )
    return 0.0;
    
    // Do simple 2-D area calculation and hope it isn't too
    // inaccurate for the 3-D refface.
  double double_area = 0.0;
  assert( point_list.prev() == point );
  prev = point_list.get_and_step();
  current = point_list.get_and_step();
    
  while( current != point )
  {
    CubitVector cross = (*prev - *point) * (*current  - *point);
    double dot_product = cross % normal;
    double_area += dot_product >= 0.0 ? cross.length() : -(cross.length());
    prev = current;
    current = point_list.get_and_step();
  }

  while( point_list.size() > 0 ) delete point_list.pop();
  return double_area / 2.0;
}
*/

struct CoEdgeDataList {
  CoEdge* coe_ptr;
  Curve* curve_ptr;
  double box_dist_sqr;
  int mark;
};


//-------------------------------------------------------------------------
// Purpose       : Find the closest point on a loop.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 08/04/98
//-------------------------------------------------------------------------
CoEdge* GeometryUtil::closest_loop_coedge( Loop* loop_ptr,
                                           const CubitVector& from_pt,
                                           CoEdge*& other_coedge,
                                           CubitVector* closest )
{
  static CoEdgeDataList* coe_array = NULL;
  static int coe_array_size = 0;
  const double ressqr = CUBIT_RESABS * CUBIT_RESABS;
  
  DLIList<CoEdge*> coedge_list;
  coedge_list.clean_out();
    
  int i;
  CubitVector closest_pt, current_pt;
  double current_dst, closest_dst;
  CoEdge* closest_coedge = 0;

  //We need the ordered coedges so that we can make sure that, when
  //we pass back other_coedge, it is the coedge following the returned
  //coedge in the loop, and not the other way around.
  loop_ptr->ordered_co_edges( coedge_list );
  
  //Make sure we have enough space in the array
  int count = coedge_list.size();
  if( coe_array_size < count )
  {
    if( coe_array ) delete [] coe_array;
    coe_array = new CoEdgeDataList[count];
    coe_array_size = count;
  }
  
  //remove from consideration any coedges for which the RefEdge occurs
  //more than once on the face.
  for( i = coedge_list.size(); i > 0; i-- )
    coedge_list.get_and_step()->get_ref_edge_ptr()->marked( 0 );
  for( i = coedge_list.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = coedge_list.get_and_step()->get_ref_edge_ptr();
    edge_ptr->marked( edge_ptr->marked() + 1 );
  }

  //Start at one of the coedges we want to remove.
  coedge_list.reset();
  for( i = coedge_list.size(); i > 0; i-- )
    if( coedge_list.get()->get_ref_edge_ptr()->marked() > 1 )
      break;
    else 
      coedge_list.step();
  
  //Populate the array with the coedges we want to keep.
  count = 0;
  int last_mark = 0;
  for( i = coedge_list.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = coedge_list.get()->get_ref_edge_ptr();
    if( edge_ptr->marked() == 1 )
    {
      coe_array[count].coe_ptr = coedge_list.get();
      coe_array[count].curve_ptr = edge_ptr->get_curve_ptr();
      coe_array[count].box_dist_sqr = 
        coe_array[count].curve_ptr->bounding_box().distance_squared(from_pt);
      coe_array[count].mark = 0;
      count++;
    }
    else if( count == 0 )
      last_mark = 1;
    else 
      coe_array[count-1].mark = 1;
    
    edge_ptr->marked(0);
    coedge_list.step();
  }
  if( last_mark == 0 )
    coe_array[count-1].mark = 1;
  
  //Do we have any coedges left?
  if( count == 0 )
    return 0;
  
  //Find the first bounding box that the point is within
  int current_coe = -1;
  for( i = 0; (i < count) && (current_coe < 0); i++ )
    if( coe_array[i].box_dist_sqr < CUBIT_RESABS )
      current_coe = i;
  
  //If the point wasn't in any bounding box, choose the closest one
  closest_dst = CUBIT_DBL_MAX;
  if( current_coe < 0 )
    for( i = 0; i < count; i++ )
      if( coe_array[i].box_dist_sqr < closest_dst ) 
      {
        closest_dst = coe_array[i].box_dist_sqr;
        current_coe = i;
      }
  
  //Start by assuming that the curve we found above is the closest
  closest_coedge = coe_array[current_coe].coe_ptr;
  coe_array[current_coe].curve_ptr->closest_point_trimmed( from_pt, closest_pt );
  closest_dst = ( from_pt - closest_pt ).length_squared();
  other_coedge = 0;
  
  int start = current_coe;
  int closest_index = start;
  int other_index = -1;
  current_coe = ( current_coe + 1 ) % count;
  
  //Now look for a closer curve
  while( start != current_coe )
  {
    if( coe_array[current_coe].box_dist_sqr <= (closest_dst+GEOMETRY_RESABS) )
    {
      coe_array[current_coe].curve_ptr
        ->closest_point_trimmed( from_pt, current_pt );
      current_dst = ( from_pt - current_pt ).length_squared();
    
      if( (closest_pt - current_pt).length_squared() < ressqr )
      {
        if( current_dst < closest_dst ) 
        {
          closest_dst = current_dst;
          closest_pt = current_pt;
          other_coedge = coe_array[current_coe].coe_ptr;
        }
        else
        {
          other_coedge = coe_array[current_coe].coe_ptr;
        }
        other_index = current_coe;
      } 
      else if( current_dst < closest_dst )
      {
        closest_coedge = coe_array[current_coe].coe_ptr;
        closest_dst = current_dst;
        closest_pt = current_pt;
        other_coedge = 0;
        closest_index = current_coe;
        other_index = -1;
      }
    }
    current_coe = ( current_coe + 1 ) % count;
  }
  
  //make sure we have things in the correct order
  if( other_coedge )
  {
    RefVertex* common = closest_coedge->get_ref_edge_ptr()->
      common_ref_vertex( other_coedge->get_ref_edge_ptr() );
    if( ! common ) 
      other_coedge = 0;
  }
  
  CubitBoolean swap = CUBIT_FALSE;
  if( other_coedge )
  {
    if( ((closest_index + 1)%count) == other_index )
      swap = CUBIT_FALSE;
    else if( ((other_index+1)%count) == closest_index )
      swap = CUBIT_TRUE;
    else
    {
      //locate the start and end of the subsection of the loop
      int begin, end;
      if( closest_index < other_index )
      { 
        begin = closest_index;
        end = other_index;
      }
      else
      {
        begin = other_index;
        end = closest_index;
      }
      
      while( (begin > 0) && (coe_array[begin-1].mark != 1) )
        begin--;
      
      while( (end < count) && (coe_array[end].mark != 1) )
        end++;
        
      if( (closest_index == begin) && (other_index == end) )
        swap = CUBIT_TRUE;
      else
        swap = CUBIT_FALSE;
    }
  }
   
  if( swap )
  {
    CoEdge* tmp = other_coedge;
    other_coedge = closest_coedge;
    closest_coedge = tmp;
  }

  if( closest ) *closest = closest_pt;
  return closest_coedge;
}
CoEdge* GeometryUtil::closest_face_coedge( RefFace* face_ptr,
                                           const CubitVector& from_pt,
                                           CoEdge*& other_coedge,
                                           CubitVector* closest )
{
  static CoEdgeDataList* coe_array = NULL;
  static int coe_array_size = 0;
  const double ressqr = GEOMETRY_RESABS * GEOMETRY_RESABS;
  
  DLIList<CoEdge*> coedge_list;
  coedge_list.clean_out();
    
  int i;
  CubitVector closest_pt, current_pt;
  double current_dst, closest_dst;
  CoEdge* closest_coedge = 0;

  face_ptr->co_edges( coedge_list );
  
  //Make sure we have enough space in the array
  int count = coedge_list.size();
  if( coe_array_size < count )
  {
    if( coe_array ) delete [] coe_array;
    coe_array = new CoEdgeDataList[count];
    coe_array_size = count;
  }
  
  //remove from consideration any coedges for which the RefEdge occurs
  //more than once on the face.
  for( i = coedge_list.size(); i > 0; i-- )
    coedge_list.get_and_step()->get_ref_edge_ptr()->marked( 0 );
  for( i = coedge_list.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = coedge_list.get_and_step()->get_ref_edge_ptr();
    edge_ptr->marked( edge_ptr->marked() + 1 );
  }

  //Start at one of the coedges we want to remove.
  coedge_list.reset();
  for( i = coedge_list.size(); i > 0; i-- )
    if( coedge_list.get()->get_ref_edge_ptr()->marked() > 1 )
      break;
    else 
      coedge_list.step();
  
  //Populate the array with the coedges we want to keep.
  count = 0;
  int last_mark = 0;
  for( i = coedge_list.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = coedge_list.get()->get_ref_edge_ptr();
    if( edge_ptr->marked() == 1 )
    {
      coe_array[count].coe_ptr = coedge_list.get();
      coe_array[count].curve_ptr = edge_ptr->get_curve_ptr();
      coe_array[count].box_dist_sqr = 
        coe_array[count].curve_ptr->bounding_box().distance_squared(from_pt);
      coe_array[count].mark = 0;
      count++;
    }
    else if( count == 0 )
      last_mark = 1;
    else 
      coe_array[count-1].mark = 1;
    
    edge_ptr->marked(0);
    coedge_list.step();
  }
  if( last_mark == 0 )
    coe_array[count-1].mark = 1;
  
  
  //Find the first bounding box that the point is within
  int current_coe = -1;
  for( i = 0; (i < count) && (current_coe < 0); i++ )
    if( coe_array[i].box_dist_sqr < ressqr )
      current_coe = i;
  
  //If the point wasn't in any bounding box, choose the closest one
  closest_dst = CUBIT_DBL_MAX;
  if( current_coe < 0 )
    for( i = 0; i < count; i++ )
      if( coe_array[i].box_dist_sqr < closest_dst ) 
      {
        closest_dst = coe_array[i].box_dist_sqr;
        current_coe = i;
      }
  
  //Start by assuming that the curve we found above is the closest
  closest_coedge = coe_array[current_coe].coe_ptr;
  coe_array[current_coe].curve_ptr->closest_point_trimmed( from_pt, closest_pt );
  closest_dst = ( from_pt - closest_pt ).length_squared();
  other_coedge = 0;
  
  int start = current_coe;
  current_coe = ( current_coe + 1 ) % count;
  
  //Now look for a closer curve
  while( start != current_coe )
  {
    if( coe_array[current_coe].box_dist_sqr <= (closest_dst+GEOMETRY_RESABS) )
    {
      coe_array[current_coe].curve_ptr
        ->closest_point_trimmed( from_pt, current_pt );
      current_dst = ( from_pt - current_pt ).length_squared();
    
      if( (closest_pt - current_pt).length_squared() < ressqr )
      {
        if( current_dst < closest_dst ) 
        {
          closest_dst = current_dst;
          closest_pt = current_pt;
          other_coedge = coe_array[current_coe].coe_ptr;
        }
        else
        {
          other_coedge = coe_array[current_coe].coe_ptr;
        }
//        other_index = current_coe;
      } 
      else if( current_dst < closest_dst )
      {
        closest_coedge = coe_array[current_coe].coe_ptr;
        closest_dst = current_dst;
        closest_pt = current_pt;
        other_coedge = 0;
//        closest_index = current_coe;
//        other_index = -1;
      }
    }
    current_coe = ( current_coe + 1 ) % count;
  }
  
  //make sure we have things in the correct order
  if( other_coedge )
  {
    RefEdge* closest_edge = closest_coedge->get_ref_edge_ptr();
    RefEdge* other_edge = other_coedge->get_ref_edge_ptr();
    RefVertex* common = closest_edge->common_ref_vertex( other_edge );
    if( ! common ) 
      other_coedge = 0;
    else if( closest_coedge->get_sense() == CUBIT_FORWARD )
    {
      if( closest_edge->start_vertex() == common )
      {
        CoEdge* tmp = closest_coedge;
        closest_coedge = other_coedge;
        other_coedge = tmp;
      }
    }
    else if( closest_coedge->get_sense() == CUBIT_REVERSED )
    {
      if( closest_edge->end_vertex() == common )
      {
        CoEdge* tmp = closest_coedge;
        closest_coedge = other_coedge;
        other_coedge = tmp;
      }
    }
  }

  if( closest ) *closest = closest_pt;
  return closest_coedge;
}

//-------------------------------------------------------------------------
// Purpose       : recursive method to make a segmented curve on a surface.
//
// Special Notes : This method is was written for use by VirtualCurve,
//                 which is why virtual curve debug flags (86) appear
//                 here.
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 08/07/98
//-------------------------------------------------------------------------
CubitStatus GeometryUtil::recursive_make_curve( Surface* surface_ptr,
                                         const CubitVector& start_pt,
                                         const CubitVector& end_pt,
                                         DLIList<CubitVector*>& segment_points,
                                         double arc_angle_tol,
                                         double midpoint_dist_tol )
{
  if( midpoint_dist_tol < GEOMETRY_RESABS )
    midpoint_dist_tol = GEOMETRY_RESABS;
  
  CubitVector seg_mid_point, mid_point, mid_pt_norm;
  CubitVector start_pt_norm, end_pt_norm, tangent;
  double angle = arc_angle_tol;
  double midpt_dist = midpoint_dist_tol;
  
  seg_mid_point = (start_pt + end_pt) / 2;
  surface_ptr->closest_point_trimmed( seg_mid_point, mid_point );
  surface_ptr->closest_point( mid_point, NULL, &mid_pt_norm   );
  surface_ptr->closest_point( start_pt,  NULL, &start_pt_norm );
  surface_ptr->closest_point( end_pt,    NULL, &end_pt_norm   );
  tangent = start_pt - end_pt;
  
  if( tangent.length() < GEOMETRY_RESABS )
  {
    PRINT_DEBUG_86( "Solution for polyline on surface not converging.\n"
                    "GeometryUtil::recursive_make_curve(..) returning "
                    "FAILURE.\n");
    return CUBIT_FAILURE;
  }
  
  CubitBoolean split = CUBIT_FALSE;
  
  if( !split && ( midpoint_dist_tol < CUBIT_DBL_MAX ) )
  {
    midpt_dist = (seg_mid_point - mid_point).length();
    if( midpt_dist > midpoint_dist_tol )
    {
      split = CUBIT_TRUE;
    }
  }

  if( !split && (arc_angle_tol < (CUBIT_PI / 2)) )
  {
    double seg_angle, half1_angle, half2_angle, sum;
    seg_angle = fabs( tangent.vector_angle( start_pt_norm, end_pt_norm ) );
    half1_angle = fabs((start_pt-mid_point).
                        vector_angle( start_pt_norm, mid_pt_norm ));
    half2_angle = fabs((end_pt-mid_point).
                        vector_angle( end_pt_norm, mid_pt_norm ));
    sum = half1_angle + half2_angle;
    angle = seg_angle > sum ? seg_angle : sum;
    
    if( angle > arc_angle_tol ) 
    {
      split = CUBIT_TRUE;
    }
  }
  
  if( split )
  {
    if( DEBUG_FLAG(86) )
    {
      GfxDebug::draw_label( ++linearized_curve_debug_count_,
        float(mid_point.x()), float(mid_point.y()), float(mid_point.z()),
        CUBIT_WHITE );
      GfxDebug::flush();
    }

    double start_len = (start_pt - mid_point).length();
    double end_len = (end_pt - mid_point).length();
    double half_tol = 0.5 * midpoint_dist_tol;
    if( ((start_len < half_tol) && (end_len < half_tol)) ||
        (start_len < GEOMETRY_RESABS) || (end_len < GEOMETRY_RESABS) )
        return CUBIT_SUCCESS;
    
    CubitStatus status;
    status = recursive_make_curve( surface_ptr, start_pt, mid_point, 
             segment_points, arc_angle_tol, midpoint_dist_tol );
    segment_points.append( new CubitVector( mid_point ) );
    
    if( status != CUBIT_SUCCESS ) return status;
    
    status = recursive_make_curve( surface_ptr, mid_point, end_pt,
             segment_points, arc_angle_tol, midpoint_dist_tol );

    return status;
  }

  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 08/03/98
//-------------------------------------------------------------------------
GeometryUtil::GeometryUtil()
{
  default_angle_tol = CUBIT_PI / 20;
  default_midpoint_tol = CUBIT_DBL_MAX;
  surface_cpu_time = 0.0;
  other_cpu_time = 0.0;
}



//-------------------------------------------------------------------------
// Purpose       : Find the closest coface in a Shell to the passed point
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/24/99
//-------------------------------------------------------------------------
CubitStatus GeometryUtil::closest_shell_coface( Shell* shell_ptr,
                                            const CubitVector& from_pt,
                                            DLIList<CoFace*>& result_set,
                                            CubitVector* closest )
{
  static DLIList<CoFace*> coface_list, shell_cofaces;
  coface_list.clean_out();
  shell_cofaces.clean_out();
  result_set.clean_out();
  
  int i;
  CubitVector closest_pt, current_pt;
  CoFace *current_cof;
  double current_dst, closest_dst;

  shell_ptr->co_faces( shell_cofaces );
  
  //remove from consideration any cofaces for which the RefFace occurs
  //more than once on the Shell.
  for( i = shell_cofaces.size(); i > 0; i-- )
    shell_cofaces.get_and_step()->get_ref_face_ptr()->marked( 0 );
  for( i = shell_cofaces.size(); i > 0; i-- )
  {
    RefFace* face_ptr = shell_cofaces.get_and_step()->get_ref_face_ptr();
    face_ptr->marked( face_ptr->marked() + 1 );
  }
  for( i = shell_cofaces.size(); i > 0; i-- )
  {
    CoFace* cof_ptr = shell_cofaces.get_and_step();
    if( cof_ptr->get_ref_face_ptr()->marked() == 1 )
      coface_list.append( cof_ptr );
  }
  for( i = shell_cofaces.size(); i > 0; i-- )
    shell_cofaces.get_and_step()->get_ref_face_ptr()->marked( 0 );
  
  if( !coface_list.size() )
    return CUBIT_FAILURE;
  
  current_cof = coface_list.get_and_step();
  current_cof->get_ref_face_ptr()->
    find_closest_point_trimmed( from_pt, closest_pt );
  closest_dst = ( from_pt - closest_pt ).length_squared();
  result_set.append( current_cof );

  for( i = coface_list.size(); i > 1; i-- )
  {
    current_cof = coface_list.get_and_step();
    current_cof->get_ref_face_ptr()->
      find_closest_point_trimmed( from_pt, current_pt );
    current_dst = ( from_pt - current_pt ).length_squared();
    
    if( fabs( current_dst - closest_dst ) < CUBIT_RESABS )
    {
      result_set.append( current_cof );
    } 
    else if( current_dst < closest_dst )
    {
      result_set.clean_out();
      result_set.append( current_cof );
      closest_dst = current_dst;
      closest_pt = current_pt;
    }
  }
  if( closest != NULL ) *closest = closest_pt;
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Test if a position is within a Shell
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/24/99
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::is_position_within_shell( 
                                            const CubitVector& position,
                                            Shell* shell_ptr,
                                            CubitVector* closest_on_shell )
{
  DLIList<CoFace*> result_cofaces;
  CubitVector closest_pt;
  closest_shell_coface( shell_ptr, position, result_cofaces, &closest_pt );
  if( closest_on_shell ) *closest_on_shell = closest_pt;
  
  //test each resulting coface
  CubitVector vector_to_surf = closest_pt - position;
  for( int i = result_cofaces.size(); i > 0; i-- )
  {
    CoFace* coface_ptr = result_cofaces.get_and_step();
    RefFace* ref_face_ptr = coface_ptr->get_ref_face_ptr();
    CubitVector normal = ref_face_ptr->normal_at( closest_pt );
    if( coface_ptr->get_sense() == CUBIT_REVERSED )
      normal *= -1.0;
    
    double sum  = (normal + vector_to_surf).length_squared();
    double diff = (normal - vector_to_surf).length_squared();
    if( diff > sum ) return CUBIT_FALSE;
  }
  return CUBIT_TRUE;
}


//-------------------------------------------------------------------------
// Purpose       : Test if a shell is a void
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/27/99
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::is_shell_a_void( Shell* shell_ptr )
{
  assert( shell_ptr != 0 );
  DLIList<CoFace*> shell_cofaces;
  
  CubitBox box = shell_ptr->bounding_box();
  double box_diag = box.diagonal().length();
  double epsilon = box_diag / 1000;
  
  CubitVector base_point( box.minimum() );
  CubitVector opposite_point;
  CubitVector point, normal;
  
  closest_shell_coface( shell_ptr, base_point, shell_cofaces, &point );
  base_point = point;
  
  if( !shell_cofaces.size() )
    return CUBIT_TRUE;
    
  CoFace* coface_ptr = shell_cofaces.get();
  normal = coface_ptr->get_ref_face_ptr()->normal_at( base_point );
  if( coface_ptr->get_sense() == CUBIT_REVERSED ) normal *= -1.0;
  normal.length( box_diag );
  
  opposite_point = base_point + normal;
  closest_shell_coface( shell_ptr, opposite_point, shell_cofaces, &point );
  
  return (base_point - point).length() < epsilon ? CUBIT_FALSE : CUBIT_TRUE;
}
  

//-------------------------------------------------------------------------
// Purpose       : Check if a point is on the inside or outside of a curve,
//                 given the surface normal, and the sense of the curve on
//                 the surface.
//
// Special Notes : The passed position is assumed to be on the surface.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/02/00
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::inside_of_curve( const CubitVector& tangent,
                                            const CubitVector& curve_pos,
                                            const CubitVector& surf_pos,
                                            const CubitVector& normal )
{
  CubitVector cross = tangent * ( surf_pos - curve_pos );
  double sum  = (cross + normal).length_squared();
  double diff = (cross - normal).length_squared();
  return (sum < diff) ? CUBIT_FALSE : CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Given a coordinate and a list of curves, determine if
//                 any of the curves has that coordinate as an end point
//                 and return that curve and whether it is start or end. 
//
// Special Notes : The passed position is assumed to be on the surface.
//
// Creator       : Steve Storm
//
// Creation Date : 05/30/00
//-------------------------------------------------------------------------
RefEdge* GeometryUtil::
find_connected_ref_edge_by_coord( CubitVector& coords, 
                                  DLIList<RefEdge*>& ref_edge_list,
                                  int& start_flg )
{
   RefEdge* ref_edge;
   DLIList<RefEdge*> tmp_list = ref_edge_list;

   ref_edge = ref_edge_list.get();
   GeometryQueryEngine* gqe_ptr = ref_edge->get_geometry_query_engine();
   double tol = gqe_ptr->get_sme_resabs_tolerance();

   tmp_list.reset();
   for( int i=0; i<tmp_list.size(); i++ )
   {
      ref_edge = tmp_list.get_and_step();

      if( ref_edge->start_coordinates().within_tolerance( coords, tol ) )
      {
         start_flg = 1;
         return ref_edge;
      }

      if( ref_edge->end_coordinates().within_tolerance( coords, tol ) )
      {
         start_flg = 0;
         return ref_edge;
      }
   }
   return NULL;
}

void 
GeometryUtil::find_connected_ref_edges_by_coord( CubitVector& coords, 
                                                 DLIList<RefEdge*>& ref_edge_list,
                                                 DLIList<RefEdge*>& connected_edges_list )
{
   RefEdge* ref_edge;
   DLIList<RefEdge*> tmp_list = ref_edge_list;

   ref_edge = ref_edge_list.get();
   GeometryQueryEngine* gqe_ptr = ref_edge->get_geometry_query_engine();
   double tol = gqe_ptr->get_sme_resabs_tolerance();

   tmp_list.reset();
   for( int i=0; i<tmp_list.size(); i++ )
   {
      ref_edge = tmp_list.get_and_step();

      if( ref_edge->start_coordinates().within_tolerance( coords, tol ) )
         connected_edges_list.append_unique( ref_edge );

      if( ref_edge->end_coordinates().within_tolerance( coords, tol ) )
         connected_edges_list.append_unique( ref_edge );
   }
}

CubitStatus 
GeometryUtil::form_ref_edge_loop_by_coord( DLIList<RefEdge*>& ref_edge_list )
{
   int i;

   if( ref_edge_list.size() == 0 )
      return CUBIT_FAILURE;

   if( ref_edge_list.size() == 1 )
   {
      double period;
      if( ref_edge_list.get()->is_periodic( period ) )
         return CUBIT_SUCCESS;
      else
      {
         // Could be a spline
         RefEdge *ref_edge_ptr = ref_edge_list.get();
         GeometryQueryEngine* gqe_ptr = ref_edge_ptr->get_geometry_query_engine();
         double tol = gqe_ptr->get_sme_resabs_tolerance();
         if( ref_edge_ptr->start_coordinates().within_tolerance( 
            ref_edge_ptr->end_coordinates(), tol ) )
            return CUBIT_SUCCESS;
      }
      return CUBIT_FAILURE;
   }

   // First form a chain (which can be a loop)
   if( form_ref_edge_chain_by_coord( ref_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   // Now check to make sure a valid loop was formed.  For a valid loop 
   // every vertex must have exactly 2 edges connected to it.  The only
   // exception is a periodic curve where there is only one vertex, which
   // was checked for before.
   DLIList<RefEdge*> connected_edges;
   RefEdge* ref_edge_ptr;
   CubitVector vector;
   for( i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
   
      connected_edges.clean_out();
      vector = ref_edge_ptr->start_coordinates();
      find_connected_ref_edges_by_coord( vector ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
         return CUBIT_FAILURE;

      connected_edges.clean_out();
      vector =  ref_edge_ptr->end_coordinates();
      find_connected_ref_edges_by_coord( vector ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
         return CUBIT_FAILURE;
   }
 
   return CUBIT_SUCCESS;
}

CubitStatus 
GeometryUtil::form_ref_edge_chain_by_coord( DLIList<RefEdge*>& ref_edge_list )
{
   if( ref_edge_list.size() == 1 )
      return CUBIT_SUCCESS;

   int i;
   int start_flg = 0;
   
   CubitVector vector_start, vector_end;
   RefEdge* curr_edge;
   RefEdge* next_edge;

   DLIList<RefEdge*> temp_list = ref_edge_list;
   DLIList<RefEdge*> ordered_list;
   
   // Start vertex must be attached to only one ref_edge in our list. 
   // Otherwise this is really a loop, which is still okay

   // Search to find the starting edge
   RefEdge* start_edge = NULL;
   for( i=0; i<temp_list.size(); i++ )
   {
      curr_edge = temp_list.get();

      // Since DLList doesn't allow a NULL item
      DLIList<RefEdge*> temp_list2;
      temp_list2 = temp_list;
      temp_list2.remove( curr_edge ); // Next item becomes current item
   
      // get the vectors to start and end
      vector_start = curr_edge->start_coordinates();
      vector_end   = curr_edge->end_coordinates();
      
      if( !find_connected_ref_edge_by_coord( vector_start ,
         temp_list2, start_flg ) ||
         !find_connected_ref_edge_by_coord( vector_end,
         temp_list2, start_flg ) )

      {
         start_edge = curr_edge;
         //temp_list.change_to( curr_edge );
         break;
      }
   }

   temp_list.reset();

   if( start_edge == NULL )
      start_edge = ref_edge_list.get();
   
   curr_edge = start_edge;
   ordered_list.append( curr_edge );
   temp_list.move_to( start_edge );
   temp_list.remove();

   next_edge = temp_list.get();

   // Attempt to use second edge in list as second edge in chain to 
   // preserve original chain direction.
   RefEdge* next_edge_tmp = NULL;
   RefEdge* next_edge_tmp2 = NULL;
   // Check start first - if this fails check end next so chain order
   // will be established properly.
   vector_start =  curr_edge->start_coordinates();
   next_edge_tmp = find_connected_ref_edge_by_coord( vector_start,
                                                     temp_list, start_flg );
   if( next_edge_tmp != next_edge )
   {
      vector_end    = curr_edge->end_coordinates();
      next_edge_tmp2 = find_connected_ref_edge_by_coord( vector_end,
                                                         temp_list, start_flg );
   }

   if( next_edge_tmp2 )
     next_edge_tmp = next_edge_tmp2;

   if( next_edge_tmp  )
   {
      curr_edge = next_edge_tmp;
      ordered_list.append( next_edge_tmp );
      temp_list.remove( next_edge_tmp );
   }
   else
   {
      return CUBIT_FAILURE;
   }

   while( temp_list.size() )
   {
      if( start_flg ) // Curr edge was found on prev start, so use end 1st
      {
	 vector_end = curr_edge->end_coordinates();
         next_edge = find_connected_ref_edge_by_coord( vector_end ,
                                                       temp_list, start_flg );
      }	 
      else
      {
	 vector_start = curr_edge->start_coordinates();
         next_edge = find_connected_ref_edge_by_coord( vector_start ,
                                                       temp_list, start_flg );
      }	 

      if( next_edge )
      {
         curr_edge = next_edge;
         ordered_list.append( next_edge );
         temp_list.remove( next_edge );
      }
      else
      {
         return CUBIT_FAILURE;
      }
   }

   ref_edge_list.clean_out();
   ref_edge_list = ordered_list;

   // Validate the chain
   DLIList<RefVertex*> problem_vertices;
   return check_valid_chain_by_coord( ref_edge_list, problem_vertices );
}

CubitStatus 
GeometryUtil::check_valid_chain_by_coord( DLIList<RefEdge*>& ref_edge_list, 
                                          DLIList<RefVertex*>& problem_vertices )
{
   CubitStatus status = CUBIT_SUCCESS;

   DLIList<RefEdge*> connected_edges;
   RefEdge* ref_edge_ptr;
   CubitVector vector_start, vector_end;
   
   for( int i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
   
      connected_edges.clean_out();
      vector_start = ref_edge_ptr->start_coordinates();
      find_connected_ref_edges_by_coord( vector_start ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 1 && connected_edges.size() != 2 )
      {
         status = CUBIT_FAILURE;
         problem_vertices.append_unique( ref_edge_ptr->start_vertex() );
      }
      
      connected_edges.clean_out();
      vector_end = ref_edge_ptr->end_coordinates();
      find_connected_ref_edges_by_coord( vector_end ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 1 && connected_edges.size() != 2 )
      {
         status = CUBIT_FAILURE;
         problem_vertices.append_unique( ref_edge_ptr->end_vertex() );
      }
   }
   return status;
}

CubitStatus 
GeometryUtil::check_valid_loop_by_coord( DLIList<RefEdge*>& ref_edge_list,
                                         DLIList<RefVertex*>& problem_vertices )
{
   CubitStatus status = CUBIT_SUCCESS;

   DLIList<RefEdge*> connected_edges;
   RefEdge* ref_edge_ptr;
   CubitVector vector_start, vector_end;
   
   for( int i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
   
      connected_edges.clean_out();
      vector_start = ref_edge_ptr->start_coordinates();
      find_connected_ref_edges_by_coord( vector_start ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
      {
         status = CUBIT_FAILURE;
         problem_vertices.append_unique( ref_edge_ptr->start_vertex() );
      }
      
      connected_edges.clean_out();
      vector_end = ref_edge_ptr->end_coordinates();
      find_connected_ref_edges_by_coord( vector_end ,
                                         ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
      {
         status = CUBIT_FAILURE;
         problem_vertices.append_unique( ref_edge_ptr->end_vertex() );
      }
   }
   return status;
}

// Added by CAT (SRS)
RefEdge* 
GeometryUtil::find_connected_ref_edge( RefVertex* ref_vertex_ptr, 
                                       DLIList<RefEdge*>& ref_edge_list,
                                       int& start_flg )
{
   RefEdge* ref_edge;
   DLIList<RefEdge*> tmp_list = ref_edge_list;

   tmp_list.reset();
   for( int i=0; i<tmp_list.size(); i++ )
   {
      ref_edge = tmp_list.get_and_step();

      if( ref_edge->start_vertex() == ref_vertex_ptr )
      {
         start_flg = 1;
         return ref_edge;
      }

      if( ref_edge->end_vertex() == ref_vertex_ptr )
      {
         start_flg = 0;
         return ref_edge;
      }
   }
   return NULL;
}

void
GeometryUtil::find_connected_ref_edges( RefVertex* ref_vertex_ptr, 
                                        DLIList<RefEdge*>& ref_edge_list,
                                        DLIList<RefEdge*>& connected_edges_list )
{
  RefEdge* ref_edge;
  DLIList<RefEdge*> tmp_list = ref_edge_list;
  
  tmp_list.reset();
  for( int i=0; i<tmp_list.size(); i++ )
  {
    ref_edge = tmp_list.get_and_step();
    
    if( ref_edge->start_vertex() == ref_vertex_ptr )
      connected_edges_list.append_unique( ref_edge );
    
    if( ref_edge->end_vertex() == ref_vertex_ptr )
      connected_edges_list.append_unique( ref_edge );
  }
}

CubitStatus 
GeometryUtil::form_ref_edge_loop( DLIList<RefEdge*>& ref_edge_list )
{
   int i;

   if( ref_edge_list.size() == 0 )
      return CUBIT_FAILURE;

   if( ref_edge_list.size() == 1 )
   {
      double period;
      if( ref_edge_list.get()->is_periodic( period ) )
         return CUBIT_SUCCESS;
   }

   // First form a chain (which can be a loop)
   if( form_ref_edge_chain( ref_edge_list ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   // Now check to make sure a valid loop was formed.  For a valid loop 
   // every vertex must have exactly 2 edges connected to it.  The only
   // exception is a periodic curve where there is only one vertex, which
   // was checked for before.
   DLIList<RefEdge*> connected_edges;
   RefEdge* ref_edge_ptr;
   for( i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
   
      connected_edges.clean_out();
      find_connected_ref_edges( ref_edge_ptr->start_vertex(), ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
         return CUBIT_FAILURE;

      connected_edges.clean_out();
      find_connected_ref_edges( ref_edge_ptr->end_vertex(), ref_edge_list, connected_edges );
      if( connected_edges.size() != 2 )
         return CUBIT_FAILURE;
   }
 
   return CUBIT_SUCCESS;
}

CubitStatus 
GeometryUtil::form_ref_edge_chain( DLIList<RefEdge*>& ref_edge_list )
{
   if( ref_edge_list.size() == 1 )
      return CUBIT_SUCCESS;

   int i;
   int start_flg = 0;
   
   RefVertex* ref_vertex_start, *ref_vertex_end;
   RefEdge* curr_edge;
   RefEdge* next_edge;

   DLIList<RefEdge*> temp_list = ref_edge_list;
   DLIList<RefEdge*> ordered_list;
   
   // Start vertex must be attached to only one ref_edge in our list. 
   // Otherwise this is really a loop, which is still okay

   // Search to find the starting edge
   RefEdge* start_edge = NULL;
   for( i=0; i<temp_list.size(); i++ )
   {
      curr_edge = temp_list.get();

      // Since DLList doesn't allow a NULL item
      DLIList<RefEdge*> temp_list2;
      temp_list2 = temp_list;
      temp_list2.remove( curr_edge ); // Next item becomes current item
   
      // get the vectors to start and end
      ref_vertex_start = curr_edge->start_vertex();
      ref_vertex_end   = curr_edge->end_vertex();
      
      if( !find_connected_ref_edge( ref_vertex_start,
         temp_list2, start_flg ) ||
         !find_connected_ref_edge( ref_vertex_end,
         temp_list2, start_flg ) )

      {
         start_edge = curr_edge;
         //temp_list.change_to( curr_edge );
         break;
      }
   }

   temp_list.reset();

   if( start_edge == NULL )
      start_edge = ref_edge_list.get();
   
   curr_edge = start_edge;
   ordered_list.append( curr_edge );
   temp_list.move_to( start_edge );
   temp_list.remove();

   next_edge = temp_list.get();

   // Attempt to use second edge in list as second edge in chain to 
   // preserve original chain direction.
   RefEdge* next_edge_tmp = NULL;
   RefEdge* next_edge_tmp2 = NULL;
   // Check start first - if this fails check end next so chain order
   // will be established properly.
   ref_vertex_start =  curr_edge->start_vertex();
   next_edge_tmp = find_connected_ref_edge( ref_vertex_start,
                                            temp_list, start_flg );
   if( next_edge_tmp != next_edge )
   {
      ref_vertex_end    = curr_edge->end_vertex();
      next_edge_tmp2 = find_connected_ref_edge( ref_vertex_end,
                                                temp_list, start_flg );
   }

   if( next_edge_tmp2 )
     next_edge_tmp = next_edge_tmp2;

   if( next_edge_tmp  )
   {
      curr_edge = next_edge_tmp;
      ordered_list.append( next_edge_tmp );
      temp_list.remove( next_edge_tmp );
   }
   else
   {
      return CUBIT_FAILURE;
   }

   while( temp_list.size() )
   {
      if( start_flg ) // Curr edge was found on prev start, so use end 1st
      {
	       ref_vertex_end = curr_edge->end_vertex();
         next_edge = find_connected_ref_edge( ref_vertex_end, temp_list, start_flg );
      }	 
      else
      {
	       ref_vertex_start = curr_edge->start_vertex();
         next_edge = find_connected_ref_edge( ref_vertex_start, temp_list, start_flg );
      }	 

      if( next_edge )
      {
         curr_edge = next_edge;
         ordered_list.append( next_edge );
         temp_list.remove( next_edge );
      }
      else
      {
         return CUBIT_FAILURE;
      }
   }

   ref_edge_list.clean_out();
   ref_edge_list = ordered_list;

   // Validate the chain
   DLIList<RefVertex*> problem_vertices;
   return check_valid_chain( ref_edge_list, problem_vertices );
}

CubitStatus 
GeometryUtil::check_valid_chain( DLIList<RefEdge*>& ref_edge_list, 
                                 DLIList<RefVertex*>& problem_vertices )
{
  CubitStatus status = CUBIT_SUCCESS;
  
  DLIList<RefEdge*> connected_edges;
  RefEdge* ref_edge_ptr;
  RefVertex* vertex_start, *vertex_end;
  
  for( int i=0; i<ref_edge_list.size(); i++ )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();
    
    connected_edges.clean_out();
    vertex_start = ref_edge_ptr->start_vertex();
    find_connected_ref_edges( vertex_start ,
      ref_edge_list, connected_edges );
    if( connected_edges.size() != 1 && connected_edges.size() != 2 )
    {
      status = CUBIT_FAILURE;
      problem_vertices.append_unique( ref_edge_ptr->start_vertex() );
    }
    
    connected_edges.clean_out();
    vertex_end = ref_edge_ptr->end_vertex();
    find_connected_ref_edges( vertex_end ,
      ref_edge_list, connected_edges );
    if( connected_edges.size() != 1 && connected_edges.size() != 2 )
    {
      status = CUBIT_FAILURE;
      problem_vertices.append_unique( ref_edge_ptr->end_vertex() );
    }
  }
  return status;
}

CubitStatus 
GeometryUtil::check_valid_loop( DLIList<RefEdge*>& ref_edge_list,
                                DLIList<RefVertex*>& problem_vertices )
{
  CubitStatus status = CUBIT_SUCCESS;
  
  DLIList<RefEdge*> connected_edges;
  RefEdge* ref_edge_ptr;
  RefVertex* vertex_start;
  RefVertex* vertex_end;
  
  for( int i=0; i<ref_edge_list.size(); i++ )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();
    
    connected_edges.clean_out();
    vertex_start = ref_edge_ptr->start_vertex();
    find_connected_ref_edges( vertex_start ,
      ref_edge_list, connected_edges );
    if( connected_edges.size() != 2 )
    {
      status = CUBIT_FAILURE;
      problem_vertices.append_unique( ref_edge_ptr->start_vertex() );
    }
    
    connected_edges.clean_out();
    vertex_end = ref_edge_ptr->end_vertex();
    find_connected_ref_edges( vertex_end ,
      ref_edge_list, connected_edges );
    if( connected_edges.size() != 2 )
    {
      status = CUBIT_FAILURE;
      problem_vertices.append_unique( ref_edge_ptr->end_vertex() );
    }
  }
  return status;
}
  
//-------------------------------------------------------------------------
// Purpose       : check vertex order and tangent direction for RefEdge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/07/00
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::valid_edge( RefEdge* edge_ptr,
                                       CubitBoolean print_error )
{
  RefVertex* start_vertex = edge_ptr->start_vertex();
  RefVertex* end_vertex = edge_ptr->end_vertex();
  CubitBoolean result = CUBIT_TRUE;
  
  if (!edge_ptr->get_curve_ptr())
  {
    if (print_error)
      PRINT_ERROR("Curve %d does not have a TopologyBridge.\n", edge_ptr->id());
    return CUBIT_FALSE;
  }
  else if(!start_vertex)
  {
    if (print_error)
      PRINT_ERROR("Curve %d has no start vertex.\n", edge_ptr->id());
    return CUBIT_FALSE;
  }
  else if(!end_vertex)
  {
    if (print_error)
      PRINT_ERROR("Curve %d has no end vertex.\n", edge_ptr->id());
    return CUBIT_FALSE;
  }
  else if (!start_vertex->get_point_ptr())
  {
    if (print_error)
      PRINT_ERROR("Vertex %d has no Point.\n", start_vertex->id());
    return CUBIT_FALSE;
  }
  else if (!end_vertex->get_point_ptr())
  {
    if (print_error)
      PRINT_ERROR("Vertex %d has no Point.\n", end_vertex->id());
    return CUBIT_FALSE;
  }
  else if (edge_ptr->get_curve_ptr()->geometry_type() == POINT_CURVE_TYPE)
  {
    if (start_vertex != end_vertex)
    {
      if (print_error) 
        PRINT_ERROR("Curve %d is a point-curve but has different "
                     "vertices (%d and %d).\n", edge_ptr->id(),
                     start_vertex->id(), end_vertex->id());
      result = CUBIT_FALSE;
    }
  }
  else
  {
    double u_start = edge_ptr->start_param();
    double u_end = edge_ptr->end_param();
    double period;
    CubitBoolean periodic = edge_ptr->is_periodic( period );
    if( start_vertex == end_vertex )
    {
      if( !periodic && print_error )
        PRINT_WARNING("Curve %d is closed but not periodic.\n",edge_ptr->id());
    }

    double u_start_vtx = edge_ptr->u_from_position( start_vertex->coordinates() );
    double u_end_vtx = edge_ptr->u_from_position( end_vertex->coordinates() );

    if (periodic)
    {
      if (period > 0.0)
      {
        while (u_start_vtx > period)
          u_start_vtx -= period;
        while (u_start > period)
          u_start -= period;
        while (u_end_vtx > period)
          u_end_vtx -= period;
        while (u_end > period)
          u_end -= period;
      }
      else if (period < 0.0)
      {
        while (u_start_vtx < period)
          u_start_vtx += period;
        while (u_start < period)
          u_start += period;
        while (u_end_vtx < period)
          u_end_vtx += period;
        while (u_end < period)
          u_end += period;
      }
    }
      


    if( (fabs( u_start_vtx - u_start ) > CUBIT_RESABS) && print_error )
      PRINT_WARNING("Vertex %d does not appear to be at the start of curve %d.\n",
          start_vertex->id(), edge_ptr->id() );
    if( (fabs( u_end_vtx - u_end ) > CUBIT_RESABS) && print_error )
      PRINT_WARNING("Vertex %d does not appear to be at the end of curve %d.\n",
          end_vertex->id(), edge_ptr->id() );

    double u_diff = u_start - u_end;
    double u_diff_vtx = u_start_vtx - u_end_vtx;
    if( (u_diff * u_diff_vtx) < 0 )
    {
      if( print_error )
      {
        PRINT_ERROR("Vertices on curve %d appear to be in the wrong order.\n",
          edge_ptr->id());
        PRINT_INFO("\tstart_param = %f, end_param = %f, u(start_vtx) = %f, "
                   "u(end_vtx) = %f\n",u_start,u_end,u_start_vtx,u_end_vtx);
      }
      result = CUBIT_FALSE;
    }

    double step = (u_end - u_start) / 100.;
    if( step > 0.0 )
    {
      if( step < CUBIT_RESABS ) step = CUBIT_RESABS;
      if( (u_start + step) > u_end ) step = u_end - u_start;
    }
    else
    {
      if( step > -CUBIT_RESABS ) step = -CUBIT_RESABS;
      if( (u_start + step) < u_end ) step = u_end - u_start;
    }

    CubitVector closest, tangent, next_pos;
    edge_ptr->closest_point( start_vertex->coordinates(), closest, &tangent );
    /* This fails a lot with spline curves, so skip it.
    double distance = (start_vertex->coordinates() - closest).length();
    if( (distance > GEOMETRY_RESABS) && print_error )
    {
      PRINT_WARNING("Vertex %d of curve %d does not appear to lie on the curve.\n",
        start_vertex->id(), edge_ptr->id() );
    }
    */
    edge_ptr->position_from_u( u_start + step, next_pos );
    next_pos -= closest;
    if( tangent.length_squared() < CUBIT_RESABS )
    {
      if( print_error )
        PRINT_ERROR("Zero-length tangent vector returned at start of curve %d\n", edge_ptr->id());
      result = CUBIT_FALSE;
    }
    else if( next_pos.length_squared() < CUBIT_RESABS )
    {
      if( print_error )
        PRINT_WARNING("Distance for delta-u of %f in [%f,%f] near start"
                      " of curve %d appears to be zero.\n",
                      step, u_start, u_end, edge_ptr->id());
    }
    else if( tangent.interior_angle( next_pos ) > (90. - CUBIT_RESABS) )
    {
      if( print_error )
        PRINT_ERROR("The tangent at the start of curve %d appears "
                    "to be in the wrong direction.\n",edge_ptr->id());
      result = CUBIT_FALSE;
    }

    edge_ptr->closest_point( end_vertex->coordinates(), closest, &tangent );
    /* This fails a lot with spline curves, so skip it.
    distance = (end_vertex->coordinates() - closest).length();
    if( (distance > GEOMETRY_RESABS) && print_error )
    {
      PRINT_WARNING("Vertex %d of curve %d does not appear to lie on the curve.\n",
        end_vertex->id(), edge_ptr->id() );
    }
    */
    edge_ptr->position_from_u( u_end - step, next_pos );
    next_pos = closest - next_pos;
    if( tangent.length_squared() < CUBIT_RESABS )
    {
      if( print_error )
        PRINT_ERROR("Zero-length tangent vector returned at end of %d\n", edge_ptr->id());
      result = CUBIT_FALSE;
    }
    else if( next_pos.length_squared() < CUBIT_RESABS )
    {
      if( print_error )
        PRINT_WARNING("Distance for delta-u of %f in [%f,%f] near end"
                      " of curve %d appears to be zero.\n",
                      step, u_start, u_end, edge_ptr->id());
    }
    else if( tangent.interior_angle( next_pos ) > (90. - CUBIT_RESABS) )
    {
      if( print_error )
        PRINT_ERROR("The tangent at the end of curve %d appears "
                    "to be in the wrong direction.\n",edge_ptr->id());
      result = CUBIT_FALSE;
    }
  }
  
  return result;
}

CubitBoolean GeometryUtil::valid_loop_coedges( Loop* loop_ptr,
                                               CubitBoolean print_error )
{
  CubitBoolean result = CUBIT_TRUE;
  RefFace* loop_owner = loop_ptr->get_ref_face_ptr();
  int surf_id = loop_owner ? loop_owner->id() : 0;

  if( !loop_ptr->get_loop_sm_ptr() && print_error )
  {
    PRINT_WARNING("Loop in Surface %d has no OSME pointer.\n",surf_id);
  }
  
  DLIList<CoEdge*> coedge_list;
  CubitStatus status = loop_ptr->ordered_co_edges( coedge_list );
  if( ! status )
  {
    if( print_error ) 
      PRINT_ERROR("Query for loop coedges failed!\n");
    return CUBIT_FALSE;
  }
  
  if( coedge_list.size() < 2 ) return CUBIT_TRUE;
  
  CoEdge* prev_coedge = coedge_list.get_and_step();
  RefEdge* prev_edge = prev_coedge->get_ref_edge_ptr();
  for( int i = coedge_list.size(); i > 0; i-- )
  {
    CoEdge* curr_coedge = coedge_list.get_and_step();
    RefEdge* curr_edge = curr_coedge->get_ref_edge_ptr();
    
    if (!curr_edge)
    {
      if (print_error)
        PRINT_ERROR("CoEdge in surface %d has no Curve.\n", surf_id);
      result = CUBIT_FAILURE;
      continue;
    }
    
    if( !curr_coedge->get_co_edge_sm_ptr() && print_error )
    {
      PRINT_WARNING("CoEdge between surface %d and curve %d has "
                    "no OSME pointer.\n", surf_id, curr_edge->id() );
    }
    
    RefVertex* common = curr_edge->common_ref_vertex( prev_edge );
    if( !common )
    {
      if( print_error )
        PRINT_ERROR("Curves %d and %d are ajacent in a Loop on surface %d, "
                    "but do not have a common vertex.\n",
          prev_edge->id(), curr_edge->id(), surf_id );
      result = CUBIT_FALSE;
    }
    
    RefVertex* start_vtx = 0;
    CubitSense sense = curr_coedge->get_sense();
    switch( sense )
    {
      case CUBIT_FORWARD:
        start_vtx = curr_edge->start_vertex();
        break;
      case CUBIT_REVERSED:
        start_vtx = curr_edge->end_vertex();
        break;
      case CUBIT_UNKNOWN:
      default:
        start_vtx = 0;
        break;
    }
    
    RefVertex* end_vtx = 0;
    switch( prev_coedge->get_sense() )
    {
      case CUBIT_FORWARD:
        end_vtx = prev_edge->end_vertex();
        break;
      case CUBIT_REVERSED:
        end_vtx = prev_edge->start_vertex();
        break;
      default:
        break;
    }
    
    if( start_vtx == NULL )
    {
      if( print_error )
        PRINT_ERROR("Curve %d has an unknown sense on surface %d.\n",
          curr_edge->id(), surf_id );
      result = CUBIT_FALSE;
    }
    else if( start_vtx != end_vtx )
    {
      if( print_error )
        PRINT_ERROR("The sense of curve %d (%s) on surface %d appears "
                    "to be incorrect with respect to the previous curve "
                    "in the loop (curve %d).\n", curr_edge->id(),
                    sense == CUBIT_FORWARD ? "forward" : "reverse",
                    surf_id, prev_edge->id() );
      result = CUBIT_FALSE;
    }
  
    prev_edge = curr_edge;
    prev_coedge = curr_coedge;
  }  
  
  return result;
}

CubitBoolean GeometryUtil::valid_shell_cofaces( Shell* shell_ptr,
                                                CubitBoolean print_error )
{
  DLIList<CoFace*> coface_list, ajacent_cofaces;
  DLIList<CoEdge*> coedge_list;
  shell_ptr->co_faces( coface_list );
  CubitBoolean result = CUBIT_TRUE;
  RefVolume* vol_ptr = shell_ptr->get_ref_volume_ptr();
  int vol_id = vol_ptr ? vol_ptr->id() : 0;
  
  if( ! shell_ptr->get_shell_sm_ptr() && print_error)
  {
    PRINT_WARNING("Shell in volume %d has no OSME pointer.\n",vol_id);
  }
      
  
  while( coface_list.size() > 0 )
  {
    CoFace* coface_ptr = coface_list.pop();
    RefFace* face_ptr = coface_ptr->get_ref_face_ptr();
    if( !face_ptr )
    {
      if( print_error )
        PRINT_ERROR("Encountered a CoFace without a RefFace in volume %d.\n",vol_id);
      result = CUBIT_FALSE;
      continue;
    }
    if (face_ptr->is_nonmanifold(shell_ptr))
    {
      continue;
    }
    CubitSense sense = coface_ptr->get_sense();
    if( sense == CUBIT_UNKNOWN )
    {
      if( print_error )
        PRINT_ERROR("Surface %d has an unknown sense with respect to volume %d.\n",
          face_ptr->id(), vol_id );
      result = CUBIT_FALSE;
      continue;
    }
 
    
    coedge_list.clean_out();
    face_ptr->co_edges( coedge_list );
    
    for( int i = coedge_list.size(); i > 0; i-- )
    {
      CoEdge* coedge_ptr = coedge_list.get_and_step();
      RefEdge* edge_ptr = coedge_ptr->get_ref_edge_ptr();
      if( !edge_ptr )
      {
        result = CUBIT_FALSE;
        if( print_error )
          PRINT_ERROR("CoEdge @ %p in Surface %d has no Curve!\n",
            coedge_ptr, face_ptr->id() );
        continue;
      }
      
      ajacent_cofaces.clean_out();
      edge_ptr->co_faces( ajacent_cofaces );
      ajacent_cofaces.intersect( coface_list );
      
      CubitSense coedge_sense = coedge_ptr->get_sense();
      if( coedge_sense == CUBIT_UNKNOWN )
      {
        if( print_error )
          PRINT_ERROR("CoEdge connecting Curve %d to surface %d "
                      "has UNKNOWN sense.\n",edge_ptr->id(), face_ptr->id() );
        result = CUBIT_FALSE;
        continue;
      }
      
      if( (ajacent_cofaces.size() > 1) && print_error )
      {
        PRINT_WARNING("Non-manifold topology at curve %d on volume %d "
                      "may result in false errors being reported.\n",
                      edge_ptr->id(), vol_id );
      }
      
      for( int j = ajacent_cofaces.size(); j > 0; j-- )
      {
        CoFace* other_coface = ajacent_cofaces.get_and_step();
        RefFace* other_face = other_coface->get_ref_face_ptr();
        if( ! other_face )
        {
          if( print_error )
            PRINT_ERROR("Encountered a CoFace on volume %d that does not "
                        "have a RefFace!\n",vol_id);
          result = CUBIT_FALSE;
          continue;
        }
        CubitSense other_coedge_sense = edge_ptr->sense( other_face );
        if( other_coedge_sense == CUBIT_UNKNOWN )
        {
          if( print_error )
            PRINT_ERROR("Curve %d has unknown sense with respect to "
                        "surface %d.\n",edge_ptr->id(), other_face->id() );
          result = CUBIT_FALSE;
          continue;
        }
        CubitSense other_sense = other_coface->get_sense();
        if( other_sense == CUBIT_UNKNOWN )
        {
          if( print_error )
            PRINT_ERROR("Surface %d has an unknown sense with respect "
                        "to volume %d.\n", other_face->id(), vol_id );
          result = CUBIT_FALSE;
          continue;
        }
        
        if( ((sense == other_sense) && (coedge_sense == other_coedge_sense))
          ||((sense != other_sense) && (coedge_sense != other_coedge_sense)) )
        {
          if( print_error )
            PRINT_ERROR("Incompatible CoFace senses found at curve %d in "
                        "volume %d.\n"
                        "  Surface  CoFace Sense  CoEdge Sense\n"
                        "  -------  ------------  ------------\n"
                        "  %7d  %12s  %12s\n  %7d  %12s  %12s\n",
                        edge_ptr->id(), vol_id,
                        face_ptr->id(), 
                        sense == CUBIT_FORWARD ? "forward" : "reverse",
                        coedge_sense == CUBIT_FORWARD ? "forward" : "reverse",
                        other_face->id(),
                        other_sense == CUBIT_FORWARD ? "forward" : "reverse",
                        other_coedge_sense == CUBIT_FORWARD ? 
                                       "forward" : "reverse");
          result = CUBIT_FALSE;
        }
      } // end for( j )
    
    } // end for( i )
  
  } // end while( coface_list.size() )
  
  return result;
}

CubitBoolean GeometryUtil::valid_topology(  TopologyEntity* topo_ptr,
                               CubitBoolean print_error,
                               DLIList<TopologyEntity*>* invalid_list )
{
  CubitBoolean result = CUBIT_TRUE;
  DLIList<Shell*> shell_list;
  DLIList<Loop*> loop_list;
  DLIList<RefEdge*> edge_list;
  Shell* shell_ptr;
  Loop* loop_ptr;
  RefEdge* edge_ptr;
  int i;  

  if( CAST_TO( topo_ptr, Body ) ||
      CAST_TO( topo_ptr, CoVolume ) ||
      CAST_TO( topo_ptr, RefVolume ) )
  {
    shell_list.clean_out();
    topo_ptr->shells( shell_list );
    for( i = shell_list.size(); i > 0; i-- )
    {
      shell_ptr = shell_list.get_and_step();
      if( ! valid_shell_cofaces( shell_ptr, print_error ) )
      {
        if( invalid_list ) invalid_list->append( shell_ptr );
        result = CUBIT_FALSE;
      }
    }
  }
  if(CAST_TO( topo_ptr, Body ) ||
     CAST_TO( topo_ptr, CoVolume ) ||
     CAST_TO( topo_ptr, RefVolume ) ||
     CAST_TO( topo_ptr, Shell ) ) 
  {
    shell_ptr = CAST_TO( topo_ptr, Shell );
    if( shell_ptr && !valid_shell_cofaces( shell_ptr, print_error ) )
    {
      if( invalid_list ) invalid_list->append( shell_ptr );
      result = CUBIT_FALSE;
    }
  }
  if(CAST_TO( topo_ptr, Body ) ||
     CAST_TO( topo_ptr, CoVolume ) ||
     CAST_TO( topo_ptr, RefVolume ) ||
     CAST_TO( topo_ptr, Shell ) ||
     CAST_TO( topo_ptr, CoFace ) ||
     CAST_TO( topo_ptr, RefFace ) )
  {
    loop_list.clean_out();
    topo_ptr->loops( loop_list );
    for( i = loop_list.size(); i > 0; i-- )
    {
      loop_ptr = loop_list.get_and_step();
      if( ! valid_loop_coedges( loop_ptr, print_error ) )
      {
        if( invalid_list ) invalid_list->append( loop_ptr );
        result = CUBIT_FALSE;
      }
    }
  }
  if(CAST_TO( topo_ptr, Body ) ||
     CAST_TO( topo_ptr, CoVolume ) ||
     CAST_TO( topo_ptr, RefVolume ) ||
     CAST_TO( topo_ptr, Shell ) ||
     CAST_TO( topo_ptr, CoFace ) ||
     CAST_TO( topo_ptr, RefFace ) ||
     CAST_TO( topo_ptr, Loop ) )
  {
    loop_ptr = CAST_TO( topo_ptr, Loop );
    if( loop_ptr && !valid_loop_coedges( loop_ptr, print_error ) )
    {
      if( invalid_list ) invalid_list->append( loop_ptr );
      result = CUBIT_FALSE;
    }
  }
  if(CAST_TO( topo_ptr, Body ) ||
     CAST_TO( topo_ptr, CoVolume ) ||
     CAST_TO( topo_ptr, RefVolume ) ||
     CAST_TO( topo_ptr, Shell ) ||
     CAST_TO( topo_ptr, CoFace ) ||
     CAST_TO( topo_ptr, RefFace ) ||
     CAST_TO( topo_ptr, Loop ) ||
     CAST_TO( topo_ptr, CoEdge ) )
  {
    edge_list.clean_out();
    topo_ptr->ref_edges( edge_list );
    for( i = edge_list.size(); i > 0; i-- )
    {
      edge_ptr = edge_list.get_and_step();
      if( !valid_edge( edge_ptr, print_error ) )
      {
        if( invalid_list ) invalid_list->append( edge_ptr );
        result = CUBIT_FALSE;
      }
    }
  }
  if(CAST_TO( topo_ptr, Body ) ||
     CAST_TO( topo_ptr, CoVolume ) ||
     CAST_TO( topo_ptr, RefVolume ) ||
     CAST_TO( topo_ptr, Shell ) ||
     CAST_TO( topo_ptr, CoFace ) ||
     CAST_TO( topo_ptr, RefFace ) ||
     CAST_TO( topo_ptr, Loop ) ||
     CAST_TO( topo_ptr, CoEdge ) ||
     CAST_TO( topo_ptr, RefEdge ) )
  {
    edge_ptr = CAST_TO( topo_ptr, RefEdge );
    if( edge_ptr && !valid_edge( edge_ptr, print_error ) )
    {
      if( invalid_list ) invalid_list->append( edge_ptr );
      result = CUBIT_FALSE;
    }
  }
  else if( CAST_TO( topo_ptr, Chain ) ||
           CAST_TO( topo_ptr, CoVertex ) ||
           CAST_TO( topo_ptr, RefVertex ) )
     ;//Do nothing
  else
  {
    if( print_error )
       PRINT_WARNING("Unknown entity type '%s' passed to GeometryUtil::"
                     "valid_topology(..).\n",topo_ptr->class_name() );
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Compare VGI topology to SM topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/19/01
//-------------------------------------------------------------------------
CubitBoolean GeometryUtil::valid_sm_topology( DLIList<RefEntity*>& entity_list,
                                              CubitBoolean print_error )
{
  CubitBoolean result = CUBIT_TRUE;
  
  DLIList<Body*>     body_list;
  DLIList<RefVolume*> vol_list, child_vols ;
  DLIList<RefFace*>  face_list, child_faces;
  DLIList<RefEdge*>  edge_list, child_edges;
  DLIList<RefVertex*> vtx_list, child_vtx;
  CAST_LIST( entity_list, body_list, Body     );
  CAST_LIST( entity_list, vol_list , RefVolume);
  CAST_LIST( entity_list, face_list, RefFace  );
  CAST_LIST( entity_list, edge_list, RefEdge  );
  CAST_LIST( entity_list, vtx_list , RefVertex);
  
  int i;
  for( i = body_list.size(); i--; )
  {
    Body* body_ptr = body_list.get_and_step();
    if( !valid_sm_topology( body_ptr, print_error ) )
      result = CUBIT_FALSE;
    child_vols.clean_out();
    body_ptr->ref_volumes( child_vols );
    vol_list.merge_unique( child_vols );
  }
  for( i = vol_list.size(); i--; )
  {
    RefVolume* vol_ptr = vol_list.get_and_step();
    if( !valid_sm_topology( vol_ptr, print_error ) )
      result = CUBIT_FALSE;
    child_faces.clean_out();
    vol_ptr->ref_faces( child_faces );
    face_list.merge_unique( child_faces );
  }
  for( i = face_list.size(); i--; )
  {
    RefFace* face_ptr = face_list.get_and_step();
    if( !valid_sm_topology( face_ptr, print_error ) )
      result = CUBIT_FALSE;
    child_edges.clean_out();
    face_ptr->ref_edges( child_edges );
    edge_list.merge_unique( child_edges );
  }
  for( i = edge_list.size(); i--; )
  {
    RefEdge* edge_ptr = edge_list.get_and_step();
    if( !valid_sm_topology( edge_ptr, print_error ) )
      result = CUBIT_FALSE;
    child_vtx.clean_out();
    edge_ptr->ref_vertices( child_vtx );
    vtx_list.merge_unique( child_vtx );
  }
  for( i = vtx_list.size(); i--; )
    if( !valid_sm_topology(vtx_list.get_and_step(), print_error ) )
      result = CUBIT_FALSE;
  
  return result;
}
CubitBoolean GeometryUtil::valid_sm_topology( Body* body_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This function checks for extraneous links in the VGI between
    // the passed body and its child volumes.  To check for missing
    // links, use valid_sm_topology(RefVolume*).
  
  DLIList<RefVolume*> vol_list;
  body_ptr->ref_volumes( vol_list );
  
    // Get all child TopologyBridges
  DLIList<TopologyBridge*> body_bridges, lump_bridges, temp_list;
  body_ptr->bridge_manager()->get_bridge_list( body_bridges );
  
  if ( body_bridges.size() != 1 )
  {
    PRINT_ERROR("%s (Body %d) has %d attached BodySM(s).\n", 
      body_ptr->entity_name().c_str(), body_ptr->id(), body_bridges.size() ); 
    result = CUBIT_FALSE;
  }

  for( i = body_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* bodysm = body_bridges.get_and_step();
    temp_list.clean_out();
    bodysm->get_children( temp_list );
    lump_bridges.merge_unique( temp_list );
  }
  
    // Check that each child VGI RefVolume has atleast one
    // TopologyBridge in the list we got from the solid modeler.
  for( i = vol_list.size(); i > 0; i-- )
  {
    RefVolume* vol_ptr = vol_list.get_and_step();
    temp_list.clean_out();
    vol_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( lump_bridges );
    if( temp_list.size() == 0 )
    {
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between %s and %s.\n",
        body_ptr->entity_name().c_str(),
        vol_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = lump_bridges.size(); i--; )
  {
    TopologyBridge* lump = lump_bridges.get_and_step();
    TopologyEntity* owner = lump->topology_entity();
    RefVolume* vol_ptr = dynamic_cast<RefVolume*>(owner);
    if( ! vol_ptr )
    {
      if( print ) PRINT_ERROR(
        "No RefVolume for Lump %p in %s\n",
        lump, body_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  } 

  return result;
}

CubitBoolean GeometryUtil::valid_sm_topology( RefVolume* vol_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the passed RefVolume and its parent Body(s)
    //  b.) Extraneous links between the RefVolume and its child shell(s)
    //  c.) Invalid Shells by calling valid_sm_topology(Shell*)
  
  DLIList<Body*> body_list;
  vol_ptr->bodies( body_list );
  DLIList<Shell*> shell_list;
  vol_ptr->shells( shell_list );
  DLIList<TopologyBridge*> vol_bridges, shell_bridges, temp_list;
  vol_ptr->bridge_manager()->get_bridge_list( vol_bridges );
  
  if ( vol_bridges.size() != 1 )
  {
    PRINT_ERROR("%s (RefVolume %d) has %d attached Lump(s).\n", 
      vol_ptr->entity_name().c_str(), vol_ptr->id(), vol_bridges.size() ); 
    result = CUBIT_FALSE;
  }
  
  for( i = vol_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* lump = vol_bridges.get_and_step();
    
    // Check for missing links between this RefVolume and parent Body(s)
    temp_list.clean_out();
    lump->get_parents( temp_list );
    for( int j = temp_list.size(); j > 0; j-- )
    {
      TopologyBridge* bodysm = temp_list.get_and_step();
      TopologyEntity* te_ptr = bodysm->topology_entity();
      Body* body_ptr = CAST_TO( te_ptr, Body );
      CubitString body_name;
      const char* body_str = "PARENT BODY";
      if( body_ptr )
      {
        body_name = body_ptr->entity_name();
        body_str = body_name.c_str();
      }
      if( (body_ptr == NULL) || !body_list.is_in_list(body_ptr) )
      {
        if( print ) PRINT_ERROR(
          "Missing VGI link between %s and %s\n",
          body_str,
          vol_ptr->entity_name().c_str() );
        result = CUBIT_FALSE;
      }
    }
    
    // Build ShellSM list for next part
    temp_list.clean_out();
    lump->get_children( temp_list );
    shell_bridges.merge_unique( temp_list );
  }
  
  // Check for extraneous VGI links between this RefVolume and child Shell(s)
  for( i = shell_list.size(); i > 0; i-- )
  {
    Shell* shell_ptr = shell_list.get_and_step();
    temp_list.clean_out();
    shell_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( shell_bridges );
    if( temp_list.size() == 0 )
    {
      DLIList<RefFace*> shell_faces;
      shell_ptr->ref_faces( shell_faces );
      CubitString face_name;
      const char* face_str = "NO REFFACE";
      if( shell_faces.size() )
      {
        face_name = shell_faces.get()->entity_name();
        face_str = face_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between %s and Shell with %s.\n",
        vol_ptr->entity_name().c_str(),
        face_str );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = shell_bridges.size(); i--; )
  {
    TopologyBridge* shell = shell_bridges.get_and_step();
    TopologyEntity* owner = shell->topology_entity();
    Shell* shell_ptr = dynamic_cast<Shell*>(owner);
    if( ! shell_ptr )
    {
      if( print ) PRINT_ERROR(
        "No Shell for ShellSM %p in %s\n",
        shell, vol_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  } 
  
  for( i = shell_list.size(); i > 0; i-- )
    if( !valid_sm_topology( shell_list.get_and_step(), print ) )
      result = CUBIT_FALSE;
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology(Shell* shell_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the passed Shell and its parent RefVolume(s)
    //  b.) Extraneous links between the Shell and its child RefFace(s)
  
  DLIList<RefVolume*> vol_list;
  shell_ptr->ref_volumes( vol_list );
  DLIList<RefFace*> face_list;
  shell_ptr->ref_faces( face_list );
  DLIList<TopologyBridge*> shell_bridges, face_bridges, temp_list;
  shell_ptr->bridge_manager()->get_bridge_list( shell_bridges );
  
  if ( shell_bridges.size() != 1 )
  {
    PRINT_ERROR("Shell %p (in Volume %d) has %d attached ShellSM(s).\n", 
      shell_ptr,  shell_ptr->get_ref_volume_ptr() ? 
      shell_ptr->get_ref_volume_ptr()->id() : 0,
      shell_bridges.size() );
    result = CUBIT_FALSE;
  }
  
  for( i = shell_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* shellsm = shell_bridges.get_and_step();
    
    // Check for missing links between this Shell and parent RefVolume(s)
    temp_list.clean_out();
    shellsm->get_parents( temp_list );
    if (temp_list.size() != 1)
    {
      PRINT_ERROR("Bad SolidModel Topology!\n");
      PRINT_ERROR("ShellSM %p in Shell %p (in Volume %d) has %d parent Lumps.\n",
       shellsm, shell_ptr,  
       shell_ptr->get_ref_volume_ptr() ? shell_ptr->get_ref_volume_ptr()->id() : 0,
       temp_list.size() );
    }
    else
    {
      TopologyBridge* lump = temp_list.get_and_step();
      TopologyEntity* te_ptr = lump->topology_entity();
      RefVolume* vol_ptr = CAST_TO( te_ptr, RefVolume );
      if( (vol_ptr == NULL) || !vol_list.is_in_list(vol_ptr) )
      {
        DLIList<RefFace*> faces;
        shell_ptr->ref_faces( faces );
        CubitString vol_name, face_name;
        const char* vol_str = "PARENT VOLUME";
        const char* face_str = "NO REFFACES";
        if( vol_ptr )
        {
          vol_name = vol_ptr->entity_name();
          vol_str = vol_name.c_str();
        }
        if( faces.size() )
        {
          face_name = faces.get()->entity_name();
          face_str = face_name.c_str();
        }
        if( print ) PRINT_ERROR(
          "Missing VGI link between %s and Shell with %s\n",
          vol_str, face_str );
        result = CUBIT_FALSE;
      }
    }
    // Build Surface list for next part
    temp_list.clean_out();
    shellsm->get_children( temp_list );
    face_bridges.merge_unique( temp_list );
  }
  
  // Check for extraneous VGI links between this Shell and child RefFaces(s)
  for( i = face_list.size(); i > 0; i-- )
  {
    RefFace* face_ptr = face_list.get_and_step();
    temp_list.clean_out();
    face_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( face_bridges );
    if( temp_list.size() == 0 )
    {
      DLIList<RefFace*> faces;
      shell_ptr->ref_faces( faces );
      CubitString face_name;
      const char* face_str = "NO REFFACES";
      if( faces.size() )
      {
        face_name = faces.get()->entity_name();
        face_str = face_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between Shell with %s and %s.\n",
        face_str,
        face_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = face_bridges.size(); i--; )
  {
    TopologyBridge* surface = face_bridges.get_and_step();
    TopologyEntity* owner = surface->topology_entity();
    RefFace* face_ptr = dynamic_cast<RefFace*>(owner);
    if( ! face_ptr )
    {
      if( print ) PRINT_ERROR(
        "No RefFace for Surface %p in Shell %p in %s\n",
        surface, shell_ptr, 
        vol_list.size() ?
          vol_list.get()->entity_name().c_str() :
          "NO REFVOLUME" );
      result = CUBIT_FALSE;
    }
  } 
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology( RefFace* face_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the passed RefFace and its parent Shell(s)
    //  b.) Extraneous links between the RefFace and its child Loop(s)
    //  c.) Invalid Loops by calling valid_sm_topology(Loop*)
  
  DLIList<Shell*> shell_list;
  face_ptr->shells( shell_list );
  DLIList<Loop*> loop_list;
  face_ptr->loops( loop_list );
  DLIList<TopologyBridge*> face_bridges, loop_bridges, temp_list;
  face_ptr->bridge_manager()->get_bridge_list( face_bridges );
  
  if ( face_bridges.size() == 0 )
  {
    PRINT_ERROR("%s (RefFace %d) has no attached Surface(s).\n", 
      face_ptr->entity_name().c_str(), face_ptr->id() ); 
    result = CUBIT_FALSE;
  }
  
  for( i = face_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* surf = face_bridges.get_and_step();
    
    // Check for missing links between this RefFace and parent Shell(s)
    temp_list.clean_out();
    surf->get_parents( temp_list );
    for( int j = temp_list.size(); j > 0; j-- )
    {
      TopologyBridge* shellsm = temp_list.get_and_step();
      TopologyEntity* te_ptr = shellsm->topology_entity();
      Shell* shell_ptr = CAST_TO( te_ptr, Shell );
      if( (shell_ptr == NULL) || !shell_list.is_in_list(shell_ptr) )
      {
        DLIList<RefVolume*> vols;
        if( shell_ptr != NULL ) 
          shell_ptr->ref_volumes( vols );
        CubitString vol_name;
        const char* vol_str = "PARENT REFVOLUME";
        if( vols.size() )
        {
          vol_name = vols.get()->entity_name();
          vol_str = vol_name.c_str();
        }
        if( print ) PRINT_ERROR(
          "Missing VGI link between Shell with %s and %s\n",
          vol_str,
          face_ptr->entity_name().c_str() );
        result = CUBIT_FALSE;
      }
    }
    
    // Build LoopSM list for next part
    temp_list.clean_out();
    surf->get_children( temp_list );
    loop_bridges.merge_unique( temp_list );
  }
  
  // Check for extraneous VGI links between this RefFace and child Loop(s)
  for( i = loop_list.size(); i > 0; i-- )
  {
    Loop* loop_ptr = loop_list.get_and_step();
    temp_list.clean_out();
    loop_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( loop_bridges );
    if( temp_list.size() == 0 )
    {
      DLIList<RefEdge*> edges;
      loop_ptr->ref_edges( edges );
      CubitString edge_name;
      const char* edge_str = "NO REFEDGE";
      if( edges.size() )
      {
        edge_name = edges.get()->entity_name();
        edge_str = edge_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between %s and Loop with %s.\n",
        face_ptr->entity_name().c_str(),
        edge_str );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = loop_bridges.size(); i--; )
  {
    TopologyBridge* loop = loop_bridges.get_and_step();
    TopologyEntity* owner = loop->topology_entity();
    Loop* loop_ptr = dynamic_cast<Loop*>(owner);
    if( ! loop_ptr )
    {
      if( print ) PRINT_ERROR(
        "No Loop for LoopSM %p in %s\n",
        loop, face_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  } 
  
  for( i = loop_list.size(); i > 0; i-- )
    if( !valid_sm_topology( loop_list.get_and_step(), print ) )
      result = CUBIT_FALSE;
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology( Loop* loop_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  CoEdge* coedge_ptr;
  int i;

    // This method checks for 
    //  a.) Missing links between the passed Loop and its parent RefFace(s)
    //  b.) Extraneous links between the Loop and its child CoEdge(s)
    //  c.) Invalid CoEdges by calling valid_sm_topology(CoEdge*)
    //  d.) CoEdgeSMs w/out CoEdges
    //  e.) Correct CoEdge order given CoEdgeSM order.
  
  DLIList<RefFace*> face_list;
  loop_ptr->ref_faces( face_list );
  DLIList<CoEdge*> coedge_list;
  loop_ptr->co_edges( coedge_list );
  DLIList<TopologyBridge*> loop_bridges, coedge_bridges, temp_list;
  loop_ptr->bridge_manager()->get_bridge_list( loop_bridges );
  
  if ( loop_bridges.size() == 0 )
  {
    PRINT_ERROR("Loop %p (in surface %d) has no attached LoopSM(s).\n", 
      loop_ptr, 
      loop_ptr->get_ref_face_ptr() ? loop_ptr->get_ref_face_ptr()->id() : 0 );
    result = CUBIT_FALSE;
  }
  
  for( i = loop_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* loopsm = loop_bridges.get_and_step();
    
    // Check for missing links between this Loop and parent RefFace(s)
    temp_list.clean_out();
    loopsm->get_parents( temp_list );
    for( int j = temp_list.size(); j > 0; j-- )
    {
      TopologyBridge* surf = temp_list.get_and_step();
      TopologyEntity* te_ptr = surf->topology_entity();
      RefFace* face_ptr = CAST_TO( te_ptr, RefFace );
      if( (face_ptr == NULL) || !face_list.is_in_list(face_ptr) )
      {
        DLIList<RefEdge*> edges;
        loop_ptr->ref_edges( edges );
        CubitString face_name;
        CubitString edge_name;
        const char* face_str = "PARENT REFFACE";
        const char* edge_str = "NO REFEDGES";
        if( face_ptr )
        {
          face_name = face_ptr->entity_name();
          face_str = face_name.c_str();
        }
        if( edges.size() )
        {
          edge_name = edges.get()->entity_name();
          edge_str = edge_name.c_str();
        }
        if( print ) PRINT_ERROR(
          "Missing VGI link between %s and Loop with %s\n", face_str, edge_str);
        result = CUBIT_FALSE;
      }
    }
    
    // Build CoEdgeSM list for next part
    temp_list.clean_out();
    loopsm->get_children( temp_list );
    coedge_bridges.merge_unique( temp_list );
  }
  
  // Check for extraneous VGI links between this Loop and child CoEdge(s)
  for( i = coedge_list.size(); i > 0; i-- )
  {
    coedge_ptr = coedge_list.get_and_step();
    temp_list.clean_out();
    coedge_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( coedge_bridges );
    if( temp_list.size() == 0 )
    {
      DLIList<RefEdge*> edges;
      DLIList<RefFace*> faces;
      loop_ptr->ref_faces( faces );
      coedge_ptr->ref_edges( edges );
      CubitString face_name;
      CubitString edge_name;
      const char* face_str = "NO REFFACE";
      const char* edge_str = "NO REFEDGE";
      if( faces.size() )
      {
        face_name = faces.get()->entity_name();
        face_str = face_name.c_str();
      }
      if( edges.size() )
      {
        edge_name = edges.get()->entity_name();
        edge_str = edge_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between Loop in %s and CoEdge with %s.\n",
        face_str, edge_str );
      result = CUBIT_FALSE;
    }
  }
  
    // Check for missing CoedgeSMs
  for( i = coedge_bridges.size(); i--; )
  {
    TopologyBridge* coedge = coedge_bridges.get_and_step();
    TopologyEntity* owner = coedge->topology_entity();
    coedge_ptr = dynamic_cast<CoEdge*>(owner);
    if( ! coedge_ptr )
    {
      CubitString face_name = "NO REFFACE";
      if( face_list.size() )
        face_name = face_list.get()->entity_name();
      if( print ) PRINT_ERROR(
        "No CoEdge for CoEdgeSM %p in Loop %p in %s\n",
        coedge, loop_ptr, face_name.c_str() );
      result = CUBIT_FALSE;
    }
  } 
  
  
  for( i = coedge_list.size(); i > 0; i-- )
    if( !valid_sm_topology( coedge_list.get_and_step(), print ) )
      result = CUBIT_FALSE;

  if( !result )
    return result;

  
  RefFace* face_ptr = loop_ptr->get_ref_face_ptr();
  CubitString face_name_str;
  const char* face_str_ptr = "NO PARENT REFFACE";
  if( face_ptr )
  {
    face_name_str = face_ptr->entity_name();
    face_str_ptr = face_name_str.c_str();
  }
  
  
    // Check correct order of CoEdges
  for( i = loop_bridges.size(); i--; )
  {
    TopologyBridge* loopsm_ptr = loop_bridges.get_and_step();
    coedge_bridges.clean_out();
    loopsm_ptr->get_children( coedge_bridges );
  
    if( coedge_bridges.size() != coedge_list.size() )
    {
      if( print ) PRINT_ERROR(
        "Loop %p in %s has %d coedges, while LoopSM %p has %d coedgesms.\n",
          loop_ptr, face_str_ptr, coedge_list.size(),
          loopsm_ptr, coedge_bridges.size() );
      result = CUBIT_FAILURE;
      continue;
    }

      // Remaining code checks the order of coedges in the loop.
      // There is no wrong order if only one or two coedges.
    if( coedge_bridges.size() <= 2 )
      continue;
  
    coedge_bridges.reset();
    coedge_ptr = dynamic_cast<CoEdge*>(coedge_bridges.get()->topology_entity());
    if( !coedge_list.move_to( coedge_ptr ) )
      continue;

    coedge_ptr = dynamic_cast<CoEdge*>(coedge_bridges.next()->topology_entity());
    if( coedge_list.next() != coedge_ptr &&
        coedge_list.prev() != coedge_ptr )
    {
      if( print ) PRINT_ERROR(
        "Order of CoEdges in Loop %p in %s is incorrect w.r.t. LoopSM %p.\n",
        loop_ptr, face_str_ptr, loopsm_ptr);
      result = CUBIT_FAILURE;
      continue;
    }

    bool loop_reversed = false;
    TopologyBridge* coedge_bridge;
    if( coedge_list.prev() == coedge_ptr )
    {
      coedge_bridge = coedge_bridges.get();
      coedge_bridges.reverse();
      coedge_bridges.move_to(coedge_bridge);
      loop_reversed = true;
    }

    for( int j = coedge_bridges.size(); j--; )
    {
      coedge_bridge = coedge_bridges.get_and_step();
      coedge_ptr = coedge_list.get_and_step();
      if( coedge_bridge->topology_entity() != coedge_ptr )
      {
        if( print ) PRINT_ERROR(
          "Order of CoEdges in Loop %p in %s is incorrect w.r.t. loopsm %p.\n",
          loop_ptr, face_str_ptr,loopsm_ptr);
        result = CUBIT_FAILURE;
        break;
      }
    }
    
    temp_list.clean_out();
    loopsm_ptr->get_parents(temp_list);
    if (temp_list.size() != 1)
      continue;
    
    bool surf_reversed = (temp_list.get()->bridge_sense() == CUBIT_REVERSED);
    if (surf_reversed != loop_reversed)
    {
      result = CUBIT_FAILURE;
      if (print) PRINT_ERROR(
        "Loop %p in %s is reversed w.r.t. loopsm %p.\n",
        loop_ptr, face_str_ptr, loopsm_ptr);
    }
    
  }       
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology(CoEdge* coedge_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the passed CoEdge and its parent Loop(s)
    //  b.) Extraneous links between the Loop and its child RefEdge(s)
  
  DLIList<Loop*> loop_list;
  coedge_ptr->loops( loop_list );
  DLIList<RefEdge*> edge_list;
  coedge_ptr->ref_edges( edge_list );
  DLIList<TopologyBridge*> coedge_bridges, edge_bridges, temp_list;
  coedge_ptr->bridge_manager()->get_bridge_list( coedge_bridges );
  
  if ( coedge_bridges.size() == 0 )
  {
    if (print) PRINT_ERROR(
      "CoEdge %p (curve %d %s in surface %d) has no attached CoEdgeSM(s).\n", 
      coedge_ptr, 
      coedge_ptr->get_ref_edge_ptr() ? coedge_ptr->get_ref_edge_ptr()->id() : 0,
      coedge_ptr->get_sense() == CUBIT_FORWARD ? "FORWARD" :
      coedge_ptr->get_sense() == CUBIT_REVERSED ? "REVERSED" : "UNKNOWN",
      coedge_ptr->get_ref_face() ? coedge_ptr->get_ref_face()->id() : 0 );
    result = CUBIT_FALSE;
  }
  
  for( i = coedge_bridges.size(); i > 0; i-- )
  {
    CoEdgeSM* coedgesm = dynamic_cast<CoEdgeSM*>(coedge_bridges.get_and_step());
    
    // Check for missing links between this CoEdge and parent Loop(s)
    temp_list.clean_out();
    coedgesm->get_parents( temp_list );
    if (temp_list.size() != 1)
    {
      PRINT_ERROR(
        "Bad SolidModel topoloy.  CoEdgeSM %p in CoEdge %p "
        "(curve %d %s in surface %d) has %d parent LoopSMs.\n", 
        coedgesm, coedge_ptr,
        coedge_ptr->get_ref_edge_ptr() ? coedge_ptr->get_ref_edge_ptr()->id() : 0,
        coedge_ptr->get_sense() == CUBIT_FORWARD ? "FORWARD" :
        coedge_ptr->get_sense() == CUBIT_REVERSED ? "REVERSED" : "UNKNOWN",
        coedge_ptr->get_ref_face() ? coedge_ptr->get_ref_face()->id() : 0,
        temp_list.size() );
      continue;
    }
                 
    TopologyBridge* loopsm = temp_list.get_and_step();
    TopologyEntity* te_ptr = loopsm->topology_entity();
    Loop* loop_ptr = CAST_TO( te_ptr, Loop );
    if( !loop_ptr )
    {
      const char* edge_str = "NO REFEDGE";
      CubitString edge_name;
      DLIList<RefEdge*> edges;
      coedge_ptr->ref_edges( edges );
      if( edges.size() )
      {
        edge_name = edges.get()->entity_name();
        edge_str = edge_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Missing VGI Entity for LoopSM %p. Discovered on CoEdge with %s\n",
        loopsm, edge_str );
      result = CUBIT_FALSE;
    }
    else if( !loop_list.is_in_list(loop_ptr) )
    {
      DLIList<RefFace*> faces;
      DLIList<RefEdge*> edges;
      loop_ptr->ref_faces( faces );
      coedge_ptr->ref_edges( edges );
      CubitString face_name, edge_name;
      const char* face_str = "NO REFFACE";
      const char* edge_str = "NO REFEDGE";
      if( faces.size() )
      {
        face_name = faces.get()->entity_name();
        face_str = face_name.c_str();
      }
      if( edges.size() )
      {
        edge_name = edges.get()->entity_name();
        edge_str = edge_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Missing VGI link between Loop in %s and CoEdge with %s\n",
        face_str, edge_str );
      result = CUBIT_FALSE;
    }
    
    
    // Build curve list for next part
    temp_list.clean_out();
    coedgesm->get_children( temp_list );
    edge_bridges.merge_unique( temp_list );
    
    // Get curve and surface connected by this CoEdgeSM
    if (temp_list.size() != 1) // should never happen!
      continue;
    TopologyBridge* curv_bridge = temp_list.get();
    temp_list.clean_out();
    loopsm->get_parents(temp_list);
    if (temp_list.size() != 1) // should never happen!
      continue;
    TopologyBridge* surf_bridge = temp_list.get();
    
    bool curv_reversed = curv_bridge->bridge_sense() == CUBIT_REVERSED;
    bool surf_reversed = surf_bridge->bridge_sense() == CUBIT_REVERSED;
    bool want_reversed = curv_reversed != surf_reversed;
    bool  are_reversed = coedgesm->sense() != coedge_ptr->get_sense();
    if (want_reversed != are_reversed)
    {
      result = CUBIT_FAILURE;
      if (print) PRINT_ERROR(
        "CoEdge %p (curve %d %s in surface %d) has wrong sense w.r.t CoEdgeSM.\n", 
        coedge_ptr, 
        coedge_ptr->get_ref_edge_ptr() ? coedge_ptr->get_ref_edge_ptr()->id() : 0,
        coedge_ptr->get_sense() == CUBIT_FORWARD ? "FORWARD" :
        coedge_ptr->get_sense() == CUBIT_REVERSED ? "REVERSED" : "UNKNOWN",
        coedge_ptr->get_ref_face() ? coedge_ptr->get_ref_face()->id() : 0 );
    }    
  }
  
  // Check for extraneous VGI links between this CoEdge and child RefEdge(s)
  for( i = edge_list.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = edge_list.get_and_step();
    temp_list.clean_out();
    edge_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( edge_bridges );
    if( temp_list.size() == 0 )
    {
      DLIList<RefEdge*> edges;
      coedge_ptr->ref_edges( edges );
      CubitString edge_name;
      const char* edge_str = "NO REFEDGE";
      if( edges.size() )
      {
        edge_name = edges.get()->entity_name();
        edge_str = edge_name.c_str();
      }
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between CoEdge with %s and %s.\n",
        edge_str,
        edge_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = edge_bridges.size(); i--; )
  {
    TopologyBridge* curve = edge_bridges.get_and_step();
    TopologyEntity* owner = curve->topology_entity();
    RefEdge* edge_ptr = dynamic_cast<RefEdge*>(owner);
    if( ! edge_ptr )
    {
      DLIList<RefFace*> faces;
      coedge_ptr->ref_faces( faces );
      CubitString face_name = "NO REFFACE";
      if( faces.size() )
        face_name = faces.get()->entity_name();
      if( print ) PRINT_ERROR(
        "No RefEdge for Curve %p in CoEdge %p in %s\n",
        curve, coedge_ptr, face_name.c_str() );
      result = CUBIT_FALSE;
    }
  } 
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology( RefEdge* edge_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the passed RefEdge and its parent CoEdge(s)
    //  b.) Extraneous links between the RefEdge and its child vert(ex/ices)
  
  DLIList<CoEdge*> coedge_list;
  edge_ptr->co_edges( coedge_list );
  DLIList<RefVertex*> vtx_list;
  edge_ptr->ref_vertices( vtx_list );
  DLIList<TopologyBridge*> edge_bridges, vtx_bridges, temp_list;
  edge_ptr->bridge_manager()->get_bridge_list( edge_bridges );
  
  if ( edge_bridges.size() == 0 )
  {
    PRINT_ERROR("%s (RefEdge %d) has no attached Curve(s).\n", 
      edge_ptr->entity_name().c_str(), edge_ptr->id());
    result = CUBIT_FALSE;
  }
  
  for( i = edge_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* curve = edge_bridges.get_and_step();
    
    // Check for missing links between this RefEdge and parent CoEdge(s)
    temp_list.clean_out();
    curve->get_parents( temp_list );
    for( int j = temp_list.size(); j > 0; j-- )
    {
      TopologyBridge* coedgesm = temp_list.get_and_step();
      TopologyEntity* te_ptr = coedgesm->topology_entity();
      CoEdge* coedge_ptr = CAST_TO( te_ptr, CoEdge );
      if( (coedge_ptr == NULL) || !coedge_list.is_in_list(coedge_ptr) )
      {
        DLIList<RefFace*> faces;
        if( coedge_ptr ) 
          coedge_ptr->ref_faces( faces );
        CubitString face_name;
        const char* face_str = "NO REFFACE";
        if( faces.size() )
        {
          face_name = faces.get()->entity_name();
          face_str = face_name.c_str();
        }
        
        if( print ) PRINT_ERROR(
          "Missing VGI link between CoEdge with %s and %s\n",
          face_str,
          edge_ptr->entity_name().c_str() );
        result = CUBIT_FALSE;
      }
    }
    
    // Build Point list for next part
    temp_list.clean_out();
    curve->get_children( temp_list );
    vtx_bridges.merge_unique( temp_list );
  }
  
  // Check for extraneous VGI links between this RefEdge and child RefVert(s)
  for( i = vtx_list.size(); i > 0; i-- )
  {
    RefVertex* vtx_ptr = vtx_list.get_and_step();
    temp_list.clean_out();
    vtx_ptr->bridge_manager()->get_bridge_list( temp_list );
    temp_list.intersect( vtx_bridges );
    if( temp_list.size() == 0 )
    {
      if( print ) PRINT_ERROR(
        "Extraneous VGI link between %s and %s.\n",
        edge_ptr->entity_name().c_str(),
        vtx_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  }
  
  for( i = vtx_bridges.size(); i--; )
  {
    TopologyBridge* point = vtx_bridges.get_and_step();
    TopologyEntity* owner = point->topology_entity();
    RefVertex* vtx_ptr = dynamic_cast<RefVertex*>(owner);
    if( ! vtx_ptr )
    {
      if( print ) PRINT_ERROR(
        "No RefVertex for Point %p in %s\n",
        point, edge_ptr->entity_name().c_str() );
      result = CUBIT_FALSE;
    }
  } 
  
  return result;
}


CubitBoolean GeometryUtil::valid_sm_topology( RefVertex* vtx_ptr, CubitBoolean print )
{
  CubitBoolean result = CUBIT_TRUE;
  int i;
  
    // This method checks for 
    //  a.) Missing links between the RefVertex and its parent RefEdge(s)
  
  DLIList<RefEdge*> edge_list;
  vtx_ptr->ref_edges( edge_list );
  DLIList<TopologyBridge*> vtx_bridges, temp_list;
  vtx_ptr->bridge_manager()->get_bridge_list( vtx_bridges );
  
  if ( vtx_bridges.size() == 0 )
  {
    PRINT_ERROR("%s (RefVertex %d) has no attached Point(s).\n", 
      vtx_ptr->entity_name().c_str(), vtx_ptr->id());
    result = CUBIT_FALSE;
  }
  
  for( i = vtx_bridges.size(); i > 0; i-- )
  {
    TopologyBridge* point = vtx_bridges.get_and_step();
    
    // Check for missing links between this RefVertex and parent RefEdge(s)
    temp_list.clean_out();
    point->get_parents( temp_list );
    for( int j = temp_list.size(); j > 0; j-- )
    {
      TopologyBridge* curve = temp_list.get_and_step();
      TopologyEntity* te_ptr = curve->topology_entity();
      RefEdge* edge_ptr = CAST_TO( te_ptr, RefEdge );
      if( !edge_ptr )
      {
        if( print ) PRINT_ERROR(
          "Missing VGI Entity for Curve %p.  Discovered as parent of %s\n",
            curve, vtx_ptr->entity_name().c_str() );
        result = CUBIT_FAILURE;
      }
      else if( !edge_list.is_in_list(edge_ptr) )
      {
        if( print ) PRINT_ERROR(
          "Missing VGI link between %s and %s\n",
          edge_ptr->entity_name().c_str(),
          vtx_ptr->entity_name().c_str() );
        result = CUBIT_FALSE;
      }
    }
  }
  
  return result;
}

void GeometryUtil::list_SM_topology( TopologyBridge* bridge, int depth )
{
  int counts[8] = {0,0,0,0,0,0,0,0};
  const char* names[8] = {"BodySM","Lump","ShellSM", "Surface", "LoopSM",
                          "CoEdgeSM", "Curve", "Point"};
  list_SM_topology( bridge, depth, 0, counts );
  
  for( int i = 0; i < 8; i++ )
    if( counts[i] == 1 )
      PRINT_INFO("%d %s, ", counts[i], names[i] );
    else if( counts[i] > 1 )
      PRINT_INFO("%d %ss, ", counts[i], names[i] );
  PRINT_INFO("\n");
}

void GeometryUtil::list_SM_topology( TopologyBridge* bridge, int depth,
                                     int indent, int counts[8] )
{
  DLIList<TopologyBridge*> relatives;

  if( depth < 0 )
  {
    bridge->get_parents( relatives );
    relatives.reset();
    for( int i = relatives.size(); i--; )
      list_SM_topology( relatives.get_and_step(), depth+1, indent+1, counts );
  }

  int index = print_topo_bridge( bridge, indent );
  if( index >= 0 && index <= 8 )
    counts[index]++;
  
  if( depth > 0 )
  {
    bridge->get_children( relatives );
    relatives.reset();
    for( int i = relatives.size(); i--; )
      list_SM_topology( relatives.get_and_step(), depth-1, indent+1, counts );
  }
}

int GeometryUtil::print_topo_bridge( TopologyBridge* bridge, int indent )
{
  const char* label = bridge ? typeid(*bridge).name() : "(null)";
#ifdef __GNUC__
  while( isdigit(*label) ) label++;
#endif

  int index = -1;
  if( dynamic_cast<BodySM*>(bridge) )
    index = 0;
  else if( dynamic_cast<Lump*>(bridge) )
    index = 1;
  else if( dynamic_cast<ShellSM*>(bridge) )
    index = 2;
  else if( dynamic_cast<Surface*>(bridge) )
    index = 3;
  else if( dynamic_cast<LoopSM*>(bridge) )
    index = 4;
  else if( dynamic_cast<CoEdgeSM*>(bridge) )
    index = 5;
  else if( dynamic_cast<Curve*>(bridge) )
    index = 6;
  else if( dynamic_cast<Point*>(bridge) )
    index = 7;
  
  CubitString name("");
  RefEntity* re = 0;
  if( bridge && (re = dynamic_cast<RefEntity*>(bridge->topology_entity())) )
    name = re->entity_name();
  PRINT_INFO("%*s%s 0x%lx (%s)\n", indent*3, "", label, (long)bridge, name.c_str());
  
  return index;
}
