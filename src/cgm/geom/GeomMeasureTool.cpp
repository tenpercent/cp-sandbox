//----------------------------------------------------------
//Class: GeomMeasureTool
//Description: Does Geometric measureing and statistical
//             gathering on various topologies.
//Creator: David R. White
//Date: 7/9/2002
//----------------------------------------------------------


#include "GeomMeasureTool.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "Body.hpp"
#include "RefGroup.hpp"
#include "RefEntityFactory.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Chain.hpp"
#include "CastTo.hpp"
#include "GMem.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeomSeg.hpp"
#include "GeomPoint.hpp"
#include "RTree.hpp"
#include "AbstractTree.hpp"
#include "IntersectionTool.hpp"
#include "CpuTimer.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "MergeTool.hpp"
#include "GeomDataObserver.hpp"
#include "Lump.hpp"
#include "ModelQueryEngine.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"


void GeomMeasureTool::get_edges_from_list(DLIList <RefEntity*> &entity_list,
                                          DLIList <RefEdge*> &ref_edges )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefEdge *ref_edge;
  
  TopologyEntity *topo_ent;
  DLIList <RefEdge*> tmp_edges;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_edges.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    if ( topo_ent == NULL )
      continue;
    topo_ent->ref_edges(tmp_edges);
    for ( jj = tmp_edges.size(); jj > 0; jj-- )
    {
      ref_edge = tmp_edges.pop();
      if ( !ref_edge->marked() )
      {
        ref_edge->marked(1);
        ref_edges.append(ref_edge);
      }
    }
  }
  for ( ii = ref_edges.size(); ii > 0; ii-- )
    ref_edges.get_and_step()->marked(0);
}

void GeomMeasureTool::get_bodies_from_list(DLIList <RefVolume*> &entity_list,
                                           DLIList <Body*> &ref_bodies )
{
    //assume the mark is zero.
  int ii, jj;
  RefVolume *entity_ptr;
  Body *ref_body;
  
  TopologyEntity *topo_ent;
  DLIList <Body*> tmp_bodies;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_bodies.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    if ( topo_ent == NULL )
      continue;
    topo_ent->bodies(tmp_bodies);
    for ( jj = tmp_bodies.size(); jj > 0; jj-- )
    {
      ref_body = tmp_bodies.pop();
      if ( !ref_body->marked() )
      {
        ref_body->marked(1);
        ref_bodies.append(ref_body);
      }
    }
  }
  for ( ii = ref_bodies.size(); ii > 0; ii-- )
    ref_bodies.get_and_step()->marked(0);
}

void GeomMeasureTool::get_faces_from_list(DLIList <RefEntity*> &entity_list,
                                          DLIList <RefFace*> &ref_faces )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefFace *ref_face;
  
  TopologyEntity *topo_ent;
  DLIList <RefFace*> tmp_faces;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_faces.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    topo_ent->ref_faces(tmp_faces);
    for ( jj = tmp_faces.size(); jj > 0; jj-- )
    {
      ref_face = tmp_faces.pop();
      if ( !ref_face->marked() )
      {
        ref_face->marked(1);
        ref_faces.append(ref_face);
      }
    }
  }
  for ( ii = ref_faces.size(); ii > 0; ii-- )
    ref_faces.get_and_step()->marked(0);
}
void GeomMeasureTool::get_volumes_from_list(DLIList <RefEntity*> &entity_list,
                                            DLIList <RefVolume*> &ref_volumes )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefVolume *ref_volume;
  
  TopologyEntity *topo_ent;
  DLIList <RefVolume*> tmp_volumes;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_volumes.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    topo_ent->ref_volumes(tmp_volumes);
    for ( jj = tmp_volumes.size(); jj > 0; jj-- )
    {
      ref_volume = tmp_volumes.pop();
      if ( !ref_volume->marked() )
      {
        ref_volume->marked(1);
        ref_volumes.append(ref_volume);
      }
    }
  }
  for ( ii = ref_volumes.size(); ii > 0; ii-- )
    ref_volumes.get_and_step()->marked(0);
}
void GeomMeasureTool::measure_curve_length(DLIList <RefEntity*> &entity_list,
                                           double &smallest,
                                           RefEdge *&smallest_edge,
                                           double &largest,
                                           RefEdge *&largest_edge,
                                           double &average,
                                           double &sum)
{
  DLIList<RefEdge*> ref_edges;
  get_edges_from_list(entity_list, ref_edges);
  measure_curve_length(ref_edges, smallest, smallest_edge,
                       largest, largest_edge, average, sum);
  
}

void GeomMeasureTool::measure_curve_length(DLIList <RefEdge*> &ref_edges,
                                           double &smallest,
                                           RefEdge *&smallest_edge,
                                           double &largest,
                                           RefEdge *&largest_edge,
                                           double &average,
                                           double &sum)
{
  int ii;
  double curr_length, total_length = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefEdge *curr_edge;
  for ( ii = ref_edges.size(); ii > 0; ii-- )
  {
    curr_edge = ref_edges.get_and_step();
    if( curr_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;
    curr_length = curr_edge->measure();
    if ( curr_length < smallest )
    {
      smallest = curr_length;
      smallest_edge = curr_edge;
    }
    if ( curr_length > largest )
    {
      largest = curr_length;
      largest_edge = curr_edge;
    }
    total_length += curr_length;
  }
  if ( ref_edges.size() != 0 )
    average = total_length/ref_edges.size();
  sum = total_length;
  
}
void GeomMeasureTool::measure_face_curves( RefFace *ref_face,
                                           double &min_curve_length,
                                           double &max_curve_ratio,
                                           RefEdge *&min_ref_edge)
{
  DLIList<DLIList<CoEdge*>*> co_edge_loops;
  DLIList <CoEdge*> *co_edge_list;
  ref_face->co_edge_loops(co_edge_loops);
  CoEdge *curr_co_edge, *next_co_edge;
  min_curve_length = CUBIT_DBL_MAX;
  double curve_length, next_curve_length, ratio;
  max_curve_ratio = -CUBIT_DBL_MAX;
  CubitBoolean ratio_set = CUBIT_FALSE;
  int ii, jj;
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list = co_edge_loops.get_and_step();
    curve_length = co_edge_list->get()->get_ref_edge_ptr()->measure();
    for ( jj = co_edge_list->size(); jj > 0; jj-- )
    {
      curr_co_edge = co_edge_list->get_and_step();
      next_co_edge = co_edge_list->get();
      next_curve_length = next_co_edge->get_ref_edge_ptr()->measure();
        //find the smallest curve
      if ( curve_length < min_curve_length )
      {
        min_curve_length = curve_length;
        min_ref_edge = curr_co_edge->get_ref_edge_ptr();
      }
      if ( co_edge_list->size() > 1 )
      {
          //find the biggest change in curve lengths.
        if ( curve_length > next_curve_length )
        {
          if ( next_curve_length < CUBIT_RESABS )
            ratio = CUBIT_DBL_MAX;
          else
            ratio = curve_length / next_curve_length;
        }
        else
        {
          if ( curve_length < CUBIT_RESABS )
            ratio = CUBIT_DBL_MAX;
          else
            ratio = next_curve_length /curve_length;
        }
        curve_length = next_curve_length;
        if ( ratio > max_curve_ratio )
        {
          ratio_set = CUBIT_TRUE;
          max_curve_ratio = ratio;
        }
      }
    }
  }
  if ( !ratio_set )
    max_curve_ratio = -1.0;
}

        
  
  
CubitStatus GeomMeasureTool::measure_face_loops( RefFace *face,
                                                 double &min_distance_between_loops,
                                                 double &min_distance_in_one_loop,
                                                 double &min_angle,
                                                 double &max_angle,
                                                 double tolerance)
{
    //This is the most tricky of these functions.
    //To do this I'm going to facet the edges, then use an AbstractTree to
    //find the minimum distance between the loops and between a single
    //loop.
  PointLoopList boundary_point_loops;
  CubitStatus stat;
  stat = get_boundary_points(face, boundary_point_loops, tolerance);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  
  SegLoopList boundary_seg_loops;
  stat = convert_to_lines( boundary_point_loops,
                           boundary_seg_loops,
                           face);
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Now add the points to the AbstractTree.
  DLIList<AbstractTree<GeomSeg*>*> atree_list;
  AbstractTree<GeomSeg*> *curr_tree;
  SegList *seg_list;
  GeomSeg  *curr_seg, *next_seg;
  int ii,jj, kk;
  double angle;
  min_angle = CUBIT_DBL_MAX;
  max_angle = 0.0;
  double creation_time=0.0;
  CpuTimer creation_timer;
  const double angle_convert = 180.0/CUBIT_PI;
  boundary_seg_loops.reset();
  int num_segs = 0;
  for (ii = 0; ii < boundary_seg_loops.size(); ii++ )
  {
    seg_list = boundary_seg_loops.get_and_step();
    creation_timer.cpu_secs();
    curr_tree = new RTree<GeomSeg*>(GEOMETRY_RESABS);
    creation_time += creation_timer.cpu_secs();
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
        //build the r-tree.
      num_segs++;
      creation_timer.cpu_secs();
      curr_seg = seg_list->get_and_step();
      curr_tree->add(curr_seg);
      creation_time += creation_timer.cpu_secs();
        //calculate the interior angles.
      next_seg = seg_list->get();

      RefEntity* temp_ent = curr_seg->get_end()->owner();
        //only measure the angle at vertices.  This
        //should make things slightly faster.
      if ( CAST_TO(temp_ent, RefVertex) )
      {
        stat = interior_angle(curr_seg, next_seg,
                              face, angle);
        if ( stat != CUBIT_SUCCESS )
        {
          min_angle = 0.0;
          max_angle = 360.0;
        }
        else if ( angle < min_angle )
        {
          min_angle = angle;
        }
        else if ( angle > max_angle )
        {
          max_angle = angle;
        }
      }
    }
    atree_list.append(curr_tree);
  }
  min_angle *= angle_convert;
  max_angle *= angle_convert;

  if ( atree_list.size() <= 1 )
    min_distance_between_loops = -1.0;
  else
  {
      //determine the minimum distance between the loops.
    CpuTimer atree_time;
    atree_time.cpu_secs();
    PointList *curr_points;
    GeomPoint *curr_point;
    DLIList <GeomSeg*> nearest_neighbors;
    
    CubitVector curr_loc;
    double closest_dist;
    min_distance_between_loops = CUBIT_DBL_MAX;
    atree_list.reset();
    for ( ii = 0; ii < atree_list.size(); ii++ )
    {
      curr_points = boundary_point_loops.get_and_step();
      for ( jj = ii+1; jj < atree_list.size(); jj++ )
      {
        curr_tree = atree_list.next(jj);
          //now for every point in curr_points, find it's closest_point
          //in the curr_tree.
        for ( kk = 0; kk < curr_points->size(); kk++ )
        {
          curr_point = curr_points->get_and_step();
          curr_loc.set(curr_point->coordinates());
          nearest_neighbors.clean_out();
          curr_tree->k_nearest_neighbor(curr_loc,
                                        1, closest_dist, nearest_neighbors,
                                        dist_sq_point_data);
          if ( nearest_neighbors.size() == 0)
          {
              //hmm, not sure what is wrong here.
            PRINT_ERROR("Can't find closest point between loops.\n");
            return CUBIT_FAILURE;
          }
          if (closest_dist < min_distance_between_loops)
            min_distance_between_loops = closest_dist;
        }
      }
    }
    double time =  atree_time.cpu_secs();
    PRINT_DEBUG_139("Time to create atree: %f\n", creation_time);
    PRINT_DEBUG_139("Time for atree nn search: %f\n", time);
  }
  if ( min_distance_between_loops != -1.0 )
  min_distance_between_loops = sqrt(min_distance_between_loops);
    //Now find the min distance inside a loop.  This is tricky
    //in that we want to find distances that are across some area
    //and not just between adjacent curves.
  atree_list.reset();
  GeomSeg *tmp_seg;
  double dist;
  min_distance_in_one_loop = CUBIT_DBL_MAX;
  int k=5;
  CubitVector curr_loc;
  DLIList <GeomSeg*> nearest_neighbors;
  RefEntity *ent_1, *ent_2;
  TopologyEntity *top_1, *top_2;
  
  
  for ( ii = atree_list.size(); ii > 0; ii-- )
  {
    curr_tree = atree_list.get_and_step();
    seg_list = boundary_seg_loops.get_and_step();
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
      nearest_neighbors.clean_out();
      curr_seg = seg_list->get_and_step();
        //get the k closest segments.
      curr_loc.set(curr_seg->get_start()->coordinates());
      curr_tree->k_nearest_neighbor(curr_loc,
                                    k, dist, nearest_neighbors,
                                    dist_sq_point_data);
        //go through these nearest_neighbors and through out
        //the segments that are attached to curr_seg.
      for ( kk = 0; kk < nearest_neighbors.size(); kk++ )
      {
        tmp_seg = nearest_neighbors.get();
        if ( tmp_seg == curr_seg ||
             tmp_seg == curr_seg->get_prev() ||
             tmp_seg == curr_seg->get_next() )
          nearest_neighbors.remove();
      }
      if ( nearest_neighbors.size() == 0 )
        continue;
        //The min distance is the distance between the curr_loc
        //and the top nearest_neighbor.
      nearest_neighbors.reset();
      tmp_seg = nearest_neighbors.get();
      ent_1 = tmp_seg->owner();
      ent_2 = curr_seg->get_start()->owner();
      top_1 = CAST_TO(ent_1, TopologyEntity);
      top_2 = CAST_TO(ent_2, TopologyEntity);
      
      if ( top_1->is_directly_related(top_2) )
        continue;
      
        //we'll need to recompute it...
      dist = dist_sq_point_data(curr_loc, tmp_seg);
      if ( dist < min_distance_in_one_loop )
        min_distance_in_one_loop = dist;
    }
  }
  if ( min_distance_in_one_loop == CUBIT_DBL_MAX )
  {
    ent_1 = (RefEntity*)face;
    DLIList<RefEntity*> entities;
    entities.append(ent_1);
    double smallest, largest, ave, sum;
    RefEdge *small_e, *large_e;
    measure_curve_length(entities, smallest, small_e,
                         largest, large_e, ave, sum);
    min_distance_in_one_loop = smallest;
  }
  else
    min_distance_in_one_loop = sqrt(min_distance_in_one_loop);
    //clean up the memory.
  int ll,mm;
  atree_list.reset();
  boundary_seg_loops.reset();
  for ( ll = boundary_seg_loops.size(); ll > 0; ll-- )
  {
    delete atree_list.pop();
    seg_list = boundary_seg_loops.pop();
    for ( mm = seg_list->size(); mm > 0; mm-- )
      delete seg_list->pop();
    delete seg_list;
  }
		
  PointList *point_list;
  for ( ll = boundary_point_loops.size(); ll > 0; ll-- )
  {
    point_list = boundary_point_loops.pop();
    for ( mm = point_list->size(); mm > 0; mm-- )
    {
      GeomPoint *tmp_point = point_list->pop();
        //delete the CubitPointData
      delete tmp_point;
    }
    delete point_list;
  }
  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------
//Function: dist_point_data
//Computes the distance between the line segment and the vector position.
//This function is used to determine the nearest neigbhor search and
//gets passed to the AbstractTree.
//-------------------------------------------------------------------
double GeomMeasureTool::dist_sq_point_data( CubitVector &curr_point,
                                            GeomSeg *&curr_seg )
{
  double node[3];
  double point_0[3], point_1[3];
  node[0] = curr_point.x();
  node[1] = curr_point.y();
  node[2] = curr_point.z();
  CubitVector start_v = curr_seg->get_start()->coordinates();
  CubitVector end_v = curr_seg->get_end()->coordinates();
  point_0[0] = start_v.x();
  point_0[1] = start_v.y();
  point_0[2] = start_v.z();
  point_1[0] = end_v.x();
  point_1[1] = end_v.y();
  point_1[2] = end_v.z();

  IntersectionTool new_tool(GEOMETRY_RESABS);
  double t;
  double dist = new_tool.distance_point_line(node, point_0,
                                             point_1, t);
  if ( dist < 0.0 )
  {
      //This happens if the point is beyond the line segment, in which
      //case we just want the closer end point.
    if ( t < 0.0 )
    {
      dist = (curr_point -start_v).length_squared();
    }
    else
    {
      dist = (curr_point - end_v).length_squared();
    }
  }
  else
  {
      //I hate to do this but everything else is in
      //terms of distance squared.
    dist = dist*dist;
  }
  return dist;
}

CubitStatus GeomMeasureTool::convert_to_lines( PointLoopList &boundary_point_loops,
                                               SegLoopList &boundary_line_loops,
                                               RefFace *ref_face )
{
    //This does more than just converts the boundary points to line segments.  It also
    //Sets up a doubily linked list through the ImprintLineSegment data structure.
    //Calling codes need to make sure all this data is deleted.
    //This function also sets the PointType for each of the boundary points.  This
    //will depend also on the boundary_1 flag.  If we are doing boundary_1 then,
    //the flag will be true, otherwise it will be false.
  int ii, jj;
  PointList *point_list;
  SegList *segment_list;
  GeomSeg *new_line, *prev = NULL, *next = NULL;
  SegList allocated_segs;
  RefEntity *entity_start, *entity_end, *seg_owner;
  boundary_point_loops.reset();
  for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
  {
    point_list = boundary_point_loops.get_and_step();
    segment_list = new SegList;
    prev = NULL;
    next = NULL;
    new_line = NULL;
    seg_owner = NULL;
    for ( jj = point_list->size(); jj > 0; jj-- )
    {
      entity_start = point_list->get()->owner();
      entity_end = point_list->next()->owner();
      seg_owner =  CAST_TO(entity_start, RefEdge);
      if ( seg_owner == NULL )
      {
        seg_owner = CAST_TO(entity_end, RefEdge);
        if ( seg_owner == NULL )
        {
          RefVertex *ref_vert_end = CAST_TO(entity_end, RefVertex);
          RefVertex *ref_vert_start = CAST_TO(entity_start, RefVertex);
          DLIList <RefEdge*> common_edges;
          ref_vert_start->common_ref_edges(ref_vert_end,
                                           common_edges,
                                           ref_face);
          if ( common_edges.size() == 0  )
          {
            PRINT_ERROR("GeomMeasureTool::convert_to_lines having problem finding owner \n" 
                        "       of boundary.  Could be indicative of corrupt geometry.\n");
            int ll;
            for ( ll = boundary_line_loops.size(); ll > 0; ll-- )
              delete boundary_line_loops.pop();
            for ( ll = allocated_segs.size(); ll > 0; ll-- )
              delete allocated_segs.pop();
            delete segment_list;
            return CUBIT_FAILURE;
          }
          if ( common_edges.size() == 1 )
            seg_owner = common_edges.get();
          else
          {
              //Now we have to decide which of these edges is the one.
              //Lets take a mid point on this segment and check the distances.
              //First find the mid_point of the line segment.
            CubitVector start = point_list->get()->coordinates();
            CubitVector end = point_list->next()->coordinates();
            CubitVector mid = (start+end)/2.0;
            CubitVector closest_point;
            
              //Now test the edges to see which is closest.
            RefEdge *closest_edge = NULL, *curr_ref_edge;
            double min_dist = CUBIT_DBL_MAX, dist;
            int kk;
            for ( kk = common_edges.size(); kk > 0; kk-- )
            {
              curr_ref_edge = common_edges.get_and_step();
              curr_ref_edge->closest_point(mid, closest_point);
              dist = (mid - closest_point).length_squared();
              if ( dist < min_dist )
              {
                min_dist = dist;
                closest_edge = curr_ref_edge;
              }
            }
            if ( closest_edge == NULL )
            {
              PRINT_ERROR("Problems determining segment ownwership.\n"
                          "Internal problem with virtual imprinting.\n");
              int ll;
              for ( ll = boundary_line_loops.size(); ll > 0; ll-- )
                delete boundary_line_loops.pop();
              for ( ll = allocated_segs.size(); ll > 0; ll-- )
                delete allocated_segs.pop();
              delete segment_list;
              return CUBIT_FAILURE;
            }
            seg_owner = closest_edge;
          }
        }
      }
      new_line = new GeomSeg( point_list->get(), point_list->next(), seg_owner);
      allocated_segs.append(new_line);
      point_list->step();
      segment_list->append(new_line);
      if ( prev != NULL )
      {
        prev->set_next(new_line);
        new_line->set_prev(prev);
      }
      prev = new_line;
    }
    if ( prev != NULL )
    {
      segment_list->reset();
      assert(prev == segment_list->prev());
      segment_list->prev()->set_next(segment_list->get());
      segment_list->get()->set_prev(segment_list->prev());
      assert(segment_list->get()->get_next() != NULL );
    }
    if ( segment_list->size() )
    {
      boundary_line_loops.append(segment_list);
    }
    else
      delete segment_list;
  }
  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::get_boundary_points( RefFace *ref_face,
                                                  PointLoopList &boundary_point_loops,
                                                  double seg_length_tol)
{
  DLIList<DLIList<CoEdge*>*> co_edge_loops;
  ref_face->co_edge_loops(co_edge_loops);
  int ii, jj;
  DLIList <CoEdge*> *co_edge_list_ptr;
  PointList *new_point_loop_ptr, tmp_point_list;
  RefEdge *ref_edge_ptr;
  CoEdge *co_edge_ptr;
  CubitStatus stat;
  CubitSense sense;
  
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list_ptr = co_edge_loops.get_and_step();
    new_point_loop_ptr = new PointList;
    for ( jj = co_edge_list_ptr->size(); jj > 0; jj-- )
    {
      co_edge_ptr = co_edge_list_ptr->get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      tmp_point_list.clean_out();
      stat = get_curve_facets( ref_edge_ptr, tmp_point_list,
                               seg_length_tol );
      if ( stat != CUBIT_SUCCESS )
        return CUBIT_FAILURE;
      if ( tmp_point_list.size() == 0 )
        continue;
      tmp_point_list.reset();
        //the points are in order from start vertex to end vertex.
        //append them now according to the loop.
      sense = co_edge_ptr->get_sense();
      if ( CUBIT_FORWARD != sense )
        tmp_point_list.reverse();
        //Now take off the last point as it is a duplicate with the
        //other list...
      tmp_point_list.reset();
      delete tmp_point_list.pop();
      (*new_point_loop_ptr) += tmp_point_list;
    }
    boundary_point_loops.append(new_point_loop_ptr);
  }
    //clean up the list memory.
  for(ii = co_edge_loops.size(); ii>0; ii-- )
    delete co_edge_loops.pop();
  co_edge_loops.clean_out();
      
  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::get_curve_facets( RefEdge* curve,
                                               PointList &segments,
                                               double seg_length_tol ) 
{
//  const double COS_ANGLE_TOL =  0.965925826289068312213715; // cos(15)
//  const double COS_ANGLE_TOL =  0.984807753012208020315654; // cos(10)
  //const double COS_ANGLE_TOL =  0.996194698091745545198705; // cos(5)
  GMem curve_graphics;
    //make this tol bigger than seg_length_tol, to
    //make sure the segments are larger than the tolerance.
  const double dist_tol = seg_length_tol;// + .05*seg_length_tol;
  //const double dist_tol_sqr = dist_tol*dist_tol;
  if ( curve->get_curve_ptr() == NULL )
  {
    PRINT_ERROR("Curve %d had no geometry, this is a bug!\n", curve->id());
    return CUBIT_FAILURE;
  }
  else if ( curve->get_curve_ptr()->get_geometry_query_engine() == NULL )
  {
    PRINT_ERROR("Curve %d had no geometry engine, this is a bug!\n",
                curve->id());
    return CUBIT_FAILURE;
  }
  curve->get_graphics( curve_graphics );
  
  GPoint* gp = curve_graphics.point_list();
  if ( gp == NULL )
  {
    if ( curve->measure() < GEOMETRY_RESABS ) 
    {
      if( curve->geometry_type() != POINT_CURVE_TYPE )
      {
        PRINT_INFO("Curve %d has a zero length, and won't be used for loop measurements.\n",
                    curve->id());
      }
    }
    else
      PRINT_ERROR("Curve %d had not faceting.  Ignoring this curve.\n",
                  curve->id());
    return CUBIT_SUCCESS;
  }
  GeomPoint *last = new GeomPoint( gp[0].x, gp[0].y, gp[0].z, curve );
  CubitVector lastv = last->coordinates();
  int num_points = curve_graphics.pointListCount;
  segments.append( last );
  int ii;
  CubitBoolean remove_second_to_end = CUBIT_FALSE;
  for ( ii = 1; ii < num_points; ii++ )
  {
    CubitVector pos(  gp[ii].x, gp[ii].y, gp[ii].z );
    CubitVector step1 = (pos - lastv);
    double len1 = step1.length();
    if( len1 < dist_tol && ii != num_points - 1) 
      continue;
    else if ( len1 < dist_tol && ii == num_points-1 )
    {
      remove_second_to_end = CUBIT_TRUE;
    }
    last = new GeomPoint(pos,curve);
    segments.append( last );
    lastv = last->coordinates();
  }
    // Now check if the segment list is reversed wrt the curve direction.
  segments.reset();
  if ( remove_second_to_end )
  {
    if ( segments.size() == 2 )
    {
        //just leave it in, let it be small...
    }
    else
    {
        //Remove the second to last one.  To do
        //this efficiently (don't do remove), pop
        //the last one, then save that and
        //re-add it after poping the second one.
      GeomPoint *temp = segments.pop();
      segments.pop();
      segments.append(temp);
    }
  }
  segments.reset();
  if( curve->start_vertex() != curve->end_vertex() )
  {
    CubitVector start_vec, end_vec;
    start_vec = curve->start_vertex()->coordinates();
    end_vec = curve->end_vertex()->coordinates();
    CubitVector start_seg = segments.get()->coordinates();
    double dist_1 = (start_seg - start_vec).length_squared();
    double dist_2 = (start_seg - end_vec).length_squared();
    if ( dist_1 > dist_2 )
      segments.reverse();
  }
  else
  {
    double u1, u2;
    u1 = curve->u_from_position( (segments.next(1)->coordinates()) );
    u2 = curve->u_from_position( (segments.next(2)->coordinates()) );    
    if( (u2 < u1) && (curve->start_param() <= curve->end_param()) )
      segments.reverse();
  }
    //clean up the periodic curve case (last seg may be too small.)
  if ( curve->start_vertex() == curve->end_vertex() )
  {
    segments.reset();
    CubitVector start_v = segments.get()->coordinates();
    CubitVector last_v = segments.prev()->coordinates();
    double dist = (start_v - last_v).length();
    if ( dist < dist_tol )
    {
        //remove the last one.
      segments.pop();
    }
  }
  if ( segments.size() < 2 )
  {
    PRINT_ERROR("Problem faceting boundary.\n");
    return CUBIT_FAILURE;
  }
  segments.get()->owner(curve->start_vertex());
  segments.prev()->owner(curve->end_vertex());
  return CUBIT_SUCCESS;
}
CubitStatus GeomMeasureTool::interior_angle(GeomSeg *first_seg,
                                            GeomSeg *next_seg,
                                            RefFace *face,
                                            double &angle )
{
  if ( first_seg->get_end() != next_seg->get_start() )
  {
    PRINT_ERROR("GeomMeasureTool::interior_angle expects segments\n"
                "to be connected (Users: this is a bug, please submit.\n");
    return CUBIT_FAILURE;
  }
  CubitVector point = first_seg->get_end()->coordinates();
  CubitVector to_prev = first_seg->get_start()->coordinates();
  to_prev -= point;
  CubitVector to_next = next_seg->get_end()->coordinates();
  to_next -= point;
  CubitVector surf_norm = face->normal_at(point);
  if ( surf_norm.length_squared() < CUBIT_RESABS )
  {
      //Try getting it at one of the other nodes...
    surf_norm = face->normal_at(first_seg->get_start()->coordinates() );
    if (surf_norm.length_squared() < CUBIT_RESABS )
    {
        //Try getting it at one of the other nodes...
      surf_norm = face->normal_at(next_seg->get_end()->coordinates() );
      if (surf_norm.length_squared() < CUBIT_RESABS )
      {
        PRINT_ERROR("Trying to get normal at surf %d.\n"
                    "       Normal length being returned equals zero.\n",
                    face->id() );
        angle = CUBIT_DBL_MAX;
        return CUBIT_FAILURE;
      }
    }
  }
  angle = surf_norm.vector_angle ( to_next, to_prev );
  return CUBIT_SUCCESS;
}

void GeomMeasureTool::measure_face_area (DLIList <RefEntity*> &entity_list,
                                         double &smallest,
                                         RefFace *&smallest_face,
                                         double &largest,
                                         RefFace *&largest_face,
                                         double &average,
                                         double &sum)
{
  DLIList<RefFace*> ref_faces;
  get_faces_from_list(entity_list, ref_faces);
  measure_face_area(ref_faces, smallest, smallest_face,
                    largest, largest_face, average, sum);
  
}

void GeomMeasureTool::measure_face_area (DLIList <RefFace*> &ref_faces,
                                         double &smallest,
                                         RefFace *&smallest_face,
                                         double &largest,
                                         RefFace *&largest_face,
                                         double &average,
                                         double &sum)
{
  int ii;
  double curr_area, total_area = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefFace *curr_face;
  for( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_area = measure_area(curr_face);
    if ( curr_area < smallest )
    {
      smallest = curr_area;
      smallest_face = curr_face;
    }
    if( curr_area > largest )
    {
      largest = curr_area;
      largest_face = curr_face;
    }
    total_area += curr_area;
  }
  if (ref_faces.size() != 0 )
    average = total_area/ref_faces.size();
  sum = total_area;
  
}

void GeomMeasureTool::measure_volume_volume (DLIList <RefEntity*> &entity_list,
                                             double &smallest,
                                             RefVolume *&smallest_volume,
                                             double &largest,
                                             RefVolume *&largest_volume,
                                             double &average,
                                             double &sum)
{
  DLIList<RefVolume*> ref_volumes;
  get_volumes_from_list(entity_list, ref_volumes);
  measure_volume_volume(ref_volumes, smallest, smallest_volume,
                        largest, largest_volume, average, sum);
  
}

void GeomMeasureTool::measure_volume_volume (DLIList <RefVolume*> &ref_volumes,
                                             double &smallest,
                                             RefVolume *&smallest_volume,
                                             double &largest,
                                             RefVolume *&largest_volume,
                                             double &average,
                                             double &sum)
{
  int ii;
  double curr_volume, total_volume = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefVolume *curr_solid;
  for ( ii = ref_volumes.size(); ii > 0; ii-- )
  {
    curr_solid = ref_volumes.get_and_step();
    curr_volume = curr_solid->measure();

    if ( curr_volume < smallest )
    {
      smallest = curr_volume;
      smallest_volume = curr_solid;
    }
    if ( curr_volume > largest )
    {
      largest = curr_volume;
      largest_volume = curr_solid;
    }
    total_volume += curr_volume;
  }
  if (ref_volumes.size() != 0 )
    average = total_volume/ref_volumes.size();
  sum = total_volume;
  
}

void GeomMeasureTool::measure_face_area_and_hydraulic_radius(RefFace *curr_face,
                                                             double &face_area,
                                                             double &face_hydraulic_radius)
{
  int ii;
  double total_edge_length = 0.0, curr_edge_length;
  RefEdge *curr_edge;
  DLIList <CoEdge*> my_co_edges;

  face_area = measure_area(curr_face);
  curr_face->co_edges(my_co_edges);
  for( ii = my_co_edges.size(); ii > 0; ii-- )
  {
    curr_edge = my_co_edges.get_and_step()->get_ref_edge_ptr();
    curr_edge_length = curr_edge->measure();
    total_edge_length += curr_edge_length;
  }
  if (total_edge_length != 0)
    face_hydraulic_radius = 4.0*(face_area / total_edge_length);
}

void GeomMeasureTool::measure_volume_volume_and_hydraulic_radius(RefVolume *curr_volume,
                                                                 double &volume_volume,
                                                                 double &volume_hydraulic_radius)
{
  int ii;
  double curr_face_area, total_surface_area = 0.0;
  DLIList <RefFace*> ref_face_list;
  RefFace *curr_face;

  volume_volume = curr_volume->measure();
  curr_volume->ref_faces(ref_face_list);
  
  for ( ii = ref_face_list.size(); ii > 0; ii-- )
  {
    curr_face = ref_face_list.get_and_step();
    curr_face_area = measure_area(curr_face);
    total_surface_area += curr_face_area;
  }
  if (total_surface_area != 0)
    volume_hydraulic_radius = 6.0*(volume_volume / total_surface_area);
 
}

void GeomMeasureTool::angles_between_volume_surfaces(RefVolume *curr_volume,
                                                     double &min_angle,
                                                     double &max_angle,
                                                     RefFace *&face_min_1,
                                                     RefFace *&face_min_2,
                                                     RefFace *&face_max_1,
                                                     RefFace *&face_max_2)
{
  int ii, jj, shared_angles = 0;
  double common_angle, sum_of_angles = 0.0;
  double smallest = CUBIT_DBL_MAX, largest = 0.0;
  DLIList <RefFace*> face_list_1, face_list_2;
  RefFace *curr_face_1, *curr_face_2;
#ifdef BOYD17
  DLIList <RefEdge*> my_ref_edges;
#endif
  RefEdge *common_edge;

  curr_volume->ref_faces(face_list_1);
  face_list_2 = face_list_1;

  for (ii = face_list_1.size(); ii > 0; ii--)
  {
    curr_face_1 = face_list_1.get_and_step();
    for(jj = face_list_2.size(); jj > 0; jj--)
    {
      curr_face_2 = face_list_2.get_and_step();

      if(curr_face_1 == curr_face_2)
        continue;
      else
      {
        common_edge = curr_face_1->common_ref_edge(curr_face_2);
        if( common_edge == NULL )
          continue;
        else
        {
          shared_angles++;
          
          if ( common_edge->get_curve_ptr()->geometry_type() != STRAIGHT_CURVE_TYPE )
          {
            int kk;
            double frac_pos = 0.00;
            for ( kk = 0; kk < 4; kk++ )
            {
              common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1, curr_face_2,
                                                                           common_edge, curr_volume,
                                                                           frac_pos);
              frac_pos += .25;
              if( common_angle < smallest )
              {
                min_angle = common_angle * (180 / CUBIT_PI);
                smallest = common_angle;
                face_min_1 = curr_face_1;
                face_min_2 = curr_face_2;
              }
              if( common_angle > largest )
              {
                max_angle = common_angle * (180 / CUBIT_PI);
                largest = common_angle;
                face_max_1 = curr_face_1;
                face_max_2 = curr_face_2;
              }
            }
          }
          else
          {
            
            common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1, curr_face_2,
                                                                         common_edge, curr_volume );
            if( common_angle < smallest )
            {
              min_angle = common_angle * (180 / CUBIT_PI);
              smallest = common_angle;
              face_min_1 = curr_face_1;
              face_min_2 = curr_face_2;
            }
            if( common_angle > largest )
            {
              max_angle = common_angle * (180 / CUBIT_PI);
              largest = common_angle;
              face_max_1 = curr_face_1;
              face_max_2 = curr_face_2;
            }
          }
          
          sum_of_angles += common_angle;
        }
      }
    }
  }
}

void GeomMeasureTool::merged_unmerged_surface_ratio(DLIList <RefVolume*> &ref_volumes,
                                                    int &merged, int &unmerged,
                                                    double &ratio)
  
{
  int ii, jj;
  RefVolume *curr_volume_1;
  DLIList <RefFace*> surface_list_1, surface_list_2;
  RefFace *curr_face_1;
  int merged_surfaces = 0, unmerged_surfaces = 0;
  
  for( ii = ref_volumes.size(); ii > 0; ii-- )
  {
    curr_volume_1 = ref_volumes.get_and_step();
    surface_list_1.clean_out();
    
    curr_volume_1->ref_faces(surface_list_1);
    for( jj = surface_list_1.size(); jj > 0; jj-- )
    {
      curr_face_1 = surface_list_1.get_and_step();
      if ( curr_face_1->marked() )
        continue;
      else
      {
        curr_face_1->marked(1);
        surface_list_2.append(curr_face_1);
      }
      if ( curr_face_1->num_ref_volumes() == 2 )
        merged_surfaces++;
      else
        unmerged_surfaces++;
    }
  }
  for ( ii = surface_list_2.size(); ii > 0; ii-- )
  {
    surface_list_2.get_and_step()->marked(0);
  }
  merged = merged_surfaces;
  unmerged = unmerged_surfaces;
  if ( unmerged_surfaces == 0 )
  {
    ratio = 100.0;
  }
  else
    ratio = (double) ((double)merged_surfaces / (double)unmerged_surfaces);
}

void GeomMeasureTool::report_intersected_volumes(DLIList <RefVolume*> &ref_vols,
                                                 DLIList <RefVolume*> &intersection_list)
{
  DLIList <RefVolume*> results;
  RefVolume *curr_vol, *curr_vol_2;
  int i, j;
  ProgressTool *progress_ptr = NULL;
  char title[29];
  int total_volumes = ref_vols.size();
  if (total_volumes > 5)
  {
    strcpy(title, "Overlapping Volumes Progress");
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, title, NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;

  for( i = 0; i < ref_vols.size(); i++ )
  {
    curr_vol = ref_vols.next(i);

    //is this body a multi-volume-body?
    Body *single_volume_body = curr_vol->body();
    DLIList<RefVolume*> body_vols;
    single_volume_body->ref_volumes( body_vols );
    if( body_vols.size() > 1 )
      single_volume_body = NULL;
      
    for( j = (i + 1); j < ref_vols.size(); j++ )
    {
      RefVolume *curr_vol2 = ref_vols.next(j);
      if ( CubitMessage::instance()->Interrupt() )
      {
          //interrpt.  We need to exit.
        if ( progress_ptr != NULL )
        {
          progress_ptr->end();
        }
          //just leave what has been calculated...
        return;
      }

      Body *single_volume_body2 = curr_vol2->body();
      DLIList<RefVolume*> body_vols2;
      single_volume_body2->ref_volumes( body_vols2 );
      if( body_vols2.size() > 1 )
        single_volume_body2 = NULL;

      //update the progress..
      if ( progress_ptr != NULL )
      {
        curr_percent = ((double)(i+1))/((double)(total_volumes));
        progress_ptr->percent(curr_percent);
      }
    
      //if both are single-volume-bodies
      if( single_volume_body && single_volume_body2 )
      {
        if( GeometryQueryTool::instance()->bodies_overlap( single_volume_body,
                                                       single_volume_body2 ) )
        {
          intersection_list.append( curr_vol );
          intersection_list.append( curr_vol2 );
        }
      }
      else if( GeometryQueryTool::instance()->volumes_overlap( curr_vol, curr_vol2 ) )
      {
       intersection_list.append( curr_vol );
       intersection_list.append( curr_vol2 );
      }
    }
  }
  
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
}

void GeomMeasureTool::report_intersected_bodies(DLIList <Body*> &ref_bodies,
                                                DLIList <Body*> &intersection_list)
{
  Body *curr_body, *curr_body_2;
  int ii, jj;
  ProgressTool *progress_ptr = NULL;
  char title[29];
  int total_bodies = ref_bodies.size();
  if (total_bodies > 5)
  {
    strcpy(title, "Overlapping Volumes Progress");
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, title, NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;
 
  
  for( ii = 0; ii < ref_bodies.size(); ii++ )
  {
    curr_body = ref_bodies.next(ii);

    for( jj = (ii + 1); jj < ref_bodies.size(); jj++ )
    {
      curr_body_2 = ref_bodies.next(jj);
      if ( CubitMessage::instance()->Interrupt() )
      {
          //interrpt.  We need to exit.
        if ( progress_ptr != NULL )
        {
          progress_ptr->end();
        }
          //just leave what has been calculated...
        return;
      }

      if (GeometryQueryTool::instance()->bodies_overlap(curr_body,
                                                        curr_body_2) )
      {
        intersection_list.append(curr_body);
        intersection_list.append(curr_body_2);
      }
    }
       //update the progress..
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(ii+1))/((double)(total_bodies));
      progress_ptr->percent(curr_percent);
    }
  }
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
}

void GeomMeasureTool::face_list_from_volume_list( DLIList <RefVolume*> &input_vols,
                                                  DLIList <RefFace*> &all_faces )
{
  DLIList <RefFace*> curr_faces;
  RefVolume *curr_vol;
  RefFace *curr_face;
  int ii, jj;

    //get the all RefFaces from input_vols.
    //put them into all_faces list.
  
  for ( ii = 0; ii < input_vols.size(); ii++ )
  {
    curr_vol = input_vols.get_and_step();
    curr_faces.clean_out();
    curr_vol->ref_faces(curr_faces);
    
    for( jj = 0; jj < curr_faces.size(); jj++ )
    {
      curr_face = curr_faces.get_and_step();
      if( curr_face->marked() )
        continue;
      else
      {
        curr_face->marked(1);
        all_faces.append(curr_face);
      }
    }
  }
  for ( ii = all_faces.size();ii >0; ii-- )
    all_faces.get_and_step()->marked(0);
}

RefFace* GeomMeasureTool::valid_start( DLIList <RefFace*> &all_faces )
{
  int ii;
  RefFace *start_face;

  for( ii = 0; ii < all_faces.size(); ii++ )
  {
    start_face = all_faces.get_and_step();

    if( !start_face->marked() && start_face->num_ref_volumes()== 1 )
      return start_face;
  }

  return NULL;
  
}

void GeomMeasureTool::find_shells( DLIList <RefVolume*> &input_vols,
                                   RefGroup *&owner_groups,
                                   int &number_of_shells)
{
  DLIList <DLIList<RefEntity*>*> list_of_lists;
  DLIList <RefEntity*> *new_list;
  DLIList <RefFace*> all_faces, stack, related_faces;
  DLIList <RefEdge*> curr_edges;
  RefFace *start_face, *curr_face, *related_face;
  RefEdge *curr_edge;
  int ii, jj;
  RefEntity *tmp_ent;
  number_of_shells = 0;
  DLIList <RefFace*> marked_faces;
  
  face_list_from_volume_list( input_vols, all_faces );
  
    // all_faces all marked - need unmarked
  for (ii = 0; ii < all_faces.size(); ii++ )
  {
    curr_face = all_faces.get_and_step();
    curr_face->marked(0);
  }
  while ( (start_face = valid_start(all_faces) ) != NULL )
  {
      //create a list.
    stack.clean_out();
    new_list = new DLIList<RefEntity*>;
    stack.append(start_face);
    start_face->marked(1);
    marked_faces.append(start_face);

    while (stack.size() > 0 )
    {
      curr_edges.clean_out();

      curr_face = stack.pop(); // get last face added
      tmp_ent = CAST_TO(curr_face, RefEntity);
      new_list->append(tmp_ent); // append that face to the newest list
      curr_face->ref_edges(curr_edges); // gets list of face's edges

      for( ii = 0; ii < curr_edges.size(); ii++ )  // loops edges
      {
        related_faces.clean_out();

        curr_edge = curr_edges.get_and_step(); // chooses next edge
        curr_edge->ref_faces(related_faces); // finds edge's related faces
        for( jj = 0; jj < related_faces.size(); jj++ ) // loops faces
        {
          related_face = related_faces.get_and_step(); // gets face
          if(!related_face->marked() && related_face->num_ref_volumes() == 1)
          {
            stack.append(related_face);
            related_face->marked(1);
            marked_faces.append(related_face);
          }
        }
      }
    }
    list_of_lists.append(new_list);
    number_of_shells++;
  }

  for (ii = 0; ii < marked_faces.size(); ii++ )
  {
    curr_face = marked_faces.get_and_step();
    curr_face->marked(0);
  }
  
  RefGroup *curr_group;
  owner_groups = RefEntityFactory::instance()->construct_RefGroup("volume_shells");
  for( ii = list_of_lists.size(); ii > 0; ii-- )
  {
    new_list = list_of_lists.pop();
    curr_group = RefEntityFactory::instance()->construct_RefGroup(*new_list);
    tmp_ent = CAST_TO(curr_group, RefEntity);
    owner_groups->add_ref_entity(tmp_ent);
    delete new_list;
  } 
}

void GeomMeasureTool::ratio_of_shells_to_volumes(int number_of_shells,
                                                 DLIList <RefVolume*> &ref_volumes,
                                                 int &number_of_volumes,
                                                 double &ratio)
{
  int ii, number = 0;

  for( ii = ref_volumes.size(); ii > 0; ii-- )
    number++;
  number_of_volumes = number;

  ratio = (double) ((double)number_of_shells / (double)number_of_volumes);
  
}
void GeomMeasureTool::print_surface_measure_summary( DLIList <RefFace*> &ref_faces )
                                                  
{
  
  double min_dist_loops=CUBIT_DBL_MAX, min_dist_loop=CUBIT_DBL_MAX;
  double curr_min_dist_loops=CUBIT_DBL_MAX, curr_min_dist_loop=CUBIT_DBL_MAX;
  double min_angle=CUBIT_DBL_MAX, max_angle=-CUBIT_DBL_MAX;
  double curr_min_angle=CUBIT_DBL_MAX, curr_max_angle=-CUBIT_DBL_MAX;
  double curr_face_area, curr_face_hydraulic_radius;
  double min_face_area=CUBIT_DBL_MAX, min_face_hydraulic_radius=CUBIT_DBL_MAX;
  double max_face_area=-CUBIT_DBL_MAX;
  int counter = 0;
  double total_area= 0.0;
  double min_curve_length = CUBIT_DBL_MAX;
  double max_curve_ratio = -CUBIT_DBL_MAX;
  double curve_length, curve_ratio;
  int ii;
       
  RefFace *curr_face;
  RefFace *face_min_dist_loops=NULL, *face_min_dist_loop=NULL;
  RefFace *face_min_angle=NULL, *face_max_angle=NULL;
  RefFace *face_min_area=NULL, *face_min_hydraulic_radius=NULL;
  RefFace *face_max_area=NULL;
  RefFace *face_min_curve_length=NULL, *face_max_curve_ratio=NULL;
  RefEdge *smallest_edge=NULL, *curr_small_edge=NULL;
       
       
  for( ii = 0; ii < ref_faces.size(); ii++ )
  {
    curr_face = ref_faces.get_and_step();
    if ( curr_face->num_ref_edges() > 0 )
    {
      GeomMeasureTool::measure_face_loops(curr_face,
                                          curr_min_dist_loops,
                                          curr_min_dist_loop,
                                          curr_min_angle,
                                          curr_max_angle,
                                          GEOMETRY_RESABS*500);
      GeomMeasureTool::measure_face_area_and_hydraulic_radius(curr_face,
                                                              curr_face_area,
                                                              curr_face_hydraulic_radius);
      GeomMeasureTool::measure_face_curves(curr_face, curve_length,
                                           curve_ratio, curr_small_edge);
    }
    else
      curr_face_area = curr_face->measure();
         
    total_area += curr_face_area;
    counter++;
    if ( curr_face->num_ref_edges() > 0 )
    {
      if ( curve_length < min_curve_length )
      {
        min_curve_length = curve_length;
        smallest_edge = curr_small_edge;
        face_min_curve_length = curr_face;
      }
      if ( curve_ratio != -1.0 && curve_ratio > max_curve_ratio )
      {
        max_curve_ratio = curve_ratio;
        face_max_curve_ratio = curr_face;
      }
      if ( curr_min_dist_loops != -1.0 && curr_min_dist_loops < min_dist_loops )
      {
        face_min_dist_loops = curr_face;
        min_dist_loops = curr_min_dist_loops;
      }
      if ( curr_min_dist_loop < min_dist_loop )
      {
        face_min_dist_loop = curr_face;
        min_dist_loop = curr_min_dist_loop;
      }
      if ( curr_min_angle < min_angle )
      {
        face_min_angle = curr_face;
        min_angle = curr_min_angle;
      }
      if ( curr_max_angle > max_angle )
      {
        face_max_angle = curr_face;
        max_angle = curr_max_angle;
      }

      if ( curr_face_hydraulic_radius < min_face_hydraulic_radius )
      {
        min_face_hydraulic_radius = curr_face_hydraulic_radius;
        face_min_hydraulic_radius = curr_face;
      }

    }
    if ( curr_face_area < min_face_area )
    {
      min_face_area = curr_face_area;
      face_min_area = curr_face;
    } 
    if ( curr_face_area > max_face_area )
    {
      max_face_area = curr_face_area;
      face_max_area = curr_face;
    }
  }
  if ( counter > 0 )
  {
    PRINT_INFO("\n******Summary of Surface Measure Info*******\n\n");
    PRINT_INFO("Measured %d surfaces: Average Area: %f \n\n", counter, total_area/counter );
    PRINT_INFO("Minimum Surface Area is: %f on surface %d\n\n",
               min_face_area, face_min_area->id());
    PRINT_INFO("Maximum Surface Area is: %f on surface %d\n\n",
               max_face_area, face_max_area->id());
    if ( face_min_hydraulic_radius )
      PRINT_INFO("Minimum Hydraulic Radius (Area/Perimeter) is: %f on surface %d\n\n",
                 min_face_hydraulic_radius, face_min_hydraulic_radius->id());
    if ( smallest_edge )
      PRINT_INFO("Minimum Curve Length is: %f for curve %d on surface %d\n\n",
                 min_curve_length, smallest_edge->id(), face_min_curve_length->id());
    if ( face_max_curve_ratio )
      PRINT_INFO("Maximum Ratio of Adjacent Curve Lengths is: %f on surface %d\n\n",
                 max_curve_ratio, face_max_curve_ratio->id());
    if ( face_min_dist_loops && face_min_dist_loops->num_loops() > 1 )
    {
      PRINT_INFO("Minimum Distance Between Loops is: %f on surface %d\n\n",
                 min_dist_loops, face_min_dist_loops->id());
    }
    if ( face_min_dist_loop )
      PRINT_INFO("Minimum Distance Between a Single Loop is: %f on surface %d\n\n",
                 min_dist_loop, face_min_dist_loop->id());
    if ( face_min_angle )
      PRINT_INFO("Minimum Interior Angle Between Curves is: %f Degrees on surface %d\n\n",
                 min_angle, face_min_angle->id());
    if ( face_max_angle )
      PRINT_INFO("Maximum Interior Angle Between Curves is: %f Degrees on surface %d\n\n",
                 max_angle, face_max_angle->id());
	center(ref_faces);
    PRINT_INFO("************End of Summary***********\n\n");
  }

}
void GeomMeasureTool::find_adjacent_face_ratios(RefVolume *curr_volume, double &max_face_ratio,
                                                RefFace *&big_face,
                                                RefFace *&small_face)
{
    //Loop over all the faces in the model.  Find the maximum face area ratio.
    //(larger face area / smaller face area)
  
  DLIList <Shell*> shells;
  curr_volume->shells(shells);
  DLIList <RefFace*> ref_face_stack, adj_faces;
  Shell* curr_shell;
  RefFace *curr_ref_face, *adj_face;
    //AreaHashTuple is defined in the .hpp file.
#ifdef BOYD17
  DLIList <AreaHashTuple*> *hash_table, hash_index;
#endif
  DLIList <AreaHashTuple*> *hash_table;
  max_face_ratio = -CUBIT_DBL_MAX;
  
  
  int ii, jj, table_size;
  int hash_val;
  double tmp_area;
  double tmp_ratio;
  
  AreaHashTuple* curr_tuple, *adj_tuple;
  
  for ( ii = shells.size(); ii > 0; ii-- )
  {
    curr_shell = shells.get_and_step();
    ref_face_stack.clean_out();
    curr_shell->ref_faces(ref_face_stack);
    table_size = ref_face_stack.size();
      //build a hash table with the most simple form
      //of hashing based on the ref face id.
      //To get this more reliable, make table size a prime
      //number...
      //Do this so that we don't have to calculate the
      //areas more than once...
    hash_table = new DLIList<AreaHashTuple*> [table_size];
    for ( jj = 0; jj < ref_face_stack.size(); jj++ )
    {
      curr_ref_face = ref_face_stack.get_and_step();
      tmp_area = curr_ref_face->measure();
      curr_tuple = new AreaHashTuple();
      curr_tuple->myFace = curr_ref_face;
      curr_tuple->myArea = tmp_area;
      hash_val = curr_ref_face->id() % table_size;
      hash_table[hash_val].append(curr_tuple);
    }
    
    while ( ref_face_stack.size() )
    {
      curr_ref_face = ref_face_stack.pop();
        //Now get the adjacent faces.
      get_adjacent_faces(curr_ref_face, adj_faces);
      find_index(hash_table, table_size, curr_ref_face, curr_tuple);
      assert(curr_tuple != NULL);
      if ( curr_tuple == NULL )
        return;
      CubitBoolean curr_is_big = CUBIT_TRUE;
      for ( jj = 0; jj < adj_faces.size(); jj++ )
      {
        adj_face = adj_faces.get_and_step();
        find_index(hash_table, table_size, adj_face, adj_tuple);
        if ( adj_tuple == NULL )
          continue;
        curr_is_big = CUBIT_TRUE;
        if ( curr_tuple->myArea < adj_tuple->myArea )
        {
          tmp_ratio = adj_tuple->myArea / curr_tuple->myArea;
          curr_is_big = CUBIT_FALSE;
        }
        else
          tmp_ratio = curr_tuple->myArea / adj_tuple->myArea;
        if ( tmp_ratio > max_face_ratio )
        {
          max_face_ratio = tmp_ratio;
          if ( curr_is_big )
          {
            big_face = curr_ref_face;
            small_face = adj_face;
          }
          else
          {
            big_face = adj_face;
            small_face = curr_ref_face;
          }
        }
      }
    }
    for ( jj = 0; jj < table_size; jj++ )
    {
        //delete the area tuples.
      while ( hash_table[jj].size() > 0 )
        delete hash_table[jj].pop();
    }
      //delete the hash_table.
    delete [] hash_table;
  }
}
void GeomMeasureTool::find_index(DLIList<AreaHashTuple*> *hash_table,
                                 int table_size,
                                 RefFace *ref_face,
                                 AreaHashTuple *&curr_tuple )
{
  int ii;
  int hash_val = ref_face->id() % table_size;
  if ( hash_table[hash_val].size() == 0 )
  {
      //this ref face isn't store here.
    curr_tuple = NULL;
    return;
  }
  int curr_size = hash_table[hash_val].size();
    //resolve collisions by separate chaining method...
  for ( ii = 0;ii < curr_size; ii++ )
  {
    if ( hash_table[hash_val].next(ii)->myFace == ref_face )
    {
      curr_tuple = hash_table[hash_val].next(ii);
      return;
    }
  }
  curr_tuple = NULL;
  return;
}
  
void GeomMeasureTool::get_adjacent_faces(RefFace* curr_face,
                                         DLIList<RefFace*> &adjacent_faces)
{
  DLIList <RefEdge*> ref_edges;
  DLIList <RefFace*> tmp_ref_faces;
  curr_face->ref_edges(ref_edges);
  RefEdge *ref_edge;
  RefFace *tmp_ref_face;
  int ii, jj;
  for ( ii = 0; ii < ref_edges.size(); ii++ )
  {
    ref_edge = ref_edges.get_and_step();
    tmp_ref_faces.clean_out();
    ref_edge->ref_faces(tmp_ref_faces);
    for ( jj = 0; jj < tmp_ref_faces.size(); jj++ )
    {
      tmp_ref_face = tmp_ref_faces.get_and_step();
      if ( !tmp_ref_face->marked() && tmp_ref_face != curr_face )
      {
        tmp_ref_face->marked(1);
        adjacent_faces.append(tmp_ref_face);
      }
      else
        continue;
    }
  }
  for ( ii = adjacent_faces.size(); ii > 0; ii-- )
    adjacent_faces.get_and_step()->marked(0);
}
void GeomMeasureTool::print_volume_measure_summary(DLIList <RefVolume*> &ref_volumes)
{
  int ii;
  RefVolume *curr_volume;
  RefFace *curr_face;
  DLIList <RefFace*> ref_faces, tmp_faces;
    //First get a list of all the ref_faces.
    //Get them uniquely.
  for ( ii = 0; ii< ref_volumes.size(); ii++ )
  {
    curr_volume = ref_volumes.get_and_step();
    tmp_faces.clean_out();
    curr_volume->ref_faces(tmp_faces);
    int jj;
    for ( jj = 0; jj < tmp_faces.size(); jj++ )
    {
      curr_face = tmp_faces.get_and_step();
      if ( curr_face->marked() )
        continue;
      else
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
    //reset the mark.
  for ( ii = 0; ii< ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);

  double volume_volume, volume_hydraulic_radius;
  double min_angle=CUBIT_DBL_MAX, max_angle=-CUBIT_DBL_MAX;
  double total_volume = 0.0, min_volume_hydraulic_radius = CUBIT_DBL_MAX;
  double min_volume=CUBIT_DBL_MAX; //, max_volume=-CUBIT_DBL_MAX;
  double curr_min_angle = 361.0, curr_max_angle =361.;
  double max_face_ratio = -CUBIT_DBL_MAX, curr_face_ratio=CUBIT_DBL_MAX;
         

  RefVolume *volume_min_volume = NULL, *volume_min_hydraulic_radius=NULL;
  RefVolume *volume_min_angle=NULL, *volume_max_angle=NULL;
  RefVolume *volume_face_ratio= NULL;
  RefFace *max_small_face=NULL, *max_big_face=NULL;
  RefFace *small_face=NULL, *big_face=NULL;
  RefFace *face_min_1=NULL, *face_min_2=NULL;
  RefFace *face_max_1=NULL, *face_max_2=NULL;
  RefFace *face_min_angle_1=NULL, *face_min_angle_2=NULL;
  RefFace *face_max_angle_1=NULL, *face_max_angle_2=NULL;
  
  
  int counter = 0;
       
  for( ii = 0; ii< ref_volumes.size(); ii++ )
  {
    curr_volume = ref_volumes.get_and_step();
    GeomMeasureTool::measure_volume_volume_and_hydraulic_radius(curr_volume,
                                                                volume_volume,
                                                                volume_hydraulic_radius);
    if ( curr_volume->num_ref_faces() > 1 )
    {
      GeomMeasureTool::angles_between_volume_surfaces(curr_volume,
                                                      curr_min_angle,
                                                      curr_max_angle,
                                                      face_min_1, face_min_2,
                                                      face_max_1, face_max_2);
      GeomMeasureTool::find_adjacent_face_ratios(curr_volume,
                                                 curr_face_ratio,
                                                 big_face, small_face);
    }
    total_volume += volume_volume;
    counter++;
    if ( volume_volume < min_volume )
    {
      volume_min_volume = curr_volume;
      min_volume = volume_volume;
    }
    if ( volume_hydraulic_radius < min_volume_hydraulic_radius )
    {
      volume_min_hydraulic_radius = curr_volume;
      min_volume_hydraulic_radius = volume_hydraulic_radius;
    }
    if ( curr_volume->num_ref_faces() > 1 &&  curr_min_angle < min_angle )
    {
      volume_min_angle = curr_volume;
      face_min_angle_1 = face_min_1;
      face_min_angle_2 = face_min_2;
      min_angle = curr_min_angle;
    }
    if ( curr_volume->num_ref_faces() > 1 && curr_max_angle > max_angle )
    {
      volume_max_angle = curr_volume;
      max_angle = curr_max_angle;
      face_max_angle_1 = face_max_1;
      face_max_angle_2 = face_max_2;
    }
    if ( curr_volume->num_ref_faces() > 1 && curr_face_ratio > max_face_ratio )
    {
      volume_face_ratio = curr_volume;
      max_face_ratio = curr_face_ratio;
      max_small_face = small_face;
      max_big_face = big_face;
    }
  }
  if ( counter > 0 )
  {
    PRINT_INFO("*************Volume Summary************\n\n");
    PRINT_INFO("Measured %d Volumes: Average Volume was %f.\n\n", 
               counter, total_volume/counter);
    PRINT_INFO("Minimum Volume was %f in volume %d\n\n",
               min_volume, volume_min_volume->id());
    PRINT_INFO("Minimum Volume Hydraulic Radius (volume/surface area) was %f in volume %d\n\n",
               min_volume_hydraulic_radius, volume_min_hydraulic_radius->id());
    if ( volume_min_angle )
      PRINT_INFO("Minimum Interior Angle Between Two Surfaces was %f,\n"
                 "\tsurfaces %d and %d in volume %d\n\n",
                 min_angle, face_min_angle_1->id(), face_min_angle_2->id(), volume_min_angle->id());
    if ( volume_max_angle )
      PRINT_INFO("Maximum Interior Angle Between Two Surfaces was %f,\n"
                 "\tsurfaces %d and %d in volume %d\n\n",
                 max_angle, face_max_angle_1->id(), face_max_angle_2->id(), volume_max_angle->id());
    if ( volume_face_ratio )
      PRINT_INFO("Maximum Area Ratios Between Adjacent Faces was %f,\n"
                 "\tsurfaces %d and %d in volume %d.\n\n",
                 max_face_ratio, max_big_face->id(), max_small_face->id(), volume_face_ratio->id());
    int merged, unmerged;
    double ratio;     
    GeomMeasureTool::merged_unmerged_surface_ratio(ref_volumes, merged, unmerged, ratio);
    PRINT_INFO("Total number of merged surfaces:\t%d\n"
               "Total number of unmerged surfaces:\t%d\n"
               "Ratio of merged to unmerged surfaces:\t%f\n\n", merged, unmerged, ratio);
    RefGroup *group_of_shells;
    int number_shells;
    GeomMeasureTool::find_shells(ref_volumes, group_of_shells, number_shells);
    if ( number_shells > 0 )
    {
      if ( number_shells == 1 )
      {
        PRINT_INFO("Found %d shell of connected surfaces.\n", number_shells);
        PRINT_INFO("These shell is contained in the group %s\n",
                   group_of_shells->entity_name().c_str());
      }
      else
      {
        PRINT_INFO("Found %d shells of connected surfaces.\n", number_shells);
        PRINT_INFO("These shells are contained in the group %s\n",
                   group_of_shells->entity_name().c_str());
      }
    }
    PRINT_INFO("\n\n");
    PRINT_INFO("Now Printing Summaries of Surfaces contained by these volumes.\n");
    GeomMeasureTool::print_surface_measure_summary(ref_faces);
  }
  return;

}

// Find all of the surfaces in the given volumes that have narrow regions.
void GeomMeasureTool::find_surfs_with_narrow_regions(DLIList <RefVolume*> &ref_vols,
                                          double tol,
                                          DLIList <RefFace*> &surfs_with_narrow_regions)
{
  int i, j;
  double tol_sq = tol*tol;
  for(i=ref_vols.size(); i--;)
  {
    RefVolume *cur_vol = ref_vols.get_and_step();
    DLIList<RefFace*> vol_faces;
    cur_vol->ref_faces(vol_faces);
    for(j=vol_faces.size(); j--;)
    {
      RefFace *cur_face = vol_faces.get_and_step();
      if(has_narrow_region(cur_face, tol, tol_sq))
      {
        surfs_with_narrow_regions.append_unique(cur_face);
      }
    }
  }
}

// Checks to see if at the given position the two edges are close together and 
// have the same tangent.
int GeomMeasureTool::is_narrow_region_at_point(RefEdge *e1,
                                               const CubitVector &pt_on_e1,
                                               RefEdge *e2,
                                               const double &tol_sq)
{
  int ret = 0;
  CubitVector closest, tan_1, tan_2;
  e2->closest_point_trimmed(pt_on_e1, closest);
  double dist = (pt_on_e1-closest).length_squared();
  if(dist < tol_sq)
  {
    e1->tangent(pt_on_e1, tan_1);
    e2->tangent(closest, tan_2);
    tan_1.normalize();
    tan_2.normalize();
    if(fabs(tan_1 % tan_2) > .9)
      ret = 1;
  }
  return ret;
}

// Checks to see if the face has a narrow region in it.  A narrow
// region here is defined as a region where two of the edges on the
// face come close enough to each other that a narrow region is 
// formed.
int GeomMeasureTool::has_narrow_region(RefFace* face, const double &tol,
                                       const double &tol_sq)
{
  int i, j, k;
  int ret = 0;
  double small_step = 5.0*tol, step;
  RefEdge *edge1, *edge2;
  CubitVector closest1, closest2, closest3;
  CubitVector p3, p2, p1, tan11, tan21, step_pos;
  RefVertex *cur_e1_vert;

  // Loop through all of the edges in the face and compare each
  // one with every other edge in the face.
  DLIList<RefEdge*> edges;
  face->ref_edges(edges);
  while(edges.size() > 1 && !ret)
  {
    // Remove the current edge each time so that we aren't
    // doing redundant comparisons.
    RefEdge *cur_edge = edges.extract();

    // Compare this edge with the remaining edges on the face.
    for(k=edges.size(); k && !ret; k--)
    {
      RefEdge *other_edge = edges.get_and_step();

      // Loop through twice to compare each edge with
      // the other.  One edge will be used as the source
      // edge and points from the source edge will be 
      // compared with points on the target edge.
      for(i=0; i<2 && !ret; i++)
      {
        if(i==0)
        {
          edge1 = cur_edge;
          edge2 = other_edge;
        }
        else
        {
          edge1 = other_edge;
          edge2 = cur_edge;
        }

        DLIList<RefVertex*> e1_verts;
        edge1->ref_vertices(e1_verts);
        for(j=e1_verts.size(); j && !ret; j--)
        {
          cur_e1_vert = e1_verts.get_and_step();
          p1 = cur_e1_vert->coordinates();

          // In trying to determine whether the face has a narrow region
          // we will simply look for one point on the source edge that is
          // close to the target edge and then step away from it a little
          // and check another point.  This is working from the ends of
          // the source edge which will not catch every case but should 
          // catch most of the cases.
          if(is_narrow_region_at_point(edge1, p1, edge2, tol_sq))
          {
            step = small_step;
            if(cur_e1_vert == edge1->end_vertex())
              step *= -1.0;
            edge1->point_from_arc_length(p1, step, p2);
            if(is_narrow_region_at_point(edge1, p2, edge2, tol_sq))
              ret = 1;
          }
        }
      }
    }
  }

  return ret;
}

void GeomMeasureTool::find_small_curves( DLIList <RefVolume*> &ref_vols,
                                         double tol,
                                         DLIList <RefEdge*> &small_curves,
                                         DLIList <double> &small_lengths)
{
  int ii, jj;
  DLIList <RefEdge*> ref_edges, temp_edges;
  RefVolume *ref_vol;
  RefEdge *curr_edge;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_edges(temp_edges);
      //uniquely add the edges.
    for ( jj = temp_edges.size(); jj > 0; jj-- )
    {
      curr_edge = temp_edges.pop();
      if ( curr_edge->marked()== 0 )
      {
        curr_edge->marked(1);
        ref_edges.append(curr_edge);
      }
    }
  }
  double length;
    //find the small curves and reset the marked flag.
  for ( ii = ref_edges.size(); ii > 0; ii-- )
  {
    curr_edge = ref_edges.get_and_step();

    //skip point curves
    if( curr_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;

    curr_edge->marked(0);
    length = curr_edge->measure();
    if ( length <= tol )
    {
      small_curves.append(curr_edge);
      small_lengths.append(length);
    }
  }
  return;
}
void GeomMeasureTool::find_small_faces( DLIList <RefVolume*> &ref_vols,
                                        double tol,
                                        DLIList <RefFace*> &small_faces)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
      //uniquely add the faces.
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      if ( curr_face->marked()== 0 )
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
  double area;
    //find the small curves and reset the marked flag.
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
      //reset the mark.
    curr_face->marked(0);
    area = measure_area(curr_face);
    if ( area <= tol )
      small_faces.append(curr_face);
  }
  return;
}
void GeomMeasureTool::find_small_faces_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                         double tol,
                                                         DLIList <RefFace*> &small_faces,
                                                         DLIList <double> &small_hyd_rad,
                                                         double &percent_planar,
                                                         double &percent_pl_co)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
      //uniquely add the faces.
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      if ( curr_face->marked()== 0 )
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
  double area;
  double length;
  int total_faces = 0;
  int total_plane = 0;
  int total_cone = 0;
  double area_plane = 0.0;
  double total_area = 0.0;
    //find the small curves and reset the marked flag.
  DLIList <CoEdge*> co_edges;
  int num_faces = ref_faces.size();
  ProgressTool *progress_ptr = NULL;
  char title[32];
  if (num_faces > 20)
  {
    strcpy(title, "Small Surface Progress");
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, title, NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_face->marked(0);
    area = measure_area(curr_face);
    total_area += area;
    total_faces++;
      //update the progress..
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_faces))/((double)(num_faces));
      progress_ptr->percent(curr_percent);
    }
    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }
    if ( curr_face->geometry_type() == PLANE_SURFACE_TYPE )
    {
      total_plane++;
      area_plane += area;
    }
    else if ( curr_face->geometry_type() == CONE_SURFACE_TYPE )
      total_cone++;
    co_edges.clean_out();
    curr_face->co_edges(co_edges);
    length = 0.0;
    for ( jj = co_edges.size(); jj > 0; jj-- )
    {
      length += co_edges.get_and_step()->get_ref_edge_ptr()->measure();
    }
      //reset the mark.
    curr_face->marked(0);
      //compute the hydraulic radius 4*(A/P).
    if ( length <= CUBIT_RESABS )
    {
      PRINT_INFO("Total Perimiter Length of Surface %d is less than tolerance.\n",
                 curr_face->id());
      
      continue;
    }
    area = 4.0*(area/length);
    if ( area <= tol )
    {
      small_faces.append(curr_face);
      small_hyd_rad.append(area);
    }
  }
  if ( total_faces > 0 )
  {
    percent_planar = 100.0*area_plane/total_area;
    percent_pl_co = 100.0*((double)total_plane+(double)total_cone)/(double)total_faces;
  }
  else
  {
    percent_planar = 0.;
    percent_pl_co = 0.;
  }
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
  
  return;
}
void GeomMeasureTool::find_small_volumes( DLIList <RefVolume*> &ref_vols,
                                          double tol,
                                          DLIList <RefVolume*> &small_volumes)
{
  int ii;
  RefVolume *ref_vol;
  double volume;
  
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    volume = ref_vol->measure();
    if ( volume <= tol )
      small_volumes.append(ref_vol);
  }
  return;
}
void GeomMeasureTool::find_small_volumes_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                           double tol,
                                                           DLIList <RefVolume*> &small_volumes,
                                                           DLIList <double> &small_hyd_rad)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  double area;
  double hyd_volume;
  
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
    area = 0.0;
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      area += measure_area(curr_face);
    }
    if ( area <= CUBIT_RESABS )
      continue;
    hyd_volume = 6.0*(ref_vol->measure() / area);
    if ( hyd_volume <= tol )
    {
      small_hyd_rad.append(hyd_volume);
      small_volumes.append(ref_vol);
    }
  }
  return;
}
///
/// Find the tangential meetings in the volume.
/// This specifically looks for surfaces that meet tangentially
/// that would be a problem.  Usually these are surfaces that
/// come into a side 180 degrees and on top there is a
/// sharpe angle.  Usually if there isn't a sharpe curve
/// angle these meetings are not bad.
/// Note that this function assumes that you are passing
/// in sets of curve edges that have small interior angles between them.
/// It also assumes that the edge pairs are ordered as they would
/// be found in the surfaces (first edge then next edge).
/// Basically, call the funciton, find_interior_curve_angles
/// before calling this function, and pass this function those results.
///
void GeomMeasureTool::find_sharp_tangential_meets( RefVolume *ref_volume,
                                                   DLIList <RefEdge*> &small_angle_edge_pairs,
                                                   DLIList <RefFace*> &tangential_surface_pairs )
{
    //Okay, given the small angle edge pairs.  See if there are surfaces that meet tangentially.
    //I'm really looking for something spefic.
  RefEdge *ref_edge_1, *ref_edge_2;
  int ii, jj;
  DLIList <Loop*> loops_1, loops_2;
  Loop *common_loop;
  const double rad_to_deg = 180.0/CUBIT_PI;
  RefFace *ref_face;
  
    //loop over for each edge pair.
  for ( ii = 0; ii < small_angle_edge_pairs.size(); ii += 2 )
  {
    ref_edge_1 = small_angle_edge_pairs.get_and_step();
    ref_edge_2 = small_angle_edge_pairs.get_and_step();
      //find the loop that these edges are both on.
    loops_1.clean_out();
    loops_2.clean_out();
    ref_edge_1->loops(loops_1);
    ref_edge_2->loops(loops_2);
    common_loop = NULL;
    for ( jj = loops_1.size(); jj > 0; jj-- )
    {
      if ( !loops_2.move_to(loops_1.get()) )
        loops_1.remove();
      else
        loops_1.step();
    }
    if ( loops_1.size() == 1 )
      common_loop = loops_1.get();
    else
    {
        //find the loop that has a co_edge from ref_edge_1
        //then a co_edge from ref_edge_2.
      DLIList <CoEdge*> co_edges;
      Loop *tmp_loop;
      int kk;
      CoEdge *co_edge;
      for ( jj = loops_1.size(); jj > 0; jj-- )
      {
        tmp_loop = loops_1.get_and_step();
        co_edges.clean_out();
        tmp_loop->ordered_co_edges(co_edges);
        for ( kk = co_edges.size(); kk > 0; kk-- )
        {
          co_edge = co_edges.get_and_step();
          if ( co_edge->get_ref_edge_ptr() == ref_edge_1 )
          {
            if ( co_edges.get()->get_ref_edge_ptr() == ref_edge_2 )
            {
              common_loop = tmp_loop;
              break;
            }
            else
                //not this loop..
              break;
          }
        }
        if ( common_loop != NULL )
          break;
      }
    }
    assert ( common_loop != NULL );
      //now find surface.
      //Then see if there are surfaces that meet at 180 degrees that
      //are 90* rotated from here.
    ref_face = common_loop->get_ref_face_ptr();
      //Okay, refface is the surface on which we have the bad angle.
      //One of these curves will have a concave angle and the other a
      //convex angle.  The surface on the other side of the convex angle
      //we need to check to see if it comes tangentially into another
      //surface (dihedral angle of 180.)
    RefFace *ref_face_1 = ref_edge_1->other_face(ref_face, ref_volume);
    if( ref_face_1 == NULL )
      continue;

    double surf_angle_1, surf_angle_2;
    RefFace *tangential_face = NULL;
    RefEdge *edge_next_to = NULL;
    surf_angle_1 =  GeometryQueryTool::instance()->surface_angle(ref_face,
                                                                 ref_face_1,
                                                                 ref_edge_1,
                                                                 ref_volume);
    surf_angle_1 *= rad_to_deg;
    RefFace *ref_face_2;
    ref_face_2 = ref_edge_2->other_face(ref_face, ref_volume);
    if( ref_face_2 == NULL )
      continue;
    surf_angle_2 = GeometryQueryTool::instance()->surface_angle(ref_face,
                                                                ref_face_2,
                                                                ref_edge_2,
                                                                ref_volume);
    surf_angle_2 *= rad_to_deg;
    if ( surf_angle_1 < 180.0 && surf_angle_2 > 180.0 )
    {
      tangential_face = ref_face_1;
      edge_next_to = ref_edge_1;
    }
    else if ( surf_angle_2 < 180.0 && surf_angle_1 > 180.0 )
    {
      tangential_face = ref_face_2;
      edge_next_to = ref_edge_2;
    }
    if ( tangential_face == NULL ) // This isn't the config we seek.
      continue;

      //Now, on this face, there is an edge, next to edge_next_to,
      //that will have an interior angle of 180.0.  If that is the
      //case then this is what we seek.
    RefVertex* ref_vert = ref_edge_1->common_ref_vertex(ref_edge_2,
                                                        ref_face);
    DLIList <RefEdge*> tmp_ref_edges;
    ref_vert->ref_edges(tmp_ref_edges);
      //Find the edge that is on tangential face and is not edge_next_to.
    RefEdge *tangential_edge=NULL, *tmp_edge;
    for ( jj = tmp_ref_edges.size(); jj > 0; jj-- )
    {
      tmp_edge = tmp_ref_edges.get_and_step();
      if ( tmp_edge != edge_next_to &&
           tmp_edge->is_directly_related(tangential_face) )
      {
        tangential_edge = tmp_edge;
        break;
      }
    }
    if ( tangential_edge == NULL )
        //not what we are looking for.
      continue;
      //Now get the other face and measure the dihedral angle.
    RefFace *other_face = tangential_edge->other_face(tangential_face,
                                                      ref_volume);
    double tan_angle;
    tan_angle = GeometryQueryTool::instance()->surface_angle(tangential_face,
                                                             other_face,
                                                             tangential_edge,
                                                             ref_volume);
    tan_angle *= rad_to_deg;
      //If they are within tolerance, then we have it...
    if ( tan_angle > 165.0 && tan_angle < 195.0 )
    {
      tangential_surface_pairs.append(tangential_face);
      tangential_surface_pairs.append(other_face);
    }
  }
}
    
      
    
void GeomMeasureTool::find_interior_curve_angles( RefVolume *ref_volume,
                                                  double upper_bound,
                                                  double lower_bound,
                                                  DLIList <RefEdge*> &large_edge_angles,
                                                  DLIList <RefEdge*> &small_edge_angles,
                                                  DLIList <double> &large_angles,
                                                  DLIList <double> &small_angles,
                                                  int &total_interior,
                                                  int &total_fuzzy)
{
  int ii, jj, kk;
#ifdef BOYD17
  DLIList <RefFace*> ref_faces, temp_faces;
#endif
  DLIList <RefFace*> ref_faces;
  RefFace *curr_face;
  total_interior = 0;
  total_fuzzy = 0;
  ref_volume->ref_faces(ref_faces);
  
    //Okay now loop over the edge loops of the reffaces
    //and find the interior angles.
  DLIList<Loop*> loops;
  Loop *curr_loop;
  DLIList <CoEdge*> co_edges;
  CoEdge *curr_edge, *next_edge;
  double angle = 0.0;
  const double RAD_TO_DEG = 180.0/CUBIT_PI;
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step(); 
    curr_face->loops(loops);
      //Now loop over these loops.
    for ( jj = loops.size(); jj > 0; jj-- )
    {
      curr_loop = loops.pop();
        //now loop over the edges in this
      co_edges.clean_out();
      curr_loop->ordered_co_edges(co_edges);
      for ( kk = co_edges.size(); kk > 0; kk-- )
      {
        curr_edge = co_edges.get_and_step();
        next_edge = co_edges.get();
	   
        if ( curr_edge == next_edge )
        {
            //This is a hardpoint
          if ( curr_edge->get_ref_edge_ptr()->geometry_type() == POINT_CURVE_TYPE )
          {
            angle = 360.0;
              //for now ignore this..
            continue;
          }
          else
          {
            RefEdge *the_edge = curr_edge->get_ref_edge_ptr();
              //Get points 1% from the ends in both directions.
            CubitVector point_start, point_end;
            CubitStatus stat = the_edge->position_from_fraction(0.01,
                                                                point_start);
            if ( stat != CUBIT_SUCCESS )
              angle = 180.0;
            else
            {
              stat = the_edge->position_from_fraction(0.99,
                                                      point_end);
              if ( stat != CUBIT_SUCCESS )
                angle = 180.0;
              else
              {
                CubitVector normal = curr_face->normal_at(point_start, NULL);
                CubitVector tangent_1, tangent_2;
                stat = the_edge->tangent( point_start, tangent_1, curr_face);
                CubitStatus stat2 = the_edge->tangent( point_end, tangent_2, curr_face);
                if ( stat != CUBIT_SUCCESS || stat2 != CUBIT_SUCCESS )
                  angle = 180.0;
                else
                {
                    //one of the tangents will be pointing the wrong
                    //direction, so reverse it based on the sense.
                  if ( curr_edge->get_sense() == CUBIT_REVERSED )
                    tangent_1 = -tangent_1;
                  else
                    tangent_2 = -tangent_2;
                  angle = normal.vector_angle(tangent_2, tangent_1);
                  angle *= RAD_TO_DEG;
                }
              }
            }
          }
        }
        else if ( curr_edge->get_ref_edge_ptr() == next_edge->get_ref_edge_ptr() )
        {
            //This is the end of a hard line. Set this angle because
            //it is tricky to measure it.
          angle = 360.0;
        }
        else if ( curr_edge != next_edge )
        {
          angle =
            GeometryQueryTool::instance()->geometric_angle(curr_edge,
                                                           next_edge);
          angle *= RAD_TO_DEG;
        }
        total_interior++;
          //count the number of angles that are fuzzy.
        if ( angle <= GEOM_END_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_END_UPPER &&
                  angle <= GEOM_SIDE_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_SIDE_UPPER &&
                  angle <= GEOM_CORNER_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_CORNER_UPPER )
          total_fuzzy++;
        if ( angle <= lower_bound )
        {
          small_edge_angles.append(curr_edge->get_ref_edge_ptr());
          small_edge_angles.append(next_edge->get_ref_edge_ptr());
          small_angles.append(angle);
        }
        if ( angle >= upper_bound )
        {
          large_edge_angles.append(curr_edge->get_ref_edge_ptr());
          large_edge_angles.append(next_edge->get_ref_edge_ptr());
          large_angles.append(angle);
        }
      }
    }
  }
  return;
}

void GeomMeasureTool::find_dihedral_angles( DLIList <RefVolume*> &ref_vols,
                                            double upper_bound,
                                            double lower_bound,
                                            DLIList <RefFace*> &large_face_angles,
                                            DLIList <RefFace*> &small_face_angles,
                                            DLIList <double> &large_angles,
                                            DLIList <double> &small_angles,
                                            int &total_interior,
                                            int &total_fuzzy,
                                            int &total_not_flat)
{
  int ii;
#ifdef BOYD17
  DLIList <RefFace*> ref_faces, temp_faces;
#endif
  RefVolume *ref_vol;
  total_interior = 0;
  total_fuzzy = 0;
  total_not_flat = 0;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    find_dihedral_angles(ref_vol, upper_bound,
                         lower_bound,
                         large_face_angles,
                         small_face_angles,
                         large_angles, small_angles,
                         total_interior, total_fuzzy, total_not_flat);
  }
  return;
}
void GeomMeasureTool::find_dihedral_angles(RefVolume *curr_volume,
                                           double upper_bound,
                                           double lower_bound,
                                           DLIList <RefFace*> &large_face_angles,
                                           DLIList <RefFace*> &small_face_angles,
                                           DLIList <double> &large_angles,
                                           DLIList <double> &small_angles,
                                           int &total_interior,
                                           int &total_fuzzy, int &total_not_flat)
{
  int ii, jj, shared_angles = 0;
  double common_angle;
  DLIList <RefFace*> face_list;
  RefFace *curr_face_1, *curr_face_2;
  RefFace *tmp_face;
  RefEdge *common_edge;
  DLIList <RefVolume*> ref_vols;
  DLIList<RefEdge*> vol_edges;
  curr_volume->ref_edges(vol_edges);

  //maybe the body is a sheet-body...don't print errors if so
  bool is_sheet_body = false;
  if( curr_volume->measure() == 0 )
    is_sheet_body = true;

  for ( ii = vol_edges.size(); ii > 0; ii-- )
  {
    common_edge = vol_edges.get_and_step();

    //skip point curves
    if( common_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;

    curr_face_1 = NULL;
    curr_face_2 = NULL;
    face_list.clean_out();
    common_edge->ref_faces(face_list);
      //find the two faces that are part of this volume.
    for ( jj = face_list.size(); jj > 0; jj-- )
    {
      tmp_face = face_list.get_and_step();
      ref_vols.clean_out();
      tmp_face->ref_volumes(ref_vols);
      if ( ref_vols.move_to(curr_volume) )
      {
        if ( curr_face_1 == NULL )
          curr_face_1 = tmp_face;
        else if ( curr_face_2 == NULL )
        {
          curr_face_2 = tmp_face;
          break;
        }
      }
    }
    if ( curr_face_1 == NULL ||
         curr_face_2 == NULL )
    {
      if( is_sheet_body == false )
        PRINT_ERROR("Problems finding connected surfaces to a curve\n"
                   "\tfor measuring dihedral angles.\n");
      continue;
    }
    shared_angles++;
    common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1,
                                                                 curr_face_2,
                                                                 common_edge,
                                                                 curr_volume);
    common_angle *= 180.0/CUBIT_PI;
    total_interior++;
      //use hard coded angles because we want these a little looser.
    if ( common_angle >= 200.0 ||
         common_angle <= 160.0 )
      total_not_flat++;
    
    if ( curr_face_1->geometry_type() != PLANE_SURFACE_TYPE ||
         curr_face_2->geometry_type() != PLANE_SURFACE_TYPE )
      total_fuzzy++;
    else if ( common_angle <= GEOM_END_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_END_UPPER &&
              common_angle <= GEOM_SIDE_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_SIDE_UPPER &&
              common_angle <= GEOM_CORNER_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_CORNER_UPPER )
      total_fuzzy++;
    
    if( common_angle <= lower_bound )
    {
      small_face_angles.append(curr_face_1);
      small_face_angles.append(curr_face_2);
      small_angles.append(common_angle);
    }
    if( common_angle >= upper_bound )
    {
      large_face_angles.append(curr_face_1);
      large_face_angles.append(curr_face_2);
      large_angles.append(common_angle);
    }
  }
    //clean up the edge marks.
  for ( ii = vol_edges.size(); ii > 0; ii-- )
    vol_edges.get_and_step()->marked(0);
}
void GeomMeasureTool::find_close_loops(DLIList <RefVolume*> &ref_vols,
                                       DLIList <RefEdge*> &close_edges,
                                       DLIList <RefFace*> &close_loop_faces,
                                       DLIList <double> &small_lengths,
                                       double tol)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
      //uniquely add the faces.
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      if ( curr_face->marked()== 0 )
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
  DLIList <RefEdge*> tmp_close_edges;
  DLIList <double> tmp_small_lengths;
  
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_face->marked(0);
    if ( curr_face->num_loops() < 2 )
      continue;
    tmp_close_edges.clean_out();
    tmp_small_lengths.clean_out();
    find_close_loops( curr_face, tmp_close_edges,
                      tmp_small_lengths, tol);
    if ( tmp_close_edges.size() > 0 )
    {
      assert(tmp_close_edges.size()%2 == 0 );
      tmp_close_edges.reset();
      tmp_small_lengths.reset();
      for ( jj = tmp_close_edges.size()/2; jj > 0; jj-- )
      {
        close_edges.append(tmp_close_edges.get_and_step());
        close_edges.append(tmp_close_edges.get_and_step());
        small_lengths.append(tmp_small_lengths.get_and_step());
        close_loop_faces.append(curr_face);
      }
    }
  }
  return;
}

  
void GeomMeasureTool::find_close_loops(RefFace *face,
                                       DLIList <RefEdge*> &close_edges,
                                       DLIList <double> &small_lengths,
                                       double tol)
{
    //only do faces with multiple loops since we are measuring the
    //distances between these loops.
  if ( face->num_loops() < 2 )
    return;
    //This is the most tricky of these functions.
    //To do this I'm going to facet the edges, then use an AbstractTree to
    //find the minimum distance between the loops and between a single
    //loop.
  PointLoopList boundary_point_loops;
  CubitStatus stat;
  stat = get_boundary_points(face, boundary_point_loops, GEOMETRY_RESABS*100);
  if ( stat != CUBIT_SUCCESS )
    return;
  
  SegLoopList boundary_seg_loops;
  stat = convert_to_lines( boundary_point_loops,
                           boundary_seg_loops,
                           face);
  if ( stat != CUBIT_SUCCESS )
    return;

    //Now add the points to the AbstractTree.
  DLIList<AbstractTree<GeomSeg*>*> atree_list;
  AbstractTree<GeomSeg*> *curr_tree;
  SegList *seg_list;
  GeomSeg  *curr_seg;
  int ii,jj, kk;
  //const double angle_convert = 180.0/CUBIT_PI;
  boundary_seg_loops.reset();
  int num_segs = 0;
  for (ii = 0; ii < boundary_seg_loops.size(); ii++ )
  {
    seg_list = boundary_seg_loops.get_and_step();
    curr_tree = new RTree<GeomSeg*>(GEOMETRY_RESABS);
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
        //build the r-tree.
      num_segs++;
      curr_seg = seg_list->get_and_step();
      curr_tree->add(curr_seg);
        //calculate the interior angles.
    }
    atree_list.append(curr_tree);
  }
    //determine the minimum distance between the loops.
  PointList *curr_points;
  GeomPoint *curr_point;
  DLIList <GeomSeg*> nearest_neighbors;
  RefEdge *closest_edge_1, *closest_edge_2;
  CubitVector curr_loc;
  double closest_dist;
  double min_for_loop;
  atree_list.reset();
  for ( ii = 0; ii < atree_list.size(); ii++ )
  {
    curr_points = boundary_point_loops.get_and_step();
    for ( jj = ii+1; jj < atree_list.size(); jj++ )
    {
      min_for_loop = CUBIT_DBL_MAX;
      closest_edge_1 = NULL;
      closest_edge_2 = NULL;
      curr_tree = atree_list.next(jj);
        //now for every point in curr_points, find it's closest_point
        //in the curr_tree.
      for ( kk = 0; kk < curr_points->size(); kk++ )
      {
        curr_point = curr_points->get_and_step();
        curr_loc.set(curr_point->coordinates());
        nearest_neighbors.clean_out();
        curr_tree->k_nearest_neighbor(curr_loc,
                                      1, closest_dist, nearest_neighbors,
                                      dist_sq_point_data);
        if ( nearest_neighbors.size() == 0)
        {
            //hmm, not sure what is wrong here.
          PRINT_ERROR("Can't find closest point between loops.\n");
        }
        if (closest_dist <= tol*tol && closest_dist < min_for_loop*min_for_loop )
        {
            //just store the smaller ones.  Don't store
          RefEntity *ref_ent = curr_point->owner();
          GeomSeg *near_seg = nearest_neighbors.get();
          RefEdge *ref_edge_1 = CAST_TO(ref_ent, RefEdge);
          if ( ref_edge_1 == NULL )
          {
              //assume this is a vertex.
            RefVertex *ref_vert = CAST_TO(ref_ent, RefVertex);
            if ( ref_vert == NULL )
            {
              PRINT_ERROR("problem with point ownership in closest loops.\n");
              continue;
            }
            DLIList <RefEdge*> ref_edges;
            ref_vert->ref_edges(ref_edges);
            int ll;
            RefEdge *the_edge = NULL;
            for ( ll = ref_edges.size(); ll > 0; ll-- )
            {
              if ( ref_edges.get()->is_directly_related(face) )
              {
                the_edge = ref_edges.get();
                break;
              }
              ref_edges.step();
            }
            if ( the_edge == NULL )
            {
              PRINT_ERROR("problem with point ownership in closest loops.\n");
              continue;
            }
            ref_edge_1 = the_edge;
          }
          ref_ent = near_seg->owner();
          RefEdge *ref_edge_2 = CAST_TO(ref_ent, RefEdge);
          if ( ref_edge_2 == NULL )
          {
            PRINT_ERROR("problem with point ownership in closest loops.\n");
            continue;
          }
          
          min_for_loop = sqrt(closest_dist);
          closest_edge_1 = ref_edge_1;
          closest_edge_2 = ref_edge_2;
        }
      }
      if ( closest_edge_1 != NULL )
      {
        close_edges.append(closest_edge_1);
        close_edges.append(closest_edge_2);
        small_lengths.append(min_for_loop);
      }
    }
  }
    //clean up the memory.
  int ll,mm;
  atree_list.reset();
  boundary_seg_loops.reset();
  for ( ll = boundary_seg_loops.size(); ll > 0; ll-- )
  {
    delete atree_list.pop();
    seg_list = boundary_seg_loops.pop();
    for ( mm = seg_list->size(); mm > 0; mm-- )
      delete seg_list->pop();
    delete seg_list;
  }
		
  PointList *point_list;
  for ( ll = boundary_point_loops.size(); ll > 0; ll-- )
  {
    point_list = boundary_point_loops.pop();
    for ( mm = point_list->size(); mm > 0; mm-- )
    {
      GeomPoint *tmp_point = point_list->pop();
        //delete the CubitPointData
      delete tmp_point;
    }
    delete point_list;
  }
  return;
}
double GeomMeasureTool::measure_area(RefFace* curr_face)
{
  double area;
  GeomDataObserver *g_d_o = GeomDataObserver::get(curr_face);
  if ( g_d_o == NULL )
    g_d_o = GeomDataObserver::create(curr_face);

  if ( g_d_o == NULL )
  {
    assert(g_d_o != NULL );
    return curr_face->measure();
  }
  if ( !g_d_o->is_measure_set() )
  {
    area  = (curr_face->get_surface_ptr())->measure();
    g_d_o->set_measure(area);
  }
  area = g_d_o->get_measure();
  return area;
}

///
///  This function simply gets the bad entities of the volume.  This
///  assumes the volume is from some solid modelar where this function
///  is defined. If not, it will just be empty...
///
void GeomMeasureTool::find_bad_geometry(RefVolume *volume,
                                        DLIList <RefEntity*> &bad_ents)
{
  Lump *lump_ptr = volume->get_lump_ptr();
  DLIList <TopologyEntity*> bad_top_ents;
  lump_ptr->validate(volume->entity_name(), bad_top_ents);
    //Now go through and get a RefEnttiy for each entity.
  RefEntity *ref_ent;
  TopologyEntity *topo_ent;
  GroupingEntity *gr_ent;
  SenseEntity *se_ent;
  
  int ii;
  for (ii = bad_top_ents.size(); ii > 0; ii-- )
  {
    topo_ent = bad_top_ents.pop();
    ref_ent = CAST_TO(topo_ent, RefEntity);
    if ( ref_ent != NULL )
    {
      bad_ents.append(ref_ent);
      continue;
    }
    gr_ent = CAST_TO(topo_ent, GroupingEntity);
    if ( gr_ent != NULL )
    {
      Chain *ch = CAST_TO(gr_ent, Chain);
      if ( ch != NULL )
      {
        RefVertex *start_v = ch->start_vertex();
        RefVertex *end_v = ch->end_vertex();
        RefEdge *ref_edge = start_v->common_ref_edge(end_v);
        bad_ents.append(ref_edge);
        continue;
      }
      Loop *lo = CAST_TO(gr_ent, Loop);
      if ( lo != NULL )
      {
        RefFace *rf = lo->get_ref_face_ptr();
        bad_ents.append(rf);
        continue;
      }
      Shell *sh = CAST_TO(gr_ent, Shell);
      if (sh != NULL )
      {
        RefVolume *rv = sh->get_ref_volume_ptr();
        bad_ents.append(rv);
        continue;
      }
      else
      {
          //this is incase there is some other entity we don't
          //know about, just get the body...
        Body *body = topo_ent->body();
        bad_ents.append(body);
        continue;
      }
    }
    se_ent = CAST_TO(topo_ent, SenseEntity);
    if ( se_ent != NULL )
    {
      BasicTopologyEntity *be = se_ent->get_basic_topology_entity_ptr();
      bad_ents.append(be);
      continue;
    }
      //If we are here we didn't get the right entity.  Just get
      //the body associated with this topology entity.
    Body *b = topo_ent->body();
    bad_ents.append(b);
  }
  return;
}

///
///  Find the irregular vertex valences.
///  Find things like vertices with valences greater than 4.
///  Assume for now that the volumes are not merged..
///
void GeomMeasureTool::find_irregular_valence( DLIList <RefVolume*> &ref_volumes,
                                              DLIList <RefVertex*> &irregular_vertices)
{
  RefVolume *ref_vol;
  int ii;
  for ( ii = 0; ii < ref_volumes.size(); ii++ )
  {
    ref_vol = ref_volumes.get_and_step();
    find_irregular_valence(ref_vol, irregular_vertices);
  }
}
///
///  Find the irregular vertex valences.
///  Find things like vertices with valences greater than 4.
///  Assume for now that the volumes are not merged..
///
void GeomMeasureTool::find_irregular_valence( RefVolume *ref_volume,
                                              DLIList <RefVertex*> &irregular_vertices)
{
  RefVertex *ref_vert;
  DLIList <RefVertex*> ref_vertices;
  ref_volume->ref_vertices(ref_vertices);
  int ii;
  for ( ii = ref_vertices.size(); ii > 0; ii-- )
  {
    ref_vert = ref_vertices.get_and_step();
    if ( ref_vert->num_ref_edges() > 4 )
      irregular_vertices.append(ref_vert);
  }
}

///
///  Find fillets and rounds.
///
void GeomMeasureTool::find_blends( RefVolume *ref_volume,
                                   DLIList <RefFace*> &blend_faces,
                                   DLIList <DLIList<RefFace*>*> &blend_groups )
{
    //Assume for now that blends have cartesian angles between
    //all the attached faces and itself.
    //Also a blend cannot be plannar.  And directly accross from
    //one blend edge, there is usually another blended edge but
    //the angle between the two normals must be orthogonal. In other
    //words their must be some sort of transition.
  
    //First go through each face and each edge on each face.
  int ii, jj;
  DLIList <RefFace*> ref_faces;
  DLIList <RefEdge*> ref_edges;
  RefFace *ref_face;
  RefEdge *ref_edge, *other_edge;
  ref_volume->ref_faces(ref_faces);
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    ref_face = ref_faces.get_and_step();
      //Test the face to see if it is a blend.
    other_edge = NULL;
    ref_edge = NULL;
    if ( is_face_blend(ref_face, ref_volume,
                       ref_edge, other_edge ) )
      blend_faces.append(ref_face);
  }
  if ( blend_faces.size() == 0 )
  {
    return;
  }
    //Find out how many different groups of surfaces there
    //are that share curves.
  DLIList <RefFace*> *blend_group;
  DLIList <RefFace*> stack;
  RefFace *other_face;
    //mark all the faces as 0.
  for ( ii =0; ii < ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);
    //mark just the blends as 1.
  for ( ii =0; ii < blend_faces.size(); ii++ )
    blend_faces.get_and_step()->marked(1);
  
  for ( ii =0; ii < blend_faces.size(); ii++ )
  {
    ref_face = blend_faces.get_and_step();
    if ( ref_face->marked() == 2 )
      continue;
    blend_group = new DLIList <RefFace*>;
    blend_groups.append(blend_group);
    stack.clean_out();
    stack.append(ref_face);
    while ( stack.size() > 0 )
    {
      ref_face = stack.pop();
      ref_face->marked(2);
      blend_group->append(ref_face);
      ref_edges.clean_out();
      ref_face->ref_edges(ref_edges);
      for ( jj = 0; jj < ref_edges.size(); jj++ )
      {
        ref_edge = ref_edges.get_and_step();
        other_face = ref_edge->other_face(ref_face, ref_volume);
        if ( other_face == NULL )
            //shouldn't happend.
          continue;
        if ( other_face->marked() == 1 )
          stack.append(other_face);
      }
    }
  }
    //cleanup marks.
  for ( ii =0; ii < ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);
  for ( ii =0; ii < blend_faces.size(); ii++ )
    blend_faces.get_and_step()->marked(0);
    //Okay we should have it now.

    //REMEMBER TO DELETE the dllists in blend_faces.
  return;
}
//--------------------------------------------------------------------
//Function: Public, is_face_blend
//Description: Determines if a face is a blend surface, returns the
//   ref_edge on one side of the blend and other_edge on the opposite side.
//   For this type of blend, only ref_edge must be tangentially meeting
//   with another surface.  Other_edge must be oriented orthogonally to
//   it and may or may not blend with another surface.  This assumes
//   a rectangular blend surface, without holes.
//---------------------------------------------------------------------
CubitBoolean GeomMeasureTool::is_face_blend(RefFace *ref_face,
                                            RefVolume *ref_volume,
                                            RefEdge *&ref_edge,
                                            RefEdge *&other_edge)
{
    //first we know that blend surfaces are not planar.
    //so remove these first.
    //Also, don't look at faces with more than 2 loops.
  if ( ref_face->geometry_type() == PLANE_SURFACE_TYPE ||
       ref_face->num_loops() > 2 )
    return CUBIT_FALSE;
  
  CubitBoolean is_cartesian;
  DLIList<RefEdge*> ref_edges;
  ref_edge = NULL;
  other_edge = NULL;
  RefFace *other_face;
  int jj;
  double angle;
  ref_face->ref_edges(ref_edges);
  for ( jj = ref_edges.size(); jj > 0; jj-- )
  {
    ref_edge = ref_edges.get_and_step();

    //Weed-out case where one edge is shared between more 
    //than 2 surfaces of the same volume
    DLIList<RefFace*> tmp_faces;
    ref_edge->ref_faces( tmp_faces );

    if( tmp_faces.size() > 2 )
    {
      int kk;
      for(kk=tmp_faces.size(); kk--;)
      {
        if( !tmp_faces.get()->is_child( ref_volume ) )
          tmp_faces.change_to(NULL);
        tmp_faces.step();
      }
      tmp_faces.remove_all_with_value( NULL );
      if( tmp_faces.size() > 2 )
        //this isn't the type of surface we are looking for...
        continue;
    }

    other_face = ref_edge->other_face(ref_face, ref_volume);
    if ( other_face == NULL )
    {
        //this isn't the type of surface we are looking for...
      break;
    }
    angle = GeometryQueryTool::instance()->surface_angle(ref_face,
                                                         other_face,
                                                         ref_edge,
                                                         ref_volume);
    angle *= 180.0/CUBIT_PI;
    is_cartesian = CUBIT_TRUE;
    if ( angle <= GEOM_SIDE_LOWER ||
         angle >= GEOM_SIDE_UPPER )
      is_cartesian = CUBIT_FALSE;
      //Okay, we have one major criteria achieved, this edge is a cartesian meet.

    if ( !is_cartesian )
      continue;
      //Now we need to check the radius of curvatures between these
      //two surfaces. I'm not totally sure about this but I think we
      //don't want the same radius of curvature.
    double k1_s1, k2_s1, k1_s2, k2_s2;
    CubitVector mid_point;
    ref_edge->mid_point(mid_point);
    ref_face->get_principal_curvatures( mid_point, k1_s1, k2_s1, ref_volume);
    other_face->get_principal_curvatures( mid_point, k1_s2, k2_s2, ref_volume);
    if (( is_equal(k1_s1, k1_s2) || is_equal(k1_s1, k2_s2) ) &&
        ( is_equal(k2_s1, k1_s2) || is_equal(k2_s1, k2_s2) ) )
        //try a different edge.
      continue;

      //Okay, this looks pretty good.
      //Now try and find an edge that is on the "opposite" side of
      // this surface. Sort of assume this is like a mapable surface.
      //check the oposite side.
    other_edge = NULL;
    CubitVector opp_point;
    if ( !find_opposite_edge( ref_edge, ref_face, other_edge, opp_point) ||
         other_edge == NULL )
        //if we can't find an opposite edge, we need to not
        //continue checking this surface.
      return CUBIT_FALSE;
    
      //Okay, we have the point opposite.  If this is a blend surface
      //the normal of this point should be pointing 90 degrees from
      //us.
    CubitVector normal_1 = ref_face->normal_at(mid_point, ref_volume);
    CubitVector normal_2 = ref_face->normal_at(opp_point, ref_volume);
    normal_1.normalize();
    normal_2.normalize();
    double dot_product = normal_1 % normal_2;
    if ( dot_product >= .5 || dot_product <= -.5  )
        //try a different edge.
      continue;
    CubitVector tangent_1, tangent_2;
    ref_edge->tangent(mid_point, tangent_1, ref_face);
    other_edge->tangent(opp_point, tangent_2, ref_face);
    tangent_1.normalize();      
    tangent_2.normalize();      
    dot_product = tangent_1 % tangent_2;
      //basically these edges should be somewhat parallel at
      //this point.
    if ( dot_product < .5 && dot_product > -.5 )
      continue;
//This was called but never checked, not sure why.
//      //check the curvature opposite.
//    double k1_s3, k2_s3;
//    ref_face->get_principal_curvatures( opp_point, k1_s3, k2_s3, ref_volume);
      //we are done with this face. It is a blend face.
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitBoolean GeomMeasureTool::find_opposite_edge( RefEdge* ref_edge,
                                                  RefFace* ref_face,
                                                  RefEdge *&other_edge,
                                                  CubitVector &closest_point)
{
    //Here we want to find the edge that is "opposite" to ref_edge.
    //Opposite is defined by assuming the surface is a square or cylinder and
    //finding the edge on the other side, so that there is a point closest
    //to the mid point on the ref_edge.
  //double need_to_turn = CUBIT_PI;
  int ii;
  DLIList <Loop*> loops;
  DLIList <CoEdge*> co_edges;
  Loop *curr_loop;
  CoEdge *co_edge, *this_co_edge;
  ref_edge->get_co_edges(co_edges, ref_face);
  if ( co_edges.size() != 1 )
    return CUBIT_FALSE;
  this_co_edge = co_edges.get();
  co_edges.clean_out();
  ref_face->loops(loops);
  if ( loops.size() > 2 )
    return CUBIT_FALSE;
  CubitBoolean found_loop = CUBIT_FALSE;
    //find the loop where this_co_edge is on.
  for(ii = loops.size(); ii > 0; ii-- )
  {
    curr_loop = loops.get_and_step();
    co_edges.clean_out();
    curr_loop->ordered_co_edges(co_edges);
    if ( co_edges.move_to(this_co_edge) )
    {
      found_loop = CUBIT_TRUE;
      break;
    }
  }
  if ( !found_loop )
    return CUBIT_FALSE;
  CoEdge *opp_co_edge = NULL;
  if ( loops.size() == 1 )
  {
      //This is the normal case, I hope.
      //this co_edge, go around the interior until
      //we turn twice.  If there are angles greater than CUBIT_PI, this isn't
      //a traditional blend.  Everything should be convex.
    if ( co_edges.size() == 4 )
    {
        //Just get the one opposite...
      opp_co_edge = co_edges.next(2);
    }
    else
    {
      //double total_angle = 0.0;
      int turn_counter = 0;
      CoEdge *next_co_edge = NULL;
      for ( ii = co_edges.size(); ii > 0; ii-- )
      {
        co_edge = co_edges.get_and_step();
        next_co_edge = co_edges.get();
        double angle = GeometryQueryTool::instance()->geometric_angle(co_edge, next_co_edge);
        angle *= 180.0/CUBIT_PI;
          //we don't deal with small or concave angles.
        if ( angle < 10.0 )
          return CUBIT_FALSE;
        else if ( angle > 185.0 )
          return CUBIT_FALSE;
        else if ( angle < 140.0 )
          turn_counter++;
        if ( turn_counter == 2 )
        {
          opp_co_edge = next_co_edge;
          break;
        }
      }
    }
  }
  else
  {
      //Just get the other loop and get one coedge in it.
    co_edges.clean_out();
    curr_loop = loops.get();
    curr_loop->ordered_co_edges(co_edges);
    opp_co_edge = co_edges.get();
  }
  if ( opp_co_edge == NULL )
    return CUBIT_FALSE;
    //This is just here to test correctness to this point, it can be taken out
    //if need be.
  if ( !co_edges.move_to(opp_co_edge) )
  {
    PRINT_ERROR("logic messed-up in find_opposite edge\n");
    assert(0);
    return CUBIT_FALSE;
  }
    //okay, we now have the co_edges in the correct position and
    //we have a canidate co_edge.
    //Search from the co_edge with the point that has closest point
    //to the mid point of this_co_edge.
  other_edge = NULL;
  CubitVector mid_point;
  ref_edge->mid_point(mid_point);
  CubitVector tmp_closest_point;
  double min_dist_sq = CUBIT_DBL_MAX;
  double dist_sq;
  RefEdge *test_ref_edge;
  CoEdge *prev_co_edge;
  CubitBoolean first = CUBIT_TRUE;
  for ( ii = co_edges.size(); ii > 0; ii-- )
  {
    co_edge = co_edges.get_and_step();
    prev_co_edge = co_edges.prev(2);
    test_ref_edge = co_edge->get_ref_edge_ptr();
      //We don't want edges that are connected to ref edge.
      //If we have them, we have gone to far around the surface.
    if ( test_ref_edge->common_ref_vertex(ref_edge) )
      break;
    if ( !first )
    {
        //break out of the for loop if we turn again. We just want the closest
        //point opposite.
      double angle = GeometryQueryTool::instance()->
        geometric_angle(prev_co_edge, co_edge);
      angle *= 180.0/CUBIT_PI;
      if ( angle > 195.0 || angle < 145.0 )
        break;
    }
    else
      first = CUBIT_FALSE;
      //now find the closest point to the mid point of ref_edge.
    test_ref_edge->closest_point_trimmed( mid_point,
                                          tmp_closest_point);
    dist_sq = (mid_point - closest_point).length_squared();
      //find and keep the closest one.
    if ( dist_sq < min_dist_sq )
    {
      other_edge = test_ref_edge;
      min_dist_sq = dist_sq;
      closest_point = tmp_closest_point;
    }
  }
  if ( other_edge == NULL )
    return CUBIT_FALSE;

  return CUBIT_TRUE;
}
CubitBoolean GeomMeasureTool::is_equal(double v1, double v2)
{
  double smallv = v1-v2;
  if ( smallv < 0.0 )
  {
    if ( smallv >= -CUBIT_RESABS )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }
  else
  {
    if ( smallv <= CUBIT_RESABS )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }
}
CubitStatus GeomMeasureTool::get_centroid( RefFace *ref_face, CubitVector &centroid, double &tot_area )
{
  int i;
  int num_tris, num_pnts, num_facets;
  GMem g_mem;
  unsigned short norm_tol = 5;
  double dist_tol = -1.0;

  ref_face->get_geometry_query_engine()->
      get_graphics(ref_face->get_surface_ptr(),num_tris, num_pnts, num_facets,
      &g_mem, norm_tol, dist_tol );

  if(num_tris < 1)
  {
      // Decrease tolerance and try again (we can get this for small features)
      norm_tol /= 2;
      ref_face->get_geometry_query_engine()->
          get_graphics(ref_face->get_surface_ptr(),num_tris, num_pnts, num_facets,
          &g_mem, norm_tol, dist_tol );
  }

  if(num_tris < 1)
  {
      // Lets give up 
      PRINT_ERROR( "Unable to find the center of a surface\n" );
      return CUBIT_FAILURE;
  }

  // Initialize
  tot_area = 0.0;
  centroid.set( 0.0, 0.0, 0.0 );

  // Loop through the triangles
  double tri_area;
  GPoint p[3];
  GPoint* plist = g_mem.point_list();
  int* facet_list = g_mem.facet_list();
  int c = 0;
  for( i=0; i<num_tris; i++ )
  {
      p[0] = plist[facet_list[++c]];
      p[2] = plist[facet_list[++c]];
      p[1] = plist[facet_list[++c]]; 
      c++;

      // Get centroid
      CubitVector p1( p[0].x, p[0].y, p[0].z );
      CubitVector p2( p[2].x, p[2].y, p[2].z );
      CubitVector p3( p[1].x, p[1].y, p[1].z );

      CubitVector center = (p1 + p2 + p3)/3.0;

      //Calculate area
      CubitVector vec1 = (p2 - p1);
	  CubitVector vec2 = (p3 - p1);

	  CubitVector cross = vec1 * vec2;
		
	  tri_area = 0.5* sqrt(cross.x() * cross.x() + cross.y() * cross.y() + cross.z() * cross.z());
		
		
	  centroid += (tri_area * center);

	  tot_area += tri_area;
        
    }
  centroid /= tot_area;
  return CUBIT_SUCCESS;
}    

CubitStatus
GeomMeasureTool::center( DLIList<RefFace*> ref_faces )
{
  int ii,id;
  double surf_area;
  double tot_area = 0.0;
  CubitVector surf_centroid;
  CubitVector centroid;
  centroid.set( 0.0, 0.0, 0.0 );

  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
	RefFace *ref_face = ref_faces.get_and_step();
	if( GeomMeasureTool::get_centroid( ref_face, surf_centroid, surf_area) != CUBIT_SUCCESS )
	  return CUBIT_FAILURE;
    centroid += (surf_area * surf_centroid);
	tot_area += surf_area;
	id = ref_face->id();
  }
  
  centroid /= tot_area;

  if( ref_faces.size() > 1 )
	PRINT_INFO("Centroid of Composite surface located at: ( %f , %f , %f )\n\n",
                 centroid.x(), centroid.y(), centroid.z());
  else
	PRINT_INFO("Centroid of Surface %i located at: ( %f , %f , %f )\n\n",
                  id, centroid.x(), centroid.y(), centroid.z());
  return CUBIT_SUCCESS;
}
