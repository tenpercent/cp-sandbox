//-----------------------------------------------------------
// Class: MedialTool2D
// Description:  Creates medials for 2D surfaces.  Uses the CMU
//               code for creating the medial.  This class
//               serves as the interface to that.
// Creator: David white
// Date: 6/21/2002
//-----------------------------------------------------------
#include "MedialTool2D.hpp"
#include "CubitMessage.hpp"
#ifdef USING_MEDIAL
#include "GfxDebug.hpp"
#include "CubitVector.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "CoEdge.hpp"
#include "Curve.hpp"
#include "RefVertex.hpp"
#include "GeometryQueryEngine.hpp"
#include "VirtualGeometryEngine.hpp"
//----------medial includes---------------------
#include "medial/medial_util/MedialDefines.hpp"
#include "medial/medial2D/Medial2DAPI.hpp"
#include "medial/medial_util/Vec3D.hpp"
#include "medial/medial2D/MedialVertex.hpp"
#include "medial/medial2D/MedialSegment.hpp"
//----------medial includes---------------------
#endif //USING_MEDIAL


MedialTool2D::MedialTool2D(RefFace *ref_face)
{
  myRefFace = ref_face;
}
MedialTool2D::~MedialTool2D()
{
    //do nothing for now.
}
#ifndef USING_MEDIAL
CubitStatus MedialTool2D::create_medial_axis(DLIList <RefEdge*>&, double, double)
{
  PRINT_ERROR("Medial Axis Code is not included with this version of CGM.\n");
  return CUBIT_FAILURE;
}
#endif //USING_MEDIAL
//---------------------------------------------------------------------------
//Functions that get compiled only when compiling in the medial stuff
//---------------------------------------------------------------------------
#ifdef USING_MEDIAL
//////////////////////////////////////////////////////////////////////
///Function Name : create_medial_axis
///Member Type   : Public
///Description   : create 2D medial axis for myRefFace. Does the following
///                1. get_boundary_points() as  DLIList <PointList*>
///                2. convert_to_segments() to get
///                      vector <vector<MedialSegment*>*>
///                3. create_medial_axis() and get results as
///                    vector<MedialVertex*> and vector<MedialSegment*>
///                4. 
///
///Date          : 09/11/2002
//////////////////////////////////////////////////////////////////////
CubitStatus MedialTool2D::create_medial_axis( DLIList <RefEdge*> &medial_axis,
                                              double cell_size,
                                              double angle_tol)
{
    //first get the loop of edges.
  PointLoopList boundary_point_loops;
  CubitStatus stat = get_boundary_points(myRefFace,
                                         boundary_point_loops);

    //Now convert these to segments.
  MedialSegmentPLoop boundary_segments;
  stat = convert_to_segments(boundary_point_loops, boundary_segments);
  
  MedialSegmentPArray output_segments;
  MedialVertexPArray output_points;
  int results = Medial2DAPI::create_medial(boundary_segments, output_points,
                                           output_segments, cell_size,
                                           angle_tol);
  if ( results != 1 )
  {
    PRINT_ERROR("Medial Extraction Failed.\n");
    return CUBIT_FAILURE;
  }

    // create geometry from the medial results
//  create_medial_geometry( output_points, output_segments, medial_axis );
  draw_results(output_points, output_segments);

  free_medial_memory(output_points, output_segments);

  return CUBIT_SUCCESS;
}

/// draw the results to the graphics output screen
void MedialTool2D::draw_results( MedialVertexPArray &output_points,
                                 MedialSegmentPArray &output_segments )
{
  int ii;
  CubitVector vec, vec_2;
  Vec3D point, start_point, end_point;
  MedialVertex *vertex, *start, *end;
  MedialSegment *seg;
  for ( ii = 0; ii < output_points.size(); ii++ )
  {
    vertex = output_points[ii];
    point = vertex->coordinates();
    convert_from_vec3d(point,vec);
    GfxDebug::draw_point(&vec,CUBIT_RED);
    GfxDebug::flush();
  }
  GfxDebug::flush();

  for ( ii = 0; ii < output_segments.size(); ii++ )
  {
    seg = output_segments[ii];
    start = seg->get_start();
    end = seg->get_end();
    start_point = start->coordinates();
    end_point = end->coordinates();
    GfxDebug::draw_line( start_point.x, start_point.y,
                         start_point.z, end_point.x,
                         end_point.y, end_point.z,
                         CUBIT_CYAN);
    GfxDebug::flush();
  }
  GfxDebug::flush();
}

//////////////////////////////////////////////////////////////////////
///Function Name : free_medial_memory
///Member Type   : Private
///Description   : free the medial memory
///Author        : Ragunath Sankaranarayanan
///Date          : 11/23/2002
//////////////////////////////////////////////////////////////////////
void MedialTool2D::free_medial_memory(MedialVertexPArray &vertices,
                                      MedialSegmentPArray &segments)
{
  int ii;
  MedialSegment *currs;
  MedialVertex *currv;
  
  for(ii=0; ii < segments.size(); ii++)
  {
    currs = segments[ii];
    delete currs;
  }
  
  for(ii=0; ii < vertices.size(); ii++)
  {
    currv = vertices[ii];
    delete currv;
  }
}


                                 
CubitStatus MedialTool2D::convert_to_segments(
  PointLoopList &boundary_point_loops,
  MedialSegmentPLoop &boundary_segments)
{
  int ii;
  PointList *point_loop;
  int jj;
  MedialVertex *start_vert, *end_vert;
  MedialSegment *new_seg;
  MedialSegmentPArray *new_vector;
  for ( ii = 0; ii < boundary_point_loops.size(); ii++ )
  {
    point_loop = boundary_point_loops.get_and_step();
    new_vector = new MedialSegmentPArray;
    for ( jj = 0; jj < point_loop->size(); jj++ )
    {
      start_vert = point_loop->get_and_step();
      end_vert = point_loop->get();
      new_seg = new MedialSegment(start_vert, end_vert);
      new_vector->push_back(new_seg);
    }
    boundary_segments.push_back(new_vector);
  }
  return CUBIT_SUCCESS;
}
 
                                               
CubitStatus MedialTool2D::get_boundary_points( RefFace *ref_face,
                                               PointLoopList &boundary_point_loops )
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
  RefVertex *temp_vert;
  double tolerance = 500*GEOMETRY_RESABS;
  
  
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list_ptr = co_edge_loops.get_and_step();
    new_point_loop_ptr = new PointList;
    for ( jj = co_edge_list_ptr->size(); jj > 0; jj-- )
    {
      co_edge_ptr = co_edge_list_ptr->get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      tmp_point_list.clean_out();
      stat = get_curve_facets( ref_edge_ptr, tmp_point_list );
      PRINT_DEBUG_129("curve %d has %d points\n",
                      ref_edge_ptr->id(),
                      tmp_point_list.size());
      if ( stat != CUBIT_SUCCESS )
        return CUBIT_FAILURE;
      tmp_point_list.reset();
        //the points are in order from start vertex to end vertex.
        //append them now according to the loop.

        //Assign the points to be part owned by the vertex, rather than the
        //curve.
      if ( ref_edge_ptr->start_vertex() !=
           ref_edge_ptr->end_vertex() )
      {
        MedialVertex *temp_point = tmp_point_list.get();
        Vec3D temp_v1 = temp_point->coordinates();
        CubitVector v1;
        convert_from_vec3d(temp_v1, v1);
        CubitVector v2 = ref_edge_ptr->start_vertex()->coordinates();
        if ( !v1.within_tolerance(v2, tolerance) )
        {
          PRINT_ERROR("Problem with surface geometry\n"
                      "Check surface %d, especially curve %d\n",
                      myRefFace->id(), ref_edge_ptr->id());
          return CUBIT_FAILURE;
        }
        temp_vert = ref_edge_ptr->start_vertex();
        temp_point = tmp_point_list.prev();
        temp_v1 = temp_point->coordinates();
        convert_from_vec3d(temp_v1, v1);
        v2 = ref_edge_ptr->end_vertex()->coordinates();    
        if (!v1.within_tolerance(v2, tolerance))
        {
          PRINT_ERROR("Problem with surface geometry\n"
                      "Check surface %d and %d, especially curve %d\n",
                      myRefFace->id(), ref_edge_ptr->id());
          return CUBIT_FAILURE;
        }
        temp_vert = ref_edge_ptr->end_vertex();
      }
      else
      {

      }
      tmp_point_list.reset();
      sense = co_edge_ptr->get_sense();
      if ( CUBIT_FORWARD != sense )
        tmp_point_list.reverse();
        //Now take off the last point as it is a duplicate with the
        //other list...
      tmp_point_list.reset();
      if ( co_edge_list_ptr->size() != 1 )
        tmp_point_list.pop();
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

CubitStatus MedialTool2D::get_curve_facets( RefEdge* curve, PointList &segments ) 
{
  const double COS_ANGLE_TOL =  0.996194698091745545198705; // cos(5)
  GMem curve_graphics;
  double tolerance = 500*GEOMETRY_RESABS;

  const double dist_tol = tolerance;
  const double dist_tol_sqr = dist_tol*dist_tol;
  int n;
  curve->get_curve_ptr()->get_geometry_query_engine()->get_graphics_edges( 
    curve, n, &curve_graphics );
  
  GPoint* gp = curve_graphics.point_list();
  CubitVector lastv(gp[0].x, gp[0].y, gp[0].z);
  Vec3D last_v3d;
  convert_to_vec3d(lastv, last_v3d);
  MedialVertex* last = new MedialVertex(last_v3d);
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
    Vec3D point;
    convert_to_vec3d(pos,point);
    last = new MedialVertex(point);
    segments.append( last );
    lastv = pos;
  }
    // Now check if the segment list is reversed wrt the curve direction.
  segments.reset();
  if ( remove_second_to_end )
  {
    if ( segments.size() == 2 )
    {
    }
    else
    {
        //Remove the second to last one.  To do
        //this efficiently (don't do remove), pop
        //the last one, then save that and
        //re-add it after poping the second one.
      MedialVertex *temp = segments.pop();
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
    Vec3D temp_pos = segments.get()->coordinates();
    CubitVector start_seg;
    convert_from_vec3d(temp_pos, start_seg);
    double dist_1 = (start_seg - start_vec).length_squared();
    double dist_2 = (start_seg - end_vec).length_squared();
    if ( dist_1 > dist_2 )
      segments.reverse();
  }
  else
  {
    double u1, u2;
    CubitVector next_1, next_2;
    Vec3D temp_next_1, temp_next_2;
    temp_next_1 = segments.next(1)->coordinates();
    temp_next_2 = segments.next(2)->coordinates();
    convert_from_vec3d(temp_next_1, next_1);
    convert_from_vec3d(temp_next_2, next_2);
    u1 = curve->u_from_position( next_1 );
    u2 = curve->u_from_position( next_2 );    
    if( (u2 < u1) && (curve->start_param() <= curve->end_param()) )
      segments.reverse();
  }
    //clean up the periodic curve case (last seg may be too small.)
  if ( curve->start_vertex() == curve->end_vertex() )
  {
    segments.reset();
    CubitVector start_v, last_v;
    Vec3D temp_start_v, temp_last_v;
    temp_start_v = segments.get()->coordinates();
    temp_last_v = segments.prev()->coordinates();
    convert_from_vec3d(temp_start_v, start_v);
    convert_from_vec3d(temp_last_v, last_v);
    double dist = (start_v - last_v).length();
    if ( dist < dist_tol )
    {
        //remove the last one.
      segments.pop();
    }
  }
  if ( segments.size() < 2 )
  {
    PRINT_ERROR("Problem with getting boundary segments.\n");
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
void MedialTool2D::convert_from_vec3d(Vec3D &vec, CubitVector &new_vec)
{
  new_vec.set(vec.x, vec.y, vec.z);
}
void MedialTool2D::convert_to_vec3d(CubitVector& vec, Vec3D &new_vec)
{
  new_vec.x = vec.x();
  new_vec.y = vec.y();
  new_vec.z = vec.z();
}

//////////////////////////////////////////////////////////////////////
///Function Name : create_medial_geometry
///Member Type   : Private
///Description   : generates geometry from the medial output.
///Creator       : Shiraj Khan                
///Date          : 11/11/2002
//////////////////////////////////////////////////////////////////////
CubitStatus MedialTool2D::create_medial_geometry(MedialVertexPArray &output_points,
                                   MedialSegmentPArray &output_segments,
								   DLIList<RefEdge*> &medial_axis)
{
	int ii;
  int jj;
  int iii;
  int jjj;
  int r=0;
  MedialSegment *seg1, *seg2, *curr_seg;
  MedialVertex  *start1, *start2, *end1, *end2;
  Vec3D start_point1, start_point2, end_point1, end_point2;
  for( ii=0; ii < output_segments.size(); ii++)//finding a segment attached to only one segment
  {
	  seg1 = output_segments[ii];
	  start1 = seg1->get_start();
	  end1 = seg1->get_end();
	  start_point1 = start1->coordinates();
	  end_point1 = end1->coordinates();
	  for ( jj = 0; jj < output_segments.size(); jj++)
	  {
		  seg2 = output_segments[jj];
		  start2 = seg2->get_start();
		  end2 = seg2->get_end();
		  start_point2 = start2->coordinates();
		  end_point2 = end2->coordinates();
		  if ( start_point1 == start_point2 || start_point1 == end_point2 || end_point1 == start_point2 || end_point1 == end_point2)
			  r=r+1;
	  }
	  if( r == 2)
		  break;
	  else 
		  r=0;
  }
   
  curr_seg = output_segments[ii];
  DLIList<DLIList<MedialSegment*>*> segment_listoflist;
  DLIList<MedialSegment*> *seg_list;
  
  r = 0;// 
  int w = 0;
  int b = 0;
  int y = 0;
  int x = 0; //no. of points at which more than one medial are emerging
  int a[50];//tracks number of medials attached to a particular point ( used as a[x] )
  int v[50];//tracks segments attached to a particular segment in the output_segments list 
  int u[50];//no. of medials emerging from a point
  int z[50][50];//first subscript represents no. of diverging points ; second represents no.of medial attached to a particular diverging point
  u[0] = -1;
  int flag = 0;
  
   seg_list  = new DLIList<MedialSegment*>;
   seg_list->append(curr_seg);
   output_segments[ii] = NULL;
   seg1 = curr_seg;
              
  for( ii=0; ii < output_segments.size()-1; ii++)
  {
      
	  start1 = seg1->get_start();
	  end1 = seg1->get_end();
	  start_point1 = start1->coordinates();
	  end_point1 = end1->coordinates();
	  
	  for ( jj = 0; jj < output_segments.size(); jj++)
	  {
		  if( output_segments[jj] == NULL)
			  continue;
		  else
		  {
			  for ( jjj = 0; jjj <= b; jjj++)
			  {
				  if ( jj == u[jjj] )
				  {
					  flag = 1;
				      break;
				  }
			  }
			  if ( flag )
				  
			  {
				  flag = 0;
			      continue;
			  }
			  else
			  {
				  seg2 = output_segments[jj];
			      start2 = seg2->get_start();
			      end2 = seg2->get_end();
			      start_point2 = start2->coordinates();
			      end_point2 = end2->coordinates();
			  }
		  }
			  
			  
		       
		  if ( start_point1 == start_point2 || start_point1 == end_point2 || end_point1 == start_point2 || end_point1 == end_point2)
		  {
			  w = r;
			  v[w] = jj;
			  curr_seg = output_segments[jj];
			  r = r+1;
			  
		  }
	  }//endof for jj

	  if( r == 1)// only one segment attached to a given segment
	  {
		  seg_list->append(curr_seg);
		  seg1 = curr_seg;
		  output_segments[v[w]] = NULL;
		  curr_seg = NULL;
		  r = 0;
		  if ( ii == output_segments.size()-2 )
		  {
			  segment_listoflist.append(seg_list);
		  }

	  }// end of if ( r == 1 )
	  
	  else
	  {
		  if( r == 0)// no segment attached to a given segment
		  {
			  segment_listoflist.append(seg_list);
			  if ( x == 0 && y == 0 )//getting a segment with only one segment attached to it
			  {
				  for( ii=0; ii < output_segments.size(); ii++)
				  {
					  if ( output_segments [ii] == NULL)
							  continue;
					  else
					  {
						  seg1 = output_segments[ii];
	                      start1 = seg1->get_start();
	                      end1 = seg1->get_end();
	                      start_point1 = start1->coordinates();
	                      end_point1 = end1->coordinates();
					  }
					  for ( jj = 0; jj < output_segments.size(); jj++)
					  {
						  if ( output_segments[jj] == NULL )
							  continue;
						  else
						  {
							  seg2 = output_segments[jj];
		                      start2 = seg2->get_start();
		                      end2 = seg2->get_end();
		                      start_point2 = start2->coordinates();
		                      end_point2 = end2->coordinates();
						  }// end of if ( output_segments[jj] == NULL )
						  if ( start_point1 == start_point2 || start_point1 == end_point2 || end_point1 == start_point2 || end_point1 == end_point2)
								  r=r+1;
					  }// end of for ( jj = 0; jj < output_segments.size(); jj++)
					  if( r == 2)
						  break;
	                  else 
						  r=0;
				  }// end of for( ii=0; ii < output_segments.size(); ii++)
				  seg1 = output_segments[ii];
				  output_segments[ii] = NULL;
				  seg_list  = new DLIList<MedialSegment*>;
				  seg_list->append(seg1);
		          r = 0;
					  
			  }
			  else
			  {
				  while( output_segments[z[x][y]] == NULL)
				  {
					  y = y-1;
			          if( y < 0)
					  {
						  x = x-1;
					      if( x < 1 )
							  break;
					       else
							   y = a[x];
					  }// end of if ( y < 0 )
			        
				  } //end of while( output_segments[z[x][y]] == NULL)
				  if ( x <= 0 && y <= -1 ) 
				  {
					  for( ii=0; ii < output_segments.size(); ii++)
					  {
						  if ( output_segments [ii] == NULL)
							  continue;
						  else
						  {
							  seg1 = output_segments[ii];
	                          start1 = seg1->get_start();
	                          end1 = seg1->get_end();
	                          start_point1 = start1->coordinates();
	                          end_point1 = end1->coordinates();
						  }// end of if ( output_segments [ii] == NULL)
						  for ( jj = 0; jj < output_segments.size(); jj++)
						  {
							  if ( output_segments[jj] == NULL )
								  continue;
							  else
							  {
								  seg2 = output_segments[jj];
		                          start2 = seg2->get_start();
		                          end2 = seg2->get_end();
		                          start_point2 = start2->coordinates();
		                          end_point2 = end2->coordinates();
							  }// end of else to  if ( output_segments[jj] == NULL )
							  if ( start_point1 == start_point2 || start_point1 == end_point2 || end_point1 == start_point2 || end_point1 == end_point2)
								        r=r+1;
						  }// end of for ( jj = 0; jj < output_segments.size(); jj++)
						  if( r == 2)
							  break;
	                      else
							  r=0;
					  }// end of for( ii=0; ii < output_segments.size(); ii++)
					  seg1 = output_segments[ii];
					  output_segments[ii] = NULL;
					  seg_list  = new DLIList<MedialSegment*>;
				      seg_list->append(seg1);
		              r = 0;
					  x = 0;
					  y = 0;
				  }
				  else
				  {
					  seg1 = output_segments[z[x][y]];
		              output_segments[z[x][y]] = NULL;
		              seg_list  = new DLIList<MedialSegment*>;
				      seg_list->append(seg1);
		              r = 0;
				          
				  }
			  }
			  
		  }// end of if ( r == 0 )
	      
		  else
		  {
			  if ( r>1 )//more than one segments attached to a given segment
			  {
				  segment_listoflist.append(seg_list);
		          x = x+1;
		          a[x] = w;
		          y = w;
				  b = w;
		          for( iii = 0; iii < r; iii++)
				  {
					  z[x][iii] = v[iii];
					  u[iii] = v[iii];
				  }
				  seg1 = output_segments[z[x][y]];
		          output_segments[z[x][y]] = NULL;
		          seg_list  = new DLIList<MedialSegment*>;
			      seg_list->append(seg1);
		          r = 0;
			  }// end of if ( r > 1)
		  }// end of else to if ( r == 0 )
	  }// end of else to if ( r == 1)
  }// end of for( ii=0; ii < output_segments.size()-1; ii++)  
	  
  segment_listoflist.reset();
  MedialSegment *first_seg, *second_seg;
  MedialVertex  *start, *end, *intermediate_start, *intermediate_end, *secondseg_start, *secondseg_end ;
  Vec3D start_point, end_point, intermediate_prev_end_point, intermediate_curr_start_point, intermediate_curr_end_point, secondseg_startpoint, secondseg_endpoint;
  DLIList<CubitVector*> curve_vectors;
  CubitVector cubit_start_vertex, cubit_end_vertex, cubit_intermediate_vertex;
  
  for ( ii = 0; ii < segment_listoflist.size(); ii++)
  {
	  seg_list  = new DLIList<MedialSegment*>;
	  seg_list = segment_listoflist.get_and_step();
	  curve_vectors.clean_out();
	  for ( jj = 0; jj < seg_list->size(); jj++)
	  {
		 first_seg = seg_list->next(jj);
		  
		  if ( seg_list->size() == 1 )//seg_list contains only one segment
		  {
			  start =first_seg->get_start();
			  end =first_seg->get_end();
			  start_point = start->coordinates();
			  end_point = end->coordinates();
			  convert_from_vec3d(start_point,cubit_start_vertex);
			  convert_from_vec3d(end_point,cubit_end_vertex);
		  }// end of if ( seg_list->size() == 1 )
		  else
		  {
			  if ( jj == 0)//for the first segment in the seg_list list
			  {
				  start =first_seg->get_start();
				  end =first_seg->get_end();
				  start_point = start->coordinates();
				  end_point = end->coordinates();
				  second_seg = seg_list->next(jj+1);
				  secondseg_start = second_seg->get_start();
				  secondseg_end = second_seg->get_end();
				  secondseg_startpoint = secondseg_start->coordinates();
				  secondseg_endpoint = secondseg_end->coordinates();
                  
				  if ( end_point == secondseg_startpoint || end_point == secondseg_endpoint )
				  {
					  convert_from_vec3d(start_point,cubit_start_vertex);
					  intermediate_prev_end_point = end_point;
				  }// end of if ( end_point == secondseg_startpoint || end_point == secondseg_endpoint )
				  else if ( start_point == secondseg_startpoint || start_point == secondseg_endpoint )
				  {
					  convert_from_vec3d(end_point,cubit_start_vertex);
					  intermediate_prev_end_point = start_point;
					  
				  }//end of if ( start_point == secondseg_startpoint || start_point == secondseg_endpoint )
				  
			  }// end of if ( jj== 0 )
			  else
			  {
				  intermediate_start =first_seg->get_start();
				  intermediate_end =first_seg->get_end();
				  intermediate_curr_start_point = intermediate_start->coordinates();
				  intermediate_curr_end_point = intermediate_end->coordinates();
				  if ( intermediate_curr_start_point == intermediate_prev_end_point )
				  {
					  convert_from_vec3d(intermediate_curr_start_point,cubit_intermediate_vertex);
					  intermediate_prev_end_point = intermediate_curr_end_point;
					  
				  }// end of if ( intermediate_curr_start_point == intermediate_prev_end_point )
				  else if ( intermediate_curr_end_point == intermediate_prev_end_point )
				  {
					  convert_from_vec3d(intermediate_curr_end_point,cubit_intermediate_vertex);
					  intermediate_prev_end_point = intermediate_curr_start_point;
					  
				  }// end of if ( intermediate_curr_end_point == intermediate_prev_end_point )
				  
				  curve_vectors.append(&cubit_intermediate_vertex);
				  
				  if ( jj == seg_list->size()-1 )//for the last segment in the seg_list
				  {
					  convert_from_vec3d(intermediate_prev_end_point,cubit_end_vertex);
					   
				  }// end of if ( jj == seg_list->size()-1 )

			  }//end of else to if ( jj == 0 )
		  }//end of else to if ( seg_list->size() == 1 ) 
	  }// end of for ( jj = 0; jj < seg_list->size(); jj++)
	  
	  RefVertex *first_ref = VirtualGeometryEngine::instance()->create_VirtualVertex(cubit_start_vertex);
	  RefVertex *last_ref = VirtualGeometryEngine::instance()->create_VirtualVertex(cubit_end_vertex);;
	  assert( first_ref && last_ref );
	  RefEdge *ref_edge = NULL;
	  if ( curve_vectors.size() == 0)//only one segment 
	  {
		  
		  ref_edge = VirtualGeometryEngine::instance()->create_VirtualEdge(first_ref, last_ref);
		  medial_axis.append( ref_edge );
		 
	  }//end of if ( curve_vectors.size() == 0)
	  else
	  {
		  
		  ref_edge = VirtualGeometryEngine::instance()->create_VirtualEdge(first_ref, last_ref, curve_vectors);
		  medial_axis.append( ref_edge );
		  
		  
	  }// end of if ( curve_vectors.size() == 0)
	  
	  seg_list->clean_out();
	  
  }//end of for ( ii = 0; ii < segment_listoflist.size(); ii++)


  return CUBIT_SUCCESS;
}

#endif //USING_MEDIAL
