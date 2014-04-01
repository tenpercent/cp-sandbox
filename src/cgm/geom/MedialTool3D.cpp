//////////////////////////////////////////////////////////////////////
/// Filename      : MedialTool3D.cpp (implementation of MedialTool3D)
///
/// Purpose       : 3-D Medial Axis extraction tool. Given a RefVolume
///                 generates the medial axis for that volume
///
/// Special Notes :
///
/// Creator       : Ragunath Sankaranarayanan
///
/// Creation Date : 08/23/2002
///
/// Owner         : Ragunath Sankaranarayanan
//////////////////////////////////////////////////////////////////////

#include "MedialTool3D.hpp"

#ifdef USING_MEDIAL
//////////////////////////////////////////////////////////////////////
/// Medial3D Includes
//////////////////////////////////////////////////////////////////////
#include "medial/medial_util/MedialDefines.hpp"
#include "medial/medial_util/Triangle3D.hpp"
#include "medial/medial_util/Vec3D.hpp"
#include "medial/medial3D/Medial3DAPI.hpp"
#include "medial/medial3D/MedialVertex3D.hpp"
#include "medial/medial3D/MedialEdge3D.hpp"
#include "medial/medial3D/MedialFace3D.hpp"

//////////////////////////////////////////////////////////////////////
/// Cubit Includes
//////////////////////////////////////////////////////////////////////
#include "GfxDebug.hpp"
#include "CubitVector.hpp"
#include "RefVertex.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#endif //USING_MEDIAL


//////////////////////////////////////////////////////////////////////
///Function Name : MedialTool3D constructor
///Member Type   : Public
///Description   : consturctor when a RefVolume ptr is given
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
MedialTool3D::MedialTool3D(RefVolume *ref_volume_ptr)
{
  myRefVolume = ref_volume_ptr;
}

//////////////////////////////////////////////////////////////////////
///Function Name : MedialTool3D destructor
///Member Type   : Public
///Description   : default constructor
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
MedialTool3D::~MedialTool3D()
{

}

//////////////////////////////////////////////////////////////////////
///Function Name : create_medial_axis_3d (when USING_MEDIAL is not defined)
///Member Type   : Public
///Description   : create the medial axis of the given RefVolume, based
///                the input parameters; Return the medial axis as a
///                list of vertices, segments or faces or a combination
///                of all the three; Refer "Medial3D.hpp" for more info
///                on the parameters;
///Date          : 08/27/2002
//////////////////////////////////////////////////////////////////////
#ifndef USING_MEDIAL
CubitStatus MedialTool3D::create_medial_axis_3d( 
                  DLIList <RefVertex*> &medial_verts,
                  DLIList <RefEdge*> &medial_edges, 
                  DLIList <RefFace*> &medial_faces,
				  bool full_search,
                  double cube_size,
                  double accuracy,
                  double min_ang_strong,
                  double min_angle_diff_components)
{
  PRINT_ERROR("Medial Axis Code is not included with this version of CGM.\n");
  return CUBIT_FAILURE;
}
#endif

//////////////////////////////////////////////////////////////////////
// Functions that gets compiled only when compiling the medial stuff
//////////////////////////////////////////////////////////////////////
#ifdef USING_MEDIAL
//////////////////////////////////////////////////////////////////////
///Function Name : create_medial_axis_3d (when USING_MEDIAL is defined)
///Member Type   : Public
///Description   : create the medial axis of the given RefVolume, based
///                the input parameters; Return the medial axis as a
///                list of vertices, segments or faces or a combination
///                of all the three; Refer "Medial3D.hpp" for more info
///                on the parameters;
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
CubitStatus MedialTool3D::create_medial_axis_3d(
  DLIList <RefVertex*> &medial_verts,
  DLIList <RefEdge*> &medial_edges,
  DLIList <RefFace*> &medial_faces,
  bool full_search,
  double cube_size,
  double accuracy,
  double min_ang_strong,
  double min_angle_diff_components)
{
  // create the array of vertices and triangles to represent objectdata3d
  Vec3DPArray my_vertices;
  Triangle3DPArray my_facets;

  MedialVertex3DPArray my_output_verts;
  MedialEdge3DPArray my_output_edges;
  MedialFace3DPArray my_output_faces;

  // get the boundary facets
  CubitStatus boundaries = get_boundary_facets(my_vertices, my_facets);

  // create the medial axis by calling the Medial3DAPI
  int results = Medial3DAPI::create_medial3d(my_vertices, my_facets,
                   full_search, cube_size, accuracy, min_ang_strong,
                   min_angle_diff_components, my_output_verts,
                   my_output_edges, my_output_faces);

  if ( results != 1 )
  {
    PRINT_ERROR("Medial Extraction Failed.\n");
    return CUBIT_FAILURE;
  }
  draw_medial_results( my_output_verts, my_output_edges, my_output_faces,
              CUBIT_RED, CUBIT_CYAN, CUBIT_RED);

  // translate the data from Medial Extractor format back to Cubit Format
//?***** Not implemented as of now.
	return CUBIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////
///Function Name : convert_vec3d_to_cubitvector
///Member Type   : Private
///Description   : convert a point from Vec3D format to CubitVector
///                format.
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
void MedialTool3D::convert_vec3d_to_cubitvector( Vec3D &my_vec3d,
                                                 CubitVector &my_cubitvec)
{
  my_cubitvec.set(my_vec3d.x, my_vec3d.y, my_vec3d.z);
}

//////////////////////////////////////////////////////////////////////
///Function Name : convert_cubitvector_to_vec3d
///Member Type   : Private
///Description   : convert a point from CubitVector format to Vec3D
///                format.
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
void MedialTool3D::convert_cubitvector_to_vec3d( CubitVector &my_cubitvec,
                                                 Vec3D &my_vec3d)
{
  my_vec3d.x = my_cubitvec.x();
  my_vec3d.y = my_cubitvec.y();
  my_vec3d.z = my_cubitvec.z();
}

//////////////////////////////////////////////////////////////////////
///Function Name : draw_medial_results
///Member Type   : Private
///Description   : Draw the extracted medial to the screen
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
void MedialTool3D::draw_medial_results(
  MedialVertex3DPArray &output_vertices,
  MedialEdge3DPArray &output_edges,
  MedialFace3DPArray &output_faces,
  int vert_color, int edge_color, int face_color)
{
  int ii;
  Vec3D curr_vec3d_1, curr_vec3d_2;

  int n1 = output_vertices.size();
  int n2 = output_edges.size();

  printf("Before drawing the results to screen\n");
  printf("Number of output vertices = %d\n",output_vertices.size());
  printf("Number of output segments = %d\n",output_edges.size());
  printf("Number of output triangles= %d\n",output_faces.size());
/*  
  // draw the vertices
  for(ii =0; ii < output_vertices.size(); ii++)
    {
      curr_vec3d_1 = output_vertices[ii]->coordinates();
      GfxDebug::draw_point(curr_vec3d_1.x, curr_vec3d_1.y, curr_vec3d_1.z, vert_color);
    }
    GfxDebug::flush();
    
  // draw the edges
  for(ii=0; ii < output_edges.size(); ii++)
  {
    // get the start, end vertices and convert them to CubitVector format
    curr_vec3d_1 = output_edges[ii]->get_start_coordinates();
    curr_vec3d_2 = output_edges[ii]->get_end_coordinates();
    GfxDebug::draw_line(curr_vec3d_1.x, curr_vec3d_1.y,
        curr_vec3d_1.z, curr_vec3d_2.x, curr_vec3d_2.y, curr_vec3d_2.z,
        edge_color);
  }
  GfxDebug::flush();
*/
  //draw faces
  for(ii=0; ii<output_faces.size(); ii++)
  {
    // Points are ordered clockwise around the tri.
    draw_medial_face(output_faces[ii], face_color);
//  void draw_tri   (GPoint p[], int color);
  }
  GfxDebug::flush();
}

//////////////////////////////////////////////////////////////////////
///Function Name : draw_medial_face
///Member Type   : Private
///Description   : Given a MedialFace3D, draw it to the cubit graphics 
///                Output window
///Date          : 10/29/2002
//////////////////////////////////////////////////////////////////////
void MedialTool3D::draw_medial_face(MedialFace3D* my_face, int tri_color)
{
  // array of three points
  GPoint my_points[3];
  Vec3D v1, v2, v3;
  v1 = my_face->get_vertex1_coordinates();
  v2 = my_face->get_vertex2_coordinates();
  v3 = my_face->get_vertex3_coordinates();
  my_points[0].x = v1.x; my_points[0].y = v1.y; my_points[0].z = v1.z;
  my_points[1].x = v2.x; my_points[1].y = v2.y; my_points[1].z = v2.z;
  my_points[2].x = v3.x; my_points[2].y = v3.y; my_points[2].z = v3.z;
    //GfxDebug::draw_tri(my_points, tri_color);
  GfxDebug::draw_polygon(my_points, 3, 8/*CUBIT_GREEN*/, CUBIT_GREEN);
}


//////////////////////////////////////////////////////////////////////
///Function Name : get_boundary_facets
///Member Type   : Private
///Description   : Get all the boundary vertices and boundary facets 
///                representing the boundary of the RefVolume.
///Note          : Do not forget to delete the DLIList pointers to facet
///                data after done with them.
///Date          : 08/23/2002
//////////////////////////////////////////////////////////////////////
CubitStatus MedialTool3D::get_boundary_facets(
  Vec3DPArray &boundary_verts,
  Triangle3DPArray &boundary_facets)
{
  // get all the ref faces of the volume
  DLIList<RefFace*> ref_faces_list;
  CubitStatus get_faces_status = myRefVolume->ref_faces(ref_faces_list);
  if(get_faces_status != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Problem in getting the faces of volume geometry\n"
                  "Check volume %d\n", myRefVolume->id());
      return CUBIT_FAILURE;
    }

  // loop through each RefFace and get its facets
  int ii;
  RefFace *my_face;
  CubitStatus face_facet_stat;
  DLIList<CubitFacet*> my_facet_list;  // delete these after done
  DLIList<CubitPoint*> my_point_list;  // delete after done

  for(ii = ref_faces_list.size(); ii > 0; ii--)
  {
    my_face = ref_faces_list.get_and_step();

    // get the facets for the current face
    // get_facets appends the new vertices and facets to the exusting list
    face_facet_stat = my_face->get_facets( my_facet_list, my_point_list);

    if(face_facet_stat != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Problem in getting the facets of face geometry\n"
                  "Check surface %d in volume %d\n",
                  my_face->id(), myRefVolume->id());
      return CUBIT_FAILURE;
    }
  }

  // translate the points to Vec3D format and flag the my_point_list
  CubitPoint *my_point;
  Vec3D *new_vec3d;

  // resize the size of "boundary_verts" and "boundary_facets" to equal that
  // of "my_point_list" and "my_facet_list"
  boundary_verts.resize(my_point_list.size());
  boundary_facets.resize(my_facet_list.size());

  // set all the pointers in vector list to NULL for safety
  for(ii=0; ii < boundary_verts.size(); ii++)
    boundary_verts[ii] = NULL;
  for(ii=0; ii < boundary_facets.size(); ii++)
    boundary_facets[ii] = NULL;

  // location of the current point
  CubitVector pt_coord;

  // go to each cubit point; set a flag and create corresponding new Vec3D
  for(ii = my_point_list.size(); ii > 0; ii--)
  {
    my_point = my_point_list.get_and_step();
    // set the mark number
    my_point->marked(ii-1);
    // create new Vec3D
    pt_coord = my_point->coordinates();
    new_vec3d = new Vec3D(pt_coord.x(), pt_coord.y(), pt_coord.z());
    boundary_verts[ii-1] = new_vec3d;
  }

  // translate the facets to Triangle3D format
  CubitFacet *my_facet;
  CubitPoint *vert0, *vert1, *vert2;
  Vec3D my_vert0, my_vert1, my_vert2;
  Triangle3D *my_tri;

  // go through each facet and create new Triangle3D
  for(ii = my_facet_list.size(); ii > 0; ii--)
  {
    my_facet = my_facet_list.get_and_step();
    // get the three points that form the facet and also their id's
    vert0 = my_facet->point(0);
    vert1 = my_facet->point(1);
    vert2 = my_facet->point(2);

    // create new Triangle3D
    my_tri = new Triangle3D(boundary_verts,
            vert0->marked(), vert1->marked(), vert2->marked());
    // assign the Triangle3D to its position in vector list
    boundary_facets[ii-1] = my_tri;
  }

  // delete the facet entities
  for(ii=my_facet_list.size(); ii >0; ii--)
    delete my_facet_list.pop();
  my_facet_list.clean_out();

  // delete the vertex entities
  for(ii=my_point_list.size(); ii>0; ii--)
    delete my_point_list.pop();
  my_point_list.clean_out();

  return CUBIT_SUCCESS;
}
#endif // USING_MEDIAL
