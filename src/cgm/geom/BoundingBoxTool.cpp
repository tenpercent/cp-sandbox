//- Class: BoundingBoxTool
//- Description: Class for bounding boxes (primarily for "tight" bounding boxes)
//- Owner: Steve Storm
//- Created: 06-October-2000

#include "CubitBox.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "RefGroup.hpp"

#include "BoundingBoxTool.hpp"
#include "AnalyticGeometryTool.hpp"
#include "GMem.hpp"

#include "GeometryQueryEngine.hpp"

#include "DLIList.hpp"

#include "SettingHandler.hpp"

CubitBoolean BoundingBoxTool::useTriangles = CUBIT_TRUE;
CubitBoolean BoundingBoxTool::useCurves = CUBIT_FALSE;
CubitBoolean BoundingBoxTool::useVertices = CUBIT_FALSE;

BoundingBoxTool::BoundingBoxTool()
{}

// Destructor
BoundingBoxTool::~BoundingBoxTool()
{}

CubitStatus 
BoundingBoxTool::get_tight_bounding_box( DLIList<RefEntity*> &ref_entity_list,
                                         CubitVector &center,
                                         CubitVector axes[3],
                                         CubitVector &extension,
                                         double ang_facet_tol,
                                         double abs_facet_tol )
{
   DLIList<CubitVector*> vec_list;

   // Expand out groups in the list
   DLIList<RefEntity*> ref_entity_list_expanded = ref_entity_list;
   expand_groups_in_list( ref_entity_list_expanded );

   // Get the facet points from all the objects in the list
   append_ref_entity_points( ref_entity_list_expanded, vec_list, ang_facet_tol, abs_facet_tol );
   
   // Get the smallest box that fits around the points making up the facetted geometry
   AnalyticGeometryTool::instance()->get_tight_bounding_box( vec_list, center, 
                                                             axes, extension );

   // Free memory
   for( int i=0; i<vec_list.size(); i++ )
   {
      CubitVector* cubit_vector_ptr =
         vec_list.get_and_step();
      delete cubit_vector_ptr;
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::get_axis_bounding_box( DLIList<RefEntity*> &ref_entity_list,
                                        CubitVector &center,
                                        CubitVector axes[3],
                                        CubitVector &extension )
{
   CubitBoolean bounding_box_found = CUBIT_FALSE;
   CubitBox bounding_box;

   // Expand out groups in the list
   DLIList<RefEntity*> ref_entity_list_expanded = ref_entity_list;
   expand_groups_in_list( ref_entity_list_expanded );

   ref_entity_list_expanded.reset();
   for( int i = ref_entity_list_expanded.size(); i>0; i-- )
   {
      RefEntity* ref_entity_ptr = ref_entity_list_expanded.get_and_step();
      if( bounding_box_found == CUBIT_FALSE ) 
      {
         bounding_box = ref_entity_ptr->bounding_box();
         bounding_box_found = CUBIT_TRUE;
      }
      else {
         bounding_box |= ref_entity_ptr->bounding_box();
      }	    
   }

   axes[0].set( 1.0, 0.0, 0.0 );
   axes[1].set( 0.0, 1.0, 0.0 );
   axes[2].set( 0.0, 0.0, 1.0 );
   
   if( bounding_box_found ) 
   {
      extension.set( bounding_box.x_range()/2.0, bounding_box.y_range()/2.0,
         bounding_box.z_range()/2.0 );
      center = bounding_box.center();
   }
   else
   {
      extension.set( 0.0, 0.0, 0.0 );
      center.set( 0.0, 0.0, 0.0 );
   }

   return CUBIT_SUCCESS;   
}

CubitStatus 
BoundingBoxTool::append_ref_entity_points( DLIList<RefEntity*> &ref_entity_list, 
                                           DLIList<CubitVector*> &vec_list,
                                           double ang_facet_tol, double abs_facet_tol )
{
   DLIList<Body*> body_list;
   CAST_LIST( ref_entity_list, body_list, Body );
   if( body_list.size() )
      append_body_points( body_list, vec_list, ang_facet_tol, abs_facet_tol );

   DLIList<RefVolume*> vol_list;
   CAST_LIST( ref_entity_list, vol_list, RefVolume );
   if( vol_list.size() )
      append_volume_points( vol_list, vec_list, ang_facet_tol, abs_facet_tol );

   DLIList<RefFace*> surface_list;
   CAST_LIST( ref_entity_list, surface_list, RefFace );
   if( surface_list.size() )
      append_surface_points( surface_list, vec_list, ang_facet_tol, abs_facet_tol );

   DLIList<RefEdge*> curve_list;
   CAST_LIST( ref_entity_list, curve_list, RefEdge );
   if( curve_list.size() )
      append_curve_points( curve_list, vec_list );

   DLIList<RefVertex*> vertex_list;
   CAST_LIST( ref_entity_list, vertex_list, RefVertex );
   if( vertex_list.size() )
      append_vertex_points( vertex_list, vec_list );

   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::append_body_points( DLIList<Body*> &body_list, 
                                     DLIList<CubitVector*> &vec_list,
                                     double ang_tol, double abs_tol )
{
   int i;
   Body *body_ptr;
   DLIList<RefFace*> ref_face_list;

   body_list.reset();
   for( i=0; i<body_list.size(); i++ )
   {
      body_ptr = body_list.get_and_step();

      DLIList<RefFace*> ref_face_list_tmp;
      body_ptr->ref_faces( ref_face_list_tmp );

      ref_face_list.merge_unique( ref_face_list_tmp );
   }

   append_surface_points( ref_face_list, vec_list, ang_tol, abs_tol);
   
   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::append_volume_points( DLIList<RefVolume*> &vol_list, 
                                       DLIList<CubitVector*> &vec_list,
                                       double ang_tol, double abs_tol )
{
   int i;
   RefVolume *vol_ptr;
   DLIList<RefFace*> ref_face_list;

   vol_list.reset();
   for( i=0; i<vol_list.size(); i++ )
   {
      vol_ptr = vol_list.get_and_step();

      DLIList<RefFace*> ref_face_list_tmp;
      vol_ptr->ref_faces( ref_face_list_tmp );

      ref_face_list.merge_unique( ref_face_list_tmp );
   }

   append_surface_points( ref_face_list, vec_list, ang_tol, abs_tol);
   
   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::append_surface_points( DLIList<RefFace*> &ref_face_list, 
                                        DLIList<CubitVector*> &vec_list,
                                        double ang_tol, double abs_tol)
{
   int i;
   RefFace *ref_face_ptr;
   DLIList<RefEdge*> ref_edge_list;
   DLIList<RefVertex*> ref_vertex_list;
   
   ref_face_list.reset();
   for( i=0; i<ref_face_list.size(); i++ )
   {
      ref_face_ptr = ref_face_list.get_and_step();

      if( useTriangles == CUBIT_TRUE )
      {
         int num_tri, num_pnt, num_facet;
         GMem *g_mem = new GMem;
      
         ref_face_ptr->get_graphics( *g_mem, (int)ang_tol, abs_tol );
         num_pnt = g_mem->pointListCount;
         num_facet = g_mem->fListCount;
         num_tri = num_facet / 4;

         GPoint* point_list = g_mem->point_list();
         int num_pnts = g_mem->point_list_size();
         int y;
         for( y=0; y<num_pnts; y++ )
         {
            CubitVector* cubit_vector_ptr = new CubitVector(
               point_list[y].x, point_list[y].y, point_list[y].z );
            vec_list.append( cubit_vector_ptr );
         }
         delete g_mem;
      }
      
      if( useCurves == CUBIT_TRUE )
      {
         // Grab the curves from the surface
         DLIList<RefEdge*> ref_edge_list_tmp;
         ref_face_ptr->ref_edges( ref_edge_list_tmp );
         ref_edge_list.merge_unique( ref_edge_list_tmp );
      }

      // Only need vertices if both surfaces and curves are false
      if( useTriangles == CUBIT_FALSE && useCurves == CUBIT_FALSE &&
         useVertices == CUBIT_TRUE )
      {
         // Grab the vertices from the surface
         DLIList<RefVertex*> ref_vertex_list_tmp;
         ref_face_ptr->ref_vertices( ref_vertex_list_tmp );
         ref_vertex_list.merge_unique( ref_vertex_list_tmp );
      }
   }

   if( useCurves == CUBIT_TRUE )
      append_curve_points( ref_edge_list, vec_list );

   if( useTriangles == CUBIT_FALSE && useCurves == CUBIT_FALSE &&
       useVertices == CUBIT_TRUE )
       append_vertex_points( ref_vertex_list, vec_list );

   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::append_curve_points( DLIList<RefEdge*> &ref_edge_list, 
                                      DLIList<CubitVector*> &vec_list )
{
   // First grab all the vertices and append their coordinates to vec_list.
   // Then do the curves, but just put inside points in.  This avoids duplicate
   // points at the vertices.

   int i;
   RefEdge *ref_edge_ptr;
   DLIList<RefVertex*> ref_vertex_list;
   ref_edge_list.reset();
   for( i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
      DLIList<RefVertex*> ref_vertex_list_temp;
      ref_edge_ptr->ref_vertices( ref_vertex_list_temp );
      ref_vertex_list.merge_unique( ref_vertex_list_temp );
   }
   append_vertex_points( ref_vertex_list, vec_list );

   // Now add the *inside* curve points

   int j, num_pnts;
   ref_edge_list.reset();
   for( i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
      GMem *g_mem = new GMem;
      
      ref_edge_ptr->get_graphics( *g_mem );
      num_pnts = g_mem->pointListCount;

      GPoint* point_list = g_mem->point_list();

      for( j=1; j<num_pnts-1; j++ )
      {
         CubitVector* cubit_vector_ptr = new CubitVector(
            point_list[j].x, point_list[j].y, point_list[j].z );
         vec_list.append( cubit_vector_ptr );
      }
      delete g_mem;
   }

   return CUBIT_SUCCESS;
}

CubitStatus 
BoundingBoxTool::append_vertex_points( DLIList<RefVertex*> &ref_vertex_list, 
                                       DLIList<CubitVector*> &vec_list )
{
   int i;
   RefVertex *ref_vertex_ptr;
   ref_vertex_list.reset();
   for( i=0; i<ref_vertex_list.size(); i++ )
   {
      ref_vertex_ptr = ref_vertex_list.get_and_step();
      CubitVector coords = ref_vertex_ptr->coordinates();
      CubitVector* cubit_vector_ptr = new CubitVector( coords );
      vec_list.append( cubit_vector_ptr );
   }
   return CUBIT_SUCCESS;
}


CubitStatus 
BoundingBoxTool::get_corner_points( CubitVector &center,
                                    CubitVector axes[3],
                                    CubitVector &extension,
                                    CubitVector& p1, CubitVector& p2,
                                    CubitVector& p3, CubitVector& p4,
                                    CubitVector& p5, CubitVector& p6,
                                    CubitVector& p7, CubitVector& p8)
{
   double x = extension.x(); double y = extension.y(); double z = extension.z();

   // Front, Bottom, Left
   center.next_point( -axes[0], x, p1 ); p1.next_point( -axes[1], y, p1 );
   p1.next_point( axes[2], z, p1 );
   
   // Front, Top, Left
   center.next_point( -axes[0], x, p2 ); p2.next_point( axes[1], y, p2 );
   p2.next_point( axes[2], z, p2 );
   
   // Front, Top, Right
   center.next_point( axes[0], x, p3 ); p3.next_point( axes[1], y, p3 );
   p3.next_point( axes[2], z, p3 );
   
   // Front, Bottom, Right
   center.next_point( axes[0], x, p4 ); p4.next_point( -axes[1], y, p4 );
   p4.next_point( axes[2], z, p4 );
   
   // Back, Bottom, Left
   center.next_point( -axes[0], x, p5 ); p5.next_point( -axes[1], y, p5 );
   p5.next_point( -axes[2], z, p5 );
   
   // Back, Top, Left
   center.next_point( -axes[0], x, p6 ); p6.next_point( axes[1], y, p6 );
   p6.next_point( -axes[2], z, p6 );
   
   // Back, Top, Right
   center.next_point( axes[0], x, p7 ); p7.next_point( axes[1], y, p7 );
   p7.next_point( -axes[2], z, p7 );
   
   // Back, Bottom, Right
   center.next_point( axes[0], x, p8 ); p8.next_point( -axes[1], y, p8 );
   p8.next_point( -axes[2], z, p8 );

   return CUBIT_SUCCESS;
}


//Initialize all settings in this class
void BoundingBoxTool::initialize_settings()
{

  SettingHandler::instance()->add_setting("Tight Surface", 
                                          BoundingBoxTool::set_use_triangles_setting, 
					  BoundingBoxTool::get_use_triangles_setting ); 

  SettingHandler::instance()->add_setting("Tight Curve", 
                                          BoundingBoxTool::set_use_curves_setting, 
					  BoundingBoxTool::get_use_curves_setting);
   
  SettingHandler::instance()->add_setting("Tight Vertex", 
					 BoundingBoxTool::set_use_vertices_setting, 
					 BoundingBoxTool::get_use_vertices_setting);  


}



CubitBoolean 
BoundingBoxTool::expand_groups_in_list( DLIList<RefEntity*> &ref_entity_list )
{
   CubitBoolean group_found = CUBIT_FALSE;
   // just have to step through original entity_list, not the expanded list
   for (int i = ref_entity_list.size(); i>0; i--)  
   {
      RefEntity *entity = ref_entity_list.get();
      RefGroup *group = CAST_TO( entity, RefGroup );      
      if ( group ) 
      {
         ref_entity_list.remove();
         group->expand_group( ref_entity_list );
         group_found = CUBIT_TRUE;
      }
      else
         ref_entity_list.step();
   }
   return group_found;
}


int BoundingBoxTool::get_use_triangles_setting()
{return useTriangles;}

void BoundingBoxTool::set_use_triangles_setting( int val )
{useTriangles = (val) ? CUBIT_TRUE : CUBIT_FALSE;}

int BoundingBoxTool::get_use_curves_setting()
{return useCurves;}

void BoundingBoxTool::set_use_curves_setting( int val )
{useCurves = (val) ? CUBIT_TRUE : CUBIT_FALSE;}

int BoundingBoxTool::get_use_vertices_setting()
{return useVertices;}

void BoundingBoxTool::set_use_vertices_setting( int val )
{useVertices = (val) ? CUBIT_TRUE : CUBIT_FALSE;}


CubitBoolean BoundingBoxTool::get_use_triangles()
{return useTriangles;}

void BoundingBoxTool::set_use_triangles( CubitBoolean val )
{useTriangles=val;}

CubitBoolean BoundingBoxTool::get_use_curves()
{return useCurves;}

void BoundingBoxTool::set_use_curves( CubitBoolean val )
{useCurves=val;}

CubitBoolean BoundingBoxTool::get_use_vertices()
{return useVertices;}

void BoundingBoxTool::set_use_vertices( CubitBoolean val )
{useVertices=val;}


