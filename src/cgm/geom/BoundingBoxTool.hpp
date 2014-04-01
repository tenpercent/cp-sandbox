//- Class: BoundingBoxTool
//- Description: Class for bounding boxes (primarily for "tight" bounding boxes)
//- Owner: Steve Storm
//- Created: 06-October-2000

#ifndef BoundingBoxTool_HPP
#define BoundingBoxTool_HPP

#include "CubitDefines.h"
#include "CubitGeomConfigure.h"

class GMem;
class CubitVector;
class RefEntity;

class CUBIT_GEOM_EXPORT BoundingBoxTool
{

public :   

   static CubitStatus get_tight_bounding_box( DLIList<RefEntity*> &ref_entity_list, 
                                       CubitVector &center,
                                       CubitVector axes[3],
                                       CubitVector &extension,
                                       double ang_facet_tol,
                                       double abs_facet_tol);
   //- Gets tightest box around the given entity list, using facets.
   //- Box is defined as center point, vector axes (unit vectors), and
   //-  extension in each axis (x extension is 1/2 width of box, for example).

   static CubitStatus get_axis_bounding_box( DLIList<RefEntity*> &ref_entity_list, 
                                      CubitVector &center,
                                      CubitVector axes[3],
                                      CubitVector &extension );
   //- Convenience function to get an axis box in the same output format as
   //- a tight box.

   static CubitStatus get_corner_points( CubitVector &center,
                                  CubitVector axes[3],
                                  CubitVector &extension,
                                  CubitVector& p1, CubitVector& p2,
                                  CubitVector& p3, CubitVector& p4,
                                  CubitVector& p5, CubitVector& p6,
                                  CubitVector& p7, CubitVector& p8);
   //- Using the center, axes and extension, gets the 8 corners of the box.

   static CubitBoolean get_use_triangles();
   static void set_use_triangles( CubitBoolean val );
   //- If triangles are used, surface facet points will be included in 
   //- the point list used to calculate the tight bounding box.  This 
   //- will include vertices and points on the curves.  This is the 
   //- default implementation.
   static int get_use_triangles_setting();
   static void set_use_triangles_setting( int val );

   static CubitBoolean get_use_curves();
   static void set_use_curves( CubitBoolean val );
   //- If curves are used, curve tesselation points will be included in
   //- the point list used to calculate the tight bounding box.  This
   //- includes the vertices on the ends of the curves.  One use for this
   //- is to find a more accurate tight bounding box, since curve 
   //- tessellations are typically more fine than surface tessellations.
   //- However, in practice, I would generally recommend just using surface
   //- tessellations.  One special case is if the user sends in a list of
   //- curves as the criteria for the tight bounding box, the curve 
   //- tessellations are always used, even if this parameter is false.
   static int get_use_curves_setting();
   static void set_use_curves_setting( int val );

   static CubitBoolean get_use_vertices();
   static void set_use_vertices( CubitBoolean val );
   //- If vertices are used, vertex points will be included in the point 
   //- list used to calculate the tight bounding box.  In extremely large
   //- models, it could be advantageous to just use vertices.  So the user
   //- would turn off both the surface and curve flags.  One special case 
   //- is if the user sends in a list of curves as the criteria for the 
   //- tight bounding box, the curve tessellations are always used, even 
   //- if the curve parameter is false and this parameter is true.
   static int get_use_vertices_setting();
   static void set_use_vertices_setting( int val );

 
   //Initialize all settings in this class
   static void initialize_settings();
    
protected :
   
   ~BoundingBoxTool();
   //- Destructor.
   
   BoundingBoxTool();
     //- Constructor for the BoundingBoxTool object
   
private :

   static CubitStatus append_ref_entity_points( DLIList<RefEntity*> &ref_entity_list, 
                                         DLIList<CubitVector*> &vec_list,
                                         double ang_tol, double abs_tol);
   static CubitStatus append_body_points( DLIList<Body*> &body_list, 
                                   DLIList<CubitVector*> &vec_list,
                                   double ang_tol, double abs_tol );
   static CubitStatus append_volume_points( DLIList<RefVolume*> &vol_list, 
                                     DLIList<CubitVector*> &vec_list,
                                     double ang_tol, double abs_tol );
   static CubitStatus append_curve_points( DLIList<RefEdge*> &ref_edge_list, 
                                    DLIList<CubitVector*> &vec_list );
   static CubitStatus append_vertex_points( DLIList<RefVertex*> &ref_vertex_list, 
                                     DLIList<CubitVector*> &vec_list );
   static CubitStatus append_surface_points( DLIList<RefFace*> &ref_face_list, 
                                      DLIList<CubitVector*> &vec_list,
                                      double ang_tol, double abs_tol );
   //- The "append" functions allocate memory for the vectors which
   //- must be freed by the calling code.

   static CubitBoolean expand_groups_in_list( DLIList<RefEntity*> &ref_entity_list );
   

  static  CubitBoolean useTriangles;
  static  CubitBoolean useCurves;
  static  CubitBoolean useVertices;
};



#endif



