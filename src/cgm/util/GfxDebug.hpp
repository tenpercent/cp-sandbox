

#ifndef GFX_DEBUG_INTERFACE_HPP
#define GFX_DEBUG_INTERFACE_HPP

// An interface for any tools, algorithms, etc... to draw things.
// NOTE: this is an "interface" in the sense that using it will not 
// create additional depedencies in your code.  If a subclass of GfxDebug
// is not instantiated, then this class acts as a stub with empty implementation.
//
// any graphics primitives drawn with this interface are temporary primitives and
// not permanent ones.  Any primitive or window created with this class are 
// owned by this class (or the user of this class in a sense) 
// and not by any other graphics model (ie GGeomModel, GMeshModel, etc...).


// We only have forward declarations and NO #includes !!! (except those already in the util folder)
template <class X> class DLIList;
class RefEntity;
class MeshEntity;
class CubitNode;
class CubitEdge;
class CubitFace;
class CubitTri;
class CubitTet;
class CubitHex;
class CubitPyramid;
//class CubitWedge;
struct GPoint;
class CubitFacet;
class CubitFacetEdge;
class GeoTet;
class GeoNode;
class CubitVector;
class RefEntity;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class Body;
class MRefEntity;
class MRefVertex;
class MRefEdge;
class MRefFace;
class MRefVolume;
class MBody;
class CubitTransformMatrix;
class CubitBox;

#include "CubitColorConstants.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT GfxDebug
{
  // the one instance of the Graphics Debug object
  // creating a second instance will delete and replace the first instance.
  static GfxDebug* mInstance;

  public:
    GfxDebug();
    virtual ~GfxDebug();

    // this will clear out all temporary data and delete any graphics windows created
    // with this interface.  It'll restore the persistent geom/mesh model.
    //static void reset();

    // this will clear out the geom/mesh model and only subsequent drawing
    // will be visible
    static void clear();

    // same as DrawingTool::display_all
    static void display_all(bool flush_highlight = true, bool import_mesh_display = false);

    // clear out the highlighted
    static void clear_highlight();

    // set the drawing mode
    static void drawing_mode(int mode);

    // rebuild all
    //static void rebuild_all();

  // Window Related Functions

    // creates a graphics window 
    // not creating a window will use an existing current window
    // all drawing goes only into all windows created here or an existing window
    static int create_window();
    // deletes a window created by this interface
    static void delete_window(int window_id);
    
    // switches a window created by this interface into a 2D drawing mode.
    static void set_2D_mode(int window_id);

    // causes the window(s) to be redrawn
    static void flush();

    // give control to the window(s) to handle mouse events
    static void mouse_xforms();

    static void display_to_world_length(double pixels, double& length);

  // Geometry Drawing Functions
    static void draw_ref_entity(RefEntity* entity, int color = CUBIT_DEFAULT_COLOR,
                                   CubitTransformMatrix* trans_mat = 0);
    static void draw_ref_vertex(RefVertex* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_edge(RefEdge* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_face(RefFace* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_face_edges(RefFace* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_volume(RefVolume* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_body(Body* entity, int color = CUBIT_DEFAULT_COLOR);
    
    //static void highlight_ref_entity(RefEntity* entity);
    static void highlight_ref_vertex(RefVertex* entity);
    static void highlight_ref_edge(RefEdge* entity);
    static void highlight_ref_face(RefFace* entity);
    static void highlight_ref_volume(RefVolume* entity);
    static void highlight_ref_body(Body* entity);
    
    static void draw_mref_entity(MRefEntity* entity, int color = CUBIT_DEFAULT_COLOR,
                                    bool draw_geom = true, bool draw_mesh = true, 
                                    CubitTransformMatrix* trans_mat = 0);
    static void draw_mref_vertex(MRefVertex* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_edge(MRefEdge* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_face(MRefFace* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_volume(MRefVolume* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_body(MBody* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_entity_mesh(MRefEntity* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_mref_volume_faces(MRefVolume* entity, int color = CUBIT_DEFAULT_COLOR);
    
    //static void highlight_mref_entity(MRefEntity* entity);
    static void highlight_mref_vertex(MRefVertex* entity);
    static void highlight_mref_edge(MRefEdge* entity);
    static void highlight_mref_face(MRefFace* entity);
    static void highlight_mref_volume(MRefVolume* entity);
    static void highlight_mref_body(MBody* entity);


  // Mesh Drawing Functions
    
    // draws mesh entities with a color from the cubit color table
    // prefer these functions to a loop drawing each entity separately
    //static void draw_mesh_entities(DLIList<MeshEntity*>* entities, int color);
    static void draw_nodes(DLIList<CubitNode*>* entities, int color);
    static void draw_edges(DLIList<CubitEdge*>* entities, int color);
    static void draw_quads(DLIList<CubitFace*>* entities, int color);
    static void draw_tris(DLIList<CubitTri*>* entities, int color);
    static void draw_tets(DLIList<CubitTet*>* entities, int color);
    static void draw_hexes(DLIList<CubitHex*>* entities, int color);
    //static void draw_pyramids(DLIList<CubitPyramid*>* entities, int color);
    //static void draw_wedges(DLIList<CubitWedge*>* entities, int color);

    static void draw_mesh_entity(MeshEntity* entity, int color);
    static void draw_node(CubitNode* entity, int color);
    static void draw_edge(CubitEdge* entity, int color);
    static void draw_quad(CubitFace* entity, int color);
    static void draw_tri(CubitTri* entity, int color);
    static void draw_tet(CubitTet* entity, int color);
    static void draw_hex(CubitHex* entity, int color);
    static void draw_pyramid(CubitPyramid* entity, int color);
    //static void draw_wedge(CubitWedge* entity, int color);
    //static void draw_node_array(CubitNode* nodes[], int color, bool ink_edges = true);

    static void highlight_node(CubitNode* entity);
    static void highlight_edge(CubitEdge* entity);
    static void highlight_quad(CubitFace* entity);
    static void highlight_tri(CubitTri* entity);
    static void highlight_tet(CubitTet* entity);
    static void highlight_hex(CubitHex* entity);
    static void highlight_pyramid(CubitPyramid* entity);
    //static void highlight_wedge(CubitWedge* entity);
    

  // Facet Drawing Functions
    static void draw_facet(CubitFacet* facet, int color, int draw_uv = 0);
    static void draw_facet_edge(CubitFacetEdge* facet_edge, int color);
    static void draw_geotet( GeoTet *tet, int color );
    static void draw_geonode( GeoNode *node, int color );

  // Generic Primitive Drawing Functions
    
    static void draw_box(CubitBox& box, int color);
    
    // draw point x,y,z of color
    static void draw_point(float x, float y, float z, int color);
    static void draw_point(const CubitVector& vector, int color);
    static void draw_point(CubitVector* vector, int color);

    // draw line of points {x,y,z}, {x,y,z} of color
    static void draw_line(float x1, float y1, float z1, 
                          float x2, float y2, float z2, int color);
    
    static void draw_line(const CubitVector& x, const CubitVector& y, int color);

    // draw a polyline with a number of points of a color
    static void draw_polyline(const GPoint* points, int num_points, int color);

    // draw a polygon with a number of points of a color and filled or not
    static void draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled = true);

    // draw a quad of a color
    static void draw_quad(const GPoint* points, int color);

    // draw a tri of a color
    static void draw_tri(const GPoint* points, int color);

    // draw a vector (tail -> head) of color
    static void draw_vector(const GPoint* tail, const GPoint* head, int color);
    static void draw_vector(const CubitVector& tail, const CubitVector& head, int color);

    // draw a shell (num points, point list, num faces, face list) of color
    static void draw_shell(int num_pts, const GPoint* points, 
                           int num_faces, const int* face_list, int color);

    // draw a label at position x,y,z of color
    static void draw_label(const char* label, float x, float y, float z, int color);
    static void draw_label(int label, float x, float y, float z, int color);

    // zoom to specified bounding box
    static void zoom(CubitBox &box);

    // highlight simple primitives
    static void highlight_points(int num_points, const float* xyzs);
    static void highlight_polylines(int num_points, const float* xyzs, int num_line_points, const int* line_list);
    static void highlight_polygons(int num_points, const float* xyzs, int num_face_points, const int* face_list );
    
  protected:

    // p is for protected in case you were wondering :)
    
    //virtual void p_reset() = 0;
    virtual void p_clear() = 0;
    virtual void p_display_all(bool,bool) = 0;
    virtual void p_clear_highlight() = 0;
    virtual void p_drawing_mode(int) = 0;
    //virtual void p_rebuild_all() = 0;
    
    virtual int p_create_window() = 0;
    virtual void p_delete_window(int) = 0;
    virtual void p_set_2D_mode(int) = 0;
    virtual void p_flush() = 0;
    virtual void p_mouse_xforms() = 0;
    virtual void p_display_to_world_length(double pixels, double& length) = 0;

    virtual void p_draw_ref_entity(RefEntity* entity, int color, CubitTransformMatrix* mat) = 0;
    virtual void p_draw_ref_vertex(RefVertex* entity, int color) = 0;
    virtual void p_draw_ref_edge(RefEdge* entity, int color) = 0;
    virtual void p_draw_ref_face(RefFace* entity, int color) = 0;
    virtual void p_draw_ref_face_edges(RefFace* entity, int color) = 0;
    virtual void p_draw_ref_volume(RefVolume* entity, int color) = 0;
    virtual void p_draw_ref_body(Body* entity, int color) = 0;
    virtual void p_draw_mref_entity(MRefEntity* entity, int color, bool geom, bool mesh, CubitTransformMatrix* mat) = 0;
    virtual void p_draw_mref_vertex(MRefVertex* entity, int color) = 0;
    virtual void p_draw_mref_edge(MRefEdge* entity, int color) = 0;
    virtual void p_draw_mref_face(MRefFace* entity, int color) = 0;
    virtual void p_draw_mref_volume(MRefVolume* entity, int color) = 0;
    virtual void p_draw_mref_body(MBody* entity, int color) = 0;
    virtual void p_draw_mref_entity_mesh(MRefEntity* entity, int color) = 0;
    virtual void p_draw_mref_volume_faces(MRefVolume* entity, int color) = 0;
    
    //virtual void p_highlight_ref_entity(RefEntity* entity) = 0;
    virtual void p_highlight_ref_vertex(RefVertex* entity) = 0;
    virtual void p_highlight_ref_edge(RefEdge* entity) = 0;
    virtual void p_highlight_ref_face(RefFace* entity) = 0;
    virtual void p_highlight_ref_volume(RefVolume* entity) = 0;
    virtual void p_highlight_ref_body(Body* entity) = 0;
    //virtual void p_highlight_mref_entity(MRefEntity* entity) = 0;
    virtual void p_highlight_mref_vertex(MRefVertex* entity) = 0;
    virtual void p_highlight_mref_edge(MRefEdge* entity) = 0;
    virtual void p_highlight_mref_face(MRefFace* entity) = 0;
    virtual void p_highlight_mref_volume(MRefVolume* entity) = 0;
    virtual void p_highlight_mref_body(MBody* entity) = 0;


    //virtual void p_draw_mesh_entities(DLIList<MeshEntity*>*, int) = 0;
    virtual void p_draw_nodes(DLIList<CubitNode*>*, int) = 0;
    virtual void p_draw_edges(DLIList<CubitEdge*>*, int) = 0;
    virtual void p_draw_quads(DLIList<CubitFace*>*, int) = 0;
    virtual void p_draw_tris(DLIList<CubitTri*>*, int) = 0;
    virtual void p_draw_tets(DLIList<CubitTet*>*, int) = 0;
    virtual void p_draw_hexes(DLIList<CubitHex*>*, int) = 0;
    //virtual void p_draw_pyramids(DLIList<CubitPyramid*>*, int) = 0;
    //virtual void p_draw_wedges(DLIList<CubitWedge*>*, int) = 0;
    virtual void p_draw_mesh_entity(MeshEntity*, int) = 0;
    virtual void p_draw_node(CubitNode*, int) = 0;
    virtual void p_draw_edge(CubitEdge*, int) = 0;
    virtual void p_draw_quad(CubitFace*, int) = 0;
    virtual void p_draw_tri(CubitTri*, int) = 0;
    virtual void p_draw_tet(CubitTet*, int) = 0;
    virtual void p_draw_hex(CubitHex*, int) = 0;
    virtual void p_draw_pyramid(CubitPyramid*, int) = 0;
    //virtual void p_draw_wedge(CubitWedge*, int) = 0;
    virtual void p_highlight_node(CubitNode* entity) = 0;
    virtual void p_highlight_edge(CubitEdge* entity) = 0;
    virtual void p_highlight_quad(CubitFace* entity) = 0;
    virtual void p_highlight_tri(CubitTri* entity) = 0;
    virtual void p_highlight_tet(CubitTet* entity) = 0;
    virtual void p_highlight_hex(CubitHex* entity) = 0;
    virtual void p_highlight_pyramid(CubitPyramid* entity) = 0;
    //virtual void p_highlight_wedge(CubitWedge* entity) = 0;

    virtual void p_draw_facet(CubitFacet*, int,int) = 0;
    virtual void p_draw_facet_edge(CubitFacetEdge*, int) = 0;
    virtual void p_draw_geotet(GeoTet*, int) = 0;
    virtual void p_draw_geonode(GeoNode *, int) = 0;
    
    virtual void p_draw_box(CubitBox& box, int color) = 0;

    virtual void p_draw_point(float, float, float, int) = 0;
    virtual void p_draw_point(const CubitVector& vector, int color) = 0;
    virtual void p_draw_line(float, float, float, float, float, float, int) = 0;
    virtual void p_draw_line(const CubitVector& x, const CubitVector& y, int color) = 0;
    virtual void p_draw_polyline(const GPoint*, int, int) = 0;
    virtual void p_draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled) = 0;
    virtual void p_draw_quad(const GPoint*, int) = 0;
    virtual void p_draw_tri(const GPoint*, int) = 0;
    virtual void p_draw_vector(const CubitVector&, const CubitVector&, int) = 0;
    virtual void p_draw_shell(int, const GPoint*, int, const int*, int) = 0;
    virtual void p_draw_label(const char*, float, float, float, int) = 0;
    virtual void p_draw_label(int, float, float, float, int) = 0;
    virtual void p_zoom(CubitBox &box) = 0;
    
    virtual void p_highlight_points(int num_points, const float* xyzs) = 0;
    virtual void p_highlight_polylines(int num_points, const float* xyzs, int num_line_points, const int* line_list) = 0;
    virtual void p_highlight_polygons(int num_points, const float* xyzs, int num_face_points, const int* face_list ) = 0;

};

#endif


