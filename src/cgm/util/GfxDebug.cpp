

#include "GfxDebug.hpp"
#include "GMem.hpp"
#include "CubitVector.hpp"
#include <assert.h>

GfxDebug* GfxDebug::mInstance = 0;

GfxDebug::GfxDebug()
{
  if(mInstance)
  {
    assert(0);
    // if you want a new instance in place of a previous one,
    // delete the previous one first
  }
  mInstance = this;
}

GfxDebug::~GfxDebug()
{
  mInstance = 0;
}



// this will clear out all temporary data and any graphics windows created
// with this interface.
/*
void GfxDebug::reset()
{
  if(!mInstance) return;
  mInstance->p_reset();
}
*/

void GfxDebug::clear()
{
  if(!mInstance) return;
  mInstance->p_clear();
}

void GfxDebug::display_all(bool flush_highlight, bool import_mesh_display)
{
  if(!mInstance) return;
  mInstance->p_display_all(flush_highlight, import_mesh_display);
}

void GfxDebug::clear_highlight()
{
  if(!mInstance) return;
  mInstance->p_clear_highlight();
}

void GfxDebug::drawing_mode(int mode)
{
  if(!mInstance) return;
  mInstance->p_drawing_mode(mode);
}
/*
void GfxDebug::rebuild_all()
{
  if(!mInstance) return;
  mInstance->p_rebuild_all();
}
*/

// Window Related Functions

// creates a graphics window and makes it the current one
// not creating a window will use an existing current window
// all drawing goes only into the current window
int GfxDebug::create_window()
{
  if(!mInstance) return -1;
  return mInstance->p_create_window();
}

// deletes a window created by this interface
void GfxDebug::delete_window(int window_id)
{
  if(!mInstance) return;
  mInstance->p_delete_window(window_id);
}

// switches a window created by this interface into a 2D drawing mode.
void GfxDebug::set_2D_mode(int window_id)
{
  if(!mInstance) return;
  mInstance->p_set_2D_mode(window_id);
}

// causes the current window to be redrawn
void GfxDebug::flush()
{
  if(!mInstance) return;
  mInstance->p_flush();
}

// give control to the current window to handle mouse events
void GfxDebug::mouse_xforms()
{
  if(!mInstance) return;
  mInstance->p_mouse_xforms();
}

// compute display to world transformation of a length
void GfxDebug::display_to_world_length(double pixels, double& length)
{
  if(!mInstance) return;
  mInstance->p_display_to_world_length(pixels, length);
}

// Geometry Drawing Functions
    
void GfxDebug::draw_ref_entity(RefEntity* entity, int color, CubitTransformMatrix* mat)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_entity(entity, color, mat);
}

void GfxDebug::draw_ref_vertex(RefVertex* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_vertex(entity, color);
}

void GfxDebug::draw_ref_edge(RefEdge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_edge(entity, color);
}

void GfxDebug::draw_ref_face(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face(entity, color);
}

void GfxDebug::draw_ref_face_edges(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face_edges(entity, color);
}

void GfxDebug::draw_ref_volume(RefVolume* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_volume(entity, color);
}

void GfxDebug::draw_ref_body(Body* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_body(entity, color);
}


void GfxDebug::draw_mref_entity(MRefEntity* entity, int color, bool draw_geom, bool draw_mesh, CubitTransformMatrix* mat)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_entity(entity, color, draw_geom, draw_mesh, mat);
}

void GfxDebug::draw_mref_vertex(MRefVertex* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_vertex(entity, color);
}

void GfxDebug::draw_mref_edge(MRefEdge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_edge(entity, color);
}

void GfxDebug::draw_mref_face(MRefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_face(entity, color);
}

void GfxDebug::draw_mref_volume(MRefVolume* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_volume(entity, color);
}

void GfxDebug::draw_mref_body(MBody* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_body(entity, color);
}

void GfxDebug::draw_mref_entity_mesh(MRefEntity* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_entity_mesh(entity, color);
}

void GfxDebug::draw_mref_volume_faces(MRefVolume* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mref_volume_faces(entity, color);
}

/*
void GfxDebug::highlight_ref_entity(RefEntity* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_entity(entity);
}
*/

void GfxDebug::highlight_ref_vertex(RefVertex* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_vertex(entity);
}

void GfxDebug::highlight_ref_edge(RefEdge* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_edge(entity);
}

void GfxDebug::highlight_ref_face(RefFace* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_face(entity);
}

void GfxDebug::highlight_ref_volume(RefVolume* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_volume(entity);
}

void GfxDebug::highlight_ref_body(Body* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_ref_body(entity);
}

/*
void GfxDebug::highlight_mref_entity(MRefEntity* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_entity(entity);
}
*/

void GfxDebug::highlight_mref_vertex(MRefVertex* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_vertex(entity);
}

void GfxDebug::highlight_mref_edge(MRefEdge* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_edge(entity);
}

void GfxDebug::highlight_mref_face(MRefFace* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_face(entity);
}

void GfxDebug::highlight_mref_volume(MRefVolume* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_volume(entity);
}

void GfxDebug::highlight_mref_body(MBody* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_mref_body(entity);
}




// Mesh Drawing Functions

// draws mesh entities with a color from the cubit color table
// prefer these functions to a loop drawing each entity separately
/*
void GfxDebug::draw_mesh_entities(DLIList<MeshEntity*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mesh_entities(entities, color);
}
*/

void GfxDebug::draw_nodes(DLIList<CubitNode*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_nodes(entities, color);
}
void GfxDebug::draw_edges(DLIList<CubitEdge*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_edges(entities, color);
}
void GfxDebug::draw_quads(DLIList<CubitFace*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_quads(entities, color);
}
void GfxDebug::draw_tris(DLIList<CubitTri*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tris(entities, color);
}
void GfxDebug::draw_tets(DLIList<CubitTet*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tets(entities, color);
}
void GfxDebug::draw_hexes(DLIList<CubitHex*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_hexes(entities, color);
}
/*
void GfxDebug::draw_pyramids(DLIList<CubitPyramid*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_pyramids(entities, color);
}
void GfxDebug::draw_wedges(DLIList<CubitWedge*>* entities, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_wedges(entities, color);
}
*/
void GfxDebug::draw_mesh_entity(MeshEntity* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_mesh_entity(entity, color);
}

void GfxDebug::draw_node(CubitNode* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_node(entity, color);
}
void GfxDebug::draw_edge(CubitEdge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_edge(entity, color);
}
void GfxDebug::draw_quad(CubitFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_quad(entity, color);
}
void GfxDebug::draw_tri(CubitTri* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tri(entity, color);
}
void GfxDebug::draw_tet(CubitTet* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tet(entity, color);
}
void GfxDebug::draw_hex(CubitHex* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_hex(entity, color);
}
void GfxDebug::draw_pyramid(CubitPyramid* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_pyramid(entity, color);
}
/*
void GfxDebug::draw_wedge(CubitWedge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_wedge(entity, color);
}
*/

void GfxDebug::highlight_node(CubitNode* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_node(entity);
}

void GfxDebug::highlight_edge(CubitEdge* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_edge(entity);
}

void GfxDebug::highlight_quad(CubitFace* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_quad(entity);
}

void GfxDebug::highlight_tri(CubitTri* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_tri(entity);
}

void GfxDebug::highlight_tet(CubitTet* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_tet(entity);
}

void GfxDebug::highlight_hex(CubitHex* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_hex(entity);
}

void GfxDebug::highlight_pyramid(CubitPyramid* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_pyramid(entity);
}

/*
void GfxDebug::highlight_wedge(CubitWedge* entity)
{
  if(!mInstance) return;
  mInstance->p_highlight_wedge(entity);
}
*/


// Facet Drawing Functions

void GfxDebug::draw_facet(CubitFacet* facet, int color, int draw_uv )
{
  if(!mInstance) return;
  mInstance->p_draw_facet(facet, color, draw_uv);
}

void GfxDebug::draw_facet_edge(CubitFacetEdge* facet_edge, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_facet_edge(facet_edge, color);
}

void GfxDebug::draw_geotet(GeoTet *tet, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_geotet(tet, color);
}

void GfxDebug::draw_geonode(GeoNode *node, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_geonode(node, color);
}

// Generic Primitive Drawing Functions

void GfxDebug::draw_box(CubitBox& box, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_box(box, color);
}

// draw point x,y,z of color
void GfxDebug::draw_point(float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(x, y, z, color);
}

void GfxDebug::draw_point(const CubitVector& vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(vector, color);
}

void GfxDebug::draw_point(CubitVector* vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(*vector, color);
}

// draw line of points {x,y,z}, {x,y,z} of color
void GfxDebug::draw_line(float x1, float y1, float z1, float x2, float y2, float z2, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x1, y1, z1, x2, y2, z2, color);
}

void GfxDebug::draw_line(const CubitVector& x, const CubitVector& y, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x, y, color);
}

// draw a polyline with a number of points of a color
void GfxDebug::draw_polyline(const GPoint* points, int num_points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_polyline(points, num_points, color);
}

// draw a polyline with a number of points of a color
void GfxDebug::draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled)
{
  if(!mInstance) return;
  mInstance->p_draw_polygon(points, num_points, color, border_color, filled);
}

// draw a quad of a color
void GfxDebug::draw_quad(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_quad(points, color);
}

// draw a tri of a color
void GfxDebug::draw_tri(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tri(points, color);
}

void GfxDebug::draw_vector(const CubitVector& tail, const CubitVector& head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(tail, head, color);
}

void GfxDebug::draw_vector(const GPoint* tail, const GPoint* head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(CubitVector(tail->x, tail->y, tail->z), CubitVector(head->x, head->y, head->z), color);
}

// draw a shell (num points, point list, num faces, face list) of color
void GfxDebug::draw_shell(int num_pts, const GPoint* points, int num_faces, const int* face_list, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_shell(num_pts, points, num_faces, face_list, color);
}

// draw a label at position x,y,z of color
void GfxDebug::draw_label(const char* label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}
void GfxDebug::draw_label(int label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}

// zoom to specfied bounding box
void GfxDebug::zoom(CubitBox &box)
{
  if(!mInstance) return;
  mInstance->p_zoom(box);
}

void GfxDebug::highlight_points(int num_points, const float* xyzs)
{
  if(!mInstance) return;
  mInstance->p_highlight_points(num_points, xyzs);
}

void GfxDebug::highlight_polylines(int num_points, const float* xyzs, int num_line_points, const int* line_list)
{
  if(!mInstance) return;
  mInstance->p_highlight_polylines(num_points, xyzs, num_line_points, line_list);
}

void GfxDebug::highlight_polygons(int num_points, const float* xyzs, int num_face_points, const int* face_list )
{
  if(!mInstance) return;
  mInstance->p_highlight_polygons(num_points, xyzs, num_face_points, face_list);
}


