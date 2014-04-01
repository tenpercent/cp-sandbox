

#include "GfxPreview.hpp"
#include "GMem.hpp"
#include "CubitVector.hpp"
#include <assert.h>

GfxPreview* GfxPreview::mInstance = 0;

GfxPreview::GfxPreview()
{
  if(mInstance)
  {
    assert(0);
    // if you want a new instance in place of a previous one,
    // delete the previous one first
  }
  mInstance = this;
}

GfxPreview::~GfxPreview()
{
  mInstance = 0;
}


// This clears preview data ONLY
void GfxPreview::clear()
{
  if(!mInstance) return;
  mInstance->p_clear();
}

// causes the current window to be redrawn
void GfxPreview::flush()
{
  if(!mInstance) return;
  mInstance->p_flush();
}

// Geometry Drawing Functions
    
void GfxPreview::draw_ref_entity(RefEntity* entity, int color, CubitTransformMatrix* mat)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_entity(entity, color, mat);
}

void GfxPreview::draw_ref_vertex(RefVertex* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_vertex(entity, color);
}

void GfxPreview::draw_ref_edge(RefEdge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_edge(entity, color);
}

void GfxPreview::draw_ref_face(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face(entity, color);
}

void GfxPreview::draw_ref_face_edges(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face_edges(entity, color);
}

void GfxPreview::draw_ref_volume(RefVolume* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_volume(entity, color);
}

void GfxPreview::draw_ref_body(Body* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_body(entity, color);
}

// Generic Primitive Drawing Functions

void GfxPreview::draw_box(CubitBox& box, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_box(box, color);
}

// draw point x,y,z of color
void GfxPreview::draw_point(float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(x, y, z, color);
}

void GfxPreview::draw_point(const CubitVector& vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(vector, color);
}

void GfxPreview::draw_point(CubitVector* vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(*vector, color);
}

// draw line of points {x,y,z}, {x,y,z} of color
void GfxPreview::draw_line(float x1, float y1, float z1, float x2, float y2, float z2, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x1, y1, z1, x2, y2, z2, color);
}

void GfxPreview::draw_line(const CubitVector& x, const CubitVector& y, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x, y, color);
}

// draw a polyline with a number of points of a color
void GfxPreview::draw_polyline(const GPoint* points, int num_points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_polyline(points, num_points, color);
}

// draw polygons given the coordinates and the connectivity.  Note that the
// connectivity (face_list) is given as the number of points in the polygon 
// followed by the point ids for that polygon.  So, the num_face_points is
// the total number of data points in the face_list such that
// num_face_points = number of points + number of polygons
void GfxPreview::draw_polygons(int num_points, const float* xyzs, int num_face_points, 
                               const int* face_list, int color )
{
  if (!mInstance) 
    return;
  mInstance->p_draw_polygons( num_points, xyzs, num_face_points, face_list, color);
}

// draw a polyline with a number of points of a color
void GfxPreview::draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled)
{
  if(!mInstance) return;
  mInstance->p_draw_polygon(points, num_points, color, border_color, filled);
}

// draw a quad of a color
void GfxPreview::draw_quad(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_quad(points, color);
}

// draw a tri of a color
void GfxPreview::draw_tri(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tri(points, color);
}

void GfxPreview::draw_vector(const CubitVector& tail, const CubitVector& head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(tail, head, color);
}

void GfxPreview::draw_vector(const GPoint* tail, const GPoint* head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(CubitVector(tail->x, tail->y, tail->z), CubitVector(head->x, head->y, head->z), color);
}

// draw a label at position x,y,z of color
void GfxPreview::draw_label(const char* label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}
void GfxPreview::draw_label(int label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}


