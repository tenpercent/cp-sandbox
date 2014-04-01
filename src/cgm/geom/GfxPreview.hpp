

#ifndef GFX_PREVIEW_INTERFACE_HPP
#define GFX_PREVIEW_INTERFACE_HPP

// An interface for any tools, algorithms, etc... to draw preview objects.
// NOTE: this is an "interface" in the sense that using it will not 
// create additional depedencies in your code.  If a subclass of GfxPreview
// is not instantiated, then this class acts as a stub with empty implementation.
// CGM users may implement preview code by deriving from this class and implementing
// graphics code in the pure virtual p_* methods.  Creating an instance of the
// derived class will result in the correct methods being called by this interface.

// PLEASE!!! use this interface to implement any graphical preview code in CGM.
//           This will identify graphics code that must be maintained for end users.

//
// any graphics primitives drawn with this interface are temporary primitives and
// not permanent ones.  All items are draw into a single drawing set owned by the 
// 


// We only have forward declarations and NO #includes !!! (except those already in the util folder)
template <class X> class DLIList;
class RefEntity;
struct GPoint;
class CubitVector;
class RefEntity;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class Body;
class CubitTransformMatrix;
class CubitBox;

#include "CubitColorConstants.hpp"
#include "CubitGeomConfigure.h"

class CUBIT_GEOM_EXPORT GfxPreview
{
  // the one instance of the Graphics Preview object
  // creating a second instance will delete and replace the first instance.
  static GfxPreview* mInstance;

  public:
    GfxPreview();
    virtual ~GfxPreview();

    // this will clear out all temporary data and delete any graphics windows created
    // with this interface.  It'll restore the persistent geom/mesh model.
    //static void reset();

    // this will clear out any objects drawn using preview graphics.
    // no other graphics will be affected so you can now do something like
    //   draw vol 4
    //   webcut vol all plane x preview
    // and then clear the preview and you will still only have vol 4 displaye.
    static void clear();  

    // causes the window(s) to be redrawn
    static void flush();

  // Geometry Drawing Functions
    static void draw_ref_entity(RefEntity* entity, int color = CUBIT_DEFAULT_COLOR,
                                   CubitTransformMatrix* trans_mat = 0);
    static void draw_ref_vertex(RefVertex* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_edge(RefEdge* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_face(RefFace* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_face_edges(RefFace* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_volume(RefVolume* entity, int color = CUBIT_DEFAULT_COLOR);
    static void draw_ref_body(Body* entity, int color = CUBIT_DEFAULT_COLOR);

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
    // draw a list of polygons (much more efficent than drawing a single polygon multiple times.)
    // pass in the coordinates and the connectivity
    // the face_list is comprised of a list num_points_in_face, face_pt1, face_pt2, face_pt3, ...
    static void draw_polygons(int num_points, const float* xyzs, int num_face_points, 
                              const int* face_list, int color);

    // draw a quad of a color
    static void draw_quad(const GPoint* points, int color);

    // draw a tri of a color
    static void draw_tri(const GPoint* points, int color);

    // draw a vector (tail -> head) of color
    static void draw_vector(const GPoint* tail, const GPoint* head, int color);
    static void draw_vector(const CubitVector& tail, const CubitVector& head, int color);

    // draw a label at position x,y,z of color
    static void draw_label(const char* label, float x, float y, float z, int color);
    static void draw_label(int label, float x, float y, float z, int color);

  protected:

    // p is for protected in case you were wondering :)
    
    virtual void p_clear() = 0;
    virtual void p_flush() = 0;

    virtual void p_draw_ref_entity(RefEntity* entity, int color, CubitTransformMatrix* mat) = 0;
    virtual void p_draw_ref_vertex(RefVertex* entity, int color) = 0;
    virtual void p_draw_ref_edge(RefEdge* entity, int color) = 0;
    virtual void p_draw_ref_face(RefFace* entity, int color) = 0;
    virtual void p_draw_ref_face_edges(RefFace* entity, int color) = 0;
    virtual void p_draw_ref_volume(RefVolume* entity, int color) = 0;
    virtual void p_draw_ref_body(Body* entity, int color) = 0;
    virtual void p_draw_box(CubitBox& box, int color) = 0;
    virtual void p_draw_point(float, float, float, int) = 0;
    virtual void p_draw_point(const CubitVector& vector, int color) = 0;
    virtual void p_draw_line(float, float, float, float, float, float, int) = 0;
    virtual void p_draw_line(const CubitVector& x, const CubitVector& y, int color) = 0;
    virtual void p_draw_polyline(const GPoint*, int, int) = 0;
    virtual void p_draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled) = 0;
    virtual void p_draw_polygons(int num_points, const float* xyzs, int num_face_points, 
                                 const int* face_list, int color) = 0;
    virtual void p_draw_quad(const GPoint*, int) = 0;
    virtual void p_draw_tri(const GPoint*, int) = 0;
    virtual void p_draw_vector(const CubitVector&, const CubitVector&, int) = 0;
    virtual void p_draw_label(const char*, float, float, float, int) = 0;
    virtual void p_draw_label(int, float, float, float, int) = 0;
};

#endif


