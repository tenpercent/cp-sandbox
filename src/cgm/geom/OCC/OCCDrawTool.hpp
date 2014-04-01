#ifndef OCC_DRAW_TOOL_HPP
#define OCC_DRAW_TOOL_HPP

#include "CubitDefines.h"

class TopoDS_Shape;
class TopoDS_Face;
class TopoDS_Edge;
class TopoDS_Vertex;
class Surface;
class Curve;

template <class X> class DLIList;

class OCCDrawTool
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  ~OCCDrawTool();

  static OCCDrawTool* instance();
  //- Gives access to the singleton object of this class

  //do not remove this function....it is handy for debugging.
  CubitStatus draw_TopoDS_Shape( TopoDS_Shape *shape, int color,
                           CubitBoolean tessellate = CUBIT_FALSE,
                           CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given TopoDS_Shape.  Handles entities of type CompSolid, Solid, Shell,
  //- Face, Loop, Edge or Vertex.  The tessellate option only applies
  //- to Faces - if true, the Face is tessellated and each triangle is drawn;
  //- otherwise, only the outer curves are drawn.

  CubitStatus draw_surface( Surface *surface, int color, CubitBoolean tessellate, CubitBoolean flush );
  CubitStatus draw_curve( Curve *curve, int color, CubitBoolean flush );

  CubitStatus draw_FACE( TopoDS_Face *face, int color, 
                         CubitBoolean tessellate = CUBIT_FALSE,
                         CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given Face.  If tessellate is true, the FACE is tessellated and
  //- each triangle is drawn; otherwise, only the outer curves are drawn.

  CubitStatus draw_EDGE( TopoDS_Edge *edge, int color, 
                         CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given EDGE

  CubitStatus draw_VERTEX( TopoDS_Vertex *vertex, int color, 
                           CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given VERTEX

protected:
   OCCDrawTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

private:
  
  static OCCDrawTool* instance_;
};

#endif
