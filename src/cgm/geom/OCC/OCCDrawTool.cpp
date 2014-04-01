//- Class:       OCCDrawTool
//- Description: Functions to draw OCC geometry (useful for debugging or
//-              in previewing operations for the user)
// 
//- Author: Jane Hu       

// ********** BEGIN OCC INCLUDES             **********
#include "gp_Pnt.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Vertex.hxx"
#include "TopoDS_CompSolid.hxx"
#include "TopoDS_Solid.hxx"
#include "TopoDS_Shell.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "BRep_Tool.hxx"
// ********** END OCC INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "OCCDrawTool.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "OCCQueryEngine.hpp"
#include "DLIList.hpp"
#include "GMem.hpp"
#include "GfxPreview.hpp"
// ********** END CUBIT INCLUDES              **********

OCCDrawTool* OCCDrawTool::instance_ = 0;

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
OCCDrawTool* OCCDrawTool::instance()
{
  if (instance_ == 0) {
    instance_ = new OCCDrawTool();
  }
  return instance_;
}

OCCDrawTool::OCCDrawTool() 
{ 
}

OCCDrawTool::~OCCDrawTool()
{
}

CubitStatus
OCCDrawTool::draw_TopoDS_Shape( TopoDS_Shape *shape, int color, 
				 CubitBoolean tessellate,
                           	 CubitBoolean flush )
{
  if( (TopoDS_CompSolid*)shape  ||
      (TopoDS_Solid*)shape      ||
      (TopoDS_Shell*)shape  )
  {
    DLIList<TopoDS_Face*> Face_list;
    TopTools_IndexedMapOfShape M;
    TopExp::MapShapes(*shape, TopAbs_FACE, M);
    int ii;
    for (ii=1; ii<=M.Extent(); ii++) 
    {
          TopologyBridge *face = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          OCCSurface *occ_face = CAST_TO(face, OCCSurface);
          Face_list.append_unique(occ_face->get_TopoDS_Face());
    } 
    int i;
    for( i=Face_list.size(); i--; )
      draw_FACE( Face_list.get_and_step(), color, tessellate );
    if( flush )
      GfxPreview::flush();
    return CUBIT_SUCCESS;
  }
  else if ((TopoDS_Face*)shape )
  {
    return draw_FACE((TopoDS_Face *)shape , color, tessellate, flush );
  }
  else if( (TopoDS_Wire*)shape) 
  {
    DLIList<TopoDS_Edge*> EDGE_list;
    TopTools_IndexedMapOfShape M;
    TopExp::MapShapes(*shape, TopAbs_EDGE, M);
    int ii;
    for (ii=1; ii<=M.Extent(); ii++)
    {
          TopologyBridge *edge = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          OCCCurve *occ_curve = CAST_TO(edge, OCCCurve);
          EDGE_list.append_unique(occ_curve->get_TopoDS_Edge());
    }
    int i;
    for( i=EDGE_list.size(); i--; )
      draw_EDGE( EDGE_list.get_and_step(), color );
    if( flush )
      GfxPreview::flush();
    return CUBIT_SUCCESS;
  }
  else if ((TopoDS_Edge*)shape )
  {
    return draw_EDGE( (TopoDS_Edge*)shape, color, flush );
  }
  else if ((TopoDS_Vertex*)shape )
  {
    return draw_VERTEX( (TopoDS_Vertex*)shape, color, flush );
  }
  else
  {
    PRINT_ERROR( "Unsupported entity type specified - cannot draw\n" );
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
OCCDrawTool::draw_EDGE( TopoDS_Edge *EDGE_ptr, int color, CubitBoolean flush )
{
  OCCQueryEngine *OQE = OCCQueryEngine::instance();

  Curve *curve_ptr = OQE->populate_topology_bridge( *EDGE_ptr );

  CubitStatus stat;
  stat = draw_curve(curve_ptr, color, flush);

  return stat;
}

CubitStatus 
OCCDrawTool::draw_surface( Surface *surface, int color, CubitBoolean tessellate ,CubitBoolean flush)
{
  OCCQueryEngine *OQE = OCCQueryEngine::instance();
  if( tessellate )
  {
    int num_tris, num_pnts, num_facets;
    GMem g_mem ;
    OQE->get_graphics( surface, num_tris, num_pnts, num_facets,
      &g_mem );

    // Draw the triangles
    GPoint p[3];
    GPoint* plist = g_mem.point_list();
    int* facet_list = g_mem.facet_list();
    int i, c = 0;
    for( i=0; i<num_tris; i++ )
    {
      p[0] = plist[facet_list[++c]];
      p[2] = plist[facet_list[++c]];
      p[1] = plist[facet_list[++c]];
      c++;
      GfxPreview::draw_tri( p, color );
    }
  }

  else
  {
    // Draw curves
    DLIList<OCCCurve*> curve_list;
    CAST_TO(surface, OCCSurface)->get_curves(curve_list);
    int i;
    for( i=curve_list.size(); i--; )
      draw_curve( (Curve*)(curve_list.get_and_step()), color, flush );
  }

  if( flush )
    GfxPreview::flush();
  return CUBIT_SUCCESS;
}

CubitStatus 
OCCDrawTool::draw_curve( Curve *curve, int color , CubitBoolean flush )
{
  GMem g_mem;
  int num_points;
  OCCQueryEngine *OQE = OCCQueryEngine::instance();
  double tol = OQE->get_sme_resabs_tolerance();
  if (curve->get_arc_length() < tol)
    return CUBIT_SUCCESS;
  // get the graphics
  CubitStatus stat;
  stat = OQE->get_graphics( curve, num_points, &g_mem );

  if (stat==CUBIT_FAILURE || num_points == 0)
  {
    PRINT_ERROR("Unable to tessellate a curve for display\n" );
    return CUBIT_FAILURE;
  }
  else
  {
    // Draw the polyline
    GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, color );
  }

  if( flush )
    GfxPreview::flush();
  return CUBIT_SUCCESS;
}

CubitStatus
OCCDrawTool::draw_FACE( TopoDS_Face *FACE_ptr, int color, 
			CubitBoolean tessellate, 
                        CubitBoolean flush )
{
  OCCQueryEngine *OQE = OCCQueryEngine::instance();

  Surface *surf_ptr = OQE->populate_topology_bridge( *FACE_ptr );

  CubitStatus stat;
  stat = draw_surface(surf_ptr, color, tessellate, flush);

  return stat;
}

CubitStatus
OCCDrawTool::draw_VERTEX( TopoDS_Vertex *VERTEX_ptr, int color, CubitBoolean flush )
{
  gp_Pnt pt = BRep_Tool::Pnt(*VERTEX_ptr);
  GfxPreview::draw_point( pt.X(), pt.Y(), pt.Z(), color );
  return CUBIT_SUCCESS;
}
