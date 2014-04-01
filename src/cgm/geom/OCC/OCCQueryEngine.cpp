//-------------------------------------------------------------------------
// Filename      : OCCQueryEngine.cpp
//
// Purpose       : Implementation of the OCCQueryEngine class.
//                 This class provides OCC-based implementations
//                 of various virtual functions in the GeometryQueryEngine
//                 hierarchy.
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 7/17/00
//
//-------------------------------------------------------------------------
#include <Standard_Stream.hxx>
//#include <Standard_SStream.hxx>
//#include <Standard_String.hxx>
//#include <stringbuf>
#include "BRep_Tool.hxx"
#include "gp_Pnt.hxx"
#include "gp_Ax1.hxx"
#include "gp_Ax2.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepBuilderAPI.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "OCCShapeAttributeSet.hpp"
//#include "OCCBinToolsShapeSet.hpp"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "Poly_Triangle.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "Handle_Poly_Triangulation.hxx"
#include "GCPnts_TangentialDeflection.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "STEPControl_Reader.hxx"
#include "IGESControl_Reader.hxx"
#include "STEPControl_Writer.hxx"
#include "IGESControl_Writer.hxx"
#include "STEPControl_StepModelType.hxx"
#include "IFSelect_ReturnStatus.hxx"
#include "BndLib_Add3dCurve.hxx"
#include "Poly_Polygon3D.hxx"
#include "Handle_Poly_Polygon3D.hxx"
#include "BRepMesh_FastDiscret.hxx"
#include "OCCQueryEngine.hpp"
#include "OCCModifyEngine.hpp"
#include "Poly_Triangulation.hxx"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "GeometryEntity.hpp"
#include "DLIList.hpp"
#include "CubitBox.hpp"
#include "CubitString.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCoFace.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "GMem.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitObserver.hpp"
#include "GfxDebug.hpp"
#include <stdio.h>
#include <errno.h>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include "TopoDS_Compound.hxx"
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <Precision.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <TDataStd_Shape.hxx>
#include <TDF_ChildIterator.hxx>
#include <BinTools_ShapeSet.hxx>
#include "Standard_Boolean.hxx"

#include "TDF_Label.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "TDocStd_Document.hxx"
#include "TCollection_ExtendedString.hxx"
#include "gp_Lin.hxx"
using namespace NCubitFile;

OCCQueryEngine* OCCQueryEngine::instance_ = NULL;

const int OCCQueryEngine::OCCQE_MAJOR_VERSION = 6;
const int OCCQueryEngine::OCCQE_MINOR_VERSION = 2;
const int OCCQueryEngine::OCCQE_SUBMINOR_VERSION = 0;

typedef std::map<int, TopologyBridge*>::value_type valType;
int OCCQueryEngine::iTotalTBCreated = 0;
int OCCQueryEngine::total_coedges = 0;
//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
OCCQueryEngine* OCCQueryEngine::instance()
{
  if (instance_ == NULL ) {
    instance_ = new OCCQueryEngine;
  }
  return instance_;
}

//================================================================================
// Description:  default constructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::OCCQueryEngine()
{
  GeometryQueryTool::instance()->add_gqe( this );
  OCCMap = new TopTools_DataMapOfShapeInteger;
  OccToCGM = new std::map<int, TopologyBridge*>;

  BodyList = new DLIList<OCCBody*>;
  WireList = new DLIList<OCCLoop*>;
  SurfaceList = new DLIList<OCCSurface*>;
  CurveList = new DLIList<OCCCurve*>;
  CubitString name("Doc");
  TCollection_ExtendedString xString((Standard_CString)name.c_str());
  MyDF = new TDocStd_Document(xString);
  mainLabel = MyDF->Main();
  EXPORT_ATTRIB = CUBIT_TRUE;
}

//================================================================================
// Description:  destructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::~OCCQueryEngine()
{
  instance_ = NULL;
  delete OCCMap;
  delete OccToCGM;
  delete BodyList;
  delete WireList;
  delete SurfaceList;
  delete CurveList;
  delete MyDF;
}

int OCCQueryEngine::get_major_version()
{
  return OCCQE_MAJOR_VERSION;
}

int OCCQueryEngine::get_minor_version()
{
  return OCCQE_MINOR_VERSION;
}

int OCCQueryEngine::get_subminor_version()
{
  return OCCQE_SUBMINOR_VERSION;
}

CubitString OCCQueryEngine::get_engine_version_string()
{
  CubitString version_string = "OCC Geometry Engine version ";
  version_string += CubitString(get_major_version());
  version_string += CubitString(".");
  version_string += CubitString(get_minor_version());
  version_string += CubitString(".");
  version_string += CubitString(get_subminor_version());
  
  return version_string;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
const type_info& OCCQueryEngine::entity_type_info() const
{
  return typeid(OCCQueryEngine);
}

//================================================================================
// Description: reflect body about an axis
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
CubitStatus OCCQueryEngine::reflect( BodySM *bodysm,
                                     const CubitVector& axis)
{
  OCCBody *body = CAST_TO(bodysm, OCCBody);
  if (!body)
    {
      PRINT_ERROR("Attempt to reflect OCC-based geometry Body.  This body is not OCCBody.");
      return CUBIT_FAILURE;
    }

  body->reflect( axis.x(), axis.y(), axis.z() );

  return CUBIT_SUCCESS;
}


//================================================================================
// Description:  This function queries OCC for the necessary facets
//               information needed in facetting a RefFace.  This
//               information is stored and output in gMem.  The
//               number of triangles, points and facets are also
//               output.
// Author     :  Jane Hu
// Date       :  10/25/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Surface* surface_ptr,
                                          int& number_triangles,
                                          int& number_points,
                                          int& number_facets,
                                          GMem* g_mem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length) const
{
  // Because this may be unnecessarily called twice,
  // say there is one triangle.
  if (!g_mem)
    {
      number_triangles = 1;
      number_points = 3;
      number_facets = 4;
      return CUBIT_SUCCESS;
    }

  OCCSurface *occ_surface_ptr = CAST_TO(surface_ptr, OCCSurface);
  TopoDS_Face * Topo_Face = occ_surface_ptr->get_TopoDS_Face();
  if (!Topo_Face)
    return CUBIT_FAILURE;

  TopLoc_Location L;
  Handle_Poly_Triangulation facets = BRep_Tool::Triangulation(*Topo_Face, L);

  gp_Trsf tf = L.Transformation();

  if(facets.IsNull() || facets->NbTriangles() == 0)
  {
    //do triangulation
    double deflection = 0.01;
    double angle  = 0.5;
    BRepAdaptor_Surface asurface(*Topo_Face);
    Bnd_Box aBox;
    BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
    BRepMesh_FastDiscret *myMesh =
    new BRepMesh_FastDiscret(deflection, *Topo_Face, aBox, angle, Standard_True, Standard_True);
    if (myMesh != NULL) delete myMesh;
    facets = BRep_Tool::Triangulation(*Topo_Face, L);
    if(facets.IsNull() || facets->NbTriangles() == 0)
    {
      PRINT_ERROR("Can't get triangulation representation for this surface.\n");
      return CUBIT_FAILURE;
    }
  }
  //if necessary, the face tolerance can be returned. now, no use.
  //double tol = BRep_Tool::Tolerance(*Topo_Face);   

  number_points = facets->NbNodes();
  number_triangles = facets->NbTriangles();
  number_facets = 4 * number_triangles; 
  
  Poly_Array1OfTriangle triangles(0, number_triangles-1);
  triangles.Assign( facets->Triangles() );
  int *facetList =  new int[number_facets];
  //needs to test that N1, N2, N3 index are starting from 0 to number_points-1
  //otherwise needs to update either facetList or gPnts to make consistent.
  //It's possible also that N's starting from 1.
  int minN = 100;
  for (int i = 0; i < triangles.Length(); i++)
    {
      Poly_Triangle triangle = triangles.Value( i );
      int N1, N2, N3;
      triangle.Get(N1, N2, N3); 
      facetList[4 * i] = 3;
      facetList[4 * i + 1] = N1;
      minN = (minN < N1 ? minN : N1);
      facetList[4 * i + 2] = N2;
      minN = (minN < N2 ? minN : N2);
      facetList[4 * i + 3] = N3;
      minN = (minN < N3 ? minN : N3);
    } 
  if(minN != 0)
  {
    //subtract the minN from the facetList for all i+1, i+2, i+3 points
    for (int i = 0; i < triangles.Length(); i++)
    {
      facetList[4 * i + 1] -= minN;
      facetList[4 * i + 2] -= minN;
      facetList[4 * i + 3] -= minN;
    }
  }
  g_mem->replace_facet_list( facetList, number_facets, number_facets); 

  TColgp_Array1OfPnt points(0,  number_points-1);
  points.Assign(facets->Nodes());
  GPoint *gPnts= new GPoint[number_points];
  for (int i = 0; i < number_points ; i ++)
    {
      gp_Pnt gp_pnt = points.Value(i);
      if( !L.IsIdentity())
        gp_pnt.Transform(tf);

      GPoint gPnt;
      gPnt.x = gp_pnt.X();
      gPnt.y = gp_pnt.Y();
      gPnt.z = gp_pnt.Z();
      gPnts[i] = gPnt;
    }
  g_mem->replace_point_list( gPnts, number_points, number_points );

  return CUBIT_SUCCESS;
}

//================================================================================
// Description: This function queries OCC for the edge information
//              needed in facetting a RefEdge.  This information is
//              stored and output in g_mem.
// Author     : Jane Hu
// Date       : 10/26/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Curve* curve_ptr,
                                          int& num_points,
                                          GMem* gMem,
                                          double /*tolerance*/ ) const
{
  //  get the OCCCurve.
  OCCCurve *occ_curve_ptr = CAST_TO(curve_ptr,OCCCurve);
  assert (gMem);

    
  TopoDS_Edge * Topo_Edge = occ_curve_ptr->get_TopoDS_Edge();
  if (!Topo_Edge)
    return CUBIT_FAILURE;

  //do tessellation
  double deflection = 0.2;
  double angle  = 0.2;
  BRepAdaptor_Curve acurve(*Topo_Edge);
  GCPnts_TangentialDeflection *myMesh = 
        new GCPnts_TangentialDeflection(acurve, angle, deflection);
  if (myMesh == NULL) 
  {
    PRINT_ERROR("Can't tessellate for this curve.\n");
    return CUBIT_FAILURE;
  }
  num_points = myMesh->NbPoints();

  //! Note: If the polygon is closed, the point of closure is 
  //! repeated at the end of its table of nodes. Thus, on a closed 
  //! triangle the function NbNodes returns 4? 
  GPoint *gPnts= new GPoint[num_points];
  for (int i = 1; i <= num_points ; i ++)
    {
      gp_Pnt gp_pnt = myMesh->Value(i);
      GPoint gPnt;
      gPnt.x = gp_pnt.X();
      gPnt.y = gp_pnt.Y();
      gPnt.z = gp_pnt.Z();
      gPnts[i-1] = gPnt;
    }
  gMem->replace_point_list( gPnts, num_points, num_points );
 
  delete myMesh;
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Given surface and number of point on u and v parametric
//              direction, find the 3-d point locations
// Author     : Jane Hu
// Date       : 10/22/07
//================================================================================
CubitStatus OCCQueryEngine::get_isoparametric_points(Surface* surface,
                                                     int &nu, int &nv,
                                                     GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
    {
      PRINT_ERROR("This surface is not OCCSurface.");
      return CUBIT_FAILURE;
    }

  the_face = *Tops_face;
  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  assert (nu > 1 && nv > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2); 
  double interval1  = (u2 - u1)/(nu-1);
  double interval2 = (v2 - v1)/(nv-1);
  
  g_mem = new GMem[nu];
  //nu and nv must be given to calculate the points.
  for (int i = 0; i < nu; i++)
    {
      Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u1 + i * interval1); 
      g_mem[i].allocate_polylines(nv-1);
      for (int j = 0; j <  nv; j++)
	{
	  gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval2);
	  g_mem[i].point_list()[j].x = pnt.X();
	  g_mem[i].point_list()[j].y = pnt.Y();
	  g_mem[i].point_list()[j].z = pnt.Z();
	}
      g_mem[i].pointListCount = nv;

    }
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_u_isoparametric_points(Surface* surface,
                                                       double v, int&n,
                                                       GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
    {
      PRINT_ERROR("This surface is not OCCSurface.");
      return CUBIT_FAILURE;
    }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  //n must be given to calculate the points.
  assert (n > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (u2 - u1)/(n -1); 
  
  Handle_Geom_Curve HGeom_curve = HGeom_surface->VIso(v);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
    {
      gp_Pnt pnt = HGeom_curve->Value(u1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
    }
  g_mem->pointListCount = n;

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_v_isoparametric_points(Surface* surface,
                                                       double u, int&n,
                                                       GMem*&g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
    {
      PRINT_ERROR("This surface is not OCCSurface.");
      return CUBIT_FAILURE;
    }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  //n must be given to calculate the points.
  assert (n > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (v2 - v1)/(n -1);

  Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
    {
      gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
    }
  g_mem->pointListCount = n;

  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus OCCQueryEngine::transform_vec_position( CubitVector const& ,
                                                    BodySM *,
                                                    CubitVector & ) const
{
  //Nobody is using this function in ACIS yet.
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description:  Calculate for intersection points between a curve and 
//               a segment defined by two points or
//               between two curves or between a curve and a surface.
// Author     :  Jane Hu
// Date       :  10/15/07
//================================================================================
CubitStatus OCCQueryEngine::get_intersections( Curve* curve, 
                                               CubitVector& point1,
                                               CubitVector& point2,
                                               DLIList<CubitVector*>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCPoint *pt1 = new OCCPoint(point1);
  OCCPoint *pt2 = new OCCPoint(point2);
  Curve *curve2 = 
    OCCModifyEngine::instance()->make_Curve(pt1, pt2);
  if (curve2 == NULL)
    {
      PRINT_ERROR( "Unable to create OCC EDGE from points\n" );
      return CUBIT_FAILURE;
    }
  
  OCCCurve *occ_curve2 = CAST_TO(curve2, OCCCurve);
  CubitStatus stat = get_intersections(occ_curve, occ_curve2, intscts, bounded, closest);
  delete_solid_model_entities(occ_curve2);
  return stat;
}

CubitStatus OCCQueryEngine::get_intersections( Curve* curve1, 
                                               Curve* curve2,
                                               DLIList<CubitVector*>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  //If this function has shortcomes in using BRepExtrema_DistShapeShape,
  //look also at IntTools_EdgeEdge.
  OCCCurve *occ_curve1 =  CAST_TO(curve1, OCCCurve);
  if (occ_curve1 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCCurve *occ_curve2 =  CAST_TO(curve2, OCCCurve);
  if (occ_curve2 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(
                                            *(occ_curve1->get_TopoDS_Edge()),
                                            *(occ_curve2->get_TopoDS_Edge()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curves.\n");
      return CUBIT_FAILURE;
    }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
    {
      int numPnt = distShapeShape.NbSolution();
      for (int i = 1; i <= numPnt; i++)
	{
	  gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
	  CubitVector* cv = new CubitVector(aPoint.X(), aPoint.Y(), aPoint.Z());
	  intscts.append(cv);
	}
    }
  return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::get_intersections( Curve* curve, Surface* surface,
                                   DLIList<CubitVector*>& intscts,
                                   CubitBoolean bounded )
{
  // There's no effect of bounded =  false. 
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCSurface *occ_surface =  CAST_TO(surface, OCCSurface);
  if (occ_surface == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }
   
  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(*(occ_curve->get_TopoDS_Edge()),
                                            *(occ_surface->get_TopoDS_Face()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curve and surface.\n");
      return CUBIT_FAILURE;
    }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
    {
      int numPnt = distShapeShape.NbSolution();
      for (int i = 1; i <= numPnt; i++)
	{
	  gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
	  CubitVector* cv = new CubitVector(aPoint.X(), aPoint.Y(), aPoint.Z());
	  intscts.append(cv);
	}
    }
 
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     : Jane Hu
// Date       : 10/30/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_extrema( DLIList<GeometryEntity*> &ref_entity_list,
				const CubitVector *dir1,
				const CubitVector *dir2,
				const CubitVector *dir3,
				CubitVector &extrema,
				GeometryEntity *&extrema_entity_ptr )
{
  //in Acis, the api_entity_extrema is used to calculate "possible 
  //self-intersecting sweeping and to align lofting sections"
  PRINT_ERROR("There's no such call in OCC correponding to Acis call."); 
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     : Jane Hu
// Date       : 10/19/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_entity_distance( GeometryEntity *entity1,
                                        GeometryEntity *entity2,
                                        CubitVector &pos1, CubitVector &pos2,
                                        double &distance )
{
  TopoDS_Shape * shape1;
  TopoDS_Shape * shape2;
  if ((shape1 = get_TopoDS_Shape_of_entity(entity1)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  if( (shape2 = get_TopoDS_Shape_of_entity( entity2 )) == NULL )
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n");
      return CUBIT_FAILURE;
    }

  BRepExtrema_DistShapeShape distShapeShape(*shape1, *shape2);
  //distShapeShape.Perform();
  
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR( "problem occured getting distance between OCC entities.\n"
		   "       Aborting.\n");
      return CUBIT_FAILURE;
    }

  distance = distShapeShape.Value();
  gp_Pnt pnt1 = distShapeShape.PointOnShape1(1);
  gp_Pnt pnt2 = distShapeShape.PointOnShape2(1);
  pos1 = CubitVector(pnt1.X(), pnt1.Y(), pnt1.Z());
  pos2 = CubitVector(pnt2.X(), pnt2.Y(), pnt2.Z());
  return CUBIT_SUCCESS;
}

TopoDS_Shape* OCCQueryEngine::get_TopoDS_Shape_of_entity(TopologyBridge *entity_ptr)
{
  if (OCCBody *body_ptr = CAST_TO( entity_ptr, OCCBody))
    {
      TopoDS_CompSolid* theShape = body_ptr->get_TopoDS_Shape();
      if (!theShape)
	{
	  PRINT_ERROR("OCCBody without TopoDS_Shape at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}
      return theShape;
    }

  else if (OCCLump * lump_ptr = CAST_TO( entity_ptr,OCCLump))
    {
      TopoDS_Solid * theSolid = lump_ptr->get_TopoDS_Solid();
      if(theSolid)
	return (TopoDS_Shape*) theSolid; 
      else
	{
	  PRINT_ERROR("OCCLump without TopoDS_Solid at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}
    }

  else if( OCCSurface *surface_ptr = CAST_TO( entity_ptr, OCCSurface))
    {
      TopoDS_Face *theFace = surface_ptr->get_TopoDS_Face();
      if(!theFace)
	{
	  PRINT_ERROR("OCCSurface without TopoDS_Face at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) theFace;
    }

  else if( OCCCurve *curve_ptr = CAST_TO( entity_ptr, OCCCurve))
    {
      TopoDS_Edge *theEdge = curve_ptr->get_TopoDS_Edge();
      if (!theEdge)
	{
	  PRINT_ERROR("OCCCurve without TopoDS_Edge at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) theEdge;
    }

  else if( OCCPoint *point_ptr = CAST_TO( entity_ptr, OCCPoint))
    {
      TopoDS_Vertex *thePoint = point_ptr->get_TopoDS_Vertex(); 
      if (!thePoint)
	{
	  PRINT_ERROR("OCCPoint without TopoDS_Point at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) thePoint;
    }
  
  PRINT_ERROR("Non-OCC TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return NULL;

}
//===========================================================================
//Function Name: save_temp_geom_file
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary FACET file
//              If file is created, CubitString 'created_file' and
//              'create_file_type' are set
//Author:       Corey Ernst
//Date:         1/18/2003
//===========================================================================

CubitStatus OCCQueryEngine::save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
						 const char *file_name,
						 const CubitString &cubit_version,
						 CubitString &created_file,
						 CubitString &created_file_type)
{
  int size_before = ref_entity_list.size();
  CubitString temp_filename(file_name);
  temp_filename += ".occ";

  if( export_solid_model( ref_entity_list, temp_filename.c_str(), "OCC",
                          cubit_version ) == CUBIT_FAILURE )
    {
      PRINT_ERROR( "Error occured while trying to save temporary OCC_BASED_GEOMETRY file\n");
      return CUBIT_FAILURE;
    }

  int size_after = ref_entity_list.size();

  if( size_before > size_after )
    {
      created_file +=  temp_filename;
      created_file_type += "OCC";
    }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:export_solid_model
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
						const char* file_name,
						const char* file_type,
						const CubitString &,
						const char*)
{
  if( strcmp( file_type, "OCC" ) != 0 && 
      strcmp( file_type, "STEP") != 0 &&
      strcmp( file_type, "IGES") != 0)
    {
      //PRINT_ERROR("The specified file type, %s, is not supported!\n", filetype );
      return CUBIT_FAILURE;
    }
 
  DLIList<OCCBody*>    OCC_bodies;
  DLIList<OCCSurface*> OCC_surfaces;
  DLIList<OCCCurve*>   OCC_curves;
  DLIList<OCCPoint*>   OCC_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
    {
      TopologyBridge* ref_entity_ptr = ref_entity_list.get();
      CubitBoolean handled = CUBIT_TRUE;

      //if it is a Vertex
      if( OCCPoint* pt = CAST_TO( ref_entity_ptr, OCCPoint) )
	OCC_points.append( pt );

      //if it is a Curve
      else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
	OCC_curves.append( curve );
    
      //if it is a surface
      else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
	OCC_surfaces.append( surf );
   
      //if it is a Body
      else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
	OCC_bodies.append( body );

      else
	handled = CUBIT_FALSE;

      if( handled == CUBIT_TRUE )
	{
	  ref_entities_handled.append( ref_entity_ptr );
	  ref_entity_list.change_to(NULL);
	}

      ref_entity_list.step();
    }

  ref_entity_list.remove_all_with_value(NULL);

  int free_body_count = OCC_bodies.size();
  int free_curve_count = OCC_curves.size();
  int free_point_count = OCC_points.size();
  int free_surface_count = OCC_surfaces.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_surface_count == 0 && 
      free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //save the facets (geometry info )
  CubitStatus status;

  //write out topology and attributes
  status = write_topology( file_name, file_type,
                           OCC_bodies, OCC_surfaces,
                           OCC_curves, OCC_points );
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_surface_count || 
      free_curve_count || free_point_count )
    PRINT_INFO( "\nExported:" );

  int flg = 0;

  if( free_body_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_body_count == 1 )
         PRINT_INFO( "%4d OCC Body to %s\n", free_body_count, file_name );
      else
         PRINT_INFO( "%4d OCC Bodies to %s\n", free_body_count, file_name );
    }

  if( free_surface_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_surface_count == 1 )
	PRINT_INFO( "%4d OCC Surface to %s\n", free_surface_count, file_name );
      else
	PRINT_INFO( "%4d OCC Surface to %s\n", free_surface_count, file_name );
    }

  if( free_curve_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
	PRINT_INFO( "%4d OCC Curve to %s\n", free_curve_count, file_name );
      else
	PRINT_INFO( "%4d OCC Curves to %s\n", free_curve_count, file_name );
    }

  if( free_point_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
	PRINT_INFO( "%4d OCC Point to %s\n", free_point_count, file_name );
      else
	PRINT_INFO( "%4d OCC Points to %s\n", free_point_count, file_name );
    }
  PRINT_INFO( "\n" );

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
						char*& p_buffer,
						int& n_buffer_size,
						bool b_export_buffer)
{
  DLIList<OCCBody*>    OCC_bodies;
  DLIList<OCCSurface*> OCC_surfaces;
  DLIList<OCCCurve*>   OCC_curves;
  DLIList<OCCPoint*>   OCC_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
    {
      TopologyBridge* ref_entity_ptr = ref_entity_list.get();
      CubitBoolean handled = CUBIT_TRUE;

      //if it is a Vertex
      if( OCCPoint* pt = CAST_TO( ref_entity_ptr, OCCPoint) )
	OCC_points.append( pt );

      //if it is a Curve
      else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
	OCC_curves.append( curve );
    
      //if it is a surface
      else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
	OCC_surfaces.append( surf );
   
      //if it is a Body
      else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
	OCC_bodies.append( body );

      else
	handled = CUBIT_FALSE;

      if( handled == CUBIT_TRUE )
	{
	  ref_entities_handled.append( ref_entity_ptr );
	  ref_entity_list.change_to(NULL);
	}

      ref_entity_list.step();
    }

  ref_entity_list.remove_all_with_value(NULL);

  int free_body_count = OCC_bodies.size();
  int free_curve_count = OCC_curves.size();
  int free_point_count = OCC_points.size();
  int free_surface_count = OCC_surfaces.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_surface_count == 0 && 
      free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //save the facets (geometry info )
  CubitStatus status;

  //write out topology and attributes
  status = write_topology( p_buffer, n_buffer_size,
			   b_export_buffer,
			   OCC_bodies, OCC_surfaces,
			   OCC_curves, OCC_points);
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_surface_count || 
      free_curve_count || free_point_count )
    if (b_export_buffer) PRINT_INFO( "\nExported:" );
    else PRINT_INFO( "\nSize checked:" );

  int flg = 0;

  if( free_body_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_body_count == 1 )
         PRINT_INFO( "%4d OCC Body to buffer\n", free_body_count );
      else
         PRINT_INFO( "%4d OCC Bodies to buffer\n", free_body_count );
    }

  if( free_surface_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_surface_count == 1 )
	PRINT_INFO( "%4d OCC Surface to buffer\n", free_surface_count );
      else
	PRINT_INFO( "%4d OCC Surface to buffer\n", free_surface_count );
    }

  if( free_curve_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
	PRINT_INFO( "%4d OCC Curve to buffer\n", free_curve_count );
      else
	PRINT_INFO( "%4d OCC Curves to buffer\n", free_curve_count );
    }

  if( free_point_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
	PRINT_INFO( "%4d OCC Point to buffer\n", free_point_count );
      else
	PRINT_INFO( "%4d OCC Points to buffer\n", free_point_count );
    }
  PRINT_INFO( "\n" );

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:export_solid_model
//Member Type:  PUBLIC
//Description:  function called for write out temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus
OCCQueryEngine::write_topology( const char* file_name,
                                const char* file_type,
                                DLIList<OCCBody*> &OCC_bodies,
                                DLIList<OCCSurface*> &OCC_surfaces,
                                DLIList<OCCCurve*> &OCC_curves,
                                DLIList<OCCPoint*> &OCC_points )
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  for (i = 0; i < OCC_bodies.size(); i++)
    {
      OCCBody* body = OCC_bodies.get_and_step();
      TopoDS_CompSolid *shape = body->get_TopoDS_Shape();

      if (shape == NULL) //sheet body or shell body
      {
         OCCSurface* surface = body->my_sheet_surface();
         OCCShell* shell = body->shell();
         if (surface == NULL && shell == NULL)
         {
           
	   PRINT_ERROR( "Wrong body structure. Internal ERROR\n" );
	   continue;
         }
         if(surface)
           B.Add(Co,*(surface->get_TopoDS_Face())); 
         else
    	   B.Add(Co,*(shell->get_TopoDS_Shell()));
         continue;
      }

      //check if this body is build backwards from lump. if so,
      //the body and its CompSolid doesn't have bounded relationship
      //established. In this case, each individual lump of the body 
      // will be exported as TopoDS_Solid without a CompSolid
      if(OCCMap->IsBound(*shape))
	B.Add(Co, *shape);
      else
	{   
	  DLIList<Lump*> lumps = body->lumps();
	  for(int i = 0; i < lumps.size(); i++)
	    {
	      OCCLump *occ_lump = (OCCLump *) lumps.get_and_step();
              TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
	      B.Add(Co, solid);
	    }
	}
    }

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  if(strcmp(file_type, "OCC") == 0)
  {
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    if(!Write(Co, const_cast<char*>(file_name),label))
      return CUBIT_FAILURE;
  } 
  else if(strcmp(file_type, "STEP") == 0)
  {
    STEPControl_Writer writer;
    writer.Model( Standard_True);
    writer.Transfer(Co, STEPControl_AsIs );
    IFSelect_ReturnStatus stat = writer.Write( (char*) file_name);
    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    }
  }

  else // IGES file
  {
    IGESControl_Writer writer;
    writer.AddShape(Co);
    writer.ComputeModel();
    Standard_Boolean  stat = writer.Write( (char*) file_name);
    if (!stat )
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}
/*
CubitStatus
OCCQueryEngine::write_topology( //const unsigned char* p_buffer,
			       //char* p_buffer,
			       ofstream& os,
				DLIList<OCCBody*> &OCC_bodies,
				DLIList<OCCSurface*> &OCC_surfaces,
				DLIList<OCCCurve*> &OCC_curves,
				DLIList<OCCPoint*> &OCC_points)
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  for (i = 0; i < OCC_bodies.size(); i++)
    {
      OCCBody* body = OCC_bodies.get_and_step();
      TopoDS_CompSolid *shape = body->get_TopoDS_Shape();

      if (shape == NULL) //sheet body or shell body
      {
         OCCSurface* surface = body->my_sheet_surface();
         OCCShell* shell = body->shell();
         if (surface == NULL && shell == NULL)
         {
           
	   PRINT_ERROR( "Wrong body structure. Internal ERROR\n" );
	   continue;
         }
         if(surface)
           B.Add(Co,*(surface->get_TopoDS_Face())); 
         else
    	   B.Add(Co,*(shell->get_TopoDS_Shell()));
         continue;
      }

      //check if this body is build backwards from lump. if so,
      //the body and its CompSolid doesn't have bounded relationship
      //established. In this case, each individual lump of the body 
      // will be exported as TopoDS_Solid without a CompSolid
      if(OCCMap->IsBound(*shape))
	B.Add(Co, *shape);
      else
	{   
	  DLIList<Lump*> lumps = body->lumps();
	  for(int i = 0; i < lumps.size(); i++)
	    {
	      OCCLump *occ_lump = (OCCLump *) lumps.get_and_step();
              TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
	      B.Add(Co, solid);
	    }
	}
    }

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  //if(strcmp(file_type, "OCC") == 0)
  //{
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    if(!Write(Co, os))
      return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::write_topology(
			       std::ostringstream& os,
				DLIList<OCCBody*> &OCC_bodies,
				DLIList<OCCSurface*> &OCC_surfaces,
				DLIList<OCCCurve*> &OCC_curves,
				DLIList<OCCPoint*> &OCC_points)
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  for (i = 0; i < OCC_bodies.size(); i++)
    {
      OCCBody* body = OCC_bodies.get_and_step();
      TopoDS_CompSolid *shape = body->get_TopoDS_Shape();

      if (shape == NULL) //sheet body or shell body
      {
         OCCSurface* surface = body->my_sheet_surface();
         OCCShell* shell = body->shell();
         if (surface == NULL && shell == NULL)
         {
           
	   PRINT_ERROR( "Wrong body structure. Internal ERROR\n" );
	   continue;
         }
         if(surface)
           B.Add(Co,*(surface->get_TopoDS_Face())); 
         else
    	   B.Add(Co,*(shell->get_TopoDS_Shell()));
         continue;
      }

      //check if this body is build backwards from lump. if so,
      //the body and its CompSolid doesn't have bounded relationship
      //established. In this case, each individual lump of the body 
      // will be exported as TopoDS_Solid without a CompSolid
      if(OCCMap->IsBound(*shape))
	B.Add(Co, *shape);
      else
	{   
	  DLIList<Lump*> lumps = body->lumps();
	  for(int i = 0; i < lumps.size(); i++)
	    {
	      OCCLump *occ_lump = (OCCLump *) lumps.get_and_step();
              TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
	      B.Add(Co, solid);
	    }
	}
    }

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  //if(strcmp(file_type, "OCC") == 0)
  //{
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    if(!Write(Co, os))
      return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}
*/
CubitStatus
OCCQueryEngine::write_topology( char*& p_buffer,
				int& n_buffer_size,
				bool b_export_buffer,
				DLIList<OCCBody*> &OCC_bodies,
				DLIList<OCCSurface*> &OCC_surfaces,
				DLIList<OCCCurve*> &OCC_curves,
				DLIList<OCCPoint*> &OCC_points)
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  for (i = 0; i < OCC_bodies.size(); i++)
    {
      OCCBody* body = OCC_bodies.get_and_step();
      TopoDS_CompSolid *shape = body->get_TopoDS_Shape();

      if (shape == NULL) //sheet body or shell body
      {
         OCCSurface* surface = body->my_sheet_surface();
         OCCShell* shell = body->shell();
         if (surface == NULL && shell == NULL)
         {
           
	   PRINT_ERROR( "Wrong body structure. Internal ERROR\n" );
	   continue;
         }
         if(surface)
           B.Add(Co,*(surface->get_TopoDS_Face())); 
         else
    	   B.Add(Co,*(shell->get_TopoDS_Shell()));
         continue;
      }

      //check if this body is build backwards from lump. if so,
      //the body and its CompSolid doesn't have bounded relationship
      //established. In this case, each individual lump of the body 
      // will be exported as TopoDS_Solid without a CompSolid
      if(OCCMap->IsBound(*shape))
	B.Add(Co, *shape);
      else
	{   
	  DLIList<Lump*> lumps = body->lumps();
	  for(int i = 0; i < lumps.size(); i++)
	    {
	      OCCLump *occ_lump = (OCCLump *) lumps.get_and_step();
              TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
	      B.Add(Co, solid);
	    }
	}
    }

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  //if(strcmp(file_type, "OCC") == 0)
  //{
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    if(!Write(Co, p_buffer, n_buffer_size, b_export_buffer))
      return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
                                   const Standard_CString File,
                                   TDF_Label label) 
{
  ofstream os;
  os.open(File, ios::out);
  if (!os.rdbuf()->is_open()) return Standard_False;
  
  CubitBoolean isGood = (os.good() && !os.eof());
  if(!isGood)
    return isGood;

  OCCShapeAttributeSet SS;
  SS.Add(Sh);

  os << "DBRep_DrawableShape\n";  // for easy Draw read
  SS.Write(os);
  isGood = os.good();
  if(isGood )
    SS.Write(Sh,os,&label);
  os.flush();
  isGood = os.good();
  os.close();
  isGood = os.good() && isGood;

  return isGood;
}
/*
CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
				   //const unsigned char* pbuffer,
				   //streambuf* rpsbuf)
				   //char* pbuffer,
				   //TDF_Label label) 
				   std::ofstream& os)
{
  os.open("test.txt", std::ios::out);
  std::filebuf* pbuf_os = os.rdbuf();
  if (!pbuf_os->is_open()) return Standard_False;
  
  CubitBoolean isGood = (os.good() && !os.eof());
  if (!isGood) return isGood;

  long os_size1 = pbuf_os->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "size os1 = " << os_size1 << std::endl;
  
  BinTools_ShapeSet SS;
  SS.Add(Sh);
  isGood = os.good();
  SS.Write(os);
  isGood = os.good();
 
  //std::streamsize size = pbuf->in_avail();
  long os_size = pbuf_os->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "size os2 = " << os_size << std::endl;
  
  char* buffer = new char[os_size];
  //os.write(buffer, os_size);
  //pbuf_os->sputn(buffer, os_size);

  //memcpy(buffer, pbuf_os, os_size);
  
  return CUBIT_TRUE;
}

CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
				   std::ostringstream& os)
{

  CubitBoolean isGood = (os.good() && !os.eof());
  if (!isGood) return isGood;

  std::stringbuf* pbuf = os.rdbuf();
  std::streamsize ssize = pbuf->in_avail(); 
  long size = pbuf->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "size os1 = " << size << " " << ssize << std::endl;
  
  BinTools_ShapeSet SS;
  SS.Add(Sh);
  isGood = os.good();
  SS.Write(os);
  isGood = os.good();
 
  ssize = pbuf->in_avail();
  size = pbuf->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "size os2 = " << size << " " << ssize << std::endl;

  char* buffer = new char[size];
  //os.write(buffer, os_size);
  pbuf->sputn(buffer, size);

  std::istringstream is;
  pbuf = is.rdbuf();

  ssize = pbuf->in_avail();
  size = pbuf->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "size io1 = " << size << " " << ssize << std::endl;

  //is.tie(&os);
  //is.read(buffer, size1);
  pbuf->pubseekpos(0, std::ios::in);
  pbuf->sgetn(buffer, size);

  ssize = pbuf->in_avail();
  size = is.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "size io2 = " << size << " " << ssize << std::endl;

  pbuf->pubsetbuf(buffer, size);

  ssize = pbuf->in_avail();
  size = is.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "size io3 = " << size << " " << ssize << std::endl;

  BinTools_ShapeSet SS2;
  SS2.Read(is);
  int nShape = SS2.NbShapes();
  std::cout << "nShpae = %d" << nShape << std::endl;

  
  // test
  std::ostringstream oss("test");
  pbuf = oss.rdbuf();
  ssize = pbuf->in_avail();
  size = pbuf->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "test size os1 = " << size << " " << ssize << std::endl;

  char* tbuffer = new char[size];
  //pbuf->sputn(tbuffer, size);
  ssize = pbuf->sgetn(tbuffer, size);

  std::istringstream iss;
  pbuf = iss.rdbuf();

  //ssize = pbuf->in_avail();
  size = pbuf->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "test size io1 = " << size << " " << ssize << std::endl;

  //is.tie(&os);
  //is.read(buffer, size1);
  pbuf->pubseekpos(0, std::ios::in);
  //pbuf->sgetn(tbuffer, size);
  //ssize = pbuf->sputn(tbuffer, size);
  pbuf->pubsetbuf(tbuffer, size);

  //ssize = pbuf->in_avail();
  size = iss.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "test size io2 = " << size << " " << ssize << std::endl;

  pbuf->pubsetbuf(tbuffer, size);

  ssize = pbuf->in_avail();
  size = iss.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "test size io3 = " << size << " " << ssize << std::endl;

  // 
  std::stringstream ssos("test2");
  ssize = ssos.rdbuf()->in_avail();
  size = ssos.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "test size os4 = " << size << " " << ssize << std::endl;
  std::string mystr = ssos.rdbuf()->str();
  std::stringstream ssis(mystr);
  ssize = ssis.rdbuf()->in_avail();
  size = ssis.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::in);
  std::cout << "test size is4 = " << size << " " << ssize << std::endl;

  std::stringbuf sbos;
  std::stringbuf sbis;
  std::iostream tos(&sbos);
  std::iostream tis(&sbis);
  BinTools_ShapeSet tSS1;
  tSS1.Add(Sh);
  isGood = tos.good();
  std::cout << "isGood=" << isGood << std::endl;
  tSS1.Write(tos);
  isGood = tos.good();
  std::cout << "isGood=" << isGood << std::endl;
  size = tos.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::out);
  std::cout << "test fos size = " << size << std::endl;
  char* tbuffer1 = new char[size];
  tos.read(tbuffer1, size);

  tis.write(tbuffer1, size);
  BinTools_ShapeSet tSS2;
  tSS2.Read(tis);
  isGood = tis.good();
  std::cout << "isGood=" << isGood << std::endl;
  nShape = tSS2.NbShapes();
  std::cout << "test nShpae = " << nShape << std::endl;

  return CUBIT_TRUE;
}
*/
CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
				   char*& pBuffer,
				   int& n_buffer_size,
				   bool b_write_buffer)
{
  std::stringbuf sb;
  std::iostream os(&sb);
  BinTools_ShapeSet SS;

  SS.Add(Sh);
  CubitBoolean isGood = os.good();
  //std::cout << "isGood=" << isGood << std::endl;
  if (!isGood) return isGood;
  SS.Write(os);
  isGood = os.good();
  //std::cout << "isGood=" << isGood << std::endl;
  if (!isGood) return isGood;

  n_buffer_size = os.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::out);
  //std::cout << "test fos size = " << n_buffer_size 
  //  << " nShape= " << SS.NbShapes() << std::endl;

  if (b_write_buffer) {
    //delete pBuffer;
    //pBuffer = new char[n_buffer_size];
    os.read(pBuffer, n_buffer_size);
  }

  return CUBIT_TRUE;
}
                                   
CubitBoolean OCCQueryEngine::Read(TopoDS_Shape& Sh,
                                  const Standard_CString File,
                                  TDF_Label label,
                                  CubitBoolean print_results)
{
  ifstream in( File );
  if (!in) {
    if (print_results) 
      PRINT_INFO("%s: Cannot open file", File );
    return CUBIT_FAILURE;
  }

  OCCShapeAttributeSet SS;
  SS.Read(in, print_results);
  int nbshapes = SS.NbShapes();
  if(!nbshapes) return CUBIT_FALSE;
  SS.Read(Sh,in,nbshapes, &label);
  return CUBIT_TRUE;
}
                                   
CubitBoolean OCCQueryEngine::Read(TopoDS_Shape& Sh,
				  const char* pBuffer,
				  const int n_buffer_size)
{
  std::stringbuf sb;
  std::iostream is(&sb);
  is.write(pBuffer, n_buffer_size);
  /*
  BinTools_ShapeSet SS;
  SS.Read(is);
  CubitBoolean isGood = is.good();
  std::cout << "isGood=" << isGood << std::endl;
  int nShape = SS.NbShapes();
  std::cout << "test nShpae = " << nShape << std::endl;
  if (!nShape) return CUBIT_FALSE;
  //Sh = SS.TopoDS_Shape();
  //SS.Read(Sh, is, nShape);
  */
  //TDF_Label label;
  //CubitBoolean print_results = false;
  //OCCBinToolsShapeSet SS;
  BinTools_ShapeSet SS;
  SS.Read(is);
  int nbshapes = SS.NbShapes();
  if (!nbshapes) return CUBIT_FALSE;
  
  //TopoDS_Shape fs = SS.Shape(1);
  //TopoDS_Shape ls = SS.Shape(nbshapes);
  /* for (int i = 1; i <= nbshapes; i++) {
    TopoDS_Shape s = SS.Shape(i);
    std::cout << "shape type" << i << "=" << s.ShapeType() << std::endl;
    }*/
  
  //SS.Read(Sh, is, nbshapes);
  Sh = SS.Shape(nbshapes);
  
  //is.seekg (0, std::ios::beg);
  //do {
  //SS.Read(Sh, is, nbshapes);
    //if (!Sh.IsNull()) {
      //AddShapes(S,SS);
      //break;
    //}
  //} while(!Sh.IsNull());
  
  /*
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);
  for (int i = 1; i <= nShape; i++) {
    B.Add(Co, SS.Shape(i));
  }

  Sh = Co;
  */
  return CUBIT_TRUE;
}

CubitStatus
OCCQueryEngine::import_temp_geom_file(FILE* file_ptr,
                                      const char* file_name,
                                      const char* file_type,
                                      DLIList<TopologyBridge*> &bridge_list )
{
  //make sure that file_type == "OCC"
  if( !strcmp( file_type,"OCC") || !strcmp( file_type,"IGES") ||
      !strcmp( file_type,"STEP") )
    return import_solid_model( file_name, file_type, bridge_list );
  else
    return CUBIT_FAILURE;
}

//===========================================================================
//Function Name:import_solid_model
//Member Type:  PUBLIC
//Description:  function called for read in temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus OCCQueryEngine::import_solid_model(
					       const char* file_name ,
					       const char* file_type,
					       DLIList<TopologyBridge*> &imported_entities,
					       CubitBoolean print_results,
					       const char* logfile_name,
					       CubitBoolean heal_step,
					       CubitBoolean import_bodies,
					       CubitBoolean import_surfaces,
					       CubitBoolean import_curves,
					       CubitBoolean import_vertices,
					       CubitBoolean free_surfaces)
{
  TopoDS_Shape *aShape = new TopoDS_Shape;
  //BRep_Builder aBuilder;
  if(strcmp(file_type ,"OCC") == 0)
  {
    Standard_Boolean result = Read(*aShape, (char*) file_name, mainLabel, print_results);
    if (result==0) return CUBIT_FAILURE;
  }
 
  else if (strcmp(file_type, "STEP") == 0)
  {
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile( (char*) file_name);
    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    } 
    reader.TransferRoots();
    *aShape = reader.OneShape(); 
  }

  else if(strcmp(file_type, "IGES") == 0)
  {
    IGESControl_Reader reader;
    const Standard_CString string1 = file_name;
    IFSelect_ReturnStatus stat = reader.ReadFile( string1);
    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    } 
    reader.TransferRoots(); 
    *aShape = reader.OneShape();
  } 
  imported_entities = populate_topology_bridge(*aShape);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::import_solid_model(DLIList<TopologyBridge*> &imported_entities,
					       const char* pBuffer,
					       const int n_buffer_size)
{
  TopoDS_Shape *aShape = new TopoDS_Shape;
  Standard_Boolean result = Read(*aShape, pBuffer, n_buffer_size);
  if (result==0) return CUBIT_FAILURE;
  
  imported_entities = populate_topology_bridge(*aShape);
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:populate_topology_bridge
//Member Type:  PUBLIC
//Description:  function called for populating topology bridge for OCC entity 
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

DLIList<TopologyBridge*> OCCQueryEngine::populate_topology_bridge(TopoDS_Shape& aShape)
{
  DLIList<TopologyBridge*> tblist;
  // suitable to populate for a TopoDS_CompSolid or TopoDS_Compound shape.
  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_COMPSOLID); Ex.More(); Ex.Next())
    tblist.append(populate_topology_bridge(TopoDS::CompSolid(Ex.Current())));

  for (Ex.Init(aShape, TopAbs_SOLID, TopAbs_COMPSOLID); Ex.More(); Ex.Next())
  {
    Lump *lump = 
    populate_topology_bridge(TopoDS::Solid(Ex.Current()), CUBIT_TRUE);
    tblist.append(CAST_TO(lump, OCCLump)->get_body());
  }

  for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
  {
    OCCShell* shell =
       populate_topology_bridge(TopoDS::Shell(Ex.Current()), CUBIT_TRUE);
    tblist.append(shell->my_body());
  }

  for (Ex.Init(aShape, TopAbs_FACE, TopAbs_SHELL); Ex.More(); Ex.Next())
  {
    Surface* face =
      populate_topology_bridge(TopoDS::Face(Ex.Current()),CUBIT_TRUE);
    if(face)
      tblist.append(CAST_TO(face, OCCSurface)->my_body());
  }

  for (Ex.Init(aShape, TopAbs_WIRE, TopAbs_FACE); Ex.More(); Ex.Next())
    populate_topology_bridge(TopoDS::Wire(Ex.Current()), CUBIT_TRUE);

  for (Ex.Init(aShape, TopAbs_EDGE, TopAbs_WIRE); Ex.More(); Ex.Next())
    tblist.append(populate_topology_bridge(TopoDS::Edge(Ex.Current())));

  for (Ex.Init(aShape, TopAbs_VERTEX, TopAbs_EDGE); Ex.More(); Ex.Next())
    tblist.append(populate_topology_bridge(TopoDS::Vertex(Ex.Current())));
  tblist.remove_all_with_value(NULL);
  return tblist;
}

BodySM* OCCQueryEngine::populate_topology_bridge(const TopoDS_CompSolid& aShape)
{
  if(aShape.IsNull())
    return (BodySM*)NULL;
  OCCBody *body;
  if (!OCCMap->IsBound(aShape))
    {
      TopoDS_CompSolid *posolid =  new TopoDS_CompSolid;
      *posolid = aShape;
      (iTotalTBCreated)++;
      body = new OCCBody(posolid);
      OCCMap->Bind(*posolid, iTotalTBCreated);
      OccToCGM->insert(valType(iTotalTBCreated,
			       (TopologyBridge*)body));
      BodyList->append(body);
    }
  else
    {
      int k = OCCMap->Find(aShape);
      body = (OCCBody*)(OccToCGM->find(k))->second;
      body->set_TopoDS_Shape(aShape);
    }

  TopExp_Explorer Ex;
  DLIList<Lump*> lumps;
  for (Ex.Init(aShape, TopAbs_SOLID); Ex.More(); Ex.Next())
  {
     Lump* lump = populate_topology_bridge(TopoDS::Solid(Ex.Current()));
     lumps.append(lump);
     CAST_TO(lump, OCCLump)->add_body(body);
  }
  body->lumps(lumps);

  return body;
}

Lump* OCCQueryEngine::populate_topology_bridge(const TopoDS_Solid& aShape,
		 			       CubitBoolean build_body)
{
  if(aShape.IsNull())
    return (Lump*)NULL;
  //one OCCBody corresponds one OCCLump
  OCCLump *lump;
  OCCBody *body;
  if (!OCCMap->IsBound(aShape))
  {
    TopoDS_Solid *posolid =  new TopoDS_Solid;
    *posolid = aShape;
    iTotalTBCreated++;
    lump = new OCCLump(posolid);
    if (build_body)
    {
      DLIList<Lump*> lumps;
      lumps.append(lump);
      body = new OCCBody(lumps);
      BodyList->append(body);
      lump->add_body(body);
    }
    OCCMap->Bind(*posolid, iTotalTBCreated);
    OccToCGM->insert(valType(iTotalTBCreated,
		       (TopologyBridge*)lump));
  }
  else 
  {
    int k = OCCMap->Find(aShape);
    lump = (OCCLump*)(OccToCGM->find(k))->second;
    lump->set_TopoDS_Solid(aShape);
    OCCBody* body = CAST_TO(lump->get_body(), OCCBody);
    DLIList<Lump*> lumps = body->lumps();
    TopoDS_CompSolid* new_top = body->make_CompSolid(lumps);
    body->set_TopoDS_Shape(*new_top);
    new_top->Nullify();
    delete new_top;
  }
 
  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_SHELL); Ex.More(); Ex.Next())
    populate_topology_bridge(TopoDS::Shell(Ex.Current()));

  return lump;
}

OCCShell* OCCQueryEngine::populate_topology_bridge(const TopoDS_Shell& aShape,
						   CubitBoolean standalone)
{
  if(aShape.IsNull())
    return (OCCShell*)NULL;
  OCCShell *shell ;
  DLIList<OCCCoFace*> cofaces_old, cofaces_new;
  if (!OCCMap->IsBound(aShape))
  {
    if(standalone)
    {
      //check if just has one Face,if so, don't make new shell.
      TopExp_Explorer Ex;
      int num_faces = 0;
      TopoDS_Face topo_face;
      for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next())
      {
        topo_face = TopoDS::Face(Ex.Current());
        num_faces++;
      }
      if(num_faces == 1)
      {
        Surface* face = populate_topology_bridge(topo_face, standalone);
        return CAST_TO(face, OCCSurface)->my_shell();
      }
    }
    TopoDS_Shell *poshell = new TopoDS_Shell;
    *poshell = aShape;
    iTotalTBCreated++;
    shell = new OCCShell(poshell);
    OCCMap->Bind(*poshell, iTotalTBCreated);
    OccToCGM->insert(valType(iTotalTBCreated,
		       (TopologyBridge*)shell));
    shell->set_body(NULL);
    shell->set_lump(NULL); 
    if(standalone)
    {
      OCCLump* lump = new OCCLump(NULL, NULL, shell);
      OCCBody* body = new OCCBody(NULL, CUBIT_FALSE, NULL, shell);
      shell->set_body(body);
      shell->set_lump(lump);
      //don't need to add body into BodyList.
    }
  }
  else
  {
    int k = OCCMap->Find(aShape);
    shell = (OCCShell*)(OccToCGM->find(k))->second;
    cofaces_old =  shell->cofaces();
    shell->set_TopoDS_Shell(aShape);
  }

  TopExp_Explorer Ex;
  DLIList<OCCCoFace*> cofaces;
  for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next())
  {
    TopoDS_Face topo_face = TopoDS::Face(Ex.Current());
    Surface* face =
      populate_topology_bridge(topo_face, CUBIT_FALSE);
    
    if(!face)
      continue;
    OCCSurface *occ_surface = CAST_TO(face, OCCSurface);
    //check if surface was a sheet surface, delete it if so
    if(occ_surface->my_shell() != NULL && occ_surface->my_shell() != shell)
    {
       //if Sheet_body, delete this sheet body
       if(occ_surface->my_body() != NULL)
       {
          delete_solid_model_entities(occ_surface->my_body());
          face =
            populate_topology_bridge(topo_face, CUBIT_FALSE);
   
          if(!face)
            continue;
          occ_surface = CAST_TO(face, OCCSurface);
       }
    }
    CubitBoolean exist = CUBIT_FALSE;
    OCCCoFace * coface = NULL;
    int size = cofaces_old.size();
    CubitSense sense ;
/*    if( aShape.Orientation() == TopAbs_REVERSED )
      sense = (topo_face.Orientation() == TopAbs_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
    else*/
      sense = (topo_face.Orientation() == TopAbs_FORWARD ? CUBIT_FORWARD : CUBIT_REVERSED); 
    for(int i = 0; i < size; i++)
    {
      coface = cofaces_old.get_and_step();
      if(coface->surface() == occ_surface)
      {
        exist = CUBIT_TRUE;
        coface->set_sense(sense);
        cofaces_new.append(coface);
        break;
      }
    }
    if(!exist)
    {
      DLIList<OCCShell*> shells = occ_surface->shells();
      for(int i = 0; i < shells.size() ; i++)
      {
        OCCShell* shell = shells.get_and_step();
        TopoDS_Shell* topo_shell = shell->get_TopoDS_Shell();
        if (!OCCMap->IsBound(*topo_shell))
        {
          DLIList<OCCCoFace*> coface_list = shell->cofaces();
          for(int j = 0; j < coface_list.size(); j++)
          {
            OCCCoFace * test_coface = coface_list.get_and_step();
            occ_surface->remove_shell(shell);
            shell->remove_coface(test_coface);
            delete test_coface;
          }
          if(!topo_shell->IsNull())
            topo_shell->Nullify();
        }
      } 
      OCCCoFace * coface = new OCCCoFace( occ_surface, shell, sense);
      cofaces_new.append(coface);
      occ_surface->add_shell(shell);
    }
 
    if(standalone)
      occ_surface->set_shell(shell);
  }
/*  if(aShape.Orientation() == TopAbs_REVERSED)
    cofaces_new.reverse();*/
  shell->cofaces(cofaces_new);
  return shell;
}

Surface* OCCQueryEngine::populate_topology_bridge(const TopoDS_Face& aShape,
                                                  CubitBoolean build_body)
{
  if(aShape.IsNull())
    return (Surface*)NULL;
  OCCSurface *surface = NULL;
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(aShape, myProps);
  double area = myProps.Mass();
  double tol = get_sme_resabs_tolerance();
  if(area < tol * tol)
    return (Surface*) NULL;

  if (!OCCMap->IsBound(aShape))
  {
    TopoDS_Face *poface = new TopoDS_Face;
    *poface = aShape;
    surface = new OCCSurface(poface);

    iTotalTBCreated++;
    OCCMap->Bind(*poface, iTotalTBCreated);
    OccToCGM->insert(valType(iTotalTBCreated,
                             (TopologyBridge*)surface));
    SurfaceList->append(surface);
    surface->set_body(NULL);
    surface->set_lump(NULL);
    surface->set_shell(NULL);
    if(build_body)
    {
      OCCShell* shell = new OCCShell(NULL, surface);
      OCCLump* lump = new OCCLump(NULL, surface);
      OCCBody* body = new OCCBody(NULL, CUBIT_TRUE, surface);
      surface->set_body(body);
      surface->set_lump(lump);
      surface->set_shell(shell);
      shell->set_body(body);
      shell->set_lump(lump);
      //Doesn't need to save sheet bodies.
    }
  } 

  else 
  {
    int k = OCCMap->Find(aShape);
    surface = (OCCSurface*)(OccToCGM->find(k))->second;
    TopoDS_Face aFace(aShape);
    surface->set_TopoDS_Face(aFace);
  }

  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_WIRE); Ex.More(); Ex.Next())
    populate_topology_bridge(TopoDS::Wire(Ex.Current()));

  return surface;
}

OCCLoop* OCCQueryEngine::populate_topology_bridge(const TopoDS_Wire& aShape,
						  CubitBoolean standalone)
{
  if(aShape.IsNull())
    return (OCCLoop*)NULL;
  OCCLoop *loop ;
  if (!OCCMap->IsBound(aShape))
    {
      TopoDS_Wire *powire = new TopoDS_Wire;
      *powire = aShape;
      iTotalTBCreated++;
      loop = new OCCLoop(powire);
      OCCMap->Bind(*powire, iTotalTBCreated);
      OccToCGM->insert(valType(iTotalTBCreated,
			       (TopologyBridge*)loop));
      if(standalone)
	WireList->append(loop);
    }
  else
    {
      int k = OCCMap->Find(aShape);
      loop = (OCCLoop*)(OccToCGM->find(k))->second;
      loop->set_TopoDS_Wire(aShape);
    }

  BRepTools_WireExplorer Ex;
  DLIList <OCCCoEdge*> coedges_old, coedges_new;
  coedges_old = loop->coedges();
  for (Ex.Init(aShape); Ex.More(); Ex.Next())
  {
    Curve* curve = populate_topology_bridge(Ex.Current());
    if(!curve)
      continue;
    OCCCurve *occ_curve = CAST_TO(curve, OCCCurve);
    DLIList<OCCLoop*> loops = occ_curve->loops();
    CubitBoolean exist = CUBIT_FALSE;
    OCCCoEdge * coedge = NULL;
    int size = coedges_old.size();
    CubitSense sense ;
/*    if( aShape.Orientation() == TopAbs_REVERSED )
      sense = (Ex.Current().Orientation() == TopAbs_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
    else*/
      sense = (Ex.Current().Orientation() == TopAbs_FORWARD ? CUBIT_FORWARD : CUBIT_REVERSED);
    for(int i = 0; i < size; i++)
    {
      coedge = coedges_old.get_and_step();
      if(coedge->curve() == curve && coedge->sense() ==  sense)
      {
        coedge->set_mark(1);
        exist = CUBIT_TRUE;
        coedge->set_sense(sense);
        coedges_new.append(coedge);
        break;
      }
    }   
    
    if(!exist)
    {
      //search through the curve loops
      for(int i = 0; i < loops.size() ; i++)
      {
        OCCLoop* occ_loop = loops.get_and_step();
        TopoDS_Wire* wire = occ_loop->get_TopoDS_Wire();
        if (wire->IsNull() || !OCCMap->IsBound(*wire))
        { 
          DLIList<OCCCoEdge*> coedge_list = occ_loop->coedges();
          for(int j = 0; j < coedge_list.size(); j++)
          {
            OCCCoEdge * test_coedge = coedge_list.get_and_step();
            occ_loop->remove_coedge(test_coedge);
            occ_curve->remove_loop(occ_loop);
            delete test_coedge;
          }
          if(!wire->IsNull())
            wire->Nullify();
          WireList->remove(occ_loop);
        }
      }
      //for the cylinder side face, there are 4 coedges, 2 of them are seam
      //edges and should have opposite sense.
      for(int i = 0; i < coedges_new.size(); i++)
      {
        coedge =  coedges_new.get_and_step();
        Curve* test_c = coedge->curve();
        if(test_c == curve)
        {
          if(sense == coedge->sense())
            sense = (sense == CUBIT_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
          break;
        }
      }

      if(occ_curve->loops().size() == 2)
      {
        //there must be a loop which doesn't have this curve anymore
        //this is been found in subtract cases while one solid becomes
        //2 solid, and one face becomes two faces. One face uses/updates
        //the old face while the other face generates face and wire from
        //new. However, in order to uses the curves, the old curve is kept
        //as possible, so curve's loops get kept, but since it's going to 
        //associate with new loop, the old loop should be removed from the 
        //loop list.
        DLIList<OCCLoop*> old_loops = occ_curve->loops();
        for (int i = 0; i < 2; i++)
        {
          OCCLoop* old_loop = old_loops.get_and_step();
          DLIList<OCCCoEdge*> test_coedges = old_loop->coedges();
          int found = 0;
          for(int j = 0; j < test_coedges.size() ; j++)
          {
            if(test_coedges.get()->curve() != curve)
              test_coedges.step();
            else
            {
              found = 1;
              break;
            }
          }
          if(!found)
            occ_curve->remove_loop(old_loop); 
        }  
      }
      //for unite case, it's possible that the a curve has 3 coedges because
      //opencascade do not perform unite on surfaces.
      coedge = new OCCCoEdge( curve, loop, sense);
      coedges_new.append(coedge);
      occ_curve->add_loop(loop);
    }
  }
/*  if(aShape.Orientation() == TopAbs_REVERSED)
    coedges_new.reverse();*/
  loop->coedges(coedges_new);

  //remove unused coedges
  for(int i = 0; i < coedges_old.size(); i++)
  {
    OCCCoEdge *coedge = coedges_old.get_and_step();
    if(coedge->get_mark() != 1)
      delete coedge;
    else
      coedge->set_mark(0);
  }
  return loop;
}

Curve* OCCQueryEngine::populate_topology_bridge(const TopoDS_Edge& aShape)
{
  if(aShape.IsNull())
    return (Curve*)NULL;
  Curve *curve;
  GProp_GProps myProps;
  BRepGProp::LinearProperties(aShape, myProps);
  double length =  myProps.Mass();
  if(length < get_sme_resabs_tolerance())
    return (Curve*) NULL;

  if (!OCCMap->IsBound(aShape)) 
    {
      TopoDS_Edge *poedge = new TopoDS_Edge;
      *poedge = aShape;
      iTotalTBCreated++;
      curve = new OCCCurve(poedge);
      
      OCCMap->Bind(*poedge, iTotalTBCreated);
      OccToCGM->insert(valType(iTotalTBCreated,
			       (TopologyBridge*)curve));
      CurveList->append((OCCCurve*)curve);
    }
  else 
    {
      int i = OCCMap->Find(aShape);
      curve = (OCCCurve*)(OccToCGM->find(i))->second;
      CAST_TO(curve, OCCCurve)->set_TopoDS_Edge(aShape);
    }

  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next())
    populate_topology_bridge(TopoDS::Vertex(Ex.Current()));

  return curve;
}

Point* OCCQueryEngine::populate_topology_bridge(const TopoDS_Vertex& aShape)
{
  if(aShape.IsNull())
    return (Point*)NULL;
  OCCPoint *point;
  if (iTotalTBCreated == 0 || !OCCMap->IsBound(aShape)) 
    {
      TopoDS_Vertex *povertex = new TopoDS_Vertex;
      *povertex = aShape;
      iTotalTBCreated++;
      point = new OCCPoint(povertex);
      OCCMap->Bind(*povertex, iTotalTBCreated);
      OccToCGM->insert(valType(iTotalTBCreated,
			       (TopologyBridge*)point));
    } 
  else 
    {
      int i = OCCMap->Find(aShape);
      point = (OCCPoint*)(OccToCGM->find(i))->second;
      point->set_TopoDS_Vertex(aShape);
    }
  return point;
}

TopologyBridge* OCCQueryEngine::occ_to_cgm(const TopoDS_Shape& shape)
{
  if(!OCCMap->IsBound(shape))
    return (TopologyBridge*) NULL;

  int k = OCCMap->Find(shape);
  return (OccToCGM->find(k))->second;
}	

//-------------------------------------------------------------------------
// Purpose       : Deletes all solid model entities associated with the
//                 Bodies in the input list.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 4/23/01
//-------------------------------------------------------------------------
void OCCQueryEngine::delete_solid_model_entities(DLIList<BodySM*>&bodyList) const
{
  BodySM* BodyPtr = NULL;
  for (int i = 0; i < bodyList.size(); i++ )
    {
      BodyPtr = bodyList.get_and_step();
      this->delete_solid_model_entities(BodyPtr);
    }

  return;
}

CubitStatus OCCQueryEngine::delete_solid_model_entities(
    GeometryEntity* ref_entity_ptr,
    bool remove_lower_entities) const
{
     //Lump
   Lump* lump = CAST_TO(ref_entity_ptr, Lump);
   if(lump != NULL)
   {
     BodySM* body = CAST_TO(lump, OCCLump)->get_body();
     DLIList<Lump*> lumps = CAST_TO(body, OCCBody)->lumps();
 
     if (remove_lower_entities)
       return delete_solid_model_entities(body);

     DLIList<TopologyBridge*> children;
     for(int i = 0; i < lumps.size(); i++)
     {
       lump = lumps.get_and_step();
       CAST_TO(lump, OCCLump)->get_children_virt(children);
     }

     CubitStatus stat = this->unhook_BodySM_from_OCC(body); 
     if(stat)
     {
       while (children.size())
          delete children.pop();
       while(lumps.size())
          delete lumps.pop();
       delete body;
     }
     return stat;
   }

     // Surface
   Surface* ref_face_ptr = CAST_TO(ref_entity_ptr, Surface);
   if (ref_face_ptr != NULL)
   {
     if (remove_lower_entities)
       return ( this->delete_solid_model_entities(ref_face_ptr) );
     CubitStatus stat = this->unhook_Surface_from_OCC(ref_face_ptr);
     if(stat)
       delete ref_face_ptr;
     return stat;
   }

     // Curve
   Curve* ref_edge_ptr = CAST_TO(ref_entity_ptr, Curve);
   if (ref_edge_ptr != NULL)
   {
      if (remove_lower_entities)
        return ( this->delete_solid_model_entities(ref_edge_ptr));
      CubitStatus stat = this->unhook_Curve_from_OCC(ref_edge_ptr);
      if(stat)
        delete ref_edge_ptr;
      return stat;
   }

     // Point
   Point* ref_vertex_ptr = CAST_TO(ref_entity_ptr, Point);
   if (ref_vertex_ptr != NULL)
   {
      return ( this->delete_solid_model_entities(ref_vertex_ptr) );
   }

     // Oops!
   PRINT_ERROR("In AcisQueryEngine::delete_solid_model_entities\n"
               "       Can only delete solid model entities underlying \n"
               "RefFaces, RefEdges and RefVertices.\n");
   return CUBIT_FAILURE;

}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCBody and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( BodySM* bodysm ) const
{
  OCCBody* occ_body = dynamic_cast<OCCBody*>(bodysm);
  if (!occ_body)
    return CUBIT_FAILURE;

  OCCSurface* occ_surface = occ_body->my_sheet_surface();
  if(occ_surface)
  { 
    delete occ_surface->my_body();
    occ_surface->set_body((OCCBody*)NULL);
    delete occ_surface->my_shell();
    delete occ_surface->my_lump();
    return delete_solid_model_entities(occ_surface);
  }

  OCCShell* occ_shell = occ_body->shell();
  if(occ_shell)
  {
    delete occ_shell->my_body();
    occ_shell->set_body((OCCBody*)NULL);
    delete occ_shell->my_lump();
    DLIList<TopologyBridge*> tb_surfaces;
    occ_shell->get_children_virt(tb_surfaces);
    CubitStatus stat = unhook_ShellSM_from_OCC(occ_shell);
    for(int k = 0; k < tb_surfaces.size(); k++)
      delete_solid_model_entities(CAST_TO(tb_surfaces.get_and_step(), Surface));
    delete occ_shell;
    return stat;
  }

  DLIList<TopologyBridge*> children;
  DLIList<Lump*> lumps = occ_body->lumps();
  DLIList<ShellSM*> shell_list;

  for(int i =0; i < lumps.size(); i++)
  {
     Lump* lump = lumps.get_and_step();
     OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     if (!occ_lump)
        continue;
     children.clean_out();
     occ_lump->get_children_virt(children);
     for(int j = 0; j < children.size(); j++)
     {
       ShellSM* shell = CAST_TO(children.get_and_step(), ShellSM);
       
       if (shell)
         shell_list.append(shell);
       DLIList<TopologyBridge*> tb_surfaces;
       shell->get_children_virt(tb_surfaces);
       for(int k = 0; k < tb_surfaces.size(); k++)
         delete_solid_model_entities(CAST_TO(tb_surfaces.get_and_step(), Surface));
     }
  }

  CubitStatus stat = unhook_BodySM_from_OCC(bodysm);

  for(int j = 0; j < shell_list.size(); j++)
     delete shell_list.get_and_step();

  for(int i =0; i < lumps.size(); i++)
     delete lumps.get_and_step(); 

  BodyList->remove(occ_body);
  delete bodysm;
  return stat;
}

CubitStatus
OCCQueryEngine::unhook_BodySM_from_OCC( BodySM* bodysm )const
{
  OCCBody* occ_body = dynamic_cast<OCCBody*>(bodysm);
  if (!occ_body)
    return CUBIT_FAILURE;

  TopoDS_Shape* shape = occ_body->get_TopoDS_Shape();

  if (!shape)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*shape) ;

  //remove the entry from the map
  int k;
  OCCBody* occ_body_find = NULL;
  if(OCCMap->IsBound(*shape))
  {
      k = OCCMap->Find(*shape);

      if(!OCCMap->UnBind(*shape))
        PRINT_ERROR("The OccBody and TopoDS_Shape pair is not in the map!");

      occ_body_find = (OCCBody*)(OccToCGM->find(k))->second;

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccBody and TopoDS_Shape pair is not in the map!");
  }

  DLIList<Lump*> lumps = occ_body->lumps();
  for(int i =0; i < lumps.size(); i++)
  {
     Lump* lump = lumps.get_and_step();
     //OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     //if(occ_lump)
     //  occ_lump->remove_body();

     unhook_Lump_from_OCC(lump);
  }

  if (!shape->IsNull())
    shape->Nullify();
  return CUBIT_SUCCESS;
} 
  
//-------------------------------------------------------------------------
// Purpose       : unhook a OCClumps and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 11/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Lump_from_OCC( Lump* lump ) const
{
  if (lump == NULL)
    return CUBIT_FAILURE;

  OCCLump* occ_lump = dynamic_cast<OCCLump*>(lump);
  if (!occ_lump)
    return CUBIT_FAILURE;

  TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();

  if(!solid)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*solid) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*solid))
    {
      k = OCCMap->Find(*solid);

      if(!OCCMap->UnBind(*solid))
	PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");

      if(!OccToCGM->erase(k))
	PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");
    }
  
  DLIList<TopologyBridge*> children;
  occ_lump->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     ShellSM* shell = CAST_TO(children.get_and_step(), ShellSM); 
     unhook_ShellSM_from_OCC(shell);
  }
  if (occ_lump->get_body() != NULL)
    BodyList->remove(CAST_TO(occ_lump->get_body(), OCCBody));

  if(!solid->IsNull())
     solid->Nullify();
  return CUBIT_SUCCESS;
} 

//-------------------------------------------------------------------------
// Purpose       : unhook a ShellSM from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_ShellSM_from_OCC( ShellSM* shell) const
{
  OCCShell* occ_shell = dynamic_cast<OCCShell*>(shell);
  if (!occ_shell)
    return CUBIT_FAILURE;

  TopoDS_Shell* Shell = occ_shell->get_TopoDS_Shell();

  if(!Shell)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*Shell))
    {
      k = OCCMap->Find(*Shell);

      if(!OCCMap->UnBind(*Shell))
        PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");
    }

  if(!Shell->IsNull())
    Shell->Nullify();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose      : unhook a list of OCCCoFaces from their underlining OCC entity.//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_CoFaces_from_OCC( DLIList<OCCCoFace*>& cofaces) const
{
  int size = cofaces.size();
  while(size > 0)
  {
     OCCCoFace* coface = cofaces.pop();

     OCCShell* shell = coface->shell();
     shell->remove_coface(coface);

     delete coface;

     size = cofaces.size();
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCSurface and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Surface* surface)const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf)
    return CUBIT_FAILURE;

  DLIList<TopologyBridge*> children;
  fsurf->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     LoopSM* loop = CAST_TO(children.get_and_step(), LoopSM);
     delete_loop(loop);
  }
  CubitStatus stat = unhook_Surface_from_OCC(surface);
  if (stat)
    delete surface;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Surface from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Surface_from_OCC( Surface* surface) const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf)
    return CUBIT_FAILURE;

  TopoDS_Face *face = fsurf->get_TopoDS_Face();

  if(!face)
    return CUBIT_FAILURE;

  unhook_cofaces_of_a_surface(fsurf);

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*face) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*face))
    {
      k = OCCMap->Find(*face);

      if(!OCCMap->UnBind(*face))
        PRINT_WARNING("The OccSurface and TopoDS_Face pair is not in the map!");

      if(!OccToCGM->erase(k))
        PRINT_WARNING("The OccSurface and TopoDS_Face pair is not in the map!");
    }
  SurfaceList->remove(fsurf);
  if(!face->IsNull())
    face->Nullify();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCLoop and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_loop( LoopSM* loopsm)const
{
  OCCLoop* occ_loop = dynamic_cast<OCCLoop*>(loopsm);
  if (!occ_loop)
    return CUBIT_FAILURE;

  DLIList<OCCCoEdge*> children;
  children = occ_loop->coedges();
  int size = children.size();
  DLIList<Curve*> curves;
  while(size > 0)
  {
     OCCCoEdge* coedge = children.pop();
     Curve* curve = coedge->curve();
     curves.append_unique(curve);
     size = children.size();
  }
   
  CubitStatus status = unhook_LoopSM_from_OCC(loopsm);

  for(int i = 0; i < curves.size(); i ++)
    delete_solid_model_entities(curves.get_and_step());

  if (status)
  {
    WireList->remove(occ_loop);
    delete loopsm;
  }
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a LoopSM from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_LoopSM_from_OCC( LoopSM* loopsm) const
{
  OCCLoop* occ_loop = dynamic_cast<OCCLoop*>(loopsm);
  if (!occ_loop)
    return CUBIT_FAILURE;

  TopoDS_Wire* wire = occ_loop->get_TopoDS_Wire();

  if(!wire)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*wire))
    {
      k = OCCMap->Find(*wire);

      if(!OCCMap->UnBind(*wire))
        PRINT_ERROR("The OccLoop and TopoDS_Wire pair is not in the map!");

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccLoop and TopoDS_Wire pair is not in the map!");
    }

  if(!wire->IsNull())
    wire->Nullify(); 
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose      : unhook a list of OCCCoEdges from their underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_CoEdges_from_OCC( DLIList<OCCCoEdge*>& coedges) const
{
  int size = coedges.size();
  while(size > 0)
  {
     OCCCoEdge* coedge = coedges.pop();

     LoopSM* loopsm = coedge->loop();
     OCCLoop* loop = CAST_TO(loopsm, OCCLoop);
     assert(loop);
     loop->remove_coedge(coedge);
 
     delete coedge;

     size = coedges.size();
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCCurve and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Curve* curve)const
{
  OCCCurve* fcurve = CAST_TO(curve, OCCCurve);
  if (!fcurve )
    return CUBIT_FAILURE;

  DLIList<TopologyBridge*> children;
  fcurve->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     Point* point = CAST_TO(children.get_and_step(), Point);
     delete_solid_model_entities(point);
  }
  
  CubitStatus stat = unhook_Curve_from_OCC(curve);
  if (!stat)
    return CUBIT_FAILURE;

  CurveList->remove(fcurve);
  delete fcurve;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Curve from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Curve_from_OCC( Curve* curve ) const
{
  OCCCurve* fcurve = dynamic_cast<OCCCurve*>(curve);
  if (!fcurve )
    return CUBIT_FAILURE;

  unhook_coedges_of_a_curve(fcurve);
    
  fcurve->clean_loops();
  TopoDS_Edge* edge = fcurve->get_TopoDS_Edge();
  if (!edge)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*edge) ;
  
  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*edge))
    {
      k = OCCMap->Find(*edge);

      if(!OCCMap->UnBind(*edge))
        PRINT_WARNING("The OccCurve and TopoDS_Edge pair is not in the map!");

      if(!OccToCGM->erase(k))
        PRINT_WARNING("The OccCurve and TopoDS_Edge pair is not in the map!");
    }
  CurveList->remove(fcurve); 
  if(!edge->IsNull())
    edge->Nullify();
  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------------
// Purpose       : Delete a OCCPoint and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Point* point) const
{
  OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);
  if (!fpoint)
    return CUBIT_FAILURE;

  CubitStatus stat = unhook_Point_from_OCC(point);
  if(stat)
    delete point;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Point from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Point_from_OCC( Point* point) const 
{
  OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);
  if (!fpoint)
    return CUBIT_FAILURE;

  TopoDS_Vertex* vertex = fpoint->get_TopoDS_Vertex();
  if (!vertex)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*vertex) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*vertex))
    {
      k = OCCMap->Find(*vertex);

      if(!OCCMap->UnBind(*vertex))
        PRINT_ERROR("The OccPoint and TopoDS_Vertex pair is not in the map!");

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccPoint and TopoDS_Vertex pair is not in the map!");
    }
  if(!vertex->IsNull())
    vertex->Nullify();
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : fire a ray at the specified body, returning the entities hit.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus OCCQueryEngine::fire_ray(BodySM * body,
                                     const CubitVector &start,
                                     const CubitVector &unit,
                                     DLIList<double>& ray_parms,
				     DLIList<GeometryEntity*> *entity_list) const
{
  CubitStatus status = CUBIT_SUCCESS;

  //- fire a ray at the specified body, returning the entities hit and
  //- the parameters along the ray; return CUBIT_FAILURE if error
  // - line body intersection. 
  gp_Pnt p(start.x(), start.y(), start.z());
  gp_Dir dir(unit.x(), unit.y(), unit.z());
  gp_Lin L(p, dir);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(L); 
  
  OCCBody *occBody = CAST_TO(body, OCCBody);
  if (occBody == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  BRepExtrema_DistShapeShape distShapeShape(edge,
					    *(occBody->get_TopoDS_Shape()));
  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input body.\n");
      return CUBIT_FAILURE;
    }

  if (distShapeShape.Value() < get_sme_resabs_tolerance())
    {
      int numPnt = distShapeShape.NbSolution();
      for (int i = 1; i <= numPnt; i++)
	{
	  double para;
	  distShapeShape.ParOnEdgeS1(i , para);
	  ray_parms.append(para);

	  TopoDS_Shape shape = distShapeShape.SupportOnShape2(i);
	  int k = OCCMap->Find(shape);
	  std::map<int,TopologyBridge*>::iterator it = OccToCGM->find(k);
	  TopologyBridge* tb = (*it).second;
	  entity_list->append((GeometryEntity*)tb);
	}
    } 
  return status;
}
double OCCQueryEngine::get_sme_resabs_tolerance() const
{
  return Precision::Confusion(); 
}
// Gets solid modeler's resolution absolute tolerance

double OCCQueryEngine::set_sme_resabs_tolerance( double p)
{
  double old_p = get_sme_resabs_tolerance();
  BRepBuilderAPI::Precision(p);
  return old_p;
}

CubitStatus OCCQueryEngine::set_int_option( const char* , int )
{
  PRINT_ERROR("OCCQueryEngine::set_int_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_dbl_option( const char* , double )
{
  PRINT_ERROR("OCCQueryEngine::set_dbl_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_str_option( const char* , const char* )
{
  PRINT_ERROR("OCCQueryEngine::set_str_option not yet implemented.\n");
  return CUBIT_FAILURE;
}
//- Set solid modeler options


//===========================================================================
//Function Name: ensure_is_ascii_stl_file
//Member Type:
//Description: returns CUBIT_TRUE in is_ascii if fp is an ascii stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================
CubitStatus OCCQueryEngine::ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii)
{

  char line[128]="";

  if (fgets(line, 128, fp)==NULL)
    {
      return CUBIT_FAILURE;
    }
  if (fgets(line, 128, fp)==NULL)
    {
      return CUBIT_FAILURE;
    }
  if (strlen(line)==127)
    {
      if (fgets(line, 128, fp)==NULL)
	{
	  return CUBIT_FAILURE;
	}
    }


  unsigned int dummy_int=0;

  while (isspace(line[dummy_int])&& dummy_int<strlen(line)) dummy_int++;

  if (strlen(line)-dummy_int>5)
    {
      if (tolower(line[dummy_int++])=='f' &&
	  tolower(line[dummy_int++])=='a' &&
	  tolower(line[dummy_int++])=='c' &&
	  tolower(line[dummy_int++])=='e' &&
	  tolower(line[dummy_int])=='t')
	{
	  if (fgets(line, 128, fp)==NULL)
	    {
	      return CUBIT_FAILURE;
	    }
	  dummy_int=0;
	  while (isspace(line[dummy_int])&& dummy_int<strlen(line))
	    {
	      dummy_int++;
	    }
	  if (strlen(line)-dummy_int>5)
	    {
	      if (tolower(line[dummy_int++])=='o' &&
		  tolower(line[dummy_int++])=='u' &&
		  tolower(line[dummy_int++])=='t' &&
		  tolower(line[dummy_int++])=='e' &&
		  tolower(line[dummy_int])=='r')
		{
		  if (fgets(line, 128, fp)==NULL)
		    {
		      return CUBIT_FAILURE;
		    }
		  dummy_int=0;
		  while (isspace(line[dummy_int])&& dummy_int<strlen(line)) {
		    dummy_int++;
		  }
		  if (strlen(line)-dummy_int>6)
		    {
		      if (tolower(line[dummy_int++])=='v' &&
			  tolower(line[dummy_int++])=='e' &&
			  tolower(line[dummy_int++])=='r' &&
			  tolower(line[dummy_int++])=='t' &&
			  tolower(line[dummy_int++])=='e'	&&
			  tolower(line[dummy_int])=='x')
			{
			  is_ascii=CUBIT_TRUE;
			}
		    }
		}
	    }
	}
    }
  return CUBIT_SUCCESS;
}


//=============================================================================
//Function:   create_super_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a list of BodySMs
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus OCCQueryEngine::create_super_bounding_box(
						      DLIList<BodySM*>& body_list,
						      CubitBox& super_box )
{
  BodySM *bodySM;
  int i;
  CubitStatus status = CUBIT_SUCCESS;

  body_list.reset();
  for ( i = 0; i < body_list.size(); i++ ) {
    bodySM = body_list.get_and_step();  
    OCCBody* occBody = CAST_TO(bodySM, OCCBody);
    super_box |= occBody->get_bounding_box();
  }

  return status;
}

CubitStatus OCCQueryEngine::restore_transform( BodySM* body )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  OCCBody* theBody = dynamic_cast<OCCBody*>(body);
  return theBody ? theBody->move( d.x(), d.y(), d.z() ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::rotate( BodySM* body, const CubitVector& v, double a )
{
  // a is in degree.
  OCCBody* occ_bod = dynamic_cast<OCCBody*>(body);
  a *= CUBIT_PI/180;
  return occ_bod ? occ_bod->rotate( v.x(), v.y(), v.z(), a ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, double factor )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( factor ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, const CubitVector& f )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( f.x(), f.y(), f.z() ) : CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Solid, Surface, Curve, or Vertex
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/23/07
//-------------------------------------------------------------------------
CubitStatus OCCQueryEngine::translate( GeometryEntity* entity,
                                       const CubitVector& v )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Vec aVec(v.x(), v.y(),v.z());
  gp_Trsf aTrsf;  
  aTrsf.SetTranslation(aVec);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::update_entity_shape(GeometryEntity* entity_ptr,
					BRepBuilderAPI_Transform* aBRepTrsf,
                                        BRepAlgoAPI_BooleanOperation *op)
{
  if (OCCBody *body_ptr = CAST_TO( entity_ptr, OCCBody))
    {
      body_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  else if( OCCSurface *surface_ptr = CAST_TO( entity_ptr, OCCSurface))
    {
      surface_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  else if( OCCCurve *curve_ptr = CAST_TO( entity_ptr, OCCCurve))
    {
       curve_ptr->update_OCC_entity(aBRepTrsf, op); 
       return CUBIT_SUCCESS;
    }

  else if( OCCPoint *point_ptr = CAST_TO( entity_ptr, OCCPoint))
    {
      point_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  PRINT_ERROR("Non-OCC TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return CUBIT_FAILURE;
}

//a is angular value of rotation in radians
CubitStatus OCCQueryEngine::rotate( GeometryEntity* entity,
                                    const CubitVector& v, double a )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax1 anAxis(aOrigin, aDir);

  //a is angular value of rotation in radians 
  gp_Trsf aTrsf;
  aTrsf.SetRotation(anAxis, a);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* entity, double f )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Trsf aTrsf;
  aTrsf.SetScaleFactor(f);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR("non-uniform scaling is not performed on OCC bodies");
  return CUBIT_FAILURE;
}

// like ACIS, here v is the normal of symmetric plane.
CubitStatus OCCQueryEngine::reflect( GeometryEntity* entity, 
                                     const CubitVector&  v)
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }
 
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax2 anAx2(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAx2);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : bodies_overlap
// Member Type: PUBLIC
// Description: determine if OCC-based bodies overlap
// Author     : Jane Hu 
// Date       : 10/07
//===============================================================================
CubitBoolean OCCQueryEngine::bodies_overlap (BodySM * body_ptr_1,
                                             BodySM * body_ptr_2 ) const
{
  OCCBody *occ_body1 = CAST_TO(body_ptr_1, OCCBody);
  if (!occ_body1)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
      return CUBIT_FALSE;
    }
 
  OCCBody *occ_body2 = CAST_TO(body_ptr_2, OCCBody);
  if (!occ_body2)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
      return CUBIT_FALSE;
    }

  CubitBox box_1 = occ_body1->get_bounding_box();
  CubitBox box_2 = occ_body2->get_bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1 = occ_body1->get_TopoDS_Shape();
  TopoDS_Shape *shape2 = occ_body2->get_TopoDS_Shape(); 
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_SOLID); Ex1.More(); Ex1.Next())
    {
      TopoDS_Solid *posolid1 =  new TopoDS_Solid;
      *posolid1 = TopoDS::Solid(Ex1.Current());
      OCCLump * lump1 = new OCCLump(posolid1); 
      for (Ex2.Init(*shape2, TopAbs_SOLID); Ex2.More(); Ex2.Next())
	{
	  TopoDS_Solid *posolid2 =  new TopoDS_Solid;
	  *posolid2 = TopoDS::Solid(Ex2.Current());
	  OCCLump * lump2 = new OCCLump(posolid2);
	  CubitBoolean is_overlap = volumes_overlap(lump1, lump2);
	  if(is_overlap)
	    return CUBIT_TRUE;
	}
    }
  return CUBIT_FALSE;
}

CubitBoolean OCCQueryEngine::volumes_overlap (Lump *lump1, Lump *lump2 ) const
{
  OCCLump *occ_lump1 = CAST_TO(lump1, OCCLump);
  if (!occ_lump1)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
      return CUBIT_FALSE;
    }

  OCCLump *occ_lump2 = CAST_TO(lump2, OCCLump);
  if (!occ_lump2)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
      return CUBIT_FALSE;
    }

  CubitBox box_1 = occ_lump1->bounding_box();
  CubitBox box_2 = occ_lump2->bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1 = (TopoDS_Shape*)(occ_lump1->get_TopoDS_Solid());
  TopoDS_Shape *shape2 = (TopoDS_Shape*)(occ_lump2->get_TopoDS_Solid());
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_FACE); Ex1.More(); Ex1.Next())  
    {
      for (Ex2.Init(*shape2, TopAbs_FACE); Ex2.More(); Ex2.Next()) 
	{
	  BRepAlgoAPI_Section section(Ex1.Current(), Ex2.Current());
	  if (section.HasGenerated())
	    return CUBIT_TRUE;
	}
    }
  return CUBIT_FALSE;
}

void OCCQueryEngine::copy_attributes(TopoDS_Shape& old_shape,
                                     TopoDS_Shape& new_shape)
{
  //update the attribute label tree
  DLIList<CubitSimpleAttrib*> list;
  OCCAttribSet::get_attributes(old_shape, list);

  for(int i = 0; !new_shape.IsNull() && i < list.size(); i ++)
  {
    CubitSimpleAttrib* s_attr = list.get_and_step();
    OCCAttribSet::append_attribute(s_attr, new_shape);
  }
}

int OCCQueryEngine::update_OCC_map(TopoDS_Shape& old_shape, 
                                   TopoDS_Shape& new_shape)
{
  if (!OCCMap->IsBound(old_shape) || old_shape.IsEqual(new_shape))
    return -1;

  //update the attribute label tree
  copy_attributes(old_shape, new_shape); 
  OCCAttribSet::remove_attribute(old_shape);

  //update CGM-OCC map
  int k = OCCMap->Find(old_shape);
  assert (k > 0 && k <= iTotalTBCreated);

  OCCMap->UnBind(old_shape);

  std::map<int, TopologyBridge*>::iterator it = OccToCGM->find(k);
  TopologyBridge* tb = NULL;
  if (it != OccToCGM->end())
     tb = (*it).second;

  else
     assert(0);

  if ((!new_shape.IsNull() && !old_shape.IsSame(new_shape)&& 
      OCCMap->IsBound(new_shape))|| new_shape.IsNull()) 
  //already has a TB built on new_shape
  {
    //delete the second TB corresponding to old_shape
    OccToCGM->erase(k);
    GeometryEntity* ge =  CAST_TO(tb, GeometryEntity);
    if(ge)
      delete_solid_model_entities( ge, CUBIT_FALSE );
    else
    {
      ShellSM * shell = CAST_TO(tb, ShellSM);
      if(shell)
      {
        DLIList<OCCCoFace*> children;
        children = CAST_TO(shell, OCCShell)->cofaces();
        while(children.size())
        {
          OCCCoFace* coface = children.pop();
          CAST_TO(coface->surface(),OCCSurface)->remove_shell(CAST_TO(shell, OCCShell));
          delete coface;
        }

        OCCLump* lump = CAST_TO(shell, OCCShell)->my_lump();
        if(lump)
        {
          delete CAST_TO(shell, OCCShell)->my_body();
          delete lump;
        }
        unhook_ShellSM_from_OCC(shell);
        delete shell;
        return k;
      }
      LoopSM* loop = CAST_TO(tb, LoopSM);
      if(loop)
      {
         DLIList<OCCCoEdge*> children;
         children = CAST_TO(loop, OCCLoop)->coedges();
         while(children.size())
         {
           OCCCoEdge* coedge = children.pop();
           CAST_TO(coedge->curve(), OCCCurve)->remove_loop(CAST_TO(loop, OCCLoop));
           delete (OCCCoEdge*)coedge;
         }
         unhook_LoopSM_from_OCC(loop);
         delete loop;
      }
    }
  }

  else
  {   
    OCCMap->Bind(new_shape, k);
    set_TopoDS_Shape(tb, new_shape);
  }
  return k;
}

void OCCQueryEngine::unhook_cofaces_of_a_surface(OCCSurface* surface)const
{
  DLIList<OCCShell*> shells;
  shells = surface->shells();
  DLIList<OCCCoFace*> cofaces;
  while ( shells.size() > 0)
  {
    DLIList<OCCCoFace*> children ;
    OCCShell* shell = shells.pop();
    children = shell->cofaces();
    for(int j = 0; j < children.size(); j++)
    {
       OCCCoFace* coface = children.get_and_step();
       if (coface->surface() == surface)
       {
          cofaces.append(coface);
          break;
       }
    }
  }
  unhook_CoFaces_from_OCC(cofaces);
}
void OCCQueryEngine::unhook_coedges_of_a_curve(OCCCurve* curve)const 
{
  DLIList<OCCLoop*> loops;
  loops = curve->loops();
  DLIList<OCCCoEdge*> coedges;
  for (int i = 0; i < loops.size(); i ++)
  {
     DLIList<OCCCoEdge *> children ;
     children = loops.get_and_step()->coedges();
     for(int j = 0; j < children.size(); j++)
     {
        OCCCoEdge* coedge = children.get_and_step();
        if (coedge->curve() == curve)
           coedges.append(coedge);
     }
  }
  unhook_CoEdges_from_OCC(coedges);
}
void OCCQueryEngine::set_TopoDS_Shape(TopologyBridge* tb,
                                      TopoDS_Shape& shape)
{
  BodySM* body = CAST_TO(tb, BodySM);
  if(body)
    return CAST_TO(body, OCCBody)->set_TopoDS_Shape(TopoDS::CompSolid(shape));

  Lump* lump = CAST_TO(tb, Lump);
  if(lump)
    return CAST_TO(lump, OCCLump)->set_TopoDS_Solid(TopoDS::Solid(shape));

  ShellSM* shell = CAST_TO(tb, ShellSM);
  if(shell)
    return CAST_TO(shell, OCCShell)->set_TopoDS_Shell(TopoDS::Shell(shape));

  Surface* surface = CAST_TO(tb, Surface);
  if (surface)
    return CAST_TO(surface, OCCSurface)->set_TopoDS_Face(TopoDS::Face(shape));

  LoopSM* loop =  CAST_TO(tb, LoopSM);
  if(loop)
    return CAST_TO(loop, OCCLoop)->set_TopoDS_Wire(TopoDS::Wire(shape));

  Curve* curve = CAST_TO(tb, Curve);
  if (curve)
    return CAST_TO(curve, OCCCurve)->set_TopoDS_Edge(TopoDS::Edge(shape));

  Point* point = CAST_TO(tb, Point);
  if(point)
    return CAST_TO(point, OCCPoint)->set_TopoDS_Vertex(TopoDS::Vertex(shape));
  
}
//EOF
