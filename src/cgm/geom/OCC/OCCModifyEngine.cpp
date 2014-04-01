//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Filename      : OCCModifyEngine.cpp
//
// Purpose       : ModifyEngine for OCC geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Author        : Jane Hu
//
// Creation Date : 1/08
//
//-------------------------------------------------------------------------
#include "gp_Pnt.hxx"
#include "gp_Ax2.hxx"
#include "gp_Dir.hxx"
#include "gp_Hypr.hxx"
#include "gp_Parab.hxx"
#include "gp_Elips.hxx"
#include "gp_Pln.hxx"
#include "gp_Circ.hxx"
#include "gp_Cylinder.hxx"
#include "gp_Cone.hxx"
#include "gp_Sphere.hxx"
#include "gp_Torus.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_MakeDraft.hxx"
#include "BRepBuilderAPI_TransitionMode.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepPrimAPI_MakeHalfSpace.hxx"
#include "BRepBuilderAPI_MakePolygon.hxx"
#include "BRepBndLib.hxx"
#include "IntersectionTool.hpp"
#include "TopoDS_Shape.hxx"
#include "TopAbs_Orientation.hxx"
#include "TopOpeBRep_Point2d.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "GC_MakeArcOfCircle.hxx"
#include "GC_MakeCircle.hxx"
#include "Geom_Circle.hxx"
#include "Geom_SurfaceOfLinearExtrusion.hxx"
#include "Geom_RectangularTrimmedSurface.hxx"
#include "Handle_Geom_RectangularTrimmedSurface.hxx"
#include "GC_MakeArcOfHyperbola.hxx"
#include "GC_MakeArcOfParabola.hxx"
#include "GC_MakeArcOfEllipse.hxx"
#include "GC_MakeSegment.hxx"
#include "GC_MakeTrimmedCone.hxx"
#include "GC_MakeTrimmedCylinder.hxx"
#include "gce_MakeElips.hxx"
#include "BRepFilletAPI_MakeFillet.hxx"
#include "BRepFilletAPI_MakeChamfer.hxx"
#include "BRepAdaptor_CompCurve.hxx"
#include "BRepFilletAPI_MakeFillet2d.hxx"
#include "ChFi2d_ConstructionError.hxx"
#include "Geom_BezierCurve.hxx"
#include "BndLib_AddSurface.hxx"
#include "Handle_Geom_Plane.hxx"
#include "Handle_Geom_OffsetCurve.hxx"
#include "Geom_OffsetCurve.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "Extrema_ExtPC.hxx"
#include "BRepPrimAPI_MakePrism.hxx"
#include "BRepSweep_Revol.hxx"
#include "BRepPrimAPI_MakeCone.hxx"
#include "BRepOffsetAPI_ThruSections.hxx"
#include "BRepLib_FuseEdges.hxx"
#include "BRepOffsetAPI_MakePipe.hxx"
#include "BRepPrimAPI_MakeTorus.hxx"
#include "BRepPrimAPI_MakeCylinder.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GC_MakeEllipse.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "ShapeExtend_Status.hxx"
#include "BRepOffsetAPI_MakeThickSolid.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepBuilderAPI_Copy.hxx"
#include "LocOpe_SplitShape.hxx"
#include "BRep_Tool.hxx"
#include "BRep_Builder.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "TopoDS.hxx"
#include "ShapeFix.hxx"
#include "TopologyBridge.hpp"
#include "ProgressTool.hpp"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepAlgoAPI_Common.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRepPrimAPI_MakeBox.hxx"
#include "BRepPrimAPI_MakeWedge.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "Handle_Geom_TrimmedCurve.hxx"
#include "Handle_ShapeBuild_ReShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "Handle_Geom_RectangularTrimmedSurface.hxx"
#include "BndLib_Add3dCurve.hxx"
#include "TopOpeBRep_EdgesIntersector.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"
#include "OCCDrawTool.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "OCCCoFace.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "TopOpeBRep_ShapeIntersector.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "CubitUtil.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"
#include "CubitFileIOWrapper.hpp"
#include "Body.hpp"
#include "GfxDebug.hpp"
#include "RefFace.hpp"
#include "CpuTimer.hpp"
#include "AppUtil.hpp"
#include "SphereEvaluator.hpp"
#include "CylinderEvaluator.hpp"
#include "GfxPreview.hpp"
#include <vector>
OCCModifyEngine* OCCModifyEngine::instance_ = 0;
#define DEBUG
double TOL = 0.0;
//===============================================================================
// Function   : OCCModifyEngine
// Member Type: PUBLIC
// Description: constructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::OCCModifyEngine()
{
//  assert( !instance_ );

    // add this modify engine to geometrymodifytool
  GeometryModifyTool::instance()->add_gme(this);
  TOL = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
}


//===============================================================================
// Function   : ~OCCModifyEngine
// Member Type: PUBLIC
// Description: destructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::~OCCModifyEngine() 
{
        instance_ = 0;
}

//===============================================================================
// Function   : make_Point
// Member Type: PUBLIC
// Description: make a geometric entity point
// Author     : Jane Hu 
// Date       : 10/07
//===============================================================================
Point* OCCModifyEngine::make_Point( CubitVector const& point) const
{
  gp_Pnt pt = gp_Pnt( point.x(), point.y(), point.z());
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);

  // Create a new PointACIS object
  return OCCQueryEngine::instance()->populate_topology_bridge( theVertex );
}

//===============================================================================
// Function   : make_Curve
//              This function creates a curve given an existing curve, copy. 
// Member Type: PUBLIC
// Description: make a curve
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve(Curve * curve_ptr) const
{
  OCCCurve* occ_curve = CAST_TO(curve_ptr, OCCCurve);
  if (!occ_curve)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
 
  TopoDS_Edge *theEdge = occ_curve->get_TopoDS_Edge();  
 
  BRepBuilderAPI_Copy api_copy(*theEdge);

  TopoDS_Shape newShape = api_copy.ModifiedShape(*theEdge);
 
  TopoDS_Edge newEdge = TopoDS::Edge(newShape);

  return OCCQueryEngine::instance()->populate_topology_bridge(newEdge);
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve by projecting a straight line defined by 
//              point1_ptr, and point2_ptr onto face_ptr, third_point
//              is used for curves that could be periodic to dertermine
//              the correct direction.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( Point const* point1_ptr,
                             Point const* point2_ptr,
                             Surface* face_ptr,
                             const CubitVector * third_point) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  GeometryType type = STRAIGHT_CURVE_TYPE;
  CubitBoolean closed = CUBIT_FALSE;
  DLIList<CubitVector*> mid_points;
  Curve* curve = NULL;
  if (point1_ptr != point2_ptr)
    curve = make_Curve(type, point1_ptr, point2_ptr, mid_points);
  else //could be a closed shape
  {
    if(third_point != NULL && face_ptr != NULL) 
    {
       closed = CUBIT_TRUE;
       Point * Pnt = make_Point(*third_point);
       curve = make_Curve(type, point1_ptr, Pnt, mid_points);
    }
    else
    {
       PRINT_ERROR("Cannot create an OCC curve from the given duplicated points.\n");
       return (Curve *)NULL;
    }
  }

  Curve* new_curve = NULL;
  if(face_ptr == NULL)
    return curve;
 
  DLIList<Point*> points;
  new_curve = 
	CAST_TO(curve, OCCCurve)->project_curve(face_ptr, points, closed, third_point); 

  delete curve;
  return new_curve;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a  spline curve by using the points on surface.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             DLIList<CubitVector*>& vector_list,
                             Surface* face_ptr) const
{
  assert(point1_ptr != NULL && point2_ptr != NULL);
  
  if (curve_type != SPLINE_CURVE_TYPE
      && curve_type != STRAIGHT_CURVE_TYPE)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve_type.\n"
                 "Candidates are SPLINE_CURVE_TYPE and STRAIGHT_CURVE_TYPE.\n");
     return (Curve *)NULL;
  }

  if (curve_type == STRAIGHT_CURVE_TYPE)
    return make_Curve(curve_type, point1_ptr, point2_ptr, NULL, CUBIT_FORWARD);

  OCCPoint* occ_point1 = CAST_TO(const_cast<Point*>(point1_ptr), OCCPoint);
  OCCPoint* occ_point2 = CAST_TO(const_cast<Point*>(point2_ptr), OCCPoint);

  if (occ_point1 == NULL || occ_point2 == NULL)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given points.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
    
  //project all points on the surface if possible
  OCCSurface* occ_face = NULL;
  if (face_ptr != NULL)
     occ_face = CAST_TO(face_ptr, OCCSurface);
 
  gp_Pnt pt;
  int size = 2+vector_list.size();
  TColgp_Array1OfPnt points(1, size);
  CubitVector* vector = NULL;
  CubitVector closest_location;
  for(int i = 1; i <= size; i++)
  {
     if (i == 1) 
     {
       TopoDS_Vertex *point = occ_point1->get_TopoDS_Vertex();
       pt = BRep_Tool::Pnt(*point);
       vector = new CubitVector(point1_ptr->coordinates());
     }
     else if (i == size)
     {
       TopoDS_Vertex *point = occ_point2->get_TopoDS_Vertex();
       pt = BRep_Tool::Pnt(*point);
       vector = new CubitVector(point2_ptr->coordinates()); 
     } 
     else
     {
       vector = vector_list.get_and_step();
       pt.SetCoord(vector->x(), vector->y(), vector->z());
     } 

     if (occ_face != NULL)
     {
       occ_face->closest_point(*vector, &closest_location);
       pt.SetCoord(closest_location.x(), closest_location.y(), closest_location.z()) ;  	 
     }

     points.SetValue(i, pt); 
  }    
     
  //make curve according to the curve type.
  if(curve_type == SPLINE_CURVE_TYPE)
  {
    Geom_BezierCurve BezierCurve(points);
    Geom_BezierCurve* curve =  new Geom_BezierCurve(BezierCurve);
    Handle(Geom_BoundedCurve) curve_ptr(curve);
    TopoDS_Vertex * vt1 = occ_point1->get_TopoDS_Vertex();
    TopoDS_Vertex * vt2 = occ_point2->get_TopoDS_Vertex(); 
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt2);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge); 
  }

  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// For STRAIGHT_CURVE_TYPE:
//    intermediate_point_ptr  is not used
//
// For PARABOLA_CURVE_TYPE
//    intermediate_point_ptr is the tip of the parabola
//
// For HYPERBOLA_CURVE_TYPE
//    intermediate_point_ptr is the center of its two foci
//
// For ELLIPSE_CURVE_TYPE
//    intermediate_point_ptr is the center of the ellipse
//    the two points are vertices, one gives the major radius, 
//    the other point gives the minor radius.
//    sense is used to determine which part of the ellipse is required
//
// For ARC_CURVE_TYPE
//    arc passes three points
//
// Author     : Jane Hu 
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             CubitVector const* intermediate_point_ptr,
                             CubitSense sense) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  DLIList<CubitVector*> mid_points;
  if (intermediate_point_ptr)
  {
    CubitVector mid_point = *intermediate_point_ptr;
    mid_points.append(&mid_point);
  }

  CubitVector v1(point1_ptr->coordinates());
  CubitVector v2(point2_ptr->coordinates());

  gp_Pnt pt1(v1.x(),v1.y(), v1.z());
  gp_Pnt pt2(v2.x(),v2.y(), v2.z());

  CubitVector v3;
  gp_Pnt pt3;

  Handle(Geom_TrimmedCurve) curve_ptr;
  if(intermediate_point_ptr != NULL)
  {
    v3 = *intermediate_point_ptr;
    pt3.SetCoord(v3.x(),v3.y(), v3.z());
  }

  if (curve_type == STRAIGHT_CURVE_TYPE)
     curve_ptr = GC_MakeSegment(pt1,pt2);

  else if (curve_type == ARC_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     curve_ptr = GC_MakeArcOfCircle(pt1, pt3, pt2);
  }

  else if (curve_type == ELLIPSE_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     
     gp_Pnt center(v3.x(), v3.y(), v3.z());

     gp_Elips ellipse;
     gce_MakeElips ellipse1(pt1	, pt2	, center);
     if(ellipse1.IsDone())
       ellipse = ellipse1.Value();
     else if(!ellipse1.IsDone() && ellipse1.Status() == gce_InvertRadius)
     {
        gce_MakeElips ellipse2(pt2, pt1, center);
        if(ellipse2.IsDone())
          ellipse = ellipse2.Value();
        else
        {
          PRINT_ERROR("Can't create an ellipse from give 3 points.\n");
          return (Curve *)NULL;
        }      
     } 
     else
     {
        PRINT_ERROR("Can't create an ellipse from give 3 points.\n");
        return (Curve *)NULL;
     }
     CubitBoolean use_sense = (sense == CUBIT_FORWARD ? CUBIT_TRUE : CUBIT_FALSE); 
     curve_ptr = GC_MakeArcOfEllipse(ellipse, pt1, pt2, use_sense);
  }

  else if(curve_type == PARABOLA_CURVE_TYPE || 
          curve_type == HYPERBOLA_CURVE_TYPE)
  {
    assert(intermediate_point_ptr != NULL);

    //find the directrix and focus of the parabola
    //or the axis, major radius and minor radius of the hyperbola
    CubitVector width_vec = v2 - v1;
    if(width_vec.length_squared() < TOL * TOL)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "2 end points are the same.\n");
       return (Curve *)NULL;
    }

    CubitVector midpoint_vec = (v1 + v2)/2.0;
    CubitVector height_vec = midpoint_vec - v3;
    gp_Pnt center(v3.x(), v3.y(), v3.z());
 
    if (height_vec.length_squared() < TOL * TOL)
    { 
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    CubitVector x = height_vec;
    x.normalize();
    gp_Dir x_dir(x.x(), x.y(), x.z());
 
    CubitVector N = x * (v2 - v1);  
    if (N.length_squared() < TOL * TOL)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    N.normalize();
    gp_Dir N_dir(N.x(), N.y(), N.z());

    gp_Ax2 axis(center, N_dir, x_dir);  

    if(curve_type == HYPERBOLA_CURVE_TYPE)
    { 
       //    (focus2) (v3) . (v2)
       //          .   .   . (midpoint = focus1)
       //                  . (v1)
       CubitVector focus2 = 2 * v3 - midpoint_vec;

       //according to the definition of hyperbola,
       //2 * a = length(v2 - focus2)-length(v2 - focus1)

       double major = (v2 - focus2).length()/2.0 - (v2 - midpoint_vec).length()/2.0;

       // if a = 1/2 length major axis, b = 1/2 length minor axis and
       // c = distance center to focus, then a*a + b*b = c*c

       double c_squared = (midpoint_vec - v3).length_squared();
       double minor = sqrt(c_squared  - major*major );
       gp_Hypr hypt(axis, major, minor);
       curve_ptr =
             GC_MakeArcOfHyperbola(hypt, pt1, pt2, CUBIT_TRUE);
    }

    else
    {
       // Find the focus of this parabola.
       // Since for a parabola with its peak at the origin, y = (1/(4*a))*x^2,
       // and since we have restricted this parabola to be symmetric (per the 
       // FastQ method, see the FastQ file getwt.f), we can use the following 
       // relationship to
       // determine "a", the distance the focus lies from the peak on the line
       // formed by the peak and the midpoint of the start and end points`
       double a = width_vec.length_squared()/(16. * height_vec.length()); 
       gp_Parab parab(axis, a);
       curve_ptr =
		GC_MakeArcOfParabola(parab, pt1, pt2, CUBIT_TRUE);
    } 
  }

  else
  {
      PRINT_ERROR("In OCCModifyEngine::make_Curve\n"
                  "       Invalid curve type.\n");
      return (Curve *)NULL;
  }

  OCCPoint* occ_pt1 = CAST_TO(const_cast<Point*>(point1_ptr),OCCPoint);
  OCCPoint* occ_pt2 = CAST_TO(const_cast<Point*>(point2_ptr),OCCPoint);
  TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
  TopoDS_Vertex * vt2 = occ_pt2->get_TopoDS_Vertex();
  TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt2);
  return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
}

//===============================================================================
// Function   : make_Surface
//              This function creates a surface given an existing surface, copy.
// Member Type: PUBLIC
// Description: make a surface, OCC allows to create a stand along surface,
//              however the CGM has a design of making all surfaces in a (sheet)
//              body. This will add complexity on all free surface related 
//              calculation and modification, and adding potential bugs too. 
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( Surface * surface_ptr,
                                 CubitBoolean extended_from) const
{
  OCCSurface* occ_surface = CAST_TO(surface_ptr, OCCSurface);
  if (!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  //Start of the codes
  double UMax, VMax, UMin, VMin;
  occ_surface->get_param_range_U(UMin, UMax);
  occ_surface->get_param_range_V(VMin, VMax);

  TopoDS_Face *theFace = occ_surface->get_TopoDS_Face();
  if( !theFace)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  TopoDS_Face newFace;
  BRepAdaptor_Surface asurface(*theFace);

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  double const height = 2*(bounding_box.diagonal()).length();
  CubitBox box = occ_surface->bounding_box();
  double ratio = height/(box.diagonal().length());

  double middleU = (UMin + UMax)/2.0;
  double middleV = (VMin + VMax)/2.0;
  double U1 = middleU - (UMax-UMin)/2.0 * ratio;
  double U2 = middleU + (UMax-UMin)/2.0 * ratio;
  double V1 = middleV - (VMax - VMin)/2.0 * ratio;
  double V2 = middleV + (VMax - VMin)/2.0 * ratio;

  if (extended_from == CUBIT_TRUE)
  {
     // We need to get the type of surface.
     GeometryType type = occ_surface->geometry_type();
     if (type  == PLANE_SURFACE_TYPE)
     {
        gp_Pln plane = asurface.Plane();
        newFace = BRepBuilderAPI_MakeFace(plane, U1, U2, V1, V2);
     }
     else if(type == CONE_SURFACE_TYPE)
     {
       //make an infinite cone.
       //Given this lets create another face that is extended from it.
       if(asurface.GetType() == GeomAbs_Cone)
       {
         gp_Cone cone = asurface.Cone();
         gp_Pnt Apex = cone.Apex();
         double semi_angle = cone.SemiAngle();
         gp_Pnt p2;
         double radius2;
         gp_XYZ xyz;
         if (semi_angle > 0)
           xyz = Apex.XYZ() + cone.Position().Direction().XYZ()*height;
             
  	 else
	   xyz = Apex.XYZ() - cone.Position().Direction().XYZ()*height;

         p2.SetXYZ(xyz);
	 radius2 = height * tan(fabs(semi_angle));
         Handle(Geom_RectangularTrimmedSurface) trimmed_cone;
         trimmed_cone = GC_MakeTrimmedCone(Apex, p2, 0, radius2); 
         newFace = BRepBuilderAPI_MakeFace(trimmed_cone);
       }
       else
       {
         gp_Cylinder cylinder = asurface.Cylinder();
         //gp_Pnt pnt = asurface.Value(UMin,V1); 
         double radius = cylinder.Radius();
         gp_Ax1 axis = cylinder.Axis(); 
         Handle(Geom_RectangularTrimmedSurface) trimmed_cyl;
         trimmed_cyl = GC_MakeTrimmedCylinder(axis, radius, height);
         newFace = BRepBuilderAPI_MakeFace(trimmed_cyl);
       } 
     }
     else if(type == SPHERE_SURFACE_TYPE)
     {
       //make a whole sphere.
       gp_Sphere sphere = asurface.Sphere();
       newFace = BRepBuilderAPI_MakeFace(sphere); 
     }
     else if(type == TORUS_SURFACE_TYPE)
     {
       //make a whole torus
       gp_Torus torus = asurface.Torus();
       newFace = BRepBuilderAPI_MakeFace(torus);
     }
     else if(type == SPLINE_SURFACE_TYPE ) 
     {
       //extend the surfaces using the equation if possible.
       Handle(Geom_BezierSurface) bezier = asurface.Bezier();
       newFace = BRepBuilderAPI_MakeFace(bezier, U1, U2, V1, V2);
     }
  }
 
  else
  {
    BRepBuilderAPI_Copy api_copy(*theFace);
    TopoDS_Shape newShape = api_copy.ModifiedShape(*theFace);
    newFace = TopoDS::Face(newShape);
  }
  
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               newFace, CUBIT_TRUE);

  return surface;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface of type surface_type, given the list of curves.
//              check edges option is done in GeometryModifyTool level, so 
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( GeometryType surface_type,
                                 DLIList<Curve*>& curve_list,
                                 Surface * old_surface_ptr,
                                 bool check_edges) const
{
  //Create TopoDS_Edge list to make a surface.
  DLIList<DLIList<TopoDS_Edge*>*> topo_edges_loops;
  curve_list.reset() ;
    
  //check no intersections of the TopoDS_Edge's.
  //need to check that no intersection in the middle of the curves, not at
  //vertices or out of boundary.

  int count = 0; //intersection point should be the same as curve_list size.
  for ( int i = 0 ; i < curve_list.size()-1 ; i++ )
  {
     for(int j = i+1; j < curve_list.size(); j ++)
     {
        DLIList<CubitVector*> intscts;
 	CubitBoolean bounded = CUBIT_TRUE;//dummy arg.
	CubitBoolean closest = CUBIT_TRUE;//dummy arg.
        CubitStatus yes_int = 
              OCCQueryEngine::instance()->get_intersections(curve_list[i],
				curve_list[j], intscts, bounded, closest);
        if(yes_int)
        {
           //check intscts point should be vertex or outside boundary.
 	   if (intscts.size() > 2 )  
	   {
	     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
             return (Surface *)NULL;
           }
           else
           {
             for(int k = 0; k < intscts.size(); k++)
             {
               CubitVector *v = intscts.get_and_step();
	       CubitPointContainment is_on = CAST_TO(curve_list[i],OCCCurve)->
					point_containment(*v);
               if (is_on == CUBIT_PNT_BOUNDARY)
               {
	 	 is_on = CAST_TO(curve_list[j],OCCCurve)->
				point_containment(*v);
		 if (is_on == CUBIT_PNT_BOUNDARY)
                   count++;
               }
               else if(is_on == CUBIT_PNT_INSIDE)
               {
                 PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
                 return (Surface *)NULL;
               }
	     }
	   }
        }
     }
  }
 
  if (count > curve_list.size()) 
  {
      PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                "       Cannot make Surface with intersecting curves.\n");
      return (Surface *)NULL;
  }

  CubitStatus stat = sort_curves(curve_list, topo_edges_loops); 
  if( stat == CUBIT_FAILURE ) //case of one disconnected curve 
     return (Surface*) NULL;
 
  // Use the topo_edges to make a topo_face
  TopoDS_Face* topo_face;
  topo_face = make_TopoDS_Face(surface_type,topo_edges_loops, old_surface_ptr);
 
  for (int i = 0; i < topo_edges_loops.size(); i++)
  {
    DLIList<TopoDS_Edge*>* topo_edges = topo_edges_loops.get_and_step();
    for(int j = 0; j < topo_edges->size(); j++)
      topo_edges->pop();
    delete topo_edges;
  }
  
  if(!topo_face)
  {
     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface object.\n");
     return (Surface *)NULL;
  }

  // make the topology bridges for the face
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               *topo_face, CUBIT_TRUE); 
  topo_face->Nullify();
  delete topo_face;
  return surface ;
}

//===============================================================================
// Function   : sort_curves
// Member Type: PROTECTED
// Description: sort the curves so they are in order and make closed loop 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::sort_curves(DLIList<Curve*> curve_list,
                        DLIList<DLIList<TopoDS_Edge*>*>& topo_edges_loops)const
{
  topo_edges_loops.clean_out();
  CubitStatus stat = CUBIT_SUCCESS;
  std::vector< DLIList<TopoDS_Edge*>* > topo_edges(curve_list.size());
  int size_in = curve_list.size();
  for(int i = 0; i < size_in; i++)
    topo_edges[i] = new DLIList<TopoDS_Edge*>;

  curve_list.reset() ;
  Curve const* curve_ptr = NULL ;
  OCCCurve* occ_curve = NULL;
  TopoDS_Edge* topo_edge = NULL;

  OCCPoint* start = NULL;
  OCCPoint* end = NULL;
  DLIList<OCCPoint*> point_list;
  CubitBoolean new_end;
  int size = curve_list.size();

  int count = 0;
  for ( int i = 0 ; i < size ; i++ )
  {
     if (i == 0)
       new_end = CUBIT_TRUE;
     for(int j = 0; j < curve_list.size(); j ++)
     {
        curve_ptr = curve_list.get() ; 
        occ_curve = CAST_TO(const_cast<Curve*>(curve_ptr), OCCCurve);

        if(occ_curve ==  NULL)
        {
           PRINT_ERROR("In OCCModifyEngine::sort_curves\n"
                       "       Got a NULL pointer to OCCCurve\n") ;
           return CUBIT_FAILURE;
        }

        point_list.clean_out();
        occ_curve->get_points(point_list);
        //assert(point_list.size()==2);

        if (i == 0)
        {
          start = point_list.get();
          end = point_list.pop();  
          break;
        }

        if(end->is_equal(*(point_list.get()), TOL) ||
           end->is_equal(*(point_list.step_and_get()),TOL)) 
        {
           end = point_list.step_and_get();
           new_end = CUBIT_TRUE;
           break;
        }
   
        else if(start->is_equal(*(point_list.get()), TOL) ||
           start->is_equal(*(point_list.step_and_get()),TOL))
        {
           start = end;
           end = point_list.step_and_get(); 
           new_end = CUBIT_TRUE;
           break;
        }
        curve_list.step();
     }

     if (new_end)//found next curve 
     {
        topo_edge = occ_curve->get_TopoDS_Edge();
        topo_edges[count]->append(topo_edge);
        curve_list.remove();
        if(start->is_equal( *end, TOL))  //formed a closed loop
        {
          i = -1;
          size = curve_list.size() ;
          topo_edges_loops.append(topo_edges[count]);
          count++;
        }
        else
          new_end = CUBIT_FALSE;
     }
     else
     {
        stat = CUBIT_FAILURE; 
        i = -1;
        size = curve_list.size();
        topo_edges_loops.append(topo_edges[count]);
        count++;
     }
  }

  if( new_end == CUBIT_FALSE ) //case of one disconnected curve
  {
    topo_edges_loops.append(topo_edges[count]); 
    stat = CUBIT_FAILURE;
  }

  for(int i = 0; i < size_in; i++)
  {
     if(topo_edges[i]->size() == 0)
       delete topo_edges[i];
  }
  return stat;
} 
//===============================================================================
// Function   : make_TopoDS_Face
// Member Type: PROTECTED
// Description: make a opoDS_Face of type surface_type, given the list of 
//              TopoDS_Edge. the TopoDS_Edge's should be in order in loops.
//              check edges option is done in GeometryModifyTool level, so
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
TopoDS_Face* OCCModifyEngine::make_TopoDS_Face(GeometryType surface_type,
			      DLIList<DLIList<TopoDS_Edge*>*> topo_edges_list,
			      Surface * old_surface_ptr)const
{
  TopoDS_Face* topo_face = NULL;
  // Make sure a supported type of surface is being requested.
  if ( surface_type != PLANE_SURFACE_TYPE  &&
       surface_type != BEST_FIT_SURFACE_TYPE)
  {
      PRINT_WARNING("In OCCGeometryEngine::make_TopoDS_Face\n"
                    "       At this time, cannot make a TopoDS_Face that isn't"
                    " planar or best fit.\n");
      return topo_face;
  }
 
  // Set the TopoDS_Face pointer, if requested.
  TopoDS_Face *fit_Face = NULL;
  Handle_Geom_Surface S;
  if ( old_surface_ptr != NULL )
  {
      OCCSurface *surf = CAST_TO(old_surface_ptr, OCCSurface );
      fit_Face = surf->get_TopoDS_Face();
      S = BRep_Tool::Surface(*fit_Face);
  }
 
  // Make a wire from the topo_edges.
  // Coincident TopoDS_Vertex will be deleted by OCC.
  if(topo_edges_list.size() == 0)
      return topo_face;

  DLIList<TopoDS_Wire*> wires;
  GProp_GProps myProps;
  double max_area  = 0.0;
  TopoDS_Wire* out_Wire = NULL;
  TopoDS_Wire test_Wire;
  DLIList<TopoDS_Edge*>* topo_edges; 
  //check and make sure the outer loop is in the first
  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    topo_edges = topo_edges_list.get_and_step();
    BRepBuilderAPI_MakeWire aWire(*(topo_edges->get()));
    for(int j = 1; j < topo_edges->size(); j++)
      aWire.Add(*(topo_edges->step_and_get()));

    test_Wire = aWire.Wire();
    wires.append(&test_Wire);
   
    if (topo_edges_list.size() == 1)
      break;

    BRepBuilderAPI_MakeFace made_face(test_Wire);

    if (!made_face.IsDone())
    {
       PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                   "   Cannot find the best fit surface for given curves.\n");
       return topo_face;
    }
    TopoDS_Face test_face = made_face.Face();
    BRepGProp::SurfaceProperties(test_face, myProps); 
    double area = myProps.Mass();
    out_Wire = max_area > area ? out_Wire : &test_Wire;
    max_area = max_area > area ? max_area : area;
  } 

  if (out_Wire)
  {
    wires.remove(out_Wire);
    wires.insert_first(out_Wire);
  }

  //create the TopoDS_Face
  CubitBoolean error = CUBIT_FALSE;

  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    TopoDS_Wire *the_wire = wires.get_and_step();
    if (i == 0)
    {
      if( old_surface_ptr != NULL )
      {
        BRepBuilderAPI_MakeFace made_face(S, *the_wire);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }
        topo_face = new TopoDS_Face(made_face.Face());
      }
      else
      {
        CubitBoolean is_planar = (surface_type == PLANE_SURFACE_TYPE ?
				  CUBIT_TRUE : CUBIT_FALSE); 
        BRepBuilderAPI_MakeFace made_face(*the_wire, is_planar);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }

        topo_face = new TopoDS_Face(made_face.Face());
      }
    }
    else
    {
      BRepBuilderAPI_MakeFace made_face(*topo_face, *the_wire);
      if (!made_face.IsDone())
      {
        error = CUBIT_TRUE;
        break;
      }
      delete topo_face;
      topo_face = new TopoDS_Face(made_face.Face());
    }
  } 

  if(error)
  {
    PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                 "   Cannot find the best fit surface for given curves.\n");
    return (TopoDS_Face*) NULL;
  }

  return topo_face;
}
//===============================================================================
// Function   : make_Lump
// Member Type: PUBLIC
// Description: make a lump of one shell
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Lump* OCCModifyEngine::make_Lump( DLIList<Surface*>& surface_list ) const
{
  if (surface_list.size() < 2) 
    return (Lump*) NULL;

  //all surfaces should be stand along surface bodies or shell bodies' surface
  DLIList<BodySM*> body_list;
  for(int i = 0; i < surface_list.size(); i++)
  {
    OCCSurface* occ_surface = CAST_TO(surface_list.get_and_step(), OCCSurface);
    if (occ_surface == NULL)
    {
       PRINT_ERROR("Cannot create an OCC lump from the given surfaces.\n"
                 "Possible incompatible geometry engines.\n");
       return (Lump *)NULL;
    }
    OCCBody* occ_body = occ_surface->my_body();
    if(occ_body && occ_body->my_sheet_surface() == NULL)
    {
      PRINT_ERROR("Cannot create an OCC lump from the given surfaces.\n"
                 "The surfaces are not free.\n");
      return (Lump *)NULL;
    }
    else if(!occ_body)
    {
      OCCShell* shell = occ_surface->my_shell();
      if(!shell)
      {
        PRINT_ERROR("This is a bug, please report it. \n");
        return (Lump*)NULL;
      }
      occ_body = shell->my_body();
    }
    body_list.append_unique(occ_body);
  }

  TopoDS_Shape aShape;
  CubitStatus stat = stitch_surfs(body_list, aShape);
  if(!stat)
  {
    PRINT_ERROR("The surfaces are not all connected, can't make a lump. \n");
    return (Lump*)NULL;
  }

  TopExp_Explorer Ex, Ex2;
  TopoDS_Shell aShell ;
  for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); Ex.More()&& stat; Ex.Next())
    aShell = TopoDS::Shell(Ex.Current());
 
  //check to make sure the aShell is closed.
  int num_edges = 0;
  int pairs = 0;
  //sometimes there's duplicate TopoDS_Edges in the shell.
  DLIList<TopoDS_Edge*> edge_list;
  for (Ex.Init(aShell, TopAbs_EDGE); Ex.More()&& stat; Ex.Next())
  {
    TopoDS_Edge edge1 = TopoDS::Edge(Ex.Current());
    TopoDS_Edge* new_edge = new TopoDS_Edge(edge1);
    edge_list.append(new_edge);
  }

  int size = edge_list.size();
  for (int i = 0; i < size; i++)
  {
    TopoDS_Edge edge1 = *edge_list[i];
    int same = 0;
    for (int j = i+1; j < edge_list.size(); j++)
    {
      TopoDS_Edge edge2 = *edge_list[j];
      if(edge1.IsEqual(edge2))
      {
           same ++;
           edge_list.remove(&edge1);
           i--;
           size--;
           break;
      }
    }
    if(same > 0)
      continue;

    else
      num_edges++;
  
    for (int j = 0; j < size; j++)  
    {
      TopoDS_Edge edge2 = *edge_list[j];    
      if (!edge1.IsEqual(edge2)&& edge1.IsSame(edge2))
      {
        pairs++;
        break;
      }
    }
  }

  for (int k = 0; k < edge_list.size(); k++)
  {
    TopoDS_Edge* edge = edge_list.get_and_step();
    edge->Nullify();
    delete edge;
  }

  if (num_edges == pairs )
    aShell.Closed(CUBIT_TRUE);

  else
    PRINT_ERROR("Surfaces must make a water-tight shape to make a lump.\n");
  
  if(aShell.Closed())
  {
    BRepBuilderAPI_MakeSolid aMakeSolid(aShell);
    if (!aMakeSolid.IsDone())
    {
       PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                   "OCC internal error.\n");
       return (Lump *)NULL;
    }

    TopoDS_Solid aSolid = aMakeSolid.Solid();

    return
      OCCQueryEngine::instance()->populate_topology_bridge(aSolid, CUBIT_TRUE); 
  }

  return (Lump*) NULL;
}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM from a surface
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( Surface *surface ) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  if(!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC body from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  OCCBody* occ_body = occ_surface->my_body();
  if(occ_body)
     return occ_body;

  TopoDS_Face* face = occ_surface->get_TopoDS_Face();
  surface = OCCQueryEngine::instance()->populate_topology_bridge(*face, CUBIT_TRUE);
   
  return CAST_TO(surface, OCCSurface)->my_body();
}



//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM given a list of Lumps.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( DLIList<Lump*>& lump_list ) const
{
  if (lump_list.size() == 0)
    return (BodySM*) NULL;

  //Create a compsolid shape, copy all BodySM's solids to create new compbody 
  DLIList<BodySM*> bodysm_list;
  TopoDS_CompSolid CS;
  BRep_Builder B;
  B.MakeCompSolid(CS);

  //Add every shape to the CompSolid
  for(int i = 0; i < lump_list.size(); i++)
  {
     Lump* lump = lump_list.get_and_step();
     OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     if(!occ_lump)
     {
        PRINT_ERROR("Cannot create an OCC BodySM from the given lumps.\n"
                    "Possible incompatible geometry engines.\n");
        return (BodySM *)NULL;
     }
     TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();
     BRepBuilderAPI_Copy api_copy(*solid);
     TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
     B.Add(CS, newShape);
  }
 
  BodySM* bodysm = OCCQueryEngine::instance()->populate_topology_bridge(CS);

  return bodysm;

}


//===============================================================================
// Function   : sphere
// Member Type: PUBLIC
// Description: build an OCC sphere
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::sphere(double radius) const
{
  if (radius <= 0)
    return (BodySM*) NULL; 

  TopoDS_Solid S = BRepPrimAPI_MakeSphere(radius);
  
  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S, 
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->get_body();
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: build an OCC brick 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( double wid, double dep, double hi ) const
{
  if (wid <= 0 || dep <=0 || hi <= 0)
    return (BodySM*)NULL;
  
  TopoDS_Solid S = BRepPrimAPI_MakeBox(wid, dep, hi);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);

  if (lump == NULL)
    return (BodySM*)NULL;

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  if(body)
    CAST_TO(body,OCCBody)->move(-wid/2.0, -dep/2.0, -hi/2.0);
  return body;
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: create an OCC brick given center axes and extension
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( const CubitVector& center, 
                                const CubitVector axes[3],
                                const CubitVector &extension) const
{
  CubitVector left_corner =  center - 0.5 * extension;
  gp_Pnt left_point(left_corner.x(), left_corner.y(), left_corner.z());
  gp_Dir main_dir(axes[2].x(), axes[2].y(), axes[2].z());
  gp_Dir x_dir(axes[0].x(), axes[0].y(), axes[0].z());
  gp_Ax2 Axis(left_point, main_dir, x_dir);
  TopoDS_Solid S = BRepPrimAPI_MakeBox( Axis, extension.x(), extension.y(),
					extension.z());

  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->get_body(); 
}

//===============================================================================
// Function   : prism
// Member Type: PUBLIC
// Description: create an OCC prism 
// Author     : Jane Hu  
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::prism( double height, int sides, double major,
                               double minor) const
{
  if(major <= 0. || minor <= 0. || (major - minor) <=  -TOL)
  {
    PRINT_ERROR("Major and minor radii must be greater than zero.\n");
    return (BodySM*)NULL;
  }

  if (sides == 4)
    return brick(2 * major, 2 * minor, height); 

  TopoDS_Wire wire ;
    make_base_for_prim_pyramid(major, minor, height, sides, wire);

  TopoDS_Face base = BRepBuilderAPI_MakeFace(wire, Standard_True);
  gp_Dir main_dir(0.0, 0.0, 1.0);
  gp_Vec norm(main_dir);
  norm *= height;
  BRepSweep_Prism swept(base, norm);
  TopoDS_Shape new_shape = swept.Shape();
  DLIList<TopologyBridge*> tbs;
  tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
  assert(tbs.size() == 1);

  BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
  return bodysm;
}

void OCCModifyEngine::make_base_for_prim_pyramid(double major,
                                                 double minor,
                                                 double height,
                                                 int sides,
                                                 TopoDS_Wire& wire)const
{
  //One of the polygon side will be perpendicular to positive x-axis.
  double y = major * sin(CUBIT_PI/sides);
  double x = sqrt(major * major - y * y);
  gp_Pnt start(x, y, -height/2.0);

  DLIList<gp_Pnt> point_list;
  double theta = 2.0/sides*CUBIT_PI;
  for(int n =1 ; n < sides; n++)
  {
    double angle = theta * (n + 0.5);
    gp_Pnt v(major * cos(angle), major * sin(angle), -height/2.0);
    point_list.append(v);
  }

  TopoDS_Edge new_edge;
  BRepBuilderAPI_MakePolygon poly_maker;
  gp_Dir main_dir(0.0, 0.0, 1.0);
  point_list.append(start);

  if (fabs(major - minor) < TOL)
    for (int i = 0; i <sides; i++)
      poly_maker.Add(point_list[i]);
  else
  {
    for (int i = 0; i <sides; i++)
    {
      x = point_list[i].X();
      if((y=point_list[i].Y()) > 0.0)
        y = sqrt((1-x*x/major/major)*minor*minor);
      else
        y = -sqrt((1-x*x/major/major)*minor*minor);
      point_list[i].SetY(y);
      poly_maker.Add(point_list[i]);
    }
  }
  poly_maker.Close();

  wire = poly_maker.Wire();
}

//===============================================================================
// Function   : pyramid
// Member Type: PUBLIC
// Description: create an OCC pyramid 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::pyramid( double height, int sides, double major,
                                 double minor, double top) const
{
  TopoDS_Solid S;

  //build the top and bottom shapes.
  TopoDS_Wire wire_bottom ; 
  make_base_for_prim_pyramid(major, minor, height, sides, wire_bottom);
  double top_minor = top * minor / major;
  TopoDS_Wire wire_top ; 
  BRepOffsetAPI_ThruSections builder(CUBIT_TRUE, CUBIT_TRUE);
  builder.AddWire(wire_bottom);
  if(top > TOL)
  {
    make_base_for_prim_pyramid(top, top_minor, -height, sides, wire_top); 
    builder.AddWire(wire_top);
  }
  else
  {
    gp_Pnt pt = gp_Pnt( 0.0, 0.0, height/2.0);
    TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
    builder.AddVertex(theVertex);
  }
  builder.Build() ;
  S = TopoDS::Solid(builder.Shape());
 
  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  return body;
  
}

//===============================================================================
// Function   : cylinder
// Member Type: PUBLIC
// Description: create an OCC cylinder, its base shape can be ellipse with
//		r1, r2 while r3 = 0; or it can be a cone with r1, 
//		r3 while r2 = 0. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::cylinder( double hi, double r1, double r2, double r3 ) const
{
  TopoDS_Solid S;
  if(r2 != r1)//elliptical based cylinder
  {
    gp_Pnt center(0.0, 0.0, 0.0);
    gp_Dir main_dir(0.0, 0.0, 1.0);
    gp_Dir x_dir(1.0, 0.0, 0.0);
    gp_Ax2 Axis(center, main_dir, x_dir); 
    Handle(Geom_Curve) curve_ptr = GC_MakeEllipse(Axis, r1, r2); 
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
    BRepBuilderAPI_MakeWire aWire(new_edge);

    TopoDS_Wire test_Wire = aWire.Wire();
    
    BRepOffsetAPI_ThruSections builder(CUBIT_TRUE, CUBIT_TRUE);
    builder.AddWire(test_Wire);
    if (r3 == 0)
    {
      gp_Pnt pt = gp_Pnt( 0.0, 0.0, hi);
      TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
      builder.AddVertex(theVertex);
    }
    else
    {
      gp_Pnt center2(0.0, 0.0,hi);
      gp_Ax2 Axis2(center2, main_dir, x_dir);
      Handle(Geom_Curve) curve_ptr = GC_MakeEllipse(Axis2, r3, r3*r2/r1);
      TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
      BRepBuilderAPI_MakeWire aWire(new_edge);
      TopoDS_Wire test_Wire = aWire.Wire();
      builder.AddWire(test_Wire);
    }
    builder.Build() ;
    S = TopoDS::Solid(builder.Shape());
  }

  else // cone
  {
    if(r1 == r3) //cylinder
      S = BRepPrimAPI_MakeCylinder(r1, hi);
    else
      S = BRepPrimAPI_MakeCone(r1, r3, hi);
  }

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::cylinder\n"
                "   Cannot create a cylinder for given radii.\n");
    return (BodySM*)NULL;
  }

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  if(body)
    CAST_TO(body, OCCBody)->move(0, 0, -hi/2.0);
  return body;
}

//===============================================================================
// Function   : torus
// Member Type: PUBLIC
// Description: create an OCC torus
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::torus( double r1, double r2 ) const
{
  if (r1 <= 0 || r2 <= 0)
    return (BodySM*) NULL;
 
  TopoDS_Solid S = BRepPrimAPI_MakeTorus(r1, r2);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::torus\n"
                "   Cannot create a torus for given radii.\n");
    return (BodySM*)NULL;
  }

  return CAST_TO(lump, OCCLump)->get_body();
}

//===============================================================================
// Function   : planar_sheet
// Member Type: PUBLIC
// Description: create an OCC planar_sheet with four vectors. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::planar_sheet ( const CubitVector& p1,
                                       const CubitVector& p2,
                                       const CubitVector& p3,
                                       const CubitVector& p4) const
{
  Point* point1 = make_Point(p1);
  Point* point2 = make_Point(p2);
  Point* point3 = make_Point(p3);
  Point* point4 = make_Point(p4);
  Curve * curve1 = make_Curve( point1, point2);
  if (curve1 == NULL)
	return (BodySM*) NULL;
  Curve * curve2 = make_Curve( point2, point3); 
  if (curve2 == NULL)
        return (BodySM*) NULL;
  Curve * curve3 = make_Curve( point3, point4);
  if (curve3 == NULL)
        return (BodySM*) NULL;
  Curve * curve4 = make_Curve( point4, point1);
  if (curve4 == NULL)
        return (BodySM*) NULL;
  DLIList<Curve*> curves;
  curves.append(curve1);
  curves.append(curve2);
  curves.append(curve3);
  curves.append(curve4);
  Surface* surface = make_Surface(PLANE_SURFACE_TYPE, curves);
  if (surface == NULL)
	return (BodySM*) NULL;
  return CAST_TO(surface,OCCSurface)->my_body();
}

//===============================================================================
// Function   : copy_body
// Member Type: PUBLIC
// Description: copy an OCC-based body
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::copy_body ( BodySM* bodyPtr ) const
{
  OCCBody* occ_body = CAST_TO(bodyPtr, OCCBody);
  if (!occ_body)
  {
     PRINT_ERROR("Cannot create an OCC bodySM from the given bodySM.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  TopoDS_CompSolid *theCS = occ_body->get_TopoDS_Shape();
  
  if (theCS == NULL) //sheet body
  {
    OCCShell* occ_shell = occ_body->shell();
    if (occ_shell != NULL)
    {
      TopoDS_Shell* shell = occ_shell->get_TopoDS_Shell();
      BRepBuilderAPI_Copy api_copy(*shell);
      TopoDS_Shape newShape = api_copy.ModifiedShape(*shell);
      TopoDS_Shell newShell = TopoDS::Shell(newShape);
      return OCCQueryEngine::instance()->populate_topology_bridge(newShell, CUBIT_TRUE)->my_body();
    }
 
    Surface* surface = make_Surface(occ_body->my_sheet_surface());
    if (surface == NULL)
    {
       PRINT_ERROR("Cannot create an OCC sheet bodySM from the given bodySM.\n");
       return (BodySM *)NULL;
    }

    OCCSurface* occ_surf = CAST_TO(surface, OCCSurface);
    return occ_surf->my_body();
  }

  if(OCCQueryEngine::instance()->OCCMap->IsBound(*theCS))
  {
    BRepBuilderAPI_Copy api_copy(*theCS);

    TopoDS_Shape newShape = api_copy.ModifiedShape(*theCS);

    TopoDS_CompSolid newCS = TopoDS::CompSolid(newShape);

    return OCCQueryEngine::instance()->populate_topology_bridge(newCS);
  }

  //single lump body
  Lump *lump = occ_body->lumps().get();
  TopoDS_Solid solid = *(CAST_TO(lump, OCCLump)->get_TopoDS_Solid());
  BRepBuilderAPI_Copy api_copy(solid);
  TopoDS_Shape newShape = api_copy.ModifiedShape(solid);
  TopoDS_Solid newSolid = TopoDS::Solid(newShape);
  lump = OCCQueryEngine::instance()->populate_topology_bridge(newSolid,
							CUBIT_TRUE);
  
  return CAST_TO(lump, OCCLump)->get_body();
}

//===============================================================================
// Function   : stitch_surfs
// Member Type: PUBLIC
// Description: stitch all surfs and try to make a shell body.
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::stitch_surfs(
                      DLIList<BodySM*>& surf_bodies,
                      BodySM*& new_body) const
{
  new_body = NULL;
  if (surf_bodies.size()==0)
    return CUBIT_SUCCESS;

  if (surf_bodies.size()==1)
  {
    new_body = surf_bodies.get();
    return CUBIT_SUCCESS;
  }

  TopoDS_Shape fuse;
  CubitStatus stat =  stitch_surfs(surf_bodies, fuse);
    
  DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->
     populate_topology_bridge(fuse);

  if (stat)
  {
    OCCBody* body = CAST_TO(tbs.get(), OCCBody);

    if (body)
      new_body = body ;
  }

  else
  {
    surf_bodies.clean_out();
    for(int i= 0; i<tbs.size(); i++)
      surf_bodies.append(CAST_TO(tbs.get_and_step(), OCCBody));
  }

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : stitch_surfs
// Member Type: PUBLIC
// Description: stitch all surfs and try to make a TopoDS_Shell .
//              called by stitch into surface body
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::stitch_surfs(
                      DLIList<BodySM*>& surf_bodies,
                      TopoDS_Shape& fuse) const
{
  if (surf_bodies.size() < 2)
    return CUBIT_SUCCESS;

  DLIList<TopoDS_Shape*> faces_to_stitch;
  for (int i = 0; i < surf_bodies.size(); i++)
  {
     BodySM * tool_body = surf_bodies.get_and_step();
     OCCBody* occ_body = CAST_TO(tool_body, OCCBody);
     OCCSurface* surface = occ_body->my_sheet_surface();
     OCCShell* shell = occ_body->shell();
     if (surface == NULL && shell == NULL)
     {
       PRINT_ERROR("Can't stitch non-sheet bodySM's. \n");
       return CUBIT_FAILURE;
     }

     delete occ_body;
     if (surface)
     {
       delete surface->my_shell();
       delete surface->my_lump();
       surface->set_shell(NULL);
       surface->set_lump(NULL);
       surface->set_body(NULL);

       TopoDS_Face* topods_face = surface->get_TopoDS_Face();
       if (topods_face != NULL)
         faces_to_stitch.append(topods_face);
     }
     else
     {
       delete shell->my_lump();
       shell->set_body(NULL);
       shell->set_lump(NULL);

       TopoDS_Shell* topods_shell = shell->get_TopoDS_Shell();
       if(topods_shell)
          faces_to_stitch.append(topods_shell);
     }
  }

  faces_to_stitch.reset();
  surf_bodies.reset();
  TopoDS_Shape* first_face  = faces_to_stitch.pop();

  TopoDS_Shape* second_face = NULL;
  for( int i = faces_to_stitch.size()-1; i >= 0; i --)
  {
     second_face = faces_to_stitch[i];
     BRepAlgoAPI_Fuse fuser(*second_face, *first_face);

     TopoDS_Shape new_shape ;
     for(int j = 0; j < 2; j++)
     {
       TopTools_IndexedMapOfShape M;
       if(j == 0)
         TopExp::MapShapes(*second_face, TopAbs_FACE, M);
       else
         TopExp::MapShapes(*first_face, TopAbs_FACE, M);
       for(int ii=1; ii<=M.Extent(); ii++)
       {
          TopoDS_Face face = TopoDS::Face(M(ii));
          TopTools_ListOfShape shapes;
          shapes.Assign(fuser.Modified(face));
          if (shapes.Extent() == 1)
          {
            new_shape = shapes.First();
            OCCSurface::update_OCC_entity(face, TopoDS::Face(new_shape), &fuser);
          }
          else if(fuser.IsDeleted(face) || shapes.Extent() > 1)
          {
            TopoDS_Face null_face;
            OCCSurface::update_OCC_entity(face, null_face, &fuser);
          }
        }
      } 

      fuse = fuser.Shape();
      first_face = &fuse;
  }

  TopExp_Explorer Ex;
  int count_shell = 0;
  for (Ex.Init(fuse, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
    count_shell++;

  if ( count_shell != 1)
  {
     PRINT_ERROR("Can't stitch all surfaces into one BodySM's. \n");
     return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : subtract
// Member Type: PUBLIC
// Description: subtract boolean operation on OCC-based bodies
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus     OCCModifyEngine::subtract(DLIList<BodySM*> &tool_body_list,
                                          DLIList<BodySM*> &from_bodies,
                                          DLIList<BodySM*> &new_bodies,
                                          bool imprint,
                                          bool keep_old) const
{
  // copy the bodies in case keep_old is true
  DLIList<TopoDS_Shape*> tool_bodies_copy;

  //for subtract function, tool-body has to be solid, 
  //otherwise it's just imprint
  DLIList<CubitBox*>* tool_boxes = new DLIList<CubitBox*>();
  DLIList<CubitBoolean> is_tool_volume;
  //keep the tool_body untouched
  CubitStatus stat = 
    get_shape_list(tool_body_list, tool_bodies_copy, is_tool_volume, CUBIT_TRUE, tool_boxes);

  if(!stat)
    return stat;

  stat = do_subtract(from_bodies, tool_bodies_copy, is_tool_volume,
                     tool_boxes, new_bodies, keep_old, imprint) ;

  //ok, we're done with all cuts, delete unnecessaries.
  CubitBoolean delete_tool_boxes = CUBIT_FALSE;
  if(tool_boxes->size() > 0)
    delete_tool_boxes = CUBIT_TRUE;
  while (tool_boxes->size())
    delete tool_boxes->pop();
  if(delete_tool_boxes)
    delete tool_boxes;
  while (tool_bodies_copy.size())
    delete tool_bodies_copy.pop();
  if(!keep_old) //delete tool_bodies
    OCCQueryEngine::instance()->delete_solid_model_entities(tool_body_list);
  return stat;
}

CubitStatus OCCModifyEngine::do_subtract(DLIList<BodySM*> &from_bodies,
                                      DLIList<TopoDS_Shape*> &tool_bodies_copy,
                                      DLIList<CubitBoolean> &is_tool_volume,
                                      DLIList<CubitBox*>* tool_boxes,
                                      DLIList<BodySM*> &new_bodies,
                                      bool keep_old,
                                      bool imprint) const
{
  DLIList<TopoDS_Shape*> from_bodies_copy;
  DLIList<CubitBoolean> is_volume;
  //get the from_bodies underling shapes
  CubitStatus stat = get_shape_list(from_bodies, from_bodies_copy, is_volume, keep_old);
  if(!stat)
  {
    for (int i = 0; i < tool_bodies_copy.size(); i++)
    {
       TopoDS_Shape* shape = tool_bodies_copy.get_and_step();
       delete shape;
    }
    tool_bodies_copy.clean_out();
    return CUBIT_FAILURE;
  } 

  //check that tool_bodies are all solid, shell and surface body can't be used
  //for subtracting solids.
  if(is_tool_volume.is_in_list(CUBIT_FALSE) && !is_volume.is_in_list(CUBIT_FALSE))
  {
     PRINT_WARNING("Surfaces or Shells can't be used to cut a solid.\n");
     while (tool_boxes->size())
       delete tool_boxes->pop();
     delete tool_boxes;
     while (tool_bodies_copy.size())
       delete tool_bodies_copy.pop();
     while (from_bodies_copy.size())
       delete from_bodies_copy.pop(); 
     return CUBIT_FAILURE;
  }

  int fraction_remaining = 100;

  // subtract the tool body from each body in the list

  CubitMessage* cmi = CubitMessage::instance(); 
  TopoDS_Shape*  from_shape = from_bodies_copy.get();
  DLIList<TopologyBridge*> tbs;
  for (int i = 0; i < from_bodies_copy.size(); i++)
  {
    CubitBoolean from_volume = is_volume.get_and_step();
    BodySM* from_body = from_bodies.get();
    CubitBox box1 = CAST_TO(from_body, OCCBody)->get_bounding_box();
    int count = 0;  //count for not preforming cut
    for(int j = 0; j < tool_bodies_copy.size(); j ++)
    {
      CubitBoolean tool_volume = is_tool_volume.get_and_step();
      if(tool_volume == CUBIT_FALSE && from_volume == CUBIT_TRUE)
      {
        PRINT_WARNING("Surfaces or Shells can't be used to cut a solid.\n");
        continue;
      }
      if (cmi->Interrupt())
      {
         PRINT_ERROR("Subtraction interrupted.  Aborting...\n");
         while (tool_boxes->size())
           delete tool_boxes->pop();
         delete tool_boxes;
         while (tool_bodies_copy.size())
            delete tool_bodies_copy.pop();
         while (from_bodies_copy.size())
            delete from_bodies_copy.pop();
         return CUBIT_FAILURE;
      }
      CubitBox tool_box = *tool_boxes->get_and_step();  
      if(!tool_box.overlap(TOL,box1))
      {
        count++;
        continue;
      } 
      TopoDS_Shape* tool_shape = tool_bodies_copy.get_and_step();

      //bodies overlap, proceed with the subtract
      BRepAlgoAPI_Cut cutter(*from_shape, *tool_shape);
      TopoDS_Shape cut_shape = cutter.Shape(); 

      //compare to see if the from_shape has gotten cut.
      CubitBoolean has_changed;
      check_operation(cut_shape, from_shape, is_volume[i], has_changed, &cutter, keep_old);

      CubitStatus stat;
      if(!has_changed && !from_shape->IsNull())
      {
        stat = CUBIT_FAILURE;
        //Add imprint code here 
        if(imprint)
        {
          stat = imprint_toposhapes(from_shape, tool_shape);
          if(!stat)
          {
            PRINT_ERROR("Can't do imprint operation on the body. \n");
            count++;
          }
          continue;
        }
      }
    }

    //ok, we're done with all cuts, construct new Body'
    if (count < tool_bodies_copy.size() )
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*from_shape);
    else
    {
      PRINT_INFO("The %d body did not change because cutting tools are not interscting with it.\n", i+1);
      from_bodies.change_to(NULL);
    }
    from_bodies.step();
    from_shape = from_bodies_copy.step_and_get();

    // done with this j iteration; write out count, if necessary
    if (from_bodies.size() * tool_bodies_copy.size() > 1)
    {
       int frac_done = (100 * (i+1)) / (from_bodies.size());
       if ((100 - frac_done) < fraction_remaining)
       {
          fraction_remaining = 100 - frac_done;
          PRINT_INFO("%d%% remaining.\n ", fraction_remaining+1);
       }
    }
  }

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }    

  //ok, we're done with all cuts, delete unnecessaries. 
  if(keep_old)
  {
    int size  = from_bodies_copy.size();
    for (int i = 0; i < size; i++)
    {
      TopoDS_Shape* shape = from_bodies_copy.pop();
      shape->Nullify();
      delete shape;
    }
  } 
  return CUBIT_SUCCESS; 
}

//===============================================================================
// Function   : imprint_toposhapes
// Member Type: PROTECTED
// Description: imprint boolean operation on OCC-based bodies.
//              from_shape must be TopoDS_Face or above, tool_shape can be
//              TopoDS_Edge or above. 
//              on_faces works only when tool_shape is an Edge, indicates that
//              those edges only imprint on the on_faces.
// Author     : Jane HU
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint_toposhapes(TopoDS_Shape*& from_shape, 
                                             TopoDS_Shape* tool_shape,
                                             DLIList<TopoDS_Face*>*on_faces)const
{
  int count = 0;   //number of imprinting
 
  //indicate if there's more faces to be imprinted
  CubitBoolean more_face = CUBIT_TRUE; 
  DLIList<TopoDS_Face*> list_for_delete;

  //list of face on from_shape that has been imprinted
  DLIList<TopoDS_Face*> from_faces; 
  while( more_face)
    {
      TopoDS_Face from_face,tool_face;
      TopoDS_Edge* common_edge = NULL;
      DLIList<TopoDS_Shape*> tool_faces_edges;
      TopTools_ListOfShape list_of_edges;
      BRepFeat_SplitShape splitor(*from_shape);
      CubitBoolean qualified = CUBIT_TRUE;
      if (tool_shape->TShape()->ShapeType() == TopAbs_EDGE)
	{
	  if(count == 1)
	    break;
      
	  DLIList<TopoDS_Face*> faces;
	  //need to delete TopoDS_Face* in faces
	  common_edge = find_imprinting_edge(*from_shape, TopoDS::Edge(*tool_shape),faces);
	  if (common_edge)
	    {
	      if (on_faces)
		qualified = CUBIT_FALSE;
	      for(int j = 0; j < faces.size(); j++)
		{
		  from_face = *faces.get();
		  for (int i = 0; on_faces && i < on_faces->size(); i++)
		    {
		      if (from_face.IsSame(*(on_faces->get_and_step())))
			{
			  qualified = CUBIT_TRUE; 
			  break;
			}
		    }
		  faces.get()->Nullify();
		  delete faces.get();
		  faces.step();
		}
	      if (qualified && (from_faces.size() == 0 || (from_faces.size() && !from_face.IsSame(*from_faces.get()))) )
		list_of_edges.Append(*common_edge);
	      else
		from_face.Nullify();
	      common_edge->Nullify();
	      delete common_edge;
	    }
	}
      else 
	{
	  TopOpeBRep_ShapeIntersector intersector;
	  intersector.InitIntersection(*from_shape, *tool_shape);

	  //find the intersecting edges and faces.
	  int max_edge = 0;
                     
	  for(; intersector.MoreIntersection(); )
	    {
	      TopoDS_Shape face1;
	      face1 = intersector.ChangeFacesIntersector().Face(1);
	      CubitBoolean has_imprinted = CUBIT_FALSE;
	      for (int j = 0; j < from_faces.size(); j++)
		{
		  TopoDS_Face* topo_face = from_faces.get_and_step();
		  if(face1.IsSame(*topo_face))
		    {
		      has_imprinted = CUBIT_TRUE;
		      break;
		    }
		}

	      if (has_imprinted == CUBIT_TRUE)
		{
		  intersector.NextIntersection();
		  continue;
		}
	      TopoDS_Shape edge_face;

	      edge_face = intersector.ChangeFacesIntersector().Face(2);
	      BRepAlgoAPI_Section section(face1, edge_face);

	      //intersection edges between face1 and edge_face
	      TopTools_ListOfShape temp_list_of_edges;
	      temp_list_of_edges.Assign(section.SectionEdges());
	      int num_edges = temp_list_of_edges.Extent();
  
	      CubitBoolean is_same = face1.IsSame(from_face);
	      CubitBoolean is_same_tool = CUBIT_FALSE;
	      for (int j = 0; j < tool_faces_edges.size(); j++)
		{
		  TopoDS_Shape* topo_face_edge = tool_faces_edges.get_and_step();
		  if(edge_face.IsSame(*topo_face_edge))
		    {
		      is_same_tool = CUBIT_TRUE;
		      break;
		    }
		}
	      if (max_edge < num_edges )
		{
		  list_of_edges.Assign(temp_list_of_edges);  
		  max_edge =  num_edges ;
		  from_face = TopoDS::Face(face1);
		  TopoDS_Shape* topo_shape = new TopoDS_Shape(edge_face);
		  DLIList<TopoDS_Shape*> shape_list;
		  for(int iii = 0; iii < tool_faces_edges.size(); iii++)
		    {
		      int size = shape_list.size();
		      shape_list.append_unique(tool_faces_edges.get_and_step());
		      if (size < shape_list.size())
			{
			  shape_list.last();
			  shape_list.get()->Nullify();
			  delete shape_list.get();
			}
		    }
		  tool_faces_edges.clean_out();
		  for(int i = 0 ; i < num_edges; i++)
		    //later has to use it num_edges times 
		    tool_faces_edges.append(topo_shape);
		}
	      else if(num_edges == max_edge && is_same && !is_same_tool) 
		//multi tool faces cut the same face
		{
		  TopTools_ListIteratorOfListOfShape Itor, temp_Itor;
		  temp_Itor.Initialize(temp_list_of_edges);
		  for(; temp_Itor.More(); temp_Itor.Next())
		    {
		      TopoDS_Edge temp_edge = TopoDS::Edge(temp_Itor.Value());
		      Itor.Initialize(list_of_edges);
		      CubitBoolean same_edge = CUBIT_FALSE;
              
		      GProp_GProps myProps1;
		      BRepGProp::LinearProperties(temp_edge, myProps1);
		      gp_Pnt center1 = myProps1.CentreOfMass();
		      for(; Itor.More(); Itor.Next())
			{
			  TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
			  GProp_GProps myProps2;
			  BRepGProp::LinearProperties(edge, myProps2);
			  gp_Pnt center2 = myProps2.CentreOfMass();
			  if(center1.IsEqual(center2, TOL))
			    {
			      same_edge = CUBIT_TRUE;
			      break;
			    }
			}
		      if(!same_edge)
			{
			  list_of_edges.Append(temp_edge);
			  TopoDS_Shape* topo_shape = new TopoDS_Shape(edge_face);
			  tool_faces_edges.append(topo_shape);
			}
		    }//end 'for'
		}//end  'else if'
              intersector.NextIntersection();
	    } //end 'for'
	}//end 'else'
      if (from_face.IsNull())
	{
	  more_face = CUBIT_FALSE;
	  DLIList<TopoDS_Shape*> shape_list;
	  for(int iii = 0; iii < tool_faces_edges.size(); iii++)
	    {
	      int size = shape_list.size();
	      shape_list.append_unique(tool_faces_edges.get_and_step());
	      if (size < shape_list.size())
		{
		  shape_list.last();
		  shape_list.get()->Nullify();
		  delete shape_list.get();
		}
	    }
	  tool_faces_edges.clean_out();

	  for (int iii=0; iii < from_faces.size(); iii++)
	    {
	      TopoDS_Face* topo_face = from_faces.get_and_step();
	      topo_face->Nullify();
	      delete topo_face;
	    }
        
	  for (int iii=0; iii < list_for_delete.size(); iii++)
	    {
	      TopoDS_Face* topo_face = list_for_delete.get_and_step();
	      topo_face->Nullify();
	      delete topo_face;
	    }
	  continue;
	}
  
      TopTools_ListIteratorOfListOfShape Itor;

      //list_of_edges is the intersection edges on tool_faces_edges 
      Itor.Initialize(list_of_edges);
      int total_edges = list_of_edges.Extent();
      DLIList<Curve*> curve_list;
      CubitBoolean topo_changed = CUBIT_FALSE;
      tool_faces_edges.reset();

      //check to see if the intersection edge is:
      //1. on from_face: if not, skip it
      //2. overlap with from_edges : if not, add the edge for splitting face
      //3. if overlap, is it the same edge:if not add it for splitting edge
      // if yes, skip it too

      Surface* face = NULL;
      if (OCCQueryEngine::instance()->OCCMap->IsBound(from_face))
	{
	  int i = OCCQueryEngine::instance()->OCCMap->Find(from_face);
	  face = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(i))->second;
	}
      else
        face = OCCQueryEngine::instance()->populate_topology_bridge(from_face);
      OCCSurface* occ_face = CAST_TO(face, OCCSurface);

      DLIList<Curve*> common_curves;
      for(; Itor.More(); Itor.Next())
	{
	  TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
	  //copy the edge for imprinting.
	  BRepBuilderAPI_Copy api_copy(edge);
	  TopoDS_Shape newShape = api_copy.ModifiedShape(edge);
	  edge = TopoDS::Edge(newShape);
	  Curve* curve = NULL;
	  curve = OCCQueryEngine::instance()->populate_topology_bridge(edge); 
	  if(curve)
	    common_curves.append(curve);
	}

      DLIList<DLIList<TopoDS_Edge*>*> temp_edge_lists;
      list_of_edges.Clear(); 
      if (common_curves.size() >= 1)
	{
	  sort_curves(common_curves, temp_edge_lists);
	  DLIList<TopoDS_Edge*>* edge_list;
	  int size = temp_edge_lists.size();
	  for(int i = 0; i < size; i++)
	    {
	      edge_list = temp_edge_lists.pop();
	      //make sure the copied edges are sharing vertices.
	      BRepBuilderAPI_MakeWire myWire;
	      edge_list->reset();
	      for(int i = 0; i < edge_list->size(); i++)
		{
		  TopoDS_Edge e = *(edge_list->get_and_step());
		  myWire.Add(e);
		}
	      TopoDS_Wire wire = myWire.Wire();
	      BRepTools_WireExplorer Ex(wire); 
	      for(; Ex.More(); Ex.Next())
		list_of_edges.Append(Ex.Current());
              edge_list->clean_out();
              delete edge_list;
	    }
	}
      for(Itor.Initialize(list_of_edges); Itor.More(); Itor.Next())
	{
	  TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
	  //check to see if the intersection edge is on from_face
	  TopExp_Explorer Ex;
	  CubitBoolean added = CUBIT_FALSE;
	  CubitBoolean skipped = CUBIT_FALSE;
	  GProp_GProps myProps1;
	  BRepGProp::LinearProperties(edge, myProps1);
	  double d1 = myProps1.Mass();
	  gp_Pnt pt = myProps1.CentreOfMass();
	  CubitVector p(pt.X(), pt.Y(), pt.Z());

	  CubitVector point_on_surf;
	  occ_face->closest_point_trimmed(p, point_on_surf);
	  if(p.distance_between(point_on_surf) > TOL) //edge not on from_face
	    {
	      skipped = CUBIT_TRUE;
	      total_edges--;
	    }

	  else 
	    {
	      for (Ex.Init(from_face, TopAbs_EDGE); Ex.More(); Ex.Next())
		{
		  //check if the edge is on from_face edges, add such edge on existing
		  //edge to split it.
		  TopoDS_Edge from_edge = TopoDS::Edge(Ex.Current());
         
		  GProp_GProps myProps2;
		  BRepGProp::LinearProperties(from_edge, myProps2);
		  double d2 = myProps2.Mass();
		  Curve* curve = NULL;
		  if (OCCQueryEngine::instance()->OCCMap->IsBound(from_edge))
		    {
		      int i = OCCQueryEngine::instance()->OCCMap->Find(from_edge);
		      curve = (OCCCurve*)(OCCQueryEngine::instance()->OccToCGM->find(i))->second;
		    }
		  else
		    curve = OCCQueryEngine::instance()->populate_topology_bridge(from_edge);
		  OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
   
		  CubitPointContainment pc = CUBIT_PNT_OFF;
		  if(occ_curve->geometry_type() == STRAIGHT_CURVE_TYPE) 
		    {
		      pc = occ_curve->point_containment(p);

		      if(pc == CUBIT_PNT_ON) //overlap
			{
			  //check if they are the same edge, so don't need to be split
			  //the overlapped edges are considered the same if they have the
			  //same length
			  if((d2 - d1) > TOL)
			    {
			      added = CUBIT_TRUE;
			      splitor.Add(edge, from_edge);
			    }
			  else
			    skipped = CUBIT_TRUE;
			  total_edges--;
			  break;
			}
		    }
		  else if(list_of_edges.Extent() == 1 && (d2 - d1) <= TOL)
		    skipped = CUBIT_TRUE;
		} 
	      if(list_of_edges.Extent() == 1 && !skipped) 
		{
		  added = CUBIT_TRUE;
		  Curve* curve = OCCQueryEngine::instance()->populate_topology_bridge(edge);
		  curve_list.append(curve); 
		}
	    } 
	  if(added)
	    {
	      topo_changed = CUBIT_TRUE;
	      continue;
	    }
	  if (!skipped)
	    {
	      //check if edge's inside from_face by checking bounding boxes  
	      BRepAdaptor_Curve acurve(edge);
	      BRepAdaptor_Surface asurface( from_face);
	      Bnd_Box aBox_edge, aBox_face;
	      BndLib_Add3dCurve::Add(acurve, Precision::Approximation(), aBox_edge);
	      BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox_face);
	      double min[3], max[3];
	      aBox_edge.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
	      CubitBox aBox_e(min, max);
	      aBox_face.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
	      CubitBox aBox_f(min, max);
	      if (aBox_e <= aBox_f)
		{
		  Curve* curve = OCCQueryEngine::instance()->populate_topology_bridge(edge);
		  curve_list.append(curve);
		}
	    }
	}

      DLIList<DLIList<TopoDS_Edge*>*> edge_lists;
      if (total_edges >= 1)
	{      
	  CubitStatus stat = CUBIT_SUCCESS;
          if(curve_list.size() > 0)
	    stat = sort_curves(curve_list, edge_lists); 
          else
	    {
	      TopoDS_Face* topo_face = new TopoDS_Face(from_face);
	      from_faces.append(topo_face);
	      continue;
	    }
	  DLIList<TopoDS_Edge*>* edge_list;
          int size = edge_lists.size();
          for (int iii = 0; iii < size; iii++)
	    {
	      edge_list = edge_lists.pop();

	      //check if the edges make a split of the surface, error out if it's
	      //just a scar on the surface
	      int count_intersection = 0;
	      if (stat == CUBIT_FAILURE) //Open wire
		{
		  count_intersection = check_intersection(edge_list, from_face);
            
		  if (count_intersection == 1 )
		    PRINT_WARNING("Cant make a scar on existing face without splitting it. \n");
		} 
	      if (stat || count_intersection == 2)
		{
		  BRepBuilderAPI_MakeWire myWire;
		  edge_list->reset(); 
		  for(int i = 0; i < edge_list->size(); i++)
		    {
		      TopoDS_Edge e = *(edge_list->get_and_step());
		      myWire.Add(e); 
		    }
		  splitor.Add(myWire.Wire(),from_face);
		  topo_changed = CUBIT_TRUE; 
                  edge_list->clean_out();
                  delete edge_list;
                  for(int i = 0; i <edge_lists.size(); i++)
                  {
                    edge_list = edge_lists.get_and_step();
                    edge_list->clean_out();
                    delete edge_list;
                  }
		  break;
		}
	      edge_list->clean_out();
	      delete edge_list;
	    }
	} 
      if(topo_changed)
	{
	  splitor.Build();
	  topo_changed = CUBIT_FALSE;
	  if(splitor.IsDone())
	    {
	      //take care of on_faces list first:after operation, the on_faces
	      // will have at least one face changed, update the pointer.
	      if (on_faces)
		{
		  for(int k = 0; k < on_faces->size(); k++)
		    {
		      TopoDS_Face* compare_face = on_faces->get_and_step();
		      if(from_face.IsSame(*compare_face))
			{
			  on_faces->remove(compare_face);
			  TopTools_ListOfShape shapes;
			  shapes.Assign(splitor.Modified(from_face)); 
			  while(shapes.Extent() > 0)
			    {
			      TopoDS_Face* face = 
				new TopoDS_Face(TopoDS::Face(shapes.First())); 
			      shapes.RemoveFirst();
			      on_faces->append(face);
			      list_for_delete.append(face);
			    }
			}
		    }
		}

	      TopoDS_Shape new_from_shape = splitor.Shape();
	      if(from_shape->TShape()->ShapeType() == TopAbs_COMPSOLID)
		{
		  TopoDS_CompSolid old_csolid = TopoDS::CompSolid(*from_shape);
		  OCCBody::update_OCC_entity(old_csolid, new_from_shape, &splitor);
		  from_shape->Nullify();
		  delete from_shape;
		  from_shape = new TopoDS_Shape(new_from_shape);
		}

	      else if(from_shape->TShape()->ShapeType() == TopAbs_SOLID)
		{
		  TopoDS_Solid old_solid = TopoDS::Solid(*from_shape);
		  OCCLump::update_OCC_entity(old_solid, new_from_shape, &splitor);
		  from_shape->Nullify();
		  delete from_shape;
		  from_shape = new TopoDS_Shape(new_from_shape);
		}
	      else if(from_shape->TShape()->ShapeType() == TopAbs_SHELL)
		{
		  TopoDS_Shell old_shell = TopoDS::Shell(*from_shape);
		  OCCShell::update_OCC_entity(old_shell,new_from_shape, &splitor);
		  from_shape->Nullify();
		  delete from_shape;
		  from_shape = new TopoDS_Shape(new_from_shape);
		}
	      else if(from_shape->TShape()->ShapeType() == TopAbs_FACE)
		{
		  TopoDS_Face old_face = TopoDS::Face(*from_shape);
		  OCCSurface::update_OCC_entity(old_face,new_from_shape, &splitor);
		  from_shape->Nullify();
		  delete from_shape;
		  from_shape = new TopoDS_Shape(new_from_shape);
		}

	      TopTools_ListOfShape shapes;
	      for(int i = 0; i < from_faces.size(); i++)
		{
		  TopoDS_Face* topo_face = from_faces.get();
		  shapes.Assign(splitor.Modified(*topo_face));
		  topo_face = new TopoDS_Face(TopoDS::Face(shapes.First()));
		  from_faces.get()->Nullify();
		  delete from_faces.get();
		  from_faces.change_to(topo_face);
		  from_faces.step();
		} 
	      count++;
	    }
	}
      else
	{
	  TopoDS_Face* topo_face = new TopoDS_Face(from_face);
	  from_faces.append(topo_face);
	} 
    }
  
  
  TopExp_Explorer Ex;
  int num_face = 0;
  for (Ex.Init(*from_shape, TopAbs_FACE); Ex.More(); Ex.Next())
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      num_face++;
    }
  
#ifdef DEBUG  
  PRINT_INFO("Total %d cuts performed, with from_shape having %d faces.\n", count, num_face); 
#endif

  if (count > 0)
    return CUBIT_SUCCESS;
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : find_imprinting_edge
// Member Type: PROTECTED
// Description: imprint boolean operation on OCC-based bodies.
//              from_shape must be TopoDS_Face or above, tool_shape must be
//              TopoDS_Edge.
// Author     : Jane HU
// Date       : 05/08
//===============================================================================
TopoDS_Edge* OCCModifyEngine::find_imprinting_edge(TopoDS_Shape& from_shape,
                                        TopoDS_Edge& tool_shape,
                                        DLIList<TopoDS_Face*>& from_faces)const
{
  TopoDS_Edge* edge = NULL;
  //list of face on from_shape that has been imprinted
  from_faces.clean_out();
  TopExp_Explorer Ex;
  for (Ex.Init(from_shape, TopAbs_FACE); Ex.More(); Ex.Next())
  {
    TopoDS_Face face = TopoDS::Face(Ex.Current());
    BRepAlgoAPI_Common intersector(face, tool_shape);
    TopTools_ListOfShape shapes;
    shapes.Assign(intersector.Modified(tool_shape));
    if (shapes.IsEmpty())
      continue;
    if ( shapes.Extent() > 1)
    {
      PRINT_ERROR("Edge has multiple intersection with the shape, make it simpler. \n");
      continue;
    }
    if (shapes.First().TShape()->ShapeType() != TopAbs_EDGE)
      continue;
    if(edge == NULL)
    {
      TopoDS_Edge common_edge = TopoDS::Edge(shapes.First());
       BRepBuilderAPI_Copy api_copy(common_edge);
       TopoDS_Shape newShape = api_copy.ModifiedShape(common_edge);
       edge = new TopoDS_Edge(TopoDS::Edge(newShape));
    }
    from_faces.append(new TopoDS_Face(face));
  }
  return edge;
}

int OCCModifyEngine::check_intersection(DLIList<TopoDS_Edge*>* edge_list,
 				        TopoDS_Face from_face)const
{
  int  count_intersection = 0;
  //Consider if edge_list has only one edge, and it intersects the from_face
  //at two different places.
  CubitBoolean double_check = CUBIT_FALSE;
  if (edge_list->size() == 1)
    double_check = CUBIT_TRUE;

  gp_Pnt pt1(0,0,0), pt2(0,0,0); 
  gp_Pnt intsec_pnt[2] = {pt1, pt2} ;
  for(int j = 0; j < edge_list->size(); j++)
  {
    TopoDS_Edge* edge = edge_list->get_and_step();
    BRepAdaptor_Curve acurve(*edge);
    double lower_bound = acurve.FirstParameter();
    double upper_bound = acurve.LastParameter();
    TopExp_Explorer Ex;
    gp_Pnt p1(0,0,0), p2(0,0,0); 
    gp_Pnt newP[2] = {p1, p2 };
    for (Ex.Init(from_face, TopAbs_EDGE); Ex.More(); Ex.Next())
    {
      TopoDS_Edge from_edge = TopoDS::Edge(Ex.Current());
      BRepAdaptor_Curve acurve2(from_edge);
      double lower_bound2 = acurve2.FirstParameter();
      double upper_bound2 = acurve2.LastParameter();
      BRepExtrema_DistShapeShape distShapeShape(*edge, from_edge);
      CubitBoolean qualified[2] = {CUBIT_FALSE, CUBIT_FALSE};
      if (distShapeShape.IsDone() && distShapeShape.Value() < TOL)
      {
        newP[0] = distShapeShape.PointOnShape1(1);
        if (double_check && distShapeShape.NbSolution() == 2)
          newP[1] = distShapeShape.PointOnShape1(2);
        double newVal[2];
        for(int j =0; j < distShapeShape.NbSolution(); j++)
        {
          Extrema_ExtPC ext(newP[j], acurve, Precision::Approximation());
          if (ext.IsDone() && (ext.NbExt() > 0)) {
            for ( int i = 1 ; i <= ext.NbExt() ; i++ ) {
              if ( ext.IsMin(i) ) {
        	newVal[j] = ext.Point(i).Parameter();
		if ((newVal[j]-lower_bound) >= -TOL && 
                    (upper_bound - newVal[j]) >= -TOL)
		{
		  qualified[j] = CUBIT_TRUE;
		  break;
		}
	      }
	    }
          }
        }
        for(int j = 0; j < distShapeShape.NbSolution(); j++)
        {
	  if (qualified[j])
	  {
	    qualified[j] = CUBIT_FALSE;
	    Extrema_ExtPC ext(newP[j], acurve2, Precision::Approximation());
	    double newVal;
	    if (ext.IsDone() && (ext.NbExt() > 0)) {
	      for ( int i = 1 ; i <= ext.NbExt() ; i++ ) {
	      	if ( ext.IsMin(i) ) {
		  newVal = ext.Point(i).Parameter();
                  if ((newVal-lower_bound2) >= -TOL &&
                      (upper_bound2 - newVal) >= -TOL)
		  {
		    qualified[j] = CUBIT_TRUE;
		    break;
		  }
		}
    	      }
	    }
	  }
        }
        for(int k = 0; count_intersection < 3 && k < distShapeShape.NbSolution(); k++)
        {
          if (qualified[k])
            count_intersection++;
          intsec_pnt[count_intersection-1] = newP[k];
          if (count_intersection == 2)
          {
            //make sure the two intersect point are not the same one 
            if (intsec_pnt[0].IsEqual(intsec_pnt[1], TOL))
              count_intersection--;
          }
        }
      }
      if (count_intersection == 2)
	break;
    } //for loop
  }
  return count_intersection;
}
//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on OCC-based bodies
// Author     : Jane Hu
// Date       : 04/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(BodySM* BodyPtr1, BodySM* BodyPtr2,
                                         BodySM*& newBody1, BodySM*& newBody2,
                                         bool  keep_old) const
{
  newBody1 = NULL;
  newBody2 = NULL;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_volume;
  
  DLIList<BodySM*> bodysm_list;
  bodysm_list.append(BodyPtr1);
  bodysm_list.append(BodyPtr2);
  
  CubitStatus stat = get_shape_list(bodysm_list,shape_list,is_volume,keep_old);

  if(!stat)
    return stat;

  TopoDS_Shape* shape1 = shape_list.get();
  TopoDS_Shape* shape2 = shape_list.step_and_get();
  DLIList<TopologyBridge*> tbs;
  stat = imprint_toposhapes(shape1, shape2);
  if(stat)
  {
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape1); 
    newBody1 = CAST_TO(tbs.get(),BodySM);
  }

  else if(!stat && keep_old)
  {
    delete shape1;
    PRINT_INFO("There's no imprint on the first body.\n");
    newBody1 = BodyPtr1;
  }

  tbs.clean_out();
  stat = imprint_toposhapes(shape2, shape1);
  if(stat)
  {
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape2);
    newBody2 = CAST_TO(tbs.get(),BodySM);     
  }
  
  else if(!stat && keep_old)
  {
    delete shape2;
    PRINT_INFO("There's no imprint on the second body.\n");
    newBody2 = BodyPtr2;
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_shape_list
// Member Type: PRIVATE
// Description: get the TopoDS_Shape list for imprinting use. 
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus OCCModifyEngine::get_shape_list(DLIList<BodySM*>& BodySM_list, 
                                         DLIList<TopoDS_Shape*>& shape_list,
                                         DLIList<CubitBoolean>& is_volume,
                                         bool  keep_old,
                                         DLIList<CubitBox*>* b_boxes) const
{
  OCCBody* occ_body = NULL;
  shape_list.clean_out();
  is_volume.clean_out();
  CubitStatus stat = CUBIT_SUCCESS;
  for(int i = 0; i <BodySM_list.size(); i++)
  {
    occ_body = CAST_TO(BodySM_list.get_and_step(), OCCBody);
    if (!occ_body)
      continue;

    OCCSurface* surface = occ_body->my_sheet_surface();
    OCCShell*   shell = occ_body->shell();
    is_volume.append( CUBIT_TRUE);

    if(b_boxes)
    {
      CubitBox *tool_box = new CubitBox(occ_body->get_bounding_box());
      b_boxes->append(tool_box);
    }

    if(surface)
    {
      TopoDS_Face* topo_face = surface->get_TopoDS_Face();
      if(!topo_face)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*topo_face);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_face);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(topo_face);
      is_volume.last();
      is_volume.change_to( CUBIT_FALSE);
    }
    else if(shell)
    {
      TopoDS_Shell* topo_shell = shell->get_TopoDS_Shell();
      if(!topo_shell)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*topo_shell);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_shell);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(topo_shell);
      is_volume.change_to( CUBIT_FALSE);
    }

    else
    {
      DLIList<Lump*> lumps = occ_body->lumps();
      if (lumps.size() > 1)
      {
        PRINT_ERROR("Can't do boolean operation on CompSolid types. \n");
        stat = CUBIT_FAILURE;
        break;
      }

      TopoDS_Solid* solid = CAST_TO(lumps.get(), OCCLump)->get_TopoDS_Solid();
      if(!solid)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*solid);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(solid);
    }
  }
  if(!stat)
  {   
    for (int i = 0; keep_old && i < shape_list.size(); i++)
    {
          TopoDS_Shape* shape = shape_list.get_and_step();
          delete shape;
    }
    shape_list.clean_out();
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : imprint multiple bodies at once
// Member Type: PUBLIC
// Description: imprint boolean operation on OCC-based bodies
// Author     : Jane HU
// Date       : 04/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint(DLIList<BodySM*> &from_body_list ,
                                     DLIList<BodySM*> &new_from_body_list,
                                     bool keep_old,
                                     DLIList<TopologyBridge*>* new_tbs,
                                     DLIList<TopologyBridge*>* att_tbs) const
{
  CubitStatus success = CUBIT_SUCCESS;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_vo;
  CubitStatus stat = get_shape_list(from_body_list, shape_list, is_vo,keep_old);

  if(!stat)
    return stat;
 
  int size = shape_list.size();
  // total number of imprints to be done
  int total_imprints = (size * (size -1))/2;

  if( size > 2 )
  {
     char message[128];
     sprintf(message, "Imprinting %d OCC Bodies", from_body_list.size() ); 
     AppUtil::instance()->progress_tool()->start(0, total_imprints, message);
  }

  for(int i = 0; i < size; i++)
  {
    TopoDS_Shape* shape1 = shape_list[i];
    CubitBoolean modified = CUBIT_FALSE;
    for(int j = i+1; j < size+i; j ++)
    {
       if (CubitMessage::instance()->Interrupt())
       {
          success = CUBIT_FAILURE;
          break;
       }

       TopoDS_Shape* shape2 = shape_list[j%size];
       DLIList<TopologyBridge*> tbs;
       CubitStatus stat = imprint_toposhapes(shape1, shape2);
       if(stat)
          modified = CUBIT_TRUE; 
    }
    if(modified)
    {
      DLIList<TopologyBridge*> tbs;
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape1);
      new_from_body_list.append(CAST_TO(tbs.get(),BodySM));
    }
    shape_list.reset();
    if( size > 2 )
      AppUtil::instance()->progress_tool()->step();
  }

  if( size > 2 )
    AppUtil::instance()->progress_tool()->end();

  if( CubitMessage::instance()->Interrupt() )
        PRINT_INFO("Imprint aborted.\n");

  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint curves onto body_list
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                           DLIList<Curve*> &ref_edge_list,
                                           DLIList<BodySM*>& new_body_list,
                                           bool keep_old,
                                           bool show_messages) const
{
  CubitStatus success = CUBIT_SUCCESS;
  DLIList<TopoDS_Shape*> shape_list, tool_shapes;
  DLIList<CubitBoolean> is_vo;
  CubitStatus stat = get_shape_list(body_list, shape_list, is_vo, keep_old);
  if (!stat)
    return stat;

  int size = ref_edge_list.size();
  // total number of imprints to be done

  if( size > 2 && show_messages)
  {
     char message[128];
     sprintf(message, "Imprinting %d OCC Bodies", body_list.size() );
     AppUtil::instance()->progress_tool()->start(0, size, message);
  }
  for (int i = 0; i < ref_edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(ref_edge_list.get_and_step(), OCCCurve) ;
    if (!curve)
      continue;

    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
    if (edge->IsNull())
      continue;
    
    if (CubitMessage::instance()->Interrupt())
    {
       success = CUBIT_FAILURE;
       break;
    }

    for(int j = 0; j < shape_list.size(); j ++)
    {
      TopoDS_Shape* shape = shape_list.get();
        
      stat = imprint_toposhapes(shape, (TopoDS_Shape*)edge);
      if (stat)
        shape_list.change_to(shape);
      shape_list.step();
      body_list.step();
    }

    if( size > 2 )
      AppUtil::instance()->progress_tool()->step();
  }   

  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
    new_body_list.append(CAST_TO(tbs.get(),BodySM));
  }

  if( size > 2 )
    AppUtil::instance()->progress_tool()->end();

  if( CubitMessage::instance()->Interrupt() )
        PRINT_INFO("Imprint aborted.\n"); 
  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: to be consistante with Acis imprint.
//              The surfaces must be part of a body, but the curves 
//              just have to be valid OCC edge.
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint( DLIList<Surface*> &ref_face_list,
                                      DLIList<Curve*> &edge_list,
                                      DLIList<BodySM*>& new_body_list,
                                      bool keep_old ) const
{
  DLIList<TopoDS_Face*> face_list;
  DLIList<TopoDS_Shape*> shape_list;
 
  face_edge_imprint(ref_face_list, edge_list, face_list, shape_list, keep_old);

  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    if (shape->TShape()->ShapeType() == TopAbs_COMPSOLID)
    {
      if(!OCCQueryEngine::instance()->OCCMap->IsBound(*shape)) 
      {
        TopExp_Explorer Ex;
        for (Ex.Init(*shape, TopAbs_SOLID);Ex.More(); Ex.Next())
        {
          TopoDS_Shape subshape = Ex.Current();
          tbs += OCCQueryEngine::instance()->populate_topology_bridge(subshape);
          new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
        }
      }
    }
    else
    {
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
      new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
    }
  }

  if (keep_old)
  {
    for(int i = 0; i < face_list.size(); i++)
    {
      TopoDS_Face* face = face_list.get();
      face->Nullify();
      delete face;
    }
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : face_edge_imprint
// Member Type: PRIVATE
// Description: to be consistante with Acis imprint.
//              The surfaces must be part of a body, but the curves
//              just have to be valid OCC edge.
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus 
OCCModifyEngine::face_edge_imprint( DLIList<Surface*> &ref_face_list,
                                    DLIList<Curve*> &edge_list,
                                    DLIList<TopoDS_Face*>& face_list,
                                    DLIList<TopoDS_Shape*>& shape_list,
                                    bool keep_old ) const
{
  for(int i = 0; i <ref_face_list.size(); i++)
  {
    OCCSurface* surface = CAST_TO(ref_face_list.get_and_step(), OCCSurface);
    if(!surface)
      continue;

    TopoDS_Face* topo_face = surface->get_TopoDS_Face();
    face_list.append(topo_face);

    if(surface->my_shell() && !surface->my_body())//shell body
      shape_list.append(surface->my_shell()->get_TopoDS_Shell());
    else if(surface->my_body()) //a sheet body
      shape_list.append(topo_face);
    else
    {
      DLIList<OCCBody*> bodies;
      surface->get_bodies(bodies);
      if(bodies.size() != 1)
      {
        PRINT_ERROR("Can't find the corresponding manifold solid body.\n");
        return CUBIT_FAILURE;
      }
      shape_list.append_unique(bodies.get()->get_TopoDS_Shape());
    }
  }

  if(keep_old)
  {
    for(int i = 0; i < shape_list.size(); i++)
    {
      TopoDS_Shape* shape = shape_list.get();
      BRepBuilderAPI_Copy api_copy(*shape);
      TopoDS_Shape newShape = api_copy.ModifiedShape(*shape);
      TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
      for(int j = 0; j < face_list.size(); j++)
      {
        TopoDS_Face* face = face_list.get();
        TopExp_Explorer Ex;
        for (Ex.Init(*shape, TopAbs_FACE); Ex.More(); Ex.Next())
        {
          if(face->IsSame(Ex.Current()))
          {
            face = new TopoDS_Face(TopoDS::Face(api_copy.ModifiedShape(*face)));
            face_list.change_to(face);
          }
        }
        face_list.step();
      }
      shape_list.change_to(Shape1);
      shape_list.step();
    }
  }

  for (int i = 0; i < edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(edge_list.get_and_step(), OCCCurve) ;
    if (!curve)
      continue;

    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
    if (edge->IsNull())
      continue;

    for(int j = 0; j < shape_list.size(); j ++)
    {
      TopoDS_Shape* shape = shape_list.get_and_step();
      imprint_toposhapes(shape, (TopoDS_Shape*)edge, &face_list);
    }
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: To be consistent with  AcisModifyEngine, althought it's hard 
//              to have a GUI interface for users to input. All surface must
//              on the same body. 
// Author     : Jane HU 
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint( DLIList<Surface*>& surface_list,
                                   DLIList<DLIList<Curve*>*>& curve_lists_list,
                                   BodySM*& new_body,
                                   bool keep_old ) const
{
  DLIList<TopoDS_Face*> face_list;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<TopoDS_Shape*> shape_list_all;
  
  assert (surface_list.size() == curve_lists_list.size());

  for(int j = 0; j < surface_list.size(); j++)
  {
    Surface* surface = surface_list.get_and_step();
    DLIList<Surface*> ref_face_list;
    ref_face_list.append(surface);
    DLIList<Curve*> *edge_list = curve_lists_list.get_and_step();
    face_edge_imprint(ref_face_list, *edge_list, face_list, shape_list, keep_old);

    for(int i = 0; i < shape_list.size(); i++)
    {
      TopoDS_Shape* shape = shape_list.get_and_step();
      shape_list_all.append_unique(shape);
    }
    shape_list.clean_out();

    if (keep_old)
    {
      for(int i = 0; i < face_list.size(); i++)
      {
        TopoDS_Face* face = face_list.get();
        face->Nullify();
        delete face;
      }
    }

    face_list.clean_out();
  }

  DLIList<BodySM*> new_body_list;
  shape_to_bodySM(shape_list, new_body_list);
  
  if (new_body_list.size() == 1)
  {
    new_body = new_body_list.get();
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : shape_to_bodySM
// Member Type: PRIVATE
// Description: After imprint, update shape list to bodySM_list
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
void OCCModifyEngine::shape_to_bodySM( DLIList<TopoDS_Shape*> shape_list,
                                       DLIList<BodySM*>& new_body_list)const
{
  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    if (shape->TShape()->ShapeType() == TopAbs_COMPSOLID)
    {
      if(!OCCQueryEngine::instance()->OCCMap->IsBound(*shape))
      {
        TopExp_Explorer Ex;
        for (Ex.Init(*shape, TopAbs_SOLID);Ex.More(); Ex.Next())
        {
          TopoDS_Shape subshape = Ex.Current();
          tbs += OCCQueryEngine::instance()->populate_topology_bridge(subshape);
          new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
        }
      }
    }
    else
    {
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
      new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
    }
  }
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: Imprints locations to bodies (for splitting curves, there's
//              no known ways to put hard points on surfaces in OCC, so I just
//              add free_vertex on OCCSurface definition, mesh should look on
//              this structure).   
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                           DLIList<CubitVector*> &vector_list,
                                           DLIList<BodySM*>& new_body_list,
                                           bool keep_old,
                                           DLIList<TopologyBridge*>*,
                                           DLIList<TopologyBridge*>* ) const
{
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_vo;
  CubitStatus stat = get_shape_list(body_list, shape_list, is_vo, keep_old);
  if(!stat)
    return stat;

  for (int i = 0; i < body_list.size(); i++)
  {
    OCCBody* body = CAST_TO(body_list.get_and_step(), OCCBody);
    TopoDS_Shape* from_shape = shape_list.get_and_step();
    if (!body)
      continue;
    DLIList<OCCSurface*> surfaces;
    DLIList<OCCCurve*> curves;

    body->get_all_surfaces(surfaces);
    body->get_all_curves(curves);
    
    CubitBoolean on_vertex = CUBIT_FALSE;
    CubitBoolean on_curve = CUBIT_FALSE;
    for (int j = 0; j < vector_list.size(); j ++)
    {
      CubitVector* v = vector_list[j];
      if(v == NULL)
        continue;
      for (int k = 0;  k < curves.size(); k ++)
      {
         OCCCurve* curve = curves.get_and_step();
         CubitPointContainment pc = curve->point_containment(*v);
         if(pc == CUBIT_PNT_BOUNDARY)
         {
           on_vertex = CUBIT_TRUE;
           break;
         }

         else if( pc == CUBIT_PNT_INSIDE)
         {
           LocOpe_SplitShape splitor(*from_shape); 
           on_curve = CUBIT_TRUE;
           TopoDS_Edge edge = *curve->get_TopoDS_Edge();
           gp_Pnt pt = gp_Pnt(v->x(), v->y(), v->z());
           TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(pt);
           double param = curve->u_from_position(*v);
           splitor.Add(vertex, param, edge);
           
           //update the curve_list
           TopTools_ListOfShape edge_shapes;
           edge_shapes.Assign(splitor.DescendantShapes(edge));
           while(edge_shapes.Extent())
           {
             TopoDS_Shape edge_shape = edge_shapes.First();
             TopoDS_Edge occ_edge = TopoDS::Edge(edge_shape);
             OCCCurve* test_curve; 
             test_curve = CAST_TO(OCCQueryEngine::instance()->populate_topology_bridge(occ_edge), OCCCurve);
             if(test_curve)
               curves.append(test_curve);
             edge_shapes.RemoveFirst();
           }
           curves.remove(curve);
           
           TopTools_ListOfShape shapes;
           shapes.Assign(splitor.DescendantShapes(*from_shape));
           if(from_shape->TShape()->ShapeType() ==TopAbs_COMPSOLID)
             OCCBody::update_OCC_entity(*from_shape, shapes.First(), 
                    (BRepBuilderAPI_MakeShape*) NULL, &splitor);
           else if(shapes.First().TShape()->ShapeType() == TopAbs_SOLID)
             OCCLump::update_OCC_entity(TopoDS::Solid(*from_shape), 
                    shapes.First(), 
                    (BRepBuilderAPI_MakeShape*) NULL, &splitor);
     
           else if(shapes.First().TShape()->ShapeType() == TopAbs_SHELL)
             OCCShell::update_OCC_entity(TopoDS::Shell(*from_shape),
                    shapes.First(), 
                    (BRepBuilderAPI_MakeShape*) NULL, &splitor);

           else if(shapes.First().TShape()->ShapeType() == TopAbs_FACE)
             OCCSurface::update_OCC_entity(TopoDS::Face(*from_shape), 
                    shapes.First(), 
                    (BRepBuilderAPI_MakeShape*) NULL, NULL, &splitor);

           from_shape->Nullify();
           delete from_shape;
           from_shape = new TopoDS_Shape(shapes.First());
           break;
         }  
       } 
       if(on_vertex || on_curve)
         continue;

       //check possible on surface
       for(int n = 0; n < surfaces.size(); n ++)
       {
          OCCSurface* surface = surfaces.get_and_step();
          if(!surface->is_position_on(*v))
            continue;
           
          CubitPointContainment ps = surface->point_containment(*v);
          if(ps == CUBIT_PNT_INSIDE)
          {
             Point* p = make_Point(*v);
             if(p)
               surface->add_hardpoint(CAST_TO(p, OCCPoint));
             break;
          }
       }
    }
  }       

  shape_to_bodySM(shape_list, new_body_list);
  
  return stat;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: Projects a list of Curves on to a list of Surfaces
//              and imprint the faces with the new Curves
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     
OCCModifyEngine::imprint_projected_edges( DLIList<Surface*> &ref_face_list,
                                          DLIList<Curve*> &ref_edge_list,
                                          DLIList<BodySM*>& new_body_list,
                                          bool keep_old_body,
                                          bool keep_free_edges) const
{
  DLIList<Curve*> projected_curves;
  CubitStatus 
     stat = project_edges(ref_face_list, ref_edge_list, projected_curves);
  if(!stat)
    return stat;

  // imprint Surface with curves
  stat = imprint(ref_face_list, projected_curves, new_body_list, keep_old_body );

  if(keep_free_edges)
     return  stat;

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< projected_curves.size();i++)
  {
    // Now delete this Curve and its underlying solid model entities

    Curve* curve = projected_curves.get_and_step();
    stat = OCCQueryEngine::instance()->delete_solid_model_entities( curve );
    if (stat == CUBIT_FAILURE)
    {
       PRINT_ERROR("In OCCQueryEngine::delete_geometry\n"
                 "       Could not delete OCCCurve.\n"
                 "       The Model database is likely corrupted "
                 "due to\n       this unsuccessful deletion.\n" );
    }
  } 
  return stat;
}
//===============================================================================
// Function   : project_edges
// Member Type: PUBLIC
// Description: Projects a list of Curves on to a list of Surfaces
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::project_edges( DLIList<Surface*> &ref_face_list,
                                            DLIList<Curve*> &ref_edge_list,
                                            DLIList<Curve*> &projected_curves,
                                            bool print_error ) const

{
  CubitVector* v = NULL;
  Curve* projected_curve = NULL;
  DLIList<Point*> points;
  //project curves onto surfaces.
  for(int i = 0; i < ref_edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(ref_edge_list.get_and_step(), OCCCurve);
    if(!curve)
       continue;

    for (int j = 0; j < ref_face_list.size(); j++)
    {
      OCCSurface* surface = CAST_TO(ref_face_list.get_and_step(), OCCSurface); 
      if(!surface)
        continue;
      if(surface->is_closed_in_U() || surface->is_closed_in_V())
      {
        if(print_error)
          PRINT_ERROR("This function can't project curves on closed surfaces.\n");
        return CUBIT_FAILURE;
      }
      
      projected_curve = NULL;
      projected_curve = curve->project_curve(surface, points, CUBIT_FALSE, v);
      if(projected_curve)
        projected_curves.append_unique(projected_curve);
    }
  }
  while(points.size() > 0)
    delete points.pop();
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: Projects a list of curves on to a list of surfaces
//              and imprint the bodies with the new curves
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus 
OCCModifyEngine::imprint_projected_edges(DLIList<Surface*> &ref_face_list,
                                         DLIList<BodySM*> &body_list,
                                         DLIList<Curve*> &ref_edge_list,
                                         DLIList<BodySM*>& new_body_list,
                                         bool keep_old,
                                         bool keep_free_edges) const
{
  DLIList<Curve*> projected_curves;
  CubitStatus
     stat = project_edges(ref_face_list, ref_edge_list, projected_curves);
  if(!stat)
    return stat; 
  return CUBIT_FAILURE;

  // imprint bodies with curves
  stat = imprint(body_list,projected_curves, new_body_list, keep_old);

  if (keep_free_edges)
        return  stat;

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< projected_curves.size();i++)
  {
    // Now delete this Curve 
    Curve* curve = projected_curves.get_and_step();
    stat = OCCQueryEngine::instance()->
          delete_solid_model_entities( curve );
    if (stat == CUBIT_FAILURE)
    {
       PRINT_ERROR("In OCCModifyEngine::delete_geometry\n"
                   "       Could not delete Curve.\n"
                   "       The Model database is likely corrupted "
                   "due to\n       this unsuccessful deletion.\n" );
    }
  }
  return stat; 
}

//===============================================================================
// Function   : intersect
// Member Type: PUBLIC
// Description: intersect boolean operation of body with list of bodies.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::intersect(BodySM*  tool_body_ptr,
                                       DLIList<BodySM*>  &from_bodies,
                                       DLIList<BodySM*>  &new_bodies,
                                       bool  keep_old) const
{
  DLIList<BodySM*> tool_bodies;
  DLIList<TopoDS_Shape*> tool_shapes;
  DLIList<CubitBoolean> is_tool_volume, is_volume;
  
  tool_bodies.append(tool_body_ptr);
  //get tool_body's underlying shape, copy it, so boolean wouldn't touch it.
  CubitStatus stat = 
       get_shape_list(tool_bodies, tool_shapes, is_tool_volume, CUBIT_TRUE); 
  if(!stat)
    return stat;

  DLIList<TopoDS_Shape*> shape_list;
  stat =  get_shape_list(from_bodies, shape_list, is_volume, keep_old);
  if(!stat)
  {
    for (int i = 0; i < tool_shapes.size(); i++)
    {
       TopoDS_Shape* shape = tool_shapes.get_and_step();
       delete shape;
    }
    tool_shapes.clean_out();
    return CUBIT_FAILURE;
  }

  TopoDS_Shape* tool_shape = tool_shapes.get();
  CubitBoolean has_changed;
  DLIList<TopologyBridge*> tbs;
  for (int i = 0; i < shape_list.size(); i++)
  { 
    TopoDS_Shape* from_shape = shape_list.get_and_step();
    BodySM* from_body = from_bodies.get_and_step();
    BRepAlgoAPI_Common intersector(*from_shape, *tool_shape);
    TopoDS_Shape common_shape = intersector.Shape();
    check_operation(common_shape, from_shape, is_volume[i], has_changed, 
                    &intersector, keep_old); 

    if(from_shape->IsNull() )
    {
      PRINT_INFO("The %d body did not have common part with the tool_body.\n", i+1);
    }
    else
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*from_shape);
  }
  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }
  
  if(tbs.size() == 0)
    stat = CUBIT_FAILURE;
    
  //ok, we're done with all cuts, delete unnecessaries.
  if(!keep_old)
    OCCQueryEngine::instance()->delete_solid_model_entities(tool_body_ptr);   

  for(int i = 0; i < tool_shapes.size(); i++)
  {
    TopoDS_Shape* shape = tool_shapes.get_and_step();
    shape->Nullify();
    delete shape;
  }

  if(keep_old)
  {
    int size  = shape_list.size();
    for (int i = 0; i < size; i++)
    {
      TopoDS_Shape* shape = shape_list.pop();
      shape->Nullify();
      delete shape;
    }
  }
  if(!stat)
    return stat;
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : check_operation
// Member Type: PRIVATE
// Description: check and update the from_shape according to type of the body.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
void OCCModifyEngine::check_operation(TopoDS_Shape& cut_shape,
                                      TopoDS_Shape*& from_shape, //output
                                      CubitBoolean  is_volume,
                                      CubitBoolean& has_changed, //output
                                      BRepAlgoAPI_BooleanOperation* op,
                                      CubitBoolean keep_old) const
{
   //compare to see if the from_shape has gotten cut.
   if(is_volume)
   {
     GProp_GProps myProps;
     BRepGProp::VolumeProperties(*from_shape, myProps);
     double orig_mass = myProps.Mass();
     TopTools_IndexedMapOfShape M;
     TopExp::MapShapes(cut_shape, TopAbs_SOLID, M);
     double after_mass = 0.0;
     CubitBoolean no_volume = CUBIT_FALSE;
     if(M.Extent() > 0)
     {
       BRepGProp::VolumeProperties(cut_shape, myProps);
       after_mass = myProps.Mass();
     }
     else
       no_volume = CUBIT_TRUE;
 
     if(fabs(-after_mass + orig_mass) <= TOL)
     {
        has_changed= CUBIT_FALSE; //common is itself
        return;
     }

     //got cut. Update the entities
     if(after_mass < TOL || no_volume) //no common section
       cut_shape.Nullify();
     has_changed = CUBIT_TRUE;
     TopExp_Explorer Ex;
     Ex.Init(*from_shape, TopAbs_SOLID);
     TopoDS_Solid old_solid = TopoDS::Solid(Ex.Current()); 
     OCCLump::update_OCC_entity(old_solid , cut_shape, op);
   }
   else
   {
     GProp_GProps myProps;
     BRepGProp::SurfaceProperties(*from_shape, myProps);
     double orig_mass = myProps.Mass();
     BRepGProp::SurfaceProperties(cut_shape, myProps);
     double after_mass = myProps.Mass();
     if(fabs(-after_mass + orig_mass) <= TOL)
     {
       has_changed= CUBIT_FALSE; //common is itself, or not cut
       return;
     }
     //got cut. Update the entities
     if(after_mass < TOL)//no common section
       cut_shape.Nullify();
     has_changed = CUBIT_TRUE;
     if(from_shape->TShape()->ShapeType() == TopAbs_SHELL)
     {
       TopoDS_Shell old_shell = TopoDS::Shell(*from_shape);
       OCCShell::update_OCC_entity(old_shell,cut_shape, op);
     }
     else
     {
       TopoDS_Face old_face = TopoDS::Face(*from_shape);
       OCCSurface::update_OCC_entity(old_face,cut_shape, op);
     }
  }
  if(keep_old)
    delete from_shape;
  from_shape = new TopoDS_Shape(cut_shape);
}

//===============================================================================
// Function   : chop
// Member Type: PUBLIC
// Description: chop boolean operation between OCC-based bodies
//              bodies has a size() = 2, a blank body and a tool body.
//              chops the blank with the  tool, returing the body formed
//              by subtracting the tool from the blank, and the body formed
//              by intersecting the tool with the blank, simultaneously.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus  OCCModifyEngine::chop(DLIList<BodySM*>& bodies, 
                                   DLIList<BodySM*> &intersectBodies, 
                                   DLIList<BodySM*> &outsideBodies,
                                   BodySM*& leftoversBody,
                                   bool keep_old ,
                                   bool nonreg) const
{
  //according to Acis chop function, leftoverBody = 0;
  leftoversBody = 0;

  //there's no effect of nonreg. keep_old mean if to keep the tool_body
  if(bodies.size() != 2)
  {
    PRINT_WARNING("Chop operation works only on two bodies. Nothing modified\n");  
    return CUBIT_FAILURE; 
  }
  
  //outsideBodies keeps the surface, curve ids if keep_old is false.
  BodySM* blank_body = bodies.get();
  
  //copy blank_body for intersect operation, because it will get changed.
  DLIList<BodySM*> tool_bodies, from_bodies;
  from_bodies.append(blank_body);
  BodySM* tool_body = bodies.step_and_get();
  tool_bodies.append(tool_body);
  
  CubitStatus stat = intersect(tool_body, from_bodies, 
                               intersectBodies, CUBIT_TRUE);

  if(!stat)
    return CUBIT_FAILURE;

  stat = subtract(tool_bodies, from_bodies, outsideBodies, 
                  CUBIT_FALSE, keep_old);
  
  return stat;
}

//===============================================================================
// Function   : unite
// Member Type: PUBLIC
// Description: unite boolean operation between OCC-based bodies
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     OCCModifyEngine::unite(DLIList<BodySM*> &bodies, 
                                       DLIList<BodySM*> &newBodies,
                                       bool keep_old) const
{
  if(bodies.size() < 2)
  {
    newBodies = bodies;
    return CUBIT_SUCCESS;
  }
  //all bodies must have only one lump in order for boolean operation to work.
  DLIList<Lump*> lumps;
  for (int i = 0; i < bodies.size(); i++)
  {
    lumps = CAST_TO(bodies.get_and_step(), OCCBody)->lumps();
    if (lumps.size() > 1)
    {
      PRINT_WARNING("All bodies must have only one lump for boolean operations to work.\n");
      return CUBIT_FAILURE;
    }
  }
 
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_volume;
  CubitStatus stat =
    get_shape_list(bodies, shape_list, is_volume, keep_old);

  if(!stat)
    return stat;

  //find a non-sheet body to be the first shape
  TopoDS_Shape* first_shape;
  TopoDS_Shape* second_shape;
  CubitBoolean first_is_volume;
  int index = 0;
  if((first_is_volume = is_volume.move_to(CUBIT_TRUE)))
    index = is_volume.get_index();

  first_shape = shape_list[index];
  shape_list.remove(first_shape);
  is_volume.step(index);
  is_volume.remove();
  bodies.step(index);
  BodySM* removed_body = bodies.remove();

  int size = shape_list.size();
  for(int i = 0; i < size; i++)
  {
    second_shape = shape_list.get_and_step();

    BRepAlgoAPI_Fuse fuser(*first_shape, *second_shape);
    TopoDS_Shape new_shape = fuser.Shape();

    CubitBoolean has_changed;
    check_operation(new_shape, first_shape, first_is_volume, has_changed, &fuser, keep_old);
 
    check_operation(new_shape,second_shape, is_volume[i], has_changed, &fuser, keep_old);
  }      

  //ok, we're done with all unites, construct new Body'
  DLIList<TopologyBridge*> tbs;
  tbs += OCCQueryEngine::instance()->populate_topology_bridge(*first_shape);

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      newBodies.append(bodysm);
  }

  //ok, we're done with all unites, delete unnecessaries.
  if(keep_old)
  {
    shape_list.append(first_shape);
    int size  = shape_list.size();
    for (int i = 0; i < size; i++)
    {
      TopoDS_Shape* shape = shape_list.pop();
      shape->Nullify();
      delete shape;
    }
  }

  bodies.append(removed_body); 
  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::thicken( DLIList<BodySM*>& bodies,
                                      DLIList<BodySM*>& new_bodies,
                                      double depth,
                                      CubitBoolean both) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : hollow
// Member Type: PUBLIC
// Description: Hollow existing solid body by remove one or several surfaces 
//              Can only take one body at a time.
//              depth > 0, thick body going outside bodies
//              depth < 0, thick body going inside bodies
// Author     : Jane Hu 
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::hollow( DLIList<BodySM*>& bodies, 
                                     DLIList<Surface*>& surfs_to_remove,
                                     DLIList<BodySM*>& new_bodies,
                                     double depth) const
{
  if(bodies.size() != 1 || surfs_to_remove.size() < 1)
  {
    PRINT_ERROR("Making thick solid in OCC will take one body and at least one surface at a time.\n"); 
    return CUBIT_FAILURE;
  }

  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_volume;
  CubitStatus stat = get_shape_list(bodies, shape_list, is_volume, CUBIT_FALSE);

  if(!stat)
    return stat;

  if(!is_volume.get())//sheet body
  {
    PRINT_ERROR("Making thick solid in OCC needs an initial solid body to hollow with.\n");
    return CUBIT_FAILURE;
  }

  //make sure the body to be hollowed has only one lump
  OCCBody* occ_body = CAST_TO(bodies.get(), OCCBody);
  DLIList<Lump*> lumps;
  lumps = occ_body->lumps();
  if(lumps.size()!=1)
  {
    PRINT_ERROR("bodies with more than one lump can't be hollowed to make a thick body.\n");
    return CUBIT_FAILURE;
  }

  //make sure surfs_to_remove are all in bodies
  DLIList<OCCSurface*> surfaces;
  TopTools_ListOfShape face_shapes;
  occ_body->get_all_surfaces(surfaces);
  for(int i = 0; i < surfs_to_remove.size(); i++)
  {
    OCCSurface* occ_surf = CAST_TO(surfs_to_remove.get(), OCCSurface);
    if(!occ_surf)
      continue;
    if(!surfaces.is_in_list(occ_surf))
      continue;
    TopoDS_Face * face = occ_surf->get_TopoDS_Face();
    face_shapes.Append(*face); 
  }

  if(face_shapes.IsEmpty())
  {
    PRINT_ERROR("The surfaces provided should be from the body to be hollowed.\n");
    return CUBIT_FAILURE;
  }
  
  double dTOL = 1.e-3; //hard coded for now, can be changed by application
  TopoDS_Shape* solid = shape_list.get();
  BRepOffsetAPI_MakeThickSolid hollower(*solid, face_shapes, depth, dTOL,
                                        BRepOffset_Skin, Standard_False,
                                        Standard_False, GeomAbs_Intersection);
  TopoDS_Shape new_shape = hollower.Shape();
  TopoDS_Solid old_solid = TopoDS::Solid(*solid);
  OCCLump::update_OCC_entity(old_solid , new_shape, &hollower); 
 
  //ok, we're done with all hollowing, construct new Body'
  DLIList<TopologyBridge*> tbs;
  tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }

  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : flip_normals
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine :: flip_normals( DLIList<Surface*>& face_list ) const
{
  DLIList<Surface*> surface_list;
  while (face_list.size())
  {
    OCCSurface* occ_surface = CAST_TO(face_list.pop(), OCCSurface);
    OCCShell* occ_shell = occ_surface->my_shell();
    DLIList<OCCSurface*> surfaces;
    surfaces.append(occ_surface);
    if(occ_shell) //find all surfaces in face_list that belong to this shell
    {
      int size = face_list.size();
      for ( int i = 0; i < size; i++)
      {
        occ_surface = CAST_TO(face_list.get(), OCCSurface); 
        if(occ_shell == occ_surface->my_shell())
          surfaces.append(CAST_TO(face_list.remove(),OCCSurface));
        else
          face_list.step();
      } 
      
      if (!occ_shell->is_sheet())
      {
        DLIList<OCCCoFace*> cofaces;
        cofaces = occ_shell->cofaces();
        for(int i = 0; i < cofaces.size(); i++)
        {
          OCCCoFace* coface =cofaces.get_and_step();
          occ_surface = coface->surface();
          if(surfaces.is_in_list(occ_surface))
          { 
            TopoDS_Face* topoface = occ_surface->get_TopoDS_Face();
            TopAbs_Orientation ori = topoface->Orientation();
            topoface->Orientation(ori == TopAbs_FORWARD ? TopAbs_REVERSED :
                                                      TopAbs_FORWARD);
            occ_surface->set_TopoDS_Face(*topoface);
            coface->set_sense(coface->sense() == CUBIT_FORWARD ? 
                                 CUBIT_REVERSED : CUBIT_FORWARD);
            surface_list.append(occ_surface);
          }
        }
      }
      else //sheet body 
      {
        TopoDS_Face* topoface = occ_surface->get_TopoDS_Face();
        TopAbs_Orientation ori = topoface->Orientation();
        topoface->Orientation(ori == TopAbs_FORWARD ? TopAbs_REVERSED :
                                                    TopAbs_FORWARD);
        occ_surface->set_TopoDS_Face(*topoface);
        surface_list.append(occ_surface);
      }
      PRINT_INFO( "Modified volume\n" );
    }
    else
      PRINT_WARNING( "Volume was not modified\n" );
  }
  face_list = surface_list;
  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : sweep_translational
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 09/08
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_translational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& sweep_vector,
  double draft_angle, //in Radius
  int draft_type, //RightCorner=1 or RoundCorner =2
  bool switchside,//not used, shell and surfaces are one sided, not like Acis
  bool rigid, //not used here, in Acis, it means whether the end surface is
              // parallel to the starting surface, or perpendicular to the path
  Surface* stop_surf,
  BodySM* to_body) const
{
  //in OCC, there's no sweep surface with draft option, this can be done by
  //creating draft shell then make solid to achieve.
  //while if draft_angle is 0, and stop_surf = 0, to_body = 0
  // directly use sweep functions.
  TopoDS_Shape *stop_shape = NULL;
  if(stop_surf)
  {
     OCCSurface* occ_surface = CAST_TO(stop_surf, OCCSurface);
     stop_shape = occ_surface->get_TopoDS_Face();
  }
  else if(to_body)
  {
     OCCBody* occ_body = CAST_TO(to_body, OCCBody);
     if(occ_body->is_sheet_body())
       stop_shape = occ_body->my_sheet_surface()->get_TopoDS_Face();
     else if(occ_body->shell())
       stop_shape = occ_body->shell()->get_TopoDS_Shell();
     else
       stop_shape = occ_body->get_TopoDS_Shape();
  }

  gp_Dir adir(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());
  gp_Vec aVec(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());
 
  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, &sweep_vector, toposhape);   
      if(!stat)
        continue;
    }
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
    }
  
    DLIList<TopologyBridge*> tbs;
    //create the draft or the sweep
    if(stop_shape == NULL && draft_angle == 0.)
    {
      BRepSweep_Prism swept(toposhape, aVec);
      TopoDS_Shape new_shape = swept.Shape();
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
      assert(tbs.size() == 1);

      BodySM* bodysm = CAST_TO(tbs.get(), BodySM); 
      if (bodysm)
        result_body_list.append(bodysm);
      continue;
    }

    BRepOffsetAPI_MakeDraft draft(toposhape, adir, draft_angle);
    BRepBuilderAPI_TransitionMode Cornertype;
    if(draft_type == 1)
      Cornertype = BRepBuilderAPI_RightCorner;
    else if(draft_type == 2)
      Cornertype = BRepBuilderAPI_RoundCorner;

    draft.SetOptions(Cornertype);
    if(stop_shape)
      draft.Perform(*stop_shape);    
    else
      draft.Perform(sweep_vector.length());
    TopoDS_Shape new_shape = draft.Shape();

    tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);

    assert(tbs.size() == 1);

    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
    if(bodysm && surface != NULL) //only gets swept side and original surfaces
    {
       //get surfaces from the shell body and create a top surface to
       //make a swept solid.
       OCCShell* occ_shell = CAST_TO(bodysm, OCCBody)->shell();
       if(!occ_shell)
       {
         PRINT_WARNING("Sweep surface failed inside OCC engine.\n");
         return CUBIT_FAILURE;
       }
       DLIList<OCCCoFace*> cofaces = occ_shell->cofaces();
       DLIList<Surface*> surface_list;
       for(int i = 0; i < cofaces.size(); i++)
         surface_list.append(cofaces.get_and_step()->surface());

       //create the top surface from edges.
       DLIList<OCCCoEdge*> coedges;
       for(int i = 0; i < surface_list.size(); i++)
         CAST_TO(surface_list.get_and_step(), OCCSurface)->get_coedges(coedges);
       for(int i = 0; i < coedges.size(); i++)
       {
         OCCCoEdge* coedge = coedges[i];
         if(coedge == NULL)
           continue;
         for(int j = i+1; j < coedges.size(); j++)
         {
            OCCCoEdge* temp_coedge = coedges[j];
            if(temp_coedge == NULL)
              continue; 
            if(coedge->curve() == temp_coedge->curve() &&
               coedge->sense() != temp_coedge->sense())
            {
              coedges.move_to(coedge);
              coedges.change_to((OCCCoEdge*)NULL);
              coedges.move_to(temp_coedge);
              coedges.change_to((OCCCoEdge*)NULL);
            }
         }
       } 
       coedges.remove_all_with_value(NULL);
       assert(coedges.size() > 0);
       DLIList<Curve*> curves;
       for(int i = 0; i < coedges.size(); i++)
         curves.append(coedges.get_and_step()->curve());

       Surface* surf = make_Surface(PLANE_SURFACE_TYPE, curves);
       if(!surf)
         surf = make_Surface(BEST_FIT_SURFACE_TYPE, curves);
       if(!surf)
       {
         PRINT_ERROR("Can't calculate for the top surface.\n");
         continue;
       }
       surface_list.append(surf);
       DLIList<BodySM*> bodies;
       create_solid_bodies_from_surfs(surface_list, bodies);

       if(bodies.size() == 1)
         bodysm = bodies.get();
       else
       {
         PRINT_WARNING("Sweep surface failed in creating solid.\n");
         return CUBIT_FAILURE;
       }
    }
    if (bodysm)
      result_body_list.append(bodysm);
  }
  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::get_sweepable_toposhape(OCCCurve*& curve,
                                              TopoDS_Shape& toposhape)const
{
  DLIList<OCCLoop*> loops;
  loops =  curve->loops();
  if( loops.size()) //not a free curve
  {
    //copy the curve
    Curve* c_curve = make_Curve(curve);
    if(c_curve)
     curve = CAST_TO(c_curve, OCCCurve);
    else
    {
      PRINT_ERROR("Can't copy the curve for sweep.\n");
      return CUBIT_FAILURE;
    }
  }
  TopoDS_Edge *edge = curve->get_TopoDS_Edge( );
  toposhape = BRepBuilderAPI_MakeWire(*edge);
  return CUBIT_SUCCESS;
}

CubitStatus OCCModifyEngine::get_sweepable_toposhape(OCCSurface*& surface,
                                                  const CubitVector* sweep_v_p,
                                                  TopoDS_Shape& toposhape)const
{
  GeometryEntity* ref_ent = NULL;
  //Make copy of the surface if it's not a sheet surface.
  Surface* c_surface = NULL;
  if(surface != NULL)
  {
    //check if the surface is sheet body, if not, copy it.
    if(surface->my_body() == NULL) //not a sheet body
    {
      c_surface = make_Surface(surface);
      if (c_surface == NULL)
      {
         PRINT_ERROR("Cannot copy surface in sweep_translational.\n");
         return CUBIT_FAILURE;
      }
      surface = CAST_TO(c_surface, OCCSurface);
    }
 
    if(sweep_v_p)
    {
      CubitVector center = surface->center_point();
      CubitVector normal;
      surface->closest_point(center,NULL,&normal);
      CubitVector sweep_vector = *sweep_v_p;
      if(normal % sweep_vector > 0)
      {
        DLIList<Surface*> surfaces;
        surfaces.append(surface);
        flip_normals(surfaces);
        surface = CAST_TO(surfaces.get(), OCCSurface);
        ref_ent = (GeometryEntity *)surface;
      }

      else if(normal % sweep_vector == 0)
      {
        PRINT_ERROR("Sweeping direction should not be on the surface.\n");
        return CUBIT_FAILURE;
      }
      else 
        ref_ent = (GeometryEntity *)surface;
    }
    else
      ref_ent = (GeometryEntity *)surface;

    if(surface->my_body() != NULL) //sheet body
    {
      delete surface->my_body();
      delete surface->my_shell();
      delete surface->my_lump();
      surface->set_shell(NULL);
      surface->set_lump(NULL);
      surface->set_body(NULL);
    }

    TopoDS_Shape* toposhape_prt = 
          OCCQueryEngine::instance()->get_TopoDS_Shape_of_entity(ref_ent);

    if(!toposhape_prt)
    {
      PRINT_WARNING("GeometryEntity without TopoDS_Shape found.\n");
      return CUBIT_FAILURE;
    }
    toposhape = *toposhape_prt;
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_perpendicular(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  double distance,
  double draft_angle,
  int draft_type,
  bool switchside, //has no effect
  bool rigid, //has no effect
  Surface* stop_surf,
  BodySM* to_body) const
{
  //find the vector perpendicular to the ref_ent normal, and sweep_translate
  //the 'distance' along this vector
  DLIList<GeometryEntity*> edge_list;
  CubitVector vec;
  for(int i = 0; i < ref_ent_list.size(); i++)
  {
     GeometryEntity *ref_ent = ref_ent_list.get_and_step();
     Surface *face = CAST_TO(ref_ent, Surface);
     Curve* edge = CAST_TO(ref_ent, Curve);
     DLIList<GeometryEntity*> face_list;
     if(face != NULL)
     {
        OCCSurface* occ_face = CAST_TO(face, OCCSurface);
        CubitVector center = occ_face->center_point();
        CubitVector closest_p, unit_normal;
        CubitStatus stat = 
                    occ_face->closest_point(center, &closest_p, &unit_normal);
        if(stat)
        {
          vec = distance * unit_normal;
          face_list.append(ref_ent);
          stat = sweep_translational(face_list, result_body_list, vec, 
                                     draft_angle, draft_type, switchside,
                                     rigid, stop_surf, to_body);
       }
     }
     else if (edge != NULL)
     {
        edge_list.append(ref_ent);
     }
  }
  if(edge_list.size())
    PRINT_ERROR("Curves cannot be swept perpendicularly, please use the vector sweep.\n");

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : sweep_rotational
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu 
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_rotational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& point,
  const CubitVector& direction,
  double angle, //in radians
  int steps,  //not used
  double draft_angle, //not used
  int draft_type,  //not used
  bool switchside, //not used
  bool make_solid,
  bool rigid,  //not used
  Surface* stop_surf,  //not used
  BodySM* to_body ) const  //not used
{
  gp_Dir adir(direction.x(), direction.y(), direction.z()); 
  gp_Pnt pt = gp_Pnt( point.x(), point.y(), point.z());
  gp_Ax1 axis = gp_Ax1(pt, adir);

  gp_Lin line = gp_Lin(axis);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(line);
  TopoDS_Edge* p_edge = new TopoDS_Edge(edge);
  OCCCurve* acurve = new OCCCurve(p_edge);
  assert(acurve);

  CubitVector start;
  CubitVector end;

  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface or curve for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, (CubitVector*)NULL, toposhape);
      if(!stat)
        continue;
      //only non-intersecting of surface and axis can be swept.
      DLIList<CubitVector*> intersect_pts;
      OCCQueryEngine::instance()->get_intersections(acurve, surface,
                                    intersect_pts, CUBIT_TRUE);
      if(intersect_pts.size() > 0)
      { 
        PRINT_ERROR("Only surfaces with no intersection point with the axis can be revolve-swept.\n");
        continue;
      } 
    }
    else if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
      //closed curve can't intersect with the axis, while open curve can only
      //intersect the axis at the end points. 
      //only curve not intersecting with axis in curve's middle locations
      //can be revolved
      DLIList<CubitVector*> intersect_pts;
      OCCQueryEngine::instance()->get_intersections(curve, acurve,
                                  intersect_pts, CUBIT_TRUE, CUBIT_TRUE);
      if(!toposhape.Closed())
      {
        //get start and end points
        DLIList<OCCPoint*> point_list;
        curve->get_points(point_list);
        assert(2 == point_list.size());
        GeometryType type = curve->geometry_type();
        start = point_list.get_and_step()->coordinates();
        end = point_list.get()->coordinates();
        CubitBoolean start_int = CUBIT_FALSE;
        CubitBoolean end_int = CUBIT_FALSE;
        if(intersect_pts.size() > 0)
        {
          CubitBoolean non_int = CUBIT_FALSE;
          for(int i = 0; i < intersect_pts.size(); i++)
          {
             CubitVector prt = *(intersect_pts.get_and_step());
             if(prt.distance_between(start) > TOL &&
                prt.distance_between(end) > TOL)
             {
                non_int = CUBIT_TRUE;
                PRINT_ERROR("Only curves with no intersection point with the axis can be revolve-swept.\n");
                break;
             }
             else if(prt.distance_between(start) <= TOL)
                start_int = CUBIT_TRUE;
             else if(prt.distance_between(end) <= TOL)
                end_int = CUBIT_TRUE;
          }
          if(non_int)
            continue;
          if(start_int && end_int && type == STRAIGHT_CURVE_TYPE)
          {
            PRINT_ERROR("Sweep along curve itself is not allowed.\n");
            continue;
          } 
        }
      }
      else
      {
        if(intersect_pts.size() > 0)
        {
          PRINT_ERROR("Only curves with no intersection point with the axis can be revolve-swept.\n");
          continue;
        }  
      }
    } 
    else
    {
      PRINT_ERROR("Only surface or curve can be revolve-swept.\n");
      continue;
    }
    TopoDS_Shape new_shape;
    DLIList<TopologyBridge*> tbs;
    if(make_solid && curve != NULL )
    //giving an open wire and want a solid
    {
      if(!toposhape.Closed())
      {
        //project the start and end points onto the axis
        CubitBoolean start_closed = CUBIT_FALSE;
        CubitBoolean end_closed = CUBIT_FALSE;
        if(acurve->point_containment(start) != CUBIT_PNT_OFF)
          start_closed = CUBIT_TRUE;
        if(acurve->point_containment(end) != CUBIT_PNT_OFF)
          end_closed = CUBIT_TRUE; 
        CubitVector start_proj, end_proj;
        TopoDS_Edge edge1, edge2;
        BRepBuilderAPI_MakeWire m_wire;
        if(!start_closed)
        {
          acurve->closest_point(start, start_proj);
          gp_Pnt pt1 = gp_Pnt( start.x(), start.y(), start.z()); 
          gp_Pnt pt2 = gp_Pnt( start_proj.x(), start_proj.y(), start_proj.z());
          edge1 = BRepBuilderAPI_MakeEdge(pt1, pt2);
          m_wire.Add(edge1);
          m_wire.Add(TopoDS::Wire(toposhape));
        }
        else
        {
          m_wire.Add(TopoDS::Wire(toposhape));
          start_proj = start;
        }
 
        if(!end_closed)
        {
          acurve->closest_point(end,end_proj);
          gp_Pnt pt1 = gp_Pnt( end.x(), end.y(), end.z());
          gp_Pnt pt2 = gp_Pnt( end_proj.x(), end_proj.y(), end_proj.z());
          edge2 = BRepBuilderAPI_MakeEdge(pt1, pt2);
          m_wire.Add(edge2);
        }
      
        else
          end_proj = end;
        
        gp_Pnt pt1 = gp_Pnt( end_proj.x(), end_proj.y(), end_proj.z());
        gp_Pnt pt2 = gp_Pnt( start_proj.x(), start_proj.y(), start_proj.z());
        TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(pt2, pt1);
        m_wire.Add(edge3);
      
        TopoDS_Wire wire = m_wire.Wire();
        toposhape = BRepBuilderAPI_MakeFace(wire);
      }
      else //closed
      {
        TopoDS_Wire wire = TopoDS::Wire(toposhape);
        toposhape = BRepBuilderAPI_MakeFace(wire);
      }
    }
    BRepSweep_Revol revol(toposhape, axis, angle);
    new_shape = revol.Shape();

    tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
    assert(tbs.size() == 1);

    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
    if (bodysm)
      result_body_list.append(bodysm);
    continue;
  }
  OCCQueryEngine::instance()->delete_solid_model_entities(acurve);
  if(result_body_list.size()>0)
    return CUBIT_SUCCESS;
  else 
    return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_along_curve
// Member Type: PUBLIC
// Description: The ref_edge_list must provide a list of curves which are
//              connected, and making G1 continuous wire.
// Author     : Jane Hu
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine::sweep_along_curve(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  DLIList<Curve*>& ref_edge_list,
  double draft_angle, //only used for straight curve case
  int draft_type, //only used for straight curve case
  bool rigid, //not used
  Surface* stop_surf, //not used
  BodySM* to_body) const //not used
{
  //make wire out of ref_edge_list
  BRepBuilderAPI_MakeWire awire;
  TopTools_ListOfShape L;
  OCCCurve* occ_curve = NULL;
  GeometryType type;
  int num_curve = 0;
  for(int i = 0; i < ref_edge_list.size(); i++)
  {
    Curve* curve = ref_edge_list.get_and_step();
    occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    TopoDS_Edge* topoedge = occ_curve->get_TopoDS_Edge( );
    BRepBuilderAPI_Copy api_copy(*topoedge);
    TopoDS_Shape newShape = api_copy.ModifiedShape(*topoedge);
    L.Append(newShape);
    type = occ_curve->geometry_type();
    num_curve++;
  }
  if(L.IsEmpty())
  {
    PRINT_ERROR("There's no valid sweeping path.\n");
    return CUBIT_FAILURE;
  }
  
  if(num_curve == 1 && type == STRAIGHT_CURVE_TYPE && draft_angle != 0.0)
  {
    DLIList<OCCPoint*> point_list;
    occ_curve->get_points(point_list);
    CubitVector v1 = point_list.get_and_step()->coordinates();
    CubitVector v2 = point_list.get()->coordinates();
    CubitVector sweep_vector = v2-v1;
    return sweep_translational(ref_ent_list,result_body_list,sweep_vector,
                               draft_angle, draft_type, CUBIT_FALSE, 
                               rigid, stop_surf, to_body); 
  }
  awire.Add(L);
  TopoDS_Wire wire;
  wire = awire.Wire();

  BRepTools_WireExplorer it(wire);
  int num_edges = 0;
  for(; it.More(); it.Next())
    num_edges++; 
  
  BRepLib_FuseEdges fuser(wire);
  fuser.SetConcatBSpl();
  fuser.Perform();
  TopoDS_Shape  spline = fuser.Shape();
  wire = TopoDS::Wire(spline);

  DLIList<TopologyBridge*> tbs;
  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface or curve for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, (CubitVector*)NULL, toposhape);
      if(!stat)
        continue;
    } 
    else if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
    }

    //sweep along the wire
    BRepOffsetAPI_MakePipe maker(wire, toposhape);
    if(!maker.IsDone())
    {
      PRINT_ERROR("Can't sweep along the provided curve(s).\n");
      continue;
    }
    TopoDS_Shape newShape = maker.Shape();
    
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(newShape);
    assert(tbs.size() == 1);

    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
    if (bodysm)
      result_body_list.append(bodysm);
    continue;
  }
  return CUBIT_SUCCESS;
}

//HEADER- Webcut-related functions

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                              const CubitVector &v1,
                              const CubitVector &v2,
                              const CubitVector &v3,
                              DLIList<BodySM*>& results_list,
                              bool imprint ) 
{
  CubitStatus stat;
  DLIList<BodySM*> new_BodySMs;
  stat = section(webcut_body_list, v1, v2, v3, new_BodySMs, true, true,false);
  if(stat == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't webcut the bodies using a plane determined by 3 points.\n");
    return stat;
  }
  
  stat = section(webcut_body_list, v1, v2, v3, new_BodySMs, false, false, false);
  if(stat == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't webcut the bodies using a plane determined by 3 points.\n");
    return stat;
  }

  if(imprint)
  {
    BodySM* new_body1, *new_body2;
    for(int i = 0; i < new_BodySMs.size()-1; i ++)
    {
      BodySM* body1 = new_BodySMs[i];
      for(int j = i+1; j < new_BodySMs.size(); j++)
      {
        BodySM* body2 = new_BodySMs[j];
        stat =  this->imprint( body1, body2, new_body1, new_body2, false);
        if(new_body1 && body1 != new_body1)
          new_BodySMs[i] = new_body1;
        if(new_body2 && body2 != new_body2)
          new_BodySMs[j] = new_body2; 
      }
    }
  }

  results_list = new_BodySMs;
  return CUBIT_SUCCESS;  
}

//===============================================================================
// Function   : webcuts a list of bodies using another Body as the tool.
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 01/09
//===============================================================================
CubitStatus    OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                                BodySM const* tool_body,
                                DLIList<BodySM*>& results_list,
                                bool imprint ) 
{
  //do intersect and subtract separately and with imprint option and keep_old
  // is true.
  //tool_body is a const pointer points to varible BodySM object
  //here trying to create a non-const pointer points to the same BodySM object.

  BodySM *body;
  *body = *tool_body;

  CubitStatus stat = intersect(body, webcut_body_list, results_list,
                               CUBIT_TRUE);
 
  if(stat)
  { 
    PRINT_ERROR("Failed to webcut the bodies.\n"); 
    return CUBIT_FAILURE;
  }

  DLIList<BodySM*> tool_bodies;
  tool_bodies.append(body);
  
  stat = subtract(tool_bodies, webcut_body_list, results_list, imprint, 
                  CUBIT_TRUE);

  return stat;
}

//===============================================================================
// Function   : webcut_across_translate
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::webcut_across_translate( DLIList<BodySM*>& /*body_list*/,
                                                          Surface* /*plane_surf1*/,
                                                          Surface* /*plane_surf2*/,
                                                          DLIList<BodySM*>& /*results_list*/,
                                                          bool /*imprint*/ ) const
{
  //from Acis : // Currently no command line to this
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : section
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 11/08
//===============================================================================
CubitStatus OCCModifyEngine::section( DLIList<BodySM*> &section_body_list,
                                      const CubitVector &point_1,
                                      const CubitVector &point_2,
                                      const CubitVector &point_3,
                                      DLIList<BodySM*>& new_body_list,
                                      bool keep_normal_side,
                                      bool keep_old,
                                      bool keep_both_sides)
{
  if (keep_both_sides == CUBIT_TRUE )
  {
     PRINT_ERROR("Please use webcut to for keep both sides option.\n");
     return CUBIT_FAILURE;
  }
 
  //Calculate normal of the section plan
  CubitVector v1, v2, normal;
  v1 = point_2 - point_1;
  v2 = point_3 - point_1; 
  normal = ~(v1 * v2); 
  if(fabs(normal.length() - 1) > TOL)
  {
     PRINT_ERROR("The three points are co-linear, and can't be used as a cutting plane.\n");
     return CUBIT_FAILURE;
  }
  
  if(keep_normal_side)
    normal *= -1;

  gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
  gp_Dir normal_dir(normal.x(), normal.y(), normal.z()); 
  gp_Pln plane(pt, normal_dir);
  gp_Vec vec(normal_dir);
  pt =  pt.Translated(vec);

  TopoDS_Face face = BRepBuilderAPI_MakeFace(plane);
  TopoDS_Solid solid = BRepPrimAPI_MakeHalfSpace(face, pt);
   
  DLIList<CubitBoolean> is_tool_volume;
  is_tool_volume.append(CUBIT_TRUE);
  DLIList<CubitBox*> tool_boxes ;
  Bnd_Box box;
  BRepBndLib::Add(solid, box);
  double min[3], max[3];
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);
  CubitBox* cBox = new CubitBox(min, max);
  
  tool_boxes.append(cBox);
  DLIList<TopoDS_Shape*> solids;
  solids.append(&solid);
  CubitStatus stat = do_subtract(section_body_list, solids, is_tool_volume,
                     &tool_boxes, new_body_list, keep_old) ;
  delete cBox;
  return stat;
}

//===============================================================================
// Function   : split_body
// Member Type: PUBLIC
// Description: Splits multiple lumps in one body into separate bodies
// Author     : Jane Hu 
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::split_body( BodySM *body_ptr,
                                         DLIList<BodySM*> &new_bodies )
{
  OCCBody* occ_body = CAST_TO(body_ptr, OCCBody);
  if(!occ_body)
  {
     PRINT_ERROR("This is not an OCC body to be split.\n");
     return CUBIT_FAILURE;
  }
  DLIList<Lump*> lumps = occ_body->lumps();
  if(lumps.size() == 1) 
  {
    new_bodies.append(body_ptr);
    return CUBIT_SUCCESS;
  }
  for(int i = 0; i < lumps.size(); i++)
  {
    Lump* lump = lumps.get_and_step();
    OCCLump* occ_lump = CAST_TO(lump, OCCLump);
    OCCSurface* occ_surface = occ_lump->my_sheet_surface();
    if(occ_surface) 
    {
      TopoDS_Face* face = occ_surface->get_TopoDS_Face();
      Surface* surface = 
        OCCQueryEngine::instance()->populate_topology_bridge(*face,
                                                             CUBIT_TRUE);
      new_bodies.append(CAST_TO(surface, OCCSurface)->my_body()); 
      continue;
    }
    OCCShell* occ_shell = occ_lump->my_shell();
    if(occ_shell) 
    {
      TopoDS_Shell* shell = occ_shell->get_TopoDS_Shell();
      OCCShell* ashell = 
         OCCQueryEngine::instance()->populate_topology_bridge(*shell,
                                                              CUBIT_TRUE);
      new_bodies.append(ashell->my_body());
      continue;
    }
    else
    {
      TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();
      Lump* alump = 
        OCCQueryEngine::instance()->populate_topology_bridge(*solid,
                                                             CUBIT_TRUE);
      new_bodies.append(CAST_TO(alump,OCCLump)->get_body());
      continue;
    }
  }
  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : reverse_body
// Member Type: PUBLIC
// Description: Turn body inside-out
// Author     : Jane Hu
// Date       : 03/03/09
//===============================================================================
CubitStatus OCCModifyEngine::reverse_body( BodySM* body_ptr )
{
  OCCBody* occ_body = CAST_TO(body_ptr, OCCBody);
  if (!occ_body)
  {
     PRINT_ERROR("Cannot reverse a non-OCC bodySM .\n"
                 "Possible incompatible geometry engines.\n");
     return CUBIT_FAILURE;
  }

  TopoDS_Shape* orig_S;
  orig_S = occ_body->get_TopoDS_Shape();
  TopoDS_Shape S;
  BRep_Builder B;
  if (orig_S->IsNull() ) //sheet body
  {
    OCCShell* occ_shell = occ_body->shell();
    OCCSurface* occ_face = occ_body->my_sheet_surface();
    if (occ_shell != NULL)
      orig_S = occ_shell->get_TopoDS_Shell();

    else
    {
       if (occ_face == NULL)
       {
         PRINT_ERROR("Cannot create an OCC sheet bodySM from the given bodySM.\n");
         return CUBIT_FAILURE;
       }
       orig_S = occ_face->get_TopoDS_Face();
    }
    S = orig_S->EmptyCopied();
    TopoDS_Iterator it(*orig_S);
    while (it.More()) {
      B.Add(S,it.Value().Reversed());
      it.Next();
    }
  }

  else
  {
    S = orig_S->EmptyCopied();
    TopoDS_Iterator it(*orig_S);
    while (it.More()) {
      B.Add(S,it.Value().Reversed());
      it.Next();
    }
    occ_body->set_TopoDS_Shape(TopoDS::CompSolid(S)); 
  }
  
  //Bind the new shape and its underlining sub-shapes.
  TopExp_Explorer Ex_orig, Ex;
  int k = -1;
  Ex.Init(S, TopAbs_COMPSOLID);
  Ex_orig.Init(*orig_S, TopAbs_COMPSOLID);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    if(OCCQueryEngine::instance()->OCCMap->IsBound(Ex.Current()))
    {
      k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());   
      OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
      OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
      TopExp_Explorer Ex_old_solid, Ex_solid;
      Ex_old_solid.Init(*orig_S,TopAbs_SOLID);
      Ex_solid.Init(S, TopAbs_SOLID);
      DLIList<Lump*> lumps = occ_body->lumps();
      for (; Ex_old_solid.More(), Ex_solid.More(); Ex_old_solid.Next(), Ex_solid.Next())
      {
        k = OCCQueryEngine::instance()->OCCMap->Find(Ex_old_solid.Current());
        OCCQueryEngine::instance()->OCCMap->UnBind(Ex_old_solid.Current());
        OCCQueryEngine::instance()->OCCMap->Bind(Ex_solid.Current(), k);
        OCCLump* occ_lump = CAST_TO(lumps.get_and_step(), OCCLump);
        occ_lump->set_TopoDS_Solid(TopoDS::Solid(Ex_solid.Current()));
      }
    } 
    
    else
    {
      Lump *lump = occ_body->lumps().get();
      OCCLump* occ_lump = CAST_TO(lump, OCCLump);
      TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
      k = OCCQueryEngine::instance()->OCCMap->Find(solid);
      OCCQueryEngine::instance()->OCCMap->UnBind(solid);
      TopExp_Explorer Ex_local;
      Ex_local.Init(S, TopAbs_SOLID);
      OCCQueryEngine::instance()->OCCMap->Bind(Ex_local.Current(), k);
      occ_lump->set_TopoDS_Solid(TopoDS::Solid(Ex_local.Current())); 
    }
  }  
      
  Ex.Init(S, TopAbs_SHELL);
  Ex_orig.Init(*orig_S, TopAbs_SHELL);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCShell *shell = (OCCShell*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    shell->set_TopoDS_Shell(TopoDS::Shell(Ex.Current())); 
  }

  Ex.Init(S, TopAbs_FACE);
  Ex_orig.Init(*orig_S, TopAbs_FACE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCSurface *surface = (OCCSurface *)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    TopoDS_Face face = TopoDS::Face(Ex.Current());
    surface->set_TopoDS_Face(face);
  }

  Ex.Init(S, TopAbs_WIRE);
  Ex_orig.Init(*orig_S, TopAbs_WIRE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCLoop* wire = (OCCLoop*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    wire->set_TopoDS_Wire(TopoDS::Wire(Ex.Current()));
  }

  Ex.Init(S, TopAbs_EDGE);
  Ex_orig.Init(*orig_S, TopAbs_EDGE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCCurve* edge = (OCCCurve*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    edge->set_TopoDS_Edge(TopoDS::Edge(Ex.Current()));
  }

  Ex.Init(S, TopAbs_VERTEX);
  Ex_orig.Init(*orig_S, TopAbs_VERTEX);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCPoint* point = (OCCPoint*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    point->set_TopoDS_Vertex(TopoDS::Vertex(Ex.Current()));
  }
  return CUBIT_SUCCESS;
}
    


//===============================================================================
// Function   : split_periodic
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::split_periodic( BodySM * /*body_ptr*/,
                                               BodySM *& /*new_body*/ )
{
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_body
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::regularize_body( BodySM * /*body_ptr*/,
                                                   BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_refentity
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus  OCCModifyEngine::regularize_entity( GeometryEntity * /*old_entity_ptr*/,  
                                                      BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : offset_curves
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::offset_curves( DLIList<Curve*>& curves, 
                                            DLIList<Curve*>& new_curves,
                                            double offset_distance,
                                            const CubitVector& offset_direction, 
                                            int gap_type )
{
  //gap_type has no effect here.
  gp_Dir offset(offset_direction.x(), offset_direction.y(), offset_direction.z()); 
  for(int i = 0 ; i < curves.size(); i++)
  {
    Curve* curve = curves.get_and_step();
    OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    TopoDS_Edge * edge = occ_curve->get_TopoDS_Edge();
    Standard_Real first;
    Standard_Real last;
    Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last); 
    Handle(Geom_OffsetCurve) new_curve = 
      new Geom_OffsetCurve(myCurve,offset_distance, offset);
    if(!new_curve)
    {
      TopologyEntity *entity = curve->topology_entity();
      BasicTopologyEntity *bte = CAST_TO(entity, BasicTopologyEntity);
      PRINT_ERROR("Can't create offset curve for curve %d.\n", bte->id());
      continue;
    }
    *edge = BRepBuilderAPI_MakeEdge(new_curve);
    Curve* offset_curve = OCCQueryEngine::instance()->populate_topology_bridge(*edge);
    new_curves.append(offset_curve);
  } 
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : trim_curve
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::trim_curve( Curve* trim_curve, 
                                    const CubitVector& trim_vector,
                                    const CubitVector& keep_vector,
                                    bool keep )
{
  OCCCurve* occ_crv = CAST_TO(trim_curve, OCCCurve);
  if(!occ_crv) 
  {
    PRINT_ERROR("This is not a OCC curve to be trimmed.\n");
    return (Curve*)NULL;
  }
  
  //Determine the trimmed curve's parameter range.
  double u1, u2;
  occ_crv->get_param_range(u1, u2);
  double trim_u = occ_crv->u_from_position(trim_vector);
  double keep_u = occ_crv->u_from_position(keep_vector);
  if(trim_u > u2+TOL || trim_u < u1 - TOL)
  {
    PRINT_ERROR("The trim_vector is outside of the curve range.\n");
    return (Curve*)NULL;
  }
 
  if(keep_u > trim_u )
     u1 =  trim_u;
  else if(keep_u < trim_u)
     u2 = trim_u;

  else
  {
    PRINT_ERROR("Can't determine which part of the curve to be kept.\n");
    return (Curve*)NULL;
  }
  //get the Geom_Curve of the OCCCurve
  TopoDS_Edge * edge = occ_crv->get_TopoDS_Edge();
  Standard_Real first;
  Standard_Real last;
  Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last);
 
  //Trim the curve
  TopoDS_Edge t_edge = BRepBuilderAPI_MakeEdge(myCurve, u1, u2);
  Curve* t_curve = OCCQueryEngine::instance()->populate_topology_bridge(t_edge);  
  if(!keep)
  {
    DLIList<OCCLoop*> loops = occ_crv->loops();
    if(loops.size() == 0)
      OCCQueryEngine::instance()->delete_solid_model_entities(trim_curve);
  }
  return t_curve;
}

//===============================================================================
// Function   : create_body_from_surfs
// Member Type: PUBLIC
// Description:
// Author     : Jane Hu
// Date       : 4/22/08
//===============================================================================
CubitStatus OCCModifyEngine::create_solid_bodies_from_surfs(DLIList<Surface*> & ref_face_list,
                                          DLIList<BodySM*>& new_bodies,
                                          bool keep_old,
                                          bool heal) const
{
  //keep_old and heal are ignored, always delete old.
  //all surfaces should be stand along surface bodies or shell bodies' surface
  Lump* lump = make_Lump(ref_face_list);
  if (!lump)
    return CUBIT_FAILURE;
  
  new_bodies.append(CAST_TO(lump, OCCLump)->get_body());
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Point* pt1, 
                                          Point* pt2,
                                          Point* pt3, 
                                          bool full )
{ 
  if(!full)
  {
    CubitVector v2(pt2->coordinates());
    return make_Curve(ARC_CURVE_TYPE,pt1,pt3, &v2, CUBIT_FORWARD);
  }
  else
  {
    CubitVector v1(pt1->coordinates());
    CubitVector v2(pt2->coordinates());
    CubitVector v3(pt3->coordinates());

    gp_Pnt gp_pt1(v1.x(),v1.y(), v1.z());
    gp_Pnt gp_pt2(v2.x(),v2.y(), v2.z());
    gp_Pnt gp_pt3(v3.x(),v3.y(), v3.z());

    Handle(Geom_Circle) curve_ptr;
    curve_ptr = GC_MakeCircle(gp_pt1,gp_pt2,gp_pt3); 

    OCCPoint* occ_pt1 = CAST_TO(const_cast<Point*>(pt1),OCCPoint);
    TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt1);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
  }
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Curve* curve1, 
                                          Curve* curve2,
                                          Curve* curve3, 
                                          bool full  )
{ 
  OCCCurve* occ_crv1 = CAST_TO(curve1, OCCCurve);
  OCCCurve* occ_crv2 = CAST_TO(curve2, OCCCurve);
  OCCCurve* occ_crv3 = CAST_TO(curve3, OCCCurve);
  GeometryType type1 = occ_crv1->geometry_type();
  GeometryType type2 = occ_crv2->geometry_type();
  GeometryType type3 = occ_crv3->geometry_type();
  if(type1 != STRAIGHT_CURVE_TYPE || type2 != STRAIGHT_CURVE_TYPE ||
     type3 != STRAIGHT_CURVE_TYPE)
  {
    PRINT_WARNING("Need three straight curves to calculate incenter.\n");
    return (Curve*) NULL;
  } 
  
  //0.check that non of the curves are parallel of each other.
  DLIList<CubitVector*> intscts;
  CubitVector vt1, vt2, vt3;
  CubitBoolean none = CUBIT_FALSE;
  OCCQueryEngine::instance()->get_intersections(curve1, curve2, intscts,none,none);
  vt1 = *intscts.get();
  intscts.clean_out();
  OCCQueryEngine::instance()->get_intersections(curve2,curve3, intscts,none,none);
  vt2 = *intscts.get();
  intscts.clean_out();
  OCCQueryEngine::instance()->get_intersections(curve3, curve1, intscts,none,none);
  vt3 = *intscts.get();

  double u11, u12, u21, u22, u31, u32;
  occ_crv1->get_param_range(u11, u12);
  occ_crv2->get_param_range(u21, u22);
  occ_crv3->get_param_range(u31, u32);

  CubitVector tangent1, tangent2, tangent3;
  occ_crv1->get_tangent(vt1, tangent1);
  occ_crv2->get_tangent(vt2, tangent2);
  occ_crv3->get_tangent(vt3, tangent3); 

  CubitVector normal1 = tangent1 * tangent2;
  CubitVector normal2 = tangent2 * tangent3;
  CubitVector normal3 = tangent3 * tangent1;
  if( normal1.length()< TOL || normal2.length()< TOL ||
      normal3.length() < TOL )
  {
    PRINT_WARNING("Three curves must be able to form a triangle.\n");
    return (Curve*) NULL;
  }

  //normals must parallel to each other, meaning all curves must be on
  //the same plane.
  normal1.normalize();
  normal2.normalize();
  normal3.normalize();
  
  CubitVector parallel1 = normal1 * normal2;
  CubitVector parallel2 = normal2 * normal3;
  CubitVector parallel3 = normal3 * normal1;
  if(parallel1.length() > TOL || parallel2.length() > TOL ||
     parallel3.length() > TOL)
  {
    PRINT_WARNING("Three curves must be able to form a triangle.\n");
    return (Curve*) NULL;
  }
  //1.find the angle between each of the two curves
  double angle1, angle2, angle3;
  angle1 = tangent1.interior_angle(tangent2);
  angle2 = tangent2.interior_angle(tangent3);
  angle3 = tangent3.interior_angle(tangent1);

  //2.create curves to bisection each of the angle passing through the
  // vertices of the triangle
  CubitVector t_curve11 = 
         vectorRotate(angle1/2.0, normal1, tangent1);  
  t_curve11.normalize();
  CubitVector p11 = vt1+t_curve11;

  CubitVector t_curve12 = 
         vectorRotate(90.0 - angle1/2.0, -normal1, tangent1);
  t_curve12.normalize();
  CubitVector p12 = vt1 + t_curve12;

  CubitVector t_curve21 =
         vectorRotate(angle2/2.0, normal2, tangent2);
  t_curve21.normalize();
  CubitVector p21 = vt2 + t_curve21;

  CubitVector t_curve22 = 
         vectorRotate(90.0 - angle2/2.0, -normal2, tangent2);
  t_curve22.normalize();
  CubitVector p22 = vt2 + t_curve22;

  CubitVector t_curve31 = 
         vectorRotate(angle3/2.0, normal3, tangent3);
  t_curve31.normalize();
  CubitVector p31 = vt3 + t_curve31;

  CubitVector t_curve32 =
         vectorRotate(90.0 - angle3/2.0, -normal3, tangent3);
  t_curve32.normalize();
  CubitVector p32 = vt3 + t_curve32;

  //3. find the three intersection points which when connected with the vertices,
  //intersect at same point.
  CubitVector line_p[4], c_p[4], c_ptemp;
  double sc, tc;

  IntersectionTool int_tool;
  for(int i = 0; i < 4; i++)
  {
    if( i == 0)
    {
      line_p[0] = vt1;
      line_p[1] = p11;
      line_p[2] = vt2;
      line_p[3] = p21;
    } 
    else if(i == 1)
      line_p[3] = p22;
    else if(i == 2)
      line_p[1] = p12;
    else
      line_p[3] = p21;  
    int_tool.closest_points_on_segments(line_p[0], line_p[1], line_p[2],
                    line_p[3], c_ptemp, c_p[i], sc, tc);
        
    //check if the closeset point from c_p[i] to three curves are on three curves
    CubitVector closest1, closest2, closest3;
    occ_crv1->closest_point(c_p[i], closest1);
    double u = occ_crv1->u_from_position(closest1);
    if(u > u11-TOL && u < u12 + TOL)
    {
      occ_crv2->closest_point(c_p[i], closest2);
      u = occ_crv2->u_from_position(closest2);
      if(u > u21 - TOL && u < u22 + TOL)
      {
        occ_crv3->closest_point(c_p[i], closest3);
        u = occ_crv3->u_from_position(closest3);
        if(u > u31 - TOL && u < u32 + TOL)
        {
        //4. use the 3 intersection points to find the arc or circle.
          OCCPoint occ_p1(closest1);
          OCCPoint occ_p2(closest2);     
          OCCPoint occ_p3(closest3);
          return create_arc_three(&occ_p1, &occ_p2, &occ_p3, full);
        } 
      }
    }
  }
  PRINT_ERROR("Can't find the tangent points to create circle.\n");
  return (Curve*) NULL;
} 

//===============================================================================
// Function   : create_arc_center_edge
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_center_edge( Point* pt1, 
                                                Point* pt2,
                                                Point* pt3,
                                                const CubitVector& normal, 
                                                double radius,
                                                bool full ) 
{ 
  CubitVector vec1 = pt1->coordinates(); // Center of arc
  CubitVector vec2 = pt2->coordinates(); // Position on arc
  CubitVector vec3 = pt3->coordinates(); // Position on arc

  CubitVector dir1( vec1, vec2 );
  CubitVector dir2( vec1, vec3 );
  // Re-adjust vec2, vec3 if radius was given
  if( radius != CUBIT_DBL_MAX )
  {
     CubitVector vec;
     vec1.next_point( dir1, radius, vec );
     if(vec.distance_between(vec2) > TOL)
     {
       vec2 = vec;
       pt2 = new OCCPoint(vec);
     }
     vec1.next_point( dir2, radius, vec );
     if(vec.distance_between(vec3) > TOL)
     {
       vec3 = vec;
       pt3 = new OCCPoint(vec);
     }
  }
 
  else
  {
    radius = vec1.distance_between(vec2);
    CubitVector vec;
    vec1.next_point( dir2, radius, vec );
    if(vec.distance_between(vec3) > TOL)
    {
       vec3 = vec;
       pt3 = new OCCPoint(vec);
    }
  }
 
  CubitVector normal_dir = normal;
  if(normal_dir.length() > TOL)
  {
    normal_dir.normalize();
    //verify sense
    if((dir1 * dir2) % normal_dir < 0.0)
    {
      Point* p = pt2;
      pt2 = pt3;
      pt3 = p;
    }
    else if((dir1 * dir2) % normal_dir == 0.0)
    {
      PRINT_ERROR("Normal can't be on the plan of the arc.\n");
      return (Curve*) NULL;
    }
  } 
  else 
  {
    normal_dir = dir1 * dir2;
    normal_dir.normalize();
  }
 
  Handle(Geom_Circle) curve_ptr;
  gp_Dir norm(normal_dir.x(), normal_dir.y(), normal_dir.z());
  gp_Pnt center = gp_Pnt( vec1.x(), vec1.y(), vec1.z());
  curve_ptr = GC_MakeCircle(center,norm,radius);

  OCCPoint* occ_pt1 = CAST_TO(const_cast<Point*>(pt2),OCCPoint);
  TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
  if(full)
  {
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt1);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge); 
  }
  else
  {
    Handle(Geom_TrimmedCurve) arc;
    gp_Pnt on_arc1 = gp_Pnt( vec2.x(), vec2.y(), vec2.z());
    gp_Pnt on_arc2 = gp_Pnt( vec3.x(), vec3.y(), vec3.z());
    arc = GC_MakeArcOfCircle(curve_ptr->Circ(), on_arc1, on_arc2, Standard_True);
    OCCPoint* occ_pt2 = CAST_TO(const_cast<Point*>(pt3),OCCPoint); 
    TopoDS_Vertex * vt2 = occ_pt2->get_TopoDS_Vertex();
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(arc, *vt1, *vt2);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
  } 
}

//===============================================================================
// Function   : create_curve_combine
// Member Type: PUBLIC
// Description: create a curve of combination of several curves.  
// Author     : Jane Hu 
// Date       : 03/09
//===============================================================================

CubitStatus 
OCCModifyEngine::create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr )
{
  int i;

  DLIList<OCCCurve*> occ_curves(curve_list.size());
  CAST_LIST( curve_list, occ_curves, OCCCurve );
  if (curve_list.size() != occ_curves.size())
  {
    PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       Not all input curves are OCC Curves.\n");
    return CUBIT_FAILURE;
  }
  
  BRepBuilderAPI_MakeWire aWire(*(occ_curves.get_and_step()->get_TopoDS_Edge()));
  for(i =1 ; i < curve_list.size(); i++)
  {
    OCCCurve* curve = occ_curves.get_and_step();
    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
      aWire.Add(*edge);
    if(!aWire.IsDone())
    {
      PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The curves are not all connected.\n");
      return CUBIT_FAILURE;
    }
  }
  TopoDS_Wire wire = aWire.Wire(); 
  BRepAdaptor_CompCurve comp_curve(wire);
  GeomAbs_CurveType type = comp_curve.GetType();
  GeomAbs_Shape cont = comp_curve.Continuity();
  if(cont < GeomAbs_G1)
  {
    PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The combined curve is not G1 continued.\n");
    return CUBIT_FAILURE;
  }

  //find the start/end vertices for the combined curve.
  double first_u = comp_curve.FirstParameter();
  double last_u = comp_curve.LastParameter();
  gp_Pnt first = comp_curve.Value(first_u);
  gp_Pnt last = comp_curve.Value(last_u);
  CubitVector first_v(first.X(), first.Y(), first.Z());
  CubitVector last_v(last.X(), last.Y(), last.Z());
  if(first_v.about_equal(last_v))
    comp_curve.SetPeriodic(Standard_True);

  DLIList<CubitVector> v_list;
  v_list.append(first_v);
  v_list.append(last_v);
  v_list.reset();

  DLIList<gp_Pnt*> V_list;
  for(int j = 0; j < 2; j++)
  {
    DLIList<TopologyBridge*> children;
    if (j == 0)
      occ_curves.reset();
    else
      occ_curves.last();
    occ_curves.get()->get_children_virt(children);
    CubitVector v = v_list.get_and_step();
    for(i = 0 ; i < children.size(); i++)
    {
      OCCPoint* vertex = CAST_TO(children.get_and_step(), OCCPoint); 
      CubitVector xyz = vertex->coordinates();
      if(xyz.about_equal(v))
      {
        gp_Pnt p ( v.x(), v.y(), v.z());
        V_list.append(&p);
        break;
      }
    }
  }
   
  V_list.reset();
  TopoDS_Edge topo_edge; 
  gp_Lin line;
  gp_Circ circle;
  gp_Elips ellip;
  gp_Hypr hypr;
  gp_Parab parab;
  Handle_Geom_BezierCurve bezier;
  Handle_Geom_BSplineCurve spline;
  switch(type)
  {
    case GeomAbs_Line:
      line = comp_curve.Line();
      topo_edge = BRepBuilderAPI_MakeEdge(line,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Circle:
      circle = comp_curve.Circle();
      topo_edge = BRepBuilderAPI_MakeEdge(circle,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Ellipse:
      ellip = comp_curve.Ellipse();
      topo_edge = BRepBuilderAPI_MakeEdge(ellip,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Hyperbola:
      hypr = comp_curve.Hyperbola();
      topo_edge = BRepBuilderAPI_MakeEdge(hypr,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Parabola:
      parab = comp_curve.Parabola();
      topo_edge = BRepBuilderAPI_MakeEdge(parab,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_BezierCurve:
      bezier = comp_curve.Bezier();
      topo_edge = BRepBuilderAPI_MakeEdge(bezier,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_BSplineCurve:
      spline = comp_curve.BSpline();
      topo_edge = BRepBuilderAPI_MakeEdge(spline,*V_list.get_and_step(), *V_list.get() );
      break;
    default:
      PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The combined curve is not G1 continued.\n");
      return CUBIT_FAILURE;
  }
  TopoDS_Edge *topo_edge_ptr = new TopoDS_Edge(topo_edge);
  OCCCurve* occ_c = new OCCCurve(topo_edge_ptr);
  new_curve_ptr = occ_c;
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_gqe
// Member Type: PUBLIC
// Description: get the facet geometry query engince instance pointer
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
GeometryQueryEngine *OCCModifyEngine::get_gqe()
{
  return OCCQueryEngine::instance();
}

//===============================================================================
// Function   : is_modify_engine
// Member Type: PUBLIC
// Description: return CUBIT_TRUE if the tb_ptr belongs to this modify engine
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitBoolean OCCModifyEngine::is_modify_engine(const TopologyBridge *tb_ptr) const 
{
  return tb_ptr->get_geometry_query_engine() == OCCQueryEngine::instance();
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* curve1, 
                                              Curve* curve2,
                                              DLIList<CubitVector*>& out_list,
                                              double offset,
                                              CubitBoolean ext_first ) 
{
  //offset the curve1 in both directions of normal direction of two curves at
  //center points.
    
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* curve1, 
                                           Surface* face_ptr,
                                           DLIList<CubitVector*> & out_list,
                                           double offset,
                                           CubitBoolean ext_surf )
{
  Surface* new_surface = face_ptr;
  if(ext_surf)
    new_surface = make_Surface(face_ptr, CUBIT_TRUE);

  BodySM* bodysm = NULL;
  CubitStatus status = CUBIT_SUCCESS;
  status = create_offset_surface(new_surface, bodysm, offset); 
  if(status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't offset surface. \n");
    return status;
  }
  OCCSurface* surface = CAST_TO(bodysm, OCCBody)->my_sheet_surface();

  status = OCCQueryEngine::instance()->get_intersections(curve1, surface, out_list);
  
  if(ext_surf || offset)
    OCCQueryEngine::instance()->delete_solid_model_entities(surface);
 
  //offset surface in opposite direction
  if(!offset)
    return status;

  status = create_offset_surface(new_surface, bodysm, -offset);
  if(status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't offset surface. \n");
    return status;
  }
  surface = CAST_TO(bodysm, OCCBody)->my_sheet_surface();
  status = OCCQueryEngine::instance()->get_intersections(curve1, surface, out_list);
  OCCQueryEngine::instance()->delete_solid_model_entities(surface);
  if(ext_surf)
    OCCQueryEngine::instance()->delete_solid_model_entities(new_surface);
  return status;
}

//===============================================================================
// Function   : surface_intersection
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu  
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::surface_intersection( Surface * surface1,
                                                   Surface * surface2,
                                                   DLIList<Curve*> &intscts,
                                                   const double tol) const
{
  OCCSurface *occ_surface1 =  CAST_TO(surface1, OCCSurface);
  if (occ_surface1 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCSurface *occ_surface2 =  CAST_TO(surface2, OCCSurface);
  if (occ_surface2 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }
  
  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(*(occ_surface1->get_TopoDS_Face()),
                                            *(occ_surface2->get_TopoDS_Face()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curve and surface.\n");
      return CUBIT_FAILURE;
    }

  if (distShapeShape.Value() < tol) //intersect
    {
      int numSol = distShapeShape.NbSolution();
      for (int i = 1; i <= numSol; i++)
        {
          TopoDS_Shape shape = distShapeShape.SupportOnShape1(i);
          if(shape.ShapeType() != TopAbs_EDGE)
            continue;

          TopoDS_Edge* edge = new TopoDS_Edge(TopoDS::Edge(shape));
          OCCCurve* cv = new OCCCurve(edge);
          intscts.append(cv);
        }
    }

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_mid_plane
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 01/09
//===============================================================================
CubitStatus OCCModifyEngine::get_mid_plane( const CubitVector & point_1,
                                            const CubitVector & point_2,
                                            const CubitVector & point_3,
                                            BodySM * body_to_trim_to,
                                            DLIList<BodySM*>& midplane_bodies ) const
{
  //Calculate normal of the mid  plane
  CubitVector v1, v2, normal;
  v1 = point_2 - point_1;
  v2 = point_3 - point_1;
  normal = ~(v1 * v2);
  if(fabs(normal.length() - 1) > TOL)
  {
     PRINT_ERROR("The three points are co-linear, and can't be used as a cutting plane.\n");
     return CUBIT_FAILURE;
  }
 
  gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
  gp_Dir normal_dir(normal.x(), normal.y(), normal.z());
  gp_Pln plane(pt, normal_dir);

  TopoDS_Face face = BRepBuilderAPI_MakeFace(plane);
  Surface *surf = OCCQueryEngine::instance()->populate_topology_bridge(face,
                                               CUBIT_TRUE);
  if(!surf)
  {
    PRINT_ERROR("Can't create cutting plane.\n");
    return CUBIT_FAILURE;
  }
  
  BodySM* tool = CAST_TO(surf, OCCSurface)->my_body();
  DLIList<BodySM*> from_bodies;
  from_bodies.append(body_to_trim_to);
   
  CubitStatus stat = intersect(tool, from_bodies, midplane_bodies, 
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);

  return stat;
}

//=============================================================================
// Function   : get_spheric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 spheric surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_spheric_mid_surface( Surface* surface_ptr1,
                                    Surface* surface_ptr2,
                                    BodySM* body_to_trim_to,
                                    DLIList<BodySM *>&midsurface_bodies ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != SPHERE_SURFACE_TYPE ||
     occ_surf2->geometry_type() != SPHERE_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be sphere type.\n");
    return CUBIT_FAILURE;
  }
 
  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  gp_Sphere sphere1 = asurface1.Sphere();
  gp_Sphere sphere2 = asurface2.Sphere();

  gp_Pnt center1 = sphere1.Location();
  gp_Pnt center2 = sphere2.Location();

  if(!center1.IsEqual(center2, TOL))
  {
    PRINT_ERROR( "Spheres need to have the same center.\n");
    return CUBIT_FAILURE;
  }

  double radius = sphere1.Radius()/2.0 + sphere2.Radius()/2.0;
  BodySM* tool = sphere(radius);
  CubitVector center(center1.X(), center1.Y(), center1.Z());
  OCCQueryEngine::instance()->translate(tool, center);

  DLIList<BodySM*> from_bodies;
  from_bodies.append(body_to_trim_to);
  
  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);

  return stat;  
}

//=============================================================================
// Function   : get_conic_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 conic surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_conic_mid_surface( Surface* surface_ptr1,
                                    Surface* surface_ptr2,
                                    BodySM* body_to_trim_to,
                                    DLIList<BodySM *>&midsurface_bodies ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != CONE_SURFACE_TYPE ||
     occ_surf2->geometry_type() != CONE_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be conic type.\n");
    return CUBIT_FAILURE;
  }  

  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  GeomAbs_SurfaceType  type1 = asurface1.GetType();
  GeomAbs_SurfaceType  type2 = asurface2.GetType();
  if(type1 != type2)
  {
    PRINT_ERROR( "Both surfaces provided should be both cylinder or cone type.\n");
    return CUBIT_FAILURE;
  }

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  double height = (bounding_box.diagonal()).length();
  OCCBody* body = CAST_TO(body_to_trim_to, OCCBody);
  CubitVector centroid;
  double volume;
  body->mass_properties(centroid, volume);
  BodySM* tool;
  if(type1 == GeomAbs_Cylinder)
  {
    gp_Cylinder cyl1 = asurface1.Cylinder();
    gp_Cylinder cyl2 = asurface2.Cylinder(); 
    gp_Ax1  axis1 = cyl1.Axis();
    gp_Ax1  axis2 = cyl2.Axis();
    if(!axis1.IsCoaxial(axis2, 0.001, TOL))
    {
      PRINT_ERROR( "Cylinders need to have the same axis of symmetry.\n");
      return CUBIT_FAILURE;
    }
    double radius = cyl1.Radius()/2.0 + cyl2.Radius()/2.0; 
    gp_Ax2 axis;
    axis.SetAxis(axis1);
    TopoDS_Solid S = BRepPrimAPI_MakeCylinder(axis, radius, height);

    Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

    if (lump == NULL)
    {
      PRINT_ERROR("In OCCModifyEngine::get_conic_mid_surface\n"
                  "   Cannot create a cylinder for given radius.\n");
      return CUBIT_FAILURE;
    }

    tool = CAST_TO(lump, OCCLump)->get_body();
    double z = centroid.z();
    z -= height/2.0;
    centroid.z(z);
    OCCQueryEngine::instance()->translate(tool, centroid);  
  }

  else  //GeomAbs_Cone
  {
    gp_Cone cone1 = asurface1.Cone();
    gp_Cone cone2 = asurface2.Cone();
    double angle1 = cone1.SemiAngle();
    double angle2 = cone2.SemiAngle();
    if(fabs(angle1 - angle2) > 0.001)
    {
      PRINT_ERROR( "Cones do not have the same semi-angle.\n");
      return CUBIT_FAILURE;
    }
    gp_Ax1  axis1 = cone1.Axis();
    gp_Ax1  axis2 = cone2.Axis();
    if(!axis1.IsCoaxial(axis2, 0.001, TOL))
    {
      PRINT_ERROR( "Cones need to have the same axis of symmetry.\n");
      return CUBIT_FAILURE;
    } 
    if(axis1.IsOpposite(axis2, 0.001))
    {
      PRINT_ERROR( "Cones need to have the same orientation of axis.\n");
      return CUBIT_FAILURE;
    }
    double r1 = cone1.RefRadius()/2.0 + cone2.RefRadius()/2.0; 
    gp_Ax3 axis;
    axis.SetAxis(axis1);
    gp_Cone cone(axis, angle1, r1);
    TopoDS_Face face = BRepBuilderAPI_MakeFace(cone);
    Surface* surface = 
      OCCQueryEngine::instance()->populate_topology_bridge(face, CUBIT_TRUE);
    tool = CAST_TO(surface,OCCSurface)->my_body();  
  } 
  DLIList<BodySM*> from_bodies;
  from_bodies.append(body_to_trim_to);

  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);

  return stat;
}

//=============================================================================
// Function   : get_toric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 toric surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_toric_mid_surface( Surface* surface_ptr1,
                                     Surface* surface_ptr2,
                                     BodySM* body_to_trim_to,
                                     DLIList<BodySM *>&midsurface_bodies ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != TORUS_SURFACE_TYPE ||
     occ_surf2->geometry_type() != TORUS_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be toric type.\n");
    return CUBIT_FAILURE;
  }

  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  gp_Torus torus1 = asurface1.Torus();
  gp_Torus torus2 = asurface2.Torus();

  gp_Pnt center1 = torus1.Location();
  gp_Pnt center2 = torus2.Location();

  if(!center1.IsEqual(center2, TOL))
  {
    PRINT_ERROR( "Torii need to have the same center.\n");
    return CUBIT_FAILURE;
  }

  double major_r1 = torus1.MajorRadius();
  double major_r2 = torus2.MajorRadius();
  if(fabs(major_r1 - major_r2) > TOL)
  {
    PRINT_ERROR( "Torii need to have the same major radius.\n");
    return CUBIT_FAILURE;
  }

  gp_Ax1 axis1 = torus1.Axis();
  gp_Ax1 axis2 = torus2.Axis();
  if(!axis1.IsCoaxial(axis2, 0.001, TOL))
  {
    PRINT_ERROR( "Torii need to have the same axis of symmetry.\n");
    return CUBIT_FAILURE;
  }

  double radius = torus1.MinorRadius()/2.0 + torus2.MinorRadius()/2.0;
  gp_Ax2 axis;
  axis.SetAxis(axis1);
  TopoDS_Solid S = BRepPrimAPI_MakeTorus(axis, major_r1, radius);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::get_toric_mid_surface\n"
                "   Cannot create a torus for given radii.\n");
    return CUBIT_FAILURE;
  }

  BodySM* tool = CAST_TO(lump, OCCLump)->get_body();

  DLIList<BodySM*> from_bodies;
  from_bodies.append(body_to_trim_to);

  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);

  return stat; 
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : Jane Hu
// Date       : 03/2009
//=============================================================================
CubitStatus OCCModifyEngine::tweak_chamfer( DLIList<Curve*> & curve_list, 
                                            double left_offset,
                                            DLIList<BodySM*> & new_bodysm_list,
                                            double right_offset,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) 
{
  CubitStatus stat;
  int count = 0;
  if(right_offset <= 0.0)
    right_offset = left_offset;

  for(int i = 0; i < curve_list.size(); i++)
  {
    BodySM * new_bodysm_ptr = NULL;
    stat = tweak_fillet(curve_list.get_and_step(), left_offset, right_offset,
                      new_bodysm_ptr , keep_old_body, CUBIT_FALSE, CUBIT_FALSE);
    if(stat && new_bodysm_ptr)
    {
      new_bodysm_list.append_unique(new_bodysm_ptr);
      count = new_bodysm_list.size();
    }
    else
      break;
  }

  if(count == 0)
    return CUBIT_FAILURE;

  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape =
          CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape();
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
      }
    }
    GfxPreview::flush();
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body 
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If 
//              multiple vertices are supplied, only one offset value is 
//              allowed and the edges are not used.
// Author     : Jane Hu
// Date       : 03/09
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_chamfer( DLIList<Point*> & point_list, 
                                double offset1,
                                DLIList<BodySM*> & new_bodysm_list,
                                Curve * edge1,
                                double offset2,
                                Curve * edge2,
                                double offset3,
                                Curve * edge3,
                                CubitBoolean keep_old_body,
                                CubitBoolean preview ) 
{
  // Sort out vertices between sheet and solid bodies
  DLIList<Point*> solid_points, sheet_points;
  DLIList<OCCSurface*> s_list;
  DLIList<OCCBody*> bodies;
  if( sort_points_by_body_type( point_list, solid_points, sheet_points, 
                                s_list, bodies ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if(solid_points.size() > 0 && solid_points.size() != bodies.size())
  {
    PRINT_ERROR( "cannot find bodies corresponding to the points.\n" );
    return CUBIT_FAILURE;
  }

  if(sheet_points.size() > 0 && sheet_points.size() != s_list.size())
  {
    PRINT_ERROR( "cannot find surfaces corresponding to the points.\n" );
    return CUBIT_FAILURE;
  }

  // Do simple forms
  if( edge1 == NULL || offset2 <= 0.0 )
  {
    if( tweak_chamfer_solid( solid_points, bodies, offset1, new_bodysm_list,
      keep_old_body, preview )== CUBIT_FAILURE )
      return CUBIT_FAILURE;
    return tweak_fillet_chamfer_sheet( sheet_points, s_list, offset1, 
           CUBIT_FALSE, new_bodysm_list, keep_old_body, preview );
  }

  if( solid_points.size() > 1 || sheet_points.size() > 1 )
  {
    PRINT_ERROR( "cannot chamfer multiple vertices with a variable radius.\n" );
    return CUBIT_FAILURE;
  }

  if( solid_points.size() )
  {
    Point *point_ptr = solid_points.get();

    if( tweak_chamfer_solid( point_ptr, bodies.get(), offset1, edge1, 
        offset2, edge2, offset3, edge3,
        new_bodysm_list, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }

  if( sheet_points.size() )
  {
    Point *point_ptr = sheet_points.get();

    if( tweak_chamfer_sheet( point_ptr, s_list.get(), offset1, edge1, offset2, 
        edge2, new_bodysm_list, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }
  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::tweak_chamfer_solid( DLIList<Point*> &point_list,
                                    DLIList<OCCBody*> &bodies,
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )
{
  for(int i = 0; i < point_list.size(); i++)
  {
    DLIList<TopologyBridge*> parents;
    OCCPoint* point = CAST_TO(point_list.get_and_step(), OCCPoint);
    OCCBody* body = bodies.get_and_step();
    if(point != NULL)
      point->get_parents_virt(parents); //OCCCurves
    assert(parents.size() == 3);
    DLIList<Curve*> curves;
    for(int j = 0; j < 3; j++)
    {
      OCCCurve* occ_curve = CAST_TO(parents.get_and_step(), OCCCurve); 
      curves.append(occ_curve);
    }
    CubitStatus stat;
    stat = tweak_chamfer_solid(point, body,radius, curves.pop(), radius, 
                        curves.pop(), radius, curves.pop(), new_bodysm_list,
                        keep_old_body, CUBIT_FALSE );
    if(!stat)
      return CUBIT_FAILURE;
  }
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape =
          CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape();
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
      }
    }
    GfxPreview::flush();
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }

  return CUBIT_SUCCESS;

}

CubitStatus
OCCModifyEngine::tweak_chamfer_solid( Point* point_ptr,
                                    OCCBody* body,
                                    double r1,
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    double r3,
                                    Curve *c3,
                                    DLIList<BodySM *> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )
{
  if(r1 <= 0.0 || r2 <= 0.0 || r3 <= 0.0)
  {
    PRINT_ERROR( "Chamfer radii must be greater than zero.\n" );
    return CUBIT_FAILURE;
  }
  
  DLIList<Curve*> curves;
  curves.append(c1);
  curves.append(c2);
  curves.append(c3);
  
  DLIList<double> radii;
  radii.append(r1);
  radii.append(r2);
  radii.append(r3);

  //check point on curves
  OCCPoint* occ_point = CAST_TO(point_ptr, OCCPoint);
  CubitVector position = occ_point->coordinates();
  DLIList<CubitVector> locations;
  for(int i = 0; i < 3; i++)
  {
    OCCCurve *occ_curve = NULL;
    occ_curve = CAST_TO(curves.get_and_step(), OCCCurve);
    double length = occ_curve->measure();

    DLIList<OCCPoint*> point_list;
    occ_curve->get_points(point_list);
    CubitBoolean in = point_list.is_in_list(occ_point); 
    if(!in)
    {
      PRINT_ERROR( "Point is not on curve.\n" );
      return CUBIT_FAILURE;
    }
    //find cutting points on curves
    double u, u1, u2;
    occ_curve->get_param_range(u1,u2);
    u = occ_curve->u_from_position(position); 
    if(fabs(u-u1) < TOL)
      u = occ_curve->u_from_arc_length(u1, radii[i]);
    else
      u = occ_curve->u_from_arc_length(u1, length-radii[i]);
    CubitVector c_p;
    occ_curve->position_from_u(u, c_p);
    locations.append(c_p);
  }

  //decide normal
  CubitVector v1, v2, normal;
  CubitVector point_1 = locations.pop();
  CubitVector point_2 = locations.pop();
  CubitVector point_3 = locations.pop();
  v1 = point_2 - point_1;
  v2 = point_1 - point_3;
  normal = ~(v1 * v2); 
  CubitVector center;
  double volume;
  body->mass_properties(center, volume); 
  CubitVector dir = ~(center - position);
  if(normal % dir > 0.0)//1, 3, 2 order
  {
    CubitVector v = point_2;
    point_2 = point_3;
    point_3 = v;
  }

  DLIList<BodySM*> bodies;
  BodySM* new_body;
  if(keep_old_body || preview)
    new_body = copy_body(body); 
  else
    new_body = body;

  bodies.append(new_body);
  const CubitVector p1 = point_1;
  const CubitVector p2 = point_2;
  const CubitVector p3 = point_3;
  CubitStatus status = section(bodies, p1, p2, p3, 
                               new_bodysm_list, true,false, false);    
  if(!status)
    return CUBIT_FAILURE;

  if(!preview)
    return CUBIT_SUCCESS;

  GfxPreview::clear();

  for(int i = 0; i < new_bodysm_list.size(); i++)
  {
     BodySM* new_bodysm = new_bodysm_list.get_and_step();
     OCCSurface* surf = CAST_TO(new_bodysm, OCCBody)->my_sheet_surface();
     TopoDS_Face* modified_shape = surf->get_TopoDS_Face();
     // Draw this face
     OCCDrawTool::instance()->draw_FACE( modified_shape, CUBIT_BLUE, CUBIT_TRUE );
  }
  GfxPreview::flush();
  OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
  new_bodysm_list.clean_out();

  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::sort_points_by_body_type( DLIList<Point*> &point_list,
                                         DLIList<Point*> &solid_points,
                                         DLIList<Point*> &sheet_points,
                                         DLIList<OCCSurface*> &s_list,
                                         DLIList<OCCBody*> &bodies )
{
  for (int i = 0; i < point_list.size(); i++)
  {
    DLIList<TopologyBridge*> parents;
    OCCPoint* point = CAST_TO(point_list.get_and_step(), OCCPoint);
    int curve_size = 0;
    if(point != NULL)
    {
      point->get_parents_virt(parents); //OCCCurves
      if(parents.size() < 2)
      {
        PRINT_ERROR( "Vertex found not attached to any surfaces.\n" );
        return CUBIT_FAILURE;
      } 
      else if(parents.size() > 3)
      {
        PRINT_ERROR( "Vertex found attached to multiple bodies.\n" );
        return CUBIT_FAILURE;
      }
      curve_size = parents.size();
    }
    
    OCCCurve* occ_curve = CAST_TO(parents.get(), OCCCurve);
    parents.clean_out();
    occ_curve->get_parents_virt(parents); //OCCCoEdges
    if(parents.size() == 0)
    {
      PRINT_ERROR( "Vertex found not attached to any surfaces.\n" );
      return CUBIT_FAILURE;
    }
    OCCCoEdge* coedge = CAST_TO(parents.get(), OCCCoEdge);
    parents.clean_out();
    coedge->get_parents_virt(parents);  //OCCLoops
    assert(parents.size() > 0);
    OCCLoop* loop = CAST_TO(parents.get(), OCCLoop);
    parents.clean_out();
    loop->get_parents_virt(parents); //OCCSurface
    assert(parents.size() > 0);
    OCCSurface* s = CAST_TO(parents.get(), OCCSurface); 
    if(s->my_body() != NULL && curve_size == 2) //sheet body
    {
      s_list.append(s);
      sheet_points.append(point);
    }
    else if(s->my_body() != NULL && curve_size == 3) //shell body
    {
      PRINT_ERROR( "Vertex found attached to multiple surfaces but not on bodies.\n" );
      return CUBIT_FAILURE;
    }
    else
    {
      solid_points.append(point);
      DLIList<OCCBody*> solid_bodies;
      s->get_bodies(solid_bodies);
      assert(solid_bodies.size() == 1);
      bodies += solid_bodies;
    }
  }
  return CUBIT_SUCCESS;
}
//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid 
//              bodies.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( DLIList<Curve*> & curve_list, 
                                           double radius,
                                           DLIList<BodySM*> & new_bodysm_list,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview ) 
{
  CubitStatus stat;
  int count = 0;
  for(int i = 0; i < curve_list.size(); i++) 
  {
    BodySM * new_bodysm_ptr = NULL;
    stat = tweak_fillet(curve_list.get_and_step(), radius, radius, 
                        new_bodysm_ptr , keep_old_body, CUBIT_FALSE);
    if(stat && new_bodysm_ptr)
    {
      new_bodysm_list.append_unique(new_bodysm_ptr);
      count = new_bodysm_list.size();
    }
    else
      break;
  }

  if(count == 0) 
    return CUBIT_FAILURE;
 
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape = 
          CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape();
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
      }
    }
    GfxPreview::flush();
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);  
    new_bodysm_list.clean_out();
  }
  
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid
//              body.  The fillet/chamfer has a variable radius from the
//              start to the end of the curve.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * curve_ptr,
                                           double start_radius,
                                           double end_radius,
                                           BodySM *& new_bodysm_ptr,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview )
{
  return tweak_fillet(curve_ptr, start_radius, end_radius, new_bodysm_ptr,
                      keep_old_body, preview, CUBIT_TRUE);
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: private
// Description: Create a round fillet (or blend) at the given curves on a solid 
//              body.  The fillet/chamfer has a variable radius from the 
//              start to the end of the curve.
// Author     : Jane Hu 
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * curve_ptr, 
                                           double start_radius,
                                           double end_radius,
                                           BodySM *& new_bodysm_ptr,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview,
                                           CubitBoolean if_fillet ) 
{
  //check if this id is valid 
  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  DLIList <OCCBody* > *bodies = oqe->BodyList;  
  DLIList<OCCCurve*> curves;
  for(int j = 0; j <  bodies->size(); j++)
  {
    OCCBody* body = bodies->get_and_step();
    body->get_all_curves(curves);
  }

  bool curve_alive = false;
  for(int j = 0; j <  curves.size(); j++)
  {
    if(curve_ptr == curves.get_and_step())
    {
      curve_alive = true; 
      break;
    }
  }

  if(!curve_alive)
  {
    PRINT_ERROR("This curve is not valid in the current model.\n");
    return CUBIT_FAILURE;
  }

  OCCCurve *occ_curve = CAST_TO(curve_ptr, OCCCurve);
  TopoDS_Edge* topo_edge = occ_curve->get_TopoDS_Edge();

  TopTools_IndexedDataMapOfShapeListOfShape M;
  DLIList<TopoDS_CompSolid*> shape_list;
  TopoDS_Face* s;
  for(int j = 0; j <  bodies->size(); j++)
  {
    OCCBody* body = bodies->get_and_step();
    TopExp_Explorer Ex;
    TopoDS_Shape ashape = *(body->get_TopoDS_Shape());
    M.Clear();
    TopExp::MapShapesAndAncestors(ashape, TopAbs_EDGE, TopAbs_COMPSOLID, M);
    if(!M.Contains(*topo_edge))
      continue;
    DLIList<Lump*> lumps = body->lumps();
    if(lumps.size() > 1)
    {
      PRINT_ERROR("Fillets must be created on single-lump solids.\n");
      return CUBIT_FAILURE;
    }
    shape_list.append_unique(body->get_TopoDS_Shape());
    
    if(!if_fillet) //for chamfer, need to know the face for the curve.
    {
      DLIList<TopologyBridge*> parents;
      occ_curve->get_parents_virt(parents); //OCCCoEdges
      assert(parents.size() > 1);
      OCCCoEdge* coedge = CAST_TO(parents.get(), OCCCoEdge);
      parents.clean_out();
      coedge->get_parents_virt(parents);  //OCCLoops
      assert(parents.size() > 0);
      OCCLoop* loop = CAST_TO(parents.get(), OCCLoop);
      parents.clean_out(); 
      loop->get_parents_virt(parents); //OCCSurface
      assert(parents.size() > 0);
      s = CAST_TO(parents.get(), OCCSurface)->get_TopoDS_Face();
    }
  }
  if(shape_list.size() != 1)
  {
    PRINT_ERROR("Fillets must be created on solids.\n");
    return CUBIT_FAILURE;
  }

  TopoDS_Shape newShape;
  TopoDS_CompSolid* shape = shape_list.get();
  if(keep_old_body)
  {
    BRepBuilderAPI_Copy api_copy(*shape);
    newShape = api_copy.ModifiedShape(*shape);
  }
  else
    newShape = *shape;

  BRepBuilderAPI_MakeShape* fillet;
  if(if_fillet)
  {
    fillet = new BRepFilletAPI_MakeFillet(newShape);
    dynamic_cast<BRepFilletAPI_MakeFillet*>(fillet)->Add(start_radius, end_radius, *topo_edge);
  }
  else
  {
    fillet = new BRepFilletAPI_MakeChamfer(newShape);
    dynamic_cast<BRepFilletAPI_MakeChamfer*>(fillet)->Add(start_radius, end_radius, *topo_edge, *s);
  }
  fillet->Build();

  if(!fillet->IsDone())
  {
    PRINT_ERROR("Can't create fillet on given curve.\n");
    return CUBIT_FAILURE;
  } 
  TopoDS_Shape modified_shape = fillet->Shape();

  if( !preview )
  {
    TopExp_Explorer Ex;
    Ex.Init(newShape, TopAbs_SOLID);
    TopoDS_Solid old_solid = TopoDS::Solid(Ex.Current());
    OCCLump::update_OCC_entity(old_solid , modified_shape, fillet);     
    DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->populate_topology_bridge(modified_shape);
    new_bodysm_ptr = CAST_TO(tbs.get(), BodySM);  
  }
  else
  {
    GfxPreview::clear();

    TopExp_Explorer Ex;
    Ex.Init(modified_shape, TopAbs_FACE); 
    for( ; Ex.More(); Ex.Next() )
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      // Draw this face
      OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
    }

    GfxPreview::flush();
  }
  delete fillet;
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Jane Hu
// Date       : 01/09 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_fillet( DLIList<Point*> & ref_vertex_list, 
                               double radius,
                               DLIList<BodySM*> & new_bodysm_list,
                               CubitBoolean keep_old_body,
                               CubitBoolean preview )
{
  DLIList<OCCSurface*> s_list;
  return tweak_fillet_chamfer_sheet(ref_vertex_list, s_list, radius, CUBIT_TRUE,
         new_bodysm_list, keep_old_body, preview);
}

CubitStatus
OCCModifyEngine::tweak_fillet_chamfer_sheet( DLIList<Point*> & ref_vertex_list,
                               DLIList<OCCSurface*> faces,
                               double radius,
                               CubitBoolean is_fillet,
                               DLIList<BodySM*> & new_bodysm_list,
                               CubitBoolean keep_old_body,
                               CubitBoolean preview )
{
  TopTools_IndexedDataMapOfShapeListOfShape M;

  for(int i = 0; i < ref_vertex_list.size(); i ++)
  {
    Point* pnt = ref_vertex_list.get_and_step();
    OCCPoint* occ_pnt = CAST_TO(pnt, OCCPoint);
    TopoDS_Vertex* vertex = occ_pnt->get_TopoDS_Vertex();
    OCCSurface* face = NULL;

    if( faces.size() == 0)
    {
      OCCQueryEngine* oqe = OCCQueryEngine::instance();
      //make sure the vertex is on sheet body, not on a volume.
      DLIList <OCCBody*> *bodies = oqe->BodyList;
      for(int k =0 ; k < bodies->size(); k++)
      {
        OCCBody* occ_body = bodies->get_and_step();
        TopExp_Explorer Ex;
        TopoDS_Shape ashape = *(occ_body->get_TopoDS_Shape());
        M.Clear();
        TopExp::MapShapesAndAncestors(ashape, TopAbs_VERTEX, TopAbs_COMPSOLID, M);
        if(M.Contains(*vertex))
        {
          PRINT_ERROR("Fillet on vertex can only be done on sheet body.\n");
          return CUBIT_FAILURE;
        }
      } 
      //find corresponding faces.
      DLIList <OCCSurface* > *surfaces = oqe->SurfaceList;
      for(int k =0 ; k < surfaces->size(); k++)
      {
        OCCSurface* occ_face = surfaces->get_and_step();
        TopoDS_Face* topo_face = occ_face->get_TopoDS_Face();
        TopExp_Explorer Ex;
        M.Clear();
        TopExp::MapShapesAndAncestors(*topo_face, TopAbs_VERTEX, TopAbs_FACE, M);
        if(!M.Contains(*vertex))
          continue;
        face = occ_face; 
        break;
      }
    }
    else
      face = faces.get_and_step();

    if(face == NULL)
    {
      PRINT_ERROR("Can't find corresponding surface for the vertex.\n");
      return CUBIT_FAILURE;
    }

    if(!is_fillet)
    {
      //find the two edges sharing the vertex.
      DLIList<OCCCurve*> curves;
      face->get_curves(curves);
      int size = curves.size();
      for(int j = 0; j < size; j ++)
      {
        DLIList<OCCPoint*> point_list;   
        OCCCurve *curve = curves.get();
        curve->get_points(point_list);
        if(!point_list.is_in_list(occ_pnt))
          curves.remove();
        else
          curves.step();
      }
      assert (curves.size()==2);
      tweak_chamfer_sheet(pnt, face, radius, curves.pop(), radius, curves.pop(),
                        new_bodysm_list, keep_old_body, CUBIT_FALSE); 
    }

    else
    {
      TopoDS_Face *shape = face->get_TopoDS_Face();
      TopoDS_Face newShape;
      if(keep_old_body)
      {
        BRepBuilderAPI_Copy api_copy(*shape);
        newShape = TopoDS::Face(api_copy.ModifiedShape(*shape));
      }
      else
        newShape = *shape;

      BRepFilletAPI_MakeFillet2d fillet(newShape);
      TopoDS_Edge fillet_edge;
      fillet_edge = fillet.AddFillet(*vertex, radius);
      fillet.Build();
      if(fillet.Status() != ChFi2d_IsDone)
      {
        PRINT_ERROR("Can't create fillet on given curve.\n");
        return CUBIT_FAILURE;
      }
      TopoDS_Shape modified_shape = fillet.Shape();

      TopExp_Explorer Ex;
      OCCSurface::update_OCC_entity(newShape , modified_shape, &fillet, vertex);
      TopoDS_Face modified_face = TopoDS::Face(modified_shape);
      Surface* surf = OCCQueryEngine::instance()->populate_topology_bridge(modified_face, CUBIT_TRUE);
      BodySM* new_bodysm_ptr = CAST_TO(surf,OCCSurface)->my_body();
      new_bodysm_list.append_unique(new_bodysm_ptr);
    }
  }

  if(!preview )
    return CUBIT_SUCCESS;

  GfxPreview::clear();

  for(int i = 0; i < new_bodysm_list.size(); i++)
  {
     BodySM* new_bodysm = new_bodysm_list.get_and_step();
     OCCSurface* surf = CAST_TO(new_bodysm, OCCBody)->my_sheet_surface();
     TopoDS_Face* modified_shape = surf->get_TopoDS_Face();
     // Draw this face
     OCCDrawTool::instance()->draw_FACE( modified_shape, CUBIT_BLUE, CUBIT_TRUE );
  }
  GfxPreview::flush();
  OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
  new_bodysm_list.clean_out();

  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::tweak_chamfer_sheet(Point* pnt,
                                     OCCSurface* face,
                                     double d1,
                                     Curve* edge1,
                                     double d2,
                                     Curve* edge2,
                                     DLIList<BodySM*> & new_bodysm_list,
                                     CubitBoolean keep_old_body,
                                     CubitBoolean preview ) 
{
  TopoDS_Face *shape = face->get_TopoDS_Face();
  TopoDS_Face newShape;
  if(keep_old_body)
  {
    BRepBuilderAPI_Copy api_copy(*shape);
    newShape = TopoDS::Face(api_copy.ModifiedShape(*shape));
  }
  else
    newShape = *shape;

  BRepFilletAPI_MakeFillet2d fillet(newShape);

  TopoDS_Edge fillet_edge;
  if(edge1 == NULL || edge2 == NULL)
  {
    PRINT_ERROR("Cannot find the two edges for the vertex.\n");
    return CUBIT_FAILURE;
  }
  TopoDS_Edge* topo_e1 = CAST_TO(edge1, OCCCurve)->get_TopoDS_Edge();
  TopoDS_Edge* topo_e2 = CAST_TO(edge2, OCCCurve)->get_TopoDS_Edge();
  TopoDS_Vertex common_v;
  TopExp::CommonVertex(*topo_e1, *topo_e2, common_v);
  fillet_edge = fillet.AddChamfer( *topo_e1, *topo_e2, d1, d2);

  fillet.Build() ;
  if(fillet.Status() != ChFi2d_IsDone)
  {
    PRINT_ERROR("Can't create chamfer on given vertex.\n");
    return CUBIT_FAILURE;
  }
  TopoDS_Shape modified_shape = fillet.Shape();

  if( !preview )
  {
    TopExp_Explorer Ex;
    Ex.Init(newShape, TopAbs_FACE);
    TopoDS_Face old_face = TopoDS::Face(Ex.Current());
    OCCSurface::update_OCC_entity(old_face , modified_shape, &fillet, &common_v);
    DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->populate_topology_bridge(modified_shape);
    BodySM* new_bodysm_ptr = CAST_TO(tbs.get(), BodySM);
    new_bodysm_list.append_unique(new_bodysm_ptr);
  }
  else
  {
    GfxPreview::clear();

    TopExp_Explorer Ex;
    Ex.Init(modified_shape, TopAbs_FACE);
    for( ; Ex.More(); Ex.Next() )
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      // Draw this face
      OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
    }

    GfxPreview::flush();
  }
  return CUBIT_SUCCESS;
}
//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Jane Hu
// Date       : 04/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Surface*> & surface_list, 
                                         const CubitVector & delta,
                                         DLIList<BodySM*> & new_bodysm_list, 
                                         CubitBoolean keep_old_body ,
                                         CubitBoolean preview) const
{
  CubitStatus stat;
  for(int i = 0 ; i < surface_list.size(); i++)
  {
    Surface* surf = surface_list.get_and_step();
    OCCSurface* occ_surf = CAST_TO(surf, OCCSurface);
    if(!occ_surf)
      continue;
    //check to make sure that the surf is not on a sheet body. 
    OCCLump* lump = occ_surf->my_lump(); 
    if(lump != NULL && (lump->my_sheet_surface() || lump->my_shell()))
    {
      PRINT_ERROR( "Cannot tweak move surfaces that are not in a solid\n");
      return CUBIT_FAILURE;
    }
    DLIList<GeometryEntity*> ref_ent_list;
    ref_ent_list.append(occ_surf);
    DLIList<BodySM*> result_bodies;
    stat = sweep_translational(ref_ent_list, result_bodies, delta, 0.0, 1,
                               false, true); 
    if(stat == CUBIT_FAILURE)
    {
      PRINT_ERROR( "Cannot tweak move the surface. \n");
      return CUBIT_FAILURE;
    }
    assert(result_bodies.size() == 1);
    
    //determine if the delta is to trim the existing body or extend it.
    CubitVector center_point;
    center_point = occ_surf->center_point();
    center_point += ~delta;
    BodySM* original_body = occ_surf->my_body();
    if(original_body == NULL)
    {
      DLIList<OCCBody*> original_bodies;
      occ_surf->get_bodies(original_bodies);
      if(original_bodies.size() > 1)
      {
        PRINT_ERROR( "Cannot tweak move the surface in non-mainfold solids. \n");
        return CUBIT_FAILURE;
      }
      else if(original_bodies.size() == 0)
      {
        PRINT_ERROR( "Interal error: Can't find associated solid. \n");       
        return CUBIT_FAILURE;
      }
      original_body = original_bodies.get();
      assert(original_body != NULL);
    } 
    CubitBoolean trim = CUBIT_FALSE;
    if(original_body->point_containment(center_point) == CUBIT_PNT_INSIDE)
       trim = CUBIT_TRUE;

    //subtract or unite the two bodies
    if(!trim)
    {
      result_bodies.insert_first(original_body);
      stat = unite(result_bodies, new_bodysm_list, keep_old_body);
    }
    else
    {
      DLIList<BodySM*> from_bodies;
      from_bodies.append(original_body);
      stat = subtract(result_bodies, from_bodies, new_bodysm_list, CUBIT_FALSE,
                      keep_old_body);
    }
    if(stat = CUBIT_FAILURE)
      return CUBIT_FAILURE;
  } 
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape =
          CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape();
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
      }
    }
    GfxPreview::flush();
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Jane Hu
// Date       : 04/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Curve*> & curves,
                                         const CubitVector & delta,
                                         DLIList<BodySM*> & new_bodysm_list, 
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview ) const
{
  gp_Dir offset_dir(delta.x(), delta.y(), delta.z());
  double length = delta.length();
  
  for(int i = 0 ; i < curves.size(); i++)
  {
    Curve* curve = curves.get_and_step();
    OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    //check to make sure that the curve is on a sheet body.
    DLIList<OCCLoop*> loops;
    loops = occ_curve->loops();
    if(loops.size() == 0)
    {
      PRINT_ERROR( "Cannot tweak move curves that are free\n");
      return CUBIT_FAILURE;
    }
    else if(loops.size() != 1)
    {
      PRINT_ERROR( "Can only tweak move curves attached to one surface\n");
      return CUBIT_FAILURE;
    }
    //determine if the delta is to trim the existing surface or extend it.
    double u_low, u_upper;
    CubitVector a_point;
    occ_curve->get_param_range(u_low, u_upper);
    occ_curve->position_from_u((u_low + u_upper)/2, a_point);
    a_point += ~delta;
    DLIList<TopologyBridge*> parents;
    loops.get()->get_parents_virt(parents);
    OCCSurface* surface = CAST_TO(parents.get(), OCCSurface);
    BodySM* original_body = surface->my_body();
    CubitBoolean trim = CUBIT_FALSE;
    if(surface->point_containment(a_point) == CUBIT_PNT_INSIDE)
       trim = CUBIT_TRUE;

    TopoDS_Edge * edge = occ_curve->get_TopoDS_Edge();
    Standard_Real first;
    Standard_Real last;
    Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last);
    Geom_SurfaceOfLinearExtrusion new_surface(myCurve, offset_dir); 
    Handle(Geom_RectangularTrimmedSurface) trimmed_surface = 
               new  Geom_RectangularTrimmedSurface(&new_surface, first, last,
                                                   0, length);    
    if(trimmed_surface == NULL)
    { 
      PRINT_ERROR( "Can not tweak move the %dth curve\n", i);
      return CUBIT_FAILURE;
    }
    TopoDS_Face FACE = BRepBuilderAPI_MakeFace(trimmed_surface);
    Surface*  extrude_surf= OCCQueryEngine::instance()->populate_topology_bridge(FACE, CUBIT_TRUE);
    BodySM* body = CAST_TO(extrude_surf, OCCSurface)->my_body();
    //subtract or unite the two surfaces
    DLIList<BodySM*> bodies;
    bodies.append(original_body);
    CubitStatus stat;
    if(!trim)
    {
      bodies.append(body);
      stat = unite(bodies, new_bodysm_list, keep_old_body);
    }
    else 
    {
      DLIList<BodySM*> tool_bodies;
      tool_bodies.append(body);
      stat = subtract(tool_bodies, bodies, new_bodysm_list, CUBIT_FALSE, 
                      keep_old_body);
    } 
    if(stat = CUBIT_FAILURE)
      return CUBIT_FAILURE;    
  }
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape =
          CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape();
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE, CUBIT_TRUE );
      }
    }
    GfxPreview::flush();
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Surface*> & /*surface_list*/, 
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Curve*> & /*curve_list*/,  
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the 
//              remaining surfaces to fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Surface*> & /*surface_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*extend_adjoining*/,
                                             CubitBoolean /*keep_surface*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove curves from a sheet body and then extend the 
//              remaining curves or fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Curve*> & /*curve_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a target 
//              surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Surface*> & /*surface_list*/,
                                           DLIList<Surface*> & ,
                                           DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*reverse_flg*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Surface*> & /*target_surfs*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean ,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a target
//              curve that is part of a sheet body.  The target is a surface 
//              created by thickening the owning surface of the target curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Curve*> & /*target_curve_ptr*/, 
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::remove_curve_slivers( BodySM* body,
                                                   double lengthlimit ) const
{
  DLIList<CubitBoolean> is_volume;
  DLIList<BodySM*> bodies;
  DLIList<TopoDS_Shape*> shapes;
  bodies.append(body);
  CubitStatus status = get_shape_list(bodies, shapes, is_volume, CUBIT_FALSE);
  if(!status)
  {
    PRINT_ERROR("Can't find underlying TopoDS_Shape for this body.\n");
    return CUBIT_FAILURE;
  }
  Handle(ShapeBuild_ReShape) context;
  TopoDS_Shape new_shape = ShapeFix::RemoveSmallEdges(*shapes.get(), 
                           lengthlimit, context);
  new_shape = context->Apply(new_shape, TopAbs_COMPSOLID);
  if(context->Status(ShapeExtend_OK))
  {
    PRINT_INFO("There's no small edges on this body.\n");
    return CUBIT_SUCCESS;
  }
     
  else if(context->Status(ShapeExtend_FAIL))
  {
    PRINT_ERROR("Small edges can't be removed from this body.\n");
    return CUBIT_FAILURE;
  }

  OCCQueryEngine::instance()->delete_solid_model_entities(body);
  DLIList<TopologyBridge*>tbs = OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
  
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Surface*>& /*ref_face_list*/, 
                                                   BodySM *& /*new_body*/,
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_u*/, 
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_v*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Curve*>& /*u_curves*/, 
                                                   DLIList<Curve*>& /*v_curves*/,
                                                   BodySM *& /*new_body*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset surface.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_offset_surface( Surface* face_ptr, 
                                                    BodySM*& new_body, 
                                                    double offset ) const
{
  //create offset surface from its center along center normal of distance 
  //"offset"
  OCCSurface *occ_surface =  CAST_TO(face_ptr, OCCSurface);
  if (occ_surface == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    } 

  Surface* c_surface = NULL;
  c_surface = make_Surface(occ_surface);
  if (c_surface == NULL)
  {
    PRINT_ERROR("Cannot copy surface in sweep_translational.\n");
    return CUBIT_FAILURE;
  }
  occ_surface = CAST_TO(c_surface, OCCSurface);

  CubitVector center = occ_surface->center_point();
  CubitVector normal;
  occ_surface->closest_point(center,NULL,&normal); 
  CubitVector v = normal * offset;
  OCCQueryEngine::instance()->translate(occ_surface, v);
  new_body = occ_surface->my_body();
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Creates an offset body.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_offset_body( BodySM* body_ptr, 
                                                 BodySM*& new_bodysm, 
                                                 double offset_distance ) const
{
  PRINT_ERROR("Function not implemented because offset_distance \n"
               "doesn't show offset direction.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a skin surface.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_skin_surface( DLIList<Curve*>& curves, 
                                                  BodySM*& new_body ) const
{
   new_body = NULL;
   Surface* surf = make_Surface(BEST_FIT_SURFACE_TYPE, curves);
   if(surf)
   {
     new_body = CAST_TO(surf, OCCSurface)->my_body();
     return CUBIT_SUCCESS;
   } 
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a shell body from lofting surfaces.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces( Surface * face1, 
                                            const double & /*takeoff1*/,
                                            Surface * face2, 
                                            const double & /*takeoff2*/,
                                            BodySM*& new_body,
                                            CubitBoolean /*arc_length_option*/, 
                                            CubitBoolean /*twist_option*/,
                                            CubitBoolean /*align_direction*/, 
                                            CubitBoolean /*perpendicular*/,
                                            CubitBoolean /*simplify_option*/ ) const
{
   BRepOffsetAPI_ThruSections loft(CUBIT_FALSE);
   CubitStatus stat = do_loft(loft, face1, face2);
   if(!stat)
     return CUBIT_FAILURE;

   TopoDS_Shape shape = loft.Shape();
   TopoDS_Shell shell = TopoDS::Shell(shape);
   OCCShell* occ_shell = OCCQueryEngine::instance()->populate_topology_bridge(shell, CUBIT_TRUE);
   if (occ_shell == NULL)
   {
     PRINT_ERROR("In OCCModifyEngine::loft_surfaces\n"
                 "   Cannot create a loft surface for given surfaces.\n");
     return CUBIT_FAILURE;
   }
   new_body = occ_shell->my_body();
   return CUBIT_SUCCESS;
}

//================================================================================
// Description: Creates a solid body by lofting surfaces between surfaces
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces_to_body( Surface * face1, 
                             const double & /*takeoff1*/,
                             Surface * face2, 
                             const double & /*takeoff2*/,
                             BodySM*& new_body,
                             CubitBoolean /*arc_length_option*/, 
                             CubitBoolean /*twist_option*/,
                             CubitBoolean /*align_direction*/, 
                             CubitBoolean /*perpendicular*/,
                             CubitBoolean /*simplify_option*/ ) const
{
   BRepOffsetAPI_ThruSections loft(CUBIT_TRUE);
   CubitStatus stat = do_loft(loft, face1, face2);
   if(!stat)
     return CUBIT_FAILURE;

   TopoDS_Shape shape = loft.Shape();
   TopoDS_Solid solid = TopoDS::Solid(shape);
   Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(solid);
   if (lump == NULL)
   {
     PRINT_ERROR("In OCCModifyEngine::loft_surfaces_to_body\n"
                 "   Cannot create a loft body for given surfaces.\n");
     return CUBIT_FAILURE;
   }
   new_body = CAST_TO(lump, OCCLump)->get_body();
   return CUBIT_SUCCESS;   
}
 
CubitStatus OCCModifyEngine::do_loft(BRepOffsetAPI_ThruSections& loft,
                                  Surface * face1,
                                  Surface * face2) const
{
   OCCSurface* surf1 = CAST_TO(face1, OCCSurface);
   OCCSurface* surf2 = CAST_TO(face2, OCCSurface);
   if(!surf1 || !surf2)
   {
     PRINT_ERROR("Surfaces are not OCC type.\n");
     return CUBIT_FAILURE;
   }
   TopoDS_Face* topo_face1 = surf1->get_TopoDS_Face();
   TopoDS_Face* topo_face2 = surf2->get_TopoDS_Face();
   TopExp_Explorer Ex;
   Ex.Init(*topo_face1, TopAbs_WIRE);
   TopoDS_Wire wire1 = TopoDS::Wire(Ex.Current());
   if(Ex.More())
   {
     PRINT_ERROR("Surface1 must have only one loop.\n");
     return CUBIT_FAILURE;
   }
   Ex.Init(*topo_face2, TopAbs_WIRE);
   TopoDS_Wire wire2 = TopoDS::Wire(Ex.Current());
   if(Ex.More())
   {
     PRINT_ERROR("Surface2 must have only one loop.\n");
     return CUBIT_FAILURE;
   }
   loft.AddWire(wire1);
   loft.AddWire(wire2);
   loft.Build();
   if(!loft.IsDone())
   {
     PRINT_ERROR("Surfaces can't be loft into a body.\n");
     return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
}  

//================================================================================
// Description: Creates a surface using a list of points, project to surface if
//              given.
// Author     : Jane Hu
// Date       : 03/02/09
//================================================================================
CubitStatus OCCModifyEngine::create_surface( DLIList<CubitVector*>& vec_list, 
                                             BodySM *& new_body, 
                                             Surface * ref_face_ptr,
                                             CubitBoolean ) const
{
   
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a weld surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_weld_surface( CubitVector & /*root*/,
                                                    Surface * /*ref_face1*/, 
                                                    double /*leg1*/, 
                                                    Surface * /*ref_face2*/, 
                                                    double /*leg2*/,
                                                    BodySM *& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::scale( BodySM *&body, const CubitVector& factors )
{
  return OCCQueryEngine::instance()->scale( body, factors );
}

CubitStatus OCCModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                               DLIList<BodySM*> &new_bodies,
                                               DLIList<TopologyBridge*>*,
                                               DLIList<TopologyBridge*>* )  const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

// EOF
