//-------------------------------------------------------------------------
// Filename      : OCCCurve.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/14/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********

// ********** END STANDARD INCLUDES        **********
// ********** BEGIN CUBIT INCLUDES         **********

#include "CastTo.hpp"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "GeometryDefines.h"
#include "OCCCurve.hpp"
#include "GeometryQueryEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "CoEdgeSM.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCPoint.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <TopExp.hxx>
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include <TopTools_IndexedMapOfShape.hxx>
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include <GCPnts_AbscissaPoint.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <Precision.hxx>
#include <Extrema_ExtPC.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "Geom_BezierCurve.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "GeomAdaptor_Curve.hxx"
#include "GCPnts_QuasiUniformAbscissa.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "TopTools_ListOfShape.hxx"
#include "BRepAlgo_NormalProjection.hxx"
#include "TopExp_Explorer.hxx"
#include "GeomLProp_CurveTool.hxx"
#include "GeomAPI_ExtremaCurveCurve.hxx"
#include "Geom_Line.hxx"
#include "Geom_Circle.hxx"
#include "Geom_Ellipse.hxx"
#include "GC_MakeLine.hxx"
#include "gp_Circ.hxx"
#include "gp_Elips.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
//#include "TopOpeBRep_ShapeIntersector.hxx"
//#include "TopOpeBRep_Point2d.hxx"
//#include "TopOpeBRep_EdgesIntersector.hxx"
//#include "TopOpeBRepTool_ShapeTool.hxx"
//#include "BRepPrimAPI_MakePrism.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor
//
// Special Notes :
//
// Creator:      : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
OCCCurve::OCCCurve( TopoDS_Edge *theEdge )
{
  myTopoDSEdge = theEdge;
  myMarked = CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
OCCCurve::~OCCCurve() 
{
  if (myTopoDSEdge)
    delete (TopoDS_Edge *)myTopoDSEdge;
}

void OCCCurve::set_TopoDS_Edge(TopoDS_Edge edge)
{
  if(edge.IsEqual(*myTopoDSEdge))
    return;

  TopoDS_Edge* the_edge = new TopoDS_Edge(edge);
  if(myTopoDSEdge)
    delete (TopoDS_Edge *)myTopoDSEdge;
  myTopoDSEdge = the_edge;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void OCCCurve::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::append_attribute(csa, *myTopoDSEdge); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void OCCCurve::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::remove_attribute(csa, *myTopoDSEdge); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void OCCCurve::remove_all_simple_attribute_virt()
  { OCCAttribSet::remove_attribute(NULL, *myTopoDSEdge); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus OCCCurve::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                               csa_list)
  { return OCCAttribSet::get_attributes(*myTopoDSEdge, csa_list); }
  
CubitStatus OCCCurve::get_simple_attribute( const CubitString& name,
                                      DLIList<CubitSimpleAttrib*>& csa_list)
  { return OCCAttribSet::get_attributes( name, *myTopoDSEdge, csa_list ); }

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCCurve::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox OCCCurve::bounding_box() const 
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  Bnd_Box aBox;
  BndLib_Add3dCurve::Add(acurve, Precision::Approximation(), aBox);
  double min[3], max[3];
  aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
  return CubitBox(min, max);
}


//-------------------------------------------------------------------------
// Purpose       : Return the length of the curve.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double OCCCurve::measure()
{
  GProp_GProps myProps;
  BRepGProp::LinearProperties(*myTopoDSEdge, myProps);
  return myProps.Mass();
}

//-------------------------------------------------------------------------
// Purpose       : Return the arc length along the Curve starting from
//                 the point represented by the parameter1 going to the 
//                 point represented by parameter2.
//
// Special Notes : The sign of the returned length value is always positive.
//                 Parameter1 and parameter2 are with respect to the EDGE.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double OCCCurve::length_from_u( double parameter1, double parameter2 )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  return GCPnts_AbscissaPoint::Length(acurve, parameter1, parameter2);
}

//-------------------------------------------------------------------------
// Purpose       : Returns CUBIT_TRUE and the associated period value. Not
//                 implemented yet
//
// Special Notes :  
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitBoolean OCCCurve::is_periodic(double& period)
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  if (acurve.IsPeriodic())
  {
    period = acurve.Period();
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//------------------------------------------------------------------
// Purpose: Returns CUBIT_TRUE and the associated parametric values, 
//          if the facet curve associated with the first EDGE is 
//          parametric.
//          Otherwise returns CUBIT_FALSE and the values of 
//          the lower and upper parametric bounds are undetermined.
//          NOT IMPLEMENTED YET
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------
CubitBoolean OCCCurve::get_param_range( double& lower_bound,
                                        double& upper_bound )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  lower_bound = acurve.FirstParameter();
  upper_bound = acurve.LastParameter();
  return CUBIT_TRUE;
}


//------------------------------------------------------------------
// Purpose:        Finds the extrema along this Curve. 
//
// Special Notes : It is the responsibility of the
//                 calling code to delete the CubitVectors added to 
//                 interior_points!
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/29/01
//-------------------------------------------------------------------
CubitStatus OCCCurve::get_interior_extrema(
  DLIList<CubitVector*>& interior_points,
  CubitSense& return_sense )
{
  // Danilov: try to use GeomAPI_ExtremaCurveCurve
  // Will do 3 primary directions seperately. 
  assert(0);
  
  DLIList<CubitVector*> point_list;
  CubitVector x(1.0, 0.0, 0.0);
  get_interior_extrema_in_direction(point_list, x);
  CubitVector y(0.0, 1.0, 0.0);
  get_interior_extrema_in_direction(point_list, y);
  CubitVector z(0.0, 0.0, 1.0);
  get_interior_extrema_in_direction(point_list, z );

  //like ACIS, return only points that aren't at an endpoint and are not 
  //close to previous point
  const double epsilon = 30.* GEOMETRY_RESABS;
  const double epsilon_squared = epsilon*epsilon;

  //get both vertices' coordinates.
  CubitVector endpoints[2];
  int i = 0;
  TopExp_Explorer aVertexExp(*myTopoDSEdge, TopAbs_VERTEX);
  while(aVertexExp.More())
  {
     TopoDS_Vertex v = TopoDS::Vertex(aVertexExp.Current());
     gp_Pnt p = BRep_Tool::Pnt(v);
     endpoints[i].x(p.X());
     endpoints[i].y(p.Y());
     endpoints[i].z(p.Z());  
     i++;
     aVertexExp.Next();
  } 

  //compare to see if the Points in point_list are interior and far apart
  int j;
  CubitVector* cubit_position = NULL;
  CubitVector * temp_position = NULL;
  point_list.sort();
  point_list.reset();
  for (j = point_list.size(); j--; )
  {
     temp_position = point_list.get_and_step();
     // save if not equal to an endpoint, or prior point
     if (temp_position->distance_between_squared(endpoints[0]) > epsilon_squared
        &&
        temp_position->distance_between_squared(endpoints[1]) > epsilon_squared)
     {
        if (!cubit_position ||
            temp_position->distance_between_squared(*cubit_position) > epsilon_squared)
        {
          cubit_position = temp_position ;
          interior_points.append( cubit_position );
        } // If point isn't close to previous point
     } // If point isn't at an endpoint
  } // for each point

  // Return sense is whatever the sense of this curve is.
  TopAbs_Orientation sense = myTopoDSEdge->Orientation();
  return_sense = (sense == TopAbs_FORWARD ? CUBIT_FORWARD : CUBIT_REVERSED);
 
  return CUBIT_SUCCESS;
}

CubitStatus OCCCurve::get_interior_extrema_in_direction(
		 	DLIList<CubitVector*>& interior_points,
			CubitVector dir)
{
  //Create a straight line.
  gp_Pnt origin(0.0, 0.0, 0.0);
  gp_Dir adir(dir.x(), dir.y(), dir.z());
  Handle(Geom_Line) line = GC_MakeLine(origin, adir);

  //get the Geom_Curve of the OCCCurve
  Standard_Real first;
  Standard_Real last;
  Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);
  
  GeomAPI_ExtremaCurveCurve extrema(myCurve, line);
  int nPnt = extrema.NbExtrema();
  for (int i = 1; i <= nPnt ; i++)
  {
    gp_Pnt P1, P2;
    extrema.Points(i, P1, P2);
    CubitVector* v = new CubitVector(P1.X(), P1.Y(), P1.Z());
    interior_points.append(v);
  }
  return CUBIT_SUCCESS;
}  
//-------------------------------------------------------------------------
// Purpose       : This function computes the point on the curve closest 
//                 to the input location.  Optionally, it can also compute
//                 the tangent and curvature on the Curve at the point on
//                 on the Curve closest to the input location.
//
// Special Notes : The tangent direction is always in the positive direction of the 
//                 owning RefEdge, regardless of the positive direction of the
//                 underlying solid model entities.
//
//                 If the calling code needs the tangent and/or the curvature,
//                 it is responsible for allocating the memory for these
//                 CubitVector(s) and sending in the relevant non-NULL
//                 pointers to this routine.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus OCCCurve::closest_point( 
  CubitVector const& location, 
  CubitVector& closest_location,
  CubitVector* tangent_ptr,
  CubitVector* curvature_ptr,
  double* param)
{  
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  Extrema_ExtPC ext(p, acurve, Precision::Approximation());
  if (!ext.IsDone())
    return CUBIT_FAILURE;

  double sqr_dist = CUBIT_DBL_MAX;
  double pparam;
  for (int i = 1; i <= ext.NbExt(); ++i) {
    if (ext.IsMin(i)) {
      double new_sqr_dist = p.SquareDistance( ext.Point(i).Value() );
      if (new_sqr_dist < sqr_dist) {
        sqr_dist = new_sqr_dist;
        newP = ext.Point(i).Value();
        pparam = ext.Point(i).Parameter();
      }
    }
  }

    // if we didn't find any minimum...
  if (sqr_dist == CUBIT_DBL_MAX)
    return CUBIT_FAILURE;

    // pass back closest point
  closest_location.set( newP.X(), newP.Y(), newP.Z() );
  if (param != NULL)
    *param = pparam;

    // pass back tangent
  if (tangent_ptr != NULL) {
    BRepLProp_CLProps CLP(acurve, 2, Precision::PConfusion());
    CLP.SetParameter( pparam );
    if (!CLP.IsTangentDefined())
      return CUBIT_FAILURE;

    gp_Dir tangent;
    CLP.Tangent(tangent);
    tangent_ptr->set( tangent.X(), tangent.Y(), tangent.Z() );
  }

    // if requested, pass back curvature
  if (curvature_ptr != NULL)
    get_curvature( closest_location, *curvature_ptr );

  return CUBIT_SUCCESS;
}


//------------------------------------------------------------------
// Purpose: This function returns the coordinate of a point in the local
//          parametric (u) space that corresponds to the input position 
//          in global (world) space.  The input point is first moved to 
//          the closest point on the Curve and the parameter value of 
//          that point is determined. 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------
CubitStatus OCCCurve::position_from_u (double u_value,
                                        CubitVector& output_position)
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  gp_Pnt p = acurve.Value(u_value);
  output_position.x(p.X());
  output_position.y(p.Y());
  output_position.z(p.Z());
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the coordinate of a point in the local
//                 parametric (u) space that corresponds to the input position 
//                 in global (world) space.  The input point is first moved to 
//                 the closest point on the Curve and the parameter value of 
//                 that point is determined. 
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
double OCCCurve::u_from_position (const CubitVector& input_position)
{
    // Get the closest point on the Curve to the input position
  CubitVector closest_point;
  double u_val;
  this->closest_point(input_position, closest_point,
                      NULL, NULL, &u_val);
    // closest_point already makes adjustments for sense and periodicity
  
  return u_val;
}

//------------------------------------------------------------------
// Purpose: This function returns the parameter value of the point 
//          that is "arc_length" away from the root point, in the
//          positive sense direction of the owning RefEdge.
//
// Special Notes : 
//   If arc_length is negative, the new point (whose parameter value
//   is being computed) is in the negative sense direction (along
//   the RefEdge) from the root point (whose parameter value is
//   root_param).
//
//   If the curve is not periodic and the new point, "arc_length"
//   away from the root point in the appropriate direction, goes
//   beyond the end point of the first EDGE, that end point is used
//   to generate the returned parameter value.
//
// If the curve is periodic and the new point, "arc_length" away
// from the root point in the appropriate direction, goes beyond
// the end point of the first EDGE, wrap around is done.  After
// wrap around, the point is treated as with other curves
//
// NOTE:
// The important assumption that is made in this routine is that
// the end points of the RefEdge that owns this CurveACIS are the
// same as the end points of the first ACIS EDGE in the list of EDGEs
// associated with this CurveACIS.
//
// Assume that the parameter root_param is with respect to the
// RefEdge as well as arc_length.  Before calling the ACIS "curve",
// we need to get them with respect to the curve.   
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/28/97
//------------------------------------------------------------------
double OCCCurve::u_from_arc_length ( double root_param,
                                       double arc_length )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  GCPnts_AbscissaPoint abs(acurve, arc_length, root_param);
  if (abs.IsDone()) return abs.Parameter();
  else return 0.0;
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying curve. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitBoolean OCCCurve::is_position_on( const CubitVector &test_position )
{
  CubitVector new_point;
  CubitStatus stat = closest_point(test_position, new_point, NULL,NULL,NULL);

  if ( !stat )
    return CUBIT_FALSE;
  CubitVector result_vec = test_position - new_point;
  if ( result_vec.length_squared() < GEOMETRY_RESABS )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the type of underlying curve. 
//
// Special Notes : It checks to see if *any* of the ACIS curves associated
//                 with the EDGEs in the list of EDGEs of this Curve is of
//                 a particular type and returns the appropriate value
//                 of the enum, CurveType.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
GeometryType OCCCurve::geometry_type()
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  if (acurve.GetType() == GeomAbs_BezierCurve)
     return SPLINE_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_BSplineCurve)
     return BSPLINE_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Line)
     return STRAIGHT_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Parabola)
     return PARABOLA_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Hyperbola)
     return HYPERBOLA_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Circle)
      return ARC_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Ellipse)
     return ELLIPSE_CURVE_TYPE;
  return UNDEFINED_CURVE_TYPE;
}

//-------------------------------------------------------------------------
// Purpose       : Return direction of point on curve
//
// Special Notes : Finds the underlying line's origin and direction unit vector
//
// Creator       : Jane Hu
//
// Creation Date : 11/14/07
//-------------------------------------------------------------------------
CubitStatus OCCCurve::get_point_direction( CubitVector& point, 
                                           CubitVector& direction )
{
  if (geometry_type() != STRAIGHT_CURVE_TYPE)
    return CUBIT_FAILURE;

  //get the underlying geometry curve
  double first,last;
  Handle(Geom_Curve) gCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);

  //get the origin and direction of the underlying curve
  Handle(Geom_Line) gLine = Handle(Geom_Line)::DownCast(gCurve);
  gp_Ax1 axis = gLine->Position();
  gp_Pnt loc = axis.Location();
  gp_Dir dir = axis.Direction();
  point.set(loc.X(), loc.Y(), loc.Z());

  //Based on the TopoDS_Edge's orientation, give the unit vector.
  if (myTopoDSEdge->Orientation() == TopAbs_FORWARD)
    direction.set(dir.X(), dir.Y(), dir.Z());
  else if(myTopoDSEdge->Orientation() == TopAbs_REVERSED)
    direction.set(-dir.X(), -dir.Y(), -dir.Z());
  return CUBIT_SUCCESS;  
}

//-------------------------------------------------------------------------
// Purpose       : Return the center and radius of an arc
//
// Special Notes : not currently implemented
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus OCCCurve::get_center_radius( CubitVector& center, 
                                           double& radius )
{
  if( geometry_type() != ELLIPSE_CURVE_TYPE &&
      geometry_type() != ARC_CURVE_TYPE )
    return CUBIT_FAILURE;
 
  //get the Geom_Curve of the OCCCurve
  Standard_Real first;
  Standard_Real last;
  Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last); 

  if (Handle(Geom_Circle) gCircle = Handle(Geom_Circle)::DownCast(myCurve))
  {
     radius = gCircle->Radius();
     gp_Circ gp_circ = gCircle->Circ();
     gp_Pnt  gp_p = gp_circ.Location();
     center.set(gp_p.X(), gp_p.Y(), gp_p.Z());
  }

  else //ellipse
  {
     Handle(Geom_Ellipse) gEllipse = Handle(Geom_Ellipse)::DownCast(myCurve);
     radius = gEllipse->MajorRadius();
     gp_Elips gp_ellip = gEllipse->Elips();
     gp_Pnt  gp_p = gp_ellip.Location();
     center.set(gp_p.X(), gp_p.Y(), gp_p.Z());
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the start parameter.
//
// Special Notes : The start param is with respect to the ref_edge.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double OCCCurve::start_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return start;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the end parameter.
//
// Special Notes : The end param is with respect to the ref_edge.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double OCCCurve::end_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return end;
}


void OCCCurve::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
{ 
   for(int i = 0; i < myLoopList.size(); i++) 
   {
      DLIList<OCCCoEdge*> coedges = myLoopList.get_and_step()->coedges();
    
      for(int j = 0; j < coedges.size(); j++)
      {
        OCCCoEdge * coedge = coedges.get_and_step();
	if(coedge->curve() == this)
	{
	  parents.append(coedge);
          break;
 	}
      }
   }
}
void OCCCurve::get_children_virt( DLIList<TopologyBridge*>& children ) 
{
	TopTools_IndexedMapOfShape M;
	TopExp::MapShapes(*myTopoDSEdge, TopAbs_VERTEX, M);
	TopologyBridge *point1, *point2;
	if (M.Extent()==1) {
		point1 = OCCQueryEngine::instance()->occ_to_cgm(M(1));
                if (point1)
		  children.append_unique(point1);
	} else if (M.Extent()==2) {
		if (  fabs(BRep_Tool::Parameter(TopoDS::Vertex(M(1)), *myTopoDSEdge)-start_param()) > 
				fabs(BRep_Tool::Parameter(TopoDS::Vertex(M(2)), *myTopoDSEdge)-start_param())  ) {
			point1 = OCCQueryEngine::instance()->occ_to_cgm(M(2));
			point2 = OCCQueryEngine::instance()->occ_to_cgm(M(1));
		} else {
			point1 = OCCQueryEngine::instance()->occ_to_cgm(M(1));
			point2 = OCCQueryEngine::instance()->occ_to_cgm(M(2));
		}
		if (point1 && point1 == point2) {
			children.append_unique(point1);
		} else {
                        if (point1)
			  children.append_unique(point1);
                        if (point2)
			  children.append_unique(point2);
		}
	}
}
 



//-------------------------------------------------------------------------
// Purpose       : Check for G1 discontinuity
//
// Special Notes : returns tangency discontinuity all along the Curve
//		   at the param, only returns minus tangent = plus tangent
//		   when it's C1 continuity.
//
// Creator       : Jane Hu
//
// Creation Date : 11/14/07
//-------------------------------------------------------------------------
CubitBoolean OCCCurve::G1_discontinuous( 
      double param, CubitVector* mtan, CubitVector* ptan )
{ 
  CubitBoolean is_discon = CUBIT_TRUE;
  double first, last;
  Handle(Geom_Curve) gCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);

  if (gCurve->Continuity() < GeomAbs_G1)
     return is_discon;

  assert(first <= param && param <= last );
  
  gp_Pnt P;
  gp_Vec V1;
  gCurve->D1(param, P, V1);
  
  mtan = new CubitVector(V1.X(), V1.Y(),V1.Z());
  ptan = new CubitVector(*mtan);
     
  return CUBIT_FALSE;
}

void OCCCurve::get_points( DLIList<OCCPoint*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSEdge, TopAbs_VERTEX, M);
  int ii;
  for (ii=M.Extent(); ii>0; ii--) {
	  TopologyBridge *point = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          if (point)
	    result_list.append_unique(dynamic_cast<OCCPoint*>(point));
  }
}

void OCCCurve::get_tangent( CubitVector const& location,
                            CubitVector& tangent)
{
    double u = u_from_position(location);
    Standard_Real first;
    Standard_Real last;
    Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);
 
    gp_Pnt p;
    gp_Vec tan;
    GeomLProp_CurveTool::D1(myCurve, u , p, tan) ;
    tangent.set(tan.X(), tan.Y(), tan.Z());
} 

void OCCCurve::get_curvature( CubitVector const& location,
                              CubitVector& curvature)
{  
    double u = u_from_position(location);
    Standard_Real first;
    Standard_Real last;
    Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);

    gp_Pnt p;
    gp_Vec tan, cur;
    GeomLProp_CurveTool::D2(myCurve, u , p, tan, cur) ;
    curvature.set(cur.X(), cur.Y(), cur.Z()); 
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********


//----------------------------------------------------------------
// Adjusts the input parameter so that it falls within the
// parameter range of this Curve, if possible.  Necessary for
// periodic curves.
//----------------------------------------------------------------
void OCCCurve::adjust_periodic_parameter(double& param)
{
  assert(0);
    // Adjustment only legal if this is a periodic curve.
  double period;
  if ( this->is_periodic(period) && (fabs(period) > CUBIT_RESABS))
  {
    double upper_bound, lower_bound;
    this->get_param_range( lower_bound, upper_bound );
    double edge_range = upper_bound - lower_bound;
		assert( edge_range > CUBIT_RESABS * 100 );

    lower_bound -= CUBIT_RESABS;
    upper_bound += CUBIT_RESABS;
    
      // Make sure period is positive
    if (period < 0.)
      period = -period;

      // Move the parameter above the low param
    while (param < lower_bound)
      param += period;
      // Move the parameter below the high param
    while (param > upper_bound)
      param -= period;
  }
}

CubitPointContainment OCCCurve::point_containment( const CubitVector &point )
{
   if (is_position_on(point) == CUBIT_TRUE)
   {
     DLIList<OCCPoint*> points;
     get_points(points);
     for (int i = 0; i < points.size(); i++)
     {
	OCCPoint* pnt = points.get_and_step();
        CubitVector v = pnt->coordinates();
        double d = v.distance_between(point); 
        if (d < GEOMETRY_RESABS)
	  return CUBIT_PNT_BOUNDARY; 
     }
     return CUBIT_PNT_ON;
   } 
   return CUBIT_PNT_OFF;
}

//----------------------------------------------------------------
// Function: to update the core Curve
//           for any movement of the body/surface/curve.
// Author: Jane Hu
//----------------------------------------------------------------
void OCCCurve::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                 BRepAlgoAPI_BooleanOperation *op)
{
  if (myMarked == 1) 
     return;

  assert(aBRepTrsf != NULL || op != NULL);
  
  TopoDS_Shape shape;
  if (aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Edge());
  else
  {
    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(*get_TopoDS_Edge()));
    if(shapes.Extent() == 0)
      shapes.Assign(op->Generated(*get_TopoDS_Edge()));
    if(shapes.Extent() == 1)
      shape = shapes.First();
    else if(shapes.Extent() > 1)
    {
      //update all attributes first.
      TopTools_ListIteratorOfListOfShape it;
      it.Initialize(shapes);
      for(; it.More(); it.Next())
      {
        shape = it.Value();
        OCCQueryEngine::instance()->copy_attributes(*get_TopoDS_Edge(), shape);
      }
      shape = shapes.First();
    }
    else if (op->IsDeleted(*get_TopoDS_Edge()))
      ;
    else
      return ;
  }
  TopoDS_Edge curve;
  if(!shape.IsNull())
    curve = TopoDS::Edge(shape);

  //set the vertices
  DLIList<TopologyBridge*> vertices;
  get_children_virt(vertices);
  for (int i = 1; i <= vertices.size(); i++)
  {
     TopologyBridge* tb = vertices.get_and_step();
     OCCPoint *point = CAST_TO(tb, OCCPoint);
     if (point)
       point->update_OCC_entity(aBRepTrsf, op);
  }
  myMarked = 1;
  OCCQueryEngine::instance()->update_OCC_map(*myTopoDSEdge, curve);
}

//===============================================================================
// Function   : project_curve
// Member Type: PUBLIC
// Description: project a curve onto a surface, if closed is true,
//              make sure it projected as two segment, then combine them
//              into a closed shape, third_point is used to determine
//              which segment to use if having two projections.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCCurve::project_curve(Surface* face_ptr, 
                               DLIList<Point*>&  normal_proj_points,
                               CubitBoolean closed,
                               const CubitVector* third_point)
{
   TopoDS_Edge* edge = get_TopoDS_Edge();
   if (edge == NULL)
   {
        PRINT_ERROR("Cannot project the curve .\n"
                 "Possible incompatible geometry engines.\n");
        return (Curve*) NULL;
   }

   TopoDS_Face* face = CAST_TO(face_ptr, OCCSurface)->get_TopoDS_Face();
   if(face == NULL)
   {
        PRINT_ERROR("Cannot project the curve to the surface.\n"
                 "Possible incompatible geometry engines.\n");
        return (Curve*) NULL;
   }

   BRepAlgo_NormalProjection aProjection;
   aProjection.Init(*face);
   aProjection.Add(*edge);
   aProjection.Build();
   if (!aProjection.IsDone())
   {
        PRINT_ERROR("Cannot project the curve to the surface.\n"
                 "OCC engine failure.\n");
        return (Curve*) NULL;
   }

   TopoDS_Shape new_shape = aProjection.Projection();//compound shape
   int num_projection = 0;
   if (new_shape.IsNull())
   {
       PRINT_ERROR("Cannot project the curve to the surface.\n");
       return (Curve*) NULL;
   }

   else
   {
     //count how many free edges and vertices the new_shape has.
     TopExp_Explorer Ex;
     for (Ex.Init(new_shape,TopAbs_EDGE); Ex.More(); Ex.Next())
       num_projection++;
     for (Ex.Init(new_shape,TopAbs_VERTEX, TopAbs_EDGE); Ex.More(); Ex.Next())
       num_projection++;
   }

   if(num_projection == 0)
   {
      PRINT_INFO("No projection on the surface.\n");
      return (Curve*) NULL;
   }

   else if ( num_projection == 1 )
   {
      if(closed == true)
       PRINT_WARNING("Cannot project the curve to create a closed projection.\n"                 "There is only one projection segment.\n");

      TopExp_Explorer Ex;
      TopoDS_Edge new_edge;
      TopoDS_Vertex new_point;
      for (Ex.Init(new_shape,TopAbs_EDGE); Ex.More(); Ex.Next())
      {
        new_edge = TopoDS::Edge(Ex.Current());
        return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
      }
      for(Ex.Init(new_shape,TopAbs_VERTEX);Ex.More(); Ex.Next())
      {
        new_point = TopoDS::Vertex(Ex.Current());
        normal_proj_points.append(OCCQueryEngine::instance()->populate_topology_bridge(new_point));
      } 
      return (Curve*) NULL;
   }

   else if (num_projection == 2)
   {
      double d;
      double first, last;
      TopExp_Explorer Ex;
      TopoDS_Edge edge1, edge2;
      TopoDS_Vertex point;

      int count = 0;
      for (Ex.Init(new_shape,TopAbs_EDGE); Ex.More(); Ex.Next())
      {
        count++;
        if(count == 1)
          edge1 = TopoDS::Edge(Ex.Current());
        if(count == 2)
          edge2 = TopoDS::Edge(Ex.Current());
      }

      for(Ex.Init(new_shape,TopAbs_VERTEX);Ex.More(); Ex.Next())
      {
        point = TopoDS::Vertex(Ex.Current());
        normal_proj_points.append(OCCQueryEngine::instance()->populate_topology_bridge(point));
      }

      if(edge1.IsNull())
        return OCCQueryEngine::instance()->populate_topology_bridge(edge2);

      if(edge2.IsNull())
        return OCCQueryEngine::instance()->populate_topology_bridge(edge1);

      if(edge1.IsNull() && edge2.IsNull())
        return (Curve*) NULL;

      Handle(Geom_Curve) myCurve1 =
                        BRep_Tool::Curve(edge1,first,last);

      Handle(Geom_Curve) myCurve2= BRep_Tool::Curve(edge2, first, last);
      //If the surface is periodic, so it has 2 projections, we just need to
      //find the segment to which the third_point is closer.
      if(closed == CUBIT_FALSE && third_point != NULL)
      {
        gp_Pnt P (third_point->x(), third_point->y(), third_point->z());
        GeomAPI_ProjectPointOnCurve projOncurve(P, myCurve1);
        if (projOncurve.NbPoints() == 0)
        {
          PRINT_ERROR("Cannot project the curve to the surface.\n"
                 "OCC engine failure.\n");
          return (Curve*) NULL;
        }
        d = projOncurve.LowerDistance();

        //Compare with the second solution
        GeomAPI_ProjectPointOnCurve projOncurve2(P, myCurve2);
        if (projOncurve2.NbPoints() == 0)
        {
           PRINT_ERROR("Cannot project the curve to the surface.\n"
                 "OCC engine failure.\n");
          return (Curve*) NULL;
        }

        double d2 = projOncurve2.LowerDistance();
        TopoDS_Edge new_edge =
                d > d2 ? edge2 : edge1 ;
        return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
      }


      else if (closed == CUBIT_TRUE)
      {
        //connect the two segment into a closed shape. Assume both segment
        // has the same curve type, create Bezier closed curve.
        GeomAdaptor_Curve acurve1(myCurve1);
        GeomAdaptor_Curve acurve2(myCurve2);
        //get 10 points of each curve, combine them to make one Bezier curve
        int NbPoints = 10;
        GCPnts_QuasiUniformAbscissa distribution1(acurve1, NbPoints);
        GCPnts_QuasiUniformAbscissa distribution2(acurve2, NbPoints);
        TColgp_Array1OfPnt points(1, 2*NbPoints-1);
        int i;
        for (i = 1; i <= NbPoints; i++)
        {
           double u = distribution1.Parameter(i);
           gp_Pnt P = myCurve1->Value(u);
           points.SetValue(i, P);
        }

        for (int j = NbPoints-1; j >= 1; j--)
        {
           double u = distribution2.Parameter(j); 
           gp_Pnt P = myCurve2->Value(u);
           points.SetValue(++i,P); 
        }    

        Geom_BezierCurve BezierCurve(points);
        Handle(Geom_Curve) curve_ptr(&BezierCurve);
        TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
        return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
      }
   }
   return (Curve*) NULL;
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
