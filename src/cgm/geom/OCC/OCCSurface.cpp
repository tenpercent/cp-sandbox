//-------------------------------------------------------------------------
// Filename      : OCCSurface.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Alexander Danilov
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

// ********** BEGIN OCC INCLUDES           **********
#include "OCCSurface.hpp"
#include "OCCQueryEngine.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

// ********** END OCC INCLUDES           **********

// ********** BEGIN CUBIT INCLUDES       **********

#include "CubitSimpleAttrib.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "GeometryQueryEngine.hpp"
#include "DLIList.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "LoopSM.hpp"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "BRepBuilderAPI_MakeShape.hxx"
#include "BRepFilletAPI_MakeFillet2d.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TopExp.hxx"
#include "BRep_Tool.hxx"
#include "LocOpe_SplitShape.hxx"
// ********** END CUBIT INCLUDES           **********


// ********** BEGIN OpenCascade INCLUDES   **********

#include <BRepAdaptor_Surface.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <Extrema_ExtPS.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "BRepClass_FaceClassifier.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
#include "ShapeAnalysis_Surface.hxx"
// ********** END OpenCascade INCLUDES      **********


// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********


// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the TopoDS_Face. 
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCSurface::OCCSurface(TopoDS_Face *theFace)
{
  myTopoDSFace = theFace;
  myShell = NULL;
  myLump = NULL;
  myBody = NULL;
}


OCCSurface::~OCCSurface() 
{
  if(myTopoDSFace)
    delete (TopoDS_Face*)myTopoDSFace;
}

void OCCSurface::set_TopoDS_Face(TopoDS_Face& face)
{
  if(face.IsEqual(*myTopoDSFace))
    return;

  TopoDS_Face* face_ptr = new TopoDS_Face(face);
  if(myTopoDSFace)
    delete (TopoDS_Face*)myTopoDSFace;
  myTopoDSFace = face_ptr ; 
}


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void OCCSurface::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::append_attribute(csa, *myTopoDSFace); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void OCCSurface::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::remove_attribute( csa , *myTopoDSFace); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void OCCSurface::remove_all_simple_attribute_virt()
  { OCCAttribSet::remove_attribute(NULL, *myTopoDSFace); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                                 csa_list)
  { return OCCAttribSet::get_attributes(*myTopoDSFace,csa_list); }

CubitStatus OCCSurface::get_simple_attribute(const CubitString& name,
                                        DLIList<CubitSimpleAttrib*>& csa_list )
  { return OCCAttribSet::get_attributes( name, *myTopoDSFace, csa_list ); }


//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCSurface::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Returns a surface type ID -- the values of these are
//                 determined by OCC.
//
// Special Notes : 
//                 This code is very ACIS-specific and could change with
//                 new versions of ACIS.  There are #defines for the
//                 various surface type ID's.  These are defined in the
//                 header files of each of the specific surface classes
//                 in ACIS.
//
//
// Creator       : Jane HU
//
// Creation Date : 12/06/07
//-------------------------------------------------------------------------
GeometryType OCCSurface::geometry_type()
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  if (asurface.GetType() == GeomAbs_BezierSurface)
     return SPLINE_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_BSplineSurface)
     return BSPLINE_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_Plane)      
     return PLANE_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_Cylinder ||
      asurface.GetType() == GeomAbs_Cone)
     return CONE_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_Sphere)
     return SPHERE_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_Torus)
      return TORUS_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_SurfaceOfRevolution)
     return REVOLUTION_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_SurfaceOfExtrusion)
     return EXTRUSION_SURFACE_TYPE;
  if (asurface.GetType() == GeomAbs_OffsetSurface)
     return OFFSET_SURFACE_TYPE;
  return UNDEFINED_SURFACE_TYPE;  
}
//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCSurface::bounding_box() const 
{
  TopoDS_Face face = *myTopoDSFace;
  BRepAdaptor_Surface asurface(face);
  Bnd_Box aBox;
  BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
  double min[3], max[3];
  aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
  return CubitBox(min, max);
}


CubitStatus OCCSurface::get_point_normal( CubitVector& location,
                                            CubitVector& normal )
{
  return closest_point( bounding_box().center(), &location, &normal );
}   

CubitStatus OCCSurface::closest_point_uv_guess(  
          CubitVector const& location,
          double& , double& ,
          CubitVector* closest_location,
          CubitVector* unit_normal )
{
  // don't use u and v guesses
 return closest_point(location, closest_location, unit_normal);
}


//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface at closest_location and the 
//                 principal curvatures(1-min, 2-max)
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::closest_point( CubitVector const& location, 
                                         CubitVector* closest_location,
                                         CubitVector* unit_normal_ptr,
                                         CubitVector* curvature_1,
                                         CubitVector* curvature_2)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0, u, v;
  int i;
  BRepLProp_SLProps SLP(asurface, 2, Precision::PConfusion());
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
			  SLP.SetParameters(u, v);
		  }
	  }
  
	if (closest_location != NULL)
 	 	*closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  	if (unit_normal_ptr != NULL) {
	  gp_Dir normal;
          //normal of a RefFace point to outside of the material
	  if (SLP.IsNormalDefined()) {
	    normal = SLP.Normal();
            CubitSense sense = get_geometry_sense();
            if(sense == CUBIT_REVERSED)
              normal.Reverse() ;
	      *unit_normal_ptr = CubitVector(normal.X(), normal.Y(), normal.Z()); 
	  }
  	}
  
        gp_Dir MaxD, MinD;
        if (SLP.IsCurvatureDefined())
        {
	   SLP.CurvatureDirections(MaxD, MinD);
           if (curvature_1 != NULL)
              *curvature_1 = CubitVector(MinD.X(), MinD.Y(), MinD.Z());
           if (curvature_2 != NULL)
              *curvature_2 = CubitVector(MaxD.X(), MaxD.Y(), MaxD.Z());
        }
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the trimmed surface to the 
//                 input location. 
//
// Special Notes : 
//-------------------------------------------------------------------------
void OCCSurface::closest_point_trimmed( CubitVector from_point, 
                                         CubitVector& point_on_surface)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(from_point.x(), from_point.y(), from_point.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0;
  int i;
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
		  }
	  }
  }
  point_on_surface = CubitVector(newP.X(), newP.Y(), newP.Z());
}

//-------------------------------------------------------------------------
// Purpose       : This functions computes the point on the surface that is 
//                 closest to the input location and then calculates the 
//                 magnitudes of the principal curvatures at this (possibly, 
//                 new) point on the surface. Specifying the RefVolume for 
//                 reference is optional.
//
// Special Notes :
//
//-------------------------------------------------------------------------

CubitStatus OCCSurface::principal_curvatures(
  CubitVector const& location, 
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0, u, v;
  int i;
  BRepLProp_SLProps SLP(asurface, 2, Precision::PConfusion());
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
			  SLP.SetParameters(u, v);
		  }
	  }
  }
  if (closest_location != NULL)
    *closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());

  if (SLP.IsCurvatureDefined())
  {
    curvature_1 = SLP.MinCurvature();
    curvature_2 = SLP.MaxCurvature();
  }
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Given values of the two parameters, get the position.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitVector OCCSurface::position_from_u_v (double u, double v)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace, Standard_False);
  gp_Pnt p = asurface.Value(u, v);
  return CubitVector (p.X(), p.Y(), p.Z());
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the {u, v} coordinates of the point 
//                 on the Surface closest to the input point (specified in 
//                 global space). The closest_location is also returned.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::u_v_from_position( CubitVector const& location,
                                             double& u,
                                             double& v,
                                             CubitVector* closest_location )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  Handle(Geom_Surface) surface = BRep_Tool::Surface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);

  GeomAPI_ProjectPointOnSurf projection( p, surface);
 
  if(projection.NbPoints() > 0)
  {               
      if(projection.LowerDistance() < Precision::Confusion() )
      {
	  projection.LowerDistanceParameters(u, v);
	  return CUBIT_SUCCESS;
      }
  }
  double minDist=0.0;
  int i;
  Extrema_ExtPS ext(p, asurface, Precision::Confusion(), Precision::Confusion());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
		  }
	  }
  }
  if (closest_location != NULL) *closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  if (ext.IsDone() && (ext.NbExt() > 0))  return CUBIT_SUCCESS;
  Handle(ShapeAnalysis_Surface) su = new ShapeAnalysis_Surface(surface);
  gp_Pnt2d suval = su->ValueOfUV(p, BRep_Tool::Tolerance(*myTopoDSFace));
  suval.Coord(u, v);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the Facet surface is periodic. Not
//                 available yet.
//                 
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic() 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUPeriodic() || asurface.IsVPeriodic())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic_in_U( double& period ) 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  if (!asurface.IsUPeriodic())
     return CUBIT_FALSE;

  period = asurface.UPeriod(); 
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic_in_V( double& period ) 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  if (!asurface.IsVPeriodic())
     return CUBIT_FALSE;

  period = asurface.VPeriod();
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction. Based on comments in SurfaceACIS: "The
//		   assumption is made that the u_param is in the
//		   bounds of the surface. 
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_singular_in_U( double u_param)
{
  //from Acis MasterIndex.htm:
  // singular_u
  // The only singularity recognized is where every value of the 
  // nonconstant parameter generates the same object-space point,
  // and these can only occur at the ends of the parameter range
  // as returned by the functions above. A plane is nonsingular 
  // in both directions. 
  double u_lower, u_upper;
  get_param_range_U( u_lower, u_upper );

  if ( u_param < u_lower - CUBIT_RESABS ||
       u_param > u_upper + CUBIT_RESABS )
  {
    PRINT_ERROR("u parameter is outside parameter bounds.\n");
    return CUBIT_FALSE;
  }

  //Currently, haven't found any singularity check in OCC.
  return CUBIT_FALSE;
}  

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction.  Not available yet.
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_singular_in_V( double v_param)
{
  double v_lower, v_upper;
  get_param_range_V( v_lower, v_upper );

  if ( v_param < v_lower - CUBIT_RESABS ||
       v_param > v_upper + CUBIT_RESABS )
  {
    PRINT_ERROR("v parameter is outside parameter bounds.\n");
    return CUBIT_FALSE;
  }

  //Currently, haven't found any singularity check in OCC.
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the U parameter.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_closed_in_U()
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUClosed())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the V parameter.
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_closed_in_V()
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsVClosed())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Calculates the derivitives at a given parameter location.
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::uv_derivitives( double u,
                                          double v,
                                          CubitVector &du,
                                          CubitVector &dv )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p;
  gp_Vec d1u, d1v;
  asurface.D1(u, v, p, d1u, d1v);
  du = CubitVector(d1u.X(), d1u.Y(), d1u.Z());
  dv = CubitVector(d1v.X(), d1v.Y(), d1v.Z());
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the surface is parametrically defined.
//                 Hopefully later this will be available.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_parametric() 
{
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in U, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::get_param_range_U( double& lower_bound,
                                            double& upper_bound )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  lower_bound = asurface.FirstUParameter();
  upper_bound = asurface.LastUParameter();
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in V, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::get_param_range_V( double& lower_bound,
                                            double& upper_bound )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  lower_bound = asurface.FirstVParameter();

  upper_bound = asurface.LastVParameter();
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Surface
//
//-------------------------------------------------------------------------
double OCCSurface::measure() 
{
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(*myTopoDSFace, myProps);
  return myProps.Mass();
}

//-------------------------------------------------------------------------
// Purpose       : Returns the center of the Surface mass
//
//-------------------------------------------------------------------------
CubitVector OCCSurface::center_point()
{
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(*myTopoDSFace, myProps);
  gp_Pnt pt = myProps.CentreOfMass();
  CubitVector v(pt.X(),pt.Y(), pt.Z());
  return v; 
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying surface.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_position_on( CubitVector &test_position )
{
  CubitVector new_point;
  CubitStatus stat = closest_point(test_position, &new_point, NULL,NULL,NULL);
  if ( !stat )
    return CUBIT_FALSE;
  CubitVector result_vec = test_position - new_point;
  if ( result_vec.length_squared() < GEOMETRY_RESABS )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

CubitPointContainment OCCSurface::point_containment( const CubitVector &point )
{
   TopoDS_Face *face = get_TopoDS_Face();
   gp_Pnt p(point.x(), point.y(), point.z());
   double tol = OCCQueryEngine::instance()->get_sme_resabs_tolerance();

   //It's checking the state of the projected point of THIS Point
   BRepClass_FaceClassifier face_classifier;
   face_classifier.Perform(*face, p, tol);
   TopAbs_State state = face_classifier.State();
   
   //if surface is part of a periodic TopoDS_Face, it'll check the point
   //againt the whole periodic Face, even it outside the occsurface 
   //boundary, if it's on its periodic extension, it'll return as in. 
   if (state == TopAbs_IN)
   {
     //double check if the point is projected on the surface
     CubitVector closest_point;
     this->closest_point_trimmed(point, closest_point);
     if(point.distance_between(closest_point) < tol) 
       return CUBIT_PNT_INSIDE;
     else
       return CUBIT_PNT_OUTSIDE;
   }
   else if (state == TopAbs_OUT)
     return CUBIT_PNT_OUTSIDE;
   else if (state == TopAbs_ON)
     return CUBIT_PNT_BOUNDARY;

   return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment OCCSurface::point_containment( double u_param, 
                                                     double v_param )
{
  CubitVector point = position_from_u_v(u_param, v_param);
  return point_containment(point);
}


CubitSense OCCSurface::get_geometry_sense()
{
  TopAbs_Orientation d = myTopoDSFace->Orientation();
  if (d == TopAbs_FORWARD) return CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) return CUBIT_REVERSED;
  else {
          printf("Check Orientation (surface)");
          return CUBIT_UNKNOWN;
  }
}

void OCCSurface::get_parents_virt( DLIList<TopologyBridge*>& parents )
{ 
  if(myShell) //shell or sheet body
  {
    parents.append(myShell);
    return;
  }

  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCBody * body = NULL;
  DLIList <OCCBody* > *bodies = oqe->BodyList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  bodies->size(); i++)
  {
     body = bodies->get_and_step();
     TopExp::MapShapesAndAncestors(*(body->get_TopoDS_Shape()),
                                   TopAbs_FACE, TopAbs_SHELL, M);
     if(!M.Contains(*(get_TopoDS_Face())))
	continue;

     const TopTools_ListOfShape& ListOfShapes =
                                M.FindFromKey(*(get_TopoDS_Face()));
     if (!ListOfShapes.IsEmpty())
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
           TopoDS_Shell Shell = TopoDS::Shell(it.Value());
           int k = oqe->OCCMap->Find(Shell);
           parents.append((OCCShell*)(oqe->OccToCGM->find(k))->second);
         }
     }
  }
}

void OCCSurface::get_children_virt( DLIList<TopologyBridge*>& children )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSFace, TopAbs_WIRE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
     TopologyBridge *loop = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
     if(loop)
       children.append_unique(loop);
  }
}

// return the sense with respect to the given shell
CubitSense OCCSurface::get_shell_sense( ShellSM* shell_ptr ) const
{
  OCCShell* shell = dynamic_cast<OCCShell*>(shell_ptr);
  if (!shell) // error
    return CUBIT_UNKNOWN;
    
  if (shell->is_sheet())  // relative sense is "both" for sheet
    return CUBIT_UNKNOWN;
  TopAbs_Orientation d = myTopoDSFace->Orientation();
  if (d == TopAbs_FORWARD) return CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) return CUBIT_REVERSED;
  else {
          printf("Check Orientation (surface)");
          return CUBIT_UNKNOWN;
  }
}


int OCCSurface::get_loops( DLIList<OCCLoop*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSFace, TopAbs_WIRE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
     TopologyBridge *loop = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
     result_list.append_unique(dynamic_cast<OCCLoop*>(loop));
  }
  return result_list.size();
}

int OCCSurface::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  DLIList<OCCLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = 0; i < loop_list.size(); i++ )
  {
    result_list += loop_list.next(i)->coedges( );
  }
  return result_list.size();
}

int OCCSurface::get_curves( DLIList<OCCCurve*>& result_list )
{
  DLIList<OCCCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    OCCCoEdge* coedge = coedge_list.get_and_step();
    OCCCurve* curve = dynamic_cast<OCCCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
  return result_list.size();
}

//----------------------------------------------------------------
// Function: to update the core Surface
//           for any movement  or Boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCSurface::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                         BRepAlgoAPI_BooleanOperation *op)
{
  assert(aBRepTrsf != NULL || op != NULL);

  TopoDS_Shape shape;
  if (aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Face());
  else
  {
    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(*get_TopoDS_Face()));
    if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(*get_TopoDS_Face()));
    if (shapes.Extent() == 1)
      shape = shapes.First();
    else if(shapes.Extent() > 1)
    {
      //update all attributes first.
      TopTools_ListIteratorOfListOfShape it;
      it.Initialize(shapes);
      for(; it.More(); it.Next())
      {
        shape = it.Value();
        OCCQueryEngine::instance()->copy_attributes(*get_TopoDS_Face(), shape);
      }
      shape = shapes.First();
    }
    else if(op->IsDeleted(*get_TopoDS_Face()))
      ;
    else
      return CUBIT_SUCCESS;
  }
 
  TopoDS_Face surface; 
  if(!shape.IsNull())
    surface = TopoDS::Face(shape);

  if (aBRepTrsf) 
  {
    //set the loops
    DLIList<OCCLoop *> loops;
    this->get_loops(loops);
    for (int i = 1; i <= loops.size(); i++)
    {
       OCCLoop *loop = loops.get_and_step();
       loop->update_OCC_entity(aBRepTrsf, op);
    }
    OCCQueryEngine::instance()->update_OCC_map(*myTopoDSFace, surface);
  }

  else if(op)
    update_OCC_entity(*myTopoDSFace, surface, op);

  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: TopoDS_Shape level function to update the core Surface
//           for any movement  or Boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCSurface::update_OCC_entity(TopoDS_Face& old_surface,
                                          TopoDS_Shape& new_surface,
                                          BRepBuilderAPI_MakeShape *op,
                                          TopoDS_Vertex* removed_vertex,
                                          LocOpe_SplitShape* sp) 
{
  //set the Wires
  TopTools_IndexedMapOfShape M, M2;
  TopoDS_Shape shape, shape2, shape_edge, shape_vertex;
  TopExp::MapShapes(old_surface, TopAbs_WIRE, M);

  TopTools_ListOfShape shapes;  
  BRepFilletAPI_MakeFillet2d* test_op = NULL;

  for (int ii=1; ii<=M.Extent(); ii++) 
  {
     TopoDS_Wire wire = TopoDS::Wire(M(ii));
     TopTools_ListOfShape shapes;
     if(op)
     {
       test_op = dynamic_cast<BRepFilletAPI_MakeFillet2d*>(op);
       if(!test_op)
         shapes.Assign(op->Modified(wire));
       if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(wire));
       if(!new_surface.IsNull())
         TopExp::MapShapes(new_surface,TopAbs_WIRE, M2);
     }
     else if(sp)
       shapes.Assign(sp->DescendantShapes(wire));

     if (shapes.Extent() == 1)
     {
       shape = shapes.First();
       if(M2.Extent() == 1)
       {
         shape2 = TopoDS::Wire(M2(1));
         if(!shape.IsSame(shape2))
           shape = shape2;
       }
       else if(M2.Extent() > 1)
         shape.Nullify();
     }
     else if(shapes.Extent() > 1)
       shape.Nullify();
     else if(op->IsDeleted(wire) || shapes.Extent() == 0)
     {
       TopTools_IndexedMapOfShape M_new;
       TopExp::MapShapes(new_surface, TopAbs_WIRE, M_new);
       if (M_new.Extent()== 1)
         shape = M_new(1);
       else
         shape.Nullify();
     }
     else
     {
       shape = wire;
       continue;
     }

     //set curves
     BRepTools_WireExplorer Ex;
      
     for(Ex.Init(wire); Ex.More();Ex.Next())
     {
       TopoDS_Edge edge = Ex.Current();
       if(op && !test_op)
       {
         shapes.Assign(op->Modified(edge));
         if(shapes.Extent() == 0)
           shapes.Assign(op->Generated(edge));
       }
         
       else if(sp)
         shapes.Assign(sp->DescendantShapes(edge));

       if (shapes.Extent() == 1)
       {
        //in fillet creating mothod, one edge could generated a face, so check
        //it here.
         TopAbs_ShapeEnum type = shapes.First().TShape()->ShapeType(); 
         if(type != TopAbs_EDGE)
           shape_edge.Nullify();
         else
           shape_edge = shapes.First();
       }
       else if (shapes.Extent() > 1)
       {
         //update all attributes first.
         TopTools_ListIteratorOfListOfShape it;
         it.Initialize(shapes);
         for(; it.More(); it.Next())
         {
           shape_edge = it.Value();
           OCCQueryEngine::instance()->copy_attributes(edge, shape_edge);
         }
         shape_edge.Nullify();
       }
       else if (op->IsDeleted(edge))
         shape_edge.Nullify(); 
       else if (test_op)
       {
         if(!test_op->IsModified(edge))
           shape_edge = edge;
         else
           shape_edge = (test_op->Modified(edge)).First();
       } 
       else
         shape_edge = edge;

       //update vertex
       TopoDS_Vertex vertex = Ex.CurrentVertex();
       shapes.Clear();
       if(test_op)
         assert(removed_vertex != NULL);

       if(op && ! test_op )
         shapes.Assign(op->Modified(vertex));
       else if(sp)
         shapes.Assign(sp->DescendantShapes(vertex));

       if (shapes.Extent() == 1)
         shape_vertex = shapes.First();

       else if(shapes.Extent() > 1)
       {
         //update all attributes first.
         TopTools_ListIteratorOfListOfShape it;
         it.Initialize(shapes);
         for(; it.More(); it.Next())
         {
           shape_vertex = it.Value();
           OCCQueryEngine::instance()->copy_attributes(vertex, shape_vertex);
         }
         shape_vertex.Nullify() ;
       }
       else if(op->IsDeleted(vertex) || (test_op && vertex.IsSame( *removed_vertex)))
         shape_vertex.Nullify() ;
         
       else
         shape_vertex = vertex;
      
       if(!vertex.IsSame(shape_vertex) )
         OCCQueryEngine::instance()->update_OCC_map(vertex, shape_vertex);

       if (!edge.IsSame(shape_edge))
         OCCQueryEngine::instance()->update_OCC_map(edge, shape_edge);
     }
     if (!wire.IsSame(shape))
       OCCQueryEngine::instance()->update_OCC_map(wire, shape);
  }
  double dTOL = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
  if (!old_surface.IsSame(new_surface))
  {
    TopAbs_ShapeEnum shapetype =  TopAbs_SHAPE;
    if(!new_surface.IsNull())
      shapetype = new_surface.TShape()->ShapeType();  
    if(shapetype == TopAbs_FACE || new_surface.IsNull())
      OCCQueryEngine::instance()->update_OCC_map(old_surface, new_surface);
    else 
    {
      TopTools_IndexedMapOfShape M;
      TopExp::MapShapes(new_surface, TopAbs_FACE, M);   
      TopoDS_Shape new_shape;
      if(M.Extent() == 1)
        new_shape = M(1);
      else if(M.Extent() > 1)
      {
        for(int i = 1; i <= M.Extent(); i++)
        {
          GProp_GProps myProps;
          BRepGProp::SurfaceProperties(old_surface, myProps);
          double orig_mass = myProps.Mass();
          gp_Pnt orig_pnt = myProps.CentreOfMass();
          BRepGProp::SurfaceProperties(M(i), myProps);
          double after_mass = myProps.Mass();
          gp_Pnt after_pnt = myProps.CentreOfMass();
          if(fabs(-after_mass + orig_mass) <= dTOL && 
             orig_pnt.IsEqual(after_pnt, dTOL))
          {
            new_shape = M(i);
            break;
          }
        }
      }
      OCCQueryEngine::instance()->update_OCC_map(old_surface, new_shape);
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus OCCSurface::get_bodies(DLIList<OCCBody*>& bodies)
{
   TopoDS_Face* topo_face = this->get_TopoDS_Face();
   OCCQueryEngine* oqe = OCCQueryEngine::instance();
   DLIList <OCCBody* > *all_bodies = oqe->BodyList;
   TopTools_IndexedDataMapOfShapeListOfShape M;
   OCCBody * body = NULL;
   for(int j = 0; j <  all_bodies->size(); j++)
   {
     body = all_bodies->get_and_step();
     TopExp_Explorer Ex;
     TopoDS_Face the_face;
     TopoDS_Shape ashape = *(body->get_TopoDS_Shape());
     M.Clear();
     TopExp::MapShapesAndAncestors(ashape, TopAbs_FACE, TopAbs_COMPSOLID, M);
     if(!M.Contains(*topo_face))
       continue;
     bodies.append_unique(body);
  }
  return CUBIT_SUCCESS;
}
// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********


// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
