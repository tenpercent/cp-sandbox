//-------------------------------------------------------------------------
// Filename      : OCCSurface.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Jane Hu  
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

#ifndef SURFACE_OCC_HPP
#define SURFACE_OCC_HPP


// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN OCC INCLUDES          **********

#include "OCCAttribSet.hpp"
#include "TopoDS_Face.hxx"

// ********** END OCC INCLUDES          **********


// ********** BEGIN CUBIT INCLUDES          **********

#include "CubitDefines.h"
#include "Surface.hpp"

// ********** END CUBIT INCLUDES          **********

class TopologyEntity;
class RefVolume;
class RefFace;
class RefVolume;
class OCCShell;
class OCCAttrib;

class OCCBody;
class OCCLump;
class OCCLoop;
class OCCCoEdge;
class OCCCurve;
class OCCPoint;
class BRepBuilderAPI_Transform;
class BRepAlgoAPI_BooleanOperation;
class BRepBuilderAPI_MakeShape;
class LocOpe_SplitShape;
class TopoDS_Vertex;
//// class CubitTransformMatrix;

class OCCSurface : public Surface
{

public :
  
  OCCSurface(TopoDS_Face *theFace);

  virtual ~OCCSurface() ;
    //- The destructor
   
  //if op is of BRepFilletAPI_MakeFillet2d type, need to input the 
  //removed_vertex for update.
  static CubitStatus update_OCC_entity(TopoDS_Face& old_surface,
                                       TopoDS_Shape& new_surface,
                                       BRepBuilderAPI_MakeShape *op,
                                       TopoDS_Vertex* removed_vertex = NULL,
                                       LocOpe_SplitShape* sp = NULL);


  void add_shell(OCCShell* shell) { myShellList.append_unique(shell);}
  DLIList<OCCShell*> shells() {return myShellList;}
  void remove_shell(OCCShell* shell) {myShellList.remove(shell);}
  void clean_shells(){myShellList.clean_out();}

  void add_hardpoint(OCCPoint* HardPoint){myHardPoints.append(HardPoint);}
  void remove_hardpoint(OCCPoint* HardPoint){myHardPoints.remove(HardPoint);}

  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
  
  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  // Added by CAT
  virtual CubitStatus get_point_normal( CubitVector& ,
                                        CubitVector& );
    //- Only valid for planar surfaces
    //- Finds the underlying plane's origin and normal vector
    //- Returns CubitFailure if not a plane.  The origin and normal 
    //- are returned directly from the underlying format for a plane.
  
  virtual void closest_point_trimmed(CubitVector from_point,
                                     CubitVector& point_on_surface);
    //R void
    //I CubitVector
    //I- point from which to find closest point on trimmed surface
    //O CubitVector
    //O- point on trimmed surface closest to passed-in point
    //- This function finds the closest point on a TRIMMED surface to the
    //- passed-in point.
  
  virtual CubitStatus closest_point_uv_guess(  
      CubitVector const& location,
      double &u, double &v,
      CubitVector* closest_location = NULL,
      CubitVector* unit_normal = NULL );

  virtual CubitSense get_shell_sense(ShellSM*) const;

  virtual CubitStatus closest_point(  
    CubitVector const& location, 
    CubitVector* closest_location = NULL,
    CubitVector* unit_normal_ptr = NULL,
    CubitVector* curvature1_ptr = NULL,
    CubitVector* curvature2_ptr = NULL);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point to which the closest point on the surface is desired.
    //O closest_location
    //O- The point on the Surface, closest to the 
    //O- input location (which might not be on the Surface).  This is
    //O- input as a reference so that the function can modify its
    //O- contents.
    //O unit_normal_ptr
    //O- The normal (represented as a unit vector) at the closest_location.
    //O- If this pointer is NULL, the normal is not returned.
    //O curvature1_ptr
    //O- The first principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //O curvature2_ptr
    //O- The second principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //- This function computes the point on the surface closest to the input 
    //- location -- i.e., closest_location. 
    //- The first Facet FACE in the list
    //- is queried.
    //-
    //- If the input pointer values of unit_normal, curvature1 and
    //- curvature2
    //- are non-NULL, the normal and principal curvatures, too, are
    //- returned.  These are computed at closest_location, not at the
    //- input location.
    //-
    //- NOTE:
    //- It is assumed that if the calling code needs the normal or the 
    //- principal curvatures, it will *allocate* space for the CubitVectors
    //- before sending in the pointers.
  
  virtual CubitStatus principal_curvatures(
    CubitVector const& location, 
    double& curvature_1,
    double& curvature_2,
    CubitVector* closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point at which the curvatures are being requested -- it is also
    //I- the point to which the closest point on the surface is returned.
    //I- curvatures.
    //O closest_location
    //O- The point on the surface, closest to the input location (this
    //O- might not be on the surface).  This is input as a reference 
    //O- so that the function can modify its contents.
    //O curvature_1/2
    //O- Returned principal curvature magnitudes.
    //- This functions computes the point on the surface that is closest
    //- to the input location and then calculates the magnitudes of the
    //- principal curvatures at this (possibly, new) point on the surface. 
  
  virtual CubitVector position_from_u_v (double u, double v);
    //R CubitVector
    //R- Returned position vector.
    //I u, v
    //I- Input point in {u.v} space
    //- This function returns the coordinates in world space of a point
    //- in the parameter space of this Surface object.
  
  virtual CubitStatus u_v_from_position (CubitVector const& location,
                                         double& u, 
                                         double& v,
                                         CubitVector*
                                         closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The input point in global space
    //O closest_point
    //O- The point on the Surface closest to the input location
    //O u, v
    //O- The returned u, v coordinate values (in local parametric space)
    //O- of the closest_point
    //I refvolume_ptr
    //- This function returns the {u, v} coordinates of the point 
    //- on the Surface closest to the input point (specified in global
    //- space). The closest_location is also returned.
  
  virtual CubitBoolean is_periodic();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- OCCSurface is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
    virtual CubitBoolean is_periodic_in_U( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the Facet surface object associated
    //- with one of the FACEs of this OCCSurface object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_periodic_in_V( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the Facet surface object associated
    //- with one of the FACEs of this OCCSurface object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_singular_in_U( double u_param );
  virtual CubitBoolean is_singular_in_V( double v_param );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I double u/v parameter value.
    //- Determines if the surface is singular in a given direction
    //- at a given parameter value.
  
  virtual CubitBoolean is_closed_in_U();  
  virtual CubitBoolean is_closed_in_V();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Determines if the surface is closed, smoothly or not in the
    //- given parameter direction.
    //- A periodic surface is always closed but a closed surface is
    //- is not always periodic.
  
  virtual CubitStatus uv_derivitives( double u_param,
                                      double v_param,
                                      CubitVector &du,
                                      CubitVector &dv );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O- du, dv
    //- Determines the u and v derivitives from the given parameter
    //- values.
  
  TopoDS_Face *get_TopoDS_Face(){return myTopoDSFace;}
  void set_TopoDS_Face(TopoDS_Face& face); 

  int get_loops(DLIList<OCCLoop*>&);

  int get_coedges(DLIList<OCCCoEdge*>&);

  int get_curves(DLIList<OCCCurve*>&); 

  void set_shell(OCCShell* shell)
  { myShell = shell;}

  void set_lump(OCCLump* lump)
  { myLump = lump;}

  void set_body(OCCBody* body)
  { myBody = body;}

  OCCShell* my_shell() {return myShell;}

  OCCLump* my_lump() {return myLump;}

  OCCBody* my_body() {return myBody;}

  CubitStatus get_bodies(DLIList<OCCBody*> &bodies);

  virtual CubitBoolean is_parametric();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
 			//- This method returns CUBIT_TRUE if a parametric representation
			//- is available for the surface
  
  virtual CubitBoolean get_param_range_U( double& lower_bound,
                                          double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in U, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean get_param_range_V( double& lower_bound,
                                          double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in V, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_position_on( CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the surface.

  virtual CubitPointContainment point_containment( const CubitVector &point );
  virtual CubitPointContainment point_containment( double u, double v );

  GeometryType geometry_type();
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point

  CubitVector center_point();
  virtual CubitSense get_geometry_sense();
    //- Return the relative surface sense. (see below)
  
  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );

  CubitBoolean is_flat();     //// Not in SurfaceACIS
  CubitBoolean is_spherical(); //// Not in SurfaceACIS
  CubitBoolean is_conical();  //// Not in SurfaceACIS

  CubitStatus update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                 BRepAlgoAPI_BooleanOperation *op = NULL);
protected: 

private:
    //CubitSense sense_;
    //- The sense of the RefFace that owns this Surface with respect
    //- to the positive sense of the first FACE in FACEPtrList_.
    //- When a Surface is first constructed, this value is arbitrarily
    //- set to CUBIT_FORWARD.
    //- In the case of Surfaces, the normal is used in determining
    //- the relative sense value.
    //- MJP NOTE:
    //- Not only does the RefFace have a sense wrt its Surface, but each
    //- Facet FACE has a sense wrt its underlying "surface" object.

  OCCAttribSet attribSet;
    //List of OCCAttrib*'s instead of CubitSimpleAttribs 

  TopoDS_Face *myTopoDSFace;

  DLIList<OCCShell*> myShellList;
  //Following 3 members are only for sheeted body.
  OCCShell* myShell;
  OCCLump* myLump;
  OCCBody* myBody;
  DLIList<OCCPoint*> myHardPoints;

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

