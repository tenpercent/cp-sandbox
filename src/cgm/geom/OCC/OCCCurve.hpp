//-------------------------------------------------------------------------
// Filename      : OCCCurve.hpp
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

#ifndef CURVE_OCC_HPP
#define CURVE_OCC_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********
// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Curve.hpp"
#include "OCCAttribSet.hpp"
#include "TopoDS_Edge.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class OCCAttrib;
class Point;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCPoint;
class BRepBuilderAPI_Transform; 
class BRepAlgoAPI_BooleanOperation;
// ********** END FORWARD DECLARATIONS     **********

class OCCCurve : public Curve
{
public :
  
  OCCCurve( TopoDS_Edge* theEdge );
  
  virtual ~OCCCurve() ;
    //- The destructor

  void set_myMarked( CubitBoolean marked) {myMarked = marked;}

  void add_loop(OCCLoop* loop) { myLoopList.append_unique(loop);}   
  DLIList<OCCLoop*> loops() {return myLoopList;}
  void remove_loop(OCCLoop* loop) {myLoopList.remove(loop);}
  void clean_loops(){myLoopList.clean_out();}

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
    //- see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
    //-
    //- If there is an error computing the length a value of -1.0 is
    //- returned.
  
  virtual double length_from_u( double parameter1,
                                double parameter2 );
    //R double
    //R- Returned length value
    //I parameter1
    //I- The first parameter value
    //I parameter2
    //I- The second parameter value
    //- This function returns the arc length along the Curve starting from
    //- the point represented by the parameter1 going to the point represented
    //- by parameter2.
    //-
    //- The sign of the returned length value is always positive.
  
  virtual CubitBoolean is_periodic( double& period);
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- Returned period value
    //- This function determines whether the underlying geometry of the
    //- Curve is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //-
    //- If it is periodic, then it returns the period in the input
    //- reference variable, "period". This value is set to 0.0 if
    //- the Curve is not periodic.
    //- 
  
  virtual CubitBoolean get_param_range( double& lower_bound,
                                        double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range.
    //O upper_bound
    //O- The upper bound of the parametric range.
    //- Returns the lower and upper parametric bounds of the 
    //- Curve.
    //-
    //- IMPORTANT NOTE:
    //- Note that the lower bound is the parameter value of the start
    //- location of the RefEdge that uses this Curve and the upper
    //- bound is the parameter value of the end location of the RefEdge
    //- that uses this Curve.  This takes into account the sense of the
    //- RefEdge with respect to the Curve (which could be REVERSED).
    //- Hence, the numerical value of the lower parameter bound could be
    //- greater than that of the upper parameter bound.
  
  virtual CubitStatus get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                           CubitSense& return_sense);
    //- Finds the extrema along this RefEdge.  An extremum is defined as
    //- a local min or max in the direction of one of the primary axial directions.
    //- O-interior_points: list of coordinates where the extrema occur.
    //- O-return_sense: Whether the interior extrema are ordered in the
    //-                 FORWARD or REVERSED direction of this RefEdge.
    //-
    //- ***IMPORTANT!!!***
    //-    This function dynamically allocates the CubitVectors appended to
    //-    interior_points.  It is the responsibility of the calling code to
    //-    delete these CubitVectors (or in the case of RefEdge, to make sure
    //-    that *it's* calling code knows that it should delete the CubitVectors)!
  
  CubitStatus get_interior_extrema_in_direction(DLIList<CubitVector*>&,
 						CubitVector);

  virtual CubitStatus closest_point( CubitVector const& location, 
                                     CubitVector& closest_location,
                                     CubitVector* tangent_ptr = NULL,
                                     CubitVector* curvature_ptr = NULL,
                                     double *param = NULL);
    //R void
    //I location
    //I- The point to which the closest point on the Curve is desired.
    //O closest_location
    //O- The point on the Curve, closest to the input location which
    //O- might not be on the Curve.  This is input as a reference 
    //O- so that the function can modify its contents.
    //O tangent_ptr
    //O- The tangent to the Curve (output as a unit vector) at the 
    //O- closest_location.
    //O curvature_ptr
    //O- The curvature of the Curve at the closest_location.
    //- This function computes the point on the Curve closest to the input 
    //- location.
    //-
    //- If the tangent and/or curvature is required, then the calling code
    //- is responsible for allocating space for the CubitVector(s) and
    //- sending in the relevant non-NULL pointers.  If either of these
    //- pointers is NULL, the related quantity is not computed.
    //-
    //- Notes:
    //- The tangent direction is always in the positive direction of the 
    //- *owning RefEdge*, regardless of the positive direction of the
    //- underlying solid model entities.
  
  void get_tangent( CubitVector const& location, 
                    CubitVector& tangent);
    //- this function returns the tangent vector at the given location
  
  void get_curvature( CubitVector const& location, 
                      CubitVector& curvature);
    //- this function returns the curvature vector at the given location
  
  virtual CubitStatus position_from_u (double u_value,
                                       CubitVector& output_position);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I u_value
    //I- The input u parameter value
    //O output_position
    //O- The output position
    //- This function returns the coordinates of a point in the global
    //- (world) space that corresponds to the input parametric position 
    //- in the local space.
    //-
    //- If the input parameter value is not defined for the Curve, then 
    //- the input CubitVector is not modified and CUBIT_FAILURE is
    //- returned. Otherwise, position is appropriately modified and
    //- CUBIT_SUCCESS is returned.
    //-
    //- If the curve is periodic, the input u_value is first "normalized"
    //- to the fundamental period of the Curve before its position
    //- in global space is determined.
  
  double u_from_position (const CubitVector& input_position);
  
  virtual double u_from_arc_length ( double root_param,
                                     double arc_length );
    //R double
    //R- Returned parameter value
    //I root_param
    //I- The parameter value of the "root point"
    //I arc_length
    //I- A distance along the Curve
    //- This function returns the parameter value of the point that is
    //- "arc_length" away from the root point in the
    //- positive sense direction of the owning RefEdge.
    //-
    //- If arc_length is negative, the new point (whose parameter value is
    //- being computed) is in the negative sense direction (along the 
    //- RefEdge) from the root point (whose parameter value is root_param).
  
  virtual CubitBoolean is_position_on( const CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the curve.
  
  GeometryType geometry_type();
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
  
  CubitStatus get_point_direction( CubitVector& origin, CubitVector& direction );
  //- Only valid for straight lines
  //- Finds the underlying line's origin and direction unit vector
  //- Returns CUBIT_FAILURE if curve is not a line

  CubitStatus get_center_radius( CubitVector& center, double& radius );
  //- Only valid for arcs
  //- Finds the underlying arc's center point and radius
  //- Returns CUBIT_FAILURE if curve is not an arc

  virtual double start_param();
    //R double parameter
    //R- start parameter of curve with respect to refEdge.
  
  virtual double end_param();
    //R double parameter
    //R- start parameter of curve with respect to refEdge.
    
  virtual CubitBoolean G1_discontinuous( double param,
                                         CubitVector* minus_tangent = NULL,
                                         CubitVector* plus_tangent = NULL );
  
  virtual CubitPointContainment point_containment( const CubitVector &point );
    //R CubitPointContainment - is the point outside, inside or on the boundary?
    //R- CUBIT_PNT_OUTSIDE, CUBIT_PNT_INSIDE, CUBIT_PNT_BOUNDARY, 
    //   CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Surface
    //I double
    //I- u coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    //I double
    //I- v coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    // NOTE: POINT MUST LIE ON THE SURFACE FOR THIS FUNCTION TO WORK PROPERLY.

  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );

  void get_points(DLIList<OCCPoint*>& point_list);
    //- Gets the list of points describing this curve.

  TopoDS_Edge *get_TopoDS_Edge( )
    { return myTopoDSEdge; } 
  void set_TopoDS_Edge(TopoDS_Edge edge);

  void update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                          BRepAlgoAPI_BooleanOperation *op = NULL );
 
  Curve* project_curve(Surface* face_ptr,
                       DLIList<Point*>&  normal_proj_points,
                       CubitBoolean closed,
                       const CubitVector* third_point);
protected: 
  
private:
  
  void adjust_periodic_parameter(double& param);
  
  OCCAttribSet attribSet;
    //List of OCCAttrib*'s instead of CubitSimpleAttribs 
  
  TopoDS_Edge *myTopoDSEdge;
  DLIList<OCCLoop*> myLoopList;
  bool periodic;
  CubitBoolean myMarked ;
};

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

