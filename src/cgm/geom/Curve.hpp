//-------------------------------------------------------------------------
// Filename      : Curve.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef CURVE_HPP
#define CURVE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "RefEdge.hpp"

#include "GeometryEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
template <class X> class DLIList;
class FacetEvalTool;
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN ENUM DEFINITIONS       **********
// ********** END ENUM DEFINITIONS         **********

class CUBIT_GEOM_EXPORT Curve : public GeometryEntity
{
public :
  
  Curve();
    //- The default constructor
  
  virtual ~Curve();
    //- The destructor  

  virtual const type_info& topology_entity_type_info() const;

  virtual CubitSense relative_sense(Curve *other_curve);
    //- given another curve, return whether the curves are "aligned" (CUBIT_FORWARD), 
    //- oppositely aligned (CUBIT_REVERSED), or inconclusive (CUBIT_UNKNOWN)
  
  virtual CubitStatus get_point_direction( CubitVector& origin, 
                                           CubitVector& direction );
    //- Only valid for straight lines
    //- Finds the underlying line's origin and direction unit vector
    //- Returns CUBIT_FAILURE if curve is not a line
  
  virtual CubitStatus get_center_radius( CubitVector& center, double& radius ) = 0;
    //- Only valid for arcs
    //- Finds the underlying arc's center point and radius
    //- Returns CUBIT_FAILURE if curve is not an arc
  
  virtual double length_from_u( double parameter1,
                                double parameter2 ) = 0;
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
  
  virtual double get_arc_length();
    //R double
    //R- Returned total length value
    //- This function returns the arc length along the entire Curve
  
  virtual double get_arc_length( const CubitVector &point1,
                                 const CubitVector &point2 );
    //R double
    //R- Returned length value
    //I point1
    //I- The first position value
    //I point2
    //I- The second position value
    //- This function returns the arc length along the Curve starting from
    //- the point1 going to point2
  
  virtual double get_arc_length( const CubitVector &point1,
                                 const int which_end );
    //R double
    //R- Returned length value
    //I point1
    //I- The first position value
    //I which_end
    //I- 0 for start, 1 for end
    //- This function returns the arc length along the Curve starting from
    //- point1 going to either the start or the end
  
  virtual CubitVector center_point();
    //R CubitVector
    //R- center point on this edge
    //- This function returns the center point on this edge, by
    //- arc length
  
  virtual CubitStatus mid_point(const CubitVector &point1,
                                const CubitVector &point2,
                                CubitVector& my_mid_point );
    //R CubitStatus
    //I CubitVector, CubitVector
    //I- points between which the mid_point is needed
    //O CubitVector
    //O- mid point on this edge
    //- This function returns the mid point between two points on this
    //-edge, by parameter
  
  virtual CubitStatus mid_point(CubitVector& my_mid_point);
    //R CubitStatus
    //O CubitVector
    //O- mid point on this edge
    //- This function returns the mid point on this edge, by parameter
  
  virtual CubitStatus position_from_fraction( const double fraction_along_curve,
                                              CubitVector &output_position );
    //R CubitStatus
    //I fraction in parameter space along curve (1/3,2/3,4/5...)
    //O- position where percent in parameter space lies.
    //-This function takes the given fraction, finds the parameter value,
    //-and calculates this position.
  
  virtual CubitBoolean is_periodic( double& period) = 0;
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
  
  virtual CubitBoolean get_param_range( double& lower_bound,
                                        double& upper_bound ) = 0;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range.
    //O upper_bound
    //O- The upper bound of the parametric range.
    //- Returns the lower and upper parametric bounds of the 
    //- Curve, based on the extent of the its first underlying entity.
    //- The boolean return value indicates whether the Curve is 
    //- parametrically defined or not -- returns CUBIT_TRUE if yes and
    //- CUBIT_FALSE, otherwise.  If the Curve is not parametrically
    //- defined, the values of the upper and lower bounds are
    //- undetermined.
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
                                           CubitSense& return_sense) = 0;
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
  
  virtual CubitStatus closest_point( CubitVector const& location, 
                                     CubitVector& closest_location,
                                     CubitVector* tangent_ptr = NULL,
                                     CubitVector* curvature_ptr = NULL,
                                     double *param = NULL) = 0;
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
    //- underlying solid model entities, if any.
  
  virtual CubitStatus closest_point_trimmed( CubitVector const& from_pt,
                                             CubitVector& result_pt );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I from_pt
    //I- The position to evaluate wrt curve.
    //I result_pt
    //I- The resulting point on the curve.
    //- This method finds the closest point to a BOUNDED curve.
    //
    // Added by Jason Kraftcheck, 07/17/98.
    // See also: CurveACIS::closest_point_trimmed(..).

  virtual CubitPointContainment point_containment( const CubitVector &point ) = 0;
    //R CubitPointContainment - is the point on bounds of the curve?
    //R- CUBIT_PNT_OFF, CUBIT_PNT_ON, CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Curve
    // NOTE: POINT MUST LIE ON THE CURVE FOR THIS FUNCTION TO WORK PROPERLY.
  
  virtual double u_from_position (const CubitVector& input_position) = 0;
    //R double
    //R- The returned "u" parameter value in local parametric space
    //I input_position
    //I- The input position for which "u" is to be computed.
    //- This function returns the coordinate of a point in the local
    //- parametric (u) space that corresponds to the input position in
    //- global (world) space.  The input point is first moved to the
    //- closest point on the Curve and the parameter value of that
    //- point is determined.
  
  virtual CubitStatus position_from_u (double u_value,
                                       CubitVector& output_position) = 0;
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
  
  virtual double u_from_arc_length ( double root_param,
                                     double arc_length ) = 0;
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
    //- A negative value for distance would force the search to go in the 
    //- negative (sense) direction of the RefEdge.
    //-
    //- NOTE:
    //- The important assumption that is made in this routine is that
    //- the end points of the RefEdge that owns this CurveACIS are the same
    //- as the end points of the first solid model entity in the list of 
    //- solid model entities associated with this Curve.
  
  virtual CubitBoolean is_position_on( const CubitVector &test_position ) = 0;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the curve.

  virtual GeometryType geometry_type()
  {return UNDEFINED_CURVE_TYPE;};
  //R GeometryType (enum)
  //R- The enumerated type of the geometric representation

  virtual double start_param() = 0;
    //R double parameter
    //R- start parameter of curve with respect to refEdge.
  
  virtual double end_param() = 0;
    //R double parameter
    //R- start parameter of curve with respect to refEdge.

  virtual CubitStatus point_from_arc_length( double root_param,
                                             const double arc_length,
                                             CubitVector &new_point);

    //R CubitStatus
    //R- CUBIT_SUCCESS, CUBIT_FAILURE
    //I- root_param-starting parameter for point on curve from which arc_length extends.
    //O- new_point, from starting parameter by arc_lenth.
  
  virtual CubitStatus point_from_arc_length( const CubitVector &root_point,
                                             double const arc_length,
                                             CubitVector& new_point );
    //R void
    //R- pointer to void
    //I- root_point-point from which arc_length extends.
    //O- new_point, from root_point by arc_lenth.
  
//=========  Add Code by SRS of Cat,  3/3/99 2:34:50 PM  =========
  virtual CubitBoolean is_tolerant(){ return CUBIT_FALSE; };
    //- This function is overloaded in the CurveACIS class only.
    //- Tolerant edges can get created by the ACIS healer if the
    //- edge cannot be healed.
//=========  Code End by SRS of Cat,  3/3/99 2:34:50 PM  =========

	virtual CubitBoolean G1_discontinuous( double u,
	                                       CubitVector* minus_tangent = NULL,
	                                       CubitVector* plus_tangent = NULL );
		//R CubitBoolean
		//R- CUBIT_TRUE/CUBIT_FALSE
		//I u
		//I- The parameter value on the curve to test at.
		//O minus_tanget
		//O- The tangent direction on the decreasing-parameter-value 
		//O- side of u.  This vector, if passed, will only be changed
		//O- when the return value is CUBIT_TRUE.
		//O plus_tanget
		//O- The tangent direction on the increasing-parameter-value
		//O- side of u.  This vector, if passed, will only be changed
		//O- when the return value is CUBIT_TRUE.
		//- Check for a G1 discontinuity in a curve.  Derived classes should
    //- override the function with a more exact implementation.  The
    //- default implementation provided in Curve is a simple numeric
    //- approximation. 

protected: 
  
private:
};


// ********** BEGIN INLINE FUNCTIONS       **********
inline
double Curve::get_arc_length()
{
  return measure();
}

inline
double Curve::get_arc_length( const CubitVector &point1,
                              const CubitVector &point2 )
{
  double param1 = u_from_position(point1);
  double param2 = u_from_position(point2);
  return length_from_u(param1, param2);
}

inline
double Curve::get_arc_length( const CubitVector &point1,
                              const int which_end )
{
  double param1 = u_from_position(point1);
  double param2;
  if (which_end == 0) param2 = start_param();
  else param2 = end_param();
  
  return length_from_u(param1, param2);
}

inline
CubitVector Curve::center_point()
{
  double param1 = start_param();
  double length = 0.5 * measure();
  double param2 = u_from_arc_length(param1, length);
  CubitVector center;
  position_from_u(param2, center);
  return center;
}

inline
CubitStatus Curve::mid_point(CubitVector &my_mid_point)
{
  double param1 = 0.5 * (start_param() + end_param());
  
  return position_from_u(param1, my_mid_point);
}

inline
CubitStatus Curve::position_from_fraction( const double fraction_along_curve,
                                           CubitVector &output_position )
{
  CubitStatus result = CUBIT_FAILURE;
  double param_1 = start_param() +
    fraction_along_curve * (end_param() - start_param() );
  result = position_from_u( param_1, output_position );
  return result;
}

inline
CubitStatus Curve::mid_point(const CubitVector &point1,
                             const CubitVector &point2,
                             CubitVector& my_mid_point )
{
  double param1 = u_from_position(point1);
  double param2 = u_from_position(point2);
  param1 = 0.5 * (param1 + param2);
  return position_from_u(param1, my_mid_point);
}

  //R CubitStatus
  //I CubitVector, CubitVector
  //I- points between which the mid_point is needed
  //O CubitVector
  //O- mid point on this edge
  //- This function returns the mid point between two points on this
  //-edge, by parameter



// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

