//-------------------------------------------------------------------------
// Filename      : Curve.cc
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

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "Curve.hpp"
#include "Point.hpp"
#include "GeometryQueryTool.hpp"

#include "CastTo.hpp"
#include "CubitMessage.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------

Curve::Curve()
{}

//-------------------------------------------------------------------------
// Purpose       : The destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------

Curve::~Curve()
{}

//-------------------------------------------------------------------------
// Purpose       : Get type of TopologyEntity this GeometryEntity
//                 should be attached to.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/14/01
//-------------------------------------------------------------------------
const type_info& Curve::topology_entity_type_info() const
{ return typeid(RefEdge); }

CubitSense Curve::relative_sense(Curve *other_curve)
{
  CubitVector center, junk, this_tangent, other_tangent;
  double start, end;
  
  if (!get_param_range(start,end)                 ||
      !position_from_u(0.5*(start+end), center)   ||
      !closest_point(center, junk, &this_tangent) ||
      !other_curve->closest_point(center, junk, &other_tangent))
    return CUBIT_UNKNOWN;
  
  double dot = this_tangent % other_tangent;
  return dot == 0.0 ? CUBIT_UNKNOWN : dot > 0. ? CUBIT_FORWARD : CUBIT_REVERSED;
}

//- When this is implemented in VirtualGeometryEngine it should be a pure vitual
CubitStatus Curve::get_point_direction( CubitVector& , CubitVector&  )
{
   return CUBIT_FAILURE;
}

//- default implementation
CubitStatus Curve::get_center_radius( CubitVector& c, double& r)
{
  CubitVector start_pos, closest_pos;
  if ( position_from_u(start_param(), start_pos) &&
       closest_point(start_pos, closest_pos, 0, &c) ) 
  {
    r = c.length();
    if (r > CUBIT_RESABS) 
    {
      r = 1.0 / r;
      c *= r * r;
      c += closest_pos;
      return CUBIT_SUCCESS;
    }
  }
  return CUBIT_FAILURE; 
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the point that is a distance
//                 arc length from the root_point.
//
// Special Notes : If the root point is on a periodic curve, it may
//                 accidentaly get the wrong parameter.  There is
//                 a special test here to make sure root point does
//                 not get set as the end parameter...
//                 If root point is at the end, we will assert later on..
//
// Creator       : David White
//
// Creation Date : 03/24/97
//-------------------------------------------------------------------------
CubitStatus Curve::point_from_arc_length ( const CubitVector &root_point,
                                           double const arc_length,
                                           CubitVector& new_point )
{
    // Get the parameter value of the root point
  double root_param = u_from_position( root_point );
  return point_from_arc_length(root_param, arc_length, new_point);
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the point that is a distance
//                 arc length from the point corresponding to root_param
//
// Special Notes : The original point_from_arc_length was split up
//                 to allow the starting parameter to be specified directly.
//                 This helps with periodic curves, where the start_param of the
//                 curve can be passed in so we're sure its starting from the
//                 start vertex.
//
// Creator       : David White (modified by Ved Vyas)
//
// Creation Date : 5/21/2006
//-------------------------------------------------------------------------

CubitStatus Curve::point_from_arc_length( double root_param,
                                          const double arc_length,
                                          CubitVector &new_point )
{
  double low_param = start_param();
  double high_param = end_param();

  if (high_param < low_param)
  {
    double temp = high_param;
    high_param = low_param;
    low_param = temp;
  }

    // Comment: The way we handle points not on the bounded curve is
    //          different for periodic and non-periodic curves!
    //          Periodic curves leave the parameter off of the curve,
    //          while non-periodic curves move the point to the closest
    //          endpoint.  Just an observation!
  
  
    // Adjust the parameter for periodic curves
  double periodic_param;

  if ( is_periodic( periodic_param ))
  {
    while (root_param < low_param)
      root_param += periodic_param;
    while (root_param > high_param)
      root_param -= periodic_param;
      // If you're moving in the positive direction and
      // you're almost at the start point...
    if (fabs( this->end_param() - root_param ) <= CUBIT_RESABS &&
        arc_length > 0.0 )
    {
        //the root param should be switched with the start param.
      root_param = this->start_param();
    }
      // If you're moving in the negative direction and
      // you're almost at the end point...
    else if (fabs( this->start_param() - root_param ) <= CUBIT_RESABS &&
             arc_length < 0.0 )
    {
      root_param = this->end_param();
    }
  }
  else if ( root_param < (low_param + CUBIT_RESABS) )
    root_param = low_param;
  else if ( root_param > (high_param - CUBIT_RESABS) )
    root_param = high_param;
 
    // Get the parameter value of the new point
  double new_param = this->u_from_arc_length( root_param,
                                              arc_length );
  
    // Now get the coordinates (in world space) representing this parameter
    // value
  return this->position_from_u(new_param, new_point);
}

//-------------------------------------------------------------------------
// Purpose       : Find the closest point on a BOUNDED curve.
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 07/17/98
//-------------------------------------------------------------------------
CubitStatus Curve::closest_point_trimmed( CubitVector const& from_pt,
                                          CubitVector& result )
{
    // Get the untrimmed point
  double param;
  if ( !closest_point(from_pt, result, NULL, NULL, &param) )
    return CUBIT_FAILURE;  
  
  double param_range = 0.0;
  double period, start_param, end_param;
  
    // Get whether periodic
  CubitBoolean is_per = this->is_periodic(period);
    // Get the parameter range
  get_param_range( start_param, end_param );
  
    // Make sure the start_param is lower than end_param.
  if (start_param > end_param)
  {
      // use param_range as temp
    param_range = start_param;
    start_param = end_param;
    end_param = param_range;
  }
  param_range = end_param - start_param;
  
    // If the Curve does not loop onto itself...
  if(  is_per == CUBIT_FALSE  ||
       param_range < period)
  {
      // If not within parameter range, return
      // the the closest endpoint
    if( (param < start_param) || (param > end_param) )
    {
      CubitVector start, end;
      position_from_u( start_param, start );
      position_from_u( end_param, end );
      result = ( (start - result).length_squared() < 
                 (end - result).length_squared() ) ? start : end ;
    }
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Check for discontinuity in the tangents of a curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/01/00
//-------------------------------------------------------------------------
CubitBoolean Curve::G1_discontinuous( double u,
	CubitVector* minus_tangent_r, CubitVector* plus_tangent_r )
{
	double start_param, end_param;
	if( ! get_param_range( start_param, end_param ) ) return CUBIT_FALSE;
	//if the curve is not parametric, just return false.
	if( start_param > end_param )
	{
		double tmp = start_param;
		start_param = end_param;
		end_param = tmp;
	}
	
	double u_minus = (2 * CUBIT_RESABS);
	double u_plus = u_minus;
	
	//can't be C1 discontinous at an end point!
	if( ((u - u_minus) < start_param) || ((u + u_plus) > end_param) )
		return CUBIT_FALSE;

	CubitVector position;
	position_from_u( u, position );
	const double res_abs_sqr = CUBIT_RESABS * CUBIT_RESABS;

	CubitVector minus_position, minus_tangent;
	position_from_u( u - u_minus, minus_position );
	double u_minus2 = u_from_position( minus_position );
	minus_tangent = position - minus_position;

	while( ((u - u_minus2) > (u - CUBIT_RESABS)) || 
	       (minus_tangent.length_squared() < res_abs_sqr) )
	{
		u_minus *= 10.;
		//can't be C1 discontinous at an end point!
		if( (u - u_minus) < start_param ) return CUBIT_FALSE;
		position_from_u( u - u_minus, minus_position );
		u_minus2 = u_from_position( minus_position );
		minus_tangent = position - minus_position;
	}
	
	CubitVector plus_position, plus_tangent;
	position_from_u( u + u_plus, plus_position );
	double u_plus2 = u_from_position( plus_position );
	plus_tangent = plus_position - position;
	
	while( ((u + u_plus2) < (u + CUBIT_RESABS)) ||
	       ((position - plus_position).length_squared() < res_abs_sqr) )
	{
		u_plus *= 10.;
		//can't be C1 discontinous at an end point!
		if( (u + u_plus) > end_param ) return CUBIT_FALSE;
		position_from_u( u + u_plus, plus_position );
		u_plus2 = u_from_position( plus_position );
		plus_tangent = plus_position - position;
	}
	
	plus_tangent.normalize();
	minus_tangent.normalize();
	if( (plus_tangent * minus_tangent).length_squared() > (2 * res_abs_sqr) )
	{
		if( plus_tangent_r ) *plus_tangent_r = plus_tangent;
		if( minus_tangent_r ) * minus_tangent_r = minus_tangent;
		return CUBIT_TRUE;
	}
	
	return CUBIT_FALSE; 
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

