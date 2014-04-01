//-------------------------------------------------------------------------
// Filename      : ParamCubitPlane.cc
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------

#include "ParamCubitPlane.hpp"
#include "CubitVector.hpp"
#include "math.h"

//-------------------------------------------------------------------------
// Purpose       : Default constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------
ParamCubitPlane::ParamCubitPlane( const CubitVector& zero_position,
                                  const CubitVector& u_direction,
																	const CubitVector& v_direction,
																	CubitStatus& result )
{
	if( mk_plane_with_points( zero_position, zero_position + u_direction, 
	    zero_position + v_direction ) == (int)CUBIT_SUCCESS )
	{
		is_plane_valid_= CUBIT_TRUE;
		result = CUBIT_SUCCESS;
	}
	else
	{
		is_plane_valid_= CUBIT_FALSE;
		result = CUBIT_FAILURE;
		return ;
	}
	
	p_ = zero_position;
	s_ = u_direction;
	t_ = v_direction;
	n_ = s_ * t_;
	
	double n_len = n_.length();
	if( 1000 * CUBIT_DBL_MIN > n_len ) n_epsilon_ = CUBIT_DBL_MIN;
	else n_epsilon_ = n_len / 1000;
}


//-------------------------------------------------------------------------
// Purpose       : closest point
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------
CubitStatus ParamCubitPlane::closest_point( const CubitVector& position,
                             CubitVector& closest_position ) const
{
	assert( is_plane_valid_);
	if( !is_plane_valid_) return CUBIT_FAILURE;

	closest_position = position - (normal() * distance(position));
	return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : closest point
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------
CubitStatus ParamCubitPlane::move_to_plane( CubitVector& position ) const
{
	const CubitVector v = position;
	CubitStatus s = closest_point( v, position );
	return s;
}		

//-------------------------------------------------------------------------
// Purpose       : make arbitrary parameterization
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------
void ParamCubitPlane::make_parameterization()
{
	//Choose the zero-point for the parameterization
	//as close to the origin as possible.
//	p_.set( 0.0, 0.0, 0.0);
//	move_to_plane( p_ );
	is_plane_valid_ = CUBIT_TRUE;
	
	const CubitVector temp_p(0.0, 0.0, 0.0);
	CubitStatus s = closest_point( temp_p, p_ );
	assert( s == CUBIT_SUCCESS );
	
	CubitVector n = normal();
	CubitVector p1;
	
	p1 = p_;
	double x = fabs( n.x() );
	double y = fabs( n.y() );
	double z = fabs( n.z() );
	
	//Choose a direction from the zero point (p_) for
	//the second point as the direction of the smallest
	//component of the normal.  The third point defining
	//the plane will be defined by the cross product of
	//the vector from the zero_point to this point and
	//the normal vector of the plane.
	if( (x <= y) && (x <= z) )
	{
		p1.x( p1.x() + 1 );
	}
	else if( (y <= x) && (y <= z) )
	{
		p1.y( p1.y() + 1 );
	}
	else
	{
		p1.z( p1.z() + 1 );
	}
	
	move_to_plane( p1 );
	s_ = p1 - p_;
	t_ = -(s_ * n);
	n_ = s_ * t_;
	
	double n_len = n_.length();
	if( 1000 * CUBIT_DBL_MIN > n_len ) n_epsilon_ = CUBIT_DBL_MIN;
	else n_epsilon_ = n_len / 1000;
	
 
 	is_plane_valid_= CUBIT_TRUE;
}
