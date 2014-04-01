//-------------------------------------------------------------------------
// Filename      : ParamCubitPlane.hpp
//
// Purpose       : A class derived from CubitPlane with a parameterization,
//                 essentially a planar surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/06/98
//-------------------------------------------------------------------------
#ifndef PARAM_CUBIT_PLANE_HPP
#define PARAM_CUBIT_PLANE_HPP

#include "CubitPlane.hpp"
#include "CubitUtilConfigure.h"

class CubitVector;


class CUBIT_UTIL_EXPORT ParamCubitPlane : public CubitPlane
{
	public:

	ParamCubitPlane( const CubitVector& zero_position,
	                 const CubitVector& u_direction,
									 const CubitVector& v_direction,
									 CubitStatus& result );
	//I zero_position
	//I- The position on the plane where (u,v) = (0,0).
	//I u_direction
	//I- A vector parallel to the plane along which u increases
	//I- and v is constant.
	//I v_direction
	//I- A vector parallel to the plane along which v increases
	//I- and u is constant.
	//O result
	//O- CUBIT_FAILURE if u_direction and v_direction are parallel, or
	//O- are of zero length.
	//- Create a parametric CubitPlane with the specified
	//- parameterization.
	//-
	//- Note:  The size of a step in the u and v directions
	//- is determined by the lenght of the respective 
	//- u and v direction vectors passed.
	
	ParamCubitPlane( const ParamCubitPlane& copy_this );
	//- Copy constructor.
	
	ParamCubitPlane( double A, double B, double C, double D );
	ParamCubitPlane( const CubitVector& normal, double D );
	ParamCubitPlane( const CubitVector& normal, const CubitVector& point );
	ParamCubitPlane( DLIList<CubitVector*>& positions );
	ParamCubitPlane( const CubitPlane& copy_this );
	//- Constructors from CubitPlane
	//- An arbitraty parameterization is generated.
	
#ifdef BOYD15
	CubitStatus u_v_from_position( const CubitVector& position,
	                               double& u, double& v ) const;
	//R CubitStatus
	//R- CUBIT_SUCCESS/CUBIT_FAILURE
	//I position
	//I- A position on the plane
	//O u, v
	//O- The parameter values at the specified position.
	//- Find the parameter values u and v at the specified
	//- position on the plane.
	
	CubitStatus position_from_u_v( double u, double v,
	                               CubitVector& position ) const;
	//R CubitStatus
	//R- CUBIT_SUCCESS/CUBIT_FAILURE
	//I u, v
	//I- parameter values
	//O position
	//O- The position on the surface at the specified parameter values.
	//- Given parameter values on the surface, pass back the 
	//- corresponding position.
#endif
	
	CubitStatus closest_point( const CubitVector& position,
	                           CubitVector& closest_position ) const;
	//R CubitStatus
	//R- CUBIT_SUCCESS/CUBIT_FAILURE
	//I position
	//I- A position to evaluate
	//O closest_position
	//O- The projection of the passed position onto the surface
	//O- along the normal.
	//- Find the closest point on the plane.
	
	CubitStatus move_to_plane( CubitVector& position ) const;
	//R CubitStatus
	//R- CUBIT_SUCCESS/CUBIT_FAILURE
	//I position
	//I- A position to evaluate
	//O position
	//O- The projection of the passed position onto the surface
	//O- along the normal.
	//- Move the passed position to the plane.
	
	private:
	
	void make_parameterization();
	//- This method is used by various constructors to generate
	//- an arbitrary parameterization after the members of the
	//- parent class CubitPlane have been initialized.
	
	CubitVector p_, s_, t_;
	CubitBoolean is_plane_valid_;
	CubitVector n_; //Non-unit normal vector.
	double n_epsilon_;
};
	  
#endif

