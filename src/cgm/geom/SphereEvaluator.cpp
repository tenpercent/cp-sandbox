//-------------------------------------------------------------------------
// Filename      : FacetSurface.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 06/06/00
//
// Owner         : David R. White
//-------------------------------------------------------------------------

#include "SphereEvaluator.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "CubitVector.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Compute and return the bounding box for this sphere.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox SphereEvaluator::bounding_box( void ) const
{
    CubitVector center = mTmatrix * mEvalData.center;
            
    CubitVector min( center.x() - mEvalData.radius,
                     center.y() - mEvalData.radius,
                     center.z() - mEvalData.radius );
    CubitVector max( center.x() + mEvalData.radius,
                     center.y() + mEvalData.radius,
                     center.z() + mEvalData.radius );

    CubitBox thebbox( min, max );

    return thebbox;
}

//-------------------------------------------------------------------------
// Purpose       : This functions computes the point on the surface that is 
//                 closest to the input location and then calculates the 
//                 magnitudes of the principal curvatures at this (possibly, 
//                 new) point on the surface.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus SphereEvaluator::principal_curvatures(
  CubitVector const& location,
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{

    curvature_1 = curvature_2 = 1.0 / mEvalData.radius;

    if ( closest_location )
    {
        closest_point( location, closest_location );
    }
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
//-------------------------------------------------------------------------
CubitStatus SphereEvaluator::closest_point( CubitVector const& location, 
                                            CubitVector* closest_location,
                                            CubitVector* unit_normal_ptr,
                                            CubitVector* curvature1_ptr,
                                            CubitVector* curvature2_ptr) const
{
    CubitTransformMatrix inverse_Tmatrix = mTmatrix;
    inverse_Tmatrix.inverse();
    CubitVector new_pt = inverse_Tmatrix * location;

    // ******************************************************
    // If requested, compute the normal at the input point.
    // ******************************************************

    if ( unit_normal_ptr != NULL )
    {
        CubitVector origin( 0.0, 0.0, 0.0 );
        CubitVector endpt = new_pt - mEvalData.center;

        origin = mTmatrix * origin;
        endpt  = mTmatrix * endpt;

        *unit_normal_ptr = endpt - origin;
        unit_normal_ptr->normalize();
    }

    // *************************************************************
    // If requested, compute the closest point to the input point.
    // *************************************************************

    if ( closest_location != NULL )
    {
        if ( location.within_tolerance( mEvalData.center, GEOMETRY_RESABS ) )
        {
            closest_location->set( mEvalData.center.x() + mEvalData.radius,
                                   mEvalData.center.y(),
                                   mEvalData.center.z() );
        }
        else
        {
            CubitVector vec = new_pt - mEvalData.center;
            vec.normalize();
            *closest_location = mEvalData.center + ( vec * mEvalData.radius );
            *closest_location = mTmatrix * (*closest_location);
        }
    }

    // ***********************************************************************
    // If requested, compute the curvature directions.
    // ***********************************************************************

    if ( curvature1_ptr &&
         curvature2_ptr )
    {
        PRINT_ERROR("Sphere's do not current support curvature pointer requests.\n");
        return CUBIT_FAILURE;
    }
  
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Given a UV parametric location, return the cooresponding
//                 XYZ location.
//-------------------------------------------------------------------------
CubitVector SphereEvaluator::position_from_u_v( double u, double v ) const
{
    double two_pi = 2.0 * CUBIT_PI;

    while ( u > two_pi ) u -= two_pi;
    while ( u < 0.0    ) u += two_pi;

    double radius_modified = mEvalData.radius * sin ( v );
    double x = radius_modified * cos( u );
    double y = radius_modified * sin( u );
    double z = mEvalData.radius * cos( v );

    CubitVector position( x, y, z );
    position = mTmatrix * position;

    return position;
}

//-------------------------------------------------------------------------
// Purpose       : Project a given XYZ position to a surface and return
//                 the UV and XYZ locations on the surface.
//-------------------------------------------------------------------------
CubitStatus SphereEvaluator::u_v_from_position
(
    CubitVector const& location,
    double& u,
    double& v,
    CubitVector* closest_location
) const
{
    CubitTransformMatrix inverse_Tmatrix = mTmatrix;
    inverse_Tmatrix.inverse();
    CubitVector transformed_pt = inverse_Tmatrix * location;

    CubitVector dir( transformed_pt );
    dir.normalize();

    v = acos( dir.z() );

    if ( v < GEOMETRY_RESABS )
    {
        u = 0.0;
    }
    else
    {
        dir.set( dir.x(), dir.y(), 0.0 );
        dir.normalize();
        u = acos( dir.x() );

        if ( dir.y() < 0.0 )
        {
            u = 2*CUBIT_PI - u;
        }
        while ( u < 0.0 ) u += 2*CUBIT_PI;
        while ( u >= 2*CUBIT_PI ) u -= 2*CUBIT_PI;
    }

    if ( closest_location )
        closest_point( location, closest_location );
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : returns a structure containing the parameters for this
//                 sphere (i.e. radius, center, etc.)
//
//-------------------------------------------------------------------------
const CubitEvaluatorData* SphereEvaluator::evaluator_data( void ) const
{
    return &mEvalData;
}

//-------------------------------------------------------------------------
// Purpose       : return parametric extremes in the U direction
//
//-------------------------------------------------------------------------
CubitBoolean SphereEvaluator::get_param_range_U
(
    double& lower_bound,
    double& upper_bound
) const
{
   lower_bound = 0.0;
   upper_bound = 2.0*CUBIT_PI;
   return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : return parametric extremes in the V direction
//
//-------------------------------------------------------------------------
CubitBoolean SphereEvaluator::get_param_range_V
(
    double& lower_bound,
    double& upper_bound
) const
{
   lower_bound = 0;
   upper_bound = CUBIT_PI;
   return CUBIT_TRUE;
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
