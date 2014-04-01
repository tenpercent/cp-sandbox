//-------------------------------------------------------------------------
// Filename      : CylinderEvaluator.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Matt Staten
//
// Creation Date : 09/15/04
//
// Owner         : Matt Staten
//-------------------------------------------------------------------------

#include "CylinderEvaluator.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "CubitVector.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Construct a CylinderEvaluator object given a structure
//                 containing the radius, center, etc. for the desired
//                 Cylinder.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CylinderEvaluator::CylinderEvaluator( const CylinderEvaluatorData *data )
{
    mEvalData = *data;

    mTmatrix.scale_about_origin( 1.0, mEvalData.base_radius_y/mEvalData.base_radius_x, 1.0 );
}

//-------------------------------------------------------------------------
// Purpose       : Compute and return the bounding box for this cylinder.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox CylinderEvaluator::bounding_box( void ) const
{
    int ipt = 0;
    CubitVector pts[8];
    CubitVector axis( 0.0, 0.0, 1.0 );
    CubitVector base_pt = base();
    CubitVector top_pt = base() + mEvalData.height * axis;
    CubitVector x_dir( 1.0, 0.0, 0.0 );
    CubitVector y_dir( 0.0, 1.0, 0.0 );
    double scale = mEvalData.base_radius_y * mEvalData.base_radius_x;

    pts[0] = base_pt + ( x_dir * mEvalData.base_radius_x );
    pts[1] = base_pt - ( x_dir * mEvalData.base_radius_x );
    pts[2] = base_pt + ( y_dir * mEvalData.base_radius_y );
    pts[3] = base_pt - ( y_dir * mEvalData.base_radius_y );
    pts[4] = top_pt  + ( x_dir * mEvalData.top_radius );
    pts[5] = top_pt  - ( x_dir * mEvalData.top_radius );
    pts[6] = top_pt  + ( y_dir * mEvalData.top_radius * scale );
    pts[7] = top_pt  - ( y_dir * mEvalData.top_radius * scale );

    for ( ipt = 0; ipt < 8; ipt++ )
    {
        pts[ipt] = mTmatrix * pts[ipt];
    }

    CubitVector min( pts[0] ),
                max( pts[0] );
    for ( ipt = 1; ipt < 8; ipt++ )
    {
        CubitVector *pt = &(pts[ipt]);
        if ( min.x() > pt->x() ) min.x( pt->x() );
        if ( min.y() > pt->y() ) min.y( pt->y() );
        if ( min.z() > pt->z() ) min.z( pt->z() );
        if ( max.x() < pt->x() ) max.x( pt->x() );
        if ( max.y() < pt->y() ) max.y( pt->y() );
        if ( max.z() < pt->z() ) max.z( pt->z() );
    }

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
CubitStatus CylinderEvaluator::principal_curvatures(
  CubitVector const& location,
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
    PRINT_ERROR("principles_curvatures not implemented for cylinders yet.\n");
    return CUBIT_FAILURE;

//    if ( closest_location )
//    {
//        return closest_point( location, closest_location );
//    }
//    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
//-------------------------------------------------------------------------
CubitStatus CylinderEvaluator::closest_point( CubitVector const& location, 
                                              CubitVector* closest_location,
                                              CubitVector* unit_normal_ptr,
                                              CubitVector* curvature1_ptr,
                                              CubitVector* curvature2_ptr) const
{
    CubitTransformMatrix inverse_Tmatrix = mTmatrix;
    inverse_Tmatrix.inverse();
    CubitVector new_pt = inverse_Tmatrix * location;
    CubitVector normal;
    double dist = sqrt( ( new_pt.x() * new_pt.x() ) + ( new_pt.y() * new_pt.y() ) );
    if ( dist < GEOMETRY_RESABS )
    {
        normal.set( 1.0, 0.0, 0.0 );
    }
    else
    {
        normal.set( new_pt.x(), new_pt.y(), 0.0 );
    }
    normal.normalize();
    if ( unit_normal_ptr )
    {
        CubitVector origin( 0.0, 0.0, 0.0 );
        CubitVector endpt( normal );

        origin = mTmatrix * origin;;
        endpt  = mTmatrix * endpt;

        *unit_normal_ptr = endpt - origin;
        unit_normal_ptr->normalize();
    }
    if ( closest_location != NULL )
    {
        CubitVector axis_pt( 0.0, 0.0, new_pt.z() );
        double v = axis_pt.z() - base().z();
        double radius = v * ( ( mEvalData.top_radius - mEvalData.base_radius_x ) / mEvalData.height ) + mEvalData.base_radius_x;

        *closest_location = axis_pt + ( radius * normal );
        *closest_location = mTmatrix * (*closest_location);
    }
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Given a UV parametric location, return the cooresponding
//                 XYZ location.
//-------------------------------------------------------------------------
CubitVector CylinderEvaluator::position_from_u_v( double u, double v ) const
{
    double radius = v * ( ( mEvalData.top_radius - mEvalData.base_radius_x ) / mEvalData.height ) + mEvalData.base_radius_x;
    CubitVector base_pt = base();
    double two_pi = 2.0 * CUBIT_PI;

    while ( u > two_pi ) u -= two_pi;
    while ( u < 0.0    ) u += two_pi;

    double x = radius * sin( u );
    double y = radius * cos( u );
    double z = base_pt.z() + v;

    CubitVector position( x, y, z );
    position = mTmatrix * position;

    return position;
}

//-------------------------------------------------------------------------
// Purpose       : Project a given XYZ position to a surface and return
//                 the UV and XYZ locations on the surface.
//-------------------------------------------------------------------------
CubitStatus CylinderEvaluator::u_v_from_position
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
    CubitVector base_pt = base();
    CubitVector dir( transformed_pt.x(), transformed_pt.y(), 0.0 );
    dir.normalize();

    v = transformed_pt.z() - base_pt.z();
    u = acos( dir.y() );
    if ( dir.x() < 0.0 )
        u = 2*CUBIT_PI - u;

    while ( u < 0.0 ) u += 2*CUBIT_PI;
    while ( u >= 2*CUBIT_PI ) u -= 2*CUBIT_PI;

    if ( closest_location )
        closest_point( location, closest_location );
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : return the dimension parameters of this cylinder.
//
//-------------------------------------------------------------------------
const CubitEvaluatorData* CylinderEvaluator::evaluator_data( void ) const
{
    return &mEvalData;
}

//-------------------------------------------------------------------------
// Purpose       : return parametric extremes in the U direction
//
//-------------------------------------------------------------------------
CubitBoolean CylinderEvaluator::get_param_range_U
(
    double& lower_bound,
    double& upper_bound
) const
{
   lower_bound = 0.0;
   upper_bound = 2*CUBIT_PI;
   return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : return parametric extremes in the V direction
//
//-------------------------------------------------------------------------
CubitBoolean CylinderEvaluator::get_param_range_V
(
    double& lower_bound,
    double& upper_bound
) const
{
   lower_bound = 0.0;
   upper_bound = mEvalData.height;
   return CUBIT_TRUE;
}

// ********** END PUBLIC FUNCTIONS         **********

CubitVector CylinderEvaluator::base() const
{
    return CubitVector( 0.0, 0.0, -0.5*mEvalData.height );
}

// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
