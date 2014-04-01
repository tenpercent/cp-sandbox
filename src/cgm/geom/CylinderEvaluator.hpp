//-------------------------------------------------------------------------
// Filename      : CylinderEvaluator.hpp
//
// Purpose       : Provide evaluators for a cylindrical surface.  Since a cylinder
//                 is an analytic primitive, its evaluators can be overwritten
//                 rather than using the ACIS, or Facet based evaluators.
//
// Special Notes :
//
// Creator       : Matthew Staten
//
// Creation Date : 09/14/04
//
// Owner         : Matthew Staten
//-------------------------------------------------------------------------

#ifndef CYLINDER_EVALUATOR_HPP
#define CYLINDER_EVALUATOR_HPP

#include "CubitEvaluator.hpp"
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "GeometryDefines.h"

class CUBIT_GEOM_EXPORT CylinderEvaluatorData : public CubitEvaluatorData
{
public:
    double height;
    double base_radius_x;
    double base_radius_y;
    double top_radius;

    virtual GeometryType ask_type() const { return CONE_SURFACE_TYPE; }
};

class CUBIT_GEOM_EXPORT CylinderEvaluator: public CubitEvaluator
{

private:

    CylinderEvaluatorData mEvalData;

protected: 

#ifdef BOYD14
    CubitVector base_transformed() const;
#endif
    CubitVector base() const;

public :
  
    CylinderEvaluator( const CylinderEvaluatorData *data );

  virtual ~CylinderEvaluator()
    {}
  
    virtual GeometryType ask_type() { return CONE_SURFACE_TYPE; }
    virtual CubitBoolean is_parametric() const { return CUBIT_TRUE; }
    virtual CubitBoolean is_periodic() const { return CUBIT_TRUE; }
    virtual CubitBoolean is_periodic_in_U( double &period ) const { period = 2*CUBIT_PI;
                                                                    return CUBIT_TRUE; }
    virtual CubitBoolean is_periodic_in_V( double &period ) const { period = mEvalData.height;
                                                                    return CUBIT_FALSE; }
    virtual CubitBoolean is_singular_in_U() const { return CUBIT_FALSE; }
    virtual CubitBoolean is_singular_in_V() const { return CUBIT_FALSE; }
    virtual CubitBoolean is_closed_in_U() const { return CUBIT_TRUE; }
    virtual CubitBoolean is_closed_in_V() const { return CUBIT_FALSE; }
    virtual CubitBoolean get_param_range_U( double& lower_bound,
                                            double& upper_bound ) const;
    virtual CubitBoolean get_param_range_V( double& lower_bound,
                                            double& upper_bound ) const;
    virtual CubitVector position_from_u_v( double u, double v ) const;
    virtual CubitStatus u_v_from_position( CubitVector const& location,
                                           double& u,
                                           double& v,
                                           CubitVector* closest_location ) const;

    virtual const CubitEvaluatorData* evaluator_data() const;

    virtual CubitBox bounding_box() const;
    //- Return the bounding box of this cylinder.

    virtual CubitStatus principal_curvatures( CubitVector const& location, 
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

    virtual CubitStatus closest_point( CubitVector const& location,
                                       CubitVector* closest_location = NULL,
                                       CubitVector* unit_normal_ptr = NULL,
                                       CubitVector* curvature1_ptr = NULL,
                                       CubitVector* curvature2_ptr = NULL ) const;
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

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

