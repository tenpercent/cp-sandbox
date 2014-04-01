//- Class: PlanarParamTool
//-------------------------------------------------------------------------
// Filename      : PlanarParamTool.hpp
//
// Purpose       : Surface Parameterization for mesh by triangulation flattening specific for Planar surfaces
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#ifndef PLANAR_PARAM_TOOL_HPP
#define PLANAR_PARAM_TOOL_HPP

#include "ParamTool.hpp"
#include "CubitVector.hpp"

class CUBIT_UTIL_EXPORT PlanarParamTool : public ParamTool
{
public:

	//- constructor
	PlanarParamTool();

	//- deconstructor
	~PlanarParamTool();
	
	CubitStatus set_up_space(CubitVector& du, CubitVector& dv, CubitVector& uv_center);

	CubitStatus transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location);

	CubitStatus transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location);


private:
	
	CubitVector Du;
	CubitVector Dv;
	CubitVector uvCenter;
   //- transformation info - used with set_up_space
};

#endif // PLANAR_PARAM_TOOL_HPP

