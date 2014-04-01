//- Class: SurfParamTool
//-------------------------------------------------------------------------
// Filename      : SurfParamTool.hpp
//
// Purpose       : This is the generic version of ParamTool when the
//                 geometry engine has a sufficient parameterization.
//                 It uses the tool's existing functions to get transform
//                 between uv and xyz spaces.
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#ifndef SURF_PARAM_TOOL_HPP
#define SURF_PARAM_TOOL_HPP

#include "ParamTool.hpp"
#include "CubitGeomConfigure.h"

class Surface;
class CubitVector;
template <class X> class DLIList;


class CUBIT_GEOM_EXPORT SurfParamTool : public ParamTool
{
public:

	//- constructor
	SurfParamTool(Surface *surf);

	//- deconstructor
	~SurfParamTool();
	
	CubitStatus set_up_space(void);

	CubitStatus transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location);

	CubitStatus transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location);


private:

	Surface *refSurf;
	//- reference surface
};

#endif // ACIS_PARAM_TOOL_HPP

