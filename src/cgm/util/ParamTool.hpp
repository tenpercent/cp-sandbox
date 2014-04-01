//- Class: ParamTool
//-------------------------------------------------------------------------
// Filename      : ParamTool.hpp
//
// Purpose       : Surface Parameterization for mesh by triangulation flattening
//
// Note:  set_up_space must be called before using the transform functions
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#ifndef PARAM_TOOL_HPP
#define PARAM_TOOL_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

class CubitVector;


class CUBIT_UTIL_EXPORT ParamTool
{
public:

    ParamTool() {}
    virtual ~ParamTool() {}

    virtual CubitStatus transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location) = 0;

    virtual CubitStatus transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location) = 0;

};

#endif // PARAM_TOOL_HPP

