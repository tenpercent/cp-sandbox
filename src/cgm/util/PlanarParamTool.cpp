//- Class: PlanarParamTool
//-------------------------------------------------------------------------
// Filename      : PlanarParamTool.cpp
//
// Purpose       : Surface Parameterization for mesh by triangulation flattening
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#include "PlanarParamTool.hpp"
#include "CastTo.hpp"
#include "GeometryDefines.h"


//-------------------------------------------------------------------------
// Function:    PlanarParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
PlanarParamTool::PlanarParamTool()
{
}

//-------------------------------------------------------------------------
// Function:    PlanarParamTool
// Description: deconstructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
PlanarParamTool::~PlanarParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PlanarParamTool::set_up_space(CubitVector& du, CubitVector& dv, CubitVector& uv_center)
{
        Du = du;
        Dv = dv;
        uvCenter = uv_center;

	return CUBIT_SUCCESS; 
}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PlanarParamTool::transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location) 
{
  // Translate to local origin at center

  CubitVector vect = xyz_location - uvCenter;

   // Multiply by transpose (inverse) of transformation vector

  uv_location.x( vect % Du );
  uv_location.y( vect % Dv );
  uv_location.z( 1.0 );

  return CUBIT_SUCCESS;
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PlanarParamTool::transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location) 
{
// Multiply by transformation matrix

  CubitVector vect;
  vect.x( uv_location.x() * Du.x() +
          uv_location.y() * Dv.x() );
  vect.y( uv_location.x() * Du.y() +
          uv_location.y() * Dv.y() );
  vect.z( uv_location.x() * Du.z() +
          uv_location.y() * Dv.z() );

   // Translate from origin

  xyz_location = vect + uvCenter;

  return CUBIT_SUCCESS;
}



//EOF
