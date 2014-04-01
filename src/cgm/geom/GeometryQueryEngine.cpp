#include "GeometryQueryEngine.hpp"
#include "CubitString.hpp"


#include "Surface.hpp"
#include "Curve.hpp"



CubitStatus GeometryQueryEngine::set_export_version(const int,
                                                    const int) 
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::set_export_allint_version(const int) 
{
  return CUBIT_FAILURE;
}

int GeometryQueryEngine::get_allint_version() 
{
  return 100*get_major_version() + get_minor_version();
}

CubitStatus GeometryQueryEngine::list_engine_versions(CubitString &versions)
{
  versions = "(none specified)";
  return CUBIT_SUCCESS;
}

//Added functions for removing the geom classes dependecy on the virtual classes

 
CubitStatus GeometryQueryEngine::get_underlying_curves(Curve * curve_ptr,
                                 DLIList<TopologyBridge*>& curve_list)
{ return CUBIT_SUCCESS; }

CubitStatus GeometryQueryEngine::get_underlying_surfaces(Surface * surf_ptr,
                                 DLIList<TopologyBridge*>& surf_list)
{ return CUBIT_SUCCESS; }
