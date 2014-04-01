

#ifndef CGM_ENGINE_DYNAMIC_LOADER_HPP
#define CGM_ENGINE_DYNAMIC_LOADER_HPP

#include "CubitDynamicLoader.hpp"
#include "CubitString.hpp"
#include "CubitGeomConfigure.h"

class GeometryQueryEngine;
class GeometryModifyEngine;
  
//! provides functionality from cgm engines residing in dynamically loadable libraries
class CUBIT_GEOM_EXPORT CGMEngineDynamicLoader
{
  public:
    //! create an instance of this engine loader with an engine name
    CGMEngineDynamicLoader(const CubitString& engine_name);
    //! delete this and unload the library if it is loaded
    virtual ~CGMEngineDynamicLoader();
    
    //! get the name of this engine
    CubitString get_engine_name();
    
    //! can set library name for the engine if it is different than
    //! the engine name
    void set_library_base_name(const CubitString& libname);
    CubitString get_library_base_name();

    //! check to see if the engine exists
    CubitBoolean engine_exists();
    
    //! load the engine, returns success or failure
    CubitStatus load_engine();
    //! unload the engine, returns success or failure
    CubitStatus unload_engine();
    
    //! get the GeometryQueryEngine instance from the engine
    //! calls load_engine if it hasn't been called yet
    //! returns NULL if library failed to load, or an engine couldn't be made
    GeometryQueryEngine* get_gqe();
    
    //! get the GeometryModifyEngine instance from the engine
    //! calls load_engine if it hasn't been called yet
    //! returns NULL if library failed to load, or an engine couldn't be made
    GeometryModifyEngine* get_gme();

  protected:

    // the name of the engine
    const CubitString mEngineName;

    // the filename of the engine library
    CubitString mEngineLibrary;

    // the handle of the engine library
    CubitDynamicLoader::LibraryHandle mLibraryHandle;
    
    // whether attempted to load the library
    CubitBoolean mLoadAttempted; 
    
    // compare versions of CGM that the engine was built with
    // and the version that is expected to be found
    CubitBoolean check_engine_version();

    GeometryQueryEngine* mQueryEngine;
    GeometryModifyEngine* mModifyEngine;
};


// macro to define symbols in CGM engine to allow 
// this dynamic loader to resolve

#if defined(WIN32)
#define CGM_EXPORT_ENGINE_SYM __declspec(dllexport)
#else
#define CGM_EXPORT_ENGINE_SYM
#endif

#define CGM_ENGINE_EXPORT_CREATE_GQE(EngineName) \
  extern "C" CGM_EXPORT_ENGINE_SYM const char* dynamic_##EngineName##_get_cgm_version() \
  { \
    return "1.0";  /* TODO -- use real versioning system, or not at all? */  \
  } \
  extern "C" CGM_EXPORT_ENGINE_SYM GeometryQueryEngine* dynamic_##EngineName##_create_gqe()

#define CGM_ENGINE_EXPORT_CREATE_GME(EngineName) \
  extern "C" CGM_EXPORT_ENGINE_SYM GeometryModifyEngine* dynamic_##EngineName##_create_gme() \


#endif // CGM_ENGINE_DYNAMIC_LOADER_HPP

