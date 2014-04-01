

#include "CGMEngineDynamicLoader.hpp"

#include "CubitMessage.hpp"

extern "C"
{
  typedef GeometryQueryEngine* (*CGMEngineCreateQueryEngine)();
  typedef GeometryModifyEngine* (*CGMEngineCreateModifyEngine)();
}

// default engine library name is something like libCGM{engine_name}.so

CGMEngineDynamicLoader::CGMEngineDynamicLoader(const CubitString& engine_name)
  : mEngineName(engine_name), 
    mEngineLibrary(CubitString("CGM") + engine_name),
    mLibraryHandle(CubitDynamicLoader::InvalidLibraryHandle),
    mLoadAttempted(CUBIT_FALSE),
    mQueryEngine(NULL),
    mModifyEngine(NULL)
{
}


CGMEngineDynamicLoader::~CGMEngineDynamicLoader()
{
  if(mLibraryHandle != CubitDynamicLoader::InvalidLibraryHandle)
  {
    unload_engine();
  }
}

CubitString CGMEngineDynamicLoader::get_engine_name()
{
  return mEngineName;
}

CubitString CGMEngineDynamicLoader::get_library_base_name()
{
  return mEngineLibrary;
}

void CGMEngineDynamicLoader::set_library_base_name(const CubitString& libname)
{
  if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
  {
    mEngineLibrary = libname;
    mLoadAttempted = CUBIT_FALSE;
  }
}

CubitBoolean CGMEngineDynamicLoader::engine_exists()
{
  if(mLibraryHandle != CubitDynamicLoader::InvalidLibraryHandle)
    return CUBIT_TRUE;

  // try for a debug version with _d on it first.
#if defined(WIN32) && defined(_DEBUG)
  CubitString debug_libname = CubitDynamicLoader::library_prefix() +
            mEngineLibrary +
            CubitString("_d") +
            CubitDynamicLoader::library_extension();
  if(CubitDynamicLoader::library_exists(debug_libname.c_str()))
    return CUBIT_TRUE;
#endif
  
  CubitString libname = CubitDynamicLoader::library_prefix() +
            mEngineLibrary +
            CubitDynamicLoader::library_extension();
  if(CubitDynamicLoader::library_exists(libname.c_str()))
    return CUBIT_TRUE;

  return CUBIT_FALSE;
}

CubitStatus CGMEngineDynamicLoader::load_engine()
{
  if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
  {
    // try for a debug version with _d on it first.
#if defined(WIN32) && defined(_DEBUG)
    CubitString debug_libname = CubitDynamicLoader::library_prefix() +
              mEngineLibrary +
              CubitString("_d") +
              CubitDynamicLoader::library_extension();
    mLibraryHandle = CubitDynamicLoader::load_library(debug_libname.c_str());
#endif

    if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
    {
      CubitString libname = CubitDynamicLoader::library_prefix() +
                            mEngineLibrary +
                            CubitDynamicLoader::library_extension();

      mLibraryHandle = CubitDynamicLoader::load_library(libname.c_str());
    }

    mLoadAttempted = CUBIT_TRUE;
    if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
    {
      CubitString error = CubitDynamicLoader::get_error();
      if(error != CubitString())
      {
        PRINT_ERROR("%s\n", error.c_str());
      }
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus CGMEngineDynamicLoader::unload_engine()
{
  if(mLibraryHandle != CubitDynamicLoader::InvalidLibraryHandle)
  {
    CubitDynamicLoader::unload_library(mLibraryHandle);
    mLoadAttempted = CUBIT_FALSE;
    mLibraryHandle = CubitDynamicLoader::InvalidLibraryHandle;
    mQueryEngine = NULL;
    mModifyEngine = NULL;
  }
  return CUBIT_SUCCESS;
}


extern "C" {
  typedef const char* (*FpVersion)();
}

CubitBoolean CGMEngineDynamicLoader::check_engine_version()
{
  // library was already loaded, so just look up the version function

  CubitString symbol = CubitString("dynamic_") + mEngineName + CubitString("_get_cgm_version");
  FpVersion fp_version = (FpVersion)CubitDynamicLoader::get_symbol_address(mLibraryHandle, symbol.c_str());

  if(!fp_version)
  {
    PRINT_ERROR("Failed to load CGM %s engine : Can't find version\n", mEngineName.c_str());
    return CUBIT_FALSE;
  }

  // TODO -- use real versioning system
  if(strcmp((*fp_version)(), "1.0") != 0)
  {
    PRINT_ERROR("Failed to load CGM %s engine : wrong version number\n", mEngineName.c_str());
    return CUBIT_FALSE;
  }

  return TRUE;
}

GeometryQueryEngine* CGMEngineDynamicLoader::get_gqe()
{
  if(mQueryEngine)
    return mQueryEngine;

  if(mLoadAttempted == CUBIT_FALSE)
    load_engine();

  if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
  {
    PRINT_ERROR("Failed to load CGM %s engine\n", mEngineName.c_str());
    return NULL;
  }

  // check for incompatible CGM version
  if(check_engine_version() == CUBIT_FALSE)
    return NULL;

  CubitString symbol = CubitString("dynamic_") + mEngineName + CubitString("_create_gqe");
  
  CGMEngineCreateQueryEngine create_query_engine = 
    (CGMEngineCreateQueryEngine)CubitDynamicLoader::get_symbol_address(mLibraryHandle, symbol.c_str());

  if(!create_query_engine)
  {
    PRINT_ERROR("Failed to get %s query engine\n", mEngineName.c_str());
    return NULL;
  }

  mQueryEngine = (*create_query_engine)();

  return mQueryEngine;
}

    
GeometryModifyEngine* CGMEngineDynamicLoader::get_gme()
{
  if(mModifyEngine)
    return mModifyEngine;

  if(mLoadAttempted == CUBIT_FALSE)
    load_engine();

  if(mLibraryHandle == CubitDynamicLoader::InvalidLibraryHandle)
  {
    PRINT_ERROR("Failed to load CGM %s engine\n", mEngineName.c_str());
    return NULL;
  }
  
  // check for incompatible CGM version
  if(check_engine_version() == CUBIT_FALSE)
    return NULL;
  
  CubitString symbol = CubitString("dynamic_") + mEngineName + CubitString("_create_gme");

  CGMEngineCreateModifyEngine create_modify_engine =
    (CGMEngineCreateModifyEngine)CubitDynamicLoader::get_symbol_address(mLibraryHandle, symbol.c_str());

  if(!create_modify_engine)
  {
    PRINT_ERROR("Failed to get %s modify engine\n", mEngineName.c_str());
    return NULL;
  }
  
  mModifyEngine = (*create_modify_engine)();

  return mModifyEngine;
}

