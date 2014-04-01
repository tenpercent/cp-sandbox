
#ifndef DYNAMIC_LOADER_HPP
#define DYNAMIC_LOADER_HPP

#if defined(__hpux)
#  include <dl.h>
#elif defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

namespace CubitDynamicLoader
{
#if defined(__hpux)
  typedef shl_t LibraryHandle;
#elif defined(_WIN32)
  typedef HMODULE LibraryHandle;
#else
  typedef void* LibraryHandle;
#endif
  
  // add a search path
  CUBIT_UTIL_EXPORT void add_search_path(const char* path);

  // check for if a library exists, only uses search paths specified
  CUBIT_UTIL_EXPORT CubitBoolean library_exists(const char* library);
  
  // open a library
  // a full path may / or may not be specified
  // use library_prefix(), library_extension() to help construct the library name
  CUBIT_UTIL_EXPORT LibraryHandle load_library(const char* library);
  
  // return the last error
  CUBIT_UTIL_EXPORT CubitString get_error();

  // close a library
  // success is returned
  CUBIT_UTIL_EXPORT CubitStatus unload_library(LibraryHandle);
  
  // get the address of a symbol in a give library
  CUBIT_UTIL_EXPORT void* get_symbol_address(LibraryHandle, const char*);
  
  // returns the library prefix for the given architecture
  CUBIT_UTIL_EXPORT const char* library_prefix();
  
  // returns the library extension for the given architecture
  CUBIT_UTIL_EXPORT const char* library_extension();

  // handle to check for load failures
  extern CUBIT_UTIL_EXPORT LibraryHandle InvalidLibraryHandle;

}

#endif

