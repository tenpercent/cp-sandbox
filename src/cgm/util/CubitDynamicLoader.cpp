
// We have 4 different implementations here
// Windows, Unix, MacOSX, HP-UX

#include "CubitDynamicLoader.hpp"

#include <vector>
#include <sys/stat.h>

#include "CubitMessage.hpp"
#include "CubitString.hpp"

static std::vector<CubitString> gSearchPaths;

void CubitDynamicLoader::add_search_path(const char* path)
{
  gSearchPaths.push_back(CubitString(path));
}

static bool absolute_path(const char* path)
{
#ifdef WIN32
  // either a '\' or a 'x:' drive letter means absolute path
  return path[0] == '\\' || (path[0] != '\0') ? (path[1] == ':') : 0;
#else
  // a '/' is always an absolute path
  return path[0] == '/';
#endif
}

CubitBoolean CubitDynamicLoader::library_exists(const char* library)
{
  struct stat buf;

  // if absolute path, test it directly
  if(absolute_path(library))
  {
    return stat(library, &buf) == 0 ? CUBIT_TRUE : CUBIT_FALSE;
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<(int)gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i] + CubitString("/") + library;
      if(stat(path.c_str(), &buf) == 0)
        return CUBIT_TRUE;
    }
  }
  
  // one more final attempt
  if (stat(library, &buf))
    return CUBIT_TRUE;

  return CUBIT_FALSE;
}

#ifdef _WIN32

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  // disable message boxes about failure to find dlls
  //UINT old_error_mode = SetErrorMode(SEM_NOOPENFILEERRORBOX | SEM_FAILCRITICALERRORS);
	UINT old_error_mode = SetErrorMode(SEM_FAILCRITICALERRORS);

  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    handle = LoadLibrary(lib);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("\\")+lib;
      handle = LoadLibrary(path.c_str());
      if(handle != InvalidLibraryHandle)
      {
        SetErrorMode(old_error_mode);
        return handle;
      }
    }
  }
  
  // one more final attempt
  handle = LoadLibrary(lib);

  SetErrorMode(old_error_mode);

  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  LPVOID lpMsgBuf;
  DWORD dw = GetLastError(); 

  FormatMessage(
    FORMAT_MESSAGE_ALLOCATE_BUFFER | 
    FORMAT_MESSAGE_FROM_SYSTEM |
    FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL,
    dw,
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
    (LPTSTR) &lpMsgBuf,
    0, NULL );

  CubitString error = (LPTSTR)lpMsgBuf;
  LocalFree(lpMsgBuf);

  return error;
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return FreeLibrary(lib) == FALSE ? CUBIT_FAILURE : CUBIT_SUCCESS;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  return GetProcAddress(lib, symbol);
}

const char* CubitDynamicLoader::library_prefix()
{
  return "";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".dll";
}

#elif defined(__hpux)

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    handle = shl_load(lib, BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("/")+lib;
      handle = shl_load(path.c_str(), BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
      if(handle != InvalidLibraryHandle)
        return handle;
    }
  }
  
  // one more final attempt
  handle = shl_load(lib, BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  return CubitString();
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return shl_unload(lib) == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  void* address;
  return shl_findsym(&lib, symbol, TYPE_PROCEDURE, &address) < 0 ? NULL : address;
}

const char* CubitDynamicLoader::library_prefix()
{
  return "lib";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".sl";
}

#elif defined(__APPLE__)

#include <mach-o/dyld.h>

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  NSObjectFileImageReturnCode rc;
  NSObjectFileImage image;
  
  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    rc = NSCreateObjectFileImageFromFile(lib, &image);
    handle = NSLinkModule(image, lib, TRUE);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("/")+lib;
      rc = NSCreateObjectFileImageFromFile(path.c_str(), &image);
      handle = NSLinkModule(image, lib, TRUE);
      if(handle != InvalidLibraryHandle)
        return handle;
    }
  }
  
  // one more final attempt
  rc = NSCreateObjectFileImageFromFile(lib, &image);
  handle = NSLinkModule(image, lib, TRUE);
  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  return CubitString();
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return NSUnLinkModule(static_cast<NSModule>(lib), FALSE) == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle, const char* sym)
{
  void *result = 0;
  // global 'C' symbols names are preceded with an underscore '_'
  char *_sym = new char[ strlen(sym) + 2 ];
  strcpy( _sym + 1, sym );
  _sym[0] = '_';
  if( NSIsSymbolNameDefined(_sym) )
  {
    NSSymbol symbol = NSLookupAndBindSymbol(_sym);
    if(symbol)
    {
      result = NSAddressOfSymbol(symbol);
    }
  }
  delete[] _sym;
  return result;
}

const char* CubitDynamicLoader::library_prefix()
{
  return "lib";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".so";
}

#else

#include <dlfcn.h>

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    handle = dlopen(lib, RTLD_LAZY);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<(int)gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("/")+lib;
      handle = dlopen(path.c_str(), RTLD_LAZY);
      if(handle != InvalidLibraryHandle)
        return handle;
    }
  }
  
  // one more final attempt
  handle = dlopen(lib, RTLD_LAZY );
  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  return dlerror();
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return dlclose(lib) == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  return dlsym(lib, symbol);
}

const char* CubitDynamicLoader::library_prefix()
{
  return "lib";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".so";
}
#endif

