PROJECT(cubit_occ)

set(cubit_occ_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "cubit_occ include dir")

SET(OCC_SRCS
    OCCShapeAttributeSet.cpp
    OCCAttribSet.cpp
    OCCBody.cpp
    OCCCoEdge.cpp
    OCCCoFace.cpp
    OCCCurve.cpp
    OCCLoop.cpp
    OCCLump.cpp
    OCCModifyEngine.cpp
    OCCPoint.cpp
    OCCQueryEngine.cpp
    OCCShell.cpp
    OCCSurface.cpp
    OCCDrawTool.cpp
   )

SET(OCC_HDRS
    OCCShapeAttributeSet.hpp
    OCCAttribSet.hpp
    OCCBody.hpp
    OCCCoEdge.hpp
    OCCCoFace.hpp
    OCCCurve.hpp
    OCCGeometryCreator.hpp
    OCCLoop.hpp
    OCCLump.hpp
    OCCModifyEngine.hpp
    OCCPoint.hpp
    OCCQueryEngine.hpp
    OCCShell.hpp
    OCCSurface.hpp
    OCCDrawTool.hpp
    )

SET(OCC_LIBS
     TKernel
     TKMath
     TKXSBase
     TKernel
     TKShHealing
     TKBRep
     TKG3d
     TKBO
     TKFeat
     TKFillet
     TKG2d
     TKOffset
     TKBool
     TKGeomAlgo
     TKGeomBase
     TKPrim
     TKTopAlgo
     TKHLR
     TKIGES
     TKMesh
     TKPrim
     TKSTEP209
     TKSTEPAttr
     TKSTEPBase
     TKSTEP
     TKSTL
     TKTopAlgo
     TKXSBase
     TKLCAF
     TKBinL
    )

if (WIN32)
	add_definitions(-DWNT)
endif (WIN32)
FIND_PATH( OCC_INCLUDES Standard_Version.hxx PATHS ${OCC_DIR} ENV CASROOT  PATH_SUFFIXES inc )
SET(OCC_DIR $ENV{CASROOT} CACHE PATH "Location where OpenCASCADE in installed")
if (WIN32)
	if (WIN64)
		FIND_LIBRARY( OCC_LIBRARY TKernel ${OCC_DIR}/win64/lib )
	else (WIN64)
		FIND_LIBRARY( OCC_LIBRARY TKernel ${OCC_DIR}/win32/lib )
	endif (WIN64)
    if (WIN64)
        set(OCC_BIN_DIR "${OCC_DIR}/win64/bin" CACHE PATH "OCC DLLs Directory")
    else (WIN64)
        set(OCC_BIN_DIR "${OCC_DIR}/win32/bin" CACHE PATH "OCC DLLs Directory")
    endif (WIN64)
else (WIN32)
	FIND_LIBRARY( OCC_LIBRARY TKernel PATHS ${OCC_DIR} ENV CASROOT  PATH_SUFFIXES lib )
endif (WIN32)
get_filename_component(OCC_LIBRARY_DIR ${OCC_LIBRARY} PATH CACHE)
message(STATUS "Using Open CASCADE Library Directory: ${OCC_LIBRARY_DIR}")
link_directories(${OCC_LIBRARY_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif()
ADD_LIBRARY( cubit_OCC ${OCC_SRCS} ${OCC_HDRS} )
INCLUDE_DIRECTORIES( ${OCC_INCLUDES} )
TARGET_LINK_LIBRARIES( cubit_OCC ${OCC_LIBS} )
if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
	SET( OCC64 "-D_OCC64")
endif(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
SET_SOURCE_FILES_PROPERTIES( ${OCC_SRCS} COMPILE_FLAGS "-DHAVE_IOSTREAM -DHAVE_IOMANIP -DHAVE_FSTREAM -DHAVE_LIMITS_H ${OCC64}" )
