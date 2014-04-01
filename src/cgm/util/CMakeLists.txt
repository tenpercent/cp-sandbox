set(cubit_util_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "cubit_util include dir")

PROJECT(cubit_util)

INCLUDE(${cubit_util_SOURCE_DIR}/UseUtil.cmake)

SET(UTIL_SRCS
    AppUtil.cpp
    ArrayBasedContainer.cpp
    CpuTimer.cpp
    Cubit2DPoint.cpp
    CubitBox.cpp
    CubitCollection.cpp
    CubitContainer.cpp
    CubitCoordinateSystem.cpp
    CubitDynamicLoader.cpp
    CubitEntity.cpp
    CubitFile.cpp
    CubitFileFEModel.cpp
    CubitFileIOWrapper.cpp
    CubitFileMetaData.cpp
    CubitFileUtil.cpp
    CubitMatrix.cpp
    CubitMessage.cpp
    CubitObservable.cpp
    CubitObserver.cpp
    CubitPlane.cpp
    CubitStack.cpp
    CubitString.cpp
    CubitTransformMatrix.cpp
    CubitUtil.cpp
    CubitVector.cpp
    DIntArray.cpp
    DLList.cpp
    DynamicArray.cpp
    GetLongOpt.cpp
    GfxDebug.cpp
    GMem.cpp
    IntersectionTool.cpp
    MemoryBlock.cpp
    MemoryManager.cpp
    ParamCubitPlane.cpp
    PlanarParamTool.cpp
    Queue.cpp
    RandomMersenne.cpp
    SDLList.cpp
    SettingHandler.cpp
    SettingHolder.cpp
    StubProgressTool.cpp
    TDUPtr.cpp
    TextProgressTool.cpp
    ToolData.cpp
    ToolDataUser.cpp
    Tree.cpp
    TtyProgressTool.cpp
    CubitEventDefines.h
    )

SET(EXTRA_UTIL_SRCS
    AllocMemManagersList.cpp
    )

IF(TEMPLATE_DEFS_INCLUDED)
SET(TEMPLATE_SRCS
    DLIList.cpp
    KDDTree.cpp
    KDDTreeNode.cpp
    OctTreeCell.cpp
    OctTree.cpp
    PriorityQueue.cpp
    RTree.cpp
    RTreeNode.cpp
    RStarTree.cpp
    RStarTreeNode.cpp
    )
ENDIF(TEMPLATE_DEFS_INCLUDED)

IF(WIN32)
  FOREACH(var ${UTIL_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(UTIL_HEADERS ${UTIL_HEADERS} ${header})
  ENDFOREACH(var ${UTIL_SRCS})
  SET(UTIL_HEADERS ${UTIL_HEADERS} CubitEquation.hpp)
ENDIF(WIN32)

ADD_DEFINITIONS(${CUBIT_UTIL_DEFINES})


INCLUDE_DIRECTORIES(${cubit_util_BINARY_DIR})

OPTION(CUBIT_UTIL_BUILD_SHARED_LIBS "Build Cubit_Util as shared library" OFF)
IF(BUILD_SHARED_LIBS)
  SET(CUBIT_UTIL_BUILD_SHARED_LIBS ON CACHE BOOL "Build Cubit_Util as shared library" FORCE)
ENDIF(BUILD_SHARED_LIBS)

IF(CUBIT_UTIL_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_util SHARED ${UTIL_SRCS} ${UTIL_HEADERS} ${TEMPLATE_SRCS} ${EXTRA_UTIL_SRCS})

  IF(CUBIT_COPY_DIR)
    # The cubit copy dir should only be set on windows.
    STRING(REGEX REPLACE "/" "\\\\" winbin "${CUBIT_COPY_DIR}")
    IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
      ADD_CUSTOM_COMMAND(TARGET cubit_util POST_BUILD
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E copy
                           $\(TargetPath\) ${winbin}\\$\(TargetName\).dll)
    ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio")
      STRING(REGEX REPLACE "/" "\\\\" winsrc "${LIBRARY_OUTPUT_PATH}")
      ADD_CUSTOM_COMMAND(TARGET cubit_util POST_BUILD
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E copy
                           ${winsrc}\\cubit_util.dll
                           ${winlib}\\cubit_util.dll)
    ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
  ENDIF(CUBIT_COPY_DIR)
ELSE(CUBIT_UTIL_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_util ${UTIL_SRCS} ${UTIL_HEADERS} ${TEMPLATE_SRCS} ${EXTRA_UTIL_SRCS})
ENDIF(CUBIT_UTIL_BUILD_SHARED_LIBS)

CONFIGURE_FILE(${cubit_util_SOURCE_DIR}/CubitUtilConfigure.h.in
               ${cubit_util_BINARY_DIR}/CubitUtilConfigure.h
               @ONLY)


SET(CMAKE_INSTALL_BINARY_DIR "bin" CACHE PATH "Install directory for binaries")
IF(CUBIT_UTIL_BUILD_SHARED_LIBS)
  INSTALL_TARGETS(/${CMAKE_INSTALL_BINARY_DIR} cubit_util)
ENDIF(CUBIT_UTIL_BUILD_SHARED_LIBS)

