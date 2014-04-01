
set(cubit_geom_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "cubit_geom include dir")

PROJECT(cubit_geom)

SUBDIRS(OCC)

#SUBDIRS(Cholla facetbool facet virtual)

# build cgm with acis support, default to on
#OPTION(CGM_ACIS "Build CGM with ACIS support" ON)
#IF(CGM_ACIS)
#SUBDIRS(ACIS)
#ENDIF(CGM_ACIS)

# build cgm with granite support, default to off
#OPTION(CGM_GRANITE "Build CGM with Granite support" OFF)
#IF(CGM_GRANITE)
#  SUBDIRS(granite)
#ENDIF(CGM_GRANITE)

# build cgm with OpenCASCADE support, default to on
#OPTION(CGM_OCC "Build CGM with OpenCASCADE support" ON)
#IF(CGM_OCC)
#  SUBDIRS(OCC)
#ENDIF(CGM_OCC)

INCLUDE(${cubit_geom_SOURCE_DIR}/../util/UseUtil.cmake)

SET(GEOM_SRCS
    AnalyticGeometryTool.cpp
    BasicTopologyEntity.cpp
    Body.cpp
    BodySM.cpp
    BoundingBoxTool.cpp
    BridgeManager.cpp
    CAActuateSet.cpp
    CADeferredAttrib.cpp
    CAEntityColor.cpp
    CAEntityId.cpp
    CAEntityName.cpp
    CAGroup.cpp
    CAMergePartner.cpp
    CAMergeStatus.cpp
    CASourceFeature.cpp
    CAUniqueId.cpp
    CGMApp.cpp
    CGMEngineDynamicLoader.cpp
    Chain.cpp
    CoEdge.cpp
    CoEdgeSM.cpp
    CoFace.cpp
    CollectionEntity.cpp
    CoVertex.cpp
    CoVolume.cpp
    CubitAttrib.cpp
    CubitAttribManager.cpp
    CubitAttribUser.cpp
    CubitPolygon.cpp
    CubitSimpleAttrib.cpp
    Curve.cpp
    CurveOverlapFacet.cpp
    CurveSM.cpp
    CylinderEvaluator.cpp
    DAG.cpp
    DagDrawingTool.cpp
    GeomDataObserver.cpp
    GeometryEntity.cpp
    GeometryFeatureEngine.cpp
    GeometryFeatureTool.cpp
    GeometryHealerEngine.cpp
    GeometryHealerTool.cpp
    GeometryModifyEngine.cpp
    GeometryModifyTool.cpp
    GeometryQueryEngine.cpp
    GeometryQueryTool.cpp
    GeometryUtil.cpp
    GeomMeasureTool.cpp
    GfxPreview.cpp
    GroupingEntity.cpp
    GSaveOpen.cpp
    Loop.cpp
    LoopSM.cpp
    Lump.cpp
    LumpSM.cpp
    MedialTool2D.cpp
    MedialTool3D.cpp
    MergeToolAssistant.cpp
    MergeTool.cpp
    MidPlaneTool.cpp
    ModelEntity.cpp
    ModelQueryEngine.cpp
    OffsetSplitTool.cpp
    OldUnmergeCode.cpp
    Point.cpp
    PointSM.cpp
    RefCollection.cpp
    RefEdge.cpp
    RefEntity.cpp
    RefEntityFactory.cpp
    RefEntityName.cpp
    RefFace.cpp
    RefGroup.cpp
    RefVertex.cpp
    RefVolume.cpp
    SenseEntity.cpp
    Shell.cpp
    ShellSM.cpp
    SphereEvaluator.cpp
    SplitSurfaceTool.cpp
    Surface.cpp
    SurfaceOverlapFacet.cpp
    SurfaceOverlapTool.cpp
    SurfParamTool.cpp
    SurfaceSM.cpp
    TDCAGE.cpp
    TDSourceFeature.cpp
    TDSplitSurface.cpp
    TDSurfaceOverlap.cpp
    TDUniqueId.cpp
    TopologyBridge.cpp
    TopologyEntity.cpp
    )

SET(EXTRA_GEOM_SRCS
    AllocMemManagersGeom.cpp
    )


IF(TEMPLATE_DEFS_INCLUDED)
SET(TEMPLATE_SRCS
    BoundaryConstrainTool.cpp
    )
ENDIF(TEMPLATE_DEFS_INCLUDED)

IF(WIN32)
  FOREACH(var ${GEOM_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(GEOM_HEADERS ${GEOM_HEADERS} ${header})
  ENDFOREACH(var ${GEOM_SRCS})
ENDIF(WIN32)

# don't add an include path to other parts of cgm
SET(THIS_INCLUDES 
          ${cubit_util_SOURCE_DIR} 
          ${cubit_util_BINARY_DIR} 
          ${cubit_geom_SOURCE_DIR} 
          ${cubit_geom_BINARY_DIR})

SET_DIRECTORY_PROPERTIES(PROPERTIES INCLUDE_DIRECTORIES "${THIS_INCLUDES}")

OPTION(CUBIT_GEOM_BUILD_SHARED_LIBS "Build cubit_geom as a shared library" OFF)

IF(BUILD_SHARED_LIBS)
  SET(CUBIT_GEOM_BUILD_SHARED_LIBS ON CACHE BOOL "Build cubit_geom as shared library" FORCE)
ENDIF(BUILD_SHARED_LIBS)

IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  # util must be shared if geom is.
  SET(CUBIT_UTIL_BUILD_SHARED_LIBS ON CACHE BOOL "Build cubit util as a shared library" FORCE)
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)


CONFIGURE_FILE(${cubit_geom_SOURCE_DIR}/CubitGeomConfigure.h.in
               ${cubit_geom_BINARY_DIR}/CubitGeomConfigure.h
               @ONLY)

IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_geom SHARED ${GEOM_SRCS} ${EXTRA_GEOM_SRCS} ${GEOM_HEADERS} ${TEMPLATE_SRCS})

  IF(CUBIT_COPY_DIR)
    # The cubit copy dir should only be set on windows.
    STRING(REGEX REPLACE "/" "\\\\" winbin "${CUBIT_COPY_DIR}")
    IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
      ADD_CUSTOM_COMMAND(TARGET cubit_geom POST_BUILD
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E copy
                           $\(TargetPath\) ${winbin}\\$\(TargetName\).dll)
    ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio")
      STRING(REGEX REPLACE "/" "\\\\" winsrc "${LIBRARY_OUTPUT_PATH}")
      ADD_CUSTOM_COMMAND(TARGET cubit_geom POST_BUILD
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E copy
                           ${winsrc}\\cubit_geom.dll
                           ${winlib}\\cubit_geom.dll)
    ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
  ENDIF(CUBIT_COPY_DIR)
ELSE(CUBIT_GEOM_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_geom ${GEOM_SRCS} ${EXTRA_GEOM_SRCS} ${GEOM_HEADERS} ${TEMPLATE_SRCS})
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)

TARGET_LINK_LIBRARIES(cubit_geom cubit_util)

IF(UNIX)
ADD_CUSTOM_TARGET(switch_acis 
                  rm -f *Acis*.o &&
                  rm -f attrib*.o &&
                  rm -f *ACIS*.o &&
                  rm -f GeometryCommands.o &&
                  rm -f Feature*.o &&
                  rm -f *App.o &&
                  cd ACIS && ${CMAKE_MAKE_PROGRAM} switch_acis)

ENDIF(UNIX)

SET(CMAKE_INSTALL_BINARY_DIR "bin" CACHE PATH "Install directory for binaries")
IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  INSTALL_TARGETS(/${CMAKE_INSTALL_BINARY_DIR} cubit_geom)
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)

