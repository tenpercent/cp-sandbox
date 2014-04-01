# Ani3D
cmake_minimum_required(VERSION 2.6)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

option(ANI3D_USE_AFT "Enable Tetrahedral Mesh Generator (aniAFT) support in Ani3D" ON)
if (ANI3D_USE_AFT)
  option(ANI3D_USE_FRT_SCG "Enable Constructive Solid Geometry support in Ani3D" ON)
  option(ANI3D_USE_FRT_MDF "Enable Discrete Surface Remeshing support in Ani3D" ON)
  option(ANI3D_USE_FRT_PRM "Enable Parametric Surface Meshing support in Ani3D" ON)
  if (ANI3D_USE_FRT_PRM)
    option(ANI3D_USE_CGM_OCC "Enable Open CASCADE Technology and CGM support in Ani3D" ON)
  endif (ANI3D_USE_FRT_PRM)
endif (ANI3D_USE_AFT)
option(ANI3D_USE_MBA "Enable Metric Based Adaptation (aniMBA) support in Ani3D" ON)
option(ANI3D_USE_FEM "Enable Finite Elements (aniFEM) support in Ani3D" ON)
option(ANI3D_USE_ILU "Enable Incomplete LU Factorization (aniILU) support in Ani3D" ON)
option(ANI3D_USE_INB "Enable Inexact Newton-Krylov Backtracking (aniINB) support in Ani3D" ON)
option(ANI3D_USE_RCB "Enable Refining/Coarsening tool (aniRCB) support in Ani3D" ON)
option(ANI3D_USE_LMR "Enable Local Metric Recovery (aniLMR) support in Ani3D" ON)
option(ANI3D_USE_VIEW "Enable Visualization (aniVIEW) support in Ani3D" ON)
option(ANI3D_BUILD_TUTORIALS "Build Tutorials for Ani3D" ON)

option(ANI3D_BUILD_BIOIMPEDANCE "Build BioImpedance project" ON)
if (ANI3D_BUILD_BIOIMPEDANCE)
  option(ANI3D_BUILD_BIOIMPEDANCE_DIFF "Build BioImpedance Diff project" ON)
  option(ANI3D_BUILD_BIOIMPEDANCE_MESH "Build BioImpedance Mesh project" ON)
endif (ANI3D_BUILD_BIOIMPEDANCE)

option(ANI3D_USE_SYSTEM_BLAS_LAPACK "Use system BLAS/LAPACK libraries" ON)

project(Ani3D C CXX Fortran)

include(FortranCInterface)

if (MSVC)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)
IF(UNIX)
   SET(MATH_LIB -lm)
ELSE(UNIX)
   SET(MATH_LIB)
ENDIF(UNIX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (ANI3D_USE_SYSTEM_BLAS_LAPACK)
    find_package(BLAS)
    find_package(LAPACK)
endif (ANI3D_USE_SYSTEM_BLAS_LAPACK)

if (NOT BLAS_FOUND)
    set(BLAS_LIBRARIES blas_ani3d)
endif (NOT BLAS_FOUND)

if (NOT LAPACK_FOUND)
    set(LAPACK_LIBRARIES lapack_ani3d)
endif (NOT LAPACK_FOUND)

add_subdirectory(src)
