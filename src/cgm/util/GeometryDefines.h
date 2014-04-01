/* Class:       GeometryDefines 
 * Description: GeometryDefines - all global geometric definitions for CUBIT.
 * Owner:       David White
 * Checked by:
 * Version:
 */

#ifndef GEOMETRY_DEFINES_HPP
#define GEOMETRY_DEFINES_HPP

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

/* CUBIT_RESABS - Values less than CUBIT_RESABS are considered to be zero. */
#ifdef __cplusplus
const double CUBIT_RESABS = 1.0E-12;
#else
#define CUBIT_RESABS 1.0E-12
#endif

/* GEOMETRY_RESABS - If the distance between two points is less */
/* than GEOMETRY_RESABS the points are considered to be identical. */
#ifdef __cplusplus
const double GEOMETRY_RESABS = 1.0E-6;
#else
#define GEOMETRY_RESABS 1.0E-6
#endif

/* Types of solid modeler engines. */
#ifdef ACIS
	#undef ACIS
#endif
enum SolidModelerType 
{ 
   NOT_A_SOLID_MODELER,
   OCC,
   ACIS,
   PROE_GEOM,                  /* Normal Pro/E model */
   PROE_FEM_MESH_SOLID,        /* Pro/Mesh models... */
   PROE_FEM_MESH_SHELL,
   PROE_FEM_MESH_MIXED,
   PROE_FEM_MESH_BOUNDARY,
   PROE_FEM_MESH_QUILT
} ;

/* Types of underlying geometric representations. */
enum GeometricRepresentationType 
{
   NONE,
   SOLID_MODEL,
   FACETTED,
   MESH_BASED,
   COMPOSITE
};

/* Types of geometry for the various GeometryEntities. */
enum GeometryType
{
    /* Point types */
  UNDEFINED_POINT_TYPE,
  
    /* Curve types */
  ARC_CURVE_TYPE,
  ELLIPSE_CURVE_TYPE,
  PARABOLA_CURVE_TYPE,
  SEGMENTED_CURVE_TYPE, /*curve defined by a chain of linear segments. */
  SPLINE_CURVE_TYPE,
  STRAIGHT_CURVE_TYPE,
  POINT_CURVE_TYPE,
  BSPLINE_CURVE_TYPE,   //OCC curve types 
  HYPERBOLA_CURVE_TYPE,
  UNDEFINED_CURVE_TYPE,
  
    /* Surface types */
  CONE_SURFACE_TYPE,
  PLANE_SURFACE_TYPE,
  SPHERE_SURFACE_TYPE,
  SPLINE_SURFACE_TYPE,
  TORUS_SURFACE_TYPE,
  BEST_FIT_SURFACE_TYPE,
  FACET_SURFACE_TYPE,
  BSPLINE_SURFACE_TYPE,     //OCC surface type
  REVOLUTION_SURFACE_TYPE,  //OCC surface type
  EXTRUSION_SURFACE_TYPE,   //OCC surface type
  OFFSET_SURFACE_TYPE,      //OCC surface type
  UNDEFINED_SURFACE_TYPE,
  
    /* Lump types */
  UNDEFINED_LUMP_TYPE
};

#endif

