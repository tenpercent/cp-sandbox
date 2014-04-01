#ifndef H_TRIA32_MESH3D
#define H_TRIA32_MESH3D


#include"struct3.h"
#include"tree32.h"


/* exported  function  function */
void  reUV( double *u, double *v, double u0, double v0 );
void rePlane(surface_mesh *pm, int vv);
void makeTria(surface_mesh *pm);


#endif

