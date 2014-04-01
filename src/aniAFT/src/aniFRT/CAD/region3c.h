#ifndef H_REGION3C_MESH3D
#define H_REGION3C_MESH3D

#include "cgmwrap.h"
#include "region3.h"

void makeAFLine(surface_mesh *pm, int *vert, int nVert, int bInverse);

/* exported  functions */
void initAFSM_(surface_mesh *pm, CGMmodel model);

#endif

