#ifndef H_TREE3_MESH3D
#define H_TREE3_MESH3D

#include "common.h"

/* exported functions */
/*double distance(double x,double y,double z,double xc,double yc,double zc);*/
double distanceS(double x,double y,double z,double xc,double yc,double zc);
PStrucFace3 addFace(surface_mesh *pm, int v1, int v2, int v3, int twin, int color);

#endif


