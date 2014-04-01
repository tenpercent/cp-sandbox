#ifndef H_TREE32_MESH3D
#define H_TREE32_MESH3D


#include"struct3.h"


/* exported  function */
PStrucEdge3 addFace32(surface_mesh *pm, int v1, int v2 );
void remFace32(surface_mesh *pm, PStrucEdge3  face);
double nearest32(surface_mesh *pm, int *vert, double x, double y, double z );
void vicinityFaces32(surface_mesh *pm, double x, double y, double z, double size);
void prepTree32(surface_mesh *pm, int nnE);


#endif


