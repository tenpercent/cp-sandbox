#ifndef H_TETRA3_MESH3D  
#define H_TETRA3_MESH3D
/* exported  function  function */
void addBadVert(surface_mesh *pm, int v);
int isBadVert(surface_mesh *pm, int v);
void addSphereVert(surface_mesh *pm, int v);
int isSphereVert(surface_mesh *pm, int v);
#endif

