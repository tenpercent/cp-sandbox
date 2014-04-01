#ifndef H_STRUCT3_MESH3D
#define H_STRUCT3_MESH3D


#include<ctype.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>



#define  S_StrucVert3     sizeof(StrucVert3)
#define  S_StrucNeigbor   sizeof(StrucNeigbor)
#define  S_StrucEdge3     sizeof(StrucEdge3)
#define  S_StrucFace3     sizeof(StrucFace3)


/* exported  function  function */
void init(surface_mesh *pm);
void addPoint(surface_mesh *pm, double x, double y, double z);


#endif

