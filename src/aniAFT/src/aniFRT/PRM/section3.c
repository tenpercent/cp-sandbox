#include "common.h"
#include "struct3.h"
#include "error3.h"
#include "section3.h"
#include "tetra3.h"
#include "tree3.h"

void makeVector(surface_mesh *pm, PStrucVert3 v, int v1, int v2) {
	v->x = pm->vert[v1].x-pm->vert[v2].x;
	v->y = pm->vert[v1].y-pm->vert[v2].y;
   v->z = pm->vert[v1].z-pm->vert[v2].z;
   return;
} /*makeVector*/


int makeNormal (PStrucVert3 v, PStrucVert3 v1, PStrucVert3 v2) {
   double p;

   v->x = v1->y*v2->z-v2->y*v1->z;
   v->y = v1->z*v2->x-v2->z*v1->x;
   v->z = v1->x*v2->y-v2->x*v1->y;
   p = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
   if (p != 0.0) {
      v->x /= p; v->y /= p; v->z /= p;
      return  1;
   }
   return  0;
} /*makeNormal*/





