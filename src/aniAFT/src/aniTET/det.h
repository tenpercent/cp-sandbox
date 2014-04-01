#include "aft.h"
void detinit(void);

REAL orient2d(REAL *pa, REAL *pb, REAL *pc);
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL orient4d(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

REAL det2i3(REAL *vertex, int v1, int v2, int v3);
int idet2i3(REAL *vertex, int v1, int v2, int v3);
REAL det3i4(REAL *vertex, int v1, int v2, int v3, int v4);
int idet3i4(REAL *vertex, int v1, int v2, int v3, int v4);
REAL det4i5(REAL *vertex, int v1, int v2, int v3, int v4, int v5);
int idet4i5(REAL *vertex, int v1, int v2, int v3, int v4, int v5);
REAL det3(REAL *v1, REAL *v2, REAL *v3);
void normvec2(REAL *v1, REAL *v2, REAL *v);
void normvec3i(REAL *vertex, int v1, int v2, int v3, REAL *v);
