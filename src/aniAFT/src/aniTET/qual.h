#include "aft.h"
REAL aft_tetra_qual(REAL *vertex, int v1, int v2, int v3, int v4);
REAL aft_tetra_qual_mba(REAL *vertex, int v1, int v2, int v3, int v4);
REAL aft_tetra_qual_delta(REAL *vertex, int v1, int v2, int v3, int v4, REAL delta);
int print_aft_stats(REAL *vertex, int nt, int *tetra);
int print_aft_stats_mba(REAL *vertex, int nt, int *tetra);
int refine3daftfunc(int *pnV, REAL *vertex, int *pnVF, int *pnT, int *tetra);
