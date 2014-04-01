#ifndef LIBAFT2DKERNEL_H
#define LIBAFT2DKERNEL_H included

int mesh2daft  (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
int mesh2daftss(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss);
int mesh2daftsslim(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss, REAL lim);
int mesh2daftfull(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss, REAL lim, int *pnC, int *color);
    
#endif //LIBAFT2DKERNEL_H
