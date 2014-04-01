#ifndef LIBAFT3DKERNEL_H
#define LIBAFT3DKERNEL_H included

int mesh3daft(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT);
int mesh3daftss(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, REAL ss);
int mesh3daftuser(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, double (*f)(double, double, double));

#endif //LIBAFT3DKERNEL_H
