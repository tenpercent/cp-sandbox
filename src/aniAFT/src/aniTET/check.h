#ifndef LIBAFT3DCHECK_H
#define LIBAFT3DCHECK_H included

int checktopology(int nV, REAL *vertex, int nF, int *face, int nT, int *tetra, int indshift);
int check2dsurface(int *pnV, double *vertex, int *pnF, int *face, int indshift, int autofix, int verbose);

#endif //LIBAFTCHECK_H

