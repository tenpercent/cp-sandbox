#include "common.h"
#include "region3.h"
#include "user3.h"

int bounSurf(surface_mesh *pm, int i, double u, double v, double *x, double *y, double *z) {
    if (i==-1) {
	if (pm->cgmfunction)  return pm->cgmfunction(i, u, v, x, y, z);
	else  return -1;
    }
    if (i==0)  return bounSurf0(pm, u, v, x, y, z);
    if (pm->bounsurf)  return pm->bounsurf(i, u, v, x, y, z);
    printf("Surface parametrization functions are not defined.\n");
    return -1;
}
double periodic(surface_mesh *pm, int dir) {
    if (pm->iSurf == 0)  return 0.0;
    if (pm->periodicfunction)  return pm->periodicfunction(pm->iSurf, dir);
    else  return 0.0;
}

void V_U(surface_mesh *pm, int i, double u, double *v) {
    if (i==0)  return V_U0(pm, u, v);
    if (pm->v_u)  return pm->v_u(i, u, v);
    printf("Curve parametrization functions are not defined.\n");
}

int bounLine(surface_mesh *pm, int i, double t, double *u, double *v) {
    if (i==0) {
	*u = t;
	*v = 0.0;
    } else if (pm->bounline) {
	return pm->bounline(i, t, u, v);
    } else if (pm->v_u) {
	*u = t;
	V_U(pm, i, t, v);
    } else {
	printf("Curve parametrization functions are not defined.\n");
	return 0;
    }
    return 1;
}

double userSizeFace(surface_mesh *pm, double x, double y, double z) {
    if (pm->fsize)  return pm->fsize(x,y,z);
    return pm->S0;
} /*userSizeFace*/
