#include "cgmwrap.h"
#include "common.h"
#include "user3c.h"

static surface_mesh *pm;

int bounSurfCGM(int i, double u, double v, double *x, double *y, double *z) {
	double xyz[3], uv[2];
	(void) i;
	uv[0] = u;  uv[1] = v;
	CGM_GetFaceCoordsFromUV(pm->cgmface, uv, xyz);
	*x = xyz[0];  *y = xyz[1];  *z = xyz[2];
	return 1; 
} /*bounSurf*/

double periodicCGM(int dummy, int dir) {
	return ((dir)?CGM_FaceIsPeriodicV:CGM_FaceIsPeriodicU)(pm->cgmface);
	(void) dummy;
}

int SurfNormCGM(double x, double y, double z, double *N) {
	double X[3];
	X[0] = x,  X[1] = y,  X[2] = z;
	CGM_GetFaceNormal(pm->cgmface, X, N);
	return 0;
}
int user3initCGM(surface_mesh *lpm) {
    pm = lpm;
    pm->periodicfunction = periodicCGM;
    pm->cgmfunction = bounSurfCGM;
    pm->cgmsurfNormal = SurfNormCGM;
    return 0;
}
