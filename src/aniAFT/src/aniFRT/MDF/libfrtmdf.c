#include <stdio.h>
#include <errno.h>
#include "surfrf.h"
#include "libfrtmdf.h"
#include "support.h"

int surface_refine(int *nV, double *vertex, int *nF, int *face, int *facematerial, int maxnV, int maxnF) {
    return aft3dsurfmeshrefiner(nV, vertex, nF, face, facematerial, maxnV, maxnF, 0);
}
int surface_refine_(int *nV, double *vertex, int *nF, int *face, int *facematerial, int *pmaxnV, int *pmaxnF) {
    return aft3dsurfmeshrefiner(nV, vertex, nF, face, facematerial, *pmaxnV, *pmaxnF, 1);
}

int surface_refine_setup_cf(double cf) {
    return libaft_internal_surfmeshrefiner_ss_setup(cf);
}
int surface_refine_setup_lim(double lim) {
    return libaft_internal_surfmeshrefiner_lim_setup(lim);
}
int surface_refine_setup_poly(double minabs, double min, double eps, double eps2/*, double dist*/) {
    return libaft_internal_surfmeshrefiner_poly_setup(minabs, min, eps, eps2/*, dist*/);
}
int surface_refine_setup_poly_extra(double minabs, double min, double eps, double eps2, double dist) {
    return libaft_internal_surfmeshrefiner_poly_setup_extra(minabs, min, eps, eps2, dist);
}
