// lib front remesh

//int surface_refine (int *nV, double *vertex, int *nF, int *face, int *facematerial, int   maxnV, int   maxnF);
int surface_refine_(int *nV, double *vertex, int *nF, int *face, int *facematerial, int *pmaxnV, int *pmaxnF);

/* controls */
int surface_refine_setup_cf(double cf);
int surface_refine_setup_lim(double lim);
int surface_refine_setup_poly(double minabs, double min, double eps, double eps2);
int surface_refine_setup_poly_extra(double minabs, double min, double eps, double eps2, double dist);

