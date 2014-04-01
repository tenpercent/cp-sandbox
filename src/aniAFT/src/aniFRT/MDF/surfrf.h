int libaft_internal_surfmeshrefiner(int *nVRT, double *vrt, int *nTRI, int *tri, int *material, int Vmax, int Tmax);
int libaft_internal_surfmeshrefiner_ss_setup(double ss);
int libaft_internal_surfmeshrefiner_lim_setup(double lim);
int libaft_internal_surfmeshrefiner_poly_setup(double minabs, double min, double eps, double eps2/*, double dist*/);
int libaft_internal_surfmeshrefiner_poly_setup_extra(double minabs, double min, double eps, double eps2, double dist);
