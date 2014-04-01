#ifndef REAL
#define REAL double
#endif

/* mesh 2D aft */
//int mesh_2d_aft_opts          (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf, int loop, int needcheck, int indexshift);
//int mesh_2d_aft               (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
//int mesh_2d_aft_loop          (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
//int mesh_2d_aft_check         (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
//int mesh_2d_aft_loop_check    (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
//int mesh_2d_aft_auto          (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT);
//int mesh_2d_aft_cf_auto       (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf);
//int mesh_2d_aft_cf            (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf);
//int mesh_2d_aft_cf_loop       (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf);
//int mesh_2d_aft_cf_check      (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf);
//int mesh_2d_aft_cf_loop_check (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf);
int mesh_2d_aft_              (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT);
int mesh_2d_aft_loop_         (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT);
int mesh_2d_aft_check_        (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT);
int mesh_2d_aft_loop_check_   (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT);
int mesh_2d_aft_auto_         (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT);
int mesh_2d_aft_cf_           (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf);
int mesh_2d_aft_cf_loop_      (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf);
int mesh_2d_aft_cf_check_     (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf);
int mesh_2d_aft_cf_loop_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf);
int mesh_2d_aft_cf_lim_loop_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf, REAL *plim);
int mesh_2d_aft_cf_auto_      (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf);
int mesh_2d_aft_full_         (int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf, int *pnC, int *color);

/* mesh 3D aft */
//int mesh_3d_aft_opts (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, REAL cf, int needcheck, int indexshift, double (*f)(double, double, double));
//int mesh_3d_aft      (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT);
//int mesh_3d_aft_cf   (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, REAL cf);
//int mesh_3d_aft_func (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, double (*f)(double, double, double));
int mesh_3d_aft_     (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT);
int mesh_3d_aft_cf_  (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT, REAL *pcf);
int mesh_3d_aft_func_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT, double (*f)(double, double, double));

/* mesh 3D ugly */
int mesh_3d_ugly (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT);
int mesh_3d_ugly_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT);

/* read/write */
int read_front     (char *fn, int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int maxnV, int maxnF, int zeroindex, int invert, int verbose);
int write_front    (char *fn, int   nV, double *vertex, int   nF, int *face, int *facematerial);
int write_front_gmv(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial);
int write_mesh_gmv (char *fn, int   nV, double *vertex, int   nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial);
int write_mesh_gmv_qual(char *fn, int   nV, double *vertex, int   nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial);
int write_mesh     (char *fn, int   nV, double *vertex, int   nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial);

/* checks */
//int check_mesh_topology (int   nV, REAL *vertex, int   nF, int *face, int   nT, int *tetra);
int check_mesh_topology_(int *pnV, REAL *vertex, int *pnF, int *face, int *pnT, int *tetra);
int check_surface_topology     (int   nV, REAL *vertex, int   nF, int *face);
int check_surface_topology_    (int *pnV, REAL *vertex, int *pnF, int *face);
int check_surface_topology_fix (int *pnV, REAL *vertex, int *pnF, int *face);
int check_surface_topology_fix_(int *pnV, REAL *vertex, int *pnF, int *face);

