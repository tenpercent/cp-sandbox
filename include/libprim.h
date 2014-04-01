int scg_make_sphere (int *pnV, int *pnF, double *Vertex, int *Index, double x, double meshsize);
int scg_make_paral (int *pnV, int *pnF, double *Vertex, int *Index, double x, double y, double z, double meshsize);
int scg_make_cylinder (int *pnV, int *pnF, double *Vertex, int *Index, double r, double h, double MS);

int scg_read_front (int *pnV, int *pnF, double *Vertex, int *Index, char *file, int if_color);
int scg_write_front (int nV, int nF, double *Vertex, int *Index, char *file);
int scg_write_front_gmv (int nV, int nF, double *Vertex, int *Index, char *file);
	
int scg_translate (int nV, double *Vertex, double x, double y, double z);
int scg_rotate (int nV, double *Vertex, double rot_x, double rot_y, double rot_z);

int scg_scale (int nV, double *Vertex, double s_x, double s_y, double s_z);

int scg_affine (int nV, double *Vertex, double *matrix);

