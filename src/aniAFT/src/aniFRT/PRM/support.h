int aft3dboundary (
	int nVVert, double *VVertxyz,
	int nLine, int *LineD, int *LineP, double *LineT,
	int nSurface, int *SurfL, int *SurfI, double *SurfT,
	void (*v_u) (int, double, double*),
	int (*bounsurf) (int, double, double, double*, double*, double*),
	int (*bounline) (int, double, double*, double*),
	double (*periodicfunction) (int, int),
	double (*fsize) (double, double, double),
	int *exportCurves,
	int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int *pnEout, int *edgeout, int *edgecolor,
	int maxnV, int maxnF, int maxnE,
	int indexshift
);
int aft3dsurfmeshrefiner(int *nV, double *vertex, int *nF, int *face, int *facematerial, int maxnV, int maxnF, int indexshift);

