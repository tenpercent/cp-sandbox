int aft3dmodel (
	CGMmodel model,
	double (*fsize) (double, double, double),
	int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int maxnV, int maxnF,
	int indexshift
);
