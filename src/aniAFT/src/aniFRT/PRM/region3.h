#ifndef H_REGION3_MESH3D
#define H_REGION3_MESH3D

#ifndef CGMWRAP_H
#define CGMface void*
#endif

typedef struct{
	int     iSurf,iV_U;
	double  tBegin,tEnd,*u,*v;
}  StrucSub;


typedef struct{
	int      nVert,*vert;
	int      vBegin,vEnd;
	int      nSub;
	StrucSub *sub;
}  StrucLine3;


typedef struct{
	int     nLine,*line,*inverse;
	int     iSurf,iNorm,iLabel,bCut;
	double  uMin,uMax,vMin,vMax;
}  StrucSurface;



#define  MAX1   50000


/* exported  functions */
int  bounSurf0(surface_mesh *pm,  double u, double v, double *x, double *y, double *z );
void V_U0(surface_mesh *pm,  double u, double *v );
double sizeFace(surface_mesh *pm,  double x, double y, double z );
int addTria(surface_mesh *pm, int v1, int v2, int v3);

void initAFS_ (
	surface_mesh *pm,
	int *pnVVert, double *VVertxyz,
	int *pnLine, int *LineD, int *LineP, double *LineT,
	int *pnSurface, int *SurfL, int *SurfI, double *SurfT
);
#endif

