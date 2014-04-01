#include <stdlib.h>
#include <stdio.h>
#include "Ani3D_FC.h"

int ANI_Metric_Eucl(double*, double*, double*, double*);
/*int ani_metric_eucl_();*/

int idet3i4(double *vertex, int v1, int v2, int v3, int v4);

int mbaFixShape(int *nv, int *nvmax, int *nb, int *nbmax, int *nt, int *ntmax, double *vrt, int *bnd, int *tet, int *labelF, int *material,
	int *nvfix, int *nbfix, int *ntfix, int *ivfix, int *ibfix, int *itfix, int *flagAuto, int *status,
	int *MaxSkipE, int *MaxQItr, int (*MetricFunction_user)(double*, double*, double*, double*), double *Quality, double *rQuality,
	int *MaxWr, int *MaxWi, double *rW, int *iW, int *iPrint, int *iERR);

int fixshape(int *pnv, double *vrt, int *pnt, int *tet, int *material, int *pnb, int *bnd, int *labelF,
	int nvfix, int ntfix, int nbfix, int nvmax, int ntmax, int nbmax) {
    int nv=*pnv, nt=*pnt, nb=*pnb;
    int MaxWr = 20000000, MaxWi = 60000000;
    int *iW;
    double *rW;
    int MaxSkipE = 300, MaxQItr = 500000;
    double Quality = 0.4, rQuality;
    int flagAuto;
    int iPrint, status, iERR;
    int *ivfix, *ibfix, *itfix;
    int i, r, d;
    int mmt, mmb = 0;

    if (!nt)  return 0;

    for (i=1, mmt = material[0]; i<nt; i++)  if (mmt > material[i])  mmt = material[i];
    for (i=0; i<nt; i++)  material[i] -= mmt-1;
    if (nb) {
	for (i=1, mmb = labelF[0]; i<nb; i++)  if (mmb > labelF[i])  mmb = labelF[i];
	for (i=0; i<nb; i++)  labelF[i] -= mmb-1;
    }

    nvmax = (8*nv < nvmax)?8*nv:nvmax;
    ntmax = (8*nt < ntmax)?8*nt:ntmax;
    nbmax = (8*nb < nbmax)?8*nb:nbmax;

    MaxWi = 5*nvmax + 10*nbmax + 25*ntmax;
    MaxWr = 17*nvmax + 3*ntmax;

    iW = (int*)malloc(sizeof(int)*MaxWi);
    rW = (double*)malloc(sizeof(double)*MaxWr);

    if (nvfix)  ivfix = (int*)malloc(sizeof(int)*nvfix);
    else  ivfix = NULL;
    if (ntfix)  itfix = (int*)malloc(sizeof(int)*ntfix);
    else  itfix = NULL;
    if (nbfix)  ibfix = (int*)malloc(sizeof(int)*nbfix);
    else  ibfix = NULL;
    for (i=0; i<nvfix; i++)  ivfix[i] = i+1;
    for (i=0; i<ntfix; i++)  itfix[i] = i+1;
    for (i=0; i<nbfix; i++)  ibfix[i] = i+1;

    flagAuto = 1;
    status = 0;
    iPrint = 1;
    r = mbaFixShape(&nv, &nvmax, &nb, &nbmax, &nt, &ntmax, vrt, bnd, tet, labelF, material,
	    &nvfix, &nbfix, &ntfix, ivfix, ibfix, itfix, &flagAuto, &status,
	    &MaxSkipE, &MaxQItr, ANI_Metric_Eucl, &Quality, &rQuality,
	    &MaxWr, &MaxWi, rW, iW, &iPrint, &iERR);
    for (i=0; i<nt; i++) {
	d = idet3i4(vrt-3, tet[4*i+0], tet[4*i+1], tet[4*i+2], tet[4*i+3]);
	if (d == 0) printf("zero volume tet\n");
	if (d<0)  d = tet[4*i+0],  tet[4*i+0] = tet[4*i+1],  tet[4*i+1] = d;
    }
    for (i=0; i<nt; i++)  material[i] += mmt-1;
    for (i=0; i<nb; i++)  labelF[i] += mmb-1;
    *pnv=nv,  *pnt=nt,  *pnb=nb;
    return r;
}
