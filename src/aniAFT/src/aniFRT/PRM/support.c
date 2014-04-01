#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"
#include "region3.h"
#include "support.h"
#include "memory3.h"

static int niceit(
	surface_mesh *pm,
	int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int maxnV, int maxnF,
	int is
) {
    int i, nf, nv;
    (void) maxnV, (void) maxnF;
    nv = *pnVout;
    nf = *pnFout;
    *pnVout += pm->nPoint;
    for (i=nv; i<*pnVout; i++) {
	vertexout[3*i+0] = pm->vert[i-nv].x;
	vertexout[3*i+1] = pm->vert[i-nv].y;
	vertexout[3*i+2] = pm->vert[i-nv].z;
    }
    *pnFout += pm->nFace;
    for (i=nf; i<*pnFout; i++) {
	faceout[3*i+0] = pm->face[i-nf].v1 + nv + is;
	faceout[3*i+1] = pm->face[i-nf].v2 + nv + is;
	faceout[3*i+2] = pm->face[i-nf].v3 + nv + is;
	facecolor[i]   = pm->face[i-nf].color;
    }
    if (pm->exportCurves) {
	for (i=0; i<pm->nexpEdge; i++) {
	    pm->expEdge[2*i+0] += is;
	    pm->expEdge[2*i+1] += is;
	}
    }
    return 0;
}

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
) {
    surface_mesh mesh;
    int r;

    memset(&mesh, 0, sizeof(mesh));
    mesh.v_u      = v_u;
    mesh.bounsurf = bounsurf;
    mesh.bounline = bounline;
    mesh.periodicfunction = periodicfunction;
    mesh.fsize    = fsize;
    mesh.exportCurves = exportCurves;
    if (pnEout) {
	if (*pnEout != 0)  printf("nE != 0 is not supported yet!");
	mesh.nexpEdge = *pnEout;
    } else mesh.nexpEdge = -1;
    mesh.expEdge = edgeout;
    mesh.expEdgeColor = edgecolor;
    mesh.nnexpEdge = maxnE;
    initAFS_(&mesh, &nVVert, VVertxyz, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT);
    r = niceit(&mesh, pnVout, vertexout, pnFout, faceout, facecolor, maxnV, maxnF, indexshift);
    if (pnEout)  *pnEout = mesh.nexpEdge;
    freeMemory(&mesh);
    return r;
}

