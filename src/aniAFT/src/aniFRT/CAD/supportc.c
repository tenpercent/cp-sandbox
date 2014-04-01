#include <stdlib.h>
#include <string.h>
#include "cgmwrap.h"
#include "common.h"
#include "region3c.h"
#include "supportc.h"
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
    return 0;
}

int aft3dmodel (
	CGMmodel model,
	double (*fsize) (double, double, double),
	int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int maxnV, int maxnF,
	int indexshift
) {
    surface_mesh mesh;
    int r;

    memset(&mesh, 0, sizeof(mesh));
    mesh.fsize    = fsize;
    initAFSM_(&mesh, model);
    r = niceit(&mesh, pnVout, vertexout, pnFout, faceout, facecolor, maxnV, maxnF, indexshift);
    freeMemory(&mesh);
    return r;
}

