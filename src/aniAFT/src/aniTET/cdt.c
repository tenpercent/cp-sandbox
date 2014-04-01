#include <math.h>
#include <string.h>
#include <stdio.h>
#include "aft.h"
#include "ugly.h"
#include "delaunay.h"
#include "helper.h"
#include "det.h"
#include "libaft.h"

//#define STRUCTCHECK
//#define DUMPFRONT

static int dump_front(char *fname, int nV, REAL *vertex, int nF, int *face, int *facecolor);

typedef struct {
    int nV, nnV, nlV;
    REAL *vertex;
    int nF, nnF;
    int *face, *facecolor, *fcolor;
    int nT, nnT, nlT;
    int *tetra, *tetracolor;
    int ncolor, *color;
    int ntV, ntF, ntT;
    REAL *tvertex;
    int *tface, *tfacecolor;
    int *ttetra, *ttetracolor;
    int *vp, *vpi;
    int *v;
    int *e;
    int *fn;
} mesh3d;

static void normvec(REAL *vertex, int v1, int v2, int v3, REAL *px) {
    REAL r;
    normvec3i(vertex, v1, v2, v3, px);
    r = sqrt(px[0]*px[0] + px[1]*px[1] + px[2]*px[2]);
    if (r > 0) r = 1.0/r;
    px[0] *= r,  px[1] *= r,  px[2] *= r;
}

static REAL angle(mesh3d *pm, int f, int g) {
    REAL nf[3], ng[3];
    normvec(pm->vertex, pm->face[3*f+0], pm->face[3*f+1], pm->face[3*f+2], nf);
    normvec(pm->vertex, pm->face[3*g+0], pm->face[3*g+1], pm->face[3*g+2], ng);
    return nf[0]*ng[0] + nf[1]*ng[1] + nf[2]*ng[2];
}

static int prepare(mesh3d *pm) {
    int i, j, k, p, errors = 0;
    REAL mc, cc;
    for (i=0; i<pm->nF; i++) {
	for (j=0; j<3; j++) {
	    pm->e[2*(3*i+j)] = pm->face[3*i+(j+1)%3];
	    pm->e[2*(3*i+j)+1] = pm->v[pm->face[3*i+j]];
	    pm->v[pm->face[3*i+j]] = 3*i+j;
	}
    }
    for (i=0; i<pm->nF; i++) {
	for (j=0; j<3; j++) {
	    p = -1, mc = -1.0;
	    k = pm->v[pm->face[3*i+j]];
	    while (k >= 0) {
		if (pm->e[2*k] == pm->face[3*i+(j+2)%3]) {
		    cc = angle(pm, i, k/3);
		    if (idet3i4(pm->vertex, pm->face[3*i+0], pm->face[3*i+1], pm->face[3*i+2], pm->face[3*(k/3) + (k%3 + 2)%3]) < 0)  cc = 2.0 - cc;
		    //printf("%d: %lf\n", idet3i4(pm->vertex, pm->face[3*i+0], pm->face[3*i+1], pm->face[3*i+2], pm->face[3*(k/3) + (k%3 + 2)%3]), cc);
		    if ((p < 0) || (mc > cc)) {
			p  = k / 3;
			mc = cc;
		    }
		}
		k = pm->e[2*k+1];
	    }
	    if (p < 0) {
		errors++;
	    }
	    pm->fn[3*i + j] = p;
	}
    }
    return errors;
}

static int pcolor(mesh3d *pm, int c) {
    int i;
    if (c==0) return 0;
    for (i=0; i<pm->ncolor; i++) {
	if (pm->color[2*i+0]==c) break;
    }
    if (i==pm->ncolor) {
	pm->ncolor++;
	pm->color[2*i+0] = c;
	pm->color[2*i+1] = c;
    }
    if (pm->color[2*i+1] != c) {
	pm->color[2*i+1] = pcolor(pm, pm->color[2*i+1]);
    }
    return pm->color[2*i+1];
}
static int cjoin(mesh3d *pm, int c1, int c2) {
    int i, p1=pcolor(pm, c1), p2=pcolor(pm, c2);
    if ((p1==0) && (p2==0)) {
	return pcolor(pm, pm->ncolor+1);
    } else if (p1==0) return p2;
    else if (p2==0) return p1;
    else if (p1 > p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p1) break;
	}
	pm->color[2*i+1] = p2;
	return pcolor(pm, c1);
    } else if (p1 < p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p2) break;
	}
	pm->color[2*i+1] = p1;
	return pcolor(pm, c2);
    } else return p1;
}

static int fill(mesh3d *pm) {
    int j, color=1;
    for (j=0; j<pm->nF; j++) {
	pm->fcolor[j] = color;
    }
    return 0;
}

static int recolor(mesh3d *pm) {
    int j, f1, f2, f3, color;
    for (j=0; j<pm->nF; j++) {
	f1 = pm->fn[3*j+0];
	f2 = pm->fn[3*j+1];
	f3 = pm->fn[3*j+2];
	color = cjoin(pm, pm->fcolor[f1], pm->fcolor[f2]);
	color = cjoin(pm, color, pm->fcolor[f3]);
	pm->fcolor[f1] = color;
	pm->fcolor[f2] = color;
	pm->fcolor[f3] = color;
	pm->fcolor[j]  = color;
    }
    for (j=0; j<pm->nF; j++) {
	pm->fcolor[j] = pcolor(pm, pm->fcolor[j]);
    }
    return 0;
}

static int bbox(int nV, REAL *vertex, REAL *dv) {
    REAL bmin[3], bmax[3];
    int t, i, j, k1, k2, k3, k;

    for (j=0; j<3; j++) {
	bmin[j] = vertex[3*0+j];
	bmax[j] = vertex[3*0+j];
    }
    for (i=1; i<nV; i++) {
	for (j=0; j<3; j++) {
	    if (bmin[j] > vertex[3*i+j])  bmin[j] = vertex[3*i+j];
	    if (bmax[j] < vertex[3*i+j])  bmax[j] = vertex[3*i+j];
	}
    }
    for (t=0; t<3; t++) {
	dv[t] = (bmin[t] + bmax[t])/2.0;
    }
    for (t=0; t<3; t++) {
	dv[3+t] = 1.0;
	while ( dv[3+t] * (bmax[t] - bmin[t]) < 1.0)  dv[3+t] *= 2.0;
	while ( dv[3+t] * (bmax[t] - bmin[t]) > 1.0)  dv[3+t] /= 2.0;
    }
    k1 = 3,  k2 = 4,  k3 = 5;
    if (dv[k1]>dv[k2])  k=k1, k1=k2, k2=k;
    if (dv[k2]>dv[k3])  k=k2, k2=k3, k3=k;
    if (dv[k1]>dv[k2])  k=k1, k1=k2, k2=k;
    if (dv[k1]<dv[k2])  dv[k1] = dv[k2]/1.0;
    if (dv[k3]>dv[k2])  dv[k3] = dv[k2]*1.0;
   
    return 0;
}
void transform(REAL *v, REAL *d) {
    int t;
    for (t=0; t<3; t++) {
	v[t] = (v[t]-d[t])*d[3+t];
    }
}
void transformi(REAL *v, REAL *d) {
    int t;
    for (t=0; t<3; t++) {
	v[t] = v[t]/d[3+t] + d[t];
    }
}

static int cna(mesh3d *pm, int i) {
    int n;
    if (pm->vp[i] < 0) {
	n = pm->ntV++;
	pm->tvertex[3*n+0] = pm->vertex[3*i+0];
	pm->tvertex[3*n+1] = pm->vertex[3*i+1];
	pm->tvertex[3*n+2] = pm->vertex[3*i+2];
	pm->vp[i]=n;
	pm->vpi[n]=i;
    }
    return pm->vp[i];
}
static int remesh(mesh3d *pm) {
    int f, i, j, k, n, m, color, nV, r, err;
    REAL dv[6];

    pm->nlV = pm->nV;
    pm->nlT = pm->nT;

    err = 0;

    for (f=0; f<pm->nF; f++) {
	if (pm->fcolor[f]<=0) continue;
	color = pm->fcolor[f];
	pm->ntV = 0;
	pm->ntF = 0;
	for (j=0; j<pm->nF; j++) {
	    if (pm->fcolor[j] == color) {
		n = pm->ntF++;
		pm->tface[3*n+0] = cna(pm, pm->face[3*j+0]);
		pm->tface[3*n+1] = cna(pm, pm->face[3*j+1]);
		pm->tface[3*n+2] = cna(pm, pm->face[3*j+2]);
		pm->tfacecolor[n] = pm->facecolor[j];
	    }
	}
	pm->ntT = 0;
	nV = pm->ntV;

	bbox(pm->ntV, pm->tvertex, dv);

	for (i=0; i<pm->ntV; i++) transform(pm->tvertex+3*i, dv);
	if (0)  dump_front("cugly.bvf", pm->ntV, pm->tvertex, pm->ntF, pm->tface, pm->tfacecolor);
	if (check_surface_topology(pm->ntV, pm->tvertex, pm->ntF, pm->tface)) {
	    libaft_3d_warn("cdt.c: remesh(): wrong topology (internal error)");
	}
	r = mesh3dugly(&pm->ntV, pm->tvertex, &pm->ntF, pm->tface, pm->tfacecolor, &pm->ntT, pm->ttetra, pm->ttetracolor, pm->nnV, pm->nnF, pm->nnT);
	for (i=0; i<pm->ntV; i++) transformi(pm->tvertex+3*i, dv);

	if (r) {
/*	    libaft_3d_warn("Meshing failed. Resulted mesh will be coruppted!");*/
	    err++;
	    for (j=0; j<pm->nF; j++) {
		if (pm->fcolor[j] == color)  pm->fcolor[j] = -color;
	    }
	} else {
	    for (j=0; j<pm->nF; j++) {
		if (pm->fcolor[j] == color)  pm->fcolor[j] = 0;
	    }
	    for (i=nV; i<pm->ntV; i++) {
		pm->vertex[3*pm->nlV + 0] = pm->tvertex[3*i + 0];
		pm->vertex[3*pm->nlV + 1] = pm->tvertex[3*i + 1];
		pm->vertex[3*pm->nlV + 2] = pm->tvertex[3*i + 2];
		pm->vpi[i] = pm->nlV;
		pm->nlV++;
	    }
	    for (m=0; m<pm->ntT; m++) {
		pm->tetra[4*pm->nlT + 0] = pm->vpi[pm->ttetra[4*m + 0]];
		pm->tetra[4*pm->nlT + 1] = pm->vpi[pm->ttetra[4*m + 1]];
		pm->tetra[4*pm->nlT + 2] = pm->vpi[pm->ttetra[4*m + 2]];
		pm->tetra[4*pm->nlT + 3] = pm->vpi[pm->ttetra[4*m + 3]];
		pm->tetracolor[pm->nlT] = pm->ttetracolor[m];
		pm->nlT++;
	    }
	}

	for (k=0; k<nV; k++) {
	    pm->vp[pm->vpi[k]] = -1;
	}
    }
    for (j=0; j<pm->nF; j++) {
	if (pm->fcolor[j])  continue;
	pm->nF--;
	n = pm->nF;
	pm->face[3*j + 0] = pm->face[3*n + 0];
	pm->face[3*j + 1] = pm->face[3*n + 1];
	pm->face[3*j + 2] = pm->face[3*n + 2];
	pm->facecolor[j] = pm->facecolor[n];
	pm->fcolor[j] = pm->fcolor[n];
	j--;
    }
    if (!err && (pm->nF > 0))  libaft_3d_warn("cdt.c: remesh(): inconsistency detected (internal error)");
    pm->nT = pm->nlT;
    pm->nV = pm->nlV;
    return err;
}

int mesh3dcdt(int *pnV, REAL *vertex,
	int *pnF, int *face, int *facecolor,
	int *pnT, int *tetra, int *tetracolor,
	int nnV, int nnF, int nnT) {
    int i, r;
    mesh3d mesh;

    detinit();

    mesh.nV = *pnV, mesh.nF = *pnF, mesh.nT = *pnT;
    mesh.nnV = nnV, mesh.nnF = nnF, mesh.nnT = nnT;
    mesh.vertex = vertex, mesh.face = face, mesh.facecolor = facecolor, mesh.tetra = tetra, mesh.tetracolor = tetracolor;

    mesh.vp     = libaft_malloc(sizeof(int) * mesh.nV);
    mesh.vpi    = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.color  = libaft_malloc(sizeof(int) * mesh.nF*2);
    mesh.fcolor = libaft_malloc(sizeof(int) * mesh.nF);
    memset(mesh.fcolor, 0, sizeof(int) * mesh.nF);

    mesh.v      = libaft_malloc(sizeof(int) * mesh.nV);
    mesh.e      = libaft_malloc(sizeof(int) * mesh.nF*6);
    mesh.fn     = libaft_malloc(sizeof(int) * mesh.nF*3);

    mesh.tvertex     = libaft_malloc(sizeof(REAL) * mesh.nnV*3);
    mesh.tface       = libaft_malloc(sizeof(int)  * mesh.nnF*3);
    mesh.tfacecolor  = libaft_malloc(sizeof(int)  * mesh.nnF);
    mesh.ttetra      = libaft_malloc(sizeof(int)  * mesh.nnT*4);
    mesh.ttetracolor = libaft_malloc(sizeof(int)  * mesh.nnT);

    mesh.ncolor = 0;

    for (i=0; i<mesh.nV; i++) mesh.vp[i] = -1;

    for (i=0; i<mesh.nV; i++) mesh.v[i] = -1;
    
    if (prepare(&mesh))  libaft_3d_warn("cdt.c: mesh3dcdt(): prepare() failed (internal error)");

    if (1)  recolor(&mesh);  else  fill(&mesh);
    r = remesh(&mesh);

    libaft_free(mesh.v);
    libaft_free(mesh.e);
    libaft_free(mesh.fn);
    libaft_free(mesh.vp);
    libaft_free(mesh.vpi);
    libaft_free(mesh.color);
    libaft_free(mesh.fcolor);
    libaft_free(mesh.tvertex);
    libaft_free(mesh.tface);
    libaft_free(mesh.tfacecolor);
    libaft_free(mesh.ttetra);
    libaft_free(mesh.ttetracolor);
    *pnV = mesh.nV, *pnF = mesh.nF, *pnT = mesh.nT;
    return r;
}

static int dump_front(char *fname, int nV, REAL *vertex, int nF, int *face, int *facecolor) {
    FILE *f;
    f = fopen(fname, "w");
    if (!f)  return perror(fname), 1;
    fwrite(&nV, sizeof(int), 1, f);
    fwrite(vertex, sizeof(REAL), 3*nV, f);
    fwrite(&nF, sizeof(int), 1, f);
    fwrite(face, sizeof(int), 3*nF, f);
    fwrite(facecolor, sizeof(int), nF, f);
    fclose(f);
    return 0;
}
