#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "aft.h"
#include "aft2d.h"
#include "aft3d.h"
#include "ugly.h"
#include "cdt.h"
#include "check.h"
#include "helper.h"
#include "libaft.h"
#include "qual.h"

static int chk2dvert(int nV, REAL *vertex, int nE, int *edge) {
	int i, j, k=0, *p;
	p = libaft_malloc(sizeof(int)*nV);
	for (i=0; i<nV; i++) {
		for (j=0; j<i; j++) {
			if ((vertex[2*i+0] == vertex[2*j+0]) && (vertex[2*i+1] == vertex[2*j+1])) break;
		}
		if (i!=j) k++;
		p[i] = j;
	}
	for (i=0; i<nE; i++) {
		edge[2*i+0] = 1+p[edge[2*i+0]-1];
		edge[2*i+1] = 1+p[edge[2*i+1]-1];
	}
	libaft_free(p);
	return k;
}

int mesh_2d_aft_opts(int *pnV, REAL *vertex,
		int *pnE, int *edge, int *edgecolor,
		int *pnT, int *tria, int *triacolor,
		int nnV, int nnE, int nnT,
		REAL cf, REAL lim, int loop, int needcheck, int indexshift) {
	int i, k, nV=*pnV, last, invert=1, ret;
	if (loop) {
		if (*pnE != 0) libaft_2d_warn("loop enabled, but nE != 0");
		last = -1, k = 0;
		for (i=0; i<nV; i++) {
			edge[2*k+0+invert] = i+1;
			edge[2*k+1-invert] = i+2;
			edgecolor[k] = 1;
			if ( (last>=0) && (vertex[2*last+0]==vertex[2*i+0]) && (vertex[2*last+1]==vertex[2*i+1]) ) {
				edge[2*(k-1)+1-invert] = last+1;
				last = -1;
			} else {
				k++;
				if (last<0) last = i;
			}
			*pnE = k;
		}
	}
	if (needcheck) chk2dvert(nV, vertex, *pnE, edge);
	if (indexshift) {
		for (i=0; i<2**pnE; edge[i++]-=indexshift);
		for (i=0; i<3**pnT; tria[i++]-=indexshift);
	}
	if (cf>0.0) ret = mesh2daftsslim(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, lim);
	else        ret = mesh2daft  (pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT);
	if (indexshift) {
		for (i=0; i<2**pnE; edge[i++]+=indexshift);
		for (i=0; i<3**pnT; tria[i++]+=indexshift);
	}
	return ret;
}

int mesh_2d_aft_full_opts(int *pnV, REAL *vertex,
		     int *pnE, int *edge, int *edgecolor,
		     int *pnT, int *tria, int *triacolor,
		     int nnV, int nnE, int nnT,
		     REAL cf, REAL lim, int loop, int needcheck, int indexshift,
		     int *pnC, int *color) {
    int i, k, nV=*pnV, last, invert=1, ret;
    if (loop) {
	if (*pnE != 0) libaft_2d_warn("loop enabled, but nE != 0");
	last = -1, k = 0;
	for (i=0; i<nV; i++) {
	    edge[2*k+0+invert] = i+1;
	    edge[2*k+1-invert] = i+2;
	    edgecolor[k] = 1;
	    if ( (last>=0) && (vertex[2*last+0]==vertex[2*i+0]) && (vertex[2*last+1]==vertex[2*i+1]) ) {
		edge[2*(k-1)+1-invert] = last+1;
		last = -1;
	    } else {
		k++;
		if (last<0) last = i;
	    }
	    *pnE = k;
	}
    }
    if (needcheck) chk2dvert(nV, vertex, *pnE, edge);
    if (indexshift) {
	for (i=0; i<2**pnE; edge[i++]-=indexshift);
	for (i=0; i<3**pnT; tria[i++]-=indexshift);
    }
    ret = mesh2daftfull(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, lim, pnC, color);
    if (indexshift) {
	for (i=0; i<2**pnE; edge[i++]+=indexshift);
	for (i=0; i<3**pnT; tria[i++]+=indexshift);
    }
    return ret;
}

int mesh_2d_aft(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 0, 0, 0);
}
int mesh_2d_aft_cf(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 0, 0, 0);
}
int mesh_2d_aft_loop(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 1, 0, 0);
}
int mesh_2d_aft_cf_loop(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 1, 0, 0);
}
int mesh_2d_aft_check(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 0, 1, 0);
}
int mesh_2d_aft_cf_check(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 0, 1, 0);
}
int mesh_2d_aft_loop_check(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 1, 1, 0);
}
int mesh_2d_aft_cf_loop_check(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 1, 1, 0);
}
int mesh_2d_aft_auto(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
	if (*pnE) return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 0, 1, 0);
	else return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 0, -1, 1, 1, 0);
}
int mesh_2d_aft_cf_auto(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL cf) {
	if (*pnE) return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 0, 1, 0);
	else return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, cf, -1, 1, 1, 0);
}


int mesh_2d_aft_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 0, 0, 1);
}
int mesh_2d_aft_cf_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 0, 0, 1);
}
int mesh_2d_aft_loop_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 1, 0, 1);
}
int mesh_2d_aft_cf_loop_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 1, 0, 1);
}
int mesh_2d_aft_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 0, 1, 1);
}
int mesh_2d_aft_cf_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 0, 1, 1);
}
int mesh_2d_aft_loop_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 1, 1, 1);
}
int mesh_2d_aft_cf_loop_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 1, 1, 1);
}
int mesh_2d_aft_cf_lim_loop_check_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf, REAL *plim) {
	return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, *plim, 1, 1, 1);
}
int mesh_2d_aft_auto_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT) {
	if (*pnE) return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 0, 1, 1);
	else return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, 0, -1, 1, 1, 1);
}
int mesh_2d_aft_cf_auto_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf) {
	if (*pnE) return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 0, 1, 1);
	else return mesh_2d_aft_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 1, 1, 1);
}
int mesh_2d_aft_full_(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int *pnnV, int *pnnE, int *pnnT, REAL *pcf, int *pnC, int *color) {
	return mesh_2d_aft_full_opts(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, *pnnV, *pnnE, *pnnT, *pcf, -1, 0, 1, 1, pnC, color);
}





static int chk3dvert(int nV, REAL *vertex, int nF, int *face) {
	int i, j, k=0, *p;
	p = libaft_malloc(sizeof(int)*nV);
	for (i=0; i<nV; i++) {
		for (j=0; j<i; j++) {
			if (  (vertex[3*i+0] == vertex[3*j+0]) &&
					(vertex[3*i+1] == vertex[3*j+1]) &&
					(vertex[3*i+2] == vertex[3*j+2])  ) break;
		}
		if (i!=j) k++;
		p[i] = j;
	}
	for (i=0; i<nF; i++) {
		face[3*i+0] = 1+p[face[3*i+0]-1];
		face[3*i+1] = 1+p[face[3*i+1]-1];
		face[3*i+2] = 1+p[face[3*i+2]-1];
	}
	libaft_free(p);
	return k;
}

int mesh_3d_ugly_opts(int *pnV, REAL *vertex,
		int *pnF, int *face, int *facecolor,
		int *pnT, int *tetra, int *tetracolor,
		int nnV, int nnF, int nnT,
		int indexshift) {
	int i, ret;
	if (indexshift) {
		for (i=0; i<3**pnF; face[i++]-=indexshift);
		for (i=0; i<4**pnT; tetra[i++]-=indexshift);
	}
//	ret = mesh3dugly(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT);
	ret = mesh3dcdt(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT);
	if (indexshift) {
		for (i=0; i<3**pnF; face[i++]+=indexshift);
		for (i=0; i<4**pnT; tetra[i++]+=indexshift);
	}
	return ret;
}

int mesh_3d_aft_opts(int *pnV, REAL *vertex,
		int *pnF, int *face, int *facecolor,
		int *pnT, int *tetra, int *tetracolor,
		int nnV, int nnF, int nnT,
		REAL cf, int needcheck, int indexshift, double (*f)(double, double, double)) {
	int i, ret;
	if (needcheck) chk3dvert(*pnV, vertex, *pnF, face);
	if (indexshift) {
		for (i=0; i<3**pnF; face[i++]-=indexshift);
		for (i=0; i<4**pnT; tetra[i++]-=indexshift);
	}
	if (cf>0.0) ret = mesh3daftss  (pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, cf);
	else if (f) ret = mesh3daftuser(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, f);
	else        ret = mesh3daft    (pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT);
	if (indexshift) {
		for (i=0; i<3**pnF; face[i++]+=indexshift);
		for (i=0; i<4**pnT; tetra[i++]+=indexshift);
	}
	return ret;
}

int mesh_3d_aft(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, 0, 0, 0, 0);
}
int mesh_3d_aft_cf(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, REAL cf) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, cf, 0, 0, 0);
}
int mesh_3d_aft_func(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, double (*f)(double, double, double)) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, 0, 0, 0, f);
}
int mesh_3d_aft_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, *pnnV, *pnnF, *pnnT, 0, 0, 1, 0);
}
int mesh_3d_aft_cf_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT, REAL *pcf) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, *pnnV, *pnnF, *pnnT, *pcf, 0, 1, 0);
}
int mesh_3d_aft_func_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT, double (*f)(double, double, double)) {
	return mesh_3d_aft_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, *pnnV, *pnnF, *pnnT, 0, 0, 1, f);
}

int mesh_3d_ugly(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT) {
	return mesh_3d_ugly_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, 0);
}
int mesh_3d_ugly_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT) {
	return mesh_3d_ugly_opts(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, *pnnV, *pnnF, *pnnT, 1);
}

#define BUFSIZE 10240
int read_front(char *filename, int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int maxnV, int maxnF, int zeroindex, int invert, int verbose) {
	int nF=*pnF, nV=*pnV;
	FILE *f;
	char buf[BUFSIZE];
	int i, j, r, nv, nf;

	if (!(f = fopen(filename, "r"))) {
		libaft_3d_warn("read_front: %s: %s", filename, strerror(errno));
		return 1;
	}
	if (verbose) printf("Reading mesh file %s.\n", filename);

	fgets(buf, BUFSIZE, f);
	i = sscanf(buf, "%d %d", &nv, &nf);
	nV += nv;
	nF += nf;
	if (nV>maxnV) {
		libaft_3d_warn("read_front: %s: maxnV", filename);
		fclose(f);
		return 2;
	}
	if (nF>maxnF) {
		libaft_3d_warn("read_front: %s: maxnF", filename);
		fclose(f);
		return 2;
	}

	for (i=nV-nv; i<nV; i++) {
		fgets(buf, BUFSIZE, f);
		sscanf(buf, "%lf %lf %lf", vertex+3*i+0, vertex+3*i+1, vertex+3*i+2);
	}
	for (i=nF-nf; i<nF; i++) {
		fgets(buf, BUFSIZE, f);
		r = sscanf(buf, "%d %d %d %d", face+3*i+0+invert, face+3*i+1-invert, face+3*i+2, facematerial+i);
		if (r < 4) facematerial[i] = 1;
		if (zeroindex) {
			face[3*i+0]++;
			face[3*i+1]++;
			face[3*i+2]++;
		}
		for (j=0; j<3; j++) face[3*i+j] += nV-nv;
	}
	fclose(f);
	if (verbose) printf("nV = %d, nF = %d\n", nv, nf);
	*pnF = nF;
	*pnV = nV;
	return 0;
}
#undef BUFSIZE

int write_front(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial) {
	int i;
	FILE *f;
	f=fopen(fn, "w");
	if (!f) {
		libaft_3d_warn("write_front: %s: %s", fn, strerror(errno));
		return 1;
	}
	fprintf(f, "%d %d\n", nV, nF);
	for (i=0; i<nV; i++)
		fprintf(f, "%20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	for (i=0; i<nF; i++)
		fprintf(f, "%d %d %d %d\n", face[3*i+0], face[3*i+1], face[3*i+2], facematerial[i]);
	fclose(f);
	return 0;
}

int write_front_gmv(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial) {
	int i,mmax=0, unkmat=0;
	FILE *f;
	f=fopen(fn, "w");
	fprintf(f, "gmvinput ascii\n\nnodev %d\n", nV);
	for (i=0; i<nV; i++) {
		fprintf(f, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	fprintf(f, "\ncells %d\n", nF);
	for (i=0; i<nF; i++) {
		fprintf(f, " tri 3\n  %d %d %d\n", face[3*i+0], face[3*i+1], face[3*i+2]);
		if (mmax<facematerial[i]) mmax = facematerial[i];
		if (facematerial[i]==0) unkmat = 1;
	}
	if (nF) {
		fprintf(f, "\nmaterial %d 0\n", mmax+unkmat);
		for (i=0; i<mmax; i++) fprintf(f,"mat%d\n", i+1);
		if (unkmat)	fprintf(f,"unknown\n");
		for (i=0; i<nF; i++) {
			fprintf(f, " %d", (facematerial[i])?facematerial[i]:mmax+unkmat);
		}
	}

	fprintf(f, "\nendgmv");
	fclose(f);
	return 0;
}

int write_mesh_gmv(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial) {
	int i,mmax=0, unkmat=0;
	FILE *f;
	f=fopen(fn, "w");
	fprintf(f, "gmvinput ascii\n\nnodev %5d\n", nV);
	for (i=0; i<nV; i++) {
		fprintf(f, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	fprintf(f, "\ncells %5d\n", nT);
	for (i=0; i<nT; i++) {
		fprintf(f, " tet 4\n  %4d %4d %4d %4d\n", tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
		if (mmax<tetramaterial[i]) mmax = tetramaterial[i];
		if (tetramaterial[i]==0) unkmat = 1;
	}
	for (i=0; i<nF; i++) {
		if (mmax<facematerial[i]) mmax = facematerial[i];
		if (facematerial[i]==0) unkmat = 1;
	}
	if (nT) {
		fprintf(f, "\nmaterial %d 0\n", mmax+unkmat);
		for (i=0; i<mmax; i++) fprintf(f,"mat%d\n", i+1);
		if (unkmat)	fprintf(f,"unknown\n");
		for (i=0; i<nT; i++) {
			fprintf(f, " %2d", (tetramaterial[i])?tetramaterial[i]:mmax+unkmat);
		}
	}

	fprintf(f, "\n\npolygons\n");
	for (i=0; i<nF; i++) {
		fprintf(f, "%3d 3", (facematerial[i])?facematerial[i]:mmax+unkmat);
		fprintf(f,      " %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+0-3], vertex[3*face[3*i+1]+0-3], vertex[3*face[3*i+2]+0-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+1-3], vertex[3*face[3*i+1]+1-3], vertex[3*face[3*i+2]+1-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+2-3], vertex[3*face[3*i+1]+2-3], vertex[3*face[3*i+2]+2-3]);
	}
	fprintf(f, "endpoly\n");
	fprintf(f, "\nendgmv");
	fclose(f);
	return 0;
}

int write_mesh_gmv_qual_mba(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial) {
	int i,mmax=0, unkmat=0;
	FILE *f;
	f=fopen(fn, "w");
	fprintf(f, "gmvinput ascii\n\nnodev %5d\n", nV);
	for (i=0; i<nV; i++) {
		fprintf(f, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	fprintf(f, "\ncells %5d\n", nT);
	for (i=0; i<nT; i++) {
		fprintf(f, " tet 4\n  %4d %4d %4d %4d\n", tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
		if (mmax<tetramaterial[i]) mmax = tetramaterial[i];
		if (tetramaterial[i]==0) unkmat = 1;
	}
	for (i=0; i<nF; i++) {
		if (mmax<facematerial[i]) mmax = facematerial[i];
		if (facematerial[i]==0) unkmat = 1;
	}
	if (nT) {
		fprintf(f, "\nmaterial %d 0\n", mmax+unkmat);
		for (i=0; i<mmax; i++) fprintf(f,"mat%d\n", i+1);
		if (unkmat)	fprintf(f,"unknown\n");
		for (i=0; i<nT; i++) {
			fprintf(f, " %2d", (tetramaterial[i])?tetramaterial[i]:mmax+unkmat);
		}
		fprintf(f, "\n\nvariable\nqual 0\n");
		for (i=0; i<nT; i++) {
			fprintf(f, " %20.15lf", aft_tetra_qual_mba(vertex, tetra[4*i+0]-1, tetra[4*i+1]-1, tetra[4*i+2]-1, tetra[4*i+3]-1));
		}
		fprintf(f, "\nendvars\n");
	}

	fprintf(f, "\n\npolygons\n");
	for (i=0; i<nF; i++) {
		fprintf(f, "%3d 3", (facematerial[i])?facematerial[i]:mmax+unkmat);
		fprintf(f,      " %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+0-3], vertex[3*face[3*i+1]+0-3], vertex[3*face[3*i+2]+0-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+1-3], vertex[3*face[3*i+1]+1-3], vertex[3*face[3*i+2]+1-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+2-3], vertex[3*face[3*i+1]+2-3], vertex[3*face[3*i+2]+2-3]);
	}
	fprintf(f, "endpoly\n");
	fprintf(f, "\nendgmv");
	fclose(f);
	return 0;
}

int write_mesh_gmv_qual(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial) {
	int i,mmax=0, unkmat=0;
	FILE *f;
	f=fopen(fn, "w");
	fprintf(f, "gmvinput ascii\n\nnodev %5d\n", nV);
	for (i=0; i<nV; i++) {
		fprintf(f, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	fprintf(f, "\ncells %5d\n", nT);
	for (i=0; i<nT; i++) {
		fprintf(f, " tet 4\n  %4d %4d %4d %4d\n", tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
		if (mmax<tetramaterial[i]) mmax = tetramaterial[i];
		if (tetramaterial[i]==0) unkmat = 1;
	}
	for (i=0; i<nF; i++) {
		if (mmax<facematerial[i]) mmax = facematerial[i];
		if (facematerial[i]==0) unkmat = 1;
	}
	if (nT) {
		fprintf(f, "\nmaterial %d 0\n", mmax+unkmat);
		for (i=0; i<mmax; i++) fprintf(f,"mat%d\n", i+1);
		if (unkmat)	fprintf(f,"unknown\n");
		for (i=0; i<nT; i++) {
			fprintf(f, " %2d", (tetramaterial[i])?tetramaterial[i]:mmax+unkmat);
		}
		fprintf(f, "\n\nvariable\nqual 0\n");
		for (i=0; i<nT; i++) {
			fprintf(f, " %20.15lf", aft_tetra_qual(vertex, tetra[4*i+0]-1, tetra[4*i+1]-1, tetra[4*i+2]-1, tetra[4*i+3]-1));
		}
		fprintf(f, "\nendvars\n");
	}

	fprintf(f, "\n\npolygons\n");
	for (i=0; i<nF; i++) {
		fprintf(f, "%3d 3", (facematerial[i])?facematerial[i]:mmax+unkmat);
		fprintf(f,      " %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+0-3], vertex[3*face[3*i+1]+0-3], vertex[3*face[3*i+2]+0-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+1-3], vertex[3*face[3*i+1]+1-3], vertex[3*face[3*i+2]+1-3]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", vertex[3*face[3*i+0]+2-3], vertex[3*face[3*i+1]+2-3], vertex[3*face[3*i+2]+2-3]);
	}
	fprintf(f, "endpoly\n");
	fprintf(f, "\nendgmv");
	fclose(f);
	return 0;
}

int write_mesh(char *fn, int nV, double *vertex, int nF, int *face, int *facematerial, int nT, int *tetra, int *tetramaterial) {
	int i;
	FILE *f;
	f=fopen(fn, "w");
	fprintf(f, "%5d\n", nV);
	for (i=0; i<nV; i++) {
		fprintf(f, "%20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	fprintf(f, "%5d\n", nT);
	for (i=0; i<nT; i++) {
		fprintf(f, "%4d %4d %4d %4d  %d\n", tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3], tetramaterial[i]);
	}
	fprintf(f, "%5d\n", nF);
	for (i=0; i<nF; i++) {
		fprintf(f, "%4d %4d %4d  %d\n", face[3*i+0], face[3*i+1], face[3*i+2], facematerial[i]);
	}
	fclose(f);
	return 0;
}


int check_mesh_topology(int nV, REAL *vertex, int nF, int *face, int nT, int *tetra) {
	return checktopology(nV, vertex, nF, face, nT, tetra, 0);
}
int check_mesh_topology_(int *pnV, REAL *vertex, int *pnF, int *face, int *pnT, int *tetra) {
	return checktopology(*pnV, vertex, *pnF, face, *pnT, tetra, 1);
}

int check_surface_topology(int nV, REAL *vertex, int nF, int *face) {
	return check2dsurface(&nV, vertex, &nF, face, 0, 0, 0);
}
int check_surface_topology_(int *pnV, REAL *vertex, int *pnF, int *face) {
	return check2dsurface(pnV, vertex, pnF, face, 1, 0, 0);
}
int check_surface_topology_fix(int *pnV, REAL *vertex, int *pnF, int *face) {
	return check2dsurface(pnV, vertex, pnF, face, 0, 1, 0);
}
int check_surface_topology_fix_(int *pnV, REAL *vertex, int *pnF, int *face) {
	return check2dsurface(pnV, vertex, pnF, face, 1, 1, 0);
}


