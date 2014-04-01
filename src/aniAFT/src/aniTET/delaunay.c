#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "aft.h"
#include "delaunay.h"
#include "helper.h"
#include "check.h"
#include "det.h"

#define STRUCTCHECK

#define TETRA 1

typedef struct {
    REAL x, y, z;
    LONGREAL r;
} sphere;

#define TMW 4
typedef struct {
    int nV, nnV, nVF;
    REAL *vertex;
    int nT, nnT;
    int *tetra;
    REAL *tetrametric;
    sphere *tetrasphere;
    int nF, nnF;
    int *face;
    int  ntpack, *tpack;
    int  npack,  *pack;
} delaunay3d;

/* debug functions */
static int savefrtraw(delaunay3d *pm);


/* Metric functions */
static LONGREAL dist3s(REAL x1, REAL y1, REAL z1,  REAL x2, REAL y2, REAL z2) {
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
}
static LONGREAL dist3si(REAL *vertex, int v1, int v2) {
    return dist3s(vertex[3*v1+0], vertex[3*v1+1], vertex[3*v1+2], vertex[3*v2+0], vertex[3*v2+1], vertex[3*v2+2]);
}
static LONGREAL dist3(REAL x1, REAL y1, REAL z1,  REAL x2, REAL y2, REAL z2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}
static LONGREAL dist3vt(delaunay3d *pm, int v, int t) {
    return dist3s(pm->vertex[3*v+0], pm->vertex[3*v+1], pm->vertex[3*v+2],
	    pm->tetrasphere[t].x, pm->tetrasphere[t].y, pm->tetrasphere[t].z);
}

/*
static int intetra(REAL *vertex, int *t, int k) {
    REAL v = det3i4(vertex, t[0], t[1], t[2], t[3]);
    REAL lim = -1e-2 * v;
    if (v <= 0.0)  libaft_3d_warn("v = %le in intetra()", v);

    if (det3i4(vertex, k, t[1], t[2], t[3]) < lim)  return printf("x... (%18le)\n", det3i4(vertex, k, t[1], t[2], t[3])/v),  0;
    if (det3i4(vertex, t[0], k, t[2], t[3]) < lim)  return printf("+x.. (%18le)\n", det3i4(vertex, t[0], k, t[2], t[3])/v),  0;
    if (det3i4(vertex, t[0], t[1], k, t[3]) < lim)  return printf("++x. (%18le)\n", det3i4(vertex, t[0], t[1], k, t[3])/v),  0;
    if (det3i4(vertex, t[0], t[1], t[2], k) < lim)  return printf("+++x (%18le)\n", det3i4(vertex, t[0], t[1], t[2], k)/v),  0;
    return printf("++++\n"),  1;
}
*/
static int intetra(REAL *vertex, int *t, int k) {
    REAL v = det3i4(vertex, t[0], t[1], t[2], t[3]);
    REAL lim = -1e-4 * v;
    if (v <= 0.0)  libaft_3d_warn("v = %le in intetra()", v);

    if (det3i4(vertex, k, t[1], t[2], t[3]) < lim)  return 0;
    if (det3i4(vertex, t[0], k, t[2], t[3]) < lim)  return 0;
    if (det3i4(vertex, t[0], t[1], k, t[3]) < lim)  return 0;
    if (det3i4(vertex, t[0], t[1], t[2], k) < lim)  return 0;
    return 1;
}
/*
static int intetra(REAL *vertex, int *t, int k) {
    (void) vertex,  (void) t,  (void) k;
    return 0;
}
*/
static int insphere(sphere *s, REAL *vertex) {
    LONGREAL d = dist3s(s->x, s->y, s->z, vertex[0], vertex[1], vertex[2]) - s->r;
    int r;
//    printf("d = %20Le ", d);
    r = (d >= -1e-12) ? 0 : 1;
//    printf("%d\n", r);
    return r;
}

/* structure functions */
static int addtetra(delaunay3d *pm, int v1, int v2, int v3, int v4) {
    if (idet3i4(pm->vertex, v1, v2, v3, v4) < 0)
	libaft_3d_warn("delaunay.c: addtetra(): bad orientation of tetra (negative), [%d %d %d %d], v = %le",
		v1, v2, v3, v4, 
		det3i4(pm->vertex, v1, v2, v3, v4));
    if (pm->nT >= pm->nnT) {
	libaft_3d_stop("Delaunay: maxTetra exceeded");
	return -1;
    } else {
	pm->tetra[4*(pm->nT)+0] = v1;
	pm->tetra[4*(pm->nT)+1] = v2;
	pm->tetra[4*(pm->nT)+2] = v3;
	pm->tetra[4*(pm->nT)+3] = v4;
	return (pm->nT)++;
    }
}
static int addvertex(delaunay3d *pm, REAL x, REAL y, REAL z) {
    if (pm->nV >= pm->nnV) {
	libaft_3d_stop("Delaunay: maxVertex exceeded");
	return -1;
    } else {
	pm->vertex[3*(pm->nV)+0] = x;
	pm->vertex[3*(pm->nV)+1] = y;
	pm->vertex[3*(pm->nV)+2] = z;
	return (pm->nV)++;
    }
}
static void tetrametric(delaunay3d *pm, int i) {
    int v0, v1, v2, v3;
    REAL v;
    REAL w0[4], w1[4], w2[4], w3[4];
    LONGREAL l1, l2, l3;


    v0 = pm->tetra[4*i+0];
    v1 = pm->tetra[4*i+1];
    v2 = pm->tetra[4*i+2];
    v3 = pm->tetra[4*i+3];

    w0[0] = pm->vertex[3*v0+0],  w1[0] = pm->vertex[3*v1+0],  w2[0] = pm->vertex[3*v2+0],  w3[0] = pm->vertex[3*v3+0];
    w0[1] = pm->vertex[3*v0+1],  w1[1] = pm->vertex[3*v1+1],  w2[1] = pm->vertex[3*v2+1],  w3[1] = pm->vertex[3*v3+1];
    w0[2] = pm->vertex[3*v0+2],  w1[2] = pm->vertex[3*v1+2],  w2[2] = pm->vertex[3*v2+2],  w3[2] = pm->vertex[3*v3+2];
    w0[3] = pm->vertex[3*v0+0],  w1[3] = pm->vertex[3*v1+0],  w2[3] = pm->vertex[3*v2+0],  w3[3] = pm->vertex[3*v3+0];

    v = det3i4(pm->vertex, v0, v1, v2, v3);
    if (v < 0.0) libaft_3d_warn("delaunay.c: tetrametric(): negative volume: %le (internal error)", v);
    v = -v;
    l1 = dist3si(pm->vertex, v1, v0);
    l2 = dist3si(pm->vertex, v2, v0);
    l3 = dist3si(pm->vertex, v3, v0);

    pm->tetrasphere[i].x = pm->vertex[3*v0+0] + (orient2d(w2+1, w3+1, w0+1)*l1 + orient2d(w3+1, w1+1, w0+1)*l2 + orient2d(w1+1, w2+1, w0+1)*l3)/(2.0*v);
    pm->tetrasphere[i].y = pm->vertex[3*v0+1] + (orient2d(w2+2, w3+2, w0+2)*l1 + orient2d(w3+2, w1+2, w0+2)*l2 + orient2d(w1+2, w2+2, w0+2)*l3)/(2.0*v);
    pm->tetrasphere[i].z = pm->vertex[3*v0+2] + (orient2d(w2+0, w3+0, w0+0)*l1 + orient2d(w3+0, w1+0, w0+0)*l2 + orient2d(w1+0, w2+0, w0+0)*l3)/(2.0*v);

    pm->tetrasphere[i].r = dist3vt(pm, v0, i);

/*    libaft_3d_warn("New tetra metric:\n (%lf, %lf, %lf)\n (%lf, %lf, %lf)\n (%lf, %lf, %lf)\n (%lf, %lf, %lf)\n %lf %lf %lf %lf",
	    w0[0], w0[1], w0[2],
	    w1[0], w1[1], w1[2],
	    w2[0], w2[1], w2[2],
	    w3[0], w3[1], w3[2],
	    dist3vt(pm, v0, i), dist3vt(pm, v1, i), dist3vt(pm, v2, i), dist3vt(pm, v3, i));*/
}
static void remtetra(delaunay3d *pm, int i) {
    if (pm->ntpack > pm->nnT) libaft_3d_stop("delaunay.c: Internal error in remtetra & packtetra");
    pm->tetra[4*i+0] = -1;
    pm->tetra[4*i+1] = -1;
    pm->tetra[4*i+2] = -1;
    pm->tetra[4*i+3] = -1;
    pm->tpack[pm->ntpack] = i;
    pm->ntpack++;
}
static void packtetra(delaunay3d *pm) {
    int i, j;
    while (pm->ntpack > 0) {
	pm->ntpack--;
	j = pm->tpack[pm->ntpack];
	pm->nT--;
	pm->tetra[4*j+0] = pm->tetra[4*(pm->nT)+0];
	pm->tetra[4*j+1] = pm->tetra[4*(pm->nT)+1];
	pm->tetra[4*j+2] = pm->tetra[4*(pm->nT)+2];
	pm->tetra[4*j+3] = pm->tetra[4*(pm->nT)+3];
	pm->tetrasphere[j].x = pm->tetrasphere[(pm->nT)].x;
	pm->tetrasphere[j].y = pm->tetrasphere[(pm->nT)].y;
	pm->tetrasphere[j].z = pm->tetrasphere[(pm->nT)].z;
	pm->tetrasphere[j].r = pm->tetrasphere[(pm->nT)].r;
	for (i=0; i<pm->ntpack; i++) if (pm->tpack[i] == pm->nT) pm->tpack[i] = j;
    }
}
static int findface(delaunay3d *pm, int v1, int v2, int v3) {
    int i, w1, w2, w3;
    for (i=0; i<pm->nF; i++) {
	w1 = pm->face[3*i+0];
	w2 = pm->face[3*i+1];
	w3 = pm->face[3*i+2];
	if ((v1==w1)&&(v2==w2)&&(v3==w3)) return i;
	if ((v2==w1)&&(v3==w2)&&(v1==w3)) return i;
	if ((v3==w1)&&(v1==w2)&&(v2==w3)) return i;
    }
    return -1;
}
static int addface(delaunay3d *pm, int v1, int v2, int v3) {
#ifdef STRUCTCHECK
    int i;
    for (i=0; i<pm->nF; i++) {
	if ((v1!=pm->face[3*i+0])&&(v1!=pm->face[3*i+1])&&(v1!=pm->face[3*i+2])) continue;
	if ((v2!=pm->face[3*i+0])&&(v2!=pm->face[3*i+1])&&(v2!=pm->face[3*i+2])) continue;
	if ((v3!=pm->face[3*i+0])&&(v3!=pm->face[3*i+1])&&(v3!=pm->face[3*i+2])) continue;
	libaft_3d_warn("delaunay.c: addface(): duplicate faces detected (%d %d %d)  (internal error)", v1, v2, v3);
    }
#endif
    if (pm->nF >= pm->nnF) {
	libaft_3d_stop("Delaunay: maxFace exceeded");
	return -1;
    } else {
	pm->face[3*(pm->nF)+0] = v1;
	pm->face[3*(pm->nF)+1] = v2;
	pm->face[3*(pm->nF)+2] = v3;
	return (pm->nF)++;
    }
}
static void remface(delaunay3d *pm, int i) {
    if (pm->npack > pm->nnF) libaft_3d_stop("delaunay.c: Internal error in remface & packface");
    pm->pack[pm->npack] = i;
    pm->npack++;
}
static void packface(delaunay3d *pm) {
    int i, j;
    while (pm->npack > 0) {
	pm->npack--;
	j = pm->pack[pm->npack];
	pm->nF--;
	pm->face[3*j+0] = pm->face[3*(pm->nF)+0];
	pm->face[3*j+1] = pm->face[3*(pm->nF)+1];
	pm->face[3*j+2] = pm->face[3*(pm->nF)+2];
	for (i=0; i<pm->npack; i++) if (pm->pack[i] == pm->nF) pm->pack[i] = j;
    }
}

static int superhex(delaunay3d *pm) {
    REAL bmin[3], bmax[3];
    REAL x, y, z, r, R;
    int t, i, j, v1, v2, v3, v4, v5, v6, v7, v8;

    for (j=0; j<3; j++) {
	bmin[j] = pm->vertex[3*0+j];
	bmax[j] = pm->vertex[3*0+j];
    }
    for (i=1; i<pm->nVF; i++) {
	for (j=0; j<3; j++) {
	    if (bmin[j] > pm->vertex[3*i+j])  bmin[j] = pm->vertex[3*i+j];
	    if (bmax[j] < pm->vertex[3*i+j])  bmax[j] = pm->vertex[3*i+j];
	}
    }
    x = (bmin[0] + bmax[0])/2.0;
    y = (bmin[1] + bmax[1])/2.0;
    z = (bmin[2] + bmax[2])/2.0;
    r = dist3(bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2])/2.0;
    R = 2.0;
    v1 = addvertex(pm, x - R*r, y - R*r, z - R*r);
    v2 = addvertex(pm, x - R*r, y + R*r, z - R*r);
    v3 = addvertex(pm, x + R*r, y + R*r, z - R*r);
    v4 = addvertex(pm, x + R*r, y - R*r, z - R*r);
    v5 = addvertex(pm, x - R*r, y - R*r, z + R*r);
    v6 = addvertex(pm, x - R*r, y + R*r, z + R*r);
    v7 = addvertex(pm, x + R*r, y + R*r, z + R*r);
    v8 = addvertex(pm, x + R*r, y - R*r, z + R*r);
    t = addtetra(pm, v1, v2, v3, v6),  tetrametric(pm, t);
    t = addtetra(pm, v1, v3, v7, v6),  tetrametric(pm, t);
    t = addtetra(pm, v5, v7, v6, v1),  tetrametric(pm, t);
    t = addtetra(pm, v1, v3, v4, v8),  tetrametric(pm, t);
    t = addtetra(pm, v1, v7, v3, v8),  tetrametric(pm, t);
    t = addtetra(pm, v5, v8, v7, v1),  tetrametric(pm, t);
    return t;
}

static int supertetra(delaunay3d *pm) {
    REAL bmin[3], bmax[3];
    REAL x, y, z, r, R;
    int t, i, j, v1, v2, v3, v4;

    for (j=0; j<3; j++) {
	bmin[j] = pm->vertex[3*0+j];
	bmax[j] = pm->vertex[3*0+j];
    }
    for (i=1; i<pm->nVF; i++) {
	for (j=0; j<3; j++) {
	    if (bmin[j] > pm->vertex[3*i+j])  bmin[j] = pm->vertex[3*i+j];
	    if (bmax[j] < pm->vertex[3*i+j])  bmax[j] = pm->vertex[3*i+j];
	}
    }
    x = (bmin[0] + bmax[0])/2.0;
    y = (bmin[1] + bmax[1])/2.0;
    z = (bmin[2] + bmax[2])/2.0;
    r = dist3(bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2])/2.0;
    R = 2.0;
    v1 = addvertex(pm, x + 4.0*R*r, y - 1.0*R*r, z - 1.0*R*r);
    v2 = addvertex(pm, x - 1.0*R*r, y + 4.0*R*r, z - 1.0*R*r);
    v3 = addvertex(pm, x - 1.0*R*r, y - 1.0*R*r, z + 4.0*R*r);
    v4 = addvertex(pm, x - 1.0*R*r, y - 1.0*R*r, z - 1.0*R*r);
    t = addtetra(pm, v1, v2, v3, v4);
    tetrametric(pm, t);
    return t;
}

static int proface(delaunay3d *pm, int v1, int v2, int v3) {
    int i;
    i = findface(pm, v3, v2, v1);
    if (i<0)  i = addface(pm, v1, v2, v3);
    else  remface(pm, i);
    return i;
}

static int checkfrt(delaunay3d *pm, int v0) {
    int v1, v2, v3, k, j, w1, w2, w3, w4, m;
    int good = 0, r = 0;
    while (!good) {
	good = 1;
	for (j=0; j<pm->nF; j++) {
	    v1 = pm->face[3*j+0],  v2 = pm->face[3*j+1],  v3 = pm->face[3*j+2];
	    if (idet3i4(pm->vertex, v0, v1, v2, v3) <= 0) {
		good = 0;
//		printf("point %d: face: %d, %d, %d. tet: ", v0, v1, v2, v3);
//		fflush(stdout);
		m = 0;
		for (k=0; k<pm->nT; k++) {
		    w1 = pm->tetra[4*k+0],  w2 = pm->tetra[4*k+1],  w3 = pm->tetra[4*k+2],  w4 = pm->tetra[4*k+3];
		    if (    ((v1==w1) || (v1==w2) || (v1==w3) || (v1==w4)) &&
			    ((v2==w1) || (v2==w2) || (v2==w3) || (v2==w4)) &&
			    ((v3==w1) || (v3==w2) || (v3==w3) || (v3==w4))    ) {
//			printf((m)?", %d":"%d", k);
			proface(pm, w2, w3, w4);
			proface(pm, w1, w4, w3);
			proface(pm, w1, w3, w2);
			proface(pm, w1, w2, w4);
			remtetra(pm, k);
			packface(pm);
			m++;
		    }
		}
//		printf(" (%d total)%s\n", m, (m==1)?"":" check this!");
		if (m==0)  return libaft_3d_warn("delaunay.c: checkfrt(): Oops, infinite loop detected! (internal warning)"),  -1;
		r++;
		break;
	    }
	}
    }
    packtetra(pm);
    return r;
}

static int incremental(delaunay3d *pm) {
    int i, j, k, v0, v1, v2, v3;

    for (i=0; i<pm->nVF; i++) {
	for (k=0; k<pm->nT; k++) {
	    if (    (intetra(pm->vertex, &(pm->tetra[4*k]), i)) ||
		    (insphere(&(pm->tetrasphere[k]), pm->vertex+3*i))    ) {
		v0 = pm->tetra[4*k+0],  v1 = pm->tetra[4*k+1],  v2 = pm->tetra[4*k+2],  v3 = pm->tetra[4*k+3];
		proface(pm, v1, v2, v3);
		proface(pm, v0, v3, v2);
		proface(pm, v0, v2, v1);
		proface(pm, v0, v1, v3);
		remtetra(pm, k);
		packface(pm);
	    }
	}
	packtetra(pm);

	checkfrt(pm, i);

	if (0) savefrtraw(pm);
	v0 = i;
	for (j=0; j<pm->nF; j++) {
	    v1 = pm->face[3*j+0],  v2 = pm->face[3*j+1],  v3 = pm->face[3*j+2];
	    k = addtetra(pm, v0, v1, v2, v3);
	    tetrametric(pm, k);
	}
	pm->nF = 0;
#ifdef SHOWPROGRESS
	printf("delaunay3D: nV = %5d, nT = %5d, %6.2lf%%\r", i+1, pm->nT, 100.0*(i+1)/pm->nVF),  fflush(stdout);
#endif
    }
    return 0;
}

static int cleanup(delaunay3d *pm) { // not used right now
    int n, k, v0, v1, v2, v3;

    n = pm->nVF;

    for (k=0; k<pm->nT; k++) {
	v0 = pm->tetra[4*k+0],  v1 = pm->tetra[4*k+1],  v2 = pm->tetra[4*k+2],  v3 = pm->tetra[4*k+3];
	if ( (v0>=n) || (v1>=n) || (v2>=n) || (v3>=n) ) remtetra(pm, k);
    }
    packtetra(pm);
    return 0;
}

static int refix(delaunay3d *pm) { // not used right now
    int i, j, k, v0, v1, v2, v3, rr, i1, i2, i3, p;

    pm->nF = 0;
    for (k=0; k<pm->nT; k++) {
	v0 = pm->tetra[4*k+0],  v1 = pm->tetra[4*k+1],  v2 = pm->tetra[4*k+2],  v3 = pm->tetra[4*k+3];
	proface(pm, v3, v2, v1);
	proface(pm, v0, v2, v3);
	proface(pm, v0, v1, v2);
	proface(pm, v0, v3, v1);
	packface(pm);
    }

    checktopology(pm->nV, pm->vertex, pm->nF, pm->face, pm->nT, pm->tetra, 0);

    rr = 0;
    while (rr) {
	rr = 0;
	for (i=0; i<3*pm->nF; i++) {
	    v0 = pm->face[i];
	    for (j=0; j<pm->nF; j++) {
		v1 = pm->face[3*j+0],  v2 = pm->face[3*j+1],  v3 = pm->face[3*j+2];
		if ((v0==v1) || (v0==v2) || (v0==v3))  continue;
		if (idet3i4(pm->vertex, v0, v1, v2, v3) > 0) {
		    p = 0;
		    i1 = findface(pm, v0, v3, v2);  if (i1>=0)  p++;
		    i2 = findface(pm, v0, v1, v3);  if (i2>=0)  p++;
		    i3 = findface(pm, v0, v2, v1);  if (i3>=0)  p++;

		    if (p<1)  continue;

		    printf("p = %d \n", p);

		    if (i1<0)  addface(pm, v0, v2, v3);  else remface(pm, i1);
		    if (i2<0)  addface(pm, v0, v3, v1);  else remface(pm, i2);
		    if (i3<0)  addface(pm, v0, v1, v2);  else remface(pm, i3);
		    remface(pm, j);
		    packface(pm);
		    addtetra(pm, v0, v1, v2, v3);
//		    printf("new tetra\n");
		    rr = 1;
		    break;
		}
	    }
	    if (rr) break;
	}
    }
    return 0;
}

static int checkvrt(delaunay3d *pm) {
    int *vc;
    int i, k, r=0;

    vc = libaft_malloc(sizeof(int) * pm->nV);
    for (i=0; i<pm->nV; i++)  vc[i] = 0;
    for (k=0; k<pm->nT; k++) {
	vc[pm->tetra[4*k+0]]++;
	vc[pm->tetra[4*k+1]]++;
	vc[pm->tetra[4*k+2]]++;
	vc[pm->tetra[4*k+3]]++;
    }
    for (i=0; i<pm->nV; i++)  if (vc[i]==0)  printf((r)?", %d":"\nDelaunay internal check.\nUnreferenced points: %d", i),  r++;
    if (r)  printf(" (%d total).\n", r);
    libaft_free(vc);
    return r;
}

int detstatsprint(char *buf);

#define MIN(x,y) (((x)<(y))?(x):(y))
/* main function */
/* NOT strictly Delaunay, and with extra points */
int delaunay3d_main(int *pnV, REAL *vertex, int *pnF, int *face, int *pnT, int *tetra, int nnV, int nnF, int nnT) {
    int i, r;
    delaunay3d mesh;
    char buf[65520];

    detinit();

    if (*pnF > 0)  libaft_3d_warn("Delaunay: nF shoud be zero in delaunay3d_main");

    mesh.nV = *pnV,  mesh.nT = *pnT,  mesh.nF = *pnF;
    mesh.nnV = MIN(nnV, mesh.nV+8);
    mesh.nnT = MIN(nnT, (mesh.nV+8)*(mesh.nV+8)*2 + 16);
    mesh.nnF = MIN(nnF, 3*nnT);
    mesh.vertex = vertex,  mesh.face = face,  mesh.tetra = tetra;

    mesh.tetrametric = libaft_malloc(sizeof(REAL) * mesh.nnT*TMW);
    mesh.tetrasphere = libaft_malloc(sizeof(sphere) * mesh.nnT);
    mesh.tpack       = libaft_malloc(sizeof(int)  * mesh.nnT);
    mesh.pack        = libaft_malloc(sizeof(int)  * mesh.nnF);

    mesh.ntpack = 0;
    mesh.npack = 0;
    mesh.nVF = mesh.nV;

    if (mesh.tetrametric && mesh.tetrasphere && mesh.tpack && mesh.face && mesh.pack) {

	//for (i=0; i<mesh.nT*4; i++) tetra[i]--;

	for (i=0; i<mesh.nT; i++) tetrametric(&mesh, i);

	if (TETRA)  supertetra(&mesh);  else  superhex(&mesh);
	incremental(&mesh);
	if (0)  cleanup(&mesh);
#ifdef SHOWPROGRESS
	printf("delaunay3D: nV = %5d, nT = %5d, fix... \r", mesh.nV, mesh.nT),  fflush(stdout);
#endif
	if (0)  refix(&mesh);

	if (0) {
	    detstatsprint(buf);
	    printf("%s\n", buf);
	}

//	mesh.nV = mesh.nVF;

#ifdef SHOWPROGRESS
	printf("delaunay3D: nV = %5d, nT = %5d, chk... \r", mesh.nV, mesh.nT),  fflush(stdout);
#endif

	r = checkvrt(&mesh);

	if (TETRA) {
	    mesh.nF=4;
	    mesh.face[0*3+0] = mesh.nVF+3,  mesh.face[0*3+1] = mesh.nVF+2,  mesh.face[0*3+2] = mesh.nVF+1;
	    mesh.face[1*3+0] = mesh.nVF+0,  mesh.face[1*3+1] = mesh.nVF+2,  mesh.face[1*3+2] = mesh.nVF+3;
	    mesh.face[2*3+0] = mesh.nVF+0,  mesh.face[2*3+1] = mesh.nVF+3,  mesh.face[2*3+2] = mesh.nVF+1;
	    mesh.face[3*3+0] = mesh.nVF+0,  mesh.face[3*3+1] = mesh.nVF+1,  mesh.face[3*3+2] = mesh.nVF+2;
	} else {
	    mesh.nF=12;
	    mesh.face[0*3+0] = mesh.nVF+0,  mesh.face[0*3+1] = mesh.nVF+1,  mesh.face[0*3+2] = mesh.nVF+2;
	    mesh.face[1*3+0] = mesh.nVF+0,  mesh.face[1*3+1] = mesh.nVF+2,  mesh.face[1*3+2] = mesh.nVF+3;
	    mesh.face[2*3+0] = mesh.nVF+4,  mesh.face[2*3+1] = mesh.nVF+6,  mesh.face[2*3+2] = mesh.nVF+5;
	    mesh.face[3*3+0] = mesh.nVF+4,  mesh.face[3*3+1] = mesh.nVF+7,  mesh.face[3*3+2] = mesh.nVF+6;
	    mesh.face[4*3+0] = mesh.nVF+0,  mesh.face[4*3+1] = mesh.nVF+5,  mesh.face[4*3+2] = mesh.nVF+1;
	    mesh.face[5*3+0] = mesh.nVF+0,  mesh.face[5*3+1] = mesh.nVF+4,  mesh.face[5*3+2] = mesh.nVF+5;
	    mesh.face[6*3+0] = mesh.nVF+7,  mesh.face[6*3+1] = mesh.nVF+3,  mesh.face[6*3+2] = mesh.nVF+2;
	    mesh.face[7*3+0] = mesh.nVF+6,  mesh.face[7*3+1] = mesh.nVF+7,  mesh.face[7*3+2] = mesh.nVF+2;
	    mesh.face[8*3+0] = mesh.nVF+5,  mesh.face[8*3+1] = mesh.nVF+2,  mesh.face[8*3+2] = mesh.nVF+1;
	    mesh.face[9*3+0] = mesh.nVF+6,  mesh.face[9*3+1] = mesh.nVF+2,  mesh.face[9*3+2] = mesh.nVF+5;
	    mesh.face[10*3+0] = mesh.nVF+3,  mesh.face[10*3+1] = mesh.nVF+7,  mesh.face[10*3+2] = mesh.nVF+0;
	    mesh.face[11*3+0] = mesh.nVF+7,  mesh.face[11*3+1] = mesh.nVF+4,  mesh.face[11*3+2] = mesh.nVF+0;
	}

	r += checktopology(mesh.nV, mesh.vertex, mesh.nF, mesh.face, mesh.nT, mesh.tetra, 0);

	//for (i=0; i<mesh.nT*4; i++) tetra[i]++;

#ifdef SHOWPROGRESS
	printf("delaunay3D: nV = %5d, nT = %5d,  done. \r", mesh.nV, mesh.nT),  fflush(stdout);
#endif
	r = 0;
    } else r = -1;

    libaft_free(mesh.tetrametric);
    libaft_free(mesh.tetrasphere);
    libaft_free(mesh.tpack);
    libaft_free(mesh.pack);

    *pnV = mesh.nV,  *pnT = mesh.nT,  *pnF = mesh.nF;
    return r;
}












static int savefrtraw(delaunay3d *pm) {
    FILE *f;
    int i;
    static int fn=0;
    char buf[1024];

    sprintf(buf, "frt%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "%d %d\n", pm->nV, pm->nF);
    for (i=0; i<pm->nV; i++) {
	fprintf(f, "%lf %lf %lf\n", pm->vertex[3*i+0], pm->vertex[3*i+1], pm->vertex[3*i+2]);
    }
    for (i=0; i<pm->nF; i++) {
	fprintf(f, "%d %d %d\n", pm->face[3*i+0]+1, pm->face[3*i+1]+1, pm->face[3*i+2]+1);
    }
    fclose(f);
    return fn++;
}
