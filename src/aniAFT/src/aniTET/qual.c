#include <stdio.h>
#include <math.h>
#include "aft.h"
#include "helper.h"
#include "det.h"

typedef int face[3];
typedef struct {
    face f;
    int n;
} flist;
typedef struct {
    int p;
    int n;
} plist;
typedef struct {
    int *i;
    face *f;
} neigh_faces;
typedef struct {
    int *i;
    int *p;
} neigh_points;

static int add_glist(int a, int b, int *png, plist *glist, int *s) {
    int c = s[a];
    while (c>=0) {
        if (glist[c].p == b)  return 0;
        c = glist[c].n;
    }
    glist[*png].p = b;
    glist[*png].n = s[a];
    s[a] = *png;
    (*png)++;
    return 1;
}

static int add_hlist(int a, face b, int *pnh, flist *hlist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (hlist[c].f[0] == b[0] && hlist[c].f[1] == b[1] && hlist[c].f[2] == b[2])  return libaft_3d_warn("hlist: entry already in list"), 0;
	c = hlist[c].n;
    }
    hlist[*pnh].f[0] = b[0];
    hlist[*pnh].f[1] = b[1];
    hlist[*pnh].f[2] = b[2];
    hlist[*pnh].n = s[a];
    s[a] = *pnh;
    (*pnh)++;
    return 1;
}

static int fill_tadj(neigh_faces *tadj, int nv, int nt, int *tetra) {
    int *s;
    int nh;
    flist *hlist;
    int i, j, n, l;
    int a;
    face b;

    hlist = (flist*)libaft_malloc(sizeof(flist)*4*nt);
    nh = 0;
    s = (int*)libaft_malloc(sizeof(int)*nv);
    for (j=0; j<nv; j++)  s[j] = -1;

    for (i=0; i<nt; i++) {
        a = tetra[4*i+0],  b[0] = tetra[4*i+1],  b[1] = tetra[4*i+2],  b[2] = tetra[4*i+3],  add_hlist(a, b, &nh, hlist, s);
        a = tetra[4*i+1],  b[0] = tetra[4*i+0],  b[1] = tetra[4*i+3],  b[2] = tetra[4*i+2],  add_hlist(a, b, &nh, hlist, s);
        a = tetra[4*i+2],  b[0] = tetra[4*i+3],  b[1] = tetra[4*i+0],  b[2] = tetra[4*i+1],  add_hlist(a, b, &nh, hlist, s);
        a = tetra[4*i+3],  b[0] = tetra[4*i+2],  b[1] = tetra[4*i+1],  b[2] = tetra[4*i+0],  add_hlist(a, b, &nh, hlist, s);
    }
    n = 0;
    for (j=0; j<nv; j++) {
        tadj->i[j] = n;
        l = s[j];
        while (l>=0) {
            tadj->f[n][0] = hlist[l].f[0];
            tadj->f[n][1] = hlist[l].f[1];
            tadj->f[n][2] = hlist[l].f[2];
            n++;
            l = hlist[l].n;
        }
    }
    tadj->i[nv] = n;

    libaft_free(s),  libaft_free(hlist);
    return n;
}

static int fill_eadj(neigh_points *eadj, int nv, int nt, int *tetra) {
    int *s;
    int ng;
    plist *glist;
    int i, j, n, l;
    int a, b, c, d;

    glist = (plist*)libaft_malloc(sizeof(plist)*12*nt);
    ng = 0;
    s = (int*)libaft_malloc(sizeof(int)*nv);
    for (j=0; j<nv; j++)  s[j] = -1;

    for (i=0; i<nt; i++) {
        a = tetra[4*i+0],  b = tetra[4*i+1],  c = tetra[4*i+2],  d = tetra[4*i+3];
        add_glist(a, b, &ng, glist, s),  add_glist(a, c, &ng, glist, s),  add_glist(a, d, &ng, glist, s);
        add_glist(b, a, &ng, glist, s),  add_glist(b, c, &ng, glist, s),  add_glist(b, d, &ng, glist, s);
        add_glist(c, a, &ng, glist, s),  add_glist(c, b, &ng, glist, s),  add_glist(c, d, &ng, glist, s);
        add_glist(d, a, &ng, glist, s),  add_glist(d, b, &ng, glist, s),  add_glist(d, c, &ng, glist, s);
    }
    n = 0;
    for (j=0; j<nv; j++) {
        eadj->i[j] = n;
        l = s[j];
        while (l>=0) {
            eadj->p[n++] = glist[l].p;
            l = glist[l].n;
        }
    }
    eadj->i[nv] = n;

    libaft_free(s),  libaft_free(glist);
    return n;
}



static REAL dist(REAL x1, REAL y1, REAL z1,  REAL x2, REAL y2, REAL z2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}
static REAL disti(REAL *vertex, int v1, int v2) {
    return dist(vertex[3*v1+0], vertex[3*v1+1], vertex[3*v1+2], vertex[3*v2+0], vertex[3*v2+1], vertex[3*v2+2]);
}
static REAL ldist(REAL x1, REAL y1, REAL z1,  REAL x2, REAL y2, REAL z2) {
    REAL s = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    return s*sqrt(s);
}
static REAL ldisti(REAL *vertex, int v1, int v2) {
    return ldist(vertex[3*v1+0], vertex[3*v1+1], vertex[3*v1+2], vertex[3*v2+0], vertex[3*v2+1], vertex[3*v2+2]);
}

#define D(a,b) (disti(vertex, a, b))

REAL aft_tetra_qual_mba(REAL *vertex, int v1, int v2, int v3, int v4) {
    REAL v, ls, q;
    v = det3i4(vertex, v1, v2, v3, v4);
    ls = D(v1,v2) + D(v1,v3) + D(v1,v4) + D(v2,v3) + D(v2,v4) + D(v3,v4);
    q = 216.0*sqrt(2.0)*v/ls/ls/ls;
    return q;
}

#define L(a,b) (ldisti(vertex, a, b))

REAL aft_tetra_qual(REAL *vertex, int v1, int v2, int v3, int v4) {
    REAL v, ls, q;
    v = det3i4(vertex, v1, v2, v3, v4);
    ls = L(v1,v2) + L(v1,v3) + L(v1,v4) + L(v2,v3) + L(v2,v4) + L(v3,v4);
    q = 6.0*sqrt(2.0)*v/ls;
    return q;
}
REAL aft_tetra_qual_delta(REAL *vertex, int v1, int v2, int v3, int v4, REAL delta) {
    REAL v, ls, q;
    v = det3i4(vertex, v1, v2, v3, v4);
    v = (v + sqrt(v*v + 4.0*delta*delta))/2.0;
    ls = L(v1,v2) + L(v1,v3) + L(v1,v4) + L(v2,v3) + L(v2,v4) + L(v3,v4);
    q = 6.0*sqrt(2.0)*v/ls;
    return q;
}

int print_aft_stats_linear(REAL *vertex, int nt, int *tetra) {
    REAL q, worst, qs;
    int stat[11], i;
    if (nt<1)  return 0;
    for (i=0; i<11; i++)  stat[i] = 0;
    worst = 1.0,  qs = 0.0;
    for (i=0; i<nt; i++) {
	q = aft_tetra_qual(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
	if (q < worst)  worst = q;
	if (q <= 0.0)  stat[10]++;
	else if (q>1.0)  stat[9]++;
	else stat[(int)ceil(q*10)-1]++;
	qs += q;
    }
    qs /= nt;
    printf("QS: q_avg: %12.8le  ", qs);
    printf(" ## %6d # ", stat[10]);
    for (i=0; i<10; i++)  printf(" %6d ", stat[i]);
    printf("\n");
    printf("    q_min: %12.8le  ", worst);
    printf(" ## %6.2lf # ", 100.0*stat[10]/nt);
    for (i=0; i<10; i++)  printf(" %6.2lf ", 100.0*stat[i]/nt);
    printf("\n");
    return 0;
}

int print_aft_stats_log(REAL *vertex, int nt, int *tetra) {
    REAL q, worst, qs;
    int stat[11], i, j;
    if (nt<1)  return 0;
    for (i=0; i<11; i++)  stat[i] = 0;
    worst = 1.0,  qs = 0.0;
    for (i=0; i<nt; i++) {
	qs += q = aft_tetra_qual(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
	if (q < worst)  worst = q;
	if (q<=0.0)  stat[10]++;
	else {
	    for (j=0; j<9; j++) {
		if (q > 0.1)  break;
		q *= 10.0;
	    }
	    stat[j]++;
	}
    }
    qs /= nt;
    printf("QS: q_avg: %15.8le  ", qs);
    printf(" ## %6d # ", stat[10]);
    for (i=0; i<10; i++)  printf(" %6d ", stat[i]);
    printf("\n");
    printf("ln  q_min: %15.8le  ", worst);
    printf(" ## %6.2lf # ", 100.0*stat[10]/nt);
    for (i=0; i<10; i++)  printf(" %6.2lf ", 100.0*stat[i]/nt);
    printf("\n");
    return 0;
}
int print_aft_stats(REAL *vertex, int nt, int *tetra) {
    return print_aft_stats_log(vertex, nt, tetra);
}

int print_aft_stats_linear_mba(REAL *vertex, int nt, int *tetra) {
    REAL q, worst, qs;
    int stat[11], i;
    if (nt<1)  return 0;
    for (i=0; i<11; i++)  stat[i] = 0;
    worst = 1.0,  qs = 0.0;
    for (i=0; i<nt; i++) {
	q = aft_tetra_qual_mba(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
	if (q < worst)  worst = q;
	if (q <= 0.0)  stat[10]++;
	else if (q>1.0)  stat[9]++;
	else stat[(int)ceil(q*10)-1]++;
	qs += q;
    }
    qs /= nt;
    printf("[mba] QS: q_avg: %12.8le  ", qs);
    printf(" ## %6d # ", stat[10]);
    for (i=0; i<10; i++)  printf(" %6d ", stat[i]);
    printf("\n");
    printf("[mba]     q_min: %12.8le  ", worst);
    printf(" ## %6.2lf # ", 100.0*stat[10]/nt);
    for (i=0; i<10; i++)  printf(" %6.2lf ", 100.0*stat[i]/nt);
    printf("\n");
    return 0;
}

int print_aft_stats_log_mba(REAL *vertex, int nt, int *tetra) {
    REAL q, worst, qs;
    int stat[11], i, j;
    if (nt<1)  return 0;
    for (i=0; i<11; i++)  stat[i] = 0;
    worst = 1.0,  qs = 0.0;
    for (i=0; i<nt; i++) {
	qs += q = aft_tetra_qual_mba(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
	if (q < worst)  worst = q;
	if (q<=0.0)  stat[10]++;
	else {
	    for (j=0; j<9; j++) {
		if (q > 0.1)  break;
		q *= 10.0;
	    }
	    stat[j]++;
	}
    }
    qs /= nt;
    printf("[mba] QS: q_avg: %15.8le  ", qs);
    printf(" ## %6d # ", stat[10]);
    for (i=0; i<10; i++)  printf(" %6d ", stat[i]);
    printf("\n");
    printf("[mba] ln  q_min: %15.8le  ", worst);
    printf(" ## %6.2lf # ", 100.0*stat[10]/nt);
    for (i=0; i<10; i++)  printf(" %6.2lf ", 100.0*stat[i]/nt);
    printf("\n");
    return 0;
}
int print_aft_stats_mba(REAL *vertex, int nt, int *tetra) {
    return print_aft_stats_log_mba(vertex, nt, tetra);
}

static LONGREAL qual_delta(REAL *vertex, int v1, int v2, int v3, int v4, REAL delta) {
    LONGREAL v, ls, q;
    v = det3i4(vertex, v1, v2, v3, v4);
    v = (v + sqrt(v*v + 4.0*delta*delta))/2.0;
    ls = L(v1,v2) + L(v1,v3) + L(v1,v4) + L(v2,v3) + L(v2,v4) + L(v3,v4);
    q = 6.0*sqrt(2.0)*v/ls;
    return q;
}
static LONGREAL qmin=1.0, qmax = 0.0, maxdf = 0.0;
static int func_xyz(REAL *vertex, int v, int v1, int v2, int v3, REAL dx[3], REAL delta, REAL ds, int n) {
    LONGREAL h = ds * 1e-6;
    LONGREAL q1, q2, f1, f2;
    int i, j;

    for (i=0; i<3; i++) {
	vertex[3*v + i] += h;
	q1 = qual_delta(vertex, v, v1, v2, v3, delta);
	vertex[3*v + i] -= h;
	for (f1 = 1.0, j = 0; j < n; j++)  f1 /= q1;
	vertex[3*v + i] -= h;
	q2 = qual_delta(vertex, v, v1, v2, v3, delta);
	vertex[3*v + i] += h;
	for (f2 = 1.0, j = 0; j < n; j++)  f2 /= q2;
//	if (fabs(f1-f2)/2.0/h > maxdf)  maxdf = fabs(f1-f2)/2.0/h;
	dx[i] = (f1 - f2)/(2.0*h);
//	if (qmin > q1)  qmin = q1;
//	if (qmin > q2)  qmin = q2;
//	if (qmax < q1)  qmax = q1;
//	if (qmax < q2)  qmax = q2;
//	printf("\t\t\tq: [%Lf, %Lf], df: %Le\r", qmin, qmax, maxdf);
    }
    return 0;
}

static REAL mins(REAL *vertex, int nt, int *tetra) {
    REAL s = 0.0, v;
    int i;

    for (i=0; i<nt; i++) {
	v = det3i4(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]);
	if ((i<1) || (s > v))  s = v;
    }
    return s;
}

static int refine_func(int nv, REAL *vertex, int nf, neigh_faces *tadj, neigh_points *eadj, int nt, int *tetra) {
    int j, l;
    int pn;
    int a, b, c, p;
    int n_iters, s;
    REAL ds, d, *dx;
    REAL dm, d1, d2, delta;
    REAL rs, r;
    REAL x[3], z[3];

    pn = 0;
    ds = 0.0;
    for (j=nf; j<nv; j++) {
        d = 0.0,  s = 0;
        for (l=eadj->i[j]; l<eadj->i[j+1]; l++) {
            c = eadj->p[l];
            d += disti(vertex, j, c);
            s++;
        }
        if (s > 0)  d /= s;
        ds += d;
	pn++;
    }
    if (pn == 0)  return 0;

    ds /= pn;

//    printf("ds = %lf\n", ds);

    dx = (REAL*)libaft_malloc(sizeof(REAL)*3*pn);

    delta = 0.1;
    n_iters = 1000; // 4000
    d = 0.05 * ds * ds / n_iters;  // 1.0
    dm = mins(vertex, nt, tetra);
    if (dm < 0.0) {
	d1 = 0.05 * -dm;
	d2 = 0.0001 * -dm;
    } else {
	d1 = 0.5 * dm;
	d2 = 0.001 * dm;
    }
//    printf("delta in [%le, %le]\n", d1, d2);
    for (s=0; s<n_iters; s++) {
//	dm = mins(vertex, nt, tetra);
//	if (dm < 0.0)  delta = 0.001*-dm;
//	else delta = 0.001*dm;
//	delta = (d1*(n_iters-s) + d2*s)/n_iters;
//	delta = 1e-6;
	delta = 1.0 * ds * ds * ds * (1.0 - 0.999*s/n_iters);
//	printf("delta = %lf\n", delta);
	rs = 0.0;
	for (j=0; j<pn; j++) {
	    p = nf+j;
	    x[0] = 0.0,  x[1] = 0.0,  x[2] = 0.0;
	    for (l=tadj->i[p]; l<tadj->i[p+1]; l++) {
		a = tadj->f[l][0];
		b = tadj->f[l][1];
		c = tadj->f[l][2];
		func_xyz(vertex, p, a, b, c, z, delta, ds, 6);
		x[0] -= z[0],  x[1] -= z[1],  x[2] -= z[2];
	    }
	    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	    if (r>1.0/ds*2) {
                r = 1.0/ds*2*(1.0 + log(r*ds/2))/r;
		x[0] *= r,  x[1] *= r,  x[2] *= r;
	    }
	    dx[3*j + 0] = x[0],  dx[3*j + 1] = x[1],  dx[3*j + 2] = x[2];
	    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	    if (rs < r)  rs = r;
	}
	for (j=0; j<pn; j++) {
	    p = nf+j;
	    vertex[3*p + 0] += d*dx[3*j + 0];
	    vertex[3*p + 1] += d*dx[3*j + 1];
	}
	printf("rs = %le\t[%6.2lf%%]\r", rs*d, 100.0*(s+1.0)/n_iters);
    }
//    printf("\n");

    qmin=1.0, qmax = 0.0, maxdf = 0.0;
    libaft_free(dx);

    return 0;
}


int refine3daftfunc(int *pnV, REAL *vertex, int *pnVF, int *pnT, int *tetra) {
    neigh_faces tadj;
    neigh_points eadj;
    REAL *bckup, r, rmax;
    int nV = *pnV, nVF = *pnVF, nT = *pnT;
    int k, i;

    bckup  = (REAL*)libaft_malloc(sizeof(REAL)*nV*3);
    for (i=0; i<3*nV; i++)  bckup[i] = vertex[i];

    tadj.i = (int* )libaft_malloc(sizeof(int )*(nV + 1));
    tadj.f = (face*)libaft_malloc(sizeof(face)*(4*nT));
    eadj.i = (int* )libaft_malloc(sizeof(int )*(nV + 1));
    eadj.p = (int* )libaft_malloc(sizeof(int )*(12*nT));

    fill_tadj(&tadj, nV, nT, tetra);
    fill_eadj(&eadj, nV, nT, tetra);

    for (k=0; k<1; k++)  {
	refine_func(nV, vertex, nVF, &tadj, &eadj, nT, tetra);
	rmax = 0.0;
	for (i=0; i<nV; i++) {
	    r =  (vertex[3*i+0]-bckup[3*i+0])*(vertex[3*i+0]-bckup[3*i+0]);
	    r += (vertex[3*i+1]-bckup[3*i+1])*(vertex[3*i+1]-bckup[3*i+1]);
	    r += (vertex[3*i+2]-bckup[3*i+2])*(vertex[3*i+2]-bckup[3*i+2]);
	    r = sqrt(r);
	    if (rmax < r)  rmax = r;
	}
	printf("dr = %le  \n", rmax);
	if (rmax < 1e2) {
	    for (i=0; i<3*nV; i++)  bckup[i] = vertex[i];
	}
	if (rmax < 1e-4) {
	    break;
	}
    }



    libaft_free(bckup);
    libaft_free(tadj.i),  libaft_free(tadj.f),  libaft_free(eadj.i),  libaft_free(eadj.p);

    return 0;
}
