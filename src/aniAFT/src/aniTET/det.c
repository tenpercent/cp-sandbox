#include <stdio.h>
#include <math.h>
#include "aft.h"

#define NODETSTATS


static LONGREAL lldet2( LONGREAL a, LONGREAL b,
			LONGREAL c, LONGREAL d ) {
    return a*d - b*c;
}
static LONGREAL lldet3( LONGREAL a, LONGREAL b, LONGREAL c,
			LONGREAL d, LONGREAL e, LONGREAL f,
			LONGREAL g, LONGREAL h, LONGREAL i ) {
    return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
}
static LONGREAL lldet4( LONGREAL a, LONGREAL b, LONGREAL c, LONGREAL d,
			LONGREAL e, LONGREAL f, LONGREAL g, LONGREAL h,
			LONGREAL i, LONGREAL j, LONGREAL k, LONGREAL l,
			LONGREAL m, LONGREAL n, LONGREAL o, LONGREAL p ) {
    return (a*f-b*e)*(k*p-l*o) + (a*g-c*e)*(l*n-j*p) + (a*h-d*e)*(j*o-k*n)
	+  (c*f-b*g)*(l*m-i*p) + (b*h-d*f)*(k*m-i*o) + (c*h-d*g)*(i*n-k*m);
}


REAL orient2d( REAL *pa, REAL *pb, REAL *pc ) {
    LONGREAL a = pa[0]-pc[0],  b = pb[0]-pc[0];
    LONGREAL c = pa[1]-pc[1],  d = pb[1]-pc[1];
    return lldet2(a, b,  c, d);
}
REAL orient3d( REAL *pa, REAL *pb, REAL *pc, REAL *pd ) {
    LONGREAL a = pa[0]-pd[0],  b = pb[0]-pd[0],  c = pc[0]-pd[0];
    LONGREAL d = pa[1]-pd[1],  e = pb[1]-pd[1],  f = pc[1]-pd[1];
    LONGREAL g = pa[2]-pd[2],  h = pb[2]-pd[2],  i = pc[2]-pd[2];
    return lldet3(a, b, c,  d, e, f,  g, h, i);
}
REAL orient4d( REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe ) {
    LONGREAL a = pa[0]-pe[0],  b = pb[0]-pe[0],  c = pc[0]-pe[0],  d = pd[0]-pe[0];
    LONGREAL e = pa[1]-pe[1],  f = pb[1]-pe[1],  g = pc[1]-pe[1],  h = pd[1]-pe[1];
    LONGREAL i = pa[2]-pe[2],  j = pb[2]-pe[2],  k = pc[2]-pe[2],  l = pd[2]-pe[2];
    LONGREAL m = pa[3]-pe[3],  n = pb[3]-pe[3],  o = pc[3]-pe[3],  p = pd[3]-pe[3];
    return lldet4(a, b, c, d,  e, f, g, h,  i, j, k, l,  m, n, o, p);
}

#ifdef DETSTATS
#define DSSIZE 256
static int detstatslog2[DSSIZE+1];
#endif
static void exactinit(void) {
#ifdef DETSTATS
    int i;
    for (i=0; i<DSSIZE; i++) detstatslog2[i] = 0;
#endif
}

void detinit(void) {
    static int initialized=0;
    int saved, cword;
    (void)saved, (void)cword;
    if (initialized) return;

    exactinit();
    initialized=1;
}

#ifdef DETSTATS
int detstatsprint(char *buf) {
    char *pc=buf;
    int i;
    pc += sprintf(pc, "Det stats:\n");
    for (i=0; i<DSSIZE; i++) {
	pc += sprintf(pc, "%2d: %d\n", i, detstatslog2[i]);
    }
    pc += sprintf(pc, " x: %d\n", detstatslog2[DSSIZE]);
    exactinit();
    return pc-buf;
}
#else
int detstatsprint(char *buf) {
    return sprintf(buf, "Det stats are disabled (define DETSTATS in det.c)\n");
}
#endif

REAL det2i3(REAL *vertex, int v1, int v2, int v3) {
    REAL r;
    r = orient2d(vertex+2*v1, vertex+2*v2, vertex+2*v3);
    return r;
}
int idet2i3(REAL *vertex, int v1, int v2, int v3) {
    REAL d = det2i3(vertex, v1, v2, v3);
    if (d>0.0) return +1;
    else if (d<0.0) return -1;
    else return 0;
}
REAL det3i4(REAL *vertex, int v1, int v2, int v3, int v4) {
    REAL r;
    r = orient3d(vertex+3*v1, vertex+3*v2, vertex+3*v3, vertex+3*v4);
    return r;
}
//#define EPS 6.938893903907228378e-18
//#define EPS 3.469446951953614189e-18
//
//#define P2 288230376151711744.0
#define P2 1125899906842624.0
//#define P2 1099511627776.0
//#define P2 1099511627776.0
int idet3i4(REAL *vertex, int v1, int v2, int v3, int v4) {
    REAL d = det3i4(vertex, v1, v2, v3, v4);
/*    REAL r = fabs(vertex[3*v1+0]) + fabs(vertex[3*v1+1]) + fabs(vertex[3*v1+2])
	+    fabs(vertex[3*v2+0]) + fabs(vertex[3*v2+1]) + fabs(vertex[3*v2+2])
	+    fabs(vertex[3*v3+0]) + fabs(vertex[3*v3+1]) + fabs(vertex[3*v3+2])
	+    fabs(vertex[3*v4+0]) + fabs(vertex[3*v4+1]) + fabs(vertex[3*v4+2]);*/
/*    if (d>0) return +1;
    else if (d<0) return -1;
    else return 0;*/
    REAL r = fabs(vertex[3*v1+0]-vertex[3*v4+0]) + fabs(vertex[3*v1+1]-vertex[3*v4+1]) + fabs(vertex[3*v1+2]-vertex[3*v4+2])
	+    fabs(vertex[3*v2+0]-vertex[3*v4+0]) + fabs(vertex[3*v2+1]-vertex[3*v4+1]) + fabs(vertex[3*v2+2]-vertex[3*v4+2])
	+    fabs(vertex[3*v3+0]-vertex[3*v4+0]) + fabs(vertex[3*v3+1]-vertex[3*v4+1]) + fabs(vertex[3*v3+2]-vertex[3*v4+2]);
    r *= r*r;
#ifdef DETSTATS
    REAL d2=fabs(d);
    int i=0;
    while ((d2<r) && (i<DSSIZE)) {
	i++;
	d2 *= 2.0;
    }
    detstatslog2[i]++;
#endif
    if (d*P2>r) return +1;
    else if (d*P2<-r) return -1;
    else return 0;
}
REAL det4i5(REAL *vertex, int v1, int v2, int v3, int v4, int v5) {
    REAL r;
    r = orient4d(vertex+4*v1, vertex+4*v2, vertex+4*v3, vertex+4*v4, vertex+4*v5);
    return r;
}
int idet4i5(REAL *vertex, int v1, int v2, int v3, int v4, int v5) {
    REAL d = det4i5(vertex, v1, v2, v3, v4, v5);
    if (d>0.0) return +1;
    else if (d<0.0) return -1;
    else return 0;
}

REAL det3(REAL *v1, REAL *v2, REAL *v3) {
    return lldet3(v1[0], v1[1], v1[2],  v2[0], v2[1], v2[2],  v3[0], v3[1], v3[2]);
}


void normvec2(REAL *v1, REAL *v2, REAL *v) {
    v[0] = lldet2(v1[1], v2[1], v1[2], v2[2]);
    v[1] = lldet2(v1[2], v2[2], v1[0], v2[0]);
    v[2] = lldet2(v1[0], v2[0], v1[1], v2[1]);
}

void normvec3i(REAL *vertex, int v1, int v2, int v3, REAL *v) {
    REAL w1[4], w2[4], w3[4];
    w1[0] = vertex[3*v1+0], w2[0] = vertex[3*v2+0], w3[0] = vertex[3*v3+0];
    w1[1] = vertex[3*v1+1], w2[1] = vertex[3*v2+1], w3[1] = vertex[3*v3+1];
    w1[2] = vertex[3*v1+2], w2[2] = vertex[3*v2+2], w3[2] = vertex[3*v3+2];
    w1[3] = vertex[3*v1+0], w2[3] = vertex[3*v2+0], w3[3] = vertex[3*v3+0];
    v[0] = orient2d(w3+1, w2+1, w1+1);
    v[1] = orient2d(w3+2, w2+2, w1+2);
    v[2] = orient2d(w3+0, w2+0, w1+0);
}
