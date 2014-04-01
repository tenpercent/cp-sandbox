/* main_prm.c
 *
 * This is example for libfrtprm and libaft usage
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "libaft.h"
#include "libfrtprm.h"
#include "fix.h"
#include "Ani3D_FC.h"

#define LAYERS 1

extern int tria_dump_front, tria_debug_front, region_dump_face; // debug flags from aftlib, could be usefull for debugging

void saveMani();

static double size = 11.0;
static double skinsize = 3.0;
static double cx[180], dux[180], dvx[180], rx[60];
static int nx = 0; /* max 60 */

static double YYSIZE = 100.0;


static double *vertexref = NULL;
static int nvref = 0;

static void addvref(double x, double y, double z, double s) {
    vertexref[4*nvref+0] = x,  vertexref[4*nvref+1] = y,  vertexref[4*nvref+2] = z,  vertexref[4*nvref+3] = s*1.25,  nvref++;
}

/* this function controls the desired size of the mesh in space */
/* x, y, z  -- coords of point in space                         */
static double fsize(double x, double y, double z) {
    double mr = 1e16, r;
    int i;

    for (i=0; i<nvref; i++) {
	r = (vertexref[4*i+0]-x)*(vertexref[4*i+0]-x) + 
	    (vertexref[4*i+1]-y)*(vertexref[4*i+1]-y) +
	    (vertexref[4*i+2]-z)*(vertexref[4*i+2]-z);
	r = sqrt(r);
	r -= vertexref[4*i+3];
	if (mr > r)  mr = r;
    }
    if (mr < 0.0)  mr = 0.0;
    r = sqrt(size*size/9.0 + mr*mr/16.0);
    if (r > size)  r = size;

    return r;
}

/* Helper functions for local minimization */
static double minimum(double (*f)(double), double x1, double x2) {
    /* 30 832040
     * 31 1346269
     * 32 2178309 */
    int n;
    double a, b, y1, y2;

    n = 32;
    a = x1,  b = x2;
    x1 = a + (b-a)*832040.0/2178309.0;
    x2 = a + (b-a)*1346269.0/2178309.0;
    y1 = f(x1),  y2 = f(x2);
    while (n > 1) {
	n--;
	if (y1 > y2)  a=x1,  x1=x2,  x2=b-(x1-a),  y1=y2,  y2=f(x2);
	else          b=x2,  x2=x1,  x1=a+(b-x2),  y2=y1,  y1=f(x1);
    }
    return  (x2+x1)/2.0;
}
static double guessx, guessy, guessr;
static double guessdist(double alpha) {
    double x = 230.0*cos(alpha)*guessr;
    double y = YYSIZE*sin(alpha)*guessr;
    return  (x-guessx)*(x-guessx) + (y-guessy)*(y-guessy);
}
static double getv(double x, double y, double ag) {
    guessx = x,  guessy = y,  guessr = 1.0;
    return minimum(guessdist, ag - M_PI/2.0, ag + M_PI/2.0);
}
static double staminx, staminy, staminz, staminu, stamingv, staminr;
static double localfuncv(double alpha) {
    double x = 230.0*cos(alpha)*staminr;
    double y = YYSIZE*sin(alpha)*staminr;
    return  (x-staminx)*(x-staminx) + (y-staminy)*(y-staminy);
}
static double localfuncu(double u) {
    staminr = (fabs(u) > 300.0) ? 1.0 - ((fabs(u)-300.0)/70.0)*((fabs(u)-300.0)/70.0) : 1.0;
    if (staminr < 0.0)  staminr = 0.0;
    staminr = sqrt(staminr);
    return (u-staminz)*(u-staminz) + localfuncv(stamingv = minimum(localfuncv, stamingv - M_PI/2.0, stamingv + M_PI/2.0));
}
static void proj(double x, double y, double z, double *u, double *v) {
    staminx = x,  staminy = y, staminz = z;
    stamingv = atan2(y, x);
    u[0] = minimum(localfuncu, z - 100.0, z + 100.0);
    staminr = (fabs(u[0]) > 300.0) ? 1.0 - ((fabs(u[0])-300.0)/70.0)*((fabs(u[0])-300.0)/70.0) : 1.0;
    if (staminr < 0.0)  staminr = 0.0;
    staminr = sqrt(staminr);
    v[0] = minimum(localfuncv, stamingv - M_PI/2.0, stamingv + M_PI/2.0);
}
static void arc(int i, double alpha, double *u, double *v) {
    double x, y, z;
    x = cx[3*i+0] + dux[3*i+0]*cos(alpha)*rx[i] + dvx[3*i+0] * sin(alpha)*rx[i];
    y = cx[3*i+1] + dux[3*i+1]*cos(alpha)*rx[i] + dvx[3*i+1] * sin(alpha)*rx[i];
    z = cx[3*i+2] + dux[3*i+2]*cos(alpha)*rx[i] + dvx[3*i+2] * sin(alpha)*rx[i];
    proj(x, y, z, u, v);
}

/* Helper functions for intersections */
static double xrc(double z) {return  110 - 50*pow(fabs(z)/280, 1.5);}
static double ar (double z) {return  90*pow(fabs(290-z)/280, 0.25);}
static double br (double z) {return  70*pow(fabs(290-z)/280, 0.25);}
static double xlc(double z) {return  -(110 - 50*pow(fabs(z)/280, 1.5));}
static double al (double z) {return  90*pow(fabs(290-z)/280, 0.25);}
static double bl (double z) {return  70*pow(fabs(290-z)/280, 0.25);}

static int sign(double x) {
    if (x >  1e-6)  return  1;
    if (x < -1e-6)  return -1;
    return 0;
}
static double zero(double (*f)(double), double x0, double x1) {
    double a, b, c;
    int sa, sb, sc, i;

    a = x0,  b = x1;
    sa = sign(f(a)),  sb = sign(f(b));
    if (sa*sb > 0) {
	if (fabs(f(a)) < 1e-2)  return a;
	if (fabs(f(b)) < 1e-2)  return b;
	printf("zero: same sign\nx0 = %lf, x1 = %lf\nf0 = %lf, f1 = %lf\n", x0, x1, f(x0), f(x1));
	return (a+b)/2.0;
    }
    if (sa == 0)  return a;
    if (sb == 0)  return b;
    for (i=0; i<30; i++) {
	c = (a+b)/2.0;
	sc = sign(f(c));
	if (sc == 0)  return c;
	if (sc == sa)  a = c;
	if (sc == sb)  b = c;
    }
    return c;
}

static double fhrz(double z) {
    return (xrc(z) - ar(z)) - (sqrt(1.0 - ((z-60.0)/60.0)*((z-60.0)/60.0))*50.0 - 20.0);
}
static double fhlz(double z) {
    return (xlc(z) + al(z)) - (sqrt(1.0 - ((z-60.0)/60.0)*((z-60.0)/60.0))*50.0 - 20.0);
}
static double fz = 0.0;
static double fhrr(double a) {
    double x = ar(fz)*cos(a) + xrc(fz);
    double y = br(fz)*sin(a);
    return ((x+20.0)/50.0)*((x+20.0)/50.0) + (y/50.0)*(y/50.0) + ((fz-60.0)/60.0)*((fz-60.0)/60.0) - 1.0;
}
static double fhrh(double a) {
    double x = 50.0*cos(a)*sqrt(1.0 - ((fz-60.0)/60.0)*((fz-60.0)/60.0)) - 20.0;
    double y = 50.0*sin(a)*sqrt(1.0 - ((fz-60.0)/60.0)*((fz-60.0)/60.0));
    return ((x-xrc(fz))/ar(fz))*((x-xrc(fz))/ar(fz)) + (y/br(fz))*(y/br(fz)) - 1.0;
}
static double fhll(double a) {
    double x = al(fz)*cos(a) + xlc(fz);
    double y = bl(fz)*sin(a);
    return ((x+20.0)/50.0)*((x+20.0)/50.0) + (y/50.0)*(y/50.0) + ((fz-60.0)/60.0)*((fz-60.0)/60.0) - 1.0;
}
static double fhlh(double a) {
    double x = 50.0*cos(a)*sqrt(1.0 - ((fz-60.0)/60.0)*((fz-60.0)/60.0)) - 20.0;
    double y = 50.0*sin(a)*sqrt(1.0 - ((fz-60.0)/60.0)*((fz-60.0)/60.0));
    return ((x-xlc(fz))/al(fz))*((x-xlc(fz))/al(fz)) + (y/bl(fz))*(y/bl(fz)) - 1.0;
}

/* Surface parameterization function
 * 1, 2 – Heart
 * 3, 4 – Right lung
 * 5, 6 – Left lung
 */
static int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
    double x = 0.0, y = 0.0, z = 0.0, r;
    if (i == 1 || i == 2) {
	r = 1.0 - ((u-60.0)/60.0)*((u-60.0)/60.0);
	if (r < 0.0)  r = 0.0;
	r = sqrt(r);
	x = 50*cos(v)*r - 20;
	y = 50*sin(v)*r;
	z = u;
    } else if (i == 3 || i == 4) {
	z = u;
	x = ar(z)*cos(v) + xrc(z);
	y = br(z)*sin(v);
    } else if (i == 5 || i == 6) {
	z = u;
	x = al(z)*cos(v) + xlc(z);
	y = bl(z)*sin(v);
    } else if (i == 7 || i == 8) {
	z = u;
	x = 230.0*cos(v);
	y = YYSIZE*sin(v);
    } else if (i == 123) {
	r = (fabs(u) > 300.0) ? 1.0 - ((fabs(u)-300.0)/70.0)*((fabs(u)-300.0)/70.0) : 1.0;
	if (r < 0.0)  r = 0.0;
	r = sqrt(r);
	x = 230.0*cos(v)*r;
	y = YYSIZE*sin(v)*r;
	z = u;
    } else {
	printf("surface_param: wrong id\n");
    }
    *px = x;
    *py = y;
    *pz = z;
    return 1;
}
static double periodic(int i, int d) {
    if ((i==123) && (d == 1))  return  2.0*M_PI;
    return  0.0;
}

/* Line parameterization function */
static int line_param(int i, double t, double *pu, double *pv) {
    double u = 0.0,  v = 0.0;
    if (i <= 0)  printf("line_param: wrong id\n");
    else if (i==1)   u =      t,  v = -M_PI;
    else if (i==2)   u =      t,  v =   0.0;
    else if (i==3)   u =      t,  v =  M_PI;
    else if (i==4)   u =    0.0,  v =     t;
    else if (i==5)   u =  280.0,  v =     t;
    else if (i==6)   u =  300.0,  v =     t;
    else if (i==7)   u = -300.0,  v =     t;
    else if (i==11)  fz = t,  u = t,  v = zero(fhrr,   0.0,   M_PI);
    else if (i==12)  fz = t,  u = t,  v = zero(fhrr, -M_PI,    0.0);
    else if (i==13)  fz = t,  u = t,  v = zero(fhrh,   0.0,   M_PI);
    else if (i==14)  fz = t,  u = t,  v = zero(fhrh, -M_PI,    0.0);
    else if (i==15)  fz = t,  u = t,  v = zero(fhll,   0.0,   M_PI);
    else if (i==16)  fz = t,  u = t,  v = zero(fhll, -M_PI,    0.0);
    else if (i==17)  fz = t,  u = t,  v = zero(fhlh,   0.0,   M_PI);
    else if (i==18)  fz = t,  u = t,  v = zero(fhlh, -M_PI,    0.0);
    else if (i==19)  fz = t,  u = t,  v = zero(fhlh,  M_PI, 2*M_PI);
    else if (i>=100)  arc(i-100, t, &u, &v);
    *pu = u,  *pv = v;
    return 1;
}

typedef double pnt[3];
typedef double vec[3];
static double dot(vec x, vec y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}
static void mkvec(pnt a, pnt b, vec x) {
    x[0] = b[0] - a[0],  x[1] = b[1] - a[1],  x[2] = b[2] - a[2];
}
static void addvec(double a, vec x, double b, vec y) {
    x[0] = a*x[0] +b*y[0],  x[1] = a*x[1] + b*y[1],  x[2] = a*x[2] + b*y[2];
}
static double dist(pnt a, pnt b) {
    vec x;
    mkvec(a, b, x);
    return sqrt(dot(x,x));
}
static double det2x2(double a, double b, double c, double d) {
    return a*d - b*c;
}
static double vp(vec x, vec y, vec z) {
    z[0] = det2x2(x[1], x[2], y[1], y[2]);
    z[1] = det2x2(x[2], x[0], y[2], y[0]);
    z[2] = det2x2(x[0], x[1], y[0], y[1]);
    return sqrt(dot(z,z));
}
static int addprism(int *pnT, int *tetra, int *tetramaterial, int a, int b, int c, int nV, int color) {
    int nT = *pnT;
    int v1, v2, v3, d;
    if (a<b)  if (a<c)  v1 = a,  v2 = (b<c)?b:c,  v3 = (b<c)?c:b,  d = (b<c)?0:1;
    else  v1 = c,  v2 = a,  v3 = b,  d = 0;
    else if (b<c)  v1 = b,  v2 = (a<c)?a:c,  v3 = (a<c)?c:a,  d = (a<c)?1:0;
    else  v1 = c,  v2 = b,  v3 = a,  d = 1;
    tetra[4*nT+0] = v1,  tetra[4*nT+1] = nV+v1,  tetra[4*nT+2+d] = v2,  tetra[4*nT+3-d] = v3,  tetramaterial[nT] = color,  nT++;
    tetra[4*nT+0] = v2,  tetra[4*nT+1] = nV+v2,  tetra[4*nT+2+d] = v3,  tetra[4*nT+3-d] = nV+v1,  tetramaterial[nT] = color,  nT++;
    tetra[4*nT+0] = v3,  tetra[4*nT+1] = nV+v3,  tetra[4*nT+2+d] = nV+v1,  tetra[4*nT+3-d] = nV+v2,  tetramaterial[nT] = color,  nT++;
    *pnT = nT;
    return 0;
}
static int fix_vertices(int *pnV, double *vertex, int nF, int *face, int nT, int *tetra, int nF2, int *face2) {
    int i, nV = *pnV, m;
    int *rename;

    rename = (int*)malloc(sizeof(int)*(nV+1));
    for (i=1; i<=nV; i++)  rename[i] = -1;
    for (i=0; i<3*nF; i++)  rename[face[i]] = 1;
    for (i=0; i<3*nF2; i++)  rename[face2[i]] = 1;
    for (i=0; i<4*nT; i++)  rename[tetra[i]] = 1;
    m = 0;
    for (i=1; i<=nV; i++) {
	if (rename[i] > 0) {
	    rename[i] = ++m;
	    addvec(0.0, vertex+3*m, 1.0, vertex+3*i);
	}
    }
    nV = m;
    for (i=0; i<3*nF; i++)  face[i] = rename[face[i]];
    for (i=0; i<3*nF2; i++)  face2[i] = rename[face2[i]];
    for (i=0; i<4*nT; i++)  tetra[i] = rename[tetra[i]];
    *pnV = m;
    free(rename);
    return 0;
}

static int keepskin(
	int *pnF, int *face, int *facematerial
	) {
    int nF = *pnF;
    int removed = 0, i;
    for (i=0; i<nF; i++) {
	if ((facematerial[i] > 1) && (facematerial[i] <= 10)) {
	    removed++;
	    nF--;
	    face[3*i+0] = face[3*nF+0];
	    face[3*i+1] = face[3*nF+1];
	    face[3*i+2] = face[3*nF+2];
	    facematerial[i] = facematerial[nF];
	    i--;
	}
    }
    for (i=0; i<nF; i++)  if (facematerial[i] > 10)  facematerial[i] -= 9;
    *pnF = nF;
    return removed;
}


static int makeskin(
	double h,
	int *pnV, double *vertex,
	int nE, int *edge, int *edgematerial,
	int *pnF, int *face, int *facematerial,
	int *pnT, int *tetra, int *tetramaterial,
	int *pnFa, int *facea, int *facemata,
	int nnV, int nnF, int nnT
) {
    int nV = *pnV,  nF = *pnF,  nT = *pnT,  nFa = *pnFa;
    int removed = 0, i, j;
    vec x, a, b;
    double *mass;
    double m;
    int v1, v2, d;
    if (4*nV+4 > nnV)  {
	printf("Please increase max nV\n");
	return -1;
    }
    mass = (double*)malloc(sizeof(double)*3*(nV+1));
    for (i=1; i<=3*nV; i++)  vertex[3*(nV+i)+0] = 0.0,  vertex[3*(nV+i)+1] = 0.0,  vertex[3*(nV+i)+2] = 0.0,  mass[i] = 0.0;
    for (i=0; i<nF; i++) {
	if (facematerial[i] == 1 || facematerial[i] > 10) {
	    mkvec(vertex+3*face[3*i+2], vertex+3*face[3*i+1], a);
	    mkvec(vertex+3*face[3*i+2], vertex+3*face[3*i+0], b);
	    m = vp(a, b, x);
	    for (j=0; j<3; j++)  addvec(1.0, vertex+3*(1*nV+face[3*i+j]), 1.0, x),  mass[0*nV+face[3*i+j]] += m;
	    for (j=0; j<3; j++)  addvec(1.0, vertex+3*(2*nV+face[3*i+j]), 1.0, x),  mass[1*nV+face[3*i+j]] += m;
	    for (j=0; j<3; j++)  addvec(1.0, vertex+3*(3*nV+face[3*i+j]), 1.0, x),  mass[2*nV+face[3*i+j]] += m;
	}
    }
    for (i=1; i<=nV; i++) {
	if (mass[i] > 0.0)
	    addvec(1.0*h/mass[i], vertex+3*(nV+i), 1.0, vertex+3*i),
	    addvec(2.0*h/mass[i], vertex+3*(2*nV+i), 1.0, vertex+3*i),
	    addvec(3.0*h/mass[i], vertex+3*(3*nV+i), 1.0, vertex+3*i);
    }
    for (i=0; i<nF; i++) {
	if (facematerial[i] == 1) {
	    addprism(&nT, tetra, tetramaterial, face[3*i+0], face[3*i+1], face[3*i+2], nV, 1);
	    facea[3*nFa+0] = nV+face[3*i+2],  facea[3*nFa+1] = nV+face[3*i+1],  facea[3*nFa+2] = nV+face[3*i+0],  facemata[nFa] = facematerial[i],  nFa++;
	    removed++;
	    nF--;
	    if (nF>0) {
		face[3*i+0] = face[3*nF+0];
		face[3*i+1] = face[3*nF+1];
		face[3*i+2] = face[3*nF+2];
		facematerial[i] = facematerial[nF];
	    }
	    i--;
	}
	if (facematerial[i] > 10) {
	    addprism(&nT, tetra, tetramaterial, face[3*i+0], face[3*i+1], face[3*i+2], nV, 1);
	    addprism(&nT, tetra, tetramaterial, nV+face[3*i+0], nV+face[3*i+1], nV+face[3*i+2], nV, 7);
	    if (LAYERS==2)  addprism(&nT, tetra, tetramaterial, 2*nV+face[3*i+0], 2*nV+face[3*i+1], 2*nV+face[3*i+2], nV, 8);
	    facea[3*nFa+0] = 2*nV+face[3*i+2],  facea[3*nFa+1] = 2*nV+face[3*i+1],  facea[3*nFa+2] = 2*nV+face[3*i+0],  facemata[nFa] = facematerial[i],  nFa++;
	    if (LAYERS==2)  facea[3*nFa+0] = 3*nV+face[3*i+2],  facea[3*nFa+1] = 3*nV+face[3*i+1],  facea[3*nFa+2] = 3*nV+face[3*i+0],  facemata[nFa] = facematerial[i],  nFa++;
	    removed++;
	    nF--;
	    if (nF>0) {
		face[3*i+0] = face[3*nF+0];
		face[3*i+1] = face[3*nF+1];
		face[3*i+2] = face[3*nF+2];
		facematerial[i] = facematerial[nF];
	    }
	    i--;
	}
    }
    for (i=0; i<nE; i++) {
	v1 = (edge[2*i+0] < edge[2*i+1]) ? edge[2*i+0] : edge[2*i+1];
	v2 = (edge[2*i+0] > edge[2*i+1]) ? edge[2*i+0] : edge[2*i+1];
	d = (edge[2*i+0] < edge[2*i+1]) ? 1 : 0;

	facea[3*nFa+d] = nV+v1,  facea[3*nFa+1-d] = 2*nV+v1,  facea[3*nFa+2] = nV+v2,  facemata[nFa] = 1,  nFa++;
	facea[3*nFa+d] = 2*nV+v2,  facea[3*nFa+1-d] = nV+v2,  facea[3*nFa+2] = 2*nV+v1,  facemata[nFa] = 1,  nFa++;
	facea[3*nFa+d] = 2*nV+v1,  facea[3*nFa+1-d] = 3*nV+v1,  facea[3*nFa+2] = 2*nV+v2,  facemata[nFa] = edgematerial[i],  nFa++;
	facea[3*nFa+d] = 3*nV+v2,  facea[3*nFa+1-d] = 2*nV+v2,  facea[3*nFa+2] = 3*nV+v1,  facemata[nFa] = edgematerial[i],  nFa++;
    }
    nV += 3*nV;
    free(mass);
    *pnV = nV;
    *pnF = nF;
    *pnT = nT;
    *pnFa = nFa;
    return removed;
    (void) nnF,  (void) nnT;
    (void) dist;
}

/* Helper functions */
static int add_edge(int *LineD, int *LineP, double *LineT, int *pnLine, int *pm, int v1, int v2, int n, ...) {
    va_list ap;
    int i;

    va_start(ap, n);
    LineD[3**pnLine+0] = v1,  LineD[3**pnLine+1] = v2,  LineD[3**pnLine+2] = n,  (*pnLine)++;
    for (i=0; i<n; i++)  LineP[2**pm+0] = va_arg(ap, int),  LineP[2**pm+1] = va_arg(ap, int),   LineT[2**pm+0] = va_arg(ap, double),  LineT[2**pm+1] = va_arg(ap, double),  (*pm)++;
    return *pnLine;
}
static int add_surf(int *SurfL, double *SurfT, int *SurfI, int *pnSurface, int *pm, int param, int c1, int c2, int dir, double umin, double umax, double vmin, double vmax, int n, ...) {
    va_list ap;
    int i;

    va_start(ap, n);
    SurfL[5**pnSurface+0] = n,  SurfL[5**pnSurface+1] = param,  SurfL[5**pnSurface+2] = c1,  SurfL[5**pnSurface+3] = c2,  SurfL[5**pnSurface+4] = dir;
    SurfT[4**pnSurface+0] = umin,  SurfT[4**pnSurface+1] = umax,  SurfT[4**pnSurface+2] = vmin,  SurfT[4**pnSurface+3] = vmax,  (*pnSurface)++;
    for (i=0; i<n; i++) {
	SurfI[2**pm+0] = va_arg(ap, int),  SurfI[2**pm+1] = va_arg(ap, int),  (*pm)++;
    }
    return *pnSurface;
}

/* Workaround 7=13 points */
static void workaround713(int *pnF, int *face, int *facematerial) {
    int nF = *pnF;
    int i;
    for (i=0; i<nF; i++) {
	if (((face[3*i+0]==7)||(face[3*i+1]==7)||(face[3*i+2]==7)) &&
		((face[3*i+0]==13)||(face[3*i+1]==13)||(face[3*i+2]==13))) {
	    nF--;
	    if (nF>0) {
		face[3*i+0] = face[3*nF+0];
		face[3*i+1] = face[3*nF+1];
		face[3*i+2] = face[3*nF+2];
		facematerial[i] = facematerial[nF];
	    }
	    i--;
	}
	if (face[3*i+0]==13)  face[3*i+0] = 7;
	if (face[3*i+1]==13)  face[3*i+1] = 7;
	if (face[3*i+2]==13)  face[3*i+2] = 7;
    }
    *pnF = nF;
}

static void addcirc(double *x, double size) {
    double u, v, r, t, q;

    if (nx>=60)  return;

    cx[3*nx+0] = x[0],  cx[3*nx+1] = x[1],  cx[3*nx+2] = x[2];
    rx[nx] = size;

    proj(x[0], x[1], x[2], &u, &v);
    t = (u - (fabs(u+300.0)-fabs(u-300.0))/2.0)/70.0;
    if ((q = 1.0 - t*t) < 0.0)  q = 0.0;
    q = sqrt(q);

//    printf("t = %lf, q = %lf\n", t, q);

    dux[3*nx+0] = -230.0*cos(v)*t,  dux[3*nx+1] = -YYSIZE*sin(v)*t,  dux[3*nx+2] = 70.0*q;
    r = dux[3*nx+0]*dux[3*nx+0] + dux[3*nx+1]*dux[3*nx+1] + dux[3*nx+2]*dux[3*nx+2],  r = sqrt(r);
    dux[3*nx+0] /= r,  dux[3*nx+1] /= r,  dux[3*nx+2] /=r;

    dvx[3*nx+0] = -230.0*sin(v)*q,  dvx[3*nx+1] = YYSIZE*cos(v)*q,  dvx[3*nx+2] = 0.0;
    r = dvx[3*nx+0]*dvx[3*nx+0] + dvx[3*nx+1]*dvx[3*nx+1] + dvx[3*nx+2]*dvx[3*nx+2],  r = sqrt(r);
    dvx[3*nx+0] /= r,  dvx[3*nx+1] /= r,  dvx[3*nx+2] /=r;

//    printf("d_u: %lf, %lf, %lf\n", dux[3*nx+0], dux[3*nx+1], dux[3*nx+2]);
//    printf("d_v: %lf, %lf, %lf\n", dvx[3*nx+0], dvx[3*nx+1], dvx[3*nx+2]);

    nx++;
}

static void addthreecirc(double *x, double size1, double size2, double size3) {
    if (nx>=57)  return;
    addcirc(x, size1);
    addcirc(x, size2);
    addcirc(x, size3);
}

static double restorex(char x, double y, double z) {
    double r;
    if ((z = fabs(z)-300.0) < 0.0)  z = 0.0;
    if ((r = 1.0 - (y/YYSIZE)*(y/YYSIZE) - (z/70.0)*(z/70.0)) < 0.0)  r = 0.0;
    return ((x=='-') ? -1 : 1)*230.0*sqrt(r);
}
static double restorey(double x, char y, double z) {
    double r;
    if ((z = fabs(z)-300.0) < 0.0)  z = 0.0;
    if ((r = 1.0 - (x/230.0)*(x/230.0) - (z/70.0)*(z/70.0)) < 0.0)  r = 0.0;
    return ((y=='-') ? -1 : 1)*YYSIZE*sqrt(r);
}
static double restorez(double x, double y, char z) {
    double r;
    if ((r = 1.0 - (x/230.0)*(x/230.0) - (y/YYSIZE)*(y/YYSIZE)) < 0.0)  r = 0.0;
    return ((z=='-') ? -1 : 1)*(300+70.0*sqrt(r));
}
/* */
#define BFS 16383
int readparams(int argc, char *argv[]) {
    char *defname="mesh.txt";
    char *fname = (argc>1)?argv[1]:defname;
    double p[3], s, s2, sw, x, y, z;
    char cx[BFS+2], cy[BFS+2], cz[BFS+2];
    FILE *f;
    char buf[BFS+1];
    int i, n;

    if (!(f = fopen(fname, "r"))) {
	perror(fname);
	return -1;
    }
    printf(" * Reading parameters from %s\n", fname);
    if (!fgets(buf, BFS, f))  return  perror(fname),  1;
    sscanf(buf, "%d %lf %lf", &n, &size, &skinsize);

    for (i=0; i<n; i++) {
	if (!fgets(buf, BFS, f))  return  perror(fname),  1;
	s = 90.0,  sw  = 5.0,  s2 = 10.0;
	if (sscanf(buf, " %lf %lf %lf %lf %lf %lf ", &x, &y, &z, &s, &sw, &s2) >= 6);
	else if (sscanf(buf, " %[+-] %lf %lf %lf %lf %lf ", cx, &y, &z, &s, &sw, &s2) >= 6)  x = restorex(cx[0], y, z);
	else if (sscanf(buf, " %lf %[+-] %lf %lf %lf %lf ", &x, cy, &z, &s, &sw, &s2) >= 6)  y = restorey(x, cy[0], z);
	else if (sscanf(buf, " %lf %lf %[+-] %lf %lf %lf ", &x, &y, cz, &s, &sw, &s2) >= 6)  z = restorez(x, y, cz[0]);
	else {
	    printf("Cannot parse string:\n%s", buf);
	    continue;
	}
	p[0] = x,  p[1] = y,  p[2] = z;
	printf("Electrode: (%lg, %lg, %lg), diams = %lg (%lg), %lg\n", p[0], p[1], p[2], s, sw, s2);
	addthreecirc(p, s/2.0, s/2.0 - sw, s2/2.0);
	addvref(x, y, z, s/2.0);
    }
    fclose(f);
    nx /= 3;
    printf("Number of electrodes: %d\n", nx);
    return 0;
}

static int colorbottom(int nV, double *vertex, int nT, int *tetra, int *tetramaterial) {
    int i, r = 0;
    double x, y, z;

    for (i=0; i<nT; i++) {
	if (tetramaterial[i] != 5)  continue;
	x = (vertex[3*tetra[4*i+0]+0] + vertex[3*tetra[4*i+1]+0] + vertex[3*tetra[4*i+2]+0] + vertex[3*tetra[4*i+3]+0]) / 4.0;
	y = (vertex[3*tetra[4*i+0]+1] + vertex[3*tetra[4*i+1]+1] + vertex[3*tetra[4*i+2]+1] + vertex[3*tetra[4*i+3]+1]) / 4.0;
	z = (vertex[3*tetra[4*i+0]+2] + vertex[3*tetra[4*i+1]+2] + vertex[3*tetra[4*i+2]+2] + vertex[3*tetra[4*i+3]+2]) / 4.0;
	if (z < 0.0)  tetramaterial[i] = 6,  r++;
    }
    return r;
}

static int recolor(int nT, int *tetramaterial) {
    int i, remap[10]={0,1,8,6,7,4,5,2,3,9};
    for (i=0; i<nT; i++)  tetramaterial[i] = remap[tetramaterial[i]];
    return 0;
}

/* Main */
int main(int argc, char* argv[]) {
    int    nnF = 2000000, nnV = 2000000, nnT = 11000000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0, *facedup = 0, *facematdup = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    i, nFdup, r, j, medge, msurf;
    int    nVVert,    nLine,  nSurface;
    int    *LineD,    *LineP;
    double *LineT;
    int    *SurfL,    *SurfI;
    double *SurfT;
    double *VVert;
    int *expCrv;
    int nE = 0, nnE = 1000000, *edge, *edgematerial;

    double zhr0, zhr1, zhl0, zhl1;
    int ehl, ehrhb, ehrht, ehr1, ehr2, ehrrb, ehrrt, ero, erb1, erb2, ert1, ert2, ehl1, ehl2, ehllb, ehllt, elo, elb1, elb2, elt1, elt2;
    int ebz1, ebz2, ebt1, ebt2, ebb1, ebb2;
    int e1, e2, e3, e4;
    int mx[256], md[256];

    int ntfix = 0,  nvfix = 0;

    int izero = 0;

    vertexref  = (double*)malloc(sizeof(double) * 4 * nnV);

    if (readparams(argc, argv))  return  0;

    // allocate memory for mesh ctructures
    vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
    face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial  = (int*)   malloc(sizeof(int)        * nnF);
    facematdup    = (int*)   malloc(sizeof(int)        * nnF);
    tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
    tetramaterial = (int*)   malloc(sizeof(int)        * nnT);
    edge          = (int*)   malloc(sizeof(int)    * 2 * nnE);
    edgematerial  = (int*)   malloc(sizeof(int)        * nnE);

    // allocate memory for boundary representation structure
    VVert = (double*)malloc(sizeof(double) * 3*100);
    nVVert = 0;
    LineD = (int*)   malloc(sizeof(int)    * 3*100);
    LineP = (int*)   malloc(sizeof(int)    * 2*2*100);
    LineT = (double*)malloc(sizeof(double) * 2*2*100);
    nLine = 0;
    SurfL = (int*)   malloc(sizeof(int)    * 5*100);
    SurfI = (int*)   malloc(sizeof(int)    * 2*2*100);
    SurfT = (double*)malloc(sizeof(double) * 4*100);
    nSurface = 0;

    expCrv = (int*)   malloc(sizeof(int)    * 100);
    memset(expCrv, 0, sizeof(int)    * 100);


    /***/
#define ADD_VERTEX(X, Y, Z);  { VVert[3*nVVert+0] = X,  VVert[3*nVVert+1] = Y,  VVert[3*nVVert+2] = Z,  nVVert++;}
#define EDGE(V1, V2, N, ...)  add_edge(LineD, LineP, LineT, &nLine, &medge, V1, V2, N, __VA_ARGS__)
#define SURF(P, C1, C2, D, U1, U2, V1, V2, N, ...); add_surf(SurfL, SurfT, SurfI, &nSurface, &msurf, P, C1, C2, D, U1, U2, V1, V2, N, __VA_ARGS__);
    /***/

    ADD_VERTEX(-230.0, -YYSIZE, -370.0);
    ADD_VERTEX( 230.0,  YYSIZE,  370.0);
    for (i=2; i<23; i++) {
	ADD_VERTEX( (double)i, 1.0,  1.0);
    }
    for (i=0; i<6*nx; i++) {
	ADD_VERTEX( 0.1,  i, 0.1);
    }

    medge = 0;

    zhr0 = zero(fhrz,  0.0,  60.0);
    zhr1 = zero(fhrz, 60.0, 120.0);

    zhl0 = zero(fhlz,  0.0,  60.0);
    zhl1 = zero(fhlz, 60.0, 120.0);

    /* EDGE(v1, v2, nSurf,  [iSurf, iLine, t_0, t1,] ...); */

    ehl   = EDGE( 2,  1,  2,   1,  3,  120.0,    0.0,   2,  1,  120.0,    0.0);
    ero   = EDGE( 3,  4,  2,   3,  2,    0.0,  280.0,   4,  2,    0.0,  280.0);
    ert1  = EDGE( 4,  5,  1,   3,  5,    0.0,   M_PI);                 
    erb1  = EDGE( 6,  3,  1,   3,  4,   M_PI,    0.0);                 
    ert2  = EDGE( 4,  5,  1,   4,  5,    0.0,  -M_PI);                 
    erb2  = EDGE( 6,  3,  1,   4,  4,  -M_PI,    0.0);                 
    ehrrt = EDGE( 5, 12,  2,   3,  3,  280.0,   zhr1,   4,  1,  280.0,   zhr1);
    ehrrb = EDGE(11,  6,  2,   3,  3,   zhr0,    0.0,   4,  1,   zhr0,    0.0);
    ehrhb = EDGE(13, 11,  2,   1,  2,   zhl0,   zhr0,   2,  2,   zhl0,   zhr0);
    ehrht = EDGE(12, 14,  2,   1,  2,   zhr1,   zhl1,   2,  2,   zhr1,   zhl1);
    ehr1  = EDGE(11, 12,  2,   3, 11,   zhr0,   zhr1,   1, 13,   zhr0,   zhr1);
    ehr2  = EDGE(11, 12,  3,   4, 12,   zhr0,   zhr1,   2, 14,   zhr0,   zhr1,   1, 14, zhr0,  zhr1);
    ehl1  = EDGE(13, 14,  2,   5, 15,   zhl0,   zhl1,   1, 17,   zhl0,   zhl1);
    ehl2  = EDGE(13, 14,  3,   6, 16,   zhl0,   zhl1,   2, 18,   zhl0,   zhl1,   1, 19, zhl0,  zhl1);
    ehllt = EDGE(14,  8,  2,   5,  2,   zhl1,  280.0,   6,  2,   zhl1,  280.0);
    ehllb = EDGE( 7, 13,  2,   5,  2,    0.0,   zhl0,   6,  2,    0.0,   zhl0);
    elo   = EDGE( 9, 10,  2,   5,  3,  280.0,    0.0,   6,  1,  280.0,    0.0);
    elt1  = EDGE( 8,  9,  1,   5,  5,    0.0,   M_PI);                 
    elb1  = EDGE(10,  7,  1,   5,  4,   M_PI,    0.0);                 
    elt2  = EDGE( 8,  9,  1,   6,  5,    0.0,  -M_PI);                 
    elb2  = EDGE(10,  7,  1,   6,  4,  -M_PI,    0.0);                 
                                                                       
    ebz1  = EDGE(15, 17,  1,   123,  4,    0.0,   M_PI);
    ebz2  = EDGE(17, 15,  1,   123,  4,  -M_PI,    0.0);


    msurf = 0;

    /* SURF(iSurf, color1, color2, direction, u0, u1, v0, v1, n,  [edge, edge_dir,] ...); */

    SURF( 1, 2, 5*1, 1,    0.0,  120.0,   0.0,   M_PI, 4,  ehrhb,1, ehr1,1, ehrht,1, ehl1,0);
    SURF( 2, 2, 5*1, 1,    0.0,  120.0, -M_PI,    0.0, 4,  ehrhb,0, ehr2,0, ehrht,0, ehl2,1);

    SURF( 3, 3, 5*1, 1,    0.0,  280.0,   0.0,   M_PI, 6,    ero,1, ert1,1, ehrrt,1, ehr1,0, ehrrb,1, erb1,1);
    SURF( 4, 3, 5*1, 1,    0.0,  280.0, -M_PI,    0.0, 6,    ero,0, ert2,0, ehrrt,0, ehr2,1, ehrrb,0, erb2,0);
    SURF( 0, 3, 5*1, 1,    0.0,    0.0,   0.0,    0.0, 2,   erb1,0, erb2,1);
    SURF( 0, 3, 5*1, 0,    0.0,    0.0,   0.0,    0.0, 2,   ert1,0, ert2,1);

    SURF( 1, 2, 3*1, 1,    0.0,  120.0, -M_PI,   M_PI, 2,   ehr1,0, ehr2,1);
                                            
    SURF( 5, 4, 5*1, 1,    0.0,  280.0,  -0.1,   M_PI, 6,  ehllb,1, ehl1,1, ehllt,1, elt1,1,   elo,1, elb1,1);
    SURF( 6, 4, 5*1, 1,    0.0,  280.0, -M_PI,    0.1, 6,  ehllb,0, ehl2,0, ehllt,0, elt2,0,   elo,0, elb2,0);
    SURF( 0, 4, 5*1, 1,    0.0,    0.0,   0.0,    0.0, 2,   elb1,0, elb2,1);
    SURF( 0, 4, 5*1, 0,    0.0,    0.0,   0.0,    0.0, 2,   elt1,0, elt2,1);

    SURF( 1, 2, 4*1, 1,    0.0,  120.0,   0.0, 2*M_PI, 2,   ehl1,1, ehl2,0);
                                            
/*    SURF( 0, 5, 6*1, 1,    0.0,    0.0,   0.0,    0.0, 6,   ebz1,1, ebz2,1,  elb1,1, elb2,0,  erb1,1, erb2,0);*/

    mx[0] = 0;

    for (i=0; i<nx; i++) {
	SURF(123, 5, 11+2*i, 1, -370.0,  370.0, -M_PI,   M_PI, 4, 
		e1 = EDGE(23+6*i, 24+6*i,  1,   123,  100+3*i, M_PI, 0.0),  0,
		e2 = EDGE(24+6*i, 23+6*i,  1,   123,  100+3*i, 0.0, -M_PI), 0,
		e3 = EDGE(25+6*i, 26+6*i,  1,   123,  101+3*i, M_PI, 0.0),  1,
		e4 = EDGE(26+6*i, 25+6*i,  1,   123,  101+3*i, 0.0, -M_PI), 1);
	mx[++mx[0]] = e1,  md[mx[0]] = 1;
	mx[++mx[0]] = e2,  md[mx[0]] = 1;
	mx[++mx[0]] = e3,  md[mx[0]] = 0;
	mx[++mx[0]] = e4,  md[mx[0]] = 0;
	expCrv[e1-1] = 1;
	expCrv[e2-1] = 1;
	expCrv[e3-1] = 1;
	expCrv[e4-1] = 1;
	SURF(123, 5, 12+2*i, 1, -370.0,  370.0, -M_PI,   M_PI, 2, 
		e1 = EDGE(27+6*i, 28+6*i,  1,   123,  102+3*i, M_PI, 0.0),  0,
		e2 = EDGE(28+6*i, 27+6*i,  1,   123,  102+3*i, 0.0, -M_PI), 0);
	mx[++mx[0]] = e1,  md[mx[0]] = 1;
	mx[++mx[0]] = e2,  md[mx[0]] = 1;
	expCrv[e1-1] = 1;
	expCrv[e2-1] = 1;
    }

    SurfL[5*nSurface+0] = (mx[0])? mx[0] : 4,  SurfL[5*nSurface+1] = 123,  SurfL[5*nSurface+2] = 5,  SurfL[5*nSurface+3] = 1,  SurfL[5*nSurface+4] = 1;
    SurfT[4*nSurface+0] = -370.0,  SurfT[4*nSurface+1] = 370.0,  SurfT[4*nSurface+2] = -M_PI,  SurfT[4*nSurface+3] = M_PI,  nSurface++;
    if (mx[0]) {
	for (i=1; i<=mx[0]; i++) {
	    SurfI[2*msurf+0] = mx[i],  SurfI[2*msurf+1] = md[i],  msurf++;
	}
    } else {
	SurfI[2*msurf+0] = ebz1,  SurfI[2*msurf+1] = 0,  msurf++;
	SurfI[2*msurf+0] = ebz2,  SurfI[2*msurf+1] = 0,  msurf++;
	SurfI[2*msurf+0] = ebz1,  SurfI[2*msurf+1] = 1,  msurf++;
	SurfI[2*msurf+0] = ebz2,  SurfI[2*msurf+1] = 1,  msurf++;
    }

//    tria_dump_front = 1;
//    tria_debug_front = 1;

    printf("\n * Generating surface mesh\n");
    i = ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
	    NULL, surface_param, line_param, periodic, fsize,
	    expCrv,
	    &nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,
	    &nnV, &nnF, &nnE
	    );
    free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);

    workaround713(&nF, face, facematerial);

    printf("INFO: nV = %d, nE = %d, nF = %d, nT = %d\n", nV, nE, nF, nT);

    /* Generate skin mesh */
    printf("\n * Generating skin mesh\n");
    nFdup = 0;
    makeskin(skinsize, &nV, vertex-3, nE, edge, edgematerial, &nF, face, facematerial, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, nnV, nnF, nnT);
    fix_vertices(&nV, vertex-3, nF, face, nT, tetra, nFdup, facedup);
    printf("INFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

    if (0) {
	write_mesh_gmv("surf.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial); // for GMV
//	write_front   ("surf.smv", nV, vertex, nF, face, facematerial); // for smv
	//		return 0; // do not mesh the volume, just exit
//	write_mesh_gmv("dups.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial); // for GMV
//	write_front   ("dups.smv", nV, vertex, nFdup, facedup, facematdup); // for smv
	return 0;
    }
    // It could be usefull to dump the triangulation of the surface in case we want to check
    // that boundary representation is correct and represents the desired region
    if (0) {
	write_mesh_gmv("surf.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial); // for GMV
	write_front   ("surf.smv", nV, vertex, nF, face, facematerial); // for smv
	//		return 0; // do not mesh the volume, just exit
	write_mesh_gmv("dups.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial); // for GMV
	write_front   ("dups.smv", nV, vertex, nFdup, facedup, facematdup); // for smv
    }
    // We will copy the front, so that it could be used in output in future.

    ntfix = nT;
    nvfix = nV;

    for (i=0; i<nF; i++) {
	facedup[3*nFdup+0] = face[3*i+0];
	facedup[3*nFdup+1] = face[3*i+1];
	facedup[3*nFdup+2] = face[3*i+2];
	facematdup[nFdup] = facematerial[i];
	nFdup++;
    }

    // Generate 3D mesh using our own size function fsize()
    printf("\n * Generating volume mesh\n");
    r = mesh_3d_aft_func_(&nV, vertex, &nF, face, facematerial, &nT, tetra, tetramaterial, &nnV, &nnF, &nnT, fsize);
    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

    if (r) {
	write_mesh_gmv("fail.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial);
    } else {
	/* Check that 3D mesh corresponds with surface mesh */
	/*printf("Checking topology: "),  fflush(stdout);*/
	if (0 && check_mesh_topology_(&nV, vertex, &nFdup, facedup, &nT, tetra)) printf("FAILED!\n");
	else {
	    /*printf("ok.\n");*/

	    if (0) {
		write_mesh         ("bfix.out", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
		write_mesh_gmv_qual("bfix.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
	    }

	    /* Improve mesh quality */
	    printf("\n * Smoothing volume mesh\n");
	    fixshape(&nV, vertex, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, 0, ntfix, 0/*nFdup*/, nnV, nnT, nnF);

	    keepskin(&nFdup, facedup, facematdup);

	    colorbottom(nV, vertex-3, nT, tetra, tetramaterial);
	    recolor(nT, tetramaterial);

	    // Write output files
	    write_mesh_gmv_qual("mesh.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
	    /*write_mesh         ("mesh.out", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);*/
	    saveMani(&nV, &nFdup, &nT,
		    vertex, facedup, tetra, facematdup, tetramaterial,
		    &izero, &izero, &izero, NULL, NULL, NULL,
		    &izero, NULL, NULL, "mesh.ani");
	}
    }

    free(vertex);
    free(face);
    free(facedup);
    free(facematerial);
    free(facematdup);
    free(tetra);
    free(tetramaterial);
    free(vertexref);
    return 0;
}

