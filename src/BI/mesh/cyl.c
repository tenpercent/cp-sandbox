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

static double size = 1.0;
static double skinsize = 1.0;
static double electrodesize = 0.5;

static double D = 10.0,  L = 60.0;

static double R;

static double h1 = 0.2,  h2 = 0.5;

static double sidelim = 30.0;
static double sideext = 30.0;


static int serie = 1;

/* this function controls the desired size of the mesh in space */
/* x, y, z  -- coords of point in space                         */
static double fsize(double x, double y, double z) {
    return size;
    (void) x,  (void) y,  (void) z;
}

/* Surface parameterization function
 */
static int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
    double x = 0.0, y = 0.0, z = 0.0;
    if (i == 1) {
	x = u,  y = R*cos(v),  z = R*sin(v);
    } else if (i == 2) {
	x = u,  y = 0.8*R*cos(v),  z = 0.8*R*sin(v);
    } else if (i == 3) {
	x = u,  y = 0.2*R*cos(v),  z = 0.2*R*sin(v);
    } else {
	printf("surface_param: wrong id\n");
    }
    *px = x,  *py = y,  *pz = z;
    return 1;
}
static double periodic(int i, int d) {
    if ((i==1) && (d == 1))  return  2.0*M_PI;
    if ((i==2) && (d == 1))  return  2.0*M_PI;
    if ((i==3) && (d == 1))  return  2.0*M_PI;
    return  0.0;
}

/* Line parameterization function */
static int line_param(int i, double t, double *pu, double *pv) {
    double u = 0.0,  v = 0.0;
    if (i==1)       u = -L/2.0,  v = t;
    else if (i==2)  u =  L/2.0,  v = t;
    else if (i==3)  u = -L/2.0+D*h1,  v = t;
    else if (i==4)  u =  L/2.0-D*h1,  v = t;
    else if (i==5)  u = -L/2.0+D*(h1+0.1),  v = t;
    else if (i==6)  u =  L/2.0-D*(h1+0.1),  v = t;
    else if (i==7)  u = -L/2.0+D*h2,  v = t;
    else if (i==8)  u =  L/2.0-D*h2,  v = t;
    else if (i==9)  u = -L/2.0+D*(h2+0.1),  v = t;
    else if (i==10)  u =  L/2.0-D*(h2+0.1),  v = t;
    else if (i==11)  u = t,  v = -0.1;
    else if (i==12)  u = t,  v =  0.1;
    else printf("line_param: wrong id\n");
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
static int fix_vertices(int *pnV, double *vertex, int nF, int *face, int nT, int *tetra, int nF2, int *face2, int nE, int *edge) {
    int i, nV = *pnV, m;
    int *rename;

    rename = (int*)malloc(sizeof(int)*(nV+1));
    for (i=1; i<=nV; i++)  rename[i] = -1;
    for (i=0; i<3*nF; i++)  rename[face[i]] = 1;
    for (i=0; i<3*nF2; i++)  rename[face2[i]] = 1;
    for (i=0; i<4*nT; i++)  rename[tetra[i]] = 1;
//    for (i=0; i<2*nE; i++)  rename[edge[i]] = 1;
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
    for (i=0; i<2*nE; i++)  edge[i] = rename[edge[i]];
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
    *pnF = nF;
    return removed;
}


static int addlayer(int type, double h,
	int *pnV, double *vertex,
	int nE, int *edge, int *edgematerial,
	int *pnF, int *face, int *facematerial,
	int *pnT, int *tetra, int *tetramaterial,
	int *pnFa, int *facea, int *facemata,
	int nnV, int nnF, int nnT
) {
    int nV = *pnV,  nF = *pnF,  nT = *pnT,  nFa = *pnFa;
    int removed = 0, i, j, k;
    vec x, a, b;
    double *mass;
    double m;
    int v1, v2, d, layers = ceil(h/size);
    int nFF = nF;
    h = h / layers;
    if ((1+layers)*nV+1+layers > nnV)  {
	printf("Please increase max nV\n");
	return -1;
    }
    mass = (double*)malloc(sizeof(double)*(nV+1));
    for (i=1; i<=layers*nV; i++)  vertex[3*(nV+i)+0] = 0.0,  vertex[3*(nV+i)+1] = 0.0,  vertex[3*(nV+i)+2] = 0.0;
    for (i=1; i<=nV; i++)  mass[i] = 0.0;
    for (i=0; i<nF; i++) {
	if ((type) ? (facematerial[i] > 10) : ((facematerial[i] == 1) || (facematerial[i] > 110))) {
	    mkvec(vertex+3*face[3*i+2], vertex+3*face[3*i+1], a);
	    mkvec(vertex+3*face[3*i+2], vertex+3*face[3*i+0], b);
	    m = vp(a, b, x);
	    for (j=0; j<3; j++)  addvec(1.0, vertex+3*(nV+face[3*i+j]), 1.0, x),  mass[face[3*i+j]] += m;
	}
    }
    for (i=1; i<=nV; i++) {
	if (mass[i] > 0.0) {
	    for (k=layers; k>1; k--)  addvec(0.0, vertex+3*(k*nV+i), 1.0, vertex+3*(nV+i));
	    for (k=layers; k>0; k--)  addvec(k*h/mass[i], vertex+3*(k*nV+i), 1.0, vertex+3*i);
	}
    }
    for (i=0; i<nFF; i++) {
	if ((type) ? ((facematerial[i] == 1) || (facematerial[i] > 10)) : ((facematerial[i] == 1) || (facematerial[i] > 100))) {
	    j = ((type==0) || (facematerial[i] > 10)) ? layers : 0;
	    if (facematerial[i]==101)  j = 0;
	    for (k=0; k<j; k++)  addprism(&nT, tetra, tetramaterial, k*nV+face[3*i+0], k*nV+face[3*i+1], k*nV+face[3*i+2], nV, (type) ? 7 : 1);
	    if (type==0)  face[3*nF+0] = j*nV+face[3*i+0],  face[3*nF+1] = j*nV+face[3*i+1],  face[3*nF+2] = j*nV+face[3*i+2],  facematerial[nF] = facematerial[i] % 100,  nF++;
	    else  facea[3*nFa+0] = j*nV+face[3*i+2],  facea[3*nFa+1] = j*nV+face[3*i+1],  facea[3*nFa+2] = j*nV+face[3*i+0],  facemata[nFa] = facematerial[i] % 100,  nFa++;
	    removed++;
	    nF--, nFF--;
	    if (nFF>0) {
		face[3*i+0] = face[3*nFF+0];
		face[3*i+1] = face[3*nFF+1];
		face[3*i+2] = face[3*nFF+2];
		facematerial[i] = facematerial[nFF];
		face[3*nFF+0] = face[3*nF+0];
		face[3*nFF+1] = face[3*nF+1];
		face[3*nFF+2] = face[3*nF+2];
		facematerial[nFF] = facematerial[nF];
	    }
	    i--;
	}
    }
    for (i=0; i<nE; i++) {
	v1 = (edge[2*i+0] < edge[2*i+1]) ? edge[2*i+0] : edge[2*i+1];
	v2 = (edge[2*i+0] > edge[2*i+1]) ? edge[2*i+0] : edge[2*i+1];
	d = (edge[2*i+0] < edge[2*i+1]) ? 1 : 0;
	if (type==0) {
	    if (edgematerial[i] % 100)  for (k=0; k<layers; k++) {
		face[3*nF+d] = k*nV+v1,  face[3*nF+1-d] = (k+1)*nV+v1,  face[3*nF+2] = k*nV+v2,  facematerial[nF] = edgematerial[i] % 100,  nF++;
		face[3*nF+d] = (k+1)*nV+v2,  face[3*nF+1-d] = k*nV+v2,  face[3*nF+2] = (k+1)*nV+v1,  facematerial[nF] = edgematerial[i] % 100,  nF++;
	    }
	} else {
	    if (edgematerial[i] / 100)  for (k=0; k<layers; k++) {
		facea[3*nFa+d] = k*nV+v1,  facea[3*nFa+1-d] = (k+1)*nV+v1,  facea[3*nFa+2] = k*nV+v2,  facemata[nFa] = edgematerial[i] / 100,  nFa++;
		facea[3*nFa+d] = (k+1)*nV+v2,  facea[3*nFa+1-d] = k*nV+v2,  facea[3*nFa+2] = (k+1)*nV+v1,  facemata[nFa] = edgematerial[i] / 100,  nFa++;
	    }
	}
	edge[2*i+0] += layers*nV,  edge[2*i+1] += layers*nV;
    }
    nV += layers*nV;
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

static int colorsides(int nV, double *vertex, int nT, int *tetra, int *tetramaterial) {
    int i, r = 0;
    double x, y, z;

    for (i=0; i<nT; i++) {
	if (tetramaterial[i] != 2)  continue;
	x = (vertex[3*tetra[4*i+0]+0] + vertex[3*tetra[4*i+1]+0] + vertex[3*tetra[4*i+2]+0] + vertex[3*tetra[4*i+3]+0]) / 4.0;
	y = (vertex[3*tetra[4*i+0]+1] + vertex[3*tetra[4*i+1]+1] + vertex[3*tetra[4*i+2]+1] + vertex[3*tetra[4*i+3]+1]) / 4.0;
	z = (vertex[3*tetra[4*i+0]+2] + vertex[3*tetra[4*i+1]+2] + vertex[3*tetra[4*i+2]+2] + vertex[3*tetra[4*i+3]+2]) / 4.0;
	if (fabs(x) > sideext)  continue;
	if (fabs(x) > sidelim)  tetramaterial[i] = 5,  r++;
    }
    return r;
}


/* Main */
int main(int argc, char* argv[]) {
    int    nnF = 600000, nnV = 600000, nnT = 3300000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0, *facedup = 0, *facematdup = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    i, nFdup, r, j, medge, msurf, m, k, p;
    int    nVVert,    nLine,  nSurface;
    int    *LineD,    *LineP;
    double *LineT;
    int    *SurfL,    *SurfI;
    double *SurfT;
    double *VVert;
    int *expCrv;
    int nE = 0, nnE = 100000, *edge, *edgematerial;

    double zhr0, zhr1, zhl0, zhl1;
    int fs1a, fs1b, fs1c, fs1d, fs2a, fs2b, fs2c, fs2d,
	fs3a, fs3b, fs3c, fs3d, fs4a, fs4b, fs4c, fs4d,
	fa3a, fa3b, fa3c, fa3d, fa4a, fa4b, fa4c, fa4d,
	fa5a, fa5b, fa5c, fa5d, fa6a, fa6b, fa6c, fa6d;
    int arc1a, arc1b, arc1c, arc1d, arc2a, arc2b, arc2c, arc2d,
	arc3a, arc3b, arc4a, arc4b, arc3c, arc3d, arc4c, arc4d,
	arc5a, arc5b, arc6a, arc6b, arc5c, arc5d, arc6c, arc6d;
    int sgm1a, sgm1b, sgm1c, sgm1d, sgm2a, sgm2b, sgm2c, sgm2d,
	sgm3a, sgm3b, sgm3c, sgm3d, sgm4a, sgm4b, sgm4c, sgm4d;
    int face1, face2, face3, face4;

    int ntfix = 0,  nvfix = 0;

    int izero = 0;

    int variant = 1;
    char type='a';

    R = D/2.0;

    if (argc>1)  sscanf(argv[1], "%d%c", &variant, &type);

    if (argc>2)  sscanf(argv[2], "%d", &serie);

    if (argc>3)  sscanf(argv[3], "%lf", &size);

    if ((variant==2) && (type!='a'))  h1 = 0.5;
    if ((variant==3) && (type!='a'))  h2 = 0.8;
    if ((variant==4) && (type!='a'))  h2 = 0.8;
    if ((variant==5) && (type!='a'))  h2 = 0.8;
    if ((variant==6) && (type!='a'))  h1 = 0.5;
    
    if (variant==2)  sidelim = L/2.0-D*(h1+0.05);
    else if (variant==3)  sidelim = L/2.0-D*(h2+0.05);
    else if (variant==4)  sidelim = L/2.0-D*(h2+0.05);
    else if (variant==5)  sidelim = L/2.0-D*(h2+0.05);
    else if (variant==6)  sidelim = L/2.0-D*(h1+0.05);
    if (variant==3)  sideext = L/2.0-D*(h1+0.05);
    else if (variant==4)  sideext = L/2.0-D*(h1+0.05);
    else if (variant==5)  sideext = L/2.0-D*(h1+0.05);

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
#define ADD_VERTEX(X, Y, Z)  (VVert[3*nVVert+0] = X,  VVert[3*nVVert+1] = Y,  VVert[3*nVVert+2] = Z,  ++nVVert)
#define EDGE(V1, V2, N, ...)  add_edge(LineD, LineP, LineT, &nLine, &medge, V1, V2, N, __VA_ARGS__)
#define SURF(P, C1, C2, D, U1, U2, V1, V2, N, ...) add_surf(SurfL, SurfT, SurfI, &nSurface, &msurf, P, C1, C2, D, U1, U2, V1, V2, N, __VA_ARGS__)
    /***/

    ADD_VERTEX(-L/2.0,  R,  R);
    ADD_VERTEX( L/2.0,  R, -R);
    ADD_VERTEX(-L/2.0, -R, -R);
    ADD_VERTEX( L/2.0, -R,  R);
    if (serie==3)  for (i=0; i<4; i++)  ADD_VERTEX(0.0, 0.0, 0.0);
    
    m = (serie>1) ? 2 : 1,  k = (serie==1) ? 1 : 1;
    p = ((variant <= 2) || (variant == 6));

    medge = 0;

    /* EDGE(v1, v2, nSurf,  [iSurf, iLine, t_0, t1,] ...); */

    arc1a = EDGE(1, 2, 1, m, 1, M_PI, 0.0);
    arc1b = EDGE(2, 1, 1, m, 1, 2.0*M_PI, M_PI);
    arc2a = EDGE(3, 4, 1, m, 2, 0.0,  M_PI);
    arc2b = EDGE(4, 3, 1, m, 2, M_PI, 2.0*M_PI);
    if (p)  expCrv[arc1a-1] = 111,  expCrv[arc1b-1] = 111,  expCrv[arc2a-1] = 112,  expCrv[arc2b-1] = 112;
    else    expCrv[arc1a-1] = 1,    expCrv[arc1b-1] = 1,    expCrv[arc2a-1] = 1,    expCrv[arc2b-1] = 1;
    if (serie==3) {
	arc1c = EDGE(nVVert-3, nVVert-2, 1, 3, 1, 0.0,  M_PI);
	arc1d = EDGE(nVVert-2, nVVert-3, 1, 3, 1, M_PI, 2.0*M_PI);
	arc2c = EDGE(nVVert-1, nVVert-0, 1, 3, 2, 0.0,  M_PI);
	arc2d = EDGE(nVVert-0, nVVert-1, 1, 3, 2, M_PI, 2.0*M_PI);
    }
    if ((variant==2) || (variant==3) || (variant==5)) {
	fa3a = ADD_VERTEX(0,0,0),  fa3b = ADD_VERTEX(0,0,0),  fa3c = ADD_VERTEX(0,0,0),  fa3d = ADD_VERTEX(0,0,0);
	fa4a = ADD_VERTEX(0,0,0),  fa4b = ADD_VERTEX(0,0,0),  fa4c = ADD_VERTEX(0,0,0),  fa4d = ADD_VERTEX(0,0,0);
	arc3a = EDGE(fa3a, fa3b, 1, m, 3,  0.0, M_PI);
	arc3b = EDGE(fa3b, fa3a, 1, m, 3,  M_PI, 2.0*M_PI);
	arc3c = EDGE(fa3c, fa3d, 1, m, 5,  M_PI, 0.0);
	arc3d = EDGE(fa3d, fa3c, 1, m, 5,  2.0*M_PI, M_PI);
	arc4a = EDGE(fa4a, fa4b, 1, m, 4,  M_PI, 0.0);
	arc4b = EDGE(fa4b, fa4a, 1, m, 4,  2.0*M_PI, M_PI);
	arc4c = EDGE(fa4c, fa4d, 1, m, 6,  0.0, M_PI);
	arc4d = EDGE(fa4d, fa4c, 1, m, 6,  M_PI, 2.0*M_PI);
	expCrv[arc3a-1] = 100,  expCrv[arc3b-1] = 100,  expCrv[arc4a-1] = 100,  expCrv[arc4b-1] = 100;
	expCrv[arc3c-1] = 100,  expCrv[arc3d-1] = 100,  expCrv[arc4c-1] = 100,  expCrv[arc4d-1] = 100;
    }
    if (variant == 3) {
	fa5a = ADD_VERTEX(0,0,0),  fa5b = ADD_VERTEX(0,0,0),  fa5c = ADD_VERTEX(0,0,0),  fa5d = ADD_VERTEX(0,0,0);
	fa6a = ADD_VERTEX(0,0,0),  fa6b = ADD_VERTEX(0,0,0),  fa6c = ADD_VERTEX(0,0,0),  fa6d = ADD_VERTEX(0,0,0);
	arc5a = EDGE(fa5a, fa5b, 1, m, 7,  0.0, M_PI);
	arc5b = EDGE(fa5b, fa5a, 1, m, 7,  M_PI, 2.0*M_PI);
	arc5c = EDGE(fa5c, fa5d, 1, m, 9,  M_PI, 0.0);
	arc5d = EDGE(fa5d, fa5c, 1, m, 9,  2.0*M_PI, M_PI);
	arc6a = EDGE(fa6a, fa6b, 1, m, 8,  M_PI, 0.0);
	arc6b = EDGE(fa6b, fa6a, 1, m, 8,  2.0*M_PI, M_PI);
	arc6c = EDGE(fa6c, fa6d, 1, m, 10, 0.0, M_PI);
	arc6d = EDGE(fa6d, fa6c, 1, m, 10, M_PI, 2.0*M_PI);
	expCrv[arc5a-1] = 100,  expCrv[arc5b-1] = 100,  expCrv[arc6a-1] = 100,  expCrv[arc6b-1] = 100;
	expCrv[arc5c-1] = 100,  expCrv[arc5d-1] = 100,  expCrv[arc6c-1] = 100,  expCrv[arc6d-1] = 100;
    }
    if ((variant == 4) || (variant == 6)) {
	fs1a = ADD_VERTEX(0,0,0),  fs1b = ADD_VERTEX(0,0,0),  fs1c = ADD_VERTEX(0,0,0),  fs1d = ADD_VERTEX(0,0,0);
	fs2a = ADD_VERTEX(0,0,0),  fs2b = ADD_VERTEX(0,0,0),  fs2c = ADD_VERTEX(0,0,0),  fs2d = ADD_VERTEX(0,0,0);
	sgm1a  = EDGE(fs1a, fs1b, 1, m, 3,  -0.1, 0.1);
	sgm1b  = EDGE(fs1b, fs1c, 1, m, 12, -L/2.0+D*h1, -L/2.0+D*(h1+0.1));
	sgm1c  = EDGE(fs1c, fs1d, 1, m, 5,  0.1, -0.1);
	sgm1d  = EDGE(fs1d, fs1a, 1, m, 11, -L/2.0+D*(h1+0.1), -L/2.0+D*h1);
	sgm2a  = EDGE(fs2b, fs2a, 1, m, 4,  0.1, -0.1);
	sgm2b  = EDGE(fs2c, fs2b, 1, m, 12,  L/2.0-D*(h1+0.1), L/2.0-D*h1);
	sgm2c  = EDGE(fs2d, fs2c, 1, m, 6,  -0.1, 0.1);
	sgm2d  = EDGE(fs2a, fs2d, 1, m, 11,  L/2.0-D*h1, L/2.0-D*(h1+0.1));
	expCrv[sgm1a-1] = 100,  expCrv[sgm1b-1] = 100,  expCrv[sgm1c-1] = 100,  expCrv[sgm1d-1] = 100;
	expCrv[sgm2a-1] = 100,  expCrv[sgm2b-1] = 100,  expCrv[sgm2c-1] = 100,  expCrv[sgm2d-1] = 100;
    }
    if ((variant == 4) || (variant == 5)) {
	fs3a = ADD_VERTEX(0,0,0),  fs3b = ADD_VERTEX(0,0,0),  fs3c = ADD_VERTEX(0,0,0),  fs3d = ADD_VERTEX(0,0,0);
	fs4a = ADD_VERTEX(0,0,0),  fs4b = ADD_VERTEX(0,0,0),  fs4c = ADD_VERTEX(0,0,0),  fs4d = ADD_VERTEX(0,0,0);
	sgm3a  = EDGE(fs3a, fs3b, 1, m, 7,  -0.1, 0.1);
	sgm3b  = EDGE(fs3b, fs3c, 1, m, 12, -L/2.0+D*h2, -L/2.0+D*(h2+0.1));
	sgm3c  = EDGE(fs3c, fs3d, 1, m, 9,  0.1, -0.1);
	sgm3d  = EDGE(fs3d, fs3a, 1, m, 11, -L/2.0+D*(h2+0.1), -L/2.0+D*h2);
	sgm4a  = EDGE(fs4b, fs4a, 1, m, 8,  0.1, -0.1);
	sgm4b  = EDGE(fs4c, fs4b, 1, m, 12,  L/2.0-D*(h2+0.1), L/2.0-D*h2);
	sgm4c = EDGE(fs4d, fs4c, 1, m, 10, -0.1, 0.1);
	sgm4d = EDGE(fs4a, fs4d, 1, m, 11,  L/2.0-D*h2, L/2.0-D*(h2+0.1));
	expCrv[sgm3a-1] = 100,  expCrv[sgm3b-1] = 100,  expCrv[sgm3c-1] = 100,  expCrv[sgm3d-1] = 100;
	expCrv[sgm4a-1] = 100,  expCrv[sgm4b-1] = 100,  expCrv[sgm4c-1] = 100,  expCrv[sgm4d-1] = 100;
    }

    msurf = 0;

    /* SURF(iSurf, color1, color2, direction, u0, u1, v0, v1, n,  [edge, edge_dir,] ...); */

    if (serie<3) {
	face1 = SURF(0, 2, (p) ? 11 : 101, 1, 0.0, 0.0, 0.0, 0.0, 2, arc1a, 0, arc1b, 0);
	face2 = SURF(0, 2, (p) ? 12 : 101, 0, 0.0, 0.0, 0.0, 0.0, 2, arc2a, 0, arc2b, 0);
    } else {
	face1 = SURF(0, 2, (p) ? 11 : 101, 1, 0.0, 0.0, 0.0, 0.0, 4, arc1a, 0, arc1b, 0, arc1c, 0, arc1d, 0);
	face2 = SURF(0, 2, (p) ? 12 : 101, 0, 0.0, 0.0, 0.0, 0.0, 4, arc2a, 0, arc2b, 0, arc2c, 1, arc2d, 1);
	SURF(0, 3, (p) ? 11 : 101, 1, 0.0, 0.0, 0.0, 0.0, 2, arc1c, 1, arc1d, 1);
	SURF(0, 3, (p) ? 12 : 101, 0, 0.0, 0.0, 0.0, 0.0, 2, arc2c, 0, arc2d, 0);
    }
    if (serie==3)  SURF(3, 3, 2,  1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc1c, 0, arc1d, 0, arc2c, 1, arc2d, 1);
    if (variant==1) {
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc1a, 1, arc1b, 1, arc2a, 1, arc2b, 1);
    } else if (variant==2) {
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc1a, 1, arc1b, 1, arc3a, 1, arc3b, 1);
	SURF(m, 2, 113, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc3a, 0, arc3b, 0, arc3c, 0, arc3d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc3c, 1, arc3d, 1, arc4c, 1, arc4d, 1);
	SURF(m, 2, 114, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4c, 0, arc4d, 0, arc4a, 0, arc4b, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4a, 1, arc4b, 1, arc2a, 1, arc2b, 1);
    } else if (variant==3) {
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc1a, 1, arc1b, 1, arc3a, 1, arc3b, 1);
	SURF(m, 2, 111, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc3a, 0, arc3b, 0, arc3c, 0, arc3d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc3c, 1, arc3d, 1, arc5a, 1, arc5b, 1);
	SURF(m, 2, 113, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc5a, 0, arc5b, 0, arc5c, 0, arc5d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc5c, 1, arc5d, 1, arc6c, 1, arc6d, 1);
	SURF(m, 2, 114, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc6c, 0, arc6d, 0, arc6a, 0, arc6b, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc6a, 1, arc6b, 1, arc4c, 1, arc4d, 1);
	SURF(m, 2, 112, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4c, 0, arc4d, 0, arc4a, 0, arc4b, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4a, 1, arc4b, 1, arc2a, 1, arc2b, 1);
    } else if (variant==4) {
	SURF(m, 2, 111, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm1a, 0, sgm1b, 0, sgm1c, 0, sgm1d, 0);
	SURF(m, 2, 112, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm2a, 0, sgm2b, 0, sgm2c, 0, sgm2d, 0);
	SURF(m, 2, 113, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm3a, 0, sgm3b, 0, sgm3c, 0, sgm3d, 0);
	SURF(m, 2, 114, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm4a, 0, sgm4b, 0, sgm4c, 0, sgm4d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 20, arc1a, 1, arc1b, 1, arc2a, 1, arc2b, 1,
		sgm1a, 1, sgm1b, 1, sgm1c, 1, sgm1d, 1,
		sgm2a, 1, sgm2b, 1, sgm2c, 1, sgm2d, 1,
		sgm3a, 1, sgm3b, 1, sgm3c, 1, sgm3d, 1,
		sgm4a, 1, sgm4b, 1, sgm4c, 1, sgm4d, 1);
    } else if (variant==5) {
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc1a, 1, arc1b, 1, arc3a, 1, arc3b, 1);
	SURF(m, 2, 111, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc3a, 0, arc3b, 0, arc3c, 0, arc3d, 0);
	SURF(m, 2, 112, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4c, 0, arc4d, 0, arc4a, 0, arc4b, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, arc4a, 1, arc4b, 1, arc2a, 1, arc2b, 1);
	SURF(m, 2, 113, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm3a, 0, sgm3b, 0, sgm3c, 0, sgm3d, 0);
	SURF(m, 2, 114, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm4a, 0, sgm4b, 0, sgm4c, 0, sgm4d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 12, arc3c, 1, arc3d, 1, arc4c, 1, arc4d, 1,
		sgm3a, 1, sgm3b, 1, sgm3c, 1, sgm3d, 1,
		sgm4a, 1, sgm4b, 1, sgm4c, 1, sgm4d, 1);
    } else if (variant==6) {
	SURF(m, 2, 113, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm1a, 0, sgm1b, 0, sgm1c, 0, sgm1d, 0);
	SURF(m, 2, 114, 1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 4, sgm2a, 0, sgm2b, 0, sgm2c, 0, sgm2d, 0);
	SURF(m, 2, k,   1, -L/2.0, L/2.0, 0.0, 2.0*M_PI, 12, arc1a, 1, arc1b, 1, arc2a, 1, arc2b, 1,
		sgm1a, 1, sgm1b, 1, sgm1c, 1, sgm1d, 1,
		sgm2a, 1, sgm2b, 1, sgm2c, 1, sgm2d, 1);
    }

    tria_dump_front = 0;
    tria_debug_front = 0;
    region_dump_face = 0;

    printf("\n * Generating surface mesh\n");
    i = ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
	    NULL, surface_param, line_param, periodic, fsize,
	    expCrv,
	    &nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,
	    &nnV, &nnF, &nnE
	    );
    free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);


    printf("INFO: nV = %d, nE = %d, nF = %d, nT = %d\n", nV, nE, nF, nT);

    /* Generate skin mesh */
    printf("\n * Generating skin mesh\n");
    nFdup = 0;
    addlayer(0, (serie==1)?0.0:skinsize, &nV, vertex-3, nE, edge, edgematerial, &nF, face, facematerial, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, nnV, nnF, nnT);
    fix_vertices(&nV, vertex-3, nF, face, nT, tetra, nFdup, facedup, nE, edge);
    addlayer(1, electrodesize, &nV, vertex-3, nE, edge, edgematerial, &nF, face, facematerial, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, nnV, nnF, nnT);
    fix_vertices(&nV, vertex-3, nF, face, nT, tetra, nFdup, facedup, nE, edge);
    printf("INFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

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

	    keepskin(&nFdup, facedup, facematdup);

	    /* Improve mesh quality */
	    printf("\n * Smoothing volume mesh\n");
	    fixshape(&nV, vertex, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, 0, ntfix, nFdup, nnV, nnT, nnF);

	    keepskin(&nFdup, facedup, facematdup);
	    
	    if (serie == 1)  colorsides(nV, vertex-3, nT, tetra, tetramaterial);

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
    return 0;
}

