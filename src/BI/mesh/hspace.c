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

static double electrodesize = 0.5;

static double R1 = 6.0,  R2 = 11.0,  R3 = 14.0;

static double Db = 400.0;

/* this function controls the desired size of the mesh in space */
/* x, y, z  -- coords of point in space                         */
static double fsize(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z) - 36.0;
    double m = size,  c = 0.6;
    if (r < 0)  r = 0.0;
    return sqrt(m*m + c*c*r*r);
}

/* Surface parameterization function
 */
static int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
    double x = 0.0, y = 0.0, z = 0.0;
    if (i == 1) {
	x = u,  y = v,  z = 0.0;
    } else {
	printf("surface_param: wrong id\n");
    }
    *px = x,  *py = y,  *pz = z;
    return 1;
}

/* Line parameterization function */
static int line_param(int i, double t, double *pu, double *pv) {
    double u = 0.0,  v = 0.0;
    if (i==1)       u = R1*cos(t),  v = R1*sin(t);
    else if (i==2)  u = R2*cos(t),  v = R2*sin(t);
    else if (i==3)  u = R3*cos(t),  v = R3*sin(t);
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


static int addlayer(double h,
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
	if (facematerial[i] > 10) {
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
	if (facematerial[i] > 10) {
	    j = layers;
	    for (k=0; k<j; k++)  addprism(&nT, tetra, tetramaterial, k*nV+face[3*i+0], k*nV+face[3*i+1], k*nV+face[3*i+2], nV, 7);
	    facea[3*nFa+0] = j*nV+face[3*i+2],  facea[3*nFa+1] = j*nV+face[3*i+1],  facea[3*nFa+2] = j*nV+face[3*i+0],  facemata[nFa] = facematerial[i],  nFa++;
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
	if (edgematerial[i])  for (k=0; k<layers; k++) {
	    facea[3*nFa+d] = k*nV+v1,  facea[3*nFa+1-d] = (k+1)*nV+v1,  facea[3*nFa+2] = k*nV+v2,  facemata[nFa] = edgematerial[i],  nFa++;
	    facea[3*nFa+d] = (k+1)*nV+v2,  facea[3*nFa+1-d] = k*nV+v2,  facea[3*nFa+2] = (k+1)*nV+v1,  facemata[nFa] = edgematerial[i],  nFa++;
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

/* Main */
int main(int argc, char* argv[]) {
    int    nnF = 200000, nnV = 200000, nnT = 1100000;
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

    int va, vb, vc, vd, ve, vf, vg, vh;
    int vc1a, vc1b, vc2a, vc2b, vc3a, vc3b;
    int arc1a, arc1b, arc2a, arc2b, arc3a, arc3b;
    int eab, ebc, ecd, eda, eef, efg, egh, ehe, eae, ebf, ecg, edh;

    int ntfix = 0,  nvfix = 0;

    int izero = 0;

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

    vc1a = ADD_VERTEX( R1, 0, 0);
    vc1b = ADD_VERTEX(-R1, 0, 0);
    vc2a = ADD_VERTEX( R2, 0, 0);
    vc2b = ADD_VERTEX(-R2, 0, 0);
    vc3a = ADD_VERTEX( R3, 0, 0);
    vc3b = ADD_VERTEX(-R3, 0, 0);
    va = ADD_VERTEX(-Db, -Db, 0);
    vb = ADD_VERTEX(-Db,  Db, 0);
    vc = ADD_VERTEX( Db,  Db, 0);
    vd = ADD_VERTEX( Db, -Db, 0);
    ve = ADD_VERTEX(-Db, -Db, -Db);
    vf = ADD_VERTEX(-Db,  Db, -Db);
    vg = ADD_VERTEX( Db,  Db, -Db);
    vh = ADD_VERTEX( Db, -Db, -Db);

    medge = 0;

    /* EDGE(v1, v2, nSurf,  [iSurf, iLine, t_0, t1,] ...); */

    arc1a = EDGE(vc1a, vc1b, 1, 1, 1, 0.0, M_PI);
    arc1b = EDGE(vc1b, vc1a, 1, 1, 1, M_PI, 2.0*M_PI);
    arc2a = EDGE(vc2a, vc2b, 1, 1, 2, M_PI, 0.0);
    arc2b = EDGE(vc2b, vc2a, 1, 1, 2, 2.0*M_PI, M_PI);
    arc3a = EDGE(vc3a, vc3b, 1, 1, 3, 0.0, M_PI);
    arc3b = EDGE(vc3b, vc3a, 1, 1, 3, M_PI, 2.0*M_PI);
    expCrv[arc1a-1] = 1,  expCrv[arc1b-1] = 1;
    expCrv[arc2a-1] = 1,  expCrv[arc2b-1] = 1;
    expCrv[arc3a-1] = 1,  expCrv[arc3b-1] = 1;
    
    eab = EDGE(va, vb, 1, 0, 0, 0.0, 0.0);
    ebc = EDGE(vb, vc, 1, 0, 0, 0.0, 0.0);
    ecd = EDGE(vc, vd, 1, 0, 0, 0.0, 0.0);
    eda = EDGE(vd, va, 1, 0, 0, 0.0, 0.0);

    eef = EDGE(ve, vf, 1, 0, 0, 0.0, 0.0);
    efg = EDGE(vf, vg, 1, 0, 0, 0.0, 0.0);
    egh = EDGE(vg, vh, 1, 0, 0, 0.0, 0.0);
    ehe = EDGE(vh, ve, 1, 0, 0, 0.0, 0.0);

    eae = EDGE(va, ve, 1, 0, 0, 0.0, 0.0);
    ebf = EDGE(vb, vf, 1, 0, 0, 0.0, 0.0);
    ecg = EDGE(vc, vg, 1, 0, 0, 0.0, 0.0);
    edh = EDGE(vd, vh, 1, 0, 0, 0.0, 0.0);

    msurf = 0;

    /* SURF(iSurf, color1, color2, direction, u0, u1, v0, v1, n,  [edge, edge_dir,] ...); */

    SURF(1, 2, 11, 0, -Db, Db, -Db, Db, 2, arc1a, 0, arc1b, 0);
    SURF(1, 1,  0, 0, -Db, Db, -Db, Db, 4, arc2a, 1, arc2b, 1, arc1a, 1, arc1b, 1);
    SURF(1, 2, 12, 0, -Db, Db, -Db, Db, 4, arc3a, 0, arc3b, 0, arc2a, 0, arc2b, 0);
    
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 6, arc3a, 1, arc3b, 1, eab, 1, ebc, 1, ecd, 1, eda, 1);
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 4, eef, 0, efg, 0, egh, 0, ehe, 0);
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 4, eab, 0, ebf, 0, eef, 1, eae, 1);
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 4, ebc, 0, ecg, 0, efg, 1, ebf, 1);
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 4, ecd, 0, edh, 0, egh, 1, ecg, 1);
    SURF(0, 1, 0, 0, 0.0, 0.0, 0.0, 0.0, 4, eda, 0, eae, 0, ehe, 1, edh, 1);

    tria_dump_front = 0;
    tria_debug_front = 0;
    region_dump_face = 0;

    printf("\n * Generating surface mesh\n");
    i = ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
	    NULL, surface_param, line_param, NULL/*periodic*/, fsize,
	    expCrv,
	    &nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,
	    &nnV, &nnF, &nnE
	    );
    free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);


    printf("INFO: nV = %d, nE = %d, nF = %d, nT = %d\n", nV, nE, nF, nT);

    /* Generate skin mesh */
    printf("\n * Generating skin mesh\n");
    nFdup = 0;
    addlayer(electrodesize, &nV, vertex-3, nE, edge, edgematerial, &nF, face, facematerial, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, nnV, nnF, nnT);
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

	    for (i=0; i<nT; i++)  tetramaterial[i] = (tetramaterial[i]==1) ? 2 : tetramaterial[i];

	    keepskin(&nFdup, facedup, facematdup);

	    /* Improve mesh quality */
	    printf("\n * Smoothing volume mesh\n");
	    fixshape(&nV, vertex, &nT, tetra, tetramaterial, &nFdup, facedup, facematdup, 0, ntfix, nFdup, nnV, nnT, nnF);

	    keepskin(&nFdup, facedup, facematdup);

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

