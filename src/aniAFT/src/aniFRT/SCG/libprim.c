#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libfrtprm.h"

#define EPS 1e-14

// this function controls the desired size of the mesh in space
static double meshsize = 0.1;
double fsize(double x, double y, double z)
{
	(void)x, (void)y, (void)z;
	return meshsize;
}


static double CX = 0.0, CY = 0.0, CZ = 0.0; // center
static double R = 1.0; // for sphere
static double CR = 1.0, CH = 1.0; // for cylinder


// We should define parametrization functions for the sphere.
// We will use two functions, first one for the top half of sphere,
// and second for bottom half.
// Parameters (u,v) will correspond to (x,y). z-coord will be calculated from (x,y).
// surface_param() gets integer i -- number of parametrization function to use
// and params (u,v). This function should return (x,y,z) = F_i(u,v).
static int surface_param_sphere(int i, double u, double v, double *px, double *py, double *pz) {
	double x = CX + R*u;
	double y = CY + R*v;
	double z = CZ;
	double r = u*u + v*v;
	double h;
	r = 1.0 - r;
	if (r>EPS) h = sqrt(r);
	else if (r<-EPS) h = -sqrt(-r);
	else h = 0.0;
	if (i==1) z += R*h;
	else if (i==2) z -= R*h;
	*px = x;
	*py = y;
	*pz = z;
	return 1;
}

// The boundary of half-sphere is a circle, we could define this circle as two arcs.
// Curves should be parametrized for each parametrized surface they belong to.
// We use V(u) parametrization, so u -> (u,V(u)) -> (x,y,z).
// v_u_param() returns v = V_i(u)
// In this example we calculate v as sqrt(1-u*u) for first arc
// and v = -sqrt(1-u*u) for second one.
static void v_u_param_sphere(int i, double u, double *pv) {
	double t = 1.0-u*u;
	double v = 0.0;
	if (t>EPS) t = sqrt(t);
	else if (t<-EPS) t = -sqrt(-t);
	else t = 0.0;
	if (i==1) v = t;
	else if (i==2) v = -t;
	*pv = v;
}

#define ADD_VERTEX(X, Y, Z); {\
	VVert[3*nVVert+0] = X;\
	VVert[3*nVVert+1] = Y;\
	VVert[3*nVVert+2] = Z;\
	nVVert++;\
}
int scg_make_sphere (int *pnV, int *pnF, double *Vertex, int *Index, double x, double MS)
{
	int nnV = 1000000, nnE = 10, nnF = 1000000, nnS = 10;
	int *face = 0, *facematerial = 0;
	double *vertex = 0;
	int i, m;

	int    nVVert,    nLine,  nSurface;
	int    *LineD,    *LineP;
	double *LineT;
	int    *SurfL,    *SurfI;
	double *SurfT;
	double *VVert;

	R = x;

	*pnV = 0, *pnF = 0;
	meshsize = MS;

	// allocating memory
	VVert = (double*)malloc(sizeof(double) * 3*nnV);
	nVVert = 0;
	LineD = (int*)   malloc(sizeof(int)    * 3*nnE);
	LineP = (int*)   malloc(sizeof(int)    * 2*2*nnE);
	LineT = (double*)malloc(sizeof(double) * 2*2*nnE);
	nLine = 0;
	SurfL = (int*)   malloc(sizeof(int)    * 5*nnS);
	SurfI = (int*)   malloc(sizeof(int)    * 2*2*nnE);
	SurfT = (double*)malloc(sizeof(double) * 4*nnS);
	nSurface = 0;

	// allocate memory for mesh ctructures
	vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
	face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facematerial  = (int*)   malloc(sizeof(int)        * nnF);

	ADD_VERTEX(CX - R, CY - R, CZ - R); // In fact coords of these points will be recalculated
	ADD_VERTEX(CX + R, CY + R, CZ + R); // This coords defines bounding box

	m = 0;

	// First arc on equator from point 1 to point 2
	LineD[3*nLine+0] = 1; LineD[3*nLine+1] = 2;
	// This arc belongs to two surfaces (halfs of the sphere)
	LineD[3*nLine+2] = 2;
	// On both halfs of the sphere this arc is parametrized by the same (u,v)
	// We will use V_1 for both surfaces, u in [-1.0, 1.0]
	LineP[2*m+0] = 1; LineP[2*m+1] = 1; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	LineP[2*m+0] = 2; LineP[2*m+1] = 1; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	nLine++; // 1

	// Second arc on equator from point 2 to point 1
	LineD[3*nLine+0] = 2; LineD[3*nLine+1] = 1;
	// This arc also belongs to two surfaces (halfs of the sphere)
	LineD[3*nLine+2] = 2;
	// On both halfs of the sphere this arc is parametrized by the same (u,v)
	// We will use V_2 for both surfaces, u in [1.0, -1.0]
	LineP[2*m+0] = 1; LineP[2*m+1] = 2; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	LineP[2*m+0] = 2; LineP[2*m+1] = 2; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	nLine++; // 2

	m = 0;

	// top half of the sphere, parametrized by F_1
	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 1;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 0;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 1; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 2; SurfI[2*m+1] = 1; m++;
	nSurface++; // 1

	// bottom half of the sphere, parametrized by F_2
	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 2;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 1;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 1; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 2; SurfI[2*m+1] = 0; m++;
	nSurface++; // 2
	// Note the difference of orientations of last two surfaces

	// Setting up our boundary parametrization functions
	// First one is for (x,y,z)=F_i(u,v) and second one is for v=V_j(u)
	// set_param_functions(surface_param_sphere, v_u_param_sphere);
	// Now we are ready to create the surface mesh
	// Setting up our own size function.
	// set_surface_size_function(fsize);
	// Making surface mesh
	i = FortranInterface_ani3d_surface_boundary(
		&nVVert, VVert,
		&nLine, LineD, LineP, LineT,
		&nSurface, SurfL, SurfI, SurfT,
		v_u_param_sphere,
		surface_param_sphere,
		NULL, NULL,
		fsize,
		pnV, vertex,
		pnF, face, facematerial,
		&nnV, &nnF
	);
	free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);

	printf("\nINFO: nV = %d, nF = %d\n", *pnV, *pnF);

	memset (Vertex, 0 , sizeof (double) * 3);
	for (i = 0; i < 3 * (*pnV); i++)
	{
		Vertex[i] = vertex[i];
	}

	for (i = 0; i < (*pnF); i++)
	{
		Index[3 * i] = face[3 * i];
		Index[3 * i + 1] = face[3 * i + 1];
		Index[3 * i + 2] = face[3 * i + 2];
	}

	printf ("nV = %d, nF = %d\n", *pnV, *pnF);
	free(vertex), free(face), free(facematerial);

	return 0;
}

// this is a macros to simplify the process of adding new flat surface with one loop
// it allocates memory for edges, fills c1
// and adds one loop (N, ...)
#define ADDFACE(C1, N, ...) \
	nE[nL] = 0; \
edge[nL] = (int*)malloc (sizeof (int) * 2 * nnE); \
c1[nL] = C1; \
ani3d_polyhedron_add_faceloop(&nE[nL], edge[nL], N, __VA_ARGS__); \
nL++;
int scg_make_paral (int *pnV, int *pnF, double *Vertex, int *Index, double x, double y, double z, double MS)
{
	enum {nnV = 1000000, nnL = 10, nnE = 10, nnF = 1000000};
	int nE[nnL];
	int *edge[nnL], c1[nnL];
	int *face = 0, *facematerial = 0;
	double *vertex = 0;
	int nL = 0;
	double xmin = -(double)x, ymin = -(double)y, zmin = -(double)z, xmax = (double)x, ymax = (double)y, zmax = (double)z;
	int i;

	*pnV = 0, *pnF = 0;
	meshsize = MS;

	// allocating memory
	vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
	face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facematerial  = (int*)   malloc(sizeof(int)        * nnF);

	// adding new 8 vertices of the box
	vertex[3*(*pnV+0)+0]=xmin, vertex[3*(*pnV+0)+1]=ymin, vertex[3*(*pnV+0)+2]=zmin;
	vertex[3*(*pnV+1)+0]=xmax, vertex[3*(*pnV+1)+1]=ymin, vertex[3*(*pnV+1)+2]=zmin;
	vertex[3*(*pnV+2)+0]=xmax, vertex[3*(*pnV+2)+1]=ymax, vertex[3*(*pnV+2)+2]=zmin;
	vertex[3*(*pnV+3)+0]=xmin, vertex[3*(*pnV+3)+1]=ymax, vertex[3*(*pnV+3)+2]=zmin;
	vertex[3*(*pnV+4)+0]=xmin, vertex[3*(*pnV+4)+1]=ymin, vertex[3*(*pnV+4)+2]=zmax;
	vertex[3*(*pnV+5)+0]=xmax, vertex[3*(*pnV+5)+1]=ymin, vertex[3*(*pnV+5)+2]=zmax;
	vertex[3*(*pnV+6)+0]=xmax, vertex[3*(*pnV+6)+1]=ymax, vertex[3*(*pnV+6)+2]=zmax;
	vertex[3*(*pnV+7)+0]=xmin, vertex[3*(*pnV+7)+1]=ymax, vertex[3*(*pnV+7)+2]=zmax;

	// adding 6 faces of the box with color 1
	ADDFACE (1, 4,  *pnV+4, *pnV+3, *pnV+2, *pnV+1);
	ADDFACE (1, 4,  *pnV+5, *pnV+6, *pnV+7, *pnV+8);
	ADDFACE (1, 4,  *pnV+1, *pnV+2, *pnV+6, *pnV+5);
	ADDFACE (1, 4,  *pnV+2, *pnV+3, *pnV+7, *pnV+6);
	ADDFACE (1, 4,  *pnV+3, *pnV+4, *pnV+8, *pnV+7);
	ADDFACE (1, 4,  *pnV+4, *pnV+1, *pnV+5, *pnV+8);
	*pnV += 8;

	// Now we are ready to create the surface mesh
	// Setting up our own size function.
	// set_surface_size_function(fsize);
	// Making surface mesh
	ani3d_polyhedron_make_front(
		pnV, vertex,
		nL, nE, edge, c1, 0,
		fsize,
		pnF, face, facematerial,
		nnV, nnF
	);

	memset (Vertex, 0 , sizeof (double) * 3);
	for (i = 0; i < 3 * (*pnV); i++)
	{
		Vertex[i] = vertex[i];
	}

	for (i = 0; i < (*pnF); i++)
	{
		Index[3 * i] = face[3 * i];
		Index[3 * i + 1] = face[3 * i + 1];
		Index[3 * i + 2] = face[3 * i + 2];
	}

	printf ("nV = %d, nF = %d\n", *pnV, *pnF);
	free(vertex), free(face), free(facematerial);

	return 0;
}

static int surface_param_cylinder(int i, double u, double v, double *px, double *py, double *pz) {
	double x = CX;
	double y = CY;
	double z = CZ;

	if (i==1) {
		x += CR*u, y += CR*v;
	} else if (i==2) {
		x += CR*u, y += CR*v, z += CH;
	} else if ((i==3) || (i==4)) {
		x += CR*cos(u+v);
		y += CR*sin(u+v);
		z += CH*v;
	}

	*px = x, *py = y, *pz = z;
	return 1;
}

static void v_u_param_cylinder(int i, double u, double *pv) {
	double v = 0.0;
	double t = 0.0;
	if ((i==1) || (i==2)) {
		t = 1.0 - u*u;
		if (t>EPS) t = sqrt(t);
		else if (t<-EPS) t = -sqrt(-t);
		else t = 0.0;
		if (i==1) v = t;
		else if (i==2) v = -t;
	} else if (i==3) {
		v = 0.0;
	} else if (i==4) {
		v = 1.0;
	} else if (i==5) {
		v = -u;
	} else if (i==6) {
		v = -u - M_PI;
	} else if (i==7) {
		v = -u + M_PI;
	}
	*pv = v;
}

// Under construction;
int scg_make_cylinder (int *pnV, int *pnF, double *Vertex, int *Index, double r, double h, double MS)
{
	int nnV = 1000000, nnE = 10, nnF = 1000000, nnS = 10;
	int *face = 0, *facematerial = 0;
	double *vertex = 0;
	int i, m;

	int    nVVert,    nLine,  nSurface;
	int    *LineD,    *LineP;
	double *LineT;
	int    *SurfL,    *SurfI;
	double *SurfT;
	double *VVert;

	CR = r, CH = h;

	*pnV = 0, *pnF = 0;
	meshsize = MS;

	// allocating memory
	VVert = (double*)malloc(sizeof(double) * 3*nnV);
	nVVert = 0;
	LineD = (int*)   malloc(sizeof(int)    * 3*nnE);
	LineP = (int*)   malloc(sizeof(int)    * 2*2*nnE);
	LineT = (double*)malloc(sizeof(double) * 2*2*nnE);
	nLine = 0;
	SurfL = (int*)   malloc(sizeof(int)    * 5*nnS);
	SurfI = (int*)   malloc(sizeof(int)    * 2*2*nnE);
	SurfT = (double*)malloc(sizeof(double) * 4*nnS);
	nSurface = 0;

	// allocate memory for mesh ctructures
	vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
	face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facematerial  = (int*)   malloc(sizeof(int)        * nnF);

	ADD_VERTEX(CX - CR, CY, CZ); // In fact coords of these points will be recalculated
	ADD_VERTEX(CX + CR, CY, CZ);
	ADD_VERTEX(CX - CR, CY - CR, CZ ); // This coords defines bounding box
	ADD_VERTEX(CX + CR, CY + CR, CZ + CH);

	m = 0;

	LineD[3*nLine+0] = 1; LineD[3*nLine+1] = 2;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 1; LineP[2*m+1] = 1; LineT[2*m+0] = -1.0; LineT[2*m+1] = 1.0; m++;
	LineP[2*m+0] = 3; LineP[2*m+1] = 3; LineT[2*m+0] = M_PI; LineT[2*m+1] = 0.0; m++;
	nLine++;

	LineD[3*nLine+0] = 2; LineD[3*nLine+1] = 1;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 1; LineP[2*m+1] = 2; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	LineP[2*m+0] = 4; LineP[2*m+1] = 3; LineT[2*m+0] =  0.0; LineT[2*m+1] = -M_PI; m++;
	nLine++;

	LineD[3*nLine+0] = 4; LineD[3*nLine+1] = 3;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 2; LineP[2*m+1] = 1; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	LineP[2*m+0] = 3; LineP[2*m+1] = 4; LineT[2*m+0] =  -1.0; LineT[2*m+1] = M_PI-1.0; m++;
	nLine++;

	LineD[3*nLine+0] = 3; LineD[3*nLine+1] = 4;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 2; LineP[2*m+1] = 2; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	LineP[2*m+0] = 4; LineP[2*m+1] = 4; LineT[2*m+0] = -M_PI-1.0; LineT[2*m+1] =  -1.0; m++;
	nLine++;

	LineD[3*nLine+0] = 2; LineD[3*nLine+1] = 4;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 3; LineP[2*m+1] = 5; LineT[2*m+0] = 0.0; LineT[2*m+1] = -1.0; m++;
	LineP[2*m+0] = 4; LineP[2*m+1] = 5; LineT[2*m+0] = 0.0; LineT[2*m+1] = -1.0; m++;
	nLine++;

	LineD[3*nLine+0] = 1; LineD[3*nLine+1] = 3;
	LineD[3*nLine+2] = 2;
	LineP[2*m+0] = 3; LineP[2*m+1] = 7; LineT[2*m+0] = M_PI; LineT[2*m+1] = M_PI-1.0; m++;
	LineP[2*m+0] = 4; LineP[2*m+1] = 6; LineT[2*m+0] = -M_PI; LineT[2*m+1] = -M_PI-1.0; m++;
	nLine++;

	m = 0;
	
	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 1;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 1;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 1; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 2; SurfI[2*m+1] = 0; m++;
	nSurface++;

	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 2;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 0;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 3; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 4; SurfI[2*m+1] = 0; m++;
	nSurface++;

	SurfL[5*nSurface+0] = 4;
	SurfL[5*nSurface+1] = 3;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 0;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  M_PI; SurfT[4*nSurface+2] = 0.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 6; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 1; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 5; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 3; SurfI[2*m+1] = 1; m++;
	nSurface++;

	SurfL[5*nSurface+0] = 4;
	SurfL[5*nSurface+1] = 4;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 0;
	SurfT[4*nSurface+0] = -M_PI-1.0; SurfT[4*nSurface+1] =  0.0; SurfT[4*nSurface+2] = 0.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 6; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 2; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 5; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 4; SurfI[2*m+1] = 1; m++;
	nSurface++;


	// Setting up our boundary parametrization functions
	// First one is for (x,y,z)=F_i(u,v) and second one is for v=V_j(u)
	// set_param_functions(surface_param_cylinder, v_u_param_cylinder);
	// Now we are ready to create the surface mesh
	// Setting up our own size function.
	// set_surface_size_function(fsize);
	// Making surface mesh
	i = FortranInterface_ani3d_surface_boundary(
		&nVVert, VVert,
		&nLine, LineD, LineP, LineT,
		&nSurface, SurfL, SurfI, SurfT,
		v_u_param_cylinder,
		surface_param_cylinder,
		NULL, NULL,
		fsize,
		pnV, vertex,
		pnF, face, facematerial,
		&nnV, &nnF
	);
	free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);

	printf("\nINFO: nV = %d, nF = %d\n", *pnV, *pnF);

	memset (Vertex, 0 , sizeof (double) * 3);
	for (i = 0; i < 3 * (*pnV); i++)
	{
		Vertex[i] = vertex[i];
	}

	for (i = 0; i < (*pnF); i++)
	{
		Index[3 * i] = face[3 * i];
		Index[3 * i + 1] = face[3 * i + 1];
		Index[3 * i + 2] = face[3 * i + 2];
	}

	printf ("nV = %d, nF = %d\n", *pnV, *pnF);
	free(vertex), free(face), free(facematerial);

	return 0;
}

int scg_read_front (int *pnV, int *pnF, double *Vertex, int *Index, char *file, int if_color)
{
    FILE *f;
    int i;

    double glf;
    int u1, u2, u3;

    if (!(f = fopen (file, "r")))
	return -1;

    if (!fscanf (f, "%d", pnV))
	return -2;

    if (!fscanf (f, "%d", pnF))
	return -3;

    if (*pnV < 0 || *pnF < 0)
	return -4;

//    Vertex[0] = Vertex[1] = Vertex[2] = 0.0f; // In case that numeration goes from 1, instead of 0;

    for (i = 0; i < 3 * (*pnV); i++)
    {
	if (!fscanf (f, "%lf", &glf))
            return -5;
	
	Vertex[i] = (double)glf;	
    }
    
    if (if_color == 1)
    {
        for (i = 0; i < (*pnF); i++)
        {
	    if (fscanf (f, "%d %d %d %lf", &u1, &u2, &u3, &glf) < 4)
                return -6;	

	    if (u1 == 0 || u2 == 0 || u3 == 0)
                printf ("Warning, there's down from 0 numeration!\n");
	
	    Index[i * 3] = u1;
	    Index[i * 3 + 1] = u2;
	    Index[i * 3 + 2] = u3;
        }
    }

    else
    {
        for (i = 0; i < (*pnF); i++)
        {
	    if (fscanf (f, "%d %d %d", &u1, &u2, &u3) < 3)
                return -6;	

	    if (u1 == 0 || u2 == 0 || u3 == 0)
                printf ("Warning, there's down from 0 numeration!\n");
	
	    Index[i * 3] = u1;
	    Index[i * 3 + 1] = u2;
	    Index[i * 3 + 2] = u3;
        }
    }

    printf ("%s read. \nnV = %d, nF = %d\n", file, *pnV, *pnF);
    
    fclose (f);
    return 0;
}

int scg_write_front (int nV, int nF, double *Vertex, int *Index, char *file)
{
    FILE *f;
    int i, c = 0;

    if (!(f = fopen (file, "w")))
	return -1;

	 fprintf(f, "%d %d\n", nV, nF);

	 for (i = 0; i < nV; i++) 
		 fprintf(f, "%20.15e %20.15e %20.15e\n", Vertex[3 * i + 0], Vertex[3 * i + 1], Vertex[3 * i + 2]);

	 for (i = 0; i < nF; i++) 
		 fprintf(f, "%d %d %d  %d\n", Index[3 * i + 0], Index[3 * i + 1], Index[3 * i + 2], c);
    
    fclose (f);
    printf("%s written.\n", file);

    return 0;
}

int scg_write_front_gmv (int nV, int nF, double *vertex, int *face, char *file)
{
	int i;
	FILE *f;
	
	f = fopen(file, "w");
	fprintf(f, "gmvinput ascii\n\nnodev %5d\n", nV);
	
	for (i = 0; i < nV; i++) 
	{
		fprintf(f, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
	}
	
	fprintf(f, "\ncells %5d\n", 0);

	fprintf(f, "\n\npolygons\n");
	for (i = 0; i < nF; i++) 
	{
		fprintf(f, "%3d 3", 0);
		fprintf(f, " %20.15lf %20.15lf %20.15lf\n", 
			vertex[3*face[3*i+0]+0], vertex[3*face[3*i+1]+0], vertex[3*face[3*i+2]+0]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", 
			vertex[3*face[3*i+0]+1], vertex[3*face[3*i+1]+1], vertex[3*face[3*i+2]+1]);
		fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", 
			vertex[3*face[3*i+0]+2], vertex[3*face[3*i+1]+2], vertex[3*face[3*i+2]+2]);
	}
	fprintf(f, "endpoly\n");
	fprintf(f, "\nendgmv");
	fclose(f);
	return 0;
}

int scg_translate (int nV, double *Vertex, double x, double y, double z)
{
    int i;

    if (x != 0)
        for (i = 0; i < nV; i++)
            Vertex[i * 3] += x;
    
    if (y != 0)
        for (i = 0; i < nV; i++)
            Vertex[i * 3 + 1] += y;
    
    if (z != 0)
        for (i = 0; i < nV; i++)
            Vertex[i * 3 + 2] += z;

    return 0;
}

int scg_rotate (int nV, double *Vertex, double rot_x, double rot_y, double rot_z)
{
    int i;

    // Rotation of rot_x degrees around (1, 0, 0) vector
    // then of rot_y degrees around (0, 1, 0) vector
    // and then of rot_z degrees around (0, 0, 1) vector
    if (rot_x != 0 || rot_y != 0 || rot_z != 0)
    {
        double r11, r12, r13, r21, r22, r23, r31, r32, r33, v1, v2, v3;
	r11 = cos (rot_y) * cos (rot_z);
	r12 = -cos (rot_x) * sin (rot_z) + sin (rot_x) * sin (rot_y) * cos (rot_z);
	r13 = sin (rot_x) * sin (rot_z) + cos (rot_x) * sin (rot_y) * cos (rot_z);
	r21 = cos (rot_y) * sin (rot_z);
	r22 = cos (rot_x) * cos (rot_z) + sin (rot_x) * sin (rot_y) * sin (rot_z);
	r23 = -sin (rot_x) * cos (rot_z) + cos (rot_x) * sin (rot_y) * sin (rot_z);
	r31 = -sin (rot_y);
	r32 = sin (rot_x) * cos (rot_y);
	r33 = cos (rot_x) * cos (rot_y);

	for (i = 0; i < nV; i++)
	{
            v1 = Vertex[i * 3];
            v2 = Vertex[i * 3 + 1];
            v3 = Vertex[i * 3 + 2];
	    
	    Vertex[i * 3] = r11 * v1 + r12 * v2 + r13 * v3;
	    Vertex[i * 3 + 1] = r21 * v1 + r22 * v2 + r23 * v3;
	    Vertex[i * 3 + 2] = r31 * v1 + r32 * v2 + r33 * v3;	
	}
    }
    
    return 0;
}


int scg_scale (int nV, double *Vertex, double s_x, double s_y, double s_z)
{
    int i;

    for (i = 0; i < nV; i++)
    {
	Vertex[3 * i] *= s_x;
	Vertex[3 * i + 1] *= s_y;
	Vertex[3 * i + 2] *= s_z;	
    }

    return 0;
}

int scg_affine (int nV, double *Vertex, double *matrix)
{
    int i;
    double v1, v2, v3;
    
    for (i = 0; i < nV; i++)
    {
	v1 = Vertex[i * 3];
	v2 = Vertex[i * 3 + 1];
	v3 = Vertex[i * 3 + 2];
	    
	Vertex[i * 3] = matrix[0] * v1 + matrix[1] * v2 + matrix[2] * v3;
	Vertex[i * 3 + 1] = matrix[3] * v1 + matrix[4] * v2 + matrix[5] * v3;
	Vertex[i * 3 + 2] = matrix[6] * v1 + matrix[7] * v2 + matrix[8] * v3;	
    }
    
    return 0;
}

