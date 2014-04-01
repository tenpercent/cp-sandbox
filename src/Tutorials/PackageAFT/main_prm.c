/* main_prm.c
 *
 * This is example for libfrtprm and libaft usage
 *
 * main_prm.exe constructs front for model "cube minus ball", where the ball is
 * inside cube and they do not intersect each other.
 * Tetra mesh is generated using mesh_3d_aft_func_
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "libaft.h"
#include "libfrtprm.h"

extern int tria_dump_front, tria_debug_front, region_dump_face; // debug flags from aftlib, could be usefull for debugging

// In this example we will create boundary representation of the following model
// A small ball will be cutted from the center of the box

double CX = 0.5, CY = 0.5, CZ = 0.5; // center of the box
double R = 0.2; // sphere radius
double L = 0.5; // box half-size


// this function controls the desired size of the mesh in space
// x, y, z  -- coords of point in space
// in this function we use size = sqrt( (cr)^2 + m^2 )
// where r is the distance to center of the small box
// c is the coefficient that controls coarseness
// m is the minimal size in the center of the small box
// near the box we will have  size = m
// far from the box we will have  size = cr
double fsize(double x, double y, double z) {
	double xc = CX;    //
	double yc = CY;    // center of the box
	double zc = CZ;    //
	double c  = 0.2;   // coarseness
	double m2 = 0.001; // minimal size in center
	double r2 = (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);
	double c2 = c*c;
	return sqrt(0.04 - r2*c2)/2.0;
	return sqrt(r2*c2 + m2);
};

// To define the surface of the sphere we should split it at least in two parts.
// We will use half-spheres.

// We should define parametrization functions for the sphere.
// We will use two functions, first one for the top half of sphere,
// and second for bottom half.
// Parameters (u,v) will correspond to (x,y). z-coord will be calculated from (x,y).
// surface_param() gets integer i -- number of parametrization function to use
// and params (u,v). This function should return (x,y,z) = F_i(u,v).
int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
	double x = CX + R*u;
	double y = CY + R*v;
	double z = CZ;
	double r = u*u + v*v;
	double h;
	if (r<=1.0) h = sqrt(1.0-r);
	else h = -sqrt(r-1.0);
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
void v_u_param(int i, double u, double *pv) {
	double t = 1.0-u*u;
	double v = 0.0;
	if (t>=0) t = sqrt(t);
	else t = -sqrt(-t);
	if (i==1) v = t;
	else if (i==2) v = -t;
	*pv = v;
}



// in this program we use indexing from one

int main(int argc, char* argv[]) {
	int    nnF = 1000000, nnV = 1000000, nnT = 1000000;
	int     nF = 0,        nV = 0,        nT = 0;
	int    nnS = 1000,    nnE = 1000;
	int    *face   = 0, *facematerial  = 0, *facedup = 0, *facematdup = 0;
	int    *tetra  = 0, *tetramaterial = 0;
	double *vertex = 0;
	int    i, nFdup, r, m;
	int    nVVert,    nLine,  nSurface;
	int    *LineD,    *LineP;
	double *LineT;
	int    *SurfL,    *SurfI;
	double *SurfT;
	double *VVert;

	(void) argc, (void) argv;


	// allocate memory for mesh ctructures
	vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
	face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facematerial  = (int*)   malloc(sizeof(int)        * nnF);
	facematdup    = (int*)   malloc(sizeof(int)        * nnF);
	tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
	tetramaterial = (int*)   malloc(sizeof(int)        * nnT);


	// allocate memory for boundary representation structure
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

//	First we will define coords of the main points.
//	Main points are the end points of each curve.
//	They are also known as V-points or V-vertices
//	We will define 8 vertices of the box
//	and two points on equator of the sphere
//	Here is an example how we can add point with coords(x,y,z):
/* 
	VVert[3*nVVert+0] = x; 
	VVert[3*nVVert+1] = y;
	VVert[3*nVVert+2] = z;
	nVVert++;
*/
// To simplify our program we can use the following macros:
#define ADD_VERTEX(X, Y, Z); {\
	VVert[3*nVVert+0] = X;\
	VVert[3*nVVert+1] = Y;\
	VVert[3*nVVert+2] = Z;\
	nVVert++;\
}
	// 8 vertices of the box
	ADD_VERTEX(CX - L, CY - L, CZ - L); // 1
	ADD_VERTEX(CX + L, CY - L, CZ - L); // 2
	ADD_VERTEX(CX + L, CY + L, CZ - L); // 3
	ADD_VERTEX(CX - L, CY + L, CZ - L); // 4
	ADD_VERTEX(CX - L, CY - L, CZ + L); // 5
	ADD_VERTEX(CX + L, CY - L, CZ + L); // 6
	ADD_VERTEX(CX + L, CY + L, CZ + L); // 7
	ADD_VERTEX(CX - L, CY + L, CZ + L); // 8
	// If the vertex is the end of parametrized curve it coords will be forcly recalculated
	// So for such points we can just say nVVert++;
	// But in this example we will still fill the coords for clarity/
	ADD_VERTEX(CX - R, CY, CZ); // 9   | In fact coords of these points will be recalculated
	ADD_VERTEX(CX + R, CY, CZ); // 10  | using parametrization of the curves
	
	m = 0;

// Now we will define our curves
// For each curve we should define indices of start point and end point
// We also should define number of parametrized surfaces this curve belongs to.
// For each such surface we define the F_i and V_j, so curve is parametrized as
// u -> (u,v) = (u,V_j(u)) -> (x,y,z) = F_i(u,v)
// Pairs (i,j) are stored in separate array in successive order
// If curve is a segment, it still needs one pair of (i,j)=(0,0)
// For each pair (i,j) we also define the value of parameter u for start and end points,
// they are stored in third array
//	Here is an example how we can add curve:
/*
	LineD[3*nLine+0] = v1; LineD[3*nLine+1] = v2; // indices of start and end points
	LineD[3*nLine+2] = 2; // this curve belongs to two surfaces
	// first one have parametrization function number 1 (F_1)
	// and in this surface curve could be parametrized by V_5(u), where u are from -1.0 to 1.0
	LineP[2*m+0] = 1; LineP[2*m+1] = 5; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	// first surface have parametrization function number 2 (F_2)
	// and in this surface curve could be parametrized by V_9(u), where u are from 0.0 to 2.0
	LineP[2*m+0] = 2; LineP[2*m+1] = 9; LineT[2*m+0] =  0.0; LineT[2*m+1] =  2.0; m++;
	nLine++;
*/
// To simplify our program we can use the following macros for segments:
#define ADD_SEGMENT(V1, V2); {\
	LineD[3*nLine+0] = V1; LineD[3*nLine+1] = V2;\
	LineD[3*nLine+2] = 1;\
	LineP[2*m+0] = 0; LineP[2*m+1] = 0; LineT[2*m+0] =  0.0; LineT[2*m+1] =  0.0; m++;\
	nLine++;\
}
	// 12 edges of the box
	ADD_SEGMENT(1, 2); // 1
	ADD_SEGMENT(2, 3); // 2
	ADD_SEGMENT(3, 4); // 3
	ADD_SEGMENT(4, 1); // 4
	ADD_SEGMENT(5, 6); // 5
	ADD_SEGMENT(6, 7); // 6
	ADD_SEGMENT(7, 8); // 7
	ADD_SEGMENT(8, 5); // 8
	ADD_SEGMENT(1, 5); // 9
	ADD_SEGMENT(2, 6); // 10
	ADD_SEGMENT(3, 7); // 11
	ADD_SEGMENT(4, 8); // 12

	// First arc on equator from point 9 to point 10
	LineD[3*nLine+0] = 9; LineD[3*nLine+1] = 10;
	// This arc belongs to two surfaces (halfs of the sphere)
	LineD[3*nLine+2] = 2;
	// On both halfs of the sphere this arc is parametrized by the same (u,v)
	// We will use V_1 for both surfaces, u in [-1.0, 1.0]
	LineP[2*m+0] = 1; LineP[2*m+1] = 1; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	LineP[2*m+0] = 2; LineP[2*m+1] = 1; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
	nLine++; // 13

	// Second arc on equator from point 10 to point 9
	LineD[3*nLine+0] = 10; LineD[3*nLine+1] = 9;
	// This arc also belongs to two surfaces (halfs of the sphere)
	LineD[3*nLine+2] = 2;
	// On both halfs of the sphere this arc is parametrized by the same (u,v)
	// We will use V_2 for both surfaces, u in [1.0, -1.0]
	LineP[2*m+0] = 1; LineP[2*m+1] = 2; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	LineP[2*m+0] = 2; LineP[2*m+1] = 2; LineT[2*m+0] =  1.0; LineT[2*m+1] = -1.0; m++;
	nLine++; // 14

	m = 0;

// Now we will define surfaces
// For each surface we define 5 integer numbers (SurfL):
// 1st: number of boundary curves
// 2nd: number of parametrization function
// 3rd: color of the face
// 4th: 0 if the face is boundary face, or the second color of the face if it is used to split volume
// 5th: 0 if orientation of the face corresponds with parametrization, 1 if it should be inverted, 0 for flat faces
// For each surface we alse define minimax values of (u,v) parametrization (SurfT)
// u_min, u_max, v_min, v_max
// For each boundary curve we define a pair of integers (SurfI):
// 1st: index of curve
// 2nd: orientation,  0 for normal, 1 for reversed
//	Here is an example how we can add surface:
/*
	SurfL[5*nSurface+0] = 2; // surface is bounded by two curves
	SurfL[5*nSurface+1] = 1; // surface is parametrized by F_1
	SurfL[5*nSurface+2] = 1; // the color of face is 1
	SurfL[5*nSurface+3] = 0; // the face is boundary
	SurfL[5*nSurface+4] = 0; // orientation of the face corresponds with parametrization
	// minimax values of parametrs
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] = -1.0; SurfT[4*nSurface+2] =  1.0; SurfT[4*nSurface+3] =  1.0;
	// first curve is curve number c1, orientation is normal
	SurfI[2*m+0] = c1; SurfI[2*m+1] = 0; m++;
	// second curve is curve number c2, orientation is reversed
	SurfI[2*m+0] = c2; SurfI[2*m+1] = 1; m++;
	nSurface++;
*/
// In common case to invert orientation of the face one should invert the 5th parameter in SurfL
// and invert orientation of all boundary curves
// To simplify our program we can use the following macros for flat quadrilaterals:
// Here for each curve (V,I) is for index of the curve V and orientation I
#define ADD_QUADRILATERAL(V1,I1, V2,I2, V3,I3, V4,I4); {\
	SurfL[5*nSurface+0] = 4;\
	SurfL[5*nSurface+1] = 0;\
	SurfL[5*nSurface+2] = 1;\
	SurfL[5*nSurface+3] = 0;\
	SurfL[5*nSurface+4] = 0;\
	SurfT[4*nSurface+0] =  0.0; SurfT[4*nSurface+1] =  0.0; SurfT[4*nSurface+2] =  0.0; SurfT[4*nSurface+3] =  0.0;\
	SurfI[2*m+0] = V1; SurfI[2*m+1] = I1; m++;\
	SurfI[2*m+0] = V2; SurfI[2*m+1] = I2; m++;\
	SurfI[2*m+0] = V3; SurfI[2*m+1] = I3; m++;\
	SurfI[2*m+0] = V4; SurfI[2*m+1] = I4; m++;\
	nSurface++;\
}
	// 6 faces of the box
	ADD_QUADRILATERAL( 4,1,  3,1,  2,1,  1,1); // 1
	ADD_QUADRILATERAL( 5,0,  6,0,  7,0,  8,0); // 2
	ADD_QUADRILATERAL( 1,0, 10,0,  5,1,  9,1); // 3
	ADD_QUADRILATERAL( 2,0, 11,0,  6,1, 10,1); // 4
	ADD_QUADRILATERAL( 3,0, 12,0,  7,1, 11,1); // 5
	ADD_QUADRILATERAL( 4,0,  9,0,  8,1, 12,1); // 6

	// top half of the sphere, parametrized by F_1
	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 1;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 1;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 13; SurfI[2*m+1] = 0; m++;
	SurfI[2*m+0] = 14; SurfI[2*m+1] = 0; m++;
	nSurface++; // 7

	// bottom half of the sphere, parametrized by F_2
	SurfL[5*nSurface+0] = 2;
	SurfL[5*nSurface+1] = 2;
	SurfL[5*nSurface+2] = 1;
	SurfL[5*nSurface+3] = 0;
	SurfL[5*nSurface+4] = 0;
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] =  1.0; SurfT[4*nSurface+2] = -1.0; SurfT[4*nSurface+3] =  1.0;
	SurfI[2*m+0] = 13; SurfI[2*m+1] = 1; m++;
	SurfI[2*m+0] = 14; SurfI[2*m+1] = 1; m++;
	nSurface++; // 8
	// Note the difference of orientations of last two surfaces

// When creating new boundary representations it could be usefull to check
// the boundary of each surface
//	tria_dump_front = 1; // enable dumping of initial front for each surface
// If triangulation of surface fails, it could be usefull to check the process
// of the triangulation
//	tria_debug_front = 1; // enable dumping front on each step. Will generate a lot of files
//	fronts will be droped in files frt_gmv.{num} for gmv
//	and frt_smv.{num} for smv or manual reference
//	{num} is 000, 001, 002, ...
//	region_dump_face = 1;

	// Setting up our boundary parametrization functions
	// First one is for (x,y,z)=F_i(u,v) and second one is for v=V_j(u)
	// Now we are ready to create the surface mesh
	// Setting up our own size function.
	// Making surface mesh
	i = ani3d_surface_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
		v_u_param, surface_param, NULL, NULL, fsize,
		&nV, vertex, &nF, face, facematerial, &nnV, &nnF);
	free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);
	
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

	// It could be usefull to dump the triangulation of the surface in case we want to check
	// that boundary representation is correct and represents the desired region
	if (0) {
		write_mesh_gmv("surf.gmv", nV, vertex, nF, face, facematerial, 0, 0, 0); // for GMV
		write_front   ("surf.smv", nV, vertex, nF, face, facematerial); // for smv
//		return 0; // do not mesh the volume, just exit
	}

	// We will copy the front, so that it could be used in output in future.
	nFdup = nF;
	memcpy(facedup, face, sizeof(int)*3*nF);
	memcpy(facematdup, facematerial, sizeof(int)*nF);

	// Generate 3D mesh using our own size function fsize()
	r = mesh_3d_aft_func_(&nV, vertex, &nF, face, facematerial, &nT, tetra, tetramaterial, &nnV, &nnF, &nnT, fsize);
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

	if (r) {
		write_mesh_gmv("fail.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial);
	} else {
		// Checking that 3D mesh corresponds with surface mesh
		printf("Cheking topology: "); fflush(stdout);
		if (check_mesh_topology_(&nV, vertex, &nFdup, facedup, &nT, tetra)) printf("FAILED!\n");
		else printf("ok.\n");

		// Write output files
		write_mesh_gmv("mesh.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
		write_mesh    ("mesh.out", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
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

