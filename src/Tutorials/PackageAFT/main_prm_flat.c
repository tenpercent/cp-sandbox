#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "libaft.h"
#include "libfrtprm.h"

#define SQR(X) ((X)*(X))

extern int tria_dump_front, tria_final_front, tria_debug_front, region_dump_face; // debug flags from aftlib, could be usefull for debugging

// all coords in flat_xy and height_z will be multiplied by these factors in the beginning of the program
// at the end the mesh will be unzoomed back
// this step is necessary, because surface_boundary() works better if the mesh is in [0,1]^3
double zoomx=0.001;
double zoomy=0.001;
double zoomz=0.005;

// x,y coordinates of key points
double flat_xy[14*2]={
	    0.0,     0.0, // O1
	27500.0,     0.0, // O2
	42500.0, 25000.0, // O3
	26000.0, 37000.0, // O4
	 8000.0, 34000.0, // O5
	    0.0, 17500.0, // O6
	33500.0, 10000.0, // O7
	25000.0, 10000.0, // O8
	25000.0, 30000.0, // O9
	35600.0, 30000.0, // O10

	32870.0, 19335.5,  //
	30500.0, 17941.25, // coords of
	30652.0, 17682.62, //  the box
	33022.0, 19076.88  //
};

// z coordinates of key levels
double height_z[9]={
	0, 135, 195, 245, 410, 516, 606,
	62.5, 72.5
};

// this function controls the desired size of the mesh in space
// x, y, z  -- coords of point in space
// in this function we use size = sqrt( (cr)^2 + m^2 )
// where r is the distance to center of the small box
// c is the coefficient that controls coarseness
// m is the minimal size in the center of the small box
// near the box we will have  size = m
// far from the box we will have  size = cr
double fsize(double x, double y, double z) {
#define XY(A,B) ((flat_xy[A*2+0]-x)*(flat_xy[B*2+1]-y) - (flat_xy[B*2+0]-x)*(flat_xy[A*2+1]-y))
	double h1=0.0, h2=0.0, h3=0.0, r2;
	double c  = 0.75; // coarseness
	double m2 = 16.0*SQR(height_z[7]-height_z[8]); // minimal size
	double c2 = c*c;
	
	if (h1 < XY(11, 10)) h1 = XY(11, 10);
	if (h2 < XY(10, 13)) h2 = XY(10, 13);
	if (h1 < XY(13, 12)) h1 = XY(13, 12);
	if (h2 < XY(12, 11)) h2 = XY(12, 11);
	if (h3 < z - height_z[8]) h3 = z - height_z[8];
	if (h3 < height_z[7] - z) h3 = height_z[7] - z;
	h1 /= sqrt(SQR(flat_xy[10*2+0]-flat_xy[11*2+0]) + SQR(flat_xy[10*2+1]-flat_xy[11*2+1]));
	h2 /= sqrt(SQR(flat_xy[10*2+0]-flat_xy[13*2+0]) + SQR(flat_xy[10*2+1]-flat_xy[13*2+1]));
	r2 = SQR(h1) + SQR(h2) + SQR(h3);
	return sqrt(sqrt(SQR(m2)+SQR(r2*c2)));
};


// this function removes faces with color `color`
// look at the end of main() to see if it could be usefull
int remove_color(int *pnF, int *face, int *facecolor, int color) {
	int nF=*pnF;
	int removed=0, i;
	for (i=0; i<nF; i++) {
		if (facecolor[i]==color) {
			removed++;
			nF--;
			face[3*i+0] = face[3*nF+0];
			face[3*i+1] = face[3*nF+1];
			face[3*i+2] = face[3*nF+2];
			facecolor[i] = facecolor[nF];
			i--;
		}
	}
	*pnF = nF;
	return removed;
}


// In this program we will define the model's boundary which is the set of flat polygons.
// We will use surface_boundary() function to mesh surface, and polyhedron_make_front() is the simple
// wrapper to this function. polyhedron_make_front() gets array of vertices (nV, vertex) and array of nS
// surfaces, each surface is defined by the number of edges -- nE[i] and the edges -- edge[i][],
// for each surface we define a face color -- color1[i], and if color2[i] is nonzero then it is
// used as a color of the opposite side of the surface. Basically color2[i] is used if surface
// is the surface between two materials.
// Resulting mesh will be stored in (nV, vertex, nF, face, facematerial)
// nnV and nnF is the maximum size of vertex and face respectively
//
// int polyhedron_make_front(int *pnV, double *vertex, int nS, int *nE, int **edge, int *color1, int *color2, int *pnF, int *face, int *facematerial, int nnV, int nnF);
//

// this function adds a loop to the edge array
// nE, edge -- array
// nP -- number of points in a loop
// ... -- points
// for example polyhedron_add_faceloop(&nE, edge, 4,  1, 2, 3, 4) will add edges
// (1,2), (2,3), (3,4) and (4,1)
// polyhedron_add_faceloop() could be used more then once for example if surface boundary consist of several disjoint loops
//
// int polyhedron_add_faceloop(int *pnE, int *edge, int nP, ...);
//


// this is a macros to simplify the process of adding new flat surface with one loop
// it allocates memory for edges, fills color1 and color2 (C1 and C2 respectively)
// and adds one loop (N, ...)
#define ADD(C1, C2, N, ...) {\
	nE[nL]=0;\
	edge[nL] = (int*)malloc(sizeof(int)*2*nnE);\
	c1[nL] = C1, c2[nL] = C2;\
	polyhedron_add_faceloop(&nE[nL], edge[nL], N, __VA_ARGS__);\
	nL++;\
}


// this function intersects line (v1,v2) with Oxz plane going through O7 and O8
// the result is stored in vertex v
void intersecty2(double *vertex, int v, int v1, int v2) {
	double y = flat_xy[2*7+1]; // y of O8
	double x1 = vertex[3*(v1-1)+0];
	double x2 = vertex[3*(v2-1)+0];
	double y1 = vertex[3*(v1-1)+1];
	double y2 = vertex[3*(v2-1)+1];
	double z1 = vertex[3*(v1-1)+2];
	double z2 = vertex[3*(v2-1)+2];
	vertex[3*(v-1)+0] = (x1*(y-y2) + x2*(y1-y))/(y1-y2);
	vertex[3*(v-1)+1] = y;
	vertex[3*(v-1)+2] = (z1*(y-y2) + z2*(y1-y))/(y1-y2);
}

// this function intersects line (v1,v2) with Oyz plane going through O8 and O9
// the result is stored in vertex v
void intersectx(double *vertex, int v, int v1, int v2) {
	double x = flat_xy[2*7+0]; // x of O8
	double x1 = vertex[3*(v1-1)+0];
	double x2 = vertex[3*(v2-1)+0];
	double y1 = vertex[3*(v1-1)+1];
	double y2 = vertex[3*(v2-1)+1];
	double z1 = vertex[3*(v1-1)+2];
	double z2 = vertex[3*(v2-1)+2];
	vertex[3*(v-1)+0] = x;
	vertex[3*(v-1)+1] = (y1*(x-x2) + y2*(x1-x))/(x1-x2);
	vertex[3*(v-1)+2] = (z1*(x-x2) + z2*(x1-x))/(x1-x2);
}

// this function intersects z-line going through O9 with plane (v1,v2,v3)
// the result is stored in vertex v
void intersectz(double *vertex, int v, int v1, int v2, int v3) {
	double x = flat_xy[2*8+0]; // x of O9
	double y = flat_xy[2*8+1]; // y of O9
	double x1 = vertex[3*(v1-1)+0];
	double x2 = vertex[3*(v2-1)+0];
	double x3 = vertex[3*(v3-1)+0];
	double y1 = vertex[3*(v1-1)+1];
	double y2 = vertex[3*(v2-1)+1];
	double y3 = vertex[3*(v3-1)+1];
	double z1 = vertex[3*(v1-1)+2];
	double z2 = vertex[3*(v2-1)+2];
	double z3 = vertex[3*(v3-1)+2];
	double a  = ((x2-x)*(y3-y)-(y2-y)*(x3-x))/((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
	double b  = ((x1-x)*(y3-y)-(y1-y)*(x3-x))/((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2));
	double c  = ((x1-x)*(y2-y)-(y1-y)*(x2-x))/((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3));
	vertex[3*(v-1)+0] = x;
	vertex[3*(v-1)+1] = y;
	vertex[3*(v-1)+2] = z1*a + z2*b + z3*c;;
}

// this function intersects z-line going through O10 with plane (v1,v2,v3)
// the result is stored in vertex v
void intersectz2(double *vertex, int v, int v1, int v2, int v3) {
	double x = flat_xy[2*9+0]; // x of O10
	double y = flat_xy[2*9+1]; // y of O10
	double x1 = vertex[3*(v1-1)+0];
	double x2 = vertex[3*(v2-1)+0];
	double x3 = vertex[3*(v3-1)+0];
	double y1 = vertex[3*(v1-1)+1];
	double y2 = vertex[3*(v2-1)+1];
	double y3 = vertex[3*(v3-1)+1];
	double z1 = vertex[3*(v1-1)+2];
	double z2 = vertex[3*(v2-1)+2];
	double z3 = vertex[3*(v3-1)+2];
	double a  = ((x2-x)*(y3-y)-(y2-y)*(x3-x))/((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
	double b  = ((x1-x)*(y3-y)-(y1-y)*(x3-x))/((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2));
	double c  = ((x1-x)*(y2-y)-(y1-y)*(x2-x))/((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3));
	vertex[3*(v-1)+0] = x;
	vertex[3*(v-1)+1] = y;
	vertex[3*(v-1)+2] = z1*a + z2*b + z3*c;;
}


// We want to construct mesh with different materials.
// All tetras in one solid connected part will have the same color, the minimal color of the boundary faces colors.
// In this example we want different color for each part of the model, and we want the same colors for each big boundary face
// Idea: we will color each face of each part in unique color, then we will mesh the surface, and create a copy of the front.
// One front will be used to construct tetra mesh, so we will recolor this front in such a way, that all triangle faces of one
// part will have the same color. In the process of 3D meshing this front will be destroyed.
// Another front will be used to write boundary faces in output. And we will recolor it in such a way, that all big faces will
// have the same color.

int main(int argc, char* argv[]) {
	int    nnF = 1000000, nnV = 1000000, nnT = 1000000;
	int     nF = 0,        nV = 0,        nT = 0;
	int    nnL = 1000,    nnE = 1000,     nL = 0;
	int    nE[nnL];
	int    *edge[nnL], c1[nnL], c2[nnL];
	int    *face   = 0, *facematerial  = 0, *facedup = 0, *facematdup = 0;
	int    *tetra  = 0, *tetramaterial = 0;
	double *vertex = 0;
	int    i, nFdup, r;

	// From now each color will be represented in form a+20*b, where b is the number of layer (0..6).
	// Layers from 0 to 5 are the real six layers, and layer 6 is in fact the small box inside layer 0.
	// fcolormap is used to recolor copy of the front for boundary face colors, we will use formula
	// newcolor = fcolormap[color % 20], so all big faces will have the same color
	// New colors 1..8 is the side colors, 9 is the bottom color, 10 is the top color.
	// Some old colors will be recolored to zero, for example old colors 9..11 will be used only on internal
	// surfaces between two solid parts. In future we can remove them from output, if we don't need them
	// colors with a=0,14 will not be used, they are mapped to '-1'
	int    fcolormap[20]={-1, 2, 3, 3, 3, 2, 2, 2, 2, 5, 5, 5, 3, 3, -1, 1, 1, 4, 4, 4};
	int    frecolormap[36]={0, 1, 2, 2, 0, 0,
				0, 4, 3, 0, 0,
				0, 6, 5, 0, 0,
				0, 8, 7, 0, 0,
				0, 9, 9, 11, 0,
				0, 10, 10, 11, 0,
				0, 0, 0, 0, 0};
	// tcolormap and tcolorbase are used to recolor front before 3D meshing
	// newcolor = tcolorbase[color/20] + tcolormap[color%20]
	// in each layer we will have one, two or three colors
	int    tcolormap[20]={-1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 3, 3, 3, -1, 1, 2, 1, 2, 3};
	int    tcolorbase[7]={0, 2, 4, 6, 8, 11, 13};

	int np=66; // it is a number of used points in first phase
	(void) argc, (void) argv;


	// allocate memory for mesh ctructures
	vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
	face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF);
	facematerial  = (int*)   malloc(sizeof(int)        * nnF);
	facematdup    = (int*)   malloc(sizeof(int)        * nnF);
	tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
	tetramaterial = (int*)   malloc(sizeof(int)        * nnT);

	// zoom in all coords in flat_xy and height_z
	for (i=0; i<14; i++) {
		flat_xy[2*i+0]*=zoomx;
		flat_xy[2*i+1]*=zoomy;
	}
	for (i=0; i<9; i++) {
		height_z[i]*=zoomz;
	}

	// fill vertex array
	// here vertices 1..10 are O1..O10 in plane z=z0
	// 11..20 are O1..O10 in plane z=z1, etc
	nV = np;
	for (i=0; i<nV; i++) {
		vertex[3*i+0] = flat_xy[(i%10)*2+0];
		vertex[3*i+1] = flat_xy[(i%10)*2+1];
		vertex[3*i+2] = height_z[i/10];
	}
	// NOTE: some points will not be used, and we will redefine them later
	// the actual values of vertex array are used only in polyhedron_make_front() function
	// before calling this function we can define vertex coords and topology in any order
	// but in this example we will first define coords, and than use vertices in topology
	// to avoid confusion

	// Now we will use a great amount of ADD() routines. Each ADD(C1, C2, ...) define one flat region
	// with colors C1 and C2 for front side and back side respectively

	// colors 15 and 16 will be mapped to bottom color at the output
	// here we define two bottom faces, one for each solid part
	// 15 will be mapped to tetracolor 1, and 16 -> 2
	// the bottom surfaces are not internal, so second color is set to zero
	// the order of points in loop defines the orientation of the surface
	// first loop will be O1, O6, O5, O4, O10, O9, O8, O7, O2 (9 points)
	ADD(15, 0, 9,    1,  6,  5,  4, 10,  9,  8,  7,  2);
	// second loop will be O7, O8, O9, O10, O3 (5 points)
	ADD(16, 0, 5,    7,  8,  9, 10,  3);

	// Now we want to construct side faces, internal faces and top faces for first layer.
	// In fact top faces will be used as bottom faces for next layer.
	// Since the topology of first 3 layers is the same we can use a for-statement
	for (i=0; i<3; i++) {
		// Colors 9, 10 and 11 in each layer will be used as internal colors for tetra materials 1, 2, 3 respectively
		// That is why if and only if the first color of surface is in [9..11] then the second is nonzero and in [9..11] too.
		// In all other cases if color1 is not in [9..11] then this surface is the boundary surface and color2 should be zero.

		// Here we define top faces of the layer, the second colors show that these faces will be used in the next layer
		// Note, that the order of these loops differs from the bottommost ones. This is because it is a top face, not bottom.
		// Color 9 is for first part, and 10 -- for second
		ADD( 9+20*i,  9+20*(i+1), 9,   11+10*i, 12+10*i, 17+10*i, 18+10*i, 19+10*i, 20+10*i, 14+10*i, 15+10*i, 16+10*i);
		ADD(10+20*i, 10+20*(i+1), 5,   17+10*i, 13+10*i, 20+10*i, 19+10*i, 18+10*i);
		
		// Now we define side faces, each with unique color
		ADD( 8+20*i, 0, 4,   11+10*i, 16+10*i,  6+10*i,  1+10*i);
		ADD( 7+20*i, 0, 4,   16+10*i, 15+10*i,  5+10*i,  6+10*i);
		ADD( 6+20*i, 0, 4,   15+10*i, 14+10*i,  4+10*i,  5+10*i);
		ADD( 5+20*i, 0, 4,   14+10*i, 20+10*i, 10+10*i,  4+10*i);
		ADD( 2+20*i, 0, 4,   17+10*i, 12+10*i,  2+10*i,  7+10*i);
		ADD( 1+20*i, 0, 4,   12+10*i, 11+10*i,  1+10*i,  2+10*i);
		ADD( 4+20*i, 0, 4,   20+10*i, 13+10*i,  3+10*i, 10+10*i);
		ADD( 3+20*i, 0, 4,   13+10*i, 17+10*i,  7+10*i,  3+10*i);

		// Here we define the internal surfaces between two parts
		// Color 9 is for first part, and 10 -- for second
		ADD( 9+20*i, 10+20*i, 4,   20+10*i, 19+10*i,  9+10*i, 10+10*i);
		ADD( 9+20*i, 10+20*i, 4,   19+10*i, 18+10*i,  8+10*i,  9+10*i);
		ADD( 9+20*i, 10+20*i, 4,   18+10*i, 17+10*i,  7+10*i,  8+10*i);
	}


	// Now we want to construct layer 3. This layer is like 0, 1 and 2.
	// But the next layer 4 will have three materials, so his bottom face should be
	// splitted additionaly. This means that the top face of the layer 3 should be also splitted.
	// Here we need an additional point: the intersection of (48,49) and (41,43),
	// ie O8O9 and O1O3 in plane z=z4.
	// Point number 51 will be never used, so we can use it for this purpose.
	intersectx(vertex, 51, 41, 43);

	// Here we contruct the top face, two loops for each part
	// Color 11 in layer 4 will be used for 3rd part.
	ADD( 9+20*3, 11+20*4, 5,   41, 42, 47, 48, 51)
	ADD( 9+20*3,  9+20*4, 7,   41, 51, 49, 50, 44, 45, 46)
	ADD(10+20*3, 11+20*4, 4,   51, 48, 47, 43);
	ADD(10+20*3, 10+20*4, 4,   43, 50, 49, 51);
	// Side faces
	ADD( 8+20*3, 0, 4,   41, 46, 36, 31);
	ADD( 7+20*3, 0, 4,   46, 45, 35, 36);
	ADD( 6+20*3, 0, 4,   45, 44, 34, 35);
	ADD( 5+20*3, 0, 4,   44, 50, 40, 34);
	ADD( 2+20*3, 0, 4,   32, 37, 47, 42);
	ADD( 1+20*3, 0, 4,   42, 41, 31, 32);
	ADD( 4+20*3, 0, 4,   50, 43, 33, 40);
	ADD( 3+20*3, 0, 4,   43, 47, 37, 33);
	// Internal faces
	ADD( 9+20*3, 10+20*3, 4,   50, 49, 39, 40);
	ADD( 9+20*3, 10+20*3, 5,   49, 51, 48, 38, 39);
	ADD( 9+20*3, 10+20*3, 4,   48, 47, 37, 38);

	
	// Layer 4
	// First we will construct the 3rd part in this layer
	// We will need two additional points: 52 and 57.
	// 52 is like 51, but shifted up a litle in such way, that 41, 52 and 53 will be on one line
	// 57 is like 47, but also shifted up in such way, that 42, 57 and 53 will be in one line
	intersectx (vertex, 52, 41, 53);
	intersecty2(vertex, 57, 42, 53);
	// The top face, color=19 is mapped to 10, in fact this is a triangle, but two sides are splitted by 52 and 57
	ADD(19+20*4, 0, 5,   41, 42, 57, 53, 52);
	// The side face, colors 12 and 13 will be mapped to 2 and 3 respectively
	// We use 12 and 13 instead of 2 and 3 because these faces belongs to the third solid part and
	// these colors should be recolored into 3rd color for tetra
	ADD(12+20*4, 0, 3,   47, 57, 42);
	ADD(13+20*4, 0, 4,   47, 43, 53, 57);
	// And the internal faces for first and second solid parts
	ADD(11+20*4, 10+20*4, 4,   43, 51, 52, 53);
	ADD(11+20*4,  9+20*4, 3,   51, 41, 52);

	// Now we will construct the 1st part in layer 4
	// We will need point 58 -- intersection of (56,53) with Oyz plane
	intersectx(vertex, 58, 56, 53);
	// Top faces, they are the bottom faces for next layer
	ADD( 9+20*4, 9+20*5, 4,   41, 52, 58, 56);
	ADD( 9+20*4, 9+20*5, 6,   56, 58, 59, 60, 54, 55);
	// Side faces
	ADD( 8+20*4, 0, 3,   56, 46, 41);
	ADD( 7+20*4, 0, 4,   45, 46, 56, 55);
	ADD( 6+20*4, 0, 4,   44, 45, 55, 54);
	ADD( 5+20*4, 0, 4,   50, 44, 54, 60);
	// Internal faces
	ADD( 9+20*4, 10+20*4, 4,   49, 50, 60, 59);
	ADD( 9+20*4, 10+20*4, 5,   59, 58, 52, 51, 49);

	// And the last, 2nd part in layer 4
	// Top faces
	ADD(10+20*4, 10+20*5, 3,   52, 53, 58);
	ADD(10+20*4, 10+20*5, 4,   58, 53, 60, 59);
	// Side face
	ADD( 4+20*4, 0, 4,   53, 43, 50, 60);

	// Layer 5
	// We will need point 61, which is a little upper than 58
	intersectx(vertex, 61, 66, 53);
	// Point 62 should be in plane (66,64,53) right on top of 59
	intersectz(vertex, 62, 66, 64, 53);
	// Point 63 should be in plane (66,64,53) a little upper than 60
	intersectz2(vertex, 63, 66, 64, 53);

	// The top faces of 1st part, color=17, mapped to 10 (top color)
	ADD(17+20*5, 0, 4,   41, 52, 61, 66);
	ADD(17+20*5, 0, 5,   66, 61, 62, 63, 64);
	ADD(17+20*5, 0, 3,   64, 65, 66);
	// Side faces
	ADD( 8+20*5, 0, 3,   66, 56, 41);
	ADD( 7+20*5, 0, 4,   55, 56, 66, 65);
	ADD( 6+20*5, 0, 4,   54, 55, 65, 64);
	ADD( 5+20*5, 0, 4,   60, 54, 64, 63);
	// Internal faces
	ADD( 9+20*5, 10+20*5, 4,   59, 60, 63, 62);
	ADD( 9+20*5, 10+20*5, 5,   62, 61, 52, 58, 59);

	// Top faces for 2nsd part, color=18, mapped to 10.
	ADD(18+20*5, 0, 3,   52, 53, 61);
	ADD(18+20*5, 0, 4,   61, 53, 63, 62);
	// Side face
	ADD( 4+20*5, 0, 3,   60, 63, 53);

	// Now we have a 6-layer domain, which could be successfully meshed into tetras.
	// But we should also add a small box in layer 0.
	// NOTE: all points from 1 to np=66 are used now.

	// Phase 2. Adding internal small box

	// Adding new 8 points
	for (i=0; i<8; i++) {
		vertex[3*(np+i)+0] = flat_xy[2*(10+i%4)+0];
		vertex[3*(np+i)+1] = flat_xy[2*(10+i%4)+1];
		vertex[3*(np+i)+2] = height_z[7+i/4];
	}
	nV += 8;
	// Adding new 6 faces. All these faces are internal.
	// Internal color is 9 for layer 6, and external is 10 for layer 0.
	// Color 9+20*6 will be mapped to tetra color 13+1=14
	ADD(9+20*6, 10, 4,   np+4, np+3, np+2, np+1);
	ADD(9+20*6, 10, 4,   np+5, np+6, np+7, np+8);
	ADD(9+20*6, 10, 4,   np+1, np+2, np+6, np+5);
	ADD(9+20*6, 10, 4,   np+2, np+3, np+7, np+6);
	ADD(9+20*6, 10, 4,   np+3, np+4, np+8, np+7);
	ADD(9+20*6, 10, 4,   np+4, np+1, np+5, np+8);

// When creating new boundary representations it could be usefull to check
// the boundary of each surface
//	tria_dump_front = 1; // enable dumping of initial front for each surface
//	tria_final_front = 1; // enable dumping of final front for each surface
// If triangulation of surface fails, it could be usefull to check the process
// of the triangulation
//	tria_debug_front = 1; // enable dumping front on each step. Will generate a lot of files
//	fronts will be droped in files frt_gmv.{num} for gmv
//	and frt_smv.{num} for smv or manual reference
//	{num} is 000, 001, 002, ...
//	region_dump_face = 1;

	// Now we are ready to create the surface mesh
	// Setting up our own size function.
	set_surface_size_function(fsize);
	// Making surface mesh
	polyhedron_make_front(&nV, vertex, nL, nE, edge, c1, c2, &nF, face, facematerial, nnV, nnF);

	// We will copy the front, so that it could be used in output in future.
	nFdup = nF;
	memcpy(facedup, face, sizeof(int)*3*nF);
	memcpy(facematdup, facematerial, sizeof(int)*nF);

	// Recoloring front to tetra colors.
	for (i=0; i<nF; i++) facematerial[i] = tcolorbase[facematerial[i]/20] + tcolormap[facematerial[i]%20];
	// Recoloring copy of the front to side faces colors
	for (i=0; i<nFdup; i++) facematdup[i] = frecolormap[5*(facematdup[i]/20) + fcolormap[facematdup[i]%20]];

	// Generate 3D mesh using our own size function fsize()
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);
	r = mesh_3d_aft_func_(&nV, vertex, &nF, face, facematerial, &nT, tetra, tetramaterial, &nnV, &nnF, &nnT, fsize);
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

	if (r) {
		write_mesh_gmv("fail.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial);
	} else {
		// Checking that 3D mesh corresponds with surface mesh
		printf("Checking topology: "); fflush(stdout);
		if (check_mesh_topology_(&nV, vertex, &nFdup, facedup, &nT, tetra)) printf("FAILED!\n");
		else printf("ok.\n");

		// Now we can remove internal faces with color=11, if we don't want them
		if (1) remove_color(&nFdup, facedup, facematdup, 0);

		// Zooming out
		for (i=0; i<nV; i++) {
			vertex[3*i+0] /= zoomx;
			vertex[3*i+1] /= zoomy;
			vertex[3*i+2] /= zoomz;
		}
		// Writing final result
		write_mesh_gmv("mesh.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
		write_mesh    ("mesh.out", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
	}
	
	
	return 0;
}

