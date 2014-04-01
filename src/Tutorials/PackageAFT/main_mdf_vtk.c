/* main_prm.c
 *
 * This is example for libfrtmdf and libaft usage
 *
 * main_mdf.exe reads front from vtk file
 * We will remesh the front to improve robustness.
 * We will mesh interior and exterior of this model.
 * Front will be doubled, one of this copy will be inverted.
 * We will add bounding box for exterior part.
 * Tetra mesh is generated using mesh_3d_aft_cf_
 * 
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaft.h"
#include "libfrtmdf.h"

// in this program we use indexing from one

// user-defined function for reading vtk format
// this function is implemented at the end of the file
int read_vtk(char *fn, int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int nnV, int nnF);

// for each face in the front make a copy and invert it
// new face colors are shifted by 'faceinc'
int double_front(int *pnF, int *face, int *facemat, int faceinc, int nnF) {
    int nF=*pnF;
    int i;
    // checking if we have enough memory
    if (2*nF>nnF) {
	printf("Increase nnF\n");
	return 1;
    }
    for (i=0; i<nF; i++) {
	face[3*(nF+i)+0] = face[3*i+2];
	face[3*(nF+i)+1] = face[3*i+1];
	face[3*(nF+i)+2] = face[3*i+0];
	facemat[nF+i] = facemat[i] + faceinc;
    }
    *pnF = 2*nF;
    return 0;
}

// compute and add bounding box surface mesh to the front
// bounding box is computed on information of vertices coords, even if vertex is not used in any face
// 'color' is the color of new faces
// 'cube'  checks if bounding box should be cube
// 'size'  determines the rate of enlargement, 1.0 -- exactly the bounding box, 1.5 -- bugger by 50%
int add_bound_box(int *pnV, double *vertex, int *pnF, int *face, int *facemat, int color, int cube, double size, int nnV, int nnF) {
    int i=0, nV=*pnV, nF=*pnF;
    double xmin, xmax, ymin, ymax, zmin, zmax, r, x, y, z;

    // checking if we have enough memory
    if (nV+8>nnV) {
	printf("Increase nnV\n");
	return 1;
    }
    if (nF+12>nnF) {
	printf("Increase nnF\n");
	return 1;
    }
    if (nV==0) return 0;

    // calculating bounding box [xmin,xmax]*[ymin,ymax]*[zmin,zmax]
    x = vertex[3*i+0], y = vertex[3*i+1], z = vertex[3*i+2];
    xmin = xmax = x;
    ymin = ymax = y;
    zmin = zmax = z;

    for (i=1; i<nV; i++) {
	x = vertex[3*i+0], y = vertex[3*i+1], z = vertex[3*i+2];
	if (xmin > x) xmin = x;
	if (xmax < x) xmax = x;
	if (ymin > y) ymin = y;
	if (ymax < y) ymax = y;
	if (zmin > z) zmin = z;
	if (zmax < z) zmax = z;
    }

    x = xmax + xmin, y = ymax + ymin, z = zmax + zmin;
    if (cube) {
	// enlarging the box to become cube
	r = xmax - xmin;
	if (r < ymax - ymin) r = ymax - ymin;
	if (r < zmax - zmin) r = zmax - zmin;
	r *= size;
	xmin = (x-r)/2.0; xmax = (x+r)/2.0;
	ymin = (y-r)/2.0; ymax = (y+r)/2.0;
	zmin = (z-r)/2.0; zmax = (z+r)/2.0;
    } else {
	// enlarging the box
	r = (xmax - xmin)*size;
	xmin = (x-r)/2.0; xmax = (x+r)/2.0;
	r = (ymax - ymin)*size;
	ymin = (y-r)/2.0; ymax = (y+r)/2.0;
	r = (zmax - zmin)*size;
	zmin = (z-r)/2.0; zmax = (z+r)/2.0;
    }

    // adding new 8 vertices
    vertex[3*(nV+0)+0]=xmin, vertex[3*(nV+0)+1]=ymin, vertex[3*(nV+0)+2]=zmin;
    vertex[3*(nV+1)+0]=xmax, vertex[3*(nV+1)+1]=ymin, vertex[3*(nV+1)+2]=zmin;
    vertex[3*(nV+2)+0]=xmin, vertex[3*(nV+2)+1]=ymax, vertex[3*(nV+2)+2]=zmin;
    vertex[3*(nV+3)+0]=xmax, vertex[3*(nV+3)+1]=ymax, vertex[3*(nV+3)+2]=zmin;
    vertex[3*(nV+4)+0]=xmin, vertex[3*(nV+4)+1]=ymin, vertex[3*(nV+4)+2]=zmax;
    vertex[3*(nV+5)+0]=xmax, vertex[3*(nV+5)+1]=ymin, vertex[3*(nV+5)+2]=zmax;
    vertex[3*(nV+6)+0]=xmin, vertex[3*(nV+6)+1]=ymax, vertex[3*(nV+6)+2]=zmax;
    vertex[3*(nV+7)+0]=xmax, vertex[3*(nV+7)+1]=ymax, vertex[3*(nV+7)+2]=zmax;

    // adding new 12 faces, very simple mesh
    face[3*(nF+ 0)+0]=nV+1, face[3*(nF+ 0)+1]=nV+3, face[3*(nF+ 0)+2]=nV+4; facemat[nF+ 0]=color;
    face[3*(nF+ 1)+0]=nV+4, face[3*(nF+ 1)+1]=nV+2, face[3*(nF+ 1)+2]=nV+1; facemat[nF+ 1]=color;
    face[3*(nF+ 2)+0]=nV+7, face[3*(nF+ 2)+1]=nV+5, face[3*(nF+ 2)+2]=nV+6; facemat[nF+ 2]=color;
    face[3*(nF+ 3)+0]=nV+6, face[3*(nF+ 3)+1]=nV+8, face[3*(nF+ 3)+2]=nV+7; facemat[nF+ 3]=color;
    face[3*(nF+ 4)+0]=nV+6, face[3*(nF+ 4)+1]=nV+2, face[3*(nF+ 4)+2]=nV+4; facemat[nF+ 4]=color;
    face[3*(nF+ 5)+0]=nV+4, face[3*(nF+ 5)+1]=nV+8, face[3*(nF+ 5)+2]=nV+6; facemat[nF+ 5]=color;
    face[3*(nF+ 6)+0]=nV+5, face[3*(nF+ 6)+1]=nV+1, face[3*(nF+ 6)+2]=nV+2; facemat[nF+ 6]=color;
    face[3*(nF+ 7)+0]=nV+2, face[3*(nF+ 7)+1]=nV+6, face[3*(nF+ 7)+2]=nV+5; facemat[nF+ 7]=color;
    face[3*(nF+ 8)+0]=nV+7, face[3*(nF+ 8)+1]=nV+3, face[3*(nF+ 8)+2]=nV+1; facemat[nF+ 8]=color;
    face[3*(nF+ 9)+0]=nV+1, face[3*(nF+ 9)+1]=nV+5, face[3*(nF+ 9)+2]=nV+7; facemat[nF+ 9]=color;
    face[3*(nF+10)+0]=nV+8, face[3*(nF+10)+1]=nV+4, face[3*(nF+10)+2]=nV+3; facemat[nF+10]=color;
    face[3*(nF+11)+0]=nV+3, face[3*(nF+11)+1]=nV+7, face[3*(nF+11)+2]=nV+8; facemat[nF+11]=color;

    *pnV = nV + 8;
    *pnF = nF + 12;
    return 0;
}


int main(int argc, char* argv[]) {
    int    nnF = 1000000, nnV = 2000000, nnT = 4000000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    nFdup,  *facedup,  *facematdup;
    int    r;
    double cf=1.2;

    // allocating memory
    vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
    face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial  = (int*)   malloc(sizeof(int)        * nnF);
    facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF); // this is for saving
    facematdup    = (int*)   malloc(sizeof(int)        * nnF); // new surface mesh
    tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
    tetramaterial = (int*)   malloc(sizeof(int)        * nnT);

    // reading surfacemesh from file
    if (argc<2) {
	printf("Usage: main_mdf.exe file.vtk\n");
	return 0;
    }
    if (read_vtk(argv[1], &nV, vertex, &nF, face, facematerial, nnV, nnF)) return 1;

    // cheking topology and merging duplicate vertices
    check_surface_topology_fix_(&nV, vertex, &nF, face);

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

    // setting parameters for surface_refine
    // Each edge could be splitted into multiple segments
    // first param bounds the mininal absolute length of these segments
    // if first param is set to zero (default value) then absolute length is not bounded below 
    // second param bounds relative length of segments
    // (value of 0.1 mean that each segment should be at least 0.1*edge_length)
    // Surface is splitted into almost flat regions, next two params control the degree of this flatness
    // 3rd and 4th params control maximum angle between normals of triangle faces
    // 3rd is for regular triangles
    // 4th is for low-quality triangles
    // if both params are nearly zero, then all regions are really flat
    // these values are something like (1 - cos \alpha)
    // if any param is set to zero or less, then the default value is used
    // you can try playing with this params
    // default values:           0, 0.05, 0.01,  0.3
    surface_refine_setup_poly(0.05, 0.05, 1e-3, 1e-2); // for bolt+nut example
    // cf is the rate at which mesh coarseness changes from the boundary to the center
    // surface_refine_setup_cf setups coarsening rate for new surface mesh
    surface_refine_setup_cf(1.3);

    // refine the surface mesh
    // this function does not require the front to be in cube [0,1]^3
    surface_refine_(&nV, vertex, &nF, face, facematerial, &nnV, &nnF);

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);
    //	write_mesh_gmv("refined.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial); // for debugging

    // if we want to mesh external part of the space
    // we should add an additional inverted copy of our front
    // each triangle will be used twice in our front with different orientation and color
    // tetra colors are defined by the colors of the boundary front
    // we will shift colors of second inverted front by 1
    double_front(&nF, face, facematerial, 1, nnF);

    // we should bound our second front, otherwise it will grow far-far-away
    // here we construct bounding cube with color 2, which is somewhat bigger than the model
    add_bound_box(&nV, vertex, &nF, face, facematerial, 2, 1, 1.5, nnV, nnF);

    // saving surface mesh before 3D meshing
    // this mesh could be used in future for checking topolgy or writing to the file
    memcpy(facedup, face, sizeof(int)*3*nnF);
    memcpy(facematdup, facematerial, sizeof(int)*nnF);
    nFdup = nF;

    // cf is the rate at which mesh coarseness changes from the boundary to the center
    // value of 1.0 should produce semi-uniform mesh
    // value of 1.2 is the default value (you can call mesh_3d_aft_ if you want to use default value)
    // value of 1.5 could be used to produce coarse mesh
    // if meshing fails try to play with this value
    cf = 1.3;
    // make 3D mesh from the front
    // this function could leave something unmeshed if advancing front fails
    // this function does not require the front to be in cube [0,1]^3
    r = mesh_3d_aft_cf_(&nV, vertex, &nF, face, facematerial, &nT, tetra, tetramaterial, &nnV, &nnF, &nnT, &cf);
    // r==0 - ok, complete meshing
    // r!=0 - something is left, (nF, face, facematerial) is the surface mesh (aka front) of the rest

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

    fprintf(stderr, "\nAFT returned: %s\n", (r)?"FAIL":"ok");

    if (r==0) {
	// Checking that 3D mesh corresponds with surface mesh
	printf("Checking topology: "); fflush(stdout);
	if ((r = check_mesh_topology_(&nV, vertex, &nFdup, facedup, &nT, tetra)))  printf("FAILED!\n"),  fprintf(stderr, "\nMESH: failed\n");
	else printf("ok.\n"),  fprintf(stderr, "\nMESH: ok\n");

	// Write output files
	write_mesh_gmv("mesh.gmv", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
	write_mesh    ("mesh.out", nV, vertex, nFdup, facedup, facematdup, nT, tetra, tetramaterial);
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nFdup, nT);
    } else {
	// Meshing failed. We save current meshed part and current front in gmv format
	printf("FAILED!\n");
	write_mesh_gmv("fail.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial);
    }
    free(vertex),  free(face),  free(facematerial),  free(facedup),  free(facematdup),  free(tetra),  free(tetramaterial);
    return r;
}



// Reads surface mesh from .vtk file
// fn      - filename
// nV      - number of vertices
// vertex  - nV*3 array of coords
// nF      - number of faces
// face    - nF*3 array of face indices (in [1..nV])
// facemat - array of face material (filled with 1 in this example)
// nnV     - max_nV
// nnF     - max_nF
#define BUFSIZE 1024
int read_vtk(char *fn, int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int nnV, int nnF) {
    FILE *f;
    char buf[BUFSIZE];
    int i, r, nv, nf;

    if (!(f = fopen(fn, "r"))) {
	perror(fn);
	return 1;
    }

    // reading vertices
    do {
	fgets(buf, BUFSIZE, f);
    } while (sscanf(buf, "POINTS %d", &nv)<1);
    // cheking if we have enough memory
    if (nv > nnV) {
	printf("nV > max_nV!\n");
	return 2;
    }
    // reading coordinates
    for (i=0; i<nv; i++) {
	fgets(buf, BUFSIZE, f);
	r = sscanf(buf, "%lf %lf %lf", vertex+3*i+0, vertex+3*i+1, vertex+3*i+2);
	if (r<3) {
	    printf("%s\nexpected \"x y z\"\n", buf);
	}
    }

    // reading triangles
    do {
	fgets(buf, BUFSIZE, f);
    } while (sscanf(buf, "POLYGONS %d", &nf)<1);
    // checking if we have enough memory
    if (nf > nnF) {
	printf("nF > max_nF!\n");
	return 2;
    }
    // reading indices
    for (i=0; i<nf; i++) {
	fgets(buf, BUFSIZE, f);
	r = sscanf(buf, "3 %d %d %d", face+3*i+0, face+3*i+1, face+3*i+2);
	if (r<3) {
	    printf("%s\nexpected \"3 i1 i2 i3\"\n", buf);
	}
	facematerial[i] = 1;

	// we increase indices here because we want them to be from 1 to nV
	// and if file they are from 0 to nV-1
	face[3*i+0]++;
	face[3*i+1]++;
	face[3*i+2]++;
    }

    // updating number of vertices and faces
    *pnV = nv;
    *pnF = nf;

    return 0;
}
#undef BUFSIZE


