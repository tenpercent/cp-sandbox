/* main_aft.c
 *
 * This is example for libaft usage
 *
 * main_aft.exe reads boundary triangulation from file specified in command line
 * Then it generates tetra mesh using mesh_3d_aft_cf_
 * 
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaft.h"

int main(int argc, char* argv[]) {
    int    zeroindex = 0, invert = 0;
    int    nnF = 1000000, nnV = 2000000, nnT = 4000000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    nFdup,  *facedup,  *facematdup;
    int    r;
    double cf=1.2;

    if (argc<2) {
	printf("Usage: FrontToTet_aft.exe file.frt  [cf]\n");
	printf("       cf - coarsening factor (default: 1.2)\n");
	return 0;
    }
    if (argc>2) {
	sscanf(argv[2], "%lf", &cf);
    }

    // Allocating memory for mesh structures
    vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
    face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial  = (int*)   malloc(sizeof(int)        * nnF);
    facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematdup    = (int*)   malloc(sizeof(int)        * nnF);
    tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
    tetramaterial = (int*)   malloc(sizeof(int)        * nnT);


    // Reading front from file argv[1]
    // File format is:
    // nV nF
    // X_1 Y_1 Z_1
    // ...
    // X_nV Y_nV Z_nV
    // i_1 j_1 k_1  color_1
    // ...
    // i_nF j_nF k_nF color_nF
    //
    // (indexing starts from 1)
    if (read_front(argv[1], &nV, vertex, &nF, face, facematerial, nnV, nnF, zeroindex, invert, 1)) {
	printf("Error in file %s\n", argv[1]);
    }

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);

    // We will copy the front, so that it could be used in output in future
    memcpy(facedup, face, sizeof(int)*3*nF);
    memcpy(facematdup, facematerial, sizeof(int)*nF);
    nFdup = nF;

    printf("cf = %lf\n", cf);

    // Generate 3D mesh
    r = mesh_3d_aft_cf_(&nV, vertex, &nF, face, facematerial, &nT, tetra, tetramaterial, &nnV, &nnF, &nnT, &cf);
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
	fprintf(stderr, "\nMESH: failed\n");
    }
    free(vertex),  free(face),  free(facematerial),  free(facedup),  free(facematdup),  free(tetra),  free(tetramaterial);
    return r;
}

