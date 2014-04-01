/* main_aft.c
 *
 * This is example for libfrtcad and libaft usage
 *
 * main_cad.exe reads brep file and generates initial front
 * Then it generates tetra mesh using mesh_3d_aft_cf_
 * 
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cgmwrap.h"
#include "libfrtcad.h"
#include "libaft.h"

extern double surfacesizeratio;
extern int region_dump_face_cgm;
extern int tria_dump_front, tria_final_front, tria_debug_front, region_dump_face; // debug flags from aftlib, could be usefull for debugging

int main(int argc, char *argv[]) {
    CGMmodel model;
    int    nnF = 1000000, nnV = 2000000, nnT = 4000000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0, *facedup = 0, *facematdup = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    nFdup, r;
    double cf = 1.2;
    char   *fext;
    (void) argc;

    // Initialize CGM
    CGM_Init();

    if (argc<2) {
	printf("Usage: CadToTet_aft.exe file.brep  [cf]  [size]\n");
	printf("       cf - coarsening factor (default: 1.2)\n");
	printf("       size - relative mesh size (default: 0.1)\n");
	return 0;
    }
    // Load CGM model from file

    // Check for file extension
    fext = argv[1] + strlen(argv[1]);
    while ((fext > argv[1]) && (fext[-1] != '.'))  fext--;
    if (strcmp(fext, "brep") == 0)  model = CGM_LoadOCCModelFromFile(argv[1], "OCC");
    else if (strcmp(fext, "rle") == 0)  model = CGM_LoadOCCModelFromFile(argv[1], "OCC");
    else if (strcmp(fext, "step") == 0)  model = CGM_LoadOCCModelFromFile(argv[1], "STEP");
    else if (strcmp(fext, "iges") == 0)  model = CGM_LoadOCCModelFromFile(argv[1], "IGES");
    else  model = CGM_LoadModelFromFile(argv[1]); // unknown file type, try default (OCC)

    // Coarsening factor
    if (argc>2)  sscanf(argv[2], "%lf", &cf);

    // Allocate memory
    vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
    face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial  = (int*)   malloc(sizeof(int)        * nnF);
    facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematdup    = (int*)   malloc(sizeof(int)        * nnF);
    tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
    tetramaterial = (int*)   malloc(sizeof(int)        * nnT);

    // this is the relative mesh size for the front
    surfacesizeratio = 0.1;
    if (argc>3)  sscanf(argv[3], "%lf", &surfacesizeratio);
    printf("surfacesizeratio=%lf\n", surfacesizeratio);

    // When creating new boundary representations it could be usefull to check
    // the boundary of each surface
    // tria_dump_front = 1; // enable dumping of initial front for each surface
    // tria_final_front = 1; // enable dumping of final front for each surface
    // If triangulation of surface fails, it could be usefull to check the process
    // of the triangulation
    // tria_debug_front = 1; // enable dumping front on each step. Will generate a lot of files
    // fronts will be droped in files frt_gmv.{num} for gmv
    // and frt_smv.{num} for smv or manual reference
    // {num} is 000, 001, 002, ...
    // region_dump_face = 1;
    region_dump_face_cgm = 1;

    // creating initial front
    ani3d_surface_cgm_model_(model, NULL, &nV, vertex, &nF, face, facematerial, &nnV, &nnF);
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

