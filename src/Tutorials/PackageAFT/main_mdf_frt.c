/* main_prm.c
 *
 * This is example for libfrtmdf and libaft usage
 *
 * main_mdf.exe reads front from smv file
 * We will remesh the front to improve robustness.
 * Tetra mesh is generated using mesh_3d_aft_cf_
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaft.h"
#include "libfrtmdf.h"

// in this program we use indexing from one

int main(int argc, char* argv[]) {
    int    zeroindex = 0, invert = 0;
    int    nnF = 1000000, nnV = 2000000, nnT = 4000000;
    int     nF = 0,        nV = 0,        nT = 0;
    int    *face   = 0, *facematerial  = 0;
    int    *tetra  = 0, *tetramaterial = 0;
    double *vertex = 0;
    int    nFdup,  *facedup,  *facematdup;
    int    r;
    // default values
    double cf=1.2, cfs=1.2;
    double pc[5] = {0.0, 0.05, 0.01,  0.3, 100.0};

    // reading surfacemesh from file
    if (argc<2) {
	printf("Usage: ParamToTet_aft.exe  file.frt  cf cfs  abs rel c1 c2 dist\n");
	printf("       cf - coarsening factor (default: %lg)\n", cf);
	printf("       cfs - coarsening factor on surface (default: %lg)\n", cfs);
	printf("       abs - absolute minimal mesh size (default: %lg)\n", pc[0]);
	printf("       rel - relative minimal mesh size (default: %lg)\n", pc[1]);
	printf("       c1, c2 - maximal angle deviation (default: %lg, %lg)\n", pc[2], pc[3]);
	printf("       dist - maximal distance deviation (default: %lg)\n", pc[4]);
	return 0;
    }
    if (argc>2)  sscanf(argv[2], "%lf", &cf);
    if (argc>3)  sscanf(argv[3], "%lf", &cfs);
    if (argc>4)  sscanf(argv[4], "%lf", pc+0);
    if (argc>5)  sscanf(argv[5], "%lf", pc+1);
    if (argc>6)  sscanf(argv[6], "%lf", pc+2);
    if (argc>7)  sscanf(argv[7], "%lf", pc+3);
    if (argc>8)  sscanf(argv[8], "%lf", pc+4);
    // Allocating memory for mesh structures
    vertex        = (double*)malloc(sizeof(double) * 3 * nnV);
    face          = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial  = (int*)   malloc(sizeof(int)        * nnF);
    facedup       = (int*)   malloc(sizeof(int)    * 3 * nnF); // this is for saving
    facematdup    = (int*)   malloc(sizeof(int)        * nnF); // new surface mesh
    tetra         = (int*)   malloc(sizeof(int)    * 4 * nnT);
    tetramaterial = (int*)   malloc(sizeof(int)        * nnT);


    printf("[cfs = %lg, cf = %lg]\n", cfs, cf);
    if (read_front(argv[1], &nV, vertex, &nF, face, facematerial, nnV, nnF, zeroindex, invert, 1)) {
	printf("Error in file %s\n", argv[1]);
    }

    // checking topology and merging duplicate vertices
    check_surface_topology_fix_(&nV, vertex, &nF, face);

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);
    // We can write initial mesh in GMV format here
    if (0)  write_mesh_gmv("initial.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial); // for debugging

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
    // Fifth parameter is the maximal distance from original plane for that regions
    // if any param is set to zero or less, then the default value is used
    // you can try playing with this params
    // default values:                  0,  0.05,  0.01,   0.3,  100
    surface_refine_setup_poly_extra(pc[0], pc[1], pc[2], pc[3], pc[4]); 

    surface_refine_setup_cf(cfs);

    // refine the surface mesh
    // this function does not require the front to be in cube [0,1]^3
    surface_refine_(&nV, vertex, &nF, face, facematerial, &nnV, &nnF);

    printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV, nF, nT);
    write_mesh_gmv("refined.gmv", nV, vertex, nF, face, facematerial, nT, tetra, tetramaterial); // for debugging
    write_front   ("refined.smv", nV, vertex, nF, face, facematerial); // for debugging

    // if we want to mesh external part of the space
    // we should add an additional inverted copy of our front
    // each triangle will be used twice in our front with different orientation and color
    // tetra colors are defined by the colors of the boundary front
    // we will shift colors of second inverted front by 1
    //double_front(&nF, face, facematerial, 1, nnF);

    // we should bound our second front, otherwise it will grow far-far-away
    // here we construct bounding cube with color 2, which is somewhat bigger than the model
    //add_bound_box(&nV, vertex, &nF, face, facematerial, 2, 1, 1.5, nnV, nnF);

    // saving surface mesh before 3D meshing
    // this mesh could be used in future for checking topolgy or writing to the file
    memcpy(facedup, face, sizeof(int)*3*nnF);
    memcpy(facematdup, facematerial, sizeof(int)*nnF);
    nFdup = nF;

    printf("cf = %lf\n", cf);
    // cf is the rate at which mesh coarseness changes from the boundary to the center
    // value of 1.0 should produce semi-uniform mesh
    // value of 1.2 is the default value (you can call mesh_3d_aft_ if you want to use default value)
    // value of 1.5 could be used to produce coarse mesh
    // if meshing fails try to play with this value
    // cf = 1.3;
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
	fprintf(stderr, "\nMESH: failed\n");
    }
    free(vertex),  free(face),  free(facematerial),  free(facedup),  free(facematdup),  free(tetra),  free(tetramaterial);
    return r;
}

