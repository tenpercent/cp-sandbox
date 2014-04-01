/* main_scg.cpp
 *
 * This is an example of libfrtscg library usage.
 * 
 * main_scg.exe generates surface meshes for primitives
 * (parallelepiped, sphere or cylinder), moves, rotates and scales them,
 * generates surface mesh for intersection/union/difference and writes it
 * to the file.
 * Then it generates tetra mesh using mesh_3d_aft_
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "libfrtscg.h"
#include "libprim.h"
#include "libaft.h"

int main ()
{
        int nnV = 1000000, nnF = 1000000, nnT = 1000000;
	int i, r;
        int nV_in1 = 0, nF_in1 = 0, nV_in2 = 0, nF_in2 = 0, nV_out = 0, nF_out = 0, nF_3d = 0, nT_out = 0;

        double *Vertex_in1 = (double*)malloc (sizeof (double) * 3 * nnV);
	int *Index_in1 = (int*)malloc (sizeof (int) * 3 * nnF);
	
	double *Vertex_in2 = (double*)malloc (sizeof (double) * 3 * nnV);
	int *Index_in2 = (int*)malloc (sizeof (int) * 3 * nnF);

	double *Vertex_out = (double*)malloc (sizeof (double) * 3 * nnV);
	int *Index_out = (int*)malloc (sizeof (int) * 3 * nnF);
	int *Tetra_out = (int*)malloc (sizeof (int) * 4 * nnT);

	int *Index_3d = (int*)malloc (sizeof (int) * 3 * nnF);
	int *facemat_3d = (int*)malloc (sizeof (int) * nnF);

	int *facematerial = (int*)malloc (sizeof (int) * nnF);
	int *tetramaterial = (int*)malloc (sizeof (int) * nnT);

	// Surface mesh generation for parallelepiped with sides 1, 2 and 1.55 and mesh step 0.3
	// scg_make_paral (&nV_in1, &nF_in1, Vertex_in1, Index_in1, 1, 2, 1.55, 0.3);

	// Surface mesh generation for sphere with radius 1.6 and mesh step 0.3
	scg_make_sphere (&nV_in1, &nF_in1, Vertex_in1, Index_in1, 1.6, 0.3);

	// Scaling x coordinates by 0.8, y coordinates by 1.2 and z coordinates by 1 
	scg_scale (nV_in1, Vertex_in1, 0.8, 1.2, 1);

	// Surface mesh generation for cylinder with radius 0.5, height 5 and mesh step 0.3
        scg_make_cylinder (&nV_in2, &nF_in2, Vertex_in2, Index_in2, 0.5, 5, 0.3);

	// Surface mesh generation for sphere with radius 1.8 and mesh step 0.3
	// scg_make_sphere (&nV_in2, &nF_in2, Vertex_in2, Index_in2, 1.8, 0.3);

	// Moving an object along (0.1, 0.5, -2.3) vector
	scg_translate (nV_in2, Vertex_in2, 0.1, 0.5, -2.3);

	// Rotation of an object by angles pitch = 0.3, yaw = 0.1 and roll = 0;
        scg_rotate (nV_in2, Vertex_in2, 0.3, 0.1, 0);
		
	// Intersection of meshec for two objects
	// the last integer parameret is an operation: 
	// 0 - union (in1 + in2)
	// 1 - difference (in1 - in2)
	// 2 - another difference (in2 - in1)
	// 3 - intersection (in1 * in2)
	scg_intersect_mesh (nV_in1, nF_in1, Vertex_in1, Index_in1,
			    nV_in2, nF_in2, Vertex_in2, Index_in2,
			    &nV_out, &nF_out, Vertex_out, Index_out, 2);

	// Fill face material array
	for (i = 0; i < nF_out; i++)
	    facematerial[i] = 1;

	// Write mesh to frt.gmv file
        write_mesh_gmv ("frt.gmv", nV_out, Vertex_out, nF_out, Index_out, facematerial, 0, NULL, NULL);
	write_front("frt.frt", nV_out, Vertex_out, nF_out, Index_out, facematerial);
	
	r = nF_3d = nT_out = 0;

	// We will copy the front, so that it could be used in output in future
	memcpy(Index_3d, Index_out, sizeof(int) * 3 * nnF);
	memcpy(facemat_3d, facematerial, sizeof(int) * nnF);
	nF_3d = nF_out;
	
	// 3D mesh gereration for resulting mesh as initial front
	r = mesh_3d_aft_(&nV_out, Vertex_out, &nF_3d, Index_3d, facemat_3d, 
			 &nT_out, Tetra_out, tetramaterial, &nnV, &nnF, &nnT);
	printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV_out, nF_3d, nT_out);

	if (r == 0) 
	{
		// Checking that 3D mesh corresponds with surface mesh
		check_mesh_topology_(&nV_out, Vertex_out, &nF_out, Index_out, &nT_out, Tetra_out);

		// Write output files
		write_mesh_gmv("mesh.gmv", nV_out, Vertex_out, nF_out, Index_out, facemat_3d, 
				nT_out, Tetra_out, tetramaterial);
		write_mesh    ("mesh.out", nV_out, Vertex_out, nF_out, Index_out, facemat_3d, 
				nT_out, Tetra_out, tetramaterial);
		printf("\nINFO: nV = %d, nF = %d, nT = %d\n", nV_out, nF_out, nT_out);
	} 
	else 
	{
		// Meshing failed. We save current meshed part and current front in gmv format
		printf("Meshing FAILED. Sorry.\n");
		write_mesh_gmv("fail.gmv", nV_out, Vertex_out, nF_3d, Index_3d, facemat_3d, 
				nT_out, Tetra_out, tetramaterial);
	}
	
	free (Vertex_in1);
	free (Index_in1);
	
	free (Vertex_in2);
	free (Index_in2);
	
	free (Vertex_out);
	free (Index_out);
	free (Tetra_out);

	free (Index_3d);
	free (facemat_3d);

	free (facematerial);
	free (tetramaterial);
	
	return 0;
}
