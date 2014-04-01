extern "C"
{
#include "libfrtprm.h"
#include "libfrtscg.h"
}
#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern "C" int scg_intersect_mesh (int nV_in1, int nF_in1, double *Vertex_in1, int *Index_in1,
			int nV_in2, int nF_in2, double *Vertex_in2, int *Index_in2,
			int *pnV_out, int *pnF_out, double *Vertex_out, int *Index_out, int op)
{
	mesh **msh;
	
	msh = (mesh**)malloc (sizeof (mesh*) * 2);

	msh[0] = new mesh ();
	msh[0]->mesh_front (nV_in1, nF_in1, Vertex_in1-3, Index_in1);

	msh[1] = new mesh ();
	msh[1]->mesh_front (nV_in2, nF_in2, Vertex_in2-3, Index_in2);

	msh[0]->in = op % 2;
	msh[1]->in = op / 2;

	msh[0]->intersect (msh[1]);

	*pnV_out = msh[0]->nV;
	*pnF_out = msh[0]->nF;

	memcpy (Vertex_out, msh[0]->Vertex+3, 3 * (msh[0]->nV) * sizeof (double));
	memcpy (Index_out, msh[0]->Index, 3 * msh[0]->nF * sizeof (int));

	free (msh);
	return 0;
}

/*
int scg_intersect_files (char *file_in1, char *file_in2, char *file_out, int op)
{
	int ret;
	mesh **msh;
	
	msh = (mesh**)malloc (sizeof (mesh*) * 2);

	msh[0] = new mesh ();
	if ((ret = msh[0]->mesh_smv (file_in1, 0)) < 0)
	{
		printf ("Error %d while reading file %s\n", ret, file_in1);
		return ret;
	}
	
	msh[1] = new mesh ();
	if ((ret = msh[1]->mesh_smv (file_in2, 0)) < 0)
	{
		printf ("Error %d while reading file %s\n", ret, file_in2);
		return ret;
	}
	
	msh[0]->in = op % 2;
	msh[1]->in = op / 2;

	msh[0]->intersect (msh[1]);
	msh[0]->mesh_write_smv("file_out");
	
	free (msh);
	return 0;
}
*/
