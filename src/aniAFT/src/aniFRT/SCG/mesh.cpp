#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
extern "C" {
#include "libprim.h"
}
#include "raytri.c"
#include "tritri1.c"


// Change vectors u[i] and v[j] their places
#define INVERSE3(u, v) \
	(u)[0] += (v)[0]; \
	(u)[1] += (v)[1]; \
	(u)[2] += (v)[2]; \
\
	(v)[0] = (u)[0] - (v)[0]; \
	(v)[1] = (u)[1] - (v)[1]; \
	(v)[2] = (u)[2] - (v)[2]; \
\
	(u)[0] -= (v)[0]; \
	(u)[1] -= (v)[1]; \
	(u)[2] -= (v)[2];


#define INVERSE2(u, v) \
	(u)[0] += (v)[0]; \
	(u)[1] += (v)[1]; \
\
	(v)[0] = (u)[0] - (v)[0]; \
	(v)[1] = (u)[1] - (v)[1]; \
\
	(u)[0] -= (v)[0]; \
	(u)[1] -= (v)[1]; \


// Scalar multiplication of (u2 - u1, v2 - v1)
#define DOT2(u1, u2, v1, v2) \
	((u2)[0] - (u1)[0]) * ((v2)[0] - (v1)[0]) + \
	((u2)[1] - (u1)[1]) * ((v2)[1] - (v1)[1]) + \
	((u2)[2] - (u1)[2]) * ((v2)[2] - (v1)[2])


// Vector multiplication of (u2 - u1, v2 - v1)
#define CROSS2(dest, v1, v2, u1, u2) \
	(dest)[0] = ((v1)[1] - (v2)[1]) * ((u1)[2] - (u2)[2]) - ((v1)[2] - (v2)[2]) * ((u1)[1] - (u2)[1]); \
	(dest)[1] = ((v1)[2] - (v2)[2]) * ((u1)[0] - (u2)[0]) - ((v1)[0] - (v2)[0]) * ((u1)[2] - (u2)[2]); \
	(dest)[2] = ((v1)[0] - (v2)[0]) * ((u1)[1] - (u2)[1]) - ((v1)[1] - (v2)[1]) * ((u1)[0] - (u2)[0]);

// Copy vector v to u
#define CPVECT(u, v) \
	(u)[0] = (v)[0]; \
	(u)[1] = (v)[1]; \
	(u)[2] = (v)[2];

// Copy vector v to u
#define CPVECT2(u, v, a) \
	(u)[0] = (v)[0] + a; \
	(u)[1] = (v)[1] + a; \
	(u)[2] = (v)[2] + a;

// Compare 2 vectors
#define EQUAL(u, v) \
	fabs ((u)[0] - (v)[0]) < EPSILON && \
	fabs ((u)[1] - (v)[1]) < EPSILON && \
	fabs ((u)[2] - (v)[2]) < EPSILON

// Check if vectors (u2 - u1) and (v2 - v1) are colinear
#define COLINEAR(u1, u2, v1, v2) \
	fabs (((u2)[1] - (u1)[1]) * ((v2)[2] - (v1)[2]) - ((u2)[2] - (u1)[2]) * ((v2)[1] - (v1)[1])) < EPSILON && \
        fabs (((u2)[2] - (u1)[2]) * ((v2)[0] - (v1)[0]) - ((u2)[0] - (u1)[0]) * ((v2)[2] - (v1)[2])) < EPSILON && \
        fabs (((u2)[0] - (u1)[0]) * ((v2)[1] - (v1)[1]) - ((u2)[1] - (u1)[1]) * ((v2)[0] - (v1)[0])) < EPSILON

// Calculates a cosine of a angle between (y, x) and (y, z)
#define COS(x, y, z) \
	(((x)[0] - (y)[0]) * ((z)[0] - (y)[0]) + \
	((x)[1] - (y)[1]) * ((z)[1] - (y)[1]) + \
	((x)[2] - (y)[2]) * ((z)[2] - (y)[2])) / \
	sqrt ((((x)[0] - (y)[0]) * ((x)[0] - (y)[0]) + \
	       ((x)[1] - (y)[1]) * ((x)[1] - (y)[1]) + \
	       ((x)[2] - (y)[2]) * ((x)[2] - (y)[2])) * \
	      (((z)[0] - (y)[0]) * ((z)[0] - (y)[0]) + \
	       ((z)[1] - (y)[1]) * ((z)[1] - (y)[1]) + \
	       ((z)[2] - (y)[2]) * ((z)[2] - (y)[2])))

// Distance between points # Index[3 * i + k] and Index[3 * i + l]
#define LENGTH(i, k, l)	sqrt ( \
	(Vertex[3 * Index[3 * i + k]] - Vertex[3 * Index[3 * i + l]]) * \
	(Vertex[3 * Index[3 * i + k]] - Vertex[3 * Index[3 * i + l]]) + \
	(Vertex[3 * Index[3 * i + k] + 1] - Vertex[3 * Index[3 * i + l] + 1]) * \
	(Vertex[3 * Index[3 * i + k] + 1] - Vertex[3 * Index[3 * i + l] + 1]) + \
	(Vertex[3 * Index[3 * i + k] + 2] - Vertex[3 * Index[3 * i + l] + 2]) * \
	(Vertex[3 * Index[3 * i + k] + 2] - Vertex[3 * Index[3 * i + l] + 2]))

// Distance between points x and y
#define DIST(x, y) sqrt ( \
	((x)[0] - (y)[0]) * ((x)[0] - (y)[0]) + \
	((x)[1] - (y)[1]) * ((x)[1] - (y)[1]) + \
	((x)[2] - (y)[2]) * ((x)[2] - (y)[2]))

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

mesh::mesh ()
{
    x = y = z = 0;
    pitch = yaw = roll = 0;

    nV = 0;
    nF = 0;

    hide = false;
    in = 1;
}

mesh::mesh (int nv, int nf)
{
    x = y = z = 0;
    pitch = yaw = roll = 0;
    nV = nv;
    nF = nf;

    hide = false;
    in = 1;

    Vertex = (double*)malloc (sizeof (double) * 3 * nV);
    Index = (int*)malloc (sizeof (int) * 3 * nF);
}
    
mesh::~mesh ()
{
    if (nV > 0)
        free (Vertex);
    if (nF > 0)
        free (Index);
}

int mesh::mesh_write_smv (char *name)
{   
    FILE *f;
    int i, c = 0;

    if (!(f = fopen (name, "w")))
	return -1;

	 fprintf(f, "%d %d\n", nV, nF);

	 for (i = 1; i <= nV; i++) 
	 {
		 fprintf(f, "%20.15e %20.15e %20.15e\n", Vertex[3 * i + 0], Vertex[3 * i + 1], Vertex[3 * i + 2]);
	 }
	 for (i = 0; i < nF; i++) 
	 {
		 fprintf(f, "%d %d %d  %d\n", Index[3 * i + 0], Index[3 * i + 1], Index[3 * i + 2], c);
//		 c = (int)(Color[i]*10) - 6;
//		 fprintf(f, "%d %d %d  %d\n", Index[3*i+((c<0)?0:2)], Index[3*i+1], Index[3*i+((c<0)?2:0)], c);
//		 fprintf(f, "%d %d %d  %d\n", Index[3*i+((c<7)?0:2)], Index[3*i+1], Index[3*i+((c<7)?2:0)], c);
	 }
    
    fclose (f);
    printf("%s\n", name);
    return 0;
}

int mesh::mesh_smv (char *name, int color)
{   
    FILE *f;
    int i;

    double glf;
    int u1, u2, u3;
    
    x = y = z = 0;
    pitch = yaw = roll = 0;

    hide = false;
    in = 1;

    if (!(f = fopen (name, "r")))
	return -1;

    if (!fscanf (f, "%d", &nV))
	return -2;

    if (!fscanf (f, "%d", &nF))
	return -3;

    if (nV < 0 || nF < 0)
	return -4;

    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1));
    Index = (int*)malloc (sizeof (int) * 3 * nF);    

    Vertex[0] = Vertex[1] = Vertex[2] = 0.0f; // In case that numeration goes from 1, instead of 0;

    for (i = 0; i < 3 * nV; i++)
    {
	if (!fscanf (f, "%lf", &glf))
            return -5;
	
	Vertex[i + 3] = (double)glf;	
    }
    
    if (color == 1)
    {
        for (i = 0; i < nF; i++)
        {
	    if (fscanf (f, "%u %u %u %lf", &u1, &u2, &u3, &glf) < 4)
                return -6;	

	    if (u1 == 0 || u2 == 0 || u3 == 0)
                printf ("Warning, there's down from 0 numeration!\n");
	
	    Index[i * 3] = (int)u1;
	    Index[i * 3 + 1] = (int)u2;
	    Index[i * 3 + 2] = (int)u3;
        }
    }

    else
    {
        for (i = 0; i < nF; i++)
        {
	    if (fscanf (f, "%u %u %u", &u1, &u2, &u3) < 3)
                return -6;	

	    if (u1 == 0 || u2 == 0 || u3 == 0)
                printf ("Warning, there's down from 0 numeration!\n");
	
	    Index[i * 3] = (int)u1;
	    Index[i * 3 + 1] = (int)u2;
	    Index[i * 3 + 2] = (int)u3;
        }
    }

    printf ("nV = %d, nF = %d\n", nV, nF);
    
    fclose (f);
    return 0;
}


int mesh::mesh_front (int nV_, int nF_, double *Vertex_, int *Index_)
{
    nV = nV_;
    nF = nF_;

    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1));
    Index = (int*)malloc (sizeof (int) * 3 * nF);

    memcpy (Vertex, Vertex_, sizeof (double) * 3 * (nV + 1));
    memcpy (Index, Index_, sizeof (int) * 3 * nF);

    return 0;
}


int mesh::primitive_paral (GLdouble x, GLdouble y, GLdouble z, double MS)
{
    double *Vertex_ = (double*)malloc (sizeof (double) * 3 * 1000000);
    int *Index_ = (int*)malloc (sizeof (int) * 3 * 1000000);
    scg_make_paral (&nV, &nF, Vertex_, Index_, x, y, z, MS);
    
    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1));
    Index = (int*)malloc (sizeof (int) * 3 * nF);

    memcpy (Vertex+3, Vertex_, sizeof (double) * 3 * (nV));
    memcpy (Index, Index_, sizeof (int) * 3 * nF);

    free (Vertex_);
    free (Index_);
    return 0;
}


int mesh::primitive_sphere (GLdouble x, double MS)
{
    double *Vertex_ = (double*)malloc (sizeof (double) * 3 * 1000000);
    int *Index_ = (int*)malloc (sizeof (int) * 3 * 1000000);
    scg_make_sphere (&nV, &nF, Vertex_, Index_, x, MS);
    
    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1));
    Index = (int*)malloc (sizeof (int) * 3 * nF);

    memcpy (Vertex+3, Vertex_, sizeof (double) * 3 * (nV));
    memcpy (Index, Index_, sizeof (int) * 3 * nF);

    free (Vertex_);
    free (Index_);

    return 0;
}


int mesh::primitive_cylinder (GLdouble r, GLdouble h, double MS)
{
    double *Vertex_ = (double*)malloc (sizeof (double) * 3 * 1000000);
    int *Index_ = (int*)malloc (sizeof (int) * 3 * 1000000);
    scg_make_cylinder (&nV, &nF, Vertex_, Index_, r, h, MS);
    
    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1));
    Index = (int*)malloc (sizeof (int) * 3 * nF);

    memcpy (Vertex+3, Vertex_, sizeof (double) * 3 * (nV));
    memcpy (Index, Index_, sizeof (int) * 3 * nF);

    free (Vertex_);
    free (Index_);

    return 0;
}


// May be rewritten into parallel...
int mesh::move (double x, double y, double z, double rot_x, double rot_y, double rot_z)
{
    scg_translate (nV, Vertex+3, x, y, z);
    
    rot_x = rot_x - pitch;
    rot_y = rot_y - yaw;
    rot_z = rot_z - roll; 

    scg_rotate (nV, Vertex+3, rot_x, rot_y, rot_z);
 
    pitch += rot_x;
    yaw += rot_y;
    roll += rot_z;
    
    return 0;
}


int remove_unused_vertices (int *pnV, double *vertex, int *pnF, int *face) {
    int nV = *pnV + 1, nF = *pnF;
    int *index, *used;
    int i, j, removed = 0;
    index = (int*)malloc (sizeof(int) * nV);
    used  = (int*)malloc (sizeof(int) * nV);
    for (i = 0; i < nV; i++) 
    {
	used[i] = 0;
	index[i] = i;
    }

    for (j = 0; j < 3 * nF; j++)
    {
	used[face[j]]++;
    }

    for (i = 0; i < nV; i++) 
    {
	if (!used[i]) 
        {
	    removed++;
	    nV--;
	    vertex[3 * i + 0] = vertex[3 * nV + 0];
	    vertex[3 * i + 1] = vertex[3 * nV + 1];
	    vertex[3 * i + 2] = vertex[3 * nV + 2];
	    used[i] = used[nV];
	    index[nV] = i;
	    i--;
	}
    }

    for (j = 0; j < 3 * nF; j++)
    {
        face[j] = index[face[j]];
    }

    for (j = 0; j < nF; j++)
    {
	if (face[3 * j] == 0 && face[3 * j + 1] == 0 && face[3 * j + 2] == 0)
	{
	    face[3 * j] = face[3 * (nF - 1)];
	    face[3 * j + 1] = face[3 * (nF - 1) + 1];
	    face[3 * j + 2] = face[3 * (nF - 1) + 2];
	    nF--;
	    j--;
	}
    }

    free(index);
    free(used);
    *pnV = nV - 1;
    *pnF = nF;
    return removed;
}

double Meshsize = 0.0;


int unify_curve (double *vertex, int l, int *loops, double k2)
{
    int k, i, start, nV = loops[l - 1];
    int *unused, *loops_uns;
    double len;
    
    if (k2 <= 0)
	return -1;
    
    unused = (int*)malloc (sizeof(int) * nV);
    loops_uns = (int*)malloc (sizeof(int) * l);
    for (i = 0; i < nV; i++) 
    {
	unused[i] = 0;
    }

    for (k = 0; k < l - 1; k++)
    {
	len = DIST (vertex + (loops[k + 1] - 1) * 3, vertex + loops[k] * 3);
	for (i = loops[k], start = loops[k], loops_uns[k] = 0; i < loops[k + 1] - 2; i++)
	{
	    if (DIST (vertex + i * 3, vertex + (i + 1) * 3) < Meshsize / k2)
	    {
		loops_uns[k]++;
		if (len < Meshsize / k2 || len < DIST (vertex + (i + 1) * 3, vertex + (i + 2) * 3))
		{
		    if (unused[i] == 1)
			loops_uns[k]--;
		    unused[i] = 1;
		    len += DIST (vertex + i * 3, vertex + (i + 1) * 3);
//		    printf ("1) point %d out. len = %f\n", i, len);
		}
		else 
		{
		    if (unused[i + 1] == 1)
			loops_uns[k]--;
		    start = i;
		    unused[i + 1] = 1;
		    len = DIST (vertex + i * 3, vertex + (i + 1) * 3);		    
//		    printf ("2) point %d out. len = %f\n", i + 1, len);
		}
	    }
	    else
		len = DIST (vertex + i * 3, vertex + (i + 1) * 3);
	}
	// Here i = loops[k + 1] - 2
	if (DIST (vertex + i * 3, vertex + (i + 1) * 3) < Meshsize / k2)
	{
	    loops_uns[k]++;
	    if ((len < Meshsize / k2 || len < DIST (vertex + (i + 1) * 3, vertex + (i + 2) * 3)) || 
		unused[loops[k]] == 1)
	    {
		if (unused[i] == 1)
		    loops_uns[k]--;
		unused[i] = 1;
		len += DIST (vertex + i * 3, vertex + (i + 1) * 3);
	    }
	    else 
	    {
		if (unused[i + 1] == 1)
		    loops_uns[k]--;
	        start = i;
	        unused[i + 1] = 1;
	        len = DIST (vertex + i * 3, vertex + (i + 1) * 3);		    
	    }
	}
	else
	    len = DIST (vertex + i * 3, vertex + (i + 1) * 3);

    }
    
    for (i = 0, start = 0; i < nV; i++)
    {
	if (unused[i] == 0)
	{
	    CPVECT (vertex + start * 3, vertex + i * 3);
	    start++;
	}
    }
    
    for (k = 0; k < l - 1; k++)
    {
	for (i = k; i < l - 1; i++)
	{
	    loops[i + 1] -= loops_uns[k];
	}
    }

    return nV - start;
}


int mesh::intersect (mesh *msh)
{
    int i, j, j2, k, l, m, nFtmp = 0, nFnew = 0, nEbound = 0, nEbound2 = 0, ret, start, current;
    double len;

    double k1 = 5, k2 = 3;
    
    if (nV < 3 || msh->nV < 3 || nF < 4 || msh->nF < 4)
    {
        printf ("Wrong or empty meshes!\n");
	return -10;
    }

//  Temporary copy of points from both meshes and indeces
    double *tmpVtx = (double*)malloc (sizeof (double) * 3 * (nV + msh->nV + 1));    
    int *tmpIdx = (int*)malloc (sizeof (int) * 3 * (nF + msh->nF));

//  An array containing information about position of a point inside or outside of the second object
    int *inside = (int*)malloc (sizeof (int) * (nV + msh->nV + 1));

//  New points added when we calculate intersections
    double *newVtx = (double*)malloc (sizeof (double) * 6 * 10 * (nF + msh->nF)); // FIXME!
    
    double *trimVtx;

    int *trimLoops;
    int *boundLoops1;
    int *boundLoops2;

    int *trimBound1;
    int *trimBound2;

    int *trimIdx;
	 
    int *used = (int*)malloc (sizeof (int) * nF * msh->nF);
    for (i = 0; i < nF * msh->nF; i++)
    {
	used[i] = -1;
    }

//  Special array of parameters for each edge of a trimming curve
//  which contains numbers of points of triangles 
    int *trimParam = (int*)malloc (sizeof (int) * 6 * nF * msh->nF);
    
    double isectpt1[3];
    double isectpt2[3];

    memcpy (tmpVtx, Vertex, 3 * (nV + 1) * sizeof (double));
    memcpy (tmpVtx + 3 * (nV + 1), msh->Vertex + 3, 3 * msh->nV * sizeof (double));

    if (in == 1 && msh->in == 0)
    {
        for (i = 0; i < nF; i++)
	{
	    j = Index[3 * i + 1];
	    Index[3 * i + 1] = Index[3 * i + 2];
	    Index[3 * i + 2] = j;
	}
    }
    
    if (in == 0 && msh->in == 1)
    {
        for (i = 0; i < msh->nF; i++)
	{
	    j = msh->Index[3 * i + 1];
	    msh->Index[3 * i + 1] = msh->Index[3 * i + 2];
	    msh->Index[3 * i + 2] = j;
	}
    }    

    for (i = 1; i <= nV; i++)
        inside[i] = msh->isInside (Vertex[3 * i], Vertex[3 * i + 1], Vertex[3 * i + 2]);
    
    for (i = 1; i <= msh->nV; i++)
        inside[i + nV] = isInside (msh->Vertex[3 * i], msh->Vertex[3 * i + 1], msh->Vertex[3 * i + 2]);

    for (i = 0; i < nF; i++)
    {
	len = MAX (DIST (Vertex + 3 * Index[3 * i], Vertex + 3 * Index[3 * i + 1]),
		   MAX (DIST (Vertex + 3 * Index[3 * i + 1], Vertex + 3 * Index[3 * i + 2]),
			DIST (Vertex + 3 * Index[3 * i + 2], Vertex + 3 * Index[3 * i])));
	if (len > Meshsize)
	    Meshsize = len;
    }
    
    for (i = 0; i < msh->nF; i++)
    {
	len = MAX (DIST (msh->Vertex + 3 * msh->Index[3 * i], msh->Vertex + 3 * msh->Index[3 * i + 1]),
		   MAX (DIST (msh->Vertex + 3 * msh->Index[3 * i + 1], msh->Vertex + 3 * msh->Index[3 * i + 2]),
			DIST (msh->Vertex + 3 * msh->Index[3 * i + 2], msh->Vertex + 3 * msh->Index[3 * i])));
	if (len > Meshsize)
	    Meshsize = len;
    }

    for (i = 0; i < nF; i++)
    {
        ret = 0;

	if (!msh->triInt (this, i, &ret, isectpt1, isectpt2))
	{
	    if (inside[Index[3 * i]] == in)
	    {
            CPVECT (tmpIdx + 3 * nFtmp, Index + 3 * i);
            nFtmp++;
	    }
	}
	
	else
	{
            if (nFnew >= nF * msh->nF)
	    {
		printf ("Oops!\n");
		return -1;
	    }
	    
            CPVECT (newVtx + 6 * nFnew, isectpt1);
            CPVECT (newVtx + 6 * nFnew + 3, isectpt2);

	    CPVECT (trimParam + 6 * nFnew, Index + 3 * i);
	    CPVECT2 (trimParam + 6 * nFnew + 3, msh->Index + 3 * (ret - 1), nV);

	    nFnew++;
	    
            while (msh->triInt (this, i, &ret, isectpt1, isectpt2))
	    {
                if (nFnew >= nF * msh->nF)
	        {
		    printf ("Ooooops!\n");
		    return -1;
	        }
		
                CPVECT (newVtx + 6 * nFnew, isectpt1);
                CPVECT (newVtx + 6 * nFnew + 3, isectpt2);

	        CPVECT (trimParam + 6 * nFnew, Index + 3 * i);
	        CPVECT2 (trimParam + 6 * nFnew + 3, msh->Index + 3 * (ret - 1), nV);
		
		nFnew++;
	    }
	}
    }     

    if (nFnew < 3)
    {
      printf ("The intersection is empty! :(\n");
        return -9;
    }

    trimVtx = (double*)malloc (sizeof (double) * 3 * nFnew);
    trimLoops = (int*)malloc (sizeof (int) * (int)(nFnew / 3 + 1));

    trimBound1 = (int*)malloc (sizeof (int) * 4 * nFnew);
    trimBound2 = (int*)malloc (sizeof (int) * 4 * nFnew);    

    CPVECT (trimVtx, newVtx);
    CPVECT (trimVtx + 3, newVtx + 3);

    trimLoops[0] = 0;
    l = 1;

    used[0] = 0;

    start = 0;

    for (i = 2; i < nFnew; i++)
    {
        for (j = start + 2; j < 2 * nFnew; j++)
	{
	    if (used[j / 2] < 0 && EQUAL (newVtx + j * 3, trimVtx + (i - 1) * 3))
	    {
                j2 = j - ((j % 2) * 2 - 1); //  2 * n      -> (2 * n + 1)
		                            // (2 * n + 1) ->  2 * n
		
		CPVECT (trimVtx + i * 3, newVtx + j2 * 3);
//		printf ("j / 2 = %d, i = %d, ", j / 2, i);
//		printPoint (trimVtx + i * 3);

		used[j / 2] *= -1;
		used[i - 1] *= j / 2;

//		printf ("used[%d] = %d   -   ", j / 2, used[j / 2]);
//		printf ("used[%d] = %d\n", i - 1, used[i - 1]);

		if (EQUAL (trimVtx + i * 3, trimVtx + start * 3))
	        {
                    for (j2 = start + 2; j2 < 2 * nFnew; j2++)
		    {
                        if (used[j2 / 2] < 0 && !(EQUAL (newVtx + j2 * 3, newVtx + (j2 + 1) * 3)))
		        {			    
                            start = j2;
			    
		            CPVECT (trimVtx + i * 3, newVtx + start * 3);
		            CPVECT (trimVtx + (i + 1) * 3, newVtx + (start + 1) * 3);
			    
			    trimLoops[l] = i;
			    l++;
		    
			    used[start / 2] *= -1;
			    used[i] *= start / 2;
			    
			    i++;
			    break;
		        }
		    }
		}		
		break;
	    }
	}
	if (j == 2 * nFnew)
	    break;
    }
    trimLoops[l] = i;
    l++;

    printf ("l = %d\n", l);

    for (i = 0; i < msh->nF; i++)
    {
        ret = 0;
	if (!triInt (msh, i, &ret, isectpt1, isectpt2))
	{
	    if (inside[nV + msh->Index[3 * i]] == msh->in)
	    {
	        CPVECT2 (tmpIdx + 3 * nFtmp, msh->Index + 3 * i, nV);
	        // tmpIdx[3 * nFtmp] = msh->Index[3 * i];
	        // tmpIdx[3 * nFtmp + 1] = msh->Index[3 * i + 2];
	        // tmpIdx[3 * nFtmp + 2] = msh->Index[3 * i + 1];
	        nFtmp++;
	    }
	}
    }

    boundLoops1 = (int*)malloc (sizeof (int) * l);
    boundLoops2 = (int*)malloc (sizeof (int) * l);

    boundLoops1[0] = 0;
    boundLoops2[0] = 0;

    for (k = 0; k < l - 1; k++)
    {
	for (i = trimLoops[k], j2 = 0; i < trimLoops[k + 1] && j2 < 2; i++)
	{
	    if (used[i] < 0)
		used[i] = -used[i];

//	    printf ("\nused[%d] = %d   ", i, used[i]);
            for (j = 0; j < 3; j++)
            {
//		printf ("tP = %d   ", trimParam[6 * used[i] + j]);
                if (inside[trimParam[6 * used[i] + j]] == in && (j2 == 0 || 
		    trimParam[6 * used[i] + j] != trimBound1[nEbound - j2]))
                {
	            trimBound1[nEbound] = trimParam[6 * used[i] + j];
	            nEbound++;
		    j2++;
		}
	    }
//	    printf ("\n");
	}

	if (j2 < 2)
	{
	    printf ("Wrong intersection!\n");
	    return -8;
	}

	if (j2 == 2)
	{
	    nEbound--;
	    if (DOT2 (tmpVtx + 3 * trimBound1[nEbound - 1], tmpVtx + 3 * trimBound1[nEbound], 
		      trimVtx + 3 * trimLoops[k], trimVtx + 3 * (trimLoops[k] + 1)) < 0)
	    {
                trimBound1[nEbound - 1] += trimBound1[nEbound];
                trimBound1[nEbound] = trimBound1[nEbound - 1] - trimBound1[nEbound];
                trimBound1[nEbound - 1] -= trimBound1[nEbound];
	    }
	}

	if (j2 == 3) // This case must not occur!!
	{
	    printf ("Achtung! 3!\n");
	    nEbound -= 2;
	    if (DOT2 (tmpVtx + 3 * trimBound1[nEbound - 1], tmpVtx + 3 * trimBound1[nEbound], 
		      trimVtx + 3 * trimLoops[k], trimVtx + 3 * (trimLoops[k] + 1)) < 0)
	    {
                trimBound1[nEbound - 1] += trimBound1[nEbound];
                trimBound1[nEbound] = trimBound1[nEbound - 1] - trimBound1[nEbound];
                trimBound1[nEbound - 1] -= trimBound1[nEbound];
	    }
	}	

	start = nEbound - 1;
	m = 0;
	trimBound1[nEbound + 1] = trimBound1[nEbound - 1];
	j2 = -1;
	while (j2 != -2 && nEbound < 100)//nFtmp)
	{
	    for (j = 0; j < nFtmp; j++)
	    {
	        if (j2 != j &&
		    ((trimBound1[nEbound] == tmpIdx[3 * j + 0] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 1]) ||
		     (trimBound1[nEbound] == tmpIdx[3 * j + 0] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 2]) ||
		     (trimBound1[nEbound] == tmpIdx[3 * j + 1] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 0]) ||
		     (trimBound1[nEbound] == tmpIdx[3 * j + 1] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 2]) ||
		     (trimBound1[nEbound] == tmpIdx[3 * j + 2] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 0]) ||
		     (trimBound1[nEbound] == tmpIdx[3 * j + 2] && trimBound1[nEbound + 1] == tmpIdx[3 * j + 1])))
		{
//		    printf ("(%d, %d), (%d, %d, %d) -> ", trimBound1[nEbound], trimBound1[nEbound + 1], 
//		    	                        tmpIdx[3 * j + 0], tmpIdx[3 * j + 1], tmpIdx[3 * j + 2]);
		    current = trimBound1[nEbound + 1];
		    trimBound1[nEbound + 1] = tmpIdx[3 * j + 0] + tmpIdx[3 * j + 1] + tmpIdx[3 * j + 2] - 
					      trimBound1[nEbound] - trimBound1[nEbound + 1];

//		    printf ("%d\n", trimBound1[nEbound + 1]);	    
		    j2 = j;

		    if ((dist_point_curve (tmpVtx + 3 * trimBound1[nEbound + 1], newVtx, nFnew, &ret) < Meshsize / k1) 
			&& nEbound > start) 
		    {
			for (i = 0; i < nFtmp; i++)
			{
			    if (trimBound1[nEbound + 1] == tmpIdx[3 * i + 0] ||
				trimBound1[nEbound + 1] == tmpIdx[3 * i + 1] ||
				trimBound1[nEbound + 1] == tmpIdx[3 * i + 2])
			    {
			        tmpIdx[3 * i + 0] = tmpIdx[3 * i + 1] = tmpIdx[3 * i + 2] = 0;
			    }
			}
			nEbound++;
			trimBound1[nEbound] = current;
			trimBound1[nEbound + 1] = trimBound1[nEbound - 1];
			j2 = -1;
		    }
		    break;
		}
	    }

	    if (j == nFtmp && j2 != -1)
	    {
	        j2 = -1;
		nEbound++;
		trimBound1[nEbound + 1] = trimBound1[nEbound - 1];
	    }

	    if (trimBound1[nEbound] == trimBound1[start + 2] && nEbound > boundLoops1[k] + 2)
	    {
		trimBound1[start] = trimBound1[nEbound - 2];
		trimBound1[start + 1] = trimBound1[nEbound - 1];
		    
		boundLoops1[k + 1] = nEbound - 1;
		nEbound--;
		j2 = -2;
		break;		
	    }
	}
    }

    /////////////////////////////////////////////////////////////////////////////////////

    for (k = 0; k < l - 1; k++)
    {
	for (i = trimLoops[k], j2 = 0; i < trimLoops[k + 1] && j2 < 2; i++)
	{
	    if (used[i] < 0)
		used[i] = -used[i];
	    
//	    printf ("used[%d] = %d   ", i, used[i]);
	    for (j = 0; j < 3; j++)
            {
//		printf ("tP = %d   ", trimParam[6 * used[i] + j + 3]);
                if (inside[trimParam[6 * used[i] + j + 3]] == msh->in && (j2 == 0 || 
		    trimParam[6 * used[i] + j + 3] != trimBound2[nEbound2 - j2]))
                {
	            trimBound2[nEbound2] = trimParam[6 * used[i] + j + 3];
	            nEbound2++;
		    j2++;
		}
	    }
//	    printf ("\n");
	}
	
	if (j2 < 2)
	{
	    printf ("Wrong intersection!\n");
	    return -8;
	}

	if (j2 == 2)
	{
	    nEbound2--;
	    if (DOT2 (tmpVtx + 3 * trimBound2[nEbound2 - 1], tmpVtx + 3 * trimBound2[nEbound2], 
		      trimVtx + 3 * trimLoops[k], trimVtx + 3 * (trimLoops[k] + 1)) < 0)
	    {
                trimBound2[nEbound2 - 1] += trimBound2[nEbound2];
                trimBound2[nEbound2] = trimBound2[nEbound2 - 1] - trimBound2[nEbound2];
                trimBound2[nEbound2 - 1] -= trimBound2[nEbound2];
	    }
	}

	if (j2 == 3) // This case must not occur!!
	{
	    printf ("Achtung! 3!\n");
	    nEbound2 -= 2;
	    if (DOT2 (tmpVtx + 3 * trimBound2[nEbound2 - 1], tmpVtx + 3 * trimBound2[nEbound2], 
		      trimVtx + 3 * trimLoops[k], trimVtx + 3 * (trimLoops[k] + 1)) < 0)
	    {
                trimBound2[nEbound2 - 1] += trimBound2[nEbound2];
                trimBound2[nEbound2] = trimBound2[nEbound2 - 1] - trimBound2[nEbound2];
                trimBound2[nEbound2 - 1] -= trimBound2[nEbound2];
	    }
	}	

	start = nEbound2 - 1;
	trimBound2[nEbound2 + 1] = trimBound2[nEbound2 - 1];
	j2 = -1;
	while (j2 != -2 && nEbound2 < nFtmp)
	{
	    for (j = 0; j < nFtmp; j++)
	    {
	        if (j2 != j &&
		    ((trimBound2[nEbound2] == tmpIdx[3 * j + 0] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 1]) ||
		     (trimBound2[nEbound2] == tmpIdx[3 * j + 0] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 2]) ||
		     (trimBound2[nEbound2] == tmpIdx[3 * j + 1] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 0]) ||
		     (trimBound2[nEbound2] == tmpIdx[3 * j + 1] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 2]) ||
		     (trimBound2[nEbound2] == tmpIdx[3 * j + 2] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 0]) ||
		     (trimBound2[nEbound2] == tmpIdx[3 * j + 2] && trimBound2[nEbound2 + 1] == tmpIdx[3 * j + 1])))
		{
//		    printf ("(%d, %d), (%d, %d, %d) -> ", trimBound2[nEbound2], trimBound2[nEbound2 + 1], 
//		    	                        tmpIdx[3 * j + 0], tmpIdx[3 * j + 1], tmpIdx[3 * j + 2]);
		    current = trimBound2[nEbound2 + 1];
		    trimBound2[nEbound2 + 1] = tmpIdx[3 * j + 0] + tmpIdx[3 * j + 1] + tmpIdx[3 * j + 2] - 
					      trimBound2[nEbound2] - trimBound2[nEbound2 + 1];

//		    printf ("%d\n", trimBound2[nEbound2 + 1]);	    
		    j2 = j;

		    if ((dist_point_curve (tmpVtx + 3 * trimBound2[nEbound2 + 1], newVtx, nFnew, &ret) < Meshsize / k1)  
			&& nEbound2 > start) 
		    {
			for (i = 0; i < nFtmp; i++)
			{
			    if (trimBound2[nEbound2 + 1] == tmpIdx[3 * i + 0] ||
				trimBound2[nEbound2 + 1] == tmpIdx[3 * i + 1] ||
				trimBound2[nEbound2 + 1] == tmpIdx[3 * i + 2])
			    {
			        tmpIdx[3 * i + 0] = tmpIdx[3 * i + 1] = tmpIdx[3 * i + 2] = 0;
			    }
			}
			nEbound2++;
			trimBound2[nEbound2] = current;
			trimBound2[nEbound2 + 1] = trimBound2[nEbound2 - 1];
			j2 = -1; 
		    }
		    break;
		}
	    }

	    if (j == nFtmp && j2 != -1)
	    {
	        j2 = -1;
		nEbound2++;
		trimBound2[nEbound2 + 1] = trimBound2[nEbound2 - 1];
	    }

	    if (trimBound2[nEbound2] == trimBound2[start + 2] && nEbound2 > boundLoops2[k] + 2)
	    {
		trimBound2[start] = trimBound2[nEbound2 - 2];
		trimBound2[start + 1] = trimBound2[nEbound2 - 1];
		    
		boundLoops2[k + 1] = nEbound2 - 1;
		nEbound2--;
		j2 = -2;
		break;		
	    }
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    i = unify_curve (trimVtx, l, trimLoops, k2);
    i += unify_curve (trimVtx, l, trimLoops, k2);
    printf ("%d points removed.\n", i);
    
    nFnew = trimLoops[l - 1];
	    
    trimIdx = (int*)malloc (sizeof (int) * 3 * (2 * nFnew + nEbound + nEbound2 - 2 * (l - 1)));

    for (k = 0, i = 0; k < l - 1; k++)
    {
//	printf ("trimLoops[k + 1] = %d, boundLoops1[k + 1] = %d, boundLoops2[k + 1] = %d\n\n",
//		trimLoops[k + 1], boundLoops1[k + 1], boundLoops2[k + 1]);
	for (j = trimLoops[k], j2 = boundLoops1[k]; 
	     i < trimLoops[k] + boundLoops2[k] + trimLoops[k + 1] + boundLoops1[k + 1] - 2 - 2 * k; i++)
	{
	    if (k % 2 == 0)
	    {
	        trimIdx[3 * i] = trimBound1[j2];
	        trimIdx[3 * i + 1] = nV + msh->nV + 1 + j;
	    }
	    else
	    {
	        trimIdx[3 * i] = nV + msh->nV + 1 + j;
	        trimIdx[3 * i + 1] = trimBound1[j2];
	    }
	    if (j2 >= boundLoops1[k + 1] - 1 || (j < trimLoops[k + 1] - 1 &&
		COS (tmpVtx + 3 * trimBound1[j2], trimVtx + 3 * j, trimVtx + 3 * (j + 1)) > 
		COS (trimVtx + 3 * j, tmpVtx + 3 * trimBound1[j2], tmpVtx + 3 * trimBound1[j2 + 1])))
	    {
		trimIdx[3 * i + 2] = nV + msh->nV + 1 + j + 1;
		j++;
	    }
	    
	    else
	    {
		trimIdx[3 * i + 2] = trimBound1[j2 + 1];
		j2++;
	    }
//	    printf("1) i = %d, j = %d, j2 = %d\n", i, j, j2);
	}

	if (k % 2 == 0)
	{
	    trimIdx[3 * i] = trimBound1[j2];
	    trimIdx[3 * i + 1] = nV + msh->nV + 1 + j;
	}
	else
	{
	    trimIdx[3 * i] = nV + msh->nV + 1 + j;
	    trimIdx[3 * i + 1] = trimBound1[j2];
	}

	trimIdx[3 * i + 2] = nV + msh->nV + 1 + trimLoops[k];
	i++;

	for (j = trimLoops[k], j2 = boundLoops2[k]; 
	     i < 2 * trimLoops[k + 1] + boundLoops1[k + 1] + boundLoops2[k + 1] - 3 - 2 * k; i++)
	{
	    if (k % 2 == 0)
	    {  
	        trimIdx[3 * i] = nV + msh->nV + 1 + j;
	        trimIdx[3 * i + 1] = trimBound2[j2];
	    }
	    else 
	    {
		trimIdx[3 * i] = trimBound2[j2];
		trimIdx[3 * i + 1] = nV + msh->nV + 1 + j;
	    }
	    if (j2 >= boundLoops2[k + 1] - 1 || (j < trimLoops[k + 1] - 1 &&
		COS (tmpVtx + 3 * trimBound2[j2], trimVtx + 3 * j, trimVtx + 3 * (j + 1)) > 
		COS (trimVtx + 3 * j, tmpVtx + 3 * trimBound2[j2], tmpVtx + 3 * trimBound2[j2 + 1])))
	    {
		trimIdx[3 * i + 2] = nV + msh->nV + 1 + j + 1;
		j++;
	    }
	    
	    else
	    {
		trimIdx[3 * i + 2] = trimBound2[j2 + 1];
		j2++;
	    }
//	    printf("2) i = %d, j = %d, j2 = %d\n", i, j, j2);
	}

	if (k % 2 == 0)
	{
  	    trimIdx[3 * i] = nV + msh->nV + 1 + j;
	    trimIdx[3 * i + 1] = trimBound2[j2];
	}
	else
	{
	    trimIdx[3 * i] = trimBound2[j2];
	    trimIdx[3 * i + 1] = nV + msh->nV + 1 + j;
	}

	trimIdx[3 * i + 2] = nV + msh->nV + 1 + trimLoops[k];
	i++;
    }
   
    msh->hide = 1;

    free (Vertex);
    free (Index);

    nV = nV + msh->nV;
    nF = nFtmp;
    
    Vertex = (double*)malloc (sizeof (double) * 3 * (nV + 1 + nFnew));
    Index = (int*)malloc (sizeof (int) * 3 * (nF + 2 * nFnew + nEbound + nEbound2 - 2 * (l - 1)));

    memcpy (Vertex, tmpVtx, 3 * (nV + 1) * sizeof (double));
    memcpy (Vertex + 3 * (nV + 1), trimVtx, 3 * nFnew * sizeof (double));
    memcpy (Index, tmpIdx, 3 * nF * sizeof (int));
    memcpy (Index + 3 * nF, trimIdx, 3 * (2 * nFnew + nEbound + nEbound2 - 2 * (l - 1)) * sizeof (int));

    for (j = 0; j < 3; j++)
    {
	for (k = nF; k < nF + 2 * nFnew + nEbound + nEbound2 - 1 - 2 * (l - 1); k++)
	{
	    if (Index[k * 3] == Index[(k + 1) * 3] && Index[k * 3 + 2] == Index[(k + 1) * 3 + 1])
	    {
		if (COS (Vertex + 3 * Index[(k + 1) * 3], 
		         Vertex + 3 * Index[(k + 1) * 3 + 2], 
		         Vertex + 3 * Index[(k + 1) * 3 + 1]) *
		    (1 - sqrt (COS (Vertex + 3 * Index[k * 3],
				    Vertex + 3 * Index[k * 3 + 1],
				    Vertex + 3 * Index[k * 3 + 2]) *
			       COS (Vertex + 3 * Index[k * 3],
				    Vertex + 3 * Index[k * 3 + 1],
				    Vertex + 3 * Index[k * 3 + 2]))) +
		    COS (Vertex + 3 * Index[k * 3],
		         Vertex + 3 * Index[k * 3 + 1],
		         Vertex + 3 * Index[k * 3 + 2]) *
		    (1 - sqrt (COS (Vertex + 3 * Index[(k + 1) * 3], 
				    Vertex + 3 * Index[(k + 1) * 3 + 2],
				    Vertex + 3 * Index[(k + 1) * 3 + 1]) *
			       COS (Vertex + 3 * Index[(k + 1) * 3], 
				    Vertex + 3 * Index[(k + 1) * 3 + 2], 
				    Vertex + 3 * Index[(k + 1) * 3 + 1]))) < 0)
	        {
		    i = Index[k* 3 + 2];
		    Index[k * 3 + 2] = Index[(k + 1) * 3 + 2];
		    Index[(k + 1) * 3 + 1] = Index[k * 3 + 1];
		    Index[(k + 1) * 3] = Index[(k + 1) * 3 + 2];
		    Index[(k + 1) * 3 + 2] = i;
	        }
	    }
	
	    if (Index[k * 3 + 1] == Index[(k + 1) * 3 + 1] && Index[k * 3 + 2] == Index[(k + 1) * 3])
	    {
	        if (COS (Vertex + 3 * Index[(k + 1) * 3 + 1], 
		         Vertex + 3 * Index[(k + 1) * 3 + 2], 
		         Vertex + 3 * Index[(k + 1) * 3]) *
		    (1 - sqrt (COS (Vertex + 3 * Index[k * 3 + 1],
				    Vertex + 3 * Index[k * 3],
				    Vertex + 3 * Index[k * 3 + 2]) *
			       COS (Vertex + 3 * Index[k * 3 + 1],
				    Vertex + 3 * Index[k * 3],
				    Vertex + 3 * Index[k * 3 + 2]))) +
		    COS (Vertex + 3 * Index[k * 3 + 1],
		         Vertex + 3 * Index[k * 3],
		         Vertex + 3 * Index[k * 3 + 2]) *
		    (1 - sqrt (COS (Vertex + 3 * Index[(k + 1) * 3 + 1], 
				    Vertex + 3 * Index[(k + 1) * 3 + 2],
				    Vertex + 3 * Index[(k + 1) * 3]) *
			       COS (Vertex + 3 * Index[(k + 1) * 3 + 1], 
				    Vertex + 3 * Index[(k + 1) * 3 + 2], 
				    Vertex + 3 * Index[(k + 1) * 3]))) < 0)
	        {
		    i = Index[k* 3 + 2];
		    Index[k * 3 + 2] = Index[(k + 1) * 3 + 2];
		    Index[(k + 1) * 3 + 1] = Index[(k + 1) * 3 + 2];		
		    Index[(k + 1) * 3] = Index[k * 3];
		    Index[(k + 1) * 3 + 2] = i;
	        }
	    }
        }
    }

    for (k = 0; k < l - 1; k++)
    {
	for (j = 0, j2 = 0; j < nF + 2 * nFnew + nEbound + nEbound2 - 2 * (l - 1); j++)
	{
	    if ((trimBound1[boundLoops1[k]] == Index[3 * j + 0] && trimBound1[boundLoops1[k] + 1] == Index[3 * j + 1]) ||
		(trimBound1[boundLoops1[k]] == Index[3 * j + 1] && trimBound1[boundLoops1[k] + 1] == Index[3 * j + 2]) ||
		(trimBound1[boundLoops1[k]] == Index[3 * j + 2] && trimBound1[boundLoops1[k] + 1] == Index[3 * j + 0]))
		j2++;
	}
	printf ("j2 = %d\n", j2);
	if (j2 != 1)
	{
	    for (j = nF + 2 * trimLoops[k] + boundLoops1[k] + boundLoops2[k] - 2 * k; 
		 j < nF + 2 * trimLoops[k + 1] + boundLoops1[k + 1] + boundLoops2[k + 1] - 2 * (k + 1); j++)
	    {
		i = Index[3 * j + 1];
		Index[3 * j + 1] = Index[3 * j + 2];
		Index[3 * j + 2] = i;
	    }
	}
    }
		    
    nV = nV + nFnew;
    nF = nF + 2 * nFnew + nEbound + nEbound2 - 2 * (l - 1);
    
    free (tmpVtx);
    free (tmpIdx);
    free (inside);
    
    free (newVtx);  
    free (trimVtx);

    free (trimLoops);
    free (boundLoops1);
    free (boundLoops2);
    
    free (trimBound1);
    free (trimBound2);
    free (trimIdx);
    
    free (used);
    free (trimParam);

    remove_unused_vertices (&nV, Vertex, &nF, Index);
    printf ("Ready! nV = %d, nF = %d, nFnew = %d, nEbound = %d, nEbound2 = %d\n", nV, nF, nFnew, nEbound, nEbound2);
    return 0;
}


// Checks if a triangle # numF from mesh msh intersects our mesh
int mesh::triInt (mesh *msh, int numF, int *ret, double pt1[3], double pt2[3])
{
    int i, ind, coplanar, total = 0;

    // Input parameters
    double v0[3];
    double v1[3];
    double v2[3];

    double u0[3];
    double u1[3];
    double u2[3];
    
    double isectpt1[3];
    double isectpt2[3];
   
//    CPVECT (u0, msh->Vertex + 3 * msh->Index[3 * numF]);
    ind = msh->Index[3 * numF];
    u0[0] = msh->Vertex[3 * ind];
    u0[1] = msh->Vertex[3 * ind + 1];
    u0[2] = msh->Vertex[3 * ind + 2];
    
//    CPVECT (u1, msh->Vertex + 3 * msh->Index[3 * numF + 1]);
    ind = msh->Index[3 * numF + 1];
    u1[0] = msh->Vertex[3 * ind];
    u1[1] = msh->Vertex[3 * ind + 1];
    u1[2] = msh->Vertex[3 * ind + 2];

//    CPVECT (u2, msh->Vertex + 3 * msh->Index[3 * numF + 2]);
    ind = msh->Index[3 * numF + 2];
    u2[0] = msh->Vertex[3 * ind];
    u2[1] = msh->Vertex[3 * ind + 1];
    u2[2] = msh->Vertex[3 * ind + 2];
    
    for (i = *ret; i < nF; i++)
    {
//	CPVECT (v0, Vertex + 3 * Index[3 * i]);
        ind = Index[3 * i];
	v0[0] = Vertex[3 * ind];
	v0[1] = Vertex[3 * ind + 1];
	v0[2] = Vertex[3 * ind + 2];
	
//	CPVECT (v1, Vertex + 3 * Index[3 * i + 1]);
        ind = Index[3 * i + 1];
	v1[0] = Vertex[3 * ind];
	v1[1] = Vertex[3 * ind + 1];
	v1[2] = Vertex[3 * ind + 2];
	
//	CPVECT (v2, Vertex + 3 * Index[3 * i + 2]);
        ind = Index[3 * i + 2];
	v2[0] = Vertex[3 * ind];
	v2[1] = Vertex[3 * ind + 1];
	v2[2] = Vertex[3 * ind + 2];

	if (tri_tri_intersect_with_isectline (v0, v1, v2, u0, u1, u2, &coplanar, isectpt1, isectpt2))
	{
            total++;
            *ret = i + 1;

//	    CPVECT (pt1, isectpt1);
	    pt1[0] = isectpt1[0];
	    pt1[1] = isectpt1[1];
	    pt1[2] = isectpt1[2];
	    
//	    CPVECT (pt2, isectpt2);
	    pt2[0] = isectpt2[0];
	    pt2[1] = isectpt2[1];
	    pt2[2] = isectpt2[2];

	    return 1;
	}
    } 

    *ret = nF;
    return 0;
}


// Checks if a point (x, y, z) is inside our mesh or outside
int mesh::isInside (double x, double y, double z)
{
    int i, ind, num = 0, res;

    // Input parameters
    double v0[3];
    double v1[3];
    double v2[3];

    double orig[3];
    double dir[3];

    double t, u, v;

    orig[0] = x;
    orig[1] = y;
    orig[2] = z;

    dir[0] = 0.5678f;
    dir[1] = 0.2954f;
    dir[2] = 0.6191f;

    for (i = 0; i < nF; i++)
    {
//	CPVECT (v0, Vertex + 3 * Index[3 * i]);
        ind = Index[3 * i];
	v0[0] = Vertex[3 * ind];
	v0[1] = Vertex[3 * ind + 1];
	v0[2] = Vertex[3 * ind + 2];
	
//	CPVECT (v1, Vertex + 3 * Index[3 * i + 1]);
        ind = Index[3 * i + 1];
	v1[0] = Vertex[3 * ind];
	v1[1] = Vertex[3 * ind + 1];
	v1[2] = Vertex[3 * ind + 2];
	
//	CPVECT (v2, Vertex + 3 * Index[3 * i + 2]);
        ind = Index[3 * i + 2];
	v2[0] = Vertex[3 * ind];
	v2[1] = Vertex[3 * ind + 1];
	v2[2] = Vertex[3 * ind + 2];

	res = intersect_triangle3 (orig, dir, v0, v1, v2, &t, &u, &v);
	
	if (res == 1 && t > 0)
            num++;
    } 

    return num % 2;
}


// Check if the (u,v) line crosses (x,y) segment
int mesh::insideTri (double *u, double *v, double *x, double *y)
{
    double dest[3], dest2[3];
    double v1[3], v2[3];

    v1[0] = x[0] - y[0];
    v1[1] = x[1] - y[1];
    v1[2] = x[2] - y[2];

    v2[0] = v[0] - u[0];
    v2[1] = v[1] - u[1];
    v2[2] = v[2] - u[2];

    CROSS (dest, v1, v2);
    CROSS (dest2, dest, v2);

    v1[0] = x[0] - u[0];
    v1[1] = x[1] - u[1];
    v1[2] = x[2] - u[2];

    v2[0] = y[0] - u[0];
    v2[1] = y[1] - u[1];
    v2[2] = y[2] - u[2];

    if (DOT (v1, dest2) * DOT (v2, dest2) < 0)
        return 1;

    return 0;
}


void mesh::printPoint (double *x)
{
    printf ("(x, y, z) = (%e, %e, %e)\n", x[0], x[1], x[2]);
}


double mesh::dist_point_curve (double *point, double *curve, int len, int *ret)
{
    int i;
    double dist = Meshsize, tmp, a, b;
    for (i = 0; i < len; i++)
    {
	a = DIST (point, curve + i * 6);
	b = DIST (point, curve + i * 6 + 3);

	tmp = MIN (a, b);
	
	if (tmp > Meshsize)
	    continue;
	
	if (tmp < dist)
	    dist = tmp;

	if (MAX (a, b) > Meshsize)
	    continue;

	tmp = DIST (point, curve + i * 6) * 
	      sqrt (1 - COS (point, curve + i * 6, curve + i * 6 + 3) * COS (point, curve + i * 6, curve + i * 6 + 3));
	if (tmp < dist)
	{
	    dist = tmp;
	    *ret = i;
	}
    }

    return dist;
}
