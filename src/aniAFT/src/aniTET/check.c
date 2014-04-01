#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aft.h"
#include "aft3d.h"
#include "helper.h"
#include "det.h"

#define V_NEGATIVE_IS_FAIL 1
#define COMPACT_OUTPUT 1

int checktopology(int nV, REAL *vertex, int nF, int *face, int nT, int *tetra, int indshift) {
    int allisok=1;
    int *vi, *vj, *tf;
    int i, j, k, m, v, v1, v2, v3, i1, i2, i3;
    REAL vol=0.0, vol2=0.0, tvol, tmin = 0.0;
    int np = 0, *ps;

    if (COMPACT_OUTPUT) {
	ps = (int*)libaft_malloc(sizeof(int)*nV);
    }

    for (i=0; i<3*nF; face[i++]-=indshift);
    for (i=0; i<4*nT; tetra[i++]-=indshift);

    for (i=0; i<nF; i++) {
	vol += det3i4(vertex, face[3*i+0], face[3*i+1], face[3*i+2], 0) / 6.0;
    }
    if (vol<=0) {
	libaft_3d_warn("Check: Inverted front, volume = %le", vol);
	allisok = 0;
    }
    for (i=0; i<nT; i++) {
	tvol = det3i4(vertex, tetra[4*i+0], tetra[4*i+1], tetra[4*i+2], tetra[4*i+3]) / 6.0;
	if (tvol<=0.0) {
	    if ((COMPACT_OUTPUT) && (ps)) {
		ps[np++] = i;
		if (tvol < tmin)  tmin = tvol;
	    } else {
		libaft_3d_warn("Check: Tetra %d volume = %le", i + indshift, tvol);
	    }
	    if (V_NEGATIVE_IS_FAIL)  allisok = 0;
	}
	vol2 += tvol;
    }
    if (fabs(vol2-vol)/vol > 1e-8) {
	libaft_3d_warn("Check: Volumes differ. Original = %lg, meshed = %lg, diff = %le", vol, vol2, vol-vol2);
	allisok = 0;
    }

    vi = (int*)libaft_malloc(sizeof(int)*(nV+1));
    vj = (int*)libaft_malloc(sizeof(int)*4*nT);
    for (i=0; i<=nV; i++) vi[i] = 0;
    for (j=0; j<4*nT; j++) {
	vi[tetra[j]+1]++;
    }
    for (i=1; i<nV; i++) vi[i] += vi[i-1];
    for (j=0; j<4*nT; j++) {
	vj[vi[tetra[j]]++] = j/4;
    }
    for (i=nV; i>0; i--) vi[i] = vi[i-1];
    vi[0] = 0;

/*    for (i=0; i<nV; i++) {
	printf("[%d]:", i);
	for (j=vi[i]; j<vi[i+1]; j++)
	    printf(" %d", vj[j]);
	printf("\n");
    }*/

    tf = (int*)libaft_malloc(sizeof(int)*4*nT);
    for (j=0; j<nT; j++) {
	for (k=0; k<4; k++) {
	    v = tetra[4*j+k];
	    v1 = tetra[4*j+(k+1)%4];
	    v2 = tetra[4*j+(k+2)%4];
	    v3 = tetra[4*j+(k+3)%4];
	    m = -1;
	    for  (i1 = vi[v1], i2 = vi[v2], i3 = vi[v3];
		    (i1<vi[v1+1])&&(i2<vi[v2+1])&&(i3<vi[v3+1]);
		 ) {
		if ((vj[i1] == vj[i2]) && (vj[i1] == vj[i3])) {
		    if (vj[i1]!=j) {
			if (m==-1) m = vj[i1];
			else {
			    libaft_3d_warn("Check: Duplicate tetras");
			    allisok = 0;
			}
		    }
		    i1++, i2++, i3++;
		} else if ((vj[i1] <= vj[i2]) && (vj[i1] <= vj[i3])) {
		    i1++;
		} else if ((vj[i2] <= vj[i1]) && (vj[i2] <= vj[i3])) {
		    i2++;
		} else if ((vj[i3] <= vj[i1]) && (vj[i3] <= vj[i2])) {
		    i3++;
		}
	    }
	    tf[4*j+k] = m;
	    if (m<0) {
		for (m=0; m<nF; m++) {
		    if (  ((v1==face[3*m+0])||(v1==face[3*m+1])||(v1==face[3*m+2])) &&
			    ((v2==face[3*m+0])||(v2==face[3*m+1])||(v2==face[3*m+2])) &&
			    ((v3==face[3*m+0])||(v3==face[3*m+1])||(v3==face[3*m+2]))  ) {
			tf[4*j+k] = nT+m;
		    }
		}
		if (tf[4*j+k]<0) {
		    libaft_3d_warn("Check: Tetra does not have neighbour tetra [%d,%d] (%d, %d %d %d)", j + indshift, k,
			    v + indshift, v1 + indshift, v2 + indshift, v3 + indshift);
		    allisok = 0;
		}
	    }
	}
    }
    if ((np) && (COMPACT_OUTPUT) && (ps)) {

	libaft_3d_warn("Check: %d tetras with negative volume (up to %le). Dumping to stderr:", np, tmin);
	for (i=0; i<np; i++)  fprintf(stderr, (i<np-1)?"%d ":"%d\n", ps[i] + indshift);
    }

    libaft_free(vi),  libaft_free(vj),  libaft_free(tf),  libaft_free(ps);

    for (i=0; i<3*nF; face[i++]+=indshift);
    for (i=0; i<4*nT; tetra[i++]+=indshift);

    return (!allisok);
}




///// 2D /////



static int cmpvxyz(double *vertex, int v1, int v2) {
    if (vertex[3*v1+0] < vertex[3*v2+0]) return -1;
    else if (vertex[3*v1+0] > vertex[3*v2+0]) return +1;
    else if (vertex[3*v1+1] < vertex[3*v2+1]) return -1;
    else if (vertex[3*v1+1] > vertex[3*v2+1]) return +1;
    else if (vertex[3*v1+2] < vertex[3*v2+2]) return -1;
    else if (vertex[3*v1+2] > vertex[3*v2+2]) return +1;
    else return 0;
}

static double *vcmp;

static int cmpxyz(const void *pv1, const void *pv2) {
    int v1 = *(int *)pv1;
    int v2 = *(int *)pv2;
    return cmpvxyz(vcmp, v1, v2);
}

static int check_xyz(int nV, double *vertex, int *index) {
    int i, j;
    int errors = 0, warnings = 0;

    vcmp = vertex;

    for (i=0; i<nV; i++) {
	index[i] = i;
    }
    qsort(index, nV, sizeof(int), cmpxyz);
    j = 0;
    for (i=1; i<nV; i++) {
	if ( !(cmpxyz(&index[i-1], &index[i])) ) {
	    if ( !(j) )	errors++;
	    j++;
	} else {
	    j = 0;
/*	    if (j) {
		printf("Check: %d duplicating vertices (%d", j+1, index[i-j-1]+1);
		while (j>0) {
		    printf(", %d", index[i-j]+1);
		    j--;
		}
		printf(")\n");
		j = 0;
	    }*/
	}
    }
    if (errors+warnings > 0)
	libaft_3d_warn("Check: %d duplicate vertices", errors);
    return errors+warnings;
}

static int fix_xyz(int *pnV, double *vertex, int *pnF, int *face, int *index) {
    int i, n = *pnV, m;
    int *rename;
    double t[3];

//    printf("Fixing xyz data.\n");
    rename = (int*)libaft_malloc(sizeof(int)*n);
    rename[index[0]] = 0; m=1;
    for (i=1; i<n; i++) {
	if ( !(cmpvxyz(vertex, index[i-1], index[i])) ) {
	    rename[index[i]] = rename[index[i-1]];
	} else {
	    rename[index[i]] = m;
	    index[m] = index[i];
	    m++;
	}
    }
    *pnV = m;
//    printf("nV: %d -> %d\n", n, m);
    for (i=0; i<(*pnF); i++) {
	face[3*i+0] = rename[face[3*i+0]];
	face[3*i+1] = rename[face[3*i+1]];
	face[3*i+2] = rename[face[3*i+2]];
    }
    for (i=0; i<m; i++) {
	if (i!=index[i]) {
	    memcpy(t, vertex+3*i, 3*sizeof(double));
	    memcpy(vertex+3*i, vertex+3*index[i], 3*sizeof(double));
	    memcpy(vertex+3*index[i], t, 3*sizeof(double));
	    index[rename[i]] = index[i];
	    rename[index[i]] = rename[i];
	    index[i] = i;
	    rename[i] = i;
	}
    }
    libaft_free(rename);
    return 0;
}

static int fix_tri(int *pnV, double *vertex, int *pnF, int *face) {  // not implemented yet
    (void)pnV,(void)vertex,(void)pnF,(void)face;
    return 0;
}

static int check_tri(int nF, int *face, int indshift, int verbose) {
    int i, j, k, p, nV = 0;
    int *vertex, *edge;

    int errors = 0, warnings = 0;

    if (verbose)  printf("Checking tri data.\n");
    for (i=0; i<nF; i++) {
	for (j=0; j<3; j++) {
	    if (face[3*i+j] > nV)
		nV = face[3*i+j];
	}
    }
    if (verbose)  printf("max vertex index used: %d\n", nV+indshift);
    nV++;
    vertex = (int*)libaft_malloc(sizeof(int)*nV);
    edge   = (int*)libaft_malloc(2*sizeof(int)*3*nF);
    for (i=0; i<nV; i++) {
	vertex[i] = -1;
    }
    for (i=0; i<nF; i++) {
	for (j=0; j<3; j++) {
	    edge[2*(3*i+j)] = face[3*i+(j+1)%3];
	    edge[2*(3*i+j)+1] = vertex[face[3*i+j]];
	    vertex[face[3*i+j]] = 3*i+j;
	}
    }
    for (i=0; i<nF; i++) {
	for (j=0; j<3; j++) {
	    k = vertex[face[3*i+j]]; p = 0;
	    while (k >= 0) {
		if (edge[2*k] == face[3*i+(j+2)%3]) p++;
		k = edge[2*k+1];
	    }
	    if (p == 0) {
		if (verbose)  libaft_3d_warn("Check: no pair for edge [%d %d]\n", face[3*i+j]+indshift, face[3*i+(j+2)%3]+indshift);
		errors++;
	    } else if (p > 1) {
		if (verbose)  printf("warning: more then one pair for edge [%d %d] (actually  %d pairs)\n", face[3*i+j]+indshift, face[3*i+(j+2)%3]+indshift, p);
		warnings++;
	    }
	}
    }
    libaft_free(vertex);
    libaft_free(edge);
    if (errors/*+warnings*/ > 0)
	libaft_3d_warn("Check: %d errors, %d warnings\n", errors, warnings);
    return errors/*+warnings*/;
}

int check2dsurface(int *pnV, double *vertex, int *pnF, int *face, int indshift, int autofix, int verbose) {
    int nV=*pnV, nF=*pnF, i;
    int *index;
    int result;

    for (i=0; i<3*nF; face[i++]-=indshift);
    index = (int*)libaft_malloc(sizeof(int)*nV);
    result=0;
    if (check_xyz(nV, vertex, index) == 0) {
	if (verbose)  libaft_3d_warn("Check: xyz coords - OK\n");
    } else if ((autofix) && (fix_xyz(pnV, vertex, pnF, face, index) == 0)) {
	if (verbose)  libaft_3d_warn("Check: updated xyz coords and topology\n");
    } else {
	result++;
    }
    nV = *pnV;
    nF = *pnF;
    if (check_tri(nF, face, indshift, verbose) == 0) {
	if (verbose)  libaft_3d_warn("Check: surface topology - OK\n");
    } else if ((autofix) && (fix_tri(pnV, vertex, pnF, face) == 0)) {
	if (verbose)  libaft_3d_warn("Check: updated surface topology\n");
    } else {
	result++;
    }
    libaft_free(index);
    for (i=0; i<3*nF; face[i++]+=indshift);
    return result;
}
