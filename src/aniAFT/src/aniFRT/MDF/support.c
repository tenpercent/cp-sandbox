#include <stdlib.h>
#include <string.h>
#include "surfrf.h"
#include "support.h"

int aft3dsurfmeshrefiner(int *nV, double *vertex, int *nF, int *face, int *facematerial, int maxnV, int maxnF, int indexshift) {
	int r, i;
	indexshift = indexshift-1;
	if (indexshift) for (i=0; i<3**nF; face[i++]-=indexshift);
	r = libaft_internal_surfmeshrefiner(nV, vertex, nF, face, facematerial, maxnV, maxnF);
	if (indexshift) for (i=0; i<3**nF; face[i++]-=indexshift);
	return r;
}

