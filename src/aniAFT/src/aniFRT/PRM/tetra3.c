#include "common.h"
#include "tetra3.h"
#include "error3.h"

void addBadVert(surface_mesh *pm, int v) {
	if (pm->nBadVert >= 10000)
		errorExit3(2, "pm->nBadVert >= 10000");
	pm->badVert[pm->nBadVert++] = v;
	return;
} /*addBadVert*/

int isBadVert(surface_mesh *pm, int v) {
	int  i;
	for (i=0; i<pm->nBadVert; i++)
		if (pm->badVert[i] == v) return 1;
	return 0;
} /*isBadVert*/

void addSphereVert(surface_mesh *pm, int v) {
	if( pm->nSphereVert >= 10000 )
		errorExit3(2, "pm->nSphereVert >= 10000");
	pm->sphereVert[pm->nSphereVert++] = v;
	return;
} /*addSphereVert*/

int isSphereVert(surface_mesh *pm, int v) {
	int  i;
	for (i=0; i<pm->nSphereVert; i++)
		if (pm->sphereVert[i] == v) return  1;
	return  0;
} /*isSphereVert*/

