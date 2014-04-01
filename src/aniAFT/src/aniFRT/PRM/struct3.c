#include <string.h>
#include "common.h"
#include "error3.h"
#include "memory3.h"

void init(surface_mesh *pm) {
    pm->nPoint = 0;
    pm->surf = 0;
    pm->nFace = 0;
    initMemory(pm);
    pm->gmesh = myAlloc(sizeof(mesh32));
    memset(pm->gmesh, 0, sizeof(mesh32));
    pm->cf = 0.0;
    pm->sizelim = -1.0;
    return;
}

void addPoint(surface_mesh *pm, double x, double y, double z) {
	if (pm->nPoint >= pm->maxPoint)
		errorExit3(2, "pm->nPoint >= pm->maxPoint");
	pm->vert[pm->nPoint].x = x;
	pm->vert[pm->nPoint].y = y;
	pm->vert[pm->nPoint].z = z;
	pm->nPoint++;
	return;
} /*addPoint*/

