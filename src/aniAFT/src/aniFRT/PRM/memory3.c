#include "common.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tree32.h"

void *myAlloc (size_t n) {
    void  *p;

    p = (void *)malloc(n);
    if (p == NULL)
	errorExit3(2,"myAlloc");

    return  p;
} /*myAlloc*/


void initMemory(surface_mesh *pm) {
    int        i,allSize=0,size[100];
    size_t     maxMemory=0,allMemory=0,sMemory[100];
    char  *    pMemory=NULL;
    int maxTetra=0;

    maxTetra = 6400000;

    //	printf("maxTetra = %7d\n", maxTetra);

    maxMemory = 48*maxTetra;


    if (maxMemory == 0)
	errorExit3(2, "Memory");

    /* init  distribution  of  memory */
    size[0] = 1 * S_StrucVert3;    /*VERT*/
    size[1] = 0; /*5 * S_StrucTetra3;*/   /*TETRA*/

    allMemory = 0;
    for (i=0; i<2; i++) allSize += size[i];
    for (i=0; i<2; i++) {
	sMemory[i] = maxMemory*size[i] / allSize;
	allMemory += sMemory[i];
    }
    /* end  distribution  of  memory */

    /* init  memory  ptr */
    pMemory = (char*)myAlloc(allMemory);

    pm->maxPoint = sMemory[0] / S_StrucVert3;
    pm->vert = (PStrucVert3)pMemory;
    pMemory += pm->maxPoint*S_StrucVert3;
    /* end  init  memory  ptr */

    pm->maxFace = 2 * pm->maxPoint;
    if (pm->maxFace < 3000) pm->maxFace = 3000;
    pm->face = myAlloc(pm->maxFace * sizeof(StrucFace3));
    pm->maxVicinityFace = 109990;
    pm->vicinityFace = myAlloc(pm->maxVicinityFace * sizeof(StrucFace4));

    pm->nnEdge = pm->maxPoint/5;
    if (pm->nnEdge < 2000) pm->nnEdge = 2000;
    pm->edge = myAlloc(pm->nnEdge * sizeof(PStrucEdge3));
    pm->tmaxVicinityFace = 4409;
    pm->tvicinityFace = myAlloc(pm->tmaxVicinityFace * sizeof(StrucFace4));

    pm->badVert = myAlloc(10000*sizeof(int));
    pm->sphereVert = myAlloc(10000*sizeof(int));

    return;
} /*initMemory*/

void freeMemory(surface_mesh *pm) {
    int i;
    free((void*)pm->vert);
    free(pm->face);
    free(pm->vicinityFace);
    for (i=pm->nEdge-1; i>=0; i--)  remFace32(pm, pm->edge[i]);
    free(pm->edge);
    free(pm->tvicinityFace);
    free(pm->badVert);
    free(pm->sphereVert);
    if (pm->gmesh) {
	free(pm->gmesh->e);
	free(pm->gmesh->heap);
	free(pm->gmesh->pack);
	free(pm->gmesh);
	pm->gmesh = NULL;
    }
    return;
} /*freeMemory*/


