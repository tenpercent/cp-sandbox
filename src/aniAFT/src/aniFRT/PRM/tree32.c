#include "common.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tree32.h"



/* Metric functions ****************************************************************************************************************************** METRIC */
static double dist3(double x1, double y1, double z1,  double x2, double y2, double z2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

/* Octree functions *********************************************************/;
#define TREEDEPTH 24
static void etreemetric(mesh32 *pm2, surface_mesh *pm, int w1, int w2, double *px, double *pr) {
    int v1, v2;
    double r, s, x, y, z;
    (void) pm2;
    
    if (w1 < w2)  v1 = w1,  v2 = w2;
    else  v1 = w2,  v2 = w1;
    x = px[0] = (pm->vert[v1].x + pm->vert[v2].x)/2.0;
    y = px[1] = (pm->vert[v1].y + pm->vert[v2].y)/2.0;
    z = px[2] = (pm->vert[v1].z + pm->vert[v2].z)/2.0;
    r = 0.0;
    s = dist3(pm->vert[v1].x, pm->vert[v1].y, pm->vert[v1].z, x, y, z);
    if (s > r) r = s;
    s = dist3(pm->vert[v2].x, pm->vert[v2].y, pm->vert[v2].z, x, y, z);
    if (s > r) r = s;
    *pr = r;
}
static void emetric(mesh32 *pm2, surface_mesh *pm, int i) {
    int v1, v2;
    v1 = pm2->e[i].v[0];
    v2 = pm2->e[i].v[1];
    etreemetric(pm2, pm, v1, v2, pm2->e[i].x, &pm2->e[i].r);
    pm2->e[i].s = pm2->e[i].r;
}
static void findbb(mesh32 *pm2, surface_mesh *pm) {
    int i;
    double d;
    if (pm->nPoint == 0) {
	pm2->bb[0] = 0.0,  pm2->bb[1] = 0.0,  pm2->bb[2] = 0.0;
	pm2->bb[3] = 0.0,  pm2->bb[4] = 0.0,  pm2->bb[5] = 0.0;
	return;
    }
    pm2->bb[0] = pm->vert[0].x,  pm2->bb[1] = pm->vert[0].y,  pm2->bb[2] = pm->vert[0].z;
    pm2->bb[3] = pm->vert[0].x,  pm2->bb[4] = pm->vert[0].y,  pm2->bb[5] = pm->vert[0].z;
    for (i=1; i<pm->nPoint; i++) {
	if (pm2->bb[0] > pm->vert[i].x)  pm2->bb[0] = pm->vert[i].x;
	if (pm2->bb[1] > pm->vert[i].y)  pm2->bb[1] = pm->vert[i].y;
	if (pm2->bb[2] > pm->vert[i].z)  pm2->bb[2] = pm->vert[i].z;
	if (pm2->bb[3] < pm->vert[i].x)  pm2->bb[3] = pm->vert[i].x;
	if (pm2->bb[4] < pm->vert[i].y)  pm2->bb[4] = pm->vert[i].y;
	if (pm2->bb[5] < pm->vert[i].z)  pm2->bb[5] = pm->vert[i].z;
    }
    d = 0.0;
    for (i=0; i<3; i++) {
	if (d < pm2->bb[3+i]-pm2->bb[i])  d = pm2->bb[3+i]-pm2->bb[i];
    }
    for (i=0; i<3; i++) {
	pm2->bb[i] = (pm2->bb[i] + pm2->bb[3+i] - d)/2.0;
	pm2->bb[3+i] = pm2->bb[i] + d;
    }
    pm2->bb[0] = pm->boxcx - pm->boxsize;
    pm2->bb[1] = pm->boxcy - pm->boxsize;
    pm2->bb[2] = pm->boxcz - pm->boxsize;
    pm2->bb[3] = pm->boxcx + pm->boxsize;
    pm2->bb[4] = pm->boxcy + pm->boxsize;
    pm2->bb[5] = pm->boxcz + pm->boxsize;
}
static int treedir(double *bb, double *x) {
    int d = 0;
    if (x[0] > (bb[0]+bb[3])/2.0) d+=1;
    if (x[1] > (bb[1]+bb[4])/2.0) d+=2;
    if (x[2] > (bb[2]+bb[5])/2.0) d+=4;
    return d;
}
static void treemove(double *bb, double s, int d) {
    if (d&1)  bb[0] += s;  else  bb[3] -= s;
    if (d&2)  bb[1] += s;  else  bb[4] -= s;
    if (d&4)  bb[2] += s;  else  bb[5] -= s;
}
static void treeback(double *bb, double s, int d) {
    if (d&1)  bb[0] -= s;  else  bb[3] += s;
    if (d&2)  bb[1] -= s;  else  bb[4] += s;
    if (d&4)  bb[2] -= s;  else  bb[5] += s;
}
static int listadd(mesh32 *pm2, int *pn, int *pnn, int **pe, int n) {
#ifdef STRUCTCHECK
    int i;
    int v1, v2, w1, w2;
    v1 = pm2->e[n].v[0];
    v2 = pm2->e[n].v[1];
    for (i=0; i<*pn; i++) {
	w1 = pm2->e[(*pe)[i]].v[0];
	w2 = pm2->e[(*pe)[i]].v[1];
	if ((v1!=w1)&&(v1!=w2)) continue;
	if ((v2!=w1)&&(v2!=w2)) continue;
	printf("pm->tc: listadd(): duplicate edges detected (%d %d)\n", v1, v2);
    }
#endif
    (void) pm2;
    if (*pn >= *pnn) {
	*pnn += 16;
	*pe = realloc(*pe, sizeof(int)* *pnn);
    }
    (*pe)[*pn] = n;
    (*pn)++;
    return 0;
}
static int listdel(int *pn, int *pnn, int **pe, int n) {
    int i;
    for (i=0; i<*pn; i++) {
	if ((*pe)[i] == n) {
	    (*pn)--;
	    (*pe)[i] = (*pe)[*pn];
	    if (*pn == 0) {
		*pnn = 0;
		free(*pe);
		*pe = 0;
	    }
	    return 0;
	}
    }
    printf("pm->tc: listdel(): element not found (internal error)\n");
    return -1;
}
static int listrep(int *pn, int **pe, int n, int m) {
    int i;
    if (n==m) return 0;
    for (i=0; i<*pn; i++) {
	if ((*pe)[i] == n) {
	    (*pe)[i] = m;
	    return 0;
	}
    }
    printf("pm->tc: listrep(): element not found (internal error)\n");
    return -1;
}
static octreenode *nodenew(void) {
    octreenode *p;
    int i;
    p = malloc(sizeof(octreenode));
    p->ne = 0;
    p->nn = 0;
    p->e = 0;
    for (i=0; i<8; i++)  p->child[i] = 0;
    return p;
}
static int nodevoid(octreenode *n) {
    int i;
    if (n->nn)  return 0;
    for (i=0; i<8; i++)  if (n->child[i])  return 0;
    if (n->e)  printf("pm->tc: nodevoid(): node void, but list is not NULL (sure memory leak)");
    return 1;
}
static octreenode *treeadd(mesh32 *pm2, int n) {
    octreenode *c, *p;
    double bb[6], x[3], s, r;
    int i, d, level;

    p = pm2->root;
    if (!p) {
	p = nodenew();
	pm2->root = p;
    }
    for (i=0; i<6; i++)  bb[i] = pm2->bb[i];
    for (i=0; i<3; i++)  x[i] = pm2->e[n].x[i];
    r = pm2->e[n].r;
    s = (pm2->bb[3]-pm2->bb[0])/2.0;
    c = p;
    level = 0;
    while ((level<TREEDEPTH) && (s > r)) {
	d = treedir(bb, x);
	treemove(bb, s, d);
	s /= 2.0;
	level++;
	p = c->child[d];
	if (!p) {
	    p = nodenew();
	    c->child[d] = p;
	}
	c = p;
    }
    listadd(pm2, &c->ne, &c->nn, &c->e, n);
    return c;
}
static int treedel(mesh32 *pm2, int n) {
    octreenode *c, *p;
    octreenode *stack[TREEDEPTH];
    int dirs[TREEDEPTH];
    double bb[6], x[3], s, r;
    int i, d, level;

    p = pm2->root;
    if (!p)  {
	printf("pm->tc: treedel(): null root in octree (internal error)");
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm2->bb[i];
    for (i=0; i<3; i++)  x[i] = pm2->e[n].x[i];
    r = pm2->e[n].r;
    s = (pm2->bb[3]-pm2->bb[0])/2.0;
    c = p;
    level = 0;
    while ((level<TREEDEPTH) && (s > r)) {
	d = treedir(bb, x);
	stack[level] = c;
	dirs[level] = d;
	treemove(bb, s, d);
	s /= 2.0;
	level++;
	p = c->child[d];
	if (!p) {
	    printf("pm->tc: treedel(): null child in octree (internal error)");
	    return -1;
	}
	c = p;
    }
    if (listdel(&c->ne, &c->nn, &c->e, n))  return -1;
    while ((level>0) && (nodevoid(c))) {
	free(c);
	level--;
	c = stack[level];
	c->child[dirs[level]] = 0;
    }
    if ((level==0) && (nodevoid(c))) {
	free(c);
	pm2->root = 0;
    }
    return 0;
}
static int wscheck(mesh32 *pm2, int i, double *x, double r) {
    return (dist3(pm2->e[i].x[0], pm2->e[i].x[1], pm2->e[i].x[2], x[0], x[1], x[2]) <= pm2->e[i].r + r) ? 1 : 0;
}
static int treewsrec(mesh32 *pm2, surface_mesh *pm, octreenode *c, double *bb, double s, double *x, double r) {
    int i, k = 0;
    if (!c)  return 0;
    for (i=0; i<3; i++) {
	if (bb[  i] - 2.0*s > x[i] + r)  return 0;
	if (bb[3+i] + 2.0*s < x[i] - r)  return 0;
    }
    for (i=0; i<c->ne; i++) {
	if (wscheck(pm2, c->e[i], x, r)) {
	    pm->tvicinityFace[pm->tnVicinityFace].face = pm->edge[pm2->e[c->e[i]].h];
	    pm->tvicinityFace[pm->tnVicinityFace].dist = 0.0;
	    pm->tnVicinityFace++;
	    if (pm->tnVicinityFace>pm->tmaxVicinityFace)
		errorExit3(3,"pm->tnVicinityFace");
	}
    }
    k += c->ne;
    for (i=0; i<8; i++) {
	treemove(bb, s, i);
	k += treewsrec(pm2, pm, c->child[i], bb, s/2.0, x, r);
	treeback(bb, s, i);
    }
    return k;
}
static int treews(mesh32 *pm2, surface_mesh *pm, double *x, double r) {
    octreenode *p;
    double bb[6], s;
    int i;

    p = pm2->root;
    if (!p)  {
	printf("pm->tc: treews(): null root in octree (internal error)");
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm2->bb[i];
    s = (pm2->bb[3]-pm2->bb[0])/2.0;
    return treewsrec(pm2, pm, p, bb, s, x, r);
}
/*static int listfindedge(mesh32 *pm2, int ne, int *e, int v1, int v2) {
    int i, w1, w2;
    for (i=0; i<ne; i++) {
	w1 = pm2->e[e[i]].v[0];
	w2 = pm2->e[e[i]].v[1];
	if ((v1==w1)&&(v2==w2)) return e[i];
	if ((v2==w1)&&(v1==w2)) return e[i];
    }
    return -1;
}
static int findedge(mesh32 *pm2, int v1, int v2) {
    octreenode *c, *p;
    double bb[6], x[3], s, r;
    int i, d, level;

    p = pm2->root;
    if (!p)  {
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm2->bb[i];
    etreemetric(pm2, pm, v1, v2, x, &r);
    s = (pm2->bb[3]-pm2->bb[0])/2.0;
    c = p;
    level = 0;
    while ((level<TREEDEPTH) && (s > r)) {
	d = treedir(bb, x);
	treemove(bb, s, d);
	s /= 2.0;
	level++;
	p = c->child[d];
	if (!p) {
	    return -1;
	}
	c = p;
    }
    return listfindedge(pm2, c->ne, c->e, v1, v2);
}*/

static int cc(surface_mesh *pm) {
//    return 0;
    int i, r;
    r = 0;
    if (pm->gmesh->nE != pm->nEdge)  printf("cc: nE != nFace\n"),  r++;
    for (i=0; i<pm->gmesh->nE; i++) {
	if ((pm->gmesh->e[pm->gmesh->heap[i]].v[0] != pm->edge[i]->v1) || (pm->gmesh->e[pm->gmesh->heap[i]].v[0] != pm->edge[i]->v1))
	    printf("cc: wrong edge\n"),  r++;
    }
    return r;
}

/* Structure functions ******************************************************/;
static void heapswap(mesh32 *pm2, surface_mesh *pm, int k, int m) {
    int t;
    PStrucEdge3 edge;
    edge = pm->edge[k],  pm->edge[k] = pm->edge[m],  pm->edge[m] = edge;
    t = pm2->heap[k],  pm2->heap[k] = pm2->heap[m],  pm2->heap[m] = t;

    pm->edge[k]->f = k;
    pm->edge[m]->f = m;
    pm2->e[pm2->heap[k]].h = k;
    pm2->e[pm2->heap[m]].h = m;
}
static int siftup(mesh32 *pm2, surface_mesh *pm, int n) {
    int i;
    while (1) {
	i = (n-1)/2;
	if ((i<0) || (pm2->e[pm2->heap[i]].s <= pm2->e[pm2->heap[n]].s)) break;
	heapswap(pm2, pm, i, n);
	n = i;
    }
    return n;
}
static int siftdown(mesh32 *pm2, surface_mesh *pm, int n) {
    int i;
    n = siftup(pm2, pm, n);
    while (1) {
	i = 0;
	if ((2*n+1 < pm2->nE) && (pm2->e[pm2->heap[n]].s > pm2->e[pm2->heap[(2*n+1)]].s)) {
	    if ((2*n+2 < pm2->nE) && (pm2->e[pm2->heap[(2*n+1)]].s > pm2->e[pm2->heap[(2*n+2)]].s))  i = 2*n+2;
	    else  i = 2*n + 1;
	} else if ((2*n+2 < pm2->nE) && (pm2->e[pm2->heap[n]].s > pm2->e[pm2->heap[(2*n+2)]].s))  i = 2*n+2;
	if (i) {
	    heapswap(pm2, pm, n, i);
	    n = i;
	} else break;
    }
    return n;
}
static int addedge(mesh32 *pm2, surface_mesh *pm, int v1, int v2) {
#ifdef STRUCTCHECK
    int i;
    for (i=0; i<pm2->nE; i++) {
	if ((v1!=pm2->e[i].v[0])&&(v1!=pm2->e[i].v[1])) continue;
	if ((v2!=pm2->e[i].v[0])&&(v2!=pm2->e[i].v[1])) continue;
	printf("pm->tc: addface(): duplicate faces detected (bad front?)\n");
    }
#endif
    if (pm2->nE >= pm2->nnE) {
	printf("AFT Front: maxFace exceeded\n");
	return -1;
    } else {
	pm2->e[(pm2->nE)].v[0] = v1;
	pm2->e[(pm2->nE)].v[1] = v2;
	emetric(pm2, pm, pm2->nE);
	pm2->heap[pm2->nE] = pm2->nE;
	pm2->e[pm2->nE].h = pm2->nE;
	pm2->e[pm2->nE].node = treeadd(pm2, pm2->nE);
	return siftup(pm2, pm, (pm2->nE)++);
    }
}
static void remedge(mesh32 *pm2, int i) {
    if (pm2->npack > pm2->nnE) {
	printf("pm->tc: Internal error in remface & packface\n");
	return;
    }
    pm2->pack[pm2->npack] = i;
    pm2->npack++;
}
static int packedge(mesh32 *pm2, surface_mesh *pm) {
    int i, j;
    int err = 0;
    octreenode *c;
    PStrucEdge3 f;
    while (pm2->npack > 0) {
	pm2->npack--;
	j = pm2->pack[pm2->npack];
	pm2->nE--;
	for (i=0; i<pm2->npack; i++) if (pm2->pack[i] == pm2->nE) pm2->pack[i] = j;
	c = pm2->e[(pm2->nE)].node;
	if (treedel(pm2, j))  err++;
	if (listrep(&c->ne, &c->e, pm2->nE, j))  err++;
	pm2->e[j].v[0] = pm2->e[(pm2->nE)].v[0];
	pm2->e[j].v[1] = pm2->e[(pm2->nE)].v[1];
	pm2->e[j].x[0] = pm2->e[(pm2->nE)].x[0];
	pm2->e[j].x[1] = pm2->e[(pm2->nE)].x[1];
	pm2->e[j].x[2] = pm2->e[(pm2->nE)].x[2];
	pm2->e[j].r = pm2->e[(pm2->nE)].r;
	pm2->e[j].s = pm2->e[(pm2->nE)].s;
	pm2->e[j].node = c;
	i = pm2->e[j].h;
	f = pm->edge[i];
	pm->edge[i] = pm->edge[pm2->e[pm2->nE].h];
	pm->edge[i]->f = i;
	pm->edge[pm2->e[pm2->nE].h] = f;
	f->f = pm2->nE;
	pm2->nE++;
	siftdown(pm2, pm, i);
	pm2->nE--;
	i = pm2->e[pm2->nE].h;
	f = pm->edge[i];
	pm->edge[i] = pm->edge[pm2->nE];
	pm->edge[i]->f = i;
	pm->edge[pm2->nE] = f;
	f->f = pm2->nE;
	pm2->heap[i] = pm2->heap[pm2->nE];
	pm2->e[pm2->heap[i]].h = i;
	siftdown(pm2, pm, i);
    }
    return err;
}

#define EPS 1e-8
PStrucEdge3 addFace32(surface_mesh *pm, int v1, int v2 ) {
    double       x,y,z;
    PStrucEdge3  face;

    face = myAlloc( S_StrucEdge3 );
    face->fail = 0;
    face->v1 = v1;   face->v2 = v2;  
    x = (pm->vert[v1].x+pm->vert[v2].x)/2.;
    y = (pm->vert[v1].y+pm->vert[v2].y)/2.;
    z = (pm->vert[v1].z+pm->vert[v2].z)/2.;
    face->x=x;   face->y=y;   face->z=z;
    face->s = dist3(pm->vert[v1].x,pm->vert[v1].y,pm->vert[v1].z,pm->vert[v2].x,pm->vert[v2].y,pm->vert[v2].z);

    if ((fabs(x-pm->boxcx)>pm->boxsize+EPS) || (fabs(y-pm->boxcy)>pm->boxsize+EPS) || (fabs(z-pm->boxcz)>pm->boxsize+EPS)) {
	printf("\nbox: x: %le, y: %le, z: %le, size: %le\n", pm->boxcx, pm->boxcy, pm->boxcz, pm->boxsize);
	printf("x: %le, y: %le, z: %le\n", x, y, z);
	errorExit3(3,"x_y_z_insert");
    }
    face->f = pm->nEdge;
    pm->edge[pm->nEdge++] = face;
    if (addedge(pm->gmesh, pm, v1, v2) != face->f)  printf("oops! recheck pm->tc\n");
    cc(pm);
    return face;
}
#undef EPS

void remFace32(surface_mesh *pm, PStrucEdge3  face) {
    remedge(pm->gmesh, pm->gmesh->heap[face->f]);
    packedge(pm->gmesh, pm);
    pm->nEdge--;
    free(face);
    cc(pm);
}

double nearest32(surface_mesh *pm, int *vert, double x, double y, double z ) {
    int          i,j,vn=0,v[2];
    double       p=0,dist=0.0;
    PStrucEdge3  face;

    for(i=0;i<pm->tnVicinityFace;i++){
	face = pm->tvicinityFace[i].face;
	v[0] = face->v1;  v[1] = face->v2; 
	for(j=0;j<2;j++){
	    p = dist3(pm->vert[v[j]].x,pm->vert[v[j]].y,pm->vert[v[j]].z,x,y,z);
	    if((p < dist) || ((i==0)&&(j==0))){
		dist = p;
		vn = v[j];
	    }
	}
    }
    vert[0] = vn;

    return( dist );
}

void vicinityFaces32(surface_mesh *pm, double x, double y, double z, double size) {
    double xyz[3];

    pm->tnVicinityFace = 0;
    xyz[0] = x,  xyz[1] = y,  xyz[2] = z;
    treews(pm->gmesh, pm, xyz, size);
}

void prepTree32(surface_mesh *pm, int nnE) {
    pm->gmesh->nnE = nnE;
    pm->gmesh->e = realloc(pm->gmesh->e, sizeof(edge3d)*nnE);
    pm->gmesh->heap = realloc(pm->gmesh->heap, sizeof(int)*nnE);
    pm->gmesh->pack = realloc(pm->gmesh->pack, sizeof(int)*nnE);
    findbb(pm->gmesh, pm);
}


