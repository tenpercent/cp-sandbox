#include <stdio.h>
#include <math.h>
#include <time.h>
#include "aft.h"
#include "aft2d.h"
#include "helper.h"
#include "det.h"
#include "refine.h"

// PROF -- enable timing
#define PROF

//#define DUMP
//#define SHOWPROGRESS

typedef int tedge[2];
typedef REAL tvertex[3];
typedef struct {
    int *ia;
    int *ja;
} neigh_edges;

typedef struct {
    int *ia;
    tedge *ja;
} neigh_trias;

typedef struct _tnode{ /* node */
    int ne, nn;
    int *e;
    struct _tnode *child[4];
} node;

typedef struct { /* edge2d */ 
    int v[2];
    int c;
    REAL x[2];
    REAL r;
    REAL s;
    int h;
    node *node;
} edge2d;

typedef struct {
    int  nV, nnV;
    REAL *vertex;
    int  nE, nnE;
    edge2d *e;
    int  *heap;
    int  nT, nnT;
    int  *tria, *triacolor;
    REAL volume, meshedvolume;
    REAL stepsize, minrho, alpha, beta, maxlength, sizelim;
    double (*fsize)(double, double);
    int  curfriend[2], curfriendlink[2];
    REAL curfriendness[2];
    int  npack, *pack;
    int  nWorkArea, nWorkSpace, nBadGuys, nGirls;
    int  *workarea, *workspace, *badguys, *girls;
    REAL *workspace_q;
    int  intedge;
    int  local;
    int  ncolor, *color;
    neigh_edges eadj;
    neigh_trias tadj;
    REAL bb[4];
    node *root;
} mesh2d;



#ifdef DUMP
/* debug functions */
static int savetmpps(mesh2d *pm);
static int savebndps(mesh2d *pm);
static int savewrkps(mesh2d *pm);
/*static int saveintps(REAL *vertex, int a, int b, int c, int u, int v);*/
#endif
static int savefrtraw(mesh2d *pm);


/* Metric functions */
static REAL dist2(REAL x1, REAL y1, REAL x2, REAL y2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
};
static REAL dist2i(REAL *vertex, int v1, int v2) {
    return dist2(vertex[2*v1+0], vertex[2*v1+1], vertex[2*v2+0], vertex[2*v2+1]);
};

/* Color functions ******************************************************************************************************************************** COLOR */
static int pcolor(mesh2d *pm, int c) {
    int i;
    for (i=0; i<pm->ncolor; i++) {
	if (pm->color[2*i+0]==c) break;
    }
    if (i==pm->ncolor) {
	pm->ncolor++;
	pm->color[2*i+0] = c;
	pm->color[2*i+1] = c;
    }
    if (pm->color[2*i+1] != c) {
	pm->color[2*i+1] = pcolor(pm, pm->color[2*i+1]);
    }
    return pm->color[2*i+1];
}
#ifdef EXTRADEBUG
static void colordump(mesh2d *pm) {
    int i;
    for (i=0; i<pm->ncolor; i++) {
	printf("%d: %d\n", pm->color[2*i+0], pm->color[2*i+1]);
    }
}
#endif
static int cjoin(mesh2d *pm, int c1, int c2) {
    int i, p1=pcolor(pm, c1), p2=pcolor(pm, c2);
    if (p1 > p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p1) break;
	}
	pm->color[2*i+1] = p2;
	return pcolor(pm, c1);
    } else if (p1 < p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p2) break;
	}
	pm->color[2*i+1] = p1;
	return pcolor(pm, c2);
    } else return p1;
}

/* Quadtree functions ***************************************************************************************************************************** OCTREE */
#define TREEDEPTH 24
static void etreemetric(mesh2d *pm, int w1, int w2, REAL *px, REAL *pr) {
    int v1, v2;
    REAL r, s, x, y;
    
    if (w1 < w2)  v1 = w1,  v2 = w2;
    else  v1 = w2,  v2 = w1;
    x = px[0] = (pm->vertex[2*v1+0] + pm->vertex[2*v2+0])/2.0;
    y = px[1] = (pm->vertex[2*v1+1] + pm->vertex[2*v2+1])/2.0;
    r = 0.0;
    s = dist2(pm->vertex[2*v1+0], pm->vertex[2*v1+1], x, y);
    if (s > r) r = s;
    s = dist2(pm->vertex[2*v2+0], pm->vertex[2*v2+1], x, y);
    if (s > r) r = s;
    *pr = r;
}
static void emetric(mesh2d *pm, int i) {
    int v1, v2;
    v1 = pm->e[i].v[0];
    v2 = pm->e[i].v[1];
    etreemetric(pm, v1, v2, pm->e[i].x, &pm->e[i].r);
    pm->e[i].s = dist2i(pm->vertex, v1, v2);
}
static void findbb(mesh2d *pm) {
    int i;
    REAL d;
    if (pm->nV == 0) {
	pm->bb[0] = 0.0,  pm->bb[1] = 0.0;
	pm->bb[2] = 0.0,  pm->bb[3] = 0.0;
	return;
    }
    pm->bb[0] = pm->vertex[0],  pm->bb[1] = pm->vertex[1];
    pm->bb[2] = pm->vertex[0],  pm->bb[3] = pm->vertex[1];
    for (i=1; i<pm->nV; i++) {
	if (pm->bb[0] > pm->vertex[2*i+0])  pm->bb[0] = pm->vertex[2*i+0];
	if (pm->bb[1] > pm->vertex[2*i+1])  pm->bb[1] = pm->vertex[2*i+1];
	if (pm->bb[2] < pm->vertex[2*i+0])  pm->bb[2] = pm->vertex[2*i+0];
	if (pm->bb[3] < pm->vertex[2*i+1])  pm->bb[3] = pm->vertex[2*i+1];
    }
    d = 0.0;
    for (i=0; i<2; i++) {
	if (d < pm->bb[2+i]-pm->bb[i])  d = pm->bb[2+i]-pm->bb[i];
    }
    for (i=0; i<2; i++) {
	pm->bb[i] = (pm->bb[i] + pm->bb[2+i] - d)/2.0;
	pm->bb[2+i] = pm->bb[i] + d;
    }
}
static int treedir(REAL *bb, REAL *x) {
    int d = 0;
    if (x[0] > (bb[0]+bb[2])/2.0) d+=1;
    if (x[1] > (bb[1]+bb[3])/2.0) d+=2;
    return d;
}
static void treemove(REAL *bb, REAL s, int d) {
    if (d&1)  bb[0] += s;  else  bb[2] -= s;
    if (d&2)  bb[1] += s;  else  bb[3] -= s;
}
static void treeback(REAL *bb, REAL s, int d) {
    if (d&1)  bb[0] -= s;  else  bb[2] += s;
    if (d&2)  bb[1] -= s;  else  bb[3] += s;
}
static int listadd(mesh2d *pm, int *pn, int *pnn, int **pe, int n) {
    #ifdef STRUCTCHECK
    int i;
    int v1, v2, w1, w2;
    v1 = pm->e[n].v[0];
    v2 = pm->e[n].v[1];
    for (i=0; i<*pn; i++) {
	w1 = pm->e[(*pe)[i]].v[0];
	w2 = pm->e[(*pe)[i]].v[1];
	if ((v1!=w1)&&(v1!=w2)) continue;
	if ((v2!=w1)&&(v2!=w2)) continue;
	libaft_2d_warn("aft2d.c: listadd(): duplicate faces detected (%d %d)", v1, v2);
    }
    #endif
    (void) pm;
    if (*pn >= *pnn) {
	*pnn += 16;
	*pe = libaft_realloc(*pe, sizeof(int)* *pnn);
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
		libaft_free(*pe);
		*pe = 0;
	    }
	    return 0;
	}
    }
    libaft_2d_warn("aft2d.c: listdel(): element not found (internal error)");
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
    libaft_2d_warn("aft2d.c: listrep(): element not found (internal error)");
    return -1;
}
static node *nodenew(void) {
    node *p;
    int i;
    p = libaft_malloc(sizeof(node));
    p->ne = 0;
    p->nn = 0;
    p->e = 0;
    for (i=0; i<4; i++)  p->child[i] = 0;
    return p;
}
static int nodevoid(node *n) {
    int i;
    if (n->nn)  return 0;
    for (i=0; i<4; i++)  if (n->child[i])  return 0;
    if (n->e)  libaft_2d_warn("aft2d.c: nodevoid(): node void, but list is not NULL (sure memory leak)");
    return 1;
}
static node *treeadd(mesh2d *pm, int n) {
    node *c, *p;
    REAL bb[4], x[2], s, r;
    int i, d, level;
    
    p = pm->root;
    if (!p) {
	p = nodenew();
	pm->root = p;
    }
    for (i=0; i<4; i++)  bb[i] = pm->bb[i];
    for (i=0; i<2; i++)  x[i] = pm->e[n].x[i];
    r = pm->e[n].r;
    s = (pm->bb[2]-pm->bb[0])/2.0;
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
    listadd(pm, &c->ne, &c->nn, &c->e, n);
    return c;
}
static int treedel(mesh2d *pm, int n) {
    node *c, *p;
    node *stack[TREEDEPTH];
    int dirs[TREEDEPTH];
    REAL bb[4], x[2], s, r;
    int i, d, level;
    
    p = pm->root;
    if (!p)  {
	libaft_2d_warn("aft2d.c: treedel(): null root in quadtree (internal error)");
	return -1;
    }
    for (i=0; i<4; i++)  bb[i] = pm->bb[i];
    for (i=0; i<2; i++)  x[i] = pm->e[n].x[i];
    r = pm->e[n].r;
    s = (pm->bb[2]-pm->bb[0])/2.0;
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
	    libaft_2d_warn("aft2d.c: treedel(): null child in octree (internal error)");
	    return -1;
	}
	c = p;
    }
    if (listdel(&c->ne, &c->nn, &c->e, n))  return -1;
    while ((level>0) && (nodevoid(c))) {
	libaft_free(c);
	level--;
	c = stack[level];
	c->child[dirs[level]] = 0;
    }
    if ((level==0) && (nodevoid(c))) {
	libaft_free(c);
	pm->root = 0;
    }
    return 0;
}
static int wscheck(mesh2d *pm, int i, REAL *x, REAL r) {
    return (dist2(pm->e[i].x[0], pm->e[i].x[1], x[0], x[1]) <= pm->e[i].r + r) ? 1 : 0;
}
static int treewsrec(mesh2d *pm, node *c, REAL *bb, REAL s, REAL *x, REAL r) {
    int i, k = 0;
    if (!c)  return 0;
    for (i=0; i<2; i++) {
	if (bb[  i] - 2.0*s > x[i] + r)  return 0;
	if (bb[2+i] + 2.0*s < x[i] - r)  return 0;
    }
    for (i=0; i<c->ne; i++) {
	if (wscheck(pm, c->e[i], x, r))  pm->workarea[pm->nWorkArea++] = c->e[i];
    }
    k += c->ne;
    for (i=0; i<4; i++) {
	treemove(bb, s, i);
	k += treewsrec(pm, c->child[i], bb, s/2.0, x, r);
	treeback(bb, s, i);
    }
    return k;
}
static int treews(mesh2d *pm, REAL *x, REAL r) {
    node *p;
    REAL bb[4], s;
    int i;
    
    p = pm->root;
    if (!p)  {
	libaft_2d_warn("aft2d.c: treews(): null root in octree (internal error)");
	return -1;
    }
    for (i=0; i<4; i++)  bb[i] = pm->bb[i];
    s = (pm->bb[2]-pm->bb[0])/2.0;
    return treewsrec(pm, p, bb, s, x, r);
}
static int listfindedge(mesh2d *pm, int ne, int *e, int v1, int v2) {
    int i, w1, w2;
    for (i=0; i<ne; i++) {
	w1 = pm->e[e[i]].v[0];
	w2 = pm->e[e[i]].v[1];
	if ((v1==w1)&&(v2==w2)) return e[i];
	if ((v2==w1)&&(v1==w2)) return e[i];
    }
    return -1;
}
static int findedge(mesh2d *pm, int v1, int v2) {
    node *c, *p;
    REAL bb[4], x[2], s, r;
    int i, d, level;
    
    p = pm->root;
    if (!p)  {
	return -1;
    }
    for (i=0; i<4; i++)  bb[i] = pm->bb[i];
    etreemetric(pm, v1, v2, x, &r);
    s = (pm->bb[3]-pm->bb[0])/2.0;
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
    return listfindedge(pm, c->ne, c->e, v1, v2);
}


/* structure functions */
static int addtria(mesh2d *pm, int v1, int v2, int v3, int color) {
    color = pcolor(pm, color);
    if (pm->nT >= pm->nnT) {
	libaft_2d_stop("maxTria exceeded");
	return -1;
    } else {
	pm->tria[3*(pm->nT)+0] = v1;
	pm->tria[3*(pm->nT)+1] = v2;
	pm->tria[3*(pm->nT)+2] = v3;
	pm->triacolor[pm->nT] = color;
	return (pm->nT)++;
    }
}
static int addvertex(mesh2d *pm, REAL x, REAL y) {
    if (pm->nV >= pm->nnV) {
	libaft_2d_stop("maxVertex exceeded");
	return -1;
    } else {
	pm->vertex[2*(pm->nV)+0] = x;
	pm->vertex[2*(pm->nV)+1] = y;
	return (pm->nV)++;
    }
}
static void heapswap(mesh2d *pm, int k, int m) {
    int t;
    t = pm->heap[k],  pm->heap[k] = pm->heap[m],  pm->heap[m] = t;
    pm->e[pm->heap[k]].h = k;
    pm->e[pm->heap[m]].h = m;
}
static int siftup(mesh2d *pm, int n) {
    int i;
    while (1) {
	i = (n-1)/2;
	if ((i<0) || (pm->e[pm->heap[i]].s <= pm->e[pm->heap[n]].s)) break;
	heapswap(pm, i, n);
	n = i;
    }
    return n;
}
static int siftdown(mesh2d *pm, int n) {
    int i;
    n = siftup(pm, n);
    while (1) {
	i = 0;
	if ((2*n+1 < pm->nE) && (pm->e[pm->heap[n]].s > pm->e[pm->heap[(2*n+1)]].s)) {
	    if ((2*n+2 < pm->nE) && (pm->e[pm->heap[(2*n+1)]].s > pm->e[pm->heap[(2*n+2)]].s))  i = 2*n+2;
	    else  i = 2*n + 1;
	} else if ((2*n+2 < pm->nE) && (pm->e[pm->heap[n]].s > pm->e[pm->heap[(2*n+2)]].s))  i = 2*n+2;
	if (i) {
	    heapswap(pm, n, i);
	    n = i;
	} else break;
    }
    return n;
}
static int addedge(mesh2d *pm, int v1, int v2, int color, REAL fs) {
    #ifdef STRUCTCHECK
    int i;
    for (i=0; i<pm->nE; i++) {
	if ((v1!=pm->e[i].v[0])&&(v1!=pm->e[i].v[1])) continue;
	if ((v2!=pm->e[i].v[0])&&(v2!=pm->e[i].v[1])) continue;
	libaft_2d_warn("aft2d.c: addface(): duplicate faces detected (bad front?)");
    }
    #endif
    color = pcolor(pm, color);
    if (pm->nE >= pm->nnE) {
	libaft_2d_stop("AFT Front: maxEdge exceeded");
	return -1;
    } else {
	pm->e[(pm->nE)].v[0] = v1;
	pm->e[(pm->nE)].v[1] = v2;
	pm->e[pm->nE].c = color;
	emetric(pm, pm->nE);
	if (fs>0.0)  pm->e[pm->nE].s = fs;
	pm->heap[pm->nE] = pm->nE;
	pm->e[pm->nE].h = pm->nE;
	pm->e[pm->nE].node = treeadd(pm, pm->nE);
	return siftup(pm, (pm->nE)++);
    }
}
static void remedge(mesh2d *pm, int i, int color) {
    cjoin(pm, pm->e[i].c, color);
    if (pm->npack > pm->nnE) libaft_2d_stop("aft2d.c: Internal error in remface & packface");
    pm->pack[pm->npack] = i;
    pm->npack++;
}
static int packedge(mesh2d *pm) {
    int i, j;
    int err = 0;
    node *c;
    while (pm->npack > 0) {
	pm->npack--;
	j = pm->pack[pm->npack];
	pm->nE--;
	for (i=0; i<pm->npack; i++) if (pm->pack[i] == pm->nE) pm->pack[i] = j;
	c = pm->e[(pm->nE)].node;
	if (treedel(pm, j))  err++;
	if (listrep(&c->ne, &c->e, pm->nE, j))  err++;
	pm->e[j].v[0] = pm->e[(pm->nE)].v[0];
	pm->e[j].v[1] = pm->e[(pm->nE)].v[1];
	pm->e[j].c = pm->e[pm->nE].c;
	pm->e[j].x[0] = pm->e[(pm->nE)].x[0];
	pm->e[j].x[1] = pm->e[(pm->nE)].x[1];
	pm->e[j].r = pm->e[(pm->nE)].r;
	pm->e[j].s = pm->e[(pm->nE)].s;
	pm->e[j].node = c;
	i = pm->e[j].h;
	pm->nE++;
	siftdown(pm, i);
	pm->nE--;
	i = pm->e[pm->nE].h;
	pm->heap[i] = pm->heap[pm->nE];
	pm->e[pm->heap[i]].h = i;
	siftdown(pm, i);
    }
    return err;
}


static int VinWS(mesh2d *pm, int v) {
    int i;
    for (i=0; i<pm->nWorkSpace; i++) if (v==pm->workspace[i]) return 1;
    return 0;
}
static int VinBG(mesh2d *pm, int v) {
    int i;
    for (i=0; i<pm->nBadGuys; i++) if (v==pm->badguys[i]) return 1;
    return 0;
}
static int VinG(mesh2d *pm, int v) {
    int i;
    for (i=0; i<pm->nGirls; i++) if (v==pm->girls[i]) return 1;
    return 0;
}


/* quality section */
static REAL func_q(REAL *vertex, int k2, int k1, int k) {
    REAL L, S;
    S = (vertex[2*k2+0]*vertex[2*k1+1] - vertex[2*k1+0]*vertex[2*k2+1] + vertex[2*k+0]*(vertex[2*k2+1]-vertex[2*k1+1]) + vertex[2*k+1]*(vertex[2*k1+0]-vertex[2*k2+0]))/2.0;
    L = 2.0*vertex[2*k+0]*vertex[2*k+0] + 2.0*vertex[2*k1+0]*vertex[2*k1+0] + 2.0*vertex[2*k2+0]*vertex[2*k2+0] - 2.0*vertex[2*k+0]*vertex[2*k1+0] - 2.0*vertex[2*k+0]*vertex[2*k2+0] - 2.0*vertex[2*k1+0]*vertex[2*k2+0] +
        2.0*vertex[2*k+1]*vertex[2*k+1] + 2.0*vertex[2*k1+1]*vertex[2*k1+1] + 2.0*vertex[2*k2+1]*vertex[2*k2+1] - 2.0*vertex[2*k+1]*vertex[2*k1+1] - 2.0*vertex[2*k+1]*vertex[2*k2+1] - 2.0*vertex[2*k1+1]*vertex[2*k2+1];
    return  4.0*sqrt(3.0)*S/L;
}
static int func_xy(REAL *vertex, int k, int k2, int k1, REAL dx[2], REAL delta, int n) {
    REAL L, S, Lx, Ly, H, Hx, Hy, f, r;
    n = 8;
    S = (vertex[2*k2+0]*vertex[2*k1+1] - vertex[2*k1+0]*vertex[2*k2+1] + vertex[2*k+0]*(vertex[2*k2+1]-vertex[2*k1+1]) + vertex[2*k+1]*(vertex[2*k1+0]-vertex[2*k2+0]))/2.0;
    H = S + sqrt(S*S + 4.0*delta*delta);
    //    H = 2.0*S;
    L = 2.0*vertex[2*k+0]*vertex[2*k+0] + 2.0*vertex[2*k1+0]*vertex[2*k1+0] + 2.0*vertex[2*k2+0]*vertex[2*k2+0] - 2.0*vertex[2*k+0]*vertex[2*k1+0] - 2.0*vertex[2*k+0]*vertex[2*k2+0] - 2.0*vertex[2*k1+0]*vertex[2*k2+0] +
    2.0*vertex[2*k+1]*vertex[2*k+1] + 2.0*vertex[2*k1+1]*vertex[2*k1+1] + 2.0*vertex[2*k2+1]*vertex[2*k2+1] - 2.0*vertex[2*k+1]*vertex[2*k1+1] - 2.0*vertex[2*k+1]*vertex[2*k2+1] - 2.0*vertex[2*k1+1]*vertex[2*k2+1];
    Lx = 4.0*vertex[2*k+0] - 2.0*(vertex[2*k1+0]+vertex[2*k2+0]);
    Ly = 4.0*vertex[2*k+1] - 2.0*(vertex[2*k1+1]+vertex[2*k2+1]);
    Hx = (vertex[2*k2+1]-vertex[2*k1+1])/2.0 + (S*(vertex[2*k2+1]-vertex[2*k1+1]))/2.0/sqrt(S*S + 4.0*delta*delta);
    Hy = (vertex[2*k1+0]-vertex[2*k2+0])/2.0 + (S*(vertex[2*k1+0]-vertex[2*k2+0]))/2.0/sqrt(S*S + 4.0*delta*delta);
    //    Hx = (vertex[2*k2+1]-vertex[2*k1+1]);
    //    Hy = (vertex[2*k1+0]-vertex[2*k2+0]);
    f = 2.0*L/H/4.0/sqrt(3.0);
    dx[0] = (2.0*Lx*H - 2.0*Hx*L)/H/H/4.0/sqrt(3.0);
    dx[1] = (2.0*Ly*H - 2.0*Hy*L)/H/H/4.0/sqrt(3.0);
    while (n>1) {
	dx[0] *= f;
	dx[1] *= f;
	n--;
    }
    /*    r = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    if (r>1.0) {
   r = (1.0 + log(r))/r;
   dx[0] *= r,  dx[1] *= r;
}*/
    return 0;
    (void) r;
}

typedef struct {
    int p;
    int n;
} plist;

static int add_glist(int a, int b, int *png, plist *glist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (glist[c].p == b) return 0;
	c = glist[c].n;
    }
    glist[*png].p = b;
    glist[*png].n = s[a];
    s[a] = *png;
    (*png)++;
    return 1;
}

static int fill_eadj(mesh2d *pm) {
    int *s;
    int ng;
    plist *glist;
    int i, j, n, l;
    int a, b, c;
    
    glist = (plist*)libaft_malloc(sizeof(plist)*6*pm->nT);
    ng = 0;
    s = (int*)libaft_malloc(sizeof(int)*pm->nV);
    for (j=0; j<pm->nV; j++)  s[j] = -1;
    
    for (i=0; i<pm->nT; i++) {
	a = pm->tria[3*i+0],  b = pm->tria[3*i+1],  c = pm->tria[3*i+2];
	add_glist(a, b, &ng, glist, s);
	add_glist(a, c, &ng, glist, s);
	add_glist(b, c, &ng, glist, s);
	add_glist(b, a, &ng, glist, s);
	add_glist(c, a, &ng, glist, s);
	add_glist(c, b, &ng, glist, s);
    }
    n = 0;
    for (j=0; j<pm->nV; j++) {
	pm->eadj.ia[j] = n;
	l = s[j];
	while (l>=0) {
	    pm->eadj.ja[n++] = glist[l].p;
	    l = glist[l].n;
	}
    }
    pm->eadj.ia[pm->nV] = n;
    
    libaft_free(s),  libaft_free(glist);
    return n;
}


typedef struct {
    tedge e;
    int n;
} elist;

static int add_hlist(int a, tedge b, int *pnh, elist *hlist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (hlist[c].e[0] == b[0] && hlist[c].e[1] == b[1]) return 0;
	c = hlist[c].n;
    }
    hlist[*pnh].e[0] = b[0];
    hlist[*pnh].e[1] = b[1];
    hlist[*pnh].n = s[a];
    s[a] = *pnh;
    (*pnh)++;
    return 1;
}
static int fill_tadj(mesh2d *pm) {
    int *s;
    int nh;
    elist *hlist;
    int i, j, n, l;
    int a;
    tedge b;
    
    hlist = (elist*)libaft_malloc(sizeof(elist)*3*pm->nT);
    nh = 0;
    s = (int*)libaft_malloc(sizeof(int)*pm->nV);
    for (j=0; j<pm->nV; j++)  s[j] = -1;
    
    for (i=0; i<pm->nT; i++) {
	a = pm->tria[3*i+0],  b[0] = pm->tria[3*i+1],  b[1] = pm->tria[3*i+2],  add_hlist(a, b, &nh, hlist, s);
	a = pm->tria[3*i+1],  b[0] = pm->tria[3*i+2],  b[1] = pm->tria[3*i+0],  add_hlist(a, b, &nh, hlist, s);
	a = pm->tria[3*i+2],  b[0] = pm->tria[3*i+0],  b[1] = pm->tria[3*i+1],  add_hlist(a, b, &nh, hlist, s);
    }
    n = 0;
    for (j=0; j<pm->nV; j++) {
	pm->tadj.ia[j] = n;
	l = s[j];
	while (l>=0) {
	    pm->tadj.ja[n][0] = hlist[l].e[0];
	    pm->tadj.ja[n][1] = hlist[l].e[1];
	    n++;
	    l = hlist[l].n;
	}
    }
    pm->tadj.ia[pm->nV] = n;
    
    libaft_free(s),  libaft_free(hlist);
    return n;
}

static int opt_func(mesh2d *pm, int nfixed) {
    REAL d, delta, r, ds, x[2], z[2], rs, as;
    int j, l, s;
    int p, b, c;
    int pn, n_iters;
    int *ps;
    tvertex *dx, *bkp;
    
    
    ps = (int*)libaft_malloc(sizeof(int)*(pm->nV - nfixed));
    pn = 0;
    ds = 0.0;
    for (j=nfixed; j<pm->nV; j++) {
	p = j;
	ps[pn++] = p;
	d = 0.0,  s = 0;
	for (l=pm->eadj.ia[p]; l<pm->eadj.ia[p+1]; l++) {
	    c = pm->eadj.ja[l];
	    d += dist2i(pm->vertex, p, c);
	    s++;
	}
	if (s > 0)  d /= s;
	ds += d;
    }
    if (pn == 0)  {
	libaft_free(ps);
	return 0;
    }
    
    ds /= pn;
    //printf("ds=%lf\n", ds);
    
    dx = (tvertex*)libaft_malloc(sizeof(tvertex)*pn);
    bkp = (tvertex*)libaft_malloc(sizeof(tvertex)*pn);
    for (j=0; j<pn; j++) {
	p = ps[j];
	bkp[j][0] = pm->vertex[2*p+0];
	bkp[j][1] = pm->vertex[2*p+1];
    }
    
    delta = 1.0;
    n_iters = 100; // 4000
    d = 0.25 * ds * ds / n_iters;  // 1.0
    while (1) {
	for (s=0; s<n_iters; s++) {
	    delta = 0.01 * ds * ds * (1.0 - 0.9*s/n_iters);
	    
	    //	delta = -as / n_iters;
	    rs = 0.0;
	    for (j=0; j<pn; j++) {
		p = ps[j];
		x[0] = 0.0,  x[1] = 0.0;
		for (l=pm->tadj.ia[p]; l<pm->tadj.ia[p+1]; l++) {
		    b = pm->tadj.ja[l][0];
		    c = pm->tadj.ja[l][1];
		    func_xy(pm->vertex, p, b, c, z, delta, 8);
		    x[0] -= z[0],  x[1] -= z[1];
		}
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
		if (r>1.0/ds) {
		    r = 1.0/ds*(1.0 + log(r*ds))/r;
		    x[0] *= r,  x[1] *= r;
		}
		dx[j][0] = x[0],  dx[j][1] = x[1];
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
		if (rs < r)  rs = r;
	    }
	    for (j=0; j<pn; j++) {
		p = ps[j];
		pm->vertex[2*p+0] += d*dx[j][0];
		pm->vertex[2*p+1] += d*dx[j][1];
	    }
	    //printf("%3d: delta = %12.8lf, d = %12.8lf\n", s, delta, d*rs/ds);
	    //savedump();
	    //scanf("%d", &j);
	    //	if (d*rs/ds < 0.01)  break;
#ifdef SHOWPROGRESS
	    printf("\t\t\t\t\t\t %3d%%\r", s+1),  fflush(stdout);
#endif
	}
#ifdef SHOWPROGRESS
	printf("\t\t\t\t\t\t      \r"),  fflush(stdout);
#endif
	break;
    }
    s = 0;
    for (j=0; j<pm->nT; j++) {
	if (func_q(pm->vertex, pm->tria[3*j+0], pm->tria[3*j+1], pm->tria[3*j+2]) <= 0.0) s++;
    }
    if (s) {
	libaft_2d_warn("Quality improvement failed, falling back to simple smoothing");
	for (j=0; j<pn; j++) {
	    p = ps[j];
	    pm->vertex[2*p+0] = bkp[j][0];
	    pm->vertex[2*p+1] = bkp[j][1];
	}
    }
    libaft_free(ps);
    libaft_free(dx);
    libaft_free(bkp);
    return s;
    (void) as;
}


/* intersection section */
static int intsect(REAL *vertex, int a, int b, int c, int u, int v) {
    int uv, dup=0;

    if ((u==a)||(u==b)||(u==c)) dup++;
    if ((v==a)||(v==b)||(v==c)) dup++;
    if (dup==0) {
	uv = idet2i3(vertex, u, v, a) + idet2i3(vertex, u, v, b) + idet2i3(vertex, u, v, c);
	if ((uv==3) || (uv==-3)) return 0;
	if (idet2i3(vertex, b, c, u) + idet2i3(vertex, b, c, v) == -2) return 0;
	if (idet2i3(vertex, c, a, u) + idet2i3(vertex, c, a, v) == -2) return 0;
	if (idet2i3(vertex, a, b, u) + idet2i3(vertex, a, b, v) == -2) return 0;
    } else if (dup==1) {
	if (idet2i3(vertex, b, c, u) + idet2i3(vertex, b, c, v) == -1) return 0;
	if (idet2i3(vertex, c, a, u) + idet2i3(vertex, c, a, v) == -1) return 0;
	if (idet2i3(vertex, a, b, u) + idet2i3(vertex, a, b, v) == -1) return 0;
    } else return 0;
    return 1;
}
static int intedge(mesh2d *pm, int a, int b, int u, int v) {
    REAL x, y, xc, yc;
    int c, i1, i2;
    x = pm->vertex[2*a+1] - pm->vertex[2*b+1], xc = (pm->vertex[2*a+0] + pm->vertex[2*b+0])/2.0;
    y = pm->vertex[2*a+0] - pm->vertex[2*b+0], yc = (pm->vertex[2*a+1] + pm->vertex[2*b+1])/2.0;
    c = addvertex(pm, xc+x, yc-y);
    if (c < 0)  return 1;
    i1 = intsect(pm->vertex, a, b, c, u, v);
    pm->nV--;
    c = addvertex(pm, xc-x, yc+y);
    if (c < 0)  return 1;
    i2 = intsect(pm->vertex, b, a, c, u, v);
    pm->nV--;
    return (i1 && i2);
}
static int intfront(mesh2d *pm) {
    int i, j, k;
    for (i=0; i<pm->nE; i++) {
	pm->nWorkArea = 0;
	if (treews(pm, pm->e[i].x, pm->e[i].r) < 0)  return -4;
	for (k=0; k<pm->nWorkArea; k++) {
	    j = pm->workarea[k];
	    if (i==j)  continue;
	    if (intedge(pm, pm->e[i].v[0], pm->e[i].v[1], pm->e[j].v[0], pm->e[j].v[1])) {
		return 1;
	    }
	}
    }
    return 0;
}


/* AFT section */
static int check(mesh2d *pm, int en, int pn) {
    int v1, v2, i, p1, p2;

    pm->intedge = -1;
    v1 = pm->e[en ].v[0];
    v2 = pm->e[en ].v[1];
    if (idet2i3(pm->vertex, v1, v2, pn) != 1) {
	libaft_2d_warn("inverted? %d (%d %d %d)", idet2i3(pm->vertex, v1, v2, pn), v1, v2, pn);
	return 1;
    }
    for (i=0; i<pm->nWorkArea; i++) {
	p1 = pm->e[pm->workarea[i]].v[0];
	p2 = pm->e[pm->workarea[i]].v[1];
	if (intsect(pm->vertex, v1, v2, pn, p1, p2)) {
	    pm->intedge = pm->workarea[i];
	    return 1;
	}
    }
    return 0;
}

static REAL height(REAL *vertex, int v0, int v1, int v) {
    REAL x, y, x0, y0, x1, y1, s, r, c;
    x = vertex[2*v+0],  y = vertex[2*v+1];
    x0 = vertex[2*v0+0],  y0 = vertex[2*v0+1];
    x1 = vertex[2*v1+0],  y1 = vertex[2*v1+1];
    c = ((x-x0)*(x1-x0) + (y-y0)*(y1-y0)) / ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    if (c<=0.0)  return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    else if (c>=1.0)  return sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));
    s = (x0-x)*(y1-y) - (y0-y)*(x1-x);
    r = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    return fabs(s)/r;
}


static int findpoint(mesh2d *pm, int en) {
    int v1, v2, vc, pn, p1, p2;
    REAL ex, ey, r, rp, x, y, rmin, rv, friendness, gx, gy;
    int i, j, k, m, v, neari, dirty;
    REAL hc, h, q;
    REAL sx[2];

    v1 = pm->e[en].v[0];
    v2 = pm->e[en].v[1];

    ex = pm->vertex[2*v1+1] - pm->vertex[2*v2+1];
    ey = pm->vertex[2*v2+0] - pm->vertex[2*v1+0];
    r = sqrt(ex*ex + ey*ey);
    if (r <= 0.0) libaft_2d_warn("zero sized edge %d [%d %d]", en, v1, v2);
    ex /= r, ey /= r;
    rp = r/2.0;
    if (pm->stepsize > 0.0) {
	r = ((pm->stepsize>=1.0) ? pm->stepsize : 1.0) * pm->e[en].s;
	if ((pm->sizelim > 0.0) && (r > pm->sizelim))  r = pm->sizelim;
    } else {
	x = (pm->vertex[2*v1+0] + pm->vertex[2*v2+0])/2.0;
	y = (pm->vertex[2*v1+1] + pm->vertex[2*v2+1])/2.0;
	rmin = pm->fsize(x, y) * 0.25 * sqrt(3.0/4.0);
	r = pm->fsize(x + rmin*ex, y + rmin*ey) * sqrt(3.0/4.0);
    }
    if (r<1.1*rp)  r = 1.1*rp;
    x = (pm->vertex[2*v1+0] + pm->vertex[2*v2+0])/2.0 + ex*sqrt(r*r-rp*rp);
    y = (pm->vertex[2*v1+1] + pm->vertex[2*v2+1])/2.0 + ey*sqrt(r*r-rp*rp);
    vc = addvertex(pm, x, y);
    if (vc < 0)  return -5;
    if (idet2i3(pm->vertex, v1, v2, vc) != 1) libaft_2d_warn("wrong edge-orientation of new vertex");
    r = 0.0;
    r += dist2i(pm->vertex, v1, vc);
    r *= 1.0001220703125;
    hc = height(pm->vertex, v1, v2, vc),  dirty = 0;
    
    if (pm->local) {
	pm->nWorkArea = 0;
	sx[0] = pm->vertex[2*vc+0],  sx[1] = pm->vertex[2*vc+1];
	i = treews(pm, sx, 2.0*r);
	if (i < 0)  return -4;
//	rmin = 0.75*r + 1.0*pm->stepsize*pm->minrho, neari = -1, pm->nWorkSpace = 0;
	rmin = (pm->beta*r + pm->minrho + sqrt((pm->beta*r - pm->minrho)*(pm->beta*r - pm->minrho) + pm->alpha))/2.0, neari = -1, pm->nWorkSpace = 0;
	for (i=0; i<pm->nWorkArea; i++) {
	    h = height(pm->vertex, pm->e[pm->workarea[i]].v[0], pm->e[pm->workarea[i]].v[1], vc);
	    if (h < 0.5*hc)  dirty++;//,  printf("dirty: %lf, %lf\n", hc, h);
	    for (j=0; j<2; j++) {
		v = pm->e[pm->workarea[i]].v[j];
		if (VinWS(pm, v)) continue;
		rv = dist2i(pm->vertex, vc, v);
		if (rv > 2.0*r) continue;
		if (idet2i3(pm->vertex, v1, v2, v) != 1) continue;
		q = func_q(pm->vertex, v1, v2, v);
		for (k=0; k<pm->nWorkSpace; k++) {
		    if (q > pm->workspace_q[k])  break;
		}
		if (k==pm->nWorkSpace) {
		    pm->workspace_q[pm->nWorkSpace] = q,  pm->workspace[pm->nWorkSpace++] = v;
		} else {
		    for (m=pm->nWorkSpace; m>k; m--)  pm->workspace[m] = pm->workspace[m-1],  pm->workspace_q[m] = pm->workspace_q[m-1];
		    pm->workspace[k] = v,  pm->workspace_q[k] = q,  pm->nWorkSpace++;
		}
		if (rv < rmin) {
		    neari = v;
		    rmin = rv;
		}
	    }
	}
    } else {
	pm->nWorkArea = 0;
	for (i=0; i<pm->nE; i++)
	    pm->workarea[pm->nWorkArea++] = i;
	neari = -1, pm->nWorkSpace = 0;
	for (i=0; i<pm->nWorkArea; i++) {
	    h = height(pm->vertex, pm->e[pm->workarea[i]].v[0], pm->e[pm->workarea[i]].v[1], vc);
	    if (h < 0.5*hc)  dirty++;//,  printf("dirty: %lf, %lf\n", hc, h);
	    for (j=0; j<2; j++) {
		v = pm->e[pm->workarea[i]].v[j];
		if (VinWS(pm, v)) continue;
		if (idet2i3(pm->vertex, v1, v2, v) != 1) continue;
		q = func_q(pm->vertex, v1, v2, v);
		for (k=0; k<pm->nWorkSpace; k++) {
		    if (q > pm->workspace_q[k])  break;
		}
		if (k==pm->nWorkSpace) {
		    pm->workspace_q[pm->nWorkSpace] = q,  pm->workspace[pm->nWorkSpace++] = v;
		} else {
		    for (m=pm->nWorkSpace; m>k; m--)  pm->workspace[m] = pm->workspace[m-1],  pm->workspace_q[m] = pm->workspace_q[m-1];
		    pm->workspace[k] = v,  pm->workspace_q[k] = q,  pm->nWorkSpace++;
		}
//		pm->workspace[pm->nWorkSpace++] = v;
	    }
	}

    }
    if (neari >= 0)  pn = neari;
    else if (dirty && (pm->nWorkSpace > 0))  pn = pm->workspace[0];
    else  pn = vc;
    /*	savewrkps(pm);*/
    /* friends */
    pm->curfriend[0] = -1; pm->curfriendness[0] = -1; pm->curfriendlink[0] = -1;
    pm->curfriend[1] = -1; pm->curfriendness[1] = -1; pm->curfriendlink[1] = -1;
    for (i=0; i<pm->nWorkArea; i++) {
	p1 = pm->e[pm->workarea[i]].v[0];
	p2 = pm->e[pm->workarea[i]].v[1];
	if ( (p2 == v1) && (idet2i3(pm->vertex, v1, v2, p1) == 1) ) {
	    gx = pm->vertex[2*p2+1] - pm->vertex[2*p1+1];
	    gy = pm->vertex[2*p1+0] - pm->vertex[2*p2+0];
	    r = sqrt(gx*gx + gy*gy);
	    gx /= r, gy /= r;
	    friendness = ex*gx + ey*gy;
	    if ((pm->curfriend[0]<0) || (pm->curfriendness[0]<friendness)) {
		pm->curfriend[0] = p1;
		pm->curfriendness[0] = friendness;
		pm->curfriendlink[0] = pm->workarea[i];
	    }
	}
	if ( (p1 == v2) && (idet2i3(pm->vertex, v1, v2, p2) == 1) ) {
	    gx = pm->vertex[2*p2+1] - pm->vertex[2*p1+1];
	    gy = pm->vertex[2*p1+0] - pm->vertex[2*p2+0];
	    r = sqrt(gx*gx + gy*gy);
	    gx /= r, gy /= r;
	    friendness = ex*gx + ey*gy;
	    if ((pm->curfriend[1]<0) || (pm->curfriendness[1]<friendness)) {
		pm->curfriend[1] = p2;
		pm->curfriendness[1] = friendness;
		pm->curfriendlink[1] = pm->workarea[i];
	    }
	}
    }
    /* searching... */
    pm->nBadGuys = 0, pm->nGirls = 0;
    while (check(pm, en, pn)) {
	pm->badguys[pm->nBadGuys++] = pn;
	if (pm->intedge >= 0) {
	    for (i=0; i<2; i++) {
		v = pm->e[pm->intedge].v[i];
		if (v==v1) continue;
		if (v==v2) continue;
		if (!VinWS(pm, v)) continue;
		if (VinBG(pm, v)) continue;
		if (VinG(pm, v)) continue;
		pm->girls[pm->nGirls++] = v;
	    }
	    if (pm->nGirls>0) {
		neari = 0; rmin = dist2i(pm->vertex, vc, pm->girls[0]);
		for (i=1; i<pm->nGirls; i++) {
		    rv = dist2i(pm->vertex, vc, pm->girls[i]);
		    if (rmin > rv) {
			neari = i;
			rmin = rv;
		    }
		}
		pn = pm->girls[neari];
		pm->nGirls--;
		pm->girls[neari] = pm->girls[pm->nGirls];
	    } else {
		if (pm->local) {
		    pm->nV--;
		    pm->local = 0;
		    vc = findpoint(pm, en);
		    pm->local = 1;
		    return vc;
		} else {
		    libaft_2d_warn("Using brute force");
		    for (i=0; i<pm->nWorkSpace; i++) {
			v = pm->workspace[i];
			if (v==v1) continue;
			if (v==v2) continue;
			if (VinBG(pm, v)) continue;
			if (VinG(pm, v)) continue;
			pm->girls[pm->nGirls++] = v;
		    }
		    if (pm->nGirls>0) {
			neari = 0; rmin = dist2i(pm->vertex, vc, pm->girls[0]);
			for (i=1; i<pm->nGirls; i++) {
			    rv = dist2i(pm->vertex, vc, pm->girls[i]);
			    if (rmin > rv) {
				neari = i;
				rmin = rv;
			    }
			}
			pn = pm->girls[neari];
			pm->nGirls--;
			pm->girls[neari] = pm->girls[pm->nGirls];
		    } else {
			if (intfront(pm)) libaft_2d_warn("Intersection occured in process");
			libaft_2d_warn("No more points left after brute force");
			return -2;
		    }
		}
	    }
	} else {
	    libaft_2d_warn("Impossible situation in findpoint, check logic ;-)");
	    return -1;
	}
    }
    if (pn!=vc) {
	pm->nV--;
    }
    return pn;
}

static int front(mesh2d *pm) {
    int en, pn, v1, v2;
    REAL ss;
    while (pm->nE > 0) {
	/* Reserve elements for one step */
	if (pm->nT >= pm->nnT-1) { /* one tria */
	    libaft_2d_stop("maxTria exceeded");
	    return 2; /* maxTria */
	}
	if (pm->nV >= pm->nnV-2) { /* one vertex + one for intersection test */
	    libaft_2d_stop("maxVertex exceeded");
	    return 3; /* maxVertex */
	}
	if (pm->nE >= pm->nnE-2) { /* two edges */
	    libaft_2d_stop("maxEdge exceeded");
	    return 4; /* maxEdge */
	}
	en = pm->heap[0];
	ss = pm->e[en].s;
	if (pm->stepsize>=1.0)  ss *= pm->stepsize;
	else ss = 0.0;
	if ((pm->sizelim > 0.0) && (ss > pm->sizelim))  ss = 0.0;
	pn = findpoint(pm, en);
	if (pn<0) {
#ifdef DUMP
	    savetmpps(pm);
	    savebndps(pm);
	    savewrkps(pm);
#endif
	    if (pn == -5)  return 3; // maxVertex
	    return pn;
	}
	v1 = pm->e[en].v[0];
	v2 = pm->e[en].v[1];
	if (addtria(pm, v1, v2, pn, pm->e[en].c) < 0)  return 2; // maxtria
	if (idet2i3(pm->vertex, v1, v2, pn) != 1) libaft_2d_warn("Wrong orientation of tria");
	pm->meshedvolume += det2i3(pm->vertex, v1, v2, pn);
	remedge(pm, en, pm->e[en].c);
	if (pm->curfriend[0] == pn) {
	    remedge(pm, pm->curfriendlink[0], pm->e[en].c);
	} else {
	    if (addedge(pm, v1, pn, pm->e[en].c, ss) < 0)  return 4; // maxedge
	}
	if (pm->curfriend[1] == pn) {
	    remedge(pm, pm->curfriendlink[1], pm->e[en].c);
	} else {
	    if (addedge(pm, pn, v2, pm->e[en].c, ss) < 0)  return 4; // maxedge
	}
	if (packedge(pm))  return -3;
#ifdef SHOWPROGRESS
	printf("nV = %5d, nE = %5d, nT = %5d, %6.2lf%%\r", pm->nV, pm->nE, pm->nT, 100.0*pm->meshedvolume/pm->volume);
	fflush(stdout);
#endif
    }
#ifdef SHOWPROGRESS
    printf("nV = %5d, nE = %5d, nT = %5d,  done. \r", pm->nV, pm->nE, pm->nT);
    fflush(stdout);
#endif
    return 0;
}

/* main function */
/* return codes: 
 *  0 -- OK
 *  1 -- possible intersection in input data
 *  2 -- max_nT exceeded
 *  3 -- max_nV exceeded
 *  4 -- max_nF exceeded
 * -1 -- internal error (logic fail)
 * -2 -- internal error (bad front)
 * -3 -- internal error (data structure fail)
 * -4 -- internal error (data structure fail)
 */
static int mesh2daftss_main(int *pnV, REAL *vertex,
	int *pnE, int *edge, int *edgecolor,
	int *pnT, int *tria, int *triacolor,
	int nnV, int nnE, int nnT,
	REAL ss, REAL sizelim, double (*f)(double, double),
	int *pnC, int *color)
{
    mesh2d mesh;
    int i, r, nVF;
    int nE = *pnE;
    int treeerr = 0;
#ifdef PROF
    time_t tr1, tr2, tr;
#endif

    detinit();

    mesh.nV = *pnV, mesh.nE = 0, mesh.nT = *pnT;
    mesh.nnV = nnV, mesh.nnE = nnE, mesh.nnT = nnT;
    mesh.vertex = vertex, mesh.tria = tria, mesh.triacolor = triacolor;

    mesh.e = libaft_malloc(sizeof(edge2d) * mesh.nnE);
    mesh.heap       = libaft_malloc(sizeof(int) * mesh.nnE);
    mesh.workarea = libaft_malloc(sizeof(int) * mesh.nnE);
    mesh.workspace = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.workspace_q = libaft_malloc(sizeof(REAL) * mesh.nnV);
    mesh.badguys = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.girls = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.pack = libaft_malloc(sizeof(int) * mesh.nnE);
    mesh.color = libaft_malloc(sizeof(int) * nE*2);
    
    mesh.root = 0;
    
    mesh.meshedvolume = 0.0;
    mesh.stepsize = ss;
    mesh.fsize = f;
    mesh.sizelim = sizelim;
    mesh.npack = 0;
    mesh.local = 1;
    mesh.ncolor = 0;
    
#ifdef PROF
    tr1 = clock();
#endif

#ifdef SHOWPROGRESS
    printf("aft warmup...\r");
    fflush(stdout);
#endif
    
    findbb(&mesh);
    
    mesh.maxlength = 0.0;
    for (i=0; i<nE; i++)  {
	if (addedge(&mesh, edge[2*i+0], edge[2*i+1], edgecolor[i], 0.0) < 0)  return 4;
	if (mesh.maxlength < mesh.e[i].s)  mesh.maxlength = mesh.e[i].s;
    }
    mesh.maxlength *= 1e+6;
    
    // Here we compute signed volume of the domain
    // Front should be oriented and closed
    // This formula could be found at http://www.cgafaq.info/wiki/Polyhedron_Volume
    for (i=0, mesh.volume = 0.0; i<mesh.nE; i++)
	mesh.volume += det2i3(mesh.vertex, 0, mesh.e[i].v[0], mesh.e[i].v[1]);
    
    
#ifdef SHOWPROGRESS
    printf("nV = %5d, nE = %5d, nT = %5d, init...\r", mesh.nV, mesh.nE, mesh.nT);
    fflush(stdout);
#endif

    for (i=0, mesh.volume = 0.0; i<mesh.nE; i++)
	mesh.volume += det2i3(mesh.vertex, 0, mesh.e[i].v[0], mesh.e[i].v[1]);

    mesh.minrho = dist2i(mesh.vertex, mesh.e[0].v[0], mesh.e[0].v[1]);

    for (i=0; i<mesh.nE; i++) {
	if (mesh.minrho > dist2i(mesh.vertex, mesh.e[i].v[0], mesh.e[i].v[1]))
	    mesh.minrho = dist2i(mesh.vertex, mesh.e[i].v[0], mesh.e[i].v[1]);
    }

    if (mesh.stepsize > 0.0)  mesh.minrho *= 0.95;
    else mesh.minrho = 0.0;
    mesh.beta = 0.5;
    mesh.minrho *= mesh.beta;
    mesh.alpha = 2.0*mesh.minrho*(1.0-mesh.beta)/mesh.beta;
    mesh.alpha *= mesh.alpha;

    nVF = mesh.nV;

    if (intfront(&mesh)) {
#ifdef PROF
    tr2 = clock();
    tr = tr2-tr1;
#endif
	r = 1;
	libaft_2d_warn("Intersection in initial boundary");
    } else {
#ifdef PROF
    tr2 = clock();
    tr = tr2-tr1;
#ifdef SHOWPROGRESS
    printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr2-tr1)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#endif
#endif
	r = front(&mesh);
#ifdef PROF
	tr1 = clock();
	tr += tr1-tr2;
#ifdef SHOWPROGRESS
	printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr1-tr2)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#else
	printf("AFT kernel time: %7.2lfs\n", ((double)tr1-tr2)/CLOCKS_PER_SEC);
	printf("nT = %d\n", mesh.nT);
	printf("time per tria = %.3le\n", ((double)tr1-tr2)/CLOCKS_PER_SEC/mesh.nT);
	printf("tria per sec = %.0lf\n", mesh.nT / (((double)tr1-tr2)/CLOCKS_PER_SEC));
#endif
#endif
    }

    if (0) savefrtraw(&mesh);

    for (i=0; i<mesh.nT; i++) {
	mesh.triacolor[i] = pcolor(&mesh, mesh.triacolor[i]);
    }
    if (pnC && color) {
	*pnC = mesh.ncolor;
	for (i=0; i<mesh.ncolor; i++)  color[2*i+0] = mesh.color[2*i+0],  color[2*i+1] = mesh.color[2*i+1];
    }
    
    for (i=0; i<mesh.nE; i++)  edge[2*i+0] = mesh.e[i].v[0],  edge[2*i+1] = mesh.e[i].v[1],  edgecolor[i]=mesh.e[i].c;
    for (i=0; i<mesh.nE; i++)  if (treedel(&mesh, i))  treeerr++;
    
    libaft_free(mesh.e);
    libaft_free(mesh.heap);
    libaft_free(mesh.workarea);
    libaft_free(mesh.workspace);
    libaft_free(mesh.workspace_q);
    libaft_free(mesh.badguys);
    libaft_free(mesh.girls);
    libaft_free(mesh.pack);
    libaft_free(mesh.color);
    if (r == 0) {
#ifdef SHOWPROGRESS
	printf("nV = %5d, nE = %5d, nT = %5d, smoothing \r", mesh.nV, mesh.nE, mesh.nT);
	fflush(stdout);
#endif
	mesh.eadj.ia   = (int* )libaft_malloc(sizeof(int )*(mesh.nV + 1));
	mesh.eadj.ja   = (int* )libaft_malloc(sizeof(int )*(6*mesh.nT));
	mesh.tadj.ia   = (int* )libaft_malloc(sizeof(int )*(mesh.nV + 1));
	mesh.tadj.ja   = (tedge*)libaft_malloc(sizeof(tedge)*(3*mesh.nT));
	fill_eadj(&mesh);
	fill_tadj(&mesh);
	if (opt_func(&mesh, nVF))  refine2daft(&mesh.nV, mesh.vertex, &nVF, &mesh.nT, mesh.tria);

#ifdef PROF
    tr1 = clock();
    tr += tr1-tr2;
#ifdef SHOWPROGRESS
    printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr1-tr2)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#else
    printf("smoothing time: %7.2lfs\n", ((double)tr1-tr2)/CLOCKS_PER_SEC);
#endif
#endif

	libaft_free(mesh.eadj.ia),  libaft_free(mesh.eadj.ja),  libaft_free(mesh.tadj.ia),  libaft_free(mesh.tadj.ja);
    }
#ifdef SHOWPROGRESS
    printf("                                              \r");
    fflush(stdout);
#endif

    *pnV = mesh.nV, *pnE = mesh.nE, *pnT = mesh.nT;
    return r;
    (void)findedge;
}
int mesh2daftss(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss) {
    return mesh2daftss_main(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, ss, -1.0, NULL, NULL, NULL);
}
int mesh2daftsslim(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss, REAL lim) {
    return mesh2daftss_main(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, ss, lim, NULL, NULL, NULL);
}
int mesh2daft(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT) {
    return mesh2daftss_main(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, 1.5, -1.0, NULL, NULL, NULL);
}
int mesh2daftfull(int *pnV, REAL *vertex, int *pnE, int *edge, int *edgecolor, int *pnT, int *tria, int *triacolor, int nnV, int nnE, int nnT, REAL ss, REAL lim, int *pnC, int *color) {
    return mesh2daftss_main(pnV, vertex, pnE, edge, edgecolor, pnT, tria, triacolor, nnV, nnE, nnT, ss, lim, NULL, pnC, color);
}












#ifdef DUMP
/* debug section */

#define zoom 1.0
#define dotw 1.0

static int savetmpps(mesh2d *pm) {
    FILE *f;
    int i;

    f=fopen("tmp.ps", "w");
    if (!f) return 1;
    fprintf(f, "10 10 translate 0 setlinewidth %lf %lf scale\n", zoom, zoom);
    fprintf(f, "/t{newpath moveto lineto lineto closepath stroke}def\n");
    fprintf(f, "/d{newpath %lf 0 360 arc closepath stroke}def\n", dotw);
    for (i=0; i<pm->nT; i++) {
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->tria[3*i+0]+0], 500*pm->vertex[2*pm->tria[3*i+0]+1]);
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->tria[3*i+1]+0], 500*pm->vertex[2*pm->tria[3*i+1]+1]);
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->tria[3*i+2]+0], 500*pm->vertex[2*pm->tria[3*i+2]+1]);
	fprintf(f, "t\n");
    }
    fprintf(f, "showpage\n");
    fclose(f);
    return 0;
}

static int savewrkps(mesh2d *pm) {
    FILE *f;
    int i;

    f=fopen("wrk.ps", "w");
    if (!f) return 1;
    fprintf(f, "10 10 translate 0 setlinewidth %lf %lf scale\n", zoom, zoom);
    fprintf(f, "/e{newpath moveto lineto stroke}def\n");
    fprintf(f, "/d{newpath %lf 0 360 arc closepath stroke}def\n", dotw);
    for (i=0; i<pm->nWorkArea; i++) {
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->e[pm->workarea[i]].v[0]+0], 500*pm->vertex[2*pm->e[pm->workarea[i]].v[0]+1]);
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->e[pm->workarea[i]].v[1]+0], 500*pm->vertex[2*pm->e[pm->workarea[i]].v[1]+1]);
	fprintf(f, "e\n");
    }
    for (i=0; i<pm->nWorkSpace; i++) {
	fprintf(f, "%lf %lf d\n", 500*pm->vertex[2*pm->workspace[i]+0], 500*pm->vertex[2*pm->workspace[i]+1]);
    }
    fprintf(f, "showpage\n");
    fclose(f);
    return 0;
}

static int savebndps(mesh2d *pm) {
    FILE *f;
    int i;

    f=fopen("bnd.ps", "w");
    if (!f) return 1;
    fprintf(f, "10 10 translate 0 setlinewidth %lf %lf scale\n", zoom, zoom);
    fprintf(f, "/e{newpath moveto lineto stroke}def\n");
    fprintf(f, "/d{newpath %lf 0 360 arc closepath stroke}def\n", dotw);
    for (i=0; i<pm->nE; i++) {
	fprintf(f, "%lf %lf d\n", 500*pm->vertex[2*pm->e[i].v[0]+0], 500*pm->vertex[2*pm->e[i].v[0]+1]);
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->e[i].v[0]+0], 500*pm->vertex[2*pm->e[i].v[0]+1]);
	fprintf(f, "%lf %lf ", 500*pm->vertex[2*pm->e[i].v[1]+0], 500*pm->vertex[2*pm->e[i].v[1]+1]);
	fprintf(f, "e\n");
    }
    fprintf(f, "showpage\n");
    fclose(f);
    return 0;
}
#endif


static int savefrtraw(mesh2d *pm) {
    FILE *f;
    int i;
    static int fn=0;
    char buf[1024];
    
    sprintf(buf, "frt%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "%d %d %d\n", pm->nV, pm->nT, pm->nE);
    for (i=0; i<pm->nV; i++) {
	fprintf(f, "%20.16le %20.16le %20.16le\n", pm->vertex[2*i+0], pm->vertex[2*i+1], 0.0);
    }
    for (i=0; i<pm->nT; i++) {
	fprintf(f, "%d %d %d\n", pm->tria[3*i+0]+1, pm->tria[3*i+1]+1, pm->tria[3*i+2]+1);
    }
    for (i=0; i<pm->nE; i++) {
	fprintf(f, "%d %d\n", pm->e[i].v[0]+1, pm->e[i].v[1]+1);
    }
    fclose(f);
    return fn++;
}
