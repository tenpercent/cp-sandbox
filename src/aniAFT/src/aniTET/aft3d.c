#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "aft.h"
#include "aft3d.h"
#include "helper.h"
#include "det.h"
#include "refine.h"
#include "qual.h"
#include "check.h"

// SHOWPROGRESS -- show progress of actions
// (included from aft.h)

// PROF -- enable timing
#define PROF

// SMOOTH_CF -- (experimental) smooth tetra size when auto sizing is used (provides good results)
#define SMOOTH_CF

// TRYUGLY -- automatically start second method
#define TRYUGLY


/* for debug */
// SHOWQUAL -- (debug) show quality stats during aft front()
//#define SHOWQUAL

// STRUCTCHECK -- (debug) check data structure consistency during aft
//#define STRUCTCHECK

// INTFRONT -- (debug) check front self-intersections during aft
//#define INTFRONT

// DUMPFRONT -- (debug) periodically dump current front
//#define DUMPFRONT

// TRACE -- (debug) dump full trace of actions to text file
//#define TRACE

// EXTRADEBUG -- (debug) dump extra debug info
//#define EXTRADEBUG

/* obsolete heuristics features  NOT SUPPORTED */ 
// CHEATPOINT -- (heuristics) insert internal point in tight spaces (not robust) NOT SUPPORTED
//#define CHEATPOINT

// MEGACHEAT -- (heuristics) use advanced version of cheatpoint (still not robust) NOT SUPPORTED
//#define MEGACHEAT

// REMOVEPOINTS -- (heuristics) delete some internal points and improve front (not robust) NOT SUPPORTED
//#define REMOVEPOINTS

// STATS -- (debug) dump heuristics stats info
//#define STATS

// INITIAL_CHECK_FATAL: 1 -- initial intersections are fatal, 0 -- try with them
#define INITIAL_CHECK_FATAL 1

/****************************************************************************/;

#ifdef TRACE
static FILE *ftrace;
static void trace(char *str, ...) {
    va_list ArgList;
    va_start(ArgList, str);
    vfprintf(ftrace, str, ArgList);
    fflush(ftrace);
}
static void traceopen(void) {
    ftrace = fopen("_trace_aft3d.txt", "w");
}
static void traceclose(void) {
    fclose(ftrace);
}
#else
static void trace(char *str, ...) {(void) str;}
static void traceopen(void) {}
static void traceclose(void) {}
#endif

#ifdef TRYUGLY
/* mesh 3D ugly */
int mesh_3d_ugly (int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT);
int mesh_3d_ugly_(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int *pnnV, int *pnnF, int *pnnT);
#endif


double RBST_MIN = 256.0;
double RBST_MAX = 1024.0;

typedef struct _tnode{ /* node */
    int ne, nn;
    int *e;
    struct _tnode *child[8];
} node;

typedef struct { /* face3d */ 
    int v[3];
    int c;
    REAL x[3];
    REAL r;
    REAL s;
    int h;
    node *node;
} face3d;

typedef struct { /* mesh3d */
    int  nV, nnV, nVF;
    REAL *vertex;
    int  nF, nnF;
    face3d *f;
    int  *heap;
    int  nT, nnT;
    int  *tetra, *tetracolor;
    REAL volume, meshedvolume;
    REAL stepsize, maxarea;
    REAL minrho, alpha, beta;
    int  curfriend[3], curfriendlink[3];
    REAL curfriendness[3];
    REAL cfnx[3], cfny[3], cfnz[3];
    int  npack,  *pack;
    int  ntpack, *tpack;
    int  ncolor, *color;
    int  nWorkArea, nWorkSpace, nBadGuys, nGirls;
    int  *workarea, *workspace, *badguys, *girls;
    int  intface;
    int  allowlocal, local, cheat, applied;
    int  statsaft, statschp8, statschpnt;
    double (*fsize)(double, double, double);
    REAL bb[6];
    node *root;
    REAL robustness;
    REAL qmin, qsum;
    int  qcnt;
} mesh3d;

int detstatsprint(char *buf);

/* forward debug functions declaration*/
#ifdef STRUCTCHECK
static int saveintraw(REAL *vertex, int a, int b, int c, int d, int u, int v, int w);
#endif
static int savefrtraw(mesh3d *pm);
static int write_front_int_gmv(char *fn, int nV, double *vertex, int nF, face3d *f);



/* Metric functions *********************************************************/;

static REAL dist3(REAL x1, REAL y1, REAL z1,  REAL x2, REAL y2, REAL z2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}
static REAL dist3i(REAL *vertex, int v1, int v2) {
    return dist3(vertex[3*v1+0], vertex[3*v1+1], vertex[3*v1+2], vertex[3*v2+0], vertex[3*v2+1], vertex[3*v2+2]);
}

/* Vector functions *********************************************************/;
/* wrapper for normvec3i() from det.c */
static void normvec(REAL *vertex, int v1, int v2, int v3, REAL *px, REAL *py, REAL *pz) {
    REAL v[3];
    normvec3i(vertex, v1, v2, v3, v);
    *px = v[0];
    *py = v[1];
    *pz = v[2];
}

/* Color functions **********************************************************/;
static int pcolor(mesh3d *pm, int c) {
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
static void colordump(mesh3d *pm) {
    int i;
    for (i=0; i<pm->ncolor; i++) {
	printf("%d: %d\n", pm->color[2*i+0], pm->color[2*i+1]);
    }
}
#endif
static int cjoin(mesh3d *pm, int c1, int c2) {
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

/* Octree functions *********************************************************/;
#define TREEDEPTH 24
static void ftreemetric(mesh3d *pm, int w1, int w2, int w3, REAL *px, REAL *pr) {
    int v1, v2, v3;
    REAL r, s, x, y, z;
    
    if      ((w1 < w2) && (w1 < w3))  v1 = w1,  v2 = w2,  v3 = w3;
    else if ((w2 < w1) && (w2 < w3))  v1 = w2,  v2 = w3,  v3 = w1;
    else if ((w3 < w1) && (w3 < w2))  v1 = w3,  v2 = w1,  v3 = w2;
    else {
	v1 = w1,  v2 = w2,  v3 = w3;
	libaft_3d_warn("aft3d.c: ftreemetric(): bad input: %d, %d, %d. Internal error\n", w1, w2, w3);
    }
    x = px[0] = (pm->vertex[3*v1+0] + pm->vertex[3*v2+0] + pm->vertex[3*v3+0])/3.0;
    y = px[1] = (pm->vertex[3*v1+1] + pm->vertex[3*v2+1] + pm->vertex[3*v3+1])/3.0;
    z = px[2] = (pm->vertex[3*v1+2] + pm->vertex[3*v2+2] + pm->vertex[3*v3+2])/3.0;
    r = 0.0;
    s = dist3(pm->vertex[3*v1+0], pm->vertex[3*v1+1], pm->vertex[3*v1+2], x, y, z);
    if (s > r) r = s;
    s = dist3(pm->vertex[3*v2+0], pm->vertex[3*v2+1], pm->vertex[3*v2+2], x, y, z);
    if (s > r) r = s;
    s = dist3(pm->vertex[3*v3+0], pm->vertex[3*v3+1], pm->vertex[3*v3+2], x, y, z);
    if (s > r) r = s;
    *pr = r;
}
static void fmetric(mesh3d *pm, int i) {
    int v1, v2, v3;
    REAL x, y, z, a, b, c;
    v1 = pm->f[i].v[0];
    v2 = pm->f[i].v[1];
    v3 = pm->f[i].v[2];
    ftreemetric(pm, v1, v2, v3, pm->f[i].x, &pm->f[i].r);
    normvec(pm->vertex, v1, v2, v3, &x, &y, &z);
    a = dist3i(pm->vertex, v2, v3);
    b = dist3i(pm->vertex, v3, v1);
    c = dist3i(pm->vertex, v1, v2);
    pm->f[i].s = a*b*c/sqrt(x*x + y*y + z*z);
/*    pm->f[i].s = 2*sqrt(x*x + y*y + z*z)/(a+b+c); */
/*    pm->f[i].s = sqrt(x*x + y*y + z*z); */
/*    pm->f[i].s = 1.0; */
}
static void findbb(mesh3d *pm) {
    int i;
    REAL d;
    if (pm->nV == 0) {
	pm->bb[0] = 0.0,  pm->bb[1] = 0.0,  pm->bb[2] = 0.0;
	pm->bb[3] = 0.0,  pm->bb[4] = 0.0,  pm->bb[5] = 0.0;
	return;
    }
    pm->bb[0] = pm->vertex[0],  pm->bb[1] = pm->vertex[1],  pm->bb[2] = pm->vertex[2];
    pm->bb[3] = pm->vertex[0],  pm->bb[4] = pm->vertex[1],  pm->bb[5] = pm->vertex[2];
    for (i=1; i<pm->nV; i++) {
	if (pm->bb[0] > pm->vertex[3*i+0])  pm->bb[0] = pm->vertex[3*i+0];
	if (pm->bb[1] > pm->vertex[3*i+1])  pm->bb[1] = pm->vertex[3*i+1];
	if (pm->bb[2] > pm->vertex[3*i+2])  pm->bb[2] = pm->vertex[3*i+2];
	if (pm->bb[3] < pm->vertex[3*i+0])  pm->bb[3] = pm->vertex[3*i+0];
	if (pm->bb[4] < pm->vertex[3*i+1])  pm->bb[4] = pm->vertex[3*i+1];
	if (pm->bb[5] < pm->vertex[3*i+2])  pm->bb[5] = pm->vertex[3*i+2];
    }
    d = 0.0;
    for (i=0; i<3; i++) {
	if (d < pm->bb[3+i]-pm->bb[i])  d = pm->bb[3+i]-pm->bb[i];
    }
    for (i=0; i<3; i++) {
	pm->bb[i] = (pm->bb[i] + pm->bb[3+i] - d)/2.0;
	pm->bb[3+i] = pm->bb[i] + d;
    }
}
static int treedir(REAL *bb, REAL *x) {
    int d = 0;
    if (x[0] > (bb[0]+bb[3])/2.0) d+=1;
    if (x[1] > (bb[1]+bb[4])/2.0) d+=2;
    if (x[2] > (bb[2]+bb[5])/2.0) d+=4;
    return d;
}
static void treemove(REAL *bb, REAL s, int d) {
    if (d&1)  bb[0] += s;  else  bb[3] -= s;
    if (d&2)  bb[1] += s;  else  bb[4] -= s;
    if (d&4)  bb[2] += s;  else  bb[5] -= s;
}
static void treeback(REAL *bb, REAL s, int d) {
    if (d&1)  bb[0] -= s;  else  bb[3] += s;
    if (d&2)  bb[1] -= s;  else  bb[4] += s;
    if (d&4)  bb[2] -= s;  else  bb[5] += s;
}
static int listadd(mesh3d *pm, int *pn, int *pnn, int **pe, int n) {
#ifdef STRUCTCHECK
    int i;
    int v1, v2, v3, w1, w2, w3;
    v1 = pm->f[n].v[0];
    v2 = pm->f[n].v[1];
    v3 = pm->f[n].v[2];
    for (i=0; i<*pn; i++) {
	w1 = pm->f[(*pe)[i]].v[0];
	w2 = pm->f[(*pe)[i]].v[1];
	w3 = pm->f[(*pe)[i]].v[2];
	if ((v1!=w1)&&(v1!=w2)&&(v1!=w3)) continue;
	if ((v2!=w1)&&(v2!=w2)&&(v2!=w3)) continue;
	if ((v3!=w1)&&(v3!=w2)&&(v3!=w3)) continue;
	libaft_3d_warn("aft3d.c: listadd(): duplicate faces detected (%d %d %d)", v1, v2, v3);
	trace("duplicate faces: new [%d, %d, %d], old [%d, %d, %d]\n", v1, v2, v3, w1, w2, w3);
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
    libaft_3d_warn("aft3d.c: listdel(): element not found (internal error)");
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
    libaft_3d_warn("aft3d.c: listrep(): element not found (internal error)");
    return -1;
}
static node *nodenew(void) {
    node *p;
    int i;
    p = libaft_malloc(sizeof(node));
    p->ne = 0;
    p->nn = 0;
    p->e = 0;
    for (i=0; i<8; i++)  p->child[i] = 0;
    return p;
}
static int nodevoid(node *n) {
    int i;
    if (n->nn)  return 0;
    for (i=0; i<8; i++)  if (n->child[i])  return 0;
    if (n->e)  libaft_3d_warn("aft3d.c: nodevoid(): node void, but list is not NULL (sure memory leak)");
    return 1;
}
static node *treeadd(mesh3d *pm, int n) {
    node *c, *p;
    REAL bb[6], x[3], s, r;
    int i, d, level;

    p = pm->root;
    if (!p) {
	p = nodenew();
	pm->root = p;
    }
    for (i=0; i<6; i++)  bb[i] = pm->bb[i];
    for (i=0; i<3; i++)  x[i] = pm->f[n].x[i];
    r = pm->f[n].r;
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
	    p = nodenew();
	    c->child[d] = p;
	}
	c = p;
    }
    listadd(pm, &c->ne, &c->nn, &c->e, n);
    return c;
}
static int treedel(mesh3d *pm, int n) {
    node *c, *p;
    node *stack[TREEDEPTH];
    int dirs[TREEDEPTH];
    REAL bb[6], x[3], s, r;
    int i, d, level;

    p = pm->root;
    if (!p)  {
	libaft_3d_warn("aft3d.c: treedel(): null root in octree (internal error)");
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm->bb[i];
    for (i=0; i<3; i++)  x[i] = pm->f[n].x[i];
    r = pm->f[n].r;
    s = (pm->bb[3]-pm->bb[0])/2.0;
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
	    libaft_3d_warn("aft3d.c: treedel(): null child in octree (internal error)");
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
static int wscheck(mesh3d *pm, int i, REAL *x, REAL r) {
    return (dist3(pm->f[i].x[0], pm->f[i].x[1], pm->f[i].x[2], x[0], x[1], x[2]) <= pm->f[i].r + r) ? 1 : 0;
}
static int treewsrec(mesh3d *pm, node *c, REAL *bb, REAL s, REAL *x, REAL r) {
    int i, k = 0;
    if (!c)  return 0;
    for (i=0; i<3; i++) {
	if (bb[  i] - 2.0*s > x[i] + r)  return 0;
	if (bb[3+i] + 2.0*s < x[i] - r)  return 0;
    }
    for (i=0; i<c->ne; i++) {
	if (wscheck(pm, c->e[i], x, r))  pm->workarea[pm->nWorkArea++] = c->e[i];
    }
    k += c->ne;
    for (i=0; i<8; i++) {
	treemove(bb, s, i);
	k += treewsrec(pm, c->child[i], bb, s/2.0, x, r);
	treeback(bb, s, i);
    }
    return k;
}
static int treews(mesh3d *pm, REAL *x, REAL r) {
    node *p;
    REAL bb[6], s;
    int i;

    p = pm->root;
    if (!p)  {
	libaft_3d_warn("aft3d.c: treews(): null root in octree (internal error)");
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm->bb[i];
    s = (pm->bb[3]-pm->bb[0])/2.0;
    return treewsrec(pm, p, bb, s, x, r);
}
static int listfindface(mesh3d *pm, int ne, int *e, int v1, int v2, int v3) {
    int i, w1, w2, w3;
    for (i=0; i<ne; i++) {
	w1 = pm->f[e[i]].v[0];
	w2 = pm->f[e[i]].v[1];
	w3 = pm->f[e[i]].v[2];
	if ((v1==w1)&&(v2==w2)&&(v3==w3)) return e[i];
	if ((v2==w1)&&(v3==w2)&&(v1==w3)) return e[i];
	if ((v3==w1)&&(v1==w2)&&(v2==w3)) return e[i];
    }
    return -1;
}
static int findface(mesh3d *pm, int v1, int v2, int v3) {
    node *c, *p;
    REAL bb[6], x[3], s, r;
    int i, d, level;

    p = pm->root;
    if (!p)  {
	return -1;
    }
    for (i=0; i<6; i++)  bb[i] = pm->bb[i];
    ftreemetric(pm, v1, v2, v3, x, &r);
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
    return listfindface(pm, c->ne, c->e, v1, v2, v3);
}

/* Structure functions ******************************************************/;
/* Return -1 if maxTetra exceeded */
static int addtetra(mesh3d *pm, int v1, int v2, int v3, int v4, int color) {
    trace("addtetra(%d, %d, %d, %d,  %d)\n", v1, v2, v3, v4, color);
    color = pcolor(pm, color);
    if (pm->nT >= pm->nnT) {
	libaft_3d_stop("AFT Front: maxTetra exceeded");
	return -1;
    } else {
	pm->tetra[4*(pm->nT)+0] = v1;
	pm->tetra[4*(pm->nT)+1] = v2;
	pm->tetra[4*(pm->nT)+2] = v3;
	pm->tetra[4*(pm->nT)+3] = v4;
	pm->tetracolor[pm->nT] = color;
	return (pm->nT)++;
    }
}
/* Return -1 if maxVertex exceeded */
static int addvertex(mesh3d *pm, REAL x, REAL y, REAL z, int crit) {
#ifdef STRUCTCHECK
    int i;
    if (crit) {
	for (i=0; i<pm->nV; i++) {
	    if (fabs(pm->vertex[3*i+0]-x) > 0.0) continue;
	    if (fabs(pm->vertex[3*i+1]-y) > 0.0) continue;
	    if (fabs(pm->vertex[3*i+2]-z) > 0.0) continue;
	    libaft_3d_warn("aft3d.c: addvertex(): duplicate points detected");
	}
    }
#endif
    if (crit==1) trace("addvertex(%.18le, %.18le, %.18le,  %d)\n", x, y, z, crit);
    if (pm->nV >= pm->nnV) { /* 3 vertices are reserved for temporary use */
	libaft_3d_stop("AFT Front: maxVertex exceeded");
	return -1;
    } else {
	pm->vertex[3*(pm->nV)+0] = x;
	pm->vertex[3*(pm->nV)+1] = y;
	pm->vertex[3*(pm->nV)+2] = z;
	return (pm->nV >= pm->nnV - 3) ? -2 : (pm->nV)++;
    }
}
static void heapswap(mesh3d *pm, int k, int m) {
    int t;
    t = pm->heap[k],  pm->heap[k] = pm->heap[m],  pm->heap[m] = t;
    pm->f[pm->heap[k]].h = k;
    pm->f[pm->heap[m]].h = m;
}
static int siftup(mesh3d *pm, int n) {
    int i;
    while (1) {
	i = (n-1)/2;
	if ((i<0) || (pm->f[pm->heap[i]].s <= pm->f[pm->heap[n]].s)) break;
	heapswap(pm, i, n);
	n = i;
    }
    return n;
}
static int siftdown(mesh3d *pm, int n) {
    int i;
    n = siftup(pm, n);
    while (1) {
	i = 0;
	if ((2*n+1 < pm->nF) && (pm->f[pm->heap[n]].s > pm->f[pm->heap[(2*n+1)]].s)) {
	    if ((2*n+2 < pm->nF) && (pm->f[pm->heap[(2*n+1)]].s > pm->f[pm->heap[(2*n+2)]].s))  i = 2*n+2;
	    else  i = 2*n + 1;
	} else if ((2*n+2 < pm->nF) && (pm->f[pm->heap[n]].s > pm->f[pm->heap[(2*n+2)]].s))  i = 2*n+2;
	if (i) {
	    heapswap(pm, n, i);
	    n = i;
	} else break;
    }
    return n;
}
/* Return -1 if maxFace exceeded */
static int addface(mesh3d *pm, int v1, int v2, int v3, int color) {
    trace("addface(%d, %d, %d,  %d)\n", v1, v2, v3, color);
#ifdef STRUCTCHECK
    int i;
    for (i=0; i<pm->nF; i++) {
	if ((v1!=pm->f[i].v[0])&&(v1!=pm->f[i].v[1])&&(v1!=pm->f[i].v[2])) continue;
	if ((v2!=pm->f[i].v[0])&&(v2!=pm->f[i].v[1])&&(v2!=pm->f[i].v[2])) continue;
	if ((v3!=pm->f[i].v[0])&&(v3!=pm->f[i].v[1])&&(v3!=pm->f[i].v[2])) continue;
	libaft_3d_warn("aft3d.c: addface(): duplicate faces detected (bad front?)");
/*	savefrtraw(pm);*/
    }
#endif
    color = pcolor(pm, color);
    if (pm->nF >= pm->nnF) {
	libaft_3d_stop("AFT Front: maxFace exceeded");
	return -1;
    } else {
	pm->f[(pm->nF)].v[0] = v1;
	pm->f[(pm->nF)].v[1] = v2;
	pm->f[(pm->nF)].v[2] = v3;
	pm->f[pm->nF].c = color;
	fmetric(pm, pm->nF);
	pm->heap[pm->nF] = pm->nF;
	pm->f[pm->nF].h = pm->nF;
	pm->f[pm->nF].node = treeadd(pm, pm->nF);
	return siftup(pm, (pm->nF)++);
    }
}
static void remface(mesh3d *pm, int i, int color) {
    cjoin(pm, pm->f[i].c, color);
    if (pm->npack > pm->nnF)  {
	libaft_3d_stop("aft3d.c: Internal error in remface & packface");
	return;
    }
    pm->pack[pm->npack] = i;
    pm->npack++;
}
static int packface(mesh3d *pm) {
    int i, j;
    int err = 0;
    node *c;
    while (pm->npack > 0) {
	pm->npack--;
	j = pm->pack[pm->npack];
	pm->nF--;
	for (i=0; i<pm->npack; i++) if (pm->pack[i] == pm->nF) pm->pack[i] = j;
	c = pm->f[(pm->nF)].node;
	if (treedel(pm, j))  err++;
	if (listrep(&c->ne, &c->e, pm->nF, j))  err++;
	pm->f[j].v[0] = pm->f[(pm->nF)].v[0];
	pm->f[j].v[1] = pm->f[(pm->nF)].v[1];
	pm->f[j].v[2] = pm->f[(pm->nF)].v[2];
	pm->f[j].c = pm->f[pm->nF].c;
	pm->f[j].x[0] = pm->f[(pm->nF)].x[0];
	pm->f[j].x[1] = pm->f[(pm->nF)].x[1];
	pm->f[j].x[2] = pm->f[(pm->nF)].x[2];
	pm->f[j].r = pm->f[(pm->nF)].r;
	pm->f[j].s = pm->f[(pm->nF)].s;
	pm->f[j].node = c;
	i = pm->f[j].h;
	pm->nF++;
	siftdown(pm, i);
	pm->nF--;
	i = pm->f[pm->nF].h;
	pm->heap[i] = pm->heap[pm->nF];
	pm->f[pm->heap[i]].h = i;
	siftdown(pm, i);
    }
    return err;
}
#ifdef REMOVEPOINTS
static void remtetra(mesh3d *pm, int i, int color) {
    cjoin(pm, pm->tetracolor[i], color);
    if (pm->ntpack > pm->nnT) libaft_3d_stop("aft3d.c: Internal error in remtetra & packtetra");
    pm->tetra[4*i+0] = -1;
    pm->tetra[4*i+1] = -1;
    pm->tetra[4*i+2] = -1;
    pm->tetra[4*i+3] = -1;
    pm->tpack[pm->ntpack] = i;
    pm->ntpack++;
}
static void packtetra(mesh3d *pm) {
    int i, j;
    while (pm->ntpack > 0) {
	pm->ntpack--;
	j = pm->tpack[pm->ntpack];
	pm->nT--;
	pm->tetra[4*j+0] = pm->tetra[4*(pm->nT)+0];
	pm->tetra[4*j+1] = pm->tetra[4*(pm->nT)+1];
	pm->tetra[4*j+2] = pm->tetra[4*(pm->nT)+2];
	pm->tetra[4*j+3] = pm->tetra[4*(pm->nT)+3];
	pm->tetracolor[j] = pm->tetracolor[pm->nT];
	for (i=0; i<pm->ntpack; i++) if (pm->tpack[i] == pm->nT) pm->tpack[i] = j;
    }
}
#endif
/* Vertex sets functions */;
static int VinWS(mesh3d *pm, int v) {
    int i;
    for (i=0; i<pm->nWorkSpace; i++)  if (v==pm->workspace[i])  return 1;
    return 0;
}
static int VinBG(mesh3d *pm, int v) {
    int i;
    for (i=0; i<pm->nBadGuys; i++)  if (v==pm->badguys[i])  return 1;
    return 0;
}
static int VinG(mesh3d *pm, int v) {
    int i;
    for (i=0; i<pm->nGirls; i++)  if (v==pm->girls[i])  return 1;
    return 0;
}

/* Intersection functions ***************************************************/;
static int intsecthelper(mesh3d *pm, int *f, int *e) {
    int i1, i2, k1, k2, k3;
    int uv, dup=0, ind=0;

    /* volume sign check */
    /*if (idet3i4(pm->vertex, f[0], f[1], f[2], f[3]) != 1)  libaft_3d_warn("intsecthelper: volume sign: idet=%d\n", idet3i4(pm->vertex, f[0], f[1], f[2], f[3]));*/

    if ((e[0]==f[0])||(e[0]==f[1])||(e[0]==f[2])) {dup++;ind=1;}
    if ((e[1]==f[0])||(e[1]==f[1])||(e[1]==f[2])) {dup++;ind=0;}
    if (dup==0) {
	i1=idet3i4(pm->vertex, f[0], f[1], f[2], e[0]);
	i2=idet3i4(pm->vertex, f[0], f[1], f[2], e[1]);
	if (i1*i2 == 1) return 0;
	else if (i1*i2 == -1) {
	    k1=idet3i4(pm->vertex, e[0], e[1], f[0], f[1]);
	    k2=idet3i4(pm->vertex, e[0], e[1], f[1], f[2]);
	    k3=idet3i4(pm->vertex, e[0], e[1], f[2], f[0]);
	    if ((k1*k2==-1) || (k2*k3==-1) || (k3*k1==-1)) return 0;
	} else if (i1 != i2) {
	    if (i1 == 0)  ind = 0;
	    if (i2 == 0)  ind = 1;
	    k1=idet3i4(pm->vertex, e[ind], f[3], f[0], f[1]);
	    k2=idet3i4(pm->vertex, e[ind], f[3], f[1], f[2]);
	    k3=idet3i4(pm->vertex, e[ind], f[3], f[2], f[0]);
	    if ((k1*k2==-1) || (k2*k3==-1) || (k3*k1==-1)) return 0;
	} else { /* i1 = i2 = 0 */
	    uv = idet3i4(pm->vertex, e[0], e[1], f[0], f[3]) 
		+idet3i4(pm->vertex, e[0], e[1], f[1], f[3])
		+idet3i4(pm->vertex, e[0], e[1], f[2], f[3]);
	    if ((uv==3) || (uv==-3)) return 0;
	    if (idet3i4(pm->vertex, f[0], f[1], e[0], f[3]) + idet3i4(pm->vertex, f[0], f[1], e[1], f[3]) == -2) return 0;
	    if (idet3i4(pm->vertex, f[1], f[2], e[0], f[3]) + idet3i4(pm->vertex, f[1], f[2], e[1], f[3]) == -2) return 0;
	    if (idet3i4(pm->vertex, f[2], f[0], e[0], f[3]) + idet3i4(pm->vertex, f[2], f[0], e[1], f[3]) == -2) return 0;
	}
    } else if (dup==1) {
	if (idet3i4(pm->vertex, f[0], f[1], f[2], e[ind])) return 0;
	else {
	    if (idet3i4(pm->vertex, f[0], f[1], e[ind], f[3]) == -1) return 0;
	    if (idet3i4(pm->vertex, f[2], f[0], e[ind], f[3]) == -1) return 0;
	    if (idet3i4(pm->vertex, f[1], f[2], e[ind], f[3]) == -1) return 0;
	}
    } else return 0;
    return 1;
}
static int intsectjerry(mesh3d *pm, int *t, int *f) {
    REAL x, y, z, xc, yc, zc;
    int d, i, j;
    int mf[4], me[4];
    int dup, ind;

    /* volume sign check */
    /*if (idet3i4(pm->vertex, t[0], t[1], t[2], t[3]) != 1)  libaft_3d_warn("intsectjerry: volume sign: idet=%d\n", idet3i4(pm->vertex, t[0], t[1], t[2], t[3]));*/

    dup = 0; ind = 0;
    if ((f[0]==t[0])||(f[0]==t[1])||(f[0]==t[2])||(f[0]==t[3])) {dup++;ind+=0;}
    if ((f[1]==t[0])||(f[1]==t[1])||(f[1]==t[2])||(f[1]==t[3])) {dup++;ind+=1;}
    if ((f[2]==t[0])||(f[2]==t[1])||(f[2]==t[2])||(f[2]==t[3])) {dup++;ind+=2;}
    if (dup==3) return 0;

    me[0] = f[0]; me[1] = f[1]; me[2] = f[2]; me[3] = f[0];
    mf[0] = t[3]; mf[1] = t[2]; mf[2] = t[1]; mf[3] = t[0];
    for (i=0; i<3; i++) if (intsecthelper(pm, mf, me+i)) return 1;
    mf[0] = t[2]; mf[1] = t[3]; mf[2] = t[0]; mf[3] = t[1];
    for (i=0; i<3; i++) if (intsecthelper(pm, mf, me+i)) return 1;
    mf[0] = t[1]; mf[1] = t[0]; mf[2] = t[3]; mf[3] = t[2];
    for (i=0; i<3; i++) if (intsecthelper(pm, mf, me+i)) return 1;
    mf[0] = t[0]; mf[1] = t[1]; mf[2] = t[2]; mf[3] = t[3];
    for (i=0; i<3; i++) if (intsecthelper(pm, mf, me+i)) return 1;

    normvec(pm->vertex, f[0], f[1], f[2], &x, &y, &z);
    xc = (pm->vertex[3*f[0]+0] + pm->vertex[3*f[1]+0] + 2*pm->vertex[3*f[2]+0])/4.0;
    yc = (pm->vertex[3*f[0]+1] + pm->vertex[3*f[1]+1] + 2*pm->vertex[3*f[2]+1])/4.0;
    zc = (pm->vertex[3*f[0]+2] + pm->vertex[3*f[1]+2] + 2*pm->vertex[3*f[2]+2])/4.0;
    d = addvertex(pm, xc+1.0*x, yc+1.0*y, zc+1.0*z, 0);
    if (d < 0)  d = pm->nV;
    else  pm->nV--;

    mf[0] = f[0]; mf[1] = f[1]; mf[2] = f[2]; mf[3] = d;
    for (i=0; i<3; i++) {
	me[0]=t[i];
	for (j=i+1; j<4; j++) {
	    me[1]=t[j];
	    if (intsecthelper(pm, mf, me)) return 1;
	}
    }

    if (dup==2) {
	i = 3-ind;
	j =  idet3i4(pm->vertex, f[i], t[1], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], f[i], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], f[i], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], t[2], f[i]);
/*	if (j==3)  libaft_3d_warn("intsectjerry: dup=2: j=3\n");
	if (j<=-3)  libaft_3d_warn("intsectjerry: dup=2: j<=-3\n");*/
	if ((j>=3) || (j<=-3)) return 1;
    } else if (dup==1) {
	if (ind==0) {i=1;j=2;}
	else if (ind==1) {i=0;j=2;}
	else {i=0;j=1;}
	i =  idet3i4(pm->vertex, f[i], t[1], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], f[i], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], f[i], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], t[2], f[i]);
	j =  idet3i4(pm->vertex, f[j], t[1], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], f[j], t[2], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], f[j], t[3])
	    +idet3i4(pm->vertex, t[0], t[1], t[2], f[j]);
/*	if (i<=-3)  libaft_3d_warn("intsectjerry: dup=1: i<=-3\n");
	if (j<=-3)  libaft_3d_warn("intsectjerry: dup=1: j<=-3\n");*/
	if (((i>=3) || (i<=-3)) && ((j>=3) || (j<=-3))) return 1;
    }

    return 0;
}
static int intsect(mesh3d *pm, int a, int b, int c, int d, int u, int v, int w) {
    int t[4]={a, b, c, d};
    int f[3]={u, v, w};
    return intsectjerry(pm, t, f);
}
static int intfacejerry(mesh3d *pm, int a, int b, int c, int u, int v, int w) {
    REAL x, y, z, xc, yc, zc, r, p;
    int d, i1, i2;
    normvec(pm->vertex, a, b, c, &x, &y, &z);
    r = sqrt(x*x + y*y + z*z);
    p = (dist3i(pm->vertex, a, b) + dist3i(pm->vertex, b, c)+ dist3i(pm->vertex, c, a))/3.0;
    if (r > 0.0)  r = p / r;
    else  libaft_3d_warn("aft3d.c: intfacejerry(): zero normal, very strange");
    x *= r,  y *= r,  z *= r;
    xc = (pm->vertex[3*a+0] + pm->vertex[3*b+0] + pm->vertex[3*c+0])/3.0;
    yc = (pm->vertex[3*a+1] + pm->vertex[3*b+1] + pm->vertex[3*c+1])/3.0;
    zc = (pm->vertex[3*a+2] + pm->vertex[3*b+2] + pm->vertex[3*c+2])/3.0;
    d = addvertex(pm, xc+x, yc+y, zc+z, 0);
    if (d < 0)  d = pm->nV++;
    i1 = intsect(pm, a, b, c, d, u, v, w);
    pm->nV--;
    if (i1==0) return 0;
    d = addvertex(pm, xc-x, yc-y, zc-z, 0);
    if (d < 0)  d = pm->nV++;
    i2 = intsect(pm, b, a, c, d, u, v, w);
    pm->nV--;
    if (i2==0) return 0;
    return 1;
}
static int intface(mesh3d *pm, int a, int b, int c, int u, int v, int w) {
    return intfacejerry(pm, a, b, c, u, v, w);
}
static int intsectfront(mesh3d *pm, int a, int b, int c) {
    int i;
    for (i=0; i<pm->nF; i++) {
	if (intface(pm, pm->f[i].v[0], pm->f[i].v[1], pm->f[i].v[2], a, b, c)) {
	    return 1;
	}
    }
    return 0;
}

static int intfront(mesh3d *pm) {
    int i, j, k;
    for (i=0; i<pm->nF; i++) {
	pm->nWorkArea = 0;
	if (treews(pm, pm->f[i].x, pm->f[i].r) < 0)  return -1;
	for (k=0; k<pm->nWorkArea; k++) {
	    j = pm->workarea[k];
	    if (i==j)  continue;
	    if (intface(pm, pm->f[i].v[0], pm->f[i].v[1], pm->f[i].v[2],
			pm->f[j].v[0], pm->f[j].v[1], pm->f[j].v[2])) {
		pm->f[i].c = -256;
		pm->f[j].c = -256;
		write_front_int_gmv("front-intersection.gmv", pm->nV, pm->vertex, pm->nF, pm->f);
		libaft_3d_warn("AFT Front: Intersection in front (dumped to front-intersection.gmv): %d",
			intface(pm, pm->f[i].v[0], pm->f[i].v[1], pm->f[i].v[2],
			    pm->f[j].v[0], pm->f[j].v[1], pm->f[j].v[2]));
		return 1;
	    }
	}
/*	printf("\t\t\t\t\t\t\t\tintfront: %4d/%4d/%6.2lf\r", i+1, pm->nF, 100.0*(i+1)/pm->nF);
	if (i % 1 == 0)  fflush(stdout);*/
    }
/*    fflush(stdout);*/
    return 0;
}


/* AFT section */;
static int checktetra(mesh3d *pm, int v1, int v2, int v3, int v4, int log) {
    static double qmin=1e10, qmax=0.0;
    double q;
    (void) qmin, (void) qmax, (void) log;

    q = aft_tetra_qual(pm->vertex, v1, v2, v3, v4);

#ifdef SHOWPROGRESS
#ifdef SHOWQUAL
    if (log) {
	if (qmin>q) qmin = q;
	if (qmax<q) qmax = q;
	printf("\r\t\t\t\t\t\t\t\tqual: %10.8lf/%10.8lf/%10.8lf\r", q, qmin, qmax),  fflush(stdout);
    }
#endif
#endif

/*    if (q<1e-8) return 1;*/
    if (q<1e-3) return 2;
    return 0;
}

static int tricheck(REAL *v, int a, int b, int c, int d, REAL robustness) {
    REAL x[3], r, alpha, beta, gamma;
    normvec(v, a, b, c, &x[0], &x[1], &x[2]);
    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    r = (dist3i(v, a, b)+dist3i(v, a, c))/2.0/r;
    x[0] *= r,  x[1] *= r,  x[2] *= r;
    x[0] += v[3*a+0],  x[1] += v[3*a+1],  x[2] += v[3*a+2];
    gamma = orient3d(v+3*b, v+3*c, v+3*d, v+3*a)/orient3d(v+3*b, v+3*c, x, v+3*a);
    if (fabs(gamma) > 1.0/robustness)  return 0;
    alpha = orient3d(v+3*d, v+3*c, x, v+3*a)/orient3d(v+3*b, v+3*c, x, v+3*a);
    if (alpha < -1.0/robustness) return 0;
    beta  = orient3d(v+3*b, v+3*d, x, v+3*a)/orient3d(v+3*b, v+3*c, x, v+3*a);
    if (beta < -1.0/robustness) return 0;
/*    y[0] = v[3*a+0] + alpha*(v[3*b+0] - v[3*a+0]) + beta*(v[3*c+0] - v[3*a+0]) + gamma*(x[0] - v[3*a+0]) - v[3*d+0];
    y[1] = v[3*a+1] + alpha*(v[3*b+1] - v[3*a+1]) + beta*(v[3*c+1] - v[3*a+1]) + gamma*(x[1] - v[3*a+1]) - v[3*d+1];
    y[2] = v[3*a+2] + alpha*(v[3*b+2] - v[3*a+2]) + beta*(v[3*c+2] - v[3*a+2]) + gamma*(x[2] - v[3*a+2]) - v[3*d+2];
    r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    if (fabs(r) > 1e-8)  libaft_3d_warn("tricheck failed: %le (%lf, %lf, %lf)\n", r, alpha, beta, gamma);
    trace("tricheck (%d, %d, %d,  %d): %lf, %lf, %lf\n", a, b, c, d, alpha, beta, gamma);
    x[0] -= v[3*a+0],  x[1] -= v[3*a+1],  x[2] -= v[3*a+2];
    trace("dists: %lf, %lf, %lf\n", dist3i(v, a, b), dist3i(v, a, c), sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));*/
    //libaft_3d_warn("tricheck: %lf, %lf, %lf\n", alpha, beta, gamma);
    return 1;
}

static int extracheck(mesh3d *pm, int fn, int pn, int w) {
    int i, j, d = 0, m, k;
    int *vf, *vw;
    REAL ex, ey, ez, gx, gy, gz, r, calpha;

    vf = pm->f[fn].v;
    vw = pm->f[w].v;
    trace("extracheck: [%d, %d, %d], %d,  [%d, %d, %d]\n", vf[0], vf[1], vf[2], pn, vw[0], vw[1], vw[2]);
    for (i=0; i<3; i++)  for (j=0; j<3; j++)  if (vw[i] == vf [j])  d++, m = i, k = j;
    if (d > 1)  return 0;
    if (d == 0) {
	for (i=0; i<3; i++) {
	    if (vw[i] != pn)  continue;
	    for (j=0; j<3; j++)  {
		if (tricheck(pm->vertex, vw[(i+0)%3], vw[(i+1)%3], vw[(i+2)%3], vf[j], pm->robustness))  return 1;
	    }
	}
	return 0;
    }
    if (vw[(m+1)%3] == pn)  normvec(pm->vertex, vw[0], vw[1], vw[2], &ex, &ey, &ez);
    else  normvec(pm->vertex, vw[2], vw[1], vw[0], &ex, &ey, &ez);
    r = sqrt(ex*ex + ey*ey + ez*ez);
    ex /= r, ey /= r, ez /= r;
    normvec(pm->vertex, vf[(k+1)%3], vf[k], pn, &gx, &gy, &gz);
    r = sqrt(gx*gx + gy*gy + gz*gz);
    gx /= r, gy /= r, gz /= r;
    calpha = ex*gx+ ey*gy + ez*gz;
    trace("checking: [%d, %d, %d] and [%d, %d, %d]: ", vw[0], vw[1], vw[2], vf[(k+1)%3], vf[k], pn);
    trace("calpha = %.12lf  \n", calpha);
    if (calpha > cos(4.0*M_PI/pm->robustness))  {
	return 1;
    }
    normvec(pm->vertex, vf[(k+2)%3], vf[k], pn, &gx, &gy, &gz);
    r = sqrt(gx*gx + gy*gy + gz*gz);
    gx /= r, gy /= r, gz /= r;
    calpha = ex*gx+ ey*gy + ez*gz;
    trace("checking: [%d, %d, %d] and [%d, %d, %d]: ", vw[0], vw[1], vw[2], vf[k], vf[(k+2)%3], pn);
    trace("calpha = %.12lf  \n", calpha);
    if (calpha > cos(4.0*M_PI/pm->robustness))  {
	return 1;
    }
    return 0;
}

static int extracheck2(mesh3d *pm, int fn, int pn, int w) {
    int i, j, d = 0, m, k;
    int *vf, *vw;

    vf = pm->f[fn].v;
    vw = pm->f[w].v;
    trace("extracheck: [%d, %d, %d], %d,  [%d, %d, %d]\n", vf[0], vf[1], vf[2], pn, vw[0], vw[1], vw[2]);
    for (i=0; i<3; i++)  for (j=0; j<3; j++)  if (vw[i] == vf [j])  d++, m = i, k = j;
    if (d != 1)  return 0;
    if (tricheck(pm->vertex, vw[(m+0)%3], vw[(m+1)%3], vw[(m+2)%3], pn, pm->robustness))  return 1;
    return 0;
}

static int check(mesh3d *pm, int fn, int pn) {
    int v1, v2, v3, i, p1, p2, p3, i_max;
    REAL ex, ey, ez, r, calpha, calpha_max;
    
    pm->intface = -1;
    v1 = pm->f[fn ].v[0];
    v2 = pm->f[fn ].v[1];
    v3 = pm->f[fn ].v[2];
    if (idet3i4(pm->vertex, v1, v2, v3, pn) != 1) {
	libaft_3d_warn("aft3d.c: check(): tetra volume is %s (%d %d %d %d)",
		       idet3i4(pm->vertex, v1, v2, v3, pn) ? "negative" : "nearly zero",
		       v1, v2, v3, pn);
	return 1;
    }
    i = checktetra(pm, v1, v2, v3, pn, 0);
    if ((i==2) && (pm->local)) return 2;
    if (i==1) return 1;
    
    trace("Check: (%d, %d, %d,  %d)\n", v1, v2, v3, pn);
    calpha_max = -1.0,  i_max = -1;
    for (i=0; i<3; i++) {
	if (pm->curfriend[i] == pn)  continue;
	normvec(pm->vertex, pm->f[fn].v[(i+2)%3], pm->f[fn].v[(i+1)%3], pn, &ex, &ey, &ez);
	r = sqrt(ex*ex + ey*ey + ez*ez);
	ex /= r, ey /= r, ez /= r;
	calpha = ex*pm->cfnx[i]+ ey*pm->cfny[i] + ez*pm->cfnz[i];
	trace("%d: (%d, %d, %d), %d: %le\n", i, pm->f[fn].v[(i+2)%3], pm->f[fn].v[(i+1)%3], pn, pm->curfriend[i], calpha);
	if (calpha > calpha_max) {
	    calpha_max = calpha;
	    i_max = i;
	}
    }
    if ((calpha_max > cos(4.0*M_PI/pm->robustness)) && (i_max >= 0)) {
	pm->intface = pm->curfriendlink[i_max];
	return 1;
    }
    for (i=0; i<pm->nWorkArea; i++) {
	p1 = pm->f[pm->workarea[i]].v[0];
	p2 = pm->f[pm->workarea[i]].v[1];
	p3 = pm->f[pm->workarea[i]].v[2];
	if ((p1==pn) || (p2==pn) || (p3==pn))  {
	    if (extracheck(pm, fn, pn, pm->workarea[i]))  return 1;
	} else {
	    if (extracheck2(pm, fn, pn, pm->workarea[i]))  return 1;
	}
	if (intsect(pm, v1, v2, v3, pn, p1, p2, p3)) {
	    pm->intface = pm->workarea[i];
	    return 1;
	}
    }
    return 0;
}

static int cheatpoint8(mesh3d *pm, int v1, int v2, int v3, int w1, int w2, int w3, int fn, int f1, int f2, int f3) {  // obsolete
    int g1, g2, g3, gn, vw, r;
    REAL p1, p2, p3, q1, q2, q3, sum, summ;
    REAL a12, a13, a21, a23, a31, a32, b1, b2, b3, c1, c2, c3, d, pp1, pp2, pp3, qq1, qq2, qq3;

    libaft_3d_warn("aft_front: cheatpoint8() is obsolete, do not use it!");

    if ((w1==w2) || (w2==w3) || (w3==w1)) return 0;

    if (intsectfront(pm, w1, w2, w3)) return 0;
    if (intsectfront(pm, v1, w2, w3)) return 0;
    if (intsectfront(pm, w1, v2, w3)) return 0;
    if (intsectfront(pm, w1, w2, v3)) return 0;
    trace("cheatpoint8(%d, %d, %d,  %d, %d, %d,  %d,  %d, %d, %d)\n", v1, v2, v3, w1, w2, w3, fn, f1, f2, f3);
    gn = findface(pm, w3, w2, w1);
    g1 = findface(pm, v1, w2, w3);
    g2 = findface(pm, v2, w3, w1);
    g3 = findface(pm, v3, w1, w2);
    if ((gn<0)||(g1<0)||(g2<0)||(g3<0)) {
	trace("linked [%d, %d %d %d]\n", gn, g1, g2, g3);
/*	libaft_3d_warn("linked [%d, %d %d %d]", gn, g1, g2, g3);
	return 0;*/
    }
    a12 = det3i4(pm->vertex, v3, v2, w1, w2);
    a13 = det3i4(pm->vertex, v3, v2, w1, w3);
    a21 = det3i4(pm->vertex, v3, w2, v1, w1);
    a23 = det3i4(pm->vertex, v3, w2, v1, w3);
    a31 = det3i4(pm->vertex, w3, v2, v1, w1);
    a32 = det3i4(pm->vertex, w3, v2, v1, w2);
    b1 = det3i4(pm->vertex, v1, v2, v3, w1);
    b2 = det3i4(pm->vertex, v1, v2, v3, w2);
    b3 = det3i4(pm->vertex, v1, v2, v3, w3);
    c1 = det3i4(pm->vertex, w3, w2, w1, v1);
    c2 = det3i4(pm->vertex, w3, w2, w1, v2);
    c3 = det3i4(pm->vertex, w3, w2, w1, v3);

    d = a12*a23*a31 + a13*a21*a32;
    p1 =   a32*a13 + a23*a12 - a13*a12;
    p2 =   a23*a31 - a23*a21 + a13*a21;
    p3 = - a32*a31 + a32*a21 + a12*a31;
    q1 =   p3*b3*a12*a23 + p2*b2*a13*a32 - p1*b1*a23*a32;
    q2 =   p3*b3*a13*a21 - p2*b2*a13*a31 + p1*b1*a23*a31;
    q3 = - p3*b3*a12*a21 + p2*b2*a12*a31 + p1*b1*a32*a21;
    if (d<0.0) {p1 *= -d, p2 *= -d, p3 *= -d;}
    qq1 =   a12*a23 + a13*a32 - a23*a32;
    qq2 =   a13*a21 - a13*a31 + a23*a31;
    qq3 = - a12*a21 + a12*a31 + a32*a21;
    pp1 =   qq3*c3*a32*a13 + qq2*c2*a23*a12 - qq1*c1*a13*a12;
    pp2 =   qq3*c3*a23*a31 - qq2*c2*a23*a21 + qq1*c1*a13*a21;
    pp3 = - qq3*c3*a32*a31 + qq2*c2*a32*a21 + qq1*c1*a12*a31;
    if (d<0.0) {qq1 *= -d, qq2 *= -d, qq3 *= -d;}

    sum = p1 + p2 + p3 + q1 + q2 + q3;
    summ = pp1 + pp2 + pp3 + qq1 + qq2 + qq3;
    p1 /= sum, p2 /= sum, p3 /= sum;
    q1 /= sum, q2 /= sum, q3 /= sum;
    pp1 /= summ, pp2 /= summ, pp3 /= summ;
    qq1 /= summ, qq2 /= summ, qq3 /= summ;
    p1 = (p1+pp1)/2.0;
    p2 = (p2+pp2)/2.0;
    p3 = (p3+pp3)/2.0;
    q1 = (q1+qq1)/2.0;
    q2 = (q2+qq2)/2.0;
    q3 = (q3+qq3)/2.0;
    if ((p1<0)||(p2<0)||(p3<0)||(q1<0)||(q2<0)||(q3<0)) {
	trace("Strange ;)\n%lf %lf %lf\n%lf %lf %lf\n", p1, p2, p3, q1, q2, q3);
/*	libaft_3d_warn("Strange ;)\n%lf %lf %lf\n%lf %lf %lf\n", p1, p2, p3, q1, q2, q3);*/
    }
    vw = addvertex(pm, 
	    (p1*pm->vertex[3*v1+0]+p2*pm->vertex[3*v2+0]+p3*pm->vertex[3*v3+0]+q1*pm->vertex[3*w1+0]+q2*pm->vertex[3*w2+0]+q3*pm->vertex[3*w3+0]),
	    (p1*pm->vertex[3*v1+1]+p2*pm->vertex[3*v2+1]+p3*pm->vertex[3*v3+1]+q1*pm->vertex[3*w1+1]+q2*pm->vertex[3*w2+1]+q3*pm->vertex[3*w3+1]),
	    (p1*pm->vertex[3*v1+2]+p2*pm->vertex[3*v2+2]+p3*pm->vertex[3*v3+2]+q1*pm->vertex[3*w1+2]+q2*pm->vertex[3*w2+2]+q3*pm->vertex[3*w3+2]), 1);
    if (vw < 0)  return -1;
    r = idet3i4(pm->vertex, v1, v2, v3, vw)
	+idet3i4(pm->vertex, v3, v2, w1, vw)
	+idet3i4(pm->vertex, v3, w2, v1, vw)
	+idet3i4(pm->vertex, w3, v2, v1, vw)
	+idet3i4(pm->vertex, v1, w2, w3, vw)
	+idet3i4(pm->vertex, w1, v2, w3, vw)
	+idet3i4(pm->vertex, w1, w2, v3, vw)
	+idet3i4(pm->vertex, w3, w2, w1, vw);
    if (r<8) {
	/*		libaft_3d_warn("cheat8: r=%d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n%20.15lf %d\n\n%lf %lf %lf\n%lf %lf %lf\n", r,
			det3i4(pm->vertex, v1, v2, v3, vw), idet3i4(pm->vertex, v1, v2, v3, vw),
			det3i4(pm->vertex, v3, v2, w1, vw), idet3i4(pm->vertex, v3, v2, w1, vw),
			det3i4(pm->vertex, v3, w2, v1, vw), idet3i4(pm->vertex, v3, w2, v1, vw),
			det3i4(pm->vertex, w3, v2, v1, vw), idet3i4(pm->vertex, w3, v2, v1, vw),
			det3i4(pm->vertex, v1, w2, w3, vw), idet3i4(pm->vertex, v1, w2, w3, vw),
			det3i4(pm->vertex, w1, v2, w3, vw), idet3i4(pm->vertex, w1, v2, w3, vw),
			det3i4(pm->vertex, w1, w2, v3, vw), idet3i4(pm->vertex, w1, w2, v3, vw),
			det3i4(pm->vertex, w3, w2, w1, vw), idet3i4(pm->vertex, w3, w2, w1, vw),
			p1, p2, p3, q1, q2, q3);*/
	pm->nV--;
	return 0;
    } else {
	pm->statschp8++;
	trace("cheat8: success\n%lf %lf %lf\n%lf %lf %lf\n", p1, p2, p3, q1, q2, q3);
    }
    if (addtetra(pm, v1, v2, v3, vw, pm->f[fn].c) < 0)  return -1;
    if (addtetra(pm, v3, v2, w1, vw, pm->f[f1].c) < 0)  return -1;
    if (addtetra(pm, v3, w2, v1, vw, pm->f[f2].c) < 0)  return -1;
    if (addtetra(pm, w3, v2, v1, vw, pm->f[f3].c) < 0)  return -1;
    if (addtetra(pm, v1, w2, w3, vw, pm->f[f1].c) < 0)  return -1;
    if (addtetra(pm, w1, v2, w3, vw, pm->f[f2].c) < 0)  return -1;
    if (addtetra(pm, w1, w2, v3, vw, pm->f[f3].c) < 0)  return -1;
    if (addtetra(pm, w3, w2, w1, vw, pm->f[fn].c) < 0)  return -1;
    remface(pm, fn, pm->f[fn].c);
    remface(pm, f1, pm->f[fn].c);
    remface(pm, f2, pm->f[fn].c);
    remface(pm, f3, pm->f[fn].c);
    if (gn>=0)  remface(pm, gn, pm->f[fn].c);
    if (g1>=0)  remface(pm, g1, pm->f[fn].c);
    if (g2>=0)  remface(pm, g2, pm->f[fn].c);
    if (g3>=0)  remface(pm, g3, pm->f[fn].c);
    packface(pm);
    if (gn<0)  addface(pm, w1, w2, w3, pm->f[fn].c);
    if (g1<0)  addface(pm, w3, w2, v1, pm->f[fn].c);
    if (g2<0)  addface(pm, w3, v2, w1, pm->f[fn].c);
    if (g3<0)  addface(pm, v3, w2, w1, pm->f[fn].c);
/*    libaft_3d_warn("cheat8");*/
    return 1;
}

/* Quality functions */;
static REAL tqual(REAL *vertex, int v1, int v2, int v3, int v4) {
#define L(a,b) (pow(dist3i(vertex,a,b),3))
    REAL v, ls, q;
    v = det3i4(vertex, v1, v2, v3, v4);
    ls = L(v1,v2) + L(v1,v3) + L(v1,v4) + L(v2,v3) + L(v2,v4) + L(v3,v4);
    q = 6.0*sqrt(2.0)*v/ls;
    return q;
#undef L
}
static int reset_qstat(mesh3d *pm) {
    pm->qmin = 1.0;
    pm->qsum = 0.0;
    pm->qcnt = 0;
    return 0;
}
static int update_qstat(mesh3d *pm, REAL q) {
    if (q < pm->qmin)  pm->qmin = q;
    pm->qsum += q;
    pm->qcnt++;
    return 0;
}
static int print_qstat(mesh3d *pm) {
    printf("Qstats: q_min = %.8le, q_avg = %.8lf\n", pm->qmin, pm->qsum/pm->qcnt);
    return 0;
}


#ifdef SMOOTH_CF
static REAL repoint(mesh3d *pm, REAL *sx, REAL r) {
    REAL s;
    int i;

    pm->nWorkArea = 0;
    if (treews(pm, sx, r*2.0) < 0)  return r;
    s = r;
    for (i=0; i<pm->nWorkArea; i++)  s += pm->f[pm->workarea[i]].r * 4.0/3.0;
    s /= (1.0 + pm->nWorkArea);
    return s;
}
#endif

static int basepoint(mesh3d *pm, int v1, int v2, int v3, REAL *px, REAL *py, REAL *pz) {
    REAL a, b, c;
    REAL ex, ey, ez, r, rmin, x, y, z;
    REAL sx[3];

    a = dist3i(pm->vertex, v2, v3);
    b = dist3i(pm->vertex, v3, v1);
    c = dist3i(pm->vertex, v1, v2);

    normvec(pm->vertex, v1, v2, v3, &ex, &ey, &ez);
    r = sqrt(ex*ex + ey*ey +ez*ez);
    if (r <= 0.0)  libaft_3d_warn("aft3d.c: basepoint(): zero sized face [%d %d %d] (bad front?)", v1, v2, v3);
    ex /= r,  ey /= r,  ez /= r;
/*    a = b = c = 1.0;*/
    x = ((b*b+c*c)*pm->vertex[3*v1+0] + (a*a+c*c)*pm->vertex[3*v2+0] + (a*a+b*b)*pm->vertex[3*v3+0])/2.0/(a*a+b*b+c*c);
    y = ((b*b+c*c)*pm->vertex[3*v1+1] + (a*a+c*c)*pm->vertex[3*v2+1] + (a*a+b*b)*pm->vertex[3*v3+1])/2.0/(a*a+b*b+c*c);
    z = ((b*b+c*c)*pm->vertex[3*v1+2] + (a*a+c*c)*pm->vertex[3*v2+2] + (a*a+b*b)*pm->vertex[3*v3+2])/2.0/(a*a+b*b+c*c);
    
    if (pm->fsize) {
	rmin = pm->fsize(x, y, z) * 0.25 * sqrt(2.0/3.0);
	r = pm->fsize(x + rmin*ex, y + rmin*ey, z + rmin*z) * sqrt(2.0/3.0);
    } else {
	a = dist3i(pm->vertex, v2, v3);
	b = dist3i(pm->vertex, v3, v1);
	c = dist3i(pm->vertex, v1, v2);
	r = pow(a*a*a + b*b*b + c*c*c, 1.0/3.0) * pm->stepsize * sqrt(2.0/3.0) / pow(3.0, 1.0/3.0);
/*	r = (a + b + c) * pm->stepsize * sqrt(2.0/3.0) / 3.0;*/
    }
    sx[0] = x + r*ex,  sx[1] = y + r*ey,  sx[2] = z + r*ez;
#ifdef SMOOTH_CF
/*    printf("\t\t\t\t\t\t\t\t%lf » ", r);*/
    r += repoint(pm, sx, r),  r /= 2.0;
/*    r = repoint(pm, sx, r);
    printf("%lf \r", r),  fflush(stdout);*/
#endif
    x += r*ex;
    y += r*ey;
    z += r*ez;
    *px = ex,  *py = ey,  *pz = ez;
    return addvertex(pm, x, y, z, 1);
}

/* Findpoint function *******************************************************/;
/* Return code
 *  0  – OK
 * -1  – max nV
 * -2  – intersection detected
 * -3  – mark as bad
 * -8  – tree structure is failed
 * -10 – cheatpoint8 (not used)
 */
static int findpoint(mesh3d *pm, int fn) {
    int v1, v2, v3, vc, pn, p1, p2, p3;
    REAL ex, ey, ez, r, rmin, rv, friendness, gx, gy, gz;
    REAL sx[3];
    int i, j, v, neari;

    v1 = pm->f[fn].v[0];
    v2 = pm->f[fn].v[1];
    v3 = pm->f[fn].v[2];
    vc = basepoint(pm, v1, v2, v3, &ex, &ey, &ez);
    if (vc < 0) {
	libaft_3d_stop("AFT Front: maxVertex exceeded");
	return -1;
    }
    if (idet3i4(pm->vertex, v1, v2, v3, vc) != 1) libaft_3d_warn("aft3d.c: findpoint(): wrong face-orientation of new vertex (internal error)");
    r = 0.0;
    if (dist3i(pm->vertex, v1, vc) > r) r = dist3i(pm->vertex, v1, vc);
    if (dist3i(pm->vertex, v2, vc) > r) r = dist3i(pm->vertex, v2, vc);
    if (dist3i(pm->vertex, v3, vc) > r) r = dist3i(pm->vertex, v3, vc);
    r *= 1.0001220703125;

    if (1 || pm->local) {  /*force this for now*/
	pm->nWorkArea = 0;
	sx[0] = pm->vertex[3*vc+0],  sx[1] = pm->vertex[3*vc+1],  sx[2] = pm->vertex[3*vc+2];
#if 0
	for (i=0; i<pm->nF; i++)
	    if (dist3(pm->f[i].x[0], pm->f[i].x[1], pm->f[i].x[2], x, y, z) <= pm->f[i].r + r)
		pm->workarea[pm->nWorkArea++] = i;
#endif
	i = treews(pm, sx, r);
	if (i < 0)  return -8;
	if (0)  printf("\t\t\t\t\t\t\t\tws: %4d, chk: %4d, [%6.2lf%%] [%6.2lf%%] \r", pm->nWorkArea, i, 100.0*pm->nWorkArea/i, 100.0*i/pm->nF),  fflush(stdout);
	rmin = 1.0*r + 0.0*pm->stepsize*pm->minrho, neari = -1, pm->nWorkSpace = 0;
/*	rmin = sqrt(0.64*r*r + pm->minrho*pm->minrho), neari = -1, pm->nWorkSpace = 0;
	rmin = (pm->beta*r + pm->minrho + sqrt((pm->beta*r - pm->minrho)*(pm->beta*r - pm->minrho) + pm->alpha))/2.0, neari = -1, pm->nWorkSpace = 0;
	rmin = pm->beta*(r + pm->minrho) + sqrt(pm->beta*pm->beta*(r - pm->minrho)*(r - pm->minrho)/4.0 + pm->alpha), neari = -1, pm->nWorkSpace = 0;*/
	for (i=0; i<pm->nWorkArea; i++) {
	    for (j=0; j<3; j++) {
		v = pm->f[pm->workarea[i]].v[j];
		if (VinWS(pm, v)) continue;
		rv = dist3i(pm->vertex, vc, v);
		if (rv > r) continue;
		if (idet3i4(pm->vertex, v1, v2, v3, v) != 1) continue;
		pm->workspace[pm->nWorkSpace++] = v;
		if (rv < rmin) {
		    neari = v;
		    rmin = rv;
		}
	    }
	}
    } else {
	pm->nWorkArea = 0;
	for (i=0; i<pm->nF; i++)
	    pm->workarea[pm->nWorkArea++] = i;
	neari = -1, pm->nWorkSpace = 0;
	for (i=0; i<pm->nWorkArea; i++) {
	    for (j=0; j<3; j++) {
		v = pm->f[pm->workarea[i]].v[j];
		if (VinWS(pm, v)) continue;
		if (idet3i4(pm->vertex, v1, v2, v3, v) != 1) continue;
		pm->workspace[pm->nWorkSpace++] = v;
	    }
	}

    }
    if (neari >= 0)  pn = neari;
    else pn = vc;
    /*	savewrkps(pm);*/
    /* friends */
    pm->curfriend[0] = -1; pm->curfriendness[0] = -1; pm->curfriendlink[0] = -1;
    pm->curfriend[1] = -1; pm->curfriendness[1] = -1; pm->curfriendlink[1] = -1;
    pm->curfriend[2] = -1; pm->curfriendness[2] = -1; pm->curfriendlink[2] = -1;
    for (i=0; i<pm->nWorkArea; i++) {
	for (j=0; j<3; j++) {
	    p1 = pm->f[pm->workarea[i]].v[(0+j)%3];
	    p2 = pm->f[pm->workarea[i]].v[(1+j)%3];
	    p3 = pm->f[pm->workarea[i]].v[(2+j)%3];
	    if ((v2==p2) && (v3==p1) /*&& (idet3i4(pm->vertex, v1, v2, v3, p3) >= 0)*/ && (p3!=v1)) {
		normvec(pm->vertex, p1, p2, p3, &gx, &gy, &gz);
		r = sqrt(gx*gx + gy*gy + gz*gz);
		gx /= r, gy /= r, gz /= r;
		friendness = ex*gx + ey*gy + ez*gz;
		if (idet3i4(pm->vertex, v1, v2, v3, p3) < 0)  friendness = 2.0 - friendness;
		if ((pm->curfriend[0]<0) || (pm->curfriendness[0]>friendness)) {
		    pm->curfriend[0] = p3;
		    pm->curfriendness[0] = friendness;
		    pm->curfriendlink[0] = pm->workarea[i];
		    pm->cfnx[0] = gx,  pm->cfny[0] = gy,  pm->cfnz[0] = gz;
		}
	    }
	    if ((v3==p2) && (v1==p1) /*&& (idet3i4(pm->vertex, v1, v2, v3, p3) >= 0)*/ && (p3!=v2)) {
		normvec(pm->vertex, p1, p2, p3, &gx, &gy, &gz);
		r = sqrt(gx*gx + gy*gy + gz*gz);
		gx /= r, gy /= r, gz /= r;
		friendness = ex*gx + ey*gy + ez*gz;
		if (idet3i4(pm->vertex, v1, v2, v3, p3) < 0)  friendness = 2.0 - friendness;
		if ((pm->curfriend[1]<0) || (pm->curfriendness[1]>friendness)) {
		    pm->curfriend[1] = p3;
		    pm->curfriendness[1] = friendness;
		    pm->curfriendlink[1] = pm->workarea[i];
		    pm->cfnx[1] = gx,  pm->cfny[1] = gy,  pm->cfnz[1] = gz;
		}
	    }
	    if ((v1==p2) && (v2==p1) /*&& (idet3i4(pm->vertex, v1, v2, v3, p3) >= 0)*/ && (p3!=v3)) {
		normvec(pm->vertex, p1, p2, p3, &gx, &gy, &gz);
		r = sqrt(gx*gx + gy*gy + gz*gz);
		gx /= r, gy /= r, gz /= r;
		friendness = ex*gx + ey*gy + ez*gz;
		if (idet3i4(pm->vertex, v1, v2, v3, p3) < 0)  friendness = 2.0 - friendness;
		if ((pm->curfriend[2]<0) || (pm->curfriendness[2]>friendness)) {
		    pm->curfriend[2] = p3;
		    pm->curfriendness[2] = friendness;
		    pm->curfriendlink[2] = pm->workarea[i];
		    pm->cfnx[2] = gx,  pm->cfny[2] = gy,  pm->cfnz[2] = gz;
		}
	    }
	}
    }
    if (pm->curfriend[0] < 0)  libaft_3d_warn("aft3d.c: findpoint(): no friend 0 (bad front!)");
    if (pm->curfriend[1] < 0)  libaft_3d_warn("aft3d.c: findpoint(): no friend 1 (bad front!)");
    if (pm->curfriend[2] < 0)  libaft_3d_warn("aft3d.c: findpoint(): no friend 2 (bad front!)");
    /*	neari = -1; rmin = r;
	for (i=0; i<3; i++) {
	v = pm->curfriend[i];
	if (v<0) continue;
	rv = dist3i(pm->vertex, vc, v);
	if (rv > r) continue;
	if (idet3i4(pm->vertex, v1, v2, v3, v) != 1) continue;
	if (rv<rmin) {
	neari = pm->curfriend[i];
	rmin = rv;
	}
	}
	if (neari>=0) {
	pn = neari;
	}*/
    /* searching... */
    pm->nBadGuys = 0, pm->nGirls = 0;
    while (check(pm, fn, pn)) {
	pm->badguys[pm->nBadGuys++] = pn;
	if (pm->intface >= 0) {
	    for (i=0; i<3; i++) {
		v = pm->f[pm->intface].v[i];
		if (v==v1) continue;
		if (v==v2) continue;
		if (v==v3) continue;
		if (!VinWS(pm, v)) continue;
		if (VinBG(pm, v)) continue;
		if (VinG(pm, v)) continue;
		pm->girls[pm->nGirls++] = v;
	    }
	}
	if (pm->nGirls>0) {
	    neari = 0; rmin = dist3i(pm->vertex, vc, pm->girls[0]);
	    for (i=1; i<pm->nGirls; i++) {
		rv = dist3i(pm->vertex, vc, pm->girls[i]);
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
		if (pm->allowlocal) {
		    pm->local = 0;
		    vc = findpoint(pm, fn);
		    pm->local = 1;
		    return vc;
		} else return -3;
	    } else {
		/*libaft_3d_warn("aft3d.c: findpoint(): starting bruteforce");*/
		for (i=0; i<pm->nWorkSpace; i++) {
		    v = pm->workspace[i];
		    if (v==v1) continue;
		    if (v==v2) continue;
		    if (v==v3) continue;
		    if (VinBG(pm, v)) continue;
		    if (VinG(pm, v)) continue;
		    pm->girls[pm->nGirls++] = v;
		}
		if (pm->nGirls>0) {
		    neari = 0; rmin = dist3i(pm->vertex, vc, pm->girls[0]);
		    for (i=1; i<pm->nGirls; i++) {
			rv = dist3i(pm->vertex, vc, pm->girls[i]);
			if (rmin > rv) {
			    neari = i;
			    rmin = rv;
			}
		    }
		    pn = pm->girls[neari];
		    pm->nGirls--;
		    pm->girls[neari] = pm->girls[pm->nGirls];
		} else {
#ifdef INTFRONT
		    if (intfront(pm)) {
			libaft_3d_warn("AFT Front: Intersection occured in process");
			return -2;
		    }
#endif
		    pm->nV--;
		    if (pm->cheat) { /*currently forced, i.e. cheatpoint8 is NOT used anymore*/
			/*pm->cheat = 0;*/
			return -3;
		    } else {
			if ((pm->curfriend[0]>=0)&&(pm->curfriend[1]>=0)&&(pm->curfriend[2]>=0)) {
			    if (cheatpoint8(pm, v1, v2, v3, pm->curfriend[0], pm->curfriend[1], pm->curfriend[2],
					fn, pm->curfriendlink[0], pm->curfriendlink[1], pm->curfriendlink[2])) {
#ifdef INTFRONT
				if (intfront(pm)) {
				    libaft_3d_warn("AFT Front: Intersection occured in process (cheat_point)");
				    return -2;
				}
#endif
				return -10;
			    }
			    return -3;
			} else {
			    return -3;
			}
		    }
		}
	    }
	}
    }
    if (pn!=vc) {
	pm->nV--;
    } else {
	pm->statsaft++;
    }
    return pn;
}

/* Cheat functions **********************************************************/;

#ifdef CHEATPOINT
static int cheat_point(mesh3d *pm) { // obsolete
    int i, j, k, m, v1, v2, v3, v4, cp, rr, color, w1, w2, fn;

    libaft_3d_warn("aft_front: cheat_point() is obsolete, do not use it!");

#ifdef MEGACHEAT
    REAL cx1, cx2, cy1, cy2, cz1, cz2, mass;
#endif
    for (i=0; i<pm->nF; i++) {
	color = pm->f[i].c;
	for (j=0; j<3; j++) {
	    v1 = pm->f[i].v[(j+0)%3];
	    v2 = pm->f[i].v[(j+1)%3];
	    v3 = pm->f[i].v[(j+2)%3];
	    pm->nWorkArea = 0; v4 = -1;
	    for (k=0; k<pm->nF; k++) {
		if (    ((pm->f[k].v[0]==v2) || (pm->f[k].v[0]==v3)) ||
			((pm->f[k].v[1]==v2) || (pm->f[k].v[1]==v3)) ||
			((pm->f[k].v[2]==v2) || (pm->f[k].v[2]==v3))  ) {
		    pm->workarea[pm->nWorkArea] = k;
		    pm->nWorkArea++;
		}
		if (v4 >= 0) continue;
		if ((pm->f[k].v[0]==v3) || (pm->f[k].v[1]==v2)) v4=pm->f[k].v[2];
		else if ((pm->f[k].v[1]==v3) || (pm->f[k].v[2]==v2)) v4=pm->f[k].v[0];
		else if ((pm->f[k].v[2]==v3) || (pm->f[k].v[0]==v2)) v4=pm->f[k].v[1];
		else continue;
	    }
	    if (v4 < 0) {
		libaft_3d_warn("aft3d.c: cheat_point(): internal check failed");
		continue;
	    }
#ifdef MEGACHEAT
	    cx1 = (pm->vertex[3*v1+0]+pm->vertex[3*v4+0])/2.0;
	    cx2 = (pm->vertex[3*v2+0]+pm->vertex[3*v3+0])/2.0;
	    cy1 = (pm->vertex[3*v1+1]+pm->vertex[3*v4+1])/2.0;
	    cy2 = (pm->vertex[3*v2+1]+pm->vertex[3*v3+1])/2.0;
	    cz1 = (pm->vertex[3*v1+2]+pm->vertex[3*v4+2])/2.0;
	    cz2 = (pm->vertex[3*v2+2]+pm->vertex[3*v3+2])/2.0;
	    mass = 2.0;
	    do {
		cp = addvertex(pm, (cx1 + (mass-1.0)*cx2)/mass, (cy1 + (mass-1.0)*cy2)/mass, (cz1 + (mass-1.0)*cz2)/mass, 2);
		if (cp < 0)  return -1;
		rr = 1;
		for (k=0; (rr)&&(k<pm->nWorkArea); k++) {
		    if (idet3i4(pm->vertex, pm->f[pm->workarea[k]].v[0], pm->f[pm->workarea[k]].v[1], pm->f[pm->workarea[k]].v[2], cp) < 1)
			rr = 0;
		}
		for (k=0; (rr)&&(k<pm->nWorkArea); k++) {
		    for (m=0; (rr)&&(m<pm->nF); m++) {
			if (intsect(pm, pm->f[pm->workarea[k]].v[0], pm->f[pm->workarea[k]].v[1], pm->f[pm->workarea[k]].v[2], cp,
				    pm->f[m].v[0], pm->f[m].v[1], pm->f[m].v[2]))
			    rr = 0;
		    }
		}
		if (rr == 0) {
		    pm->nV--;
		    mass *= 2.0;
		    if (mass>1e30) {
			//libaft_3d_warn("mass");
			break;
		    }
		}
	    } while (rr == 0);
#else
	    cp = addvertex(pm, 
		    (pm->vertex[3*v1+0]+pm->vertex[3*v2+0]+pm->vertex[3*v3+0]+pm->vertex[3*v4+0])/4.0,
		    (pm->vertex[3*v1+1]+pm->vertex[3*v2+1]+pm->vertex[3*v3+1]+pm->vertex[3*v4+1])/4.0,
		    (pm->vertex[3*v1+2]+pm->vertex[3*v2+2]+pm->vertex[3*v3+2]+pm->vertex[3*v4+2])/4.0, 1);
	    if (cp < 0)  return -1;
	    rr = 1;
	    for (k=0; (rr)&&(k<pm->nWorkArea); k++) {
		if (idet3i4(pm->vertex, pm->f[pm->workarea[k]].v[0], pm->f[pm->workarea[k]].v[1], pm->f[pm->workarea[k]].v[2], cp) < 1)
		    rr = 0;
	    }
	    for (k=0; (rr)&&(k<pm->nWorkArea); k++) {
		for (m=0; (rr)&&(m<pm->nF); m++) {
		    if (intsect(pm, pm->f[pm->workarea[k]].v[0], pm->f[3*pm->workarea[k]].v[1], pm->f[pm->workarea[k]].v[2], cp,
				pm->f[m].v[0], pm->f[m].v[1], pm->f[m].v[2]))
			rr = 0;
		}
	    }
#endif
	    if (rr) {
		pm->statschpnt++;
/*#ifdef MEGACHEAT
		libaft_3d_warn("Inserting cheat point (mass=%lf)", mass);
#endif*/
		for (k=0; k<pm->nWorkArea; k++) {
		    if (addtetra(pm, pm->f[pm->workarea[k]].v[0], pm->f[pm->workarea[k]].v[1], pm->f[pm->workarea[k]].v[2], cp, color) < 0)  return -1;
		    if ((pm->f[pm->workarea[k]].v[0]==v2)||(pm->f[pm->workarea[k]].v[0]==v3)) {
			w1 = pm->f[pm->workarea[k]].v[1];
			w2 = pm->f[pm->workarea[k]].v[2];
		    } else if ((pm->f[pm->workarea[k]].v[1]==v2)||(pm->f[pm->workarea[k]].v[1]==v3)) {
			w1 = pm->f[pm->workarea[k]].v[2];
			w2 = pm->f[pm->workarea[k]].v[0];
		    } else if ((pm->f[pm->workarea[k]].v[2]==v2)||(pm->f[pm->workarea[k]].v[2]==v3)) {
			w1 = pm->f[pm->workarea[k]].v[0];
			w2 = pm->f[pm->workarea[k]].v[1];
		    } else continue;
		    if ( (w1==v2) || (w1==v3) || (w2==v2) || (w2==v3) ) continue;
		    fn = findface(pm, cp, w2, w1);
		    if (fn>=0)  remface(pm, fn, color);
		    else  addface(pm, cp, w1, w2, color);
		}
		for (k=0; k<pm->nWorkArea; k++)
		    remface(pm, pm->workarea[k], color);
		packface(pm);
		return 1;
	    } else {
#ifndef MEGACHEAT
		pm->nV--;
#endif
	    }
	}
    }
    return 0;
}
#endif // CHEATPOINT
#ifdef REMOVEPOINTS
static int remove_points(mesh3d *pm) { // obsolete
    int i, j, k, m, n;
    int v1, v2, p, q1, q2, q3, w1, w2, w3;
    REAL x, y, z;

    libaft_3d_warn("aft_front: remove_points() is obsolete, do not use it!");

    for (i=0; i<pm->nF; i++) {
	for (j=0; j<3; j++) {
	    v1 = pm->f[i].v[(j+1)%3];
	    v2 = pm->f[i].v[(j+2)%3];
	    if (v2 < pm->nVF) continue;
	    x = pm->vertex[3*v2+0];
	    y = pm->vertex[3*v2+1];
	    z = pm->vertex[3*v2+2];
	    pm->vertex[3*v2+0] = pm->vertex[3*v1+0];
	    pm->vertex[3*v2+1] = pm->vertex[3*v1+1];
	    pm->vertex[3*v2+2] = pm->vertex[3*v1+2];
	    /*libaft_3d_warn("try");*/
	    p = 1;
	    for (k=0; k<pm->nT; k++) {
		if (   ((pm->tetra[4*k+0]==v2) || (pm->tetra[4*k+1]==v2) || (pm->tetra[4*k+2]==v2) || (pm->tetra[4*k+3]==v2)) &&
			!((pm->tetra[4*k+0]==v1) || (pm->tetra[4*k+1]==v1) || (pm->tetra[4*k+2]==v1) || (pm->tetra[4*k+3]==v1))   ) {
		    if (idet3i4(pm->vertex, pm->tetra[4*k+0], pm->tetra[4*k+1], pm->tetra[4*k+2], pm->tetra[4*k+3])<1) {
			p = 0;
			break;
		    }
		}
	    }
	    if (p) {
		for (k=0; k<pm->nF; k++) {
		    if (     ((pm->f[k].v[0]==v2) || (pm->f[k].v[1]==v2) || (pm->f[k].v[2]==v2)) &&
			    !((pm->f[k].v[0]==v1) || (pm->f[k].v[1]==v1) || (pm->f[k].v[2]==v1))  ) {
			q1 = pm->f[k].v[0];
			q2 = pm->f[k].v[1];
			q3 = pm->f[k].v[2];
			for (m=0; m<pm->nF; m++) {
			    w1 = pm->f[m].v[0];
			    w2 = pm->f[m].v[1];
			    w3 = pm->f[m].v[2];
			    if ((w1==v1) || (w2==v1) || (w3==v1))  continue;
			    if (q1==v2)  q1=v1;
			    if (q2==v2)  q2=v1;
			    if (q3==v2)  q3=v1;
			    if (w1==v2)  w1=v1;
			    if (w2==v2)  w2=v1;
			    if (w3==v2)  w3=v1;
			    if (intface(pm, q1, q2, q3, w1, w2, w3)) {
				p = 0;
				break;
			    }
			}
		    }
		}
	    }
	    if (p) {
		for (k=0; k<pm->nF; k++) {
		    if (  ((pm->f[k].v[0]==v1) || (pm->f[k].v[1]==v1) || (pm->f[k].v[2]==v1)) &&
			    ((pm->f[k].v[0]==v2) || (pm->f[k].v[1]==v2) || (pm->f[k].v[2]==v2))  ) {
			remface(pm, k, pm->f[k].c);
		    }
		}
		for (k=0; k<pm->nT; k++) {
		    if (  ((pm->tetra[4*k+0]==v1) || (pm->tetra[4*k+1]==v1) || (pm->tetra[4*k+2]==v1) || (pm->tetra[4*k+3]==v1)) &&
			    ((pm->tetra[4*k+0]==v2) || (pm->tetra[4*k+1]==v2) || (pm->tetra[4*k+2]==v2) || (pm->tetra[4*k+3]==v2))  ) {
			remtetra(pm, k, pm->tetracolor[k]);
		    }
		}
		packface(pm);
		packtetra(pm);
		for (k=0; k<pm->nF; k++) {
		    for (m=0; m<3; m++) {
			if (pm->f[k].v[m]!=v2) continue;
			pm->f[k].v[m] = v1;
			q1 = v1;
			q2 = pm->f[k].v[(m+1)%3];
			q3 = pm->f[k].v[(m+2)%3];
			for (n=0; n<pm->nF; n++) {
			    w1 = pm->f[n].v[0];
			    w2 = pm->f[n].v[1];
			    w3 = pm->f[n].v[2];
			    if (  ((q3==w1)&&(q2==w2)&&(q1==w3)) ||
				    ((q2==w1)&&(q1==w2)&&(q3==w3)) ||
				    ((q1==w1)&&(q3==w2)&&(q2==w3))  ) {
				remface(pm, k, pm->f[k].c);
				remface(pm, n, pm->f[n].c);
				break;
			    }
			}
		    }
		    fmetric(pm, k);
		    siftdown(pm, pm->f[k].h);
		}
		packface(pm);
		/*for (k=0; k<pm->nF; k++) libaft_3d_warn("face %d: %d %d %d", k, pm->f[k].v[0], pm->f[k].v[1], pm->f[k].v[2]);*/
		for (k=0; k<4*pm->nT; k++)
		    if (pm->tetra[k]==v2) pm->tetra[k]=v1;
		/*libaft_3d_warn("Experimental: Point removed ;-)");*/
		return 1;
	    } else {
		pm->vertex[3*v2+0] = x;
		pm->vertex[3*v2+1] = y;
		pm->vertex[3*v2+2] = z;
	    }
	}

    }
    return 0;
}
#endif // REMOVEPOINTS

/* Front function ***********************************************************/;
/* Return code
 *  0 – OK
 *  1 – max nV
 *  2 – max nF
 *  3 – max nT
 * -2 – intersection detected
 * -3 – try ugly
 * -4 – volume check failed
 * -8 – tree structure is failed
 */
static int front(mesh3d *pm) {
    int i, j, fn, pn, v[3], nsteps = 0;
    REAL ss;
    if (0) savefrtraw(pm);
    pm->allowlocal = 0; pm->applied = 0;
    while (pm->nF > 0) {
	/* Reserve elements for one step */
	if (pm->nV >= pm->nnV-3) { /* one vertex + two for intersection test */
	    libaft_3d_stop("AFT Front: maxVertex exceeded");
	    return 1; /* max nV */
	}
	if (pm->nF >= pm->nnF-3) { /* three faces */
	    libaft_3d_stop("AFT Front: maxFace exceeded");
	    return 2; /* max nF */
	}
	if (pm->nT >= pm->nnT-1) { /* one tetrahedron */
	    libaft_3d_stop("AFT Front: maxTetra exceeded");
	    return 3; /* maxEdge */
	}
	if (pm->meshedvolume > (1025.0/1024.0)*pm->volume) {
	    libaft_3d_warn("AFT Front: meshed volume is greater then initial volume (internal error or bad front)");
	    return -4; /* volume check failed */
	}
	trace("front step\n");
	/*trace("initial topology: %d\n", check_surface_topology(pm->nV, pm->vertex, pm->nF, pm->face));*/
	fn = pm->heap[0];
	ss = pm->f[fn].s;
	if (ss >= pm->maxarea)  fn = -1;
	if (fn == -1) {
#ifdef INTFRONT
	    if (intfront(pm)) {
		libaft_3d_warn("AFT Front: Intersection occured in process (before)");
		return -2; /* intersection detected */
	    }
#endif
	    if (pm->allowlocal) {
		/*pm->applied    = 0; // FIXME
		pm->allowlocal = 0; // FIXME*/
		if (pm->applied) {
		    pm->applied = 0;
		    pm->robustness *= 2.0;
		    if (pm->robustness > RBST_MAX)  pm->robustness = RBST_MAX;
		    /*libaft_3d_warn("front restart");*/
		    if (0)  savefrtraw(pm);
		    for (i=0; i<pm->nF; i++) {
			fmetric(pm, i);
			siftup(pm, pm->f[i].h);
		    }
		    continue;
		} else {
#ifdef CHEATPOINT
		    if (cheat_point(pm))  {
#ifdef INTFRONT
			if (intfront(pm)) {
			    libaft_3d_warn("AFT Front: Intersection occured in process (cheat_point)");
			    return -2; /* intersection detected */
			}
#endif
			continue;
		    }
#endif
#ifdef REMOVEPOINTS
		    if (remove_points(pm))  {
#ifdef INTFRONT
			if (intfront(pm)) {
			    libaft_3d_warn("AFT Front: Intersection occured in process (remove_points)");
			    return -2; /* intersection detected */
			}
#endif
			continue;
		    }
#endif
		    return -3;
		}
	    } else {
		pm->allowlocal = 1;
		for (i=0; i<pm->nF; i++) {
		    fmetric(pm, i);
		    siftup(pm, pm->f[i].h);
		}
		continue;
	    }
	}
	trace("face %d\n", fn);
	pn = findpoint(pm, fn);
	trace("point %d\n", pn);
	if (pn==-3) {
	    pm->f[fn].s += pm->maxarea;
	    siftdown(pm, pm->f[fn].h);
	    continue;
	} else if (pn==-10) {
#ifdef SHOWPROGRESS
	    printf("aft kernel: nV = %5d, nF = %5d, nT = %5d, %6.2lf%%  \r", pm->nV, pm->nF, pm->nT, 100.0*pm->meshedvolume/pm->volume),  fflush(stdout);
#endif
	    continue;
	} else if (pn==-1) {
	    return 1; /* max nV */
	} else if (pn<0) {
	    return pn; /* report back findpoint() error code */
	}
	v[0] = pm->f[fn].v[0];
	v[1] = pm->f[fn].v[1];
	v[2] = pm->f[fn].v[2];
	checktetra(pm, v[0], v[1], v[2], pn, 1);
	if (addtetra(pm, v[0], v[1], v[2], pn, pm->f[fn].c) < 0)  return 3; /* max nT */
	update_qstat(pm, tqual(pm->vertex, v[0], v[1], v[2], pn));
	pm->applied = 1;
	if (idet3i4(pm->vertex, v[0], v[1], v[2], pn) != 1) libaft_3d_warn("aft3d.c: front(): bad orientation of tetra (%s), [%d %d %d %d], v = %le",
		(idet3i4(pm->vertex, v[0], v[1], v[2], pn))?"negative":"nearly zero", 
		v[0], v[1], v[2], pn, 
		det3i4(pm->vertex, v[0], v[1], v[2], pn));
	pm->meshedvolume += det3i4(pm->vertex, v[0], v[1], v[2], pn) / 6.0;
	remface(pm, fn, pm->f[fn].c);
	for (i=0; i<3; i++) {
	    if ((j = findface(pm, v[(i+2)%3], v[(i+1)%3], pn))>=0) {
		/*if (pm->curfriend[i] == pn) {
		remface(pm, pm->curfriendlink[i], pm->f[fn].c);*/
		remface(pm, j, pm->f[fn].c);
	    } else {
		if (addface(pm, v[(i+1)%3], v[(i+2)%3], pn, pm->f[fn].c) < 0)  return 2; /* max nF */
#ifdef STRUCTCHECK
		if (pm->curfriend[i]>=0) {
		    if (intsect(pm, v[0], v[1], v[2], pn, v[(i+1)%3], v[(i+2)%3], pm->curfriend[i])) {
			libaft_3d_warn("aft3d.c: front(): Ooooops! ;-) (internal error)");
		    }
		    if (0) saveintraw(pm->vertex, v[0], v[1], v[2], pn, v[(i+1)%3], v[(i+2)%3], pm->curfriend[i]);
		    /*savefrtraw(pm);*/
		}
#endif
	    }
	}
	if (packface(pm))  return -8; /* octree failed */
#ifdef SHOWPROGRESS
	printf("aft kernel: nV = %5d, nF = %5d, nT = %5d, %6.2lf%%  \r", pm->nV, pm->nF, pm->nT, 100.0*pm->meshedvolume/pm->volume);
	/*if (nsteps % 16 == 0)  fflush(stdout);*/
	fflush(stdout);
#endif
	trace("normal step\n");
	/*trace("front topology: %d\n", check_surface_topology(pm->nV, pm->vertex, pm->nF, pm->face));*/
	nsteps++;
#ifdef DUMPFRONT
	if (nsteps % 1 == 0)  trace("front dump %d\n", savefrtraw(pm));
#endif
#ifdef INTFRONT
	if (nsteps % 100000 == 0) {
	    if (intfront(pm)) {
		trace("front dump %d\n", savefrtraw(pm));
		libaft_3d_warn("AFT Front: Intersection occured in process (periodic check)");
		return -2; /* intersection detected */
	    }
	}
#endif
	/*for (i=0; i<pm->nF; i++)  if (pm->f[pm->heap[i]].s < pm->f[pm->heap[((i-1)/2)]].s)  printf("\nheap failed\n");
	printf("\nheap:\n");
	for (i=0; i<pm->nF; i++)  printf("%3d [%d]: %.10lf (%d) %c\n", i+1, pm->heap[i],  pm->f[pm->heap[i]].s, (i-1)/2+1, (pm->f[pm->heap[i]].s < pm->f[pm->heap[((i-1)/2)]].s)?'!':' ');*/
    }
#ifdef SHOWPROGRESS
    printf("aft kernel: nV = %5d, nF = %5d, nT = %5d,  done.                                                    \r", pm->nV, pm->nF, pm->nT),  fflush(stdout);
#endif
    return 0; /* OK */
}

/* Main function ************************************************************/;
/* return codes: 
*  0 -- OK
*  1 -- max_nV exceeded
*  2 -- max_nF exceeded
*  3 -- max_nT exceeded
*  4 -- wrong topology
*  8 -- front intersection
* -1 -- internal error (memory allocation fail)
* -2 –- intersection detected
* -4 –- volume check failed
* -8 -- tree structure is failed
*/
static int mesh3daftss_main(int *pnV, REAL *vertex,
	int *pnF, int *face, int *facecolor,
	int *pnT, int *tetra, int *tetracolor,
	int nnV, int nnF, int nnT,
	REAL ss, double (*f)(double, double, double))
{
    mesh3d mesh;
    int i, j, r = 0;
    int err = 0;
    int nF = *pnF;
    char buf[65536];
#ifdef PROF
    time_t tr1, tr2, tr;
#endif

    memset(&mesh, 0, sizeof(mesh));

    traceopen();

    detinit();

    mesh.nV = *pnV, mesh.nF = 0, mesh.nT = *pnT;
    mesh.nnV = nnV, mesh.nnF = nnF, mesh.nnT = nnT;
    mesh.vertex = vertex, mesh.tetra = tetra, mesh.tetracolor = tetracolor;

    mesh.f          = libaft_malloc(sizeof(face3d) * mesh.nnF);
    mesh.heap       = libaft_malloc(sizeof(int) * mesh.nnF);
    mesh.workarea   = libaft_malloc(sizeof(int) * mesh.nnF);
    mesh.workspace  = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.badguys    = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.girls      = libaft_malloc(sizeof(int) * mesh.nnV);
    mesh.pack       = libaft_malloc(sizeof(int) * mesh.nnF);
    mesh.tpack      = libaft_malloc(sizeof(int) * mesh.nnT);
    mesh.color      = libaft_malloc(sizeof(int) * nF*2);
    
    if (!(mesh.f && mesh.heap && mesh.workarea && mesh.workspace && mesh.badguys && mesh.girls && mesh.pack && mesh.tpack && mesh.color)) {
	libaft_free(mesh.f);
	libaft_free(mesh.heap);
	libaft_free(mesh.workarea);
	libaft_free(mesh.workspace);
	libaft_free(mesh.badguys);
	libaft_free(mesh.girls);
	libaft_free(mesh.pack);
	libaft_free(mesh.tpack);
	libaft_free(mesh.color);
	return -1; /* memory allocation fail */
    }

    mesh.npack  = 0;
    mesh.ntpack = 0;
    mesh.ncolor = 0;
    mesh.local  = 1;
    mesh.cheat  = 1;

    mesh.statsaft   = 0;
    mesh.statschp8  = 0;
    mesh.statschpnt = 0;

    mesh.root = 0;

#ifdef PROF
    tr1 = clock();
#endif

#ifdef SHOWPROGRESS
    printf("aft warmup...\r");
    fflush(stdout);
#endif

    findbb(&mesh);
    trace("bbox: [%lf, %lf] x [%lf, %lf] x [%lf, %lf]\n", mesh.bb[0], mesh.bb[3], mesh.bb[1], mesh.bb[4], mesh.bb[2], mesh.bb[5]);

    mesh.maxarea = 0.0;
    for (i=0; i<nF; i++)  {
	if (addface(&mesh, face[3*i+0], face[3*i+1], face[3*i+2], facecolor[i]) < 0) {
	    r = 2; /* max nF */
	    break;
	}
	if (mesh.maxarea < mesh.f[i].s)  mesh.maxarea = mesh.f[i].s;
    }
    mesh.maxarea *= 1e+6;

    /*Here we compute signed volume of the domain
    Front should be oriented and closed
    This formula could be found at http://www.cgafaq.info/wiki/Polyhedron_Volume*/
    for (i=0, mesh.volume = 0.0; i<mesh.nF; i++)
	mesh.volume += det3i4(mesh.vertex, mesh.f[i].v[0], mesh.f[i].v[1], mesh.f[i].v[2], 0) / 6.0;

    mesh.meshedvolume = 0.0;

    mesh.minrho = dist3i(mesh.vertex, mesh.f[0].v[0], mesh.f[0].v[1]);

    for (i=0; i<mesh.nF; i++) {
	for (j=0; j<3; j++) {
	    if (mesh.minrho > dist3i(mesh.vertex, mesh.f[i].v[j%3], mesh.f[i].v[(j+1)%3]))
		mesh.minrho = dist3i(mesh.vertex, mesh.f[i].v[j%3], mesh.f[i].v[(j+1)%3]);
	}
    }
    mesh.beta = 0.9;
    mesh.alpha = mesh.minrho*(1.0-mesh.beta)/2.0;
    mesh.alpha *= mesh.alpha;

#ifdef SHOWPROGRESS
    printf("mesh.volume = %lf\n", mesh.volume);
    printf("mesh.minrho = %lf\n", mesh.minrho);
    printf("aft kernel: nV = %5d, nF = %5d, nT = %5d, init...\r", mesh.nV, mesh.nF, mesh.nT);
    fflush(stdout);
#endif

    trace("nV = %5d, nF = %5d, nT = %5d\n", mesh.nV, mesh.nF, mesh.nT);

    mesh.stepsize = ss;
    mesh.fsize    = f;
    mesh.robustness = RBST_MIN;

    trace("stepsize = %lf\n", mesh.stepsize);

#ifdef EXTRADEBUG
    printf("\nheap:\n");
    for (i=0; i<mesh.nF; i++)  printf("%3d: %lf (%d) %c\n", i+1,  mesh.f[i].s, (i-1)/2+1, (mesh.f[i].s < mesh.f[((i-1)/2)].s)?'!':' ');
#endif

    mesh.nVF = mesh.nV;
    
    reset_qstat(&mesh);

    if (check2dsurface(&mesh.nV, mesh.vertex, &mesh.nF, face, 0, 0, 0)) {
	libaft_3d_warn("AFT Front: Bad initial front topology");
	r |= 0x4;
    }
    if (intfront(&mesh)) {
	libaft_3d_warn("AFT Front: Intersection in initial boundary");
	if (INITIAL_CHECK_FATAL)  r |= 0x8;
    }
#ifdef PROF
    tr2 = clock();
    tr = tr2-tr1;
#ifdef SHOWPROGRESS
    printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr2-tr1)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#endif
#endif
    if (r) {
	libaft_3d_warn("AFT Front: Bad initial front");
	for (i=0; i<mesh.nF; i++)  treedel(&mesh, i);
    } else {
	r = front(&mesh);
	trace("front() returned %d\n", r);
	
#ifdef PROF
	tr1 = clock();
	tr += tr1-tr2;
#ifdef SHOWPROGRESS
	printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr1-tr2)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#else
	printf("AFT kernel time: %7.2lfs\n", ((double)tr1-tr2)/CLOCKS_PER_SEC);
	printf("nT = %d\n", mesh.nT);
	printf("time per tetra = %.3le\n", ((double)tr1-tr2)/CLOCKS_PER_SEC/mesh.nT);
#endif
#endif
	
	for (i=0; i<mesh.nF; i++)  face[3*i+0] = mesh.f[i].v[0],  face[3*i+1] = mesh.f[i].v[1],  face[3*i+2] = mesh.f[i].v[2],  facecolor[i]=mesh.f[i].c;
	for (i=0; i<mesh.nF; i++)  if (treedel(&mesh, i))  err++;
	if (err)  r = -8; /* octree fail */
	fprintf(stderr, "\nAFTINFO: nf = %d\n", mesh.nF);
	printf("\n");
	print_qstat(&mesh);
	/*print_aft_stats(mesh.vertex, mesh.nT, mesh.tetra);*/
	print_aft_stats_mba(mesh.vertex, mesh.nT, mesh.tetra);
	if (r==-3) {
#ifdef TRYUGLY
	    libaft_3d_warn("AFT Front: stopped, starting second method");
	    /*savefrtraw(&mesh);
	    refine3daft(&mesh.nV, mesh.vertex, &mesh.nVF, &mesh.nT, mesh.tetra);
	    savefrtraw(&mesh);*/
	    r = mesh_3d_ugly(&mesh.nV, mesh.vertex, &mesh.nF, face, facecolor, &mesh.nT, mesh.tetra, mesh.tetracolor, mesh.nnV, mesh.nnF, mesh.nnT);
#ifdef PROF
	    tr2 = clock();
	    tr += tr2-tr1;
#ifdef SHOWPROGRESS
	    printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr2-tr1)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#else
	    printf("second method time: %7.2lfs\n", ((double)tr2-tr1)/CLOCKS_PER_SEC);
#endif
#endif
	    if (r)  libaft_3d_warn("AFT Front: second method failed");
#else
	    libaft_3d_warn("AFT Front: stopped");
#endif
	}
    }
    
#ifdef STATS
    libaft_3d_warn("AFT Front: Vertex stats: AFT: %d, CHP8: %d, CHPNT: %d.",
	    mesh.statsaft, mesh.statschp8, mesh.statschpnt);
#endif
    trace("Vertex stats: AFT: %d, CHP8: %d, CHPNT: %d.\n", mesh.statsaft, mesh.statschp8, mesh.statschpnt);

    /*	colordump(&mesh);
	for (i=0; i<mesh.ncolor; i++) {
	printf(" %d", pcolor(&mesh, mesh.color[2*i+0]));
	}
	printf("\n");
	*/
    for (i=0; i<mesh.nT; i++) {
	mesh.tetracolor[i] = pcolor(&mesh, mesh.tetracolor[i]);
    }

    libaft_free(mesh.f);
    libaft_free(mesh.heap);
    libaft_free(mesh.workarea);
    libaft_free(mesh.workspace);
    libaft_free(mesh.badguys);
    libaft_free(mesh.girls);
    libaft_free(mesh.pack);
    libaft_free(mesh.tpack);
    libaft_free(mesh.color);

    detstatsprint(buf);
    trace("%s\n", buf);

    printf("\n");
    /*print_aft_stats(mesh.vertex, mesh.nT, mesh.tetra);*/
    print_aft_stats_mba(mesh.vertex, mesh.nT, mesh.tetra);

    refine3daft(&mesh.nV, mesh.vertex, &mesh.nVF, &mesh.nT, mesh.tetra);
    /*refine3daftfunc(&mesh.nV, mesh.vertex, &mesh.nVF, &mesh.nT, mesh.tetra);*/
#ifdef PROF
    tr1 = clock();
    tr += tr1-tr2;
#ifdef SHOWPROGRESS
    printf("\r\t\t\t\t\t\t\t\t[%7.2lfs, %7.2lfs]\n", ((double)tr1-tr2)/CLOCKS_PER_SEC, (double)tr/CLOCKS_PER_SEC);
#else
    printf("smoothing time: %7.2lfs\n", ((double)tr1-tr2)/CLOCKS_PER_SEC);
#endif
#endif

    traceclose();

#ifdef SHOWPROGRESS
    printf("                                                        \r");
    fflush(stdout);
#endif
    printf("\n");
    print_aft_stats(mesh.vertex, mesh.nT, mesh.tetra);
    print_aft_stats_mba(mesh.vertex, mesh.nT, mesh.tetra);
    
    *pnV = mesh.nV, *pnF = mesh.nF, *pnT = mesh.nT;
    return r;
}
/* Wrappers *****************************************************************/;
int mesh3daftss(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, REAL ss) {
    return mesh3daftss_main(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, ss, 0);
}
int mesh3daft(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT) {
    return mesh3daftss_main(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, 1.125, 0);
}
int mesh3daftuser(int *pnV, REAL *vertex, int *pnF, int *face, int *facecolor, int *pnT, int *tetra, int *tetracolor, int nnV, int nnF, int nnT, double (*f)(double, double, double)) {
    return mesh3daftss_main(pnV, vertex, pnF, face, facecolor, pnT, tetra, tetracolor, nnV, nnF, nnT, 1.125, f);
}


/* Debug functions **********************************************************/;
#ifdef STRUCTCHECK
static int saveintraw(REAL *vertex, int a, int b, int c, int d, int u, int v, int w) {
    FILE *f;
    static int fn=0;
    char buf[1024];

    sprintf(buf, "int%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "7 5\n");
    fprintf(f, "%lf %lf %lf\n", vertex[3*a+0], vertex[3*a+1], vertex[3*a+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*b+0], vertex[3*b+1], vertex[3*b+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*c+0], vertex[3*c+1], vertex[3*c+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*d+0], vertex[3*d+1], vertex[3*d+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*u+0], vertex[3*u+1], vertex[3*u+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*v+0], vertex[3*v+1], vertex[3*v+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*w+0], vertex[3*w+1], vertex[3*w+2]);
    fprintf(f, "1 2 3\n");
    fprintf(f, "2 3 4\n");
    fprintf(f, "3 4 1\n");
    fprintf(f, "4 1 2\n");
    fprintf(f, "5 6 7\n");
    fclose(f);
    return fn++;
}
#endif
static int savefrtraw(mesh3d *pm) {
    FILE *f;
    int i;
    static int fn=0;
    char buf[1024];

    sprintf(buf, "frt%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "%d %d\n", pm->nV, pm->nF);
    for (i=0; i<pm->nV; i++) {
	fprintf(f, "%20.16le %20.16le %20.16le\n", pm->vertex[3*i+0], pm->vertex[3*i+1], pm->vertex[3*i+2]);
    }
    for (i=0; i<pm->nF; i++) {
	fprintf(f, "%d %d %d\n", pm->f[i].v[0]+1, pm->f[i].v[1]+1, pm->f[i].v[2]+1);
    }
    fclose(f);
    return fn++;
}
static int write_front_int_gmv(char *fn, int nV, double *vertex, int nF, face3d *f) {
    int i,mmax=0, unkmat=0;
    FILE *file;
    file=fopen(fn, "w");
    fprintf(file, "gmvinput ascii\n\nnodev %5d\n", nV);
    for (i=0; i<nV; i++) {
	fprintf(file, "  %20.15lf %20.15lf %20.15lf\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2]);
    }
    fprintf(file, "\ncells %5d\n", nF);
    for (i=0; i<nF; i++) {
	fprintf(file, " tri 3\n  %4d %4d %4d\n", f[i].v[0]+1, f[i].v[1]+1, f[i].v[2]+1);
	if (mmax<f[i].c) mmax = f[i].c;
	if (f[i].c==0) unkmat = 1;
    }
    if (nF) {
	fprintf(file, "\nmaterial %d 0\n", mmax+unkmat+1);
	for (i=0; i<mmax; i++) fprintf(file,"mat%d\n", i+1);
	if (unkmat)	fprintf(file,"unknown\n");
	fprintf(file,"intersect\n");
	for (i=0; i<nF; i++) {
	    fprintf(file, " %2d", (f[i].c>0)?f[i].c:(f[i].c==-256)?mmax+unkmat+1:mmax+unkmat);
	}
    }

    fprintf(file, "\nendgmv");
    fclose(file);
    return 0;
}
