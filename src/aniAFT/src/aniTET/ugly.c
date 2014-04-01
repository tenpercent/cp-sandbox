#include <math.h>
#include <string.h>
#include <stdio.h>
#include "aft.h"
#include "ugly.h"
#include "delaunay.h"
#include "helper.h"
#include "check.h"
#include "det.h"
#include "libaft.h"

//#define STRUCTCHECK
//#define DUMPFRONT
//#define FULLDEBUG

#define DUMP 0

typedef struct {
    int nV, nnV, lnV;
    REAL *vertex;
    int nF, nnF, lnF;
    int *face, *facecolor;
    int *dupface, ndupF;
    int nE, nnE;
    int *edge;
    int nT, nnT, lnT, mnT;
    int *tetra, *tetracolor;
    int npack, *pack;
    int ntpack, *tpack;
    int ncolor, *color;
    int nP, nnP;
    int *pair;
    int *vtrace;
    int nD, nnD, lnD;
    int *deck, *dused;
    int nbnd, *bnd;
    int nH, *hide, nHedge, nHface;
    int *wingstart, *wingnext, nW, nnW;
    char *wingdir;
    REAL *wingnorm;
} mesh3d;



/* debug functions (at the end of the file) */
static int gmvwrite(mesh3d *pm, int flag);
static int save_ugly_int(REAL *vertex, int a, int b, int c, int d, int x);
static int savefrtraw(mesh3d *pm);



/* Metric functions */
//#define SIGN_EPS 1.0/1073741824.0/1.0
// previous value #define SIGN_EPS 1e-10
#define SIGN_EPS 1e-12
static double sign_eps_default = SIGN_EPS;
static double sign_eps_face = SIGN_EPS;
static double sign_eps_big = SIGN_EPS;
static double sign_eps = SIGN_EPS;
//static double sign_eps_flat = 1.0/1073741824.0/1.0;
static double sign_eps_flat = 1e-10;
static double edge_eps = 1e-10;
static int sign(REAL x) {
    if (x>sign_eps) return 1;
    else if (x<-sign_eps) return -1;
    else return 0;
}

static int sign_flat(REAL x) {
    if (x>sign_eps_flat) return 1;
    else if (x<-sign_eps_flat) return -1;
    else return 0;
}
static REAL dot(REAL *v1, REAL *v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}
static void vnorm(REAL *v) {
    REAL d;
    d = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    d = sqrt(d);
    if (d != 0) d = 1.0/d;
    v[0] *= d,  v[1] *= d,  v[2] *= d;
}
static REAL vcos(REAL *vertex, int v1, int v2, int v3) {
    REAL e1[3], e2[3];
    e1[0] = vertex[3*v1+0] - vertex[3*v2+0];
    e1[1] = vertex[3*v1+1] - vertex[3*v2+1];
    e1[2] = vertex[3*v1+2] - vertex[3*v2+2];
    e2[0] = vertex[3*v3+0] - vertex[3*v2+0];
    e2[1] = vertex[3*v3+1] - vertex[3*v2+1];
    e2[2] = vertex[3*v3+2] - vertex[3*v2+2];
    vnorm(e1),  vnorm(e2);
    return dot(e1, e2);
}


/* colorify functions */
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
static int cjoin(mesh3d *pm, int c1, int c2) {
    int i, p1=pcolor(pm, c1), p2=pcolor(pm, c2);
    if (p1 > p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p1)  break;
	}
	pm->color[2*i+1] = p2;
	return pcolor(pm, c1);
    } else  if (p1 < p2) {
	for (i=0; i<pm->ncolor; i++) {
	    if (pm->color[2*i+0]==p2)  break;
	}
	pm->color[2*i+1] = p1;
	return pcolor(pm, c2);
    } else  return p1;
}

/* structure functions */
static int addtetra(mesh3d *pm, int v1, int v2, int v3, int v4, int color) {
    if (idet3i4(pm->vertex, v1, v2, v3, v4) < 0)
	libaft_3d_warn("ugly.c: addtetra(): Bad orientation of new tetra: (negative), [%d %d %d %d] (%le)",
		v1, v2, v3, v4, det3i4(pm->vertex, v1, v2, v3, v4));
    color = pcolor(pm, color);
    if (pm->nT >= pm->nnT) {
	libaft_3d_stop("3D meshing: maxTetra exceeded");
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
static int addvertex(mesh3d *pm, REAL x, REAL y, REAL z) {
#ifdef STRUCTCHECK
    int i;
    for (i=0; i<pm->nV; i++) {
	if (fabs(pm->vertex[3*i+0]-x) > sign_eps_flat)  continue;
	if (fabs(pm->vertex[3*i+1]-y) > sign_eps_flat)  continue;
	if (fabs(pm->vertex[3*i+2]-z) > sign_eps_flat)  continue;
	libaft_3d_warn("ugly.c: addvertex(): duplicate points detected");
    }
#endif
    if (pm->nV >= pm->nnV) {
	libaft_3d_stop("3D meshing: maxVertex exceeded");
	return -1;
    } else {
	pm->vertex[3*(pm->nV)+0] = x;
	pm->vertex[3*(pm->nV)+1] = y;
	pm->vertex[3*(pm->nV)+2] = z;
	return (pm->nV)++;
    }
}
static int addface(mesh3d *pm, int v1, int v2, int v3, int color) {
    
#ifdef STRUCTCHECK
    int i;
    for (i=pm->lnF; i<pm->nF; i++) {
	if ((v1!=pm->face[3*i+0])&&(v1!=pm->face[3*i+1])&&(v1!=pm->face[3*i+2])) continue;
	if ((v2!=pm->face[3*i+0])&&(v2!=pm->face[3*i+1])&&(v2!=pm->face[3*i+2])) continue;
	if ((v3!=pm->face[3*i+0])&&(v3!=pm->face[3*i+1])&&(v3!=pm->face[3*i+2])) continue;
	libaft_3d_warn("ugly.c: addface(): duplicate faces detected (internal error)");
    }
#endif
    
    color = pcolor(pm, color);
    if (pm->nF >= pm->nnF) {
	libaft_3d_stop("3D meshing: maxFace exceeded");
	return -1;
    } else {
	pm->face[3*(pm->nF)+0] = v1;
	pm->face[3*(pm->nF)+1] = v2;
	pm->face[3*(pm->nF)+2] = v3;
	pm->facecolor[pm->nF] = color;
	return (pm->nF)++;
    }
}
static void remface(mesh3d *pm, int i, int color) {
    cjoin(pm, pm->facecolor[i], color);
    if (pm->npack > pm->nnF) libaft_3d_stop("ugly.c: Internal error in remface & packface");
    pm->pack[pm->npack] = i;
    pm->npack++;
}
static int findface(mesh3d *pm, int v1, int v2, int v3) {
    int i, w1, w2, w3;
    for (i=pm->lnF; i<pm->nF; i++) {
	w1 = pm->face[3*i+0];
	w2 = pm->face[3*i+1];
	w3 = pm->face[3*i+2];
	if ((v1==w1)&&(v2==w2)&&(v3==w3)) return i;
	if ((v2==w1)&&(v3==w2)&&(v1==w3)) return i;
	if ((v3==w1)&&(v1==w2)&&(v2==w3)) return i;
    }
    return -1;
}
static int findlface(mesh3d *pm, int v1, int v2, int v3) {
    int i, w1, w2, w3;
    for (i=0; i<pm->lnF; i++) {
	w1 = pm->face[3*i+0];
	w2 = pm->face[3*i+1];
	w3 = pm->face[3*i+2];
	if (    ((w1==v1)||(w1==v2)||(w1==v3)) &&
		((w2==v1)||(w2==v2)||(w2==v3)) &&
		((w3==v1)||(w3==v2)||(w3==v3))    ) return i;
    }
    return -1;
}
static void packface(mesh3d *pm) {
    int i, j;
    while (pm->npack > 0) {
	pm->npack--;
	j = pm->pack[pm->npack];
	pm->nF--;
	pm->face[3*j+0] = pm->face[3*(pm->nF)+0];
	pm->face[3*j+1] = pm->face[3*(pm->nF)+1];
	pm->face[3*j+2] = pm->face[3*(pm->nF)+2];
	pm->facecolor[j] = pm->facecolor[pm->nF];
	for (i=0; i<pm->npack; i++) if (pm->pack[i] == pm->nF) pm->pack[i] = j;
    }
}
static void remtetra(mesh3d *pm, int i, int color) {
    cjoin(pm, pm->tetracolor[i], color);
    if (pm->ntpack > pm->nnT) libaft_3d_stop("ugly.c: Internal error in remtetra & packtetra");
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

static int caddedge(mesh3d *pm, int v1, int v2) {
    int i;
    for (i=0; i<pm->nE; i++) {
	if (  (v1==pm->edge[2*i+0]) && (v2==pm->edge[2*i+1])  )  return i;
	if (  (v2==pm->edge[2*i+0]) && (v1==pm->edge[2*i+1])  )  return i;
    }
    pm->edge[2*pm->nE+0] = v1;
    pm->edge[2*pm->nE+1] = v2;
    return pm->nE++;
}

/**********************************************************************/
static int edge2edge(mesh3d *pm, int i0, int i1, int i2, int i3) {
    int i;
    REAL alpha[3], beta[3], gamma[3], phi0[3], tmp[3];
    REAL D;
    REAL eps = edge_eps;
    REAL X;

    for (i=0; i<3; i++)  alpha[i] = pm->vertex[3*i0+i] - pm->vertex[3*i1+i];
    for (i=0; i<3; i++)  beta[i]  = pm->vertex[3*i3+i] - pm->vertex[3*i2+i];
    for (i=0; i<3; i++)  phi0[i]  = pm->vertex[3*i0+i] + pm->vertex[3*i1+i] - pm->vertex[3*i2+i] - pm->vertex[3*i3+i];
    X = sqrt((dot(alpha, alpha) + dot(beta, beta) + dot(phi0, phi0))/3.0);
//    X = exp2((int)round(log2(X)));
    for (i=0; i<3; i++)  alpha[i] /= X;
    for (i=0; i<3; i++)  beta[i]  /= X;
    for (i=0; i<3; i++)  phi0[i]  /= X;
   
    normvec2(alpha, beta, gamma);
    
    D = dot(gamma, gamma);
    
    if (D < eps) {
	normvec2(phi0, alpha, tmp);
	if (dot(tmp, tmp) > eps)  return 0;
	normvec2(phi0, beta, tmp);
	if (dot(tmp, tmp) > eps)  return 0;
	tmp[0] = alpha[0] + beta[0],  tmp[1] = alpha[1] + beta[1],  tmp[2] = alpha[2] + beta[2];
	if (fabs(dot(phi0, tmp)) > dot(tmp, tmp) + eps)  return 0;
	tmp[0] = alpha[0] - beta[0],  tmp[1] = alpha[1] - beta[1],  tmp[2] = alpha[2] - beta[2];
	if (fabs(dot(phi0, tmp)) > dot(tmp, tmp) + eps)  return 0;
	return 1;
    } else {
	if (fabs(dot(phi0, gamma)) > eps)  return 0;
	if (fabs(det3(phi0, beta, gamma)) > D + eps)  return 0;
	if (fabs(det3(alpha, phi0, gamma)) > D + eps)  return 0;
	return 1;
    }
}
/**********************************************************************/
static int bufadd(int *pnx, int *x, int y) {
    int n=*pnx;
    int i;
    if (y<0)  return -1;
    for (i=0; i<n; i++) {
	if (x[i] == y)  return i;
    }
    if (n>=10)  libaft_3d_warn("ugly.c: bufadd(): nmax exceeded (internal error)");
    x[n++] = y;
    *pnx = n;
    return n;
}
static int check_face_vertex(mesh3d *pm, int v, int f) {
    int i, j, x, k;
    k = 0;
    for (i=0; i<3; i++) {
	x = pm->vtrace[3*v+i];
	if (x<0) break;
	for (j=0; j<3; j++) {
	    if (x == pm->face[3*f+j]) {
		k++;
		break;
	    }
	}
    }
    return (k < i) ? 0 : 1;
}
static int force_face(mesh3d *pm, int v1, int v2) {
    int i, j, r1, r2, p1, p2, p3, t1, t2, t3, k, kk, kkk, v, w1, w2, w3, w0;
    REAL m1, m2, m3, m, rl1, rl2;

    sign_eps = sign_eps_face;
    for (i=0; i<pm->nF; i++) {
	if (check_face_vertex(pm, v1, i))  continue;
	if (check_face_vertex(pm, v2, i))  continue;
	t1 = pm->face[3*i+0];
	t2 = pm->face[3*i+1];
	t3 = pm->face[3*i+2];
	
	m = det3i4(pm->vertex, v1, t1, t2, t3) - det3i4(pm->vertex, v2, t1, t2, t3);
	if (sign_flat(m)==0) continue;
	rl1 = det3i4(pm->vertex, v1, t1, t2, t3)/m,  r1 = sign_flat(rl1);
	rl2 = det3i4(pm->vertex, v2, t1, t2, t3)/m,  r2 = sign_flat(rl2);
	if (r1*r2 != -1) continue;
	m1 = det3i4(pm->vertex, v1, v2, t2, t3)/m,  p1 = sign(m1);
	m2 = det3i4(pm->vertex, v1, v2, t3, t1)/m,  p2 = sign(m2);
	m3 = det3i4(pm->vertex, v1, v2, t1, t2)/m,  p3 = sign(m3);
	if (r1<0) {
	    m1 = -m1,  m2 = -m2,  m3 = -m3;
	    p1 = -p1,  p2 = -p2,  p3 = -p3;
	}
	if (p1+p2+p3<-1)  {
	    if (0) save_ugly_int(pm->vertex, v1, t1, t2, t3, v2);
	    libaft_3d_warn("3D meshing: internal error");
	    libaft_3d_warn("p1=%d, p2=%d, p3=%d,  m1=%le, m2=%le, m3=%le.\n", p1, p2, p3, m1, m2, m3);
	    libaft_3d_warn("v1=%d, v2=%d,  t1=%d, t2=%d, t3=%d.\n", v1, v2, t1, t2, t3);
	    libaft_3d_warn("r1=%d (%le), r2=%d (%le).\n", r1, det3i4(pm->vertex, v1, t1, t2, t3), r2, det3i4(pm->vertex, v2, t1, t2, t3));
	    libaft_3d_warn("force_face: p1+p2+p3 < -1");
	}
	if ((p1<0) || (p2<0) || (p3<0)) continue;
	if (p1+p2+p3==3) {
	    m = m1+m2+m3;
	    v = addvertex(pm,
		    (m1*pm->vertex[3*t1+0] + m2*pm->vertex[3*t2+0] + m3*pm->vertex[3*t3+0])/m,
		    (m1*pm->vertex[3*t1+1] + m2*pm->vertex[3*t2+1] + m3*pm->vertex[3*t3+1])/m,
		    (m1*pm->vertex[3*t1+2] + m2*pm->vertex[3*t2+2] + m3*pm->vertex[3*t3+2])/m);
	    pm->vtrace[3*v+0] = t1;
	    pm->vtrace[3*v+1] = t2;
	    pm->vtrace[3*v+2] = t3;
	    pm->nP = 0;
	    for (k=pm->lnT; k<pm->nT; k++) {
		for (kk=0; kk<4; kk++) {
		    w0 = pm->tetra[4*k + kk];
		    if (w0!=v1) continue;
		    for (kkk=1; kkk<4; kkk++) {
			w1 = pm->tetra[4*k + (kk+kkk)%4];
			if (w1!=v2) continue;
			w2 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?1:0))%3)%4];
			w3 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?0:1))%3)%4];
			pm->pair[2*pm->nP + 0] = w2;
			pm->pair[2*pm->nP + 1] = w3;
			pm->nP++;
			remtetra(pm, k, 0);
			break;
		    }
		}
	    }
	    packtetra(pm);
	    if (pm->nP<2) libaft_3d_warn("ugly.c: force_face(): case3: something is wrong... (internal error)");
	    /*printf("%d pairs\n", pm->nP);*/
	    for (j=0; j<pm->nP; j++) {
		w2 = pm->pair[2*j+0];
		w3 = pm->pair[2*j+1];
		/*printf("%3d %3d\n", w2, w3);*/
		addtetra(pm, v1, v, w2, w3, 0);
		addtetra(pm, v, v2, w2, w3, 0);
//		caddedge(pm, v, v1);
//		caddedge(pm, v, v2);
//		caddedge(pm, v, w2);
//		caddedge(pm, v, w3);
	    }
//	    libaft_3d_warn("force_face: case3 found!");
	    return 1 + force_face(pm, v1, v) + force_face(pm, v, v2);
	} else {
	    libaft_3d_warn("3D meshing: internal error");
	    libaft_3d_warn("p1=%d, p2=%d, p3=%d,  m1=%le, m2=%le, m3=%le.\n", p1, p2, p3, m1, m2, m3);
	    libaft_3d_warn("v1=%d, v2=%d,  t1=%d, t2=%d, t3=%d.\n", v1, v2, t1, t2, t3);
	    libaft_3d_warn("r1=%d (%le), r2=%d (%le).\n", r1, rl1, r2, rl2);
	    libaft_3d_warn("force_face: do not trust me");
	    return -99999;
	}
    }
    return 0;
}
/**********************************************************************/
#ifdef STRUCTCHECK
static int check_valid_tetra(mesh3d *pm, int k) {
    int i, kk, v1, v2, v3, r=0;
    int buf[10], nbuf;
    for (kk=0; kk<4; kk++) {
	v1 = pm->tetra[4*k + (kk^3)];
	v2 = pm->tetra[4*k + (kk^2)];
	v3 = pm->tetra[4*k + (kk^1)];
	nbuf = 0;
	for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v1+i]);
	for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v2+i]);
	for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v3+i]);
	if (nbuf<3)  {
	    libaft_3d_warn("3D meshing: internal error");
	    libaft_3d_warn("check_valid_tetra(): tetra %d: ntrace < 3", k);
	    libaft_3d_warn("additional info:\nface: %d %d %d\nvertex %d: %d %d %d\nvertex %d: %d %d %d\nvertex %d: %d %d %d\n",
			    v1, v2, v3,
			    v1, pm->vtrace[3*v1+0],  pm->vtrace[3*v1+1], pm->vtrace[3*v1+2],
			    v2, pm->vtrace[3*v2+0],  pm->vtrace[3*v2+1], pm->vtrace[3*v2+2],
			    v3, pm->vtrace[3*v3+0],  pm->vtrace[3*v3+1], pm->vtrace[3*v3+2] );
			    r++;
	}
    }
    return r;
}
static int check_tetra_faces(mesh3d *pm) {
    int k, r = 0;
    for (k=pm->lnT; k<pm->nT; k++) {
	r += check_valid_tetra(pm, k);
    }
    return r;
}
#endif

static int force_pre_edge(mesh3d *pm, int v1, int v2, int c1, int c2);
static int force_edge_split_face(mesh3d *pm, int v1, int v2, int c1, int c2, int t1, int t2, int t3, REAL m1, REAL m2, REAL m3) {
    int v, k, kk, w0, w1, w2, w3;
    REAL m;
    m = m1+m2+m3;
    v = addvertex( pm,
		   (m1*pm->vertex[3*t1+0] + m2*pm->vertex[3*t2+0] + m3*pm->vertex[3*t3+0])/m,
		   (m1*pm->vertex[3*t1+1] + m2*pm->vertex[3*t2+1] + m3*pm->vertex[3*t3+1])/m,
		   (m1*pm->vertex[3*t1+2] + m2*pm->vertex[3*t2+2] + m3*pm->vertex[3*t3+2])/m);
    pm->vtrace[3*v+0] = c1;
    pm->vtrace[3*v+1] = c2;
    pm->vtrace[3*v+2] = -1;
    for (k=pm->lnT; k<pm->nT; k++) {
	for (kk=0; kk<4; kk++) {
	    w1 = pm->tetra[4*k + (kk+((kk%2)?3:1))%4];
	    w2 = pm->tetra[4*k + (kk+2)%4];
	    w3 = pm->tetra[4*k + (kk+((kk%2)?1:3))%4];
	    if (((t1==w3) && (t2==w2) && (t3==w1)) ||
		((t2==w3) && (t3==w2) && (t1==w1)) ||
		((t3==w3) && (t1==w2) && (t2==w1))) {
		w0 = pm->tetra[4*k + kk];
		remtetra(pm, k, 0);
		packtetra(pm);
#ifdef STRUCTCHECK
		check_valid_tetra(pm, addtetra(pm, v1, t1, t2, v, 0));
		check_valid_tetra(pm, addtetra(pm, v1, t2, t3, v, 0));
		check_valid_tetra(pm, addtetra(pm, v1, t3, t1, v, 0));
		check_valid_tetra(pm, addtetra(pm, w0, t3, t2, v, 0));
		check_valid_tetra(pm, addtetra(pm, w0, t2, t1, v, 0));
		check_valid_tetra(pm, addtetra(pm, w0, t1, t3, v, 0));
#else
		addtetra(pm, v1, t1, t2, v, 0);
		addtetra(pm, v1, t2, t3, v, 0);
		addtetra(pm, v1, t3, t1, v, 0);
		addtetra(pm, w0, t3, t2, v, 0);
		addtetra(pm, w0, t2, t1, v, 0);
		addtetra(pm, w0, t1, t3, v, 0);
#endif
		return 1 + force_pre_edge(pm, v, v2, c1, c2);
	    }
	}
    }
    libaft_3d_warn("ugly.c: force_edge_split_face(): split_face case: failed (internal error)");
    return -999999;
}
static int force_edge_split_edge(mesh3d *pm, int v1, int v2, int c1, int c2, int t1, int t2, REAL m1, REAL m2) {
    int i, k, kk, kkk, v, w0, w1, w2, w3;
    REAL m;
    (void) v1;
    
    m = m1+m2;
    v = addvertex( pm,
		   (m1*pm->vertex[3*t1+0] + m2*pm->vertex[3*t2+0])/m,
		   (m1*pm->vertex[3*t1+1] + m2*pm->vertex[3*t2+1])/m,
		   (m1*pm->vertex[3*t1+2] + m2*pm->vertex[3*t2+2])/m );
    pm->vtrace[3*v+0] = c1;
    pm->vtrace[3*v+1] = c2;
    pm->vtrace[3*v+2] = -1;
    pm->nP = 0;
    for (k=pm->lnT; k<pm->nT; k++) {
	for (kk=0; kk<4; kk++) {
	    w0 = pm->tetra[4*k + kk];
	    if (w0!=t1)  continue;
	    for (kkk=1; kkk<4; kkk++) {
		w1 = pm->tetra[4*k + (kk+kkk)%4];
		if (w1!=t2)  continue;
		w2 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?1:0))%3)%4];
		w3 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?0:1))%3)%4];
		pm->pair[2*pm->nP + 0] = w2;
		pm->pair[2*pm->nP + 1] = w3;
		pm->nP++;
		remtetra(pm, k, 0);
		break;
	    }
	}
    }
    packtetra(pm);
    if (pm->nP<2) libaft_3d_warn("ugly.c: force_edge_split_edge(): edge-to-edge case: number of tetras < 2 (internal error)");
    for (i=0; i<pm->nP; i++) {
	w2 = pm->pair[2*i+0];
	w3 = pm->pair[2*i+1];
#ifdef STRUCTCHECK
	check_valid_tetra(pm, addtetra(pm, t1, v, w2, w3, 0));
	check_valid_tetra(pm, addtetra(pm, v, t2, w2, w3, 0));
#else
	addtetra(pm, t1, v, w2, w3, 0);
	addtetra(pm, v, t2, w2, w3, 0);
#endif
    }
    return 1 + force_pre_edge(pm, v, v2, c1, c2);
}
static int force_pre_edge(mesh3d *pm, int v1, int v2, int c1, int c2) {
    int i, j, jj, r1, r2;
    int t[3];
    REAL m[3], mm, rl1, rl2;
    sign_eps = sign_eps_default;
    for (j=pm->lnT; j<pm->nT; j++) {
	for (jj=0; jj<4; jj++) {
	    if (pm->tetra[4*j+jj] == v1) {
		t[0] = pm->tetra[4*j + (jj+((jj%2)?3:1))%4];
		t[1] = pm->tetra[4*j + (jj+2)%4];
		t[2] = pm->tetra[4*j + (jj+((jj%2)?1:3))%4];
		for (i=0; i<3; i++) {
		    if (t[i]==v2)  return 0;
		}
		rl1 = det3i4(pm->vertex, v1, t[0], t[1], t[2]);
		rl2 = det3i4(pm->vertex, v2, t[0], t[1], t[2]);
		mm = rl1 - rl2;
		if (sign_flat(mm)==0)  continue;
		rl1 /= mm,  r1 = sign_flat(rl1);
		if ((r1 < 0) && (mm > 0)) libaft_3d_warn("ugly.c: force_pre_edge(): r1<=0 (v = %le)  (internal error)", rl1);
		rl2 /= mm,  r2 = sign_flat(rl2);
		if (r1*r2 > 0)  continue;
		for (i=0; i<3; i++) {
		    m[i] = det3i4(pm->vertex, v1, v2, t[(i+1)%3], t[(i+2)%3])/mm;
		}
		for (i=0; i<3; i++) {
		    if (edge2edge(pm, v1, v2, t[i], t[(i+1)%3]))  return force_edge_split_edge(pm, v1, v2, c1, c2, t[i], t[(i+1)%3], m[i], m[(i+1)%3]);
		}
		if ((m[0] > 0.0) && (m[1] > 0.0) && (m[2] > 0.0)) {
		    remtetra(pm, j, 0);
		    return force_edge_split_face(pm, v1, v2, c1, c2, t[0], t[1], t[2], m[0], m[1], m[2]);
		}
		if (m[0] < -sign_eps)  continue;
		if (m[1] < -sign_eps)  continue;
		if (m[2] < -sign_eps)  continue;
		if (0) save_ugly_int(pm->vertex, v1, t[0], t[1], t[2], v2);
		libaft_3d_warn("3D meshing: internal error");
		libaft_3d_warn("m1=%le, m2=%le, m3=%le.\n", m[0], m[1], m[2]);
		libaft_3d_warn("v1=%d, v2=%d,  t1=%d, t2=%d, t3=%d.\n", v1, v2, t[0], t[1], t[2]);
		libaft_3d_warn("r1=%d (%le), r2=%d (%le).\n", r1, rl1, r2, rl2);
		libaft_3d_warn("force_pre_edge: case1, the impossible one :-O");
		}
	    }
	}
    return 0;
}
/**********************************************************************/
static int force_edge(mesh3d *pm, int v1, int v2, int c1, int c2) {
    int i, j, jj, k, kk, kkk, t1, t2, t3, p1, p2, p3, v, w1, w2, w3, w0, r1, r2;
    REAL m1, m2, m3, m, rl1, rl2;
    sign_eps = sign_eps_default;
    libaft_3d_warn("ugly.c: force_edge() is deprecated! Please, use force_pre_edge()");
    while (1) {
	for (j=pm->lnT; j<pm->nT; j++) {
	    for (jj=0; jj<4; jj++) {
		if (pm->tetra[4*j+jj] == v1) {
		    t1 = pm->tetra[4*j + (jj+((jj%2)?3:1))%4];
		    t2 = pm->tetra[4*j + (jj+2)%4];
		    t3 = pm->tetra[4*j + (jj+((jj%2)?1:3))%4];
		    if ((t1==v2) || (t2==v2) || (t3==v2))  return 0;
		    m = det3i4(pm->vertex, v1, t1, t2, t3) - det3i4(pm->vertex, v2, t1, t2, t3);
		    if (sign_flat(m)==0)  continue;
		    rl1 = det3i4(pm->vertex, v1, t1, t2, t3)/m,  r1 = sign(rl1);
		    if ((r1 < 0) && (m > 0)) libaft_3d_warn("ugly.c: force_edge(): r1<=0 (v = %le)  (internal error)", rl1);
		    rl2 = det3i4(pm->vertex, v2, t1, t2, t3)/m,  r2 = sign(rl2);
		    if (r1*r2>0)  continue;
		    m1 = det3i4(pm->vertex, v1, v2, t2, t3)/m,  p1 = sign(m1);
		    m2 = det3i4(pm->vertex, v1, v2, t3, t1)/m,  p2 = sign(m2);
		    m3 = det3i4(pm->vertex, v1, v2, t1, t2)/m,  p3 = sign(m3);
		    if (r1<0) {
			m1 = -m1,  m2 = -m2,  m3 = -m3;
			p1 = -p1,  p2 = -p2,  p3 = -p3;
		    }
		    if (p1+p2+p3<-1) libaft_3d_warn("ugly.c: force_edge(): p1+p2+p3 < -1  (internal error)");
		    if ((p1<0) || (p2<0) || (p3<0))  continue;
		    if (p1+p2+p3==3) {
			m = m1+m2+m3;
			v = addvertex(pm,
				(m1*pm->vertex[3*t1+0] + m2*pm->vertex[3*t2+0] + m3*pm->vertex[3*t3+0])/m,
				(m1*pm->vertex[3*t1+1] + m2*pm->vertex[3*t2+1] + m3*pm->vertex[3*t3+1])/m,
				(m1*pm->vertex[3*t1+2] + m2*pm->vertex[3*t2+2] + m3*pm->vertex[3*t3+2])/m);
			pm->vtrace[3*v+0] = c1;
			pm->vtrace[3*v+1] = c2;
			pm->vtrace[3*v+2] = -1;
			remtetra(pm, j, 0);
			for (k=pm->lnT; k<pm->nT; k++) {
			    for (kk=0; kk<4; kk++) {
				w1 = pm->tetra[4*k + (kk+((kk%2)?3:1))%4];
				w2 = pm->tetra[4*k + (kk+2)%4];
				w3 = pm->tetra[4*k + (kk+((kk%2)?1:3))%4];
				if (((t1==w3) && (t2==w2) && (t3==w1)) ||
					((t2==w3) && (t3==w2) && (t1==w1)) ||
					((t3==w3) && (t1==w2) && (t2==w1))) {
				    w0 = pm->tetra[4*k + kk];
				    remtetra(pm, k, 0);
				    packtetra(pm);
				    addtetra(pm, v1, t1, t2, v, 0);
				    addtetra(pm, v1, t2, t3, v, 0);
				    addtetra(pm, v1, t3, t1, v, 0);
				    addtetra(pm, w0, t3, t2, v, 0);
				    addtetra(pm, w0, t2, t1, v, 0);
				    addtetra(pm, w0, t1, t3, v, 0);
				    return 1 + force_edge(pm, v, v2, c1, c2);
				}
			    }
			}
			libaft_3d_warn("ugly.c: force_edge(): case3 failed (internal error)");
		    } else if (p1+p2+p3==2) {
			if (p1==0)  {t1 = t3; m1 = m3;}
			if (p2==0)  {t2 = t3; m2 = m3;}
			m = m1+m2;
			v = addvertex(pm,
				(m1*pm->vertex[3*t1+0] + m2*pm->vertex[3*t2+0])/m,
				(m1*pm->vertex[3*t1+1] + m2*pm->vertex[3*t2+1])/m,
				(m1*pm->vertex[3*t1+2] + m2*pm->vertex[3*t2+2])/m);
			pm->vtrace[3*v+0] = c1;
			pm->vtrace[3*v+1] = c2;
			pm->vtrace[3*v+2] = -1;
			pm->nP = 0;
			for (k=pm->lnT; k<pm->nT; k++) {
			    for (kk=0; kk<4; kk++) {
				w0 = pm->tetra[4*k + kk];
				if (w0!=t1)  continue;
				for (kkk=1; kkk<4; kkk++) {
				    w1 = pm->tetra[4*k + (kk+kkk)%4];
				    if (w1!=t2)  continue;
				    w2 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?1:0))%3)%4];
				    w3 = pm->tetra[4*k + (kk+1+(kkk+((kk%2)?0:1))%3)%4];
				    pm->pair[2*pm->nP + 0] = w2;
				    pm->pair[2*pm->nP + 1] = w3;
				    pm->nP++;
				    remtetra(pm, k, 0);
				    break;
				}
			    }
			}
			packtetra(pm);
			if (pm->nP<2) libaft_3d_warn("ugly.c: force_edge(): case2: something is wrong...");
			/*printf("%d pairs\n", pm->nP);*/
			for (i=0; i<pm->nP; i++) {
			    w2 = pm->pair[2*i+0];
			    w3 = pm->pair[2*i+1];
			    /*printf("%3d %3d\n", w2, w3);*/
			    addtetra(pm, t1, v, w2, w3, 0);
			    addtetra(pm, v, t2, w2, w3, 0);
			}
			return 1 + force_edge(pm, v, v2, c1, c2);
		    }
		    if (0) save_ugly_int(pm->vertex, v1, t1, t2, t3, v2);
		    libaft_3d_warn("3D meshing: internal error");
		    libaft_3d_warn("p1=%d, p2=%d, p3=%d,  m1=%le, m2=%le, m3=%le.\n", p1, p2, p3, m1, m2, m3);
		    libaft_3d_warn("v1=%d, v2=%d,  t1=%d, t2=%d, t3=%d.\n", v1, v2, t1, t2, t3);
		    libaft_3d_warn("r1=%d (%le), r2=%d (%le).\n", r1, rl1, r2, rl2);
		    libaft_3d_warn("force_edge: case1, the impossible one :-O");
		}
	    }
	}
	if (sign_eps < sign_eps_big)  sign_eps *= 2.0;
	else  break;
    }
    libaft_3d_warn("3D meshing: internal error");
    libaft_3d_warn("force_edge failed!!!\nedge: <%d>, <%d>, (eps=%le)", v1+1, v2+1, sign_eps);
    sign_eps = sign_eps_default;
    return 0;
}
/**********************************************************************/
static int step1(mesh3d *pm) {
    int i, v1, v2;

    for (i=0; i<pm->nE; i++) {
	v1 = pm->edge[2*i+0];
	v2 = pm->edge[2*i+1];
	if (0)  force_edge(pm, v1, v2, v1, v2);
	if (1)  force_pre_edge(pm, v1, v2, v1, v2);
#ifdef SHOWPROGRESS
	printf("force edge: %6.2lf%%\r", 100.0*(i+1)/pm->nE),  fflush(stdout);
#endif
    }
    packtetra(pm);
    return 0;
}

static int step2(mesh3d *pm) {
    int i, v1, v2, r, k;
    if (pm->lnF != pm->nF) libaft_3d_warn("ugly.c: step2(): lnF != nF  (internal error)");

    r = 0;
    for (i=0; i<pm->nE; i++) {
	v1 = pm->edge[2*i+0];
	v2 = pm->edge[2*i+1];
	k = force_face(pm, v1, v2);
	if (k<0)  return -1;
	r += k;
#ifdef SHOWPROGRESS
	printf("force face: %6.2lf%%\r", 100.0*(i+1)/pm->nE),  fflush(stdout);
#endif
    }
    packtetra(pm);
    return r;
}

static int faceedges(mesh3d *pm) {
    int i, v1, v2, v3;
    pm->nE = 0;
    for (i=0; i<pm->nF; i++) {
	v1 = pm->face[3*i+0];
	v2 = pm->face[3*i+1];
	v3 = pm->face[3*i+2];
	caddedge(pm, v1, v2);
	caddedge(pm, v2, v3);
	caddedge(pm, v3, v1);
    }
    return pm->nE;
}
static int tetraedges(mesh3d *pm) {
    int i, v1, v2, v3, v4;
    pm->nE = 0;
    for (i=0; i<pm->nT; i++) {
	v1 = pm->tetra[4*i+0];
	v2 = pm->tetra[4*i+1];
	v3 = pm->tetra[4*i+2];
	v4 = pm->tetra[4*i+3];
	caddedge(pm, v1, v2);
	caddedge(pm, v2, v3);
	caddedge(pm, v3, v1);
	caddedge(pm, v1, v4);
	caddedge(pm, v2, v4);
	caddedge(pm, v3, v4);
    }
    return pm->nE;
}
static int selectfaces(mesh3d *pm) {
    int i, k, kk, v1, v2, v3, face;
    int buf[10], nbuf;
    REAL n1[3], n2[3], s;
    for (k=pm->lnT; k<pm->mnT; k++) {
	for (kk=0; kk<4; kk++) {
	    v1 = pm->tetra[4*k + (kk^3)];
	    v2 = pm->tetra[4*k + (kk^2)];
	    v3 = pm->tetra[4*k + (kk^1)];
	    nbuf = 0;
	    for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v1+i]);
	    for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v2+i]);
	    for (i=0; i<3; i++)  bufadd(&nbuf, buf, pm->vtrace[3*v3+i]);
	    if (nbuf<3)  {
		libaft_3d_warn("3D meshing: internal error");
		libaft_3d_warn("ugly.c: selectfaces(): ntrace < 3");
		libaft_3d_warn("additional info:\nface: %d %d %d\nvertex %d: %d %d %d\nvertex %d: %d %d %d\nvertex %d: %d %d %d\n",
				v1, v2, v3,
				v1, pm->vtrace[3*v1+0],  pm->vtrace[3*v1+1], pm->vtrace[3*v1+2],
				v2, pm->vtrace[3*v2+0],  pm->vtrace[3*v2+1], pm->vtrace[3*v2+2],
				v3, pm->vtrace[3*v3+0],  pm->vtrace[3*v3+1], pm->vtrace[3*v3+2] );
		continue;
	    }
	    if (nbuf>3)  continue;
//	    for (i=0; i<nbuf; i++)  printf(" %d", buf[i]);
	    face = findlface(pm, buf[0], buf[1], buf[2]);
//	    printf("  face: %d\n", face);
	    if (face < 0)  continue;
	    normvec3i(pm->vertex, v1, v2, v3, n1);
	    normvec3i(pm->vertex, pm->face[3*face+0], pm->face[3*face+1], pm->face[3*face+2], n2);
	    vnorm(n1),  vnorm(n2);
	    s = dot(n1, n2);
//	    printf("%11.8lf\n", s);
	    if (s>0)  addface(pm, v1, v2, v3, pm->facecolor[face]);
/*	    else  addface(pm, v3, v2, v1, pm->facecolor[face]);*/
	}
    }
    return 0;
}
static int findpoint(mesh3d *pm, int fn) {
    int v1, v2, v3, k, kk, pn, w1, w2, w3;

    v1 = pm->face[3*fn+0];
    v2 = pm->face[3*fn+1];
    v3 = pm->face[3*fn+2];
    for (k=pm->lnT; k<pm->mnT; k++) {
	for (kk=0; kk<4; kk++) {
	    w1 = pm->tetra[4*k + (kk^3)];
	    w2 = pm->tetra[4*k + (kk^2)];
	    w3 = pm->tetra[4*k + (kk^1)];
	    pn = pm->tetra[4*k + (kk^0)];
	    if ((v1==w1)&&(v2==w2)&&(v3==w3))  return pn;
	    if ((v2==w1)&&(v3==w2)&&(v1==w3))  return pn;
	    if ((v3==w1)&&(v1==w2)&&(v2==w3))  return pn;
	}
    }
    return -1;
}
static int front(mesh3d *pm) {
    int i, j, pn, fn, v[3], t;
    if (0)  savefrtraw(pm);
    while (pm->nF > pm->lnF) {
	fn = pm->nF-1;
	pn = findpoint(pm, fn);
	if (pn<0) {
	    libaft_3d_warn("ugly.c: front(): findpoint failed (internal error)");
	    return pn;
	}
	v[0] = pm->face[3*fn+0];
	v[1] = pm->face[3*fn+1];
	v[2] = pm->face[3*fn+2];
	addtetra(pm, v[0], v[1], v[2], pn, pm->facecolor[fn]);
	if (idet3i4(pm->vertex, v[0], v[1], v[2], pn) != 1)
	    libaft_3d_warn("ugly.c: front(): bad volume of new tetra [%d %d %d %d]: %le",
		    v[0], v[1], v[2], pn, det3i4(pm->vertex, v[0], v[1], v[2], pn));
	remface(pm, fn, pm->facecolor[fn]);
	for (i=0; i<3; i++) {
	    if ((j = findface(pm, v[(i+2)%3], v[(i+1)%3], pn))>=0)  remface(pm, j, pm->facecolor[fn]);
	    else  addface(pm, v[(i+1)%3], v[(i+2)%3], pn, pm->facecolor[fn]);
	}
	packface(pm);
#ifdef SHOWPROGRESS
	printf("ugly front: nF = % 4d, nT = % 4d  \r", pm->nF-pm->lnF, pm->nT-pm->mnT),  fflush(stdout);
#endif
	t = pm->nF - pm->lnF;
	if (0)  savefrtraw(pm);
	if (check2dsurface(&pm->nV, pm->vertex, &t, pm->face + 3*pm->lnF, 0, 0, 0)) {
	    libaft_3d_warn("ugly.c: front(): surface topology is broken (internal error)");
	}
    }
#ifdef SHOWPROGRESS
    printf("ugly front: done.                 \r"),  fflush(stdout);
#endif
    return 0;
}

static int deckchkadd2(mesh3d *pm, int v1, int v2, int v3) {
    int i;
    
    for (i=pm->nD; i<pm->lnD; i++) {
	if (    ((pm->deck[3*i+0] == v3) && (pm->deck[3*i+1] == v2) && (pm->deck[3*i+2] == v1)) ||
	    ((pm->deck[3*i+0] == v2) && (pm->deck[3*i+1] == v1) && (pm->deck[3*i+2] == v3)) ||
	    ((pm->deck[3*i+0] == v1) && (pm->deck[3*i+1] == v3) && (pm->deck[3*i+2] == v2))    ) {
	    pm->lnD--;
	pm->deck[3*i+0] = pm->deck[3*pm->lnD+0];
	pm->deck[3*i+1] = pm->deck[3*pm->lnD+1];
	pm->deck[3*i+2] = pm->deck[3*pm->lnD+2];
	return -1;
	}
    }
    if (pm->lnD >= pm->nnD) {
	libaft_3d_stop("ugly.c: maxD exceeded (internal problem)");
	return -1;
    }
    pm->deck[3*pm->lnD + 0] = v1;
    pm->deck[3*pm->lnD + 1] = v2;
    pm->deck[3*pm->lnD + 2] = v3;
    return pm->lnD++;
}
static int deckchkadd(mesh3d *pm, int v1, int v2, int v3) {
    int i;

    for (i=0; i<pm->nD; i++) {
	if (    ((pm->deck[3*i+0] == v3) && (pm->deck[3*i+1] == v2) && (pm->deck[3*i+2] == v1)) ||
		((pm->deck[3*i+0] == v2) && (pm->deck[3*i+1] == v1) && (pm->deck[3*i+2] == v3)) ||
		((pm->deck[3*i+0] == v1) && (pm->deck[3*i+1] == v3) && (pm->deck[3*i+2] == v2))    ) {
	    pm->nD--;
	    pm->deck[3*i+0] = pm->deck[3*pm->nD+0];
	    pm->deck[3*i+1] = pm->deck[3*pm->nD+1];
	    pm->deck[3*i+2] = pm->deck[3*pm->nD+2];
	    return -1;
	}
    }
    if (pm->nD >= pm->nnD) {
	libaft_3d_stop("ugly.c: maxD exceeded (internal problem)");
	return -1;
    }
    pm->deck[3*pm->nD + 0] = v1;
    pm->deck[3*pm->nD + 1] = v2;
    pm->deck[3*pm->nD + 2] = v3;
    return pm->nD++;
}
static void addmass(mesh3d *pm, int k, REAL *V) {
    int i, j, m;
    for (j=0; j<4; j++) {
	m = pm->tetra[4*k + j];
	for (i=0; i<3; i++)  V[i] += 0.25*pm->vertex[3*m+i];
    }
}
/*static int inv(mesh3d *pm) {
    int i, k;
    for (i=0; i<pm->ntpack; i++) {
	k = pm->tpack[i];
	if (det3i4(pm->vertex, pm->tetra[4*k+0], pm->tetra[4*k+1], pm->tetra[4*k+2], pm->tetra[4*k+3]) <= 0.0) return 1;
    }
    return 0;
}
static int trymove(mesh3d *pm, int i, REAL *V0, REAL *V1) {
    int t = 0;
    REAL M=1.0;
    while (inv(pm) && (t<500)) {
	//M /= 2.0;
	//M = (50.0-t)/50.0;
	M = exp(-(t/10.0));
	t++;
	pm->vertex[3*i+0] = M*V1[0] + (1.0-M)*V0[0];
	pm->vertex[3*i+1] = M*V1[1] + (1.0-M)*V0[1];
	pm->vertex[3*i+2] = M*V1[2] + (1.0-M)*V0[2];
    }
    return (t>=500);
}*/
static REAL minv(mesh3d *pm, int x) {
    REAL rmin = -1.0, r;
    int i, k;
    for (i=0; i<pm->ntpack; i++) {
	k = pm->tpack[i];
	if ((pm->tetra[4*k+0] != x) && (pm->tetra[4*k+1] != x) && (pm->tetra[4*k+2] != x) && (pm->tetra[4*k+3] != x))  continue;
	r = det3i4(pm->vertex, pm->tetra[4*k+0], pm->tetra[4*k+1], pm->tetra[4*k+2], pm->tetra[4*k+3]);
	if ((i) && (r >= rmin))  continue;
	rmin = r;
    }
    return rmin;
}
static int trymove(mesh3d *pm, int i, REAL *V0, REAL *V1) { // minimax
    int t = 0, tmax = -1;
    REAL M=1.0, r, rmax = -1.0;

    for (t=0; t<500; t++) {
	M = exp(-(t/10.0));
	pm->vertex[3*i+0] = M*V1[0] + (1.0-M)*V0[0];
	pm->vertex[3*i+1] = M*V1[1] + (1.0-M)*V0[1];
	pm->vertex[3*i+2] = M*V1[2] + (1.0-M)*V0[2];
	r = minv(pm, i);
	if (r > rmax) {
	    rmax = r;
	    tmax = t;
	}
    }
    if (tmax >= 0) {
	M = exp(-(tmax/10.0));
	pm->vertex[3*i+0] = M*V1[0] + (1.0-M)*V0[0];
	pm->vertex[3*i+1] = M*V1[1] + (1.0-M)*V0[1];
	pm->vertex[3*i+2] = M*V1[2] + (1.0-M)*V0[2];
    }
#ifdef FULLDEBUG
    if (rmax < 0.0)  libaft_3d_warn("ugly.c: trymove(): debug: %d %le", tmax, rmax);
#endif
    return ((tmax < 0) || (rmax < 0.0));
}
static int findbedge(mesh3d *pm, int c1, int c2, int *w) {
    int i, k, j, m, *v;
    k = 0;
    for (i=0; i<pm->nD; i++) {
	v = pm->deck + 3*i+1;
	for (j=0; j<2; j++) {
	    for (m=0; m<k; m++)  if (v[j] == w[m]) break;
	    if (m<k) continue;
	    if (    (v[j]==c1) || (v[j]==c2) ||
		    ((pm->vtrace[3*v[j]+0]==c1)&&(pm->vtrace[3*v[j]+1]==c2)&&(pm->vtrace[3*v[j]+2]<0)) ||
		    ((pm->vtrace[3*v[j]+0]==c2)&&(pm->vtrace[3*v[j]+1]==c1)&&(pm->vtrace[3*v[j]+2]<0))    )  {
		if (k>=2)  {libaft_3d_warn("ugly.c: findbedge(): k>2  (internal error)"); return -1;}
		w[k] = v[j];
		k++;
	    }
	}
    }
    return k;
}
static int splitdeck(mesh3d *pm, int pn, int *pw) {
    int k, wk, w, we, nused = 0;
    
    for (k=0; k<pm->nD; k++)  pm->dused[k] = 0;
    
    pm->nW = 0;
    
    while (nused < pm->nD) {
	for (wk=0; wk<2; wk++) {
	    if (pm->nW >= pm->nnW) {
		pm->nnW += 16;
		pm->wingstart = libaft_realloc(pm->wingstart, sizeof(int)*pm->nnW);
		pm->wingdir   = libaft_realloc(pm->wingstart, sizeof(char)*pm->nnW);
		pm->wingnorm  = libaft_realloc(pm->wingstart, sizeof(REAL)*pm->nnW*4);
	    }
	    pm->wingstart[pm->nW] = -1;
	    pm->wingdir[pm->nW] = wk;
	    w  = pw[wk];
	    we = pw[1-wk];
	    while (we != w) {
		for (k=0; k<pm->nD; k++) {
		    if (pm->dused[k])  continue;
		    if (pm->deck[3*k+0] != pn)  libaft_3d_warn("ugly.c: splitdeck(): deck corrupted  (internal error)");
		    if (pm->deck[3*k+2] != we)  continue;
		    pm->wingnext[k] = pm->wingstart[pm->nW];
		    pm->wingstart[pm->nW] = k;
		    we = pm->deck[3*k+1];
		    pm->dused[k] = 1;
		    nused++;
		    break;
		}
		if (k >= pm->nD) {
		    libaft_3d_warn("ugly.c: splitdeck(): search failed  (internal error)");
		    return 1;
		}
	    }
	    pm->nW++;
	}
    }
#ifdef FULLDEBUG
    if (pm->nW > 2)  libaft_3d_warn("Non-manifold case found! :)\n(nW = %d)", pm->nW);
#endif
    return 0;
}
static int tricheck(mesh3d *pm, int v1, int v2, int v3, int pn) {
    int v, i;
    for (i=0; i<pm->nbnd; i++) {
	v = pm->bnd[i];
	if ((v==v1) || (v==v2) || (v==v3))  continue;
	if (idet3i4(pm->vertex, v, v2, v3, pn) < 0) continue;
	if (idet3i4(pm->vertex, v1, v, v3, pn) < 0) continue;
	if (idet3i4(pm->vertex, v1, v2, v, pn) < 0) continue;
	return 1;
    }
    return 0;
}
static int reconstruct(mesh3d *pm, int pn, int color) {
    int i, j, v1, v2, v3, mk;
    REAL malpha, alpha;

    while (pm->nbnd > 2) {
	mk = -1,  malpha = 0.0;
	v1 = -1,  v2 = -1,  v3 = -1;
	for (i=0; i<pm->nbnd; i++) {
	    v1 = pm->bnd[(i+pm->nbnd-1) % pm->nbnd];
	    v2 = pm->bnd[(i) % pm->nbnd];
	    v3 = pm->bnd[(i+1) % pm->nbnd];
	    if (idet3i4(pm->vertex, v1, v2, v3, pn) < 1)  continue;
	    if (tricheck(pm, v1, v2, v3, pn))  continue;
	    alpha = vcos(pm->vertex, v1, v2, v3);
	    if ((mk<0) || (alpha > malpha))  malpha = alpha,  mk = i;
	}
	if (mk<0) {
	    for (i=0; i<pm->nbnd; i++) {
		v1 = pm->bnd[(i+pm->nbnd-1) % pm->nbnd];
		v2 = pm->bnd[(i) % pm->nbnd];
		v3 = pm->bnd[(i+1) % pm->nbnd];
		if (idet3i4(pm->vertex, v1, v2, v3, pn) < 0)  continue;
		if (tricheck(pm, v1, v2, v3, pn))  continue;
		alpha = vcos(pm->vertex, v1, v2, v3);
		if ((mk<0) || (alpha > malpha))  malpha = alpha,  mk = i;
	    }
	    if (mk<0) {
		libaft_3d_warn("ugly.c: reconstruct(): no candidates, trying first point  (warning)");
		mk = 0;
	    }
	}
	i = mk;
	v1 = pm->bnd[(i+pm->nbnd-1) % pm->nbnd];
	v2 = pm->bnd[(i) % pm->nbnd];
	v3 = pm->bnd[(i+1) % pm->nbnd];
	pm->tpack[pm->ntpack++] = addtetra(pm, v1, v2, v3, pn, color);
	pm->nbnd--;
	for (j=i; j<pm->nbnd; j++)  pm->bnd[j] = pm->bnd[j+1];
    }
    return 0;
}

static int facereconstruct(mesh3d *pm, int pn, int color) {
    int k, w, we;

    for (k=0; k<pm->nD; k++)  pm->dused[k] = 0;
    pm->nbnd = 0;
    w  = pm->deck[3*0+2];
    we = pm->deck[3*0+1];
    pm->bnd[pm->nbnd++] = we;
    while (w != we) {
	for (k=0; k<pm->nD; k++) {
	    if (pm->dused[k])  continue;
	    if (pm->deck[3*k+1] != w)  continue;
	    pm->bnd[pm->nbnd++] = w;
	    w = pm->deck[3*k+2];
	    pm->dused[k] = 1;
	    break;
	}
	if (k >= pm->nD) {
	    libaft_3d_warn("ugly.c: facereconstruct(): failed  (internal error)");
	    return 1;
	}
    }
    return reconstruct(pm, pn, color);
}
static int selectwing(mesh3d *pm, int orig, int curr, int shift) {
    int i, k, kk, kkk, w, v[3], ff;
    pm->lnD = pm->nD;
    for (w=shift; w<shift+2; w++) {
	k = pm->wingstart[w];
	if (k<0)  libaft_3d_warn("ugly.c: edgewingreconstruct(): void wing  (internal error)");
	while (k>=0) {
	    deckchkadd2(pm, pm->deck[3*k+0], pm->deck[3*k+1], pm->deck[3*k+2]);
	    k = pm->wingnext[k];
	}
    }
    while (pm->lnD > pm->nD) {
	ff = 0;
	for (i=0; i<pm->ntpack; i++) {
	    k = pm->tpack[i];
	    for (kk=0; kk<4; kk++) {
		if (pm->tetra[4*k+kk] == orig) {
		    v[0] = pm->tetra[4*k + (kk^1)];
		    v[1] = pm->tetra[4*k + (kk^2)];
		    v[2] = pm->tetra[4*k + (kk^3)];
		    for (kkk=0; kkk<3; kkk++) {
			if ((v[kkk]==pm->deck[3*pm->nD+1]) && (v[(kkk+1)%3]==pm->deck[3*pm->nD+2])) {
			    deckchkadd2(pm, orig, v[1], v[0]);
			    deckchkadd2(pm, orig, v[2], v[1]);
			    deckchkadd2(pm, orig, v[0], v[2]);
			    pm->tetra[4*k+kk] = curr;
			    ff++;
			    break;
			}
		    }
		}
	    }
	    if (ff)  break;
	}
	if (ff == 0) {
	    libaft_3d_warn("ugly.c: selectwing(): failed to select tetras (ff=%d)  (internal error)", ff);
	    break;
	}
    }
    
    return 0;
}
static int edgewingreconstruct(mesh3d *pm, int pn, int shift, int color) {
    int k, w, r[2];
    
    for (w=shift; w<shift+2; w++) {
	pm->nbnd = 0;
	k = pm->wingstart[w];
	if (k<0)  libaft_3d_warn("ugly.c: edgewingreconstruct(): void wing  (internal error)");
	pm->bnd[pm->nbnd++] = pm->deck[3*k+1];
	while (k>=0) {
	    pm->bnd[pm->nbnd++] = pm->deck[3*k+2];
	    k = pm->wingnext[k];
	}
	r[w-shift] = reconstruct(pm, pn, color);
    }
    return r[0] + r[1];
}
/*
static int edgereconstruct(mesh3d *pm, int pn, int *pw, int color) {
    int k, wk, w, we, r[2];

    for (k=0; k<pm->nD; k++)  pm->dused[k] = 0;
    for (wk=0; wk<2; wk++) {
	pm->nbnd = 0;
	w  = pw[wk];
	we = pw[1-wk];
	pm->bnd[pm->nbnd++] = we;
	while (w != we) {
	    for (k=0; k<pm->nD; k++) {
		if (pm->dused[k])  continue;
		if (pm->deck[3*k+1] != w)  continue;
		pm->bnd[pm->nbnd++] = w;
		w = pm->deck[3*k+2];
		pm->dused[k] = 1;
		break;
	    }
	    if (k >= pm->nD) {
		libaft_3d_warn("ugly.c: edgereconstruct(): failed  (internal error)");
		return 1;
	    }
	}
	r[wk] = reconstruct(pm, pn, color);
    }
    return r[0] + r[1];
}*/
static int normdeck2(mesh3d *pm, int pn, REAL *d) {
    int i, j, v1, v2, v3;
    REAL e[3];
    for (j=0; j<3; j++)  d[j] = 0.0;
    for (i=0; i<pm->nbnd-1; i++) {
	v1 = pn;
	v2 = pm->bnd[i];
	v3 = pm->bnd[i+1];
	normvec3i(pm->vertex, v2, v1, v3, e);
	for (j=0; j<3; j++)  d[j] += e[j];
    }
    return 0;
}
static int getwingnorms(mesh3d *pm, int pn) {
    int i, j, k;
    REAL norm[3], r;
    for (i=0; i<pm->nW; i++) {
	pm->nbnd = 0;
	k = pm->wingstart[i];
	if (k<0)  libaft_3d_warn("ugly.c: getwingnorms(): void wing  (internal error)");
	pm->bnd[pm->nbnd++] = pm->deck[3*k+1];
	while (k>=0) {
	    pm->bnd[pm->nbnd++] = pm->deck[3*k+2];
	    k = pm->wingnext[k];
	}
	normdeck2(pm, pn, norm);
	r = sqrt(dot(norm, norm));
	pm->wingnorm[4*i+3] = 32.0 * sqrt(r/M_PI);
	for (j=0; j<3; j++)  pm->wingnorm[4*i+j] = norm[j] / r * ((pm->wingdir[i])? 1 : -1);
    }
    return 0;
}
static REAL wingangle(REAL x[3], REAL y[3], REAL z[3]) {
    if (det3(x, y, z) <= 0.0)  return 1.0 - dot(x, y);
    else return 3.0 + dot(x, y);
}
static int sortwings(mesh3d *pm, int *pw) {
    int i, j, k, mink, ti;
    REAL z[3], mina, a, tr;

    for (j=0; j<3; j++)  z[j] = pm->vertex[3*pw[0]+j] - pm->vertex[3*pw[1]+j];
    for (i=1; i<pm->nW; i++) {
	mink = i;
	mina = wingangle(pm->wingnorm + 4*(i-1), pm->wingnorm + 4*mink, z);
	for (k=i+1; k<pm->nW; k++) {
	    a = wingangle(pm->wingnorm + 4*(i-1), pm->wingnorm + 4*k, z);
	    if (a < mina) {
		mink = k;
		mina = a;
	    }
	}
	if (mink > i) {
	    ti = pm->wingdir[i],    pm->wingdir[i]   = pm->wingdir[mink],    pm->wingdir[mink]   = ti;
	    ti = pm->wingstart[i],  pm->wingstart[i] = pm->wingstart[mink],  pm->wingstart[mink] = ti;
	    for (j=0; j<4; j++)  tr = pm->wingnorm[4*i+j],  pm->wingnorm[4*i+j] = pm->wingnorm[4*mink+j],  pm->wingnorm[4*mink+j] = tr;
	}
    }
    for (i=1; i<pm->nW; i++) {
	if (pm->wingdir[i-1] + pm->wingdir[i] != 1)  libaft_3d_warn("ugly.c: sortwings(): wrong order after sort (internal error)");
    }
#if FULLDEBUG
    if (pm->nW > 2) {
	printf("\nSorted:\n");
	for (i=0; i<pm->nW; i++)  printf((i)?", %lf":"%lf", wingangle(pm->wingnorm + 0, pm->wingnorm + 4*i, z));
	printf("\nSorted start:\n");
	for (i=0; i<pm->nW; i++)  printf((i)?", %d":"%d", pm->wingstart[i]);
	printf("\n");
    }
#endif
    return 0;
}
/*
static int normdeckedge(mesh3d *pm, int pn, int *pw, REAL *d) {
    int k, wk, w, we, j, nused = 0;
    REAL f[6];
    REAL r, r1, r2;
    
    for (k=0; k<pm->nD; k++)  pm->dused[k] = 0;

    for (wk=0; wk<2; wk++) {
	pm->nbnd = 0;
	w  = pw[wk];
	we = pw[1-wk];
	while (w != we) {
	    for (k=0; k<pm->nD; k++) {
		if (pm->dused[k])  continue;
		if (pm->deck[3*k+0] != pn)  libaft_3d_warn("ugly.c: normdeckedge(): deck corrupted  (internal error)");
		if (pm->deck[3*k+1] != w)  continue;
		pm->bnd[pm->nbnd++] = w;
		w = pm->deck[3*k+2];
		pm->dused[k] = 1;
		nused++;
		break;
	    }
	    if (k >= pm->nD) {
		libaft_3d_warn("ugly.c: normdeckedge(): failed  (internal error)");
		return 1;
	    }
	}
	pm->bnd[pm->nbnd++] = we;
	normdeck2(pm, pn, f+wk*3);
    }
    if (pm->nD != nused) {
	libaft_3d_warn("ugly.c: normdeckedge(): unused tria in deck, non-manifold region?  (internal error)");
    }
    r1 = sqrt(sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]) / M_PI),  r1 = 32.0 * r1;
    r2 = sqrt(sqrt(f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) / M_PI),  r2 = 32.0 * r2;
    r = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    for (j=0; j<3; j++)  f[j] /= r;
    r = sqrt(f[3]*f[3] + f[4]*f[4] + f[5]*f[5]);
    for (j=3; j<6; j++)  f[j] /= r;
    for (j=0; j<3; j++)  d[j] = f[j] + f[3+j];
    r = (r1+r2)/2.0;
    for (j=0; j<3; j++)  d[j] *= r;
    return 0;
}
*/
static REAL vq(mesh3d *pm, int i) {
    REAL rmin = 1e6, r;
    int k; 
    for (k=pm->mnT; k<pm->nT; k++) {
	if ((pm->tetra[4*k+0] == i) || (pm->tetra[4*k+1] == i) ||  (pm->tetra[4*k+2] == i) || (pm->tetra[4*k+3] == i)) {
	    r = det3i4(pm->vertex, pm->tetra[4*k+0], pm->tetra[4*k+1], pm->tetra[4*k+2], pm->tetra[4*k+3]);
	    if (r < rmin)  rmin = r;
	}
    }
    return rmin;
}
static int hide_prep_edge(mesh3d *pm) {
    int i;
    pm->nH = 0;
    for (i=0; i<pm->nV; i++) {
	if (pm->vtrace[3*i+1] < 0)  continue;
	if (pm->vtrace[3*i+2] >= 0)  continue;
	pm->hide[pm->nH++] = i;
    }
    return pm->nH;
}
static int hide_prep_face(mesh3d *pm) {
    int i;
    pm->nH = 0;
    for (i=0; i<pm->nV; i++) {
	if (pm->vtrace[3*i+1] < 0)  continue;
	if (pm->vtrace[3*i+2] < 0)  continue;
	pm->hide[pm->nH++] = i;
    }
    return pm->nH;
}
static int hide_edge_vertex(mesh3d *pm, int i) {
    int j, mi, c1, c2;
    int k, kk, v1, v2, v3;
    int w[2];
    int color = 0;
    int r = 0;
    REAL V0[3], V1[3], VN[3];

    c1 = pm->vtrace[3*i+0],  c2 = pm->vtrace[3*i+1];
    pm->ntpack = 0,  pm->nD = 0;
    for (j=0; j<3; j++)  V0[j] = pm->vertex[3*i+j];
    for (k=pm->mnT; k<pm->nT; k++) {
	for (kk=0; kk<4; kk++) {
	    if (pm->tetra[4*k+kk] == i) {
		v1 = pm->tetra[4*k + (kk^1)];
		v2 = pm->tetra[4*k + (kk^2)];
		v3 = pm->tetra[4*k + (kk^3)];
		color = pm->tetracolor[k];
		deckchkadd(pm, i, v1, v2);
		deckchkadd(pm, i, v2, v3);
		deckchkadd(pm, i, v3, v1);
		pm->tpack[pm->ntpack] = k;
		pm->ntpack++;
	    }
	}
    }
    if (pm->ntpack==0)  return 0;
#ifdef FULLDEBUG
    printf("edge: [%d %d]\n", c1, c2);
    printf("%d tetras, %d indeck\n", pm->ntpack, pm->nD);
#endif
    if (findbedge(pm, c1, c2, w)!=2) {
	libaft_3d_warn("ugly.c: hide_edge_vertex(): findbedge() failed!  (internal error)\nvertex: %d, edge: %d, %d,  w: %d, %d",
		i+1, c1+1, c2+1, w[0]+1, w[1]+1);
	return 1;
    }
    splitdeck(pm, i, w);
    getwingnorms(pm, i);
    sortwings(pm, w);
    if (pm->nW % 2)  libaft_3d_warn("ugly.c: hide_edge_vertex(): odd number of wings (internal error)");
    for (k=pm->nW/2-1; k>=0; k--) {
	if (k)  mi = addvertex(pm, V0[0], V0[1], V0[2]);
	else mi = i;
	for (j=0; j<3; j++) {
	    VN[j] = (pm->wingnorm[4*2*k+3] + pm->wingnorm[4*(2*k+1)+3])/2.0 * (pm->wingnorm[4*(2*k+1)+j] - pm->wingnorm[4*2*k+j]);
	    V1[j] = V0[j] - VN[j];
	    pm->vertex[3*mi+j] = V1[j];
	}
	if (k)  selectwing(pm, i, mi, 2*k);
	if (edgewingreconstruct(pm, mi, 2*k, color)) {
	    libaft_3d_warn("ugly.c: hide_edge_vertex(): edge reconstruct failed  (internal error)");
	    r++;
	}
	if (trymove(pm, mi, V0, V1))  pm->nHedge++;
	for (j=0; j<3; j++) pm->vtrace[3*mi+j] = -1;
    }
    return r;
}
static int hide_face_vertex(mesh3d *pm, int i) {
    int j, c1, c2, c3;
    int k, kk, v1, v2, v3;
    int color = 0;
    REAL V0[3], V1[3];

    c1 = pm->vtrace[3*i+0],  c2 = pm->vtrace[3*i+1],  c3 = pm->vtrace[3*i+2];
    pm->ntpack = 0,  pm->nD = 0;
    for (j=0; j<3; j++)  V0[j] = pm->vertex[3*i+j];
    for (j=0; j<3; j++)  V1[j] = 0.0;
    for (k=pm->mnT; k<pm->nT; k++) {
	for (kk=0; kk<4; kk++) {
	    if (pm->tetra[4*k+kk] == i) {
		v1 = pm->tetra[4*k + (kk^1)];
		v2 = pm->tetra[4*k + (kk^2)];
		v3 = pm->tetra[4*k + (kk^3)];
		color = pm->tetracolor[k];
		deckchkadd(pm, i, v1, v2);
		deckchkadd(pm, i, v2, v3);
		deckchkadd(pm, i, v3, v1);
		pm->tpack[pm->ntpack] = k;
		addmass(pm, k, V1),  pm->ntpack++;
	    }
	}
    }
    if (pm->ntpack==0)  return 0;
    if (pm->nD==0) {
	libaft_3d_warn("ugly.c: hide_face_vertex(): point is already inside region (?) (warning)");
	return 0;
    }
    for (j=0; j<3; j++)  V1[j] /= pm->ntpack;
    for (j=0; j<3; j++)  V1[j] = V0[j] + 32.0*(V1[j]-V0[j]);
#ifdef FULLDEBUG
    printf("face: [%d %d %d]\n", c1, c2, c3);
    printf("%d tetras, %d indeck\n", pm->ntpack, pm->nD);
#endif
    for (j=0; j<3; j++)  pm->vertex[3*i+j] = V1[j];
    if (facereconstruct(pm, i, color)) {
	libaft_3d_warn("ugly.c: facereconstruct() failed!  (internal error)");
	return 1;
    }
    if (trymove(pm, i, V0, V1))  pm->nHface++;
    for (j=0; j<3; j++) pm->vtrace[3*i+j] = -1;
    return 0;
}
static int hide_edge(mesh3d *pm) {
    int k = 0, i, j, r = 0;
    REAL rmax = 0.0, rc;
    hide_prep_edge(pm);
    while (pm->nH > 0) {
	for (j=0; j<pm->nH; j++) {
	    rc = vq(pm, pm->hide[j]);
	    if ((j) && (rc <= rmax))  continue;
	    rmax = rc;
	    k = j;
	}
	i = pm->hide[k];
	if (hide_edge_vertex(pm, i))  r++;
	pm->nH--;
	pm->hide[k] = pm->hide[pm->nH];
#ifdef SHOWPROGRESS
	printf("edge hide:  nH = % 4d \r", pm->nH),  fflush(stdout);
#endif
    }
    return r;
}
static int hide_face(mesh3d *pm) {
    int k = 0, i, j, r = 0;
    REAL rmax = 0.0, rc;
    hide_prep_face(pm);
    while (pm->nH > 0) {
	for (j=0; j<pm->nH; j++) {
	    rc = vq(pm, pm->hide[j]);
	    if ((j) && (rc <= rmax))  continue;
	    rmax = rc;
	    k = j;
	}
	i = pm->hide[k];
	if (hide_face_vertex(pm, i))  r++;
	pm->nH--;
	pm->hide[k] = pm->hide[pm->nH];
#ifdef SHOWPROGRESS
	printf("face hide:  nH = % 4d \r", pm->nH),  fflush(stdout);
#endif
    }
    return r;
}
static int hide(mesh3d *pm) {
    int r;
    pm->nHedge = 0;
    pm->nHface = 0;
    r = hide_edge(pm) + hide_face(pm);
#ifdef SHOWPROGRESS
    printf("                                                \r"),  fflush(stdout);
#endif
    if (pm->nHedge + pm->nHface > 0) {
	libaft_3d_warn("3D meshing: Failed to correctly move %d points (%d from edges, %d from faces)", pm->nHedge + pm->nHface, pm->nHedge, pm->nHface);
    }
    return r;
}

int mesh3dugly(int *pnV, REAL *vertex,
	int *pnF, int *face, int *facecolor,
	int *pnT, int *tetra, int *tetracolor,
	int nnV, int nnF, int nnT) {
    mesh3d mesh;
    int i, j, r, st, t, nev;
    (void) st;

    detinit();

    mesh.nV  = *pnV,  mesh.nF  = *pnF,  mesh.nT  = *pnT;
    mesh.nnV =  nnV,  mesh.nnF =  nnF,  mesh.nnT =  nnT;

    mesh.vertex = vertex;
    mesh.face   = face,    mesh.facecolor  = facecolor;
    mesh.tetra  = tetra,   mesh.tetracolor = tetracolor;

    mesh.nE     = 0;
    mesh.nnE    = 4 * mesh.nnT;
    mesh.edge   = libaft_malloc(sizeof(int) * mesh.nnE*2);

    mesh.nP     = 0;
    mesh.nnP    = 3 * mesh.nF;
    mesh.pair   = libaft_malloc(sizeof(int) * mesh.nnP*2);

    mesh.nD     = 0;
    mesh.nnD    = 3 * mesh.nnT;
    mesh.deck   = libaft_malloc(sizeof(int) * mesh.nnD*3);
    mesh.dused  = libaft_malloc(sizeof(int) * mesh.nnD);

    mesh.nW        = 0;
    mesh.nnW       = 16;
    mesh.wingstart = libaft_malloc(sizeof(int)  * mesh.nnW);
    mesh.wingnorm  = libaft_malloc(sizeof(REAL) * mesh.nnW*4);
    mesh.wingdir   = libaft_malloc(sizeof(char) * mesh.nnW);
    mesh.wingnext  = libaft_malloc(sizeof(int)  * mesh.nnD);
    
    mesh.nbnd   = 0;
    mesh.bnd    = libaft_malloc(sizeof(int) * mesh.nnV);

    mesh.nH     = 0;
    mesh.hide   = libaft_malloc(sizeof(int) * mesh.nnV);

    mesh.vtrace = libaft_malloc(sizeof(int) * mesh.nnV*3);

    mesh.npack  = 0;
    mesh.pack   = libaft_malloc(sizeof(int) * mesh.nnF);
    mesh.ntpack = 0;
    mesh.tpack  = libaft_malloc(sizeof(int) * mesh.nnT);
    mesh.ncolor = 0;
    mesh.color  = libaft_malloc(sizeof(int) * (mesh.nF+1)*2);
    
    mesh.ndupF   = 0;
    mesh.dupface = libaft_malloc(sizeof(int) * mesh.nnF*3);

    mesh.lnV = mesh.nV;
    mesh.lnF = mesh.nF;
    mesh.lnT = mesh.nT;

    if (mesh.edge && mesh.pair && mesh.deck && mesh.dused && mesh.wingnext && mesh.wingstart && mesh.wingnorm && mesh.wingdir && mesh.bnd && mesh.hide && mesh.vtrace && mesh.pack && mesh.tpack && mesh.color) {

	if (mesh.nT>0)  libaft_3d_warn("ugly.c: mesh3dugly(): nT>0  (bad input)");

	for (i=0; i<mesh.lnV; i++) {
	    mesh.vtrace[3*i+0] = i;
	    mesh.vtrace[3*i+1] = -1;
	    mesh.vtrace[3*i+2] = -1;
	}

	if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/

	delaunay3d_main(&mesh.nV, mesh.vertex, &mesh.ndupF, mesh.dupface, &mesh.nT, mesh.tetra, nnV, nnF, nnT);	
	for (i=mesh.lnT; i<mesh.nT; i++) mesh.tetracolor[i] = 0;

	if (mesh.nV < mesh.lnV)  libaft_3d_warn("ugly.c: mesh3dugly(): Delaunay removed some points! (internal error)");
	nev = mesh.nV - mesh.lnV;

	for (i=mesh.lnV; i<mesh.nV; i++) {
	    mesh.vtrace[3*i+0] = i;
	    mesh.vtrace[3*i+1] = -1;
	    mesh.vtrace[3*i+2] = -1;
	}
	if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/

	if (checktopology(mesh.nV, mesh.vertex, mesh.ndupF, mesh.dupface, mesh.nT, mesh.tetra, 0))
	    libaft_3d_warn("ugly.c: mesh3dugly(): Wrong topology after Delaunay (internal error)");

#ifdef SHOWPROGRESS
	printf("force edge: init...                          \r"),  fflush(stdout);
#endif
	faceedges(&mesh);
	step1(&mesh);
#ifdef STRUCTCHECK	
	check_tetra_faces(&mesh);
#endif
	if (checktopology(mesh.nV, mesh.vertex, mesh.ndupF, mesh.dupface, mesh.nT, mesh.tetra, 0))
	    libaft_3d_warn("ugly.c: mesh3dugly(): Wrong topology after Step 1 (internal error)");

	if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/

	i = 0;
	while (1) {
#ifdef SHOWPROGRESS
	    printf("force face: init...                          \r"),  fflush(stdout);
#endif
	    tetraedges(&mesh);
	    i++;
	    r = step2(&mesh);
#ifdef STRUCTCHECK	
	    check_tetra_faces(&mesh);
#endif
	    if (r == 0)  break;
	    if ((i>5) || (r<0)) {
		if (DUMP)  gmvwrite(&mesh, 0);
		libaft_free(mesh.vtrace);
		libaft_free(mesh.pack);
		libaft_free(mesh.tpack);
		libaft_free(mesh.color);
		libaft_free(mesh.bnd);
		libaft_free(mesh.hide);
		libaft_free(mesh.deck);
		libaft_free(mesh.dused);
		libaft_free(mesh.wingnext);
		libaft_free(mesh.wingstart);
		libaft_free(mesh.wingnorm);
		libaft_free(mesh.wingdir);
		libaft_free(mesh.pair);
		libaft_free(mesh.edge);
		return 1;
	    }
	}

	if (checktopology(mesh.nV, mesh.vertex, mesh.ndupF, mesh.dupface, mesh.nT, mesh.tetra, 0))
	    libaft_3d_warn("ugly.c: mesh3dugly(): Wrong topology after Step 2 (internal error)");

	if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/

#ifdef SHOWPROGRESS
	printf("ugly front: init...                          \r"),  fflush(stdout);
#endif
	mesh.mnT = mesh.nT;
	if (mesh.nF != mesh.lnF)  libaft_3d_warn("ugly.c: mesh3dugly(): Unexpected new faces! (internal error)");
	selectfaces(&mesh);

	if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	t = mesh.nF - mesh.lnF;
	if (check2dsurface(&mesh.nV, mesh.vertex, &t, mesh.face + 3*mesh.lnF, 0, 0, 0)) {
//	if (check_surface_topology(mesh.nV, mesh.vertex, mesh.nF - mesh.lnF, mesh.face + 3*mesh.lnF)) {
//   gmvwrite(&mesh, 0); /******************************************************/
//   gmvwrite(&mesh, 1); /******************************************************/
	    libaft_free(mesh.vtrace);
	    libaft_free(mesh.pack);
	    libaft_free(mesh.tpack);
	    libaft_free(mesh.color);
	    libaft_free(mesh.bnd);
	    libaft_free(mesh.hide);
	    libaft_free(mesh.deck);
	    libaft_free(mesh.dused);
	    libaft_free(mesh.wingnext);
	    libaft_free(mesh.wingstart);
	    libaft_free(mesh.wingnorm);
	    libaft_free(mesh.wingdir);
	    libaft_free(mesh.pair);
	    libaft_free(mesh.edge);
	    return 1;
	}

	mesh.ndupF = mesh.nF - mesh.lnF;
	for (i=0; i<mesh.ndupF; i++) {
	    mesh.dupface[3*i + 0] = mesh.face[3*(i+mesh.lnF) + 0];
	    mesh.dupface[3*i + 1] = mesh.face[3*(i+mesh.lnF) + 1];
	    mesh.dupface[3*i + 2] = mesh.face[3*(i+mesh.lnF) + 2];
	}

	r = front(&mesh);
	if (r) {
	    if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/
	    if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	}
	for (i=mesh.mnT; i<mesh.nT; i++)  mesh.tetracolor[i] = pcolor(&mesh, mesh.tetracolor[i]);
	if (mesh.ntpack > 0)  libaft_3d_warn("ugly.c: mesh3dugly(): Tetras not packed! (internal error)");

	if (checktopology(mesh.nV, mesh.vertex, mesh.ndupF, mesh.dupface, mesh.nT - mesh.mnT, mesh.tetra + 4*mesh.mnT, 0))
	    libaft_3d_warn("ugly.c: mesh3dugly(): Wrong topology after front (internal error)");

	if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/

	if (!r) {
	    if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	    r = hide(&mesh);
#ifdef SHOWPROGRESS
	    printf("3D meshing: done.                               \r"),  fflush(stdout);
#endif
	    if (r) {
		if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/
		if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	    }
	    if (checktopology(mesh.nV, mesh.vertex, mesh.lnF, mesh.face, mesh.nT - mesh.mnT, mesh.tetra + 4*mesh.mnT, 0)) {
		libaft_3d_warn("ugly.c: mesh3dugly(): Wrong topology after hide (internal error)");
		if (DUMP)  gmvwrite(&mesh, 0); /******************************************************/
		if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	    }

	    if (DUMP)  gmvwrite(&mesh, 1); /******************************************************/
	    for (i=mesh.mnT, j=mesh.lnT; i<mesh.nT; i++, j++) {
		mesh.tetra[4*j+0]  = mesh.tetra[4*i+0];
		mesh.tetra[4*j+1]  = mesh.tetra[4*i+1];
		mesh.tetra[4*j+2]  = mesh.tetra[4*i+2];
		mesh.tetra[4*j+3]  = mesh.tetra[4*i+3];
		mesh.tetracolor[j] = mesh.tetracolor[i];
	    }
	    mesh.nT = j;
	    mesh.nF = 0;
	    if (nev > 0) {
		for (j=mesh.lnT; j<mesh.nT; j++) {
		    for (i=0; i<4; i++) {
			if ((mesh.tetra[4*j+i] >= mesh.lnV) && (mesh.tetra[4*j+i] < mesh.lnV+nev))
			    libaft_3d_warn("ugly.c: mesh3dugly(): Extra point still in use (internal error)");
		    }
		}
		for (i=mesh.lnV; i<mesh.nV-nev; i++) {
		    mesh.vertex[3*i+0] = mesh.vertex[3*(i+nev)+0];
		    mesh.vertex[3*i+1] = mesh.vertex[3*(i+nev)+1];
		    mesh.vertex[3*i+2] = mesh.vertex[3*(i+nev)+2];
		}
		mesh.nV -= nev;
		for (j=mesh.lnT; j<mesh.nT; j++) {
		    for (i=0; i<4; i++) {
			if (mesh.tetra[4*j+i] >= mesh.lnV)  mesh.tetra[4*j+i] -= nev;
		    }
		}
	    }
	}
    } else  r = -1;

    libaft_free(mesh.dupface);
    libaft_free(mesh.vtrace);
    libaft_free(mesh.pack);
    libaft_free(mesh.tpack);
    libaft_free(mesh.color);
    libaft_free(mesh.bnd);
    libaft_free(mesh.hide);
    libaft_free(mesh.deck);
    libaft_free(mesh.dused);
    libaft_free(mesh.wingnext);
    libaft_free(mesh.wingstart);
    libaft_free(mesh.wingnorm);
    libaft_free(mesh.wingdir);
    libaft_free(mesh.pair);
    libaft_free(mesh.edge);
    *pnV = mesh.nV, *pnF = mesh.nF, *pnT = mesh.nT;
    return r;
}





static int gmvwrite(mesh3d *pm, int flag) {
    int i, j;
    static int num=0;
    FILE *f;
    char fn[1024];
    sprintf(fn, "dump%04d.gmv", num);
    f=fopen(fn, "w");
    fprintf(f, "gmvinput ascii\n\n nodes %5d\n", pm->nV);
    for (j=0; j<3; j++) {
	for (i=0; i<pm->nV; i++) {
	    fprintf(f, " %20.15lf", pm->vertex[3*i+j]);
	}
	fprintf(f, "\n");
    }
    fprintf(f, "\n cells %5d\n", pm->nT-((flag)?pm->mnT:0));
    for (i=(flag)?pm->mnT:0; i<pm->nT; i++) {
	fprintf(f, "  tet 4\n  %4d %4d %4d %4d\n", pm->tetra[4*i+0]+1, pm->tetra[4*i+1]+1, pm->tetra[4*i+2]+1, pm->tetra[4*i+3]+1);
    }
    fprintf(f, "\npolygons\n");
    for (i=(flag)?pm->lnF:0; i<pm->nF; i++) {
	if (pm->face[3*i] < 0) continue;
	fprintf(f, "%3d 3", pm->facecolor[i]);
	fprintf(f,      " %20.15lf %20.15lf %20.15lf\n", pm->vertex[3*pm->face[3*i+0]+0], pm->vertex[3*pm->face[3*i+1]+0], pm->vertex[3*pm->face[3*i+2]+0]);
	fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", pm->vertex[3*pm->face[3*i+0]+1], pm->vertex[3*pm->face[3*i+1]+1], pm->vertex[3*pm->face[3*i+2]+1]);
	fprintf(f, "      %20.15lf %20.15lf %20.15lf\n", pm->vertex[3*pm->face[3*i+0]+2], pm->vertex[3*pm->face[3*i+1]+2], pm->vertex[3*pm->face[3*i+2]+2]);
    }
    fprintf(f, "endpoly\n");
    fprintf(f, "\nendgmv");
    fclose(f);
    printf("%04d: nV=%d, nF=%d, nT=%d\n", num, pm->nV, pm->nF-((flag)?pm->lnF:0), pm->nT-((flag)?pm->mnT:0));
    return num++;
}


static int save_ugly_int(REAL *vertex, int a, int b, int c, int d, int x) {
    FILE *f;
    static int fn=0;
    char buf[1024];

    sprintf(buf, "ugly_int%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "5 4 1 5\n");
    fprintf(f, "%lf %lf %lf\n", vertex[3*a+0], vertex[3*a+1], vertex[3*a+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*b+0], vertex[3*b+1], vertex[3*b+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*c+0], vertex[3*c+1], vertex[3*c+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*d+0], vertex[3*d+1], vertex[3*d+2]);
    fprintf(f, "%lf %lf %lf\n", vertex[3*x+0], vertex[3*x+1], vertex[3*x+2]);
    fprintf(f, "1 2 3\n");
    fprintf(f, "2 3 4\n");
    fprintf(f, "3 4 1\n");
    fprintf(f, "4 1 2\n");
    fprintf(f, "1 5\n");
    fprintf(f, "1 A\n");
    fprintf(f, "2 B\n");
    fprintf(f, "3 C\n");
    fprintf(f, "4 D\n");
    fprintf(f, "5 X\n");
    fclose(f);
    return fn++;
}
static int savefrtraw(mesh3d *pm) {
    FILE *f;
    int i;
    static int fn=0;
    char buf[1024];

    sprintf(buf, "frt%d.raw", fn);
    f=fopen(buf, "w");
    if (!f) return 1;
    fprintf(f, "%d %d\n", pm->nV, pm->nF-pm->lnF);
    for (i=0; i<pm->nV; i++) {
	fprintf(f, "%20.16le %20.16le %20.16le\n", pm->vertex[3*i+0], pm->vertex[3*i+1], pm->vertex[3*i+2]);
    }
    for (i=pm->lnF; i<pm->nF; i++) {
	fprintf(f, "%d %d %d\n", pm->face[3*i+0]+1, pm->face[3*i+1]+1, pm->face[3*i+2]+1);
    }
    fclose(f);
    return fn++;
}
