#include <time.h>
#include "aft.h"
#include "common.h"
#include "tetra3.h"
#include "tria32.h"
#include "region3.h"
#include "error3.h"
#include "tree3.h"
#include "section3.h"
#include "user3.h"

// PROF -- enable timing
#define PROF


int tria_dump_front = 0;
int tria_debug_front = 0;
int tria_final_front = 0;

static double minNear = 0.5;

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

#define N_SPLIT 16
#define N_DEEP 20
#define N_ANGLE 16

/*static int temp_hist[32];*/

static int dump();
static int dump3();

static int goodPlane32(surface_mesh *pm, int v1, int v2, int v);
static void surfNormal(surface_mesh *pm, int iSurf, int iNorm, double u, double v, StrucVert3 *norm);
static int choseWorkVert(surface_mesh *pm);

static int shiftperiodiclocal(surface_mesh *pm, int k, int i);

static double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/

static double ddet3(double a, double b, double c,
		    double d, double e, double f,
		    double g, double h, double i) {
    return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
}
static int idet4(surface_mesh *pm, int a, int b, int c, int d) {
    double len = (distance(pm->vert[a].x, pm->vert[a].y, pm->vert[a].z,  pm->vert[d].x, pm->vert[d].y, pm->vert[d].z) +
		  distance(pm->vert[b].x, pm->vert[b].y, pm->vert[b].z,  pm->vert[d].x, pm->vert[d].y, pm->vert[d].z) +
		  distance(pm->vert[c].x, pm->vert[c].y, pm->vert[c].z,  pm->vert[d].x, pm->vert[d].y, pm->vert[d].z)) / 3.0;
    double ddet = ddet3(pm->vert[a].x-pm->vert[d].x, pm->vert[b].x-pm->vert[d].x, pm->vert[c].x-pm->vert[d].x,
	    		pm->vert[a].y-pm->vert[d].y, pm->vert[b].y-pm->vert[d].y, pm->vert[c].y-pm->vert[d].y,
			pm->vert[a].z-pm->vert[d].z, pm->vert[b].z-pm->vert[d].z, pm->vert[c].z-pm->vert[d].z);
    ddet /= len*len*len;
    if (ddet > 1e-5)  return +1;
    else if (ddet < -1e-5)  return -1;
    else return 0;
}

static int faceIntersectNew(surface_mesh *pm, int a, int b, int c, int u, int v, int d) {
    int uv, dup=0;

    if ((u==a)||(u==b)||(u==c)) dup++;
    if ((v==a)||(v==b)||(v==c)) dup++;
    if (dup==0) {
	uv = idet4(pm, u, v, a, d) + idet4(pm, u, v, b, d) + idet4(pm, u, v, c, d);
	if ((uv==3) || (uv==-3)) return 0;
	if (idet4(pm, b, c, u, d) + idet4(pm, b, c, v, d) == -2) return 0;
	if (idet4(pm, c, a, u, d) + idet4(pm, c, a, v, d) == -2) return 0;
	if (idet4(pm, a, b, u, d) + idet4(pm, a, b, v, d) == -2) return 0;
    } else if (dup==1) {
	if (idet4(pm, b, c, u, d) + idet4(pm, b, c, v, d) == -1) return 0;
	if (idet4(pm, c, a, u, d) + idet4(pm, c, a, v, d) == -1) return 0;
	if (idet4(pm, a, b, u, d) + idet4(pm, a, b, v, d) == -1) return 0;
    } else return 0;
    return 1;
}

static int checkIntersect(surface_mesh *pm, PStrucEdge3  face, int v) {
    int i, iFace, v1, v2, workVert, bln;
    PStrucEdge3 f;

    v1 = face->v1;  v2 = face->v2;

    for (i=0; i<3; i++) {
	pm->sneigBool[i]= 0;
    }
    iFace = 0;
    while (iFace < pm->tnVicinityFace) {
	f = pm->tvicinityFace[iFace].face;
	if ( (v1==f->v2) && (v==f->v1) ) {
	    pm->sneigBool[0] = 1;
	    pm->sneigFace[0] = f;
	}

	if ( (v2==f->v1) && (v==f->v2) ) {
	    pm->sneigBool[1] = 1;
	    pm->sneigFace[1] = f;
	}
	iFace++;
    }

    iFace = bln = 0;
    while (iFace < pm->tnVicinityFace) {
	if (1) {
	    f = pm->tvicinityFace[iFace].face;
	    {/* proection on perpendicular to surfNormal */
//		StrucVert3  norm,e1,e2;
		double  x,y,z,p1,p2,p,len;

		shiftperiodiclocal(pm, v1, v);

		x = (pm->vert[f->v1].x + pm->vert[f->v2].x)/2.0;
		y = (pm->vert[f->v1].y + pm->vert[f->v2].y)/2.0;
		z = (pm->vert[f->v1].z + pm->vert[f->v2].z)/2.0;
		p1 = 1.0/(1e-6 + distance(x,y,z, pm->vert[v1].x,pm->vert[v1].y,pm->vert[v1].z));
//		surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v1].u, pm->vert[v1].v, &e1);
		p2 = 1.0/(1e-6 + distance(x,y,z, pm->vert[v2].x,pm->vert[v2].y,pm->vert[v2].z));
//		surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v2].u, pm->vert[v2].v, &e2);
		p =  1.0/(1e-6 + distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z));
//		surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v].u, pm->vert[v].v, &norm);
		len = p + p1 + p2;
//		pm->sv[0].x = (p1*e1.x+p2*e2.x+p*norm.x)/len;
//		pm->sv[0].y = (p1*e1.y+p2*e2.y+p*norm.y)/len;
//		pm->sv[0].z = (p1*e1.z+p2*e2.z+p*norm.z)/len;
		surfNormal(pm, pm->iSurf, pm->iNorm,
			(p1*pm->vert[v1].u+p2*pm->vert[v2].u+p*pm->vert[v].u)/len,
			(p1*pm->vert[v1].v+p2*pm->vert[v2].v+p*pm->vert[v].v)/len, &pm->sv[0]);
		x += pm->sv[0].x*3.0/len,  y += pm->sv[0].y*3.0/len,  z += pm->sv[0].z*3.0/len;
		pm->nPoint++;
		addPoint(pm, x,y,z);
		pm->nPoint--;
		pm->nPoint--;
	    }

	    if (1 && faceIntersectNew(pm, v1, v2, v, f->v1, f->v2, pm->nPoint+1)) {
		bln++;
		pm->ssect = pm->tvicinityFace[iFace].face;
		workVert = choseWorkVert(pm);
		if (workVert >= 0) {
		    return  bln;
		}
	    }
	}
	iFace++;
    }/* while  iFace */

    return  bln;
}/*checkIntersect*/


static int choseWorkVert(surface_mesh *pm) {
    int i, best=-1, v[2], b[2]={1,1};
    int v1, v2;
    double p, len=-1.0;

    v[0] = pm->ssect->v1;
    v[1] = pm->ssect->v2;

    v1 = pm->swork->v1;   v2 = pm->swork->v2;
    for (i=0; i<2; i++) {
	if ( (v[i] == v1) || (v[i] == v2) || isBadVert(pm, v[i]) || !goodPlane32(pm, v1, v2, v[i]) )  b[i] = 0;
    }
    for (i=0; i<2; i++) {
	if (b[i]) {
	    p = distanceS(pm->vert[v1].x,pm->vert[v1].y,pm->vert[v1].z,  pm->vert[v[i]].x,pm->vert[v[i]].y,pm->vert[v[i]].z) +
		distanceS(pm->vert[v2].x,pm->vert[v2].y,pm->vert[v2].z,  pm->vert[v[i]].x,pm->vert[v[i]].y,pm->vert[v[i]].z);
	    if ((best<0) || (p < len)) {
		len = p;
		best = v[i];
	    }
	}
    }
    return best;
} /*choseWorkVert*/


static void checkNearEdge(surface_mesh *pm) {
    int i, iFace, v1, v2, v, bestFace=-1, best;
    PStrucEdge3 f;
    StrucVert3  e1,e2;
    double c, p, dist=-1.0, x, y, z, len;
    int vs[2], b[2]={1,1};

    v = pm->nPoint;

    for (iFace = 0; iFace < pm->tnVicinityFace; iFace++) {
	f = pm->tvicinityFace[iFace].face;
	v1 = f->v1,  v2 = f->v2;
	makeVector(pm, &e1, v,  v1);
	makeVector(pm, &e2, v2, v1);
	c = e1.x*e2.x + e1.y*e2.y + e1.z*e2.z;
	c /= e2.x*e2.x + e2.y*e2.y + e2.z*e2.z;
	if ((c>=0.0) && (c<=1.0)) {
	    x = (1.0-c)*pm->vert[v1].x + c*pm->vert[v2].x;
	    y = (1.0-c)*pm->vert[v1].y + c*pm->vert[v2].y;
	    z = (1.0-c)*pm->vert[v1].z + c*pm->vert[v2].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	} else if (c<0.5) {
	    x = pm->vert[v1].x;
	    y = pm->vert[v1].y;
	    z = pm->vert[v1].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	} else {
	    x = pm->vert[v2].x;
	    y = pm->vert[v2].y;
	    z = pm->vert[v2].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	}
	if ((bestFace < 0) || (p < dist)) {
	    dist = p;
	    bestFace = iFace;
	}
    }
    if (bestFace >= 0) {
	f = pm->tvicinityFace[bestFace].face;
	vs[0] = f->v1;
	vs[1] = f->v2;
	v1 = pm->swork->v1;   v2 = pm->swork->v2;
	for (i=0; i<2; i++) {
	    if ( (vs[i] == v1) || (vs[i] == v2) || isBadVert(pm, vs[i]) || !goodPlane32(pm, v1, v2, vs[i]) )  b[i] = 0;
	}
	best = -1,  len = -1.0;
	for (i=0; i<2; i++) {
	    if (b[i]) {
		p = distance(pm->vert[v].x,pm->vert[v].y,pm->vert[v].z,  pm->vert[vs[i]].x,pm->vert[vs[i]].y,pm->vert[vs[i]].z);
		if ((best<0) || (p < len)) {
		    len = p;
		    best = vs[i];
		}
	    }
	}
	if (best>=0)  pm->sworkVert = best;
    }
} /*checkNearEdge*/


static void checkNearEdgeVik(surface_mesh *pm, double threshold) {
    int i, iFace, v1, v2, v, bestFace=-1, best;
    PStrucEdge3 f;
    StrucVert3  e1,e2;
    double c, p, dist=-1.0, x, y, z, len;
    int vs[2], b[2]={1,1};

    v = pm->nPoint;

    for (iFace = 0; iFace < pm->tnVicinityFace; iFace++) {
	f = pm->tvicinityFace[iFace].face;
	v1 = f->v1,  v2 = f->v2;
	makeVector(pm, &e1, v,  v1);
	makeVector(pm, &e2, v2, v1);
	c = e1.x*e2.x + e1.y*e2.y + e1.z*e2.z;
	c /= e2.x*e2.x + e2.y*e2.y + e2.z*e2.z;
	if ((c>=0.0) && (c<=1.0)) {
	    x = (1.0-c)*pm->vert[v1].x + c*pm->vert[v2].x;
	    y = (1.0-c)*pm->vert[v1].y + c*pm->vert[v2].y;
	    z = (1.0-c)*pm->vert[v1].z + c*pm->vert[v2].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	} else if (c<0.5) {
	    x = pm->vert[v1].x;
	    y = pm->vert[v1].y;
	    z = pm->vert[v1].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	} else {
	    x = pm->vert[v2].x;
	    y = pm->vert[v2].y;
	    z = pm->vert[v2].z;
	    p = distance(x,y,z, pm->vert[v].x,pm->vert[v].y,pm->vert[v].z);
	}
	if ((bestFace < 0) || (p < dist)) {
	    dist = p;
	    bestFace = iFace;
	}
    }
    if ((bestFace >= 0) && (dist < threshold)) {
	f = pm->tvicinityFace[bestFace].face;
	vs[0] = f->v1;
	vs[1] = f->v2;
	v1 = pm->swork->v1;   v2 = pm->swork->v2;
	for (i=0; i<2; i++) {
	    if ( (vs[i] == v1) || (vs[i] == v2) || isBadVert(pm, vs[i]) || !goodPlane32(pm, v1, v2, vs[i]) )  b[i] = 0;
	}
	best = -1,  len = -1.0;
	for (i=0; i<2; i++) {
	    if (b[i]) {
		p = distance(pm->vert[v].x,pm->vert[v].y,pm->vert[v].z,  pm->vert[vs[i]].x,pm->vert[vs[i]].y,pm->vert[vs[i]].z);
		if ((best<0) || (p < len)) {
		    len = p;
		    best = vs[i];
		}
	    }
	}
	if (best>=0)  pm->sworkVert = best;
    }
} /*checkNearEdge*/


double determ(double a11, double a12, double a21, double a22) {
    return a11*a22 - a12*a21;
} /*determ*/


void surfNormal(surface_mesh *pm, int iSurf, int iNorm, double u, double v, StrucVert3 *norm) {
    double x0, y0, z0,  x, y, z;
    double dxdu, dydu, dzdu, dxdv, dydv, dzdv;
    double du, dv, er;
    double nxyz[3];

    bounSurf(pm, iSurf, u, v, &x0, &y0, &z0);
    er = 0.1*sizeFace(pm, x0, y0, z0);

    if (0 && (iSurf == -1) && (pm->cgmsurfNormal)) {
	pm->cgmsurfNormal(x0, y0, z0, nxyz);
	norm->x = nxyz[0],  norm->y = nxyz[1],  norm->z = nxyz[2];
	du = -sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z);
	norm->x /= du;
	norm->y /= du;
	norm->z /= du;
/*	if (iNorm > 0) {
	    norm->x = -norm->x;
	    norm->y = -norm->y;
	    norm->z = -norm->z;
	}*/
	return;
    }

    du = (u < (pm->uMax + pm->uMin)/2.0) ? pm->uMax - pm->uMin : pm->uMin - pm->uMax;
    do {
	du /= 2.0;
	bounSurf(pm, iSurf, u+du, v, &x, &y, &z);
    } while (distance(x0,y0,z0, x,y,z) > er);
    dxdu = (x - x0)/du;
    dydu = (y - y0)/du;
    dzdu = (z - z0)/du;

    dv = (v < (pm->vMax + pm->vMin)/2.0) ? pm->vMax - pm->vMin : pm->vMin - pm->vMax;
    do {
	dv /= 2.0;
	bounSurf(pm, iSurf, u, v+dv, &x, &y, &z);
    } while (distance(x0,y0,z0, x,y,z) > er);
    dxdv = (x - x0)/dv;
    dydv = (y - y0)/dv;
    dzdv = (z - z0)/dv;

    norm->x = -determ(dydu, dzdu, dydv, dzdv);
    norm->y = -determ(dzdu, dxdu, dzdv, dxdv);
    norm->z = -determ(dxdu, dydu, dxdv, dydv);

    du = sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z);
    if (du < 1e-12) {
	if (dxdu*dxdu + dydu*dydu + dzdu*dzdu > dxdv*dxdv + dydv*dydv + dzdv*dzdv) {
	    printf("				du %c %c ", (du>0)?'+':'-', (dv>0)?'+':'-');
	    du = (u < (pm->uMax + pm->uMin)/2.0) ? pm->uMax - pm->uMin : pm->uMin - pm->uMax;
	    dv = (v < (pm->vMax + pm->vMin)/2.0) ? pm->vMax - pm->vMin : pm->vMin - pm->vMax;
	    do {
		du /= 2.0;
		dv /= 2.0;
		bounSurf(pm, iSurf, u+du, v+dv, &x, &y, &z);
	    } while (distance(x0,y0,z0, x,y,z) > er);
	    if (du*dv > 0.0)  iNorm = (iNorm > 0) ? 0 : 1;
	    dxdv = (x - x0)/sqrt(du*du + dv*dv);
	    dydv = (y - y0)/sqrt(du*du + dv*dv);
	    dzdv = (z - z0)/sqrt(du*du + dv*dv);
	    if (v-dv >= pm->vMin && v-dv <= pm->vMax)  bounSurf(pm, iSurf, u+du, v-dv, &x, &y, &z),  dv = sqrt(du*du + dv*dv),  printf("1\n");
	    else /*if (u-du >= pm->uMin && u-du <= pm->uMax) */ bounSurf(pm, iSurf, u-du, v+dv, &x, &y, &z),  dv = -sqrt(du*du + dv*dv),  printf("2\n");
	    //	else bounSurf(pm, iSurf, u+du, v, &x, &y, &z),  dv = du;
	    dxdu = (x - x0)/dv;
	    dydu = (y - y0)/dv;
	    dzdu = (z - z0)/dv;
	} else {
	    du = (u < (pm->uMax + pm->uMin)/2.0) ? pm->uMax - pm->uMin : pm->uMin - pm->uMax;
	    dv = (v < (pm->vMax + pm->vMin)/2.0) ? pm->vMax - pm->vMin : pm->vMin - pm->vMax;
	    do {
		du /= 2.0;
		dv /= 2.0;
		bounSurf(pm, iSurf, u+du, v+dv, &x, &y, &z);
	    } while (distance(x0,y0,z0, x,y,z) > er);
	    if (du*dv > 0.0)  iNorm = (iNorm > 0) ? 0 : 1;
	    printf("				dv %c %c ", (du>0)?'+':'-', (dv>0)?'+':'-');
	    dxdu = (x - x0)/sqrt(du*du + dv*dv);
	    dydu = (y - y0)/sqrt(du*du + dv*dv);
	    dzdu = (z - z0)/sqrt(du*du + dv*dv);
	    if (u-du >= pm->uMin && u-du <= pm->uMax)  bounSurf(pm, iSurf, u-du, v+dv, &x, &y, &z),  du = -sqrt(du*du + dv*dv),  printf("1\n");
	    else /*if (u-du >= pm->uMin && u-du <= pm->uMax) */ bounSurf(pm, iSurf, u+du, v-dv, &x, &y, &z),  du = sqrt(du*du + dv*dv),  printf("2\n");;
	    //	else bounSurf(pm, iSurf, u+du, v, &x, &y, &z),  dv = du;
	    dxdv = (x - x0)/du;
	    dydv = (y - y0)/du;
	    dzdv = (z - z0)/du;
	}
	norm->x = -determ(dydu, dzdu, dydv, dzdv);
	norm->y = -determ(dzdu, dxdu, dzdv, dxdv);
	norm->z = -determ(dxdu, dydu, dxdv, dydv);

	du = sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z);
	if (du < 1e-12) errorExit3(3," du == 0.   in  surfNormal ");
    }
    norm->x /= du;
    norm->y /= du;
    norm->z /= du;

    if (iNorm > 0) {
	norm->x = -norm->x;
	norm->y = -norm->y;
	norm->z = -norm->z;
    }

    return;
} /*surfNormal*/


void rePlane(surface_mesh *pm, int vv) {
    double x, y, z, u, v;

    x = pm->vert[vv].x - pm->x0;
    y = pm->vert[vv].y - pm->y0;
    z = pm->vert[vv].z - pm->z0;

    u = x*pm->x1 + y*pm->y1 + z*pm->z1;
    v = x*pm->x2 + y*pm->y2 + z*pm->z2;

    pm->vert[vv].u = u;
    pm->vert[vv].v = v;

    return;
} /*rePlane*/


static void goodPlane32prep(surface_mesh *pm, int v1, int v2) {
    return surfNormal(pm, pm->iSurf, pm->iNorm, (pm->vert[v1].u+pm->vert[v2].u)/2.0, (pm->vert[v1].v+pm->vert[v2].v)/2.0, &pm->gp);
}
static int goodPlane32(surface_mesh *pm, int v1, int v2, int v) {
    double p;
    StrucVert3 norm, e1, e2;

    makeVector(pm, &e1, v, v2);
    makeVector(pm, &e2, v2, v1);
    makeNormal(&norm, &e1, &e2);
/*    surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v1].u, pm->vert[v1].v, &e1);
    surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v2].u, pm->vert[v2].v, &e2);
    e1.x += e2.x;
    e1.y += e2.y;
    e1.z += e2.z;*/
//    surfNormal(pm, pm->iSurf, pm->iNorm, (pm->vert[v1].u+pm->vert[v2].u)/2.0, (pm->vert[v1].v+pm->vert[v2].v)/2.0, &e1);
//    p = norm.x*e1.x + norm.y*e1.y + norm.z*e1.z;
    p = norm.x*pm->gp.x + norm.y*pm->gp.y + norm.z*pm->gp.z;
    if (p < 1.e-7) return 0;
    return 1;
}

static int shiftperiodic(surface_mesh *pm, int k) {
    double pu, pv;
    double x, y, z;
    int i, q, j;
//    if (pm->iSurf >= 0)  return 0;
    pu = periodic(pm, 0),  pv = periodic(pm, 1);
    if (pu > 0.0) {
	pm->uMax = pm->vert[k].u + pu/2.0;
	pm->uMin = pm->vert[k].u - pu/2.0;
	for (j=0; j<pm->nEdge; j++) {
	    i = pm->edge[j]->v1;
	    if (i==k)  continue;
	    q = 0;
	    while (pm->vert[i].u > pm->uMax)  pm->vert[i].u -= pu,  q--;
	    while (pm->vert[i].u < pm->uMin)  pm->vert[i].u += pu,  q++;
	    if (q) {
		bounSurf(pm, pm->iSurf, pm->vert[i].u, pm->vert[i].v, &x, &y, &z);
		if (distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z) > 1e-3)
		    printf("U period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			    pu, q, distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z),
			    pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z, pm->vert[i].u, pm->vert[i].v);
//			pm->vert[i].u -= pu*q;
	    }
	}
    }
    if (pv > 0.0) {
	pm->vMax = pm->vert[k].v + pv/2.0;
	pm->vMin = pm->vert[k].v - pv/2.0;
	for (j=0; j<pm->nEdge; j++) {
	    i = pm->edge[j]->v1;
	    if (i==k)  continue;
	    q = 0;
	    while (pm->vert[i].v > pm->vMax)  pm->vert[i].v -= pv,  q--;
	    while (pm->vert[i].v < pm->vMin)  pm->vert[i].v += pv,  q++;
	    if (q) {
		bounSurf(pm, pm->iSurf, pm->vert[i].u, pm->vert[i].v, &x, &y, &z);
		if (distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z) > 1e-3)
		    printf("V period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			    pv, q, distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z),
			    pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z, pm->vert[i].u, pm->vert[i].v);
//			pm->vert[i].v -= pv*q;
	    }
	}
    }
    return 0;
}

static int shiftperiodiclocal(surface_mesh *pm, int k, int i) {
    double pu, pv;
    double x, y, z;
    int q;
    pu = periodic(pm, 0),  pv = periodic(pm, 1);
    if (pu > 0.0) {
	pm->uMax = pm->vert[k].u + pu/2.0;
	pm->uMin = pm->vert[k].u - pu/2.0;
	q = 0;
	while (pm->vert[i].u > pm->uMax)  pm->vert[i].u -= pu,  q--;
	while (pm->vert[i].u < pm->uMin)  pm->vert[i].u += pu,  q++;
	if (q) {
	    bounSurf(pm, pm->iSurf, pm->vert[i].u, pm->vert[i].v, &x, &y, &z);
	    if (distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z) > 1e-3)
		printf("U period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			pu, q, distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z),
			pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z, pm->vert[i].u, pm->vert[i].v);
	    //		pm->vert[i].u -= pu*q;
	}
    }
    if (pv > 0.0) {
	pm->vMax = pm->vert[k].v + pv/2.0;
	pm->vMin = pm->vert[k].v - pv/2.0;
	q = 0;
	while (pm->vert[i].v > pm->vMax)  pm->vert[i].v -= pv,  q--;
	while (pm->vert[i].v < pm->vMin)  pm->vert[i].v += pv,  q++;
	if (q) {
	    bounSurf(pm, pm->iSurf, pm->vert[i].u, pm->vert[i].v, &x, &y, &z);
	    if (distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z) > 1e-3)
		printf("V period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			pv, q, distance(pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z),
			pm->vert[i].x, pm->vert[i].y, pm->vert[i].z,  x, y, z, pm->vert[i].u, pm->vert[i].v);
	    //		pm->vert[i].v -= pv*q;
	}
    }
    return 0;
}

#define NQ 256
static int newPoint(surface_mesh *pm, PStrucEdge3 face) {
    int i, j, v1, v2, nearVert;
    int nTest=0, sec=0;
    int bound=0;
    double weakness = 1.0, realsize;
    double x, y, z,  x1, y1, z1,  x2, y2, z2, size, dist, d1, d2;
    StrucVert3 norm, e3, e4;

    double fi1, p;
    double u, v, du, dv;

    double fi, fiMin, fiMax;
    double ra, raMin, raMax;
    double sin0, cos0, si, co;

    StrucVert3  surfNorm;

    v1 = face->v1;  v2 = face->v2;
    shiftperiodiclocal(pm, v1, v2);
    goodPlane32prep(pm, v1, v2);
    if (0)  shiftperiodic(pm, v1);
    x1 = pm->vert[v1].x;  y1 = pm->vert[v1].y;  z1 = pm->vert[v1].z;
    x2 = pm->vert[v2].x;  y2 = pm->vert[v2].y;  z2 = pm->vert[v2].z;

    u = pm->vert[v1].u;
    v = pm->vert[v1].v;
    bounSurf(pm, pm->iSurf, u, v, &x, &y, &z);
    surfNormal(pm, pm->iSurf, pm->iNorm, (pm->vert[v1].u + pm->vert[v2].u)/2.0, (pm->vert[v1].v + pm->vert[v2].v)/2.0, &surfNorm);
    pm->uSave = pm->vert[v1].u;
    pm->vSave = pm->vert[v1].v;
    if ((pm->uSave > pm->uMax) || (pm->uSave < pm->uMin) || (pm->vSave > pm->vMax) || (pm->vSave < pm->vMin)) {
	fprintf(stderr, "\nparams (u,v) are out of param-box\n");
    }
    dist = distance(x1, y1, z1,  x2, y2, z2);
    if (pm->cf <= 0.0)  size = sizeFace(pm, 0.5*(x1+x2), 0.5*(y1+y2), 0.5*(z1+z2));
    else  size = pm->cf * dist;
    if ((pm->sizelim > 0.0) && (size > pm->sizelim))  size = pm->sizelim;
    //	size = dist*1.3;
    if (size < 0.5*dist)  {
	size = 0.6*dist;
	printf("\n...increasing local mesh size...\n");
    }
    realsize = size;

    du = pm->vert[v2].u - pm->vert[v1].u;
    dv = pm->vert[v2].v - pm->vert[v1].v;
    p = sqrt(du*du + dv*dv);
    if (p == 0.)  {
	printf("\nv1=%d, v2=%d, dist=%lf\n", v1, v2, dist);
	printf("%lf, %lf  :: %lf, %lf\n", pm->vert[v1].u, pm->vert[v1].v, pm->vert[v2].u, pm->vert[v2].v);
	printf("%lf, %lf, %lf\n", x1, y1, z1);
	printf("%lf, %lf, %lf\n", x2, y2, z2);
	printf("%lf, %lf, %lf\n", x, y, z);
	if (0)  dump(pm);
	printf("p == 0. in newPoint \n");
	return -21;
//	errorExit3(2, "p == 0. in newPoint ");
    }
    cos0 = du/p;
    sin0 = dv/p;
    sec = -1;
    while (1) {
	sec++;
	if (sec>2*N_SPLIT) {
	    weakness *= 2,  sec = -1,  bound++;
	    if (weakness > 20.0) {
		bound = -1;
		co = cos0;
		si = sin0;
		ra = dist*0.5;
		u = pm->uSave;
		v = pm->vSave;
		u += ra*co;
		v += ra*si;
		x = (x1+x2)/2.0,  y = (y1+y2)/2.0,  z = (z1+z2)/2.0;
		pm->vert[pm->nPoint].u = u;
		pm->vert[pm->nPoint].v = v;
		addPoint(pm, x,y,z);
		pm->nPoint--;
		break;
		size = (dist/2.0 + 3.0*size)/4.0;
		weakness = 1;
		if (fabs(dist - 2.0*size)/dist < 1e-3) {
		    return -20;
		    co = cos0;
		    si = sin0;
		    ra = dist*0.5;
		    break;
		}
	    }
	    continue;
	}
	u = pm->uSave;
	v = pm->vSave;
	if (sec % 2 == 0) {
	    fiMin = M_PI/N_SPLIT * (sec/2);
	    fiMax = M_PI/N_SPLIT * (sec/2 + 1);
	} else {
	    fiMin = 2*M_PI - M_PI/N_SPLIT * (sec/2);
	    fiMax = 2*M_PI - M_PI/N_SPLIT * (sec/2 + 1);
	}
	for (i=0; i<N_ANGLE; i++) {
	    fi = 0.5*(fiMin + fiMax);
	    co = cos0*cos(fi) - sin0*sin(fi);
	    si = sin0*cos(fi) + cos0*sin(fi);
	    raMin = 0;
	    raMax = 2.0*sqrt((pm->uMax-pm->uMin)*(pm->uMax-pm->uMin) + (pm->vMax-pm->vMin)*(pm->vMax-pm->vMin));
	    for (j=0; j<N_DEEP; j++) {/*calculate ra for fi*/
		ra = 0.5*(raMin + raMax);
		if ((u+ra*co > pm->uMax) || (u+ra*co < pm->uMin) || (v+ra*si > pm->vMax) || (v+ra*si < pm->vMin)) {
		    raMax = ra;
		    continue;
		}
		bounSurf(pm, pm->iSurf,u+ra*co,v+ra*si,&x,&y,&z);
		d1 = distance(x1,y1,z1,x,y,z);
		if (fabs(d1-size)/size < 0.0001)
		    break;
		else if (d1 > size)
		    raMax = ra;
		else
		    raMin = ra;
	    }
	    if ((u+ra*co > pm->uMax) || (u+ra*co < pm->uMin) || (v+ra*si > pm->vMax) || (v+ra*si < pm->vMin)) {
		fiMax = fi;
	    }
	    d1 = distance(x1,y1,z1,x,y,z);
	    d2 = distance(x2,y2,z2,x,y,z);
	    if (fabs(d2-d1)/d1 < 0.0001)
		break;
	    else if (d1 < d2)
		fiMax = fi;
	    else
		fiMin = fi;
	}
	if ((u+ra*co > pm->uMax) || (u+ra*co < pm->uMin) || (v+ra*si > pm->vMax) || (v+ra*si < pm->vMin))  continue;
//	if ((!weak)&&(( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3 ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3 )))  continue;
	if (( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3 * weakness ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3 * weakness ))  continue;
	u += ra*co;
	v += ra*si;
	if (bounSurf(pm, pm->iSurf,u,v,&x,&y,&z) == -1)  continue;
	pm->vert[pm->nPoint].u = u;
	pm->vert[pm->nPoint].v = v;
	addPoint(pm, x,y,z);
	pm->nPoint--;

	/*   printf("v= %5d   %5.5lf  %5.5lf  x=%5.5lf  y=%5.5lf  z=%5.5lf \n",pm->nPoint,u,v,x0,y0,z0); */
	makeVector(pm, &e3, pm->nPoint, v2);
	makeVector(pm, &e4, v2, v1);
	makeNormal(&norm, &e3, &e4);
//	surfNormal(pm, pm->iSurf, pm->iNorm, pm->vert[v1].u, pm->vert[v1].v, &e3);
	surfNormal(pm, pm->iSurf, pm->iNorm, (pm->vert[v1].u+pm->vert[v2].u)/2.0, (pm->vert[v1].v+pm->vert[v2].v)/2.0, &e3);
	fi1 = norm.x*e3.x + norm.y*e3.y + norm.z*e3.z;
	if (fabs(fi1) < 1e-3) {
	    printf("small cos!\nfi1 = %le\n", fi1);
	}
	if (fi1 < 1.e-7)  continue;
	/*if (!bound)  temp_hist[sec]++;*/
	break;
    }

    if (bounSurf(pm, pm->iSurf,pm->vert[pm->nPoint].u,pm->vert[pm->nPoint].v,&x,&y,&z) == -1)  printf("\nbounsurf failed!\n");
//    if ((fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3))
    if ((bound>=0) && (( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-1 ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-1 ))) {
	printf("	size = %lf, s1 = %lf (%lf), s2 = %lf (%lf)\n", size, distance(x1,y1,z1,x,y,z), (size-distance(x1,y1,z1,x,y,z))/size,
		distance(x2,y2,z2,x,y,z), (size-distance(x2,y2,z2,x,y,z))/size);
	if (0)  dump(pm);
	return -10;
    }
    pm->sworkVert = pm->nPoint;
    pm->swork = face;

    /******************  TEST  ****************/
//    x = (pm->vert[v1].x+pm->vert[v2].x+pm->vert[pm->nPoint].x)/3.;
//    y = (pm->vert[v1].y+pm->vert[v2].y+pm->vert[pm->nPoint].y)/3.;
//    z = (pm->vert[v1].z+pm->vert[v2].z+pm->vert[pm->nPoint].z)/3.;
    x = pm->vert[pm->nPoint].x;
    y = pm->vert[pm->nPoint].y;
    z = pm->vert[pm->nPoint].z;
    size = max(distance(x,y,z, pm->vert[v1].x,pm->vert[v1].y,pm->vert[v1].z),
	    distance(x,y,z, pm->vert[v2].x,pm->vert[v2].y,pm->vert[v2].z));
    if (bound < 0) {
	if (realsize > sqrt(2.0)*size)  size = sqrt(realsize*realsize - size*size);
    }
    if (size > realsize*1.2) {
	fprintf(stderr, "\nvery big size\n");
	if (bound < 0)  size = realsize*1.2;
    }
//    size = max(size, distance(x,y,z,pm->vert[pm->nPoint].x,pm->vert[pm->nPoint].y,pm->vert[pm->nPoint].z));
    vicinityFaces32(pm, x, y, z, 2.0 * size * 1025.0/1024.0);

    if (1) {  /* pm->sphereVert */
	int          i,iFace,v[2];
	double       min,p,rad;
	PStrucEdge3  face;

/*	x = 0.5*(pm->vert[v1].x+pm->vert[v2].x);
	y = 0.5*(pm->vert[v1].y+pm->vert[v2].y);
	z = 0.5*(pm->vert[v1].z+pm->vert[v2].z);
	radS = distanceS(x,y,z,pm->vert[pm->sworkVert].x,pm->vert[pm->sworkVert].y,pm->vert[pm->sworkVert].z);  
*/
//	rad = size*sqrt(2.0);
	rad = size;
	if (bound<0)  rad = rad*2.0;
	pm->nSphereVert = 0;
	for(iFace=0;iFace<pm->tnVicinityFace;iFace++){
	    face = pm->tvicinityFace[iFace].face;
	    v[0] = face->v1;
	    v[1] = face->v2;
	    for (i=0; i<2; i++) {
		if ((v[i]==v1) || (v[i]==v2))  continue;
		if (isSphereVert(pm, v[i]))  continue;
		if (!goodPlane32(pm, v1,v2,v[i]))  continue;
		p = distance(x,y,z, pm->vert[v[i]].x,pm->vert[v[i]].y,pm->vert[v[i]].z);  
		if (p < rad)  addSphereVert(pm, v[i]);
	    }
	}
	min = -1.0;
	v[0] = -1;
	for (i=0; i<pm->nSphereVert; i++) {
	    p = distance(x,y,z, pm->vert[pm->sphereVert[i]].x,pm->vert[pm->sphereVert[i]].y,pm->vert[pm->sphereVert[i]].z);
	    if ((v[0]<0) || (p < min)) {
		min = p;
		v[0] = pm->sphereVert[i];
	    }  
	}
	if ( (v[0] != -1) && (min < 0.9*size) )  pm->sworkVert = v[0];
	if ((bound < 0) && (pm->nSphereVert == 1))  pm->sworkVert = pm->sphereVert[0]; /* Last triangle workaround */
    }  /* pm->sphereVert */     

    x = pm->vert[pm->nPoint].x;
    y = pm->vert[pm->nPoint].y;
    z = pm->vert[pm->nPoint].z;
    if (pm->sworkVert == pm->nPoint) {
	dist = nearest32(pm, &nearVert,x,y,z);
	if ( (dist < minNear*size) && (nearVert != v1) && (nearVert != v2) ) {
	    if (goodPlane32(pm, v1,v2,nearVert))  pm->sworkVert = nearVert;
	}
    }
    if ( (pm->sworkVert == pm->nPoint) && (bound) )  checkNearEdge(pm);
    if (pm->sworkVert == pm->nPoint)  checkNearEdgeVik(pm, 0.2*size);
    if ( (pm->sworkVert == pm->nPoint) && (bound<0) )  return  -20;
    pm->nBadVert = 0;
    nTest = 0;

    x = (pm->vert[v1].x+pm->vert[v2].x+pm->vert[pm->nPoint].x)/3.;
    y = (pm->vert[v1].y+pm->vert[v2].y+pm->vert[pm->nPoint].y)/3.;
    z = (pm->vert[v1].z+pm->vert[v2].z+pm->vert[pm->nPoint].z)/3.;

    while (1) {
	nTest++;
	if( nTest > NQ )
		errorExit3(2," nTest > NQ ");

	if (size < distance(x,y,z,pm->vert[pm->sworkVert].x,pm->vert[pm->sworkVert].y,pm->vert[pm->sworkVert].z)) {
	    size = distance(x,y,z,pm->vert[pm->sworkVert].x,pm->vert[pm->sworkVert].y,pm->vert[pm->sworkVert].z);
	    vicinityFaces32(pm, x, y, z, size*1025.0/1024.0);
	}
	if (checkIntersect(pm, face, pm->sworkVert) == 0)  break;

	addBadVert(pm, pm->sworkVert);
	pm->sworkVert = choseWorkVert(pm);

	if (pm->sworkVert < 0) {
	    double  p, len=-1.0;
	    int best=-1;

	    for (i=0; i<pm->nSphereVert; i++) {
		if (isBadVert(pm, pm->sphereVert[i]))  continue;
		p = distanceS(pm->vert[v1].x,pm->vert[v1].y,pm->vert[v1].z,
			pm->vert[pm->sphereVert[i]].x,pm->vert[pm->sphereVert[i]].y,pm->vert[pm->sphereVert[i]].z) +
		    distanceS(pm->vert[v2].x,pm->vert[v2].y,pm->vert[v2].z,
			    pm->vert[pm->sphereVert[i]].x,pm->vert[pm->sphereVert[i]].y,pm->vert[pm->sphereVert[i]].z);
		if ((best<0) || (p<len)) {
		    len = p;
		    best = pm->sphereVert[i];
		}  
	    }
	    if (best >= 0)  pm->sworkVert = best;
	    else {
		printf("\nposible intersection in TRIA32\n");
		printf("Please, try with smaller mesh size\n");
		if (0)  dump3(pm);
		return -1;
	    }
	}
    }
    return 1;
} /*newPoint*/


/* newTria return codes:
 *  0 – OK
 *  2 – max_nT exceeded
 * -2 – box (check front orientation, check bounding box)
 */
static int newTria(surface_mesh *pm) {
    int i, r;
    int v[3];
    PStrucEdge3 face;

    for (i=0; i<pm->nEdge; i++) {
	face = pm->edge[i];
	if (face->fail)  continue;
	if ((r = newPoint(pm, face)) > 0)  break;
	face->fail = 1;
    }
    if (i>=pm->nEdge)  {
	if (!pm->fresh) {
	    printf("\nFront restart\n");
	    for (i=0; i<pm->nEdge; i++)  pm->edge[i]->fail = 0;
	    for (i=0; i<pm->nEdge; i++) {
		face = pm->edge[i];
		if (face->fail)  continue;
		if ((r = newPoint(pm, face)) > 0)  break;
		face->fail = 1;
	    }
	    if (i>=pm->nEdge)  {
		if (0)  dump(pm);
		return r;
	    }
	    pm->fresh = 1;
	} else {
	    if (0) dump(pm);
	    return r;
	}
    }

    v[0] = face->v1;
    v[1] = pm->sworkVert;
    v[2] = face->v2;
    /*printf("%5d: %5d %5d %5d\n",pm->nTria,v[0],v[2],v[1]);*/
    if (addTria(pm, v[0], v[2], v[1]) < 0)  return 2;
    remFace32(pm, face);

    if (pm->sworkVert == pm->nPoint) {
	if (    (fabs(pm->vert[pm->nPoint].x - pm->boxcx) > pm->boxsize*1.0000001) ||
		(fabs(pm->vert[pm->nPoint].y - pm->boxcy) > pm->boxsize*1.0000001) ||
		(fabs(pm->vert[pm->nPoint].z - pm->boxcz) > pm->boxsize*1.0000001)    ) {
	    printf("\nbounding box: x: %le, y: %le, z: %le, size: %le\n", pm->boxcx, pm->boxcy, pm->boxcz, pm->boxsize);
	    printf("new point:    x: %le, y: %le, z: %le\n", pm->vert[pm->nPoint].x, pm->vert[pm->nPoint].y, pm->vert[pm->nPoint].z);
	    return -2;
	}
	pm->nPoint++;
    }
    for (i=0; i<2; i++) {
	if (pm->sneigBool[i]) {
	    remFace32(pm, pm->sneigFace[i]);
	} else {
	    addFace32(pm, v[i], v[i+1]);
	}
    }
    pm->fresh = 0;
    return  0;
} /*newTria*/


static int dump(surface_mesh *pm) {
    static int num=0;
    char fname[1024];
    FILE *f;
    int i;
    sprintf(fname, "_dump_frt2_smv.%03d", num);
    f = fopen(fname, "w");

    fprintf(f, "%d %d %d %d 1\n", pm->nPoint, pm->nTria, pm->nEdge, (pm->nEdge)?pm->nEdge+2:0);
    for (i=0; i<pm->nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", pm->vert[i].x, pm->vert[i].y, pm->vert[i].z);
    for (i=0; i<pm->nTria; i++) fprintf(f, "%d %d %d\n", pm->v1[i]+1, pm->v2[i]+1, pm->v3[i]+1);
    for (i=0; i<pm->nEdge; i++) fprintf(f, "%d %d\n", pm->edge[i]->v1+1, pm->edge[i]->v2+1);
    for (i=0; i<pm->nEdge; i++) fprintf(f, "%d\n", pm->edge[i]->v1+1);
    fprintf(f, "%d A\n", pm->edge[0]->v1+1);
    fprintf(f, "%d B\n", pm->edge[0]->v2+1);
    fclose(f);
    return num++;
}

static int dump3(surface_mesh *pm) {
    static int num=0;
    char fname[1024];
    FILE *f;
    int i;
    sprintf(fname, "_dump_frt3_smv.%03d", num);
    f = fopen(fname, "w");

    fprintf(f, "%d %d %d %d 1\n", pm->nPoint+1, pm->nTria, pm->nEdge, pm->nSphereVert+3);
    for (i=0; i<=pm->nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", pm->vert[i].x, pm->vert[i].y, pm->vert[i].z);
    for (i=0; i<pm->nTria; i++) fprintf(f, "%d %d %d\n", pm->v1[i]+1, pm->v2[i]+1, pm->v3[i]+1);
    for (i=0; i<pm->nEdge; i++) fprintf(f, "%d %d\n", pm->edge[i]->v1+1, pm->edge[i]->v2+1);
    for (i=0; i<pm->nSphereVert; i++) fprintf(f, "%d\n", pm->sphereVert[i]+1);
    fprintf(f, "%d A\n", pm->edge[0]->v1+1);
    fprintf(f, "%d B\n", pm->edge[0]->v2+1);
    fprintf(f, "%d C\n", (pm->sworkVert >= 0) ? pm->sworkVert+1 : pm->nPoint+1);
    fclose(f);
    return num++;
}


void makeTria(surface_mesh *pm) {
    int r = 0;
    /*int temp_i;*/
#ifdef PROF
    time_t tr1, tr2, tr;
#endif

    pm->fresh = 0;
#ifdef PROF
    tr1 = clock();
#endif
    if (tria_dump_front) dump(pm);
#ifdef PROF
    tr2 = clock();
    tr = tr2-tr1;
#endif
    pm->suMax = pm->uMax,  pm->suMin = pm->uMin;
    pm->svMax = pm->vMax,  pm->svMin = pm->vMin;
    /*for (temp_i=0; temp_i<32; temp_i++)  temp_hist[temp_i] = 0;*/
    while (pm->nEdge > 0) {

	r = newTria(pm);
	if (r)  break;
#ifdef SHOWPROGRESS
	printf("\r nP = %5d  nT = %5d  nE = %5d", pm->nPoint, pm->nTria, pm->nEdge);
	fflush(stdout);
#endif
	if (tria_debug_front) dump(pm);
    }
    if (r==-20) {
	while (pm->nEdge > 0) {
	    remFace32(pm, pm->edge[pm->nEdge-1]);
	}
	r = 0;
    }
    if (!r) {
#ifdef PROF
	tr1 = clock();
	tr += tr1-tr2;
#ifdef SHOWPROGRESS
	printf("\r\t\t\t\t\t\t[%7.2lfs, tria/sec %5.0lf]\r", ((double)tr1-tr2)/CLOCKS_PER_SEC, pm->nTria / (((double)tr1-tr2)/CLOCKS_PER_SEC));
#endif
#endif
#ifdef SHOWPROGRESS
	printf("\r nP = %5d  nT = %5d  done.        \n", pm->nPoint, pm->nTria);
#endif
	/*for (temp_i=0; temp_i<32; temp_i++)  printf("%d ", temp_hist[temp_i]);  printf("\n");*/
	if (tria_final_front) dump(pm);
	pm->uMax = pm->suMax,  pm->uMin = pm->suMin;
	pm->vMax = pm->svMax,  pm->vMin = pm->svMin;
    }
    return;
} /*makeTria*/
