#include "common.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tria32.h"
#include "refine3.h"
#include "region3.h"
#include "user3.h"
#include "det.h"

int region_dump_face = 0;

static double precision = 1e-3;

double  surfacesizeratio = 0.05;

double sizeFace(surface_mesh *pm, double x, double y, double z) {
    return userSizeFace(pm, x,y,z);
}


int bounSurf0(surface_mesh *pm, double u, double v, double *x, double *y, double *z) {
    x[0] = pm->x0 + u*pm->x1 + v*pm->x2;
    y[0] = pm->y0 + u*pm->y1 + v*pm->y2;
    z[0] = pm->z0 + u*pm->z1 + v*pm->z2;
    return 1;
}


void V_U0(surface_mesh *pm, double u, double *v) {
    (void)u,  (void) pm;
    v[0] = 0.;
    return;
}


static double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/

static double normvec(double x1, double y1, double z1, double x2, double y2, double z2, double *px, double *py, double *pz) {
	double w1[4], w2[4], w3[4], x, y, z;
	w1[0] = x1, w2[0] = x2, w3[0] = 0.0;
	w1[1] = y1, w2[1] = y2, w3[1] = 0.0;
	w1[2] = z1, w2[2] = z2, w3[2] = 0.0;
	w1[3] = x1, w2[3] = x2, w3[3] = 0.0;
	*px = x = orient2d(w3+1, w2+1, w1+1);
	*py = y = orient2d(w3+2, w2+2, w1+2);
	*pz = z = orient2d(w3+0, w2+0, w1+0);
	return x*x + y*y + z*z;
}


static double nextT(surface_mesh *pm, int iSurf, int iV_U, double t) {
    int i=0;
    double u, v, t0, t1, x0, y0, z0, x1, y1, z1, x, y, z, s1=0., s, er;

    t0 = t;
    t1 = pm->tEnd;
    bounLine(pm, iV_U, t0, &u, &v);
    bounSurf(pm, iSurf, u, v, &x0, &y0, &z0);
    bounLine(pm, iV_U, t1, &u, &v);
    bounSurf(pm, iSurf, u, v, &x1, &y1, &z1);

    s = sizeFace(pm, 0.5*(x0+x1), 0.5*(y0+y1), 0.5*(z0+z1));
    er = 0.01*s;
    if (distance(x0, y0, z0, x1, y1, z1) <= s)
	return  -1.0e11;
    while (fabs(s1-s) >= er) {
	t = 0.5*(t0+t1);
	bounLine(pm, iV_U, t, &u, &v);
	bounSurf(pm, iSurf, u, v, &x, &y, &z);
	s1 = distance(x0, y0, z0, x, y, z);
	s = sizeFace(pm, 0.5*(x0+x), 0.5*(y0+y), 0.5*(z0+z));
	er = 0.01*s;
	if (s1>s) t1=t; else t0=t;
	i++;
	if (i > 100) {
	    if (fabs(s1-s) >= 10.*er)
		errorExit3(3, "i>100 in  nextT");
	    break;
	}
    }

    return t;
}

static void smoothingLine(surface_mesh *pm, int *vert, double *pT, double *pU, double *pV, int iSurf, int iV_U, int n) {
    int i, j=0, jj;
    double t, tp, tn, tu, tv, sp, sn;
    double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2;
    double du, k=0.0, k0, k1, fi0_, fi1_, f1, p;
    int v, vp, vn;

    while (j < 5) {
	for (i=n-2; i>0; i--) {
	    v  = vert[i];
	    vp = vert[i-1];  tp = pT[i-1];
	    vn = vert[i+1];  tn = pT[i+1];
	    if (tp > tn) {
		t  = tp;
		tp = tn;
		tn = t;
	    }
	    x1 = pm->vert[vp].x;
	    y1 = pm->vert[vp].y;
	    z1 = pm->vert[vp].z;
	    x2 = pm->vert[vn].x;
	    y2 = pm->vert[vn].y;
	    z2 = pm->vert[vn].z;

	    sp = sizeFace(pm, 0.5*(pm->vert[v].x+x1), 0.5*(pm->vert[v].y+y1), 0.5*(pm->vert[v].z+z1));
	    sn = sizeFace(pm, 0.5*(pm->vert[v].x+x2), 0.5*(pm->vert[v].y+y2), 0.5*(pm->vert[v].z+z2));

	    /* find  point  on  surface */
	    t = pT[i];
	    bounLine(pm, iV_U, t, &tu, &tv);
	    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
	    f1 = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
		(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
	    for (jj=0; jj<222; jj++) { /* NEW METH */
		if (f1 < 0.01) break;
		du = 0.2; /*???*/
		do {
		    du /= 2.0;
		    if (iSurf>=0) {
			bounLine(pm, iV_U, t+du, &tu, &tv);
			bounSurf(pm, iSurf, tu, tv, &x0, &y0, &z0);
		    } else {
#ifdef CGM
			CGM_GetEdgeCoordsFromU(edge, t+du, xyz);
			x0 = xyz[0]; y0 = xyz[1]; z0 = xyz[2];
#endif
		    }
		    p = distanceS(x, y, z, x0, y0, z0);
		} while (p > 0.1*sp*sn);

		k0 = -0.9;
		k1 = 1.0;
		if (du == 0.0)
		    errorExit3(3, "du == 0.0 ");
		if (t + k0*du < tp)
		    k0 = (tp-t)/du;
		if (t + k1*du > tn)
		    k1 = (tn-t)/du;
		while (k1 - k0 > 0.00001) {/*min of funk*/
		    k = 0.5*(k0+k1);
		    bounLine(pm, iV_U, t+k*du, &tu, &tv);
		    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		    fi0_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
			(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
		    bounLine(pm, iV_U, t+(k+fabs(k)*0.01)*du, &tu, &tv);
		    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		    fi1_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
			(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
		    if (fi1_ - fi0_ > 0.)
			k1 = k;
		    else
			k0 = k;
		}
		t += k*du;
		bounLine(pm, iV_U, t, &tu, &tv);
		bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		f1 = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
		    (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
	    } /*NEW METH*/
	    bounLine(pm, iV_U, t, &tu, &tv);
	    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
	    pU[i] = tu,  pV[i] = tv,  pT[i] = t;
	    pm->vert[v].x = x;
	    pm->vert[v].y = y;
	    pm->vert[v].z = z;
	}
	j++;
    }
    return;
}

static void adjustingLine(surface_mesh *pm, int *vert, double *pT, double *pU, double *pV, int iSurf, int iV_U, int n) {
    int i, j;
    double t, tp, tn, ta, tb, tu, tv;
    double x, y, z, x1, y1, z1, x2, y2, z2, xa, ya, za, xb, yb, zb, xt, yt, zt;
    int v, vp, vn;

    vp = vert[0];   tp = pT[0];
    vn = vert[n-1]; tn = pT[n-1];
    x1 = pm->vert[vp].x;
    y1 = pm->vert[vp].y;
    z1 = pm->vert[vp].z;
    x2 = pm->vert[vn].x;
    y2 = pm->vert[vn].y;
    z2 = pm->vert[vn].z;

    if (distance(x2, y2, z2, x1, y1, z1) == 0)
	errorExit3(2, " bounds for V_U parametrization point to the same node ");

    t = pT[0];
    bounLine(pm, iV_U, t, &tu, &tv);
    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
    pU[0] = tu,  pV[0] = tv;
    if (distance(x1, y1, z1, x, y, z) > precision) {
	printf("x1 = %12.11lf y1 = %12.11lf z1 = %12.11lf \n", x1, y1, z1);
	printf("x2 = %12.11lf y2 = %12.11lf z2 = %12.11lf \n", x2, y2, z2);
	printf("x = %12.11lf y = %12.11lf z = %12.11lf \n", x, y, z);
	errorExit3(2, " bounds for V_U parametrization do not correspond, pU[0] "); 
    }
    t = pT[n-1];
    bounLine(pm, iV_U, t, &tu, &tv);
    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
    pU[n-1] = tu,  pV[n-1] = tv;
    if (distance(x2, y2, z2, x, y, z) > precision) {
	printf("x1 = %12.11lf y1 = %12.11lf z1 = %12.11lf \n", x1, y1, z1);
	printf("x2 = %12.11lf y2 = %12.11lf z2 = %12.11lf \n", x2, y2, z2);
	printf("x = %12.11lf y = %12.11lf z = %12.11lf \n", x, y, z);
	errorExit3(2, " bounds for V_U parametrization do not correspond, pU[last] ");
    } 

    for (i=n-2; i>0; i--) {
	v  = vert[i];
	x  = pm->vert[v].x;
	y  = pm->vert[v].y;
	z  = pm->vert[v].z;
	xa = x1;
	ya = y1;
	za = z1;
	xb = x2;
	yb = y2;
	zb = z2;
	ta = tp;
	tb = tn;

	for (j=0; j<40; j++) {
	    if (distance(xa, ya, za, x, y, z) < precision) {
		t = ta; break;
	    }
	    if (distance(xb, yb, zb, x, y, z) < precision) {
		t = tb; break;
	    }
	    t = (ta + tb)/2.0;
	    bounLine(pm, iV_U, t, &tu, &tv);
	    bounSurf(pm, iSurf, tu, tv, &xt, &yt, &zt);
	    if (distance(xa, ya, za, x, y, z) < distance(xb, yb, zb, x, y, z)) {
		if (  (distance(xa, ya, za, xt, yt, zt) > distance(xa, ya, za, x, y, z)) &&
			(distance(xa, ya, za, xt, yt, zt) > distance(xt, yt, zt, x, y, z))  ) {
		    tb = t; xb = xt; yb = yt; zb = zt;
		} else {
		    ta = t; xa = xt; ya = yt; za = zt;
		}
	    } else {
		if (  (distance(xt, yt, zt, xb, yb, zb) > distance(xt, yt, zt, x, y, z)) &&
			(distance(xt, yt, zt, xb, yb, zb) > distance(xb, yb, zb, x, y, z))  ) {
		    ta = t; xa = xt; ya = yt; za = zt;
		} else {
		    tb = t; xb = xt; yb = yt; zb = zt;
		}
	    }
	}
	bounLine(pm, iV_U, t, &tu, &tv);
	pU[i] = tu,  pV[i] = tv,  pT[i] = t;
    }
    return;
}


void makeAFLine(surface_mesh *pm, int *vert, int nVert, int bInverse) {
	int i;

	for (i=0; i<nVert-1; i++) {
		if (vert[i]==vert[i+1]) continue; // cheating ;-)
		if (bInverse!=0) {
			addFace32(pm, vert[i+1], vert[i]);
//			printf("edge %d %d\n", vert[i+1], vert[i]);
		}
		if (bInverse!=1) {
			addFace32(pm, vert[i], vert[i+1]);
//			printf("edge %d %d\n", vert[i], vert[i+1]);
		}
	}
	return;
} /*makeAFLine*/

static int dump(surface_mesh *pm) {
	static int num=0;
	char fname[1024];
	FILE *f;
	int i;
	sprintf(fname, "_dump_region_smv.%03d", num);
	f = fopen(fname, "w");

	fprintf(f, "%d %d 0 0\n", pm->nPoint, pm->nTria);
	for (i=0; i<pm->nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", pm->vert[i].x, pm->vert[i].y, pm->vert[i].z);
	for (i=0; i<pm->nTria; i++) fprintf(f, "%d %d %d\n", pm->v1[i]+1, pm->v2[i]+1, pm->v3[i]+1);
	fclose(f);
	return num++;
}

static int writeBound(surface_mesh *pm) {
    if (region_dump_face) dump(pm);
    return pm->nTria;
}


void initAFS_ (
	surface_mesh *pm,
	int *pnVVert, double *VVertxyz,
	int *pnLine, int *LineD, int *LineP, double *LineT,
	int *pnSurface, int *SurfL, int *SurfI, double *SurfT
) {
	int i, j, k, iLine, iSub, m;

	int nVVert = *pnVVert;
	int nLine = *pnLine;
	int nSurface = *pnSurface;

	int vBegin, vEnd, iSurf, iNorm, iV_U;
	double tBegin, tEnd;
	double uMin, uMax, vMin, vMax;
	double pmin[3]={0,0,0}, pmax[3]={0,0,0};

	int bInverse, bCut, iLabel, nBoundTria=0;

	int nSub, iSurfLine;
	int nVert;
	double x, y, z, x1, y1, z1, u, v, t;
	int *boundVert, *vVert;
	double *boundT, *boundU, *boundV;
	StrucLine3 *line;
	StrucSurface *surface;

	init(pm);
	detinit();

	pm->boxcx = 0.5;
	pm->boxcy = 0.5;
	pm->boxcz = 0.5;
	pm->boxsize = 0.5;

	boundVert = (int *)myAlloc(MAX1 * sizeof(int));
	boundT = (double *)myAlloc(MAX1 * sizeof(double));
	boundU = (double *)myAlloc(MAX1 * sizeof(double));
	boundV = (double *)myAlloc(MAX1 * sizeof(double));

	vVert = (int *)myAlloc(nVVert * sizeof(int));

	for (i=0; i<nVVert; i++) {
		x = VVertxyz[3*i+0];
		y = VVertxyz[3*i+1];
		z = VVertxyz[3*i+2];
		if (i==0) {
			pmin[0] = pmax[0] = x;
			pmin[1] = pmax[1] = y;
			pmin[2] = pmax[2] = z;
		} else {
			if (pmin[0]>x) pmin[0]=x;
			if (pmax[0]<x) pmax[0]=x;
			if (pmin[1]>y) pmin[1]=y;
			if (pmax[1]<y) pmax[1]=y;
			if (pmin[2]>z) pmin[2]=z;
			if (pmax[2]<z) pmax[2]=z;
		}
		vVert[i] = pm->nPoint;
		addPoint(pm, x, y, z);
	}
	
//	printf("BOX:%lf %lf %lf -- %lf %lf %lf\n", pmin[0], pmin[1], pmin[2], pmax[0], pmax[1], pmax[2]);
	pm->S0 = (pmax[0]-pmin[0]);
	if (pmax[1]-pmin[1] > pm->S0)  pm->S0 = pmax[1]-pmin[1];
	if (pmax[2]-pmin[2] > pm->S0)  pm->S0 = pmax[2]-pmin[2];
	pm->S0 *= surfacesizeratio;
//	printf("S0: %lf\n", S0);
	
	pm->boxcx = 0.5*(pmin[0]+pmax[0]);
	pm->boxcy = 0.5*(pmin[1]+pmax[1]);
	pm->boxcz = 0.5*(pmin[2]+pmax[2]);
	pm->boxsize = pmax[0]-pmin[0];
	if (pmax[1]-pmin[1]>pm->boxsize) pm->boxsize = pmax[1]-pmin[1];
	if (pmax[2]-pmin[2]>pm->boxsize) pm->boxsize = pmax[2]-pmin[2];


	line = (StrucLine3 *)myAlloc(2*nLine * sizeof(StrucLine3) );

	for (iSub=0, iLine=0; iLine<nLine; iLine++) {
		vBegin = LineD[3*iLine+0];
		vEnd   = LineD[3*iLine+1];
		nSub   = LineD[3*iLine+2];
		line[iLine].vBegin = vBegin - 1;
		line[iLine].vEnd   = vEnd - 1;
		line[iLine].nSub   = nSub;
		line[iLine].sub = (StrucSub *)myAlloc(nSub * sizeof(StrucSub));
		for (i=0; i<nSub; i++, iSub++) {
			iSurf = LineP[2*iSub+0];
			iV_U  = LineP[2*iSub+1];
			tBegin = LineT[2*iSub+0];
			tEnd   = LineT[2*iSub+1];
			line[iLine].sub[i].iSurf = iSurf;
			line[iLine].sub[i].iV_U = iV_U;
			line[iLine].sub[i].tBegin = tBegin;
			line[iLine].sub[i].tEnd = tEnd;
		}
	}

	/* for  mooving  vVert ... */
	for (iLine=0; iLine<nLine; iLine++) {
		iSurf = line[iLine].sub[0].iSurf;
		if (iSurf == 0)
			continue;
		vBegin = line[iLine].vBegin;
		vEnd = line[iLine].vEnd;
		tBegin = line[iLine].sub[0].tBegin;
		tEnd = line[iLine].sub[0].tEnd;
		iV_U = line[iLine].sub[0].iV_U;

		bounLine(pm, iV_U, tBegin, &u, &v);
		bounSurf(pm, iSurf, u, v, &x, &y, &z);
		pm->vert[vBegin].x = x;
		pm->vert[vBegin].y = y;
		pm->vert[vBegin].z = z;
		pm->vert[vBegin].u = u;
		pm->vert[vBegin].v = v;

		bounLine(pm, iV_U, tEnd, &u, &v);
		bounSurf(pm, iSurf, u, v, &x, &y, &z);
		pm->vert[vEnd].x = x;
		pm->vert[vEnd].y = y;
		pm->vert[vEnd].z = z;
		pm->vert[vEnd].u = u;
		pm->vert[vEnd].v = v;
	}

	surface = (StrucSurface *)myAlloc(nSurface * sizeof(StrucSurface));

	for (iSurfLine=0, iLine=0; iLine<nSurface; iLine++) {
		i      = SurfL[5*iLine+0];
		iSurf  = SurfL[5*iLine+1];
		iLabel = SurfL[5*iLine+2];
		bCut   = SurfL[5*iLine+3];
		iNorm  = SurfL[5*iLine+4];
		uMin = SurfT[4*iLine+0];
		uMax = SurfT[4*iLine+1];
		vMin = SurfT[4*iLine+2];
		vMax = SurfT[4*iLine+3];
		surface[iLine].nLine = i;
		surface[iLine].line = (int *)myAlloc(i * sizeof(int));
		surface[iLine].inverse = (int *)myAlloc(i * sizeof(int));

		for (j=0; j<i; j++, iSurfLine++) {
			k        = SurfI[2*iSurfLine+0];
			bInverse = SurfI[2*iSurfLine+1];
			surface[iLine].line[j] = k - 1;
			surface[iLine].inverse[j] = bInverse;
		}
		surface[iLine].iSurf = iSurf;
		surface[iLine].iLabel = iLabel;
		surface[iLine].bCut = bCut;
		surface[iLine].iNorm = iNorm;
		surface[iLine].uMin = uMin;
		surface[iLine].uMax = uMax;
		surface[iLine].vMin = vMin;
		surface[iLine].vMax = vMax;
	}


	/* lines into  sides */
	for (iLine=0; iLine<nLine; iLine++) {
		vBegin = line[iLine].vBegin;
		vEnd = line[iLine].vEnd;
		iSurf = line[iLine].sub[0].iSurf;
		if (iSurf > 0) {
			pm->tBegin = tBegin = line[iLine].sub[0].tBegin;
			pm->tEnd = tEnd = line[iLine].sub[0].tEnd;
			iV_U = line[iLine].sub[0].iV_U;
		} else {/*iSurf == 0*/
			pm->tBegin = tBegin = 0.;
			pm->tEnd = tEnd = 1.;
			iV_U = 0;
			pm->x0 = pm->vert[vBegin].x;
			pm->y0 = pm->vert[vBegin].y;
			pm->z0 = pm->vert[vBegin].z;
			pm->x1 = pm->vert[vEnd].x - pm->x0;
			pm->y1 = pm->vert[vEnd].y - pm->y0;
			pm->z1 = pm->vert[vEnd].z - pm->z0;
			pm->x2 = 0.;
			pm->y2 = 0.;
			pm->z2 = 0.;
		}

		nVert = 0;
		bounLine(pm, iV_U, tBegin, &u, &v);
		boundT[nVert] = tBegin;
		boundU[nVert] = u;
		boundV[nVert] = v;
		boundVert[nVert++] = vVert[vBegin];
		t = tBegin;
		for (j=0; ; j++) {
			if (nVert >= MAX1)
				errorExit3(4, "MAX1");
			t = nextT(pm, iSurf, iV_U, t);
			if (t < -1.e10)
				break;

			bounLine(pm, iV_U, t, &u, &v);
			pm->vert[pm->nPoint].u = u;
			pm->vert[pm->nPoint].v = v;
			bounSurf(pm, iSurf, u, v, &x, &y, &z);
			boundT[nVert] = t;
			boundU[nVert] = u;
			boundV[nVert] = v;
			boundVert[nVert++] = pm->nPoint;
			addPoint(pm, x, y, z);
		}/* for(j) */
		bounLine(pm, iV_U, tEnd, &u, &v);
		boundT[nVert] = tEnd;
		boundU[nVert] = u;
		boundV[nVert] = v;
		boundVert[nVert++] = vVert[vEnd];
		smoothingLine(pm, boundVert, boundT, boundU, boundV, iSurf, iV_U, nVert); 

		line[iLine].nVert = nVert;
		line[iLine].vert = (int *)myAlloc(nVert * sizeof(int));
		line[iLine].sub[0].u = (double *)myAlloc(nVert * sizeof(double));
		line[iLine].sub[0].v = (double *)myAlloc(nVert * sizeof(double));
		for (j=0; j<nVert; j++) {
			line[iLine].vert[j] = boundVert[j];
			line[iLine].sub[0].u[j] = boundU[j];
			line[iLine].sub[0].v[j] = boundV[j];
		}

		if (pm->exportCurves) {
		    if (pm->exportCurves[iLine] > 0) {
			for (j=0; j<nVert-1; j++) {
			    if (boundVert[j]==boundVert[j+1])  continue;
			    if (pm->nexpEdge >= pm->nnexpEdge) {
				printf("max nE is small\n");
				break;
			    }
			    pm->expEdge[2*pm->nexpEdge+0] = boundVert[j];
			    pm->expEdge[2*pm->nexpEdge+1] = boundVert[j+1];
			    pm->expEdgeColor[pm->nexpEdge] = pm->exportCurves[iLine];
			    pm->nexpEdge++;
			}
		    }
		}

		for (iSub=1; iSub<line[iLine].nSub; iSub++) {
			vBegin = line[iLine].vBegin;
			vEnd = line[iLine].vEnd;
			pm->tBegin = tBegin = line[iLine].sub[iSub].tBegin;
			pm->tEnd = tEnd = line[iLine].sub[iSub].tEnd;
			iV_U = line[iLine].sub[iSub].iV_U;
			iSurf = line[iLine].sub[iSub].iSurf;


			boundT[0] = tBegin;
			boundT[nVert-1] = tEnd;

			adjustingLine(pm, boundVert, boundT, boundU, boundV, iSurf, iV_U, nVert);
			line[iLine].sub[iSub].u = (double *)myAlloc(nVert * sizeof(double));
			line[iLine].sub[iSub].v = (double *)myAlloc(nVert * sizeof(double));
			for (j=0; j<nVert; j++) {
				line[iLine].sub[iSub].u[j] = boundU[j];
				line[iLine].sub[iSub].v[j] = boundV[j];
			}
		}/*for iSub < nSub*/
	}/* for iLine */
	pm->nLinePoint = pm->nPoint;

	pm->maxTria = pm->maxFace;
	pm->v1 = (int *)myAlloc(pm->maxTria * sizeof(int));
	pm->v2 = (int *)myAlloc(pm->maxTria * sizeof(int));
	pm->v3 = (int *)myAlloc(pm->maxTria * sizeof(int));

	pm->countCrvT = 0;
	pm->cbnd = (int *)myAlloc(pm->maxTria * sizeof(int));
	pm->ciSurf = (int *)myAlloc(pm->maxTria * sizeof(int));
	pm->cu1 = (double *)myAlloc(pm->maxTria * sizeof(double));
	pm->cu2 = (double *)myAlloc(pm->maxTria * sizeof(double));
	pm->cu3 = (double *)myAlloc(pm->maxTria * sizeof(double));
	pm->cv1 = (double *)myAlloc(pm->maxTria * sizeof(double));
	pm->cv2 = (double *)myAlloc(pm->maxTria * sizeof(double));
	pm->cv3 = (double *)myAlloc(pm->maxTria * sizeof(double));


	for (iSurf=0; iSurf<nSurface; iSurf++) {
		pm->nEdge = 0;
		nLine = surface[iSurf].nLine;
		iLabel = surface[iSurf].iLabel;
		bCut = surface[iSurf].bCut;
		prepTree32(pm, pm->nnEdge);

		pm->nTria = 0;
		pm->uMin = surface[iSurf].uMin;
		pm->uMax = surface[iSurf].uMax;
		pm->vMin = surface[iSurf].vMin;
		pm->vMax = surface[iSurf].vMax;
		pm->iSurf = surface[iSurf].iSurf;
		if (pm->iSurf > 0)
			pm->iNorm = surface[iSurf].iNorm;
		else
			pm->iNorm = 0;
		for (iLine=0; iLine<nLine; iLine++) {
			i = surface[iSurf].line[iLine];
			j = line[i].nSub;
			if ( (pm->iSurf>0) && (j>=1) ) {/*reload u&v*/
				iSub = -1;
				for (k=0; k<j; k++) {
					if (line[i].sub[k].iSurf == pm->iSurf) {
						iSub = k;
						break;
					}
				}
				if (iSub == -1)
					errorExit3(3, " iSub == 0");
				for (k=0; k<line[i].nVert; k++) {
					pm->vert[ line[i].vert[k] ].u = line[i].sub[iSub].u[k];
					pm->vert[ line[i].vert[k] ].v = line[i].sub[iSub].v[k];
				}
			}
			bInverse = surface[iSurf].inverse[iLine];
			makeAFLine(pm, line[i].vert, line[i].nVert, bInverse);
		}/*for(iLine<nLine)*/

		if (pm->iSurf == 0) {/*init  for  plane*/
			pm->uMin = 0.;
			pm->uMax = 1.;
			pm->vMin = 0.;
			pm->vMax = 1.;
			i = surface[iSurf].line[0];
			x = 0.0;
			y = 0.0;
			z = 0.0;
			for (k=0; k<pm->nEdge; k++) {
			    pm->vv1 = pm->edge[k]->v1;
			    pm->vv2 = pm->edge[k]->v2;
			    normvec(pm->vert[pm->vv1].x, pm->vert[pm->vv1].y, pm->vert[pm->vv1].z,  pm->vert[pm->vv2].x, pm->vert[pm->vv2].y, pm->vert[pm->vv2].z, &x1, &y1, &z1);
			    x+=x1, y+=y1, z+=z1;
			}

/*			for (iLine=0; iLine<nLine; iLine++) {
				i = surface[iSurf].line[iLine];
				x1 = 0.0,  y1 = 0.0,  z1 = 0.0;
				for (k=0; k<line[i].vert[line[i].nVert]; k++) {
				}
				pm->vv1 = line[i].vert[0];
				pm->vv2 = line[i].vert[ line[i].nVert - 1 ];
				normvec(pm->vert[pm->vv1].x, pm->vert[pm->vv1].y, pm->vert[pm->vv1].z,  pm->vert[pm->vv2].x, pm->vert[pm->vv2].y, pm->vert[pm->vv2].z, &x1, &y1, &z1);
				if (surface[iSurf].inverse[iLine]) {
					x-=x1, y-=y1, z-=z1;
				} else {
					x+=x1, y+=y1, z+=z1;
				}
			}*/
			u = sqrt(x*x + y*y + z*z);
			if (u == 0.0)
				errorExit3(3, " u == 0.0   in  surfNormal ");
			x /= u;
			y /= u;
			z /= u;
			pm->vv0 = line[i].vert[0];
			pm->vv1 = line[i].vert[ line[i].nVert - 1 ];
			if (surface[iSurf].inverse[0] > 0) {/*for  right  normal*/
				k = pm->vv1;
				pm->vv1 = pm->vv0;
				pm->vv0 = k;
			}
			x1 = pm->vert[pm->vv1].x - pm->vert[pm->vv0].x;
			y1 = pm->vert[pm->vv1].y - pm->vert[pm->vv0].y;
			z1 = pm->vert[pm->vv1].z - pm->vert[pm->vv0].z;
			u = sqrt(x1*x1 + y1*y1 + z1*z1);
			if (u == 0.0)
				errorExit3(3, " u == 0.0   in  surfNormal ");
			x1 /= u;
			y1 /= u;
			z1 /= u;

			pm->x0 = pm->vert[pm->vv0].x;
			pm->y0 = pm->vert[pm->vv0].y;
			pm->z0 = pm->vert[pm->vv0].z;
			pm->x1 = x1;
			pm->y1 = y1;
			pm->z1 = z1;
			normvec(x, y, z,  x1, y1, z1,  &pm->x2, &pm->y2, &pm->z2);

			for (iLine=0; iLine<nLine; iLine++) {/*reload u&v*/
				i = surface[iSurf].line[iLine];
				for (j=0; j<line[i].nVert; j++) {
					m = line[i].vert[j];
					rePlane(pm, m);
					if (pm->vert[m].u > pm->uMax) pm->uMax = pm->vert[m].u;
					if (pm->vert[m].u < pm->uMin) pm->uMin = pm->vert[m].u;
					if (pm->vert[m].v > pm->vMax) pm->vMax = pm->vert[m].v;
					if (pm->vert[m].v < pm->vMin) pm->vMin = pm->vert[m].v;
				}
			}
		}/*if( pm->iSurf == 0 )*/

		makeTria(pm);
		for (i=0; i<10; i++) {
			smoothingSurf(pm); 
		}
		nBoundTria += writeBound(pm);

		/* make AF for surface  */
		for (i=0; i<pm->nTria; i++) {
			addFace(pm, pm->v1[i], pm->v2[i], pm->v3[i], 0, iLabel);
			if (bCut)
				addFace(pm, pm->v1[i], pm->v3[i], pm->v2[i], 1, bCut);
		}
	}/*for(iSurf<nSurface)*/

	pm->nBoundPoint = pm->nPoint;
//	outBound(nBoundTria);
	free(boundVert);
	free(boundT);
	free(boundU);
	free(boundV);
	free(vVert);
	nLine = *pnLine;
	for (iLine=0; iLine<nLine; iLine++) {
		for (i=0; i<line[iLine].nSub; i++) {
			free(line[iLine].sub[i].u);
			free(line[iLine].sub[i].v);
		}
		free(line[iLine].vert);
		free(line[iLine].sub);
	}
	free(line);
	for (iLine=0; iLine<nSurface; iLine++) {
		free(surface[iLine].line);
		free(surface[iLine].inverse);
	}
	free(surface);
	free(pm->v1);
	free(pm->v2);
	free(pm->v3);
	free(pm->cbnd);
	free(pm->ciSurf);
	free(pm->cu1);
	free(pm->cu2);
	free(pm->cu3);
	free(pm->cv1);
	free(pm->cv2);
	free(pm->cv3);
	return;
} /* makeAFSurface */


int addTria(surface_mesh *pm, int v1, int v2, int v3) {
    if (pm->nTria >= pm->maxTria)  return printf("nTria >= maxTria\n"), -1;
    pm->v1[pm->nTria] = v1;
    pm->v2[pm->nTria] = v2;
    pm->v3[pm->nTria] = v3;
    return pm->nTria++;
} /*addTria*/

