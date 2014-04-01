#include "cgmwrap.h"
#include "common.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tria32.h"
#include "refine3.h"
#include "region3c.h"
#include "user3c.h"

#define VERBOSE_CONSTRUCT 0

int region_dump_face_cgm = 0;

extern double  surfacesizeratio;

static int dump(surface_mesh *pm) {
    static int num=0;
    char fname[1024];
    FILE *f;
    int i;
    sprintf(fname, "_dump_region_cgm_smv.%03d", num);
    f = fopen(fname, "w");
    fprintf(f, "%d %d 0 0\n", pm->nPoint, pm->nTria);
    for (i=0; i<pm->nPoint; i++)
	fprintf(f, "%20.15lf %20.15lf %20.15lf\n", pm->vert[i].x, pm->vert[i].y, pm->vert[i].z);
    for (i=0; i<pm->nTria; i++) fprintf(f, "%d %d %d\n", pm->v1[i]+1, pm->v2[i]+1, pm->v3[i]+1);
    fclose(f);
    return num++;
}

static double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/










static int nextUinEdge(surface_mesh *pm, CGMedge edge, double t, double *tout) {
    int i=0;
    double t0, t1, s1, s, er, torg;
    double p[3], p0[3], p1[3];

    torg = t0 = t;
    t1 = pm->tEnd;
    CGM_GetEdgeCoordsFromU(edge, t0, p0);
    CGM_GetEdgeCoordsFromU(edge, t1, p1);

    s = sizeFace(pm, 0.5*(p0[0]+p1[0]), 0.5*(p0[1]+p1[1]), 0.5*(p0[2]+p1[2]));
    er = 0.01*s;
    s1 = CGM_EdgeArcLengthParam(edge, t0, t1);
    if (s1 <= s) return 1;
    s1 = 0.0;
    while (fabs(s1-s) >= er) {
	t = 0.5*(t0+t1);
	CGM_GetEdgeCoordsFromU(edge, t, p);
	s1 = CGM_EdgeArcLengthParam(edge, torg, t);
	/*s1 = distance(p0[0], p0[1], p0[2], p[0], p[1], p[2]);*/
	s = sizeFace(pm, 0.5*(p0[0]+p[0]), 0.5*(p0[1]+p[1]), 0.5*(p0[2]+p[2]));
	er = 0.01*s;
	if (s1>s) t1=t; else t0=t;
	i++;
	if (i > 100) {
	    if (fabs(s1-s) >= 10.*er)
		errorExit3(3, "i>100 in  nextUinEdge");
	    break;
	}
    }
    *tout = t;

    return 0;
}

void smoothingLineCGM(surface_mesh *pm, int *vert, double *pT, double *pU, double *pV, int iSurf, int iV_U, int n, CGMedge edge) {
    int i, j=0, jj;
    double t, tp, tn, tu, tv, sp, sn;
    double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2;
    double du, k=0.0, k0, k1, fi0_, fi1_, f1, p;
    int v, vp, vn;
    double xyz[3];

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
	    if (iSurf>=0) {
		bounLine(pm, iV_U, t, &tu, &tv);
		bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
	    } else {
		CGM_GetEdgeCoordsFromU(edge, t, xyz);
		x = xyz[0]; y = xyz[1]; z = xyz[2];
	    }
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
			CGM_GetEdgeCoordsFromU(edge, t+du, xyz);
			x0 = xyz[0]; y0 = xyz[1]; z0 = xyz[2];
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
		    if (iSurf>=0) {
			bounLine(pm, iV_U, t+k*du, &tu, &tv);
			bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		    } else {
			CGM_GetEdgeCoordsFromU(edge, t+k*du, xyz);
			x = xyz[0]; y = xyz[1]; z = xyz[2];
		    }
		    fi0_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
			(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
		    if (iSurf>=0) {
			bounLine(pm, iV_U, t+(k+fabs(k)*0.01)*du, &tu, &tv);
			bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		    } else {
			CGM_GetEdgeCoordsFromU(edge, t+(k+fabs(k)*0.01)*du, xyz);
			x = xyz[0]; y = xyz[1]; z = xyz[2];
		    }
		    fi1_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
			(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
		    if (fi1_ - fi0_ > 0.)
			k1 = k;
		    else
			k0 = k;
		}
		t += k*du;
		if (iSurf>=0) {
		    bounLine(pm, iV_U, t, &tu, &tv);
		    bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		} else {
		    CGM_GetEdgeCoordsFromU(edge, t, xyz);
		    x = xyz[0]; y = xyz[1]; z = xyz[2];
		}
		f1 = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
		    (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
	    } /*NEW METH*/
	    if (iSurf>=0) {
		bounLine(pm, iV_U, t, &tu, &tv);
		bounSurf(pm, iSurf, tu, tv, &x, &y, &z);
		pU[i] = tu,  pV[i] = tv;
	    } else {
		CGM_GetEdgeCoordsFromU(edge, t, xyz);
		x = xyz[0]; y = xyz[1]; z = xyz[2];
	    }
	    pT[i] = t;
	    pm->vert[v].x = x;
	    pm->vert[v].y = y;
	    pm->vert[v].z = z;
	}
	j++;
    }
    return;
}





static int writeBoundCGM(surface_mesh *pm) {
    if (region_dump_face_cgm) dump(pm);
    return pm->nTria;
}


void initAFSM_ (
	surface_mesh *pm,
	CGMmodel model
) {
	int i, j, k, iLine, iSurface, iVert;

	int nVVert = CGM_NumVertices(model);
	int nLine = CGM_NumEdges(model);
	int nSurface = CGM_NumFaces(model);
	
	int vBegin, vEnd, iSurf, iNorm, iV_U;
	double tBegin, tEnd;
	double uMin, uMax, vMin, vMax;
	

	int bInverse, bCut, iLabel, nBoundTria=0;

	int nSub;
	int nVert;
	double x, y, z, u, p[3], uv[2], pmin[3], pmax[3];
	int *boundVert, *vVert;
	double *boundT;
	StrucLine3 *line;
	StrucSurface *surface;
	CGMvertex vertex;
	CGMedge edge;
	CGMface face;
	
	init(pm);
	user3initCGM(pm);

	CGM_ModelBoundingBox(model, pmin, pmax);
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

	vVert = (int*)myAlloc(nVVert*sizeof(int));
	for (i=0; i<nVVert; i++) {
		vertex = CGM_ithVertex(model, i);
		CGM_GetVertexCoords(vertex, p);
		vVert[i] = pm->nPoint;
		addPoint(pm, p[0], p[1], p[2]);
	}
	
	line = (StrucLine3 *)myAlloc(2*nLine * sizeof(StrucLine3) );

	for (iLine=0; iLine<nLine; iLine++) {
		edge = CGM_ithEdge(model, iLine);
		vBegin = CGM_VertexIndex(model, CGM_EdgeStartVertex(edge));
		vEnd   = CGM_VertexIndex(model, CGM_EdgeEndVertex(edge));
		nSub   = 1;
		line[iLine].vBegin = vBegin;
		line[iLine].vEnd   = vEnd;
		line[iLine].nSub   = nSub;
		line[iLine].sub = (StrucSub *)myAlloc(nSub * sizeof(StrucSub));
		iSurf = -1;
		iV_U  = iLine;
		tBegin = CGM_EdgeStartParam(edge);
		tEnd   = CGM_EdgeEndParam(edge);
		line[iLine].sub[0].iSurf = iSurf;
		line[iLine].sub[0].iV_U = iV_U;
		line[iLine].sub[0].tBegin = tBegin;
		line[iLine].sub[0].tEnd = tEnd;
	}

	surface = (StrucSurface *)myAlloc(nSurface * sizeof(StrucSurface));

	for (iSurface=0; iSurface<nSurface; iSurface++) {
	        if (VERBOSE_CONSTRUCT)  printf("Surface %d:\n", iSurface);
		face = CGM_ithFace(model, iSurface);
		i      = CGM_NumEdgesInFace(face);
		iSurf  = -1;
		iLabel = iSurface + 1; // 
		bCut = -1;
		iNorm = 0;
		for (j=0; j<CGM_NumVolumes(model); j++) {
			k = CGM_VolumeFaceOrientation(CGM_ithVolume(model, j), face);
			if (VERBOSE_CONSTRUCT)  printf("Surface %d, volume %d, orientation: %d\n", iSurface, j, k);
			if (k<0)  {
//			    printf("AniFrtCad: CGM_VolumeFaceOrientation returned %d\n", k);
			    continue;
			}
			bCut++;
			iNorm = k;
		}
		if ((bCut<0) || (bCut>1)) {
			printf("Face #%d share %d volumes\n", iSurface, bCut+1);
			bCut = 0;
		}
		CGM_GetFaceParamRange(face, &uMin, &uMax, &vMin, &vMax);
		surface[iSurface].nLine = i;
		surface[iSurface].line = (int *)myAlloc(i * sizeof(int));
		surface[iSurface].inverse = (int *)myAlloc(i * sizeof(int));
//		printf("surface %d:\n", iSurface);
		for (j=0; j<i; j++) {
			edge = CGM_ithEdgeInFace(face, j);
			k        = CGM_EdgeIndex(model, edge);
			bInverse = CGM_FaceEdgeOrientation(face, edge);
			if (VERBOSE_CONSTRUCT)  printf("Surface %d, edge %d, orientation: %d\n", iSurface, j, bInverse);
//			if (bInverse < 0) {
//			    printf("AniFrtCad: CGM_FaceEdgeOrientation returned %d\nAssuming internal edge\n", bInverse);
//			}
//			if (iNorm && bInverse==0)  bInverse = 1;
//			if (iNorm && bInverse==1)  bInverse = 0;
			surface[iSurface].line[j] = k;
			surface[iSurface].inverse[j] = bInverse;
//			printf("\tline %d [%c]\n", k, (bInverse==1)?'*':((bInverse==0)?' ':'?'));
		}
		surface[iSurface].iSurf = iSurf;
		surface[iSurface].iLabel = iLabel;
		surface[iSurface].bCut = bCut;
		surface[iSurface].iNorm = iNorm;
		surface[iSurface].uMin = uMin;
		surface[iSurface].uMax = uMax;
		surface[iSurface].vMin = vMin;
		surface[iSurface].vMax = vMax;
	}
	
	boundVert = (int *)myAlloc(MAX1 * sizeof(int));
	boundT = (double *)myAlloc(MAX1 * sizeof(double));
	for (iLine=0; iLine<nLine; iLine++) {
		edge = CGM_ithEdge(model, iLine);
		vBegin = line[iLine].vBegin;
		vEnd = line[iLine].vEnd;
		iSurf = line[iLine].sub[0].iSurf;
		pm->tBegin = tBegin = line[iLine].sub[0].tBegin;
		pm->tEnd = tEnd = line[iLine].sub[0].tEnd;
		iV_U = line[iLine].sub[0].iV_U;
		nVert = 0;
		boundT[nVert] = tBegin;
		boundVert[nVert++] = vVert[vBegin];
		u = tBegin;
		for (j=0; ; j++) {
			if (nVert >= MAX1)
				errorExit3(4, "MAX1");
			if (nextUinEdge(pm, edge, u, &u)) break;

			pm->vert[pm->nPoint].u = u;
			CGM_GetEdgeCoordsFromU(edge, u, p);
			x = p[0];  y = p[1];  z = p[2];
			boundT[nVert] = u;
			boundVert[nVert++] = pm->nPoint;
			addPoint(pm, x, y, z);
		}/* for(j) */
		boundT[nVert] = tEnd;
		boundVert[nVert++] = vVert[vEnd];
		smoothingLineCGM(pm, boundVert, boundT, 0, 0, iSurf, iV_U, nVert, edge); 

		line[iLine].nVert = nVert;
		line[iLine].vert = (int *)myAlloc(nVert * sizeof(int));
		line[iLine].sub[0].u = (double *)myAlloc(nVert * sizeof(double));
		line[iLine].sub[0].v = 0;
		for (j=0; j<nVert; j++) {
			line[iLine].vert[j] = boundVert[j];
			line[iLine].sub[0].u[j] = boundT[j];
		}
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
		face = CGM_ithFace(model, iSurf);
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
		pm->iNorm = surface[iSurf].iNorm;
		pm->cgmface = face;
		for (iLine=0; iLine<nLine; iLine++) {
			i = surface[iSurf].line[iLine];
			j = line[i].nSub;
			for (k=0; k<line[i].nVert; k++) {
				iVert = line[i].vert[k];
				p[0] = pm->vert[iVert].x;
				p[1] = pm->vert[iVert].y;
				p[2] = pm->vert[iVert].z;
				CGM_GetFaceUVFromCoords(face, p, uv);
				pm->vert[iVert].u = uv[0];
				pm->vert[iVert].v = uv[1];
//				printf("Vertex %d: x=%lf, y=%lf, z=%lf, \tu=%lf, v=%lf\n", iVert,  pm->vert[iVert].x,  pm->vert[iVert].y,  pm->vert[iVert].z,  pm->vert[iVert].u, pm->vert[iVert].v);
			}
			bInverse = surface[iSurf].inverse[iLine];
			makeAFLine(pm, line[i].vert, line[i].nVert, bInverse);
		}/*for(iLine<nLine)*/

		makeTria(pm);
		for (i=0; i<5; i++) {
			smoothingSurf(pm); 
		}
		nBoundTria += writeBoundCGM(pm);

		/* make AF for surface  */
		for (i=0; i<pm->nTria; i++) {
			addFace(pm, pm->v1[i], pm->v2[i], pm->v3[i], 0, iLabel);
			if (bCut)
				addFace(pm, pm->v1[i], pm->v3[i], pm->v2[i], 1, iLabel);
		}
	}/*for(iSurf<nSurface)*/

	pm->nBoundPoint = pm->nPoint;
//	outBound(nBoundTria);
	free(boundVert);
	free(boundT);
	free(vVert);
	nLine = CGM_NumEdges(model);
	for (iLine=0; iLine<nLine; iLine++) {
		free(line[iLine].vert);
		free(line[iLine].sub[0].u);
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
}

