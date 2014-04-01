#include "common.h"
#include "memory3.h"
#include "tree3.h"
#include "user3.h"
#include "error3.h"
#include "region3.h"
#include "refine3.h"

static double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/

static void calcNeigTetra2(surface_mesh *pm) {
	int i, j, v1;

	for (i=0; i<pm->nPoint; i++)
		pm->neigTetra[i].n = 0;
	for (i=0; i<pm->nTria; i++) {
		v1 = pm->v1[i];
		j = pm->neigTetra[v1].n;
		if (j >= MAX_NEIGBOR)
			errorExit3(3, "j >= NAX_NEIGBOR");
		pm->neigTetra[v1].n++;
      pm->neigTetra[v1].neig[j] = i;

      v1 = pm->v2[i];
      j = pm->neigTetra[v1].n;
      if (j >= MAX_NEIGBOR)
        errorExit3(3,"j >= NAX_NEIGBOR");
      pm->neigTetra[v1].n++;
      pm->neigTetra[v1].neig[j] = i;

      v1 = pm->v3[i];
      j = pm->neigTetra[v1].n;
      if (j >= MAX_NEIGBOR)
        errorExit3(3,"j >= NAX_NEIGBOR");
      pm->neigTetra[v1].n++;
      pm->neigTetra[v1].neig[j] = i;
   }
   return;
} /*calcNeigTetra2*/


static void calcNeigbor2(surface_mesh *pm) {
	int i, j, k, n, iTria, vert[3*MAX_NEIGBOR];

	for (i=0; i<pm->nPoint; i++) {
		if (i < pm->nLinePoint)
			pm->neigbor[i].n = -1;
		else
			pm->neigbor[i].n = 0;
	}
	for (i=0; i<pm->nPoint; i++) {
		n = 0;
		for (j=0; j<pm->neigTetra[i].n; j++) {
			iTria = pm->neigTetra[i].neig[j];
			vert[3*j+0] = pm->v1[iTria];
			vert[3*j+1] = pm->v2[iTria];
			vert[3*j+2] = pm->v3[iTria];
		}
		for (j=0; j<3*pm->neigTetra[i].n; j++) {
			if (vert[j] == i)
				continue;
			iTria = 0;
			for (k=0; k<n; k++)
				if (vert[j] == pm->neigbor[i].neig[k]) {
					iTria = 1;
					break;
				}
			if (iTria) continue;
			pm->neigbor[i].neig[n] = vert[j];
			n++;
		}
		if (pm->neigbor[i].n < 0)
			pm->neigbor[i].n = -n;
		else
			pm->neigbor[i].n = n;
	}
   return;
} /*calcNeigbor2*/


void smoothingSurf(surface_mesh *pm) {
	int    i, j, k, n, nn;
	double x0, y0, z0, xx, yy, zz, x, y, z, u, v;
	double f1, du, dv, du1, dv1, du2, dv2, maxdu, maxdv, d1u, d1v, size;

	/* !!! alloc only for  2*pm->nPoint */
	pm->neigTetra = (PStrucNeigbor)myAlloc(2*pm->nPoint*S_StrucNeigbor);   // cbD
	pm->neigbor   = (PStrucNeigbor)myAlloc(2*pm->nPoint*S_StrucNeigbor);   // cbD
	calcNeigTetra2(pm);
	calcNeigbor2(pm);

	for (i=0; i<pm->nPoint; i++) {
		n = pm->neigbor[i].n;
		if (n > 0) {
			xx = yy = zz = 0.0;
			for (j=0; j<n; j++) {
				k = pm->neigbor[i].neig[j];
				nn = pm->neigbor[k].n;
				x = pm->vert[k].x;
				y = pm->vert[k].y;
				z = pm->vert[k].z;
				xx += x;  yy += y;  zz += z;
			}
			xx /= n;  yy /= n;  zz /= n;

			u = pm->vert[i].u;
			v = pm->vert[i].v;
			bounSurf(pm, pm->iSurf, u, v, &x, &y, &z);
			size = sizeFace(pm, x, y, z);
			f1 = distance(xx, yy, zz, x, y, z);
			for (j=0; j<15; j++) { /* NEW METH */
				du = 0.2*(pm->uMax - pm->uMin);
				do {
					du /= 2.0;
					bounSurf(pm, pm->iSurf, u+du, v, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.01*size);
				d1u = (distance(xx, yy, zz, x0, y0, z0) - f1)/du;

				dv = 0.2*(pm->vMax - pm->vMin);
				do {
					dv /= 2.0;
					bounSurf(pm, pm->iSurf, u, v+dv, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.01*size);
				d1v = (distance(xx, yy, zz, x0, y0, z0) - f1)/dv;

				du = 0.6*(pm->uMax - pm->uMin);
				do {
					du /= 2.0;
					bounSurf(pm, pm->iSurf, u+du, v, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

				dv = 0.6*(pm->vMax - pm->vMin);
				do {
					dv /= 2.0;
					bounSurf(pm, pm->iSurf, u, v+dv, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

				maxdu = du;
				maxdv = dv;
				if (d1u == 0. && d1v == 0)
					errorExit3(3, " der == 0.   in  smooth ");
				if (fabs(d1u) > fabs(d1v)) {
					dv1 = maxdv;
					du1 = (-f1 - d1v*dv)/d1u;
					dv2 = -maxdv;
					du2 = (-f1 - d1v*dv)/d1u;
					if (fabs(du1) < fabs(du2)) {
						du = du1;
						dv = dv1;
					} else {
						du = du2;
						dv = dv2;
					}
					if (du > maxdu) {
						dv *= (maxdu/du);
						du = maxdu;
					}
					if (du < -maxdu) {
						dv *= (-maxdu/du);
						du = -maxdu;
					}
				} else {
					du1 = maxdu;
					dv1 = (-f1 - d1u*du)/d1v;
					du2 = maxdu;
					dv2 = (-f1 - d1u*du)/d1v;
					if (fabs(dv1) < fabs(dv2)) {
						du = du1;
						dv = dv1;
					} else {
						du = du2;
						dv = dv2;
					}
					if (dv > maxdv) {
						du *= (maxdv/dv);
						dv = maxdv;
					}
					if (dv < -maxdv) {
						du *= (-maxdv/dv);
						dv = -maxdv;
					}
				}
				u += du;
				v += dv;
				bounSurf(pm, pm->iSurf, u, v, &x, &y, &z);
				f1 = distance(xx, yy, zz, x, y, z);
			}/*NEW METH*/
			pm->vert[i].u = u;
			pm->vert[i].v = v;
			pm->vert[i].x = x;
			pm->vert[i].y = y;
			pm->vert[i].z = z;
		}
	}/* for i */

	free(pm->neigbor);
	free(pm->neigTetra);

	return;
} /*smoothingSurf*/


