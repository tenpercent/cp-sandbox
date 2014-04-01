#include <stdio.h>
#include "aft.h"
#include "helper.h"
#include "det.h"

//#define SHOWPROGRESS

int refine2daft(int *pnV, REAL *vertex, int *pnVF, int *pnT, int *tria) {
    int nV = *pnV, nT = *pnT, nVF = *pnVF;
    int *vi, *vj;
    int i, j, k, m, v, t, good, vm;
    REAL step, x, y;
    REAL *dv;

    dv = (REAL*)libaft_malloc(sizeof(REAL)*nV*2);
    vi = (int*)libaft_malloc(sizeof(int)*(nV+1));
    vj = (int*)libaft_malloc(sizeof(int)*3*nT);
    for (i=0; i<=nV; i++) vi[i] = 0;
    for (j=0; j<3*nT; j++) {
	vi[tria[j]+1]++;
    }
    for (i=1; i<nV; i++) vi[i] += vi[i-1];
    for (j=0; j<3*nT; j++) {
	vj[vi[tria[j]]++] = j/3;
    }
    for (i=nV; i>0; i--) vi[i] = vi[i-1];
    vi[0] = 0;

    for (k=0; k<10; k++) {
	step = 1.0/(k+2.0);
	for (i=0; i<nV; i++) {
	    dv[2*i+0] = vertex[2*i+0];
	    dv[2*i+1] = vertex[2*i+1];
	}
	for (i=nVF; i<nV; i++) {
	    m = 0;
	    x = 0.0, y = 0.0;
	    for (j=vi[i]; j<vi[i+1]; j++) {
		for (t=0; t<3; t++) {
		    v = tria[3*vj[j]+t];
		    if (v==i) continue;
		    x += dv[2*v+0];
		    y += dv[2*v+1];
		    m++;
		}
	    }
	    if (m) {
		x /= m, y /= m;
		vertex[2*i+0] = x*step + dv[2*i+0]*(1.0-step);
		vertex[2*i+1] = y*step + dv[2*i+1]*(1.0-step);
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, moving %6.2lf%%            \r", k+1, 100.0*(1.0+i)/nV);
	    fflush(stdout);
#endif
	}
	good = 1;
	for (j=0; j<nT; j++) {
	    if (idet2i3(vertex, tria[3*j+0], tria[3*j+1], tria[3*j+2])<1) {
		good = 0;
		break;
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, checking %6.2lf%%          \r", k+1, 100.0*(1.0+j)/nT);
	    fflush(stdout);
#endif
	}
	if (!good) {
	    for (i=0; i<nV; i++) {
		vertex[2*i+0] = dv[2*i+0];
		vertex[2*i+1] = dv[2*i+1];
	    }
	    break;
	}
    }
    if (!good) {
	for (k=0; k<10; k++) {
	    vm = 0;
	    step = 1.0/(k+2.0);
	    for (i=nVF; i<nV; i++) {
		m = 0;
		x = 0.0, y = 0.0;
		for (j=vi[i]; j<vi[i+1]; j++) {
		    for (t=0; t<3; t++) {
			v = tria[3*vj[j]+t];
			if (v==i) continue;
			x += dv[2*v+0];
			y += dv[2*v+1];
			m++;
		    }
		}
		if (m) {
		    dv[2*i+0] = vertex[2*i+0];
		    dv[2*i+1] = vertex[2*i+1];
		    x /= m, y /= m;
		    vertex[2*i+0] = x*step + dv[2*i+0]*(1.0-step);
		    vertex[2*i+1] = y*step + dv[2*i+1]*(1.0-step);
		    good = 1;
		    for (t=vi[i]; t<vi[i+1]; t++) {
			j = vj[t];
			if (idet2i3(vertex, tria[3*j+0], tria[3*j+1], tria[3*j+2])<1) {
			    good = 0;
			    break;
			}
		    }
		    if (!good) {
			vertex[2*i+0] = dv[2*i+0];
			vertex[2*i+1] = dv[2*i+1];
		    } else {
			vm++;
		    }
		}
#ifdef SHOWPROGRESS
		printf("S2: %2d, moving %6.2lf%%            \r", k+1, 100.0*(1.0+i)/nV);
		fflush(stdout);
#endif
	    }
	}
    }
#ifdef SHOWPROGRESS
    printf("                                   \r");
    fflush(stdout);
#endif
    libaft_free(dv),  libaft_free(vi),  libaft_free(vj);
    return 0;
}

static int add_glist(int a, int b, int *png, int *glist, int *s) {
    int c = s[a];
    while (c>=0) {
        if (glist[2*c+0] == b) return 0;
        c = glist[2*c+1];
    }
    glist[2**png+0] = b;
    glist[2**png+1] = s[a];
    s[a] = *png;
    (*png)++;
    return 1;
}



int refine3daft_new(int *pnV, REAL *vertex, int *pnVF, int *pnT, int *tetra) {
    int nV = *pnV, nVF = *pnVF, nT = *pnT;
    int *s, ng, *glist, j, i, a, b, c, d, n, *ia, *ja, l;
    int k, m, v, t, good;
    REAL step, x, y, z;
    REAL *dv;

    glist = (int*)libaft_malloc(2*sizeof(int)*12*nT);
    ng = 0;
    s = (int*)libaft_malloc(sizeof(int)*nV);
    for (j=0; j<nV; j++)  s[j] = -1;
    for (i=0; i<nT; i++) {
	a = tetra[4*i + 0];
	b = tetra[4*i + 1];
	c = tetra[4*i + 2];
	d = tetra[4*i + 3];
	add_glist(a, b, &ng, glist, s);
	add_glist(a, c, &ng, glist, s);
	add_glist(a, d, &ng, glist, s);
	add_glist(b, a, &ng, glist, s);
	add_glist(b, c, &ng, glist, s);
	add_glist(b, d, &ng, glist, s);
	add_glist(c, a, &ng, glist, s);
	add_glist(c, b, &ng, glist, s);
	add_glist(c, d, &ng, glist, s);
	add_glist(d, a, &ng, glist, s);
	add_glist(d, b, &ng, glist, s);
	add_glist(d, c, &ng, glist, s);
    }
    ia = (int*)libaft_malloc(sizeof(int)*(nV+1));
    ja = (int*)libaft_malloc(sizeof(int)*(12*nT));
    n = 0;
    for (j=0; j<nV; j++) {
	ia[j] = n;
	l = s[j];
	while (l>=0) {
	    ja[n++] = glist[2*l+0];
	    l = glist[2*l+1];
	}
    }
    ia[nV] = n;

    dv = (REAL*)libaft_malloc(sizeof(REAL)*nV*3);

    for (k=0; k<10; k++) {
	step = 1.0/20.0;
	for (i=0; i<nV; i++) {
	    dv[3*i+0] = vertex[3*i+0];
	    dv[3*i+1] = vertex[3*i+1];
	    dv[3*i+2] = vertex[3*i+2];
	}
	for (i=nVF; i<nV; i++) {
	    m = 0;
	    x = 0.0, y = 0.0, z = 0.0;
	    for (j=ia[i]; j<ia[i+1]; j++) {
		for (t=0; t<4; t++) {
		    v = tetra[4*ja[j]+t];
		    if (v==i) continue;
		    x += dv[3*v+0];
		    y += dv[3*v+1];
		    z += dv[3*v+2];
		    m++;
		}
	    }
	    if (m) {
		x /= m, y /= m, z /=m;
		vertex[3*i+0] = x*step + dv[3*i+0]*(1.0-step);
		vertex[3*i+1] = y*step + dv[3*i+1]*(1.0-step);
		vertex[3*i+2] = z*step + dv[3*i+2]*(1.0-step);
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, moving %6.2lf%%            \r", k+1, 100.0*(1.0+i)/nV);
	    fflush(stdout);
#endif
	}
	good = 1;
	for (j=0; j<nT; j++) {
	    if (idet3i4(vertex, tetra[4*j+0], tetra[4*j+1], tetra[4*j+2], tetra[4*j+3])<1) {
		good = 0;
		break;
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, checking %6.2lf%%          \r", k+1, 100.0*(1.0+j)/nT);
	    fflush(stdout);
#endif
	}
	if (!good) {
	    for (i=nVF; i<nV; i++) {
		vertex[3*i+0] = dv[3*i+0];
		vertex[3*i+1] = dv[3*i+1];
		vertex[3*i+2] = dv[3*i+2];
	    }
	    break;
	}
    }
    libaft_free(dv);
    libaft_free(ia),  libaft_free(ja);
    libaft_free(glist),  libaft_free(s);
    return 0;
}

int refine3daft(int *pnV, REAL *vertex, int *pnVF, int *pnT, int *tetra) {
    int nV = *pnV, nVF = *pnVF, nT = *pnT;
    int *vi, *vj;
    int i, j, k, m, v, t, good, vm;
    REAL step, x, y, z;
    REAL *dv;

    dv = (REAL*)libaft_malloc(sizeof(REAL)*nV*3);
    vi = (int*)libaft_malloc(sizeof(int)*(nV+1));
    vj = (int*)libaft_malloc(sizeof(int)*4*nT);
    for (i=0; i<=nV; i++) vi[i] = 0;
    for (j=0; j<4*nT; j++) {
	vi[tetra[j]+1]++;
    }
    for (i=1; i<nV; i++) vi[i] += vi[i-1];
    for (j=0; j<4*nT; j++) {
	vj[vi[tetra[j]]++] = j/4;
    }
    for (i=nV; i>0; i--) vi[i] = vi[i-1];
    vi[0] = 0;

    for (k=0; k<10; k++) {
//	step = 1.0/20.0;
	step = 1.0/(k+2.0);
	for (i=0; i<nV; i++) {
	    dv[3*i+0] = vertex[3*i+0];
	    dv[3*i+1] = vertex[3*i+1];
	    dv[3*i+2] = vertex[3*i+2];
	}
	for (i=nVF; i<nV; i++) {
	    m = 0;
	    x = 0.0, y = 0.0, z = 0.0;
	    for (j=vi[i]; j<vi[i+1]; j++) {
		for (t=0; t<4; t++) {
		    v = tetra[4*vj[j]+t];
		    if (v==i) continue;
		    x += dv[3*v+0];
		    y += dv[3*v+1];
		    z += dv[3*v+2];
		    m++;
		}
	    }
	    if (m) {
		x /= m, y /= m, z /=m;
		vertex[3*i+0] = x*step + dv[3*i+0]*(1.0-step);
		vertex[3*i+1] = y*step + dv[3*i+1]*(1.0-step);
		vertex[3*i+2] = z*step + dv[3*i+2]*(1.0-step);
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, moving %6.2lf%%            \r", k+1, 100.0*(1.0+i)/nV);
	    fflush(stdout);
#endif
	}
	good = 1;
	for (j=0; j<nT; j++) {
	    if (idet3i4(vertex, tetra[4*j+0], tetra[4*j+1], tetra[4*j+2], tetra[4*j+3])<1) {
		good = 0;
		break;
	    }
#ifdef SHOWPROGRESS
	    printf("S1: %2d, checking %6.2lf%%          \r", k+1, 100.0*(1.0+j)/nT);
	    fflush(stdout);
#endif
	}
	if (!good) {
//	    printf("\nS1: failed! (v = %le, v_o = %le)\n", det3i4(vertex, tetra[4*j+0], tetra[4*j+1], tetra[4*j+2], tetra[4*j+3]),
//		    det3i4(dv, tetra[4*j+0], tetra[4*j+1], tetra[4*j+2], tetra[4*j+3]));
	    for (i=0; i<nV; i++) {
		vertex[3*i+0] = dv[3*i+0];
		vertex[3*i+1] = dv[3*i+1];
		vertex[3*i+2] = dv[3*i+2];
	    }
	    break;
	}
    }
    if (!good) {
	for (k=0; k<10; k++) {
	    vm = 0;
//	    step = 1.0/20.0;
	    step = 1.0/(k+2.0);
	    for (i=nVF; i<nV; i++) {
		m = 0;
		x = 0.0, y = 0.0, z = 0.0;
		for (j=vi[i]; j<vi[i+1]; j++) {
		    for (t=0; t<4; t++) {
			v = tetra[4*vj[j]+t];
			if (v==i) continue;
			x += dv[3*v+0];
			y += dv[3*v+1];
			z += dv[3*v+2];
			m++;
		    }
		}
		if (m) {
		    dv[3*i+0] = vertex[3*i+0];
		    dv[3*i+1] = vertex[3*i+1];
		    dv[3*i+2] = vertex[3*i+2];
		    x /= m, y /= m, z /=m;
		    vertex[3*i+0] = x*step + dv[3*i+0]*(1.0-step);
		    vertex[3*i+1] = y*step + dv[3*i+1]*(1.0-step);
		    vertex[3*i+2] = z*step + dv[3*i+2]*(1.0-step);
		    good = 1;
		    for (t=vi[i]; t<vi[i+1]; t++) {
			j = vj[t];
			if (idet3i4(vertex, tetra[4*j+0], tetra[4*j+1], tetra[4*j+2], tetra[4*j+3])<1) {
			    good = 0;
			    break;
			}
		    }
		    if (!good) {
			vertex[3*i+0] = dv[3*i+0];
			vertex[3*i+1] = dv[3*i+1];
			vertex[3*i+2] = dv[3*i+2];
		    } else {
			vm++;
		    }
		}
#ifdef SHOWPROGRESS
		printf("S2: %2d, moving %6.2lf%% [% 4d]     \r", k+1, 100.0*(1.0+i)/nV, vm);
		fflush(stdout);
#endif
	    }
//	    printf("\n");
	}
    }
#ifdef SHOWPROGRESS
    printf("                                   \r");
    fflush(stdout);
#endif

    libaft_free(dv),  libaft_free(vi),  libaft_free(vj);
    return 0;
}
