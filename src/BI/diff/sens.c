#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double s;
    int m;
    double v;
    double x[3];
    int id;
} cell;

typedef struct {
    int nc, nnc;
    cell *c;
    double st;
    double *s;
} bar;

typedef struct {
    int nc, nnc;
    bar p, m;
    double *s;
} smesh;

static smesh mesh = {0, 0};
static smesh *pm = &mesh;

static int sensadd(bar *pb, int id, double s, int m, double v, double x0, double x1, double x2){
    if (pb->nc+1 >= pb->nnc) {
	pb->nnc += 65536;
	pb->c = realloc(pb->c, sizeof(cell)*pb->nnc);
	pb->s = realloc(pb->s, sizeof(double)*pb->nnc);
	if (!pb->c)  return  perror("realloc"), 2;
	if (!pb->s)  return  perror("realloc"), 2;
    }
    pb->c[pb->nc].s = s,  pb->c[pb->nc].m = m,  pb->c[pb->nc].v = v,  pb->c[pb->nc].x[0] = x0,  pb->c[pb->nc].x[1] = x1,  pb->c[pb->nc].x[2] = x2;
    pb->c[pb->nc++].id = id;
    return 0;
}

int addsens_(double *ps, int *pma, double *pv, double *px, double *py, double *pz) {
    double s = *ps, v = *pv, x0 = *px, x1 = *py, x2 = *pz;
    int m = *pma;

    if (pm->nc+1 >= pm->nnc) {
	pm->nnc += 65536;
	pm->s = realloc(pm->s, sizeof(double)*pm->nnc);
	if (!pm->s)  return  perror("realloc"), 2;
    }
    if (s > 0.0)  sensadd(&pm->p, pm->nc, s, m, v, x0, x1, x2);
    else if (s < 0.0)  sensadd(&pm->m, pm->nc, -s, m, v, x0, x1, x2);
    else  pm->s[pm->nc] = 0.0;
    pm->nc++;
    return 0;
}

static int cellsort(const void *p1, const void *p2) {
    cell *c1 = (cell*)p1,  *c2 = (cell*)p2;

    if (c1->s > c2->s)  return -1;
    else if (c1->s < c2->s)  return 1;
    else return 0;
}

static double summ(bar *pb) {
    int i;
    double s;
    for (i=0; i<pb->nc; i++) {
	pb->s[i] = pb->st;
	s = pb->c[i].s * pb->c[i].v;
	pb->st += s;
    }
    pb->s[pb->nc] = pb->st;
    return pb->st;
}

static int remark(bar *pb, double *s, double d) {
    int i;
    for (i=0; i<pb->nc; i++) {
	s[pb->c[i].id] = d * (pb->st - pb->s[i])/pb->st;
    }
    return 0;
}

static double ffind(bar *pb, double x) {
    int a = 0, b = pb->nc, c;
    double r;

    while (b>a+1) {
	c = (a+b)/2;
	if (pb->c[c].s < x)  b = c;
	else a = c;
    }
    c = (fabs(x - pb->c[a].s) < fabs(x - pb->c[a+1].s)) ? a : a+1;
    r = (pb->st - pb->s[c])/pb->st * 100.0;
    if (r>100.0)  r = 100.0;
    return r;
}

int sensremap_(int *pn, double *sns, double *snm) {
    int i, n = *pn;
    double x, y;

    qsort(mesh.p.c, mesh.p.nc, sizeof(cell), cellsort);
    qsort(mesh.m.c, mesh.m.nc, sizeof(cell), cellsort);

    summ(&mesh.p);
    summ(&mesh.m);

    remark(&mesh.p, mesh.s, 100.0);
    remark(&mesh.m, mesh.s, -100.0);

    for (i=0; i<n; i++) {
	x = sns[i];
	if (x > 0.0)  y = ffind(&pm->p, x);
	else if (x < 0.0)  y = -ffind(&pm->m, -x);
	else  y = 0.0;
	snm[i] = y;
    }
    return 0;
}

