#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

int savebin_(int *pd, char *labels, ...) {
    FILE *f;
    va_list ap;
    int nc = *pd, i, n1, n2, *data;
    char *fmt = strdup(labels), *pc, *pc2;
    char fname[16384];

    va_start(ap, labels);

    pc = strchr(fmt, ':');
    *pc = 0;

    for (i=0; i<nc; i++) {
	n1 = ((int*)va_arg(ap, int*))[0];
	n2 = ((int*)va_arg(ap, int*))[0];
	data = (int*)va_arg(ap, int*);
	pc++;
	pc2 = strchr(pc, '.');
	*pc2 = 0;
	sprintf(fname, "%s.%s", fmt, pc);
	if (!(f = fopen(fname, "w")))  perror(fname);
	fwrite(data, sizeof(double), n1*n2, f);
	fclose(f);
	if (0) { // what?
	    fprintf(stderr, "	<DataItem Dimensions=\"");
	    fprintf(stderr, (n2>1)?"%d %d":"%d", n1, n2);
	    fprintf(stderr, "\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\">%s.%s</DataItem>\n", fmt, pc);
	}
	pc = pc2;
    }
    free(fmt);
    return 0;
}

int loadbin_(int *pnnv, int *pnnf, int *pnnt, int *pnv, int *pnf, int *pnt, double *vertex, int *tria, int *tetra, int *labelf, int *labelt, char *name) {
    FILE *f;
    char *bname = strdup(name);
    char fname[16384];
    int nv, nf, nt, i;

    strchr(bname, '*')[0] = 0;

    fprintf(stderr, "name: \"%s\"\n", bname);

    sprintf(fname, "%s.txt", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fscanf(f, "%d %d %d\n", &nv, &nt, &nf);
    fclose(f);
    if (nv>=*pnnv)  fprintf(stderr, "maxnv\n");
    if (nf>=*pnnf)  fprintf(stderr, "maxnf\n");
    if (nt>=*pnnt)  fprintf(stderr, "maxnt\n");

    sprintf(fname, "%s.crd", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fread(vertex, sizeof(double), 3*nv, f);
    fclose(f);

    sprintf(fname, "%s.tet", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fread(tetra, sizeof(int), 4*nt, f);
    fclose(f);

    sprintf(fname, "%s.lbt", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fread(labelt, sizeof(int), nt, f);
    fclose(f);

    sprintf(fname, "%s.tri", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fread(tria, sizeof(int), 4*nf, f);
    fclose(f);

    sprintf(fname, "%s.lbf", bname);
    if (!(f = fopen(fname, "r")))  perror(fname);
    fread(labelf, sizeof(int), nf, f);
    fclose(f);

    for (i=0; i<4*nt; i++)  tetra[i]++;
    for (i=0; i<3*nf; i++)  tria[i]++;

    *pnv = nv,  *pnf = nf,  *pnt = nt;

    free(bname);

    return 0;
}
