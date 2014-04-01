#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "def.h"

static double def_crude = 1.5;
static double lim_size = -1.0;

extern "C" int libaft_internal_surfmeshrefiner_ss_setup(double ss) {
	def_crude = ss;
	return 0;
}
extern "C" void libaft_internal_surfmeshrefiner_poly_setup_extralim(double lim);
extern "C" int libaft_internal_surfmeshrefiner_lim_setup(double lim) {
	libaft_internal_surfmeshrefiner_poly_setup_extralim(lim);
	lim_size = lim;
	return 0;
}

static double M_EPS = GL_EPS;

extern double TriangleValue(POINT vertices[3]);

extern "C" {
#include "libaft.h"
}

extern int *Rtri, *Rmaterial, RVmax, RTmax;
extern double *Rvrt;
extern int nver, ntr;
static double crude;
extern POINT *tmppoints;
extern FILE *fver, *ftri;
extern double thresholds[];

const int nvmax = 150000;       
const int ntmax = 2*nvmax;    
const int nbmax = 100000;
double vrt[2*nvmax], h[nvmax], crv[2*nvmax];
int iFNC[nvmax], material[ntmax], tri[3*ntmax], bnd[2*nbmax], bndc[nbmax];
int qual_mas_res[6]={0, 0, 0, 0, 0, 0};

void file_combiner (char *a, char *b)
{
	char *tmp;
	tmp = new char[200];
	FILE *f, *g;
	f = fopen(a, "a");
	g = fopen(b, "r");
	while(fgets(tmp, 200, g))
		fprintf(f, "%s", tmp);
	fclose(f);
	fclose(g);
	delete [] tmp;
}

inline int sign (double a)                         
{
	return(a>0 ? 1 : (a == 0 ? 0 : -1));
}

double** matrix_mult (double **a, double **b, double **c, int n) 
{
	int i, j, k;
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
		{
		   c[i][j] = 0;
		   for(k=0; k<n; k++)
			   c[i][j] += a[i][k]*b[k][j];
		}
	return(c);
}

double** trans_A (double **A, int n)           
{
	int i, j;
	double tmp;
	for(i=0; i<n; i++)
		for(j=n-1; j>i; j--)
		{
			tmp = A[j][i];
			A[j][i] = A[i][j];
			A[i][j] = tmp;
		}
	return(A);
}

void def_turn (double **a, double fi, int axe)  
{
	int i, j;
	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			if(i == j)
				i == axe ? a[i][j] = 1 : a[i][j] = cos(fi);
			else
				if(i == axe || j == axe)
					a[i][j] = 0;
				else
					i<j ? a[i][j] = sin(fi) : a[i][j] = -sin(fi);
		}
}

void transform_determing (TRIANGLE XYZ, POINT *points,double **A, double **A1, POINT& d)
{
	int i, j;
	POINT p = points[XYZ.vertex[1]-1]-points[XYZ.vertex[0]-1];
	POINT q = points[XYZ.vertex[2]-1]-points[XYZ.vertex[0]-1];
	double **tmpXYZ, ***turn, **res, **tmp, fi, psi, eta;
	tmpXYZ = new double*[3];
	res = new double*[3];
	turn = new double**[3];
	for(i=0; i<3; i++)
	{
		tmpXYZ[i] = new double[3];
		turn[i] = new double*[3];
		res[i] = new double[3];
		for(j=0; j<3; j++)
			turn[i][j] = new double[3];
	}
	d = -points[XYZ.vertex[0]-1];
	d.min_width = 100;
	for(i=0; i<3; i++)
		tmpXYZ[0][i] = 0;
	tmpXYZ[1][0] = p.x;
	tmpXYZ[1][1] = p.y;
	tmpXYZ[1][2] = p.z;
	tmpXYZ[2][0] = q.x;
	tmpXYZ[2][1] = q.y;
	tmpXYZ[2][2] = q.z;
	(tmpXYZ[1][1] == 0 && tmpXYZ[1][0] == 0) ? eta = 0 : eta = atan2(tmpXYZ[1][1], tmpXYZ[1][0]);
	def_turn(turn[2], eta, 2);
	matrix_mult(turn[2], trans_A(tmpXYZ, 3), res, 3);
	tmp = trans_A(res, 3);
	res = tmpXYZ;
	tmpXYZ = tmp;
	psi = atan2(tmpXYZ[1][2], tmpXYZ[1][0]);
	def_turn (turn[1], psi, 1);
	matrix_mult(turn[1], trans_A(tmpXYZ, 3), res, 3);
	tmp=trans_A(res, 3);
	res = tmpXYZ;
	tmpXYZ = tmp;
	fi = atan2(tmpXYZ[2][2], tmpXYZ[2][1]);
	def_turn(turn[0], fi, 0);
	matrix_mult(turn[0], matrix_mult(turn[1], turn[2], res, 3), A, 3);
	matrix_mult(trans_A(turn[2], 3), matrix_mult(trans_A(turn[1], 3), trans_A(turn[0], 3), res ,3), A1, 3);
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
			delete [] turn[i][j];
		delete [] turn[i];
		delete [] tmpXYZ[i];
		delete [] res[i];
	}
	delete [] turn;
	delete [] res;
	delete [] tmpXYZ;
}

POINT transform (double **A, const POINT d, const POINT x, int dir)
{
	int i, j;
	POINT t = x;
	double tmp[3], res[3];
	if(dir == 1)
		t = t+d;
	tmp[0] = t.x;
	tmp[1] = t.y;
	tmp[2] = t.z;
	for(i=0; i<3; i++)
	{
		res[i] = 0;
		for(j=0; j<3; j++)
			res[i] += A[i][j]*tmp[j];
	}
	t.x = res[0];
	t.y = res[1];
	t.z = res[2];
	if(dir == -1)
	t = t-d;
	return(t);
}

double getz (POINT X, POINT Y, POINT Z, POINT x)     
{
	POINT p = Y-X, q = Z-X, normal;
	double d;
	normal.x = p.y*q.z-p.z*q.y;
	normal.y = p.z*q.x-p.x*q.z;
	normal.z = p.x*q.y-p.y*q.x;
	d = normal*X;
	return ((d-normal*x)/normal.z);
}

void scale (POINT& pt, POINT shift, double stretch, int direction)
{
//	direction == 1 ? pt = pt/stretch-shift : pt = (pt+shift)*stretch;
}

double coef_detection (POINT *boundary, ITYPE boundary_size, POINT& shift)
{
	int i/*, j, k*/;
	double res;
	POINT min, max; 
	max = boundary[0];
	min = boundary[0];
	for(i=0; i<boundary_size; i++)
	{
		if(min.x>boundary[i].x)
			min.x = boundary[i].x;
		if(min.y>boundary[i].y)
			min.y = boundary[i].y;
		if(max.x<boundary[i].x)
			max.x = boundary[i].x;
		if(max.y<boundary[i].y)
			max.y = boundary[i].y;
	}
	(max.x-min.x-max.y+min.y)>0 ? res = fabs(max.x-min.x) : res = fabs(max.y-min.y);
	res /= 0.8;
	shift.x = min.x/res-0.1;
	shift.y = min.y/res-0.1;
	return(res);
}

int dot_vector_classify (POINT a, POINT b, POINT c)           
{
	POINT p = b-a, q = c-a, r = c-b;
	double s = p.x*q.y-p.y*q.x;
	int f = 0;
	p.z = 0;
	q.z = 0;
	r.z = 0;
	if(p*q<0 || (p*q)>(p*p))
		f = -2;
	if((f != -2) && fabs(s/(2*p.R()))<M_EPS)
		return(0);
	if(q.R()<M_EPS || r.R()<M_EPS)
		return(0);
	if(s<0)
		return(1);
	if(s>0)
		return(-1);
	return(f == -2 ? 2 : 0);
}

int dot_vector_classify2 (const POINT& a, const POINT& b, const POINT& c)           
{
	double pr1, d = sqrt((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y));

	pr1 = (c.x-a.x)*(b.y-a.y) - (c.y-a.y)*(b.x-a.x);
    
    if(pr1/d >= 1.e-16) 
        return 1;
    if(pr1 <= -1.e-16)
        return -1;
	if(sqrt((c.x-a.x)*(c.x-a.x)+(c.y-a.y)*(c.y-a.y))<1.e-16 ||sqrt((c.x-b.x)*(c.x-b.x)+(c.y-b.y)*(c.y-b.y))<1.e-16)
		return 0;
	return 2;
}

int vector_vector_classify (const POINT& a, const POINT& b, const POINT& c, const POINT& d)           
{
	int ca = dot_vector_classify2(c, d, a), cb = dot_vector_classify2(c, d, b);
	int cc = dot_vector_classify2(a, b, c), cd = dot_vector_classify2(a, b, d);
	if((ca*cb == -1 && cc*cd == -1) || abs(ca*cb) == 2 ||  abs(cc*cd) == 2)
		return (1);
	return (0);
}

int dot_tryangle_classify (const POINT a, const POINT b, const POINT c, const POINT x)  
{
	double flag[3];
	if((flag[0] = dot_vector_classify (a,b,x)) == 0)
		return (2);
	if((flag[1] = dot_vector_classify (b,c,x)) == 0)
		return (2);
	if((flag[2] = dot_vector_classify (c,a,x)) == 0)
		return (2);
	if(flag[2] == flag[1] && flag[2] == flag[0])
		return (1);
	return (0);
}

int orientation (POINT *pt, int lenght)
{
	int i, num = 0;
	double minx = pt[0].x, sign;
	POINT p, q;
	p.z = 0;
	q.z = 0;
	for(i=0; i<lenght-1; i++)
		if(pt[i].x<minx)
		{
			minx = pt[i].x;
			num = i;
		}
	p = pt[num+1];
	if(num == 0)
		q = pt[lenght-2];
	else
		q = pt[num-1];
	q = q-pt[num];
	p = p-pt[num];
	sign = q.x*p.y-q.y*p.x;
	if(sign>0 || (sign == 0 && p.y>q.y))
		return -1;
	else
		return 1;
}

void reorder(POINT** bounds, ITYPE* bd_size, ITYPE bound_num)
{
	int i, j, flag, int_tmp;
	bool *min = new bool[bound_num], *max = new bool[bound_num];
	POINT tmp, *tm;
	double minx = bounds[0][0].x, maxx = bounds[0][0].x;
	for(j=0; j<bound_num; j++)
		for(i=0; i<bd_size[j]; i++)
		{
			if(minx>bounds[j][i].x)
				minx = bounds[j][i].x;
			if(maxx<bounds[j][i].x)
				maxx = bounds[j][i].x;
		}
	for(j=0; j<bound_num; j++)
	{
		min[j] = 0;
		max[j] = 0;
		for(i=0; i<bd_size[j]; i++)
		{
			if(fabs(minx-bounds[j][i].x)<GL_EPS)
				min[j] = 1;
			if(fabs(maxx-bounds[j][i].x)<GL_EPS)
				max[j] = 1;
		}
	}
	for(j=0; j<bound_num; j++)
		if(max[j] == 1 && min[j] == 1)
		{
			tm = bounds[0];
			bounds[0] = bounds[j];
			bounds[j] = tm;
			int_tmp = bd_size[0];
			bd_size[0] = bd_size[j];
			bd_size[j] = int_tmp;
		}
	for(j=0; j<bound_num; j++)
	{
		flag = orientation(bounds[j], bd_size[j]);
		if((flag == 1 && j == 0)||(flag == -1 && j != 0))
			for(i=0; i<(int)bd_size[j]/2; i++)
			{
				tmp = bounds[j][bd_size[j]-1-i];
				bounds[j][bd_size[j]-1-i] = bounds[j][i];
				bounds[j][i] = tmp;
			}
	}
	delete [] min;
	delete [] max;
}

int MeshPolyreg3D(MESH_DATA& mesh, POLYREGION& polyreg, int bound_num)
{
	double stretch, check_res/*, try_num=0*/;
	POINT shift, **bounds, *pr, tmp, check[3];
	TRIANGLE index;
	int i, j, flag, error, repeat, retry = 1;
	int nv = 0, nb = 0, nt = 0/*, nc = 0*/, *bd_size, beg;
	int nnv = nvmax, nnb = nbmax, nnt = ntmax;
	shift.z = 0;
	bounds = new POINT*[bound_num];
	bd_size = new int[bound_num];
	tmp = polyreg.boundary[0];
	beg = 0;
	j = 0;
	for(i = 2; i<polyreg.boundary_size; i++)
	{
		if((tmp-polyreg.boundary[i]).R()<GL_EPS)
		{
			bd_size[j] = i-beg+1;
			//printf("%d\n",bd_size[j]);
			bounds[j] = new POINT[bd_size[j]];
			memcpy(bounds[j], polyreg.boundary+beg, bd_size[j]*sizeof(POINT));
			j++;
			beg = i+1;
			i++;
			if(beg<polyreg.boundary_size)
				tmp = polyreg.boundary[beg];
		}
	}	
	//printf("%d\n",bound_num);
	for(j=0; j<bound_num; j++)
		for(i=0; i<bd_size[j]; i++)
			bounds[j][i] = transform(polyreg.m_rotate, polyreg.m_move, bounds[j][i], 1);

	/*FILE *o;
				//o = fopen("bad_border", "w");
				//for(i=0;i<bd_size[0];i++)
				//	fprintf(o,"%lf %lf\n",bounds[0][i].x,bounds[0][i].y);
				//fclose(o);
				o = fopen("bor", "w");
				fprintf(o,"%d\n",polyreg.boundary_size-bound_num);
				for(j=0;j<bound_num;j++)
					for(i=0;i<bd_size[j]-1;i++)
					{
						fprintf(o,"%lf %lf ",bounds[j][i].x,bounds[j][i].y);
						fprintf(o,"%lf %lf\n",bounds[j][i+1].x,bounds[j][i+1].y);   
					}
				fclose(o);
	
	*/
	if(bound_num>1)
		reorder(bounds, bd_size, bound_num);
	else
	{
		if(orientation(bounds[0], bd_size[0]) == 1)
			for(i=0; i<(int)bd_size[0]/2; i++)
			{
				tmp = bounds[0][bd_size[0]-1-i];
				bounds[0][bd_size[0]-1-i] = bounds[0][i];
				bounds[0][i] = tmp;
			}
	}
	stretch = coef_detection(bounds[0], bd_size[0], shift);
	for(j=0; j<bound_num; j++)
		for(i=0; i<bd_size[j]; i++)
			scale(bounds[j][i], shift,stretch, 1);
	nb = 0; nv = 0;
	for(j=0; j<bound_num; j++) {
		for(i=0; i<bd_size[j]; i++)
		{
			vrt[2*nv+0] = bounds[j][i].x;
			vrt[2*nv+1] = bounds[j][i].y;
			nv++;
		}
	}
//	nv = polyreg.boundary_size;
	crude = def_crude;
	while(retry)
	{
		retry = 0;
		error = 1;
		while (crude<2.51 && error)
		{
			/*error = 0;*/
			nv = polyreg.boundary_size;
			nb = 0;
			nt = 0;
			error = mesh_2d_aft_cf_lim_loop_check_(&nv, vrt, &nb, bnd, bndc, &nt, tri, material, &nnv, &nnb, &nnt, &crude, &lim_size);
			if (error) {
/*					FILE *o;
					char buf[1024];
					o = fopen("bad_border", "w");
					for(i=0;i<bd_size[0];i++)
						fprintf(o,"%lf %lf\n",bounds[0][i].x,bounds[0][i].y);
					fclose(o);
					o = fopen("bor", "w");
					fprintf(o,"%d\n",polyreg.boundary_size-bound_num);
					for(j=0;j<bound_num;j++)
						for(i=0;i<bd_size[j]-1;i++)
						{
							fprintf(o,"%lf %lf ",bounds[j][i].x,bounds[j][i].y);
							fprintf(o,"%lf %lf\n",bounds[j][i+1].x,bounds[j][i+1].y);   
						}
					fclose(o);
					scanf("%s", buf);
				};
			}
			catch(...)
			{*/
				/*error = 1;*/
				if (error<0) {
					crude += 0.1;
					printf("WARNING!!! crude had been increased to %lf\n", crude);
				}

				/*FILE *o;
				o = fopen("bad_border", "w");
				for(i=0;i<bd_size[0];i++)
					fprintf(o,"%lf %lf\n",bounds[0][i].x,bounds[0][i].y);
				fclose(o);
				o = fopen("bor", "w");
				fprintf(o,"%d\n",polyreg.boundary_size-bound_num);
				for(j=0;j<bound_num;j++)
					for(i=0;i<bd_size[j]-1;i++)
					{
						fprintf(o,"%lf %lf ",bounds[j][i].x,bounds[j][i].y);
						fprintf(o,"%lf %lf\n",bounds[j][i+1].x,bounds[j][i+1].y);   
					}
				fclose(o);*/

				if ((crude>2.5) || (error>0))
				{	
					for(j=0; j<bound_num; j++)
						delete [] bounds[j];
					delete [] bounds;
					delete [] bd_size;
					printf("ERROR!!! rebuilding polyregion now\n");
					return (0);	
				}
			}
		}
		pr = new POINT[nv];
		for(i=0; i<nv && !retry; i++)
		{
			pr[i].x = vrt[2*i];
			pr[i].y = vrt[2*i+1];
			pr[i].z = 0;
			scale(pr[i], shift, stretch, -1);
			M_EPS = GL_EPS;
			flag = 0;
			repeat = 1;
			while(repeat && !retry)
			{
				j = 0;
				while((j<polyreg.triangle_size) && (flag == 0) && !retry)
				{
					index = mesh.triangles[polyreg.triangles[j]];
					if(dot_tryangle_classify(tmppoints[index.vertex[0]-1], tmppoints[index.vertex[1]-1], tmppoints[index.vertex[2]-1], pr[i]) == 0)
   						j++;
					else
					{
						flag = 1;
						repeat = 0;
						pr[i].z = getz(tmppoints[index.vertex[0]-1], tmppoints[index.vertex[1]-1], tmppoints[index.vertex[2]-1], pr[i]);
						pr[i] = transform(polyreg.m_rotate_back, polyreg.m_move, pr[i], -1);
					}
				} 
				if(j == polyreg.triangle_size)
				{
					M_EPS *= 10;			
					printf("WARNING!!! M_EPS had been increased to %e\n", M_EPS);
					if(M_EPS/GL_EPS > 1.e6)
					{
						crude+=0.1;
						retry = 1;
						delete [] pr;
						pr = NULL;
						printf("WARNING 2!!! crude had been increased to %lf\n", crude);
						if (crude>2.51)
						{
							printf("ERROR2!!! rebuilding polyregion now\n");
							for(j=0; j<bound_num; j++)
								delete [] bounds[j];
							delete [] bounds;
							delete [] bd_size;
							return (0);
						}
					}
				}
			}
		}
	}

	/*ftri = fopen("rest.txt", "w");
	fprintf(ftri, "  %d  %d\n", nv, nt);
	for(i=0; i<nv; i++)
		fprintf(ftri, "  %lf  %lf  %lf\n", pr[i].x, pr[i].y , pr[i].z);*/
	if(nver+nv > RVmax)
    {
        printf("It is necessary to increase Vmax parameter!!!");
        exit(1);
    }
	if(ntr+nt > RTmax)
    {
        printf("It is necessary to increase Tmax parameter!!!");
        exit(1);
    }
	for(i=0; i<nv; i++)
	{
		Rvrt[3*(nver+i)] = pr[i].x;
		Rvrt[3*(nver+i)+1] = pr[i].y; 
		Rvrt[3*(nver+i)+2] = pr[i].z;
	}
	for(j=0; j<nt; j++)
	{
		for(int k=0; k<3; k++)
		{
			check[k] = pr[tri[j*3+k]-1];
			Rtri[3*(ntr+j)+k] = tri[j*3+k]+nver;
			//fprintf(ftri, " %d", tri[j*3+k]);
		}
		check_res = TriangleValue(check);
		for(int k=0; k<6; k++)
            if(check_res < thresholds[k])
            {
                qual_mas_res[k]++;
                break;
            }
		Rmaterial[ntr+j] = mesh.triangles[polyreg.triangles[0]].color;
		//fprintf(ftri, " %d\n", mesh.triangles[polyreg.triangles[0]].color);
	}
	//fclose(ftri);
	nver+=nv;
	ntr+=nt;

	delete [] pr;
	for(j=0; j<bound_num; j++)
		delete [] bounds[j];
	delete [] bounds;
	delete [] bd_size;
	return (1);
}
