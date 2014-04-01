#include "def.h"

//extern "C" {
//#include "helper.h"
//}

POINT *tmppoints;
int *Rtri, *Rmaterial, RVmax, RTmax;
double *Rvrt;

extern void PolyregionConstruct(MESH_DATA& mesh);
extern void file_combiner (char *a,char *b);
extern int nver, ntr;
extern int qual_mas_res[];
extern double thresholds[];

static int compare_index(const void* a, const void* b) 
{
	Point1 a1=*(Point1*)a, b1=*(Point1*)b;
	if ((fabs(a1.point[0]-b1.point[0])+fabs(a1.point[1]-b1.point[1])+fabs(a1.point[2]-b1.point[2]))<3*GL_EPS)
		return 0;
	else
		if(((a1.point[0]-b1.point[0])>GL_EPS)||((fabs(a1.point[0]-b1.point[0])<GL_EPS)&&((a1.point[1]-b1.point[1])>GL_EPS))||((fabs(a1.point[0]-b1.point[0])+fabs(a1.point[1]-b1.point[1])<2*GL_EPS&&(a1.point[2]-b1.point[2])>GL_EPS)))
			return 1;
		else
			return -1;
}

static int compare_index1(const void* a, const void* b) 
{
	Point1 a1=*(Point1*)a, b1=*(Point1*)b;
	return a1.before-b1.before;
}

extern "C" int libaft_internal_surfmeshrefiner(int *nVRT, double *vrt, int *nTRI, int *tri, int *material, int Vmax, int Tmax)
{
	int i, j;
	MESH_DATA mesh;
	Point1 *pt;
	mesh.points_size = *nVRT;
	mesh.triangle_size = *nTRI;

	double *tmp;
	tmppoints=new POINT[mesh.points_size];
	mesh.points=new POINT[mesh.points_size];
	mesh.triangles=new TRIANGLE[mesh.triangle_size];

	for(i=0; i<mesh.points_size; i++)
	{   
		mesh.points[i].x = vrt[3*i];
		mesh.points[i].y = vrt[3*i+1];
		mesh.points[i].z = vrt[3*i+2];
	}

	for(i=0; i<mesh.triangle_size; i++)
	{
		mesh.triangles[i].vertex[0] = tri[3*i];
		mesh.triangles[i].vertex[1] = tri[3*i+1]; 
		mesh.triangles[i].vertex[2] = tri[3*i+2];
		mesh.triangles[i].color = material[i];
		mesh.triangles[i].color = max(0, mesh.triangles[i].color);
	}

	RVmax = Vmax;
	RTmax = Tmax;


	Rvrt = vrt;
	Rtri = tri;
	Rmaterial = material;

	PolyregionConstruct(mesh);

	delete [] tmppoints;
	delete [] mesh.points;
	delete [] mesh.triangles;
	pt = new Point1[nver];
	for(i=0;i<nver;i++)
	{
		pt[i].point=&vrt[3*i];
		pt[i].before=i+1;
		pt[i].after=i+1;
	}
	qsort(pt, nver, sizeof(Point1), compare_index);
	pt[0].after=1;
	for(i=1;i<nver;i++)
		if((fabs(pt[i].point[0]-pt[i-1].point[0])+fabs(pt[i].point[1]-pt[i-1].point[1])+fabs(pt[i].point[2]-pt[i-1].point[2]))<3*GL_EPS)
			pt[i].after=pt[i-1].after;
		else 
		{
			pt[i].after=pt[i-1].after+1;

		}
	*nVRT = pt[nver-1].after;
	if (*nVRT > Vmax) {
//		TODO warning
//		libaft_3d_warn("Vmax");
		*nVRT = Vmax;
	}
	tmp = new double[3*(*nVRT)];
	for(j=0;j<3;j++)
		tmp[j]=pt[0].point[j];
	for(i=1;i<nver;i++)
		if(pt[i].after != pt[i-1].after)
			for(j=0;j<3;j++)
				tmp[3*(pt[i].after-1)+j]=pt[i].point[j];
	/*	delete [] vrt;
		vrt=tmp;*/
	for(i=0; i<3**nVRT; i++)
		vrt[i]=tmp[i];

	delete [] tmp;


	//for(i=1;i<nver;i++)
	//	printf("%lf %lf %lf\n",pt[i].point[0],pt[i].point[1],pt[i].point[2]);
	qsort(pt, nver, sizeof(Point1), compare_index1);
	*nTRI = ntr;
	if (*nTRI > Tmax) {
//		TODO warning
//		libaft_3d_warn("Tmax");
		*nTRI = Tmax;
	}
	for(i=0; i<3*(*nTRI); i++)
		tri[i]=pt[tri[i]-1].after;
	delete [] pt;
	return 0;
}
