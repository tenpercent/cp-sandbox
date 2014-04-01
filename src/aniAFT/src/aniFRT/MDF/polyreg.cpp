#include "def.h"
#include "stack.h"

static double MINABS_WIDTH = 0.0;
static double MIN_WIDTH = 0.05; ///defines minimal relative triangle's width
///relative width = absolute width/length of a boundary edge

static double PLANE_EPS = 0.01; ///define maximal polyregion's plane deviation
static double PLANE_EPS2 = 0.3; ///define maximal polyregion's plane deviation for "bad" triangles
static double MAX_DIST = 100; ///defines maximal relative distance between 
///a polyregion's triangle and the polyregion's plane
///absolute distance = relative distance * lenght of the longest edge

static double MAX_WIDTH = -1.0;

extern "C" int libaft_internal_surfmeshrefiner_poly_setup(double minabs, double min, double eps, double eps2/*, double dist*/) {
    if (minabs>0.0) MINABS_WIDTH = minabs;
    if (min>0.0)    MIN_WIDTH    = min;
    if (eps>0.0)    PLANE_EPS    = eps;
    if (eps2>0.0)   PLANE_EPS2   = eps2;
//    if (dist>0.0)   MAX_DIST     = dist;
    return 0;
}

extern "C" int libaft_internal_surfmeshrefiner_poly_setup_extra(double minabs, double min, double eps, double eps2, double dist) {
    if (minabs>0.0) MINABS_WIDTH = minabs;
    if (min>0.0)    MIN_WIDTH    = min;
    if (eps>0.0)    PLANE_EPS    = eps;
    if (eps2>0.0)   PLANE_EPS2   = eps2;
    if (dist>0.0)   MAX_DIST     = dist;
    return 0;
}

extern "C" void libaft_internal_surfmeshrefiner_poly_setup_extralim(double lim) {
    MAX_WIDTH = lim;
}

extern int vector_vector_classify (const POINT& a, const POINT& b, const POINT& c, const POINT& d);
extern void transform_determing (TRIANGLE XYZ, POINT *points,double **A, double **A1, POINT& d);
extern int MeshPolyreg3D(MESH_DATA& mesh, POLYREGION& polyreg, int bound_num);
extern POINT transform (double **A, POINT d, POINT x, int dir);
extern POINT *tmppoints;

double global_eps = PLANE_EPS;
int nver=0;
int ntr=0;
double max_length = 0;
double thresholds[6] = {1.e-5, 1.e-4, 1.e-3, 1.e-2, 0.1, 1};
EDGE* global_edges; ///contain all edges of the task

///Calculates triangle quality value
double TriangleValue(POINT vertices[3])
{
    double sides[3], max_side=0, p=0;

    for(int i=0; i<3; i++)
    {
	int i1 = i+1;
	if(i1>=3)
	    i1=0;

	sides[i]= sqrt(pow(vertices[i1].x-vertices[i].x, 2) + pow(vertices[i1].y-vertices[i].y, 2) + pow(vertices[i1].z-vertices[i].z, 2));

	if(max_side < sides[i])
	    max_side = sides[i];

	p+=sides[i];
    }

    p/=2;

    double r = sqrt((p-sides[0])*(p-sides[1])*(p-sides[2])/p);

    return r/max_side*2*sqrt(3.0);
}

///Calculate the width of the triangle (square/max. side)
double TriangleArea(POINT vertices[3])
{
    double sides[3], p=0.0;

    for(int i=0; i<3; i++)
    {
	int i1 = (i+1) % 3;

	sides[i]= sqrt(pow(vertices[i1].x-vertices[i].x, 2) + pow(vertices[i1].y-vertices[i].y, 2) + pow(vertices[i1].z-vertices[i].z, 2));

	p+=sides[i];
    }

    p/=2;

    ///Heron formula
    double Square = sqrt(p*(p-sides[0])*(p-sides[1])*(p-sides[2]));

    return Square;
}

///Calculate the width of the triangle (square/max. side)
double TriangleWidth(POINT vertices[3])
{
    double sides[3], max_side=0, p=0;
    double min_side = 1000;;

    for(int i=0; i<3; i++)
    {
	int i1 = i+1;
	if(i1>=3)
	    i1=0;

	sides[i]= sqrt(pow(vertices[i1].x-vertices[i].x, 2) + pow(vertices[i1].y-vertices[i].y, 2) + pow(vertices[i1].z-vertices[i].z, 2));

	if(max_side < sides[i])
	    max_side = sides[i];

	if(min_side > sides[i])
	    min_side = sides[i];

	p+=sides[i];
    }

    p/=2;

    ///Heron formula
    double Square = sqrt(p*(p-sides[0])*(p-sides[1])*(p-sides[2]));
    double value = min(min_side, 2*Square/max_side);

    return value;
}

///calculate normal vector for the plane containing pt
POINT CalcPlaneNormal(POINT pt[3])
{
    POINT pt21=pt[1]-pt[0], pt31=pt[2]-pt[0];
    POINT normal;

    normal.x = pt21.y*pt31.z - pt21.z*pt31.y;
    normal.y = -(pt21.x*pt31.z - pt21.z*pt31.x);
    normal.z = pt21.x*pt31.y - pt21.y*pt31.x;

    return normal/normal.R();
}

///check if triangles trig1 and trig2 are neighbours and trig2 almost belongs to the plane with normal normal1
bool CheckForNeighbour(const TRIANGLE& trig1, const TRIANGLE& trig2, const MESH_DATA& mesh, const POLYREGION* polyreg)
{
    POINT vertices1[3], vertices2[3];
    int num_equal = 0, i, j;
    POINT equal_pt[2], not_equal_pt[2];

    if(trig1.color!=trig2.color)
	return false;

    for(i=0; i<3; i++)
    {
	vertices1[i] = mesh.points[trig1.vertex[i]-1];
	vertices2[i] = mesh.points[trig2.vertex[i]-1];
    }

    double norm_eps = global_eps;
    double quality = TriangleValue(vertices2);
    if(quality < 0.1)
	norm_eps = PLANE_EPS2;

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    if(vertices1[i] == vertices2[j])
	    {
		equal_pt[num_equal++] = vertices1[i];
		break;
	    }

    if(num_equal<2)
	return false;

    num_equal = 0;
    for(i=0; i<3; i++)
    {
	for(j=0; j<2; j++)
	{
	    POINT pt = equal_pt[j]-vertices1[i];
	    if(pt.R()<0.00001)
		break;
	}
	if(j == 2)
	{
	    not_equal_pt[num_equal++] = vertices1[i];
	    break;
	}
    }

    for(i=0; i<3; i++)
    {
	for(j=0; j<2; j++)
	{
	    POINT pt = equal_pt[j]-vertices2[i];
	    if(pt.R()<0.00001)
		break;
	}
	if(j == 2)
	{
	    not_equal_pt[num_equal++] = vertices2[i];
	    break;
	}
    }

    memcpy(vertices1, equal_pt, sizeof(equal_pt));
    vertices1[2] = not_equal_pt[0];

    memcpy(vertices2, equal_pt, sizeof(equal_pt));
    vertices2[2] = not_equal_pt[1];

    POINT normal2 = CalcPlaneNormal(vertices2);
    double scalar_product = min(fabs(1-normal2*polyreg->normal), fabs(1+normal2*polyreg->normal));
    if(scalar_product > norm_eps)
	return false;
    else
    {
	POINT normal = CalcPlaneNormal(vertices1);
	if(normal2*normal < -0.000001)
	{
	    for(j=0; j<3; j++)
		if(fabs(polyreg->a*vertices2[j].x+polyreg->b*vertices2[j].y+polyreg->c*vertices2[j].z-polyreg->p) > MAX_DIST*max_length)
		    return false;
	    return true;
	}
	else
	    return false;
    }
}

/**method calculates new polyregions' boundaries and call function of constructing new triangulation;
\param mesh - initial data(loaded from file), \param cur_polyreg - all information about the current polyregion,
\param boundary - boundary(s) of the current polyregion (indices of the boundary points),
\param global_edges - all edges of the task.*/  
void PolyregProcessing(MESH_DATA& mesh, POLYREGION* cur_polyreg, ITYPE* boundary, EDGE* global_edges)
{
    POINT* new_boundary = new POINT[MAX_BOUNDARY];

    ITYPE boundary_size = 0;
    int repeat_index = -1, start_ind = 0; ///for defining subboundaries

    for(int j=0; j<cur_polyreg->boundary_size-1; j++)
    {
	ITYPE left = boundary[j]-1, right = boundary[j+1]-1;
	POINT left_pt = mesh.points[left], right_pt = mesh.points[right];
	POINT difference = right_pt-left_pt;
	double length = difference.R(); ///edge's length
	double min_width[2] = {100, 100};

	if(repeat_index == -1)
	{
	    if(j>0) ///the current boundary is finished. Go to the next one
	    {
		repeat_index = right;
		continue;
	    }
	    else
		repeat_index = left;
	}

	///find min_width for the edge's vertices
	for(int m=0; m<2; m++)
	{
	    POINT* vertice;
	    if(!m)
		vertice = &left_pt;
	    else
		vertice = &right_pt;

	    if(vertice->min_width<0)
	    {
		vertice->min_width = 100;
		///consider all edges which are incedent to the given point
		for(int i=0; i<vertice->edge_size; i++)
		{
		    EDGE& vrem_edge = global_edges[vertice->edges[i]];

		    if(vrem_edge.min_width<0)
		    {
			vrem_edge.min_width = 100;
			///define minimal width among all triangles containing vrem_edge
			for(int k=0; k<vrem_edge.triangle_size; k++)
			{
			    POINT vertices[3];
			    TRIANGLE& triangle = mesh.triangles[vrem_edge.triangles[k]];
			    for(int l=0; l<3; l++)
				vertices[l] = mesh.points[triangle.vertex[l]-1];
			    vrem_edge.min_width = min(vrem_edge.min_width, TriangleWidth(vertices));
			}
		    }
		    vertice->min_width = min(vertice->min_width, vrem_edge.min_width);
		}
	    }

	    min_width[m] = vertice->min_width;
	}

	for(int w=0; w<2; w++) {
	    if(min_width[w]/length < MIN_WIDTH)
		min_width[w] = length*MIN_WIDTH;
	    if(min_width[w] < MINABS_WIDTH)
		min_width[w] = MINABS_WIDTH;
	    if((MAX_WIDTH>0.0) && (min_width[w] > MAX_WIDTH))  min_width[w] = MAX_WIDTH;
	}

	///define interpolated function describing grid step on the current edge
	double part_count = 2*length/(min_width[0]+min_width[1])-1; ///number of parts that the current edge must be divided into
	int part_num = (int)(part_count+0.5);
	min_width[0] *= (double(part_count+1)/(part_num+1));
	min_width[1] *= (double(part_count+1)/(part_num+1));

	part_count = part_num;

	double coef = (min_width[1]-min_width[0])/part_count;
	difference /= length;

	double coef_ind = part_count/part_num;

	if(repeat_index == right)
	{
	    //part_count += 1;
	    //part_num++;
	    repeat_index = -1;
	}

	if(boundary_size+part_num >= MAX_BOUNDARY)
	{
	    printf("It is necessary to increase MAX_BOUNDARY parameter!!!");
	    exit(1);
	}

	///divide initial and construct new boundary
	new_boundary[boundary_size++] = left_pt;
	for(int k=0; k<part_num; k++)
	{
	    double ind = coef_ind*k;
	    new_boundary[boundary_size++] = left_pt + difference*(min_width[0]+coef*ind/2)*(ind+1);
	}

	if(repeat_index == -1)
	{
	    new_boundary[boundary_size++] = new_boundary[start_ind];
	    start_ind = boundary_size;
	}
    }

    cur_polyreg->boundary_size = boundary_size;
    cur_polyreg->boundary = new POINT[boundary_size];
    memcpy(cur_polyreg->boundary, new_boundary, sizeof(POINT)*boundary_size);
    delete [] new_boundary;
}

///build polyregions consisting of points and triangles from mesh
void PolyregionConstruct(MESH_DATA& mesh)
{
    int i, k, j, m;
    int wrong = 0, repair = 0;
    ITYPE poly_size=1;
    ITYPE *triangles_tmp = new ITYPE[mesh.triangle_size];
    ITYPE *triangles = new ITYPE[mesh.triangle_size];
    ITYPE *points = new ITYPE[mesh.points_size], first_pt, last_pt;
    bool *marks = new bool[mesh.triangle_size], *reg_marks = new bool[mesh.triangle_size];
    int left_to_mark = mesh.triangle_size, begin_tr=0, begin_tr_prev, cur_index=0, polyreg_size=0;
    memset(marks, 0, sizeof(bool)*mesh.triangle_size);
    typedef PAIR<TRIANGLE, ITYPE> TRIG_TYPE;
    CStack<TRIG_TYPE> stack(mesh.triangle_size);

    EDGE boundary[MAX_BOUNDARY];
    ITYPE real_boundary[MAX_BOUNDARY], vrem_boundary[MAX_BOUNDARY];
    int edge_size=0;

    global_edges = new EDGE[MAX_EDGES];
    int glob_edge_size=0;

    int qual_mas[6], poly_sizes[5];
    //double thresholds[6] = {1.e-2, 1.e-1, 0.5, 0.7, 0.9, 1};
    int max_sizes[5] = {1, 10, 100, 1000, 10000};

    memset(qual_mas, 0, sizeof(qual_mas));
    memset(poly_sizes, 0, sizeof(poly_sizes));

    ///create surface edges
    for(i=0; i<mesh.triangle_size; i++)
    {
	ITYPE* vertices = mesh.triangles[i].vertex;

	POINT v_tr[3] = {mesh.points[vertices[0]-1], mesh.points[vertices[1]-1], mesh.points[vertices[2]-1]};
	double res = TriangleValue(v_tr);
	for(k=0; k<6; k++)
	    if(res < thresholds[k])
	    {
		qual_mas[k]++;
		break;
	    }

	for(k=0; k<3; k++)
	{
	    bool flag = true;

	    int k1 = k+1;
	    if(k1>2)
		k1 = 0;

	    PAIR<ITYPE, ITYPE> pair = PAIR<ITYPE, ITYPE>(vertices[k], vertices[k1]);
	    POINT* point = &v_tr[k];

	    for(j=0; j<point->edge_size; j++)
	    {
		int index = point->edges[j];
		EDGE* edge = &global_edges[index];
		if(edge->boundary == pair)
		{
		    mesh.triangles[i].edges[k] = index;
		    int tr_esize = global_edges[index].triangle_size++;
		    global_edges[index].triangles[tr_esize] = i;

		    if(tr_esize >= MAX_EDGE_TRIANGLE)
		    {
			printf("It is necessary to increase MAX_EDGE_TRIANGLE parameter!!!");
			exit(1);
		    }

		    flag = false;
		    break;
		}
	    }

	    if(flag)
	    {
		mesh.triangles[i].edges[k] = glob_edge_size;

		global_edges[glob_edge_size].index = glob_edge_size;
		global_edges[glob_edge_size].triangles[0] = i;
		global_edges[glob_edge_size].triangle_size++;
		global_edges[glob_edge_size].boundary = pair;

		POINT dist = mesh.points[pair.elem2-1]-mesh.points[pair.elem1-1];
		if(max_length < dist.R())
		    max_length = dist.R();

		int pt_esize1 = mesh.points[pair.elem1-1].edge_size++, pt_esize2 = mesh.points[pair.elem2-1].edge_size++;
		mesh.points[pair.elem1-1].edges[pt_esize1] = glob_edge_size;
		mesh.points[pair.elem2-1].edges[pt_esize2] = glob_edge_size;
		glob_edge_size++;

		if(pt_esize1 >= MAX_DEGREE || pt_esize2 >= MAX_DEGREE)
		{
		    printf("It is necessary to increase MAX_DEGREE parameter!!!");
		    exit(1);
		}

		if(glob_edge_size>=MAX_EDGES)
		{
		    printf("It is necessary to increase MAX_EDGES parameter!!!");
		    exit(1);
		}
	    }
	}
    }

#ifdef DEBUGGER
    ///Find maximal degree among all vertices and index of the found vertice
    int max_sz=0 ,ind;
    for(k=0; k<mesh.points_size; k++)
	if(max_sz<mesh.points[k].edge_size)
	{
	    max_sz = mesh.points[k].edge_size;
	    ind = k;
	}

    ///Find maximal number of triangles with a common edge
    int min_tr=100, max_tr=0;
    for(k=0; k<glob_edge_size; k++)
	if(max_tr < global_edges[k].triangle_size)
	{
	    max_tr = global_edges[k].triangle_size;
	    ind = k;
	}
	else if(min_tr > global_edges[k].triangle_size)
	    min_tr = global_edges[k].triangle_size;
#endif ///DEBUGGER

    ///until all triangles become considered
    while(left_to_mark>0)
    {
	poly_size=1;

	short repeat = 0;
	do
	{
	    POLYREGION* cur_polyreg = new POLYREGION;
	    double area, maxarea = 0.0;

	    ///find the first unconsidered triangle
	    cur_index = -1;
	    for(i=0; i<mesh.triangle_size; i++)
		if(!marks[i])
		{
		    POINT vertices[3];
		    int l;
		    for(l=0; l<3; l++)
			vertices[l] = mesh.points[mesh.triangles[i].vertex[l]-1];
		    area = TriangleArea(vertices);
		    if ((cur_index<0) || (area > maxarea)) {
			maxarea = area;
			cur_index = i;
		    }
		}
	    begin_tr_prev = begin_tr;
	    //begin_tr = cur_index+1;

	    ///stack for recursion for an advancing front algorithm
	    stack.Push(TRIG_TYPE(mesh.triangles[cur_index], cur_index));

	    triangles_tmp[0]=cur_index;
	    transform_determing (mesh.triangles[cur_index], mesh.points, cur_polyreg->m_rotate, cur_polyreg->m_rotate_back, cur_polyreg->m_move);

	    ///define a current polyregion
	    POINT vertices[3];
	    for(i=0; i<3; i++)
		vertices[i] = mesh.points[mesh.triangles[cur_index].vertex[i]-1];
	    cur_polyreg->normal = CalcPlaneNormal(vertices);

	    for(int vert=0; vert<3; vert++)
		tmppoints[mesh.triangles[cur_index].vertex[vert]-1] =
		    transform(cur_polyreg->m_rotate, cur_polyreg->m_move, vertices[vert], 1);

	    cur_polyreg->p = -(vertices[0].x*cur_polyreg->normal.x+vertices[0].y*cur_polyreg->normal.y+vertices[0].z*cur_polyreg->normal.z);
	    double sign = 1;
	    if(cur_polyreg->p > 0)
		sign = -1;
	    cur_polyreg->a = cur_polyreg->normal.x*sign;
	    cur_polyreg->b = cur_polyreg->normal.y*sign;
	    cur_polyreg->c = cur_polyreg->normal.z*sign;
	    cur_polyreg->p *= sign;

	    cur_polyreg->color = mesh.triangles[cur_index].color;

	    memset(reg_marks, 0, sizeof(bool)*mesh.triangle_size);
	    marks[cur_index] = true;
	    reg_marks[cur_index] = true;

	    edge_size = 0;
	    //first_pt = MAX_POINTS;
	    first_pt = mesh.points_size;
	    last_pt = 0;
	    memset(points, 0, sizeof(ITYPE)*mesh.points_size);

	    while(stack.Size())
	    {
		TRIG_TYPE trig;
		stack.Pop(&trig);

		//fill current polyregion
		cur_index = trig.elem2;
		triangles[cur_polyreg->triangle_size] = cur_index;
		cur_polyreg->triangle_size++;
		left_to_mark--;

		///find all neighbours of the given triangle and save its vertices
		TRIANGLE* triangle = &trig.elem1;
		for(k=0; k<3; k++)
		{
		    int cur_vertex = triangle->vertex[k]-1;
		    if(cur_vertex < first_pt)
			first_pt = cur_vertex;
		    if(cur_vertex > last_pt)
			last_pt = cur_vertex;
		    points[cur_vertex]++;

		    EDGE* edge = &global_edges[triangle->edges[k]];
		    int neighbour_count=0;
		    if(edge->triangle_size > 2)
		    {
			for(i=0; i<edge_size; i++)
			    if(edge->boundary == boundary[i].boundary)
				break;
			if(i == edge_size)
			    boundary[edge_size++] = *edge;
			continue;
		    }

		    for(i=0; i<edge->triangle_size; i++)
		    {
			int ind = edge->triangles[i];
			if(ind != cur_index)
			    neighbour_count++;
			else
			    continue;

			if(!marks[ind])
			{
			    bool intersect = true;
			    POINT add_point;
			    int vrem_ind;

			    for(int v=0; v<3; v++)
			    {
				vrem_ind = mesh.triangles[ind].vertex[v];
				if(vrem_ind != edge->boundary.elem1 && vrem_ind != edge->boundary.elem2)
				{
				    add_point = mesh.points[vrem_ind-1];
				    break;
				}
			    }

			    tmppoints[vrem_ind-1] =
				transform(cur_polyreg->m_rotate, cur_polyreg->m_move, add_point, 1);

			    if(CheckForNeighbour(*triangle, mesh.triangles[ind], mesh, cur_polyreg))
			    {
				if(repeat)
				{
				    TRIANGLE* all_trig = mesh.triangles;
				    TRIANGLE* vrem_tr = &all_trig[ind];
				    for(j=0; j<3 && intersect; j++)
				    {
					EDGE* vrem_ed = &global_edges[vrem_tr->edges[j]];
					if(*vrem_ed == *edge)
					    continue;

					POINT pt1 = tmppoints[vrem_ed->boundary.elem1-1], pt2 = tmppoints[vrem_ed->boundary.elem2-1];
					for(int l=0; l<poly_size && intersect; l++)
					{
					    TRIANGLE* vrem_tr1 = &all_trig[triangles_tmp[l]];
					    for(m=0; m<3 && intersect; m++)
					    {
						EDGE* vrem_ed1 = &global_edges[vrem_tr1->edges[m]];
						if(vector_vector_classify(pt1, pt2,
							    tmppoints[vrem_ed1->boundary.elem1-1],
							    tmppoints[vrem_ed1->boundary.elem2-1]) == 1 && intersect)
						    intersect=false;
					    }
					}
				    }
				}
			    }
			    else
				intersect=false;

			    if(intersect)
			    {
				stack.Push(TRIG_TYPE(mesh.triangles[ind], ind));
				marks[ind] = true;
				reg_marks[ind] = true;
				triangles_tmp[poly_size]=ind;
				poly_size++;
			    }
			    else
				boundary[edge_size++] = *edge; ///add new edge to polyregion's border
			}
			else if(!reg_marks[ind])
			    boundary[edge_size++] = *edge;
		    }

		    ///the edge is a part of the boundary of a hole
		    if(!neighbour_count)
			boundary[edge_size++] = *edge;
		}
	    }

	    ///remove boundary edges which belong to more than one triangles
	    for(i=0; i<edge_size; i++)
	    {
		int edge = boundary[i].index;
		int count = 0;
		for(j=0; j<cur_polyreg->triangle_size; j++)
		{
		    TRIANGLE* trig = &mesh.triangles[triangles[j]];
		    if(edge == trig->edges[0] || edge == trig->edges[1] || edge == trig->edges[2])
			count++;
		}

		if(count != 1)
		{
		    memmove(&boundary[i], &boundary[i+1], sizeof(EDGE)*(edge_size-i-1));
		    edge_size--;
		    i--;
		}
	    }

	    ///save all polyregion's points
	    for(i = first_pt; i<=last_pt; i++)
		if(points[i])
		    points[cur_polyreg->points_size++] = i;

#ifdef DEBUGGER
	    ///check that all polygone's edges belong to less than 3 triangles
	    for(i=0; i<cur_polyreg->triangle_size; i++)
		for(k=0; k<3; k++)
		{
		    int edge = mesh.triangles[triangles[i]].edges[k];
		    int edge_belong=1;
		    for(j=i+1; j<cur_polyreg->triangle_size; j++)
			for(m=0; m<3; m++)
			    if(edge == mesh.triangles[triangles[j]].edges[m])
				edge_belong++;

		    if(edge_belong>2)
		    {
			printf("The polygone's edge belongs to more than 2 triangles!!! Decrease PLANE_EPS parameter.");
			exit(1);
		    }
		}

	    int count = 0;
	    for(j=0; j<edge_size; j++)
	    {
		count = 0;
		int pt_num = boundary[j].boundary.elem1;
		for(i=0; i<edge_size; i++)
		    if(boundary[i].boundary.elem1 == pt_num || boundary[i].boundary.elem2 == pt_num)
			count++;
		if(count%2 != 0)
		    int stop = 1;

		count = 0;
		pt_num = boundary[j].boundary.elem2;
		for(i=0; i<edge_size; i++)
		    if(boundary[i].boundary.elem1 == pt_num || boundary[i].boundary.elem2 == pt_num)
			count++;
		if(count%2 != 0)
		    int stop = 1;
	    }
#endif

	    real_boundary[0] = boundary[0].boundary.elem1;
	    int start_pos = boundary[0].boundary.elem1;
	    boundary[0].boundary.elem1 = -1;
	    cur_polyreg->boundary_size++;
	    ITYPE last_pos = boundary[0].boundary.elem2;
	    int cycle = 0, begin=1, edge_consider = 1;

	    ///put boundary's points in a right order
	    int number_of_boundaries = 1; ///full number of boundaries (external and internal) 
	    while(edge_consider<edge_size)
	    {
		if(cycle == edge_size)
		    for(i=1; i<edge_size; i++)
			if(boundary[i].boundary.elem1>=0)
			{
			    number_of_boundaries++;
			    if(boundary[i].boundary.elem1 != start_pos)
				real_boundary[cur_polyreg->boundary_size++] = start_pos;

			    real_boundary[cur_polyreg->boundary_size++] = boundary[i].boundary.elem1;
			    start_pos = boundary[i].boundary.elem1;
			    last_pos = boundary[i].boundary.elem2;
			    boundary[i].boundary.elem1 = -1;
			    edge_consider++;
			    begin = i+1;
			    break;
			}

		///find a point which coincides with the current boundary's end (last_pos)
		for(cycle=begin; cycle<edge_size; cycle++)
		    if(boundary[cycle].boundary.elem1>=0)
		    {
			if(boundary[cycle].boundary.elem1 == last_pos)
			{
			    real_boundary[cur_polyreg->boundary_size++] = boundary[cycle].boundary.elem1;
			    if(start_pos == boundary[cycle].boundary.elem1)
				real_boundary[cur_polyreg->boundary_size++] = boundary[cycle].boundary.elem1;
			    last_pos = boundary[cycle].boundary.elem2;
			    boundary[cycle].boundary.elem1 = -1;
			    edge_consider++;
			    break;
			}
			else if(boundary[cycle].boundary.elem2 == last_pos)
			{
			    real_boundary[cur_polyreg->boundary_size++] = boundary[cycle].boundary.elem2;
			    if(start_pos == boundary[cycle].boundary.elem2)
				real_boundary[cur_polyreg->boundary_size++] = boundary[cycle].boundary.elem2;
			    last_pos = boundary[cycle].boundary.elem1;
			    boundary[cycle].boundary.elem1 = -1;
			    edge_consider++;
			    break;
			}
		    }
		begin=0;
	    }
	    real_boundary[cur_polyreg->boundary_size++] = start_pos;

	    number_of_boundaries=1;
	    int start_bd = real_boundary[0];

	    ///extract embedded boundaries
	    for(i=1; i<cur_polyreg->boundary_size-1; i++)
	    {
		for(j=i+1; j < cur_polyreg->boundary_size-1; j++)
		{
		    if(real_boundary[j] == start_bd)
			break;
		    if((real_boundary[i] == real_boundary[j]) && (real_boundary[i] != start_bd))
		    {
			memcpy(vrem_boundary, &real_boundary[i], sizeof(ITYPE)*(j-i+1));
			memmove(&real_boundary[i+1], &real_boundary[j+1], sizeof(ITYPE)*(cur_polyreg->boundary_size-j-1));
			memcpy(&real_boundary[i+cur_polyreg->boundary_size-j], vrem_boundary, sizeof(ITYPE)*(j-i+1));
			cur_polyreg->boundary_size++;
			break;
		    }
		}

		if(real_boundary[i] == start_bd)
		{
		    number_of_boundaries++;
		    start_bd = real_boundary[i+1];
		    i++;
		}
	    }

	    //        printf("\rCurrent polyregion %d   Triangles left to consider %d", polyreg_size, left_to_mark);
	    //		  fprintf(stderr, "\rCurrent polyregion %6d   Triangles left to consider %6d", polyreg_size, left_to_mark); // added by Alexander Danilov
#ifdef SHOWPROGRESS
	    printf("\t\t\t\t\t\tPr: %6d, Tl: %6d\r", polyreg_size, left_to_mark);
	    fflush(stdout);
#endif

	    ///(then replace by cur_polyreg->triangles = triangles and cur_polyreg->points = points)
	    cur_polyreg->triangles = new ITYPE[cur_polyreg->triangle_size];
	    memcpy(cur_polyreg->triangles, triangles, sizeof(ITYPE)*cur_polyreg->triangle_size);

	    cur_polyreg->points = new ITYPE[cur_polyreg->points_size];
	    memcpy(cur_polyreg->points, points, sizeof(ITYPE)*cur_polyreg->points_size);

#ifdef DEBUGGER
	    //////TEMPORARY!!! Just for checking constructed polyregions!!!//////////////////////////////////////////////////////
	    FILE* u ;
	    u = fopen("boundary_init", "w");

	    fprintf(u,"%d\n",cur_polyreg->boundary_size-number_of_boundaries);

	    double st_x = mesh.points[real_boundary[0]-1].x, st_y = mesh.points[real_boundary[0]-1].y;
	    for(m=0; m<cur_polyreg->boundary_size-1; m++)
	    {
		fprintf(u,"%lf %lf ",mesh.points[real_boundary[m]-1].x,mesh.points[real_boundary[m]-1].y);
		fprintf(u,"%lf %lf\n",mesh.points[real_boundary[m+1]-1].x,mesh.points[real_boundary[m+1]-1].y);

		if(fabs(mesh.points[real_boundary[m+1]-1].x - st_x)<0.00001 && fabs(mesh.points[real_boundary[m+1]-1].y - st_y)<0.00001)
		{
		    if(m+2<cur_polyreg->boundary_size)
		    {
			st_x = mesh.points[real_boundary[m+2]-1].x;
			st_y = mesh.points[real_boundary[m+2]-1].y;
		    }
		    m++;
		}
	    }
	    fclose(u);
#endif
	    /*		int prev_size = cur_polyreg->boundary_size; */
	    PolyregProcessing(mesh, cur_polyreg, real_boundary, global_edges);

	    for(k=0; k<5; k++)
		if(cur_polyreg->triangle_size <= max_sizes[k])
		{
		    poly_sizes[k]++;
		    break;
		}

	    int tr_size = cur_polyreg->triangle_size;
//	    mesh.triangles[cur_polyreg->triangles[0]].color = polyreg_size; /// FIXME! REMOVE ME !
	    if(!MeshPolyreg3D(mesh,*cur_polyreg,number_of_boundaries))
	    {
		poly_size = 1;
		left_to_mark += tr_size;
		begin_tr = begin_tr_prev;
		for(i=0; i<tr_size; i++)
		    marks[triangles[i]] = false;

		if(repeat)
		{
		    ///output bad polygon
		    if(fabs(global_eps - PLANE_EPS) < 0.0001)
		    {
			repair++;
			fprintf(stderr, "anifrtmdf: bad_polygon!\n");
			/*                    char str[30] = "";
					      sprintf(str, "bad_polygon%d", repair);
					      int size = strlen(str);
					      char* vrem_str = new char[size+1];
					      memcpy(vrem_str, str, sizeof(char)*size);
					      vrem_str[size] = '\0';

			///Write results
			FILE* f = fopen(vrem_str, "w");
			delete [] vrem_str;

			fprintf(f, " %d", cur_polyreg->points_size);
			fprintf(f, " %d", cur_polyreg->triangle_size);
			fprintf(f, "    #nodes and #faces\n");

			for(j=0; j<cur_polyreg->points_size; j++)
			{
			fprintf(f, "  %lf", mesh.points[cur_polyreg->points[j]].x);
			fprintf(f, "  %lf", mesh.points[cur_polyreg->points[j]].y);
			fprintf(f, "  %lf\n", mesh.points[cur_polyreg->points[j]].z);
			}

			for(j=0; j<cur_polyreg->triangle_size; j++)
			{
			int index = cur_polyreg->triangles[j];
			for(k=0; k<3; k++)
			for(i=0; i<cur_polyreg->points_size; i++)
			if(mesh.triangles[index].vertex[k] == cur_polyreg->points[i]+1)
			{
			fprintf(f, " %d", i+1);
			break;
			}
			fprintf(f, " %d\n", mesh.triangles[index].color);
			}

			fclose(f);*/
		    }
		    global_eps /= 2;

		    if(global_eps > PLANE_EPS/1000)
			repeat = 0;
		}
		repeat++;
	    }
	    else
		repeat = 0;

#ifdef DEBUGGER
	    //////TEMPORARY!!! Just for checking constructed polyregions!!!//////////////////////////////////////////////////////

	    ///Write results
	    FILE* f = fopen("temp", "w");

	    fprintf(f, " %d", cur_polyreg->points_size);
	    fprintf(f, " %d", cur_polyreg->triangle_size);
	    fprintf(f, "    #nodes and #faces\n");

	    for(j=0; j<cur_polyreg->points_size; j++)
	    {
		fprintf(f, "  %lf", mesh.points[cur_polyreg->points[j]].x);
		fprintf(f, "  %lf", mesh.points[cur_polyreg->points[j]].y);
		fprintf(f, "  %lf\n", mesh.points[cur_polyreg->points[j]].z);
	    }

	    for(j=0; j<cur_polyreg->triangle_size; j++)
	    {
		int index = cur_polyreg->triangles[j];
		for(k=0; k<3; k++)
		    for(i=0; i<cur_polyreg->points_size; i++)
			if(mesh.triangles[index].vertex[k] == cur_polyreg->points[i]+1)
			{
			    fprintf(f, " %d", i+1);
			    break;
			}
		fprintf(f, " %d\n", mesh.triangles[index].color);
	    }

	    fclose(f);


	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif ///DEBUGGER

	    delete cur_polyreg;        
	}
	while(repeat == 1);

	global_eps = PLANE_EPS;

	if(repeat == 2)
	    wrong++;

	polyreg_size++;
    }
#ifdef SHOWPROGRESS
    printf("\t\t\t\t\t\t                      \r");
    fflush(stdout);
#endif

    //	printf("\nNumber of wrong polygones %d\n", wrong);

    //statistics output
    /*    FILE* statist = fopen("init_stat.txt", "w");

	  fprintf(statist, "Initial quality statistics\n");
	  fprintf(statist, "<%.0e    <%.0e    <%.0e    <%.0e    <%.0e    <%.0e\n", 
	  thresholds[0], thresholds[1], thresholds[2], thresholds[3], thresholds[4], thresholds[5]);
	  fprintf(statist, "%d           %d         %d           %d          %d         %d\n\n",
	  qual_mas[0], qual_mas[1], qual_mas[2], qual_mas[3], qual_mas[4], qual_mas[5]);

	  fprintf(statist, "Polyregions statistics (number of triangles)\n");
	  fprintf(statist, "<=%d    <=%d    <=%d    <=%d    <=%d \n", 
	  max_sizes[0], max_sizes[1], max_sizes[2], max_sizes[3], max_sizes[4]);
	  fprintf(statist, "%d       %d       %d       %d     %d",
	  poly_sizes[0], poly_sizes[1], poly_sizes[2], poly_sizes[3], poly_sizes[4]);

	  fclose(statist);*/

    delete [] points;
    delete [] marks;
    delete [] reg_marks;
    delete [] triangles;
    delete [] triangles_tmp;
    delete [] global_edges;
}
