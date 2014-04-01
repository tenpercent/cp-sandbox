#ifndef _DEF_H_
#define _DEF_H_

#define SHOWPROGRESS

#ifndef DEBUG
  #define DEBUG_ ///only for debug version
#endif ///DEBUG

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "stack.h"

const int MAX_BOUNDARY = 100000;       ///maximal number of points in a polyregion's boundary
const int MAX_EDGES   = 2000000;       ///maximal number of edges
const int MAX_EDGE_TRIANGLE = 9;       ///maximal number of triangles containing the common edge
const int MAX_DEGREE       = 50;       ///maximal degree of a vertice

//const double GL_EPS = 1.e-11;
const double GL_EPS = 1.e-11;

//#define DEBUGGER ///only for debug version

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

typedef int ITYPE;

struct Point1
{
    double *point;
    ITYPE before;
    ITYPE after;
};

///contains information about a vertice 
struct POINT
{
    POINT()
    {
        edge_size = 0;
        min_width = -9999;
    }
    
    bool operator==(const POINT& pt) const
    {
        if(pt.x == x && pt.y == y && pt.z == z)
            return true;
        else
            return false;
    }

    POINT operator+(const POINT& pt) const
    {
        POINT rez;
        
        rez.x = x + pt.x;
        rez.y = y + pt.y;
        rez.z = z + pt.z;

        return rez;
    }
    
    POINT operator-(const POINT& pt) const
    {
        POINT rez;
        
        rez.x = x - pt.x;
        rez.y = y - pt.y;
        rez.z = z - pt.z;

        return rez;
    }

    POINT operator-(void) const
    {
        POINT rez;
        
        rez.x = -x;
        rez.y = -y;
        rez.z = -z;

        return rez;
    }

    ///scalar product
    double operator*(const POINT& pt) const
    {
        return x*pt.x + y*pt.y + z*pt.z;
    }

    POINT operator/(const double& coef) const
    {
        POINT rez;
        
        rez.x = x/coef;
        rez.y = y/coef;
        rez.z = z/coef;

        return rez;
    }

    POINT operator*(const double& coef) const
    {
        POINT rez;
        
        rez.x = x*coef;
        rez.y = y*coef;
        rez.z = z*coef;

        return rez;
    }

    void operator*=(const double& coef)
    {
        POINT rez;
        
        x *= coef;
        y *= coef;
        z *= coef;
    }

	void operator/=(const double& coef)
    {
        POINT rez;
        
        x /= coef;
        y /= coef;
        z /= coef;
    }

	double R() const
    {
        return sqrt(x*x + y*y + z*z);
    }

    double x, y, z;

    ITYPE edge_size;
    ITYPE edges[MAX_DEGREE]; ///edges which are incident to the given point
    double min_width; ///width referred to the given point
};

///describes triangle
struct TRIANGLE
{
    ITYPE edges[3];
    ITYPE vertex[3];
    ITYPE color;
};

///contains all input data
struct MESH_DATA
{
    ITYPE points_size;
    POINT *points;

    ITYPE triangle_size;
    TRIANGLE *triangles;
};

///describes edge
struct EDGE
{
    EDGE()
    {
        triangle_size = 0;
        min_width = -9999;
    }

    bool operator==(const EDGE& edge) const
    {
        if((boundary.elem1 == edge.boundary.elem1 && boundary.elem2 == edge.boundary.elem2) ||
            (boundary.elem1 == edge.boundary.elem2 && boundary.elem1 == edge.boundary.elem2))
            return true;
        else
            return false;
    }
    
    ITYPE index;                 ///index of the given edge
    PAIR<ITYPE, ITYPE> boundary; ///boundary points' indices
    char triangle_size;
    ITYPE triangles[MAX_EDGE_TRIANGLE]; ///numbers of triangles which the given edge belongs to
    double min_width; ///minimal width among all triangles which contain the given edge
};

///describes full information about the given polyregion (points, triangles, plane's parameters)
struct POLYREGION
{
    POLYREGION()
    {
        boundary_size = 0;
        triangle_size = 0;
        points_size = 0;
        color = 0;
        boundary = NULL;
        triangles = NULL;
        points = NULL;

		m_rotate = new double*[3];
		m_rotate_back = new double*[3];
		for(int i=0; i<3; i++)
		{
			m_rotate[i] = new double[3];
			m_rotate_back[i] = new double[3];
		}
    }

    ~POLYREGION()
    {
        if(boundary)
            delete [] boundary;
        if(triangles)
            delete [] triangles;
        if(points)
            delete [] points;

		for(int i=0; i<3; i++)
		{
			delete [] m_rotate[i];
			delete [] m_rotate_back[i];
		}
		delete [] m_rotate;
		delete [] m_rotate_back;
    }

    POINT normal;    ///polyregion plane's normal
    double a, b, c, p;
    ITYPE color;

	POINT m_move;
	double **m_rotate, **m_rotate_back;

    ITYPE boundary_size;
    POINT* boundary; ///new boundary's points

    ITYPE points_size;
    ITYPE* points;  ///points's indices

    ITYPE triangle_size;
    ITYPE* triangles; ///triangles' indices
};

#endif //_DEF_H_
