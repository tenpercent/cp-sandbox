#ifndef MESH_H
#define MESH_H

#define PI 3.14159265358979f  
#define EPS 10e-6

#include "GL/gl.h"

class mesh
{
public:
// Coordinates of the mesh:
    GLdouble x, y, z;
    
// Angles:
    double pitch, yaw, roll;
    
    int nV, nF;
    GLdouble *Vertex;
    GLint *Index;
    
    mesh ();
    mesh (int, int);
    ~mesh ();

    bool hide;
    int in;

    int mesh_smv (char *name, int color);
    int mesh_front (int nV_, int nF_, double *Vertex_, int *Index_);
    int mesh_write_smv (char *name);

    int primitive_paral (GLdouble x, GLdouble y, GLdouble z, double MS);
    int primitive_sphere (GLdouble x, double MS);
    int primitive_cylinder (GLdouble r, GLdouble h, double MS);
    int move (GLdouble x, GLdouble y, GLdouble z, double rot_x, double rot_y, double rot_z);

    int intersect (mesh *msh);
    int triInt (mesh *msh, int numF, int *ret, GLdouble pt1[3], GLdouble pt2[3]);
    int isInside (GLdouble x, GLdouble y, GLdouble z);
    int insideTri (GLdouble *u, GLdouble *v, GLdouble *x, GLdouble *y);

    void printPoint (GLdouble *x);
    double dist_point_curve (double *point, double *curve, int len, int *ret);
};

#endif // MESH_H
