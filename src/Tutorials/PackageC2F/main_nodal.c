#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include  "ani3D.h"



int main() 
{
/* local variables */
  int     i, n;
  int     nP, nF, nE;
  int     face[2], tet[3], label, max_iters;
  double  xy[3], Q, Lp;

/* create a pointer to the main structure */
  ani3D   aniReal, *ani;

/* mesh control parameters */
  int     nEStar;
  double  *metric, mem_factor;
  char    *fname;


/* initialize a mesh with nEStar triangles */
  ani = &aniReal;

  nEStar = 20000;
  mem_factor = 1.5;
  ani3D_INIT( ani, nEStar, mem_factor );


/* load a mesh from a file */
  fname = "../data/ramp.ani";
  ani3D_load_mesh( ani, fname );


/* allocate more memory than set by default */
  ani3D_set_max_points(   ani, 10000 );
  ani3D_set_max_elements( ani, 50000 );


/* example of accessing mesh points */
  nP = ani3D_number_points( ani ); 
  printf("\nNumber of points is%5d  (we print two)\n", nP);

  for( i=0; i<2; i++ ) {
     ani3D_get_point( ani, i, xy );

     printf("%4d %6.3f %6.3f %6.3f\n", i, xy[0], xy[1], xy[2]);
  }


/* example of accessing boundary edges */
  nF = ani3D_number_faces( ani ); 
  printf("\nNumber of boundary faces is%5d  (we print two)\n", nF);

  for( i=0; i<2; i++ ) {
     ani3D_get_face( ani, i, face, &label );

     printf("%4d  pts: %5d %5d %5d  label =%2d\n", i, face[0],face[1],face[2], label);
  }


/* example of accessing elements */
  nE = ani3D_number_elements( ani ); 
  printf("\nNumber of triangles is%5d  (we print two)\n", nE);
  
  for( i=0; i<2; i++ ) {
     ani3D_get_element( ani, i, tet, &label );
     printf("%4d  pts: %5d %5d %5d %5d  label =%2d\n", i, tet[0],tet[1],tet[2],tet[3], label);
  }


/* control parameter */
  Q = 0.4;
  ani3D_set_quality(   ani, Q );
  ani3D_get_norm(      ani, &Lp );
  ani3D_get_max_iters( ani, &max_iters );

  printf("\nFinal mesh quality be %5.2f  after at most %5d iterations\n", Q, max_iters);
  if( Lp == 0 ) printf("The mesh will be optimal in the maximum norm\n\n");
  else          printf("The mesh will be optimal in Lp norm, Lp = %5.2f\n\n", Lp);


/* create an artificial diagonal metric (1, 2, 3) */
  nP = ani3D_number_points( ani );
  metric = malloc( 6*nP * sizeof(double) );

  for( n=0,i=0; i<nP; i++ ) {
    metric[n++] = 1.0; 
    metric[n++] = 2.0; 
    metric[n++] = 3.0; 
    metric[n++] = 0.0; 
    metric[n++] = 0.0; 
    metric[n++] = 0.0; 
  }
  

/* adapt the mesh to metric */
  ani3D_nodal( ani, metric );


/* save the mesh in a file */
  fname = "save.ani";
  ani3D_save_mesh( ani, fname );


/* draw final mesh. The name must terminate with .ps */
  fname = "save.gmv";
  ani3D_draw_mesh( ani, fname );


/* kill the mesh */
  ani3D_KILL( ani );

  exit(0);
}


