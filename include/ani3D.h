#ifndef   __ANI_GEN_H
#define   __ANI_GEN_H


#ifdef   __cplusplus
extern "C" {
#endif


typedef long int   f77i;
typedef double     f77r;
typedef char       f77c;
typedef long int   f77b;


typedef struct{
/* variables for main calls */
  f77i      nP,    nF,    nE;
  f77i      nPv,   nFv,   nEv;
  f77i      MaxP,  MaxF,  MaxE;
  f77i      MaxPV, MaxFV, MaxEV;

  f77r      *XYP;
  f77i      *IPF, *IPE;

  f77i      nEStar;

  f77i      *IPV, *IFV, *IEV, *lbF, *lbE;
  f77b      flagAuto;
  f77i      status;

  f77i      MaxSkipE, MaxQItr;
  f77r      *Metric, Quality, rQuality, Lp;

  f77i      MaxWr, MaxWi, *iW;
  f77r      *rW;
  f77i      iPrint;

/* other variables */
  double    mem_factor;
} ani3D;



/* C INTERFACES TO FORTRAN ROUTINES */
void loadmani_(
/* M */ f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77r*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77c*);

void savemani_(
/* M */ f77i*, f77i*, f77i*, f77r*, f77i*, f77i*,
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77i*, f77b*, f77i*, f77i*,
        f77c*);


void loads_(f77i*, f77r*, char*);
void saves_(f77i*, f77r*, char*);

void mbaanalytic_(
/* M */ f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, 
        f77r*, f77i*, f77i*, f77i*, f77i*,
        f77i*, 
/* D */ f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77b*, f77i*,
/* Q */ f77i*, f77i*,
        void*, f77r*, f77r*, 
/* W */ f77i*, f77i*, f77r*, f77i*,
        f77i*, f77i*);

void mbanodal_(
/* M */ f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, 
        f77r*, f77i*, f77i*, f77i*, f77i*,
        f77i*, 
/* D */ f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77b*, f77i*,
/* Q */ f77i*, f77i*,
        f77r*, f77r*, f77r*, 
/* W */ f77i*, f77i*, f77r*, f77i*,
        f77i*, f77i*);

void savemgmv_(f77i*, f77i*, f77i*, 
               f77r*, f77i*, f77i*, f77i*, f77i*, f77c*, f77i*);



/* BASIC ROUTINES */
int  ani3D_INIT( ani3D* ani, f77i nEStar, f77r mem_factor );
int  ani3D_KILL( ani3D* ani );

int  ani3D_load_mesh( ani3D* ani, char* fname );
int  ani3D_save_mesh( ani3D* ani, char* fname );

int  ani3D_analytic( ani3D* ani, void* metric );
int  ani3D_nodal(    ani3D* ani, double* metric );

int  ani3D_draw_mesh(ani3D* ani, char* fname );


/* ELEMENTAL ROUTINES */
int  ani3D_number_points(   ani3D* ani );
int  ani3D_number_edges(    ani3D* ani );
int  ani3D_number_elements( ani3D* ani );

void ani3D_set_max_points(  ani3D* ani, int maxP );
void ani3D_set_max_edges(   ani3D* ani, int maxF );
void ani3D_set_max_elements(ani3D* ani, int maxE );

void ani3D_get_point( ani3D* ani, int i, double* xy );
void ani3D_set_point( ani3D* ani, int i, double* xy );
void ani3D_fix_point( ani3D* ani, int i );

void ani3D_get_edge( ani3D* ani, int i, int* edge, int* icrv, int* label );
void ani3D_set_edge( ani3D* ani, int i, int* edge, int  icrv, int  label );
void ani3D_fix_edge( ani3D* ani, int i );

void ani3D_get_element( ani3D* ani, int i, int* tri, int* label );
void ani3D_set_element( ani3D* ani, int i, int* tri, int  label );
void ani3D_fix_element( ani3D* ani, int i );

void ani3D_get_quality( ani3D* ani, double* rQ );
void ani3D_set_quality( ani3D* ani, double   Q );

void ani3D_get_norm(    ani3D* ani, double* Lp );
void ani3D_set_norm(    ani3D* ani, double  Lp );

void ani3D_get_status(  ani3D* ani, int* status );
void ani3D_set_status(  ani3D* ani, int  status );

void ani3D_get_max_iters(  ani3D* ani, int* max_iters );
void ani3D_set_max_iters(  ani3D* ani, int  max_iters );

void ani3D_get_max_basket( ani3D* ani, int* max_basket );
void ani3D_set_max_basket( ani3D* ani, int  max_basket );


#ifdef   __cplusplus
           }
#endif

#endif

