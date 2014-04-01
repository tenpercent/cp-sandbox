#include <stdlib.h>

#include "ani3D.h"


int ani3D_INIT( ani3D* ani, f77i nEStar, f77r mem_factor )
{
  int mem;

  /* default parameters */
  (*ani).nP = 0;
  (*ani).nF = 0;
  (*ani).nE = 0;

  (*ani).nPv = 0;
  (*ani).nFv = 0;
  (*ani).nEv = 0;

  (*ani).MaxPV = 100;
  (*ani).MaxFV = 100;
  (*ani).MaxEV = 100;

  (*ani).flagAuto = 1;
  (*ani).status = 0;

  (*ani).MaxSkipE = 500;
  (*ani).MaxQItr = 15000;
  (*ani).Quality = 0.4;
  (*ani).Lp = 0;

  (*ani).mem_factor = mem_factor;

  (*ani).iPrint = 1;


  /* default mesh memory */
  (*ani).nEStar = nEStar;

  (*ani).MaxE = nEStar * mem_factor;
  (*ani).MaxF = nEStar * mem_factor / 4;
  (*ani).MaxP = nEStar * mem_factor / 2;
}


int ani3D_KILL( ani3D* ani )
{
  free(ani->XYP);
  free(ani->IPF);
  free(ani->IPE);

  free(ani->lbF);
  free(ani->lbE);

  free(ani->IPV);
  free(ani->IFV);
  free(ani->IEV);

  return 0;
}


/* C INTERFACE TO FORTRAN ROUTINES */
int ani3D_load_mesh( ani3D* ani, const char* fname )
{
  int mem, *IPP, *IFF;
 
  loadmani_header_(
        &ani->nP, &ani->nF, &ani->nE, 
        &ani->nPv, &ani->nFv, &ani->nEv, 
        fname);

  mem = ani->MaxP = C2Fmax(ani->nP, ani->MaxP);
  ani->XYP = malloc( 3*mem * sizeof(f77r) );

  mem = ani->MaxF = C2Fmax(ani->nF, ani->MaxF);
  ani->lbF = malloc(   mem * sizeof(f77i) );
  ani->IPF = malloc( 3*mem * sizeof(f77i) );

  mem = ani->MaxE = C2Fmax(ani->nE, ani->MaxE);
  ani->lbE = malloc(   mem * sizeof(f77i) );
  ani->IPE = malloc( 4*mem * sizeof(f77i) );
 
  mem = ani->MaxPV = C2Fmax(ani->nPv, ani->MaxPV);
  ani->IPV = malloc( mem * sizeof(f77i) );

  mem = ani->MaxFV = C2Fmax(ani->nFv, ani->MaxFV);
  ani->IFV = malloc( mem * sizeof(f77i) );

  mem = ani->MaxEV = C2Fmax(ani->nEv, ani->MaxEV);
  ani->IEV = malloc( mem * sizeof(f77i) );
 
  loadmani_(
/* M */ &ani->MaxP, &ani->MaxF, &ani->MaxE, &ani->nP, &ani->nF, &ani->nE,
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE,
        &ani->nPv, &ani->nFv, &ani->nEv, ani->IPV, ani->IFV, ani->IEV,
        IPP, IFF,
        fname);
}


int ani3D_save_mesh( ani3D* ani, const char* fname )
{
  int flagI = 0, *IPP, *IFF;

  savemani_(
/* M */ &ani->nP, &ani->nF, &ani->nE, 
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE,
        &ani->nPv, &ani->nFv, &ani->nEv, ani->IPV, ani->IFV, ani->IEV,
        &flagI, IPP, IFF,
        fname);
}


int ani3D_save_mesh_gmv( ani3D* ani, const char* fname )
{
  int mem, *iW;

  mem = ani->nE;
  iW = malloc( mem * sizeof(f77i) );

  savemgmv_(
/* M */ &ani->nP, &ani->nF, &ani->nE, 
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE,
        fname, iW);

  free(iW);
}


int ani3D_analytic( ani3D* ani, void* metric )
{
  int  mem;
  f77i iERR;

  /* working memory */
  mem =  7 * ani->MaxP + ani->nP 
      +  7 * ani->MaxF 
      + 21 * ani->MaxE + 13 * ani->nE + 10000;

  ani->MaxWi = mem;
  ani->iW    = malloc( mem * sizeof(f77i) );

  mem = 14 * ani->MaxP + 14 * ani->nP 
      + 19 * ani->MaxE + 10000; 

  ani->MaxWr = mem;
  ani->rW    = malloc( mem * sizeof(f77r) );

  /* run the code */
  mbaanalytic_(
/* M */ &ani->nP, &ani->MaxP, &ani->nF, &ani->MaxF, &ani->nE, &ani->MaxE, 
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE,
        &ani->nEStar, 
/* D */ &ani->nFv, &ani->nFv, &ani->nEv, ani->IPV, ani->IFV, ani->IEV,
        &ani->flagAuto, &ani->status,
/* Q */ &ani->MaxSkipE, &ani->MaxQItr,
        metric, &ani->Quality, &ani->rQuality, 
/* W */ &ani->MaxWr, &ani->MaxWi, ani->rW, ani->iW,
        &ani->iPrint, &iERR);


  /* free working memory */
  free(ani->iW);
  free(ani->rW);

  return iERR;
}


int ani3D_nodal( ani3D* ani, double* metric )
{
  int  i, mem;
  f77i iERR;

  /* working memory */
  mem =  7 * ani->MaxP + ani->nP 
      +  7 * ani->MaxF 
      + 21 * ani->MaxE + 13 * ani->nE + 10000;

  ani->MaxWi = mem;
  ani->iW    = malloc( mem * sizeof(f77i) );

  mem = 14 * ani->MaxP + 14 * ani->nP 
      + 19 * ani->MaxE + 10000; 

  ani->MaxWr = mem;
  ani->rW    = malloc( mem * sizeof(f77r) );


  /* copy solution */
  mem = 6 * ani->MaxP;
  ani->Metric = malloc( mem * sizeof(f77r) );

  for( i=0; i<6*ani->nP; i++ ) ani->Metric[i] = metric[i];


  /* run the code */
  mbanodal_(
/* M */ &ani->nP, &ani->MaxP, &ani->nF, &ani->MaxF, &ani->nE, &ani->MaxE, 
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE,
        &ani->nEStar, 
/* D */ &ani->nFv, &ani->nFv, &ani->nEv, ani->IPV, ani->IFV, ani->IEV,
        &ani->flagAuto, &ani->status,
/* Q */ &ani->MaxSkipE, &ani->MaxQItr,
        ani->Metric, &ani->Quality, &ani->rQuality, 
/* W */ &ani->MaxWr, &ani->MaxWi, ani->rW, ani->iW,
        &ani->iPrint, &iERR);


  /* free working memory */
  free(ani->Metric);

  free(ani->iW);
  free(ani->rW);

  return iERR;
}


int  ani3D_draw_mesh(ani3D* ani, char* fname )
{
  int mem;

  mem = ani->nE;
  ani->iW = calloc( mem, sizeof(f77i) );

  savemgmv_(
/* M */ &ani->nP, &ani->nF, &ani->nE, 
        ani->XYP, ani->IPF, ani->IPE, ani->lbF, ani->lbE, 
        fname, ani->iW);

  free(ani->iW);
}


/* ELEMENTAL ROUTINES */
inline int ani3D_number_points(   ani3D* ani ) { return ani->nP; }
inline int ani3D_number_faces(    ani3D* ani ) { return ani->nF; }
inline int ani3D_number_elements( ani3D* ani ) { return ani->nE; }

inline int C2Fmax(int a, int b) { return (a>b)? a : b; }


/* set maximal number of points, edges and elements */
void ani3D_set_max_points( ani3D* ani, int maxP ) 
{ 
  if( maxP > ani->MaxP ) {
  ani->MaxP = maxP; 
  ani->XYP = realloc( ani->XYP, 3*maxP * sizeof(f77r) );
  }
}

void ani3D_set_max_faces( ani3D* ani, int maxF ) 
{ 
  if( maxF > ani->MaxF ) {
  ani->MaxF = maxF; 
  ani->IPF  = realloc( ani->IPF, 3*maxF * sizeof(f77i) );
  }
}

void ani3D_set_max_elements( ani3D* ani, int maxE ) 
{ 
  if( maxE > ani->MaxE ) {
  ani->MaxE = maxE; 
  ani->lbE = realloc( ani->lbE,   maxE * sizeof(f77i) );
  ani->IPE = realloc( ani->IPE, 4*maxE * sizeof(f77i) );
  }
}


/* operations with points */
inline void ani3D_get_point( ani3D* ani, int i, double* xy )
{
  int k;
  for( k=0; k<3; k++ ) xy[k] = ani->XYP[3*i+k];   
}

void ani3D_set_point( ani3D* ani, int i, double* xy )
{
  int k;
  for( k=0; k<3; k++ ) ani->XYP[3*i+k] = xy[k];   
}

void ani3D_fix_point( ani3D* ani, int i )
{
  int k, nPv, flag, mem;

  flag = 0;
  nPv = ani->nPv;

  for( k=0; k<nPv; k++ ) 
    if( ani->IPV[k] = i ) { flag = 1; break; }

  if( !flag ) { 
    nPv++; 
    if( nPv == ani->MaxPV ) {
       mem = ani->MaxPV + 100;
       ani->IPV = realloc( ani->IPV, mem * sizeof(f77i) );
    }

    ani->IPV[nPv] = i; 
  };
}


/* operations with boundary edges */
void ani3D_get_face( ani3D* ani, int i, int* face, int* label )
{
  int k;
  for( k=0; k<3; k++ ) face[k] = ani->IPF[3*i+k];
  *label = ani->lbF[i];
}

void ani3D_set_face( ani3D* ani, int i, int* face, int label )
{
  int k;
  for( k=0; k<3; k++ ) ani->IPF[3*i+k] = face[k];
  ani->lbF[i] = label;
}

void ani3D_fix_face( ani3D* ani, int i )
{
  int k, nFv, flag, mem;

  flag = 0;
  nFv = ani->nFv;

  for( k=0; k<nFv; k++ ) 
    if( ani->IFV[k] = i ) { flag = 1; break; }

  if( !flag ) { 
    nFv++; 
    if( nFv == ani->MaxFV ) {
       mem = ani->MaxFV + 100;
       ani->IFV = realloc( ani->IFV, mem * sizeof(f77i) );
    }

    ani->IFV[nFv] = i; 
  };
}


/* operations with elements */
void ani3D_get_element( ani3D* ani, int i, int* tet, int* label )
{
  int k;
  for( k=0; k<4; k++ ) tet[k] = ani->IPE[4*i+k];
  *label = ani->lbE[i];
}

void ani3D_set_element( ani3D* ani, int i, int* tet, int label )
{
  int k;
  for( k=0; k<4; k++ ) ani->IPE[3*i+k] = tet[k];
  ani->lbE[i] = label;
}

void ani3D_fix_element( ani3D* ani, int i )
{
  int k, nEv, flag, mem;

  flag = 0;
  nEv = ani->nEv;

  for( k=0; k<nEv; k++ ) 
    if( ani->IEV[k] = i ) { flag = 1; break; }

  if( !flag ) { 
    nEv++; 
    if( nEv == ani->MaxEV ) {
       mem = ani->MaxEV + 100;
       ani->IEV = realloc( ani->IEV, mem * sizeof(f77i) );
    }

    ani->IEV[nEv] = i; 
  };
}


/* control parameers */
inline void ani3D_get_quality( ani3D* ani, double* Q ) { *Q = ani->rQuality; }
inline void ani3D_set_quality( ani3D* ani, double  Q ) { ani->Quality = Q; }

inline void ani3D_get_norm( ani3D* ani, double* Lp ) { *Lp = ani->Lp; }
inline void ani3D_set_norm( ani3D* ani, double  Lp ) { ani->Lp = Lp; }

inline void ani3D_get_status( ani3D* ani, int* status ) { *status = ani->status; }
inline void ani3D_set_status( ani3D* ani, int  status ) { ani->status = status; }

inline void ani3D_get_max_iters( ani3D* ani, int* max_iters ) { *max_iters = ani->MaxQItr; }
inline void ani3D_set_max_iters( ani3D* ani, int  max_iters ) { ani->MaxQItr = max_iters; }

inline void ani3D_get_max_basket( ani3D* ani, int* max_basket ) { *max_basket = ani->MaxSkipE; }
inline void ani3D_set_max_basket( ani3D* ani, int  max_basket ) { ani->MaxSkipE = max_basket; }

