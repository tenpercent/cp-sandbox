//-Class: TDSplitSurface.hpp

#ifndef TD_SPLIT_SURFACE_HPP
#define TD_SPLIT_SURFACE_HPP

#include "ToolData.hpp"
#include "CastTo.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "Cubit2DPoint.hpp"
#include "CubitGeomConfigure.h"

class RefFace;
class RefEdge;
class RefVertex;
class CoEdge;

//================================================================================
// Description: A trivial class to hold two parameter values so that they 
//              can be stored in a DLIList.  They are the min and max parameter
//              space along the composite curve. (Split Surface Param)
// Author     : Steve Storm
// Date       : 2/3/2004
//================================================================================
class CUBIT_GEOM_EXPORT SSParam
{
public:

  SSParam( double min, double max );
  ~SSParam();

  double umin();
  double umax();

private:

  double uMin;
  double uMax;
};

inline double SSParam::umin()
{ return uMin; }

inline double SSParam::umax()
{ return uMax; }

//================================================================================
// Description: This class (Split Surface Side) holds a chain of curves on one 
//              side of the surface.  It is needed to handle queries using a 
//              composite curve concept.
// Author     : Steve Storm
// Date       : 2/3/2004
//================================================================================
class CUBIT_GEOM_EXPORT SSSide
{
public:

  SSSide( RefFace *ref_face_ptr, DLIList<CoEdge*> &co_edges,
          const CubitVector *collapsed_loc_ptr = NULL );
  ~SSSide();

  CubitStatus position_from_u( double u_value,
                               CubitVector& output_position );

  CubitStatus u_from_position ( const CubitVector& input_position,
                                double &u );
  CubitStatus u_from_position ( const CubitVector& input_position,
                                CoEdge *co_edge_ptr, SSParam *param,
                                double &u );

  CubitBoolean is_vertex_on( RefVertex *ref_vertex_ptr );
  //- Determines if the given vertex lies on this side spatially.

  double length();
  //- Return curve length of this side

  CubitStatus build_param_list_from_facets( double tolerance = 0.1 );
  //- This function builds paramList using the graphics facets.  The
  //- coordList will be filled during syncronization - syncronize_lists.
  //- The tolerance is used by the geometry engine as the allowable deviation 
  //- of the facets from the actual curve.

  CubitStatus build_param_list( double fraction, double distance, int num_segs,
                                DLIList<RefVertex*> &through_vertex_list );
  //- This function fills the paramList and coordList with the 
  //- appropriate number of points for the A and C sides of the surface.
  //- This is typically one point at the 50% location; however it could 
  //- be at a specified fraction value, at a vertex location or multiple 
  //- points evenly spaced.

  CubitStatus syncronize_lists( SSSide *other_side,
                                double param_tol=.1 );
  //- Update lists in "this" and other_side.  Parameter values from
  //- both lists will be merged so that paramList will be identical
  //- on both SSSides.  This also populates coordList.  Use this 
  //- function to syncronize the B-D sides.

  DLIList<CoEdge*> *co_edges();
  //- Returns a pointer to this side's CoEdge list

  int coord_list_size();
  //- Return size of coordList
  void coord_list_reset();
  //- Set coordList to start of list
  void coord_list_last();
  //- Set coordList to end of list
  CubitVector *coord_list_get();
  //- Return value of current item in coordList
  CubitVector *coord_list_get_and_step();
  //- Get and step on the coordList
  CubitVector *coord_list_get_and_back();
  //- Get and back on the coordList

  int param_list_size();
  //- Return size of paramList
  void param_list_reset();
  //- Set paramList to start of list
  void param_list_last();
  //- Set paramList to end of list
  double param_list_get();
  //- Return value of current item in paramList
  double param_list_get_and_step();
  //- Get and step on the paramList
  double param_list_get_and_back();
  //- Get and back on the paramList

  CubitBoolean is_collapsed();
  //- Return whether this side is collapsed or not

private:

  RefFace *refFacePtr;
  DLIList<CoEdge*> coEdgeChain;
  double paramHigh; // Parameter space of chain goes from 0.0 to curve length
  DLIList<SSParam*> coEdgeParamList; // Holds low and high parameter on each CoEdge
  
  DLIList<CubitVector*> coordList; // Coordinates for determining split locations
  DLIList<double> paramList; // Parameters for determining split locations
  // Note: sideA and sideC just contain the interior locs, while sideB and sideD 
  //       contain interior plus end locs built from curve facetting

  CubitBoolean isCollapsed; // CUBIT_TRUE if edge is collapsed (just a point)
};

inline double SSSide::length()
{ return paramHigh; }

inline DLIList<CoEdge*> *SSSide::co_edges()
{ return &coEdgeChain; }

inline int SSSide::coord_list_size()
{ return coordList.size(); }
inline void SSSide::coord_list_reset()
{ coordList.reset(); }
inline void SSSide::coord_list_last()
{ coordList.last(); }
inline CubitVector *SSSide::coord_list_get()
{ return coordList.get(); }
inline CubitVector *SSSide::coord_list_get_and_step()
{ return coordList.get_and_step(); }
inline CubitVector *SSSide::coord_list_get_and_back()
{ return coordList.get_and_back(); }

inline int SSSide::param_list_size()
{ return paramList.size(); }
inline void SSSide::param_list_reset()
{ paramList.reset(); }
inline void SSSide::param_list_last()
{ paramList.last(); }
inline double SSSide::param_list_get()
{ return paramList.get(); }
inline double SSSide::param_list_get_and_step()
{ return paramList.get_and_step(); }
inline double SSSide::param_list_get_and_back()
{ return paramList.get_and_back(); }

inline CubitBoolean SSSide::is_collapsed()
{ return isCollapsed;}

//================================================================================
// Description: This class holds data on the surface defining the sides of
//              the surface to be split.
// Author     : Steve Storm
// Date       : 2/3/2004
//================================================================================
class CUBIT_GEOM_EXPORT TDSplitSurface : public ToolData
{
public:

  TDSplitSurface( int vertex_type );
  //- Constructor if placed on a CoEdge

  TDSplitSurface( RefFace *ref_face_ptr );
  //- Constructor if placed on a RefFace

  ~TDSplitSurface();
  //- Destructor

  RefFace *ref_face_ptr();
  //- Return this TDSplitSurface's RefFace ptr
  
  void add_type( int type );
  int get_type();
  //- For CoEdges

  CubitStatus add_coedges( DLIList<CoEdge*> &co_edge_list,
                           int side_interval[] );
  CubitStatus add_a_coedges( DLIList<CoEdge*> &a_coedges,
                             RefVertex *start_vertex_ptr = NULL );
  CubitStatus add_b_coedges( DLIList<CoEdge*> &b_coedges, 
                             RefVertex *start_vertex_ptr = NULL );
  CubitStatus add_c_coedges( DLIList<CoEdge*> &c_coedges,
                             RefVertex *start_vertex_ptr = NULL );
  CubitStatus add_d_coedges( DLIList<CoEdge*> &d_coedges, 
                             RefVertex *start_vertex_ptr = NULL );
  //- Add the curves into SSSide.  Note any side could be collapsed because of
  //- a triangle, if that is true a vertex needs to be provided to give the 
  //- coordinate of the side.

  CubitStatus tessellate_sides( double tol, double fraction, double distance,
                                int num_segs,
                                DLIList<RefVertex*> &through_vertex_list );
  //- Tessellates curves on the sides to find boundary coordinates in preparation
  //- for the mapping algirithm.  To give a smooth mapping, the tessellations
  //- are matched across the sides (ie., all locations on either side are at the
  //- same parameter value along the side).  The through_vertex_list if optional -
  //- if specified, the split will pass through vertices in the through_vertex_list 
  //- that are on sides A or C spatially.
  
  RefVertex *start_vertex( CoEdge *co_edge_ptr );
  //- Get the starting vertex of the given CoEdge.

  DLIList<CoEdge*> *get_a_coedges();
  DLIList<CoEdge*> *get_b_coedges();
  DLIList<CoEdge*> *get_c_coedges();
  DLIList<CoEdge*> *get_d_coedges();
  //- Just return a pointer to the curves

  double length_a();
  double length_b();
  double length_c();
  double length_d();
  //- Return the length of the given side

  int coord_list_size_a();
  void coord_list_reset_a();
  void coord_list_last_a();
  CubitVector *coord_list_get_a();
  CubitVector *coord_list_get_and_step_a();
  CubitVector *coord_list_get_and_back_a();
  //- Used to traverse the coordinates on side A

  int coord_list_size_b();
  void coord_list_reset_b();
  void coord_list_last_b();
  CubitVector *coord_list_get_b();
  CubitVector *coord_list_get_and_step_b();
  CubitVector *coord_list_get_and_back_b();
  //- Used to traverse the coordinates on side B

  int coord_list_size_c();
  void coord_list_reset_c();
  void coord_list_last_c();
  CubitVector *coord_list_get_c();
  CubitVector *coord_list_get_and_step_c();
  CubitVector *coord_list_get_and_back_c();
  //- Used to traverse the coordinates on side C

  int coord_list_size_d();
  void coord_list_reset_d();
  void coord_list_last_d();
  CubitVector *coord_list_get_d();
  CubitVector *coord_list_get_and_step_d();
  CubitVector *coord_list_get_and_back_d();
  //- Used to traverse the coordinates on side D

  int param_list_size_a();
  void param_list_reset_a();
  void param_list_last_a();
  double param_list_get_a();
  double param_list_get_and_step_a();
  double param_list_get_and_back_a();
  //- Used to traverse the parameters on side A

  int param_list_size_b();
  void param_list_reset_b();
  void param_list_last_b();
  double param_list_get_b();
  double param_list_get_and_step_b();
  double param_list_get_and_back_b();
  //- Used to traverse the parameters on side B

  int param_list_size_c();
  void param_list_reset_c();
  void param_list_last_c();
  double param_list_get_c();
  double param_list_get_and_step_c();
  double param_list_get_and_back_c();
  //- Used to traverse the parameters on side C

  int param_list_size_d();
  void param_list_reset_d();
  void param_list_last_d();
  double param_list_get_d();
  double param_list_get_and_step_d();
  double param_list_get_and_back_d();
  //- Used to traverse the parameters on side D

  CubitBoolean is_a_collapsed();
  CubitBoolean is_b_collapsed();
  CubitBoolean is_c_collapsed();
  CubitBoolean is_d_collapsed();
  //- Return whether sides are collapsed or not

  static int is_split_surface(const ToolData* td)
     {return (CAST_TO(td, const TDSplitSurface) != NULL);}

private:

  RefFace *refFacePtr;

  SSSide *sideA;
  SSSide *sideB;
  SSSide *sideC;
  SSSide *sideD;

  // For vertices
  int vertexType;
};

inline RefFace *TDSplitSurface::ref_face_ptr()
{ return refFacePtr; }

inline int TDSplitSurface::get_type()
{ return vertexType; }

inline double TDSplitSurface::length_a()
{ return sideA->length(); }
inline double TDSplitSurface::length_b()
{ return sideB->length(); }
inline double TDSplitSurface::length_c()
{ return sideC->length(); }
inline double TDSplitSurface::length_d()
{ return sideD->length(); }

inline int TDSplitSurface::coord_list_size_a()
{ return sideA->coord_list_size(); }
inline void TDSplitSurface::coord_list_reset_a() 
{ sideA->coord_list_reset(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_step_a()
{ return sideA->coord_list_get_and_step(); }

inline int TDSplitSurface::coord_list_size_b()
{ return sideB->coord_list_size(); }
inline void TDSplitSurface::coord_list_reset_b() 
{ sideB->coord_list_reset(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_step_b()
{ return sideB->coord_list_get_and_step(); }

inline int TDSplitSurface::coord_list_size_c()
{ return sideC->coord_list_size(); }
inline void TDSplitSurface::coord_list_reset_c() 
{ sideC->coord_list_reset(); }
inline void TDSplitSurface::coord_list_last_c()
{ sideC->coord_list_last(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_step_c()
{ return sideC->coord_list_get_and_step(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_back_c()
{ return sideC->coord_list_get_and_back(); }

inline int TDSplitSurface::coord_list_size_d()
{ return sideD->coord_list_size(); }
inline void TDSplitSurface::coord_list_reset_d() 
{ sideD->coord_list_reset(); }
inline void TDSplitSurface::coord_list_last_d()
{ sideD->coord_list_last(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_step_d()
{ return sideD->coord_list_get_and_step(); }
inline CubitVector *TDSplitSurface::coord_list_get_and_back_d()
{ return sideD->coord_list_get_and_back(); }

inline void TDSplitSurface::param_list_reset_a() 
{ sideA->param_list_reset(); }
inline void TDSplitSurface::param_list_last_a()
{ sideA->param_list_last(); }
inline double TDSplitSurface::param_list_get_and_step_a()
{ return sideA->param_list_get_and_step(); }

inline void TDSplitSurface::param_list_reset_b() 
{ sideB->param_list_reset(); }
inline double TDSplitSurface::param_list_get_and_step_b()
{ return sideB->param_list_get_and_step(); }

inline void TDSplitSurface::param_list_last_c()
{ sideC->param_list_last(); }
inline double TDSplitSurface::param_list_get_and_back_c()
{ return sideC->param_list_get_and_back(); }

inline void TDSplitSurface::param_list_last_d()
{ sideD->param_list_last(); }
inline double TDSplitSurface::param_list_get_and_back_d()
{ return sideD->param_list_get_and_back(); }

inline CubitBoolean TDSplitSurface::is_a_collapsed()
{ return sideA->is_collapsed(); }

inline CubitBoolean TDSplitSurface::is_b_collapsed()
{ return sideB->is_collapsed(); }

inline CubitBoolean TDSplitSurface::is_c_collapsed()
{ return sideC->is_collapsed(); }

inline CubitBoolean TDSplitSurface::is_d_collapsed()
{ return sideD->is_collapsed(); }

#endif // TD_SPLIT_SURFACE_HPP
