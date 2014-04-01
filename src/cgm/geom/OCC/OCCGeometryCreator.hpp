//-------------------------------------------------------------------------
// Filename      : OCCGeometryCreator.hpp 
//
// Purpose       : For a list of facets, builds a topology of
//                 surfaces, curves and vertices 
//
// Special Notes : 
//
// Creator       : Steven J. Owen
//
// Date          : 4/29/2001
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef MESHGEOMETRYCREATOR_HPP
#define MESHGEOMETRYCREATOR_HPP

#include "DLIList.hpp"

class FacetEntity;
class FacetCurveMesh;
class FacetSurfaceMesh;
class FacetVolumeMesh;
class LoopSM;
class Body;
class Surface;
class FacetPointMesh;
class CubitPoint;
class CubitFacetEdge;
class CubitFacet;

class FacetGeometryCreator
{
public:
  FacetGeometryCreator( );
  FacetGeometryCreator(DLIList<FacetEntity*> &face_list,
                       DLIList<FacetEntity*> &edge_list,
                       DLIList<FacetEntity*> &point_list );
  //- Constructor

  ~FacetGeometryCreator();
  //- Destructor

#ifdef BOYD14
  CubitStatus create_geometry(CubitBoolean use_feature_angle,
                              double angle,
                              int interp_order,
                              CubitBoolean smooth_non_manifold,
                              CubitBoolean split_surfaces,
                              DLIList<Surface *> &surface_list);
    //- Creates the geometry infrastructure based on
    //- the given mesh data.
#endif
  void print_me();
   
private:

  DLIList<FacetEntity*> faceList;
  DLIList<FacetEntity*> edgeList;
  DLIList<FacetEntity*> pointList;
  DLIList<FacetSurfaceMesh*> facetSurfaceList;
  DLIList<FacetCurveMesh*> facetCurveList;
  DLIList<FacetPointMesh*> facetPointList;

  DLIList<FacetCurveMesh*> *hashCurveArray;
  int hashCurveSize;
  DLIList<FacetPointMesh*> *hashPointArray;
  int hashPointSize;

  void set_up_tool_datas( );
  void delete_tool_datas( );

#ifdef BOYD14
  int facet_dimension(FacetEntity *facet_ptr);
    //- returns the dimension of the facet entity.
#endif

  CubitStatus create_volume_boundaries( DLIList<FacetSurfaceMesh*> &facet_surface_sheets,
                                        CubitBoolean use_feature_angle,
                                        double min_dot,
                                        CubitBoolean split_surfaces );
    //- creates the correct facetsurfaces for all the element blocks.

  CubitStatus create_surface_boundaries( DLIList<FacetSurfaceMesh*> &facet_surface_list,
                                         DLIList<FacetCurveMesh*> &facet_curve_list,
                                         CubitBoolean use_feature_angle,
                                         double min_dot );
    //- creates the correct blockcurves for all the surfaces

  CubitStatus create_curve_boundaries( DLIList<FacetCurveMesh*> &facet_curve_list,
                                       DLIList<FacetPointMesh*> &facet_point_list );
    //- creates the correct blockpoints for all the curves

  CubitStatus classify_edge( FacetEntity *edge_ptr,
                             DLIList<FacetCurveMesh*> &facet_curve_list,
                             FacetSurfaceMesh *fsm_ptr );
    //- sorts a edge into its correct curve based on its associated
    //- blocks and sidesets.  Creates a new facet curve if necessary.

  CubitStatus classify_point(CubitPoint *point_ptr,
                             DLIList<FacetPointMesh*> &facet_point_list,
                             FacetCurveMesh *fcm_ptr ); 
    //- sorts a node into correct point based on its associated
    //- curve.  Creates a new facet point if necessary
   
  CubitStatus build_geometry( DLIList<FacetSurfaceMesh*> &facet_surface_list,
                              DLIList<FacetCurveMesh*> &facet_curve_list,
                              DLIList<FacetPointMesh*> &facet_point_list,
                              int interp_order,
                              CubitBoolean use_feature_angle,
                              double min_dot,
                              CubitBoolean smooth_non_manifold ,
                              CubitBoolean split_surfaces );
    //- build the CUBIT geometry based on the Facet entity class lists

  CubitStatus build_point_geometry( DLIList<FacetPointMesh*> &facet_point_list );
    //- From the facet point list, create geometric points for each

  CubitStatus build_curve_geometry( DLIList<FacetCurveMesh*> &facet_curve_list );
    //- From the facet curve list, create geometric curves for each

  CubitStatus build_loop_geometry( DLIList<FacetCurveMesh*> &facet_curve_list,
                                   FacetSurfaceMesh *fsm_ptr,
                                   DLIList<LoopSM*> &loop_list,
                                   int &ncurves );
    //- From the facet curve list of a surface, create geometric loops

  CubitStatus build_surface_geometry( DLIList<FacetSurfaceMesh*> &facet_surface_list,
                                      int interp_order, double min_dot );
    //- From the facet surface list, create geometric surface,
    //- loops and coedges for each surface in the list

  CubitStatus init_hash_curves( );
  void delete_hash_curves( );
  int get_curve_hash_key( DLIList<FacetSurfaceMesh*> *bsm_list_ptr );
    // functions for hashing curves - to speed up edge classification

  CubitStatus init_hash_points( );
  void delete_hash_points( );
  int get_point_hash_key( DLIList<FacetCurveMesh*> *bsm_list_ptr );
    // functions for hashing points - to speed up node classification

  CubitStatus clean_geometry( CubitBoolean smooth_non_manifold,
                            CubitBoolean split_surfaces,
                            CubitBoolean use_feature_angle,
                            double min_dot,
                            DLIList <FacetCurveMesh *> &facet_curve_list );

    //- fix the edge control points and the normals so they are conforming
    //- (or non-conforming) accross curves


};

#endif
