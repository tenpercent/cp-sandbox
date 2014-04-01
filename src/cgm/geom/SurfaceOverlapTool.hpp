//- Class:       SurfaceOverlapTool
//- Description: Utilities to debug imprinting/merging problems
//- Owner:       Steve Storm
//- Created:     22 October 1999

#ifndef SurfaceOverlapTool_HPP
#define SurfaceOverlapTool_HPP

#include <map>
#include "GMem.hpp"
#include "AnalyticGeometryTool.hpp"
#include "CubitGeomConfigure.h"
#include "CurveOverlapFacet.hpp"
#include "SurfaceOverlapFacet.hpp"
#include "AbstractTree.hpp"

template <class X> class DLIList;
class RefFace;

class CUBIT_GEOM_EXPORT SurfaceOverlapTool
{

public :   

   CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list, 
                                          DLIList<RefEntity*> &faces_to_draw );
   CubitStatus find_overlapping_surfaces( DLIList<Body*> &body_list,
                                          DLIList<RefEntity*> &faces_to_draw );
   CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                          DLIList<RefFace*> &ref_face_list1,
                                          DLIList<RefFace*> &ref_face_list2,
                                          DLIList<RefEntity*> &faces_to_draw,
                                          CubitBoolean show_messages = CUBIT_FALSE);
   CubitStatus find_overlapping_surfaces( DLIList<Body*> &body_list,
                                          DLIList<RefFace*> &ref_face_list1,
                                          DLIList<RefFace*> &ref_face_list2,
                                          DLIList<RefEntity*> &faces_to_draw,
                                          CubitBoolean show_messages = CUBIT_FALSE);

   CubitStatus find_overlapping_surfaces( DLIList<BodySM*> &body_list,
                                          DLIList<Surface*> &surface_list1,
                                          DLIList<Surface*> &surface_list2 );

   // Searches for surfaces that overlap each other and are good
   // candidates for imprinting.  It can even find those with gaps
   // between them (according to settings).
   // First method is for interactive use.  They print messages and
   // draw the surfaces as it finds them, as dictated by settings.
   // The second method is for non-interactive use.  It simply returns the
   // surface pairs in lists.  Note that the same surface can occur 
   // multiple times in the lists.  The first list will contain sequential
   // surfaces - ie., if surface 1 overlaps with 2 other surfaces, list1
   // will contain surface 1, then surface 1, then the next surface, 
   // etc..  If show_messages is CUBIT_TRUE, you will see "Finding surface
   // overlap...", and groupResults, listPairs and displayPairs will be
   // used, otherwise they will be false (default).

   CubitStatus find_overlapping_curves( DLIList<Body*> &bodies,
                   std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map );
   CubitStatus find_overlapping_curves( DLIList<BodySM*> &bodies,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map );

   CubitBoolean check_overlap( RefEdge *edge1, RefEdge *edge2, 
               std::map<RefEdge*, DLIList<CurveOverlapFacet*>* > *facet_map );

   CubitBoolean check_overlap( Curve *curve1, Curve *curve2, 
               std::map<Curve*, DLIList<CurveOverlapFacet*>* > *facet_map );

   void list_settings();
     // List the settings used for surface overlap detection

   unsigned short get_facet_ang_tol();
   void set_facet_ang_tol( unsigned short angle_tol );


   //Initialize all settings in this class
   static void initialize_settings();

   static int get_facet_ang_tol_setting();
   static void set_facet_ang_tol_setting( int angle_tol );
   
   static double get_facet_abs_tol();
   static void set_facet_abs_tol( double absolute_tol );

   static double get_gap_min();
   static void set_gap_min( double val );

   static double get_gap_max();
   static void set_gap_max( double val );

   static double get_angle_min();
   static void set_angle_min( double val );

   static double get_angle_max();
   static void set_angle_max( double val );

   static double get_overlap_tolerance();
   static void set_overlap_tolerance( double val );

   static int get_normal_type(); //
   static void set_normal_type( int type ); //

   static CubitString get_normal_type_setting();
   static void set_normal_type_setting( CubitString type );

   static CubitBoolean get_group_results(); //
   static void set_group_results( CubitBoolean setting ); //

   static CubitBoolean get_list_pairs(); //
   static void set_list_pairs( CubitBoolean setting ); //

   static CubitBoolean get_display_pairs(); //
   static void set_display_pairs( CubitBoolean setting ); //

   static int get_group_results_setting();
   static void set_group_results_setting( int setting );

   static int get_list_pairs_setting();
   static  void set_list_pairs_setting( int setting );

   static  int get_display_pairs_setting();
   static  void set_display_pairs_setting( int setting );

   static CubitBoolean get_imprint();
   static void set_imprint( CubitBoolean setting );

   static CubitBoolean get_check_within_bodies();
   static void set_check_within_bodies( CubitBoolean setting );

   static CubitBoolean get_check_across_bodies();
   static void set_check_across_bodies( CubitBoolean setting );

   static SurfaceOverlapTool* instance();
     // Returns a static pointer to unique instance of this class.
   
   ~SurfaceOverlapTool();
     //- Destructor.

protected :
   
private :

  CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                         DLIList<RefFace*> &ref_face_list1,
                                         DLIList<RefFace*> &ref_face_list2,
                                         DLIList<RefEntity*> &pair_list,
                                         CubitBoolean list_pairs,
                                         int prog_step );
  // Private workhorse function for finding the overlapping surfaces. Here 
  // pair_list is a list of all the surfaces that were found that overlap.
  // If prog_step == -1 it does not step the progress bar. 

   CubitBoolean check_overlap( RefFace *ref_face_ptr1,
                               RefFace *ref_face_ptr2,
                               CubitBoolean abort );
   CubitBoolean check_overlap( Surface *surface1,
                               Surface *surface2,
              std::map<Surface*, DLIList<SurfaceOverlapFacet*>* > *facet_map,
              std::map<Surface*, double > *area_map,
              std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* > *a_tree_map ); 

                               
   //CubitStatus draw_facets( GMem* gMem, int color = 2 );
   //void draw_triangle( Triangle3 &T, int color = 2 );

   CubitStatus imprint(DLIList<RefFace*> &ref_face_list1,
                       DLIList<RefFace*> &ref_face_list2);
   //- Imprints the bodies containing the lists of surfaces.  Also
   //- uses the edges imprinting rather than the normal imprint.

   int num_descendants(Body *body);
    //- Counts the total number of descendants of the body.
    //- Only counts Ref* entities.
  
   CubitStatus copy_edges_in_list(DLIList<RefEdge*> &old_list,
                                  DLIList<RefEdge*> &new_list );
    //- copy's each of the edges in the old list, and puts the
    //- stand alone edges in the new list.
    //- NOTE: THESE ARE NEWLY CREATED EDGES AND SHOULD BE DELETED
    //-       BY THE CALLING FUNCTION.

   CubitStatus delete_edges_in_list(DLIList<RefEdge*> &edge_list);
    //- Deletes the RefEdges in the list via GeometryQueryTool.

   static SurfaceOverlapTool* instance_;
     // Static pointer to unique instance of this class
   
   SurfaceOverlapTool();
     //- Constructor for the SurfaceOverlapTool object

   static  double facetAbsTol;
   static  unsigned short facetAngTol;
   static  double gapMin;
   static  double gapMax;
   static  double angleMin;
   static  double angleMax;
   static  int normalType; // 1=any, 2=opposite, 3=same
   static  CubitBoolean groupResults;
   static  CubitBoolean listPairs;
   static  CubitBoolean displayPairs;
   static  CubitBoolean imprintResults;
   static  double overlapTolerance;
   static  CubitBoolean checkWithinBodies;
   static  CubitBoolean checkAcrossBodies;
};

inline
void SurfaceOverlapTool::set_facet_ang_tol( unsigned short val )
{facetAngTol=val;}

inline
double SurfaceOverlapTool::get_facet_abs_tol()
{return facetAbsTol;}
inline
void SurfaceOverlapTool::set_facet_abs_tol( double val )
{facetAbsTol=val;}

inline
double SurfaceOverlapTool::get_gap_min()
{return gapMin;}
inline
void SurfaceOverlapTool::set_gap_min( double val )
{gapMin=val;}

inline
double SurfaceOverlapTool::get_gap_max()
{return gapMax;}
inline
void SurfaceOverlapTool::set_gap_max( double val )
{gapMax=val;}

inline
double SurfaceOverlapTool::get_angle_min()
{return angleMin;}
inline
void SurfaceOverlapTool::set_angle_min( double val )
{angleMin=val;}

inline
double SurfaceOverlapTool::get_angle_max()
{return angleMax;}
inline
void SurfaceOverlapTool::set_angle_max( double val )
{angleMax=val;}

inline
double SurfaceOverlapTool::get_overlap_tolerance()
{return overlapTolerance;}
inline
void SurfaceOverlapTool::set_overlap_tolerance( double val )
{overlapTolerance=val;}

#ifdef CAT
inline
int SurfaceOverlapTool::get_normal_type()
{return normalType;}
#endif
inline
void SurfaceOverlapTool::set_normal_type( int type )
{normalType=type;}

#ifdef CAT
inline
CubitBoolean SurfaceOverlapTool::get_group_results()
{return groupResults;}
#endif
inline
void SurfaceOverlapTool::set_group_results( CubitBoolean setting )
{groupResults=setting;}

#ifdef CAT
inline
CubitBoolean SurfaceOverlapTool::get_list_pairs()
{return listPairs;}
#endif
inline
void SurfaceOverlapTool::set_list_pairs( CubitBoolean setting )
{listPairs=setting;}

#ifdef CAT
inline
CubitBoolean SurfaceOverlapTool::get_display_pairs()
{return displayPairs;}
#endif
inline
void SurfaceOverlapTool::set_display_pairs( CubitBoolean setting )
{displayPairs=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_imprint()
{return imprintResults;}
inline
void SurfaceOverlapTool::set_imprint( CubitBoolean setting )
{imprintResults=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_check_within_bodies()
{return checkWithinBodies;}
inline
void SurfaceOverlapTool::set_check_within_bodies( CubitBoolean setting )
{checkWithinBodies=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_check_across_bodies()
{return checkAcrossBodies;}
inline
void SurfaceOverlapTool::set_check_across_bodies( CubitBoolean setting )
{checkAcrossBodies=setting;}  

#endif



