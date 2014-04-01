//-----------------------------------------------------------------------------
// Filename      : GeometryHealerTool.hpp
//
// Purpose       : Define the healer interface for all solid modeling engines.
//
// Special Notes : This is a singleton pattern class for the healer functions.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 7/21/03
//
// Owner         : 
//-----------------------------------------------------------------------------

#ifndef _GEOMETRYHEALERTOOL_HPP
#define _GEOMETRYHEALERTOOL_HPP

// *** BEGIN INCLUDES *** //

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

// *** END INCLUDES *** //

// *** BEGIN FORWARD DECLARATIONS *** //

class GeometryHealerEngine;

class TopologyBridge;
class TopologyEntity;

class RefEntity;
class RefFace;
class RefEdge;

class Body;

// *** END FORWARD DECLARATIONS *** //

class CUBIT_GEOM_EXPORT GeometryHealerTool
{
   public:

      // singleton pattern class instance interface
      static GeometryHealerTool* instance( GeometryHealerEngine* GHEPtr = NULL );

      // destructor
      ~GeometryHealerTool();

      // *** BEGIN ENGINE OPERATIONS *** //

      // set the default healer engine
//      void set_default_engine( GeometryHealerEngine *ghe_ptr );

      // add a healer engine to the list
      void add_ghe( GeometryHealerEngine *ghe_ptr );

      // remove a healer engine from the list
//      CubitStatus remove_ghe( GeometryHealerEngine *ghe_ptr );

      // get the list of healer engines
//      void get_ghe_list( DLIList<GeometryHealerEngine*> &ghe_list );

   private:
      // returns the healer engine of an entity
      GeometryHealerEngine* get_engine( TopologyBridge *tb_ptr ) const;
      GeometryHealerEngine* get_engine( TopologyEntity *te_ptr ) const;

      // determines if entities belong to the same engine
      CubitBoolean same_healer_engine( DLIList<RefEntity*> &ref_entity_list,
						                     CubitBoolean check_children ) const;
      CubitBoolean same_healer_engine( DLIList<TopologyEntity*> &topo_list ) const;

      // *** END ENGINE OPERATIONS *** //

      // *** BEGIN HEALER FUNCTIONS *** //
   public:
      CubitStatus auto_heal_bodies( DLIList<Body*> &body_list, 
                                    DLIList<Body*> &new_body_list,
                                    DLIList<TopologyEntity*> &bad_geometry,
                                    CubitBoolean rebuild = CUBIT_FALSE, 
                                    CubitBoolean keep_old = CUBIT_FALSE,
                                    CubitBoolean make_tolerant = CUBIT_FALSE,
                                    FILE* logfile_ptr = NULL );
      // Uses the autohealer to heal the given body list.  The rebuild option
      // can be used for more rigorous healing, where each surface is pulled off,
      // healed, and then stitched back into a body.

      CubitStatus heal_bodies( DLIList<Body*> &body_list, 
                               DLIList<Body*> &new_body_list,
                               DLIList<TopologyEntity*> &bad_geometry,
                               CubitBoolean rebuild = CUBIT_FALSE,
                               CubitBoolean keep_old = CUBIT_FALSE,
                               CubitBoolean make_tolerant = CUBIT_FALSE,
                               FILE* logfile_ptr = NULL );
      // Heals the input bodies, using either autoheal or incremental,
      // depending on the switches that are setup.  Writes to the logfile
      // if it's been opened.

      CubitStatus analyze_badgeom( DLIList<Body*> &body_list, 
                                   DLIList<TopologyEntity*> &bad_geometry,
                                   FILE* logfile = NULL );
      // Uses the healing husk to find bad geometry and give feedback to 
      // the user.  A lot more could be done here - this is just a quick 
      // overview for the user (the user will have no idea why the identified 
      // geometry is bad).  A logfile can give a little more information.

      CubitStatus get_badgeom( DLIList<Body*> &body_list, 
                               DLIList<TopologyEntity*> &bad_geometry );
      // Shows the bad geometry.  The geometry must have been analyzed first.
      // If body_list is empty, shows only for those bodies that have been
      // analyzed.

      CubitStatus get_tcurves( DLIList<Body*> &body_list,
                               DLIList<RefEdge*> &t_curves );
      // Get tolerant curves.  If the body list is empty, all tolerant 
      // curves are retrieved.

      CubitStatus heal_incremental( DLIList<Body*> &body_list, 
                                    DLIList<Body*> &new_bodies,
                                    DLIList<TopologyEntity*> &bad_geometry,
                                    double simplify_tol = CUBIT_DBL_MAX,
                                    double stitch_min_tol = CUBIT_DBL_MAX,
                                    double stitch_max_tol = CUBIT_DBL_MAX,
                                    double geombuild_tol = CUBIT_DBL_MAX,
                                    double analytic_tol = CUBIT_DBL_MAX,
                                    double isospline_tol = CUBIT_DBL_MAX,
                                    double reblend_classify_tol = CUBIT_DBL_MAX,
                                    double reblend_tol = CUBIT_DBL_MAX,
                                    CubitBoolean keep_old = CUBIT_FALSE,
                                    CubitBoolean make_tolerant = CUBIT_FALSE,
                                    FILE* logfile_ptr = NULL);
      // Uses the healing husk to perform one or more healing operations.
      // Which operations are performed is determined by switches in the
      // healer.  The user can set default tolerances (separate from this 
      // function, override the defaults with this function, or just use the 
      // tolerances calculated by the healer (preferred).
      //
      // The possible healing steps are:
      //    1) preprocess - removes zero-length edges, sliver faces, duplicate 
      //                    vertices. This is the first step which is always 
      //                    done (user shouldn't suppress this).  Tolerance is resabs.
      //    2) simplify - simplify NURBS into analytic.   
      //                     Default simplify_tol = .0001
      //    3) stitch - stitch geometry together.  Iterative from min to max 
      //                tolerance.
      //                   Default min tol = 10e-5
      //                   Default max tol = 1
      //    4) geombuild - geometry-related healing.  
      //                      Default geombuild_tol = .01
      //                      Default analytic_tangency_tol = .01
      //                      Default isolspline_solver_tol = .01      
      //          Individual geombuild steps can be (instead of doing all):
      //             analytic - performs all of the stages of the analytic solver 
      //                        subphase of the geometry building phase. The analytic 
      //                        solver subphase attempts to heal all edges and 
      //                        vertices shared by analytic surfaces. 
      //             isospline - performs all of the stages of the isospline solver
      //                         subphase of the geometry building phase. The 
      //                         isospline solver attempts to heal all edges shared 
      //                         by tangential isoparametric surfaces (e.g., the 
      //                         intersection curve is an isoparametric curve of 
      //                         both splines in the intersection). 
      //             reblend - future option
      //             sharpedge - performs all of the stages of the sharp edge solver 
      //                         subphase of the geometry building phase. The sharp 
      //                         edge solver attempts to heal all edges and vertices 
      //                         that are shared by surfaces that intersect sharply. 
      //                         This includes nontangential surface junctions. 
      //             genericspline - performs all of the stages of the generic spline 
      //                             solver subphase of the geometry building phase. 
      //                             The generic spline solver attempts to heal 
      //                             generic tangential spline junctions, (e.g., the 
      //                             intersection curve is not an isoparametric curve 
      //                             of both splines in the intersection). 
      //             wrapup - handles remaining pcurves, wraps up
      //                      geometry buiilding phase (user shouldn't suppress)
      //    5) postprocess - correction of neg-area faces, duplicate vertices, edge groups;
      //                     last step - always done (user shouldn't suppress)

      void list_incremental();
      // Lists the current settings for incremental healing

      void list_tolerances( DLIList<Body*> &body_list );
      // Lists out the tolerances of each body, separately

      double get_default_simplify_tol();
      void set_default_simplify_tol( double tol );
      double get_default_stitch_min_tol();
      void set_default_stitch_min_tol( double tol );
      double get_default_stitch_max_tol();
      void set_default_stitch_max_tol( double tol );
      double get_default_geombuild_tol();
      void set_default_geombuild_tol( double tol );
      double get_default_analytic_tol();
      void set_default_analytic_tol( double tol );
      double get_default_isospline_tol();
      void set_default_isospline_tol( double tol );
      double get_default_reblend_classify_tol();
      void set_default_reblend_classify_tol( double tol );
      double get_default_reblend_tol();
      void set_default_reblend_tol( double tol );
      void reset_default_tolerances();
      void list_default_tolerances();
      // Functions to set the default tolerances used.  The healer calculates
      // the default tolerance per body.  These functions allow the user to override
      // these defaults for all bodies healed.  In incremental healing, the user can 
      // override these defaults by sending-in tolerances.  For autoheal, these
      // defaults are used if they are set, otherwise the healer calculates
      // intelligent defaults.

      void clean_attributes( DLIList<Body*>& body_list );
      // Cleans healer attributes from the bodies.  These can be left if the 
      // CleanAtt setting was used when doing analysis or healing.

      CubitBoolean get_cleanatt_flg();
      void set_cleanatt_flg( CubitBoolean flg );
      // Get/set flags which determine whether to clean attributes after
      // analysis and healing.

      int get_show_method(); // 0-none, 1-highlight, 2-draw
      void set_show_method( int method );
      // Method determines how bad geometry is displayed

      CubitBoolean get_show_summary_flg();
      void set_show_summary_flg( CubitBoolean flg );
      // Flag to determine whether to list a summary when showing bad geometry

      CubitBoolean get_show_details_flg();
      void set_show_details_flg( CubitBoolean flg );
      // Flag to determine whether to list details when showing bad geometry

      CubitBoolean get_show_on_heal_flg();
      void set_show_on_heal_flg( CubitBoolean flg );
      // Flag to determine whether to show bad geometry after healing

      CubitBoolean get_check_vol_on_heal_flg();
      void set_check_vol_on_heal_flg( CubitBoolean flg );
      double get_vol_on_heal_limit();
      void set_vol_on_heal_limit( double limit );
      // Allow for checking of short curves after healing

      CubitBoolean get_check_surf_on_heal_flg();
      void set_check_surf_on_heal_flg( CubitBoolean flg );
      double get_surf_on_heal_limit();
      void set_surf_on_heal_limit( double limit );
      // Allow for checking of small surfaces after healing

      CubitBoolean get_check_curve_on_heal_flg();
      void set_check_curve_on_heal_flg( CubitBoolean flg );
      double get_curve_on_heal_limit();
      void set_curve_on_heal_limit( double limit );
      // Allow for checking of short curves after healing

      CubitBoolean get_show_bad_vertices_flg();
      void set_show_bad_vertices_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_curves_flg();
      void set_show_bad_curves_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_coedges_flg();
      void set_show_bad_coedges_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_loops_flg();
      void set_show_bad_loops_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_surfaces_flg();
      void set_show_bad_surfaces_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_shells_flg();
      void set_show_bad_shells_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_volumes_flg();
      void set_show_bad_volumes_flg( CubitBoolean flg );
      CubitBoolean get_show_bad_bodies_flg();
      void set_show_bad_bodies_flg( CubitBoolean flg );
      // Get/set for flags for what to show during analysis/show

      void list_onshow_flgs();
      // Function to list to user what the current onshow flags are

      CubitBoolean get_inc_preprocess_flg();
      void set_inc_preprocess_flg( CubitBoolean flg );
      CubitBoolean get_inc_simplify_flg();
      void set_inc_simplify_flg( CubitBoolean flg );
      CubitBoolean get_inc_stitch_flg();
      void set_inc_stitch_flg( CubitBoolean flg );
      CubitBoolean get_inc_geombuild_flg();
      void set_inc_geombuild_flg( CubitBoolean flg );
      CubitBoolean get_inc_analytic_flg();
      void set_inc_analytic_flg( CubitBoolean flg );
      CubitBoolean get_inc_isospline_flg();
      void set_inc_isospline_flg( CubitBoolean flg );
      CubitBoolean get_inc_reblend_flg();
      void set_inc_reblend_flg( CubitBoolean flg );
      CubitBoolean get_inc_sharpedge_flg();
      void set_inc_sharpedge_flg( CubitBoolean flg );
      CubitBoolean get_inc_genericspline_flg();
      void set_inc_genericspline_flg( CubitBoolean flg );
      CubitBoolean get_inc_wrapup_flg();
      void set_inc_wrapup_flg( CubitBoolean flg );
      CubitBoolean get_inc_postprocess_flg();
      void set_inc_postprocess_flg( CubitBoolean flg );
      // Functions for controlling incremental healing

      CubitStatus force_simplify_to_plane( DLIList<RefFace*> &ref_face_list, DLIList<Body*>& new_body_list, 
                                           CubitBoolean keep = CUBIT_FALSE );
      CubitStatus force_simplify_to_cylinder( DLIList<RefFace*> &ref_face_list, DLIList<Body*>& new_body_list, 
                                              CubitBoolean keep = CUBIT_FALSE );
      CubitStatus force_simplify_to_cone( DLIList<RefFace*> &ref_face_list, DLIList<Body*>& new_body_list, 
                                          CubitBoolean keep = CUBIT_FALSE );
      CubitStatus force_simplify_to_sphere( DLIList<RefFace*> &ref_face_list, DLIList<Body*>& new_body_list, 
                                            CubitBoolean keep = CUBIT_FALSE );
      CubitStatus force_simplify_to_torus( DLIList<RefFace*> &ref_face_list, DLIList<Body*>& new_body_list, 
                                           CubitBoolean keep = CUBIT_FALSE );
      // Forces a spline surface to be an analytical of the type specified.

      // *** END HEALER FUNCTIONS *** //

   protected:

      // constructor
      GeometryHealerTool( GeometryHealerEngine *GHEPtr );

   private:

      // singleton pattern class instance interface
      static GeometryHealerTool* instance_;

      // list of healer engines
      DLIList<GeometryHealerEngine*> gheList;
};

#endif

