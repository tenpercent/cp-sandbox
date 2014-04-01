//-----------------------------------------------------------------------------
// Filename      : GeometryHealerTool.cpp
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

// *** BEGIN INCLUDES *** //

#include "GeometryHealerTool.hpp"

#include "GeometryHealerEngine.hpp"

#include "Body.hpp"

#include "RefEntityFactory.hpp"
#include "TopologyEntity.hpp"

#include "CastTo.hpp"
#include "RefFace.hpp"

// *** END INCLUDES *** //

GeometryHealerTool* GeometryHealerTool::instance_ = 0;

// *** BEGIN PUBLIC FUNCTIONS *** //

//-----------------------------------------------------------------------------
// Purpose       : Controls access/creation of the sole instance of this class.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
GeometryHealerTool* GeometryHealerTool::instance( GeometryHealerEngine *GHEPtr )
{
   // Check to see if we have created an instance of the class;
   // if not, proceed to create one.
   if (instance_ == 0) 
   {
      // When creating the instance, we should always have a valid
      // ghePtr. If not, complain.
      instance_ = new GeometryHealerTool(GHEPtr);

      // Check to make sure there's a ref entity factory extant:
      // RefEntityFactory *factory = 
      RefEntityFactory::instance();
   }
   // If there is an existing instance of the class, check if there
   // is a request to set the default engine. If so, do so.
   else if ( GHEPtr != NULL && !instance_->gheList.move_to(GHEPtr) )
   {
      delete instance_->gheList.remove();
      instance_->gheList.insert(GHEPtr);
   }

   // Return the a pointer to the instance of the class.
   return instance_;
}

//-----------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
GeometryHealerTool::~GeometryHealerTool()
{
   for (int i = gheList.size(); i > 0; i--) 
      delete gheList.get_and_step();

   gheList.clean_out();
   instance_ = NULL;
}

// *** BEGIN ENGINE OPERATIONS *** //
/*
//-----------------------------------------------------------------------------
// Purpose       : Sets the default engine.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
void GeometryHealerTool::set_default_engine( GeometryHealerEngine *ghe_ptr )
{
   default_ghe = ghe_ptr;
}
*/
//-----------------------------------------------------------------------------
// Purpose       : Adds a healer engine to the list.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
void GeometryHealerTool::add_ghe( GeometryHealerEngine *ghe_ptr )
{
   assert(ghe_ptr != 0);

   // for now, GeometryHealerTool is only set up for a single healer engine
   // so if more than one healer is being added all the code for this class needs
   // to be reviewed - Byron 08/12/2003
   assert(gheList.size() == 0);

   if (!gheList.move_to(ghe_ptr))
      gheList.append(ghe_ptr);
}
/*
//-----------------------------------------------------------------------------
// Purpose       : Removes a healer engine from the list.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
CubitStatus GeometryHealerTool::remove_ghe( GeometryHealerEngine *ghe_ptr )
{
   assert(ghe_ptr != 0);
   CubitStatus status = CUBIT_FAILURE;

   if (gheList.move_to(ghe_ptr))
   {
      gheList.remove();
      status = CUBIT_SUCCESS;
   }

   return status;
}

//-----------------------------------------------------------------------------
// Purpose       : Gets the list of healer engines.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
void GeometryHealerTool::get_ghe_list( DLIList<GeometryHealerEngine*> &ghe_list )
{
   ghe_list += gheList;
}

//-----------------------------------------------------------------------------
// Purpose       : Gets the current healer engine (first in the list).
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
GeometryHealerEngine *GeometryHealerTool::get_ghe()
{
   GeometryHealerEngine *ghe = NULL;

   if (gheList.size())
   {
      gheList.reset();
      ghe = gheList.get();
   }

   return ghe;
}
*/
//-----------------------------------------------------------------------------
// Purpose       : Returns the healer engine of an entity.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 08/01/03
//-----------------------------------------------------------------------------
GeometryHealerEngine* GeometryHealerTool::get_engine( TopologyBridge *tb_ptr ) const
{
  GeometryHealerEngine *ghe;

  for (int i = 0; i < gheList.size(); i++)
  {
    ghe = gheList.next(i);
    if (ghe->is_healer_engine(tb_ptr))
       return ghe;
  }

  return NULL;
}

//-----------------------------------------------------------------------------
// Purpose       : Returns the healer engine of an entity.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 08/01/03
//-----------------------------------------------------------------------------
GeometryHealerEngine* GeometryHealerTool::get_engine( TopologyEntity *te_ptr ) const
{
  GeometryHealerEngine *ghe;

  TopologyBridge *tb_ptr = te_ptr->bridge_manager()->topology_bridge();
  
  for (int i = 0; i < gheList.size(); i++)
  {
    ghe = gheList.next(i);
    if (ghe->is_healer_engine(tb_ptr))
       return ghe;
  }
  
  return NULL;
}

//-----------------------------------------------------------------------------
// Purpose       : Determines if entities are from the same engine.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 08/01/03
//-----------------------------------------------------------------------------
CubitBoolean GeometryHealerTool::same_healer_engine( DLIList<RefEntity*> &ref_entity_list,
						                                   CubitBoolean check_children ) const
{
   DLIList<RefEntity*> complete_entity_list;

   //Check the check_children option and check all the children if necessary
   if (check_children)
   {
      //Make a complete list of all the RefEntities and their children
      DLIList<RefEntity*> temp = ref_entity_list;
      RefEntity* ref_entity_ptr;

      for (int i = 0; i < ref_entity_list.size(); i++)
      {
         ref_entity_ptr = ref_entity_list.get_and_step();
         complete_entity_list.clean_out();
         ref_entity_ptr->get_all_child_ref_entities(complete_entity_list);
         temp += complete_entity_list;
      }

      complete_entity_list.clean_out();
      complete_entity_list.merge_unique(temp);
   }

   //Now make sure all the RefEntities are from the same geometry engine
   DLIList<TopologyEntity*> te_list;
   CAST_LIST(complete_entity_list, te_list, TopologyEntity);
   return same_healer_engine(te_list);
}

//-----------------------------------------------------------------------------
// Purpose       : Determines if entities are from the same engine.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 08/01/03
//-----------------------------------------------------------------------------
CubitBoolean GeometryHealerTool::same_healer_engine( DLIList<TopologyEntity*> &topo_list ) const
{
   GeometryHealerEngine *gePtr1 = get_engine(topo_list.get_and_step());
   GeometryHealerEngine *gePtr2;

   for (int i = 1; i < topo_list.size(); i++)
   {
      gePtr2 = get_engine(topo_list.get_and_step());
      if (gePtr1 != gePtr2)
      {
         return CUBIT_FALSE;   
      }
   }
   return CUBIT_TRUE;
}

// *** END ENGINE OPERATIONS *** //

// *** BEGIN HEALER FUNCTIONS *** //

CubitStatus GeometryHealerTool::auto_heal_bodies( DLIList<Body*> &body_list, 
                                                  DLIList<Body*> &new_body_list,
                                                  DLIList<TopologyEntity*> &bad_geometry,
                                                  CubitBoolean rebuild, CubitBoolean keep_old,
                                                  CubitBoolean make_tolerant, FILE* logfile_ptr )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
   {
      CubitStatus healer_status = GHEPtr->auto_heal_bodies(body_list, new_body_list, bad_geometry,
                                      rebuild, keep_old, make_tolerant, logfile_ptr);

      if( healer_status == CUBIT_SUCCESS )
        CubitObserver::notify_static_observers(NULL, HEALER_COMPLETED);

      return healer_status;
   }
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::heal_bodies( DLIList<Body*> &body_list, 
                                             DLIList<Body*> &new_body_list,
                                             DLIList<TopologyEntity*> &bad_geometry,
                                             CubitBoolean rebuild, CubitBoolean keep_old,
                                             CubitBoolean make_tolerant, FILE* logfile_ptr )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      return GHEPtr->heal_bodies(body_list, new_body_list, bad_geometry,
                                 rebuild, keep_old, make_tolerant, logfile_ptr);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::analyze_badgeom( DLIList<Body*> &body_list, 
                                                 DLIList<TopologyEntity*> &bad_geometry,
                                                 FILE* logfile)
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      return GHEPtr->analyze_badgeom(body_list, bad_geometry, logfile);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::get_badgeom( DLIList<Body*> &body_list,
                                             DLIList<TopologyEntity*> &bad_geometry )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      return GHEPtr->get_badgeom(body_list, bad_geometry);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::get_tcurves( DLIList<Body*> &body_list,
                                              DLIList<RefEdge*> &t_curves )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      return GHEPtr->get_tcurves(body_list, t_curves);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::heal_incremental( DLIList<Body*> &body_list, 
                                                  DLIList<Body*> &new_bodies,
                                                  DLIList<TopologyEntity*> &bad_geometry,
                                                  double simplify_tol, double stitch_min_tol,
                                                  double stitch_max_tol, double geombuild_tol,
                                                  double analytic_tol, double isospline_tol,
                                                  double reblend_classify_tol, double reblend_tol,
                                                  CubitBoolean keep_old, CubitBoolean make_tolerant,
                                                  FILE* logfile_ptr )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      return GHEPtr->heal_incremental(body_list, new_bodies, bad_geometry, simplify_tol, 
                                          stitch_min_tol, stitch_max_tol, geombuild_tol,
                                          analytic_tol, isospline_tol, reblend_classify_tol,
                                          reblend_tol, keep_old, make_tolerant, logfile_ptr);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

void GeometryHealerTool::list_incremental()
{
   gheList.get()->list_incremental();
}

void GeometryHealerTool::list_tolerances( DLIList<Body*> &body_list )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      GHEPtr->list_tolerances(body_list);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
}

double GeometryHealerTool::get_default_simplify_tol()
{
   return gheList.get()->get_default_simplify_tol();
}

void GeometryHealerTool::set_default_simplify_tol( double tol )
{
   gheList.get()->set_default_simplify_tol(tol);
}

double GeometryHealerTool::get_default_stitch_min_tol()
{
   return gheList.get()->get_default_stitch_min_tol();
}

void GeometryHealerTool::set_default_stitch_min_tol( double tol )
{
   gheList.get()->set_default_stitch_min_tol(tol);
}

double GeometryHealerTool::get_default_stitch_max_tol()
{
   return gheList.get()->get_default_stitch_max_tol();
}

void GeometryHealerTool::set_default_stitch_max_tol( double tol )
{
   gheList.get()->set_default_stitch_max_tol(tol);
}

double GeometryHealerTool::get_default_geombuild_tol()
{
   return gheList.get()->get_default_geombuild_tol();
}

void GeometryHealerTool::set_default_geombuild_tol( double tol )
{
   gheList.get()->set_default_geombuild_tol(tol);
}

double GeometryHealerTool::get_default_analytic_tol()
{
   return gheList.get()->get_default_analytic_tol();
}

void GeometryHealerTool::set_default_analytic_tol( double tol )
{
   gheList.get()->set_default_analytic_tol(tol);
}

double GeometryHealerTool::get_default_isospline_tol()
{
   return gheList.get()->get_default_isospline_tol();
}

void GeometryHealerTool::set_default_isospline_tol( double tol )
{
   gheList.get()->set_default_isospline_tol(tol);
}

double GeometryHealerTool::get_default_reblend_classify_tol()
{
   return gheList.get()->get_default_reblend_classify_tol();
}

void GeometryHealerTool::set_default_reblend_classify_tol( double tol )
{
   gheList.get()->set_default_reblend_classify_tol(tol);
}

double GeometryHealerTool::get_default_reblend_tol()
{
   return gheList.get()->get_default_reblend_tol();
}

void GeometryHealerTool::set_default_reblend_tol( double tol )
{
   gheList.get()->set_default_reblend_tol(tol);
}

void GeometryHealerTool::reset_default_tolerances()
{
   gheList.get()->reset_default_tolerances();
}

void GeometryHealerTool::list_default_tolerances()
{
   gheList.get()->list_default_tolerances();
}

void GeometryHealerTool::clean_attributes( DLIList<Body*>& body_list )
{
   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST_TO_PARENT(body_list, ref_entity_list);

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING bodies from different\n"
                  "       geometry engines is not allowed.\n");
   }

   GeometryHealerEngine* GHEPtr = get_engine(body_list.get());
   if (GHEPtr)
      GHEPtr->clean_attributes(body_list);
   else
      PRINT_ERROR( "Bodies are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
}

CubitBoolean GeometryHealerTool::get_cleanatt_flg()
{
   return gheList.get()->get_cleanatt_flg();
}

void GeometryHealerTool::set_cleanatt_flg( CubitBoolean flg )
{
   gheList.get()->set_cleanatt_flg(flg);
}

int GeometryHealerTool::get_show_method()
{
   return gheList.get()->get_show_method();
}

void GeometryHealerTool::set_show_method( int method )
{
   gheList.get()->set_show_method(method);
}

CubitBoolean GeometryHealerTool::get_show_summary_flg()
{
   return gheList.get()->get_show_summary_flg();
}

void GeometryHealerTool::set_show_summary_flg( CubitBoolean flg )
{
   gheList.get()->set_show_summary_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_details_flg()
{
   return gheList.get()->get_show_details_flg();
}

void GeometryHealerTool::set_show_details_flg( CubitBoolean flg )
{
   gheList.get()->set_show_details_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_on_heal_flg()
{
   return gheList.get()->get_show_on_heal_flg();
}

void GeometryHealerTool::set_show_on_heal_flg( CubitBoolean flg )
{
   gheList.get()->set_show_on_heal_flg(flg);
}

CubitBoolean GeometryHealerTool::get_check_vol_on_heal_flg()
{
   return gheList.get()->get_check_vol_on_heal_flg();
}

void GeometryHealerTool::set_check_vol_on_heal_flg( CubitBoolean flg )
{
   gheList.get()->set_check_vol_on_heal_flg(flg);
}

double GeometryHealerTool::get_vol_on_heal_limit()
{
   return gheList.get()->get_vol_on_heal_limit();
}

void GeometryHealerTool::set_vol_on_heal_limit( double limit )
{
   gheList.get()->set_vol_on_heal_limit(limit);
}

CubitBoolean GeometryHealerTool::get_check_surf_on_heal_flg()
{
   return gheList.get()->get_check_surf_on_heal_flg();
}

void GeometryHealerTool::set_check_surf_on_heal_flg( CubitBoolean flg )
{
   gheList.get()->set_check_surf_on_heal_flg(flg);
}

double GeometryHealerTool::get_surf_on_heal_limit()
{
   return gheList.get()->get_surf_on_heal_limit();
}

void GeometryHealerTool::set_surf_on_heal_limit( double limit )
{
   gheList.get()->set_surf_on_heal_limit(limit);
}

CubitBoolean GeometryHealerTool::get_check_curve_on_heal_flg()
{
   return gheList.get()->get_check_curve_on_heal_flg();
}

void GeometryHealerTool::set_check_curve_on_heal_flg( CubitBoolean flg )
{
   gheList.get()->set_check_curve_on_heal_flg(flg);
}

double GeometryHealerTool::get_curve_on_heal_limit()
{
   return gheList.get()->get_curve_on_heal_limit();
}

void GeometryHealerTool::set_curve_on_heal_limit( double limit )
{
   gheList.get()->set_curve_on_heal_limit(limit);
}

CubitBoolean GeometryHealerTool::get_show_bad_vertices_flg()
{
   return gheList.get()->get_show_bad_vertices_flg();
}

void GeometryHealerTool::set_show_bad_vertices_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_vertices_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_curves_flg()
{
   return gheList.get()->get_show_bad_curves_flg();
}

void GeometryHealerTool::set_show_bad_curves_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_curves_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_coedges_flg()
{
   return gheList.get()->get_show_bad_coedges_flg();
}

void GeometryHealerTool::set_show_bad_coedges_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_coedges_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_loops_flg()
{
   return gheList.get()->get_show_bad_loops_flg();
}

void GeometryHealerTool::set_show_bad_loops_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_loops_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_surfaces_flg()
{
   return gheList.get()->get_show_bad_surfaces_flg();
}

void GeometryHealerTool::set_show_bad_surfaces_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_surfaces_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_shells_flg()
{
   return gheList.get()->get_show_bad_shells_flg();
}

void GeometryHealerTool::set_show_bad_shells_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_shells_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_volumes_flg()
{
   return gheList.get()->get_show_bad_volumes_flg();
}

void GeometryHealerTool::set_show_bad_volumes_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_volumes_flg(flg);
}

CubitBoolean GeometryHealerTool::get_show_bad_bodies_flg()
{
   return gheList.get()->get_show_bad_bodies_flg();
}

void GeometryHealerTool::set_show_bad_bodies_flg( CubitBoolean flg )
{
   gheList.get()->set_show_bad_bodies_flg(flg);
}

void GeometryHealerTool::list_onshow_flgs()
{
   gheList.get()->list_onshow_flgs();
}

CubitBoolean GeometryHealerTool::get_inc_preprocess_flg()
{
   return gheList.get()->get_inc_preprocess_flg();
}

void GeometryHealerTool::set_inc_preprocess_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_preprocess_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_simplify_flg()
{
   return gheList.get()->get_inc_simplify_flg();
}

void GeometryHealerTool::set_inc_simplify_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_simplify_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_stitch_flg()
{
   return gheList.get()->get_inc_stitch_flg();
}

void GeometryHealerTool::set_inc_stitch_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_stitch_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_geombuild_flg()
{
   return gheList.get()->get_inc_geombuild_flg();
}

void GeometryHealerTool::set_inc_geombuild_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_geombuild_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_analytic_flg()
{
   return gheList.get()->get_inc_analytic_flg();
}

void GeometryHealerTool::set_inc_analytic_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_analytic_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_isospline_flg()
{
   return gheList.get()->get_inc_isospline_flg();
}

void GeometryHealerTool::set_inc_isospline_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_isospline_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_reblend_flg()
{
   return gheList.get()->get_inc_reblend_flg();
}

void GeometryHealerTool::set_inc_reblend_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_reblend_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_sharpedge_flg()
{
   return gheList.get()->get_inc_sharpedge_flg();
}

void GeometryHealerTool::set_inc_sharpedge_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_sharpedge_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_genericspline_flg()
{
   return gheList.get()->get_inc_genericspline_flg();
}

void GeometryHealerTool::set_inc_genericspline_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_genericspline_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_wrapup_flg()
{
   return gheList.get()->get_inc_wrapup_flg();
}

void GeometryHealerTool::set_inc_wrapup_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_wrapup_flg(flg);
}

CubitBoolean GeometryHealerTool::get_inc_postprocess_flg()
{
   return gheList.get()->get_inc_postprocess_flg();
}

void GeometryHealerTool::set_inc_postprocess_flg( CubitBoolean flg )
{
   gheList.get()->set_inc_postprocess_flg(flg);
}

CubitStatus GeometryHealerTool::force_simplify_to_plane( DLIList<RefFace*> &ref_face_list, 
                                                         DLIList<Body*>& new_body_list, 
                                                         CubitBoolean keep )
{
   ref_face_list.reset();
   DLIList<RefEntity*> ref_entity_list(ref_face_list.size());
   CAST_LIST_TO_PARENT( ref_face_list, ref_entity_list );

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING faces from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   ref_face_list.reset();
   GeometryHealerEngine* GHEPtr = get_engine((TopologyEntity*)(ref_face_list.get()));
   if (GHEPtr)
      return GHEPtr->force_simplify_to_plane(ref_face_list, new_body_list, keep);
   else
      PRINT_ERROR( "Faces are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::force_simplify_to_cylinder( DLIList<RefFace*> &ref_face_list, 
                                                            DLIList<Body*>& new_body_list, 
                                                            CubitBoolean keep )
{
   ref_face_list.reset();
   DLIList<RefEntity*> ref_entity_list(ref_face_list.size());
   CAST_LIST_TO_PARENT( ref_face_list, ref_entity_list );

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING faces from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   ref_face_list.reset();
   GeometryHealerEngine* GHEPtr = get_engine((TopologyEntity*)(ref_face_list.get()));
   if (GHEPtr)
      return GHEPtr->force_simplify_to_cylinder(ref_face_list, new_body_list, keep);
   else
      PRINT_ERROR( "Faces are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::force_simplify_to_cone( DLIList<RefFace*> &ref_face_list, 
                                                        DLIList<Body*>& new_body_list, 
                                                        CubitBoolean keep )
{
   ref_face_list.reset();
   DLIList<RefEntity*> ref_entity_list(ref_face_list.size());
   CAST_LIST_TO_PARENT( ref_face_list, ref_entity_list );

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING faces from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   ref_face_list.reset();
   GeometryHealerEngine* GHEPtr = get_engine((TopologyEntity*)(ref_face_list.get()));
   if (GHEPtr)
      return GHEPtr->force_simplify_to_cone(ref_face_list, new_body_list, keep);
   else
      PRINT_ERROR( "Faces are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::force_simplify_to_sphere( DLIList<RefFace*> &ref_face_list, 
                                                          DLIList<Body*>& new_body_list, 
                                                          CubitBoolean keep )
{
   ref_face_list.reset();
   DLIList<RefEntity*> ref_entity_list(ref_face_list.size());
   CAST_LIST_TO_PARENT( ref_face_list, ref_entity_list );

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING faces from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   ref_face_list.reset();
   GeometryHealerEngine* GHEPtr = get_engine((TopologyEntity*)(ref_face_list.get()));
   if (GHEPtr)
      return GHEPtr->force_simplify_to_sphere(ref_face_list, new_body_list, keep);
   else
      PRINT_ERROR( "Faces are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

CubitStatus GeometryHealerTool::force_simplify_to_torus( DLIList<RefFace*> &ref_face_list, 
                                                         DLIList<Body*>& new_body_list, 
                                                         CubitBoolean keep )
{
   ref_face_list.reset();
   DLIList<RefEntity*> ref_entity_list(ref_face_list.size());
   CAST_LIST_TO_PARENT( ref_face_list, ref_entity_list );

   if (!same_healer_engine(ref_entity_list, CUBIT_TRUE))
   {
      PRINT_ERROR("HEALING faces from different\n"
                  "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
   }

   ref_face_list.reset();
   GeometryHealerEngine* GHEPtr = get_engine((TopologyEntity*)(ref_face_list.get()));
   if (GHEPtr)
      return GHEPtr->force_simplify_to_torus(ref_face_list, new_body_list, keep);
   else
      PRINT_ERROR( "Faces are of a geometry engine without a healer\n"
                   "         and cannot be healed.\n");
   return CUBIT_FAILURE;
}

// *** END HEALER FUNCTIONS *** //

// *** END PUBLIC FUNCTIONS *** //

// *** BEGIN PROTECTED FUNCTIONS *** //

//-----------------------------------------------------------------------------
// Purpose       : Constructor.
//
// Creator       : Tyronne Lim (CAT)
//
// Creation Date : 07/21/03
//-----------------------------------------------------------------------------
GeometryHealerTool::GeometryHealerTool( GeometryHealerEngine* GHEPtr )
{
  if (GHEPtr != NULL)
     add_ghe(GHEPtr);
}

// *** END PROTECTED FUNCTIONS *** //

// *** BEGIN PRIVATE FUNCTIONS *** //

// *** END PRIVATE FUNCTIONS *** //
