// *** BEGIN INCLUDES *** //
#include "GeometryFeatureTool.hpp"
#include "Body.hpp"
#include "RefEntityFactory.hpp"
#include "TopologyEntity.hpp"
#include "CastTo.hpp"
#include "RefFace.hpp"
#include "BasicTopologyEntity.hpp"
#include "TDSourceFeature.hpp"
// *** END INCLUDES *** //

GeometryFeatureTool* GeometryFeatureTool::instance_ = 0;

// *** BEGIN PUBLIC FUNCTIONS *** //
GeometryFeatureTool* GeometryFeatureTool::instance( GeometryFeatureEngine *GFEPtr )
{
   // Check to see if we have created an instance of the class;
   // if not, proceed to create one.
   if (instance_ == 0)
   {
      // When creating the instance, we should always have a valid
      // gfePtr. If not, complain.
      instance_ = new GeometryFeatureTool(GFEPtr);

      // Check to make sure there's a ref entity factory extant:
      // RefEntityFactory *factory =
      RefEntityFactory::instance();
   }
   // If there is an existing instance of the class, check if there
   // is a request to set the default engine. If so, do so.
   else if ( GFEPtr != NULL && !instance_->gfeList.move_to(GFEPtr) )
   {
      delete instance_->gfeList.remove();
      instance_->gfeList.insert(GFEPtr);
   }

   // Return the a pointer to the instance of the class.
   return instance_;
}

GeometryFeatureTool::~GeometryFeatureTool()
{
   for (int i = gfeList.size(); i > 0; i--)
      delete gfeList.get_and_step();

   gfeList.clean_out();
   instance_ = NULL;
}

// *** BEGIN ENGINE OPERATIONS *** //
void GeometryFeatureTool::add_gfe( GeometryFeatureEngine *gfe_ptr )
{
   assert(gfe_ptr != 0);

   if (!gfeList.move_to(gfe_ptr))
      gfeList.append(gfe_ptr);
}

GeometryFeatureEngine* GeometryFeatureTool::get_engine( TopologyBridge *tb_ptr ) const
{
  GeometryFeatureEngine *gfe;

  for (int i = 0; i < gfeList.size(); i++)
  {
    gfe = gfeList.next(i);
    if (gfe->is_feature_engine(tb_ptr))
       return gfe;
  }

  return NULL;
}

GeometryFeatureEngine* GeometryFeatureTool::get_engine( TopologyEntity *te_ptr ) const
{
  GeometryFeatureEngine *gfe;

  TopologyBridge *tb_ptr = te_ptr->bridge_manager()->topology_bridge();

  for (int i = 0; i < gfeList.size(); i++)
  {
    gfe = gfeList.next(i);
    if (gfe->is_feature_engine(tb_ptr))
       return gfe;
  }

  return NULL;
}

CubitBoolean GeometryFeatureTool::same_feature_engine( DLIList<RefEntity*> &ref_entity_list,
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
   return same_feature_engine(te_list);
}

CubitBoolean GeometryFeatureTool::same_feature_engine( DLIList<TopologyEntity*> &topo_list ) const
{
   GeometryFeatureEngine *gePtr1 = get_engine(topo_list.get_and_step());
   GeometryFeatureEngine *gePtr2;

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

// *** BEGIN FEATURE FUNCTIONS *** //

GeometryFeatureEngine::FeatureType
GeometryFeatureTool::feature_type(BasicTopologyEntity* bte)
{
    // first check to see if the feature information is in the tool data
    TDSourceFeature *tdss = 0;
    tdss = (TDSourceFeature *)bte->
        get_TD(&TDSourceFeature::is_source_feature);

    if(tdss)
        return tdss->source_feature();

    // if the feature information is not in the tool data 
    // then call the engine specific FeatureEngine
    TopologyEntity* topo_ptr = CAST_TO(bte, TopologyEntity);

    if (!topo_ptr)
        return GeometryFeatureEngine::FEATURE_UNDEFINED;

    GeometryFeatureEngine* gfe = get_engine(topo_ptr);

    if(!gfe)
        return GeometryFeatureEngine::FEATURE_UNDEFINED;

    return gfe->feature_type(bte->get_geometry_entity_ptr());
}

CubitStatus GeometryFeatureTool::get_related_by_feature(
    BasicTopologyEntity* bte,
    DLIList<GeometryEntity*>& return_list)
{
    TopologyEntity* topo_ptr = CAST_TO(bte, TopologyEntity);

    if (!topo_ptr)
        return CUBIT_FAILURE;

    GeometryFeatureEngine* gfe = get_engine(topo_ptr);

    if(!gfe)
        return CUBIT_FAILURE;

    return gfe->get_related_by_feature(bte->get_geometry_entity_ptr(),return_list);
}

DLIList<GeometryFeatureEngine::FeatureType> 
GeometryFeatureTool::available_feature_types(BasicTopologyEntity* bte)
{
    TopologyEntity* topo_ptr = CAST_TO(bte, TopologyEntity);
    DLIList<GeometryFeatureEngine::FeatureType> null_list;
    if (!topo_ptr)
        return null_list;

    GeometryFeatureEngine* gfe = get_engine(topo_ptr);

    if(!gfe)
        return null_list;

    return gfe->available_feature_types(bte->get_geometry_entity_ptr());
}

// *** END PUBLIC FUNCTIONS *** //

// *** BEGIN PROTECTED FUNCTIONS *** //

GeometryFeatureTool::GeometryFeatureTool( GeometryFeatureEngine* GFEPtr )
{
  if (GFEPtr != NULL)
     add_gfe(GFEPtr);
}

// *** END PROTECTED FUNCTIONS *** //

// *** BEGIN PRIVATE FUNCTIONS *** //

// *** END PRIVATE FUNCTIONS *** //
