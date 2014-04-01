
//-------------------------------------------------------------------------
// Filename      : GeometryModifyTool.cpp
//
// Purpose       : The tool used to create and modify geometry
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
//-------------------------------------------------------------------------


#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "GeometryModifyTool.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "GeometryEntity.hpp"
#include "GeometryQueryTool.hpp"
#include "MergeTool.hpp"
#include "SurfaceOverlapTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeometryModifyEngine.hpp"
#include "AnalyticGeometryTool.hpp"
#include "DAG.hpp"
#include "TopologyBridge.hpp"
#include "ModelQueryEngine.hpp"
#include "CADefines.hpp"

#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "RefEntityName.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"

#include "CoVertex.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "CoVolume.hpp"

#include "Chain.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Body.hpp"

#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "BodySM.hpp"
#include "ShellSM.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"

#include "CubitAttrib.hpp"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"
#include "CubitTransformMatrix.hpp"

#include "DLIList.hpp"

#include "DRefFaceArray.hpp"
#include "DRefEdgeArray.hpp"
#include "DRefVertexArray.hpp"


#include "CubitMessage.hpp"
#include "SettingHandler.hpp"

#include "CastTo.hpp"
#include "CpuTimer.hpp"

#include "BridgeManager.hpp"
#include "SplitSurfaceTool.hpp"
#include "OffsetSplitTool.hpp"




/* Work around stupid #define hz equal to HZ on IBM */
#ifdef hz
#  undef hz
#endif

GeometryModifyTool* GeometryModifyTool::instance_ = 0;
CubitBoolean GeometryModifyTool::allEdgesImprint = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::groupImprint = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::newIds = CUBIT_FALSE;
CubitBoolean GeometryModifyTool::sepAfterWebcut = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::booleanRegularize= CUBIT_TRUE;
CubitBoolean GeometryModifyTool::oldNames = CUBIT_FALSE;
RefEntity* GeometryModifyTool::copyingEntity = NULL;

//-------------------------------------------------------------------------
// Purpose       : Controls access and creation of the sole instance of this
//                 class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------

GeometryModifyTool* GeometryModifyTool::instance(GeometryModifyEngine *GMEPtr)
{
     // Check to see if we have created an instance of the class
     // If not, proceed to create one.

   if (instance_ == 0)
   {
        // When creating the instance, we should always have a valid
        // SMEPtr. If not, complain.

     instance_ = new GeometryModifyTool (GMEPtr) ;

        // check to make sure there's a ref entity factory extant
     //RefEntityFactory *factory =
     RefEntityFactory::instance();
   }

     // If there is an existing instance of the class, check if there
     // was a request to set default solid modeling engine. If so, be nice
     // to the calling routine :) :) and kindly set the default solid
     // modeling engine.

   else if ( GMEPtr != NULL && instance_->gmeList.move_to(GMEPtr)) {
     delete instance_->gmeList.remove();
     instance_->gmeList.insert(GMEPtr);
   }

     // Return the a pointer to the instance of the class.

   return instance_ ;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryModifyTool::~GeometryModifyTool ()
{
   instance_ = NULL;
     //Kill the geometry modify engine(s)
   int i;
   for (i = gmeList.size(); i > 0; i--)
     delete gmeList.get_and_step();

   gmeList.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a sphere with the
//                 given radius
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::sphere(double radius)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure the radius is not zero or less
   if ( radius <= 0.0 )
   {
      PRINT_ERROR("In GeometryModifyTool::sphere\n"
                  "       Cannot make a sphere of radius %f\n",
                  radius);
      return NULL ;
   }

   BodySM* body_sm = gmeList.get()->sphere(radius);
   if (!body_sm)
   {
      PRINT_ERROR("In GeometryModifyTool::sphere\n"
                  "       Problems building a volume from the sphere.\n");
      return NULL ;
   }

   return GeometryQueryTool::instance()->make_Body(body_sm);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a brick with given width,
//                 depth, and height.
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::brick(double width, double depth, double height)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that none of the input values are zero or less
   if ( width <= GEOMETRY_RESABS || depth <= GEOMETRY_RESABS || height <= GEOMETRY_RESABS )
   {
      PRINT_ERROR("In GeometryModifyTool::brick\n"
                  "    Cannot make a cuboid of size %f x %f x %f\n"
                  "       All dimensions must be > 0.0\n",
                  width, depth, height);
      return NULL;
   }

   BodySM* bodyPtr = gmeList.get()->brick(width, depth, height) ;

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::brick\n"
                  "       Problem creating a brick.\n") ;
      return NULL ;
   }

   return GeometryQueryTool::instance()->make_Body(bodyPtr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a brick at given location,
//                 orientation, and size.
//
// Special Notes : If one of the dimensions is zero, a planar sheet is created.
//
// Creator       : Steve Storm
//
// Creation Date : 10/09/2000
//-------------------------------------------------------------------------
Body* GeometryModifyTool::brick(const CubitVector &center,
                                const CubitVector axes[3],
                                const CubitVector &extension)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   double width = 2.0*extension.x();
   double height = 2.0*extension.y();
   double depth = 2.0*extension.z();
     // First make sure that entity creation is possible.  Allow at most
     // only one of the dimensions to be zero - in which case a planar
     // sheet is created.
   if ( width < 0.0 || height < 0.0 || depth < 0.0 )
   {
      PRINT_ERROR( "Cannot make a brick of size %f x %f x %f\n"
                   "     Negative dimensions are not allowed.\n",
                   width, height, depth );
      return NULL ;
   }

   int wz = width < GEOMETRY_RESABS;
   int hz = height < GEOMETRY_RESABS;
   int dz = depth < GEOMETRY_RESABS;
   if( wz || hz || dz )
   {
      if( (wz + hz + dz) > 1 )
      {
         PRINT_ERROR( "Cannot make a sheet of size %f x %f x %f\n"
            "     At least two dimensions must be nonzero.\n",
            width, height, depth );
         return NULL ;
      }

      // Make a sheet body instead of a cuboid
      CubitVector sheet_axes[2];
      if( wz )
      {
         sheet_axes[0] = axes[1];
         sheet_axes[1] = axes[2];
         return planar_sheet( center, sheet_axes, height, depth );
      }
      else if( hz )
      {
         sheet_axes[0] = axes[2];
         sheet_axes[1] = axes[0];
         return planar_sheet( center, sheet_axes, depth, width );
      }
      else
      {
         sheet_axes[0] = axes[0];
         sheet_axes[1] = axes[1];
         return planar_sheet( center, sheet_axes, width, height );
      }
   }

   else
   {
      BodySM* bodyPtr = gmeList.get()->brick(center, axes, extension) ;

      if (bodyPtr == NULL)
      {
         PRINT_ERROR("In GeometryTool::brick\n"
            "       Problem creating a brick.\n") ;
         return NULL;
      }

      return GeometryQueryTool::instance()->make_Body(bodyPtr);
   }
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a prism with height,
//                 sides, major, minor radii
//
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::prism( double height, int sides, double major, double minor)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( major < minor ||
        major <= GEOMETRY_RESABS || minor <= GEOMETRY_RESABS ||
        height <= GEOMETRY_RESABS || sides < 3 )
   {
      PRINT_ERROR("In GeometryModifyTool::prism\n"
                  "   Cannot make a prism of major-radius = %f,"
                  " minor-radius = %f, height = %f, and"
                  " number of sides = %d\n"
                  "   All dimensions must be > 0.0, the major-radius\n"
                  " should be greater than or equal to the minor-"
                  "radius,\n and the number of sides should be greater"
                  " than 2.\n", major, minor, height, sides);

      return NULL;
   }

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->prism(height,  sides, major, minor) ;

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::prism\n"
                  "       Problems building a volume from the prism.\n") ;
      return NULL;
   }

   return GeometryQueryTool::instance()->make_Body(bodyPtr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a pyramid with given
//                 height, major, minor,top, sides
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::pyramid  ( double height, int sides, double major,
                             double minor,  double top)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( major < minor ||
        major <= GEOMETRY_RESABS || minor <= GEOMETRY_RESABS ||
        height <= GEOMETRY_RESABS || top < 0.0 || sides < 3  )
   {
      PRINT_ERROR("In GeometryModifyTool::pyramid\n"
                  "    Cannot make a pyramid of major-radius = %f,"
                  " minor-radius = %f, height = %f,"
                  " top %f, and number of sides = %d\n"
                  "     All dimensions must be > 0.0, the major-radius"
                  " should be greater than the minor-radius, the "
                  "number of sides should be greater than 2.\n",
                  major, minor, height, top, sides);
      return NULL;
   }

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->pyramid ( height, sides, major, minor, top);
   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::pyramid\n"
                  "      Problems building a volume from the pyramid.\n");
      return NULL;
   }

   return GeometryQueryTool::instance()->make_Body(bodyPtr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a frustum of a cone with
//                 height,
//                 input radius in x-direction at base
//                 input radius in y-direction at base
//                 input radius in x-direction at top
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::cylinder ( double hi, double r1, double r2, double r3 )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( hi <= GEOMETRY_RESABS || r1 <= GEOMETRY_RESABS ||
        r2 <= GEOMETRY_RESABS || r3 < 0.0 )
   {
      PRINT_ERROR("In GeometryModifyTool::cylinder\n"
                  " Cannot make a frustum of a cone with height = %f,"
                  " lower x-direction radius = %f,"
                  " lower y-direction radius = %f, and "
                  "top radius = %f\n"
                  "       All dimensions must be > 0.0\n",
                  hi, r1, r2, r3);
      return NULL;
   }

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->cylinder( hi, r1, r2, r3);

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::cylinder\n"
                  "       Problems building a volume from the conical frustum.\n");
      return NULL;
   }

   return GeometryQueryTool::instance()->make_Body(bodyPtr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a torus with major_radius
//                 and minor_radius
//
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::torus( double r1, double r2 )
{   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }


     // First make sure that the input values make sense
   if ( r1 <= r2 || r1 <= GEOMETRY_RESABS || r2 <= 0.0)
   {
      PRINT_ERROR("In GeometryModifyTool::torus\n"
                  "  Cannot make a torus of major-radius = %f and"
                  " minor-radius = %f\n",
                  r1, r2);
      return NULL;
   }

     // Create a Body that represents the torus
   BodySM* bodyPtr = gmeList.get()->torus(r1, r2) ;

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::torus\n"
                  "       Problems building a volume from the torus.\n") ;
      return NULL ;
   }

   return GeometryQueryTool::instance()->make_Body(bodyPtr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the given corner locations.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 05/10/05
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet( const CubitVector& p1,
                                        const CubitVector& p2,
                                        const CubitVector& p3,
                                        const CubitVector& p4 )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

  // Create a Body that represents the sheet
  BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;

  if( body_ptr == NULL )
  {
    PRINT_ERROR("In GeometryTool::planar_sheet\n"
      "       Problems building a volume from the sheet.\n") ;
    return NULL;
  }

  return GeometryQueryTool::instance()->make_Body(body_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the given input location,
//                 orientation and size.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 10/09/00
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet ( const CubitVector &center,
                                         const CubitVector axes[2],
                                         double width, double height )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   CubitVector p1, p2, p3, p4;

   // Get the corners of the sheet
   center.next_point( axes[0], width/2.0, p1 );
   p1.next_point( axes[1], -height/2.0, p1 );
   p1.next_point( axes[1], height, p2 );
   p2.next_point( axes[0], -width, p3 );
   p3.next_point( axes[1], -height, p4 );

   // Create a Body that represents the sheet
   BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;

   if( body_ptr == NULL )
   {
      PRINT_ERROR("In GeometryTool::planar_sheet\n"
         "       Problems building a volume from the sheet.\n") ;
      return NULL;
   }

   return GeometryQueryTool::instance()->make_Body(body_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the minimal amount of area
//                 to cut through the given 3D bounding box.  The sheet
//                 can be extended optionally to make it larger than this
//                 minimal size.
//
// Special Notes : extension_types: 0-none, 1-percentage, 2-absolute
//
// Creator       : Steve Storm
//
// Creation Date : 10/19/98
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet ( const CubitPlane& plane,
                                         const CubitBox& bounding_box,
                                         int extension_type,
                                         double extension )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   CubitVector p1, p2, p3, p4;

   // Get the corners of the sheet
   if( AnalyticGeometryTool::instance()->
      min_pln_box_int_corners( plane, bounding_box, extension_type,
                               extension, p1, p2, p3, p4 ) == CUBIT_FAILURE )
      return NULL;

   // Create a Body that represents the sheet
   BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;

   if( body_ptr == NULL )
   {
      PRINT_ERROR("In GeometryModifyTool::planar_sheet\n"
         "       Problems building a volume from the sheet.\n") ;
      return NULL ;
   }

   return GeometryQueryTool::instance()->make_Body(body_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a location
//                 to create a RefVertex. The underlying representation
//                 of the RefVertex is determined by the default
//                 GeometryModifyEngine.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------

RefVertex* GeometryModifyTool::make_RefVertex(
                                        CubitVector const& point, int color) const
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // Call the default GeometryModifyEngine to create a new Point
   Point* point_ptr = gmeList.get()->make_Point(point);

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( point_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefVertex\n"
                    "         Got a NULL pointer to a Point.\n"
                    "         Cannot make a RefVertex.\n") ;
      return (RefVertex *)NULL ;
   }

     // Use the Point to create a RefVertex
   RefVertex* ref_vertex_ptr = RefEntityFactory::instance()->construct_RefVertex(point_ptr) ;

   ref_vertex_ptr->color(color);

     // Send a message to the model indicating the vertex was created
   CubitObserver::notify_static_observers(ref_vertex_ptr, FREE_REF_ENTITY_GENERATED);

     // Return the newly created RefVertex
   return ref_vertex_ptr ;
}

Body *GeometryModifyTool::make_Body(Surface *surface) const
{
  // make a Body from a Surface
  BodySM *body_sm = surface->bodysm();

  GeometryModifyEngine* gme = get_engine(surface);

  if (gme && !body_sm)
    body_sm = gme->make_BodySM(surface);

  assert(body_sm != 0);

  return GeometryQueryTool::instance()->make_Body(body_sm);
}

//-------------------------------------------------------------------------
// Purpose       : This function takes RefEdge type information, too
//
// Creator       : David White
//
// Creation Date : 10/9/97
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(RefEdge *ref_edge_ptr, bool copy_attribs ) const
{
  if ( ref_edge_ptr == NULL )
  {
    PRINT_ERROR("curve is NULL\n");
    return (RefEdge*)NULL;
  }

  TopologyBridge* bridge = 0;
  GeometryModifyEngine* engine = get_engine(ref_edge_ptr, &bridge);
  Curve *old_curve = dynamic_cast<Curve*>(bridge);
  if (engine == NULL)
  {
     PRINT_ERROR( "%s (curve %d) does not have a modify engine.\n",
                  ref_edge_ptr->entity_name().c_str(),
                  ref_edge_ptr->id() );
     return 0;
  }

  Curve *tmp_curve;
  if( copy_attribs )
  {
    DLIList<RefEntity*> tmp_list;
    tmp_list.append( ref_edge_ptr );
    TopologyBridge *curve_bridge;
    prepare_for_copy( ref_edge_ptr, curve_bridge );
    tmp_curve = CAST_TO( curve_bridge, Curve);
  }
  else
    tmp_curve = old_curve;

  Curve *new_curve = engine->make_Curve( tmp_curve );
  if (!new_curve)
  {
    if( copy_attribs )
      clean_up_from_copy_failure( old_curve );
    return (RefEdge *)NULL;
  }

  TopologyBridge *new_curve_bridge;
  if( copy_attribs )
  {
    new_curve_bridge = new_curve;
    finish_copy( new_curve_bridge, old_curve );
    new_curve = CAST_TO( new_curve_bridge, Curve);
  }

  // Complete the task of linking this new Curve into the rest of the
  // geometry datastructures and return the new RefEdge.
  RefEdge *new_ref_edge = GeometryQueryTool::instance()->make_free_RefEdge(new_curve);

  return new_ref_edge;
}

CubitStatus GeometryModifyTool::prepare_for_copy( RefEntity *ref_ent,
                                                  TopologyBridge *&top_bridge )
{
  //save attribute settings
  CGMApp::instance()->save_current_attribute_states();

  CGMApp::instance()->attrib_manager()->auto_flag(1);
  // Groups are saved directly in the Cubit file, not with attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GROUP, CUBIT_FALSE);

  // The mesh container is only used for the Exodus save/resore method
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_MESH_CONTAINER, CUBIT_FALSE);

  // Genesis Entities are saved directly in the Cubit file, not with attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);

  //Don't save out entity ids
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_ENTITY_ID, CUBIT_FALSE);

  //Don't save out graphics
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);

  //Flag to say "we're copying" attributes
  copyingEntity = ref_ent;

  //Get all children
  DLIList<RefEntity*> child_list;
  ref_ent->get_all_child_ref_entities( child_list );

  //Puts attributes on all the entities
  child_list.append( ref_ent );
  CubitAttribUser::auto_update_cubit_attrib( child_list );

  //Prepare virtual
  TopologyEntity *topo_ptr = dynamic_cast<TopologyEntity*>( ref_ent );
  top_bridge = topo_ptr->bridge_manager()->topology_bridge();
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( top_bridge );
  GeometryQueryTool::instance()->ige_export_geom( bridge_list );

  //Should only have 1 bridge in the list
  assert( bridge_list.size() == 1 );
  TopologyBridge *tmp_bridge_after = bridge_list.get();
  if( top_bridge != tmp_bridge_after)
    top_bridge = tmp_bridge_after;

  //Done copying attributes
  copyingEntity = NULL;

  return CUBIT_SUCCESS;

}

CubitStatus GeometryModifyTool::clean_up_from_copy_failure( TopologyBridge *old_bridge )
{
  //Remove attributes on underlying entities of virtual geometry
  //entities which do not have corresponding RefEntities.
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( old_bridge );
  GeometryQueryTool::instance()->ige_remove_attributes( bridge_list );

  //Remove attributes on original entity and children
  DLIList<RefEntity*> child_list;
  RefEntity *ref_ent = CAST_TO( old_bridge->topology_entity(), RefEntity );
  ref_ent->get_all_child_ref_entities( child_list );
  //child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );

  CGMApp::instance()->restore_previous_attribute_states();

return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::finish_copy( TopologyBridge *&new_bridge,
                                             TopologyBridge *old_bridge )
{
  //Remove attributes on underlying entities of virtual geometry
  //entities which do not have corresponding RefEntities.
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( old_bridge );
  GeometryQueryTool::instance()->ige_remove_attributes( bridge_list );

  //Remove attributes on original entity and children
  DLIList<RefEntity*> child_list;
  RefEntity *ref_ent = CAST_TO( old_bridge->topology_entity(), RefEntity );
  ref_ent->get_all_child_ref_entities( child_list );
  //child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );


  //Restore virtual
  //Could create virtual geometry here so the Topology Bridge can change
  bridge_list.clean_out();
  bridge_list.append( new_bridge );
  TopologyBridge *tmp_bridge_before = new_bridge;
  GeometryQueryTool::instance()->ige_import_geom( bridge_list );
  assert( bridge_list.size() == 1 );
  if( tmp_bridge_before != bridge_list.get() )
    new_bridge = bridge_list.get();

  //make the RefEntities
  Curve *curve = NULL;
  Surface *surface = NULL;
  Lump *lump = NULL;
  BodySM *body = NULL;

  RefEntity *new_ref_ent = NULL;
  if( (curve = CAST_TO( new_bridge, Curve ) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_RefEdge( curve );
  else if( (surface = CAST_TO( new_bridge, Surface) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_RefFace( surface );
  else if( (lump = CAST_TO( new_bridge, Lump) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_Body( lump->bodysm() );
  else if( (body = CAST_TO( new_bridge, BodySM) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_Body( body );

  //actuate the attributes on everything
  child_list.clean_out();
  new_ref_ent->get_all_child_ref_entities( child_list );
  child_list.append( new_ref_ent );
  GeometryQueryTool::import_actuate( child_list );

  //Remove attributes on new entity and children
  child_list.clean_out();
  new_ref_ent->get_all_child_ref_entities( child_list );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( new_ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );

  CGMApp::instance()->restore_previous_attribute_states();
  return CUBIT_SUCCESS;
}

GeometryModifyEngine*
GeometryModifyTool::make_RefEdge_common( RefVertex* start_vertex,
                                         RefVertex* end_vertex,
                                         Point*& start_point,
                                         Point*& end_point,
                                         RefFace* ref_face,
                                         Surface** surface ) const
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme = 0;
  start_point = end_point = 0;

  bool new_start_point = start_vertex->get_parents() > 0;
  bool new_end_point = end_vertex->get_parents() > 0;

  if (ref_face)
    entity_list.append( ref_face );
  if (!new_start_point)
    entity_list.append( start_vertex );
  if (!new_end_point)
    entity_list.append( end_vertex );

  if (entity_list.size())
    gme = common_modify_engine( entity_list, bridge_list );

  if (gme)
  {
    bridge_list.reset();
    if (ref_face)
      *surface = dynamic_cast<Surface*>(bridge_list.get_and_step());
    if (!new_start_point)
      start_point = dynamic_cast<Point*>(bridge_list.get_and_step());
    if (!new_end_point)
      end_point = dynamic_cast<Point*>(bridge_list.get_and_step());
  }
  else if (ref_face)
  {
    gme = get_engine( ref_face, &bridge );
    if (!gme)
    {
      PRINT_ERROR("No modify engine for surface %d\n", ref_face->id()) ;
      return 0;
    }
    *surface = dynamic_cast<Surface*>(bridge);
    GeometryQueryEngine* gqe = bridge->get_geometry_query_engine();
    if (!new_start_point)
      start_point = dynamic_cast<Point*>( start_vertex->
                                  bridge_manager()->topology_bridge( gqe ) );
    if (!new_end_point)
      end_point = dynamic_cast<Point*>( end_vertex->
                                  bridge_manager()->topology_bridge( gqe ) );
  }
  else if (!new_start_point && (gme = get_engine( start_vertex, &bridge )))
  {
    start_point = dynamic_cast<Point*>(bridge);
    if (!new_end_point)
      end_point = dynamic_cast<Point*>( end_vertex->bridge_manager()->
                    topology_bridge( bridge->get_geometry_query_engine() ) );
  }
  else if (!new_end_point && (gme = get_engine( end_vertex, &bridge )))
  {
    end_point = dynamic_cast<Point*>(bridge);
  }
  else
  {
    gme = get_gme();
  }


  if (!start_point)
    start_point = gme->make_Point( start_vertex->coordinates() );
  if (!end_point)
    end_point = gme->make_Point( end_vertex->coordinates());

  return gme;
}

//-------------------------------------------------------------------------
// Purpose       : This function uses the vertices and surface to create
//                 a curve along this surface.  The solid modeler will actually
//                 do most of the work.
// Note          : The optional third vertex will not be part of the curve,
//                 it is used for interpolation purposes.
//
// Creator       : David White
//
// Creation Date : 10/9/97
//
// Rewriting for new GeometryQueryTool interface - J.Kraftcheck 9/03
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    RefFace* ref_face_ptr,
                                    RefVertex const* ref_vertex_3 ) const
{
	  //make sure that the vertices are on the ref-face.
   CubitVector vert_1 = ref_vertex_1->coordinates();
   CubitVector vert_2 = ref_vertex_2->coordinates();

   ref_face_ptr->move_to_surface(vert_1);
   ref_face_ptr->move_to_surface(vert_2);

   if (!ref_vertex_1->coordinates().within_tolerance(vert_1, GEOMETRY_RESABS))
   {
      PRINT_ERROR("vertices must lie within tolerance to the given"
                  " surface.\n"
                  "%s (Vertex %d) does not lie on %s (Surface %d).\n",
                  ref_vertex_1->entity_name().c_str(),
                  ref_vertex_1->id(),
                  ref_face_ptr->entity_name().c_str(),
                  ref_face_ptr->id() );
      return (RefEdge*)NULL;
   }
   else if (!ref_vertex_2->coordinates().within_tolerance(vert_2, GEOMETRY_RESABS))
   {
      PRINT_ERROR("vertices must lie within tolerance to the given"
                  " surface.\n"
                  "%s (Vertex %d) does not lie on %s (Surface %d).\n",
                  ref_vertex_2->entity_name().c_str(),
                  ref_vertex_2->id(),
                  ref_face_ptr->entity_name().c_str(),
                  ref_face_ptr->id() );
      return (RefEdge*)NULL;
   }
     //Now let us find the points on the surface.  We want to
     //create them as we go for accuracy.

	  // Get the GME of the first RefVertex.
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the end Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;
   Surface* surface_ptr = NULL;

   // Look for a common GeometryModifyEngine
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2,
                                ref_face_ptr, &surface_ptr );

  //If we did not find a common GeometryModifyEngine, fail
  if( ! GMEPtr )
  {
    PRINT_ERROR("Surface %d, vertex %d and vertex %d do not "
      "belong to the same geometric modeling engine.\n", ref_face_ptr->id(),
      ref_vertex_1->id(), ref_vertex_2->id() );

    return 0;
  }

   CubitVector *third_vector_ptr = NULL;
   CubitVector third_vector;
   if ( ref_vertex_3 != NULL )
   {
      third_vector = ref_vertex_3->coordinates();
      third_vector_ptr = &third_vector;
      ref_face_ptr->move_to_surface( third_vector );
   }

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL );

     // Request the GME to create a Curve using the Points
   Curve *curve_ptr;
   curve_ptr = GMEPtr->make_Curve( point_ptr1,
                          point_ptr2,
                          surface_ptr,
                          third_vector_ptr);


     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "\tProblems making a spline curve from Vertex %d and %d\n"
                    "\tand the input list of positions.\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   return GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function takes RefEdge type information, two
//                 RefVertices, and a list of positions in space (represented
//                 by CubitVectors) to create a RefEdge. The RefVertices must
//                 be associated with the same GeometryModifyEngine.
//                 The return value can be a NULL pointer, if the RefEdge
//                 cannot be succesfully made for some reason.
//
// Special Notes : The interpolated curve is a spline curve.
//                 If the input refface_ptr is not NULL, the points are first
//                 moved to the surface before interpolation.
//
// Creator       : Malcolm Panthaki
//
// Creation Date : 05/15/97
//
// Rewriting for new GeometryQueryTool interface - J.Kraftcheck 9/03
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(GeometryType ref_edge_type,
                                    RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    DLIList<CubitVector*>& vector_list,
                                    RefFace* refface_ptr) const
{
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the end Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;
   Surface* surface_ptr = NULL;



   // Look for a common GeometryModifyEngine
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2,
                                refface_ptr, &surface_ptr );
  if (!GMEPtr)
    return 0;

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL );

     // Request the GME to create a Curve using the Points
   Curve* curve_ptr = GMEPtr->make_Curve(ref_edge_type,
                                         point_ptr1,
                                         point_ptr2,
                                         vector_list,
                                         refface_ptr ? refface_ptr->get_surface_ptr() : 0);

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "         Got a NULL pointer to a Curve.\n"
                    "         Problems making a spline RefEdge from RefVertex %d and %d\n"
                    "         and the input list of positions.\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   return GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information, two
//                 RefVertices, an intermediate position, and a sense
//                 information to create a RefEdge. The RefVertices must
//                 be associated with the same GeometryModifyEngine.
//                 The return value can be a NULL pointer, if the RefEdge
//                 cannot be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(GeometryType ref_edge_type,
                                    RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    CubitVector const* intermediate_point,
                                    CubitSense sense) const
{
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;

    // Look for a common geometric modeling engine to use
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2 );
  if (!GMEPtr)
    return 0;

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL ) ;

     // Request the GME to create a Curve using the Points
   Curve* curve_ptr = GMEPtr->make_Curve(ref_edge_type,
                                         point_ptr1,
                                         point_ptr2,
                                         intermediate_point,
                                         sense) ;

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "         Got a NULL pointer to a Curve.\n"
                    "         Problems making RefEdge from RefVertex %d and %d\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   return GeometryQueryTool::instance()->make_free_RefEdge( curve_ptr );
}
//-------------------------------------------------------------------------
// Purpose       : This function creates a RefFace from an existing RefFace.
//                 The new face is a sheet body, double sided face with no
//                 volume to the body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------

RefFace* GeometryModifyTool::make_RefFace(RefFace *from_ref_face,
                                    CubitBoolean extended_from) const
{
   GeometryModifyEngine* GME_ptr = get_engine(from_ref_face);

   if (!GME_ptr)
   {
      PRINT_ERROR("Cannot create surface from another virtual surface.\n");
      return NULL;
   }
     //From the surface create a new surface.
   TopologyBridge* bridge = 0;
   GeometryModifyEngine* engine = get_engine(from_ref_face, &bridge);
   Surface* old_surface_ptr = dynamic_cast<Surface*>(bridge);
   if ( engine == NULL )
   {
      PRINT_ERROR("%s (surface %d) does not have a modify engine.\n",
                  from_ref_face->entity_name().c_str(),
                  from_ref_face->id() );
      return (RefFace*)NULL;
   }

  //this list will get all the TB's what we'll be copying
#ifdef BOYD17
  DLIList<TopologyBridge*> original_bridges;
#endif
  TopologyBridge *top_bridge = old_surface_ptr;
  if( !extended_from )
  {
    prepare_for_copy( from_ref_face, top_bridge );
  }


  Surface *tmp_surface = CAST_TO( top_bridge, Surface );
  Surface* new_surface_ptr = engine->make_Surface( tmp_surface,
                                                    extended_from );

  if (!new_surface_ptr)
  {
    PRINT_ERROR("Surface copy failed.\n");
    if( !extended_from )
      clean_up_from_copy_failure( old_surface_ptr );
    return 0;
  }

  if( !extended_from  )
  {
    TopologyBridge *top_bridge_new = new_surface_ptr;
    finish_copy( top_bridge_new, old_surface_ptr );
    new_surface_ptr = CAST_TO( top_bridge_new, Surface );
  }

  Body *new_Body = make_Body(new_surface_ptr);
  DLIList<RefFace*> ref_faces;
  new_Body->ref_faces(ref_faces);
  assert(ref_faces.size() > 0);

  return ref_faces.get();
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a list of
//                 RefEdges to create a Body with just one RefFace. The
//                 underlying representation of the RefFace is determined
//                 by the GeometryModifyEngine of the RefEdges. All the
//                 RefEdges in the list must be associated with the same
//                 GeometryModifyEngine. The return value can be a
//                 NULL pointer, if the RefFace cannot be succesfully
//                 made for some reason.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 10/23/97
//-------------------------------------------------------------------------

RefFace* GeometryModifyTool::make_RefFace(GeometryType ref_face_type,
                                    DLIList<RefEdge*>& ref_edge_list,
                                    RefFace *ref_face_ptr,
                                    bool check_edges ) const
{
    //Look for a common GeometryModifyEngine for all of
    //the RefEdges.

  const int count = ref_edge_list.size() + (ref_face_ptr ? 1 : 0);
  DLIList<TopologyBridge*> bridge_list(count);
  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( ref_edge_list, entity_list );
  if( ref_face_ptr ) entity_list.append( ref_face_ptr );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return 0;
  }

   Surface* old_surface_ptr = 0;
   if (ref_face_ptr)
    old_surface_ptr = dynamic_cast<Surface*>(bridge_list.pop());

   //perform edge checks
   DLIList<RefEdge*> copied_ref_edges;
   if( check_edges ) 
   { 
      DLIList<RefEdge*> vtx_edges;
      DLIList<RefVertex*> vtx_list; 
      for ( int i = ref_edge_list.size(); i > 0; i-- )
      {
         RefEdge *ref_edge = ref_edge_list.get();

         CubitBoolean other_edge = CUBIT_FALSE;
         ref_edge->ref_vertices(vtx_list);
         while( vtx_list.size() )
         {
            vtx_edges.clean_out();
            vtx_list.pop()->ref_edges( vtx_edges );
            while( vtx_edges.size() )
            {
               if( ! ref_edge_list.is_in_list( vtx_edges.pop() ) )
               {
                   other_edge = CUBIT_TRUE;
                   vtx_list.clean_out();
                   break;
               }
            }
            if(other_edge)
	       break;
         }

         if ( other_edge || (ref_edge->get_parents() > 0) )
         {
            RefEdge *replacement_edge = GeometryModifyTool::instance()->make_RefEdge( ref_edge );
            if (!replacement_edge)
            {
                PRINT_WARNING("Creation of Surface Unsuccessful\n");
                return (RefFace *)NULL;
            }
            ref_edge_list.change_to( replacement_edge );
            copied_ref_edges.append( replacement_edge );
         }
         ref_edge_list.step();
      }//end 'for' loop
   }

   DLIList<Curve*> curve_list(ref_edge_list.size());
   if (copied_ref_edges.size() > 0)
   {
     entity_list.clean_out();
     bridge_list.clean_out();
     CAST_LIST_TO_PARENT( ref_edge_list, entity_list );

     GME_ptr =
        common_modify_engine( entity_list, bridge_list ); 
   }
   CAST_LIST( bridge_list, curve_list, Curve );

     // Use the Curves to create a Surface
   Surface* surface_ptr = GME_ptr->make_Surface(ref_face_type, curve_list,
                                                old_surface_ptr, CUBIT_FALSE) ;

   if (surface_ptr == NULL) {
     PRINT_ERROR("Couldn't make new RefFace.\n");
     for(int i=copied_ref_edges.size(); i--; )
       GeometryQueryTool::instance()->delete_RefEdge( copied_ref_edges.get_and_step() );
     return NULL;
   }

   GeometryQueryTool* const gqt = GeometryQueryTool::instance();
   RefFace* result_face = gqt->make_free_RefFace(surface_ptr);
   gqt->cleanout_deactivated_geometry();
   return result_face;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a list of
//                 RefFaces to create a RefVolume. The underlying
//                 representation of the RefVolume is determined by the
//                 GeometryModifyEngine of the RefFaces. All the
//                 RefFaces in the list must be associated with the same
//                 GeometryModifyEngine. The return value can be a
//                 NULL pointer, if the RefVolume cannot be succesfully
//                 made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------
RefVolume*
GeometryModifyTool::make_RefVolume( DLIList<RefFace*>& ref_face_list) const
{
  GeometryModifyEngine* GME_ptr = 0;

  DLIList<TopologyEntity*> entity_list(ref_face_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_face_list.size());

    //Look for a common GeometryModifyEngine
 CAST_LIST_TO_PARENT( ref_face_list, entity_list );
 GME_ptr = common_modify_engine( entity_list, bridge_list );

  if( !GME_ptr )
  {
    PRINT_ERROR("The surfaces are not owned by the same geometry "
      "engine.  Cannot create a volume from them.\n");
    return 0;
  }

    //Get the list of surfaces.
  DLIList<Surface*> surface_list(ref_face_list.size());
  CAST_LIST( bridge_list, surface_list, Surface );
  surface_list.remove_all_with_value(0);
  assert(surface_list.size() == ref_face_list.size());

     // Use the Surfaces to create a Lump
   Lump* lump_ptr = GME_ptr->make_Lump(surface_list) ;

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if (lump_ptr == NULL)
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefVolume\n"
                    "         Got a NULL pointer to a Lump.\n"
                    "         Cannot make a RefVolume.\n") ;
      return (RefVolume *)NULL;
   }

     // Use the Lump to make a RefVolume
    CubitBoolean modified;
    RefVolume* ref_volume_ptr =
      GeometryQueryTool::instance()->make_RefVolume(lump_ptr, modified);

     // Return the new RefVolume
   return ref_volume_ptr ;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a list of RefVolumes to create a
//                 Body. The underlying representation of the Body is
//                 determined by the GeometryModifyEngine of the
//                 RefVolumes. All the RefVolumes in the list must be
//                 associated with the same GeometryModifyEngine. The
//                 return value can be a NULL pointer, if the RefFace
//                 cannot be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(DLIList<RefVolume*>& ref_volume_list) const
{
  DLIList<TopologyEntity*> topo_list;
  DLIList<Lump*> lump_list;
  DLIList<TopologyBridge*> bridge_list;
  CAST_LIST_TO_PARENT( ref_volume_list, topo_list );
  GeometryModifyEngine* GME_ptr =
    common_modify_engine( topo_list, bridge_list, CUBIT_FALSE );
  if( ! GME_ptr )
  {
    PRINT_ERROR("The specified volumes do not belong to the same "
      "geometry engine, and therefore cannot be used to form a "
      "single body.");
    return 0;
  }
  CAST_LIST( bridge_list, lump_list, Lump );
  assert( ref_volume_list.size() == lump_list.size() );


     // Use the Lumps to create a BodySM
   BodySM* bodySM_ptr = GME_ptr->make_BodySM(lump_list) ;

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( bodySM_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_Body\n"
                    "         Got a NULL pointer to a BodySM.\n"
                    "         Cannot make a Body.\n") ;
      return (Body *)NULL;
   }

   DLIList<BodySM*> bodysm_list;
   bodysm_list.append(bodySM_ptr);
  
   DLIList<Body*> input_bodies;
   for (int i =0; i < ref_volume_list.size(); i++)
     input_bodies.append(ref_volume_list.get_and_step()->get_body_ptr());

   DLIList<Body*> results;

   finish_sm_op(input_bodies, bodysm_list, results, CUBIT_TRUE);
   
   return results.get();
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a sheet body from an existing
//                 RefFace.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(RefFace *from_ref_face,
                              CubitBoolean extended_from) const
{
     // Given the arguments, make a RefFace.
   RefFace* new_ref_face = this->make_RefFace(from_ref_face,
                                              extended_from);
   if (!new_ref_face)
   {
      return (Body *)NULL;
   }

   DLIList<Body*> bodies;
   new_ref_face->bodies(bodies);

   if (!bodies.size()) {
     GeometryEntity *ge_ptr = new_ref_face->get_geometry_entity_ptr();
     Surface *surf_ptr = CAST_TO(ge_ptr, Surface);
     assert(surf_ptr != 0);

     BodySM *body_sm = gmeList.get()->make_BodySM(surf_ptr);
     return GeometryQueryTool::instance()->make_Body(body_sm);
   }

   else return bodies.get();
}
//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a list of
//                 RefEdges to create a Body that has just one RefFace.
//                 The underlying representation of the Body is
//                 determined by the GeometryModifyEngine of the
//                 RefEdges. All the RefEdges in the list must be
//                 associated with the same GeometryModifyEngine. The
//                 return value can be a NULL pointer, if the Body cannot
//                 be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(GeometryType ref_face_type,
                              DLIList<RefEdge*>& ref_edge_list,
                              RefFace *ref_face_ptr) const
{
     // Given the arguments, make a RefFace.
   RefFace* new_ref_face = this->make_RefFace(ref_face_type,
                                              ref_edge_list,
                                              ref_face_ptr) ;

   if( new_ref_face == NULL )
      return NULL;

      // If new_ref_face doesn't have a body, create one.
   DLIList<Body*> bodies;
   new_ref_face->bodies(bodies);

   if (!bodies.size()) {
     Surface *surf_ptr = new_ref_face->get_surface_ptr();
     assert(surf_ptr != 0);
     GeometryModifyEngine* engine = get_engine(surf_ptr);
     assert(engine != 0);

     BodySM *body_sm = engine->make_BodySM(surf_ptr);
     Body* body = GeometryQueryTool::instance()->make_Body(body_sm);
     bodies.append( body );
   }

     // Return the body containing the new RefFace.
   return bodies.get();
}

//-------------------------------------------------------------------------
// Purpose       : The following functions copy the input Body to create
//                 a new one.  Some transformations may be applied to the
//                 newly copied Body before returning it:
//                    copy
//                    copy_and_move
//                    copy_and_rotate
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/1/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::copy_body ( Body* bodyPtr )
{
   BodySM* body_sm = bodyPtr->get_body_sm_ptr();
   if (!body_sm)
   {
     PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", bodyPtr->id());
     return 0;
   }

   //this list will get all the TB's what we'll be copying
#ifdef BOYD17
   DLIList<TopologyBridge*> original_bridges;
#endif
   DLIList<RefEntity*> tmp_list;
   tmp_list.append( bodyPtr );
   TopologyBridge *top_bridge;
   prepare_for_copy( bodyPtr, top_bridge );

   BodySM *tmp_body_sm = CAST_TO( top_bridge, BodySM );

   GeometryModifyEngine* GMEPtr = get_engine(tmp_body_sm);
   if ( GMEPtr == NULL )
   {
     clean_up_from_copy_failure( top_bridge );
     PRINT_ERROR("Cannot copy volume %d\n"
                 "     Copy is supported for bodies based on solid "
                 "models only\n", bodyPtr->ref_volume()->id() );
     return NULL ;
   }

   BodySM* new_body = GMEPtr->copy_body(tmp_body_sm);

   if (!new_body)
   {
     clean_up_from_copy_failure( top_bridge );
     PRINT_ERROR("Failed to copy volume %d\n", bodyPtr->ref_volume()->id());
     return 0;
   }

   TopologyBridge *top_bridge_new = new_body;
   finish_copy( top_bridge_new, body_sm );

   new_body = CAST_TO( top_bridge_new, BodySM );
   return GeometryQueryTool::instance()->make_Body(new_body);
}

//-------------------------------------------------------------------------
// Purpose       : Check if bodies can be webcut
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::okay_to_modify(
                                 DLIList<Body*>& webcut_body_list,
                                 const char* op ) const
{
  CubitStatus ret = CUBIT_SUCCESS;

  if(strcmp(op, "WEBCUT") &&
     strcmp(op, "CHOP") &&
     strcmp(op, "IMPRINT") &&
     strcmp(op, "SPLIT"))
  {
    if ( contains_intermediate_geom(webcut_body_list) )
    {
      PRINT_ERROR("Performing %s on volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n",
                  op);
      ret = CUBIT_FAILURE;
    }
  }
  else
  {
    if(contains_partitions(webcut_body_list))
    {
      PRINT_ERROR("Performing %s on volumes containing virtual partitions.\n"
                  "This is currently not supported.\n\n",
                  op);
      ret = CUBIT_FAILURE;
    }
  }

  return ret;
}


//-------------------------------------------------------------------------
// Purpose       : Common code for finishing webcut operations.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::finish_webcut( DLIList<Body*>& webcut_body_list,
                                               DLIList<BodySM*>& result_body_sms,
                                               CubitBoolean merge,
                                               CubitStatus status,
                                               DLIList<Body*>& result_list,
                                               CubitBoolean print_info) const
{
   if (!finish_sm_op(webcut_body_list, result_body_sms, result_list,(bool)print_info))
     status = CUBIT_FAILURE;

   DLIList<Body*> temp_result_list;
   if (status)
   {
     status = separate_body_after_webcut( result_list, temp_result_list);
     result_list = temp_result_list;
   }

   if (merge && status)
   {
     DLIList<Body*> temp_results(result_list);
     temp_results += webcut_body_list;
     status = MergeTool::instance()->merge_bodies( temp_results );
   }

   return status;
}


//-------------------------------------------------------------------------
// Purpose       : Complete solid modeling operation on bodies
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::finish_sm_op( DLIList<Body*>& input_bodies,
                                              DLIList<BodySM*>& new_bodies,
                                              DLIList<Body*>& result_bodies,
                                              bool print_info ) const
{
  int i;
  GeometryQueryTool* gqt = GeometryQueryTool::instance();

  DLIList<int> updated_ids, created_ids, destroyed_ids;
  DLIList<int> updated_vol_ids, created_vol_ids, destroyed_vol_ids;

    // Remove dead bodies
  input_bodies.reset();
  for (i = input_bodies.size(); i--;)
  {
    Body* body = input_bodies.step_and_get();

    BodySM* bodysm = body->get_body_sm_ptr();
    if (!bodysm)
    {
        // If the body was destroyed, update the list
      destroyed_ids.append( body->id() );

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        if( !temp_vols.get()->get_lump_ptr() )
          destroyed_vol_ids.append( temp_vols.get_and_step()->id() );
      }

      input_bodies.change_to(0);
      gqt->destroy_dead_entity(body);
    }
    else
    {
      remove_dead_entity_names( body );
    }
  }
  gqt->cleanout_deactivated_geometry();

    // Create new bodies
  new_bodies.last();
  for (i = new_bodies.size(); i--; )
  {
    BodySM* bodysm = new_bodies.step_and_get();
    bool newbody = bodysm->owner() == 0;
    Body* body = gqt->make_Body(bodysm);
    result_bodies.append(body);

    if (newbody)
    {
      created_ids.append(body->id());

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        created_vol_ids.append( temp_vols.get_and_step()->id() );
      }
    }
    else
    {
      updated_ids.append(body->id());

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        updated_vol_ids.append( temp_vols.get_and_step()->id() );
      }
    }
  }
  gqt->cleanout_deactivated_geometry();

  if (print_info)
  {
    if( DEBUG_FLAG( 153 ) )
    {
      if (created_ids.size())
         CubitUtil::list_entity_ids( "Created body(s): ", created_ids );
      if (updated_ids.size())
         CubitUtil::list_entity_ids( "Updated body(s): ", updated_ids );
      if (destroyed_ids.size())
         CubitUtil::list_entity_ids( "Destroyed body(s): ", destroyed_ids );
    }

    if (created_vol_ids.size())
       CubitUtil::list_entity_ids( "Created volume(s): ", created_vol_ids );
    if (updated_vol_ids.size())
       CubitUtil::list_entity_ids( "Updated volume(s): ", updated_vol_ids );
    if (destroyed_vol_ids.size())
       CubitUtil::list_entity_ids( "Destroyed volume(s): ", destroyed_vol_ids );
  }

  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : Webcut a body with a cylinder given the input parameters.
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 10/28/97
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::webcut_with_cylinder(
                                       DLIList<Body*>& webcut_body_list,
                                       double radius,
                                       const CubitVector &axis,
                                       const CubitVector &center,
                                       DLIList<Body*>& results_list,
                                       CubitBoolean imprint,
                                       CubitBoolean merge)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  CubitStatus rval = CUBIT_SUCCESS;

  const int count = webcut_body_list.size();
  DLIList<BodySM*> result_sm_list;
  DLIList<Body*> body_list(webcut_body_list);
  DLIList<BodySM*> engine_body_sms(count);
  DLIList<Body*> engine_bodies(count);
  GeometryModifyEngine* gme = 0;

  do_attribute_setup();

  while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
  {
    push_vg_attributes_before_modify(engine_body_sms);

    CubitStatus status = webcut_w_cylinder( engine_body_sms, radius, axis,
              center, result_sm_list, imprint );

    restore_vg_after_modify(result_sm_list, engine_bodies);

    if ( status != CUBIT_FAILURE )
      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list );

    engine_bodies.clean_out();
    engine_body_sms.clean_out();
    result_sm_list.clean_out();

    if ( status == CUBIT_FAILURE )
    {
      rval = CUBIT_FAILURE;
      break;
    }
  }

  do_attribute_cleanup();

  return rval;
}

CubitStatus GeometryModifyTool::webcut_w_cylinder(
                                        DLIList<BodySM*> &webcut_body_list,
                                        double radius,
                                        const CubitVector &axis,
                                        const CubitVector &center,
                                        DLIList<BodySM*>& results_list,
                                        bool imprint )
{
  GeometryModifyEngine* gme = 0;
  gme = get_engine(webcut_body_list.get()); 
  assert(gme);

  double max_size =  0.;
  //lets find the distance to the center for each body and take
  //the max.
  double curr;
  CubitVector cent_bod;
  CubitBox bounding_box;
  BodySM *body_ptr;
  bounding_box = webcut_body_list[0]->bounding_box();
  cent_bod =  bounding_box.center();
  cent_bod = cent_bod - center;
  curr = cent_bod.length();
  if ( curr > max_size )
     max_size = curr;


  for ( int ii = webcut_body_list.size()-1; ii > 0; ii-- )
    {
      body_ptr = webcut_body_list[ii];
      bounding_box |= body_ptr->bounding_box();
      cent_bod = body_ptr->bounding_box().center();
      cent_bod = cent_bod - center;
      curr = cent_bod.length();
      if ( curr > max_size )
        max_size = curr;
    }

  curr = bounding_box.diagonal().length();

  if ( curr > max_size )
     max_size = curr;

  double height = 0.0;
  if ( center.x() > max_size )
    {
      height = 500.0 * center.x();
    }
  else if ( center.y() > max_size )
    {
      height = 500.0 * center.y();
    }
  else if ( center.z() > max_size )
    {
      height = 500.0 * center.z();
    }
  else
    {
      height = 500.0 * max_size;
    }

  //lets make certain we have a valid height..
  if ( height < GEOMETRY_RESABS )
    {
      height = 500.0;
    }

  BodySM *cutting_tool_ptr = gme->cylinder( height, radius, radius, radius );

  if( cutting_tool_ptr == NULL )
    return CUBIT_FAILURE;

  //transform the cyclinder to cernter and axis
  // The current frustum is centered on the z axis.
  CubitVector axis2(0., 0., 1.0 );
  //now find the normal to the current axis and axis we want to be
  //at. This normal is where we will rotate about.
  CubitVector normal_axis = axis2 * axis;
  if ( normal_axis.length() > CUBIT_RESABS )
    {
       //angle in degrees.
       double angle = normal_axis.vector_angle( axis2, axis );
       gme->get_gqe()->rotate(cutting_tool_ptr, normal_axis, angle);
    }
  gme->get_gqe()->translate(cutting_tool_ptr, center);

  CubitStatus stat =
    gme->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

  // Delete the BodySM that was created to be used as a tool
  gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

  return stat;
}

void GeometryModifyTool::do_attribute_setup(void)
{
  //save attribute settings
  CGMApp::instance()->save_current_attribute_states();
  // enable update, actuate, write, read for composite attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  // enable update, actuate, write, read for partition attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_PARTITION_VG, CUBIT_TRUE);
}

void GeometryModifyTool::do_attribute_cleanup(void)
{
  CGMApp::instance()->restore_previous_attribute_states();
}

// Push the virtual geometry attributes down to the underlying solid model
// topology.
void GeometryModifyTool::push_vg_attributes_before_modify(DLIList<BodySM*> &old_sms)
{
  DLIList<TopologyBridge*> top_bridges;
  CAST_LIST_TO_PARENT(old_sms, top_bridges);
  GeometryQueryTool::instance()->ige_export_geom(top_bridges);
}

// Push imprint specific attributes down to the underlying solid model
// topology prior to imprinting.
void GeometryModifyTool::push_imprint_attributes_before_modify(DLIList<BodySM*> &body_sms)
{
  GeometryQueryTool::instance()->ige_push_imprint_attributes_before_modify(body_sms);
}

// Cleanup imprint attributes that were pushed onto the underlying solid
// model topology.
void GeometryModifyTool::remove_imprint_attributes_after_modify(DLIList<BodySM*> &old_sms,
                                                                DLIList<BodySM*> &new_sms)
{
  GeometryQueryTool::instance()->ige_remove_imprint_attributes_after_modify(old_sms, new_sms);
}

#ifdef CAT
CubitStatus GeometryModifyTool::webcut_across_translate(
                                          DLIList<Body*> &webcut_body_list,
                                          RefFace *ref_face_top,
                                          RefFace *ref_face_bottom,
                                          DLIList<Body*> &results_list,
                                          CubitBoolean imprint,
                                          CubitBoolean merge )
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  DLIList<Body*> original_body_list = webcut_body_list;
  const int count = webcut_body_list.size() + 2;
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(webcut_body_list, entity_list);
  entity_list.append( ref_face_top );
  entity_list.append( ref_face_bottom );
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if ( !gme )
  {
     PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                 "different modeling engines is not allowed.\n"
                 "Delete uncommon geometry on these volumes before operation.\n\n");
     return CUBIT_FAILURE;
  }

  Surface* surface_bottom = dynamic_cast<Surface*>(bridge_list.pop());
  Surface* surface_top    = dynamic_cast<Surface*>(bridge_list.pop());
  DLIList<BodySM*> result_sm_list, webcut_sm_list(bridge_list.size());
  CAST_LIST(bridge_list, webcut_sm_list, BodySM);

  do_attribute_setup();
  push_vg_attributes_before_modify(webcut_sm_list);

  CubitStatus result_val = gme->webcut_across_translate (
    webcut_sm_list, surface_top, surface_bottom, result_sm_list, imprint );

  restore_vg_after_modify(result_sm_list, original_body_list);

  result_val = finish_webcut( webcut_body_list, result_sm_list, merge,
                              result_val, results_list);
  do_attribute_cleanup();

  return result_val;
}
#endif

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies using a loop
//                 of curves
//
// Special Notes :
//
// Creator       : Eric W. Nielsen
//
// Creation Date : 12/20/99
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::webcut_with_curve_loop(
                                         DLIList<Body*>& webcut_body_list,
                                         DLIList<RefEdge*>& refedge_list,
                                         DLIList<Body*>& results_list,
                                         CubitBoolean imprint)

{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

   int i;
   DLIList<Body*> original_body_list = webcut_body_list;
   const int count = webcut_body_list.size() + refedge_list.size();
   DLIList<TopologyEntity*> entity_list(count);
   DLIList<TopologyBridge*> bridge_list(count);
   CAST_LIST_TO_PARENT(webcut_body_list, entity_list);
   refedge_list.reset();
   for (i = refedge_list.size(); i--; )
     entity_list.append(refedge_list.get_and_step());
   GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

   if ( !gme )
   {
      PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   DLIList<BodySM*> webcut_sm_list(webcut_body_list.size()), result_sm_list;
   DLIList<Curve*> curve_list(refedge_list.size());
   CAST_LIST(bridge_list, webcut_sm_list, BodySM);
   CAST_LIST(bridge_list, curve_list, Curve);
   assert(webcut_sm_list.size() == webcut_body_list.size());
   assert(curve_list.size() == refedge_list.size());

   GeometryType surface_type = PLANE_SURFACE_TYPE;
   CubitStatus stat;

   //make copies of the curves.
   DLIList<Curve*> copied_curves;
   Curve* temp_curve = NULL;
   for (int i = curve_list.size(); i--;)
     {
       temp_curve = gme->make_Curve(curve_list.get_and_step());
       if(temp_curve != NULL)
         copied_curves.append(temp_curve);
     }

   //make a face out of the curves
   Surface * surf = gme->make_Surface(surface_type, copied_curves, NULL, false );
   if (surf == NULL)
     {
       PRINT_ERROR("webcut tool surface is not created from acis.\n");
       return CUBIT_FAILURE;
     }


   //get cutting tool BodySM.
   BodySM* cutting_tool_ptr = gme->make_BodySM(surf);
   assert(cutting_tool_ptr );

   do_attribute_setup();
   push_vg_attributes_before_modify(webcut_sm_list);

   stat = gme->webcut(
                     webcut_sm_list, cutting_tool_ptr, result_sm_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   restore_vg_after_modify(result_sm_list, original_body_list);

   stat = finish_webcut( webcut_body_list, result_sm_list, CUBIT_FALSE,
                               stat, results_list );
   do_attribute_cleanup();

   return stat;
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies through a plane
//                 defined by a refFace. The newly created bodies are
//                 merged and imprinted depeding on the respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/16/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_surface(
                                      DLIList<Body*>& webcut_body_list,
                                      RefFace* refFace,
                                      DLIList<Body*>& results_list,
                                      CubitBoolean imprint,
                                      CubitBoolean merge)
{
  //make sure that this refface is planar, or you'll get unexpected results
  if( refFace->geometry_type() != PLANE_SURFACE_TYPE )
  {
    PRINT_ERROR("Surface %d is not planar.\n", refFace->id() );
    return CUBIT_FAILURE;
  }

     // Using the face, get three positions (CubitVectors :) :) )
     // The positions can be used by a GeometryModifyEngine to generate its
     // own webcutting tool and perform the webcut.
   CubitVector vector1 = refFace->position_from_u_v(0, 1) ;
   CubitVector vector3 = refFace->position_from_u_v(1, 0) ;
   CubitVector vector2 = refFace->position_from_u_v(0, 0) ;

   CubitStatus result_val = this->webcut_with_plane(webcut_body_list, vector1,
                                  vector2, vector3, results_list, imprint, merge) ;

   return result_val;
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies using a plane
//                 defined by a refFace. The newly created bodies are
//                 merged and imprinted depending on the respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/16/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_body(
                                   DLIList<Body*>& webcut_body_list,
                                   Body* tool_body,
                                   DLIList<Body*>& results_list,
                                   CubitBoolean imprint,
                                   CubitBoolean merge)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  if( tool_body->is_sheet_body() )
    return webcut_with_sheet( webcut_body_list, tool_body, results_list, imprint );

  DLIList<Body*> original_body_list = webcut_body_list;
  DLIList<BodySM*> body_sm_list(webcut_body_list.size()), result_sm_list;
  BodySM* tool_sm = tool_body->get_body_sm_ptr();
  GeometryModifyEngine* gme = 0;
  if (!tool_sm || !(gme = get_engine(tool_sm)))
  {
    PRINT_DEBUG_153("No modify engine available for body %d.\n", tool_body->id());
    PRINT_ERROR("No modify engine available for volume %d.\n", tool_body->ref_volume()->id());
    return CUBIT_FAILURE;
  }

  webcut_body_list.reset();
  for (int i = webcut_body_list.size(); i--; )
  {
    Body* body_ptr = webcut_body_list.get_and_step();
    BodySM* sm_ptr = body_ptr->get_body_sm_ptr();
    if (!sm_ptr || get_engine(sm_ptr) != gme)
    {
      PRINT_DEBUG_153("ERROR: Body %d does not belong to the same geometric engine "
                  "as body %d\n", body_ptr->id(), tool_body->id());
      PRINT_ERROR("Volume %d does not belong to the same geometric engine "
                  "as volume %d\n", body_ptr->ref_volume()->id(), tool_body->ref_volume()->id());
      return CUBIT_FAILURE;
    }

    body_sm_list.append(sm_ptr);
  }

  do_attribute_setup();

  push_vg_attributes_before_modify(body_sm_list);

  int count = gme->webcut(body_sm_list, tool_sm, result_sm_list, imprint);

  restore_vg_after_modify(result_sm_list, original_body_list);

  CubitStatus ret = finish_webcut(webcut_body_list, result_sm_list, merge,
                        count > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE,
                        results_list);

  do_attribute_cleanup();

  return ret;
}

CubitStatus GeometryModifyTool::section(
                                   DLIList<Body*> &section_body_list,
                                   const CubitVector &point_1,
                                   const CubitVector &point_2,
                                   const CubitVector &point_3,
                                   DLIList<Body*> &new_body_list,
                                   CubitBoolean keep_normal_side,
                                   CubitBoolean keep_old )
{
  if (!okay_to_modify( section_body_list, "SECTION" ))
    return CUBIT_FAILURE;
   CubitStatus rval = CUBIT_SUCCESS;

   const int count = section_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(section_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;
   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
      CubitStatus result = gme->section( engine_body_sms,
                                         point_1, point_2, point_3,
                                         result_sm_list,
                                         keep_normal_side, keep_old );
      if (!finish_sm_op( engine_bodies, result_sm_list, new_body_list ))
        result = CUBIT_FAILURE;
      if (!result)
        rval = CUBIT_FAILURE;

      engine_body_sms.clean_out();
      engine_bodies.clean_out();
      result_sm_list.clean_out();
   }

   return rval;
}

CubitStatus
GeometryModifyTool::offset_curves( DLIList<RefEdge*>& ref_edge_list,
                                   double offset_distance,
                                   const CubitVector& offset_direction,
                                   int gap_type )
{
  // Make sure all curves are from same geometry engine
  DLIList<TopologyEntity*> entity_list(ref_edge_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_edge_list.size());
  CAST_LIST_TO_PARENT(ref_edge_list, entity_list);
  GeometryModifyEngine* gme_ptr = common_modify_engine(entity_list, bridge_list);
  if ( !gme_ptr )
  {
    PRINT_ERROR("Can't create offset curves that don't all originate from same\n"
      "       modeling engine.\n");
    return CUBIT_FAILURE;
  }

  DLIList<Curve*> curve_list(bridge_list.size()), result_list;
  CAST_LIST(bridge_list, curve_list, Curve);
  assert(curve_list.size() == ref_edge_list.size());

  CubitStatus rval = gme_ptr->offset_curves( curve_list, result_list,
                        offset_distance, offset_direction, gap_type );
  assert( rval || !result_list.size() );
  result_list.reset();
  for (int i = result_list.size(); i--; )
    GeometryQueryTool::instance()->make_free_RefEdge(result_list.get_and_step());

  return rval;
}

CubitStatus GeometryModifyTool::align_body( Body *body_ptr,
                                            RefFace* my_face,
                                            RefFace *target_face,
                                            CubitVector &my_center,
                                            CubitVector &axis,
                                            CubitVector &target_center,
                                            double &angle)
{
     // Initialize angle to zero, so we are passing back a valid
     // value if we never need to rotate the body
   angle = 0.0;

     //Align the body, which contains my_face, with the target_face.

     //Essentially we want to align the body at my_face with the target_face.
   my_center = my_face->center_point();
   target_center = target_face->center_point();

     //Get the correct volume from my_volume.
   DLIList<RefVolume*> ref_vols;
   body_ptr->ref_volumes(ref_vols);
   DLIList<RefFace*> ref_face_list;
   RefVolume *ref_vol = NULL;

   CubitBoolean found_vol = CUBIT_FALSE;

   for ( int ii = ref_vols.size(); ii > 0; ii-- )
   {
      ref_vol = ref_vols.get_and_step();
      ref_face_list.clean_out();
      ref_vol->ref_faces( ref_face_list );
      if ( ref_face_list.move_to(my_face) )
      {
         found_vol = CUBIT_TRUE;
         break;
      }
   }
   if ( !found_vol )
   {
      PRINT_ERROR( "Cannot find appropriate volume that %s (surface%d) is on.\n",
                   my_face->entity_name().c_str(),
                   my_face->id() );
      PRINT_DEBUG_153( "%s (surface%d) is not on %s (body %d)\n",
                       my_face->entity_name().c_str(),
                       my_face->id(),
                       body_ptr->entity_name().c_str(),
                       body_ptr->id() );
      return CUBIT_FAILURE;
   }
   CubitVector my_vector = my_face->normal_at( my_center, ref_vol );
   CubitVector target_vector = target_face->normal_at( target_center );

   axis = my_vector * target_vector;
     //Move the body so that the center of my_face
     //is at the origin.
   GeometryQueryTool::instance()->translate( body_ptr, -my_center );

     //If the axis length is zero we dont want to do a rotation,
     //they are already parallel..
   if (axis.length() > CUBIT_RESABS )
   {
      angle = axis.vector_angle( my_vector, target_vector );

        //Now we have the angle and the axis to rotate about...
      angle = 180.0*angle/CUBIT_PI;

        //Now rotate the body about the axis.
      GeometryQueryTool::instance()->rotate( body_ptr, axis, angle );
   }
     //Now move the body to the location of the target_center.
   GeometryQueryTool::instance()->translate( body_ptr, target_center );

     //That should do it.
   return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::align_body(Body *body_ptr,
                                           DLIList<RefEntity*>& ref_ent_list,
                                           CubitVector first_vector,
                                           CubitVector second_vector,
                                           CubitVector &my_center_1,
                                           CubitVector &my_center_2,
                                           CubitVector &axis_of_rot,
                                           double &angle,
                                           double &angle_2)
{
   DLIList<RefFace*> ref_face_list;
   DLIList<RefVertex*> ref_vertex_list;

   CAST_LIST(ref_ent_list, ref_face_list, RefFace);
   CAST_LIST(ref_ent_list, ref_vertex_list, RefVertex);

   angle = 0.0;
   angle_2 = 0.0;

   if( ref_face_list.size() == 1 &&
      ref_vertex_list.size() == 0 )
   {

      RefFace *my_face = ref_face_list.get_and_step();
      my_center_1 = my_face->center_point();
      my_center_2 = my_center_1;

      //Get the correct volume from my_volume.
      DLIList<RefVolume*> ref_vols;
      body_ptr->ref_volumes(ref_vols);
      DLIList<RefFace*> face_list;
      RefVolume *ref_vol = NULL;

      CubitBoolean found_vol = CUBIT_FALSE;

      for ( int ii = ref_vols.size(); ii > 0; ii-- )
      {
         ref_vol = ref_vols.get_and_step();
         face_list.clean_out();
         ref_vol->ref_faces( face_list );

         if ( face_list.move_to(my_face) )
         {
            found_vol = CUBIT_TRUE;
            break;
         }
      }

      if ( !found_vol )
      {
         PRINT_ERROR( "Can't find appropriate volume that surface %d is on.\n",
                      my_face->id() );

         PRINT_DEBUG_153( "surface%d is not on body %d\n",
            my_face->id(),
            body_ptr->id() );
         return CUBIT_FAILURE;
      }

      //translate body so surface's centroid is at 0,0,0
      GeometryQueryTool::instance()->translate( body_ptr, -my_center_1 );

      //get axis of rotation
      CubitVector my_vector = my_face->normal_at( my_center_1, ref_vol );
      my_vector.normalize();
      axis_of_rot = first_vector * my_vector;

      //get angle of rotation and rotate
      if (axis_of_rot.length() > CUBIT_RESABS )
      {
         angle = axis_of_rot.vector_angle( my_vector, first_vector );
         //Now we have the angle and the axis to rotate about...
         angle = 180.0*angle/CUBIT_PI;
         //Now rotate the body about the axis.
         GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, angle );
      }
      else  //is aligned already----flip 180 degrees
      {
        axis_of_rot = second_vector;
        angle = 180.0;
        GeometryQueryTool::instance()->rotate( body_ptr, second_vector, angle );
      }

      //translate body back
      GeometryQueryTool::instance()->translate( body_ptr, my_center_1 );

      //That should do it.
      return CUBIT_SUCCESS;
  }
  else if( ref_face_list.size() == 0 &&
     ref_vertex_list.size() == 2 )
  {

     RefVertex *my_vertex_1 = ref_vertex_list.get_and_step();
     RefVertex *my_vertex_2 = ref_vertex_list.get_and_step();
     my_center_1 = my_vertex_1->center_point();
     CubitVector my_center_ver_2 = my_vertex_2->center_point();

     //Get the correct volume from my_volume.
     DLIList<RefVolume*> ref_vols;
     body_ptr->ref_volumes(ref_vols);
     DLIList<RefVertex*> vertex_list;
     RefVolume *ref_vol;

     CubitBoolean found_vol = CUBIT_FALSE;

     for ( int ii = ref_vols.size(); ii > 0; ii-- )
     {
        ref_vol = ref_vols.get_and_step();
        vertex_list.clean_out();
        ref_vol->ref_vertices( vertex_list );

        if ( vertex_list.move_to(my_vertex_1)&&
           vertex_list.move_to(my_vertex_2))
        {
           found_vol = CUBIT_TRUE;
           break;
        }
     }

     if ( !found_vol )
     {
        PRINT_ERROR( "Can't find appropriate volume that vertex %d or %d is on.\n",
                     my_vertex_1->id(),
                     my_vertex_2->id() );

        PRINT_DEBUG_153( "vertex %d or %d is not on body %d\n",
           my_vertex_1->id(),
           my_vertex_2->id(),
           body_ptr->id() );
        return CUBIT_FAILURE;
     }

     CubitVector my_vector = my_center_ver_2 - my_center_1;
     CubitVector target_vector = first_vector;
     axis_of_rot = my_vector * target_vector;

     //Move the body so that the center of the vertex
     //is at the origin.
     GeometryQueryTool::instance()->translate( body_ptr, -my_center_1 );

     my_center_2 = my_vertex_1->center_point();
     CubitVector temp_center_ver_2 = my_vertex_2->center_point();
     my_vector = temp_center_ver_2 - my_center_2;

     my_vector.normalize();
     target_vector.normalize();
     CubitVector test_vec = my_vector - target_vector;

     if (axis_of_rot.length() > CUBIT_RESABS )
     {
        angle = axis_of_rot.vector_angle( my_vector, target_vector );
        //Now we have the angle and the axis to rotate about...
        angle = 180.0*angle/CUBIT_PI;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, angle );
     }
     else if ( test_vec.normalize() > CUBIT_RESABS)
     {
        axis_of_rot = my_center_1 * target_vector;
        angle += 180.0;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, 180.0 );
     }

     CubitVector body_center = body_ptr->center_point();
     body_center.normalize();
     CubitVector temp_vec = first_vector * body_center;

     // room for improvment
     temp_vec.set( temp_vec.x()*temp_vec.x(),
        temp_vec.y()*temp_vec.y(),
        temp_vec.z()*temp_vec.z() );

     body_center += temp_vec;

     if( body_center.x() < -CUBIT_RESABS ||
        body_center.y() < -CUBIT_RESABS ||
        body_center.z() < -CUBIT_RESABS )
     {
        GeometryModifyEngine *gePtr1 = get_engine(body_ptr);
        DLIList<Surface*> bad_face_list;
        body_ptr->get_body_sm_ptr()->surfaces(bad_face_list);
        angle_2= 180.0;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, first_vector, angle_2 );
        if (bad_face_list.size() && body_ptr->is_sheet_body())
        {
           CubitStatus flip_result = gePtr1->flip_normals(bad_face_list);
           PRINT_WARNING("Reversing all surface normals within sheet body %d \n"
                         "to remain in the positive plane\n",body_ptr->id());
           if ( flip_result == CUBIT_FAILURE )
           {
              return CUBIT_FAILURE;
           }
        }
     }
     return CUBIT_SUCCESS;
  }
  else
  {
     PRINT_ERROR("Incorrect syntax type 'help align'. Nothing aligned\n");
     return CUBIT_FAILURE;
  }
}


CubitStatus GeometryModifyTool::sweep_setup(
                         const char* name,
                         DLIList<RefEntity*>& input_entity_list,
                         DLIList<Body*>& output_body_list,
                         GeometryModifyEngine*& output_engine,
                         CubitBoolean& changed_new_ids,
                         DLIList<GeometryEntity*>& output_geom_list,
                         DLIList<RefEdge*>* input_edge_list,
                         DLIList<Curve*>* output_curve_list )
{
  int i;

  if (input_entity_list.size() == 0)
    return CUBIT_FAILURE;

  DLIList<RefFace*> ref_face_list;
  DLIList<RefEdge*> ref_edge_list;
  CAST_LIST( input_entity_list, ref_face_list, RefFace );
  CAST_LIST( input_entity_list, ref_edge_list, RefEdge );
  if( ref_face_list.size() && ref_edge_list.size() )
  {
    PRINT_ERROR( "Currently cannot sweep curves and surfaces at the same time\n" );
    return CUBIT_FAILURE;
  }

  if( (ref_face_list.size()+ref_edge_list.size()) != input_entity_list.size() )
  {
    PRINT_ERROR( "Can only sweep surfaces or curves\n" );
    return CUBIT_FAILURE;
  }

  // Currently faces to be swept have to belong to different bodies.  It
  // doesn't work well to sweep multiple faces from the same body (also
  // current implementation in AGE crashes).
  if( ref_face_list.size() )
  {
    int free_surf_cnt = 0; // Number of free surfaces being swept
    DLIList<Body*> body_list;
    RefFace *ref_face_ptr;
    for( i=ref_face_list.size(); i--; )
    {
      ref_face_ptr = ref_face_list.get_and_step();
      DLIList<Body*> tmp_body_list;
      ref_face_ptr->bodies( tmp_body_list );
      if( tmp_body_list.size() > 1 )
      {
        PRINT_ERROR( "Surface %d belongs to more than one volume; cannot sweep.\n",
          ref_face_ptr->id() );
        return CUBIT_FAILURE;
      }

      if( tmp_body_list.size() > 0 )
        body_list.append_unique( tmp_body_list.get() );
      else
        free_surf_cnt++;
    }

    if( body_list.size()+free_surf_cnt != ref_face_list.size() )
    {
      PRINT_ERROR( "The surfaces to be swept cannot belong to the same volume\n" );
      return CUBIT_FAILURE;
    }

    output_body_list = body_list;
  }

  // Make sure all entities are from same modify engine
  DLIList<TopologyEntity*> te_list;
  CAST_LIST(input_entity_list, te_list, TopologyEntity);

  // This is for sweeping along a curve list
  if (input_edge_list)
  {
    input_edge_list->reset();
    for( i=input_edge_list->size(); i--; )
      te_list.append( input_edge_list->get_and_step() );
  }

  DLIList<TopologyBridge*> bridge_list( te_list.size() );
  output_engine = common_modify_engine( te_list, bridge_list );
  if( output_engine == NULL )
  {
     PRINT_ERROR("Can't sweep with entities from different modeling engines.\n");
     return CUBIT_FAILURE;
  }

  // Check for virtual anywhere in participating bodies.  This will catch cases where
  // the result from sweeping would interact with some virtual somewhere. bwc 12/22/05.
  for(i=output_body_list.size(); i--;)
  {
     if(GeometryQueryTool::instance()->contains_intermediate_geometry(
        output_body_list.get_and_step()))
     {
        PRINT_ERROR("Can't sweep faces of bodies containing virtual geometry.\n");
        return CUBIT_FAILURE;
     }
  }

  changed_new_ids = CUBIT_FALSE;

  TopologyBridge* bridge_ptr;
  bridge_list.reset();
  for( i=input_entity_list.size(); i--; )
  {
    bridge_ptr = bridge_list.get_and_step();
    GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge_ptr);
    output_geom_list.append(geom_ptr);
  }

  if( input_edge_list )
  {
    for( i=input_edge_list->size(); i--; )
    {
      bridge_ptr = bridge_list.get_and_step();
      Curve* curve_ptr = dynamic_cast<Curve*>(bridge_ptr);
      output_curve_list->append(curve_ptr);
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::sweep_finish(
                                     const char* const,
                                     DLIList<Body*>& input_body_list,
                                     DLIList<BodySM*>& new_body_list,
                                     CubitBoolean changed_new_ids )
{
  int i;
  DLIList<BodySM*> regen_list( new_body_list );
  GeometryQueryTool *gqt = GeometryQueryTool::instance();

  input_body_list.reset();
  for (i = input_body_list.size(); i--; )
  {
    Body* body = input_body_list.get_and_step();
    BodySM* bodysm = body->get_body_sm_ptr();
    if (bodysm)
    {
      regen_list.append_unique(bodysm);
      remove_dead_entity_names(body);
    }
    else
    {
      gqt->destroy_dead_entity(body);
    }
  }
  gqt->cleanout_deactivated_geometry();

  regen_list.reset();
  for (i = regen_list.size(); i--; )
  {
    BodySM* bodysm = regen_list.get_and_step();
    Body* body = gqt->make_Body(bodysm);
    PRINT_INFO("%s volume %d\n",
      new_body_list.is_in_list(bodysm) ? "Created swept" : "Updated",
      body->ref_volume()->id());
  }

  if (changed_new_ids) set_new_ids(CUBIT_FALSE);
  gqt->cleanout_deactivated_geometry();

  return CUBIT_SUCCESS;
}





CubitStatus GeometryModifyTool::sweep_rotational(
                                        DLIList<RefEntity*>& ref_ent_list,
                                        const CubitVector& point,
                                        const CubitVector& sweep_axis,
                                        double angle,
                                        int steps,
                                        double draft_angle,
                                        int draft_type,
                                        CubitBoolean switchside,
                                        CubitBoolean make_solid,
                                        CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("rotational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1-> sweep_rotational( geom_list,
                                                   result_list,
                                                   point,
                                                   sweep_axis,
                                                   angle,
                                                   steps,
                                                   draft_angle,
                                                   draft_type,
                                                   switchside,
                                                   make_solid,
                                                   rigid);

  if (!sweep_finish("rotational", body_list, result_list, change_newids))
    status = CUBIT_FAILURE;

  body_list.clean_out();
  for(int i = 0; i < result_list.size(); i++)
  {
    Body* body = CAST_TO(result_list.get_and_step()->topology_entity(),Body );
    if(body)
      body_list.append(body);
  }
  CAST_LIST( body_list, ref_ent_list, RefEntity);
  return status;
}

CubitStatus GeometryModifyTool::sweep_translational(
                                         DLIList<RefEntity*>& ref_ent_list,
                                         const CubitVector& sweep_vector,
                                         double draft_angle,
                                         int draft_type,
                                         CubitBoolean switchside,
                                         CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1->
    sweep_translational( geom_list,
                         result_list,
                         sweep_vector,
                         draft_angle,
                         draft_type,
                         switchside,
                         rigid);

  if (!sweep_finish("translational", body_list, result_list, change_newids))
    status = CUBIT_FAILURE;

  body_list.clean_out();
  for(int i = 0; i < result_list.size(); i++)
  {
    Body* body = CAST_TO(result_list.get_and_step()->topology_entity(),Body ); 
    if(body)
      body_list.append(body);
  }
  CAST_LIST( body_list, ref_ent_list, RefEntity);
  return status;
}

CubitStatus GeometryModifyTool::sweep_target(CubitPlane ref_plane,
											 DLIList<RefEntity*>& ref_ent_list)
{
	double distance1;
	double distance2;
	double distance3;
	double temp_counter=0;
	CubitVector begin_point;
	CubitVector target_begin_point;
	CubitVector get_mid_point;
	CubitVector target_mid_point;
	CubitVector end_point;
	CubitVector target_end_point;
	DLIList<RefEdge*> edge_list;
	CAST_LIST(ref_ent_list, edge_list, RefEdge);	
	CubitVector result;
	double max_distance_for_edge=0;
	CubitVector max_result;
	//make sure that there are actually edges in list
	if(edge_list.size() > 0)
	{
		//this part of the code steps through all edges that could be
		//selected by user and finds the largest distance of all edges from three
		//points: midpoint, endpoint, and beginpoint
		for (int i=0;i<edge_list.size();i++)
		{
			//find the midpoint vector of the edge
			edge_list[i]->mid_point(get_mid_point);

			//Project the midpoint onto the specified plane
			target_mid_point = ref_plane.project(get_mid_point);

			//Calculate the distance between the mid_point, and target_point
			distance1 = target_mid_point.distance_between(get_mid_point);

			//the next two blocks are just copies of the above three steps
			double fraction_along_curve = 1;
			edge_list[i]->position_from_fraction(fraction_along_curve, end_point);
			target_end_point = ref_plane.project(end_point);
			distance2 = target_end_point.distance_between(end_point);

			fraction_along_curve = 0;
			edge_list[i]->position_from_fraction(fraction_along_curve, begin_point);
			target_begin_point = ref_plane.project(begin_point);
			distance3 = target_begin_point.distance_between(begin_point);

			//see which of the three distances is greater for that edge
			//and compare its distance to the other edges that have already been calculated
			//saving the vector of the largest distance
			if (distance1>distance2 && distance1>distance3)
			{
				result = 2*(target_mid_point - get_mid_point);
				max_distance_for_edge=distance1;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
			else if (distance2>distance3)
			{
				result = 2*(target_end_point - end_point);
				max_distance_for_edge=distance2;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
			else
			{
				result = 2*(target_begin_point - begin_point);
				max_distance_for_edge=distance3;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
		}

		//using the 'facet edges' defined by the drawing geometry take the leading
		//and trailing edge projection vectors of this small edge.  With that,
		//calc the dot product and if negative, the vectors are opposite of each
		//other meaning the curve travels through the plane and the sweep function will fail
		for (int ii=0;ii<edge_list.size();ii++)
		{
			Curve *facet_curve;
			facet_curve=edge_list[ii]->get_curve_ptr();
			int num_points;
			//int color = 2;
			CubitStatus response;
			GMem g_mem;

			//get number of points and their locations
			//on the curve as defined by the drawing geometry algorithm
			response = facet_curve->get_geometry_query_engine()->
				get_graphics( facet_curve, num_points, &g_mem );

			if (response==CUBIT_FAILURE || num_points == 0)
			{
				PRINT_WARNING("Facet Curve calling function failed\n" );
			}

			GPoint *point_data = g_mem.point_list();

			for (int jj=0; jj<= (num_points-1); jj++)
			{
				//get first points vectors
				CubitVector point_1;
				point_1.x(point_data[jj].x);
				point_1.y(point_data[jj].y);
				point_1.z(point_data[jj].z);
				//get second points vectors
				CubitVector point_2;
				point_2.x(point_data[jj+1].x);
				point_2.y(point_data[jj+1].y);
				point_2.z(point_data[jj+1].z);
				//project the two points onto target plane
				CubitVector target_point_1;
				target_point_1 = ref_plane.project(point_1);
				CubitVector target_point_2;
				target_point_2 = ref_plane.project(point_2);
				//calc vector from point on curve to point on surface
				CubitVector vec_1 = point_1 - target_point_1;
				CubitVector vec_2 = point_2 - target_point_2;
				//make them unit vectors
				vec_1.normalize();
				vec_2.normalize();
				//calculate dot product
				double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();
				//check to see if dot product sign is zero
				//the and statement checks for if the first or last vertex of the edge
				//which may be sitting on the surface, this is alright
				//and should still be able to sweep
				if (dot<0 && (jj+1!=num_points || jj==0))
				{
					PRINT_ERROR( "The edge is traveling through the plane\n" );
					PRINT_ERROR( "and thus forces the sweep direction to switch\n" );
					PRINT_ERROR( "direction. Try splitting the edge.\n" );
					return CUBIT_FAILURE;
				}
			}
		}
		DLIList<BodySM*> webcut_results_list;
		DLIList<Body*> body_list;
		DLIList<BodySM*> sweep_result_list;
		DLIList<GeometryEntity*> geom_list;
		GeometryModifyEngine* gePtr1 = 0;
		CubitBoolean change_newids;

		if (!sweep_setup("translational", ref_ent_list , body_list, gePtr1,
			change_newids, geom_list))
			return CUBIT_FAILURE;

		//below block is default settings to be fed into sweep_translational
		CubitBoolean rigid = CUBIT_FALSE;
		CubitBoolean switchside = CUBIT_FALSE;
		int draft_type = 0;
		double draft_angle=0.0;

		//sweep the curve through and hopefully past the target surface
		CubitStatus status = gePtr1->sweep_translational( geom_list,sweep_result_list,
			max_result,draft_angle, draft_type,switchside,rigid);

		if (status == 0)
		{
			PRINT_ERROR( "Sweep operation failed!\n" );
			//delete sweep_result_list memory
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			return CUBIT_FAILURE;
		}

		//below lines are used to define the ref_plane with three points for webcut input
		//find normal vector to user specified plane
		CubitVector vec1 = ref_plane.normal();
		CubitVector vec2;
		CubitVector vec3;
		//get orthogonal vectors of the normal vector
		vec1.orthogonal_vectors(vec2, vec3);
		//place the orthogonal vectors at a point on the plane as opposed to the origin
		vec2=vec2+target_mid_point;
		vec3=vec3+target_mid_point;

		//do a webcut with a plane created from the three projected points above
		CubitStatus status2 = gePtr1->webcut(sweep_result_list,target_mid_point,
			vec2,vec3,webcut_results_list);

		if (status2 == 0)
		{
			PRINT_ERROR( "Sweep operation worked; however, webcut operation failed.\n" );
			//delete memory since it failed
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);
			return CUBIT_FAILURE;
		}

		if (webcut_results_list.size()<=0)
		{
			PRINT_ERROR( "Number of bodies from webcut is zero, unable to perform rest of sweep operation\n" );
			//delete memory since it failed
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);
			return CUBIT_FAILURE;
		}

		CubitVector created_curve_mid;
		CubitVector user_mid;
		DLIList<BodySM*> keep_bodies_list;
		DLIList<Curve*> curve_list;

		//generate a list of curves from the webcut_results_list
		for (int counter = 0;counter < webcut_results_list.size(); counter++)
		{
			BodySM* OO = webcut_results_list[counter];
			OO->curves(curve_list);
		}

		//step through each of the user specified edges
		for (int edge_counter = 0; edge_counter < edge_list.size(); edge_counter++)
		{
			//find the midpoint of a user specified edge
			edge_list[edge_counter]->position_from_fraction(.5,user_mid);
			double min_dist=DBL_MAX;
			double distance;
			Curve* closest_curve=0;

			//step through all of the curves
			for (int counter2 = 0; counter2 < curve_list.size(); counter2++)
			{
				//find the midpoint of the created curve
				curve_list[counter2]->position_from_fraction(.5,created_curve_mid);
				//calculate the distance between the two midpoints
				distance=created_curve_mid.distance_between(user_mid);

				//find the minimum distance between all the curves
				if (distance<min_dist)
				{
					//reset the min_distance
					min_dist=distance;
					//keep track of which curve is associated with the shortest distance
					closest_curve=curve_list[counter2];
				}
			}
			//find the parent body of the minimum distance curve
			BodySM* body_sm = closest_curve->bodysm();
			//append that body to a keep list
			keep_bodies_list.append(body_sm);
		}
		//delete bodies which repeat (useful only when in granite engine)
		keep_bodies_list.uniquify_unordered();
		//subtract the keep_list from the webcut_results_list
		webcut_results_list -=keep_bodies_list;
		//delete webcut_results_list memory
		gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);

		//builds ref bodies
		if (!sweep_finish("translational", body_list, keep_bodies_list, change_newids))
		{
			gePtr1->get_gqe()->delete_solid_model_entities(keep_bodies_list);
			status = CUBIT_FAILURE;
		}
	}
	else
	{
		PRINT_ERROR( "No edge(s) found - sweep creation failed.\n" );
		return CUBIT_FAILURE;
	}

	return CUBIT_SUCCESS; 
}

CubitStatus GeometryModifyTool::sweep_surface_target(CubitPlane ref_plane,
													 DLIList<RefEntity*>& ref_ent_list)
{
	DLIList<RefFace*> surface_list;
	CAST_LIST(ref_ent_list, surface_list, RefFace);
	double distance1;
	double distance2;
	double distance3;
	double temp_counter=0;
	CubitVector begin_point;
	CubitVector target_begin_point;
	CubitVector get_mid_point;
	CubitVector target_mid_point;
	CubitVector end_point;
	CubitVector target_end_point;
	CubitVector result;
	double max_distance_for_edge=0;
	CubitVector max_result;
	DLIList<RefEdge*> surf_edge_list;

	//make sure that only one surface has been selected
	if(surface_list.size() == 0)
	{
		PRINT_ERROR( "No edge(s) found - sweep surface to target failed.\n" );
		return CUBIT_FAILURE;
	}

	//make sure that there are actually surfaces in list
	else if(surface_list.size() > 1)
	{
		PRINT_ERROR( "You can only use this operation \n with one surface selected at a time\n" );
		return CUBIT_FAILURE;
	}
	else
	{
		//this part of the code steps through all edges that could be
		//selected by user and finds the largest distance of all edges from three
		//points: midpoint, endpoint, and beginpoint
		for (int i=0;i<surface_list.size();i++)
		{
			surface_list[i]->ref_edges(surf_edge_list);
			for (int j=0;j<surf_edge_list.size();j++)
			{
				//get midpoint of edge on surface
				surf_edge_list[j]->mid_point(get_mid_point);
				//Project the midpoint of each surface edge onto the specified plane
				target_mid_point = ref_plane.project(get_mid_point);		

				//Calculate the distance between the mid_point, and target_point
				distance1 = target_mid_point.distance_between(get_mid_point);

				//the next two blocks are just copies of the above three steps

				int fraction_along_curve = 1;
				surf_edge_list[j]->position_from_fraction(fraction_along_curve, end_point);
				target_end_point = ref_plane.project(end_point);
				distance2 = target_end_point.distance_between(end_point);

				fraction_along_curve = 0;
				surf_edge_list[j]->position_from_fraction(fraction_along_curve, begin_point);
				target_begin_point = ref_plane.project(begin_point);
				distance3 = target_begin_point.distance_between(begin_point);

				//see which of the three distances is greater of the edge
				//and compare its distance to the other edges
				if (distance1>distance2 && distance1>distance3)
				{
					result = 2*(target_mid_point - get_mid_point);
					max_distance_for_edge=distance1;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				else if (distance2>distance3)
				{
					result = 2*(target_end_point - end_point);
					max_distance_for_edge=distance2;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				else
				{
					result = 2*(target_begin_point - begin_point);
					max_distance_for_edge=distance3;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				//just checking to make sure the user didn't specify a surface as both
				//a target surface and sweeping surface which would result in a CUBIT crash
				if (i==surface_list.size()-1 && j==surf_edge_list.size()-1 &&
					max_distance_for_edge < 0.000000000000001)
				{
					PRINT_ERROR( "The sweep distance is less than the geometry tolerance\n" );
					PRINT_ERROR( "This may be caused by selecting the same\n");
					PRINT_ERROR( "sweep surface and target surface\n" );
					return CUBIT_FAILURE;
				}
			}
		}
		//using the facet edges defined by the drawing geometry code; take the leading
		//and trailing edge vectors of each facet edge;  With that, calc the dot product
		//and if negative (meaning the vector directions switched) the curve travels
		//through the plane and the sweep function will fail
		surface_list.reset();
		surf_edge_list.clean_out();
		for (int kk=0;kk<surface_list.size();kk++)
		{
			surface_list[kk]->ref_edges(surf_edge_list);
			for (int ii=0;ii<surf_edge_list.size();ii++)
			{
				Curve *facet_curve;
				facet_curve=surf_edge_list[ii]->get_curve_ptr();
				int num_points;
				//int color = 2;
				CubitStatus response;
				GMem g_mem;

				//get number of points and their locations
				//on the curve as defined by the drawing geometry algorithm
				response = facet_curve->get_geometry_query_engine()->
					get_graphics( facet_curve, num_points, &g_mem );

				if (response==CUBIT_FAILURE || num_points == 0)
				{
					PRINT_WARNING("Unable to preview a curve\n" );
				}

				GPoint *point_data = g_mem.point_list();

				for (int jj=0; jj<= (num_points-1); jj++)
				{
					//get first points vectors
					CubitVector point_1;
					point_1.x(point_data[jj].x);
					point_1.y(point_data[jj].y);
					point_1.z(point_data[jj].z);
					//get second points vectors
					CubitVector point_2;
					point_2.x(point_data[jj+1].x);
					point_2.y(point_data[jj+1].y);
					point_2.z(point_data[jj+1].z);
					//project the two points onto target plane
					CubitVector target_point_1;
					target_point_1 = ref_plane.project(point_1);
					CubitVector target_point_2;
					target_point_2 = ref_plane.project(point_2);
					//calc vector from point on curve to point on surface
					CubitVector vec_1 = point_1 - target_point_1;
					CubitVector vec_2 = point_2 - target_point_2;
					//make them unit vectors
					vec_1.normalize();
					vec_2.normalize();
					//double dot = DotProduct(vec_1,vec_2);
					//calculate dot product
					double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();

					//check to see if dot product sign is zero
					//the and statement checks for if the first or last vertex of the edge
					//which may be sitting on the surface, this is alright
					//and should still be able to sweep
					if (dot<0 && (jj+1!=num_points || jj==0))
					{
						PRINT_ERROR( "The surface is traveling through the plane\n" );
						PRINT_ERROR( "and thus forces the sweep direction to switch\n" );
						PRINT_ERROR( "direction. Try splitting the surface.\n" );
						return CUBIT_FAILURE;
					}
				}
			}
		}

		DLIList<Body*> body_list;
		DLIList<GeometryEntity*> geom_list;
		GeometryModifyEngine* gePtr1 = 0;
		CubitBoolean change_newids;

		if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
			change_newids, geom_list))
			return CUBIT_FAILURE;

		DLIList<BodySM*> sweep_result_list;

		//below block is default settings to be fed into sweep_translational
		CubitBoolean rigid = CUBIT_FALSE;
		CubitBoolean switchside = CUBIT_FALSE;
		int draft_type = 0;
		double draft_angle=0.0;

		//find normal vector to user specified plane
		CubitVector vec1 = ref_plane.normal();
		CubitVector vec2;
		CubitVector vec3;
		//get orthogonal vectors of the normal vector
		vec1.orthogonal_vectors(vec2, vec3);
		//place the orthogonal vectors at a point on the plane as opposed to the origin
		vec2=vec2+target_mid_point;
		vec3=vec3+target_mid_point;
		DLIList<BodySM*> webcut_results_list;
		//get a point on the first user specified surface before it gets consumed in the sweep_translational function
		CubitVector mid_point_surface;
		mid_point_surface = surface_list[0]->center_point();

		//sweep the curve down through and hopefully past the target surface
		CubitStatus status = gePtr1->sweep_translational( geom_list,sweep_result_list,
			max_result,draft_angle, draft_type,switchside,rigid);

		if (status == 0)
		{
			//If in here, sweep_translational failed so delete result_list and
			//print an error to the screen for the user
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			PRINT_ERROR( "Error occured in the sweep operation.\n" );
			return CUBIT_FAILURE;
		}

		//do a webcut with a plane created from the three projected points
		CubitStatus status2 = gePtr1->webcut(sweep_result_list,target_mid_point,
			vec2,vec3,webcut_results_list);

		if (status2 == 0)
		{
			//If in here, webcut operation failed so delete result_list and
			//print an error to the screen for the user
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);
			PRINT_ERROR( "Error occured in the webcut operation.\n" );
			return CUBIT_FAILURE;
		}

		DLIList<BodySM*> keep_bodies_list = webcut_results_list;
		int bodies_deleted_counter=0;
		CubitVector target_mid_point_volume;
		CubitVector mid_point_volume;
		CubitVector target_mid_point_surface;
		double volume_size;

		//project the mid_point of user specified surface onto target plane
		target_mid_point_surface = ref_plane.project(mid_point_surface);
		//calculate vector between the two points and then normalize
		CubitVector vec_1 = mid_point_surface - target_mid_point_surface;
		vec_1.normalize();

		//step through the 
		for (int counter=0;counter<webcut_results_list.size();counter++)
		{
			//find the geometric midpoint of the body and project that point on the target plane
			webcut_results_list[counter]->mass_properties(mid_point_volume,volume_size);
			target_mid_point_volume = ref_plane.project(mid_point_volume);
			//generate a vector between the two points and then normalize
			CubitVector vec_2 = mid_point_volume - target_mid_point_volume;
			vec_2.normalize();

			//calculate the dot product
			double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();

			//if a negative dot product delete it because it is on the opposite side of the target plane
			if (dot < 0)
			{
				keep_bodies_list.remove(webcut_results_list[counter]);
				bodies_deleted_counter = bodies_deleted_counter + 1;
			}

		}

		//test to see if all bodies have been deleted and if so let the user know
		if (bodies_deleted_counter == webcut_results_list.size())
		{
			PRINT_ERROR( "All sweeped surfaces deleted - sweep_target failed.\n" );
			PRINT_ERROR( "This may be due to granite engine limitations and/or\n" );
			PRINT_ERROR( "angle between curve and target surface\n" );
			return CUBIT_FAILURE;
		}

		//delete webcut_results_list since it is no longer of use
		webcut_results_list -= keep_bodies_list;
		gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);

		//builds ref bodies
		if (!sweep_finish("translational", body_list, keep_bodies_list, change_newids))
			status = CUBIT_FAILURE;

	}
	return CUBIT_SUCCESS;
}


CubitStatus GeometryModifyTool::sweep_perpendicular( DLIList<RefEntity*>& ref_ent_list,
                                                     double distance,
                                                     double draft_angle,
                                                     int draft_type,
                                                     CubitBoolean switchside,
                                                     CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1->
    sweep_perpendicular( geom_list,
                         result_list,
                         distance,
                         draft_angle,
                         draft_type,
                         switchside,
                         rigid);

  if (!sweep_finish("perpendicular", body_list, result_list, change_newids))
    status = CUBIT_FAILURE;

  body_list.clean_out();
  for(int i = 0; i < result_list.size(); i++)
  {
    Body* body = CAST_TO(result_list.get_and_step()->topology_entity(),Body );
    if(body)
      body_list.append(body);
  }
  CAST_LIST( body_list, ref_ent_list, RefEntity);
  return status;
}
CubitStatus GeometryModifyTool::sweep_along_curve(DLIList<RefEntity*>& ref_ent_list,
                                                  DLIList<RefEdge*>& ref_edge_list,
                                                  double draft_angle,
                                                  int draft_type,
                                                  CubitBoolean rigid)

{
   DLIList<GeometryEntity*> geom_list(ref_ent_list.size());
   DLIList<Curve*> curve_list(ref_edge_list.size());
   DLIList<Body*> body_list(ref_ent_list.size());
   GeometryModifyEngine* engine_ptr = 0;
   CubitBoolean changed_new_ids = CUBIT_FALSE;
   CubitStatus status = sweep_setup( "along_curve",
                                     ref_ent_list,
                                     body_list,
                                     engine_ptr,
                                     changed_new_ids,
                                     geom_list,
                                     &ref_edge_list,
                                     &curve_list );
   if (status != CUBIT_SUCCESS)
    return status;



   DLIList<BodySM*> result_list;
   status = engine_ptr->sweep_along_curve( geom_list,
                                           result_list,
                                           curve_list,
                                           draft_angle,
                                           draft_type,
                                           rigid);

   if (!sweep_finish("along_curve", body_list, result_list, changed_new_ids))
    status = CUBIT_FAILURE;

   body_list.clean_out();
   for(int i = 0; i < result_list.size(); i++)
   {
     Body* body = CAST_TO(result_list.get_and_step()->topology_entity(),Body );
     if(body)
       body_list.append(body);
   }
   CAST_LIST( body_list, ref_ent_list, RefEntity);
   return status;

}
void GeometryModifyTool::initialize_settings() {


  SettingHandler::instance()->add_setting("Group Imprint",
					  GeometryModifyTool::set_group_imprint,
					  GeometryModifyTool::get_group_imprint);

  SettingHandler::instance()->add_setting("NonRegImprint",
					  GeometryModifyTool::set_all_edges_imprint,
					  GeometryModifyTool::get_all_edges_imprint);

  SettingHandler::instance()->add_setting("New Ids",
					  GeometryModifyTool::set_new_ids,
					  GeometryModifyTool::get_new_ids);

  SettingHandler::instance()->add_setting("Separate After Webcut",
					  GeometryModifyTool::set_sep_after_webcut_setting,
					  GeometryModifyTool::get_sep_after_webcut_setting);

  SettingHandler::instance()->add_setting("old names",
            GeometryModifyTool::set_old_names,
            GeometryModifyTool::get_old_names);

}



// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : Constructor of the GeometryModifyTool class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryModifyTool::GeometryModifyTool(GeometryModifyEngine*gme_ptr)
{
  if (gme_ptr) gmeList.append(gme_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies through a plane.
//                 The newly created bodies are and merged depeding on the
//                 respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_brick(
                                           DLIList<Body*>& webcut_body_list,
                                           const CubitVector &center,
                                           const CubitVector axes[3],
                                           const CubitVector &extension,
                                           DLIList<Body*> &results_list,
                                           CubitBoolean imprint,
                                           CubitBoolean merge )
{
   // Make sure that entity creation is possible.  Allow at most
   // only one of the dimensions to be zero - in which case a planar
   // sheet is used to do the cutting.
   double width = 2.0*extension.x();
   double height = 2.0*extension.y();
   double depth = 2.0*extension.z();
   if ( width < 0.0 || height < 0.0 || depth < 0.0 )
   {
      PRINT_ERROR( "Cannot make a brick of size %f x %f x %f\n"
                   "     Negative dimensions are not allowed.\n",
                   width, height, depth );
      return CUBIT_FAILURE ;
   }
   
   BodySM *cutting_tool_ptr = NULL;
   CubitBoolean is_sheet_body = CUBIT_FALSE;
   CubitVector p1, p2, p3, p4;
   int wz = width < GEOMETRY_RESABS;
   int hz = height < GEOMETRY_RESABS;
   int dz = depth < GEOMETRY_RESABS;
   int num_zero_dim = wz + hz + dz;
   if( num_zero_dim > 0 )
   {
      if( num_zero_dim > 1 )
      {
         PRINT_ERROR( "Cannot make a sheet of size %f x %f x %f\n"
            "     At least two dimensions must be nonzero.\n",
            width, height, depth );
         return CUBIT_FAILURE ;
      }

      // Make a sheet body instead of a cuboid
      is_sheet_body = CUBIT_TRUE;
      CubitVector sheet_axes[2];
      if( wz )
      {
         sheet_axes[0] = axes[1];
         sheet_axes[1] = axes[2];
         width = depth;
      }
      else if( hz )
      {
         sheet_axes[0] = axes[2];
         sheet_axes[1] = axes[0];
         height = depth;
      }
      else
      {
         sheet_axes[0] = axes[0];
         sheet_axes[1] = axes[1];
      }

      // Create the planar sheet to cut with
      // Get the corners of the sheet
      center.next_point( axes[0], width/2.0, p1 );
      p1.next_point( axes[1], -height/2.0, p1 );
      p1.next_point( axes[1], height, p2 );
      p2.next_point( axes[0], -width, p3 );
      p3.next_point( axes[1], -height, p4 );
   }

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   CubitStatus rval = CUBIT_SUCCESS;

   const int count = webcut_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(webcut_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;

   do_attribute_setup();

   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
      push_vg_attributes_before_modify(engine_body_sms);

      // Create the brick to cut with
      if (is_sheet_body)
	cutting_tool_ptr = gme->planar_sheet(p1,p2,p3,p4);
      else
        cutting_tool_ptr = gme->brick( center, axes, extension );
      if( cutting_tool_ptr == NULL )
         return CUBIT_FAILURE;

      CubitStatus status = gme->webcut (
        engine_body_sms, cutting_tool_ptr, result_sm_list, imprint );

      // Delete the BodySM that was created to be used as a tool
      gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

      restore_vg_after_modify(result_sm_list, engine_bodies);

      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list );
      if (!status)
        rval = CUBIT_FAILURE;

      engine_bodies.clean_out();
      engine_body_sms.clean_out();
      result_sm_list.clean_out();
   }

   do_attribute_cleanup();

   return rval;
}

CubitStatus GeometryModifyTool::webcut_with_planar_sheet(
                                           DLIList<Body*>& webcut_body_list,
                                           const CubitVector &center,
                                           const CubitVector axes[2],
                                           double width, double height,
                                           DLIList<Body*> &results_list,
                                           CubitBoolean imprint,
                                           CubitBoolean merge )
{
   if ( width <= GEOMETRY_RESABS || height <= GEOMETRY_RESABS )
   {
      PRINT_ERROR( "Cannot webcut with a sheet of size %f x %f\n"
                   "     Negative or zero dimensions are not allowed.\n",
                   width, height );
      return CUBIT_FAILURE ;
   }

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   CubitStatus rval = CUBIT_SUCCESS;

   const int count = webcut_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(webcut_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;
   do_attribute_setup();
   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
      push_vg_attributes_before_modify(engine_body_sms);

      // Create the planar sheet to cut with
      CubitVector p1, p2, p3, p4;

      // Get the corners of the sheet
      center.next_point( axes[0], width/2.0, p1 );
      p1.next_point( axes[1], -height/2.0, p1 );
      p1.next_point( axes[1], height, p2 );
      p2.next_point( axes[0], -width, p3 );
      p3.next_point( axes[1], -height, p4 );

      BodySM *cutting_tool_ptr = gme->planar_sheet(p1,p2,p3,p4);
      if( cutting_tool_ptr == NULL )
         return CUBIT_FAILURE;

      CubitStatus status = gme->webcut (
        engine_body_sms, cutting_tool_ptr, result_sm_list, imprint );
 
      // Delete the BodySM that was created to be used as a tool
      gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

      restore_vg_after_modify(result_sm_list, engine_bodies);

      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list );
      if (!status)
        rval = CUBIT_FAILURE;

      engine_bodies.clean_out();
      engine_body_sms.clean_out();
      result_sm_list.clean_out();
   }

   do_attribute_cleanup();

   return rval;
}

CubitStatus GeometryModifyTool::webcut_with_plane(
                                    DLIList<Body*>& webcut_body_list,
                                    const CubitVector &vector1,
                                    const CubitVector &vector2,
                                    const CubitVector &vector3,
                                    DLIList<Body*>& results_list,
                                    CubitBoolean imprint,
                                    CubitBoolean merge)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  CubitStatus rval = CUBIT_SUCCESS;

  const int count = webcut_body_list.size();
  DLIList<BodySM*> temp_sm_list(webcut_body_list.size());
  DLIList<BodySM*> result_sm_list;
  DLIList<Body*> body_list(webcut_body_list);
  DLIList<BodySM*> engine_body_sms(count);
  DLIList<Body*> engine_bodies(count);
  GeometryModifyEngine* gme = 0;

  do_attribute_setup();

  while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
  {
    push_vg_attributes_before_modify(engine_body_sms);

    CubitStatus status = gme->webcut(engine_body_sms, vector1, vector2,
              vector3, result_sm_list, imprint );

    restore_vg_after_modify(result_sm_list, engine_bodies);

    if ( status != CUBIT_FAILURE )
      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list );

    engine_bodies.clean_out();
    engine_body_sms.clean_out();
    result_sm_list.clean_out();

    if ( status == CUBIT_FAILURE )
    {
      rval = CUBIT_FAILURE;
      break;
    }
  }

  do_attribute_cleanup();

  return rval;
}

CubitStatus GeometryModifyTool::restore_vg_after_modify(DLIList<BodySM*> &new_sms,
                                                        DLIList<Body*> &old_bodies)
{
  DLIList<TopologyBridge*> old_bridges(old_bodies.size());
  DLIList<TopologyBridge*> new_bridges(new_sms.size());
  CAST_LIST(new_sms, new_bridges, TopologyBridge);

  int k;
  for(k = old_bodies.size(); k>0; k--)
  {
    Body *body = old_bodies.get_and_step();
    TopologyBridge *tb = body->bridge_manager()->topology_bridge();
    if(tb)
    {
      old_bridges.append(tb);
    }
  }

  DLIList<TopologyBridge*> all_bridges;
  all_bridges = new_bridges;
  for(k=old_bridges.size(); k--;)
    all_bridges.append_unique(old_bridges.get_and_step());

  // After a real geometry operation some of the virtual topology bridges may
  // have been blown away and some of them may have just had underlying
  // topology modified.  These calls will try to recognize virtual geometry
  // that has been modified so that it can be deactivated and rebuilt.
  // I am calling it on both of these lists because it is not always clear
  // which one will have the virtual we are interested in.
  GeometryQueryTool::instance()->ige_remove_modified(all_bridges);

  // Now that we have removed any virtual that was affected by the operation we
  // can examine the remaining virtual and remove any COMPOSITE_GEOM attributes
  // that are unneeded since the virtual already exists and wasn't modified by 
  // the operation.  Another way of saying this is that we don't need to process the
  // COMPOSITE_GEOM attributes to rebuild the virtual layer if the virtual layer is 
  // already there.
  GeometryQueryTool::instance()->ige_remove_attributes_from_unmodifed_virtual(all_bridges);

  //Restore virtual
  GeometryQueryTool::instance()->ige_import_geom( all_bridges );

  // At this point we don't need any more attributes on the underlying
  // entities so make sure they are cleaned up.
  GeometryQueryTool::instance()->ige_remove_attributes( all_bridges );

  return CUBIT_SUCCESS;
}

GeometryModifyEngine* GeometryModifyTool::group_bodies_by_engine(
                                  DLIList<Body*>& remaining_bodies,
                                  DLIList<Body*>& engine_bodies,
                                  DLIList<BodySM*>& engine_body_sms ) const
{
  int i = remaining_bodies.size();
  remaining_bodies.reset();
  GeometryModifyEngine* engine = 0;

  if (i == 0)
    return 0;

    // Skip over any bodies that don't have a modify engine.
  while (i--)
  {
    Body* body = remaining_bodies.get();
    TopologyBridge* bridge = 0;
    engine = get_engine(body, &bridge);
    if (engine)
    {
      remaining_bodies.change_to(0);
      engine_bodies.append(body);
      engine_body_sms.append(dynamic_cast<BodySM*>(bridge));
      remaining_bodies.step();
      break;
    }
    remaining_bodies.step();
  }

  // catch case where no engine was found
  if (0 == engine)
  {
    PRINT_WARNING("No geometry modify engine found for this operation.");
    return engine;
  }


    // Get remaining bodies with same modify engine.
  while (i--)
  {
    Body* body = remaining_bodies.get();
    TopologyBridge* bridge = 0;
    if (get_engine(body, &bridge) == engine)
    {
      remaining_bodies.change_to(0);
      engine_bodies.append(body);
      engine_body_sms.append(dynamic_cast<BodySM*>(bridge));
    }
    remaining_bodies.step();
  }

  remaining_bodies.remove_all_with_value(0);
  return engine;
}




// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

CubitStatus GeometryModifyTool::unite( DLIList<BodySM*> &body_sm_list,
                                       DLIList<BodySM*> &new_body_sm_list,
                                       bool keep_old )
{
  //this assumes that all bodies have the same modify engine
  GeometryModifyEngine *gme = get_engine( body_sm_list.get() );
  CubitStatus result = gme->unite(body_sm_list, new_body_sm_list, keep_old);
  return result;
}


CubitStatus GeometryModifyTool::unite( DLIList<Body*> &bodies,
                                       DLIList<Body*> &newBodies,
                                       bool keep_old )
{
   if (bodies.size() <= 1)
   {
      PRINT_WARNING("There is only one body in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }
   if (!okay_to_modify( bodies, "UNITE" ))
     return CUBIT_FAILURE;

  int i;
  const int count = bodies.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  bodies.reset();
  for (i = bodies.size(); i--; )
    entity_list.append_unique(bodies.get_and_step());
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if ( !gme )
  {
      PRINT_ERROR("Performing UNITE with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
  }

    // give 1st body in "old_bodies" list all the
    //  names of all bodies being united

  std::list<CubitString> names_list;
  DLIList<CubitString*> entity_names;

  //Loop through second till end body
  for (i=0; i < bodies.size(); i++)
  {

      //See if body has names
    if(bodies.get()->num_names())
    {
      //Put the names in a list
      bodies.get()->entity_names( entity_names );
      entity_names.reset();

      //Loop through names
      for (int ij = 0; ij < entity_names.size(); ij++)
        names_list.push_back( *entity_names.get_and_step() );

      entity_names.clean_out();
      bodies.get()->remove_entity_names();
    }
    bodies.step();
  }

  DLIList<BodySM*> body_sm_list(bridge_list.size());
  DLIList<BodySM*> new_bodies;
  CAST_LIST(bridge_list, body_sm_list, BodySM);
  CubitStatus result = unite(body_sm_list, new_bodies, keep_old);

  DLIList<Body*> result_list;
  if (!finish_sm_op(bodies, new_bodies, result_list))
    result = CUBIT_FAILURE;

  if (result)
  {
    newBodies += result_list;

   int i;
   for( i=result_list.size(); i--; )
   {
     //Add names to 1st body
     std::list<CubitString>::iterator iter, end = names_list.end();
     for (iter = names_list.begin(); iter != end; ++iter)
       result_list.get_and_step()->entity_name( *iter );
   }
  }
  else
  {
     PRINT_ERROR("UNITE failed\n");
  }

  return result;
}

CubitStatus GeometryModifyTool::chop( DLIList<Body*>& bodies,
                                      DLIList<Body*> &intersectBodies,
                                      DLIList<Body*> &outsideBodies,
                                      Body*& leftoversBody,
                                      bool keep_old,
                                      bool nonreg )
{
   leftoversBody = 0;
   if (bodies.size() <= 1)
   {
      PRINT_WARNING("There is only one body in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }
   if (!okay_to_modify( bodies, "CHOP" ))
     return CUBIT_FAILURE;

  int i;
  const int count = bodies.size();
  DLIList<Body*> original_body_list = bodies;
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  bodies.reset();
  for (i = bodies.size(); i--; )
    entity_list.append_unique(bodies.get_and_step());
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if ( !gme )
  {
      PRINT_ERROR("Performing CHOP with volumes containing geometry\n"
                  " from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
  }

   DLIList<BodySM*> body_sm_list;
   CAST_LIST( bridge_list, body_sm_list, BodySM );
   DLIList<BodySM*> intersect_bodies, outside_bodies;
   BodySM *leftovers_body = 0;

   do_attribute_setup();
   push_vg_attributes_before_modify(body_sm_list);

   CubitStatus result = gme->chop( body_sm_list, intersect_bodies,
                          outside_bodies, leftovers_body, keep_old, nonreg );

   if( result == CUBIT_FAILURE )
   {
     PRINT_ERROR("CHOP failed\n");
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }

   restore_vg_after_modify(intersect_bodies, original_body_list);
   restore_vg_after_modify(outside_bodies, original_body_list);

   DLIList<Body*> result_bodies;
   body_sm_list.clean_out();
   body_sm_list += intersect_bodies;

   if (!finish_sm_op(bodies, body_sm_list, result_bodies))
   {
     PRINT_ERROR("CHOP failed\n");
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }
   intersectBodies += result_bodies;
   bodies.clean_out();

   body_sm_list.clean_out();
   body_sm_list += outside_bodies;
   result_bodies.clean_out();
   if (!finish_sm_op(bodies, body_sm_list, result_bodies))
   {
     PRINT_ERROR("CHOP failed\n");
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }
   outsideBodies += result_bodies;

   if( leftovers_body )
   {
     body_sm_list.clean_out();
     body_sm_list.append( leftovers_body );
     result_bodies.clean_out();
     if (!finish_sm_op(bodies, body_sm_list, result_bodies))
     {
       PRINT_ERROR("CHOP failed\n");
       do_attribute_cleanup();
       return CUBIT_FAILURE;
     }
     leftoversBody = result_bodies.get();
   }

  do_attribute_cleanup();
  return CUBIT_SUCCESS;
}
CubitStatus GeometryModifyTool::hollow( DLIList<Body*>& bodies,
                                        DLIList<RefFace*> faces_to_remove,
                                        DLIList<Body*>& new_bodies,
                                        double depth)
{
  if (bodies.size() <= 0 || faces_to_remove.size() <= 0)
  {
     PRINT_WARNING("Needs at least one body and one face. Nothing modified\n");
     return CUBIT_FAILURE;
  }

  if (!okay_to_modify( bodies, "HOLLOW" ))
    return CUBIT_FAILURE;

  // Get the GeometryEngine for each Body of the list to check
  // if they are the same and if they are GeometryModifyEngine

  const int count = bodies.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(bodies, entity_list);
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if (!gme)
  {
     PRINT_ERROR("Performing THICKEN with volumes containing geometry\n"
                 " from different modeling engines is not allowed.\n"
                 "Delete uncommon geometry on these volumes before operation.\n\n");
     return CUBIT_FAILURE;
  }

  DLIList<BodySM*> new_sms(count);
  DLIList<BodySM*> body_sms(count);
  CAST_LIST(bridge_list, body_sms, BodySM);

  DLIList <Surface*> surfs_to_remove;
  for(int i = 0 ; i < faces_to_remove.size(); i++)
  {
    Surface* surf = faces_to_remove.get_and_step()->get_surface_ptr();
    if(surf)
      surfs_to_remove.append(surf); 
  }

  CubitStatus result = gme->hollow( body_sms, surfs_to_remove, new_sms, depth);

  // check for resued entities, they have been moved and we need to notify observers
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_sms.size(); i++)
  {
    BodySM* bodysm = new_sms.get_and_step();
    DLIList<TopologyBridge*> to_check;
    DLIList<TopologyBridge*> tmp;
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    DLIList<Curve*> curves;
    bodysm->curves(curves);
    DLIList<Point*> points;
    bodysm->points(points);
    to_check.append(bodysm);
    to_check.append(bodysm->lump());
    CAST_LIST_TO_PARENT(surfs, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(curves, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(points, tmp);
    to_check += tmp;

    int k;
    for(k=0; k<to_check.size(); k++)
      if(BridgeManager* m = to_check.get_and_step()->bridge_manager())
        if(TopologyEntity* t = m->topology_entity())
          entities_to_update.append(CAST_TO(t, RefEntity));

  }

  if (!finish_sm_op(bodies, new_sms, new_bodies))
    result = CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return result; 
}
CubitStatus GeometryModifyTool::thicken( DLIList<Body*>& bodies,
                                         DLIList<Body*>& new_bodies,
                                         double depth,
                                         bool both )
{
  if (bodies.size() <= 0)
  {
     PRINT_WARNING("There is only one body in the list. Nothing modified\n");
     return CUBIT_FAILURE;
  }

  if (!okay_to_modify( bodies, "THICKEN" ))
    return CUBIT_FAILURE;

  // Get the GeometryEngine for each Body of the list to check
  // if they are the same and if they are GeometryModifyEngine

  const int count = bodies.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(bodies, entity_list);
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if (!gme)
  {
     PRINT_ERROR("Performing THICKEN with volumes containing geometry\n"
                 " from different modeling engines is not allowed.\n"
                 "Delete uncommon geometry on these volumes before operation.\n\n");
     return CUBIT_FAILURE;
  }

  DLIList<BodySM*> new_sms(count);
  DLIList<BodySM*> body_sms(count);
  CAST_LIST(bridge_list, body_sms, BodySM);
  
  CubitStatus result = gme->thicken( body_sms, new_sms, depth, both);

  // check for resued entities, they have been moved and we need to notify observers
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_sms.size(); i++)
  {
    BodySM* bodysm = new_sms.get_and_step();
    DLIList<TopologyBridge*> to_check;
    DLIList<TopologyBridge*> tmp;
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    DLIList<Curve*> curves;
    bodysm->curves(curves);
    DLIList<Point*> points;
    bodysm->points(points);
    to_check.append(bodysm);
    to_check.append(bodysm->lump());
    CAST_LIST_TO_PARENT(surfs, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(curves, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(points, tmp);
    to_check += tmp;

    int k;
    for(k=0; k<to_check.size(); k++)
      if(BridgeManager* m = to_check.get_and_step()->bridge_manager())
        if(TopologyEntity* t = m->topology_entity())
          entities_to_update.append(CAST_TO(t, RefEntity));

  }

  if (!finish_sm_op(bodies, new_sms, new_bodies))
    result = CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return result;
}


CubitStatus GeometryModifyTool::validate_normals(DLIList<Body*>& bodies,
                                                 RefFace *surf_ref,
                                                 bool reverse)
{
   if (bodies.size() <= 0)
   {
      PRINT_WARNING("There are no entities in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }

   DLIList<RefEntity*> temp;
   CAST_LIST_TO_PARENT(bodies, temp);
   if ( !same_modify_engine(temp, CUBIT_TRUE))
   {
      PRINT_ERROR("Performing VALIDATE NORMALS with volumes containing geometry\n"
            "      from different modeling engines is not allowed. Delete uncommon\n"
            "      geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   // Get the GeometryEngine for each Body of the list to check
   // if they are the same and if they are GeometryModifyEngine
   GeometryModifyEngine* gePtr1 = get_engine(bodies.get());
   GeometryModifyEngine* gePtr2;
   Body* Body_ptr = NULL;
   bodies.reset();
   for( int i = bodies.size(); i > 0; i--)
   {
      Body_ptr = bodies.get_and_step();
      gePtr2 = get_engine(Body_ptr);
      if (gePtr1 != gePtr2)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "  Volumes are associated with different GMEs. \n");
         return CUBIT_FAILURE;
      }
      if ( gePtr2 == NULL)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "       Body %d is not associated with a SME.\n",
            Body_ptr->id());
         return CUBIT_FAILURE;
      }
   }


   DLIList<RefFace*> face_list;
   DLIList<RefFace*> ref_face_list;
   DLIList<RefFace*> free_face_list;
   DLIList<RefFace*> bad_face_list;

   bodies.reset();

   Body* BodyPtr = bodies.get_and_step();
   BodyPtr->ref_faces(face_list);
   BodyPtr->ref_faces(free_face_list);


   RefFace *ref_face_ptr;
   RefFace *inter_face_ptr;

   if(surf_ref != NULL)   // getting the starting surface
   {
     ref_face_ptr = surf_ref;
   }
   else
   {
     ref_face_ptr = face_list.get_and_step();
   }

   ref_face_list.append(ref_face_ptr);
   free_face_list.remove(ref_face_ptr);
   free_face_list.reset();

   while(free_face_list.size())
   {
      DLIList<RefEdge*> curve_list;
      ref_face_ptr->ref_edges(curve_list);
      free_face_list.reset();
      for(int jj=free_face_list.size(); jj > 0; jj--)  // getting a new searching surface
      {
         inter_face_ptr = free_face_list.get_and_step();
         DLIList<RefEdge*> inter_curve_list;
         inter_face_ptr->ref_edges(inter_curve_list);
         curve_list.reset();
        // PRINT_INFO("base face %d working on face %d\n", ref_face_ptr->id(),inter_face_ptr->id());

         for (int k= curve_list.size(); k > 0; k--)  // looping through all of the surface curves
         {
            RefEdge *ref_check_curve = curve_list.step_and_get();
            inter_curve_list.reset();

            for (int kk = inter_curve_list.size(); kk > 0; kk--)
            {
               RefEdge *check_curve = inter_curve_list.step_and_get();

               if(ref_check_curve == check_curve)  // finding if a surface is connected
               {
                  DLIList<CoEdge*> coedge_list;
                  check_curve->co_edges(coedge_list);

                  CoEdge* first_coedge = coedge_list.get_and_step();
                  CoEdge* second_coedge = coedge_list.get_and_step();

                  if((first_coedge->get_sense() == second_coedge->get_sense() &&
                     !bad_face_list.is_in_list(ref_face_ptr) ) ||
                     (first_coedge->get_sense() != second_coedge->get_sense() &&
                     bad_face_list.is_in_list(ref_face_ptr) ))    // finding if a surface has a fliped normal
                  {
                     bad_face_list.append(inter_face_ptr);
                     PRINT_INFO("Surface %d is not consistent\n", inter_face_ptr->id());
                  }

                  ref_face_list.append(inter_face_ptr);  // adding to the searched list
                  free_face_list.remove(inter_face_ptr); // removing from the unsearched list
               }

            }

         }

      }

      ref_face_list.remove(ref_face_ptr);     // removeing from the searched list
      ref_face_list.last();
      if(ref_face_list.size() <= 0)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "       all surfaces must be connected\n");
         return CUBIT_FAILURE;
      }
      ref_face_ptr = ref_face_list.get();
  }

  if (reverse && bad_face_list.size())
  {
//     CubitStatus result = gePtr1->flip_normals(bad_face_list);
//     if ( result == CUBIT_FAILURE )
//     {
        return CUBIT_FAILURE;
//    }
  }
  else if(!bad_face_list.size())
  {
    PRINT_INFO("All surfaces are consistent\n");
  }
  return CUBIT_SUCCESS;
}



CubitStatus GeometryModifyTool::subtract ( Body* tool_body, DLIList<Body*> &from_bodies,
                                           DLIList<Body*> &new_bodies,
                                           bool imprint,
                                           bool keep_old )
{
  DLIList<Body*> temp_body_list;
  temp_body_list.append(tool_body);
  return subtract(temp_body_list, from_bodies, new_bodies, imprint, keep_old);
}


CubitStatus GeometryModifyTool::subtract( DLIList<Body*>  &tool_body_list,
                                          DLIList<Body*> &from_bodies,
                                          DLIList<Body*> &new_bodies,
                                          bool imprint,
                                          bool keep_old )
{
   if(tool_body_list.size() == 0 )
       return CUBIT_FAILURE;
   if(from_bodies.size() == 0 )
       return CUBIT_FAILURE;
   DLIList<Body*> tem_bodies(tool_body_list);

   // cannot subtract from self
   int old_size = from_bodies.size();
   from_bodies -= tool_body_list;
   if (from_bodies.size() != old_size)
   {
     PRINT_WARNING("Cannot subtract body from itself.  Ignoring \"from\" body.\n");
     if (!from_bodies.size())
       return CUBIT_FAILURE;
   }

   tem_bodies += from_bodies;
   if (!okay_to_modify( tem_bodies, "SUBTRACT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> tool_sms(tool_body_list.size());
   DLIList<BodySM*> from_sms(from_bodies.size());
   GeometryModifyEngine* gme = common_modify_engine( tool_body_list, tool_sms );
   GeometryModifyEngine* gme2 = common_modify_engine( from_bodies, from_sms );
   if (!gme || gme != gme2)
   {
      PRINT_ERROR("Performing SUBTRACTION with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

     // Do the subtract operation
   DLIList<BodySM*> new_sms;
   CubitStatus result = gme->subtract(tool_sms, from_sms, new_sms, imprint, keep_old );

   if (!finish_sm_op(tem_bodies, new_sms, new_bodies))
      result = CUBIT_FAILURE;

   if ( result == CUBIT_FAILURE )
   {
     PRINT_ERROR("Subtract FAILED\n" );
     return CUBIT_FAILURE;
   }

   return CUBIT_SUCCESS;
}


CubitStatus GeometryModifyTool::intersect ( Body *tool_body_ptr,
                                            DLIList<Body*> &from_bodies,
                                            DLIList<Body*> &new_bodies,
                                            bool keep_old )
{
   if(tool_body_ptr == NULL )
       return CUBIT_FAILURE;
   if(from_bodies.size() == 0 || from_bodies.get() == NULL )
       return CUBIT_FAILURE;

   DLIList<Body*> tem_bodies = from_bodies;
   tem_bodies.append( tool_body_ptr );
   if (!okay_to_modify( tem_bodies, "INTERSECT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> from_sm_list(tem_bodies.size());
   GeometryModifyEngine* engine = common_modify_engine(tem_bodies, from_sm_list);
   if ( NULL == engine )
   {
      PRINT_ERROR("Performing INTERSECTION with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   BodySM* tool_sm = from_sm_list.pop();


   //cannot intersect tool with itself
   from_sm_list.remove_all_with_value( tool_sm );
   if( from_sm_list.size() == 0 )
   {
     PRINT_ERROR("Cannot intersect volume %d from itself\n",
                  tool_body_ptr->ref_volume()->id() );
     return CUBIT_FAILURE;
   }

     // Do the intersect operation
   DLIList<BodySM*> new_sms;
   CubitStatus result =
       engine->intersect(tool_sm, from_sm_list, new_sms, keep_old );
   if(!finish_sm_op(tem_bodies, new_sms, new_bodies))
      result = CUBIT_FAILURE;

   if ( result == CUBIT_FAILURE )
   {
     PRINT_ERROR("Intersect FAILED\n" );

      return CUBIT_FAILURE;
   }

   return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &from_body_list,
                                         DLIList<Body*> &new_body_list,
                                         CubitBoolean keep_old )
{
   //if (get_group_imprint() == CUBIT_FALSE)
   //  return imprint_singly( from_body_list, new_body_list, keep_old );

     // Check the GeometryEngine for each of the Body's; check to
     // make sure they're all the same
   from_body_list.reset();
   if (!okay_to_modify( from_body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

     //Check for repeats in each individual list and for overlap
     //between the two lists.
   from_body_list.uniquify_ordered();

   DLIList<BodySM*> from_sms(from_body_list.size()), new_sms;
   GeometryModifyEngine* gePtr1 = common_modify_engine(from_body_list, from_sms);
   if ( !gePtr1 )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

#ifdef BOYD17
   DLIList<Body*> dummy_body_list;
#endif

   int process_composites = 0;
   if(contains_composites(from_body_list))
     process_composites = 1;

   if(process_composites)
   {
      // Push virtual attributes down to solid model topology before
      // doing the imprint.
      do_attribute_setup();
      push_vg_attributes_before_modify(from_sms);
      // This must be done after pushing the vg atts because it uses them.
      push_imprint_attributes_before_modify(from_sms);
   }

   DLIList<TopologyBridge*> new_tbs, att_tbs;
   CubitStatus result =  gePtr1->imprint(from_sms, new_sms, keep_old, &new_tbs,
     &att_tbs);

   int i, j;
   if(process_composites)
   {
      // Analyze the results and adjust virtual attributes as necessary.
      GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
        new_sms, from_body_list);

      // Clean up attributes.
      remove_imprint_attributes_after_modify(from_sms, new_sms);

      // Restore the virtual geometry.
      restore_vg_after_modify(new_sms, from_body_list);
   }

   if (get_old_names() == CUBIT_FALSE)
   {
     if (!finish_sm_op(from_body_list, new_sms, new_body_list))
       result = CUBIT_FAILURE;

     if(process_composites)
       do_attribute_cleanup();

     return result;
   }

   if(process_composites)
     do_attribute_cleanup();

   // If old_names is true, need to make sure things are deleted in
   // the correct order so that entities get the same @A type extension
   // on their names.

   // Update existing bodies.
   from_body_list.reset();
   for (i = from_body_list.size(); i--; )
   {
     Body* body = from_body_list.get();
     BodySM* body_sm = body->get_body_sm_ptr();
     if (!body_sm)
     {
       GeometryQueryTool::instance()->destroy_dead_entity(body);
       from_body_list.change_to(0);
     }
     else
     {
       remove_dead_entity_names(body);
       GeometryQueryTool::instance()->make_Body(body_sm);
     }
     from_body_list.step();
   }

      // Construct new bodies
   new_sms.reset();
   for (j = new_sms.size(); j--; )
   {
     BodySM* body_sm = new_sms.get_and_step();
     Body* body = GeometryQueryTool::instance()->make_Body(body_sm);
     new_body_list.append(body);
   }
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   return result;
}

CubitStatus GeometryModifyTool::scale( Body *&body,
                                       const CubitVector& factors, bool check_to_transform )
{
  if( check_to_transform )
    if (!GeometryQueryTool::instance()->okay_to_transform( body ))
      return CUBIT_FAILURE;

  BodySM* bodysm = body->get_body_sm_ptr();
  GeometryModifyEngine* engine = get_engine( bodysm );
  CubitStatus result;
  if( !engine )
  {
    GeometryQueryEngine* tmp_engine = bodysm->get_geometry_query_engine();
    result = tmp_engine->scale( bodysm, factors );
  }
  else
    result = engine->scale( bodysm, factors );

  //for non-uniform scaling, topology can change...need to update stuff
  if( factors.x() != factors.y() ||
      factors.y() != factors.z() ||
      factors.z() != factors.z() )
    body = GeometryQueryTool::instance()->make_Body(bodysm);

  if (result)
  {
    CubitTransformMatrix xform;
    xform.scale_about_origin( factors );
    GeometryQueryTool::instance()->notify_intermediate_of_transform( body, xform );
    GeometryQueryTool::instance()->notify_observers_of_transform( body );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );

  return result;
}


/**********************************************************************************
CubitStatus GeometryModifyTool::imprint_singly( DLIList<Body*> &from_body_list,
                                                DLIList<Body*> &new_body_list,
                                                CubitBoolean keep_old )
{

     // Check the GeometryEngine for each of the Body's; check to
     // make sure they're all the same
   from_body_list.reset();
   if (!okay_to_modify( from_body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

     //Check for repeats in each individual list and for overlap
     //between the two lists.
   from_body_list.uniquify_ordered();

#ifdef BOYD17
   DLIList<BodySM*> from_sms(from_body_list.size()), new_sms;
#endif
   DLIList<BodySM*> from_sms(from_body_list.size());
   GeometryModifyEngine* gePtr1 = common_modify_engine(from_body_list, from_sms);
   if ( !gePtr1 )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   CubitBoolean reset_new_ids = CUBIT_FALSE;
   if (get_new_ids() == CUBIT_FALSE) {
     PRINT_WARNING("New ids must be TRUE when group imprint FALSE; "
                   "setting new ids TRUE for this operation only.\n");
     set_new_ids(CUBIT_TRUE);
     reset_new_ids = CUBIT_TRUE;
   }

     // do the imprinting; bounding box should be checked in
     // SMEEngine function
   new_body_list = from_body_list;
   GeometryQueryTool* gqt = GeometryQueryTool::instance();

   int i;
   for (i = 0; i < new_body_list.size(); i++) {
     for (int j = 1; j < new_body_list.size()-i; j++) {
       new_body_list.reset();
       new_body_list.step(i);
       Body *body_1 = new_body_list.get();
       Body *body_2 = new_body_list.next(j);
       Body *newBodyPtr1 = NULL;
       Body *newBodyPtr2 = NULL;

       if (body_1 == body_2) {
         PRINT_WARNING("Can't imprint a volume with itself.\n");

         if (reset_new_ids == CUBIT_TRUE) set_new_ids(CUBIT_FALSE);

         return CUBIT_FAILURE;
       }

       BodySM *new_sm_1, *new_sm_2;
       CubitStatus status = gePtr1->imprint(body_1->get_body_sm_ptr(),
                                            body_2->get_body_sm_ptr(),
                                            new_sm_1, new_sm_2,
                                            keep_old);

       if ( status != CUBIT_FAILURE &&
            (new_sm_1 != NULL ||
             new_sm_2 != NULL) )
       {
         from_body_list.reset();
         from_body_list.step(i);

           //put the new ones in the new list and
           //remove the olds ones.
         if ( new_sm_1 != NULL )
         {
           if (!body_1->get_body_sm_ptr())
             gqt->destroy_dead_entity(body_1);
           newBodyPtr1 = gqt->make_Body(new_sm_1);
           new_body_list.change_to(newBodyPtr1);
           from_body_list.change_to(0);
         }
         else
           gqt->make_Body(body_1->get_body_sm_ptr());

         new_body_list.step(j);
         from_body_list.step(j);
         if ( new_sm_2 != NULL )
         {
           if (!body_2->get_body_sm_ptr())
             gqt->destroy_dead_entity(body_2);
           newBodyPtr2 = gqt->make_Body(new_sm_2);
           new_body_list.change_to(newBodyPtr2);
           from_body_list.change_to(NULL);
         }
         else
           gqt->make_Body(body_2->get_body_sm_ptr());
       }
       gqt->cleanout_deactivated_geometry();
     }
   }

   from_body_list.remove_all_with_value(NULL);
   Body *temp_body;
   for (i = from_body_list.size(); i > 0; i--)
   {
     temp_body = from_body_list.get_and_step();
     while (new_body_list.move_to(temp_body)) new_body_list.remove();
   }
   PRINT_INFO("\n");
   if( DEBUG_FLAG( 153 ) )
   {
     PRINT_INFO( "  New Body(ies) created:");
     new_body_list.reset();
     for (i = 0; i < new_body_list.size(); i++)
     {
       if (i != 0) PRINT_INFO( ",");
       PRINT_INFO( " %d", new_body_list.get_and_step()->id());
     }
     PRINT_INFO("\n");
     PRINT_INFO( "  Original Body(ies) retained:");
     from_body_list.reset();
     for (i = 0; i < from_body_list.size(); i++)
     {
       if (i != 0) PRINT_INFO( ",");
       PRINT_INFO( " %d", from_body_list.get_and_step()->id());
     }
     PRINT_INFO("\n");
   }

   PRINT_INFO( "  New Volume(s) created:");
   new_body_list.reset();
   DLIList<RefVolume*> new_vol_list;
   for (i = 0; i < new_body_list.size(); i++)
   {
     DLIList<RefVolume*> t2;
     new_body_list.get_and_step()->ref_volumes( t2 );
     new_vol_list += t2;
   }

   for( i = 0; i < new_vol_list.size(); i++ )
   {
     if (i != 0) PRINT_INFO( ",");
     PRINT_INFO( " %d", new_vol_list.get_and_step()->id());
   }
   PRINT_INFO("\n");

   PRINT_INFO( "  Original Volume(s) retained:");
   from_body_list.reset();
   DLIList<RefVolume*> from_vol_list;
   for (i = 0; i < from_body_list.size(); i++)
   {
     DLIList<RefVolume*> t2;
     from_body_list.get_and_step()->ref_volumes( t2 );
     from_vol_list += t2;
   }

   for( i = 0; i < from_vol_list.size(); i++ )
   {
     if (i != 0) PRINT_INFO( ",");
     PRINT_INFO( " %d", from_vol_list.get_and_step()->id());
   }
   PRINT_INFO("\n");


   if (reset_new_ids) set_new_ids(CUBIT_FALSE);

   return CUBIT_SUCCESS;
}
*******************************************************************************/

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &body_list,
                                         DLIList<RefEdge*> &ref_edge_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean show_messages)
{
   // Check the GeometryEngine for each of the bodies; check to
   // make sure they're all the same
   body_list.reset();
   int i;

   if (!okay_to_modify( body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

   const int count = body_list.size() + ref_edge_list.size();
   DLIList<TopologyEntity*> entity_list(count);
   DLIList<TopologyBridge*> bridge_list(count);
   CAST_LIST_TO_PARENT(body_list, entity_list);
   ref_edge_list.reset();
   for (i = ref_edge_list.size(); i--;)
     entity_list.append(ref_edge_list.get_and_step());

   GeometryModifyEngine* gePtr1 = common_modify_engine(entity_list, bridge_list);
   DLIList<BodySM*> body_sm_list(body_list.size());
   DLIList<Curve*> curve_list(ref_edge_list.size());
   CAST_LIST(bridge_list, body_sm_list, BodySM);
   CAST_LIST(bridge_list, curve_list, Curve);

   if ( !gePtr1 ||
        body_sm_list.size() != body_list.size() ||
        curve_list.size() != ref_edge_list.size() )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   DLIList<BodySM*> new_sm_list;
   CubitStatus status = gePtr1->imprint( body_sm_list, curve_list,
                                         new_sm_list, keep_old_body, show_messages );
   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   return status;
}

CubitStatus GeometryModifyTool::imprint( DLIList<RefFace*> &ref_face_list,
                                         DLIList<RefEdge*> &ref_edge_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body )
{
  //get the owning bodies of the faces and edges
    DLIList<Body*> body_list;
    int j;
    for(j=ref_face_list.size(); j--;)
      ref_face_list.get_and_step()->bodies( body_list );
    for(j=ref_edge_list.size(); j--;)
      ref_edge_list.get_and_step()->bodies( body_list );
    body_list.uniquify_unordered();
   if (!okay_to_modify( body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

   DLIList<ModelEntity*> temp_list, temp_list_2, body_me_list;
   CAST_LIST_TO_PARENT(ref_face_list, temp_list);
   CAST_LIST_TO_PARENT(ref_edge_list, temp_list_2);
   temp_list += temp_list_2;
   ModelQueryEngine::instance()->query_model(
    temp_list, DagType::body_type(), body_me_list );

   DLIList<Surface*> surf_list(ref_face_list.size());
   DLIList<Curve*> curve_list(ref_edge_list.size());
   GeometryModifyEngine* gePtr1 = common_modify_engine( ref_face_list,
                                                        ref_edge_list,
                                                        surf_list,
                                                        curve_list );

   if ( !gePtr1 )
     {
        PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                    "different modeling engines is not allowed.\n"
                    "Delete uncommon geometry on these volumes before operation.\n\n");
        return CUBIT_FAILURE;
     }

   body_list.clean_out();
   CAST_LIST(body_me_list, body_list, Body);
   DLIList<BodySM*> new_sm_list;

   CubitStatus status = gePtr1->imprint( surf_list, curve_list,
                                         new_sm_list, keep_old_body );
   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   body_list.reset();
   for (int i = body_list.size(); i--; )
   {
     Body* body = body_list.get_and_step();
     BodySM* bodysm = body->get_body_sm_ptr();
     assert(!!bodysm);
     GeometryQueryTool::instance()->make_Body(bodysm);
   }
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   return status;

}

CubitStatus GeometryModifyTool::imprint( DLIList<Surface*> &surface_list,
                                         DLIList<DLIList<Curve*>*> &curve_lists_list,
                                         Body*& /*new_body*/,
                                         CubitBoolean keep_old_body )
{
  // Get parent bodies
  int i;
  DLIList<Body*> old_body_list;
  surface_list.reset();
  for( i=surface_list.size(); i--; )
  {
    Surface *surf_ptr = surface_list.get_and_step();
    RefEntity* ref_ent = dynamic_cast<RefEntity*>(surf_ptr->topology_entity());
    RefFace *ref_face_ptr = CAST_TO( ref_ent, RefFace );

    if( ref_face_ptr )
    {
      DLIList<Body*> body_list;
      ref_face_ptr->bodies( body_list );
      old_body_list.merge_unique( body_list );
    }
  }

  if( old_body_list.size() > 1 )
  {
    PRINT_ERROR( "This operation requires all surfaces to be from the same volume\n" );
    // Note: this restriction could be pretty easily lifted by sorting the
    //       input lists and calling the GeometryModifyEngine for each body
    //       separately, or having engines handle this.
    return CUBIT_FAILURE;
  }

  // Check engines - must all be the same
  GeometryModifyEngine* gme;
  surface_list.reset();
  gme = get_engine( surface_list.get() );
  for( i=surface_list.size(); i--; )
  {
    Surface *surf_ptr = surface_list.get_and_step();
    GeometryModifyEngine* gme2 = get_engine( surf_ptr );
    if( gme != gme2 )
    {
      PRINT_ERROR( "All surfaces being imprinted must be from the same geometry engine\n" );
      return CUBIT_FAILURE;
    }
  }

  int j;
  for( i=curve_lists_list.size(); i--; )
  {
    DLIList<Curve*> *curve_list_ptr = curve_lists_list.get_and_step();
    for( j=curve_list_ptr->size(); j--; )
    {
      Curve *curve_ptr = curve_list_ptr->get_and_step();
      GeometryModifyEngine* gme2 = get_engine( curve_ptr );
      if( gme != gme2 )
      {
        PRINT_ERROR( "Curves used to imprint must be from same geometry engine as Surface\n" );
        return CUBIT_FAILURE;
      }
    }
  }

  BodySM* new_body_sm = 0;

  CubitStatus status = gme->imprint( surface_list, curve_lists_list,
    new_body_sm, keep_old_body );


  DLIList<Body*> new_body_list;
  DLIList<BodySM*> new_sm_list;
  new_sm_list.append( new_body_sm );

  if (!finish_sm_op(old_body_list, new_sm_list, new_body_list))
    status = CUBIT_FAILURE;

  return status;
}

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &body_list,
                                         DLIList<CubitVector*> &vector_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body )
{
   // Check the GeometryEngine for each of the RefEdges; check to
   // make sure they're all the same
   body_list.reset();

   if (!okay_to_modify( body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> body_sm_list(body_list.size()), new_sm_list;
   GeometryModifyEngine* gePtr1 = common_modify_engine(body_list, body_sm_list);
   if ( !gePtr1 )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   int process_composites = 0;
   if(contains_composites(body_list))
     process_composites = 1;

   if(process_composites)
   {
      // Push virtual attributes down to solid model topology before
      // doing the imprint.
      do_attribute_setup();
      push_vg_attributes_before_modify(body_sm_list);
      // This must be done after pushing the vg atts because it uses them.
      push_imprint_attributes_before_modify(body_sm_list);
   }
   DLIList<TopologyBridge*> new_tbs, att_tbs;

   CubitStatus status = gePtr1->imprint( body_sm_list, vector_list,
                                         new_sm_list, keep_old_body, &new_tbs,
                                            &att_tbs );

   if(process_composites)
   {
      // Analyze the results and adjust virtual attributes as necessary.
      GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
                      new_sm_list, body_list);

      // Clean up attributes.
      remove_imprint_attributes_after_modify(body_sm_list, new_sm_list);

      // Restore the virtual geometry.
      restore_vg_after_modify(new_sm_list, body_list);
   }

   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   if(process_composites)
     do_attribute_cleanup();

   return status;
}
CubitStatus GeometryModifyTool::project_edges( DLIList<RefFace*> &ref_face_list,
                                               DLIList<RefEdge*> &ref_edge_list_in,
                                               DLIList<RefEdge*> &ref_edge_list_new)
{
   // Check the GeometryEngine for each of the RefEdges; check to
   // make sure they're all the same
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Curve*> curve_list_in(ref_edge_list_in.size()), curve_list_new;
  GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                    ref_edge_list_in,
                                                    surface_list,
                                                    curve_list_in );

  if ( !gme ) {
    PRINT_ERROR("In GeometryTool::create_blend\n"
                "       Curves have different modify engines.\n");
    return CUBIT_FAILURE;
  }

   CubitStatus status = gme->
     project_edges( surface_list, curve_list_in, curve_list_new);

   curve_list_new.reset();
   for (int i = curve_list_new.size(); i--; )
   {
     Curve* curve = curve_list_new.get_and_step();
     RefEdge* new_edge = GeometryQueryTool::instance()->make_free_RefEdge(curve);
     PRINT_INFO("Created Curve %d\n", new_edge->id());
     ref_edge_list_new.append(new_edge);
   }

   return status;

}

CubitStatus
GeometryModifyTool::imprint_projected_edges(DLIList<RefFace*> &ref_face_list,
                                            DLIList<RefEdge*> &ref_edge_list_in,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean keep_free_edges)
{

   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Curve*> curve_list(ref_edge_list_in.size());
   GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                     ref_edge_list_in,
                                                     surface_list,
                                                     curve_list );
   if ( !gme )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
   CAST_LIST_TO_PARENT(ref_face_list, query_input);
   ModelQueryEngine::instance()->query_model(query_input, DagType::body_type(), query_output);
   DLIList<Body*> body_list(query_output.size());
   CAST_LIST(query_output, body_list, Body);
   DLIList<BodySM*> new_sm_list;

   CubitStatus status = gme->imprint_projected_edges( surface_list, curve_list,
                                                      new_sm_list, keep_old_body,keep_free_edges);

   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   return status;
}

CubitStatus
GeometryModifyTool::imprint_projected_edges(DLIList<RefFace*> &ref_face_list,
                                            DLIList<Body*> &body_list,
                                            DLIList<RefEdge*> &ref_edge_list_in,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean keep_free_edges)
{

   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Curve*>     curve_list(ref_edge_list_in.size());
   DLIList<BodySM*>  body_sm_list(body_list.size()), new_sm_list;
   GeometryModifyEngine* gme = common_modify_engine(ref_face_list,
                                                    ref_edge_list_in,
                                                    surface_list,
                                                    curve_list);
   if (!gme || gme != common_modify_engine(body_list, body_sm_list))
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

     // Get RefFace bodies
   DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
   CAST_LIST_TO_PARENT(ref_face_list, query_input);
   ModelQueryEngine::instance()->
     query_model( query_input, DagType::body_type(), query_output );
   DLIList<Body*> face_body_list;
   CAST_LIST(query_output, face_body_list, Body);

   CubitStatus status = gme->imprint_projected_edges( surface_list,
                                                      body_sm_list,
                                                      curve_list,
                                                      new_sm_list,
                                                      keep_old_body,
                                                      keep_free_edges);
   face_body_list += body_list;
   if (!finish_sm_op(face_body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   return status;

}


RefEntity* GeometryModifyTool::copy_refentity( RefEntity *old_entity )
{
   RefEntity *new_entity = NULL;

   if( old_entity == NULL )
       return (RefEntity*) NULL;

   if( CAST_TO( old_entity, RefVolume ) )
       return (RefEntity*) NULL;

   RefFace *old_face = CAST_TO( old_entity, RefFace );
   if( old_face )
   {
      RefFace *new_face = instance()->make_RefFace( old_face );
      new_entity = CAST_TO( new_face, RefEntity );
   }

   RefEdge *old_edge = CAST_TO( old_entity, RefEdge );
   if( old_edge )
   {
      RefEdge *new_edge = instance()->make_RefEdge( old_edge );
      new_entity = CAST_TO( new_edge, RefEntity );
   }

   RefVertex *old_vert = CAST_TO( old_entity, RefVertex );
   if( old_vert )
   {
      RefVertex *new_vert = instance()->make_RefVertex( old_vert->coordinates() );
      new_entity = CAST_TO( new_vert, RefEntity );
   }

   return new_entity;
}

//  Calls solid modeling engine to webcut with a sheet body.
CubitStatus GeometryModifyTool::webcut_with_sheet( DLIList<Body*> &webcut_body_list,
                                                   Body *sheet_body,
                                                   DLIList<Body*> &new_bodies,
                                                   CubitBoolean imprint )
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   webcut_body_list.append(sheet_body);
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   if ( !gme )
   {
      PRINT_ERROR("Performing WEBCUTS with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }
   BodySM* tool_sm = body_sm_list.pop();

   do_attribute_setup();

   push_vg_attributes_before_modify(body_sm_list);

   CubitStatus stat = gme->webcut (
       body_sm_list, tool_sm, new_sms, imprint  );

   restore_vg_after_modify(new_sms, original_body_list);

   stat = finish_webcut( webcut_body_list, new_sms, false, stat, new_bodies );
    // leave webcut_body_list as we found it -- remove appended tool body
   webcut_body_list.pop();

   do_attribute_cleanup();

   return stat;
}
//  Calls solid modeling engine to webcut with a sheet body.
CubitStatus GeometryModifyTool::webcut_with_extended_surf( DLIList<Body*> &webcut_body_list,
                                                           RefFace *extend_from,
                                                           DLIList<Body*> &new_bodies,
                                                           int &num_cut,
                                                           CubitBoolean imprint )
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   Surface* surf_ptr = 0;
   if (gme) {
     TopologyBridge* bridge = extend_from->bridge_manager()->topology_bridge(gme->get_gqe());
     surf_ptr = dynamic_cast<Surface*>(bridge);
   }

   if ( !surf_ptr )
   {
      PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   //make the extended face
   Surface * surf = gme->make_Surface(surf_ptr, true);
   if (surf == NULL)
     {
       PRINT_ERROR("webcut tool surface is not created from acis.\n");
       return CUBIT_FAILURE;
     }
 
   //get cutting tool BodySM.
   BodySM* cutting_tool_ptr = surf->bodysm();
   assert(cutting_tool_ptr );

   do_attribute_setup();

   push_vg_attributes_before_modify(body_sm_list);

   //change sjs@cat 1/27/04
   
   CubitStatus stat = gme->webcut ( //gmeList.get()->webcut_with_extended_surf (
       body_sm_list, cutting_tool_ptr, new_sms, imprint  );

   num_cut = new_sms.size();

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   restore_vg_after_modify(new_sms, original_body_list);

   CubitStatus ret = finish_webcut(webcut_body_list, new_sms, false, stat, new_bodies );

   do_attribute_cleanup();

   return ret;
}


CubitStatus GeometryModifyTool::webcut_with_sweep_surfaces_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            bool up_to_next,
                            DLIList<Body*> &new_bodies,
                            CubitBoolean imprint,
                            CubitBoolean merge )
{

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Surface*> surfaces_to_sweep;
   if (gme)
   {
     GeometryModifyEngine* surf_gme = common_modify_engine( tool_faces, surfaces_to_sweep );
     if( gme != surf_gme )
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   BodySM* cutting_tool_ptr = NULL;
   CubitStatus stat = prepare_surface_sweep(body_sm_list,surfaces_to_sweep,
                           sweep_axis,false,false,false,
                           up_to_next,stop_surface, NULL, cutting_tool_ptr, &point, &angle);
   if(stat == CUBIT_FAILURE )
        return stat;

   do_attribute_setup();

   push_vg_attributes_before_modify(body_sm_list);

   stat = gme->webcut(
       body_sm_list, cutting_tool_ptr, new_sms, imprint );

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   restore_vg_after_modify(new_sms, original_body_list);

   stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies );

   do_attribute_cleanup();

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_curves_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            DLIList<Body*> &new_bodies,
                            CubitBoolean imprint,
                            CubitBoolean merge )
{

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Curve*> curves_to_sweep;
   if (gme)
   {
     DLIList<RefFace*> dummy1;
     DLIList<Surface*> dummy2;
     GeometryModifyEngine* surf_gme = common_modify_engine( dummy1, tool_curves,
                                                            dummy2, curves_to_sweep );
     if( gme != surf_gme )
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   //sweep the curves.
   DLIList<GeometryEntity*> ref_ent_list;
   for(int i = 0; i < curves_to_sweep.size(); i++)
     ref_ent_list.append((GeometryEntity*)(curves_to_sweep.get_and_step())); 

   DLIList<BodySM*> swept_bodies;
   CubitStatus stat = gme->sweep_rotational(ref_ent_list,swept_bodies,point,
                                      sweep_axis,angle,0, 0.0, 0,false,false,
                                      false,stop_surface);
   if(stat == CUBIT_FAILURE  && swept_bodies.size() == 0)
     return stat;

   //stitch faces together
   BodySM* cutting_tool_ptr = NULL;
   stat = gme->stitch_surfs(swept_bodies, cutting_tool_ptr);

   if(stat == CUBIT_FAILURE)
     {
       //delete all swept faces
       for(int i = 0; i < swept_bodies.size(); i++)
         gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
       return stat;
     }

   do_attribute_setup();
   push_vg_attributes_before_modify(body_sm_list);

   stat = gme->webcut(
       body_sm_list, cutting_tool_ptr, new_sms, imprint);

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   //delete all swept faces
   for(int i = 0; i < swept_bodies.size(); i++)
      gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;

   restore_vg_after_modify(new_sms, original_body_list);

   stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies );

   do_attribute_cleanup();

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_curves(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            const CubitVector sweep_vector,
                            bool through_all,
                            RefFace *stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            CubitBoolean imprint,
                            CubitBoolean merge )
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Curve*> curves_to_sweep;
   if (gme)
   {
     DLIList<RefFace*> dummy1;
     DLIList<Surface*> dummy2;
     GeometryModifyEngine* surf_gme = common_modify_engine( dummy1, tool_curves,
                                                            dummy2, curves_to_sweep );
     if( gme != surf_gme ){
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
       return CUBIT_FAILURE;
     }
   }
   else{
     PRINT_ERROR("All volumes not from the same modeling engine.\n");
     return CUBIT_FAILURE;
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   CubitVector tmp_sweep_vector = sweep_vector;
   //get model bbox info...will scale sweep vector by its diagonal
   //so that we go far enough
   if( through_all || stop_surf )
     {
       CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
       tmp_sweep_vector.normalize();
       tmp_sweep_vector*=(2*bounding_box.diagonal().length());
     }

   DLIList<GeometryEntity*> ref_ent_list;
   for(int i = 0; i < curves_to_sweep.size(); i++)
     ref_ent_list.append((GeometryEntity*)(curves_to_sweep.get_and_step()));
   DLIList<BodySM*> swept_bodies;
   CubitStatus stat;

   Curve *curve_to_sweep_along = NULL;
   if( edge_to_sweep_along )
   {
     GeometryModifyEngine *tmp_gme = get_engine(edge_to_sweep_along);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Curve %d does is not of same modeling as other input geometry.\n",
                    edge_to_sweep_along->id() );
       return CUBIT_FAILURE;
     }

     curve_to_sweep_along = edge_to_sweep_along->get_curve_ptr();
     DLIList<Curve*> curves_to_sweep_along;
     curves_to_sweep_along.append(curve_to_sweep_along);
     stat = gme->sweep_along_curve(ref_ent_list, swept_bodies,
                              curves_to_sweep_along,0.0,0,false,stop_surface);

     if (!stat && swept_bodies.size() == 0)
       return stat;
   }

   else
   {
     stat = gme->sweep_translational(ref_ent_list, swept_bodies,
                            tmp_sweep_vector,0.0,0, false, false, stop_surface);

     if (!stat && swept_bodies.size() == 0)
       return stat;
   }

   //stitch faces together
   BodySM* cutting_tool_ptr = NULL;
   stat = gme->stitch_surfs(swept_bodies, cutting_tool_ptr);

   if(stat == CUBIT_FAILURE)
   {
       //delete all swept faces
       for(int i = 0; i < swept_bodies.size(); i++)
         gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
       return stat;
   }

   do_attribute_setup();
   push_vg_attributes_before_modify(body_sm_list);

   stat = gme->webcut(
       body_sm_list, cutting_tool_ptr, new_sms, imprint );

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   //delete all swept faces
   for(int i = 0; i < swept_bodies.size(); i++)
      gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;

   restore_vg_after_modify(new_sms, original_body_list);

   stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies );

   do_attribute_cleanup();

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_surfaces(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            const CubitVector sweep_vector,
                            bool sweep_perp,
                            bool through_all,
                            bool outward,
                            bool up_to_next,
                            RefFace* stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            CubitBoolean imprint,
                            CubitBoolean merge )
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Surface*> surfaces_to_sweep;
   if (gme)
   {
     GeometryModifyEngine* surf_gme = common_modify_engine( tool_faces, surfaces_to_sweep );
     if( gme != surf_gme ){
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
       return CUBIT_FAILURE;
     }
   }
   else{
     PRINT_ERROR("All volumes not from the same modeling engine.\n");
     return CUBIT_FAILURE;
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   Curve *curve_to_sweep_along = NULL;
   if( edge_to_sweep_along )
   {
     GeometryModifyEngine *tmp_gme = get_engine(edge_to_sweep_along);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Curve %d does is not of same modeling as other input geometry.\n",
                    edge_to_sweep_along->id() );
       return CUBIT_FAILURE;
     }

     curve_to_sweep_along = edge_to_sweep_along->get_curve_ptr();
   }

   do_attribute_setup();
   push_vg_attributes_before_modify(body_sm_list);

   BodySM* cutting_tool_ptr = NULL;
   CubitStatus stat = prepare_surface_sweep(
      body_sm_list, surfaces_to_sweep, sweep_vector, sweep_perp, through_all, outward,
      up_to_next, stop_surface, curve_to_sweep_along, cutting_tool_ptr );
   if (stat == CUBIT_FAILURE)
     return stat;

   stat = gme->webcut(
       body_sm_list, cutting_tool_ptr, new_sms, imprint );

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   restore_vg_after_modify(new_sms, original_body_list);

   stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies );

   do_attribute_cleanup();

   return stat;
}

//-------------------------------------------------------------------------
// Purpose       : split a multiple volume body into several bodies
//                 each containing a single volume.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::split_body( Body *body_ptr,
                                            DLIList<Body*> &new_bodies ) const
{
   int i;
   DLIList<RefVolume*> ref_vols;
   body_ptr->ref_volumes(ref_vols);
   if ( ref_vols.size() < 2 )
   {
       //no need to split...
     new_bodies.append(body_ptr);
     return CUBIT_SUCCESS;
   }

   DLIList<Body*> b_list;
   b_list.append(body_ptr);
   if (!okay_to_modify( b_list, "SPLIT" ))
     return CUBIT_FAILURE;

     // Call the owning GeometryModifyEngine to split the body
     // so that there is one volume per body.
   BodySM* bodysm_ptr = body_ptr->get_body_sm_ptr();
   GeometryModifyEngine* engine_ptr = get_engine(bodysm_ptr);
   if (!engine_ptr) {
     PRINT_ERROR("There is no modify engine available for this volume."
                 " Volume cannot be split.\n");
     return CUBIT_FAILURE;
   }

   DLIList<BodySM*> new_sm_list;
   CubitStatus stat = engine_ptr->split_body( bodysm_ptr, new_sm_list );
   if ( new_sm_list.size() == 0 )
   {
     PRINT_ERROR("failed in splitting volumes, orginal was lost.\n");
     GeometryQueryTool::instance()->destroy_dead_entity(body_ptr);
     return CUBIT_FAILURE;
   }

   bodysm_ptr = body_ptr->get_body_sm_ptr();
   if (bodysm_ptr)
   {
     remove_dead_entity_names(body_ptr);
     body_ptr = GeometryQueryTool::instance()->make_Body(bodysm_ptr);
   }
   else
   {
     GeometryQueryTool::instance()->destroy_dead_entity(body_ptr);
   }

   new_sm_list.reset();
   for (i = new_sm_list.size(); i--; )
   {
     bodysm_ptr = new_sm_list.get_and_step();
     Body* body = GeometryQueryTool::instance()->make_Body(bodysm_ptr);
     new_bodies.append(body);
   }

   if ( !GeometryModifyTool::instance()->get_new_ids())
   {
       //Now reuse the body ids.
     DLIList<RefVolume*> new_ref_vols;
     for ( int ii = new_bodies.size(); ii > 0; ii-- )
     {
       Body *temp_body = new_bodies.get_and_step();
       new_ref_vols.clean_out();
       temp_body->ref_volumes(new_ref_vols);
       int vol_id = new_ref_vols.get()->id();
       if ( temp_body->id() != vol_id )
       {
           //Check to see if this id is being used by
           //some other body.
         if ( RefEntityFactory::instance()->get_body(vol_id) == NULL )
         {
           temp_body->set_id(vol_id);
         }
           //if it is in use, then we shouldn't mess around with it...
       }

     }
   }

   GeometryQueryTool::instance()->cleanout_deactivated_geometry();
   return stat;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse a body
//
// Special Notes : Moved from Body
//
// Creator       : Jason Kraftcheck
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::reverse( Body* body )
{
   BodySM* body_sm = body->get_body_sm_ptr();
   GeometryModifyEngine* gme = get_engine( body_sm );
   if (!gme) {
      PRINT_ERROR("Body %d does not have a modify engine.\n", body->id());
      return CUBIT_FAILURE;
   }

   CubitStatus stat = gme->reverse_body( body_sm );
   if ( CUBIT_SUCCESS != stat )
   {
     PRINT_ERROR("Reverse failed.\n");
     return CUBIT_FAILURE;
   }
   else
   {
     GeometryQueryTool::instance()->make_Body( body_sm );
     body->notify_all_observers( GEOMETRY_MODIFIED );
     PRINT_INFO("Reversed body %d.\n", body->id());
     return CUBIT_SUCCESS;
   }
}

//-------------------------------------------------------------------------
// Purpose       : split a periodic surface.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::split_periodic(Body *body_ptr,
                                               Body *&new_body_ptr )
{
   DLIList<Body*> b_list(1);
   b_list.append(body_ptr);
   if (!okay_to_modify( b_list, "SPLIT" ))
     return CUBIT_FAILURE;

   BodySM* body_sm = body_ptr->get_body_sm_ptr();
   GeometryModifyEngine* gme = get_engine(body_sm);
   if (!gme) {
      PRINT_ERROR("Volume %d does not have a modify engine.\n", body_ptr->id());
      return CUBIT_FAILURE;
   }

   BodySM* new_sm = 0;

     // Call the default GeometryModifyEngine to create a new Point
   CubitStatus stat = gme->split_periodic( body_sm, new_sm );

   update_body(body_ptr);

   new_body_ptr = 0;
   if (new_sm)
     new_body_ptr = GeometryQueryTool::instance()->make_Body(new_sm);

   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   return stat;
}

//===============================================================================
// Function   : split_surface
// Member Type: PUBLIC
// Description: Split a single surface into one or more pieces
// Author     : Steve Storm (CAT)
// Date       : 07/04
//===============================================================================
CubitStatus
GeometryModifyTool::split_surface( RefFace *ref_face_ptr,
                                   DLIList<CubitVector*> &locations,
                                   DLIList<DLIList<CubitVector*>*> &vec_lists,
                                   CubitBoolean preview_flg,
                                   CubitBoolean create_ref_edges_flg )
{
  // Check for virtual geometry
  DLIList<RefFace*> ref_face_list;
  ref_face_list.append( ref_face_ptr );

  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on surface %d before\n"
      "       operation.\n", ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  int i;
  DLIList<Body*> old_body_list;
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  //bad geom with no body -- dont try to imprint this...
  //quick and dirty fix by (aga@cat|1/7/04)
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "Surface %d is not contained within a parent body\n."
      "      It cannot be split.\n", ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  SplitSurfaceTool sst;
  if( preview_flg == CUBIT_TRUE )
    return sst.preview( ref_face_ptr, locations, vec_lists, create_ref_edges_flg );
  else
   return sst.split_surface( ref_face_ptr, locations, vec_lists );
}

//===============================================================================
// Function   : split_surfaces
// Member Type: PUBLIC
// Description: Split a chain of surfaces into one or more pieces
// Author     : Steve Storm (CAT)
// Date       : 01/04
//===============================================================================
CubitStatus
GeometryModifyTool::split_surfaces( DLIList<RefFace*> &ref_face_list,
                                    int num_segs,
                                    double fraction,
                                    double distance,
                                    RefEdge *from_curve_ptr,
                                    DLIList<RefVertex*> &corner_vertex_list,
                                    DLIList<RefVertex*> &through_vertex_list,
                                    RefEdge *curve_dir_ptr,
                                    CubitBoolean preview_flg,
                                    CubitBoolean create_ref_edges_flg )
{
  // Get parent bodies - all surfs must be from same body
  int i;
  DLIList<Body*> old_body_list;
  RefFace *ref_face_ptr;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  if( old_body_list.size() > 1 )
  {
    PRINT_ERROR( "This operation requires all surfaces to be from the same volume\n" );
    // Note: this restriction could be pretty easily lifted by sorting the
    //       input lists and calling the SplitSurfaceTool separately for each set of
    //       surfaces on each body.
    return CUBIT_FAILURE;
  }

  //bad geom with no body -- dont try to imprint this...
  //quick and dirty fix by (aga@cat|1/7/04)
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "A surface is not contained within a parent body.\n"
      "       It cannot be split.\n");
    return CUBIT_FAILURE;
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Make sure all surfaces are from same geometry engine
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
  {
    PRINT_ERROR("Performing SPLIT with surfaces containing geometry from\n"
      "different modeling engines is not allowed.\n"
      "Delete uncommon geometry on these surfaces before operation.\n\n");
    return CUBIT_FAILURE;
  }

  SplitSurfaceTool split_tool;
  return split_tool.split_surfaces( ref_face_list, num_segs, fraction,
    distance, from_curve_ptr, corner_vertex_list, through_vertex_list,
    curve_dir_ptr, preview_flg, create_ref_edges_flg );
}

//===============================================================================
// Function   : split_surfaces_offset
// Member Type: PUBLIC
// Description: Split a list of surface offset from a curve
// Author     : Sam Showman (CAT)
// Date       : 05/05
//===============================================================================
CubitStatus
GeometryModifyTool::split_surfaces_offset(DLIList<RefFace*> &ref_face_list,
                                          DLIList<RefEdge*> &edge_list,
                                          int num_segs,
                                          double distance,
                                          CubitBoolean partition_flg,
                                          CubitBoolean blunt_flg,
                                          CubitBoolean preview_flg,
                                          CubitBoolean create_ref_edges_flg)
{
  int i;
  DLIList<Body*> old_body_list;
  RefFace *ref_face_ptr;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  //check for geometry with no body
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "A surface is not contained within a parent body.\n"
      "       It cannot be split.\n");
    return CUBIT_FAILURE;
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Make sure all surfaces are from same geometry engine
  // Make sure all curves are from same geometry engine
  CAST_LIST_TO_PARENT(edge_list, ref_ent_list);
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
  {
    PRINT_ERROR("Performing SPLIT with geometry from\n"
      "different modeling engines is not allowed.\n");
    return CUBIT_FAILURE;
  }

  OffsetSplitTool split_tool;
  return split_tool.split_surfaces_offset(
      ref_face_list,
      edge_list,
      num_segs,
      distance,
      partition_flg,
      blunt_flg,
      preview_flg,
      create_ref_edges_flg);
}

//-------------------------------------------------------------------------
// Purpose       : clean RefEntity
//
// Special Notes :
//
// Creator       : Lingyun Pan (CAT)
//
// Creation Date : 07/15/01
//-------------------------------------------------------------------------

CubitStatus
GeometryModifyTool::regularize_refentity(RefEntity *old_entity_ptr,
													               Body *&new_body_ptr)
{
   DLIList<RefEntity*> tmp_ref_ent_list(1);
   tmp_ref_ent_list.append( old_entity_ptr );
   if( GeometryQueryTool::instance()->contains_intermediate_geometry(tmp_ref_ent_list) )
   {
     PRINT_ERROR("%s %d or owning parent(s) is virtual.  Regularizing virtual\n"
       "       geometry is not allowed. Delete virtual geometry first.\n",
               old_entity_ptr->class_name(), old_entity_ptr->id() );
     return CUBIT_FAILURE;
   }

   BasicTopologyEntity* bte_ptr = dynamic_cast<BasicTopologyEntity*>(old_entity_ptr);
   if (!bte_ptr)
   {
     PRINT_ERROR("Invalid entity passed to regularize_refentity(..)\n");
     return CUBIT_FAILURE;
   }

   DLIList<Body*> body_list;
   bte_ptr->bodies(body_list);

   DLIList<TopologyBridge*> bridge_list;
   bte_ptr->bridge_manager()->get_bridge_list(bridge_list);
   bridge_list.reset();

   DLIList<BodySM*> new_sm_list;
   DLIList<Body*> new_body_list;
   CubitStatus stat = CUBIT_SUCCESS;

   for (int i = bridge_list.size(); i--; )
   {
     TopologyBridge* bridge = bridge_list.get_and_step();
     GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge);
     GeometryModifyEngine* gme = get_engine(geom_ptr);
     if (!gme) continue;

     BodySM* new_body_sm = 0;
     if (!gme->regularize_entity( geom_ptr, new_body_sm ))
       stat = CUBIT_FAILURE;

     if (new_body_sm)
       new_sm_list.append(new_body_sm);
   }

   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     stat = CUBIT_FAILURE;

   new_body_ptr = new_body_list.size() ? new_body_list.get() : 0;
   return stat;
}

//-------------------------------------------------------------------------
// Purpose       : split a multiple volume body into several bodies
//                 each containing a single volume.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::regularize_body(Body *body_ptr,
                                                Body *&new_body )
{
   DLIList<Body*> b_list;
   b_list.append(body_ptr);
   if (!okay_to_modify( b_list, "REGULARIZE" ))
     return CUBIT_FAILURE;

   BodySM* body_sm = body_ptr->get_body_sm_ptr(), *new_sm = 0;
   GeometryModifyEngine* gme = get_engine(body_sm);
   if (!gme) {
     PRINT_ERROR("Volume does not have a modify engine.\n" );
     return CUBIT_FAILURE;
  }

   CubitStatus stat = gme->regularize_body( body_sm, new_sm );
   if ( new_sm == NULL )
   {
      PRINT_ERROR("REGULARIZATION failure.\n");
      return CUBIT_FAILURE;
   }

   body_sm = body_ptr->get_body_sm_ptr();
   update_body(body_ptr);

   new_body = GeometryQueryTool::instance()->make_Body(new_sm);
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   return stat;
}

CubitStatus
GeometryModifyTool::create_solid_bodies_from_surfs( DLIList<RefFace*> &ref_face_list,
                                            DLIList<Body*> &new_bodies,
                                            CubitBoolean keep_old,
                                            CubitBoolean heal ) const
{
    //First check to make sure the data is all here.
  for ( int ii = ref_face_list.size(); ii > 0; ii-- )
  {
    RefFace *ref_face = ref_face_list.get_and_step();
    DLIList<Body*> bodies;
    ref_face->bodies(bodies);
    if ( bodies.size() > 1 )
    {
      PRINT_ERROR("Can't create a volume with %s (surface %d) is attached to more\n"
                  "than one volume, or if the attached body is not a sheet body.\n",
                  ref_face->entity_name().c_str(),
                  ref_face->id());
      return CUBIT_FAILURE;
    }
    else if ( bodies.size() == 1 )
    {
      if (!bodies.get()->is_sheet_body())
      {
        PRINT_ERROR("Can't create a volume with %s (surface %d), it is\n"
                    "attached to a body that is not a sheet body.\n",
                    ref_face->entity_name().c_str(),
                    ref_face->id());
        return CUBIT_FAILURE;
      }
    }
  }

  DLIList<TopologyEntity*> entity_list(ref_face_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_face_list.size());
  DLIList<Surface*>       surface_list(ref_face_list.size());
  GeometryModifyEngine* gme;

  CAST_LIST_TO_PARENT(ref_face_list, entity_list);
  gme = common_modify_engine(entity_list, bridge_list);
  CAST_LIST(bridge_list, surface_list, Surface);

  if (!gme) {
    PRINT_ERROR("Surfaces do not share a common modify engine.\n");
    return CUBIT_FAILURE;
  }

  DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
  DLIList<Body*> body_list;
  

  CAST_LIST_TO_PARENT(ref_face_list, query_input);
  ModelQueryEngine::instance()->
    query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST( query_output, body_list, Body );

  int i;
  DLIList<RefFace*> free_face_list;
  for ( i=ref_face_list.size(); i--; )
  {
    RefFace* ref_face = ref_face_list.get_and_step();
    query_output.clean_out();
    ModelQueryEngine::instance()->
      query_model( *ref_face, DagType::body_type(), query_output );
    if (!query_output.size())
      free_face_list.append(ref_face);
  }

  DLIList<BodySM*> new_bodies_sm;
  CubitStatus stat = gme->create_solid_bodies_from_surfs( surface_list, new_bodies_sm, keep_old, heal );
  DLIList<BodySM*> body_sm_list;
  for ( i=new_bodies_sm.size(); i--; )
    body_sm_list.append( new_bodies_sm.get_and_step() );

  if (!finish_sm_op( body_list, body_sm_list, new_bodies))
    stat = CUBIT_FAILURE;

  DLIList<int> id_list (free_face_list.size());
  while (free_face_list.size())
  {
    RefFace* ref_face = free_face_list.pop();
    if (!ref_face->get_surface_ptr())
    {
      id_list.append(ref_face->id());
      GeometryQueryTool::instance()->
        destroy_dead_entity( ref_face );
    }
  }
  GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  if (id_list.size())
    CubitUtil::list_entity_ids( "Destroyed surface(s) ", id_list );

  //new_body = new_body_list.size() ? new_body_list.get() : 0;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : Check for virtual geometry
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
bool GeometryModifyTool::contains_intermediate_geom( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_intermediate_geom(bridges);
}

//-------------------------------------------------------------------------
// Purpose       : Check for virtual geometry
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
bool GeometryModifyTool::contains_intermediate_geom(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}

bool GeometryModifyTool::contains_partitions( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_partitions(bridges);
}

bool GeometryModifyTool::contains_composites( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_composites(bridges);
}

bool GeometryModifyTool::contains_composites(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner) &&
      GeometryQueryTool::instance()->ige_is_composite(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}

bool GeometryModifyTool::contains_partitions(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner) &&
      GeometryQueryTool::instance()->ige_is_partition(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}


CubitBoolean GeometryModifyTool::same_modify_engine(DLIList<TopologyEntity*> &topo_list) const
{
  GeometryModifyEngine *gePtr1 = get_engine(topo_list.get_and_step());
  GeometryModifyEngine *gePtr2;
  for( int i = 1; i < topo_list.size(); i++)
  {
    gePtr2 = get_engine(topo_list.get_and_step());
    if(gePtr1 != gePtr2)
    {
      return CUBIT_FALSE;
    }
  }
  return CUBIT_TRUE;
}

CubitBoolean GeometryModifyTool::same_modify_engine(DLIList<RefEntity*> &ref_entity_list,
						    CubitBoolean check_children) const
{

  DLIList<RefEntity*> complete_entity_list;

  //Check the check_children option and check all the children if necessary
  if (check_children)
  {
    //Make a complete list of all the RefEntitys and their children
    DLIList<RefEntity*> temp = ref_entity_list;
    RefEntity* ref_entity_ptr;
    int i;
    for (i = 0; i < ref_entity_list.size(); i++)
    {
      ref_entity_ptr = ref_entity_list.get_and_step();
      complete_entity_list.clean_out();
      ref_entity_ptr->get_all_child_ref_entities(complete_entity_list);
      temp += complete_entity_list;
    }
    complete_entity_list.clean_out();
    complete_entity_list.merge_unique(temp);
  }

  //Now make sure all the RefEntitys are from the same geometry engine
  DLIList<TopologyEntity*> te_list;
  CAST_LIST(complete_entity_list, te_list, TopologyEntity);
  return same_modify_engine(te_list);

}

CubitStatus
GeometryModifyTool::trim_curve( RefEdge* trim_curve,
                                const CubitVector& trim_vector,
                                const CubitVector& keep_vector )
{
  // Use geometry engine to find intersections
  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme_ptr = get_engine(trim_curve, &bridge);
  Curve *new_curve = 0, *curve = dynamic_cast<Curve*>(bridge);

  if (!gme_ptr || !curve)
    return CUBIT_FAILURE;

  new_curve = gme_ptr->trim_curve( curve, trim_vector, keep_vector );
  if (!new_curve)
    return CUBIT_FAILURE;

  GeometryQueryTool::instance()->destroy_dead_entity( trim_curve );

  RefEdge* new_edge = GeometryQueryTool::instance()->make_free_RefEdge(new_curve);
  PRINT_INFO("Created curve %d\n", new_edge->id());

  return CUBIT_SUCCESS;
}

CubitStatus
GeometryModifyTool::surface_intersection( RefFace *ref_face1,
                                          RefFace *ref_face2,
                                          DLIList<RefEdge*> &ref_edge_list )
{
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  entity_list.append(ref_face1);
  entity_list.append(ref_face2);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Surfaces %d and %d do not share a common solid modeller.\n",
      ref_face1->id(), ref_face2->id() );
    return CUBIT_FAILURE;
  }

  bridge_list.reset();
  Surface* surf0 = dynamic_cast<Surface*>(bridge_list.get_and_step());
  Surface* surf1 = dynamic_cast<Surface*>(bridge_list.get());

  GeometryQueryEngine* gqe = surf0->get_geometry_query_engine();
  // Note the user can set the following value through
  //  set geometry accuracy <val>
  double resabs = gqe->get_sme_resabs_tolerance();

  DLIList<Curve*> curve_list;
  if( gme->surface_intersection( surf0, surf1, curve_list, resabs ) ==
    CUBIT_FAILURE )
    return CUBIT_FAILURE;

  int i;
  Curve *curve_ptr;
  RefEdge *ref_edge_ptr;
  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    curve_ptr = curve_list.get_and_step();
    ref_edge_ptr = GeometryQueryTool::instance()->
      make_free_RefEdge( curve_ptr );
    if( ref_edge_ptr )
    {
      ref_edge_list.append( ref_edge_ptr );
    }
    else
      delete curve_ptr;
  }

  return CUBIT_SUCCESS;
}

RefEdge*
GeometryModifyTool::create_arc_three( RefVertex* ref_vertex1,
                                      RefVertex* ref_vertex2,
                                      RefVertex* ref_vertex3,
                                      CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_vertex1);
  entity_list.append(ref_vertex2);
  entity_list.append(ref_vertex3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Vertices do not share a common solid modeller.\n");
    return 0;
  }

  //if we can reuse vertices, we decide here
  if( full )
  {
    bool need_new_start_point = ref_vertex1->get_parents() > 0;
    if (need_new_start_point)
    {
      bridge_list.reset();
      Point *start_point = gme->make_Point( ref_vertex1->coordinates() );
      bridge_list.change_to( start_point );
    }
  }
  else
  {
    bool need_new_start_point = ref_vertex1->get_parents() > 0;
    bool need_new_end_point = ref_vertex3->get_parents() > 0;

    if (need_new_start_point)
    {
      bridge_list.reset();
      Point *start_point = gme->make_Point( ref_vertex1->coordinates() );
      bridge_list.change_to( start_point );
    }
    if (need_new_end_point)
    {
      bridge_list.reset();
      bridge_list.last();
      Point *end_point = gme->make_Point( ref_vertex3->coordinates());
      bridge_list.change_to( end_point );
    }
  }

  bridge_list.reset();
  Point* point0 = dynamic_cast<Point*>(bridge_list.next(0));
  Point* point1 = dynamic_cast<Point*>(bridge_list.next(1));
  Point* point2 = dynamic_cast<Point*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_three( point0, point1, point2, full );
  if (!curve)
    return 0;

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);
  PRINT_INFO("Created curve %d\n", result->id());
  return result;
}

RefEdge*
GeometryModifyTool::create_arc_three( RefEdge* ref_edge1,
                                      RefEdge* ref_edge2,
                                      RefEdge* ref_edge3,
                                      CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_edge1);
  entity_list.append(ref_edge2);
  entity_list.append(ref_edge3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("The input curves must be associated with a common solid modeller.\n");
    return NULL;
  }

  bridge_list.reset();
  Curve* curve0 = dynamic_cast<Curve*>(bridge_list.next(0));
  Curve* curve1 = dynamic_cast<Curve*>(bridge_list.next(1));
  Curve* curve2 = dynamic_cast<Curve*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_three( curve0, curve1, curve2, full );
  if (!curve)
    return 0;

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);
  PRINT_INFO("Created curve %d\n", result->id());
  return result;
}

RefEdge*
GeometryModifyTool::create_arc_center_edge( RefVertex* ref_vertex1,
                                            RefVertex* ref_vertex2,
                                            RefVertex* ref_vertex3,
                                            const CubitVector &normal,
                                            double radius,
                                            CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_vertex1);
  entity_list.append(ref_vertex2);
  entity_list.append(ref_vertex3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Vertices do not share a common modify engine.\n");
    return 0;
  }

  //if we can reuse vertices, we decide here
  if( full )
  {
    bool need_new_start_point = ref_vertex3->get_parents() > 0;
    if (need_new_start_point)
    {
      bridge_list.reset();
      Point *start_point = gme->make_Point( ref_vertex3->coordinates() );
      bridge_list.change_to( start_point );
    }
  }
  else
  {
    bool need_new_start_point = ref_vertex2->get_parents() > 0;
    bool need_new_end_point = ref_vertex3->get_parents() > 0;

    if (need_new_start_point)
    {
      bridge_list.reset();
      bridge_list.step();
      Point *start_point = gme->make_Point( ref_vertex2->coordinates() );
      bridge_list.change_to( start_point );
    }
    if (need_new_end_point)
    {
      bridge_list.last();
      Point *end_point = gme->make_Point( ref_vertex3->coordinates());
      bridge_list.change_to( end_point );
    }
  }


  bridge_list.reset();
  Point* point0 = dynamic_cast<Point*>(bridge_list.next(0));
  Point* point1 = dynamic_cast<Point*>(bridge_list.next(1));
  Point* point2 = dynamic_cast<Point*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_center_edge( point0, point1, point2,
                                              normal, radius, full );
  if (!curve)
    return 0;

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);
  PRINT_INFO("Created curve %d\n", result->id());
  return result;
}

CubitStatus
GeometryModifyTool::create_curve_combine( DLIList<RefEdge*>& ref_edge_list,
                                    RefEdge *&new_ref_edge_ptr )
{
  DLIList<TopologyBridge*> bridge_list(ref_edge_list.size());
  CubitStatus result = CUBIT_FAILURE;

  int count = ref_edge_list.size();
  if (count == 0)
  {
    PRINT_ERROR("No edges to combine.\n");
    return result;
  }

  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( ref_edge_list, entity_list );
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Edges do not share a common modify engine.\n");
    return result;
  }

  Curve* new_curve_ptr;
  DLIList<Curve*> curve_list(count);
  CAST_LIST( bridge_list, curve_list, Curve );
  result = gme->create_curve_combine(curve_list, new_curve_ptr);

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Find a common geometry engine
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/00
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<TopologyEntity*>& topology_list,
                                          DLIList<TopologyBridge*>& engine_bridges,
                                          CubitBoolean /*allow_virtual_engine*/ ) const
{
  topology_list.reset();

  TopologyEntity* topo_ptr = topology_list.get_and_step();
  if (!topo_ptr)
     return (GeometryModifyEngine*)NULL;
  DLIList<TopologyBridge*> first_bridge_list;
  topo_ptr->bridge_manager()->get_bridge_list( first_bridge_list );

  first_bridge_list.reset();
  GeometryModifyEngine* gme_ptr = 0;
  for( int i = first_bridge_list.size(); i > 0; i-- )
  {
    TopologyBridge* bridge_ptr = first_bridge_list.get_and_step();
    engine_bridges.clean_out();
    engine_bridges.append( bridge_ptr );
    gme_ptr = get_engine(bridge_ptr);

    if( !gme_ptr )
     return (GeometryModifyEngine*)NULL;

    topology_list.reset();
    topology_list.step(); //skip first entry
    for( int j = topology_list.size(); j > 1; j-- )
    {
      topo_ptr = topology_list.get_and_step();
      bridge_ptr = topo_ptr->bridge_manager()->topology_bridge(gme_ptr->get_gqe());
      if( bridge_ptr ) engine_bridges.append( bridge_ptr );
      else break;
    }

    if( engine_bridges.size() == topology_list.size() )
      break;

    gme_ptr = 0;
  }

  if( !gme_ptr )
  {
    engine_bridges.clean_out();
    PRINT_ERROR("Entities do not belong to the same geometry engine.\n");
  }

  return gme_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Find common modify engine for Bodies
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/26/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<Body*>& input,
                                          DLIList<BodySM*>& output ) const
{
  input.reset();
  Body* body_ptr = input.get();
  BodySM* body_sm_ptr = body_ptr->get_body_sm_ptr();
  GeometryModifyEngine* engine = get_engine(body_sm_ptr);

  for (int i = input.size(); i--; )
  {
    body_ptr = input.get_and_step();
    body_sm_ptr = body_ptr->get_body_sm_ptr();
    output.append(body_sm_ptr);

    if (!body_sm_ptr)
    {
      PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", body_ptr->id());
      output.clean_out();
      return 0;
    }

    if (get_engine(body_sm_ptr) != engine)
    {
      output.clean_out();
      return 0;
    }
  }

  return engine;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefFaces and RefEdges.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/26/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefFace*>& face_list,
                                          DLIList<RefEdge*>& edge_list,
                                          DLIList<Surface*>& surf_list,
                                          DLIList<Curve*>& curv_list ) const
{
  int i;
  const int count = face_list.size() + edge_list.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);

  face_list.reset();
  edge_list.reset();
  for (i = face_list.size(); i--; )
    entity_list.append(face_list.get_and_step());
  for (i = edge_list.size(); i--; )
    entity_list.append(edge_list.get_and_step());

  GeometryModifyEngine* engine = common_modify_engine(entity_list, bridge_list);
  if (!engine)
    return 0;

  entity_list.reset();
  CAST_LIST(bridge_list, surf_list, Surface);
  CAST_LIST(bridge_list, curv_list, Curve  );
  if (surf_list.size() != face_list.size() || curv_list.size() != edge_list.size())
    return 0;

  return engine;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefFaces.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/05/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefFace*>& face_list,
                                          DLIList<Surface*>& surf_list ) const
{
  const int size = face_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( face_list, topo_list );
  result = common_modify_engine( topo_list, geom_list );

  CAST_LIST( geom_list, surf_list, Surface );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefEdges.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Steve Storm
//
// Creation Date : 03/09/05
//-------------------------------------------------------------------------
GeometryModifyEngine* GeometryModifyTool::common_modify_engine(
                                      DLIList<RefEdge*>& edge_list,
                                      DLIList<Curve*>& curve_list ) const
{
  const int size = edge_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( edge_list, topo_list );
  result = common_modify_engine( topo_list, geom_list );

  CAST_LIST( geom_list, curve_list, Curve );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefVertices.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Steve Storm
//
// Creation Date : 03/25/05
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefVertex*>& vertex_list,
                                          DLIList<Point*>& point_list ) const
{
  const int size = vertex_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( vertex_list, topo_list );
  result = common_modify_engine( topo_list, geom_list );

  CAST_LIST( geom_list, point_list, Point );
  return result;
}

// Separates the list of bodies so that there is one body per volume
// after a webcut.  Checks the sepAfterWebcut flag.
CubitStatus GeometryModifyTool::separate_body_after_webcut( DLIList<Body*> &input_list,
                                              DLIList<Body*> &output_list) const
{
    //First see if we should spearate.
  if ( !sepAfterWebcut )
  {
    output_list = input_list;
    return CUBIT_SUCCESS;
  }
  DLIList<Body*> temp_body_list;
  for ( int ii = input_list.size(); ii > 0; ii-- )
  {
    Body *body_ptr = input_list.get_and_step();
    temp_body_list.clean_out();
    split_body(body_ptr, temp_body_list);
    output_list += temp_body_list;
  }
  return CUBIT_SUCCESS;
}

GeometryModifyEngine *GeometryModifyTool::get_engine(TopologyBridge *tb_ptr) const
{
  int i;
  GeometryModifyEngine *gme;
  for (i = 0; i < gmeList.size(); i++) {
    gme = gmeList.next(i);
    if (gme->is_modify_engine(tb_ptr)) return gme;
  }
  return NULL;
}

GeometryModifyEngine *GeometryModifyTool::get_engine(TopologyEntity *te_ptr,
                                                     TopologyBridge **bridge) const
{
  int i;
  GeometryModifyEngine *gme = 0;
  TopologyBridge* tb_ptr = NULL;
  BridgeManager* bm = te_ptr->bridge_manager();
  DLIList<TopologyBridge*> bridges(bm->number_of_bridges());
  bm->get_bridge_list(bridges);
  bridges.reset();
  for (i = bridges.size(); !gme && i--; )
  {
    tb_ptr = bridges.get_and_step();
    gme = get_engine(tb_ptr);
  }

  if (bridge && gme)
    *bridge = tb_ptr;

  return gme;
}

CubitStatus GeometryModifyTool::get_offset_intersections( RefEdge* ref_edge1, RefEdge* ref_edge2,
                                                          DLIList<CubitVector*>& intersection_list,
                                                          double offset, CubitBoolean ext_first )
{
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  entity_list.append(ref_edge1);
  entity_list.append(ref_edge2);
  GeometryModifyEngine* gme_ptr1 = common_modify_engine(entity_list,bridge_list);

  if( gme_ptr1 == 0 )
  {
      PRINT_ERROR( "Curves %d and %d do not have the same underlying geometry modeling engine\n"
                   "      For intersection calculations, they must be the same\n",
         ref_edge1->id(), ref_edge2->id() );
      return CUBIT_FAILURE;
   }

   bridge_list.reset();
   Curve* curve0 = dynamic_cast<Curve*>(bridge_list.next(0));
   Curve* curve1 = dynamic_cast<Curve*>(bridge_list.next(1));

   return gme_ptr1->get_offset_intersections( curve0, curve1, intersection_list,
                                              offset, ext_first );
}

CubitStatus GeometryModifyTool::get_offset_intersections( RefEdge* ref_edge_ptr,
                                                          RefFace* ref_face_ptr,
                                                          DLIList<CubitVector*>& intersection_list,
                                                          double offset,
                                                          CubitBoolean ext_surf )
{
  // If curve is straight and surface is planar, compute their intersection;
  // otherwise use the geometry engine to do it.

  Curve* curve_ptr = ref_edge_ptr->get_curve_ptr();
  Surface* surface_ptr = ref_face_ptr->get_surface_ptr();

  if( curve_ptr == NULL )
  {
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Curve %d\n"
      "       This is a bug - please report it\n", ref_edge_ptr->id() );
    return CUBIT_FAILURE;
  }

  if( surface_ptr == NULL )
  {
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Surface %d\n"
      "       This is a bug - please report it\n", ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  // If straight line and plane, find location right here analytically.
  CubitVector pln_origin, pln_normal;
  CubitVector crv_origin, crv_direction;
  if( ref_face_ptr->get_point_normal( pln_origin, pln_normal ) &&
     ref_edge_ptr->get_point_direction( crv_origin, crv_direction ) )
  {
     double pln_orig[3], pln_norm[3];
     pln_orig[0]=pln_origin.x(); pln_orig[1]=pln_origin.y(); pln_orig[2]=pln_origin.z();
     pln_norm[0]=pln_normal.x(); pln_norm[1]=pln_normal.y(); pln_norm[2]=pln_normal.z();

     double crv_orig[3], crv_dir[3];
     crv_orig[0]=crv_origin.x(); crv_orig[1]=crv_origin.y(); crv_orig[2]=crv_origin.z();
     crv_dir[0]=crv_direction.x(); crv_dir[1]=crv_direction.y(); crv_dir[2]=crv_direction.z();

     AnalyticGeometryTool *agt = AnalyticGeometryTool::instance();

     if( agt->is_vec_perp( pln_norm, crv_dir ) )
     {
        PRINT_ERROR( "Line is parallel to the plane - intersection not possible\n" );
        return CUBIT_FAILURE;
     }

     double ang = agt->angle_vec_vec( pln_norm, crv_dir );
     if( ang > AGT_PI_DIV_2 )
        ang = AGT_PI - ang;

     double int_pnt[3];
     agt->int_ln_pln( crv_orig, crv_dir, pln_orig, pln_norm, int_pnt );

     // ang2 is angle between line and plane
     double ang2 = AGT_PI_DIV_2 - ang;

     double hypotenuse = offset/sin(ang2);

     double final_pnt[3];

     agt->next_pnt( int_pnt, crv_dir, hypotenuse, final_pnt );

     double end1[3], end2[3];
     CubitVector start, end;
     start = ref_edge_ptr->start_coordinates();
     end = ref_edge_ptr->end_coordinates();
     end1[0]=start.x(); end1[1]=start.y(); end1[2]=start.z();
     end2[0]=end.x(); end2[1]=end.y(); end2[2]=end.z();

     CubitVector curve_position;
     if( agt->is_pnt_on_ln_seg( final_pnt, end1, end2 ) )
     {
        curve_position.x( final_pnt[0] );
        curve_position.y( final_pnt[1] );
        curve_position.z( final_pnt[2] );
     }
     else
     {
        agt->reverse_vec( crv_dir, crv_dir );

        agt->next_pnt( int_pnt, crv_dir, hypotenuse, final_pnt );

        if( agt->is_pnt_on_ln_seg( final_pnt, end1, end2 ) )
        {
           curve_position.x( final_pnt[0] );
           curve_position.y( final_pnt[1] );
           curve_position.z( final_pnt[2] );
        }
        else
        {
           PRINT_ERROR( "Resultant point does not lie on bounded curve\n" );
           return CUBIT_FAILURE;
        }
     }
     CubitVector *vec_ptr = new CubitVector( curve_position );
     intersection_list.append( vec_ptr );
  }
  else
  {
     DLIList<TopologyEntity*> entity_list(2);
     DLIList<TopologyBridge*> bridge_list(2);
     entity_list.append(ref_edge_ptr);
     entity_list.append(ref_face_ptr);
     GeometryModifyEngine* gme_ptr1 = common_modify_engine(entity_list,bridge_list);

     if( gme_ptr1 == 0 )
     {
        PRINT_ERROR( "Curve %d and Surface %d do not have the same underlying geometry modeling engine\n"
           "       For intersection calculations, they must be the same\n",
           ref_edge_ptr->id(), ref_face_ptr->id() );
        return CUBIT_FAILURE;
     }

     bridge_list.reset();
     curve_ptr = dynamic_cast<Curve*>(bridge_list.next(0));
     surface_ptr = dynamic_cast<Surface*>(bridge_list.next(1));

     // Use geometry engine to find intersections
     return gme_ptr1->get_offset_intersections( curve_ptr, surface_ptr,
                                                intersection_list, offset, ext_surf );

  }

  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::get_mid_plane( RefFace *ref_face1,
                                               RefFace *ref_face2,
                                               Body *body_to_trim_to,
                                               DLIList<RefFace*> &mid_plane_surfs ) const
{
  if( ref_face1 == ref_face2 )
  {
    PRINT_ERROR("Cannot create midplane between the same surface.\n"
                "       Surface %d was entered twice\n",  ref_face1->id() ); 
    return CUBIT_FAILURE;
  }

  BodySM* body_sm_to_trim_to = body_to_trim_to->get_body_sm_ptr();
  GeometryModifyEngine *gme1_ptr = get_engine(body_sm_to_trim_to);
  if ( !gme1_ptr )
  {
    PRINT_ERROR("Geometry can't be modified, no associated modify engine.\n");
    return CUBIT_FAILURE;
  }

  CubitVector normal_1, normal_2, point_1, point_2, point_3;
  CubitPlane plane_1, plane_2;
  CubitVector p_mid, n_mid;

  point_1 = ref_face1->center_point();
  point_2 = ref_face2->center_point();

  normal_1 = ref_face1->normal_at(point_1);
  normal_2 = ref_face2->normal_at(point_2);

  plane_1 = CubitPlane(normal_1,point_1);
  plane_2 = CubitPlane(normal_2,point_2);

  if(point_1 == point_2)
  {
    PRINT_ERROR("Since both surfaces share the same point, the midplane is not well-defined\n");
    return CUBIT_FAILURE;
  }
  else
  {
    CubitVector temp1 = point_2;
    temp1 = plane_1.project(temp1);
    temp1 -= point_2;
    if ( temp1.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
    {
      PRINT_ERROR("Since both planes are the same, the mid-plane is not well-defined.\n");
      return CUBIT_FAILURE;
    }
  }

  if(plane_1.normal()==plane_2.normal() || plane_1.normal()==-plane_2.normal())
  {
    p_mid = (point_1+point_2)/2;
    n_mid = plane_1.normal();
  }
  else
  {

    CubitVector direction_of_line;
    plane_1.intersect(plane_2,p_mid,direction_of_line);
    direction_of_line.normalize();

    // Find if point_1 and point_2 are on the line of intersection
    // If they are, then the mid-plane is not well-defined
    CubitVector p1 = point_1-p_mid;
    CubitVector p2 = point_2-p_mid;
    p1.normalize();
    p2.normalize();

    if(p1==direction_of_line || p1==-direction_of_line)
    {
      PRINT_ERROR("P1 is on the line of intersection.\n");
      return CUBIT_FAILURE;
    }

    if(p2==direction_of_line || p2==-direction_of_line)
    {
      PRINT_ERROR("P2 is on the line of intersection.\n");
      return CUBIT_FAILURE;
    }

    CubitVector v1 = p1 - (p1%direction_of_line)*direction_of_line;
    v1.normalize();

    CubitVector v2 = p2 - (p2%direction_of_line)*direction_of_line;
    v2.normalize();

    n_mid = v1 - v2;
    n_mid.normalize();
  }

  CubitPlane mid_plane(n_mid, p_mid);
  point_1 = p_mid;
    //find three points that will define the infinite plane from the
    //mid plane.
  CubitVector test1(1,0,0), test1n(-1,0,0),test2(0,1,0);
  CubitVector direction1;
    //through the point in any direction just not along the
    //normal direction.
  if(n_mid != test1 && n_mid != test1n )
    direction1 = test1 + n_mid;
  else
    direction1 = test2 + n_mid;

  point_2 = p_mid + direction1;
  point_2 = mid_plane.project(point_2);

  direction1 = point_2-point_1;
  CubitVector direction2 = direction1*n_mid;
  point_3 = point_1 + direction2;

  DLIList<BodySM*> midplane_bodysm_list ;
  CubitStatus ret = gme1_ptr->get_mid_plane(point_1, point_2, point_3,
                                    body_sm_to_trim_to, midplane_bodysm_list );

  if (midplane_bodysm_list.size() > 0)
  {
#ifdef BOYD17
    DLIList<Body*> bodies;
    DLIList<Surface*> surfs;
#endif

    Body *midplane_body;

    for (int j = 0; j < midplane_bodysm_list.size(); j++)
    {
      BodySM* midplane_body_sm = midplane_bodysm_list.get_and_step();
      midplane_body = GeometryQueryTool::instance()->make_Body(midplane_body_sm);

      DLIList<RefFace*> ref_faces;
      midplane_body->ref_faces( ref_faces );

      //make each surface of the body into its own body
      int i;
      for( i=0; i<ref_faces.size(); i++ )
      {
        RefEntity *new_entity_ptr;
        new_entity_ptr = GeometryModifyTool::instance()->copy_refentity(ref_faces.get_and_step());
        RefFace *ref_face_ptr = CAST_TO(new_entity_ptr, RefFace);
        mid_plane_surfs.append( ref_face_ptr );
      }
      GeometryQueryTool::instance()->delete_Body( midplane_body );
    }
  }
  else
    return CUBIT_FAILURE;

  return ret;
}

//=============================================================================
// Function   : get_planar_mid_surface
// Member Type: LOCAL
// Description: Calculates a mid-surface between 2 planar surfaces.
// Author     : Philippe Pebay
// Date       : 05/15/06
//=============================================================================
CubitStatus get_planar_mid_surface( RefFace* ref_face1,
				    RefFace* ref_face2,
				    BodySM* body_sm_to_trim_to,
				    DLIList<BodySM*>& midsurface_bodysm_list,
				    GeometryModifyEngine *gme_ptr )
{
    CubitVector normal_1, normal_2, point_1, point_2, point_3;
    CubitPlane plane_1, plane_2;
    CubitVector p_mid, n_mid;

    point_1 = ref_face1->center_point();
    point_2 = ref_face2->center_point();

    normal_1 = ref_face1->normal_at(point_1);
    normal_2 = ref_face2->normal_at(point_2);

    plane_1 = CubitPlane(normal_1,point_1);
    plane_2 = CubitPlane(normal_2,point_2);

    if(point_1 == point_2)
    {
      PRINT_ERROR( "In GeometryModifyTool:: get_planar_mid_surface\n"
		   "       Since both surfaces share the same point, the midsurface is not well-defined\n");
      return CUBIT_FAILURE;
    }
    else
    {
      CubitVector temp1 = point_2;
      temp1 = plane_1.project(temp1);
      temp1 -= point_2;
      if ( temp1.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
	{
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       Since both planes are the same, the midsurface is not well-defined.\n");
	return CUBIT_FAILURE;
	}
    }

    if ( ( normal_1.about_equal( normal_2 ) ) || ( (-normal_1).about_equal( normal_2 ) ) )
    {
      p_mid = (point_1+point_2)/2;
      n_mid = plane_1.normal();
    }
    else
    {
      CubitVector direction_of_line;
      plane_1.intersect(plane_2,p_mid,direction_of_line);
      direction_of_line.normalize();

      // Find if point_1 and point_2 are on the line of intersection
      // If they are, then the mid-plane is not well-defined
      CubitVector p1 = point_1-p_mid;
      CubitVector p2 = point_2-p_mid;
      p1.normalize();
      p2.normalize();

      if(p1==direction_of_line || p1==-direction_of_line)
      {
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       P1 is on the line of intersection.\n");
	return CUBIT_FAILURE;
      }

      if(p2==direction_of_line || p2==-direction_of_line)
      {
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       P2 is on the line of intersection.\n");
	return CUBIT_FAILURE;
      }

      CubitVector v1 = p1 - (p1%direction_of_line)*direction_of_line;
      v1.normalize();

      CubitVector v2 = p2 - (p2%direction_of_line)*direction_of_line;
      v2.normalize();

      n_mid = v1 - v2;
      n_mid.normalize();
    }

    CubitPlane mid_plane(n_mid, p_mid);
    point_1 = p_mid;

    //find three points that will define the infinite plane from the
    //mid plane.through the point in any direction just not along the
    //normal direction
    CubitVector Xdir(1,0,0), Ydir(0,1,0);
    CubitVector direction1;

    if ( ( ! n_mid.about_equal( Xdir ) ) && ( ! (-n_mid).about_equal( Xdir ) ) )
      direction1 = Xdir + n_mid;
    else
      direction1 = Ydir + n_mid;

    point_2 = p_mid + direction1;
    point_2 = mid_plane.project(point_2);

    direction1 = point_2-point_1;
    CubitVector direction2 = direction1*n_mid;
    point_3 = point_1 + direction2;

    CubitStatus ret = gme_ptr->get_mid_plane(point_1, point_2, point_3,
					      body_sm_to_trim_to, midsurface_bodysm_list );
    return ret;
}

//=============================================================================
// Function   : get_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 surfaces (calls appropriate
//              specialized methods if needed).
// Author     : Philippe Pebay
// Date       : 03/07/06
//=============================================================================
CubitStatus GeometryModifyTool::get_mid_surface( RefFace *ref_face1,
                                               RefFace *ref_face2,
                                               Body *body_to_trim_to,
                                               DLIList<RefFace*> &mid_surface_surfs ) const
{
  if( ref_face1 == ref_face2 )
  {
    PRINT_ERROR("Cannot create midplane between the same surface.\n"
                "       Surface %d was entered twice\n",  ref_face1->id() ); 
    return CUBIT_FAILURE;
  }

  BodySM *body_sm_to_trim_to = body_to_trim_to->get_body_sm_ptr();
  GeometryModifyEngine *gme1_ptr = get_engine(body_sm_to_trim_to);
  if ( !gme1_ptr )
  {
    PRINT_ERROR("In GeometryModifyTool::get_mid_surface\n"
		"       Geometry can't be modified, no associated modify engine.\n");
    return CUBIT_FAILURE;
  }

  bool found_case = false;
  CubitStatus ret;
  DLIList<BodySM*> midsurfaces;
  // Plane to plane case
  if ( ( ref_face1->geometry_type() == PLANE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == PLANE_SURFACE_TYPE ) )
  {
    found_case = true;
    ret = get_planar_mid_surface( ref_face1, ref_face2, body_sm_to_trim_to, midsurfaces, gme1_ptr );
  }

  // Quadric to quadric cases
  if ( ( ( ref_face1->geometry_type() == SPHERE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == SPHERE_SURFACE_TYPE ) )
       ||
       ( ( ref_face1->geometry_type() == CONE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == CONE_SURFACE_TYPE ) )
       ||
       ( ( ref_face1->geometry_type() == TORUS_SURFACE_TYPE ) && ( ref_face2->geometry_type() == TORUS_SURFACE_TYPE ) ) )
  {
    found_case = true;
    DLIList<TopologyEntity*> entity_list(2);
    DLIList<TopologyBridge*> bridge_list(2);

    entity_list.append(ref_face1);
    entity_list.append(ref_face2);
    GeometryModifyEngine* gme2_ptr = common_modify_engine(entity_list,bridge_list);

    if( gme2_ptr == 0 )
    {
      PRINT_ERROR( "In GeometryModifyTool::get_mid_surface\n"
		   "       Surfaces %d and %d do not have the same underlying geometry modeling engine.\n"
		   "       For midsurface calculations, they must be the same\n",
		   ref_face1->id(), ref_face2->id() );
      return CUBIT_FAILURE;
    }

    if( gme1_ptr != gme2_ptr )
    {
      PRINT_ERROR( "In GeometryModifyTool::get_mid_surface\n"
		   "       Body and surfaces do not have the same underlying geometry modeling engine.\n"
		   "       For midsurface calculations, they must be the same\n");
      return CUBIT_FAILURE;
    }

    bridge_list.reset();
    Surface* surface1_ptr = dynamic_cast<Surface*>(bridge_list.next(0));
    Surface* surface2_ptr = dynamic_cast<Surface*>(bridge_list.next(1));

    // Sphere to sphere case
    if ( ( ref_face1->geometry_type() == SPHERE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == SPHERE_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_spheric_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurfaces );
    }

    // Cone to cone case
    if ( ( ref_face1->geometry_type() == CONE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == CONE_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_conic_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurfaces );
    }

    // Torus to torus case
    if ( ( ref_face1->geometry_type() == TORUS_SURFACE_TYPE ) && ( ref_face2->geometry_type() == TORUS_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_toric_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurfaces );
    }
  }

  // Unsupported pair of surfaces
  if ( ! found_case )
  {
    PRINT_ERROR("In GeometryModifyTool::get_mid_surface\n"
		"       Midsurface calculation not yet supported for such a pair of surfaces.\n");
    return CUBIT_FAILURE;
  }

  if ( midsurfaces.size() > 0)
  {
    Body *midsurface_body;
    for(int j = 0; j < midsurfaces.size()  ;j++)
    {
      BodySM* midsurface_body_sm = midsurfaces.get_and_step();
      midsurface_body = GeometryQueryTool::instance()->make_Body(midsurface_body_sm);

      DLIList<RefFace*> ref_faces;
      midsurface_body->ref_faces( ref_faces );

      //make each surface of the body into its own body
      int i;
      for( i=0; i<ref_faces.size(); i++ )
      {
        RefEntity *new_entity_ptr;
        new_entity_ptr = GeometryModifyTool::instance()->copy_refentity(ref_faces.get_and_step());
        RefFace *ref_face_ptr = CAST_TO(new_entity_ptr, RefFace);
        mid_surface_surfs.append( ref_face_ptr );
      }
      GeometryQueryTool::instance()->delete_Body( midsurface_body );
    }
    return ret;
  }
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::set_default_gme(GeometryModifyEngine* GMEPtr)
{
    if(!GMEPtr)
        return CUBIT_FAILURE;

    int i;
    for (i = 0; i < gmeList.size(); i++)
    {
        if(GMEPtr == gmeList.get())
            break;

        gmeList.step();
    }

    if(i == gmeList.size())
        return CUBIT_FAILURE;

    GeometryModifyEngine* temp_ptr = gmeList.get();
    gmeList.remove();
    gmeList.insert_first(temp_ptr);
    return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Consolidate common setup code for tweak surface operations
// Author     : Jason Kraftcheck
// Date       : 11/05/03
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefFace*> &input_faces,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Surface*> &output_surfaces )
{
  if( input_faces.size() == 0 )
  {
    PRINT_ERROR( "No surfaces specified\n" );
    return 0;
  }

  //are any of the input faces merged? ...disallow that
  int i;
  for( i=input_faces.size(); i--; )
  {
    RefFace *tmp_face = input_faces.get();
    if( tmp_face->is_merged() )
    {
      input_faces.change_to( NULL );
      PRINT_ERROR("Surface %d is a merged surface.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_face->id(), name );
    }
    input_faces.step();
  }
  input_faces.remove_all_with_value( NULL );
  if( input_faces.size() == 0 )
    return 0;

  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_faces.size()), query_output;
  CAST_LIST_TO_PARENT(input_faces, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  // Check for virtual geometry
  if ( contains_intermediate_geom(output_bodies))
  {
    PRINT_ERROR("%s surfaces on volumes containing virtual geometry\n"
      "       is not allowed.\n"
      "       Delete virtual geometry on these volumes before operation.\n",
      name);
    return 0;
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_faces, output_surfaces );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s surfaces on volumes containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: Consolidate common setup code for tweak curve operations
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefEdge*> &input_curves,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Curve*> &output_curves )
{
  if( input_curves.size() == 0 )
  {
    PRINT_ERROR( "No curves specified\n" );
    return 0;
  }

  //are any of the input curves merged? ...disallow that
  int i;
  for( i=input_curves.size(); i--; )
  {
    RefEdge *tmp_edge = input_curves.get();
    if( tmp_edge->is_merged() )
    {
      input_curves.change_to( NULL );
      PRINT_ERROR("Curve %d is a merged curve.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_edge->id(), name );
    }
    input_curves.step();
  }
  input_curves.remove_all_with_value( NULL );
  if( input_curves.size() == 0 )
    return 0;

  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_curves.size()), query_output;
  CAST_LIST_TO_PARENT(input_curves, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  // Check for virtual geometry
  if ( contains_intermediate_geom(output_bodies))
  {
    PRINT_ERROR("%s curves on surfaces containing virtual geometry\n"
      "       is not allowed.\n"
      "       Delete virtual geometry on these entities before operation.\n",
      name);
    return 0;
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_curves, output_curves );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s curves on entities containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: Consolidate common setup code for tweak vertex operations
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefVertex*> &input_vertices,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Point*> &output_points )
{
  if( input_vertices.size() == 0 )
  {
    PRINT_ERROR( "No vertices specified\n" );
    return 0;
  }

  //are any of the input vertices merged? ...disallow that
  int i;
  for( i=input_vertices.size(); i--; )
  {
    RefVertex *tmp_vert = input_vertices.get();
    if( tmp_vert->is_merged() )
    {
      input_vertices.change_to( NULL );
      PRINT_ERROR("Vertex %d is a merged vertex.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_vert->id(), name );
    }
    input_vertices.step();
  }
  input_vertices.remove_all_with_value( NULL );
  if( input_vertices.size() == 0 )
    return 0;




  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_vertices.size()), query_output;
  CAST_LIST_TO_PARENT(input_vertices, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  // Check for virtual geometry
  if ( contains_intermediate_geom(output_bodies))
  {
    PRINT_ERROR("%s vertices on surfaces containing virtual geometry\n"
      "       is not allowed.\n"
      "       Delete virtual geometry on these entities before operation.\n",
      name);
    return 0;
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_vertices, output_points );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s vertices on entities containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_chamfer( DLIList<RefEdge*> &ref_edge_list,
                                               double left_offset,
                                               DLIList<Body*> &new_body_list,
                                               double right_offset,
                                               CubitBoolean keep_old_body,
                                               CubitBoolean preview )
{
  // Make sure curves are not part of a sheet body
  // Get unique volumes that the curves are attached to
  DLIList<RefVolume*> ref_volume_list;
  int i;
  RefEdge *ref_edge_ptr;
  for( i=ref_edge_list.size(); i--; )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();
    DLIList<RefVolume*> tmp_ref_volume_list;
    ref_edge_ptr->ref_volumes( tmp_ref_volume_list );
    ref_volume_list.merge_unique( tmp_ref_volume_list );
  }

  RefVolume *ref_volume_ptr;
  for( i=ref_volume_list.size(); i--; )
  {
    ref_volume_ptr = ref_volume_list.get_and_step();
    if( ref_volume_ptr->is_sheet() )
    {
      PRINT_ERROR( "Cannot chamfer curves on sheet bodies\n" );
      return CUBIT_FAILURE;
    }
  }

  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Chamfering", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Do chamfering
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_chamfer( curve_list, left_offset, new_bodysm_list,
    right_offset, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If
//              multiple vertices are supplied, only one offset value is
//              allowed and the edges are not used.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus
GeometryModifyTool::tweak_chamfer( DLIList<RefVertex*> &ref_vertex_list,
                                   double offset1,
                                   DLIList<Body*> &new_body_list,
                                   RefEdge *edge1,
                                   double offset2,
                                   RefEdge *edge2,
                                   double offset3,
                                   RefEdge *edge3,
                                   CubitBoolean keep_old_body,
                                   CubitBoolean preview )
{

  if( offset1 <= 0.0 )
  {
    PRINT_ERROR( "Chamfer radius not specified.\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Point*> point_list(ref_vertex_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_vertex_list, "Chamfering", old_body_list, point_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  RefVertex *ref_vertex_ptr = ref_vertex_list.get();
  if( ref_vertex_list.size() == 1 && (edge1 || edge2 || edge3) )
  {
    // Make sure input edges are attached to the input vertex
    DLIList<RefEdge*> ref_edge_list;
    ref_vertex_ptr->ref_edges( ref_edge_list );

    if( edge1 && !ref_edge_list.is_in_list( edge1 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    if( edge2 && !ref_edge_list.is_in_list( edge2 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    if( edge3 && !ref_edge_list.is_in_list( edge3 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    // Make sure offsets supplied
    if( edge1 && offset1 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge1->id() );
      return CUBIT_FAILURE;
    }

    if( edge2 && offset2 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge2->id() );
      return CUBIT_FAILURE;
    }

    if( edge3 && offset3 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge3->id() );
      return CUBIT_FAILURE;
    }
  }

  if( point_list.size() > 1 && offset2 > 0.0 )
  {
    PRINT_ERROR( "Cannot supply multiple radii when chamfering multiple vertices.\n" );
    return CUBIT_FAILURE;
  }

  Curve *curve1 = NULL;
  Curve *curve2 = NULL;
  Curve *curve3 = NULL;

  if( edge1 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge1, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge1->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve1 = dynamic_cast<Curve*>(bridge);
  }

  if( edge2 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge2, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge2->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve2 = dynamic_cast<Curve*>(bridge);
  }

  if( edge3 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge3, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge3->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve3 = dynamic_cast<Curve*>(bridge);
  }

  // Do chamfering
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_chamfer( point_list, offset1, new_bodysm_list, curve1, offset2,
    curve2, offset3, curve3, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given curves on solid
//              bodies.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_fillet( DLIList<RefEdge*> &ref_edge_list,
                                              double radius,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  // Make sure curves are not part of a sheet body
  // Get unique volumes that the curves are attached to
  DLIList<RefVolume*> ref_volume_list;
  int i;
  RefEdge *ref_edge_ptr;
  for( i=ref_edge_list.size(); i--; )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();
    DLIList<RefVolume*> tmp_ref_volume_list;
    ref_edge_ptr->ref_volumes( tmp_ref_volume_list );
    ref_volume_list.merge_unique( tmp_ref_volume_list );
  }

  RefVolume *ref_volume_ptr;
  for( i=ref_volume_list.size(); i--; )
  {
    ref_volume_ptr = ref_volume_list.get_and_step();
    if( ref_volume_ptr->is_sheet() )
    {
      PRINT_ERROR( "Cannot fillet curves on sheet bodies\n" );
      return CUBIT_FAILURE;
    }
  }

  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Filleting", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Do filleting
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_fillet(curve_list, radius, new_bodysm_list, keep_old_body,
    preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given curves on a solid
//              body.  The fillet has a variable radius from the start to the
//              end of the curve.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_fillet( RefEdge *ref_edge_ptr,
                                              double start_radius,
                                              double end_radius,
                                              Body *&new_body_ptr,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  // Make sure curve is not part of a sheet body
  // Get unique volumes that the curves are attached to
  DLIList<RefVolume*> ref_volume_list;
  ref_edge_ptr->ref_volumes( ref_volume_list );

  int i;
  RefVolume *ref_volume_ptr;
  for( i=ref_volume_list.size(); i--; )
  {
    ref_volume_ptr = ref_volume_list.get_and_step();
    if( ref_volume_ptr->is_sheet() )
    {
      PRINT_ERROR( "Cannot fillet curves on sheet bodies\n" );
      return CUBIT_FAILURE;
    }
  }

  DLIList<Curve*> curve_list(1);
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  DLIList<RefEdge*> ref_edge_list(1);
  ref_edge_list.append( ref_edge_ptr );

  gme_ptr = tweak_setup( ref_edge_list, "Filleting", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  Curve *curve_ptr = curve_list.get();

  // Do filleting
  BodySM *new_bodysm_ptr;
  if( gme_ptr->tweak_fillet( curve_ptr, start_radius, end_radius, new_bodysm_ptr,
    keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    DLIList<BodySM*> new_bodysm_list;
    new_bodysm_list.append( new_bodysm_ptr );
    DLIList<Body*> new_body_list;
    if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
      return CUBIT_FAILURE;

    new_body_ptr = new_body_list.get();
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus
GeometryModifyTool::tweak_fillet( DLIList<RefVertex*> &ref_vertex_list,
                                  double radius,
                                  DLIList<Body*> &new_body_list,
                                  CubitBoolean keep_old_body,
                                  CubitBoolean preview )
{
  DLIList<Point*> point_list(ref_vertex_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_vertex_list, "Filleting", old_body_list, point_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Do filleting
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_fillet( point_list, radius, new_bodysm_list, keep_old_body,
    preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_move( DLIList<RefFace*>& ref_face_list,
                                            const CubitVector &delta,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Moving", old_body_list, surface_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Do move
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_move( surface_list, delta, new_bodysm_list, keep_old_body,
    preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && ref_face_list.is_in_list(ref_face))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_move( DLIList<RefEdge*>& ref_edge_list,
                                            const CubitVector &delta,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Moving", old_body_list, curve_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Do move
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_move( curve_list, delta, new_bodysm_list, keep_old_body,
    preview )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_offset( DLIList<RefFace*>& ref_face_list,
                                              double offset_distance,
                                              DLIList<Body*>& new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Offsetting", old_body_list, surface_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Do offset
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_offset( surface_list, offset_distance, new_bodysm_list,
    keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && ref_face_list.is_in_list(ref_face))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_offset( DLIList<RefEdge*>& ref_edge_list,
                                              double offset_distance,
                                              DLIList<Body*>& new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Offsetting", old_body_list, curve_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Do offset
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_offset( curve_list, offset_distance, new_bodysm_list,
    keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Function to remove surfaces from a body and then extend the
//              remaining surfaces to fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_remove( DLIList<RefFace*> &ref_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_adjoining,
                                              CubitBoolean keep_surface,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean individual,
                                              CubitBoolean preview )
{
  // Split things up if individual
  if (individual && extend_adjoining )
  {
    //build a surfae to volume map
    std::map<RefFace*, RefVolume*> surface_to_volume_map;

    int i;
    for (i = ref_face_list.size(); i--; )
    {
      RefFace *tmp_face = ref_face_list.get_and_step();
      RefVolume *tmp_vol = tmp_face->ref_volume();
      surface_to_volume_map.insert( std::map<RefFace*, RefVolume*>::value_type( tmp_face, tmp_vol));
    }

    DLIList<RefFace*> one_ref_face;
    CubitStatus total_rv = CUBIT_FAILURE;

    // Succeed if any one surface succeeds.
    for (i = ref_face_list.size(); i--; )
    {
      //make sure that the surface to remove is still in the body...
      //that it hasn't been removed from a previous tweak operation
      RefFace *tmp_face = ref_face_list.get_and_step();
      std::map<RefFace*, RefVolume*>::iterator tmp_iter;
      tmp_iter = surface_to_volume_map.find( tmp_face ); 
      RefVolume *tmp_vol = tmp_iter->second;
      DLIList<RefFace*> ref_face_list;
      tmp_vol->ref_faces( ref_face_list );
      if( !ref_face_list.move_to( tmp_face ) )
        continue;

      one_ref_face.clean_out();
      one_ref_face.append( tmp_face );
      int id = one_ref_face.get()->id();
      CubitStatus rv = this->tweak_remove(one_ref_face, new_body_list,
        extend_adjoining, keep_surface, keep_old_body, false, preview );
      if (rv)
      {
        new_body_list.uniquify_unordered();
        total_rv = CUBIT_SUCCESS;
        if( !preview )
          PRINT_INFO("Successfully removed Surface %d\n\n", id);
        else
          PRINT_INFO("Successfully removed Surface %d in preview\n\n", id);
        
        //see if we have a multi-volume body or multiple bodies
        //if so, we know the original volume was destroyed, so we 
        //cannot remove any more surfaces because the check above is 
        //highly likely to crash
        bool volume_destroyed = false;
        if( new_body_list.size() > 1 ) 
          volume_destroyed = true;
        else if( new_body_list.size() )
        {
          if( new_body_list.get()->num_ref_volumes() > 1 )
            volume_destroyed = true;
        }

        if( volume_destroyed == true  && i )
        {
          PRINT_WARNING("Unable to remove more surfaces because multiple bodies\n"
                        "       have been produced from removing surfaces individually\n" );
          return total_rv; 
        }
      }
      else
      {
        if( !preview )
          PRINT_INFO("Unable to remove Surface %d\n\n", id);
        else
          PRINT_INFO("Unable to remove Surface %d in preview\n\n", id);
      }
    }
    return total_rv;
  }

   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Body*> old_body_list;
   GeometryModifyEngine* gme_ptr;

   //collect all neighboring surfaces to those in the list
   int i,j;
   DLIList<RefFace*> neighboring_surfaces;
   for( i=ref_face_list.size(); i--; )
   {
     RefFace *tmp_face = ref_face_list.get_and_step();
     DLIList<RefEdge*> ref_edge_list;
     tmp_face->ref_edges( ref_edge_list );
     for( j=ref_edge_list.size(); j--; )
       ref_edge_list.get_and_step()->ref_faces( neighboring_surfaces );
   }

  //uniquify and add other surfaces
  neighboring_surfaces.uniquify_unordered();
  neighboring_surfaces += ref_face_list;

   gme_ptr = tweak_setup( ref_face_list, "Removing", old_body_list, surface_list );
   if (!gme_ptr)
     return CUBIT_FAILURE;

   // Do remove
   DLIList<BodySM*> new_bodysm_list;
   if( gme_ptr->tweak_remove( surface_list, new_bodysm_list, extend_adjoining,
     keep_surface, keep_old_body, preview ) == CUBIT_FAILURE )
     return CUBIT_FAILURE;

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if( ref_face && neighboring_surfaces.is_in_list(ref_face) )
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
          {
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
            {
              if(TopologyEntity* t = m->topology_entity())
              {
                entities_to_update.append(CAST_TO(t, RefEntity));
                //RefEntity *ref_ent = CAST_TO(t, RefEntity );
              }
            }
          }
        }
      }
    }
  }


  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

   return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Function to remove curves from a sheet body and then extend the
//              remaining curves or fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_remove( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Removing", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Do remove
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_remove( curve_list, new_bodysm_list, keep_old_body, preview )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes up to a target
//              surface or set of connected target surfaces.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefFace*> &ref_face_list,
                                              DLIList<RefFace*> &target_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Surface*> target_surf_list(target_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_face_list, "Tweaking", old_body_list, surface_list );
  if (!gme_ptr1)
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_face_list, "Tweaking", old_body_list2, target_surf_list );
  if (!gme_ptr2)
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target surfaces must belong to same geometry engine as tweaked surfaces.\n" );
    return CUBIT_FAILURE;
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr1->tweak_target( surface_list, target_surf_list, new_bodysm_list,
    reverse_flg, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    // some times the target surface gets modified
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && (ref_face_list.is_in_list(ref_face) ||
            target_face_list.is_in_list(ref_face)))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes up to a target
//              plane.
// Author     : Steve Storm
// Date       : 05/10/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefFace*> &ref_face_list,
                                              CubitPlane &plane,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size()+1);
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Tweaking", old_body_list, surface_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Create a target Surface from the plane (arbitrarily 10.0 X 10.0 in size)

  // Get corners of the surface
  CubitVector normal = plane.normal();
  CubitVector x, y;
  normal.orthogonal_vectors( x, y );
  CubitVector p1 = plane.point_on_plane();
  CubitVector p2, p3, p4;
  p1.next_point( x, 5.0, p1 );
  p1.next_point( y, 5.0, p1 );
  p1.next_point( -x, 10.0, p2 );
  p2.next_point( -y, 10.0, p3 );
  p3.next_point( x, 10.0, p4 );

  BodySM* bodysm_ptr = gme_ptr->planar_sheet( p1, p2, p3, p4 );
  if( !bodysm_ptr )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Surface*> target_surf_list;
  bodysm_ptr->surfaces( target_surf_list );
  if( !target_surf_list.size() )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_target( surface_list, target_surf_list, new_bodysm_list,
    reverse_flg, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Delete temporary sheet body
  bodysm_ptr->get_geometry_query_engine()->delete_solid_model_entities( bodysm_ptr );

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface or set of target surfaces.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<RefFace*> &target_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Surface*> target_surf_list(target_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr1 )
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_face_list, "Tweaking", old_body_list2, target_surf_list );
  if( !gme_ptr2 )
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target surface(s) must belong to same geometry engine as tweaked curves.\n" );
    return CUBIT_FAILURE;
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr1->tweak_target( curve_list, target_surf_list, new_bodysm_list,
    reverse_flg, keep_old, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target plane.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              CubitPlane &plane,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Create a target Surface from the plane (arbitrarily 10.0 X 10.0 in size)

  // Get corners of the surface
  CubitVector normal = plane.normal();
  CubitVector x, y;
  normal.orthogonal_vectors( x, y );
  CubitVector p1 = plane.point_on_plane();
  CubitVector p2, p3, p4;
  p1.next_point( x, 5.0, p1 );
  p1.next_point( y, 5.0, p1 );
  p1.next_point( -x, 10.0, p2 );
  p2.next_point( -y, 10.0, p3 );
  p3.next_point( x, 10.0, p4 );

  BodySM* bodysm_ptr = gme_ptr->planar_sheet( p1, p2, p3, p4 );
  if( !bodysm_ptr )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Surface*> target_surf_list;
  bodysm_ptr->surfaces( target_surf_list );
  if( !target_surf_list.size() )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_target( curve_list, target_surf_list, new_bodysm_list,
    reverse_flg, keep_old, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Delete temporary sheet body
  bodysm_ptr->get_geometry_query_engine()->delete_solid_model_entities( bodysm_ptr );

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a sheet body or bodies up to target
//              curves that are part of a sheet body.  The target is a surface
//              created by thickening the owning surface of the target curves.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<RefEdge*> &target_edge_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Curve*> target_curve_list(target_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr1 )
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_edge_list, "Tweaking", old_body_list2, target_curve_list );
  if( !gme_ptr2 )
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target curves must belong to same geometry engine as tweaked curves.\n" );
    return CUBIT_FAILURE;
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr1->tweak_target( curve_list, target_curve_list, new_bodysm_list,
    reverse_flg, keep_old, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if (!finish_sm_op( old_body_list, new_bodysm_list, new_body_list ))
    return CUBIT_FAILURE;

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}


// KGM
#if 0
bool GeometryModifyTool::contains_intermediate_geometry(DLIList<RefFace*>& ref_face_list) const
{
  for (int i = 0; i < ref_face_list.size(); i++)
    if (GeometryQueryTool::instance()->
        contains_intermediate_geometry(ref_face_list.next(i)))
      return true;

  return false;
}
#endif

//the following surface tool operations added by Tyronne Lim (CAT) ********************
CubitStatus GeometryModifyTool::create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body,
                                                    DLIList<DLIList<CubitVector*>*> &vec_lists_u,
                                                    DLIList<DLIList<CubitVector*>*> &vec_lists_v,
                                                    double net_tol, CubitBoolean heal )
{
   GeometryModifyEngine* GMEPtr = get_engine(ref_face_list.get());
   return GMEPtr->create_net_surface( ref_face_list, new_body, vec_lists_u, vec_lists_v, net_tol, heal );
}

CubitStatus GeometryModifyTool::create_net_surface( DLIList<RefEdge*>& u_curves, DLIList<RefEdge*>& v_curves,
                                                    double net_tol, CubitBoolean heal )
{
  DLIList<TopologyBridge*> bridge_list;
  DLIList<TopologyEntity*> entity_list;
  CAST_LIST_TO_PARENT( u_curves, entity_list );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  DLIList<Curve*> curves_in_u(bridge_list.size());
  CAST_LIST( bridge_list, curves_in_u, Curve );

  bridge_list.clean_out();
  entity_list.clean_out();
  CAST_LIST_TO_PARENT( v_curves, entity_list );

  GeometryModifyEngine* dummy_GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! dummy_GME_ptr || dummy_GME_ptr != GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  DLIList<Curve*> curves_in_v(bridge_list.size());
  CAST_LIST( bridge_list, curves_in_v, Curve );

  BodySM *new_body = NULL;
  if( !GME_ptr->create_net_surface( curves_in_u, curves_in_v, new_body, net_tol, heal ) )
    return CUBIT_FAILURE;

  if( GeometryQueryTool::instance()->make_Body( new_body ) )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::create_offset_surface( RefFace* ref_face_ptr, double offset_distance )
{
  TopologyBridge *bridge_ptr = NULL;
  GeometryModifyEngine* GMEPtr = get_engine(ref_face_ptr, &bridge_ptr );

  Surface *tmp_surf = NULL;
  tmp_surf = CAST_TO( bridge_ptr, Surface );

  BodySM *new_body;
  if( !GMEPtr->create_offset_surface( tmp_surf, new_body, offset_distance ) )
    return CUBIT_FAILURE;

  if( GeometryQueryTool::instance()->make_Body( new_body ) )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::create_offset_body( Body *body_ptr, Body *&new_body, double offset_distance )
{
  GeometryModifyEngine* GMEPtr = get_engine(body_ptr);

  BodySM *body_sm = body_ptr->get_body_sm_ptr();
  if (!body_sm)
  {
    PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", body_ptr->id());
    return CUBIT_FAILURE;
  }

  BodySM *new_body_sm;
  if( !GMEPtr->create_offset_body( body_sm, new_body_sm, offset_distance ) )
    return CUBIT_FAILURE;

  new_body = GeometryQueryTool::instance()->make_Body( new_body_sm );

  if( new_body )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::create_skin_surface( DLIList<RefEdge*>& ref_edges, Body*& new_body )
{
  DLIList<TopologyBridge*> bridge_list;
  DLIList<TopologyEntity*> entity_list;
  CAST_LIST_TO_PARENT( ref_edges, entity_list );

  if (ref_edges.size() < 2)
  {
     PRINT_ERROR("Must specify at least 2 curves to create a skinned surface.\n");
     return CUBIT_FAILURE;
  }

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  DLIList<Curve*> curves_to_skin(bridge_list.size());
  CAST_LIST( bridge_list, curves_to_skin, Curve );

  BodySM *new_body_sm = NULL;
  if( !GME_ptr->create_skin_surface( curves_to_skin, new_body_sm ) )
    return CUBIT_FAILURE;

  new_body = NULL;
  new_body = GeometryQueryTool::instance()->make_Body( new_body_sm );

  if( new_body )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}


CubitStatus GeometryModifyTool::loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
                                                      RefFace *face2, const double &takeoff2,
                                                      Body*& new_body,
                                                      CubitBoolean arc_length_option, CubitBoolean twist_option,
                                                      CubitBoolean align_direction, CubitBoolean perpendicular,
                                                      CubitBoolean simplify_option)
{

    DLIList<RefFace*> loft_faces;
    loft_faces.append(face1);
    loft_faces.append(face2);
    DLIList<Surface*> loft_surfaces;

    // Get engine and correspoding geom entities
    GeometryModifyEngine* result_ptr;
    result_ptr = common_modify_engine( loft_faces, loft_surfaces );
    if (!result_ptr)
    {
        PRINT_ERROR("Loft surfaces on volumes containing surfaces from different\n"
            "       geometry engines is not allowed.\n");
        return CUBIT_FAILURE;
    }

    if(2!=loft_surfaces.size())
        return CUBIT_FAILURE;

    loft_surfaces.reset();
    BodySM* new_body_sm = 0;

    CubitStatus result = result_ptr->loft_surfaces_to_body(
        loft_surfaces.get_and_step(),
        takeoff1,
        loft_surfaces.get_and_step(),
        takeoff2,
        new_body_sm,
        arc_length_option,
        twist_option,
        align_direction,
        perpendicular,
        simplify_option);

    if(result && new_body_sm)
    {
        new_body = GeometryQueryTool::instance()->make_Body(new_body_sm);
    }

   return result;
}

CubitStatus GeometryModifyTool::create_surface( DLIList<CubitVector*>& vec_list, Body *&new_body,
                                                RefFace *ref_face_ptr,CubitBoolean project_points )
{
   GeometryModifyEngine* GMEPtr = gmeList.get();
   if( ref_face_ptr )
   {
     if( GMEPtr != get_engine(ref_face_ptr) )
     {
       PRINT_ERROR("Geometry engine of Surface %d is not the active geometry engine.\n", ref_face_ptr->id() );
       PRINT_INFO("      Use command \"Set Geometry Engine ...\" to set to correct engine.\n");
       return CUBIT_FAILURE;
     }
   }
   BodySM* body_sm = NULL;
   Surface *project_to_surface = NULL;
   if( ref_face_ptr )
     project_to_surface = ref_face_ptr->get_surface_ptr();
   CubitStatus stat = GMEPtr->create_surface( vec_list, body_sm, project_to_surface, project_points );

   if( stat == CUBIT_FAILURE )
     return stat;
   if( body_sm )
     new_body = GeometryQueryTool::instance()->make_Body(body_sm);
   if( !new_body )
     return CUBIT_FAILURE;
   else
     return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::create_weld_surface( CubitVector &root,
                                                    RefFace *ref_face1,
                                                    double leg1,
                                                    RefFace *ref_face2,
                                                    double leg2,
                                                    Body *&new_body )
{
    GeometryModifyEngine* GMEPtr;
    DLIList<RefFace*> ref_faces;
    ref_faces.append(ref_face1);
    ref_faces.append(ref_face2);
    DLIList<Surface*> surfaces;

    GMEPtr = common_modify_engine(ref_faces, surfaces);
    if (!GMEPtr)
    {
        PRINT_ERROR("Loft surfaces on volumes containing surfaces from different\n"
            "       geometry engines is not allowed.\n");
        return CUBIT_FAILURE;
    }

    surfaces.reset();
    BodySM* new_body_sm = 0;
    CubitStatus result = GMEPtr->create_weld_surface(
        root,
        surfaces.get_and_step(),
        leg1,
        surfaces.get_and_step(),
        leg2,
        new_body_sm );

    if(result && new_body_sm)
    {
        new_body = GeometryQueryTool::instance()->make_Body(new_body_sm);
    }
    return result;
}
//end of surface tool operations ******************************************************

//-------------------------------------------------------------------------
// Purpose       : Remove entity names from dead entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/08/03
//-------------------------------------------------------------------------
void GeometryModifyTool::remove_dead_entity_names( RefEntity* entity ) const
{
  TopologyEntity* topo_ent = dynamic_cast<TopologyEntity*>(entity);
  if (topo_ent->bridge_manager()->topology_bridge() == NULL)
    RefEntityName::instance()->remove_refentity_name( entity, CUBIT_TRUE );

  DLIList<RefEntity*> children;
  entity->get_child_ref_entities(children);
  children.last();
  for (int i = children.size(); i--; )
  {
    //PRINT_INFO("Removing dead entity on %s %d\n", children.get()->class_name(), children.get()->id() );
    remove_dead_entity_names( children.step_and_get() );
  }
}

//-------------------------------------------------------------------------
// Purpose       : Destroy or update modified body, as appropriate.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/08/03
//-------------------------------------------------------------------------
Body* GeometryModifyTool::update_body( Body* body ) const
{
  BodySM* body_sm = body->get_body_sm_ptr();
  if (body_sm)
    return GeometryQueryTool::instance()->make_Body(body_sm);

  GeometryQueryTool::instance()->destroy_dead_entity(body);
  return 0;
}

CubitStatus GeometryModifyTool::tolerant_imprint( DLIList<Body*> &bodies,
                                      DLIList<Body*> &new_bodies, bool merge )
{
  //make sure all bodies are from the same modify engine
  DLIList<BodySM*> body_sm_list;
  GeometryModifyEngine* gme = common_modify_engine(bodies, body_sm_list);
  if ( !gme )
  {
    PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                "from different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }


  //make sure that merge tolerance is not inapproiate for model
  int i;
  CubitBox bounding_box( CubitVector(0,0,0),
                         CubitVector(CUBIT_DBL_MAX, 
                                     CUBIT_DBL_MAX,
                                     CUBIT_DBL_MAX ) );
  for( i=bodies.size(); i--; )
  {
    CubitBox tmp_box = bodies.get_and_step()->bounding_box();
    if(bounding_box.max_x() == CUBIT_DBL_MAX)
      bounding_box = tmp_box;
    else if( tmp_box.diagonal().length_squared() < 
        bounding_box.diagonal().length_squared() )
      bounding_box = tmp_box;
  }
   
  //get the merge tolerance 
  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;
 
  //if the merge tolerance is greater than 1/10th the length of the 
  //diagonal of the bounding box of the smallest volume, fail!
  double tenth_smallest_bbox = 0.1*(bounding_box.diagonal().length()); 
  if( tolerance > tenth_smallest_bbox ) 
  {
    PRINT_ERROR("Merge tolerance is set excessively high.  Must be lowner than %lf\n", 
                 tenth_smallest_bbox );
    PRINT_INFO("       (Merge tolerance must be less than than 1/10th of the diagonal\n"
                        "of the bounding box of the smallest volume)\n");
    return CUBIT_FAILURE;
  }

  body_sm_list.clean_out();
  DLIList<Body*> old_body_list;
  old_body_list += bodies;
  bodies.reset();
  for( i=bodies.size(); i--; )
    body_sm_list.append( bodies.get_and_step()->get_body_sm_ptr() );

  int process_composites = 0;
  if(contains_composites(bodies))
    process_composites = 1;

  if(process_composites)
  {
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sm_list);
    // This must be done after pushing the vg atts because it uses them.
    push_imprint_attributes_before_modify(body_sm_list);
  }

  DLIList<BodySM*> new_body_list;
  DLIList<TopologyBridge*> new_tbs, att_tbs;
  CubitStatus result = gme->tolerant_imprint( body_sm_list, new_body_list, &new_tbs, &att_tbs ) ;

  if(process_composites)
  {
    // Analyze the results and adjust virtual attributes as necessary.
    GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
      new_body_list, bodies);

    // Clean up attributes.
    remove_imprint_attributes_after_modify(body_sm_list, new_body_list);

    // Restore the virtual geometry.
    restore_vg_after_modify(new_body_list, bodies);
  }

  if(result == CUBIT_FAILURE)
  {
    if(process_composites)
      do_attribute_cleanup();
    return result;
  }

  result = finish_sm_op( bodies, body_sm_list, new_bodies );

  if(process_composites)
    do_attribute_cleanup();

  if(!result)
    return CUBIT_FAILURE;

  if( merge )
    MergeTool::instance()->merge_bodies( bodies );

  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::remove_curve_slivers( DLIList<Body*> &bodies,
                                              double lengthlimit )
{
  DLIList<BodySM*> body_sm_list;
  GeometryModifyEngine* gme = common_modify_engine(bodies, body_sm_list);
  if ( !gme )
  {
    PRINT_ERROR("Curve sliver removal only supported on ACIS geometry\n");
    return CUBIT_FAILURE;
  }

  CubitStatus status = CUBIT_FAILURE;
  DLIList<BodySM*> modified_bodies;
  int i;
  for( i=body_sm_list.size(); i--; )
  {
    BodySM *tmp_body_sm = body_sm_list.get_and_step();
    if( gme->remove_curve_slivers( tmp_body_sm, lengthlimit ) == CUBIT_SUCCESS )
    {
      modified_bodies.append( tmp_body_sm );
      status = CUBIT_SUCCESS;
    }
  }

  DLIList<Body*> dummy_list;
  if(!finish_sm_op( bodies, modified_bodies, dummy_list))
    return CUBIT_FAILURE;

  return status;
}

void GeometryModifyTool::determine_solutions_for_eliminating_small_surface(RefFace *face,
                                                                           DLIList<CubitString> &display_strings,
                                                                           DLIList<CubitString> &command_strings)
{
  // Build strings for potential composite solutions.
  // Build strings for potential remove_face solutions.
  // Build strings for potential remove_topology solutions.
  // Build strings for potential regularize solutions.
}

CubitStatus GeometryModifyTool::prepare_surface_sweep(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& sweep_vector,
                              bool sweep_perp,
                              bool through_all,
                              bool outward,
                              bool up_to_next,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              BodySM* &cutting_tool_ptr ,
                              const CubitVector* point,
                              double *angle)
{
  GeometryModifyEngine* gme = get_engine(blank_bodies.get());

  if(surfaces.size() == 0 )
    return CUBIT_FAILURE;

  DLIList<GeometryEntity*> ref_ent_list;
  Surface * temp_face = NULL;
  for(int i = 0; i < surfaces.size(); i++)
    {
      //copy the faces before sweep
      temp_face = gme->make_Surface(surfaces.get_and_step());
      if (temp_face)
        ref_ent_list.append((GeometryEntity*)temp_face);
    }

  BodySM* to_body = NULL;
  CubitStatus stat = CUBIT_SUCCESS;
  if(up_to_next && blank_bodies.size() > 1) //unite all bland_bodies
    {
       DLIList<BodySM*> newBodies;
       DLIList<BodySM*> copied_bodies;
       for(int i = 0; i < blank_bodies.size(); i++)
         copied_bodies.append(gme->copy_body(blank_bodies.get_and_step()));

       stat = gme->unite(copied_bodies, newBodies);
       if(stat == CUBIT_FAILURE)
         {
           PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
           PRINT_INFO("Try the 'stop surface <id>' option instead\n");
           return stat;
         }
       to_body = newBodies.get();
    }

  else if(up_to_next && blank_bodies.size() == 1)
    to_body = gme->copy_body(blank_bodies.get());

  DLIList<BodySM*> swept_bodies;
  if (point && angle) //sweep_surface_rotated
    stat = gme->sweep_rotational(ref_ent_list,swept_bodies,*point,
                           sweep_vector, *angle,0, 0.0,0,false,false,
                           false,stop_surf, to_body);  
	
  else
  {
    CubitVector tmp_sweep_vector = sweep_vector;

    //get model bbox info...will scale sweep vector by its diagonal
    //so that we go far enough
    if( through_all || stop_surf || up_to_next)
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.normalize();
      tmp_sweep_vector*=(2*bounding_box.diagonal().length());
    }
    
    //see if we're sweeping along a specified curve
    if( curve_to_sweep_along )
      {
        DLIList<Curve*> curves_to_sweep_along;
        curves_to_sweep_along.append(curve_to_sweep_along);
        stat = gme->sweep_along_curve(ref_ent_list, swept_bodies,
                                 curves_to_sweep_along, 0.0,0,false,stop_surf,
                                 to_body);
      }

    else if (sweep_perp )
      stat = gme->sweep_perpendicular(ref_ent_list, swept_bodies,
                                 tmp_sweep_vector.length(),0.0,0,outward,false,
                                 stop_surf, to_body);
    else
      stat = gme->sweep_translational(ref_ent_list, swept_bodies,
                                 tmp_sweep_vector,0.0,0, false, false, stop_surf,
                                 to_body);
  }

  if(stat == CUBIT_FAILURE && swept_bodies.size() == 0)
    {
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              gme->get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }

       return stat;
    }

  //if there are more than 1, unite them all
  DLIList<BodySM*>  newBodies;
  if (swept_bodies.size() > 1)
     stat = gme->unite(swept_bodies, newBodies);
  else
     newBodies = swept_bodies;

  if(stat == CUBIT_FAILURE || newBodies.size()!= 1)
    {
       PRINT_ERROR("webcut tool body is not created from acis.\n");
       //delete the swept_bodies
       BodySM* tmp_body = NULL;
       for (int i = swept_bodies.size(); i--;)
         {
           tmp_body= swept_bodies.get_and_step();
           if (tmp_body)
             gme->get_gqe()->delete_solid_model_entities(tmp_body);
         }

       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              gme->get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }
       return CUBIT_FAILURE;
    }
  
    cutting_tool_ptr = newBodies.get();
    return stat;
}
