//-------------------------------------------------------------------------
// Filename      : TopologyEntity.cpp
//
// Purpose       : Implementation of the TopologyEntity class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/31/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitMessage.hpp"

#include "CastTo.hpp"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"

#include "ModelQueryEngine.hpp"

#include "Body.hpp"
#include "DLIList.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoVolume.hpp"
#include "CoFace.hpp"
#include "CoEdge.hpp"
#include "CoVertex.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Casting query functions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
#define DECLARE_TOPO_ENT_QUERY_FUNC( TYPE, NAME, QUERYTYPE ) \
 CubitStatus TopologyEntity::NAME(DLIList<TYPE*>& list) \
 { \
   DLIList<ModelEntity*> temp_list ; \
   CubitStatus result; \
   ModelQueryEngine *const mqe = ModelQueryEngine::instance(); \
   \
   result = mqe->query_model( *this, DagType::QUERYTYPE(), temp_list ); \
   if (result == CUBIT_FAILURE) \
   { \
      PRINT_ERROR("In TopologyEntity::" #NAME "\n"); \
      PRINT_ERROR("       Query failed for unknown reason.\n"); \
      return CUBIT_FAILURE; \
   } \
   \
   temp_list.reset(); \
   for (int i = temp_list.size(); i--; ) \
   { \
     TYPE* ptr = dynamic_cast<TYPE*>(temp_list.get_and_step()); \
     assert(!!ptr); \
     list.append(ptr); \
   } \
   \
   return CUBIT_SUCCESS; \
}



//-------------------------------------------------------------------------
// Purpose       : Default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------

TopologyEntity::TopologyEntity() 
    : bridgeMan(this)
{
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/10/96
//-------------------------------------------------------------------------

TopologyEntity::~TopologyEntity() 
{
}


GeometryQueryEngine*
TopologyEntity::get_geometry_query_engine() const
{
  TopologyBridge* bridge = 
    const_cast<TopologyEntity*>(this)->bridgeMan.topology_bridge();
  if (bridge)
    return bridge->get_geometry_query_engine();
  else
    return NULL;
}


DECLARE_TOPO_ENT_QUERY_FUNC(      Body,       bodies,       body_type )
DECLARE_TOPO_ENT_QUERY_FUNC(  CoVolume,   co_volumes,  co_volume_type )
DECLARE_TOPO_ENT_QUERY_FUNC( RefVolume,  ref_volumes, ref_volume_type )
DECLARE_TOPO_ENT_QUERY_FUNC(     Shell,       shells,      shell_type )
DECLARE_TOPO_ENT_QUERY_FUNC(    CoFace,     co_faces,    co_face_type )
DECLARE_TOPO_ENT_QUERY_FUNC(   RefFace,    ref_faces,   ref_face_type )
DECLARE_TOPO_ENT_QUERY_FUNC(      Loop,        loops,       loop_type )
DECLARE_TOPO_ENT_QUERY_FUNC(    CoEdge,     co_edges,    co_edge_type )
DECLARE_TOPO_ENT_QUERY_FUNC(   RefEdge,    ref_edges,   ref_edge_type )
DECLARE_TOPO_ENT_QUERY_FUNC(     Chain,       chains,      chain_type )
DECLARE_TOPO_ENT_QUERY_FUNC(  CoVertex,  co_vertices,  co_vertex_type )
DECLARE_TOPO_ENT_QUERY_FUNC( RefVertex, ref_vertices, ref_vertex_type )


RefVertex *TopologyEntity::ref_vertex()
{
  DLIList<RefVertex*> verts;
  ref_vertices(verts);
  if (verts.size() > 0) return verts.get();
  else return NULL;
}

RefEdge *TopologyEntity::ref_edge()
{
  DLIList<RefEdge*> edges;
  ref_edges(edges);
  if (edges.size() > 0) return edges.get();
  else return NULL;
}

RefFace *TopologyEntity::ref_face()
{
  DLIList<RefFace*> faces;
  ref_faces(faces);
  if (faces.size() > 0) return faces.get();
  else return NULL;
}

RefVolume *TopologyEntity::ref_volume()
{
  DLIList<RefVolume*> volumes;
  ref_volumes(volumes);
  if (volumes.size() > 0) return volumes.get();
  else return NULL;
}

Body *TopologyEntity::body()
{
  DLIList<Body*> these_bodies;
  this->bodies(these_bodies);
  if (these_bodies.size() > 0) return these_bodies.get();
  else return NULL;
}


CoEdge *TopologyEntity::co_edge() 
{
  DLIList<CoEdge*> these_co_edges;
  co_edges(these_co_edges);
  if (these_co_edges.size() > 0) return these_co_edges.get();
  else return NULL;
}

int TopologyEntity::num_loops() 
{
  DLIList<Loop*> these_loops;
  loops(these_loops);
  return these_loops.size();
}


int TopologyEntity::num_ref_volumes() 
{
  DLIList<RefVolume*> these_ref_volumes;
  ref_volumes(these_ref_volumes);
  return these_ref_volumes.size();
}

int TopologyEntity::num_ref_faces() 
{
  DLIList<RefFace*> these_ref_faces;
  ref_faces(these_ref_faces);
  return these_ref_faces.size();
}

int TopologyEntity::num_ref_edges() 
{
  DLIList<RefEdge*> these_ref_edges;
  ref_edges(these_ref_edges);
  return these_ref_edges.size();
}

int TopologyEntity::num_ref_vertices() 
{
  DLIList<RefVertex*> these_ref_vertices;
  ref_vertices(these_ref_vertices);
  return these_ref_vertices.size();
}


  //- return the number of connected entities of the given type

CubitBoolean TopologyEntity::is_directly_related( TopologyEntity* entity )
{
    // NULL?
  if ( !entity )
    return CUBIT_FALSE;
  
    // self?
  if (entity == this)
    return CUBIT_TRUE;  

    // same type but not self?
  if ( dag_type() == entity->dag_type() )  
    return CUBIT_FALSE;
  
    // Get the entities of the right type.
  DLIList<ModelEntity*> kin_folk;
  CubitStatus result = ModelQueryEngine::instance()->
    query_model( *this, entity->dag_type(), kin_folk );
    //only fails if types are the same, caught above. 
  assert(result != CUBIT_FAILURE); 

    // Search entities for passed in entity.
  ModelEntity* model_entity = CAST_TO(entity, ModelEntity) ;
  return kin_folk.is_in_list( model_entity );
}

CubitStatus TopologyEntity::set_topology_bridge(TopologyBridge* TB_ptr)
{
  assert (TB_ptr != NULL);
  CubitStatus rv = CUBIT_SUCCESS;
  
  BridgeManager* bridges = bridge_manager();
  bridges->remove_all_bridges();
  bridges->add_bridge(TB_ptr);

  return rv;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

