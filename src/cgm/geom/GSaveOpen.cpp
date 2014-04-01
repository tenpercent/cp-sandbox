#include "GSaveOpen.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "RefEntityFactory.hpp"

// Initialize Globals
int GSaveOpen::gsoSetsIds = 0;
int GSaveOpen::gsoIncBodyId = 0;
int GSaveOpen::gsoIncRefVolumeId = 0;
int GSaveOpen::gsoIncRefFaceId = 0;
int GSaveOpen::gsoIncRefEdgeId = 0;
int GSaveOpen::gsoIncRefVertexId = 0;
int GSaveOpen::gsoErrorCount = 0;
DLIList<int> GSaveOpen::gsoErrorIdList;

GSaveOpen::GSaveOpen()
{
  RefEntityFactory *ref = RefEntityFactory::instance();

  gsoSetsIds = 1;
  gsoIncBodyId = ref->current_body_id();
  gsoIncRefVolumeId = ref->current_volume_id();
  gsoIncRefFaceId = ref->current_face_id();
  gsoIncRefEdgeId = ref->current_edge_id();
  gsoIncRefVertexId = ref->current_vertex_id();
  gsoErrorCount = 0;
}

GSaveOpen::~GSaveOpen() 
{
  // Reset globals
  gsoSetsIds = 0;
  gsoIncBodyId = 0;
  gsoIncRefVolumeId = 0;
  gsoIncRefFaceId = 0;
  gsoIncRefEdgeId = 0;
  gsoIncRefVertexId = 0;
  gsoErrorCount = 0;
  gsoErrorIdList.clean_out();
}

int 
GSaveOpen::get_id_inc( RefEntity *entity )
{
  const type_info* my_type = &entity->entity_type_info();

  if( *my_type == typeid(Body) )
    return gsoIncBodyId;
  else if( *my_type == typeid(RefVolume) )
    return gsoIncRefVolumeId;
  else if( *my_type == typeid(RefFace) )
    return gsoIncRefFaceId;
  else if( *my_type == typeid(RefEdge) )
    return gsoIncRefEdgeId;
  else if( *my_type == typeid(RefVertex) )
    return gsoIncRefVertexId;

  return -1; // This function doesn't understand the type
}


