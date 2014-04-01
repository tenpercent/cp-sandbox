//- Class:          CubitAttrib
//- Owner:          Greg Nielson
//- Description:    implementation of the CubitAttrib class.
//- Checked By:
//- Version:

#include "CastTo.hpp"
#include "CubitAttrib.hpp"
#include "CubitAttribUser.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CADeferredAttrib.hpp"
#include "RefEntity.hpp"
#include "DLIList.hpp"
#include "RefEntityFactory.hpp"
#include "MergeTool.hpp"
#include "ModelQueryEngine.hpp"
#include "GeometryQueryTool.hpp"


    

CubitAttrib::CubitAttrib(RefEntity *attrib_owner)
{
  attribOwnerEntity = attrib_owner;
  hasActuated = CUBIT_FALSE;
  hasUpdated = CUBIT_FALSE;
  hasWritten = CUBIT_FALSE;
  deleteAttrib = CUBIT_FALSE;
  nextAttrib = NULL;
  
    // add this to the owner
  if (attrib_owner) attrib_owner->add_cubit_attrib(this);
}

CubitAttrib::~CubitAttrib() 
{
  if( !hasActuated )
    CADeferredAttrib::remove_unactuated_ca( this );
}

CubitStatus CubitAttrib::actuate_list(DLIList<RefEntity*> entity_list)
{
  RefEntity * ref_ent;
  for(int i = entity_list.size(); i > 0; i--)
  {
    ref_ent = entity_list.get_and_step();
    ref_ent->actuate_cubit_attrib(CA_ENTITY_NAME);
    ref_ent->actuate_cubit_attrib ( CA_UNIQUE_ID );
    ref_ent->actuate_cubit_attrib(CA_SIZING_FUNCTION_SKELETON);    
    ref_ent->actuate_cubit_attrib(CA_MESH_INTERVAL);
    ref_ent->actuate_cubit_attrib(CA_GROUP);
    ref_ent->actuate_cubit_attrib(CA_GENESIS_ENTITY);
//    ref_ent->actuate_cubit_attrib ( CA_ENTITY_ID );
    ref_ent->actuate_cubit_attrib ( CA_MESH_SCHEME );
    ref_ent->actuate_cubit_attrib ( CA_SMOOTH_SCHEME );
    ref_ent->actuate_cubit_attrib ( CA_PARTITION_VG );
    ref_ent->actuate_cubit_attrib ( CA_COMPOSITE_VG );
    ref_ent->actuate_cubit_attrib ( CA_VIRTUAL_VG );
    ref_ent->actuate_cubit_attrib(CA_MERGE_PARTNER);
    ref_ent->actuate_cubit_attrib(CA_DEFERRED_ATTRIB);
    ref_ent->actuate_cubit_attrib(CA_MESH_CONTAINER);
    ref_ent->actuate_cubit_attrib(CA_BODIES);
    ref_ent->actuate_cubit_attrib ( CA_ENTITY_ID );
    ref_ent->actuate_cubit_attrib(CA_ENTITY_COLOR);
//#ifdef CAT
    ref_ent->actuate_cubit_attrib(CA_VERTEX_FORCE);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_FORCE);
    ref_ent->actuate_cubit_attrib(CA_CURVE_FORCE);
    ref_ent->actuate_cubit_attrib(CA_VERTEX_DISPLACEMENT);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_DISPLACEMENT);
    ref_ent->actuate_cubit_attrib(CA_CURVE_DISPLACEMENT);
    ref_ent->actuate_cubit_attrib(CA_VOLUME_DISPLACEMENT);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_PRESSURE);
    ref_ent->actuate_cubit_attrib(CA_CURVE_PRESSURE);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_TEMPERATURE);
    ref_ent->actuate_cubit_attrib(CA_CURVE_TEMPERATURE);
    ref_ent->actuate_cubit_attrib(CA_VERTEX_TEMPERATURE);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_HEATFLUX);
    ref_ent->actuate_cubit_attrib(CA_CURVE_HEATFLUX);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_CONVECTION);
    ref_ent->actuate_cubit_attrib(CA_CURVE_CONVECTION);
    ref_ent->actuate_cubit_attrib(CA_SURFACE_CONTACT);
    ref_ent->actuate_cubit_attrib(CA_CURVE_CONTACT);
    ref_ent->actuate_cubit_attrib(CA_COORD_SYS);
    ref_ent->actuate_cubit_attrib(CA_PROPERTY_BLOCK);
    ref_ent->actuate_cubit_attrib(CA_MATERIAL_BLOCK);
//#endif
    ref_ent->actuate_cubit_attrib(CA_MERGE_STATUS);
  } 
  return CUBIT_SUCCESS;
}
  
void CubitAttrib::has_written(CubitBoolean set_has_written)
{
  hasWritten = set_has_written;

    // if the written flag is being set to true, reset the hasUpdated flag
  if (CUBIT_TRUE == hasWritten)
    hasUpdated = CUBIT_FALSE;
}

CubitBoolean CubitAttrib::has_written() const
{return hasWritten;}

void CubitAttrib::remove_attribute()
{
  if (has_written())
    attribOwnerEntity->remove_attrib_geometry_entity(this);
}

void CubitAttrib::add_attribute()
{
  attribOwnerEntity->add_cubit_attrib(this);
}

int CubitAttrib::equivalent(CubitSimpleAttrib* csa_ptr)
{
    //- return true if the csa and this are equivalent
   CubitSimpleAttrib* this_csa_ptr = cubit_simple_attrib();
       
   CubitBoolean equivalent =
       CubitSimpleAttrib::equivalent(csa_ptr, this_csa_ptr);
   delete this_csa_ptr;
   return equivalent;
}

void CubitAttrib::print()
{
    // print some details about this attrib
  PRINT_INFO("Attrib type %s, Owner = %s %d, Actuated=%d, Updated=%d, "
             "Written=%d, Delete=%d\n",
             att_internal_name(),
             (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
             (attribOwnerEntity ? attribOwnerEntity->id() : 0),
             hasActuated, hasUpdated,
             hasWritten, deleteAttrib);

}

