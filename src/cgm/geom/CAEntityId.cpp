//- File:           CAEntityId.cpp
//- Owner:          Dong Zhu
//- Description:    Cubit Attribute for entity ids.
//- Checked By:
//- Version:


#include "CAEntityId.hpp"
#include "CAMergePartner.hpp"
#include "TDUniqueId.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "CastTo.hpp"
#include "MergeTool.hpp"
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GSaveOpen.hpp"
#include "CADeferredAttrib.hpp"

CubitAttrib* CAEntityId_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAEntityId *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAEntityId(entity);
  }
  else
  {
    new_attrib = new CAEntityId(entity, p_csa);
  }

  return new_attrib;
}

CAEntityId::CAEntityId(RefEntity* new_attrib_owner,
                       CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  assert ( csa_ptr != 0 );
   PRINT_DEBUG_94( "Creating ENTITY_ID attribute from CSA for %s %d\n",
      (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
      (attribOwnerEntity ? attribOwnerEntity->id() : 0));
   
   DLIList<int*> *i_list = csa_ptr->int_data_list();

   assert(i_list && i_list->size() > 0);
   i_list->reset();
   entityId = *( i_list->get_and_step() );

   if (i_list->size() > 1) {
     boundingUid = *(i_list->get_and_step());
     boundingSense = (CubitSense) *(i_list->get_and_step());
   }
   else {
     boundingUid = CUBIT_INT_MIN;
     boundingSense = CUBIT_UNKNOWN;
   }

   DLIList<double*> *d_list = csa_ptr->double_data_list();
   if (d_list && d_list->size() > 0) {
     RefEdge *edge = CAST_TO(new_attrib_owner, RefEdge);
     assert(d_list->size() == 3 && 
            edge && edge->start_vertex() == edge->end_vertex());
     d_list->reset();
     boundingXYZ = new CubitVector(*d_list->get_and_step(),
                                   *d_list->get_and_step(),
                                   *d_list->get());
   }
   else boundingXYZ = NULL;
}

CAEntityId::CAEntityId(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  entityId = 0;
  boundingUid = CUBIT_INT_MIN;
  boundingSense = CUBIT_UNKNOWN;
  boundingXYZ = NULL;
  
  PRINT_DEBUG_94( "Creating ENTITY_ID attribute for %s %d\n",
              (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
              (attribOwnerEntity ? attribOwnerEntity->id() : 0));
}

CAEntityId::~CAEntityId()
{
}


CubitStatus CAEntityId::actuate()
{
  if ( hasActuated)
    return CUBIT_SUCCESS;
  
  if ( !attribOwnerEntity )
    return CUBIT_FAILURE;

  int this_owner_id = attribOwnerEntity->id();

  deleteAttrib = CUBIT_FALSE;

  // First, check to see if we are importing a Cubit file, in which
  // case we may need to add an increment to the CAEntityId id, for
  // consistency.  In the Cubit file restore, mesh gets reattached to
  // geometry based on the geometry id.  When importing into an 
  // existing session of Cubit (where a model already exists), we 
  // need to use a consistent id increment so that the mesh gets 
  // restored to the proper geometry entity.
  int id_inc = 0;
  if( GSaveOpen::gso_sets_ids() )
  {
    id_inc = GSaveOpen::get_id_inc( attribOwnerEntity );
    if( id_inc == -1 )
      id_inc = 0;
  }
    
    // three possibilities:

    // 1) owner entity id is the same as this attrib id
  if (entityId == this_owner_id+id_inc) {
    hasActuated = CUBIT_TRUE;
    PRINT_DEBUG_102("CAEntityId::actuate: already have id for %s %d\n",
                    attribOwnerEntity->class_name(), this_owner_id+id_inc);
    return CUBIT_SUCCESS;
  }

  RefEntity *other_entity = 
    GeometryQueryTool::instance()->get_ref_entity(attribOwnerEntity->class_name(), 
                                                  entityId+id_inc);
    // 2) already an entity with the new id;
  if (other_entity) {
      // 2a) if other entity has a CAMP attribute, this owner has one too,
      //     and they both have the same unique id, these entities will
      //     get merged together; need to make sure we assign the lower id
      //     to the entity with the correct sense or direction
    DLIList<CubitAttrib*> att_list;
    other_entity->find_cubit_attrib_type(CA_MERGE_PARTNER, att_list);
    CAMergePartner *camp_ptr_other = NULL;
    CAMergePartner *camp_ptr_this = NULL;
    if (att_list.size() > 0) camp_ptr_other = CAST_TO(att_list.get(), CAMergePartner);
    att_list.clean_out();
    attribOwnerEntity->find_cubit_attrib_type(CA_MERGE_PARTNER, att_list);
    if (att_list.size() > 0) camp_ptr_this = CAST_TO(att_list.get(), CAMergePartner);
    
    if (camp_ptr_other && camp_ptr_this &&
        camp_ptr_other->merge_id()+id_inc == camp_ptr_this->merge_id()+id_inc) {
      PRINT_DEBUG_102("CAEntityId::actuate: another entity with same id & merge partner,"
                      " this = %s %d, other = %s %d\n",
                      attribOwnerEntity->class_name(), this_owner_id+id_inc, other_entity->class_name(),
                      entityId+id_inc);

      if (boundingUid != CUBIT_INT_MIN) {
          // check the sense of the *other* entity; can't check sense of this
          // entity yet, because the topology hasn't been filled yet
        RefEdge *edge = CAST_TO(other_entity, RefEdge);
        RefFace *face = CAST_TO(other_entity, RefFace);
        CubitBoolean switch_ids = CUBIT_FALSE;
        
        if (edge) {
          
          if (edge->start_vertex() != edge->end_vertex()) {
              // check the start vertex of the other entity
            ToolDataUser *tdu = TDUniqueId::find_td_unique_id(boundingUid, other_entity);
            RefVertex *vert = CAST_TO(tdu, RefVertex);
            assert(vert != 0);
            if (vert != edge->start_vertex()) 
                // other entity doesn't have the right sense, so we'll have to switch
                // (we don't know whether this one will have the right sense, and can't
                // check 'cuz this one's topology hasn't been completed yet)

                // only switch if other entity's id is lower, so that it's not the one kept
              if (entityId+id_inc > this_owner_id+id_inc) switch_ids = CUBIT_TRUE;
          }
          else if (boundingXYZ) {
              // else we have a single-vertex curve and an xyz value; 
              // check position of 1/3 parameter 
            const double ONE_THIRD = 1.0/3.0;
            CubitVector test_vec;
            /*CubitStatus result = */
            edge->position_from_fraction( ONE_THIRD, test_vec);
                                                               
            if (!GeometryQueryTool::instance()->about_spatially_equal(*boundingXYZ, test_vec)) {
              
                // only switch if other entity's id is lower, so that it's not the one kept
              if (entityId+id_inc > this_owner_id+id_inc) switch_ids = CUBIT_TRUE;
            }
          }
        }
        
        else if (face) {
            // check sense of other face wrt bounding uid entity; other face should be
            // fully constructed by now
          ToolDataUser *tdu = TDUniqueId::find_td_unique_id(boundingUid, other_entity);

          edge = CAST_TO(tdu, RefEdge);
          if (!edge) 
            PRINT_WARNING("CAEntityId::actuate: didn't find edge with uid %d for face %d.\n",
                          boundingUid, face->id());
          else if (edge->sense(face) != boundingSense && entityId+id_inc > this_owner_id+id_inc) 
            switch_ids = CUBIT_TRUE;
        }

        if (switch_ids) {
          attribOwnerEntity->set_id(0, CUBIT_FALSE);
          other_entity->set_id(this_owner_id+id_inc);
          PRINT_DEBUG_102("CAEntityId::actuate: other entity with same id, merge partner, and correct sense;"
                          " switching id's, this = %s %d, other = %s %d\n",
                          attribOwnerEntity->class_name(), this_owner_id+id_inc, other_entity->class_name(),
                          entityId+id_inc);
          attribOwnerEntity->set_id (entityId+id_inc, CUBIT_FALSE);
          attribOwnerEntity->color(CUBIT_DEFAULT_COLOR);
        }
      }
      
      hasActuated = CUBIT_TRUE;
      return CUBIT_SUCCESS;
    }
    
      // 2b) if other entity has a CAEID attribute, check to make sure it's not
      //     the same, and if not, switch real ids with it (it will get changed later);
      //     otherwise, print a warning
    other_entity->find_cubit_attrib_type(CA_ENTITY_ID, att_list);
    CAEntityId *other_caeid = (att_list.size() ?  CAST_TO(att_list.get(), CAEntityId)
                               : NULL);
    if ( other_caeid && other_caeid->id()+id_inc != entityId+id_inc) {
        // need to reset owner entity id first, so that we don't have
        // two identical ids active at the same time (messes up the
        // graphics)
      
      attribOwnerEntity->set_id(0, CUBIT_FALSE);
      other_entity->set_id(this_owner_id+id_inc);
      PRINT_DEBUG_102("CAEntityId::actuate: other entity with same id, NO merge partner;"
                      " switching id's, this = %s %d, other = %s %d\n",
                      attribOwnerEntity->class_name(), this_owner_id+id_inc, other_entity->class_name(),
                      entityId+id_inc);
    }
  
    else if( camp_ptr_other )
    {
      DLIList<RefEntity*> merge_partners;
      camp_ptr_other->merge_prepare( merge_partners );
      int lowest_id = entityId + id_inc;
      for( int i = merge_partners.size(); i--; )
        if( merge_partners.step_and_get()->id() < lowest_id )
          lowest_id = merge_partners.get()->id();
      
      if( lowest_id < entityId+id_inc )
        // Entity will merge later and loose its current Id,
        // so we can take the Id for this entity.
      {
        attribOwnerEntity->set_id(0, CUBIT_FALSE);
        other_entity->set_id( this_owner_id + id_inc );
      }
    
    }
  }

      // 2b) else print an error and don't actuate
  if( other_entity && other_entity->id() == (entityId+id_inc) )
  {
//    PRINT_WARNING ( "Duplicate entity id attribute on %s %d; duplicated id = %d.\n"
//                    "This sometimes happens when id attributes are exported without exporting"
//                    " merge attribute.\n",
//                    attribOwnerEntity->class_name(), this_owner_id+id_inc, entityId+id_inc);
//    hasActuated = CUBIT_TRUE;
    CADeferredAttrib::add_unactuated_ca(this);
    return CUBIT_FAILURE;
  }
  
  else PRINT_DEBUG_102("CAEntityId::actuate: setting id of %s %d; new id = %d\n",
                       attribOwnerEntity->class_name(), this_owner_id+id_inc, entityId+id_inc);
  
    // ok, now set the id and return
  attribOwnerEntity->set_id (entityId+id_inc, CUBIT_FALSE);
  attribOwnerEntity->color(CUBIT_DEFAULT_COLOR);

  hasActuated = CUBIT_TRUE;

  return CUBIT_SUCCESS;
}

CubitStatus CAEntityId::update()
{
  if ( hasUpdated ) 
     return CUBIT_SUCCESS;

  if ( !attribOwnerEntity)
    return CUBIT_FAILURE;

  entityId = attribOwnerEntity->id ( );

    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // set the sense data
  TopologyEntity *topo_entity = CAST_TO(attribOwnerEntity, TopologyEntity);
  assert(topo_entity != 0);
  if (MergeTool::instance()->entity_merged(topo_entity) == CUBIT_FALSE)
    return CUBIT_SUCCESS;
  
    // get the uid of the first child entity
  RefEdge *edge;
  RefFace *face;
  if ((edge = CAST_TO(attribOwnerEntity, RefEdge)) != NULL) {
    RefVertex *vert = edge->start_vertex();
    boundingUid = TDUniqueId::get_unique_id(vert);
    boundingSense = CUBIT_FORWARD;

      // need to check for single-vertex curves
    if (vert == edge->end_vertex()) {
        // if it is one, store the xyz of 1/3 down the curve
      const double ONE_THIRD = 1.0/3.0;
        // Find the point 1/3 along *this*
      if (!boundingXYZ) boundingXYZ = new CubitVector();
      CubitStatus result = edge->position_from_fraction( ONE_THIRD,
                                                         *boundingXYZ);
      if ( result != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Error in CAEntityId::update;"
                    "Can't find position 1/3 along curve.\n");
        return CUBIT_FAILURE;
      }
    }
  }

  else if ((face = CAST_TO(attribOwnerEntity, RefFace)) != NULL) {
    DLIList<RefEdge*> edges;
    face->ref_edges(edges);
    if (edges.size()) {
      edge = edges.get_and_step();
      boundingUid = TDUniqueId::get_unique_id(edge);
      boundingSense = edge->sense(face); 
      return CUBIT_SUCCESS;
    }
  }
    
  return CUBIT_SUCCESS;
}

void CAEntityId::merge_owner(CubitAttrib *deletable_attrib)
{
    // take the id with the lowest value
  CAEntityId *other_caeid = CAST_TO(deletable_attrib, CAEntityId);
  
  if (other_caeid &&
      (entityId == 0 || entityId > other_caeid->id()))
    entityId = other_caeid->id();
}

CubitSimpleAttrib* CAEntityId::cubit_simple_attrib()
{
  DLIList<CubitString*> cs_list;
  DLIList<double> d_list;
  DLIList<int> i_list;

  i_list.append ( entityId );
  i_list.append ( boundingUid );
  i_list.append ( boundingSense );
  if (boundingXYZ) {
    d_list.append ( boundingXYZ->x() );
    d_list.append ( boundingXYZ->y() );
    d_list.append ( boundingXYZ->z() );
  }
    
  cs_list.append(new CubitString(att_internal_name()));

  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list,
                                                          &d_list,
                                                          &i_list);
  int i;
  for ( i = cs_list.size(); i--;) delete cs_list.get_and_step();
  
  return csattrib_ptr;
}

void CAEntityId::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAEntityId: owner = %s %d:   id=%d, color=%d\n",
             attribOwnerEntity->class_name(), attribOwnerEntity->id(),
             entityId, 0);
}
