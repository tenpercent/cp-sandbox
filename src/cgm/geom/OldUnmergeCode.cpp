#include "OldUnmergeCode.hpp"

#include "Body.hpp"
#include "CoVolume.hpp"
#include "RefVolume.hpp"
#include "Shell.hpp"
#include "CoFace.hpp"
#include "RefFace.hpp"
#include "Loop.hpp"
#include "CoEdge.hpp"
#include "RefEdge.hpp"
#include "Chain.hpp"
#include "CoVertex.hpp"
#include "RefVertex.hpp"

#include "BodySM.hpp"
#include "Lump.hpp"
#include "ShellSM.hpp"
#include "Surface.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"
#include "Curve.hpp"
#include "Point.hpp"

#include "GeometryQueryTool.hpp"
#include "RefEntityFactory.hpp"
#include "CpuTimer.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"

#include "MergeTool.hpp"
#include "MergeToolAssistant.hpp"
#include "SettingHandler.hpp"


void OldUnmergeCode::initialize_settings()
{
  SettingHandler::instance()->add_setting("unmerge new ids",
                                          OldUnmergeCode::set_use_old_unmerge_code,
                                          OldUnmergeCode::get_use_old_unmerge_code);
}

bool OldUnmergeCode::useOldUnmergeCode = false;

   
bool OldUnmergeCode::get_use_old_unmerge_code()
    { return useOldUnmergeCode; }
    
void OldUnmergeCode::set_use_old_unmerge_code( bool value )
    { useOldUnmergeCode = value; }




OldUnmergeCode& OldUnmergeCode::instance()
{
  static OldUnmergeCode instance_;
  return instance_;
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge everything
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/27/01
//-------------------------------------------------------------------------
CubitStatus OldUnmergeCode::unmerge_all()
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge_all();
  
  int i;
  CubitStatus result = CUBIT_SUCCESS;
  CubitBoolean top = start_unmerge();
  
  
  DLIList<RefFace*>  face_list;
  DLIList<RefEdge*>  edge_list;
  DLIList<RefVertex*> vtx_list;
  
  GeometryQueryTool::instance()->ref_faces( face_list );
  GeometryQueryTool::instance()->ref_edges( edge_list );
  GeometryQueryTool::instance()->ref_vertices( vtx_list );
  
  for( i = face_list.size(); (i > 0) && !CubitMessage::instance()->Interrupt(); i-- )
    if( ! unmerge(face_list.get_and_step(),CUBIT_FALSE) )
      result = CUBIT_FAILURE;
 
  for( i = edge_list.size(); (i > 0) && !CubitMessage::instance()->Interrupt(); i-- )
    if( ! unmerge(edge_list.get_and_step(),CUBIT_FALSE) ) 
      result = CUBIT_FAILURE;
  
  for( i = vtx_list.size(); (i > 0) && !CubitMessage::instance()->Interrupt(); i-- )
    if( ! unmerge(vtx_list.get_and_step()) )
      result = CUBIT_FAILURE;
  
  end_unmerge(top);
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Unmerge RefEntities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
CubitStatus OldUnmergeCode::unmerge( DLIList<RefEntity*> &entity_list,
                                CubitBoolean descend )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(entity_list, descend);
  
  CubitBoolean top = start_unmerge();
  
  for( int i = entity_list.size(); (i > 0) && !CubitMessage::instance()->Interrupt(); i-- )
    unmerge( entity_list.get_and_step(), descend );
  
  end_unmerge(top);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge a RefEntity
//
// Special Notes : All parents must be unmerged.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
CubitStatus OldUnmergeCode::unmerge( RefEntity* entity_ptr, CubitBoolean descend )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(entity_ptr, descend);
  
  if( CAST_TO( entity_ptr, Body ) )
     return descend ? unmerge(CAST_TO(entity_ptr,Body)) : CUBIT_FAILURE;
  else if( CAST_TO( entity_ptr, RefVolume ) )
     return descend ? unmerge(CAST_TO(entity_ptr,RefVolume)) : CUBIT_FAILURE;
  else if( CAST_TO( entity_ptr, RefFace ) )
     return unmerge( CAST_TO(entity_ptr,RefFace), descend );
  else if( CAST_TO( entity_ptr, RefEdge ) )
     return unmerge( CAST_TO(entity_ptr,RefEdge), descend );
  else if( CAST_TO( entity_ptr, RefVertex ) )
     return unmerge( CAST_TO(entity_ptr,RefVertex) );
  else
  {
    PRINT_ERROR("Bad Entity \"%s\" in "
                "OldUnmergeCode::unmerge(RefEntity*,CubitBoolean)\n",
                entity_ptr->class_name());
      return CUBIT_FAILURE;
  }
}

CubitStatus OldUnmergeCode::unmerge( Body* body_ptr )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(body_ptr);
  
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  DLIList<RefVolume*> vol_list;
  body_ptr->ref_volumes(vol_list);
  for( int i = vol_list.size(); (i > 0) && !CubitMessage::instance()->Interrupt(); i-- )
    if( !unmerge(vol_list.get_and_step()) )
      result = CUBIT_FAILURE;
  end_unmerge(top);
  return result;
}

CubitStatus OldUnmergeCode::unmerge( RefVolume* vol_ptr )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(vol_ptr);
  
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  DLIList<RefFace*> face_list, new_faces;
  DLIList<RefEdge*> edge_list, new_edges;
#ifdef BOYD17 
  DLIList<RefVertex*> vtx_list, new_vertices;
#endif
  DLIList<RefVertex*> vtx_list;
  vol_ptr->ref_faces(face_list);
  vol_ptr->ref_edges(edge_list);
  vol_ptr->ref_vertices(vtx_list);
  face_list.reset();
  while( face_list.size() && !CubitMessage::instance()->Interrupt() )
  {
    RefFace* old_face = face_list.extract();
    RefFace* new_face = unmerge(old_face,vol_ptr);
    if( new_face ) 
      new_faces.append(new_face);
    else
      result = CUBIT_FAILURE;
  }
  
  while( edge_list.size() && !CubitMessage::instance()->Interrupt() )
  {
    RefEdge* old_edge = edge_list.extract();
    face_list.clean_out();
    old_edge->ref_faces(face_list);
    while( face_list.size() )
    {
      RefFace* face_ptr = face_list.extract();
      if( new_faces.is_in_list(face_ptr) )
      {
        RefEdge* new_edge = unmerge(old_edge,face_ptr);
        if( new_edge )
          new_edges.append(new_edge);
        else
          result = CUBIT_FAILURE;
        break;
      }
    }
  }
  
  
  while( vtx_list.size() && !CubitMessage::instance()->Interrupt() )
  {
    RefVertex* old_vtx = vtx_list.extract();
    edge_list.clean_out();
    old_vtx->ref_edges(edge_list);
    while( edge_list.size() )
    {
      RefEdge* edge_ptr = edge_list.extract();
      if( new_edges.is_in_list(edge_ptr) )
      {
        RefVertex* new_vtx = unmerge(old_vtx,edge_ptr);
        if( !new_vtx )
          result = CUBIT_FAILURE;
        break;
      }
    }
  }
  
  end_unmerge(top);
  return result;
}

CubitStatus OldUnmergeCode::unmerge( RefFace* face_ptr, CubitBoolean descend )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(face_ptr);
  
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  CubitBoolean reversed;
  int i;
  
//  if( !face_ptr->can_modify() )
//    return CUBIT_FAILURE;
  
  DLIList<RefFace*> unmerged_faces;
  DLIList<RefVolume*> vol_list;
  face_ptr->ref_volumes(vol_list);
  if( vol_list.size() )
  {
    for( i = vol_list.size(); i > 0; i-- )
    {
      RefFace* new_face = unmerge(face_ptr,vol_list.get_and_step());
      if( new_face )
        unmerged_faces.append_unique(new_face);
      else
        result = CUBIT_FAILURE;
    }
  }
  else
  {
    DLIList<TopologyBridge*> bridge_list;
    face_ptr->bridge_manager()->get_bridge_list( bridge_list );
    
    //Try top remove each bridge
    for( i = bridge_list.size(); i > 0; i-- )
    {
      //but stop if there is only one left
      if( face_ptr->bridge_manager()->number_of_bridges() == 1 )
      {
        unmerged_faces.append(face_ptr);
        break;
      }
      
      TopologyBridge* bridge_ptr = bridge_list.get_and_step();
      Surface* surf_ptr = CAST_TO(bridge_ptr,Surface);
      assert(surf_ptr != 0);
      RefFace* new_face = split_out_Surface(surf_ptr, reversed);
      if( new_face )
      {
        unmerged_faces.append(new_face);
  
          //Notify merge assistants of unmerge
        DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
        for( int a = assistant_list_.size(); a > 0; a-- )
          assistant_list_.get_and_step()->
            unmerged( face_ptr, new_face, reversed );
      }
      else
        result = CUBIT_FAILURE;
    }
    if( !unmerged_faces.is_in_list(face_ptr) )
      result = CUBIT_FAILURE;
  }
  
  if( !descend ) 
  {
    end_unmerge(top);
    return result;
  }
  
  DLIList<RefEdge*> edge_list, unmerged_edges;
  for( i = unmerged_faces.size(); i > 0; i-- )
  {
    face_ptr = unmerged_faces.get_and_step();
    edge_list.clean_out();
    face_ptr->ref_edges(edge_list);
    for( int j = edge_list.size(); j > 0; j-- )
    {
      RefEdge* new_edge = unmerge(edge_list.get_and_step(),face_ptr);
      if( new_edge )
        unmerged_edges.append_unique(new_edge);
      else
        result = CUBIT_FAILURE;
    }
  }
  
  DLIList<RefVertex*> vtx_list;
  for( i = unmerged_edges.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = unmerged_edges.get_and_step();
    vtx_list.clean_out();
    edge_ptr->ref_vertices(vtx_list);
    for( int j = vtx_list.size(); j > 0; j-- )
    {
      RefVertex* new_vtx = unmerge(vtx_list.get_and_step(),edge_ptr);
      if( !new_vtx ) result = CUBIT_FAILURE;
    }
  }
       
  end_unmerge(top);
  return result;
}


CubitStatus OldUnmergeCode::unmerge( RefEdge* edge_ptr, CubitBoolean descend )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(edge_ptr);
  
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  CubitBoolean reversed;
  int i;
  
//  if( !edge_ptr->can_modify() )
//    return CUBIT_FAILURE;
  
  DLIList<RefEdge*> unmerged_edges;
  DLIList<RefFace*> face_list;
  edge_ptr->ref_faces(face_list);
  if( face_list.size() )
  {
    for( i = face_list.size(); i > 0; i-- )
    {
      RefEdge* new_edge = unmerge(edge_ptr,face_list.get_and_step());
      if( new_edge )
        unmerged_edges.append_unique(new_edge);
      else
        result = CUBIT_FAILURE;
    }
  }
  else
  {
    DLIList<TopologyBridge*> bridge_list;
    DLIList<Curve*> curve_list;
    edge_ptr->bridge_manager()->get_bridge_list( bridge_list );
    
    //Try top remove each bridge
    for( i = bridge_list.size(); i > 0; i-- )
    {
      //but stop if there is only one left
      if( edge_ptr->bridge_manager()->number_of_bridges() == 1 )
      {
        unmerged_edges.append(edge_ptr);
        break;
      }
      
      TopologyBridge* bridge_ptr = bridge_list.get_and_step();
      curve_list.clean_out();
      curve_list.append( CAST_TO(bridge_ptr,Curve) );
      RefEdge* new_edge = split_out_Curves(curve_list, reversed);
      if( new_edge )
      {
        unmerged_edges.append(new_edge);
  
          //Notify merge assistants of unmerge
        DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
        for( int a = assistant_list_.size(); a > 0; a-- )
          assistant_list_.get_and_step()->
            unmerged( edge_ptr, new_edge, reversed );
      }
      else
        result = CUBIT_FAILURE;
    }
    if( !unmerged_edges.is_in_list(edge_ptr) )
      result = CUBIT_FAILURE;
  }
  
  if( !descend ) 
  {
    end_unmerge(top);
    return result;
  }
  
  DLIList<RefVertex*> vtx_list;
  for( i = unmerged_edges.size(); i > 0; i-- )
  {
    edge_ptr = unmerged_edges.get_and_step();
    vtx_list.clean_out();
    edge_ptr->ref_vertices(vtx_list);
    for( int j = vtx_list.size(); j > 0; j-- )
    {
      RefVertex* new_vtx = unmerge(vtx_list.get_and_step(),edge_ptr);
      if( !new_vtx )
        result = CUBIT_FAILURE;
    }
  }
  
  end_unmerge(top);
  return result;
}


CubitStatus OldUnmergeCode::unmerge( RefVertex* vtx_ptr )
{
  if (!get_use_old_unmerge_code())
    return MergeTool::instance()->unmerge(vtx_ptr);
  
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  int i;
  
//  if( !vtx_ptr->can_modify() )
//    return CUBIT_FAILURE;
  
  DLIList<RefVertex*> unmerged_vtxs;
  DLIList<RefEdge*> edge_list;
  vtx_ptr->ref_edges(edge_list);
  if( edge_list.size() )
  {
    for( i = edge_list.size(); i > 0; i-- )
    {
      RefVertex* new_vtx = unmerge(vtx_ptr,edge_list.get_and_step());
      if( new_vtx )
        unmerged_vtxs.append_unique(new_vtx);
      else
        result = CUBIT_FAILURE;
    }
  }
  else
  {
    DLIList<TopologyBridge*> bridge_list;
    DLIList<Point*> point_list;
    vtx_ptr->bridge_manager()->get_bridge_list( bridge_list );
    
    //Try top remove each bridge
    for( i = bridge_list.size(); i > 0; i-- )
    {
      //but stop if there is only one left
      if( vtx_ptr->bridge_manager()->number_of_bridges() == 1 )
      {
        unmerged_vtxs.append(vtx_ptr);
        break;
      }
      
      TopologyBridge* bridge_ptr = bridge_list.get_and_step();
      point_list.clean_out();
      point_list.append( CAST_TO(bridge_ptr,Point) );
      RefVertex* new_vtx = split_out_Points(point_list);
      if( new_vtx )
      {
        unmerged_vtxs.append(new_vtx);
  
          //Notify merge assistants of unmerge
        DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
        for( int a = assistant_list_.size(); a > 0; a-- )
          assistant_list_.get_and_step()->
            unmerged( vtx_ptr, new_vtx, CUBIT_FALSE );
      }
      else
        result = CUBIT_FAILURE;
    }
    if( !unmerged_vtxs.is_in_list(vtx_ptr) )
      result = CUBIT_FAILURE;
  }
  
  end_unmerge(top);
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Given an unmerged parent, and a merged child, 
//                 unmerge the topology bridge of the child that 
//                 corresponds to the passed parent.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
RefFace* OldUnmergeCode::unmerge( RefFace* face_ptr, RefVolume* vol_ptr )
{  
  assert( face_ptr && vol_ptr );
  CubitBoolean top = start_unmerge();
  CubitBoolean reversed;
  int i;
  
//  if( !face_ptr->can_modify() )
//  {
//    end_unmerge(top);
//    return 0;
//  }
  
  //If passed entity is not merged, just return it.
  if( face_ptr->bridge_manager()->number_of_bridges() < 2 )
  {
    end_unmerge(top);
    return face_ptr;
  }
  
  //Find the Surfaces of face_ptr that go with vol_ptr
  Lump* lump_ptr = CAST_TO( vol_ptr->get_geometry_entity_ptr(), Lump );
  assert(lump_ptr != 0);
  DLIList<Lump*> surf_lumps;
  DLIList<Surface*> surface_list;
  DLIList<TopologyBridge*> bridge_list;
  face_ptr->bridge_manager()->get_bridge_list(bridge_list);
  for( i = bridge_list.size(); i > 0; i-- )
  {
    TopologyBridge* bridge_ptr = bridge_list.get_and_step();
    surf_lumps.clean_out();
    bridge_ptr->lumps(surf_lumps);
    for( int j = surf_lumps.size(); j> 0; j-- )
    {
      if( surf_lumps.get_and_step()->topology_entity() == vol_ptr )
      {
        surface_list.append(CAST_TO(bridge_ptr,Surface));
        break;
      }
    }
  }
  //face_ptr and vol_ptr have no association in SM topology!!
  if( ! surface_list.size() )
  {
    end_unmerge(top);
    return 0;
  }
  
  //Assuming this is true, until this assert fails
  assert(surface_list.size() == 1);
  Surface* surface_ptr = surface_list.get();

  // Unmerge any child GroupingEntities and SenseEntities,
  // and reconstruct a new owner.
  RefFace* new_entity 
    = split_out_Surface( surface_ptr, reversed );
  if( !new_entity )
  {
    end_unmerge(top);
    return 0;
  }

  // Find any links from parent to old entity, and move them
  // to the new, unmerged entity.
  DLIList<CoFace*> coface_list;
  DLIList<ShellSM*> shellsm_list;
  surface_ptr->shellsms( shellsm_list );
  for( i = shellsm_list.size(); i > 0; i-- )
  {
    ShellSM* shellsm = shellsm_list.get_and_step();
    Shell* shell_ptr = CAST_TO(shellsm->topology_entity(),Shell);
    if( !shell_ptr ) continue;
    
    coface_list.clean_out();
    shell_ptr->co_faces(coface_list);
    for( int j = coface_list.size(); j > 0; j-- )
    {
      CoFace* cof_ptr = coface_list.get_and_step();
      if( cof_ptr->get_ref_face_ptr() == face_ptr )
      {
        cof_ptr->switch_basic_topology_entity(new_entity);
        if( reversed ) cof_ptr->reverse_sense();
      }
    }
  }
  
  //We have changed the topology of the passed RefVolume.
  //Note this for later.
  unmerge_modified.append( face_ptr );
  
  //Notify merge assistants of unmerge
  DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
  for( int a = assistant_list_.size(); a > 0; a-- )
    assistant_list_.get_and_step()->
      unmerged( face_ptr, new_entity, reversed );
  
  end_unmerge(top);  
  return new_entity;
}


RefEdge* OldUnmergeCode::unmerge( RefEdge* edge_ptr, RefFace* face_ptr )
{  
  assert( face_ptr && edge_ptr );
  CubitBoolean top = start_unmerge();
  CubitBoolean reversed;
  int i;
  
//  if( !edge_ptr->can_modify() )
//  {
//    end_unmerge(top);
//    return 0;
//  }
  
  //If passed entity is not merged, just return it.
  if( edge_ptr->bridge_manager()->number_of_bridges() < 2 )
  {
    end_unmerge(top);
    return edge_ptr;
  }
  
  //If parent is merged, cannot proceed.
  if( face_ptr->bridge_manager()->number_of_bridges() != 1 )
  {
    end_unmerge(top);
    return 0;
  }
  
  //Parent is not merged.  This is the only Surface
  Surface* surf_ptr = CAST_TO( face_ptr->get_geometry_entity_ptr(), Surface );
  assert(surf_ptr != 0);
  
  //Find curves in edge_ptr associated with surf_ptr.
  //We need to check links in both directions becuase with virtual
  //geometry, links may only exist in one direction.
  //
  //Check downward links.
  DLIList<Curve*> curve_list, surf_curves;
  surf_ptr->curves(surf_curves);
  for( i = surf_curves.size(); i > 0; i-- )
  {
    Curve* curve_ptr = surf_curves.get_and_step();
    if( curve_ptr->topology_entity() == edge_ptr )
      curve_list.append(curve_ptr);
//    else if( curve_ptr->topology_entity()->is_parasite(edge_ptr) )
//      curve_list.append(curve_ptr);
  }
  //Check upward links.
  DLIList<TopologyBridge*> bridge_list;
  DLIList<Surface*> curve_surfs;
  edge_ptr->bridge_manager()->get_bridge_list(bridge_list);
  for( i = bridge_list.size(); i > 0; i-- )
  {
    TopologyBridge* bridge_ptr = bridge_list.get_and_step();
    curve_surfs.clean_out();
    bridge_ptr->surfaces(curve_surfs);
    for( int j = curve_surfs.size(); j > 0; j-- )
    {
      TopologyEntity* te_ptr = curve_surfs.get_and_step()->topology_entity();
      //if( te_ptr == face_ptr )
      if( (te_ptr == face_ptr) /*|| face_ptr->is_host(te_ptr)*/ )
      {
        curve_list.append_unique(CAST_TO(bridge_ptr,Curve));
        break;
      }
    }
  }
  
  //If there is no association between the passed edge_ptr and
  //face_ptr in the SM topology
  if( ! curve_list.size() )
  {
    end_unmerge(top);
    return 0;
  }
  
  // Unmerge any child GroupingEntities and SenseEntities,
  // and reconstruct a new owner.
  RefEdge* new_entity 
    = split_out_Curves( curve_list, reversed );
  if( !new_entity )
  {
    end_unmerge(top);
    return 0;
  }
  
  //We have changed the topology of the passed RefFace.
  //Note this for later.
  unmerge_modified.append( edge_ptr );

  // Find any links from parent to old entity, and move them
  // to the new, unmerged entity.
  DLIList<CoEdge*> coedge_list;
  DLIList<Curve*> coe_curves, tmp_list;
  DLIList<TopologyBridge*> tb_list;
  edge_ptr->get_co_edges (coedge_list);
  for( i = coedge_list.size(); i > 0; i-- )
  {
    tb_list.clean_out();
    coe_curves.clean_out();
    CoEdge* coe_ptr = coedge_list.get_and_step();
    coe_ptr->bridge_manager()->get_bridge_list(tb_list);
    for( int j = tb_list.size(); j> 0; j-- )
    {
      tmp_list.clean_out();
      tb_list.get_and_step()->curves(tmp_list);
      coe_curves.merge_unique(tmp_list);
    }
    coe_curves.intersect(curve_list);
    
    if( coe_curves.size() )
    {
      coe_ptr->switch_basic_topology_entity(new_entity);
      if( reversed ) 
        coe_ptr->reverse_sense();
    }
  }
  
  //Notify merge assistants of unmerge
  DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
  for( int a = assistant_list_.size(); a > 0; a-- )
    assistant_list_.get_and_step()->
      unmerged( edge_ptr, new_entity, reversed );
  
  
  end_unmerge(top);  
  return new_entity;
}

RefVertex* OldUnmergeCode::unmerge( RefVertex* vtx_ptr, RefEdge* edge_ptr )
{  
  assert( vtx_ptr && edge_ptr );
  CubitBoolean top = start_unmerge();
  int i;
  
//  if( !vtx_ptr->can_modify() )
//  {
//    end_unmerge(top);
//    return 0;
//  }
  
  //If passed entity is not merged, just return it.
  if( vtx_ptr->bridge_manager()->number_of_bridges() < 2 )
  {
    end_unmerge(top);
    return vtx_ptr;
  }
  
  //If parent is merged, cannot proceed.
  if( edge_ptr->bridge_manager()->number_of_bridges() != 1 )
  {
    end_unmerge(top);
    return 0;
  }
  
  //Parent is not merged,  This is the only Curve.
  Curve* curve_ptr = CAST_TO( edge_ptr->get_geometry_entity_ptr(), Curve );
  assert(curve_ptr != 0);

  //Find points in vtx_ptr associated with curve_ptr.
  DLIList<Point*> curve_points, point_list;
  curve_ptr->points(curve_points);
  for( i = curve_points.size(); i > 0; i-- )
  {
    Point* point_ptr = curve_points.get_and_step();
    if( point_ptr->topology_entity() == vtx_ptr )
      point_list.append(point_ptr);
  }

  //If parent and child are not related in the SolidModelingEngine,
  //return failure.
  if( ! point_list.size() )
  {
    end_unmerge(top);
    return 0;
  }

  // Find any links from parent to old entity, and move them
  // to the new, unmerged entity.
  DLIList<RefEdge*> edge_list; 
  DLIList<Point*> curve_pts;
  DLIList<CoVertex*> cvtx_to_change;
  
  vtx_ptr->ref_edges( edge_list );
  for( i = edge_list.size(); i > 0; i-- )
  {
    edge_ptr = edge_list.get_and_step();
    if( edge_ptr->bridge_manager()->number_of_bridges() > 1 )
      continue;
    
    curve_ptr = edge_ptr->get_curve_ptr();
    assert( curve_ptr != 0 );
    curve_pts.clean_out();
    curve_ptr->points(curve_pts);
    curve_pts.intersect(point_list);
    if( ! curve_pts.size() )
      continue;
      
    CoVertex* start_cvtx = edge_ptr->get_chain_ptr()->start_co_vertex();
    CoVertex*   end_cvtx = edge_ptr->get_chain_ptr()->  end_co_vertex();
    if (start_cvtx && start_cvtx->get_ref_vertex_ptr() == vtx_ptr)
      cvtx_to_change.append( start_cvtx );
    if (end_cvtx && end_cvtx->get_ref_vertex_ptr() == vtx_ptr)
      cvtx_to_change.append( end_cvtx );
      
  }

  // Make new RefVertex with points.
  RefVertex* new_entity = split_out_Points( point_list );
  if( !new_entity )
  {
    end_unmerge(top);
    return 0;
  }

  for( i = cvtx_to_change.size(); i--; )
  {
    cvtx_to_change.get_and_step()->switch_basic_topology_entity(new_entity);
  }
  
  //We have changed the topology of the passed RefEdge.
  //Note this for later.
  unmerge_modified.append( vtx_ptr );
  
  //Notify merge assistants of unmerge
  DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
  for( int a = assistant_list_.size(); a > 0; a-- )
    assistant_list_.get_and_step()->
      unmerged( vtx_ptr, new_entity, CUBIT_FALSE );
  
  
  end_unmerge(top);  
  return new_entity;
}


void OldUnmergeCode::remove_CAEntityId_attrib( TopologyBridge* tb_ptr )
{
  DLIList<CubitSimpleAttrib*> attrib_list;
  tb_ptr->get_simple_attribute( attrib_list );
  for( int i = attrib_list.size(); i > 0; i-- )
  {
    CubitSimpleAttrib* attrib = attrib_list.get_and_step();
    if( attrib->character_type() == "ENTITY_ID" )
    {
      tb_ptr->remove_simple_attribute_virt( attrib );
    }
  }
}  


//-------------------------------------------------------------------------
// Purpose       : Remove a surface from its owning RefFace, and make a new
//                 RefFace with the Surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
RefFace* OldUnmergeCode::split_out_Surface( Surface* surface, CubitBoolean& reversed )
{
  // Remove the child from the old owner.
  reversed = false;
  TopologyEntity* te_ptr = surface->topology_entity();
  RefFace* old_entity = CAST_TO( te_ptr, RefFace );
  if( te_ptr )
  {
    assert( old_entity != 0 );
    
    te_ptr->bridge_manager()->remove_bridge( surface );
    remove_CAEntityId_attrib( surface );
    old_unmerged.append( old_entity );
    
    if (te_ptr->bridge_manager()->topology_bridge()->bridge_sense() == CUBIT_REVERSED)
    { 
      te_ptr->reverse_topology();
      te_ptr->bridge_manager()->reverse_bridge_senses();
      reversed = !reversed;
    }
    if (surface->bridge_sense() == CUBIT_REVERSED)
    {
      surface->reverse_bridge_sense();
      reversed = !reversed;
    }
  }
  
  surface->set_saved_id(0);
  
  // make new refface
  RefFace* new_entity = RefEntityFactory::instance()->construct_RefFace( surface );

  // Unmerge Loops
  DLIList<LoopSM*> loopsms;
  surface->loopsms( loopsms );
  for( int i = loopsms.size(); i > 0; i-- )
  {
    Loop* new_loop = split_out_Loop(loopsms.get_and_step(), new_entity, reversed);
    new_entity->add_grouping_entity( new_loop );
  }
  
  new_unmerged.append(new_entity);
  event_list.append( new UnMergeEvent( old_entity, new_entity ) );
  
  return new_entity;
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge Loop and CoEdges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/23/01
//-------------------------------------------------------------------------
Loop* OldUnmergeCode::split_out_Loop( LoopSM* loopsm, RefFace* /*new_entity*/, CubitBoolean reverse )
{
  // Remove from existing TopologyEntity
  TopologyEntity* topo = loopsm->topology_entity();
  loopsm->bridge_manager()->remove_bridge( loopsm );

  int i;
  Loop* old_loop = CAST_TO( topo, Loop );  
  assert( old_loop != 0 );
  //RefFace* old_face = old_loop->get_ref_face_ptr();
  Loop* new_loop = new Loop( loopsm );
  DLIList<TopologyBridge*> tb_list;
  DLIList<CoEdge*> coedges;
#ifdef BOYD17 
  DLIList<Surface*> coesm_surfs;
#endif
  old_loop->ordered_co_edges( coedges );
  if( !reverse ) coedges.reverse();
  
  CoEdge* prev = 0;
  while( coedges.size() > 0 )
  {
    coedges.last();
    CoEdge* coedge_ptr = coedges.remove();
    
    // Find the corresponding CoEdgeSM to unmerge.
    tb_list.clean_out();
    coedge_ptr->bridge_manager()->get_bridge_list(tb_list);
    CoEdgeSM* coedgesm = NULL;  //Match to LoopSM
    for( i = tb_list.size(); i > 0; i-- )
    {
      TopologyBridge* tb_ptr = tb_list.get_and_step();
      if( tb_ptr->loopsm() == loopsm )
      {
        coedgesm = CAST_TO(tb_ptr,CoEdgeSM);
        break;
      }
    }
    
    assert(coedgesm != 0);
    
    // Unmerge CoEdge
    coedge_ptr->bridge_manager()->remove_bridge( coedgesm );
    RefEdge* edge_ptr = coedge_ptr->get_ref_edge_ptr();
    CubitSense sense = coedge_ptr->get_sense();
    if( reverse ) 
      sense = CubitUtil::opposite_sense( sense );
    CoEdge* new_coedge = new CoEdge( edge_ptr, sense );
    new_coedge->set_co_edge_sm_ptr( coedgesm );
    new_loop->add_sense_entity( new_coedge, prev );
    prev = new_coedge;
  }
  
  return new_loop;
}


    
//-------------------------------------------------------------------------
// Purpose       : VG provides downward query of TopologyBridges, but
//                 not always upward queries.  So to do upward queries,
//                 do downward queries and search for source object.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/26/01
//-------------------------------------------------------------------------
void OldUnmergeCode::find_curves( Point* point_ptr, DLIList<Curve*>& result_set )
{
  //Do normal TB query for real geometry
  point_ptr->curves(result_set);
  
  //Now account for missing TB links due to virtual geometry
  RefVertex* vtx_ptr = CAST_TO(point_ptr->topology_entity(),RefVertex);
  if( vtx_ptr )
  {
    DLIList<RefEdge*> vtx_edges;
    DLIList<Point*> curve_pts;
    DLIList<TopologyBridge*> edge_bridges;
    vtx_ptr->ref_edges(vtx_edges);
    for( int i = vtx_edges.size(); i > 0; i-- )
    {
      RefEdge* edge_ptr = vtx_edges.get_and_step();
      edge_bridges.clean_out();
      edge_ptr->bridge_manager()->get_bridge_list( edge_bridges );
      for( int j = edge_bridges.size(); j > 0; j-- )
      {
        TopologyBridge* bridge_ptr = edge_bridges.get_and_step();
        curve_pts.clean_out();
        bridge_ptr->points(curve_pts);
        if( curve_pts.is_in_list( point_ptr ) )
        {
          Curve* curve_ptr = CAST_TO(bridge_ptr,Curve);
          assert(curve_ptr != 0);
          result_set.append_unique(curve_ptr);
        }
      }
    }
  }
}
    
//-------------------------------------------------------------------------
// Purpose       : VG provides downward query of TopologyBridges, but
//                 not always upward queries.  So to do upward queries,
//                 do downward queries and search for source object.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/26/01
//-------------------------------------------------------------------------
void OldUnmergeCode::find_surfaces( Curve* curve_ptr, DLIList<Surface*>& result_set )
{
  //Do normal TB query for real geometry
  curve_ptr->surfaces(result_set);
}


//-------------------------------------------------------------------------
// Purpose       : Unmerge a Curve from its owning RefEdge, and make a new
//                 RefEdge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
RefEdge* OldUnmergeCode::split_out_Curves( DLIList<Curve*>& curve_list, CubitBoolean& reversed )
{
  int i;
  
  // If any parent Surfaces are still merged, we
  // cannot continue.
  DLIList<Surface*> surface_list;
  for( i = curve_list.size(); i > 0; i-- )
  {
    surface_list.clean_out();
    find_surfaces( curve_list.get_and_step(), surface_list );
    for( int j = surface_list.size(); j > 0; j-- )
      if( surface_list.get_and_step()->bridge_manager()
        ->number_of_bridges() != curve_list.size() )
        return NULL;
  }
  
  // Remove the curves from the old owner.
  curve_list.reset();
  Curve* curve = curve_list.get();
  TopologyEntity* te_ptr = curve->topology_entity();
  for( i = curve_list.size(); i > 0; i-- )
    if( curve_list.get_and_step()->topology_entity() != te_ptr )
      return 0;
      
  assert( te_ptr->bridge_manager()->number_of_bridges() > curve_list.size() );
  
  reversed = false;
  RefEdge* old_entity = CAST_TO( te_ptr, RefEdge );
  if( te_ptr )
  {
    assert( old_entity != 0 );
    
    curve_list.reset();
    for( i = curve_list.size(); i > 0; i-- )
    {
      Curve* ptr = curve_list.get_and_step();
      te_ptr->bridge_manager()->remove_bridge( ptr );
      remove_CAEntityId_attrib( ptr );
    }
    old_unmerged.append( old_entity );

    curve_list.reset();
    if (te_ptr->bridge_manager()->topology_bridge()->bridge_sense() == CUBIT_REVERSED)
    {
      te_ptr->reverse_topology();
      te_ptr->bridge_manager()->reverse_bridge_senses();
      reversed = !reversed;
    }
    
    if (curve->bridge_sense() == CUBIT_REVERSED)
    {
      for (i = curve_list.size(); i--; )
        curve_list.get_and_step()->reverse_bridge_sense();
      reversed = !reversed;
    }
  }
  
  curve_list.reset();
  for ( i = curve_list.size(); i--; )
    curve_list.get_and_step()->set_saved_id(0);
  
  RefEdge* new_entity = RefEntityFactory::instance()->construct_RefEdge( curve );
  Chain* chain = new Chain;
  new_entity->add_grouping_entity(chain);

  DLIList<TopologyBridge*> points(2);
  curve->get_children(points);
  points.reset();
  for (i = 0; i < 2; i++)
  {
    CoVertex* cvtx = new CoVertex;
    TopologyEntity* owner = points.get_and_step()->topology_entity();
    cvtx->attach_basic_topology_entity( dynamic_cast<RefVertex*>(owner) );
    chain->add_sense_entity(cvtx);
  }
  
  for( i = curve_list.size(); i > 0; i-- )
  {
    Curve* a_curve = curve_list.get_and_step();
    if( curve != a_curve )
      new_entity->bridge_manager()->add_bridge(a_curve);
  }  

  PRINT_DEBUG_19("  Unmerged %s from %s.\n", 
             new_entity->entity_name().c_str(), 
             old_entity->entity_name().c_str() );
  
  
  
  new_unmerged.append( new_entity );
  event_list.append( new UnMergeEvent( old_entity, new_entity ) );

  return new_entity;
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge a Point from its Vertex and make new Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
RefVertex* OldUnmergeCode::split_out_Points( DLIList<Point*>& point_list )
{
  int i;
  
  // If any parent Curves are still merged, we 
  // cannot continue.
  DLIList<Curve*> curve_list;
  for( i = point_list.size(); i > 0; i-- )
  {
    curve_list.clean_out();
    find_curves(point_list.get_and_step(),curve_list);
    for( int j = curve_list.size(); j > 0; j-- )
      if( curve_list.get_and_step()->bridge_manager()
        ->number_of_bridges() != point_list.size() )
        return NULL;
  }
  
  // Remove the child from the old owner.
  TopologyEntity* te_ptr = point_list.get()->topology_entity();
  
  // All points must have same bridge manager.
  for( i = point_list.size(); i > 0; i-- )
    if( point_list.get_and_step()->topology_entity() != te_ptr )
      return NULL;
  
  RefVertex* old_entity = CAST_TO( te_ptr, RefVertex );
  if( te_ptr )
  {
    assert(old_entity != 0);
    for( i = point_list.size(); i > 0; i-- )
    {
      Point* ptr = point_list.get_and_step();
      te_ptr->bridge_manager()->remove_bridge( ptr );
      remove_CAEntityId_attrib( ptr );
    }
    old_unmerged.append( old_entity );
  }

  point_list.reset();
  for (i = point_list.size(); i--; )
    point_list.get_and_step()->set_saved_id(0);

  RefVertex* new_entity = RefEntityFactory::instance()
    ->construct_RefVertex(point_list.get_and_step());
  for( i = point_list.size(); i > 1; i-- )
    new_entity->bridge_manager()->add_bridge(point_list.get_and_step());
    
  PRINT_DEBUG_19("    Unmerged %s from %s.\n", 
             new_entity->entity_name().c_str(), 
             old_entity->entity_name().c_str() );
  
  new_unmerged.append( new_entity );
  event_list.append( new UnMergeEvent( old_entity, new_entity ) );

  return new_entity;
}

//-------------------------------------------------------------------------
// Purpose       : Handle sending various events as a result of unmerging.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
void OldUnmergeCode::cleanup_unmerge()
{
  DLIList<RefEntity*> temp_list, marked_list;
  int i, j;
  CpuTimer timer;
  
  //Remove entities from unmerge_modified that
  // a) are duplicate entries in the list
  // b) also exist in new_unmerged
  // c) have parent entities already in unmerge_modified.
  
  //Change unmerge_modified to be a list of the immediate
  //parents of the entities currently in the list
  marked_list = unmerge_modified;
  unmerge_modified.clean_out();
  for( i = marked_list.size(); i--; )
  {
    temp_list.clean_out();
    marked_list.get_and_step()->get_parent_ref_entities( temp_list );
    unmerge_modified += temp_list;
  }
  for( i = new_unmerged.size(); i > 0; i-- )
  {
    temp_list.clean_out();
    new_unmerged.get_and_step()->get_parent_ref_entities( temp_list );
    unmerge_modified += temp_list;
  }
  
  //Reset marks on all entities in list
  for( i = unmerge_modified.size(); i > 0; i-- )
    unmerge_modified.get_and_step()->marked(0);

  //Set instance count (mark) to one for all in new_unmerged
  for( i = new_unmerged.size(); i > 0; i-- )
    new_unmerged.get_and_step()->marked(1);
    
  
    
  //Increment instance count (mark) for each occurance
  //in unmerge_modified.
  for( i = unmerge_modified.size(); i > 0; i-- )
  {
    unmerge_modified.get()->marked( unmerge_modified.get()->marked()+1 );
    temp_list.clean_out();
    unmerge_modified.get()->get_all_child_ref_entities(temp_list);
    marked_list += temp_list; 
    unmerge_modified.step();
  }
  
  for( i = marked_list.size(); i > 0; i-- )
  {
    marked_list.get()->marked( marked_list.get()->marked() + 1 );
    marked_list.step();
  }
  
  //If the entity's mark is non-zero after decrementing
  //the mark, remove it because it occurs again later in
  //the list.
  unmerge_modified.reset();
  for( i = unmerge_modified.size(); i > 0; i-- )
  {
    unmerge_modified.get()->marked( unmerge_modified.get()->marked() - 1 );
    if( unmerge_modified.get()->marked() )
      unmerge_modified.extract();
    else
      unmerge_modified.step();
  }

  //Clear out other marks we set.
  for( i = new_unmerged.size(); i > 0; i-- )
    new_unmerged.get_and_step()->marked(0);
  for( i = marked_list.size(); i > 0; i-- )
    marked_list.get_and_step()->marked(0);
  
  
  //Now remove duplicates from old_unmerged. 
  for( i = old_unmerged.size(); i > 0; i-- )
    old_unmerged.get_and_step()->marked(1);
  old_unmerged.reset();
  for( i = old_unmerged.size(); i > 0; i-- )
  {
    if( old_unmerged.get()->marked() )
      old_unmerged.get_and_step()->marked(0);
    else
      old_unmerged.extract();
  }
  
  
  // Finally, done rearranging lists.  Now actually do
  // the important stuff.
  
  // Count occurances of each type of entity, and 
  // output counts.
  int face_count = 0, edge_count = 0, vtx_count = 0;
  new_unmerged.reset();
  for( i = new_unmerged.size(); i > 0; i-- )
  {
    if( CAST_TO( new_unmerged.get(), RefFace ) )
       face_count++;
    else if( CAST_TO( new_unmerged.get(), RefEdge ) )
       edge_count++;
    else if( CAST_TO( new_unmerged.get(), RefVertex ) )
       vtx_count++;
    else
      assert( 0 /*Bad EntityType*/ );

    new_unmerged.step();
  }
  
  if( CubitMessage::instance()->Interrupt() ) PRINT_WARNING("Unmerge aborted.\n");
  if( face_count ) PRINT_INFO("%4d surfaces unmerged.\n",face_count);
  if( edge_count ) PRINT_INFO("%4d  curves  unmerged.\n",edge_count);
  if( vtx_count )  PRINT_INFO("%4d vertices unmerged.\n",vtx_count);
  PRINT_DEBUG_19("Unmerge completed in %0.2f seconds.\n",
                                         timer.cpu_secs() );
  
  //Let the user know if we actually unmerged anything.
  if( !new_unmerged.size() )
    PRINT_INFO("No entities unmerged.\n");
 
  int event_count = new_unmerged.size() + event_list.size() + unmerge_modified.size();
  ProgressTool* progress = 0;
  if (event_count > 20)
    progress = AppUtil::instance()->progress_tool();
  if (progress)
    progress->start( 0, event_count, "Updating Graphics" );
 
  //Info for DEBUG output
  const char* entity_names[] = { "Vertices", "Curves", "Surfaces" };
  timer.cpu_secs();
  
  //Add all new entities to the graphics.  Do lowest-dimension
  //entities first, so that they exist in the graphics when the
  //higher-dimension entities owning them are added.
  for( j = 0; j < 3; j++ )
  {
    int count = 0; //for PRINT_DEBUG statement below.
      
    for( i = new_unmerged.size(); i > 0; i-- )
    {
      RefEntity* entity = new_unmerged.get_and_step();
      if( entity->dimension() != j ) 
        continue;
      
      count++;
      
      //Is this a free RefEntity?
      temp_list.clean_out();
      entity->get_parent_ref_entities( temp_list );
      
      if( temp_list.size() == 0 )
        CubitObserver::notify_static_observers( entity, FREE_REF_ENTITY_GENERATED );

      if (progress)
        progress->step();
      PRINT_DEBUG_19("\t\tAdded %s to graphics\n", entity->entity_name().c_str() );

      //entity->notify_all_observers( NEW_ENTITY_UNMERGED );
    }
    PRINT_DEBUG_19("\tAdded %d new %s to graphics in %0.2f seconds.\n",
      count, entity_names[j], timer.cpu_secs() );
  }
  
  while( event_list.size() )
  {
    UnMergeEvent* event = event_list.pop();
    event->new_entity->notify_all_observers( *event );
    event->old_entity->notify_all_observers( *event );
    delete event;
    if (progress)
      progress->step();
  }
  
  //Update for other non-new entities that got modified
  //as a part of unmerging.
  for( i = unmerge_modified.size(); i > 0; i-- )
  {
    PRINT_DEBUG_19("\t\tUpdated graphics for %s.\n", unmerge_modified.get()->entity_name().c_str() );
    CubitObserver::notify_static_observers( unmerge_modified.get_and_step(), TOPOLOGY_MODIFIED );
    if (progress)
      progress->step();  
  }
  if (progress)
    progress->end();
  
  PRINT_DEBUG_19("Sent TOPOLOGY_MODIFIED event for %d entities in %0.2f "
                 "seconds\n",unmerge_modified.size(), timer.cpu_secs());
  PRINT_DEBUG_19("Total time to update after unmerge: %0.2f seconds.\n",
    timer.elapsed() );

  DLIList<MergeToolAssistant*>& assistant_list_ = MergeTool::instance()->assistant_list_;
  for( int a = assistant_list_.size(); a--; )
    assistant_list_.get_and_step()->finish_unmerge();
}

    

