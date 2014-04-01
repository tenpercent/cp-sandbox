//- Class: MergeTool
//- Description: Location of all Merge and Unmerge functionality
//- Owner: Steven Jankovich
//- Created: 27 April 1998
//- Checked by:

#include <assert.h>

#include "GeometryEntity.hpp"
#include "MergeTool.hpp"
#include "MergeToolAssistant.hpp"
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "ModelQueryEngine.hpp"

#include "CubitObserver.hpp"
#include "MergeEvent.hpp"
#include "CubitUtil.hpp"
#include "RefEntity.hpp"
#include "BasicTopologyEntity.hpp"
#include "GroupingEntity.hpp"
#include "SenseEntity.hpp"

#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "RefGroup.hpp"

#include "CoVertex.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"

#include "Chain.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Body.hpp"
#include "BodySM.hpp"

#include "ToolData.hpp"
#include "TDCompare.hpp"
#include "UnMergeEvent.hpp"

#include "DLIList.hpp"
#include "DRefFaceArray.hpp"
#include "DRefEdgeArray.hpp"
#include "DRefVertexArray.hpp"

#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"
#include "ShellSM.hpp"

#include "CpuTimer.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"
#include "CastTo.hpp"
#include "RTree.hpp"
#include "AbstractTree.hpp"
#include "SettingHandler.hpp"

MergeTool* MergeTool::instance_ = NULL;
CubitBoolean MergeTool::mergeCalled = CUBIT_FALSE;
CubitBoolean MergeTool::groupResults = CUBIT_FALSE;
CubitBoolean MergeTool::destroyDeadGeometry = CUBIT_TRUE;

// Constructor
MergeTool* MergeTool::instance()
{
   if( instance_ == NULL )
       instance_ = new MergeTool;
   return instance_;
}

void MergeTool::initialize_settings( )
{
}

MergeTool::MergeTool()
{
     //This is private.  It can only be accessed through
     //the instance() function above.
   mergeCalled = CUBIT_FALSE;
   
   unmerged_list_in_use = CUBIT_FALSE;
   displayProgress = false;
   destroyDeadGeometry = true;
     //These shouldn't get deleted, they
     //get deleted when other groups get deleted (reset).
   lastSurfsMerged = NULL;
   lastCurvsMerged = NULL;
   lastVertsMerged = NULL;
}

// Destructor
MergeTool::~MergeTool()
{
  while( assistant_list_.size() )
    delete assistant_list_.pop();
  instance_ = 0;
}

//Public Functions:
CubitBoolean MergeTool::contains_merged_entities( DLIList<RefEntity*> &ref_entities)
{
     //Loop through the entities and their children to find
     //merged entities.  For now, just do it quickly, so
     //if we find a merged entity, return true.

  DLIList<RefEntity*> all_entities, temp_entities;

  ref_entities.reset();
  ref_entities.get()->get_all_child_ref_entities(ref_entities, temp_entities);
  temp_entities += ref_entities;

  all_entities.merge_unique(temp_entities);
  
  for (int i = all_entities.size(); i > 0; i--) {
    RefEntity *temp_entity = all_entities.get_and_step();
    if (entity_merged(CAST_TO(temp_entity, TopologyEntity)))
      return CUBIT_TRUE;
  }
  
  return CUBIT_FALSE;
}

CubitBoolean MergeTool::contains_merged_children( Body *body, 
                                                  DLIList<RefEntity*> &merged_children)
{
  
  RefEntity *ref_ent = NULL;
  ref_ent = CAST_TO(body, RefEntity);
  DLIList<RefEntity*> ref_ent_list;
  ref_ent->get_all_child_ref_entities( ref_ent_list );
  int i;
  for( i=ref_ent_list.size(); i--;)
  {
    ref_ent = ref_ent_list.get_and_step();
    if(entity_merged(CAST_TO(ref_ent, TopologyEntity)))
      merged_children.append( ref_ent );
  }

  if( merged_children.size() )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}


CubitBoolean MergeTool::contains_merged_entities( DLIList<Body*> &bodies )
{
     //Loop through the bodies and their entities to find
     //merged entities.  For now, just do it quickly, so
     //if we find a merged entity, return true.
  DLIList<RefEntity*> ref_entities;
  CAST_LIST_TO_PARENT(bodies, ref_entities);
  return contains_merged_entities(ref_entities);
}

CubitBoolean MergeTool::entity_merged( TopologyEntity *entity )
{
  if (entity->bridge_manager()->number_of_bridges() > 1)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

CubitStatus MergeTool::merge_all_bodies()
{
   PRINT_INFO( "\n...Merging all features in the model\n" );
   set_merge_occurance( CUBIT_TRUE );
   
   if( merge_all_reffaces() == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Surface merging failed\n" );
      return CUBIT_FAILURE;
   }
   
   if( merge_all_refedges() == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Curve merging failed\n" );
      return CUBIT_FAILURE;
   }
   
   if( merge_all_refvertices() == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Vertex merging failed\n" );
      return CUBIT_FAILURE;
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_bodies( DLIList<Body*>& body_list )
{
   set_merge_occurance( CUBIT_TRUE );
   
   DLIList<ModelEntity*> query_input(body_list.size()), query_output;
   CAST_LIST_TO_PARENT(body_list, query_input);
   ModelQueryEngine *const mqe = ModelQueryEngine::instance();
   
     // Get the RefFaces
   query_output.clean_out();
   mqe->query_model( query_input, DagType::ref_face_type(), query_output );
   DLIList<RefFace*> refface_list(query_output.size());
   CAST_LIST(query_output, refface_list, RefFace);
   
     // Merge the RefFaces
   if( merge_reffaces( refface_list ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Surface merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
     // Get the RefEdges
   query_output.clean_out();
   mqe->query_model( query_input, DagType::ref_edge_type(), query_output );
   DLIList<RefEdge*> refedge_list(query_output.size());
   CAST_LIST(query_output, refedge_list, RefEdge);
   
     // Merge the RefEdges
   if( merge_refedges( refedge_list ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Curve merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
     // Get the RefVertices
   query_output.clean_out();
   mqe->query_model( query_input, DagType::ref_vertex_type(), query_output );
   DLIList<RefVertex*> refvertex_list(query_output.size());   
   CAST_LIST(query_output, refvertex_list, RefVertex);
   
     // Merge the RefVertices
   if( merge_refvertices( refvertex_list ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Vertex merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_volumes( DLIList<RefVolume*>& vol_list,
                                      CubitBoolean print_info )
{
   int i;
   set_merge_occurance( CUBIT_TRUE );
   vol_list.reset();
   
     // Get the RefFaces
   DLIList<RefFace*> refface_list;
   for( i = vol_list.size(); i > 0; i-- )
       vol_list.get_and_step()->ref_faces( refface_list );
   
     // Merge the RefFaces
   if( merge_reffaces( refface_list, print_info ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Surface merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
     // Get the RefEdges
   DLIList<RefEdge*> refedge_list;
   for( i = vol_list.size(); i > 0; i-- )
       vol_list.get_and_step()->ref_edges( refedge_list );
   
     // Merge the RefEdges
   if( merge_refedges( refedge_list, CUBIT_TRUE, print_info ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Curve merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
     // Get the RefVertices
   DLIList<RefVertex*> refvertex_list;
   for( i = vol_list.size(); i > 0; i-- )
       vol_list.get_and_step()->ref_vertices( refvertex_list );
   
     // Merge the RefVertices
   if( merge_refvertices( refvertex_list, print_info ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Vertex merging failed\n" );
      return CUBIT_FAILURE;      
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_all_reffaces()
{
//   if( !displayProgress )
     PRINT_INFO( "\n...Merging all Surfaces in the model\n" );
   set_merge_occurance( CUBIT_TRUE );
   
     // Get the list of all the RefFaces in the model
   DLIList<RefFace*> refface_list;
   GeometryQueryTool::instance()->ref_faces(refface_list);
   
     // Merge the RefFaces in this list
   if ( DEBUG_FLAG(139) )
   {
     test_r_tree(refface_list);
//     test_r_star_tree(refface_list);
     test_no_tree(refface_list);
   }
   if( merge_reffaces( refface_list ) == CUBIT_FAILURE )
   {
     PRINT_ERROR( "Surface merging failed\n" );
     return CUBIT_FAILURE;
   }

   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_reffaces_old( DLIList<RefFace*>& refface_list,
                                           CubitBoolean print_info)
{
    // The basic algorithm is to step through the input list,
    // compare every entity with every other entity, and
    // merge the entities that are spatially equivalent.
    // In the merge, the entity with the lowest ID is retained.
    // This is complicated by the fact that once an entity has
    // been merged out, we can no longer access that pointer
    // since it is invalid. So...we need to remove it from the
    // list, but still keep track of where we are in the list.
    // If the entity that is deleted is before the retained item
    // in the list, then we do not want to step or we will....
  set_merge_occurance( CUBIT_TRUE );
  
  CpuTimer timer;
  int merge_count = 0;
  
  DRefFaceArray refface_array( refface_list.size() );
  refface_list.reset();
  
    // Remove entities that should not be automatically merged
  int i = 0;
  for( i = refface_list.size(); i > 0; i-- )
  {
    RefFace *curr_face = refface_list.get_and_step();
    if( curr_face->is_mergeable() )
      refface_array.append( curr_face );
  }
  
  ProgressTool* progress = 0;
  if( displayProgress )
  {
    char info[64];
    sprintf(info, "Comparing %d Surfaces for Merge", refface_array.size() );
    progress = AppUtil::instance()->progress_tool();
    if (progress)
    {
      progress->start( 0, refface_array.size(), "Comparing Surfaces:", 
                       info, CUBIT_TRUE, CUBIT_TRUE );
    }
  }
  
    // Now find overlapping RefFaces and merge them.
    // Make sure that the operation is not performed on
    // RefFaces that are deactivated.
  int j = 0;
  RefFace* refface_ptr = NULL;
  RefFace* compare_refface_ptr = NULL;
  int array_size = refface_array.size();
  DLIList<RefFace*> faces_merged;
  
  for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
  {
    if( progress )
       progress->step();
    
    refface_ptr = refface_array[i];
    if( refface_ptr == NULL )
      continue;
    
      // There is no need to compare if the first RefFace is
      // deactivated.
    if( refface_ptr->deactivated() == CUBIT_TRUE )
      continue ;
    
    j = i+1;
    for( j = i+1; (j < array_size) && !CubitMessage::instance()->Interrupt(); j++ )
    {
      compare_refface_ptr = refface_array[j];
      if( compare_refface_ptr == NULL )
        continue;
      
        // Make sure we are not comparing the same entity
        // and that they are not deactivated.
      if( ( refface_ptr == compare_refface_ptr ) ||
          ( compare_refface_ptr->deactivated() == CUBIT_TRUE ) )
        continue;
      
        // IMPORTANT: this compare is for merging, so we set the
        //            last parameter to CUBIT_TRUE so that the
        //            compare sets the TDCompare data which will
        //            be used in the RefEntity compare and merge.
        // By Jihong
      
      double geom_factor = GeometryQueryTool::get_geometry_factor();
      CubitBoolean status =
        refface_ptr->about_spatially_equal( compare_refface_ptr,
                                            geom_factor,
                                            CUBIT_TRUE,
                                            GeometryQueryTool::instance()->get_merge_test_bbox(),
                                            GeometryQueryTool::instance()->get_merge_test_internal() );
      if( status == CUBIT_FALSE )
        continue;
      
        // If we are in this block, we want to merge the 
        // two entities. If they do not merge, there was
        // an error somewhere.
      
        //First test to see if all the children of these
        // reffaces are mergeable
      if( !compare_refface_ptr->children_mergeable() )
      {
        PRINT_WARNING( "Cannot merge surfaces %d and %d\n"
                       "     Make sure all merge flags are on.\n",
                       refface_ptr->id(), compare_refface_ptr->id() );
        continue;
      }
      
        //Now let us test to see if we are merging two faces that
        //belong to the same volume.
      DLIList<RefVolume*> ref_vols_1, ref_vols_2;
      refface_ptr->ref_volumes( ref_vols_1 );
      compare_refface_ptr->ref_volumes( ref_vols_2 );
      ref_vols_2.intersect( ref_vols_1 );
      if( ref_vols_2.size() > 0 )
      {
          PRINT_ERROR( "Tolerance problems, trying to merge "
                       "surfaces\non the same volume.\n"
                       "  %s (surface %d) and %s (surface %d) on\n"
                       "  %s (volume %d)\n",
                       refface_ptr->entity_name().c_str(),
                       refface_ptr->id(),
                       compare_refface_ptr->entity_name().c_str(),
                       compare_refface_ptr->id(),
                       ref_vols_2.get()->entity_name().c_str(),
                       ref_vols_2.get()->id() );
          if (print_info) PRINT_INFO( "Try changing the merge tolerance.\n" );
          continue;
      }

        // Always retain the entity with the lowest id.
      int nullify = j;
      if( refface_ptr->id() > compare_refface_ptr->id() )
      {
        std::swap(refface_ptr, compare_refface_ptr);
        nullify  = i;
      }
      if (groupResults )
      {
        faces_merged.append(refface_ptr);
      }

        // Now check if merge is okay with all assistants.
      CubitBoolean assistant_says_no = CUBIT_FALSE;
      for( int k = assistant_list_.size(); k > 0; k-- )
      {
        if( ! assistant_list_.get_and_step()
          ->can_merge( refface_ptr, compare_refface_ptr ) )
        {
          assistant_says_no = CUBIT_TRUE;
          break;
        }
      }
      if( assistant_says_no )
        continue; 

        // Need to retain these so that the pointers are not
        // accessed after a merge operation when the 'deleted'
        // pointer may be invalid.
      int retained_id = refface_ptr->id();
      int deleted_id  = compare_refface_ptr->id();

      PRINT_DEBUG_19( "Consolidating RefFace %d and "
                      "%d...\n", retained_id, deleted_id );
      if( merge_BTE( refface_ptr, compare_refface_ptr ) )
      {
        merge_count++;
        if (print_info && !progress) 
          PRINT_INFO( "Surface %d and %d consolidated\n",
                                    retained_id, deleted_id );
        
          // The 'deleted' RefFace is now gone. It is an
          // error to access that pointer, so we need to
          // get it out of the list.
        refface_array[nullify] = NULL;
      }
      else
      {
        PRINT_ERROR( "Surface %d and %d NOT consolidated\n"
                     "       The use counts were updated first and"
                     " are now out of date.\n",
                     retained_id, deleted_id);
      }
      
      if (nullify == i)
        break;
    }
  }
  
    // Remove the crud that accumulated during the merge
    // operations, from the geometry database.

  complete_merge();
  if (progress)
      progress->end();
  PRINT_DEBUG_3( "Merge Reffaces time: %f secs\n",
                 timer.cpu_secs() );
  
  PRINT_DEBUG_19("Cleaning out TDCompare data from RefFaces...\n");
   
  if ( groupResults && faces_merged.size() )
  {
    DLIList<RefEntity*> refentity_list;
    RefFace *tmp_face;
    for (int iii = faces_merged.size(); iii > 0; iii-- )
    {
      tmp_face = faces_merged.get_and_step();
      if ( !tmp_face->deactivated() )
        refentity_list.append(tmp_face);
    }
    RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_surfs_merged");
    new_group->add_ref_entity( refentity_list );
    if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                               "  Group contains surfaces that were merged during\n"
                               "  current merge operation\n",
                               new_group->entity_name().c_str(),
                               new_group->id());
    lastSurfsMerged = new_group;
  }
    //since surface mergeing has occured and we either don't want the
    //group or didn't store it, set the lastsurfs to null.
  else
    lastSurfsMerged = NULL;
  
  if( destroyDeadGeometry )
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  
  PRINT_DEBUG_3( "cleanout time: %f secs\n", timer.cpu_secs() );
  if (print_info)
  {
    PRINT_INFO( "Consolidated %d pair", merge_count);
    if(merge_count > 1)
       PRINT_INFO("s");
    PRINT_INFO( " of surfaces\n");
  }
  
  if( CubitMessage::instance()->Interrupt() )
  {
     PRINT_WARNING("Surface merging aborted.\n");
     return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_reffaces( DLIList<RefFace*>& refface_list,
                                       CubitBoolean print_info)
{
  if (refface_list.size() < 20)
  {
    if ( merge_reffaces_old(refface_list, print_info) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
    else
      return CUBIT_SUCCESS;
  }
    // The basic algorithm is to step through the input list,
    // compare every entity with every entity within range of
    // its bounding box and merge the spatially equivellant entities.
    // An R-Tree is used to efficiently find the entities within range.
    // In the merge, the entity with the lowest ID is retained.
    // This is complicated by the fact that once an entity has
    // been merged out, we can no longer access that pointer
    // since it is invalid. So...we need to remove it from the
    // tree, but still keep track of where we are in the list.
  set_merge_occurance( CUBIT_TRUE );
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  AbstractTree <RefFace*> *a_tree = new RTree<RefFace*> (GEOMETRY_RESABS*geom_factor);
  CpuTimer timer;
  int merge_count = 0;
  
  DRefFaceArray refface_array( refface_list.size() );
  refface_list.reset();
  
    // Remove entities that should not be automatically merged
  int i = 0;
  int loop_size = refface_list.size();
  CpuTimer time_to_build;
  for( i = 0; i < loop_size; i++ )
  {
    RefFace *curr_face = refface_list.get_and_step();
    if( curr_face->is_mergeable() )
    {
      refface_array.append( curr_face );
      a_tree->add(curr_face);
    }
  }
  PRINT_DEBUG_3( "Time to build r_tree %f secs, with %d entries\n",
                 time_to_build.cpu_secs(), refface_array.size() );

    //initialize the marked flag for fast nulification...
  int array_size = refface_array.size();
  for ( i = 0; i < array_size; i++)
    refface_array[i]->marked(i);
  
  ProgressTool* progress = 0;
  if( displayProgress )
  {
    char info[64];
    sprintf(info, "Comparing %d Surfaces for Merge", refface_array.size() );
    progress = AppUtil::instance()->progress_tool();
    if (progress)
    {
      progress->start( 0, refface_array.size(), "Comparing Surfaces:", 
                       info, CUBIT_TRUE, CUBIT_TRUE );
    }
  }
  
    // Now find overlapping RefFaces and merge them.
    // Make sure that the operation is not performed on
    // RefFaces that are deactivated.
  int j = 0;
  RefFace* refface_ptr = NULL;
  RefFace* compare_refface_ptr = NULL;
  DLIList<RefFace*> faces_merged, faces_in_range;
  CubitBox temp_box;
  for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
  {
    if( progress ) progress->step();
    
    refface_ptr = refface_array[i];
    if( refface_ptr == NULL )
      continue;
    
      // There is no need to compare if the first RefFace is
      // deactivated.
    if( refface_ptr->deactivated() == CUBIT_TRUE )
      continue ;

      //Now from the atree, get the surfaces that are close.
    temp_box = refface_ptr->bounding_box();
    faces_in_range.clean_out();
    a_tree->find(temp_box, faces_in_range);
    
    for( j = 0; (j < faces_in_range.size()) && !CubitMessage::instance()->Interrupt(); j++ )
    {
      compare_refface_ptr = faces_in_range.get_and_step();
      if( compare_refface_ptr == NULL )
        continue;
      
        // Make sure we are not comparing the same entity
        // and that they are not deactivated.
      if( ( refface_ptr == compare_refface_ptr ) ||
          ( compare_refface_ptr->deactivated() == CUBIT_TRUE ) )
        continue;
      
        // IMPORTANT: this compare is for merging, so we set the
        //            last parameter to CUBIT_TRUE so that the
        //            compare sets the TDCompare data which will
        //            be used in the RefEntity compare and merge.
        // By Jihong
      
      geom_factor = GeometryQueryTool::get_geometry_factor();
      CubitBoolean status =
        refface_ptr->about_spatially_equal( compare_refface_ptr,
                                            geom_factor,
                                            CUBIT_TRUE,
                                            GeometryQueryTool::instance()->get_merge_test_bbox(),
                                            GeometryQueryTool::instance()->get_merge_test_internal() );
      if( status == CUBIT_FALSE )
        continue;
      
        // If we are in this block, we want to merge the 
        // two entities. If they do not merge, there was
        // an error somewhere.
      
        //First test to see if all the children of these
        // reffaces are mergeable
      if( !compare_refface_ptr->children_mergeable() )
      {
        PRINT_WARNING( "Cannot merge surfaces %d and %d\n"
                       "     Make sure all merge flags are on.\n",
                       refface_ptr->id(), compare_refface_ptr->id() );
        continue;
      }
      
        //Now let us test to see if we are merging two faces that
        //belong to the same volume.
      DLIList<RefVolume*> ref_vols_1, ref_vols_2;
      refface_ptr->ref_volumes( ref_vols_1 );
      compare_refface_ptr->ref_volumes( ref_vols_2 );
      ref_vols_2.intersect( ref_vols_1 );
      if( ref_vols_2.size() > 0 )
      {
          PRINT_ERROR( "Tolerance problems, trying to merge "
                       "surfaces\non the same volume.\n"
                       "  %s (surface %d) and %s (surface %d) on\n"
                       "  %s (volume %d)\n",
                       refface_ptr->entity_name().c_str(),
                       refface_ptr->id(),
                       compare_refface_ptr->entity_name().c_str(),
                       compare_refface_ptr->id(),
                       ref_vols_2.get()->entity_name().c_str(),
                       ref_vols_2.get()->id() );
          if (print_info) PRINT_INFO( "Try changing the merge tolerance.\n" );
          continue;
      }
 
        // Always retain the entity with the lowest id.
      int nullify = compare_refface_ptr->marked();
      if( refface_ptr->id() > compare_refface_ptr->id() )
      {
        std::swap(refface_ptr, compare_refface_ptr);
        nullify = i;
      }
      if (groupResults )
      {
        faces_merged.append(refface_ptr);
      }

        // Now check if merge is okay with all assistants.
      CubitBoolean assistant_says_no = CUBIT_FALSE;
      for( int k = assistant_list_.size(); k > 0; k-- )
      {
        if( ! assistant_list_.get_and_step()
          ->can_merge( refface_ptr, compare_refface_ptr ) )
        {
          assistant_says_no = CUBIT_TRUE;
          break;
        }
      }
      if( assistant_says_no )
        continue; 

        // Need to retain these so that the pointers are not
        // accessed after a merge operation when the 'deleted'
        // pointer may be invalid.
      int retained_id = refface_ptr->id();
      int deleted_id  = compare_refface_ptr->id();
      
      PRINT_DEBUG_19( "Consolidating RefFace %d and "
                      "%d...\n", retained_id, deleted_id );
      
      a_tree->remove(compare_refface_ptr);
      if( merge_BTE( refface_ptr, compare_refface_ptr ) )
      {
        merge_count++;
        if (print_info && !progress) 
          PRINT_INFO( "Surface %d and %d consolidated\n",
                      retained_id, deleted_id );
        
          // The 'deleted' RefFace is now gone. It is an
          // error to access that pointer, so we need to
          // get it out of the list.
        refface_array[nullify] = NULL;
      }
      else
      {
        PRINT_ERROR( "Surface %d and %d NOT consolidated\n"
                     "       The use counts were updated first and"
                     " are now out of date.\n",
                     retained_id, deleted_id);
        a_tree->add(compare_refface_ptr);
      }
      
      if (nullify == i)
        break;
      
    }
  }
    //clean the marks.
  for( i = 0; i < array_size; i++ )
  {
    if ( refface_array[i] != NULL )
      refface_array[i]->marked(0);
  }
  complete_merge();
  if (progress)
      progress->end();
  CpuTimer time_to_destroy;
  delete a_tree;
  PRINT_DEBUG_3( "Time to destroy r_tree %f secs.\n",
                  time_to_destroy.cpu_secs());

      // Remove the crud that accumulated during the merge
    // operations, from the geometry database.
  PRINT_DEBUG_3( "Merge Reffaces time: %f secs\n",
                 timer.cpu_secs() );
  
  PRINT_DEBUG_19("Cleaning out TDCompare data from RefFaces...\n");
  if ( groupResults && faces_merged.size() )
  {
    DLIList<RefEntity*> refentity_list;
    RefFace *tmp_face;
    for (int iii = faces_merged.size(); iii > 0; iii-- )
    {
      tmp_face = faces_merged.get_and_step();
      if ( !tmp_face->deactivated() )
        refentity_list.append(tmp_face);
    }
    RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_surfs_merged");
    new_group->add_ref_entity( refentity_list );
    if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                               "  Group contains surfaces that were merged during\n"
                               "  current merge operation\n",
                               new_group->entity_name().c_str(),
                               new_group->id());
    lastSurfsMerged = new_group;
  }
    //since surface mergeing has occured and we either don't want the
    //group or didn't store it, set the lastsurfs to null.
  else
    lastSurfsMerged = NULL;
  
  if( destroyDeadGeometry )
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();

  //for updating GUI tree
  CubitObserver::notify_static_observers(NULL, MERGE_COMPLETED);

  PRINT_DEBUG_3( "cleanout time: %f secs\n", timer.cpu_secs() );
  if (print_info)
  {
    PRINT_INFO( "Consolidated %d pair", merge_count);
    if(merge_count > 1)
       PRINT_INFO("s");
    PRINT_INFO( " of surfaces\n");
  }
  if( CubitMessage::instance()->Interrupt() )
  {
     PRINT_WARNING("Surface merging aborted.\n");
     return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}
 
CubitStatus MergeTool::merge_all_refedges()
{
//   if( !displayProgress )
     PRINT_INFO( "\n...Merging all Curves in the model\n" );
   set_merge_occurance( CUBIT_TRUE );
   
     // Get the list of all the RefEdges in the model
   DLIList<RefEdge*> refedge_list_model;
   GeometryQueryTool::instance()->ref_edges(refedge_list_model);
   
     // Merge the RefEdges in this list
   if( merge_refedges( refedge_list_model ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Curve merging failed\n" );
      return CUBIT_FAILURE;
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_mergeable_refentities( DLIList<RefEntity*> &entities,
                                                   DLIList< DLIList<RefFace*>*> &lists_of_mergeable_ref_faces,
                                                   DLIList< DLIList<RefEdge*>*> &lists_of_mergeable_ref_edges,
                                                   DLIList< DLIList<RefVertex*>*> &lists_of_mergeable_ref_vertices)
{
  DLIList<RefEntity*> tmp_ents;

  int i;
  
  //only want surfaces and above 
  for( i=entities.size(); i--; )
  {
    RefEntity *ref_ent = entities.get_and_step();
    if( ref_ent->dimension() > 1 || ref_ent->dimension() == -1 )
      tmp_ents.append( ref_ent );
  }
  find_mergeable_reffaces( tmp_ents, lists_of_mergeable_ref_faces, false );

  //only want curves and above 
  tmp_ents.clean_out();
  for( i=entities.size(); i--; )
  {
    RefEntity *ref_ent = entities.get_and_step();
    if( ref_ent->dimension() > 0 || ref_ent->dimension() == -1 )
      tmp_ents.append( ref_ent );
  }
  find_mergeable_refedges( tmp_ents, lists_of_mergeable_ref_edges, false );

  find_mergeable_refvertices( entities, lists_of_mergeable_ref_vertices, false);
  
/* 
  //remove the compare partners off everything
  DLIList<RefEntity*> all_ents;
  RefEntity::get_all_child_ref_entities( entities, all_ents );
  PRINT_INFO("Got %d children\n", all_ents.size() );
  all_ents += entities;

  int i;
  for( i=all_ents.size(); i--; )
    all_ents.get_and_step()->remove_compare_data(); 
*/
  remove_compare_data();

  return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_mergeable_reffaces( DLIList<RefEntity*> &entities,
                                                DLIList< DLIList<RefFace*>*> &lists_of_mergeable_ref_faces,
                                                bool clean_up_compare_data )
{
  DLIList<TopologyEntity*> t_ents;
  CAST_LIST(entities, t_ents, TopologyEntity);

  DLIList<RefFace*> refface_list;
  int i;
  for( i=t_ents.size(); i--; )
  {
    DLIList<RefFace*> tmp_faces;
    t_ents.get_and_step()->ref_faces( tmp_faces );
    refface_list += tmp_faces;
  }

  double geom_factor = GeometryQueryTool::get_geometry_factor();
  AbstractTree <RefFace*> *a_tree = new RTree<RefFace*> (GEOMETRY_RESABS*geom_factor);
  
  DRefFaceArray refface_array( refface_list.size() );
  refface_list.reset();
  
    // Remove entities that should not be automatically merged
  int loop_size = refface_list.size();
  CpuTimer time_to_build;
  for( i = 0; i < loop_size; i++ )
  {
    RefFace *curr_face = refface_list.get_and_step();
    if( curr_face->is_mergeable() )
    {
      refface_array.append( curr_face );
      a_tree->add(curr_face);
    }
  }
  PRINT_DEBUG_3( "Time to build r_tree %f secs, with %d entries\n",
                 time_to_build.cpu_secs(), refface_array.size() );

    //initialize the marked flag for fast nulification...
  int array_size = refface_array.size();
  for ( i = 0; i < array_size; i++)
    refface_array[i]->marked(i);
  
  ProgressTool* progress = 0;
  if( displayProgress )
  {
    char info[64];
    sprintf(info, "Comparing %d Surfaces for Merge", refface_array.size() );
    progress = AppUtil::instance()->progress_tool();
    if (progress)
    {
      progress->start( 0, refface_array.size(), "Comparing Surfaces:", 
                       info, CUBIT_TRUE, CUBIT_TRUE );
    }
  }
  
    // Now find overlapping RefFaces and merge them.
    // Make sure that the operation is not performed on
    // RefFaces that are deactivated.
  int j = 0;
  RefFace* refface_ptr = NULL;
  RefFace* compare_refface_ptr = NULL;
#ifdef BOYD17 
  DLIList<RefFace*> faces_merged, faces_in_range;
#endif
  DLIList<RefFace*> faces_in_range;
  
  CubitBox temp_box;
  for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
  {
    if( progress ) progress->step();
    
    refface_ptr = refface_array[i];
    if( refface_ptr == NULL )
      continue;
    
      // There is no need to compare if the first RefFace is
      // deactivated.
    if( refface_ptr->deactivated() == CUBIT_TRUE )
      continue ;

      //Now from the atree, get the surfaces that are close.
    temp_box = refface_ptr->bounding_box();
    faces_in_range.clean_out();
    a_tree->find(temp_box, faces_in_range);

    DLIList<RefVolume*> tmp_vols;
    refface_ptr->ref_volumes( tmp_vols);
    RefVolume* tmp_vol= tmp_vols.get();

    DLIList<RefFace*> *new_list = NULL;
    for( j = 0; (j < faces_in_range.size()) && !CubitMessage::instance()->Interrupt(); j++ )
    {
      compare_refface_ptr = faces_in_range.get_and_step();
      if( compare_refface_ptr == NULL )
        continue;
    
      if( compare_refface_ptr == refface_ptr )
        continue;

      //surfaces in the same bodies won't merge
      tmp_vols.clean_out();
      compare_refface_ptr->ref_volumes( tmp_vols );
      RefVolume* compare_vol= tmp_vols.get();
      if( tmp_vol == compare_vol )
        continue;

        // Make sure we are not comparing the same entity
        // and that they are not deactivated.
      if( ( refface_ptr == compare_refface_ptr ) ||
          ( compare_refface_ptr->deactivated() == CUBIT_TRUE ) )
        continue;
      
      geom_factor = GeometryQueryTool::get_geometry_factor();
      CubitBoolean status =
        refface_ptr->about_spatially_equal( compare_refface_ptr,
                                            geom_factor,
                                            CUBIT_TRUE,
                                            GeometryQueryTool::instance()->get_merge_test_bbox(),
                                            GeometryQueryTool::instance()->get_merge_test_internal() );
      if( status == CUBIT_FALSE )
        continue;
      
        //First test to see if all the children of these
        // reffaces are mergeable
      if( !compare_refface_ptr->children_mergeable() )
        continue;
      
        //Now let us test to see if we are merging two faces that
        //belong to the same volume.
      DLIList<RefVolume*> ref_vols_1, ref_vols_2;
      refface_ptr->ref_volumes( ref_vols_1 );
      compare_refface_ptr->ref_volumes( ref_vols_2 );
      ref_vols_2.intersect( ref_vols_1 );
      if( ref_vols_2.size() > 0 )
      {
          PRINT_ERROR( "Tolerance problems, trying to find mergeable "
                       "surfaces\non the same volume.\n"
                       "  %s (surface %d) and %s (surface %d) on\n"
                       "  %s (volume %d)\n",
                       refface_ptr->entity_name().c_str(),
                       refface_ptr->id(),
                       compare_refface_ptr->entity_name().c_str(),
                       compare_refface_ptr->id(),
                       ref_vols_2.get()->entity_name().c_str(),
                       ref_vols_2.get()->id() );
          continue;
      }
 
        // Always retain the entity with the lowest id.
      int nullify = compare_refface_ptr->marked();

        // Now check if merge is okay with all assistants.
      CubitBoolean assistant_says_no = CUBIT_FALSE;
      for( int k = assistant_list_.size(); k > 0; k-- )
      {
        if( ! assistant_list_.get_and_step()
          ->can_merge( refface_ptr, compare_refface_ptr ) )
        {
          assistant_says_no = CUBIT_TRUE;
          break;
        }
      }
      if( assistant_says_no )
        continue; 

      a_tree->remove(compare_refface_ptr);

      refface_array[nullify] = NULL;

      if( new_list == NULL ) 
      {
        new_list = new DLIList<RefFace*>; 
        new_list->append( refface_ptr );
        new_list->append( compare_refface_ptr );
      }
      else
        new_list->append( compare_refface_ptr );
    }

    if( new_list )
      lists_of_mergeable_ref_faces.append( new_list );
  }
    //clean the marks.
  for( i = 0; i < array_size; i++ )
  {
    if ( refface_array[i] != NULL )
      refface_array[i]->marked(0);
  }

  if( clean_up_compare_data )
    remove_compare_data();

  if (progress)
      progress->end();

  delete a_tree;

  if( CubitMessage::instance()->Interrupt() )
  {
     PRINT_WARNING("Finding mergeable surfaces aborted.\n");
     return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_mergeable_refedges( DLIList<RefEntity*> &entities,
                                                DLIList< DLIList<RefEdge*>*> &lists_of_mergeable_ref_edges,
                                                bool clean_up_compare_data )
{
  
  DLIList<TopologyEntity*> t_ents;
  CAST_LIST(entities, t_ents, TopologyEntity);

  DLIList<RefEdge*> refedge_list;
  int i;
  for( i=t_ents.size(); i--; )
  {
    DLIList<RefEdge*> tmp_edges;
    t_ents.get_and_step()->ref_edges( tmp_edges );
    refedge_list += tmp_edges;
  }

  if( refedge_list.size() == 0 )
    return CUBIT_SUCCESS;

   DRefEdgeArray refedge_array( refedge_list.size() );
   refedge_list.reset();

   DLIList<RefEdge*> edges_on_two_surfs; 
  
   for( i = refedge_list.size(); i > 0; i-- )
   {
      // Remove entities that should not be automatically merged
      RefEdge *curr_edge = refedge_list.get_and_step();
      if( curr_edge->is_mergeable() )
      {
        //put all edges in the list first that are attached 
        //to only one or zero surface.  This prevents curves on sheet
        //bodies from being excluded. 
        DLIList<RefFace*> tmp_faces;
        curr_edge->ref_faces( tmp_faces );
        if( tmp_faces.size() < 2 )
          refedge_array.append( curr_edge );
        else
          edges_on_two_surfs.append( curr_edge );
      }
   }
    
   for( i = edges_on_two_surfs.size(); i--; )
     refedge_array.append( edges_on_two_surfs.get_and_step() ); 

#ifdef BOYD17
   DLIList<RefEdge*> edges_merged;
#endif
  
   ProgressTool* progress = 0;

   int j = 0;
   RefEdge* refedge_ptr = NULL;
   RefEdge* compare_refedge_ptr = NULL;
   int array_size = refedge_array.size();
   for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
   {
      if( progress ) progress->step();
      
      refedge_ptr = refedge_array[i];
      if( refedge_ptr == NULL )
         continue;
      
        // There is not need to compare if the first RefEdge is
        // dactivated.
      if ( refedge_ptr->deactivated() == CUBIT_TRUE )
          continue ;

      if( refedge_ptr->get_compare_partner() )
      {
        //see if owning reface has compare partner...if so
        //edge is merged already.
        DLIList<RefFace*> tmp_faces;
        refedge_ptr->ref_faces( tmp_faces );

        int kk;
        bool owning_surface_has_compare_partner = false;
        for( kk=tmp_faces.size(); kk--; )
          if( tmp_faces.get_and_step()->get_compare_partner() )
            owning_surface_has_compare_partner = true;

        if( owning_surface_has_compare_partner )
          continue;
      }

      DLIList<RefEdge*> *new_list = NULL;
      
      for( j = i+1; (j < array_size) && !CubitMessage::instance()->Interrupt(); j++ )
      {
        compare_refedge_ptr = refedge_array[j];
         if( compare_refedge_ptr == NULL )
             continue;
         
           // Make sure we are not comparing the same entity
           // and that they are not deactivated.
         if( ( refedge_ptr == compare_refedge_ptr ) ||
             ( compare_refedge_ptr->deactivated() == CUBIT_TRUE ) )
             continue;
         
         double geom_factor = GeometryQueryTool::get_geometry_factor();
         CubitBoolean status =
             refedge_ptr->about_spatially_equal(compare_refedge_ptr,
                                                geom_factor,
                                                (CubitSense*)NULL,
                                                CUBIT_TRUE);
         if( status == CUBIT_FALSE )
             continue;

         
           //Make sure the children of the refedges are mergeable
         if( !refedge_ptr->children_mergeable() )
            continue;
         
           //Now let us test to see if we are merging two edges that
           //belong to the same face
         DLIList<RefFace*> ref_faces_1, ref_faces_2;
         refedge_ptr->ref_faces( ref_faces_1 );
         compare_refedge_ptr->ref_faces( ref_faces_2 );
         ref_faces_2.intersect( ref_faces_1 );
         if( ref_faces_2.size() > 0 )
         {
               PRINT_ERROR( "Tolerance problems, trying to find mergeable"
                            " curves\non the same surface.\n"
                            "%s (curve %d) and %s (curve %d) on\n"
                            "%s (surface %d)\n",
                            refedge_ptr->entity_name().c_str(),
                            refedge_ptr->id(),
                            compare_refedge_ptr->entity_name().c_str(),
                            compare_refedge_ptr->id(),
                            ref_faces_2.get()->entity_name().c_str(),
                            ref_faces_2.get()->id() );
               PRINT_INFO( "Try changing the merge tolerance.\n" );
               continue;
         }
           // If we are in this block, we want to merge the 
           // two entities. If they do not merge, there was
           // an error somewhere.
      
           // Always retain the entity with the lowest id.
         int nullify = j;
         if( refedge_ptr->id() > compare_refedge_ptr->id() )
         {
            std::swap(refedge_ptr, compare_refedge_ptr);
            nullify  = i;
         }

          // Now check if merge is okay with all assistants.
        CubitBoolean assistant_says_no = CUBIT_FALSE;
        for( int k = assistant_list_.size(); k > 0; k-- )
        {
          if( ! assistant_list_.get_and_step()
            ->can_merge( refedge_ptr, compare_refedge_ptr ) )
          {
            assistant_says_no = CUBIT_TRUE;
            break;
          }
        }
        if( assistant_says_no )
          continue; 

           // Need to retain these so that the pointers are not
           // accessed after a merge operation when the 'deleted'
           // pointer may be invalid.
         int retained_id = refedge_ptr->id();
         int deleted_id  = compare_refedge_ptr->id();
         
         if( new_list == NULL ) 
         {
           new_list = new DLIList<RefEdge*>; 
           new_list->append( refedge_ptr );
           new_list->append( compare_refedge_ptr );
         }
         else
           new_list->append( compare_refedge_ptr );
         
         refedge_array[nullify] = NULL;
         
         if (nullify == i)
           break;
      }

      if( new_list )
        lists_of_mergeable_ref_edges.append( new_list );
   }
   
  if( clean_up_compare_data )
    remove_compare_data();

  if(progress)
    progress->end();

   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_only_mergeable_surfaces ( DLIList<BodySM*> &body_list, 
                  DLIList< DLIList<Surface*>*> &lists_of_mergeable_surfaces )
{
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double merge_tolerance = geom_factor*GEOMETRY_RESABS;
  int i,j,k,l;

  std::map< Surface*, DLIList<Surface*>*> surf_to_list_map;
  std::map< Surface*, DLIList<Surface*>*>::iterator list_iter; 

  DLIList<BodySM*> tmp_body_list = body_list;
  tmp_body_list.reset();
  for(i=tmp_body_list.size(); i--; )
  {
    BodySM *tmp_body1 = tmp_body_list.pop();
    CubitBox body1_box = tmp_body1->bounding_box(); 
    for(j=tmp_body_list.size(); j--; )
    {
      BodySM *tmp_body2 = tmp_body_list.get_and_step();
      CubitBox body2_box = tmp_body2->bounding_box(); 
    
      if( body1_box.overlap( merge_tolerance, body2_box ) )
      {
        DLIList<Surface*> body1_surfs;
        DLIList<Surface*> body2_surfs;

        tmp_body1->surfaces( body1_surfs );
        tmp_body2->surfaces( body2_surfs );
        
        for( k=body1_surfs.size(); k--; )
        {
          Surface *body1_surf = body1_surfs.pop();
          CubitBox surf1_box = body1_surf->bounding_box();

          for(l=body2_surfs.size(); l--; )
          {
            Surface *body2_surf = body2_surfs.get();

            if( body2_surf == NULL )
            {
              body2_surfs.step();
              continue;
            }
  
            CubitBox surf2_box = body2_surf->bounding_box();
           
            if( surf1_box.overlap( merge_tolerance, surf2_box ) )
            {

              if( about_spatially_equal( body1_surf, body2_surf, geom_factor) )
              {
                //check to see if surfs have already been inserted into lists
                DLIList<Surface*> *surf1_list = NULL;
                DLIList<Surface*> *surf2_list = NULL;

                list_iter = surf_to_list_map.find( body1_surf);
                if( list_iter != surf_to_list_map.end() )
                  surf1_list = list_iter->second;

                list_iter = surf_to_list_map.find( body2_surf);
                if( list_iter != surf_to_list_map.end() )
                  surf2_list = list_iter->second;

                if( surf1_list == NULL && surf2_list == NULL )
                {
                  surf1_list = new DLIList<Surface*>;
                  surf1_list->append( body1_surf );
                  surf1_list->append( body2_surf );
                  surf_to_list_map.insert( std::map<Surface*,
                        DLIList<Surface*>*>::value_type( body1_surf, surf1_list )); 
                  surf_to_list_map.insert( std::map<Surface*,
                        DLIList<Surface*>*>::value_type( body2_surf, surf1_list )); 
                  lists_of_mergeable_surfaces.append( surf1_list );
                  break;
                }
                else if( surf1_list == NULL )
                {
                  PRINT_ERROR("A surface cannot be merged with more than 1 other surface\n");
                }
                else if( surf2_list == NULL )
                {
                  PRINT_ERROR("A surface cannot be merged with more than 1 other surface\n");
                }

                body2_surfs.change_to( NULL );
              }
            }
            body2_surfs.step();
          }
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_only_mergeable_curves( DLIList<BodySM*> &body_list, 
                  DLIList< DLIList<Curve*>*> &lists_of_mergeable_curves )
{
  //collect all the curves from off the bodies
  DLIList<Curve*> all_curves;
  body_list.reset();
  int i;
  for(i=body_list.size(); i--; )
  {
    BodySM* tmp_body = body_list.get_and_step();
    DLIList<Curve*> tmp_curves;
    tmp_body->curves( tmp_curves);
    all_curves += tmp_curves;
  }
 
  double geom_factor = GeometryQueryTool::get_geometry_factor();

  //build up a tree for speed purposes
  AbstractTree <Curve*> *a_tree = new RTree<Curve*> (GEOMETRY_RESABS*geom_factor);
  for( i=all_curves.size(); i--; )
    a_tree->add( all_curves.get_and_step() );

  std::map< Curve*, DLIList<Curve*>*> curve_to_list_map;
  std::map< Curve*, DLIList<Curve*>*>::iterator list_iter; 

  for( i=all_curves.size(); i--; )
  {
    Curve *curr_curve = all_curves.get_and_step();
    if( curr_curve == NULL )
      continue;
    
    //get close curves
    DLIList<Curve*> close_curves;
    a_tree->find(curr_curve->bounding_box(), close_curves);
    
    int j;
    for( j=close_curves.size(); j--; )
    {
      Curve *other_curve = close_curves.get_and_step();
      if( curr_curve == other_curve )
        continue;
  
      CubitSense rel_sense;
      CubitBoolean abs = about_spatially_equal( curr_curve, other_curve, rel_sense,
                                             geom_factor ); 
      if( abs )
      {
        //check to see if curves have already been inserted into lists
        DLIList<Curve*> *curve1_list = NULL;

        list_iter = curve_to_list_map.find( curr_curve );
        if( list_iter != curve_to_list_map.end() )
          curve1_list = list_iter->second;

        if( curve1_list == NULL ) 
        {
          curve1_list = new DLIList<Curve*>;
          curve1_list->append( curr_curve );
          curve1_list->append( other_curve );
          curve_to_list_map.insert( std::map<Curve*,
                DLIList<Curve*>*>::value_type( curr_curve , curve1_list )); 
          curve_to_list_map.insert( std::map<Curve*,
                DLIList<Curve*>*>::value_type( other_curve, curve1_list )); 
          lists_of_mergeable_curves.append( curve1_list );
        }
        else 
        {
          curve1_list->append( other_curve );
          curve_to_list_map.insert( std::map<Curve*,
                DLIList<Curve*>*>::value_type( other_curve, curve1_list )); 
        }
      
        //remove mergeable curves from list and reset list
        int item_index = all_curves.where_is_item( other_curve );
        if( item_index > 0 ) 
        {
          int curr_index = all_curves.get_index();
          all_curves.reset();
          all_curves.step( item_index );
          all_curves.change_to( NULL );
          all_curves.reset();
          all_curves.step( curr_index );
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus MergeTool::find_mergeable_refvertices( DLIList<RefEntity*> &entities,
                                                DLIList< DLIList<RefVertex*>*> &lists_of_mergeable_ref_vertices,
                                                bool clean_up_compare_data )
{
  
  DLIList<TopologyEntity*> t_ents;
  CAST_LIST(entities, t_ents, TopologyEntity);

  DLIList<RefVertex*> refvertex_list;
  int i;
  for( i=t_ents.size(); i--; )
  {
    DLIList<RefVertex*> tmp_verts;
    t_ents.get_and_step()->ref_vertices( tmp_verts);
    refvertex_list += tmp_verts;
  }

  if( refvertex_list.size() == 0 )
    return CUBIT_SUCCESS;

   DRefVertexArray refvertex_array( refvertex_list.size() );
   refvertex_list.reset();
   
     // Remove entities that should not be automatically merged
   for( i = refvertex_list.size(); i > 0; i-- )
   {
      RefVertex *curr_vert = refvertex_list.get_and_step();
      if( curr_vert->is_mergeable() )
          refvertex_array.append( curr_vert );
   }
  
   ProgressTool* progress = 0;

     // Now find overlapping RefVertices and merge them.
     // Make sure that the operation is not performed on
     // RefVertices that are deactivated.
   int j = 0;
   RefVertex* refvertex_ptr = NULL;
   RefVertex* compare_refvertex_ptr = NULL;
   int array_size = refvertex_array.size();
   for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
   {
      if( progress ) progress->step();
      
      refvertex_ptr = refvertex_array[i];
      if( refvertex_ptr == NULL )
          continue;

      if( refvertex_ptr->get_compare_partner() )
        continue;

        // There is not need to compare if the first RefVertex is
        // dactivated.
      if( refvertex_ptr->deactivated() == CUBIT_TRUE )
          continue ;
      
      DLIList<RefVertex*> *new_list = NULL;

      for( j = i+1; (j < array_size) && !CubitMessage::instance()->Interrupt(); j++ )
      {
         compare_refvertex_ptr = refvertex_array[j];
         if( compare_refvertex_ptr == NULL )
             continue;
         
           // Make sure we are not comparing the same entity
           // and that they are not deactivated.
         if( ( refvertex_ptr == compare_refvertex_ptr ) ||
             ( compare_refvertex_ptr->deactivated() == CUBIT_TRUE ) )
             continue;
         
         double geom_factor = GeometryQueryTool::get_geometry_factor();
         CubitBoolean status =
             refvertex_ptr->about_spatially_equal( compare_refvertex_ptr,
                                                   geom_factor );
         if( status == CUBIT_FALSE )
             continue;
         
           //Make sure we arn't merging two vertices on a
           //curve.
         DLIList<RefEdge*> edges_1, edges_2;
         refvertex_ptr->ref_edges( edges_1 );
         compare_refvertex_ptr->ref_edges( edges_2 );
         edges_2.intersect( edges_1 );
         if( edges_2.size() > 0 )
         {
               PRINT_ERROR( "Tolerance problems, trying to merge"
                            " vertices\non the same curve.\n"
                            "%s (vertex %d) and %s (vertex %d) on\n"
                            "%s (curve %d)\n",
                            refvertex_ptr->entity_name().c_str(),
                            refvertex_ptr->id(),
                            compare_refvertex_ptr->entity_name().c_str(),
                            compare_refvertex_ptr->id(),
                            edges_2.get()->entity_name().c_str(),
                            edges_2.get()->id() );
               PRINT_INFO( "Try changing the merge tolerance.\n" );
               continue;
         }

           // Always retain the entity with the lowest id.
         int nullify = j;
         if( refvertex_ptr->id() > compare_refvertex_ptr->id() )
         {
            std::swap(refvertex_ptr, compare_refvertex_ptr);
            nullify  = i;
         }

          // Now check if merge is okay with all assistants.
        CubitBoolean assistant_says_no = CUBIT_FALSE;
        for( int k = assistant_list_.size(); k > 0; k-- )
        {
          if( ! assistant_list_.get_and_step()
            ->can_merge( refvertex_ptr, compare_refvertex_ptr ) )
          {
            assistant_says_no = CUBIT_TRUE;
            break;
          }
        }
        if( assistant_says_no )
          continue; 

           // Need to retain these so that the pointers are not
           // accessed after a merge operation when the 'deleted'
           // pointer may be invalid.
         int retained_id = refvertex_ptr->id();
         int deleted_id  = compare_refvertex_ptr->id();

         refvertex_array[nullify] = NULL;

         if( new_list == NULL ) 
         {
           new_list = new DLIList<RefVertex*>; 
           new_list->append( refvertex_ptr );
           new_list->append( compare_refvertex_ptr );
         }
         else
           new_list->append( compare_refvertex_ptr );

         if( nullify == i )
             break;
      }

      if( new_list )
        lists_of_mergeable_ref_vertices.append( new_list );
   }
   
   if(progress)
     progress->end();

   if( clean_up_compare_data )
     remove_compare_data();

   if( CubitMessage::instance()->Interrupt() )
   {
     PRINT_WARNING("Vertex merging aborted.\n");
     return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
}


CubitStatus MergeTool::merge_refedges( DLIList<RefEdge*>& refedge_list,
                                       CubitBoolean should_clean_out,
                                       CubitBoolean print_info)
{
  if( refedge_list.size() < 20 )
    return old_merge_refedges( refedge_list, should_clean_out, print_info );

  set_merge_occurance( CUBIT_TRUE );
  int merge_count = 0;
  CpuTimer timer;

  double geom_factor = GeometryQueryTool::get_geometry_factor();
  AbstractTree <RefEdge*> *a_tree = new RTree<RefEdge*> (GEOMETRY_RESABS*geom_factor);

  DRefEdgeArray refedge_array( refedge_list.size() );
  refedge_list.reset();
  
    // Remove entities that should not be automatically merged
  int i = 0;
  int loop_size = refedge_list.size();
  CpuTimer time_to_build;
  for( i = 0; i < loop_size; i++ )
  {
    RefEdge *curr_edge = refedge_list.get_and_step();
    if( curr_edge->is_mergeable() )
    {
      refedge_array.append( curr_edge );
      a_tree->add(curr_edge);
    }
  }

    //initialize the marked flag for fast nulification...
  int array_size = refedge_array.size();
  for ( i = 0; i < array_size; i++)
    refedge_array[i]->marked(i);

  ProgressTool* progress = 0;
  if( displayProgress )
  {
    char info[64];
    sprintf(info, "Comparing %d Curves for Merge", refedge_array.size() );
    progress = AppUtil::instance()->progress_tool();
    if (progress)
    {
      progress->start( 0, refedge_array.size(), "Comparing Curves:", 
                       info, CUBIT_TRUE, CUBIT_TRUE );
    }
  }

    // Now find overlapping RefEdges and merge them.
    // Make sure that the operation is not performed on
    // RefEdges that are deactivated.
  int j = 0;
  RefEdge* refedge_ptr = NULL;
  RefEdge* compare_refedge_ptr = NULL;
  DLIList<RefEdge*> edges_merged, edges_in_range;
  CubitBox temp_box;
  for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
  {
    if( progress ) progress->step();
    
    refedge_ptr = refedge_array[i];
    if( refedge_ptr == NULL )
      continue;
    
      // There is no need to compare if the first RefEdge is
      // deactivated.
    if( refedge_ptr->deactivated() == CUBIT_TRUE )
      continue ;

      //Now from the atree, get the edges that are close.
    temp_box = refedge_ptr->bounding_box();
    edges_in_range.clean_out();
    a_tree->find(temp_box, edges_in_range);
    
    for( j = 0; (j < edges_in_range.size()) && !CubitMessage::instance()->Interrupt(); j++ )
    {
      compare_refedge_ptr = edges_in_range.get_and_step();
      if( compare_refedge_ptr == NULL )
        continue;
      
        // Make sure we are not comparing the same entity
        // and that they are not deactivated.
      if( ( refedge_ptr == compare_refedge_ptr ) ||
          ( compare_refedge_ptr->deactivated() == CUBIT_TRUE ) )
        continue;
      
        // IMPORTANT: this compare is for merging, so we set the
        //            last parameter to CUBIT_TRUE so that the
        //            compare sets the TDCompare data which will
        //            be used in the RefEntity compare and merge.
        // By Jihong
      
      geom_factor = GeometryQueryTool::get_geometry_factor();
      CubitBoolean status =
        refedge_ptr->about_spatially_equal( compare_refedge_ptr,
                                            geom_factor,
                                            (CubitSense*)NULL,
                                            CUBIT_TRUE );
      if( status == CUBIT_FALSE )
        continue;
      
        // If we are in this block, we want to merge the 
        // two entities. If they do not merge, there was
        // an error somewhere.
      
        //First test to see if all the children of these
        // refedges are mergeable
      if( !compare_refedge_ptr->children_mergeable() )
      {
        PRINT_WARNING( "Cannot merge curve %d and %d\n"
                       "     Make sure all merge flags are on.\n",
                       refedge_ptr->id(), compare_refedge_ptr->id() );
        continue;
      }
      
        // Always retain the entity with the lowest id.
      int nullify = compare_refedge_ptr->marked();
      if( refedge_ptr->id() > compare_refedge_ptr->id() )
      {
        std::swap(refedge_ptr, compare_refedge_ptr);
        nullify = i;
      }
      if (groupResults )
      {
        edges_merged.append(refedge_ptr);
      }

        // Now check if merge is okay with all assistants.
      CubitBoolean assistant_says_no = CUBIT_FALSE;
      for( int k = assistant_list_.size(); k > 0; k-- )
      {
        if( ! assistant_list_.get_and_step()
          ->can_merge( refedge_ptr, compare_refedge_ptr ) )
        {
          assistant_says_no = CUBIT_TRUE;
          break;
        }
      }
      if( assistant_says_no )
        continue; 

        // Need to retain these so that the pointers are not
        // accessed after a merge operation when the 'deleted'
        // pointer may be invalid.
      int retained_id = refedge_ptr->id();
      int deleted_id  = compare_refedge_ptr->id();
      
      PRINT_DEBUG_19( "Consolidating RefEdge %d and "
                      "%d...\n", retained_id, deleted_id );
      
      a_tree->remove(compare_refedge_ptr);
      if( merge_BTE( refedge_ptr, compare_refedge_ptr ) )
      {
        merge_count++;
        if (print_info && !progress) 
          PRINT_INFO( "Curve %d and %d consolidated\n",
                      retained_id, deleted_id );
        
          // The 'deleted' RefEdge is now gone. It is an
          // error to access that pointer, so we need to
          // get it out of the list.
        refedge_array[nullify] = NULL;
      }
      else
      {
        PRINT_ERROR( "Curve %d and %d NOT consolidated\n"
                     "       The use counts were updated first and"
                     " are now out of date.\n",
                     retained_id, deleted_id);
        a_tree->add(compare_refedge_ptr);
      }
      
      if (nullify == i)
        break;
      
    }
  }

    //clean the marks.
  for( i = 0; i < array_size; i++ )
  {
    if ( refedge_array[i] != NULL )
      refedge_array[i]->marked(0);
  }
  complete_merge();
  if (progress)
      progress->end();
  CpuTimer time_to_destroy;
  delete a_tree;
  PRINT_DEBUG_3( "Time to destroy r_tree %f secs.\n",
                  time_to_destroy.cpu_secs());

      // Remove the crud that accumulated during the merge
    // operations, from the geometry database.
  PRINT_DEBUG_3( "Merge RefEdge time: %f secs\n",
                 timer.cpu_secs() );
  
  PRINT_DEBUG_19("Cleaning out TDCompare data from RefEdges...\n");
  if ( groupResults && edges_merged.size() )
  {
    DLIList<RefEntity*> refentity_list;
    RefEdge *tmp_edge;
    for (int iii = edges_merged.size(); iii > 0; iii-- )
    {
      tmp_edge = edges_merged.get_and_step();
      if ( !tmp_edge->deactivated() )
        refentity_list.append(tmp_edge);
    }
    RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_curvs_merged");
    new_group->add_ref_entity( refentity_list );
    if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                               "  Group contains curves that were merged during\n"
                               "  current merge operation\n",
                               new_group->entity_name().c_str(),
                               new_group->id());
    lastCurvsMerged = new_group;
  }
  else
    lastCurvsMerged = NULL;
  
  if( destroyDeadGeometry )
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  
  //for updating GUI tree
  CubitObserver::notify_static_observers(NULL, MERGE_COMPLETED);
  
  PRINT_DEBUG_3( "cleanout time: %f secs\n", timer.cpu_secs() );
  if (print_info)
  {
    PRINT_INFO( "Consolidated %d pair", merge_count);
    if(merge_count > 1)
       PRINT_INFO("s");
    PRINT_INFO( " of curves \n");
  }
  if( CubitMessage::instance()->Interrupt() )
  {
     PRINT_WARNING("Curve merging aborted.\n");
     return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}


CubitStatus MergeTool::old_merge_refedges( DLIList<RefEdge*>& refedge_list,
                                       CubitBoolean should_clean_out,
                                       CubitBoolean print_info)
{
     // The basic algorithm is to step through the input list,
     // compare every entity with every other entity, and
     // merge the entities that are spatially equivalent.
     // In the merge, the entity with the lowest ID is retained.
     // NOTE: RefEdges can participate in multiple merges.
   set_merge_occurance( CUBIT_TRUE );
   int merge_count = 0;
   CpuTimer timer;
   DRefEdgeArray refedge_array( refedge_list.size() );
   refedge_list.reset();
   
     // Remove entities that should not be automatically merged
   int i = 0;
   for( i = refedge_list.size(); i > 0; i-- )
   {
      RefEdge *curr_edge = refedge_list.get_and_step();
      if( curr_edge->is_mergeable() )
          refedge_array.append( curr_edge );
   }
   DLIList<RefEdge*> edges_merged;
  
  
  ProgressTool* progress = 0;
  //if( displayProgress )
  //{
  //  progress = AppUtil::instance()->progress_tool();
  //  if (progress)
  //  {
  //    progress->start( 0, refedge_array.size(), "Comparing Curves:", 
  //                     0, CUBIT_TRUE, CUBIT_TRUE );
  //  }
  //}
   
     // Now find overlapping Refedges and merge them.
     // Make sure that the operation is not performed
     // on Refedges that are deactivated.
   int j = 0;
   RefEdge* refedge_ptr = NULL;
   RefEdge* compare_refedge_ptr = NULL;
   int array_size = refedge_array.size();
   for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
   {
      if( progress ) progress->step();
      
      refedge_ptr = refedge_array[i];
      if( refedge_ptr == NULL )
         continue;
      
        // There is not need to compare if the first RefEdge is
        // dactivated.
      if ( refedge_ptr->deactivated() == CUBIT_TRUE )
          continue ;
      
//         // Get the GeometryQueryEngine of the Refedges
//         // Try to merge this RefEdge only if it is associated
//         // with a SolidModelingEngine.
//       GeometryQueryEngine* firstGMEPtr = 
//           refedge_ptr->get_geometry_query_engine();
//       SolidModelingEngine* SMEPtr = CAST_TO( firstGMEPtr, SolidModelingEngine );
      
//         // If the RefEdge is not associated with a
//         // SolidModelingEngine, go on to the next RefEdge.
//       if ( SMEPtr == NULL )
//           continue ;
      
      for( j = i+1; (j < array_size) && !CubitMessage::instance()->Interrupt(); j++ )
      {
        compare_refedge_ptr = refedge_array[j];
         if( compare_refedge_ptr == NULL )
             continue;
         
           // Make sure we are not comparing the same entity
           // and that they are not deactivated.
         if( ( refedge_ptr == compare_refedge_ptr ) ||
             ( compare_refedge_ptr->deactivated() == CUBIT_TRUE ) )
             continue;
         
//            // Get the GeometryQueryEngine of the second RefEdge.
//            // Make sure that both engines are same before proceeding
//            // with the merge.
//          GeometryQueryEngine* secondGMEPtr = 
//              compare_refedge_ptr->get_geometry_query_engine() ;
         
//            // If the two engines are different, move on to the 
//            // next RefEdge.
//          if( firstGMEPtr != secondGMEPtr ) 
//              continue;
         
           // IMPORTANT: this compare is for merging, so we set the
           //            last parameter to CUBIT_TRUE so that the
           //            compare sets the TDCompare data which will
           //            be used in the RefEntity compare and merge.
           // By Jihong 
         
         double geom_factor = GeometryQueryTool::get_geometry_factor();
         CubitBoolean status =
             refedge_ptr->about_spatially_equal(compare_refedge_ptr,
                                                geom_factor,
                                                (CubitSense*)NULL,
                                                CUBIT_TRUE );
         if( status == CUBIT_FALSE )
             continue;
         
           //Make sure the children of the refedges are mergeable
         if( !refedge_ptr->children_mergeable() )
         {
            PRINT_WARNING( "Cannot merge curves %d and %d\n"
                           "     Make sure merge flags are on\n",
                           refedge_ptr->id(), compare_refedge_ptr->id() );
            continue;
         }
         
           //Now let us test to see if we are merging two edges that
           //belong to the same face
         DLIList<RefFace*> ref_faces_1, ref_faces_2;
         refedge_ptr->ref_faces( ref_faces_1 );
         compare_refedge_ptr->ref_faces( ref_faces_2 );
         ref_faces_2.intersect( ref_faces_1 );
         if( ref_faces_2.size() > 0 )
         {
               PRINT_ERROR( "Tolerance problems, trying to merge"
                            " curves\non the same surface.\n"
                            "%s (curve %d) and %s (curve %d) on\n"
                            "%s (surface %d)\n",
                            refedge_ptr->entity_name().c_str(),
                            refedge_ptr->id(),
                            compare_refedge_ptr->entity_name().c_str(),
                            compare_refedge_ptr->id(),
                            ref_faces_2.get()->entity_name().c_str(),
                            ref_faces_2.get()->id() );
               PRINT_INFO( "Try changing the merge tolerance.\n" );
               continue;
         }
           // If we are in this block, we want to merge the 
           // two entities. If they do not merge, there was
           // an error somewhere.
      
           // Always retain the entity with the lowest id.
         int nullify = j;
         if( refedge_ptr->id() > compare_refedge_ptr->id() )
         {
            std::swap(refedge_ptr, compare_refedge_ptr);
            nullify  = i;
         }

          // Now check if merge is okay with all assistants.
        CubitBoolean assistant_says_no = CUBIT_FALSE;
        for( int k = assistant_list_.size(); k > 0; k-- )
        {
          if( ! assistant_list_.get_and_step()
            ->can_merge( refedge_ptr, compare_refedge_ptr ) )
          {
            assistant_says_no = CUBIT_TRUE;
            break;
          }
        }
        if( assistant_says_no )
          continue; 

           // Need to retain these so that the pointers are not
           // accessed after a merge operation when the 'deleted'
           // pointer may be invalid.
         int retained_id = refedge_ptr->id();
         int deleted_id  = compare_refedge_ptr->id();
         if ( groupResults )
         {
           edges_merged.append(refedge_ptr);
         }
         if( merge_BTE( refedge_ptr, compare_refedge_ptr ) )
         {
            merge_count++;
            if (print_info && !progress) 
              PRINT_INFO( "Curve %d and %d consolidated\n",
                                        retained_id, deleted_id);
            
              // The 'deleted' RefEdge is now gone. It is an
              // error to access that pointer, so we need to
              // get it out of the list.
            refedge_array[nullify] = NULL;
         }
         else
         {
            PRINT_ERROR( "Curve %d and %d NOT consolidated\n"
                         "       The use counts were updated first and"
                         " are now out of date.\n",
                         retained_id, deleted_id);
         }
         
         if (nullify == i)
           break;
      }
   }
   
     // Remove the crud that accumulated during the
     // merge operations, from the geometry database.
   PRINT_DEBUG_3( "Merge RefEdges time: %f secs\n",
                timer.cpu_secs() );
   
   PRINT_DEBUG_19(
                "Cleaning out TDCompare data from RefEdges...\n" );
   complete_merge();
   if (progress)
       progress->end();
   
   if ( groupResults && edges_merged.size()  )
   {
     DLIList<RefEntity*> refentity_list;
     RefEdge *tmp_edge;
     for (int iii = edges_merged.size(); iii > 0; iii-- )
     {
       tmp_edge = edges_merged.get_and_step();
       if ( !tmp_edge->deactivated() )
         refentity_list.append(tmp_edge);
     }
     RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_curvs_merged");
     new_group->add_ref_entity( refentity_list );
     if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                                "  Group contains curves that were seperatly merged during\n"
                                "  current merge operation (ie, not during surface merge)\n",
                                new_group->entity_name().c_str(),
                                new_group->id());
     lastCurvsMerged = new_group;
   }
     //set this to null otherwise.
   else
     lastCurvsMerged = NULL;
   
   if( should_clean_out == CUBIT_TRUE && destroyDeadGeometry )
   {
      GeometryQueryTool::instance()->cleanout_deactivated_geometry();
      PRINT_DEBUG_3( "cleanout time: %f secs\n",
                   timer.cpu_secs() );
   }
   if (print_info) PRINT_INFO( "Consolidated %d curves\n", merge_count );
   if( CubitMessage::instance()->Interrupt() )
   {
     PRINT_WARNING("Curve merging aborted.\n");
     return CUBIT_FAILURE;
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus MergeTool::merge_all_refvertices()
{
//   if( !displayProgress )
     PRINT_INFO( "\n...Merging all Vertices in the model\n" );
   
   set_merge_occurance( CUBIT_TRUE );
   
     // Get the list of all the RefVertices in the model
   DLIList<RefVertex*> refvertex_list_model;
   GeometryQueryTool::instance()->ref_vertices(refvertex_list_model);
   
     // Merge the RefVerties in this list
   if( merge_refvertices( refvertex_list_model ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Vertex merging failed\n" );
      return CUBIT_FAILURE;
   }
   
   return CUBIT_SUCCESS;
}


CubitStatus MergeTool::merge_refvertices( DLIList<RefVertex*>& refvertex_list,
                                          CubitBoolean print_info)
{
  if( refvertex_list.size() < 20 )
    return old_merge_refvertices( refvertex_list, print_info );

  set_merge_occurance( CUBIT_TRUE );
   
  CpuTimer timer;
  int merge_count = 0;
  
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double tol = GEOMETRY_RESABS*geom_factor;
  AbstractTree <RefVertex*> *a_tree = new RTree<RefVertex*> ( tol );
  DRefVertexArray refvertex_array( refvertex_list.size() );
  refvertex_list.reset();
   
    // Remove entities that should not be automatically merged
  int i = 0;
  for( i = refvertex_list.size(); i > 0; i-- )
  {
    RefVertex *curr_vert = refvertex_list.get_and_step();
    if( curr_vert->is_mergeable() )
    {
      refvertex_array.append( curr_vert );
      a_tree->add( curr_vert );
    }
  }
  
    //initialize the marked flag for fast nulification...
  int array_size = refvertex_array.size();
  for ( i = 0; i < array_size; i++)
    refvertex_array[i]->marked(i);

  DLIList<RefVertex*> vertices_merged;
  
  ProgressTool* progress = 0;

   // Now find overlapping RefVertices and merge them.
   // Make sure that the operation is not performed on
   // RefVertices that are deactivated.
  int j = 0;
  RefVertex* refvertex_ptr = NULL;
  RefVertex* compare_refvertex_ptr = NULL;
  for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
  {
    if( progress ) progress->step();
      
    refvertex_ptr = refvertex_array[i];
    if( refvertex_ptr == NULL )
      continue;

    // There is not need to compare if the first RefVertex is
    // dactivated.
    if( refvertex_ptr->deactivated() == CUBIT_TRUE )
      continue ;
     
    DLIList<RefVertex*> close_verts;
    a_tree->find(refvertex_ptr->bounding_box(), close_verts);
      
    for( j = 0; (j < close_verts.size()) && !CubitMessage::instance()->Interrupt(); j++ )
    {
      compare_refvertex_ptr = close_verts.get_and_step(); 
     
      //skip vertices already handled
      if( refvertex_array[ compare_refvertex_ptr->marked() ] == NULL )
        continue;

      if( compare_refvertex_ptr == NULL )
        continue;

      // Make sure we are not comparing the same entity
      // and that they are not deactivated.
      if( ( refvertex_ptr == compare_refvertex_ptr ) ||
        ( compare_refvertex_ptr->deactivated() == CUBIT_TRUE ) )
        continue;
         
      // IMPORTANT: this compare is for merging, so we set the
      //            last parameter to CUBIT_TRUE so that the
      //            compare sets the TDCompare data which will
      //            be used in the RefEntity compare and merge.
      // By Jihong 
      double geom_factor = GeometryQueryTool::get_geometry_factor();
      CubitBoolean status =
        refvertex_ptr->about_spatially_equal( compare_refvertex_ptr,
                                                   geom_factor,
                                                   CUBIT_TRUE );
      if( status == CUBIT_FALSE )
        continue;
         
      //Make sure we arn't merging two vertices on a
      //curve.
      DLIList<RefEdge*> edges_1, edges_2;
      refvertex_ptr->ref_edges( edges_1 );
      compare_refvertex_ptr->ref_edges( edges_2 );
      edges_2.intersect( edges_1 );
      if( edges_2.size() > 0 )
      {
        PRINT_ERROR( "Tolerance problems, trying to merge"
                     " vertices\non the same curve.\n"
                     "%s (vertex %d) and %s (vertex %d) on\n"
                     "%s (curve %d)\n",
                     refvertex_ptr->entity_name().c_str(),
                     refvertex_ptr->id(),
                     compare_refvertex_ptr->entity_name().c_str(),
                     compare_refvertex_ptr->id(),
                     edges_2.get()->entity_name().c_str(),
                     edges_2.get()->id() );
        PRINT_INFO( "Try changing the merge tolerance.\n" );
        continue;
      }

      // Always retain the entity with the lowest id.
      int nullify = compare_refvertex_ptr->marked();
      if( refvertex_ptr->id() > compare_refvertex_ptr->id() )
      {
        std::swap(refvertex_ptr, compare_refvertex_ptr);
        nullify  = i;
      }

      // Now check if merge is okay with all assistants.
      CubitBoolean assistant_says_no = CUBIT_FALSE;
      for( int k = assistant_list_.size(); k > 0; k-- )
      {
        if( ! assistant_list_.get_and_step()
          ->can_merge( refvertex_ptr, compare_refvertex_ptr ) )
        {
          assistant_says_no = CUBIT_TRUE;
          break;
        }
      }
      if( assistant_says_no )
        continue; 

      // Need to retain these so that the pointers are not
      // accessed after a merge operation when the 'deleted'
      // pointer may be invalid.
      int retained_id = refvertex_ptr->id();
      int deleted_id  = compare_refvertex_ptr->id();
      if (groupResults)
      {
        vertices_merged.append(refvertex_ptr);
      }

      a_tree->remove( compare_refvertex_ptr );
      if( merge_BTE( refvertex_ptr, compare_refvertex_ptr ) )
      {
        merge_count++;
        if(print_info && !progress) 
          PRINT_INFO( "Vertex %d and %d consolidated\n",
                       retained_id, deleted_id);
            
        // The 'deleted' RefVertex is now gone. It is an
        // error to access that pointer, so we need to
        // get it out of the list.
        refvertex_array[nullify] = NULL;
      }
      else
      {
        PRINT_ERROR( "Vertex %d and %d NOT consolidated\n"
                     "       The use counts were updated first and"
                     " are now out of date.\n",
                     retained_id, deleted_id );
        a_tree->add( compare_refvertex_ptr );
      }
      if( nullify == i )
      break;
    }
  }
   
  //clean the marks.
  for( i = 0; i < array_size; i++ )
  {
    if ( refvertex_array[i] != NULL )
      refvertex_array[i]->marked(0);
  }

     // Remove the crud that accumulated during the merge
     // operations, from the geometry database.
  PRINT_DEBUG_3( "Merge RefVertexs time: %f secs\n",
                timer.cpu_secs() );
   
  PRINT_DEBUG_19( "Cleaning out TDCompare data from RefVertices...\n");
  complete_merge();
  if(progress)
    progress->end();
  delete a_tree;
  if ( groupResults && vertices_merged.size() )
  {
    DLIList<RefEntity*> refentity_list;
    RefVertex *tmp_vertex;
    for (int iii = vertices_merged.size(); iii > 0; iii-- )
    {
      tmp_vertex = vertices_merged.get_and_step();
      if ( !tmp_vertex->deactivated() )
        refentity_list.append(tmp_vertex);
    }
    RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_verts_merged");
    new_group->add_ref_entity( refentity_list );
    if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                               "  Group contains curves that were seperatly merged during\n"
                               "  current merge operation (ie, not during surface merge)\n",
                               new_group->entity_name().c_str(),
                               new_group->id());
    lastVertsMerged = new_group;
  }   
     //set this to null otherwise.
  else
    lastVertsMerged = NULL;

  if( destroyDeadGeometry )
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  PRINT_DEBUG_3( "cleanout time: %f secs\n",
                timer.cpu_secs() );

  //for updating GUI tree
  CubitObserver::notify_static_observers(NULL, MERGE_COMPLETED);

  if (print_info) PRINT_INFO( "Consolidated %d pairs of vertices\n", merge_count );
  if( CubitMessage::instance()->Interrupt() )
  {
    PRINT_WARNING("Vertex merging aborted.\n");
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Purpose       : merge all items in the RefVertex list if
//                 they can be merged.
//
// Special Notes : 
//
// Creator       : Jihong Ma
//
// Creation Date : 11/27/96
//----------------------------------------------------------------
CubitStatus MergeTool::old_merge_refvertices( DLIList<RefVertex*>& refvertex_list,
                                          CubitBoolean print_info)
{
     // The basic algorithm is to step through the input list,
     // compare every entity with every other entity, and
     // merge the entities that are spatially equivalent.
     // In the merge, the entity with the lowest ID is retained.
     // NOTE: RefVertices can participate in multiple merges.
   set_merge_occurance( CUBIT_TRUE );
   
   CpuTimer timer;
   int merge_count = 0;
   
   DRefVertexArray refvertex_array( refvertex_list.size() );
   refvertex_list.reset();
   
     // Remove entities that should not be automatically merged
   int i = 0;
   for( i = refvertex_list.size(); i > 0; i-- )
   {
      RefVertex *curr_vert = refvertex_list.get_and_step();
      if( curr_vert->is_mergeable() )
          refvertex_array.append( curr_vert );
   }
   DLIList<RefVertex*> vertices_merged;
  
  ProgressTool* progress = 0;
//  if( displayProgress )
//  {
//    progress = AppUtil::instance()->progress_tool();
//    if (progress)
//    {
//      progress->start( 0, refvertex_array.size(), 
//                       "Comparing Vertices:", 0, CUBIT_TRUE, CUBIT_TRUE );
//    }
//  }

     // Now find overlapping RefVertices and merge them.
     // Make sure that the operation is not performed on
     // RefVertices that are deactivated.
   int j = 0;
   RefVertex* refvertex_ptr = NULL;
   RefVertex* compare_refvertex_ptr = NULL;
   int array_size = refvertex_array.size();
   for( i = 0; (i < array_size) && !CubitMessage::instance()->Interrupt(); i++ )
   {
      if( progress ) progress->step();
      
      refvertex_ptr = refvertex_array[i];
      if( refvertex_ptr == NULL )
          continue;
      
        // There is not need to compare if the first RefVertex is
        // dactivated.
      if( refvertex_ptr->deactivated() == CUBIT_TRUE )
          continue ;
      
//         // Get the GeometryQueryEngine of the Refedges
//         // Try to merge this RefVertex only if it is associated
//         // with a SolidModelingEngine.
//       GeometryQueryEngine* firstGMEPtr = 
//           refvertex_ptr->get_geometry_query_engine();
//       SolidModelingEngine* SMEPtr = CAST_TO( firstGMEPtr, SolidModelingEngine );
      
//         // If the RefVertex is not associated with a
//         // SolidModelingEngine, go on to the next RefVertex.
//       if( SMEPtr == NULL )
//           continue ;
      
      for( j = i+1; (j < array_size) && !CubitMessage::instance()->Interrupt(); j++ )
      {
         compare_refvertex_ptr = refvertex_array[j];
         if( compare_refvertex_ptr == NULL )
             continue;
         
           // Make sure we are not comparing the same entity
           // and that they are not deactivated.
         if( ( refvertex_ptr == compare_refvertex_ptr ) ||
             ( compare_refvertex_ptr->deactivated() == CUBIT_TRUE ) )
             continue;
         
//            // Get the GeometryQueryEngine of the second
//            // RefVertex.  Make sure that both engines are same
//            // before proceeding with the merge.
//          GeometryQueryEngine* secondGMEPtr = 
//              compare_refvertex_ptr->get_geometry_query_engine() ;
         
//            // If the two engines are different, move on to the 
//            // next RefVertex.
//          if( firstGMEPtr != secondGMEPtr ) 
//              continue ; 
         
           // IMPORTANT: this compare is for merging, so we set the
           //            last parameter to CUBIT_TRUE so that the
           //            compare sets the TDCompare data which will
           //            be used in the RefEntity compare and merge.
           // By Jihong 
         double geom_factor = GeometryQueryTool::get_geometry_factor();
         CubitBoolean status =
             refvertex_ptr->about_spatially_equal( compare_refvertex_ptr,
                                                   geom_factor,
                                                   CUBIT_TRUE );
         if( status == CUBIT_FALSE )
             continue;
         
           //Make sure we arn't merging two vertices on a
           //curve.
         DLIList<RefEdge*> edges_1, edges_2;
         refvertex_ptr->ref_edges( edges_1 );
         compare_refvertex_ptr->ref_edges( edges_2 );
         edges_2.intersect( edges_1 );
         if( edges_2.size() > 0 )
         {
               PRINT_ERROR( "Tolerance problems, trying to merge"
                            " vertices\non the same curve.\n"
                            "%s (vertex %d) and %s (vertex %d) on\n"
                            "%s (curve %d)\n",
                            refvertex_ptr->entity_name().c_str(),
                            refvertex_ptr->id(),
                            compare_refvertex_ptr->entity_name().c_str(),
                            compare_refvertex_ptr->id(),
                            edges_2.get()->entity_name().c_str(),
                            edges_2.get()->id() );
               PRINT_INFO( "Try changing the merge tolerance.\n" );
               continue;
         }

           // Always retain the entity with the lowest id.
         int nullify = j;
         if( refvertex_ptr->id() > compare_refvertex_ptr->id() )
         {
            std::swap(refvertex_ptr, compare_refvertex_ptr);
            nullify  = i;
         }

          // Now check if merge is okay with all assistants.
        CubitBoolean assistant_says_no = CUBIT_FALSE;
        for( int k = assistant_list_.size(); k > 0; k-- )
        {
          if( ! assistant_list_.get_and_step()
            ->can_merge( refvertex_ptr, compare_refvertex_ptr ) )
          {
            assistant_says_no = CUBIT_TRUE;
            break;
          }
        }
        if( assistant_says_no )
          continue; 

           // Need to retain these so that the pointers are not
           // accessed after a merge operation when the 'deleted'
           // pointer may be invalid.
         int retained_id = refvertex_ptr->id();
         int deleted_id  = compare_refvertex_ptr->id();
         if (groupResults)
         {
           vertices_merged.append(refvertex_ptr);
         }
         if( merge_BTE( refvertex_ptr, compare_refvertex_ptr ) )
         {
            merge_count++;
            if (print_info && !progress) 
              PRINT_INFO( "Vertex %d and %d consolidated\n",
                                        retained_id, deleted_id);
            
              // The 'deleted' RefVertex is now gone. It is an
              // error to access that pointer, so we need to
              // get it out of the list.
            refvertex_array[nullify] = NULL;
         }
         else
         {
            PRINT_ERROR( "Vertex %d and %d NOT consolidated\n"
                         "       The use counts were updated first and"
                         " are now out of date.\n",
                         retained_id, deleted_id );
         }

         if( nullify == i )
             break;
      }
   }
   
     // Remove the crud that accumulated during the merge
     // operations, from the geometry database.
   PRINT_DEBUG_3( "Merge RefVertexs time: %f secs\n",
                timer.cpu_secs() );
   
   PRINT_DEBUG_19(
                "Cleaning out TDCompare data from RefVertices...\n");
   complete_merge();
   if (progress)
       progress->end();
   if ( groupResults && vertices_merged.size() )
   {
     DLIList<RefEntity*> refentity_list;
     RefVertex *tmp_vertex;
     for (int iii = vertices_merged.size(); iii > 0; iii-- )
     {
       tmp_vertex = vertices_merged.get_and_step();
       if ( !tmp_vertex->deactivated() )
         refentity_list.append(tmp_vertex);
     }
     RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup("gr_verts_merged");
     new_group->add_ref_entity( refentity_list );
     if (print_info) PRINT_INFO("Created new group %s (Group %d)\n"
                                "  Group contains curves that were seperatly merged during\n"
                                "  current merge operation (ie, not during surface merge)\n",
                                new_group->entity_name().c_str(),
                                new_group->id());
     lastVertsMerged = new_group;
   }   
     //set this to null otherwise.
   else
     lastVertsMerged = NULL;

   if( destroyDeadGeometry )
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
   PRINT_DEBUG_3( "cleanout time: %f secs\n",
                timer.cpu_secs() );
   
   if (print_info) PRINT_INFO( "Consolidated %d pairs of vertices\n", merge_count );
   if( CubitMessage::instance()->Interrupt() )
   {
     PRINT_WARNING("Vertex merging aborted.\n");
     return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
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
CubitStatus MergeTool::unmerge_all()
{
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
  if(result == CUBIT_SUCCESS)
     set_merge_occurance(CUBIT_FALSE);
  
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
CubitStatus MergeTool::unmerge( DLIList<RefEntity*> &entity_list,
                                CubitBoolean descend )
{
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
CubitStatus MergeTool::unmerge( RefEntity* entity_ptr, CubitBoolean descend )
{
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
                "MergeTool::unmerge(RefEntity*,CubitBoolean)\n",
                entity_ptr->class_name());
      return CUBIT_FAILURE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::unmerge( Body* body_ptr )
{
  DLIList<Body*> list(1);
  list.append( body_ptr );
  return separate_bodies( list );
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::unmerge( RefVolume* vol_ptr )
{
  DLIList<RefVolume*> list(1);
  list.append( vol_ptr );
  return separate_volumes( list );
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::unmerge( RefFace* face_ptr, CubitBoolean descend )
{
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  int i;
  
  DLIList<TopologyBridge*> bridge_list;
  face_ptr->bridge_manager()->get_bridge_list(bridge_list);
  if (bridge_list.size() < 2)
    return CUBIT_SUCCESS;
  
  DLIList<Surface*> surf_list;
  bridge_list.reset();
  for (i = bridge_list.size(); i > 1; i--)
  {
    Surface* surf = dynamic_cast<Surface*>(bridge_list.step_and_get());
    surf_list.clean_out();
    surf_list.append( surf );
    if (0 == separate_face( surf_list, descend ))
    {
      result = CUBIT_FAILURE;
      break;
    }
  }
       
  end_unmerge(top);
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::unmerge( RefEdge* edge_ptr, CubitBoolean descend )
{
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  int i;
  
  DLIList<TopologyBridge*> bridge_list;
  edge_ptr->bridge_manager()->get_bridge_list(bridge_list);
  if (bridge_list.size() < 2)
    return CUBIT_SUCCESS;
  
  DLIList<Curve*> curve_list;
  bridge_list.reset();
  for (i = bridge_list.size(); i > 1; i--)
  {
    Curve* curve = dynamic_cast<Curve*>(bridge_list.step_and_get());
    curve_list.clean_out();
    curve_list.append( curve );
    if (0 == separate_edge( curve_list, descend ))
    {
      result = CUBIT_FAILURE;
      break;
    }
  }
  
  end_unmerge(top);
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::unmerge( RefVertex* vtx_ptr )
{
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  int i;
  
  DLIList<TopologyBridge*> bridge_list;
  vtx_ptr->bridge_manager()->get_bridge_list(bridge_list);
  if (bridge_list.size() < 2)
    return CUBIT_SUCCESS;
  
  DLIList<Point*> point_list;
  bridge_list.reset();
  for (i = bridge_list.size(); i > 1; i--)
  {
    Point* point = dynamic_cast<Point*>(bridge_list.step_and_get());
    point_list.clean_out();
    point_list.append( point );
    if (0 == separate_vertex( point_list ))
    {
      result = CUBIT_FAILURE;
      break;
    }
  }

  end_unmerge(top);
  return result;
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
void MergeTool::cleanup_unmerge()
{
  std::set<CubitObservable*> modified_list, new_list;
  int i;
//   CpuTimer timer;
  int vtx_count = 0, curve_count = 0, surf_count = 0;
  
  DLIList<RefEntity*> parents;
  assert( new_unmerged.size() == old_unmerged.size() );
  new_unmerged.reset();
  old_unmerged.reset();
  for (i = new_unmerged.size(); i--; )
  {
    RefEntity* new_ptr = new_unmerged.get_and_step();
    RefEntity* old_ptr = old_unmerged.get_and_step();
    if (dynamic_cast<RefVertex*>(new_ptr))
      vtx_count++;
    else if(dynamic_cast<RefEdge*>(new_ptr))
      curve_count++;
    if (dynamic_cast<RefFace*>(new_ptr))
      surf_count++;
      

    new_ptr->get_parent_ref_entities( parents );
    
    if (parents.size() == 0)
    {
      if (new_list.insert(new_ptr).second)
        CubitObserver::notify_static_observers( new_ptr, FREE_REF_ENTITY_GENERATED );
    }
    
    while (parents.size())
      modified_list.insert( parents.pop() );

    old_ptr->get_parent_ref_entities( parents );
    while (parents.size())
      modified_list.insert( parents.pop() );
    modified_list.insert( old_ptr );
    
    UnMergeEvent event( old_ptr, new_ptr );
    old_ptr->notify_all_observers( event );
  }
  
  std::set<CubitObservable*>::iterator iter;
  for (iter = modified_list.begin(); iter != modified_list.end(); ++iter)
    (*iter)->notify_all_observers( TOPOLOGY_MODIFIED );

  for( int a = assistant_list_.size(); a--; )
    assistant_list_.get_and_step()->finish_unmerge();
    
  if (vtx_count)
    PRINT_INFO("Unmerged %d vertices.\n", vtx_count);
  if (curve_count)
    PRINT_INFO("Unmerged %d curves.\n", curve_count);
  if (surf_count)
    PRINT_INFO("Unmerged %d surfaces.\n", surf_count);
  if (vtx_count + curve_count + surf_count == 0)
    PRINT_INFO("No entities unmerged.\n");
}

    

void MergeTool::compare_notify(RefEntity *entity, CubitEventType event)
{
    //- notifies MergeTool about comparisons found and put on ref entities
  assert(event == COMPARISON_FOUND);
  compareEntityList.append_unique(entity);
}

CubitStatus MergeTool::merge_BTE( BasicTopologyEntity* keeper_entity,
                                  BasicTopologyEntity* dead_entity )
{
  CubitStatus result = CUBIT_FAILURE;
  DLIList<ModelEntity*> query_results, query_results_2;
  
    // Make sure that the 2 BTE's are of the same type
  if( keeper_entity->dag_type() != dead_entity->dag_type() )
  {
    PRINT_ERROR( "In MergeTool::merge_BTE(), \n"
                 "      two different types of entities "
                 "are merged. \n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n");
    assert(0);
  }
  
     // Debug info
  CubitString keeper_entity_name(""), dead_entity_name("");
  if (DEBUG_FLAG(19))
  {
    keeper_entity_name = keeper_entity->entity_name();
    dead_entity_name = dead_entity->entity_name();
  
    PRINT_DEBUG_19("...In MergeTool::merge_BTE - "
                  "Merging %s and %s\n",
                  keeper_entity_name.c_str(),
                  dead_entity_name.c_str() );
  }
  
    // Sanity Check: don't merge the entities if either of them
    // is already deactivated
  if( keeper_entity->deactivated() )
  {
    PRINT_ERROR( "%s has already been merged\n",
                 keeper_entity->entity_name().c_str() );
    return CUBIT_FAILURE;
  }
  if( dead_entity->deactivated() )
  {
    PRINT_ERROR( "%s has already been merged\n",
                 dead_entity->entity_name().c_str() );
    return CUBIT_FAILURE;
  }
  
    // sanity check - don't merge entity with itself
  if( keeper_entity == dead_entity )
  {
      // Debug info
    PRINT_DEBUG_19( "   Did not merge %s and %s. "
                    "Cannot merge entity with itself.\n",
                    keeper_entity_name.c_str(),
                    dead_entity_name.c_str() );
    return CUBIT_SUCCESS;
  }
  
    // sanity check - don't merge entities w/ same parent
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  DagType parent_type = keeper_entity->get_parent_ref_entity_type();
  mqe->query_model( *keeper_entity, parent_type, query_results );
  mqe->query_model( *dead_entity, parent_type, query_results_2 );
  query_results.intersect( query_results_2 );
  if (query_results.size())
  {
    PRINT_ERROR( "In MergeTool::merge_BTE()\n"
                 "  Attempt to merge two entities with same parent.\n"
                 "  %s (%s %d) and %s (%s %d)\n",
                 keeper_entity->entity_name().c_str(),
                 keeper_entity->class_name(), keeper_entity->id(),
                 dead_entity->entity_name().c_str(),
                 dead_entity->class_name(), dead_entity->id());
    return CUBIT_FAILURE;
  }

    // Merge the GroupingEntitys of the BasicTopologyEntity
  DLIList<GroupingEntity*> keeper_GE_list;
  DLIList<GroupingEntity*> dead_GE_list;
  
    // First get the GroupingEntities of the BTE's
  keeper_entity->get_grouping_entity_list( keeper_GE_list );
  dead_entity->get_grouping_entity_list( dead_GE_list );
  
    // Make sure they have the same number of GroupingEntities
  if( keeper_GE_list.size() != dead_GE_list.size() )
  {
    PRINT_ERROR( "In MergeTool::merge_BTE()\n"
                 "    the two entities have different "
                 "numbers of GroupingEntities.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0);
  }
  
    // Merge all child BTEs
  BasicTopologyEntity *bte_ptr_1, *bte_ptr_2;
  DagType child_type = keeper_entity->get_child_ref_entity_type();
  if (child_type.is_valid())
  {
    query_results.clean_out();
    mqe->query_model( *keeper_entity, child_type, query_results );
    while (query_results.size())
    {
      bte_ptr_1 = dynamic_cast<BasicTopologyEntity*>(query_results.pop());
      bte_ptr_2 = dynamic_cast<BasicTopologyEntity*>(bte_ptr_1->get_compare_partner());
      if (!bte_ptr_2 || bte_ptr_2->deactivated())
        continue;
      
      if (bte_ptr_2->get_compare_partner() != bte_ptr_1)  
      {
        PRINT_ERROR("Bad compare partner TDs encountered at %s:%d.\n"
                    "This is a bug.  Please report it.\n",
                    __FILE__, __LINE__ );
        return CUBIT_FAILURE;
      }
      
      CubitStatus merge_status;
      if (bte_ptr_1->id() < bte_ptr_2->id())
        merge_status = merge_BTE(bte_ptr_1, bte_ptr_2);
      else
        merge_status = merge_BTE(bte_ptr_2, bte_ptr_1);
        
      if (!merge_status)
        return CUBIT_FAILURE;
    }
  }
      
    // If RefFace or RefEdge, adjust sense of parent sense entities,
    // if necessary.  This was previously handled by the 
    // switch_child_notify() callback in DAGNode/ModelEntity.
    // However, with virtual geometry, there may be CoEdges/CoFaces
    // that do not get merged with anything and still need to be
    // updated.  Also, by taking this out of switch_child_notify,
    // compare_alignment() does not need to be called for every
    // SenseEntitiy merge.  We only need to call it once.
  CubitBoolean switch_sense = CUBIT_FALSE;
  if(CAST_TO( keeper_entity, RefFace ) )
  {
    RefFace* keep_face = CAST_TO(keeper_entity,RefFace);
    RefFace* dead_face = CAST_TO(  dead_entity,RefFace);
    if( keep_face->compare_alignment(dead_face) == CUBIT_REVERSED )
    {
      switch_sense = CUBIT_TRUE;
    }
    warn_about_refface_sense( keep_face, dead_face, switch_sense );
  }
  else if( CAST_TO( keeper_entity, RefEdge ) )
  {
    RefEdge* keep_edge = CAST_TO(keeper_entity,RefEdge);
    RefEdge* dead_edge = CAST_TO(  dead_entity,RefEdge);
    CubitSense sense;
    CubitBoolean junk;
    keep_edge->relative_sense( dead_edge, 
      GeometryQueryTool::get_geometry_factor(),
      &sense, junk, CUBIT_TRUE );
    if( sense == CUBIT_REVERSED )
      switch_sense = CUBIT_TRUE;
  }

  
    // Let any assistants know that we are about to merge
    // these two entities.
  for( int a = assistant_list_.size(); a > 0; a-- )
    assistant_list_.get_and_step()
      ->merging( keeper_entity, dead_entity, switch_sense );
  
  
    // Now find the matching pairs of GroupingEntities to merge
  GroupingEntity *keeper_GE = NULL;
  GroupingEntity *dead_GE = NULL;
  
  CubitStatus found_flag;
  keeper_GE_list.reset();
  for( int i = 0; i< keeper_GE_list.size(); i++ )
  {  
    found_flag = CUBIT_FAILURE;
    keeper_GE = keeper_GE_list.get_and_step(); 
    
      // if keeper_GE is deactivated, 
      // then skip this one and come to the next one
    if( keeper_GE->deactivated() == CUBIT_TRUE )
      continue;
    
    dead_GE_list.reset();
    for( int j = 0; j < dead_GE_list.size(); j++ )
    {
      dead_GE = dead_GE_list.get_and_step();
      
        // Before doing any compares, check to see if this
        // GroupingEntity has already been deactivated
      if( dead_GE->deactivated() == CUBIT_TRUE )
        continue;
      
        // if find any item in dead_GE_list matches the item in
        // keeper_GE_list, then merge these two items, remove
        // dead_GE from the dead_GE_list, continue to the next
        // item in keeper_GE_list.
      if( compare_GE( keeper_GE, dead_GE ) )
      {
        if( merge_GE( keeper_GE, dead_GE ) == CUBIT_FAILURE )
        {
          PRINT_ERROR("In MergeTool::merge_BTE()\n"
                      "    Cannot merge the GroupingEntities\n"
                      "  THIS IS A BUG - PLEASE REPORT IT!\n" );
          return CUBIT_FAILURE;
        }
        else
        {
            // A match was found and merged. Deactivation of
            // dead_GE is done in ModelEntity::merge_links().
          found_flag = CUBIT_SUCCESS;
          dead_GE_list.reset();
          break;
        }
      }
    }
    
    if( found_flag == CUBIT_FAILURE )
    {
      PRINT_ERROR( "In MergeTool::merge_BTE()\n"
                   "       cannot find the matching GroupingEntities.\n"
                   "  This may be due to curves smaller than the\n"
                   "  merge tolerance.  You will probably want to\n"
                   "  modify the geometry.\n");
      return CUBIT_FAILURE;
    }
  }
  
    // Merge the name(s) of dead_entity to those of keeper_entity
  keeper_entity->merge_entity_names( dead_entity );
 
  bool is_free_entity = false;
  if( dead_entity->num_parent_ref_entities() == 0 )
    is_free_entity = true;

    // Next, merge the links of these two BTEs
  SenseEntity* co_edge_ptr;
  result = CUBIT_SUCCESS;
  while( (co_edge_ptr = dead_entity->get_first_sense_entity_ptr()) )  
  {
    if (switch_sense)
      co_edge_ptr->reverse_sense();
    if (!dead_entity->remove_sense_entity(co_edge_ptr))
      result = CUBIT_FAILURE;
    if (!keeper_entity->add_sense_entity(co_edge_ptr))
      result = CUBIT_FAILURE;
  }
  
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_BTE()\n"
                 "    Could not merge the links of %s and %s.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n",
                 keeper_entity->entity_name().c_str(),
                 dead_entity->entity_name().c_str() );
    assert(0);
  }
  
    // Save IDs for geometry from dead entity.
    // Only do this for the first geometry if the dead
    // entity is already a merge of several.  Others should
    // already have IDs saved from whatever they where 
    // before they were merged.
  GeometryEntity* geom = dead_entity->get_geometry_entity_ptr();
  geom->set_saved_id( dead_entity->id() );


  // TODO -- Suggestion to make this merge code more friendly for observers:
  //         1. emit only one merge event after the merge actually happened
  //         2. delete entities after the merge event was emitted 
  
    // Destroy the old entity
   
  dead_entity->notify_all_observers( MergeEvent(dead_entity, keeper_entity) );

  dead_entity->notify_all_observers( MODEL_ENTITY_DESTRUCTED );
  if( is_free_entity ) //is free entity...top level
    dead_entity->notify_all_observers( TOP_LEVEL_ENTITY_DESTRUCTED );

  dead_entity->deactivated(CUBIT_TRUE);
  
  BridgeManager* keeper_GEs = keeper_entity->bridge_manager();
  BridgeManager* dead_GEs = dead_entity->bridge_manager();
  result = keeper_GEs->merge( dead_GEs, 
                              switch_sense ? CUBIT_REVERSED : CUBIT_FORWARD );
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_BTE()\n"
                 "    Could not merge the GeometryEntities "
                 "of %s and %s.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n",
                 keeper_entity->entity_name().c_str(),
                 dead_entity->entity_name().c_str() );
    assert(0);
  }
  
    // Debug info
  PRINT_DEBUG_19(
    "\n...Merging of %s (retained) and %s (deleted) "
    "successful.\n\n\n", keeper_entity_name.c_str(),
    dead_entity_name.c_str() );
  
    // Tell the model about the merge
    // It must be a RefFace, a RefEdge or a RefVertex
  MergeEvent merge_event(dead_entity, keeper_entity);
  merge_event.set_event_type(ENTITY_SURVIVED_MERGE);
  CubitObserver::notify_static_observers(keeper_entity, merge_event);
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : When merging RefFaces with the same sense, check
//                 CoFace senses to see if the adjacent RefVolumes overlap.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/04/04
//-------------------------------------------------------------------------
void MergeTool::warn_about_refface_sense( RefFace* face1, RefFace* face2,
                                          bool faces_reversed )
{
  DLIList<CoFace*> coface_list_1, coface_list_2;
  face1->co_faces(coface_list_1);
  face2->co_faces(coface_list_2);
  bool non_manifold_1 = false, non_manifold_2 = false;
  while (coface_list_1.size())
  {
    CoFace* coface1 = coface_list_1.pop();
    if (face1->is_nonmanifold(coface1->get_shell_ptr()))
    {
      non_manifold_1 = true;
      continue;
    }
    
    for (int i = coface_list_2.size(); i--; )
    {
      CoFace* coface2 = coface_list_2.pop();
      if (face2->is_nonmanifold(coface2->get_shell_ptr()))
      {
        non_manifold_2 = true;
        continue;
      }
      
      bool cofaces_reversed = (coface1->get_sense() == coface2->get_sense());
      if (faces_reversed != cofaces_reversed)
      {
        RefVolume* vol1 = coface1->get_shell_ptr()->get_ref_volume_ptr();
        RefVolume* vol2 = coface2->get_shell_ptr()->get_ref_volume_ptr();
        PRINT_WARNING("Merging %s (surface %d) and %s (surface %d): "
                      "%s (volume %d) and %s (volume %d) appear to be "
                      "on the same side of the surfaces.  This may "
                      "indicate bad geometry.\n",
                      face1->entity_name().c_str(), face1->id(),
                      face2->entity_name().c_str(), face2->id(),
                       vol1->entity_name().c_str(),  vol1->id(),
                       vol2->entity_name().c_str(),  vol2->id());
      }
    }
  }
      
  if (non_manifold_1)
    PRINT_WARNING("Merging non-manifold surface %d (%s).  Sheet body?\n",
      face1->id(), face1->entity_name().c_str());
  if (non_manifold_2)
    PRINT_WARNING("Merging non-manifold surface %d (%s).  Sheet body?\n",
      face2->id(), face2->entity_name().c_str());
}


//----------------------------------------------------------------
// Purpose       : Merge "this" GroupingEntity with the input
//                 GroupingEntity
//                 
// Special Notes : 
//
// Creator       : jihong Ma
//
// Creation Date : 11/26/96
//----------------------------------------------------------------
CubitStatus MergeTool::merge_GE( GroupingEntity *keeper_entity,
                                 GroupingEntity *dead_entity )
{
    // first check if merging keeper_entity with itself
  if( keeper_entity == dead_entity )
    return CUBIT_SUCCESS;
  
  CubitStatus result = CUBIT_FAILURE;
  
    // Merge is done whenever we find matching SenseEntities from
    // keeper_entity and dead_entity in compare_and_merge() with a 
    // CUBIT_TRUE input parameter.
  result = compare_and_merge( CUBIT_TRUE, keeper_entity, dead_entity );
  
    // When merging SenseEntities fails, print error and return
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_GE()\n"
                 "  Could not merge the SenseEntities\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    return CUBIT_FAILURE;
  }
  
    // After succussfully merging SenseEntities, merge the links
    // of the GeometryEntities.
  assert(dead_entity->get_children() == 0);
  dead_entity->deactivated(CUBIT_TRUE);
  BasicTopologyEntity* parent = dead_entity->get_basic_topology_entity_ptr();
  result = parent->remove_grouping_entity(dead_entity);
  
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_GE()\n"
                 "  Could not merge the GroupingEntity links.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // merge the OSME's.
  BridgeManager* keeper_OSMEs = keeper_entity->bridge_manager();
  BridgeManager* dead_OSMEs = dead_entity->bridge_manager();
  
  result = keeper_OSMEs->merge(dead_OSMEs,CUBIT_UNKNOWN);
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_GE()\n"
                 "  Could not merge the OSME's.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0);
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;      
}


//----------------------------------------------------------------
// Purpose       : Compare this GroupingEntity to the input
//                 GroupingEntity to compare for spatial equality.
//                 
// Special Notes : 
//
// Creator       : jihong Ma
//
// Creation Date : 11/26/96
//----------------------------------------------------------------
CubitBoolean MergeTool::compare_GE( GroupingEntity* keeper_entity,
                                    GroupingEntity* dead_entity )
{
    // Comparison is done by comparing the SenseEntity list of
    // keeper_entity with the SenseEntity list of dead_entity.
    // If we find that all the SenseEntities of keeper_entity
    // match with all the SenseEntities of dead_entity
    // then the GroupingEntities are spatially equal.
  
  if( compare_and_merge( CUBIT_FALSE, keeper_entity, dead_entity ) )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

//---------------------------------------------------------------
// Purpose       : spatially compare this GroupingEntity with the
//                 input GroupingEntity. But if the merge_flag is
//                 CUBIT_TRUE, merge the matching SenseEntity
//                 pairs whenever found.
//                 
// Special Notes : Merging the OSME's and merging links are done
//                 in merge_GE() function. In this function, ONLY
//                 the matching SenseEntities can be merged when
//                 asked to merge them ( i.e., merge_flag is
//                 set to CUBIT_TRUE ).
//
//                 Grouping entities compare successfully if their
//                 associated lists of SenseEntities compare
//                 successfully.  For example, 2 Loops compare if
//                 their lists of CoEdges contain matching (ones
//                 that compare) pairs.
//
//                 The first trial check that is done is to make
//                 sure that the GroupingEntities have the same
//                 number of SenseEntities.  If they do not, then
//                 we ASSUME that the GroupingEntities are not
//                 spatially equal.  This is not strictly a good
//                 assumption, but we don't have a more exact
//                 algorithm that takes care of the case where,
//                 for example, the number of RefEdges associated
//                 with two Loops are different, but the Loops
//                 themselves are spatialy equal.  Consider the
//                 case of 2 Loops representing the exact same
//                 square -- one can have 4 RefEdges and the other
//                 could have 5 (just bisect one of the previous
//                 ones...)   
//
// Creator       : jihong Ma
//
// Creation Date : 12/05/96
//---------------------------------------------------------------
CubitStatus MergeTool::compare_and_merge( CubitBoolean merge_flag,
                                          GroupingEntity* keeper_entity,
                                          GroupingEntity* dead_entity )
{
    // check if they are equal
  if( keeper_entity == dead_entity )
    return CUBIT_SUCCESS;
  
    // Make sure that the 2 GroupingEntities are of the same type
  if( keeper_entity->dag_type() !=  dead_entity->dag_type() )
  {
    PRINT_ERROR( "In MergeTool::compare_and_merge()\n"
                 "  Entities of different type are being compared.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    return CUBIT_FAILURE;
  }
  
    // Now compare, and merge if asked, the SenseEntities   
    // Get the SenseEntity lists
  DLIList<SenseEntity*> keeper_SE_list;
  DLIList<SenseEntity*> dead_SE_list;
  
  keeper_entity->get_sense_entity_list( keeper_SE_list );
  dead_entity->get_sense_entity_list( dead_SE_list );
  
    // Make sure the lists are of the same size. If not, then
    // there is a bug here as the 2 GroupingEntities should not
    // have passed the "compare" test.
  if( keeper_SE_list.size() != dead_SE_list.size() )
  {  
    if( merge_flag == CUBIT_TRUE )
    {
      PRINT_ERROR("  THIS IS A BUG - PLEASE REPORT IT!\n");
      assert(0);
    }
    return CUBIT_FAILURE;
  }
  
  SenseEntity* keeper_SE = NULL;
  SenseEntity* dead_SE = NULL;
  CubitStatus found_flag;
  
    // Now find the matching pairs of SenseEntities
  keeper_SE_list.reset();
  for( int i = 0; i < keeper_SE_list.size(); i++ )
  {
    found_flag = CUBIT_FAILURE;
    keeper_SE = keeper_SE_list.get_and_step(); 
    
      // If this SenseEntity has already been deactivated,
      // then skip to the next one
    if( keeper_SE->deactivated() == CUBIT_TRUE )
      continue;
    
    dead_SE_list.reset();
    for( int j = 0; j < dead_SE_list.size(); j++ ) 
    {
      dead_SE = dead_SE_list.get_and_step();
      
        // if dead_SE is deactivated, then skip to next one
      if( dead_SE->deactivated() == CUBIT_TRUE )
        continue;
      
        // look for an item in dead_SE_list that matches an
        // item in keeper_SE_list
      if( compare_SE( keeper_SE, dead_SE ) )
      {
          // if asked to merge, then merge the matching
          // SenseEntities remove the dead_SE from the
          // dead_SE_list.  Continue to the next item in
          // keeper_SE_list
        if( merge_flag == CUBIT_TRUE )
        {
          if( merge_SE( keeper_SE, dead_SE ) == CUBIT_FAILURE )
          {
            PRINT_ERROR( "In MergeTool::compare_and_merge()\n"
                         "   Cannot merge the matched SenseEntity\n"
                         "  THIS IS A BUG - PLEASE REPORT IT!\n" );
            return CUBIT_FAILURE;
          }
        }
        
          // if not asked to merge, just break to outer loop to
          // compare the next SenseEntity.
        found_flag = CUBIT_SUCCESS;
        dead_SE_list.reset();
        break;
      }
    }

      // If we were unable to find a match for one of the SEs,
      // return failure.
    if (found_flag == CUBIT_FAILURE)
      return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Purpose       : merge two SenseEntities
//                 
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 11/20/96
//----------------------------------------------------------------
CubitStatus MergeTool::merge_SE( SenseEntity* keeper_entity,
                                 SenseEntity* dead_entity )
{
    // first check if merging the SenseEntity with itself
  if( keeper_entity == dead_entity )
    return CUBIT_SUCCESS;
  
  CubitStatus result = CUBIT_FAILURE;
  
    // Make sure that the 2 SenseEntities are of the same type
  if( keeper_entity->dag_type() !=  dead_entity->dag_type() )
  {
    PRINT_ERROR( "In MergeTool::merge_SE()\n"
                 "   Merging SenseEntities of different types\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    return CUBIT_FAILURE;
  }
  
  BasicTopologyEntity* keeper_BTE;
  BasicTopologyEntity* dead_BTE;
  
    // Get the BTE's of the SenseEntities being merged
  keeper_BTE = keeper_entity->get_basic_topology_entity_ptr();
  dead_BTE = dead_entity->get_basic_topology_entity_ptr();

    // Compare them to make sure that these SenseEntities can
    // really be merged.
  if( keeper_BTE != dead_BTE)
  {
    PRINT_ERROR( "In MergeTool::merge_SE()\n"
                 "   Two SenseEntities have incompatible BTEs\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0);
    return CUBIT_FAILURE;
  }
  
  
  dead_entity->notify_observers( MODEL_ENTITY_DESTRUCTED );
  
    // Now that the BTE's have been successfully merged, merge the 
    // the links of the SenseEntities.
  dead_entity->deactivated(CUBIT_TRUE);
  GroupingEntity* parent = dead_entity->get_grouping_entity_ptr();
  BasicTopologyEntity* child = dead_entity->get_basic_topology_entity_ptr();
  result = (CubitStatus)(parent->remove_sense_entity(dead_entity) &&
                         child->remove_sense_entity(dead_entity));
  if( result == CUBIT_FAILURE )
  {
    PRINT_ERROR( "In MergeTool::merge_SE()\n"
                 "  Couldn't merge the links of the SenseEntities.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // Now that the BTE's have been successfully merged, merge the
    // OSME's of the SenseEntities being merged
  BridgeManager* keeper_OSMEs = keeper_entity->bridge_manager();
  BridgeManager* dead_OSMEs = dead_entity->bridge_manager();
  
  result = keeper_OSMEs->merge(dead_OSMEs,CUBIT_UNKNOWN);
  if (result == CUBIT_FAILURE)
  {
    PRINT_ERROR( "In MergeTool::merge_SE()\n"
                 "   Could not merge the OSME's.\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    assert(0) ;
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Purpose       : Compare two SenseEntities
//                 
// Special Notes : Sense entities compare successfully if their
//                 associated BasicTopologyEntities compare
//                 successfully.  For example, 2 CoEdges compare
//                 if their RefEdges compare.
//
// Creator       : jihong Ma
//
// Creation Date : 11/20/96
//----------------------------------------------------------------
CubitBoolean MergeTool::compare_SE( SenseEntity* keeper_entity,
                                    SenseEntity* dead_entity )
{
    // check if they are the same entity
  if( keeper_entity == dead_entity )
    return CUBIT_TRUE;
  
    // first check if they are the same type
  if( keeper_entity->dag_type() != dead_entity->dag_type() )
  {
    PRINT_ERROR( "In MergeTool::compare_SE()\n"
                 "   SenseEntities of different type are compared\n"
                 "  THIS IS A BUG - PLEASE REPORT IT!\n" );
    return CUBIT_TRUE;
  }
  
    // Get their BTE's and compare them
  BasicTopologyEntity *keeper_BTE = NULL;
  BasicTopologyEntity *dead_BTE = NULL;
  
  keeper_BTE = keeper_entity->get_basic_topology_entity_ptr();
  dead_BTE = dead_entity->get_basic_topology_entity_ptr();
  
    // Return the status of comparing the BTE's
  return keeper_BTE == dead_BTE;
}

void MergeTool::complete_merge()
{
     // Merge operation was completed successfully or was aborted midway.
     // Clean up the temporary compare data that was added to the RefEntities
     // being compared (during the merge) and cleanout the local (Model) lists
     // that store pointers to such entities. 

  remove_compare_data();
      
        // Now clear the lists 
  compareEntityList.clean_out() ;
  mergeSurvivorEntityList.clean_out();
  
    // Notify assistants
  for( int j = assistant_list_.size(); j--; )
    assistant_list_.get_and_step()->finish_merge();
}

void MergeTool::remove_compare_data()
{
  compareEntityList.reset();
  for(int i = compareEntityList.size(); i > 0 ; i-- )
  {
    RefEntity* ref_ent = compareEntityList.get();
    if ( CAST_TO( ref_ent, RefVertex ) ||
         CAST_TO( ref_ent, RefEdge ) ||
         CAST_TO( ref_ent, RefFace ) )
    {
        // Remove the TDCompare data attached to the entity
      PRINT_DEBUG_19( 
          "Model::notify Removing compare_TD from %s %d\n",
          ref_ent->class_name(),
          ref_ent->id());
      ref_ent->remove_compare_data();
    }
    else
    {
      PRINT_WARNING("WARNING:Something went wrong with the merging data.\n");
      compareEntityList.remove();
    }
    compareEntityList.step();
  }
  compareEntityList.clean_out();
}


void MergeTool::remove_merge_tool_assistant( MergeToolAssistant* mta_ptr )
{
  if( assistant_list_.move_to( mta_ptr ) )
    assistant_list_.remove();
}

void MergeTool::add_merge_tool_assistant( MergeToolAssistant* mta_ptr )
{
  if( !assistant_list_.is_in_list( mta_ptr ) )
    assistant_list_.append( mta_ptr );
}
MergeToolAssistant* MergeTool::find_merge_tool_assistant( const type_info& type )
{
  for( int i = assistant_list_.size(); i > 0; i-- )
  {
    if( typeid( *(assistant_list_.step_and_get()) ) == type )
      return assistant_list_.get();
  }
  return 0;
}
void MergeTool::test_r_tree(DLIList <RefFace*> &refface_list)
{
  CpuTimer timer;
  timer.cpu_secs();
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double tol = GEOMETRY_RESABS*geom_factor;
  AbstractTree <RefFace*> *a_tree = new RTree<RefFace*> (tol);
  
  DRefFaceArray refface_array( refface_list.size() );
  refface_list.reset();
  
    // Remove entities that should not be automatically merged
  int i = 0;
  int j;
  int loop_size = refface_list.size();
  for( i = 0; i < loop_size; i++ )
  {
    RefFace *curr_face = refface_list.get_and_step();
    if( curr_face->is_mergeable() )
    {
      refface_array.append( curr_face );
      a_tree->add(curr_face);
    }
  }
  double time_to_build = timer.cpu_secs();
  
    //initialize the marked flag for fast nulification...
  int array_size = refface_array.size();
  RefFace *ref_face, *ref_face1;
  CubitBox temp_box, temp_box2;
  DLIList<RefFace*> faces_in_range;
  int hit = 0;
  for ( i = 0; i < array_size; i++ )
  {
    ref_face = refface_array[i];
    temp_box = ref_face->bounding_box();
    faces_in_range.clean_out();
    a_tree->find(temp_box, faces_in_range);
    for ( j = 0; j<faces_in_range.size(); j++)
    {
      ref_face1 = faces_in_range.get_and_step();
      temp_box2 = ref_face1->bounding_box();
      if ( temp_box.overlap(tol, temp_box2) )
        hit++;
    }
  }
  PRINT_INFO( "TREE: Total Merge Reffaces time: %f secs\n"
              "\tTime to build: %f secs\n",
              timer.cpu_secs()+time_to_build, time_to_build );
}

// void MergeTool::test_r_star_tree(DLIList <RefFace*> &refface_list)
// {
//   CpuTimer timer;
//   timer.cpu_secs();
//   double geom_factor = GeometryQueryTool::get_geometry_factor();
//   double tol = GEOMETRY_RESABS*geom_factor;
//   RStarTree <RefFace*> *r_tree = new RStarTree<RefFace*> (tol);
  
//   DRefFaceArray refface_array( refface_list.size() );
//   refface_list.reset();
  
//     // Remove entities that should not be automatically merged
//   int i = 0;
//   int j;
//   int loop_size = refface_list.size();
//   for( i = 0; i < loop_size; i++ )
//   {
//     RefFace *curr_face = refface_list.get_and_step();
//     if( curr_face->is_mergeable() )
//     {
//       refface_array.append( curr_face );
//       r_tree->add(curr_face);
//     }
//   }
//   double time_to_build = timer.cpu_secs();
//     //initialize the marked flag for fast nulification...
//   int array_size = refface_array.size();
//   RefFace *ref_face, *ref_face1;
//   CubitBox temp_box, temp_box2;
//   DLIList<RefFace*> faces_in_range;
//   int hit = 0;
//   for ( i = 0; i < array_size; i++ )
//   {
//     ref_face = refface_array[i];
//     temp_box = ref_face->bounding_box();
//     faces_in_range.clean_out();
//     r_tree->find(temp_box, faces_in_range);
//     for ( j = 0; j<faces_in_range.size(); j++)
//     {
//       ref_face1 = faces_in_range.get_and_step();
//       temp_box2 = ref_face1->bounding_box();
//       if ( temp_box.overlap(tol, temp_box2) )
//         hit++;
//     }
//   }
//   PRINT_INFO( "RSTARTREE: Total Merge Reffaces time: %f secs\n"
//                  "\tTime to build: %f secs\n",
//                  timer.cpu_secs()+time_to_build, time_to_build );
// }
void MergeTool::test_no_tree(DLIList <RefFace*> &refface_list)
{
  CpuTimer timer;
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double tol = GEOMETRY_RESABS*geom_factor;
  DRefFaceArray refface_array( refface_list.size() );
  refface_list.reset();
  
    // Remove entities that should not be automatically merged
  int i = 0;
  int j;
  int loop_size = refface_list.size();
  for( i = 0; i < loop_size; i++ )
  {
    RefFace *curr_face = refface_list.get_and_step();
    if( curr_face->is_mergeable() )
      refface_array.append( curr_face );
  }
    //initialize the marked flag for fast nulification...
  int array_size = refface_array.size();
  RefFace *ref_face, *ref_face1;
  CubitBox temp_box, temp_box2;
  int hit = 0;
  for ( i = 0; i < array_size; i++ )
  {
    ref_face = refface_array[i];
    temp_box = ref_face->bounding_box();
    for ( j = i+1; j<array_size; j++)
    {
      ref_face1 = refface_array[j];
      temp_box2 = ref_face1->bounding_box();
      if ( temp_box.overlap(tol, temp_box2) )
        hit++;
    }
  }
  PRINT_INFO( "NO TREE: Merge Reffaces time: %f secs\n",
                timer.cpu_secs() );
}


//-------------------------------------------------------------------------
// Purpose       : Force-merge vertices
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
RefVertex* MergeTool::force_merge( RefVertex* vtx1, RefVertex* vtx2 )
{
  if (vtx1 == vtx2)
    return vtx1;
  
  if (vtx1->id() > vtx2->id())
    std::swap(vtx1, vtx2);
  
  DLIList<RefFace*> faces1, faces2;
  DLIList<RefEdge*> edges1, edges2;
  
  vtx1->ref_faces( faces1 );
  vtx2->ref_faces( faces2 );
  faces1.intersect( faces2 );
  if (faces1.size())
  {
    PRINT_ERROR("Vertices %d and %d are both in Surface %d.  "
                "Cannot merge vertices in the same surface.\n",
                vtx1->id(), vtx2->id(), faces1.get()->id());
    return NULL;
  }
  
  vtx1->ref_edges( edges1 );
  vtx2->ref_edges( edges2 );
  edges1.intersect( edges2 );
  if (edges1.size())
  {
    PRINT_ERROR("Vertices %d and %d are both in Curve %d.  "
                "Cannot merge vertices in the same curve.\n",
                vtx1->id(), vtx2->id(), edges1.get()->id());
    return NULL;
  }
  
  int k;
  for( k = assistant_list_.size(); k > 0; k-- )
    if( ! assistant_list_.get_and_step()->can_merge( vtx1, vtx2 ) )
      break;
  
  if (k)
    return NULL;
  
  if (!merge_BTE( vtx1, vtx2 ))
    return NULL;
  
  if (destroyDeadGeometry)
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  
  return vtx1;
}


//-------------------------------------------------------------------------
// Purpose       : Force-merge curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
RefEdge* MergeTool::force_merge( RefEdge* edge1, RefEdge* edge2 )
{
  if (edge1 == edge2)
    return edge1;
  
  if (edge1->id() > edge2->id())
    std::swap(edge1, edge2);
  
  DLIList<RefFace*> faces1, faces2;
  
  edge1->ref_faces( faces1 );
  edge2->ref_faces( faces2 );
  faces1.intersect( faces2 );
  if (faces1.size())
  {
    PRINT_ERROR("Curves %d and %d are both in Surface %d.  "
                "Cannot merge curves in the same surface.\n",
                edge1->id(), edge2->id(), faces1.get()->id());
    return NULL;
  }
  
  bool closed1 = edge1->start_vertex() == edge1->end_vertex();
  bool closed2 = edge2->start_vertex() == edge2->end_vertex();
  if (closed1 != closed2)
  {
    PRINT_ERROR("Curves %d and %d do not have the same number of vertices.\n",
      edge1->id(), edge2->id());
    return NULL;
  }
  
  int k;
  for( k = assistant_list_.size(); k > 0; k-- )
    if( ! assistant_list_.get_and_step()->can_merge( edge1, edge2 ) )
      break;
  
  if (k)
    return NULL;
    
  
  CubitSense sense;
  CubitBoolean equal;
  if (!edge1->relative_sense( edge2, GEOMETRY_RESABS, &sense, equal, true))
    return NULL;

  if (closed1)
  {
    if (edge1->start_vertex() != edge2->start_vertex())
      edge1->start_vertex()->notify( edge2->start_vertex(), COMPARISON_FOUND );
  }
  else
  {
    RefVertex *start2 = edge2->start_vertex();
    RefVertex *end2   = edge2->end_vertex();
    if (sense == CUBIT_REVERSED)
      std::swap(start2, end2);
      
    if (edge1->start_vertex() == end2 || edge1->end_vertex() == start2)
    {
      PRINT_ERROR("Error merging curves:  invalid relative sense calculation.\n");
      return NULL;
    }
    
    if (edge1->start_vertex() != start2)
      edge1->start_vertex()->notify( start2, COMPARISON_FOUND );
    if (edge1->end_vertex() != end2)
      edge1->end_vertex()->notify( end2, COMPARISON_FOUND );
  }  
  
  if (!merge_BTE( edge1, edge2 ))
  {
    PRINT_ERROR("Merge failed.\n");
    remove_compare_data();
    return NULL;
  }
  
  remove_compare_data();
  if (destroyDeadGeometry)
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  return edge1;
}

//-------------------------------------------------------------------------
// Purpose       : Force-merge RefFaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
RefFace* MergeTool::force_merge( RefFace* face1, RefFace* face2 )
{
  if (face1 == face2)
    return face1;
  
  int i, j;
  if (face1->id() > face2->id())
    std::swap(face1, face2);
  
  DLIList<RefVolume*> vols1, vols2;
  
  face1->ref_volumes( vols1 );
  face2->ref_volumes( vols2 );
  vols1.intersect( vols2 );
  if (vols1.size())
  {
    PRINT_ERROR("Surfaces %d and %d are both in Volume %d.  "
                "Cannot merge surfaces in the same volume.\n",
                face1->id(), face2->id(), vols1.get()->id());
    return NULL;
  }
  
  DLIList<RefVertex*> vertices1, vertices2;
  DLIList<Loop*> loops1, loops2, tmp_loops;
  
  face1->loops( loops1 );
  face2->loops( loops2 );
  face1->ref_vertices( vertices1 );
  face2->ref_vertices( vertices2 );
  
  if (loops1.size() != loops2.size() || vertices1.size() != vertices2.size())
  {
    PRINT_ERROR("Surfaces %d and %d do not have equivalent topology.\n",
      face1->id(), face2->id() );
    return 0;
  }
  
    // Need relative sense so we know in which order to compare
    // coedges of a loop and the expected sense of the coedges.
  CubitSense rel_sense = face1->compare_alignment( face2 );
  if (rel_sense == CUBIT_UNKNOWN)
  {
    PRINT_ERROR("Cound not calculate relative sense of surfaces.\n");
    return 0;
  }
  
    // Remove from consideration any vertices that are already merged.
  DLIList<RefVertex*> common( vertices1 );
  common.intersect( vertices2 );
  vertices1 -= common;
  vertices2 -= common;
  
    // Match vertex pairs
  vertices1.reset();
  for (i = vertices1.size(); i--; )
  {
    RefVertex* vtx1 = vertices1.get_and_step();
    
    RefVertex* closest = 0;
    double shortest = CUBIT_DBL_MAX;
    vertices2.reset();
    for (j = vertices2.size(); j--; )
    {
      RefVertex* vtx2 = vertices2.get_and_step();
      double d = (vtx1->coordinates() - vtx2->coordinates()).length_squared();
      if (d < shortest)
      {
        shortest = d;
        closest = vtx2;
      }
    }
    
    vertices2.move_to( closest );
    vertices2.extract();
    if (vtx1 != closest)
      vtx1->notify( closest, COMPARISON_FOUND );
  }
  
    // Compare all loops, find RefEdge merge pairs
  DLIList<CoEdge*> coedges1, coedges2;
  loops1.reset();
  for (i = loops1.size(); i--; )
  {
    Loop* loop1 = loops1.get_and_step();
    
      // Given a vertex in a loop in the first surface,
      // find the loop connected to that vertex and in
      // the second surface.
    vertices1.clean_out();
    loop1->ref_vertices( vertices1 );
    RefVertex* vtx = vertices1.get();
    if (vtx->get_compare_partner())
      vtx = dynamic_cast<RefVertex*>(vtx->get_compare_partner());
    tmp_loops.clean_out();
    vtx->loops( tmp_loops );
    Loop* loop2 = 0;
    for (j = tmp_loops.size(); j--; )
    {
      if (tmp_loops.step_and_get()->get_ref_face_ptr() == face2)
      {
        loop2 = tmp_loops.get();
        break;
      }
    }
    
    if (!loop2 || !loops2.move_to(loop2))
    {
      PRINT_ERROR("RefFace topology does not match.  Cannot merge.\n");
      remove_compare_data();
      return 0;
    }
    loops2.extract();
    
      // Compare loop coedges
    coedges1.clean_out();
    coedges2.clean_out();
    loop1->co_edges( coedges1 );
    loop2->co_edges( coedges2 );
    if (coedges1.size() != coedges2.size())
    {
      PRINT_ERROR("RefFace topology does not match.  Cannot merge.\n");
      remove_compare_data();
      return 0;
    }
    if (CUBIT_REVERSED == rel_sense)
      coedges2.reverse();
    
      // Given a coedge in the first loop, find the 
      // matching coedge in the second loop
    coedges1.reset();
    coedges2.reset();
    RefVertex* start1 = coedges1.get()->start_vertex();
    RefVertex* end1 = coedges1.get()->end_vertex();
    if (CUBIT_REVERSED == rel_sense)
      std::swap(start1, end1);
    if (start1->get_compare_partner())
      start1 = dynamic_cast<RefVertex*>(start1->get_compare_partner());
    if (end1->get_compare_partner())
      end1 = dynamic_cast<RefVertex*>(end1->get_compare_partner());  
      
    for (j = coedges2.size(); j > 0; j-- )
    {
      if (coedges2.get()->start_vertex() == start1 &&
          coedges2.get()->end_vertex() == end1)
        break;
      
      coedges2.step();
    }
    
      // No matching coedge
    if (!j)
    {
      PRINT_ERROR("RefFace topology does not match.  Cannot merge.\n");
      remove_compare_data();
      return 0;
    }
    
      // Check that remaining coedges match, and mark RefEdges accordingly
    for (j = coedges1.size(); j--; )
    {
      CoEdge* coedge1 = coedges1.get_and_step();
      CoEdge* coedge2 = coedges2.get_and_step();
      if (coedge1->get_ref_edge_ptr() == coedge2->get_ref_edge_ptr())
        continue;
      
      RefVertex* start1 = coedge1->start_vertex();
      RefVertex* end1 = coedge1->end_vertex();
      if (CUBIT_REVERSED == rel_sense)
        std::swap(start1, end1);
      if (start1->get_compare_partner())
        start1 = dynamic_cast<RefVertex*>(start1->get_compare_partner());
      if (end1->get_compare_partner())
        end1 = dynamic_cast<RefVertex*>(end1->get_compare_partner());  
      
      if (coedge2->start_vertex() != start1 ||
          coedge2->end_vertex() != end1)
      {
        PRINT_ERROR("RefFace topology does not match.  Merge aborted.\n");
        remove_compare_data();
        return 0;
      }
      
      coedge1->get_ref_edge_ptr()->notify( coedge2->get_ref_edge_ptr(), 
                                           COMPARISON_FOUND );
    }
  } // for(loops1)
  
    // check if mesh can be merged, etc.
  for( i = assistant_list_.size(); i > 0; i-- )
    if( ! assistant_list_.get_and_step()->can_merge( face1, face2 ) )
      { remove_compare_data();  return 0; }
  
    // merge
  if (!merge_BTE( face1, face2 ))
  {
    PRINT_ERROR("Merge Failed.\n");
    remove_compare_data();
    return 0;
  }
  
  remove_compare_data();
  if (destroyDeadGeometry)
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  return face1;
}

//-------------------------------------------------------------------------
// Purpose       : Force-merge RefEntities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
RefEntity* MergeTool::force_merge( RefEntity* ent1, RefEntity* ent2 )
{
  if (RefFace* face1 = dynamic_cast<RefFace*>(ent1))
  {
    if (RefFace* face2 = dynamic_cast<RefFace*>(ent2))
      return force_merge( face1, face2 );
  }
  else if (RefEdge* edge1 = dynamic_cast<RefEdge*>(ent1))
  {
    if (RefEdge* edge2 = dynamic_cast<RefEdge*>(ent2))
      return force_merge( edge1, edge2 );
  }
  else if (RefVertex* vtx1 = dynamic_cast<RefVertex*>(ent1))
  {
    if (RefVertex* vtx2 = dynamic_cast<RefVertex*>(ent2))
      return force_merge( vtx1, vtx2 );
  }
  
  PRINT_ERROR("Invalid entities passed to MergeTool::force_merge\n");
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Force-merge RefEntities
//
// Special Notes : Provied for use by CAMergeAttribute only
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
RefEntity* MergeTool::force_merge( const DLIList<RefEntity*>& list )
{
  RefEntity* result = list.get();
  for (int i = 1; result && i < list.size(); i++ )
    result = force_merge( list.next(i), result );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Un-merge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::separate_bodies( DLIList<Body*>& separate_list,
                                        DLIList<Body*>* from_list )
{
  CubitBoolean top = start_unmerge();
  DLIList<RefVolume*> volume_list, tmp_list;
  separate_list.reset();
  for (int i = separate_list.size(); i--; )
  {
    Body* body = separate_list.get_and_step();
    tmp_list.clean_out();
    body->ref_volumes( tmp_list );
    volume_list += tmp_list;
  }
  
  CubitStatus result;
  if (from_list == NULL)
  {
    result = separate_volumes( volume_list, NULL );
  }
  else
  {
    DLIList<RefVolume*> from_vols;
    from_list->reset();
    for (int i = from_list->size(); i--; )
    {
      Body* body = from_list->get_and_step();
      tmp_list.clean_out();
      body->ref_volumes( tmp_list );
      from_vols += tmp_list;
    }
   
    result = separate_volumes( volume_list, &from_vols );
  }
  
  end_unmerge(top);
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::separate_volumes( DLIList<RefVolume*>& split_list,
                                         DLIList<RefVolume*>* from_list )
{
  DLIList<TopologyEntity*> entity_list( split_list.size() );
  CAST_LIST_TO_PARENT( split_list, entity_list );
  if (from_list)
  {
    DLIList<TopologyEntity*> from_list2( from_list->size() );
    CAST_LIST_TO_PARENT( *from_list, from_list2 );
    return separate_entities( entity_list, &from_list2 );
  }
  else
  { 
    return separate_entities( entity_list );
  }
}
CubitStatus MergeTool::separate_faces( DLIList<RefFace*>& split_list,
                                       DLIList<RefFace*>* from_list )
{
  DLIList<TopologyEntity*> entity_list( split_list.size() );
  CAST_LIST_TO_PARENT( split_list, entity_list );
  if (from_list)
  {
    DLIList<TopologyEntity*> from_list2( from_list->size() );
    CAST_LIST_TO_PARENT( *from_list, from_list2 );
    return separate_entities( entity_list, &from_list2 );
  }
  else
  { 
    return separate_entities( entity_list );
  }
}
CubitStatus MergeTool::separate_edges( DLIList<RefEdge*>& split_list,
                                       DLIList<RefEdge*>* from_list )
{
  DLIList<TopologyEntity*> entity_list( split_list.size() );
  CAST_LIST_TO_PARENT( split_list, entity_list );
  if (from_list)
  {
    DLIList<TopologyEntity*> from_list2( from_list->size() );
    CAST_LIST_TO_PARENT( *from_list, from_list2 );
    return separate_entities( entity_list, &from_list2 );
  }
  else
  { 
    return separate_entities( entity_list );
  }
}


//-------------------------------------------------------------------------
// Purpose       : Unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/02/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::separate_entities( DLIList<TopologyEntity*>& volume_list,
                                          DLIList<TopologyEntity*>* from_list )
{
  int i, j, k, l;
  CubitBoolean top = start_unmerge();
  CubitStatus result = CUBIT_SUCCESS;
  
  if (!volume_list.size())
    return CUBIT_FAILURE;
  
  DagType parent_type = volume_list.get()->dag_type();
  for (i = volume_list.size(); i--;)
    if (volume_list.get_and_step()->dag_type() != parent_type)
      return CUBIT_FAILURE;
  if (from_list)
    for (i = from_list->size(); i--; )
      if (from_list->get_and_step()->dag_type() != parent_type)
        return CUBIT_FAILURE;
    
  
  DLIList<ModelEntity*> query_input( volume_list.size() ), query_output;
  CAST_LIST_TO_PARENT( volume_list, query_input );
  DLIList<TopologyBridge*> bridge_list, split_list, parent_bridges;
  DLIList<TopologyEntity*> entity_vol_list;
  DLIList<BasicTopologyEntity*> entity_list;
  DLIList<Lump*> lump_list;
  DLIList<Surface*> surface_list;
  DLIList<Curve*> curve_list;
  DLIList<Point*> point_list;
  
  
    // Loop once for each type, 
  DagType types[] = { DagType::ref_face_type(), 
                      DagType::ref_edge_type(),
                      DagType::ref_vertex_type() };
  for (j = 0; j < 3; j++)
  {
    if (types[j] >= parent_type)
      continue;
      
        // Get RefFaces/RefEdges/RefVertices from RefVolumes
    query_output.clean_out();
    entity_list.clean_out();
    ModelQueryEngine::instance()->query_model( query_input, types[j], query_output );
    CAST_LIST( query_output, entity_list, BasicTopologyEntity );

    entity_list.reset();
    for (i = entity_list.size(); i--; )
    {
      BasicTopologyEntity* entity = entity_list.get_and_step();
      split_list.clean_out();
    
        // Get parent volumes to unmerge from 
      entity_vol_list.clean_out();
      query_output.clean_out();
      ModelQueryEngine::instance()->query_model( *entity, parent_type, query_output );
      CAST_LIST( query_output, entity_vol_list, TopologyEntity );
      entity_vol_list -= volume_list;
      if (from_list)
        entity_vol_list.intersect( *from_list );
      if (!entity_vol_list.size())
        continue;
    
        // For each merged bridge in the entity
      bridge_list.clean_out();
      entity->bridge_manager()->get_bridge_list( bridge_list );
      bridge_list.reset();
      for (k = bridge_list.size(); k--; )
      { 
        TopologyBridge* bridge = bridge_list.get_and_step();
      
        parent_bridges.clean_out();
        if (parent_type == DagType::ref_volume_type())
        {
          lump_list.clean_out();
          bridge->lumps( lump_list );
          CAST_LIST_TO_PARENT( lump_list, parent_bridges );
        }
        else if(parent_type == DagType::ref_face_type())
        {
          surface_list.clean_out();
          bridge->surfaces( surface_list );
          CAST_LIST_TO_PARENT( surface_list, parent_bridges );
        }
        else if(parent_type == DagType::ref_edge_type())
        {
          curve_list.clean_out();
          bridge->curves( curve_list );
          CAST_LIST_TO_PARENT( curve_list, parent_bridges );
        }
        else
        {
          assert(0);
          return CUBIT_FAILURE;
        }
        
          // Check if bridge is in RefVolumes to unmerge from
          
        for (l = parent_bridges.size(); l--; )
        {
          TopologyBridge* lump = parent_bridges.get_and_step();
          if (entity_vol_list.is_in_list(lump->topology_entity()))
            split_list.append( bridge );
        }
      }
    
        // If there are bridges to unmerge...
      if (split_list.size() != 0 && split_list.size() != bridge_list.size())
      {
        if (types[j] == DagType::ref_face_type())
        {
          surface_list.clean_out();
          CAST_LIST( split_list, surface_list, Surface );
          RefFace* face = separate_face( surface_list, false );
          if (0 == face)
            result = CUBIT_FAILURE;
        }
        else if(types[j] == DagType::ref_edge_type())
        {
          curve_list.clean_out();
          CAST_LIST( split_list, curve_list, Curve );
          RefEdge* edge = separate_edge( curve_list, false );
          if (0 == edge)
            result = CUBIT_FAILURE;
        }
        else 
        {
          point_list.clean_out();
          CAST_LIST( split_list, point_list, Point );
          assert( split_list.size() == point_list.size() );
          RefVertex* vtx = separate_vertex( point_list );
          if (0 == vtx)
            result = CUBIT_FAILURE;
        }
      }
    } // for(i in entity_list)
  } // for(j in type)

  end_unmerge(top);
  return result;
}



//-------------------------------------------------------------------------
// Purpose       : Common code for separate functions.
//                 Check if bridges have same owner (are merged together),
//                 and if check_parents == true, check if unmerging the 
//                 bridges will invalidate the parent topology.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
BasicTopologyEntity*  MergeTool::can_separate( DLIList<TopologyBridge*>& bridges,
                                               bool check_parents )
{
  int i;
  BasicTopologyEntity* bte = 0;
  
  TBOwner* owner = bridges.get_and_step()->owner();
  for (i = bridges.size(); i > 1; i--)
    if (bridges.get_and_step()->owner() != owner)
      return 0;
  
  BridgeManager* bmanager = dynamic_cast<BridgeManager*>(owner);
  if (!bmanager)
    { assert(0); return 0; }
  
  bte = dynamic_cast<BasicTopologyEntity*>(bmanager->topology_entity());
  if (!bte)
    { assert(0); return 0; }
  
  if (bmanager->number_of_bridges() == bridges.size())
    return 0;
  
  if (!check_parents) 
    return bte;
  
  DLIList<TopologyBridge*> parent_bridges, bte_bridges;
  DLIList<TBOwner*> parent_owners;
  bte->bridge_manager()->get_bridge_list( bte_bridges );
  bte_bridges -= bridges;
  
  for (i = bte_bridges.size(); i--;)
  {
    bte_bridges.step_and_get()->get_parents( parent_bridges );
    while (parent_bridges.size())
      parent_owners.append( parent_bridges.pop()->owner() );
  }
  for (i = bridges.size(); i--; )
  {
    bridges.step_and_get()->get_parents( parent_bridges );
    while (parent_bridges.size())
      if (parent_owners.is_in_list( parent_bridges.pop()->owner() ))
        return 0;
  }
  
  return bte;
}


//-------------------------------------------------------------------------
// Purpose       : Split a RefFace
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
RefFace* MergeTool::separate_face( DLIList<Surface*>& surfaces,
                                   bool unmerge_curves )
{
  int i, j, k;
  
  DLIList<TopologyBridge*> bridge_list( surfaces.size() );
  CAST_LIST_TO_PARENT( surfaces, bridge_list );
  RefFace* old_face = dynamic_cast<RefFace*>( can_separate( bridge_list, false ) );
  if (!old_face)
    return 0;
  
  CubitBoolean top = start_unmerge();
  
    // Split the RefFace
    
    // Remove surfaces from old entity
  for (i = surfaces.size(); i--; )
    old_face->bridge_manager()->remove_bridge( surfaces.get_and_step() );
  check_saved_id( old_face );
  surfaces.reset();
    // Use surface with smallest saved ID as primary
  for (j = 0, i = 1; i < surfaces.size(); i++)
    if (surfaces.next(i)->get_saved_id() < surfaces.next(j)->get_saved_id())
      j = i;
  surfaces.step( j );
    // Create new face
  RefFace* new_face = RefEntityFactory::instance()->construct_RefFace( surfaces.get() );
  for (i = surfaces.size(); i > 1; i-- )
    new_face->bridge_manager()->add_bridge( surfaces.step_and_get() );
  
    // Move CoFaces
  bridge_list.clean_out();
  DLIList<TopologyBridge*> bridge_parents;
  DLIList<TopologyEntity*> shells;
  DLIList<CoFace*> co_faces;
  old_face->co_faces( co_faces );
  for (i = surfaces.size(); i > 0; i--)
  {
    surfaces.get_and_step()->get_parents( bridge_list );
    while (bridge_list.size())
      shells.append( bridge_list.pop()->topology_entity() );
  }
  co_faces.reset();
  for (i = co_faces.size(); i--; )
  {
    CoFace* co_face = co_faces.get_and_step();
    if (shells.is_in_list( co_face->get_grouping_entity_ptr() ))
    {
      co_face->switch_basic_topology_entity( new_face );
    }
  }
  
    // Split loops and coedges
  DLIList<Loop*> loops;
  DLIList<CoEdge*> coedges;
  DLIList<SenseEntity*> new_coedges, junk_list;
  DLIList<TopologyBridge*> loop_bridges;
  old_face->loops( loops );
  loops.reset();
  for (i = loops.size(); i--; )
  {
    Loop* old_loop = loops.get_and_step();
    Loop* new_loop = new Loop;
    new_face->add_grouping_entity( new_loop );
    
    loop_bridges.clean_out();
    bridge_list.clean_out();
    old_loop->bridge_manager()->get_bridge_list( bridge_list );
    bridge_list.reset();
    for (j = bridge_list.size(); j--; )
    {
      TopologyBridge* loopsm = bridge_list.get_and_step();
      loopsm->get_parents( bridge_parents );
      Surface* loop_surf = dynamic_cast<Surface*>(bridge_parents.pop());
      assert (loop_surf && !bridge_parents.size());
      if (surfaces.is_in_list(loop_surf))
      {
        old_loop->bridge_manager()->remove_bridge( loopsm );
        new_loop->bridge_manager()->add_bridge( loopsm );
        loop_bridges.append( loopsm );
      }
    }
    
    new_coedges.clean_out();
    coedges.clean_out();
    old_loop->co_edges( coedges );
    coedges.reset();
    for (j = coedges.size(); j--; )
    {
      CoEdge* old_coedge = coedges.get_and_step();
      CoEdge* new_coedge = new CoEdge( old_coedge->get_ref_edge_ptr(),
                                       old_coedge->get_sense() );
      new_coedges.append( new_coedge );
      
      bridge_list.clean_out();
      old_coedge->bridge_manager()->get_bridge_list( bridge_list );
      bridge_list.reset();
      for (k = bridge_list.size(); k--; )
      {
        TopologyBridge* coedgesm = bridge_list.get_and_step();
        coedgesm->get_parents( bridge_parents );
        LoopSM* coedgesm_loop = dynamic_cast<LoopSM*>(bridge_parents.pop());
        assert (coedgesm_loop && !bridge_parents.size());
        if (loop_bridges.is_in_list(coedgesm_loop))
        {
          old_coedge->bridge_manager()->remove_bridge( coedgesm );
          new_coedge->bridge_manager()->add_bridge( coedgesm );
        }
      }
    }
    
    new_loop->set_sense_entity_list( new_coedges, junk_list );
  }
  
    // Check if should reverse new face
  bool reverse_new = true;
  for (i = surfaces.size(); i--; )
    if (surfaces.get_and_step()->bridge_sense() != CUBIT_REVERSED)
      reverse_new = false;
  if (reverse_new)
  {
    new_face->bridge_manager()->reverse_bridge_senses();
    new_face->reverse_topology();
  }
  
    // Check if should reverse old face
  bool reverse_old = true;
  bridge_list.clean_out();
  old_face->bridge_manager()->get_bridge_list( bridge_list );
  for (i = bridge_list.size(); i--; )
    if (bridge_list.get_and_step()->bridge_sense() != CUBIT_REVERSED)
      reverse_old = false;
  if (reverse_old)
  {
    old_face->bridge_manager()->reverse_bridge_senses();
    old_face->reverse_topology();
  }
  
    // Misc stuff for updating other code for changed topology
  bool reversed = reverse_old != reverse_new;
  new_unmerged.append( new_face );
  old_unmerged.append( old_face );
  for (i = assistant_list_.size(); i--; )
    assistant_list_.get_and_step()->unmerged( old_face, new_face, reversed );
  
  if (unmerge_curves)
  {  
    DLIList<RefEdge*> edge_list;
    DLIList<Curve*> curve_list;
    DLIList<TopologyBridge*> bridge_children;
    new_face->ref_edges( edge_list );
    edge_list.reset();
    
    for (i = edge_list.size(); i--;)
    {
      curve_list.clean_out();
      
      RefEdge* edge = edge_list.get_and_step();
      bridge_list.clean_out();
      edge->bridge_manager()->get_bridge_list( bridge_list );
      if (bridge_list.size() < 2)
        continue;
      
        // Find the curve(s) in the just-unmerged refface
      bridge_list.reset();
      for (j = bridge_list.size(); j--; )
      {
        TopologyBridge* bridge = bridge_list.get_and_step();
        bridge_parents.clean_out();
        bridge->get_parents( bridge_parents );
        bridge_parents.reset();
        bool in_old = false, in_new = false;
        while (bridge_parents.size())
        {
          CoEdge* coedge = dynamic_cast<CoEdge*>(bridge_parents.pop()->topology_entity());
          BasicTopologyEntity* bte = coedge->get_parent_basic_topology_entity_ptr();
          if (bte == new_face)
            in_new = true;
          else if(bte == old_face)
            in_old = true;
        }
        
        if (in_old && in_new)
        {
          curve_list.clean_out();
          break;
        }
        else if(in_new)
        {
          curve_list.append( dynamic_cast<Curve*>(bridge) );
        }
      } // for( j in bridge_list )
      
        // Find curve(s) that must remain merged.
      curve_list.reset();
      for (j = 0; j < curve_list.size(); j++)
      {
        bridge_parents.clean_out();
        curve_list.get_and_step()->get_parents( bridge_parents );
        bridge_parents.reset();
        for (k = 0; k< bridge_parents.size(); k++)
        {
          bridge_list.clean_out();
          BridgeManager* bm = bridge_parents.get_and_step()->bridge_manager();
          bm->get_bridge_list( bridge_list );
          
          bridge_list.reset();
          for (int l = bridge_list.size(); l--; )
          {
            bridge_children.clean_out();
            bridge_list.get_and_step()->get_children( bridge_children );
            assert (bridge_children.size() == 1); // bridges are coedges, must have 1 curve
            Curve* curve = dynamic_cast<Curve*>(bridge_children.get());
            assert (curve->owner() == edge->bridge_manager());
            curve_list.append_unique( curve );
          }
        }
      } // for( j in curve_list )
      
      if (curve_list.size() != 0 && 
          curve_list.size() != edge->bridge_manager()->number_of_bridges())
      {
        separate_edge( curve_list, true );
      }
    } // for (i in edge_list)
  } // if (unmerge_curves)
  
  end_unmerge(top);
  return new_face;
}


//-------------------------------------------------------------------------
// Purpose       : Split a RefEdge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
RefEdge* MergeTool::separate_edge( DLIList<Curve*>& curves,
                                   bool unmerge_vertices )
{
  int i, j, k;
  
  DLIList<TopologyBridge*> bridge_list( curves.size() );
  CAST_LIST_TO_PARENT( curves, bridge_list );
  RefEdge* old_edge = dynamic_cast<RefEdge*>( can_separate( bridge_list, true ) );
  if (!old_edge)
    return 0;
    
  CubitBoolean top = start_unmerge();
     
    // Split the RefEdge
    
    // Remove curves from old edge
  for (i = curves.size(); i--; )
    old_edge->bridge_manager()->remove_bridge( curves.get_and_step() );
  check_saved_id( old_edge );
  curves.reset();
    // Use curve with smallest saved ID as primary
  for (j = 0, i = 1; i < curves.size(); i++)
    if (curves.next(i)->get_saved_id() < curves.next(j)->get_saved_id())
      j = i;
  curves.step( j );
    // Create new edge
  RefEdge* new_edge = RefEntityFactory::instance()->construct_RefEdge( curves.get() );
  for (i = curves.size(); i > 1; i-- )
    new_edge->bridge_manager()->add_bridge( curves.step_and_get() );
  
    // Move CoEdges
  bridge_list.clean_out();
  DLIList<TopologyEntity*> new_coedges;
  DLIList<CoEdge*> old_coedges;
  old_edge->co_edges( old_coedges );
  for (i = curves.size(); i > 0; i--)
  {
    curves.get_and_step()->get_parents( bridge_list );
    while (bridge_list.size())
      new_coedges.append( bridge_list.pop()->topology_entity() );
  }
  old_coedges.reset();
  for (i = old_coedges.size(); i--; )
  {
    CoEdge* co_edge = old_coedges.get_and_step();
    if (new_coedges.is_in_list( co_edge ))
    {
      co_edge->switch_basic_topology_entity( new_edge );
    }
  }
  
    // Construct chain and co-vertices
  Chain* new_chain = new Chain;
  new_edge->add_grouping_entity( new_chain );
  CoVertex* start = new CoVertex( old_edge->start_vertex() );
  new_chain->add_sense_entity( start );
  CoVertex* end   = new CoVertex( old_edge->end_vertex() );
  new_chain->add_sense_entity( end, start );
  
    // Check if should reverse new edge
  bool reverse_new = true;
  for (i = curves.size(); i--; )
    if (curves.get_and_step()->bridge_sense() != CUBIT_REVERSED)
      reverse_new = false;
  if (reverse_new)
  {
    new_edge->bridge_manager()->reverse_bridge_senses();
    new_edge->reverse_topology();
  }
  
    // Check if should reverse old edge
  bool reverse_old = true;
  bridge_list.clean_out();
  old_edge->bridge_manager()->get_bridge_list( bridge_list );
  for (i = bridge_list.size(); i--; )
    if (bridge_list.get_and_step()->bridge_sense() != CUBIT_REVERSED)
      reverse_old = false;
  if (reverse_old)
  {
    old_edge->bridge_manager()->reverse_bridge_senses();
    old_edge->reverse_topology();
  }
  
    // Misc stuff for updating other code for changed topology
  bool reversed = reverse_old != reverse_new;
  new_unmerged.append( new_edge );
  old_unmerged.append( old_edge );
  for (i = assistant_list_.size(); i--; )
    assistant_list_.get_and_step()->unmerged( old_edge, new_edge, reversed );
  
  if (unmerge_vertices)
  {
      // Split vertices
    DLIList<Point*> point_list;
    DLIList<TopologyBridge*> parent_curves, curve_points;

    int n_vert = old_edge->num_ref_vertices();
    assert(n_vert < 3);

    for (i = 0; i < n_vert; i++)
    {
      point_list.clean_out();
      RefVertex* vtx = i ? old_edge->end_vertex() : old_edge->start_vertex();
      CubitSense sense = i ? CUBIT_REVERSED : CUBIT_FORWARD;

      curves.reset();
      for (j = curves.size(); j--; )
      {
        Curve* curve = curves.get_and_step();
        bridge_list.clean_out();
        curve->get_children( bridge_list );
        bridge_list.reset();
        if (reversed != (curve->bridge_sense() != sense))
          bridge_list.step();
        Point* point = dynamic_cast<Point*>(bridge_list.get());
        assert (point->owner() == vtx->bridge_manager());
        point_list.append( point );
      }
      
      point_list.reset();
      for (j = 0; j < point_list.size(); j++)
      {
        Point* point = point_list.get_and_step();
        parent_curves.clean_out();
        point->get_parents( parent_curves );
        
        parent_curves.reset();
        for (k = parent_curves.size(); k--; )
        {
          TopologyBridge* curve = parent_curves.get_and_step();
          BridgeManager* bm = curve->bridge_manager();
          bridge_list.clean_out();
          bm->get_bridge_list( bridge_list );
          
          bridge_list.reset();
          for (int l = bridge_list.size(); l--; )
          {
            curve_points.clean_out();
            bridge_list.get_and_step()->get_children( curve_points );
            assert(curve_points.size() < 3);
            if (curve_points.get()->owner() == point->owner())
              point_list.append_unique( dynamic_cast<Point*>(curve_points.get()) );
            else if (curve_points.next()->owner() == point->owner())
              point_list.append_unique( dynamic_cast<Point*>(curve_points.next()) );
          }
        }
      }
 
      separate_vertex( point_list );
    }
  }
  
  end_unmerge(top);  
  return new_edge;
}
  


//-------------------------------------------------------------------------
// Purpose       : Split a RefVertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/01/04
//-------------------------------------------------------------------------
RefVertex* MergeTool::separate_vertex( DLIList<Point*>& points )
{
  int i, j;
  
  DLIList<TopologyBridge*> bridge_list( points.size() );
  CAST_LIST_TO_PARENT( points, bridge_list );
  RefVertex* old_vtx = dynamic_cast<RefVertex*>( can_separate( bridge_list, true ) );
  if (!old_vtx)
    return 0;
  
  CubitBoolean top = start_unmerge();
  
  
    // Split the RefVertex
  
    // Remove points from old vertex
  for (i = points.size(); i--; )
    old_vtx->bridge_manager()->remove_bridge( points.get_and_step() );
  check_saved_id( old_vtx );
  points.reset();
    // Use point with smallest saved ID as primary
  for (j = 0, i = 1; i < points.size(); i++)
    if (points.next(i)->get_saved_id() < points.next(j)->get_saved_id())
      j = i;
  points.step( j );
    // Create new vertex
  RefVertex* new_vtx = RefEntityFactory::instance()->construct_RefVertex( points.get() );
  for (i = points.size(); i > 1; i-- )
    new_vtx->bridge_manager()->add_bridge( points.step_and_get() );
  
    // Move CoVertices
  bridge_list.clean_out();
  DLIList<TopologyEntity*> edges;
  DLIList<CoVertex*> co_vertices;
  old_vtx->co_vertices( co_vertices );
  for (i = points.size(); i > 0; i--)
  {
    points.get_and_step()->get_parents( bridge_list );
    while (bridge_list.size())
      edges.append( bridge_list.pop()->topology_entity() );
  }
  co_vertices.reset();
  for (i = co_vertices.size(); i--; )
  {
    CoVertex* covtx = co_vertices.get_and_step();
    if (edges.is_in_list( covtx->get_parent_basic_topology_entity_ptr() ))
    {
      covtx->switch_basic_topology_entity( new_vtx );
    }
  }
  
    // Misc stuff for updating other code for changed topology
  new_unmerged.append( new_vtx );
  old_unmerged.append( old_vtx );
  for (i = assistant_list_.size(); i--; )
    assistant_list_.get_and_step()->unmerged( old_vtx, new_vtx, false );

  end_unmerge(top);
  return new_vtx;
}

    
//-------------------------------------------------------------------------
// Purpose       : Set entity ID to smallest saved bridge ID
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/03/04
//-------------------------------------------------------------------------
CubitStatus MergeTool::check_saved_id( BasicTopologyEntity* bte )
{
  DLIList<TopologyBridge*> bridge_list;
  bte->bridge_manager()->get_bridge_list( bridge_list );
  
  bridge_list.reset();
  int smallest = dynamic_cast<GeometryEntity*>(bridge_list.get())->get_saved_id();
  if (smallest == bte->id())
    return CUBIT_SUCCESS;
  
  for (int i = 1; i < bridge_list.size(); i++)
  {
    int id = dynamic_cast<GeometryEntity*>(bridge_list.next(i))->get_saved_id();
//    if (!id)
//      ;
    if (id == bte->id())
      return CUBIT_SUCCESS;
    else if(id < smallest)
      smallest = id;
  }
  
  if (smallest)
    bte->set_id(smallest);
  
  return CUBIT_SUCCESS;
}


CubitBoolean MergeTool::about_spatially_equal( Surface *surf_1, Surface *surf_2,
                                               double tolerance_factor )
{
  if( surf_1 == surf_2 )
    return CUBIT_TRUE;

  //make sure number of edges between surfaces is identical 
  DLIList<Curve*> surf_1_edges;
  DLIList<Curve*> surf_2_edges;
  surf_1->curves( surf_1_edges );
  surf_2->curves( surf_2_edges );
  if( surf_1_edges.size() != surf_2_edges.size() )
    return CUBIT_FALSE;

  int test_internal = GeometryQueryTool::get_merge_test_internal();

  CubitBoolean test_bbox = GeometryQueryTool::get_merge_test_bbox(); 

  if( test_internal == 2 )
  {
    const GeometryType surf_1_type = surf_1->geometry_type();
    const GeometryType surf_2_type = surf_2->geometry_type();

    if (surf_1_type != SPLINE_SURFACE_TYPE    &&
        surf_1_type != BEST_FIT_SURFACE_TYPE  &&
        surf_1_type != UNDEFINED_SURFACE_TYPE &&
        surf_2_type != SPLINE_SURFACE_TYPE    &&
        surf_2_type != BEST_FIT_SURFACE_TYPE  &&
        surf_2_type != UNDEFINED_SURFACE_TYPE   )
      test_internal = 0;
     else
      test_bbox = CUBIT_FALSE;
  }

  //make sure number of loops between surfaces is identical
  DLIList<LoopSM*> surf_1_loops;
  DLIList<LoopSM*> surf_2_loops;
  surf_1->loopsms( surf_1_loops );
  surf_2->loopsms( surf_2_loops );
  if( surf_1_loops.size() != surf_2_loops.size() )
    return CUBIT_FALSE;

  //compare algnment (sense)
  CubitSense relative_sense = CUBIT_FORWARD;
  CubitVector bbox_center = surf_1->bounding_box().center();
  CubitVector center_point_1;
  surf_1->closest_point( bbox_center, &center_point_1 ); 
  
  bbox_center = surf_2->bounding_box().center();
  CubitVector center_point_2;
  surf_2->closest_point( bbox_center, &center_point_2 ); 

  CubitVector normal_1, normal_2;
  surf_1->closest_point( center_point_1, NULL, &normal_1 );
  surf_2->closest_point( center_point_2, NULL, &normal_2 );
  double dot = normal_1 % normal_2;
  if ( dot < 0 )
     relative_sense = CUBIT_REVERSED;
 
  //compare loops 
  int i,j;
  for( i=surf_1_loops.size(); i--; )
  {
    LoopSM *surf_1_loop = surf_1_loops.get_and_step();
    bool loop_match = false;
    
    //check every loop in surf_2_loops to see if it matches
    //the current loop from surf_1_loops 
    for( j=surf_2_loops.size(); j-- && !loop_match; )
    {
      LoopSM *surf_2_loop = surf_2_loops.get_and_step();
      loop_match = about_spatially_equal( surf_1_loop, surf_2_loop, relative_sense, 
                                          tolerance_factor );

     if( ! loop_match )
       return CUBIT_FALSE;

    }
  }
  
  return CUBIT_TRUE;

}

CubitBoolean MergeTool::about_spatially_equal( LoopSM *loop_1, LoopSM *loop_2, 
                                    CubitSense relative_sense, double tolerance_factor )
{
  
  //get the coedges from the loops
  DLIList<CoEdgeSM*> loop_1_coedges;
  DLIList<CoEdgeSM*> loop_2_coedges;
  loop_1->coedgesms( loop_1_coedges );
  loop_2->coedgesms( loop_2_coedges );

  if( loop_1_coedges.size() != loop_2_coedges.size() )
    return CUBIT_FALSE;

  // Want to compare coedges in order, so make sure we have
  // them in the correct order.
  if (relative_sense == CUBIT_REVERSED)
    loop_1_coedges.reverse();

  // Try to match all coedges.  Begin with the first coedge
  // in this loop.  For each coedge in the other loop that 
  // it matches, check if all the other coedges match in the
  // correct order.
  int other_loop_index = 0;
  loop_1_coedges.reset();
  loop_2_coedges.reset();
  CoEdgeSM* coedge_1 = loop_1_coedges.get_and_step();
  int i;
  for (i = loop_2_coedges.size(); i--; )
  {
    // Loop until we find a matching CoEdge
    CoEdgeSM* coedge_2 = loop_2_coedges.get_and_step();
  
    if (!about_spatially_equal( coedge_1, coedge_2, relative_sense,tolerance_factor ))
      continue;
    
    // Found a matching coedge.  Now try to match all the
    // others in the correct order.
    bool match = true;
    other_loop_index = loop_2_coedges.get_index();
    for (int j = loop_2_coedges.size() - 1; j-- && match; )
    {
      coedge_1 = loop_1_coedges.get_and_step();
      coedge_2 = loop_2_coedges.get_and_step();
      match = about_spatially_equal( coedge_1, coedge_2, relative_sense, tolerance_factor );
    }

    // Matched all coedges, in order.  Done.
    if (match)
      return CUBIT_TRUE;

    // Try again, as perhaps the first coedge of this loop
    // also matches some other one in the second loop and
    // if we start with that one, the remaining coedges will
    // also match.
    loop_1_coedges.reset();
    coedge_1 = loop_1_coedges.get_and_step();
    loop_2_coedges.reset();
    loop_2_coedges.step( other_loop_index );
  }

  // If here, loops didn't match.
  return CUBIT_FALSE;

}
                                        
CubitBoolean MergeTool::about_spatially_equal( CoEdgeSM *coedge_1, CoEdgeSM *coedge_2,
                                               CubitSense relative_sense,
                                               double tolerance_factor )
{
  DLIList<Curve*> curves_1; 
  DLIList<Curve*> curves_2; 
  coedge_1->curves( curves_1 );
  coedge_2->curves( curves_2 );

  Curve* curve_1 = curves_1.get();
  Curve* curve_2 = curves_2.get();

  CubitSense edge_sense;

  if( !about_spatially_equal( curve_1, curve_2, edge_sense, tolerance_factor ) ) 
    return CUBIT_FALSE;

  if( curve_1->geometry_type() == POINT_CURVE_TYPE ||
      curve_2->geometry_type() == POINT_CURVE_TYPE )
    return CUBIT_TRUE;

  if (edge_sense == CUBIT_UNKNOWN)
  {
    PRINT_WARNING("Failed to determine relative sense of curves.\n");
    return CUBIT_TRUE;
  }

  bool coedges_reversed = coedge_1->sense() != coedge_2->sense();
  bool want_reversed = edge_sense != relative_sense;
  if (coedges_reversed == want_reversed)
    return CUBIT_TRUE;
  
  return CUBIT_FALSE;
}

CubitBoolean MergeTool::about_spatially_equal( Curve *curve_1, Curve *curve_2,
                                               CubitSense &relative_sense, 
                                               double tolerance_factor )
{
  if( curve_1 == curve_2 )
    return CUBIT_TRUE;

  relative_sense = CUBIT_FORWARD;
  const double ONE_THIRD = 1.0/3.0;

  // Find the point 1/3 along curve_1 
  CubitVector test_point_1, test_point_2;
  if( curve_1->position_from_fraction( ONE_THIRD, test_point_1 ) != CUBIT_SUCCESS )
    return CUBIT_FALSE;

  // See if the 1/3 point on curve_1 lies on curve_2
  if ( curve_2->closest_point_trimmed(test_point_1, test_point_2)
       != CUBIT_SUCCESS )
  {
    return CUBIT_FALSE;
  }

  if ( GeometryQueryTool::instance()->
       about_spatially_equal(test_point_1, test_point_2, tolerance_factor )
       != CUBIT_SUCCESS )
  {
    return CUBIT_FALSE;
  }
  
  CubitVector tangent_1, tangent_2;
  if( curve_1->closest_point(test_point_2, test_point_1, &tangent_1) != CUBIT_SUCCESS )
    return CUBIT_FALSE;
  if( curve_2->closest_point(test_point_1, test_point_2, &tangent_2) != CUBIT_SUCCESS )
    return CUBIT_FALSE;

  //If one of the curves is zero-length, it will have a zero
  //tangent vector.
  double len_product = tangent_1.length() * tangent_2.length();
  if( len_product > CUBIT_DBL_MIN )
  {
    double dot_product = (tangent_1 % tangent_2);
    if (dot_product < 0)
    relative_sense = CUBIT_REVERSED;
  }
  else
  {
    //If one of the tangents is zero-length, one of the curves had
    //better be as well.
    assert( (curve_1->measure() * curve_2->measure()) < CUBIT_RESABS );
  }

  //compare the start and end vertices to be spatially equal.
  DLIList<Point*> curve_1_points(2), curve_2_points(2);
  curve_1->points( curve_1_points );
  curve_2->points( curve_2_points );

  if( curve_1->bridge_sense() == CUBIT_REVERSED )
    curve_1_points.reverse();
  if( curve_2->bridge_sense() == CUBIT_REVERSED )
    curve_2_points.reverse();

  Point* curve_1_start = curve_1_points.get(); 
  curve_1_points.last();
  Point* curve_1_end =  curve_1_points.get();

  Point* curve_2_start = curve_2_points.get(); 
  curve_2_points.last();
  Point* curve_2_end =  curve_2_points.get();

  if (relative_sense == CUBIT_REVERSED)
    std::swap(curve_2_start, curve_2_end);

  if (curve_1_start == curve_1_end ||curve_2_start == curve_2_end)
  {
    if ((curve_1_start != curve_1_end) ||
        (curve_2_start != curve_2_end) ||
        !about_spatially_equal(curve_1_start, curve_2_start, tolerance_factor))
      return CUBIT_FALSE;
  }
  else
  {
    if ((curve_1_start == curve_2_end) ||
        (curve_1_end == curve_2_start) ||
        !about_spatially_equal(curve_1_start, curve_2_start, tolerance_factor ) ||
        !about_spatially_equal(curve_1_end, curve_2_end, tolerance_factor ))
      return CUBIT_FALSE;
  }

  return CUBIT_TRUE;

}

CubitBoolean MergeTool::about_spatially_equal( Point *point_1, Point *point_2,
                                               double tolerance_factor )
{
  if( point_1 == point_2 )
    return CUBIT_TRUE;

  CubitVector point_1_pos = point_1->coordinates();
  CubitVector point_2_pos = point_2->coordinates();
  
  if (!GeometryQueryTool::instance()-> about_spatially_equal(
     point_1_pos, point_2_pos, tolerance_factor))
    return CUBIT_FALSE;

   return CUBIT_TRUE;
}


