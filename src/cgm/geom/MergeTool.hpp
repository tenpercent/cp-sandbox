//- Class:       MergeTool
//- Description: Location of all Merge and Unmerge functionality
//- Owner:       Steven Jankovich
//- Created: 27 April 1998
//- Checked By:  
//- Version: $Id: 

#ifndef MERGETOOL_HPP
#define MERGETOOL_HPP 

#include "CubitDefines.h"
#include "CubitEventDefines.h"
#include "GeometryDefines.h"
#include "CpuTimer.hpp"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif

class MergeToolAssistant;

class RefEntity;
class TopologyEntity;
class BasicTopologyEntity;
class GroupingEntity;
class SenseEntity;
class TopologyBridge;
class GeometryEntity;

class RefFace;
class RefEdge;
class RefVertex;
class Loop;
class CoEdge;

class BodySM;
class Body;
class RefEntity;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class RefGroup;
//template <class X> class DLIList;
class GeometryEntity;
class CAMergePartner;
class UnMergeEvent;

class Surface;
class Curve;
class Point;
class LoopSM;
class CoEdgeSM;



class CUBIT_GEOM_EXPORT MergeTool
{
   friend class CAMergePartner;
   friend class OldUnmergeCode;
   public :
    
   bool displayProgress;
    // Show progress meter when merging
   
   static MergeTool* instance();
     // Returns a static pointer to unique instance of this class.
   
   ~MergeTool();
     //- Destructor.
     
   CubitBoolean contains_merged_entities( DLIList<Body*> &bodies );
     //- Tests the entities in the list of bodies to see if they 
     //- have been merged with other entities and result in
     //- of multiple volumes.  This returns CUBIT_TRUE at the first
     //- lower entity that shares two volumes.
   
  CubitBoolean contains_merged_children( Body *body,
                                         DLIList<RefEntity*> &merged_children );


  CubitBoolean contains_merged_entities( DLIList<RefEntity*> &ref_entities);
     //- Tests the entities in the list to see if they or their descendents
     //- have been merged with other entities
  
  CubitBoolean parents_contain_merged_entities( DLIList<RefEntity*> &ref_entities);
     //- Tests the entities in the list to see if their ANY children of ancestors 
     //- have been merged with other entities
   
   CubitBoolean entity_merged( TopologyEntity *entity );
     //- Takes a RefEntity, makes sure it is not a volume, and
     //- checks to see if it has more than one volume.
   
   CubitStatus merge_all_bodies();
     //- Compare all RefFaces, RefEdges, and RefVertices in the
     //- model and merge the matches
   
   CubitStatus merge_bodies( DLIList<Body*>& refbody_list );
     //- Compare all RefFaces, RefEdges, and RefVertices in 
     //- refbody_list and merge the matches
   
   CubitStatus merge_volumes( DLIList<RefVolume*>& vol_list,
                              CubitBoolean print_info = CUBIT_TRUE );
     //- Compare all RefFaces, RefEdges, and RefVertices in 
     //- vol_list and merge the matches
   
   CubitStatus merge_all_reffaces    ();
     //- Compare all RefFaces in the model and merge the matches
   
   CubitStatus merge_reffaces_old( DLIList<RefFace*>& refface_list,
                               CubitBoolean print_info = CUBIT_TRUE);
     //- Compare all input RefFaces and merge the matches
   
   CubitStatus merge_reffaces( DLIList<RefFace*>& refface_list,
                               CubitBoolean print_info = CUBIT_TRUE);
     //- Compare all input RefFaces and merge the matches
     //- Uses AbstractTree rather than O(nlogn) comparisons.

   CubitStatus merge_all_refedges();
     //- Compare all RefEdges in the model and merge the matches
   
   CubitStatus old_merge_refedges( DLIList<RefEdge*>& refedge_list,
                               CubitBoolean should_clean_out = CUBIT_TRUE,
                               CubitBoolean print_info = CUBIT_TRUE);
   CubitStatus merge_refedges( DLIList<RefEdge*>& refedge_list,
                               CubitBoolean should_clean_out = CUBIT_TRUE,
                               CubitBoolean print_info = CUBIT_TRUE);
     //- Compare all input RefEdges and merge the matches
     //- BE CAREFUL with the should_clean_out flag.  If you set
     //- it to false, then YOU (the caller) are responsible for
     //- cleaning out the deactivated geometry.
   
   CubitStatus merge_all_refvertices();
     //- Compare all RefVertices in the model and merge the matches
   
   CubitStatus old_merge_refvertices( DLIList<RefVertex*>& refvertex_list,
                                  CubitBoolean print_info = CUBIT_TRUE );
   CubitStatus merge_refvertices( DLIList<RefVertex*>& refvertex_list,
                                  CubitBoolean print_info = CUBIT_TRUE );
     //- Compare all input RefVertices and merge the matches
   
   CubitStatus merge_entities( DLIList<RefEntity*>& entity_list,
                               CubitBoolean should_clean_out = CUBIT_TRUE,
                               CubitBoolean print_info = CUBIT_TRUE);
     //- merge the entities in list; asserts if they're not all the same type
   
   CubitStatus unmerge_all();
     //- Unmerge everything.
   
   CubitStatus 
   unmerge( DLIList<RefEntity*>& entity_list, CubitBoolean descend = CUBIT_TRUE );
     //- Unmerge entities in list.
     //- If decend is true, will decend topology graph, unmerging child topology
     //- of the passed topology.  If decend is false, a.) passing bodies or volumes
     //- in entity list will have no effect and b.) when a surface is unmerged
     //- its child curves will not be unmerged, and the child vertices will not be 
     //- unmerged when an edge is unmerged.
     
   CubitStatus unmerge( RefEntity* entity_ptr, CubitBoolean descend = CUBIT_TRUE );
   CubitStatus unmerge( RefFace*   face_ptr,   CubitBoolean descend = CUBIT_TRUE );
   CubitStatus unmerge( RefEdge*   edge_ptr,   CubitBoolean descend = CUBIT_TRUE );
   CubitStatus unmerge( Body*      body_ptr   );
   CubitStatus unmerge( RefVolume* vol_ptr    );
   CubitStatus unmerge( RefVertex* vertex_ptr );
     //- Unmerge the passed entity.  All parents must already be
     //- unmerged.  
     //- If decend is true, will decend topology graph, unmerging child topology
     //- of the passed topology.  If decend is false, a.) passing bodies or volumes
     //- in entity list will have no effect and b.) when a surface is unmerged
     //- its child curves will not be unmerged, and the child vertices will not be 
     //- unmerged when an edge is unmerged.
     
   int unmerged_entities( DLIList<RefEntity*>* entities = NULL ) const;
     //- Return number of entities unmerged in last call to one
     //- of the unmerge methods.  If an entity list is passed,
     //- it will be populated with unmerged entities.
     
    CubitStatus separate_bodies( DLIList<Body*>& separate_list,
                                 DLIList<Body*>* from_list = NULL );
      //- Unmerge such that the group of entities in
      //- separate_list share no topology with the entities
      //- in from_list, or if from_list is is null, any other
      //- entities.
   
    CubitStatus separate_volumes( DLIList<RefVolume*>& separate_list,
                                  DLIList<RefVolume*>* from_list = NULL );
      //- Unmerge such that the group of entities in
      //- separate_list share no topology with the entities
      //- in from_list, or if from_list is is null, any other
      //- entities.
   
    CubitStatus separate_faces( DLIList<RefFace*>& separate_list,
                                DLIList<RefFace*>* from_list = NULL );
      //- Unmerge such that the group of entities in
      //- separate_list share no topology with the entities
      //- in from_list, or if from_list is is null, any other
      //- entities.
  
    CubitStatus separate_edges( DLIList<RefEdge*>& separate_list,
                                DLIList<RefEdge*>* from_list = NULL );
      //- Unmerge such that the group of entities in
      //- separate_list share no topology with the entities
      //- in from_list, or if from_list is is null, any other
      //- entities.
   
   RefFace* force_merge( RefFace* face1, RefFace* face2 );
   RefEdge* force_merge( RefEdge* edge1, RefEdge* edge2 );
   RefVertex* force_merge( RefVertex* vtx1, RefVertex* vtx2 );
   RefEntity* force_merge( RefEntity* ent1, RefEntity* ent2 );
   RefEntity* force_merge( const DLIList<RefEntity*>& list );
   
   static CubitBoolean merge_has_occured();
   static void set_merge_occurance( CubitBoolean t_or_f );
     //- Sets/Gets the mergeCalled flag.  This signifies
     //- that currently a merge operation has been called for 
     //- the model.

   static void group_results( CubitBoolean t_or_f );
    //- tells us to group the results or not.
    
   static void set_new_ids_on_unmerge( CubitBoolean value );
   static CubitBoolean get_new_ids_on_unmerge();
   
   static void initialize_settings();
  
  void compare_notify(RefEntity *entity, CubitEventType event);
    //- notifies MergeTool about comparisons found and put on ref entities
  
  void remove_compare_data();
    //- Remove TDCompares from RefEntities.
    
    void add_merge_tool_assistant( MergeToolAssistant* mta_ptr );
    void remove_merge_tool_assistant( MergeToolAssistant* mta_ptr );
    MergeToolAssistant* find_merge_tool_assistant( const type_info& type );
  
  static void destroy_dead_geometry( CubitBoolean yes_no )
    { destroyDeadGeometry = yes_no; }

  RefGroup* get_group_last_merged_surfs()
    {
        //clear this out once this function has been called.
      RefGroup *temp = lastSurfsMerged;
      lastSurfsMerged = NULL;
      return temp;
    }
    ///
    /// This function is to be used only right
    /// after merging is called.  It is a way to
    /// access the groups that are created during the
    /// merge.  Note the groups can be destroyed
    /// and these pointers can be stale if the
    /// function is not called immediatly following merging...
    ///
  RefGroup* get_group_last_merged_curvs()
    {
        //clear this out once this function has been called.
      RefGroup *temp = lastCurvsMerged;
      lastCurvsMerged = NULL;
      return temp;
    }
    ///
    /// This function is to be used only right
    /// after merging is called.  It is a way to
    /// access the groups that are created during the
    /// merge.  Note the groups can be destroyed
    /// and these pointers can be stale if the
    /// function is not called immediatly following merging...
    ///
  RefGroup* get_group_last_merged_verts()
    {
        //clear this out once this function has been called.
      RefGroup *temp = lastVertsMerged;
      lastVertsMerged = NULL;
      return temp;
    }
    ///
    /// This function is to be used only right
    /// after merging is called.  It is a way to
    /// access the groups that are created during the
    /// merge.  Note the groups can be destroyed
    /// and these pointers can be stale if the
    /// function is not called immediatly following merging...
    ///

  //Note:  The caller of the following 4 functions is responsible to 
  //delete the DLIList*s that are returned.
   CubitStatus find_mergeable_refentities( DLIList<RefEntity*> &entities,
              DLIList< DLIList<RefFace*>*> &lists_of_mergeable_ref_faces,
              DLIList< DLIList<RefEdge*>*> &lists_of_mergeable_ref_edges,
              DLIList< DLIList<RefVertex*>*> &lists_of_mergeable_ref_vertices);
   CubitStatus find_mergeable_reffaces( DLIList<RefEntity*> &entities,
              DLIList< DLIList<RefFace*>*> &lists_of_mergeable_ref_faces,
              bool clean_up_compare_data = true );
   CubitStatus find_mergeable_refedges( DLIList<RefEntity*> &entities,
              DLIList< DLIList<RefEdge*>*> &lists_of_mergeable_ref_edges,
              bool clean_up_compare_data = true );
   CubitStatus find_mergeable_refvertices( DLIList<RefEntity*> &entities,
              DLIList< DLIList<RefVertex*>*> &lists_of_mergeable_ref_vertices,
              bool clean_up_compare_data = true );

   //Faster comparison that only check for mergeable refedges between volumes.
   //It reports all edges, even when the owning face is mergeable
   CubitStatus find_only_mergeable_curves( DLIList<BodySM*> &body_list, 
                  DLIList< DLIList<Curve*>*> &lists_of_mergeable_curves );

   CubitStatus find_only_mergeable_surfaces( DLIList<BodySM*> &body_list, 
                  DLIList< DLIList<Surface*>*> &lists_of_mergeable_surfaces);

   CubitStatus find_only_mergeable_refedges( DLIList<Body*> &body_list, 
                  DLIList< DLIList<RefEdge*>*> &lists_of_mergeable_ref_edges );
                                             
   CubitBoolean about_spatially_equal( Surface *surf_1, Surface *surf_2,
                                       double tolerance_factor );
   CubitBoolean about_spatially_equal( LoopSM *loop_1, LoopSM *loop_2, 
                                       CubitSense relative_sense, double tolerance_factor );
   CubitBoolean about_spatially_equal( CoEdgeSM *coedge_1, CoEdgeSM *coedge_2, 
                                       CubitSense relative_sense, double tolerance_factor );
   CubitBoolean about_spatially_equal( Curve *curve_1, Curve *curve_2, 
                                       CubitSense &relative_sense, double tolerance_factor );
   CubitBoolean about_spatially_equal( Point *point_1, Point *point_2, 
                                               double tolerance_factor );
  
   protected :
   
   CubitStatus separate_entities( DLIList<TopologyEntity*>& separate_list,
                                  DLIList<TopologyEntity*>* from_list = NULL );
    //- Common implementation for public separate() functions.
    //- All passed entities must be of the same type.
    //- Unmerge such that the group of entities in
    //- separate_list share no topology with the entities
    //- in from_list, or if from_list is is null, any other
    //- entities.
   
   BasicTopologyEntity* can_separate( DLIList<TopologyBridge*>& bridge_list,
                                      bool check_parents );
    //- Check if the passed list of bridges can be unmerged and if
    //- so, return their owning BTE.  If check_parents is false, skip
    //- check for merged parents.
  
   RefFace* separate_face( DLIList<Surface*>& bridges, bool descend );
    //- Split a merged entity into two such that the returned, new
    //- entity contains the passed list of bridges.  If descend
    //- is false, skip attempt to unmerge child entities.
  
   RefEdge* separate_edge( DLIList<Curve*>& bridges, bool descend );
    //- Split a merged entity into two such that the returned, new
    //- entity contains the passed list of bridges.  If descend
    //- is false, skip attempt to unmerge child entities.
   
   RefVertex* separate_vertex( DLIList<Point*>& bridges );
    //- Split a merged entity into two such that the returned, new
    //- entity contains the passed list of bridges. 
    
   void cleanup_unmerge();
    //- Does post-processing for unmerge (sending events and such.)
    
   CubitStatus check_saved_id( BasicTopologyEntity* bte );
    //- Used as part of unmerging.  If the original geometry
    //- entity corresponding to the ID of the passed BTE has
    //- been removed/unmerged from the BTE, assign the BTE a
    //- new ID from it's current entities so that the new,
    //- unmerged entity can be assigned the current ID.

    
   private :
   
   static void warn_about_refface_sense( RefFace* face_ptr_1,
                                         RefFace* face_ptr_2,
                                         bool faces_reversed );
   
   static CubitBoolean destroyDeadGeometry;
   
   static MergeTool* instance_;
     // Static pointer to unique instance of this class
   
   static CubitBoolean mergeCalled;
     //- This static flag tells us if a merge operation on the
     //- current model has been performed.  After a FULL unmerge,
     //- this should be reset, as well after a "reset" to be FALSE.

  RefGroup *lastSurfsMerged;
  RefGroup *lastCurvsMerged;
  RefGroup *lastVertsMerged;
    ///
    /// Contains pointers to the groups storing the most
    /// recently merged surfaces, curves and vertices.
    ///
  
  DLIList<MergeToolAssistant*> assistant_list_;
    //- List of helper objects to handle updates outside the realm
    //- of CGM, like updating mesh.

   static CubitBoolean groupResults;
    //- Tells the tool to group the merging results.
   
  DLIList<RefEntity*>  compareEntityList;
    //- A list containing entities that have found partners with
    //- whom they compare.

  DLIList<RefEntity*>  mergeSurvivorEntityList;
    //- A list containing entities that were involved in a merge
    //- operation and survived.  
    
  DLIList<RefEntity*> new_unmerged;
  DLIList<RefEntity*> old_unmerged;
  CubitBoolean unmerged_list_in_use;
     //- A list new entities created by unmerging, and a list of
     //- existing entities that were unmerged, and a flag set by the highest
     //- level method called to begin an unmerge.  The flag signifies
     //- that lower-level methods should not clear the list when
     //- starting/finishing, and should not try to update the graphics.
     //-
     //- Any method that sets unmerged_list_in_use to true MUST
     //- set it to false before returning.
     
  CubitBoolean start_unmerge();
  void end_unmerge( CubitBoolean top );
     //- Handle setting unmerged_list_in_use, clearing lists,
     //- and calling cleanup_unmerge().

   MergeTool();
     //- Constructor for the MergeTool object
   
   CubitStatus merge_BTE( BasicTopologyEntity* keeper_entity,
                          BasicTopologyEntity* dead_entity );
     //R CubitStatus
     //R- CUBIT_SUCCESS/FAILURE
     //I keeper_entity, dead_entity
     //I- The BasicTopologyEntities that are to be merged.
     //- This function merges the two BasicTopologyEntities

   CubitBoolean compare_BTE( BasicTopologyEntity * keeper_entity,
                             BasicTopologyEntity * dead_entity ) const;
     //R CubitBoolean
     //R- CUBIT_TRUE/CUBIT_FALSE
     //I basicTopoEntityPtr
     //I- The BTE that "this" BTE must be "compared" to.
     //- This function compares "this" BasicTopologyEntity with
     //- the input BasicTopologyEntity.  This is a spatial comparison.
     //- IMPORTANT NOTE:
     //- The actual spatial comparison of the underlying GeometricEntities 
     //- is NOT done in this routine.  Before this routine is called,
     //- the SolidModelingEngine compare routines would have to be called.
     //- If matches are found there, the matching pairs of BTE's are tagged
     //- with TDCompare objects (containing doubly-linked pointers).  This
     //- function merely checks to see if these TD objects exist on the
     //- "this" and basicTopoEntityPtr objects.  If they do, the function
     //- checks their pointers to make sure they point to each other.
     //- Returns CUBIT_SUCCESS if the comparison was successful.
   
   CubitStatus merge_GE( GroupingEntity* keeper_entity,
                         GroupingEntity* dead_entity );
     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //I keeper_entity, dead_entity
     //I- The GroupingEntities that are to be merged.
     //- This function merges the two GroupingEntities

   CubitBoolean compare_GE( GroupingEntity* keeper_entity,
                            GroupingEntity* dead_entity );
     //R CubitBoolean
     //R- CUBIT_TRUE/CUBIT_FALSE
     //I dead_GroupingEntityPtr
     //I- This GroupingEntity with which comparision is to be done.
     //- This function patially compares this GroupingEntity to the input 
     //- GroupingEntity. The function returns CUBIT_TRUE if the two
     //- objects are spatially equal, CUBIT_FALSE otherwise.
     //- MJP NOTE:
     //- The first trial check that is done is to make sure that the
     //- GroupingEntities have the same number of SenseEntities.
     //- If they do not, then we ASSUME that the GroupingEntities
     //- are not spatially equal.  This is not strictly a good
     //- assumption, but we don't have a more exact algorithm
     //- that takes care of the case where, for example, the number 
     //- of RefEdges associated with two Loops are different, but
     //- the Loops themselves are spatialy equal.  Consider the case
     //- of 2 Loops representing the exact same square -- one can
     //- have 4 RefEdges and the other could have 5 (just bisect one
     //- of the previous ones...)   
   
   CubitStatus compare_and_merge( CubitBoolean merge_flag,
                                  GroupingEntity* keeper_entity,
                                  GroupingEntity* dead_entity );
     //R CubitStatus
     //R- the result of spatial comparison if if_merge is CUBIT_FALSE,
     //R- or the result of merge if if_merge is CUBIT_TRUE
     //I if_merge
     //I- If this flag is CUBIT_TRUE, whenever find the matching SenseEntities,
     //I do the merge. otherwise, skip the merge block.
     //- This private function always do the spatially comparison of "this"
     //- GroupingEntity with the input GroupingEntity, by comparing the
     //- associated SenseEntities of "this" and input GroupingEntity. If all
     //- the SenseEntites associated with "this" GroupingEntity have matching
     //- SenseEntity associated with the input GroupingEntity, then "this" and
     //- the input GroupingEntity are spatially equal. 
     //- SPECIAL NOTES : 
     //  Merging the OSME's and merging links are done in merge()
     //  function. In this function, ONLY the matching SenseEntities
     //  can be merged when asked to merge them ( i.e., if_merge
     //  flag is set to CUBIT_TRUE ).
   
   CubitStatus merge_SE( SenseEntity* keeper_entity,
                         SenseEntity* dead_entity );
     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //I dead_SenseEntityPtr
     //I- A pointer to a SenseEntity which is about to be merged 
     //I- with "this".
     //- This function merges the input dead_SenseEntityPtr with this 
     //- SenseEntity. If the operation can fail if the two entitys are
     //- of different types (derived types.) At the end of a successful
     //- operation, the dead_SenseEntityPtr is an entity that is 
     //- dangling in space with no connection to any other entity. The
     //- return value is CUBIT_SUCCESS when a successful merge occurs,
     //- CUBIT_FAILURE otherwise.

   CubitBoolean compare_SE( SenseEntity* keeper_entity,
                            SenseEntity* dead_entity );
     //R CubitBoolean
     //R- CUBIT_TRUE/CUBIT_FALSE
     //I dead_SenseEntityPtr
     //I- This SenseEntity with which comparision is to be done.
     //- This function patially compares this SenseEntity to the input 
     //- SenseEntity. The function returns CUBIT_TRUE if the two
     //- objects are spatially equal, CUBIT_FALSE otherwise.

  void complete_merge();
    //- completed or abort a merge; clean up TDs and lists

  void test_r_tree(DLIList<RefFace*> &refface_list);
  void test_r_star_tree(DLIList<RefFace*> &refface_list);
  void test_no_tree(DLIList<RefFace*> &refface_list);
    //- temp timing functions for testing with and without and r-tree.
  
};
inline void MergeTool::group_results(CubitBoolean t_or_f)
{groupResults = t_or_f;}
inline CubitBoolean MergeTool::merge_has_occured()
{ return mergeCalled;}
inline void MergeTool::set_merge_occurance( CubitBoolean t_or_f )
{ mergeCalled = t_or_f;}


//-------------------------------------------------------------------------
// Purpose       : Start/stop unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
inline CubitBoolean MergeTool::start_unmerge()
{ 
  if( unmerged_list_in_use ) return CUBIT_FALSE;
  new_unmerged.clean_out();
  old_unmerged.clean_out();
  unmerged_list_in_use = CUBIT_TRUE;
  return CUBIT_TRUE;
}
inline void MergeTool::end_unmerge( CubitBoolean top )
{
  if( top )
  {
    unmerged_list_in_use = CUBIT_FALSE;
    cleanup_unmerge();
  }
} 

#endif

