#ifndef OLD_UNMERGE_CODE_HPP
#define OLD_UNMERGE_CODE_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "UnMergeEvent.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;
class Body;
class RefVolume;
class RefFace;
class Loop;
class RefEdge;
class RefVertex;

class TopologyBridge;
class Surface;
class LoopSM;
class Curve;
class Point;


class CUBIT_GEOM_EXPORT OldUnmergeCode
{
  public:
  
   static OldUnmergeCode& instance();

   static void initialize_settings();

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
   
   
   static bool get_use_old_unmerge_code();
   static void set_use_old_unmerge_code( bool value );
    // If false (default), functions in this class
    // will just return the result of the new code
    // in MergeTool.  If true, the old unmerge code
    // in this class will be used.
   
   private:
     
   RefFace*   unmerge( RefFace*  face_ptr, RefVolume* parent_ptr );
   RefEdge*   unmerge( RefEdge*  edge_ptr, RefFace*   parent_ptr );
   RefVertex* unmerge( RefVertex* vtx_ptr, RefEdge*   parent_ptr );
     //- Split out the merged entity in entity_ptr that is associated
     //- with the passed parent, remove entity_ptr from parent, attach
     //- new, unmerged entity to parent.
     //-
     //- If the passed entity is not a merged entity, no change is made
     //- and the passed entity is returned.
     //-
     //- If the passed entity is a merged entity, the new entity created
     //- by the unmerge is returned.
     //-
     //- NULL is returned if an error was encountered.  It is an error
     //- if the passed parent_ptr is a merged entity.
   
   RefFace*    split_out_Surface( Surface* surface_ptr, CubitBoolean& reversed );
   RefEdge*    split_out_Curves ( DLIList<Curve*>& curve_list, CubitBoolean& reversed );
   RefVertex*  split_out_Points ( DLIList<Point*>& point_list );
    //- Split out merged entities from owning BridgeManager, and
    //- construct new topology, including new child GroupingEntities
    //- and child SenseEntities. 
    //-
    //- NOTE: All passed TopologyBridges MUST belong to the same BridgeManager.                        

   Loop*       split_out_Loop   ( LoopSM* loopsm, RefFace* new_loop_owner,
                                  CubitBoolean reverse );
    //- Unmerge a loop and its child coedges.  Loop and coedge directions
    //- for new entities are reversed if reverse == CUBIT_TRUE.
    
   void cleanup_unmerge();
    //- Does post-processing for unmerge (sending events and such.)
    
  static void remove_CAEntityId_attrib( TopologyBridge* tb_ptr );
    
  
  static void find_curves( Point* point_ptr, DLIList<Curve*>& result_set );
  static void find_surfaces( Curve* curve_ptr, DLIList<Surface*>& result_set);
    //- Virtual Geometry will provide all downward topology bridge queries,
    //- but not always the corresponding upwards queries of the topology
    //- bridges.  These methods are an implementation of the upward queries
    //- from the downward ones.  While somewhat more expensive, this works
    //- around the shortcomming of the VG topology bridge traversals.
    
  DLIList<RefEntity*> new_unmerged;
  DLIList<RefEntity*> old_unmerged;
  DLIList<RefEntity*> unmerge_modified;
  DLIList<UnMergeEvent*> event_list;
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

  OldUnmergeCode() : unmerged_list_in_use(false) {}
  
  static bool useOldUnmergeCode;
  
};

//-------------------------------------------------------------------------
// Purpose       : Start/stop unmerge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/18/01
//-------------------------------------------------------------------------
inline CubitBoolean OldUnmergeCode::start_unmerge()
{ 
  if( unmerged_list_in_use ) return CUBIT_FALSE;
  new_unmerged.clean_out();
  old_unmerged.clean_out();
  unmerge_modified.clean_out();
  unmerged_list_in_use = CUBIT_TRUE;
  return CUBIT_TRUE;
}
inline void OldUnmergeCode::end_unmerge( CubitBoolean top )
{
  if( top )
  {
    unmerged_list_in_use = CUBIT_FALSE;
    cleanup_unmerge();
  }
} 


#endif
