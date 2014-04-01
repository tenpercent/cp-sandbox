//-------------------------------------------------------------------------
// Filename      : MergeToolAssistant.hpp
//
// Purpose       : Abstract class defining methods for doing other 
//                 other operations as a part of merging, for example
//                 handling mesh.
//
// Special Notes : Creating an instance of this class automagically
//                 registers it with MergeTool.  Similarly destroying
//                 the instance will unregister it, and destroying
//                 the MergeTool instance will destroy all registered
//                 MergeToolAssistants.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/26/01
//-------------------------------------------------------------------------

#ifndef MERGE_TOOL_ASSISTANT_HPP
#define MERGE_TOOL_ASSISTANT_HPP

#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;
class RefFace;
class RefEdge;
class RefVertex;

class CUBIT_GEOM_EXPORT MergeToolAssistant
{

  public:
  
    MergeToolAssistant();
    //- Constructor
    //- Registers this assistant with MergeTool.
    
    
    virtual ~MergeToolAssistant();
    //- Destructor
    //- Unregister this assistant with MergeTool.
    
    
    virtual CubitString class_identifier() = 0;
    //- Return a CubitString for MergeTool to use in error messages.
    
    virtual CubitBoolean 
    can_merge( RefFace* keep_face_ptr, RefFace* dead_face_ptr ) = 0;
    //- Allow MergeToolAssistant to abort the merge.  The
    //- implementation of MergeToolAssistant is expected
    //- to check the mergeability of all child RefEdges and
    //- RefVertices as well.  The merge partner for each 
    //- child entity can be determined by calling 
    //- RefEntity::get_compare_partner().  If get_compare_partner()
    //- returns NULL, this should indicate that the child entities
    //- are already merged.
    
    virtual CubitBoolean 
    can_merge( RefEdge* keep_edge_ptr, RefEdge* dead_edge_ptr ) = 0;
    //- Allow MergeToolAssistant to abort the merge.  The
    //- implementation of MergeToolAssistant is expected
    //- to check the mergeability of all child RefVertices. 
    //- The merge partner for each child entity can be determined 
    //- by calling RefEntity::get_compare_partner().  If 
    //- get_compare_partner() returns NULL, this should indicate 
    //- that the child entities are already merged.
    
    virtual CubitBoolean 
    can_merge( RefVertex* keep_vtx_ptr, RefVertex* dead_vtx_ptr ) = 0;
    //- Allow MergeToolAssistant to abort the merge. 
    
    virtual void
    merging( RefEntity* keep_entity_ptr, 
             RefEntity* dead_entity_ptr,
             CubitBoolean opposite_sense ) = 0;
    //- This method is called to notify an assistant of a merge in
    //- progress.  This method is called for each RefEntity merged.
    //- For example when merging two RefEdges, this method will
    //- first be called for each RefVertex as it is merged, and
    //- then for the RefEdge.  This method is called after all 
    //- child entities have been merged, but just prior to merging
    //- the passed entites.  For example if two RefFaces are being
    //- merged, when this method is called for those RefFaces, all
    //- child RefEdges, CoEdges and Loops have already been merged.  
    //- The RefFaces will SHARE A COMMON SET OF LOOPS.
    
    virtual void
    unmerged( RefEntity* old_entity_ptr, 
              RefEntity* new_entity_ptr,
              CubitBoolean reversed ) = 0;
    //I old_entity_ptr 
    //I- A pointer to the RefEntity that existed prior to the unmerge.
    //I new_entity_ptr 
    //I- A pointer to the new RefEntity created as a
    //I- GeometryEntity was split off of old_entity_ptr.
    //I old_sense_switched
    //I- True of the geometric sense of the old entity was changed
    //I- as a part of the unmerge.
    //I new_sense_switched
    //I- True if the geometric sense of the old entity was not changed
    //I- and the new entity has the opposite sense of the old.
    //- Notifies MergeToolAssistant that an unmerge has occured.
    //- This method will be called for all unmerged RefEntities.  For
    //- example, when a RefEdge is unmerged, the method will first be
    //- called for the RefEdge while the child vertices are still 
    //- merged.  It will then be called for the child RefVertices 
    //- if/when they are unmerged.
    
    virtual void finish_merge()   = 0;
    //- Let assisitant know that merge operation is completed
    //- in case any temporary data neds to be cleaned up.
    
    virtual void finish_unmerge() = 0;
    //- Final call to clean up any data.  Called just before
    //- MergeTool returns control to its caller.
    
};

#endif

