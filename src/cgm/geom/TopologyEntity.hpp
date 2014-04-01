//-------------------------------------------------------------------------
// Filename      : TopologyEntity.hpp
//
// Purpose       : This file contains the declarations of the class 
//                 TopologyEntity.
//                 This class is the base class of all the topological 
//                 entities in the Virtual Geometry Interface subsystem.
//
// Special Notes : This is a pure virtual class.
//
// Creator       : Malcolm Panthaki 
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef TOPOLOGY_ENTITY_HPP
#define TOPOLOGY_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "ModelEntity.hpp"
#include "BridgeManager.hpp"

#include "CubitEventDefines.h"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefEntity ;
template <class X> class DLIList ;
class RefVolume ;
class RefFace ;
class RefEdge ;
class RefVertex ;
class Shell;
class Loop;
class Chain;
class CoVolume;
class CoFace;
class CoEdge;
class CoVertex;
class GeometryQueryEngine ;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class Body;
class RefGroup;

// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

class CUBIT_GEOM_EXPORT TopologyEntity: public ModelEntity
{
public :
  
  TopologyEntity() ;
    //- The default constructor
  
  virtual ~TopologyEntity();
    //- A pure virtual destructor - this ensures that the class is
    //- non-instantiable.

  static const char* get_class_name()
     { return "TopologyEntity"; }
  
  virtual const char* class_name() const
     { return get_class_name(); }
  
  GeometryQueryEngine* get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometry query engine associated with
    //R- the object.
    //- This function returns a pointer to the geometry query engine
    //- associated with the object.  If this TE has been merged, the
    //- merged TEs may have different GQEs.  In that case, the GQE used
    //- for evaluations is returned.
  
  CubitBoolean is_directly_related( TopologyEntity* entity );
    //- Return TRUE if entity is this or a direct parent or child of this.
  
  CubitStatus bodies(DLIList<Body*>& body_list);
  CubitStatus shells(DLIList<Shell*>& shell_list);
  CubitStatus loops( DLIList<Loop*>& loop_list );
  CubitStatus chains(DLIList<Chain*>& chain_list);
  CubitStatus ref_volumes(DLIList<RefVolume*>& ref_volume_list);
  CubitStatus ref_faces(DLIList<RefFace*>& ref_face_list);
  CubitStatus ref_edges(DLIList<RefEdge*>& ref_edge_list);
  CubitStatus ref_vertices(DLIList<RefVertex*>& ref_vertex_list);
  CubitStatus co_volumes(DLIList<CoVolume*>& co_volume_list);
  CubitStatus co_faces(DLIList<CoFace*>& co_face_list);
  CubitStatus co_edges(DLIList<CoEdge*>& co_edge_list);
  CubitStatus co_vertices(DLIList<CoVertex*>& co_vertex_list);

  RefVertex *ref_vertex();
  RefEdge *ref_edge();
  RefFace *ref_face();
  RefVolume *ref_volume();
  Body *body();
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //O body/ref_volume/ref_face/ref_edge/ref_vertex/shell/loop/chain/covolume/
    //O coface/coedge/covertex _list
    //O- List of Bodies/RefVolumes/RefFaces/RefEdges/RefVertex'es/Loops/
    //O- Shells/Chains/CoVolumes/CoFaces/CoEdge/CoVertices
    //- Append the appropriate list of TopologyEntities associated with
    //- this TopologyEntity.  
    //- There are 2 types of queries involved:
    //- a) Return all "owned" entities such as all the RefFaces owned by
    //-    a RefVolume.
    //- b) Return all "owning" entities such as all the Refvolumes that
    //-    own this RefFace
    //- CUBIT_FAILURE is returned if this entity has the same type as
    //- the requested list - e.g., if ref_edges is called on a RefEdge
    //- object. 
    //-

#ifdef BOYD14
  Shell *shell();
  Loop *loop();
  Chain *chain();
  CoVolume *co_volume();
  CoFace *co_face();
  CoVertex *co_vertex();
    //- get a single connected entity of the given type

  int num_bodies();
  int num_shells();
  int num_chains();
#endif
  CoEdge *co_edge();
  int num_loops();
  int num_ref_volumes();
  int num_ref_faces();
  int num_ref_edges();
  int num_ref_vertices();
#ifdef BOYD14
  int num_co_volumes();
  int num_co_faces();
  int num_co_edges();
  int num_co_vertices();
    //- return the number of connected entities of the given type
#endif

  const BridgeManager* bridge_manager() const
    { return &bridgeMan; }
  BridgeManager* bridge_manager()
    { return &bridgeMan; }
/*  
  virtual EntityType get_topology_bridge_type() const
    { return TopologyBridge_TYPE; }
    // Returns the generic TB type associated with this TE.
*/  
  virtual void reverse_topology() {}
    // This method is used by owned BridgeManagers to
    // notify their owning TopologyEntity when the sense
    // of the geometry has changed (due to a different
    // TopologyBridge being used as the primary rep.)

protected:

  CubitStatus set_topology_bridge(TopologyBridge* TB_ptr);
    // Sets the TB for this TE to be TB_ptr ONLY (removes other TBs).
  
private:
  BridgeManager bridgeMan;
  
  TopologyEntity( const TopologyEntity& );
  void operator=( const TopologyEntity& );
};

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

#endif

