//-------------------------------------------------------------------------
// Filename      : BridgeManager.hpp
//
// Purpose       : Manages the TopologyBridges being used by
//                 a single TopologyEntity.  Encapsulates the
//                 merging/unmerging of TopologyBridges into/from
//                 a TopologyEntity.
//
// Creator       : Darryl Melander
//
// Creation Date : 02/20/99
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------
#ifndef BRIDGE_MANAGER_HPP
#define BRIDGE_MANAGER_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "TBOwner.hpp"
#include "TopologyBridge.hpp"

#include <vector>

class TopologyEntity;
class GeometryQueryEngine;

class CUBIT_GEOM_EXPORT BridgeManager : public TBOwner
{
public:
    // Constructor/Destructor
  BridgeManager(TopologyEntity* owner);
  virtual ~BridgeManager();
  
    // Bridge Management Functions
  virtual CubitStatus add_bridge(TopologyBridge* bridge);
  CubitStatus add_bridge_as_representation(TopologyBridge* bridge);
  virtual CubitStatus remove_bridge(TopologyBridge* bridge );
  virtual CubitStatus bridge_destroyed( TopologyBridge* bridge );
  CubitStatus remove_all_bridges();
  CubitStatus merge (BridgeManager* dead_manager, 
                     CubitSense relative_sense );
  virtual CubitStatus swap_bridge( TopologyBridge* old_tb, 
                                   TopologyBridge* new_tb,
                                   bool reversed );
  
    // Bridge Query Functions
  virtual CubitBoolean contains_bridge(TopologyBridge* bridge) const;
  int get_bridge_list(DLIList<TopologyBridge*>& bridge_list) const;
  TopologyBridge* topology_bridge() const
    {
      if(!mergeList.empty())
        return *(mergeList.begin());
      else
        return 0;
    }
		
    //Bridge Query Functions, querying for bridges owned by
    //a particular GeometryQueryEngine.
  int get_bridge_list(DLIList<TopologyBridge*>& bridge_list,
                      GeometryQueryEngine* gme_ptr ) const;
  TopologyBridge* topology_bridge( GeometryQueryEngine* gme_ptr ) const;
  int number_of_bridges() const;
  
    // TE Query Function (Can't change TE after creation)
  TopologyEntity* topology_entity() const
    {return topologyEntity;}

  void reverse_bridge_senses();
  // Switch the relative sense of all TopologyBridges.
  // If update_first is true, change the principal bridge
  // to be one that has a forward relative sense, if there
  // is one.
  
  virtual void notify_reversed( TopologyBridge* bridge );

private:

  TopologyEntity* const topologyEntity;
  std::vector<TopologyBridge*> mergeList;
};

#endif

