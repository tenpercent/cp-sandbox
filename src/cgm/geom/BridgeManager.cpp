//-------------------------------------------------------------------------
// Filename      : BridgeManager.cpp
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
#include "BridgeManager.hpp"
#include "TopologyBridge.hpp"
#include "TopologyEntity.hpp"

#include <algorithm>

// Constructor.
// Attaches this BridgeManager to 'owner'.
BridgeManager::BridgeManager(TopologyEntity* owner)
    : topologyEntity(owner)
{
}

// Destructor
BridgeManager::~BridgeManager()
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=mergeList.begin(); iter!=mergeList.end(); iter++)
  {
    delete (*iter);
  }
  mergeList.clear();
}

// Adds this TB to the child list.  Tells 'bridge'
// that it belongs to this BridgeManager.  If bridge is
// already in this BridgeManager, we return FAILURE.
// This doesn't mean 'bridge' won't be in this manager,
// it just means it won't be in the list twice.
CubitStatus BridgeManager::add_bridge(TopologyBridge* bridge)
{
  if( !bridge || bridge->owner() )
  {
    assert( bridge && !bridge->owner() );
    return CUBIT_FAILURE;
  }

  mergeList.push_back(bridge);
  
  bridge->owner(this);
  return CUBIT_SUCCESS;
}

// This function sets the 
CubitStatus BridgeManager::add_bridge_as_representation(TopologyBridge* bridge)
{
    // Make sure bridge isn't already owned by someone else
  if( !bridge || bridge->owner() ) 
  {
    assert( bridge && !bridge->owner() ); 
    return CUBIT_FAILURE;
  }

  // Add the bridge to the beginning of the list
  mergeList.insert(mergeList.begin(), bridge );
    
  bridge->owner(this);
  return CUBIT_SUCCESS;
}  

// Removes 'bridge' from this.  If it wasn't in this,
// then return FAILURE.
CubitStatus BridgeManager::remove_bridge(TopologyBridge* bridge )
{
  std::vector<TopologyBridge*>::iterator iter;
  iter = std::find(mergeList.begin(), mergeList.end(), bridge);
  if (iter != mergeList.end())
    mergeList.erase(iter);
  else
    return CUBIT_FAILURE;
  
  bridge->owner(0);
  return CUBIT_SUCCESS;
}

CubitStatus BridgeManager::bridge_destroyed( TopologyBridge* bridge )
{
//  if( firstBridge == bridge &&
//      bridge->next_bridge() &&
//      bridge->next_bridge()->bridge_sense() != bridge->bridge_sense() )
//    topologyEntity->reverse_topology();

  return remove_bridge( bridge );
}
  
CubitStatus BridgeManager::remove_all_bridges()
{
  if (mergeList.empty())
    return CUBIT_SUCCESS;

  std::vector<TopologyBridge*>::iterator iter;
  std::vector<TopologyBridge*> temp_list;

  // make a copy of the list so we can modify the real one
  temp_list = mergeList;

  for (iter=temp_list.begin(); iter!=temp_list.end(); iter++)
  {
    if (!remove_bridge( *iter ))
      break;
  }

  return mergeList.empty() ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

// Transfer all bridges from 'dead' to 'this'.
// Set correct pointers.
CubitStatus BridgeManager::merge(BridgeManager* dead_manager, 
                                 CubitSense relative_sense )
{
  if (dead_manager->mergeList.empty())
    return CUBIT_SUCCESS;
  
    // Remove bridge list from dead BridgeManager
  std::vector<TopologyBridge*> dead_list = dead_manager->mergeList;
  dead_manager->mergeList.clear();

    // Set owner and reverse sense if necessary
    // for all bridges in merge list.
  std::vector<TopologyBridge*>::iterator iter;
  for (iter = dead_list.begin(); iter != dead_list.end(); iter++)
  {
    (*iter)->owner(this);
    if( relative_sense == CUBIT_REVERSED )
      (*iter)->reverse_bridge_sense();
  }
  
    // append merge list to end of this manager's bridge list
  mergeList.insert(mergeList.end(), dead_list.begin(), dead_list.end());
  
  return CUBIT_SUCCESS;
}

// Indicates whether this bridge is in this manager.
CubitBoolean BridgeManager::contains_bridge(TopologyBridge* bridge) const
{
  return bridge->owner() == (TBOwner*)this;
}

// Appends all children TBs to bridge_list.
// There is no check for duplicates in the list.
// The number of TBs added is returned.
int BridgeManager::get_bridge_list(
  DLIList<TopologyBridge*>& bridge_list) const
{
  std::vector<TopologyBridge*>::const_iterator iter;
  for (iter = mergeList.begin(); iter != mergeList.end(); iter++)
  {
    bridge_list.append(*iter);
  }
  return mergeList.size();
}


// Get all the child TB's that are owned by the passed
// BridgeManager.
int BridgeManager::get_bridge_list( 
	DLIList<TopologyBridge*>& bridge_list,
  GeometryQueryEngine* gqe_ptr ) const
{
  std::vector<TopologyBridge*>::const_iterator iter;
  for (iter=mergeList.begin(); iter!=mergeList.end(); iter++)
  {
    if ( (*iter)->get_geometry_query_engine() == gqe_ptr )
      bridge_list.append(*iter);
  }
  return bridge_list.size();
  // TODO - BWH - return on this function should be consistent with other bridge list function
}

TopologyBridge* BridgeManager::topology_bridge(
  GeometryQueryEngine* gqe_ptr ) const
{
  std::vector<TopologyBridge*>::const_iterator iter;
  for (iter=mergeList.begin(); iter!=mergeList.end(); iter++)
  {
    if ( (*iter)->get_geometry_query_engine() == gqe_ptr )
      return *iter;
  }
  
  return 0;
}


void BridgeManager::reverse_bridge_senses ()
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=mergeList.begin(); iter!=mergeList.end(); iter++)
  {
    (*iter)->reverse_bridge_sense();
  }
}  

CubitStatus BridgeManager::swap_bridge( TopologyBridge* old_tb,
                                        TopologyBridge* new_tb,
                                        bool reversed )
{
    // make sure new_tb isn't already owned by someone else
  if( new_tb->owner() != NULL )
  {
    assert( new_tb->owner() == 0 );
    return CUBIT_FAILURE;
  }

    // Replace old bridge with new bridge at same location
    // in linked list.

  std::vector<TopologyBridge*>::iterator iter;
  iter = std::find(mergeList.begin(), mergeList.end(), old_tb);
  if (iter != mergeList.end())
  {
    iter = mergeList.erase(iter);
    mergeList.insert(iter, new_tb);
  }
  else
  {
    assert(false);
    return CUBIT_FAILURE;
  }
  
    // Update owner pointers
  
  new_tb->owner(this);
  old_tb->owner(0);
  
    // If 'reversed == true' then make sure bridges have
    // opposite sense.
  
  bool same_sense = new_tb->bridge_sense() == old_tb->bridge_sense();
  if (reversed == same_sense) // if (reversed XOR same_sense) 
    new_tb->reverse_bridge_sense();
    
  return CUBIT_SUCCESS;
}

int BridgeManager::number_of_bridges() const
{
  return mergeList.size();
}

void BridgeManager::notify_reversed( TopologyBridge* bridge )
{
  bridge->reverse_bridge_sense();
}
