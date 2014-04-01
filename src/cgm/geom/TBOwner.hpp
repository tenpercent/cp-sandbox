//-------------------------------------------------------------------------
// Filename      : TBOwner.hpp
//
// Purpose       : Interface for an object that owns a TopologyBridge
//
// Special Notes : Interface implemented by BridgeManager and virtual classes.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/28/01
//-------------------------------------------------------------------------
#ifndef TB_OWNER_HPP
#define TB_OWNER_HPP

#ifndef CUBITOBJECT_HPP
#  include "CubitDefines.h"
#endif

#include "TopologyBridge.hpp"

class CUBIT_GEOM_EXPORT TBOwner
{

  public:

    virtual ~TBOwner() {}
  
    virtual CubitStatus remove_bridge( TopologyBridge* bridge ) = 0;
    
    virtual CubitStatus bridge_destroyed( TopologyBridge* bridge )
      { return remove_bridge( bridge ); }
    
    virtual CubitStatus swap_bridge( TopologyBridge* old_tb, 
                                     TopologyBridge* new_tb,
                                     bool reversed ) = 0;
                                     
    virtual CubitBoolean contains_bridge( TopologyBridge* bridge ) const
      { return static_cast<CubitBoolean>(bridge->owner() == this); }
    
    virtual void notify_reversed( TopologyBridge* ) = 0;

    virtual void notify_joined( TopologyBridge* /* dead_bridge      */,
                                TopologyBridge* /* combined_bridge  */) {}
    
    virtual void notify_split ( TopologyBridge* /* new_bridge       */,
                                TopologyBridge* /* split_from       */) {}
    
    virtual void notify_merged( TopologyBridge*   dead_bridge         ,
                                TopologyBridge* /* coincident_bridge*/) 
      { remove_bridge( dead_bridge ); }
    
    virtual void notify_copied( TopologyBridge* /* new_bridge       */,
                                TopologyBridge* /* split_from       */) {}
    virtual void notify_topology_modified( TopologyBridge* ) {}

    
  protected:
};

#endif

