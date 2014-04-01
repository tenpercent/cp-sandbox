//-------------------------------------------------------------------------
// Filename      : TBOwnerSet.hpp
//
// Purpose       : Interface for an adapter allowing a TopologyBridge
//                 to have multiple owners.  I.e. An object implementing
//                 this interface is the owner the bridge keeps a pointer
//                 to.  The object implementing this interface maintains   
//                 links to a set of TopologyBridges (partitions) of the
//                 actual geometry.
//
// Special Notes : This interface is provided to abstract the one facet
//                 of partition geometry that TopologyBridge needs to be
//                 aware of.  
//
//                 The owners contained in a TBOwnerSet need not (and
//                 probably will not) implement the TBOwner interface.
//                 They must be TopologyBridges.  (This doesn't make a whole
//                 lot of sense, beyond the fact that that is what is
//                 required for partition geometry.)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/10/02
//-------------------------------------------------------------------------
#ifndef TB_OWNER_SET
#define TB_OWNER_SET

#ifndef TB_OWNER_HPP
#  include "TBOwner.hpp"
#endif

class TopologyBridge;
template <class X> class DLIList;

class TBOwnerSet : public TBOwner
{

  public:
  
    virtual void get_owners( DLIList<TopologyBridge*>& owner_list ) const = 0;
    
    virtual int get_owner_layer() const = 0;
};

#endif

