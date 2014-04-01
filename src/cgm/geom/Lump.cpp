//-------------------------------------------------------------------------
// Filename      : Lump.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#include "Lump.hpp"
#include "RefVolume.hpp"


//-------------------------------------------------------------------------
// Purpose       : The default constructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------

Lump::Lump()
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------

Lump::~Lump()
{}

//-------------------------------------------------------------------------
// Purpose       : Get type of TopologyEntity this GeometryEntity
//                 should be attached to.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/14/01
//-------------------------------------------------------------------------
const type_info& Lump::topology_entity_type_info() const
{ return typeid(RefVolume); }


