//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Chain.C
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

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "Chain.hpp"

#include "SenseEntity.hpp"
#include "CoVertex.hpp"

#include "BasicTopologyEntity.hpp"

#include "DLIList.hpp"

#include "CastTo.hpp"

#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "ModelEntity.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
Chain::Chain()
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
Chain::~Chain()
{
}

//-------------------------------------------------------------------------
// Purpose       : These functions return the RefVertex associated with 
//                 the first/last CoVertex'es of this Chain.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
RefVertex* Chain::start_vertex()
{
    // Get first CoVertex in chain
  SenseEntity* se_ptr = get_first_sense_entity_ptr();
  if (!se_ptr)
    return 0;
  
    // Get RefVertex from CoVertex
  BasicTopologyEntity* bte_ptr = se_ptr->get_basic_topology_entity_ptr();
  return dynamic_cast<RefVertex*>(bte_ptr);
}

RefVertex* Chain::end_vertex()
{
    // Get last CoVertex in chain
  SenseEntity* se_ptr = get_last_sense_entity_ptr();
  if (!se_ptr)
    return 0;
  
    // Get RefVertex from CoVertex
  BasicTopologyEntity* bte_ptr = se_ptr->get_basic_topology_entity_ptr();
  return dynamic_cast<RefVertex*>(bte_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : These functions return the 
//                 first/last CoVertex'es of this Chain.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
CoVertex* Chain::start_co_vertex()
{
    // Get first CoVertex in chain
  SenseEntity* se_ptr = get_first_sense_entity_ptr();
  return dynamic_cast<CoVertex*>(se_ptr);
}

CoVertex* Chain::end_co_vertex()
{
    // Get last CoVertex in chain
  SenseEntity* se_ptr = get_last_sense_entity_ptr();
  return dynamic_cast<CoVertex*>(se_ptr);
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

// ********** END PUBLIC FUNCTIONS         **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

