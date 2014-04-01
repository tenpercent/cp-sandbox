//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : BodySM.cc
//
// Purpose       : To implement the Body solid model class.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//
// Owner         : Stephen J. Verzi
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

#include "BodySM.hpp"
#include "Body.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
BodySM::BodySM()
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
BodySM::~BodySM()
{
}

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
const type_info& BodySM::topology_entity_type_info() const
{ return typeid(Body); }

CubitBox BodySM::bounding_box()
{
  DLIList<Lump*> body_lumps;
  this->lumps( body_lumps );

  CubitBox bbox;
  bbox = body_lumps.get_and_step()->bounding_box();

  int i;
  for( i=body_lumps.size()-1; i--; )
  {
    Lump *tmp_lump = body_lumps.get_and_step();
    bbox |= tmp_lump->bounding_box();
  }

  return bbox;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

