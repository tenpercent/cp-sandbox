//-------------------------------------------------------------------------
// Filename      : MergeToolAssistant.cpp
//
// Purpose       : Implementation of concrete members of MergeToolAssistant
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/26/01
//-------------------------------------------------------------------------

#include "MergeToolAssistant.hpp"
#include "MergeTool.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor - register assistant with MergeTool
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/26/01
//-------------------------------------------------------------------------
MergeToolAssistant::MergeToolAssistant()
{
}

//-------------------------------------------------------------------------
// Purpose       : Destructor - unregister assistant with MergeTool
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/26/01
//-------------------------------------------------------------------------
MergeToolAssistant::~MergeToolAssistant()
{
  MergeTool::instance()->remove_merge_tool_assistant( this );
}

