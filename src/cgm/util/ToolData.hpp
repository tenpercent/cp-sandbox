//- Class: ToolData
//- Owner: Scott Mitchell
//- Description: This is tool specific data that should chain off of CubitNode.
//-              This is the generic form, all instances should be
//-              more specific.
//- Checked By: 
//- Version: $Id: 

#ifndef TOOL_DATA_HPP
#define TOOL_DATA_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"
#include <assert.h>

class CUBIT_UTIL_EXPORT ToolData
{
  private:

    ToolData *nextToolData;
    //- next tool_data in the chain. Last one has this member NULL.

  protected:

    ToolData() {nextToolData = NULL;}
    //- protected constructor to eliminate creation by other than
    //- derived classes

  public:

    virtual ~ToolData();
    //- destructor
 
    ToolData* next_tool_data()               { return nextToolData; }
    void      next_tool_data(ToolData *data) { assert( data != this );
    nextToolData = data; }
    //- access to next tooldata in chain
  
    // handy for ToolDataUser::delete_TD, see e.g. DoubletPillower.cc
  
};
// ********** BEGIN INLINE FUNCTIONS       **********

// ********** END INLINE FUNCTIONS       **********


#endif // TOOL_DATA_HPP








