//- Class: ToolDataUser
//- Owner: Scott Mitchell
//- Description: This is the tool data interface that every class that uses a
//-              tool data needs to specify. E.g. CubitNode and CubitFace
//-              are derived from this class.
//- Checked By: 
//- Version: $Id: 

#ifndef TOOL_DATA_USER_HPP
#define TOOL_DATA_USER_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

#ifdef NT
class __single_inheritance ToolData;
#else
class ToolData;
#endif


typedef int (*IdentityFn)(const ToolData* td);
template <class X> class DLIList;

class CUBIT_UTIL_EXPORT ToolDataUser
{

private:
    ToolDataUser( const ToolDataUser& );
    void operator=( const ToolDataUser&);

  ToolData *toolData;
  //- generic pointer to extra data needed by a particular VolSmoothTool.
  //- see the _TD functions below for how to access this

  void      tool_data(ToolData *set_data) {toolData = set_data;}
  ToolData* tool_data()    const          {return toolData;}
  //- all external access/setting should be done with the _TD functions.
  
public:

  ToolDataUser () { toolData = NULL; }

  virtual ~ToolDataUser();
  //- automatically deletes all the chained ToolDatas
  
  CubitBoolean delete_TD(IdentityFn specified_type);
  //- delete the specific type of tool data from the chain.
  //- return true if something was actually found and deleted.

  CubitBoolean delete_TD(ToolData *td);
  //- delete the specific tool data from the chain.
  //- return true if something was actually found and deleted.
  
  ToolData *remove_TD(IdentityFn specified_type);
  //- remove the specific type of tool data from the chain, and
  //- return it.

  ToolData *remove_TD(ToolData *td);
  //- remove the specific tool data from the chain, and
  //- return it.
  
  void add_TD(ToolData *new_td);
  //- add the new_td to the beginning of the tool_data chain
  
  virtual ToolData *get_TD(IdentityFn specified_type);
  virtual ToolData const *get_TD(IdentityFn specified_type) const;
  virtual void get_all_TDs(IdentityFn specified_type, 
                           DLIList <ToolData *> *all_tds) const;
  //- get the specific type of ToolData in the chain.
  //- returns null if not found.
  
};


#endif // TOOL_DATA_USER_HPP





