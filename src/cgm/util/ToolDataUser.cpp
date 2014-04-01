// Class: ToolDataUser

#include <assert.h>
#include "ToolDataUser.hpp"
#include "ToolData.hpp"
#include "DLIList.hpp"

ToolData *ToolDataUser::remove_TD( IdentityFn specified_type )
{
  ToolData *td = tool_data();
  ToolData *td_prev = NULL;
  while (td) {
    if ( (*specified_type)(td) ) {
      if (td_prev)
        td_prev->next_tool_data( td->next_tool_data() );
      else
        toolData = td->next_tool_data();
      td->next_tool_data( NULL );
      return td;
    }
    td_prev = td;
    td = td->next_tool_data();
  }
  return NULL;
}

ToolData *ToolDataUser::remove_TD( ToolData *td_remove )
{
  ToolData *td = tool_data();
  ToolData *td_prev = NULL;
  while (td) {
    if ( td_remove == td ) {
      if (td_prev)
        td_prev->next_tool_data( td->next_tool_data() );
      else
        toolData = td->next_tool_data();
      td->next_tool_data( NULL );
      return td;
    }
    td_prev = td;
    td = td->next_tool_data();
  }
  return NULL;
}

CubitBoolean ToolDataUser::delete_TD( IdentityFn specified_type )
{
  ToolData *td = remove_TD( specified_type );
  if (td)
  {
    delete td;
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitBoolean ToolDataUser::delete_TD( ToolData *td_remove )
{
  ToolData *td = remove_TD( td_remove );
  if (td)
  {
    delete td;
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

void ToolDataUser::add_TD( ToolData *new_td )
{
  assert( new_td != NULL );
  new_td->next_tool_data( toolData );
  toolData = new_td;
}

ToolData *ToolDataUser::get_TD( IdentityFn specified_type )
{
  ToolData *td = tool_data();
  while (td) {
    if ( (*specified_type)(td) ) 
      return td;
    td = td->next_tool_data();
  }
  return NULL;
}

ToolData const* ToolDataUser::get_TD( IdentityFn specified_type ) const
{
  ToolData *td = tool_data();
  while (td) {
    if ( (*specified_type)(td) ) 
      return td;
    td = td->next_tool_data();
  }
  return NULL;
}

void ToolDataUser::get_all_TDs(IdentityFn specified_type, 
                               DLIList <ToolData *> *all_tds) const
{
  ToolData *td = tool_data();
  while (td) 
  {
    if ( (*specified_type)(td) ) 
      all_tds->append(td);
    td = td->next_tool_data();
  }  
}


ToolDataUser::~ToolDataUser()
{
  //delete all ToolData's chained off this user.
  ToolData *td = tool_data();
  while ( td ) {
    ToolData *next = td->next_tool_data();
    delete td;
    td = next;
  }

    // set the first TD to NULL
  tool_data(NULL);
}



