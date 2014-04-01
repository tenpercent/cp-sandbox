//- Class: TDUniqueId
//- Owner: Tim Tautges
//- Description: This tool data generates unique ids using a random number
//-    generator, and keeps an id in this object
//- A NOTE ABOUT RANDOM NUMBER GENERATION USED IN THIS CLASS:
//- A 3rd party random number generator is used in this class, to overcome limitations
//- of the (16-bit) random number generator on windows systems.  This generator is in
//- the TRandomMersenne class; see that class for documentation of the actual RNG source.
//- In this class (TDUniqueId), the random number generator is seeded DURING THE FIRST CALL
//- TO GET A UNIQUE ID with the output of time(NULL).  This "function-static" method is 
//- used to further randomize the seed to decrease the liklihood of getting duplicate uids.
//- 
//- Checked By: 
//- Version:

#ifndef TD_UNIQUE_ID
#define TD_UNIQUE_ID

#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "CubitDefines.h"
#include "CubitGeomConfigure.h"
#include <map>

template <class X> class DLIList;
class ToolDataUser;
class RefEntity;
class TDUniqueId;

typedef std::multimap<long, TDUniqueId*> TDUIDList;

//This map is used when copying an entity. It maps the unique 
//id of the original to the unique id of the copy
typedef std::map<long, long> COPYUIDMap;

class CUBIT_GEOM_EXPORT TDUniqueId : public ToolData
{
private:

  int uniqueId;
    //- unique id of the owning entity

  ToolDataUser *ownerEntity;
    //- back pointer to the owning entity (needed for sorted lists of
    //- TDUniqueId's)
  
  static MemoryManager memoryManager;
    //- memory management object

  static TDUIDList uniqueIdList;
    //- static list of all entities containing unique ids

  static COPYUIDMap mapForCopying;
    //- maps original unique id to that on copy of entity 

  static int initialize();
    //- initialize the random number generator

  static TDUIDList &unique_id_list();
    //- get a reference to the unique id map

public:

  static int generate_unique_id();
  
  TDUniqueId(ToolDataUser *owner, const int id = 0);
    
  virtual ~TDUniqueId();
    //-constructor and destructor

#ifdef BOYD14
    //- clear out the uid list
  static void clear_uniqueid_list();
#endif
    
    //- clear copy map
  static void clear_copy_map();

  static int is_unique_id(const ToolData* td);
  
  static int get_unique_id(ToolDataUser *owner,
                           const CubitBoolean create_new = CUBIT_TRUE);
    //- get the unique id for owner, or create a new one

  int unique_id(); /*{return uniqueId;};*/
  void unique_id(const int id) {uniqueId = id;};
    //- get/set functions for unique id

  static int get_unique_id_for_copy( int original_id );
    //- when copying an entity that has unique id, gets another
    //unique id for the copy

  ToolDataUser *owner_entity() {return ownerEntity;};
  void owner_entity(ToolDataUser *owner) {ownerEntity = owner;};
    //- get/set functions for ownerEntity
  
  static ToolDataUser *find_td_unique_id(const int temp_id,
                                         const RefEntity *related_entity = NULL);
    //- find the tdu with id temp_id (sorts the list if not sorted)

  static int find_td_unique_id(const int temp_id,
                               DLIList<ToolDataUser*> &td_list,
                               const RefEntity *related_entity = NULL);
    //- find all tdus with id temp_id (sorts the list if not sorted); returns num found

  SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators

  static void set_memory_allocation_increment(int increment = 0)
      {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment

  static void destroy_memory()
      {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

};

inline TDUIDList &TDUniqueId::unique_id_list() 
{
  return uniqueIdList;
}
#endif 


