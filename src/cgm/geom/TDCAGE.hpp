//- Class:          TDCAGE
//- Descpription:   Stores data used in the CAGenesisEntity attribute
//- Owner:          Tim Tautges
//- Checked by:
//- Version: 

#ifndef TDCAGE_HPP
#define TDCAGE_HPP

#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

class RefGroup;
class RefEntity;

class CUBIT_GEOM_EXPORT TDCAGE : public ToolData 
{
private:

  static MemoryManager memoryManager;
  //- memory management object

  int uniqueId;

  DLIList<RefGroup *> groupList;
    //- keeps group pointers for group/sequence # pairs

  DLIList<int> sequenceList;
    //- keeps sequence #'s for group/sequence # pairs

public:

  TDCAGE(int unique_id);

  ~TDCAGE();

  static int is_cage(const ToolData* td);
  
  int unique_id() {return uniqueId;}
  void unique_id(int id) {uniqueId = id;}
    //- get/set unique id

  SetDynamicMemoryAllocation(memoryManager)
//- class specific new and delete operators
    
  static void set_memory_allocation_increment( int increment = 0 )
  {memoryManager.set_memory_allocation_increment(increment);}
  //- set block memory size increment

  static void destroy_memory()
  {memoryManager.destroy_memory();}
  //- destroy all memory allocated to this object

  static int group_sequence_number(RefGroup *group, const RefEntity *entity);
    //- return the sequence number of the given entity in the given group

  int td_sequence_number(const RefGroup *group);
    //- given a group, return the sequence number paired with that group

  void insert_group(RefGroup *group, const int seq_num);
    //- add the group/sequence number pair to the lists

  static void insert_entity(RefEntity *entity, const int seq_num,
                            RefGroup *into_group);
    //- insert the entity into the group following the sequence number
    //- given

  void initialize_group_sequence_list(RefEntity *entity);
    //- get the groups owning entity and build the sequence lists
};

inline TDCAGE::TDCAGE(int unique_id) : ToolData()
{
  uniqueId = unique_id;
}

inline TDCAGE::~TDCAGE()
{
  //- empty
}

inline void TDCAGE::insert_group(RefGroup *group, const int seq_num) 
{
  groupList.append(group);
  sequenceList.append(seq_num);
}

#endif // TD_CAGE_HPP
    
