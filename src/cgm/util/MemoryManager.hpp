//-       Class: MemoryManager
//-       Owner: Jim Hipp
//- Description: MemoryManager provides object information and stack
//-              storage pointer management used in association with
//-              MemoryBlock and MemoryAllocation classes.
//- Checked By: 
//-    Version:
 
#ifndef MEMORY_MANAGER_HPP
#define MEMORY_MANAGER_HPP

#define USE_DYNAMIC_MEMORY_ALLOCATION
//- comment out this #define to use global system level dynamic memory
//- allocation ... Otherwise overloaded operators new and delete will
//- be set to call operator_new and operator_delete defined in this module

const  int DEFAULT_MEMORY_ALLOC_SIZE = 1024;
//- Default MemoryBlock allocation size

const  int STATIC_MEMORY_MANAGER = 1;
//- Static Memory Manager definition argument

#include <stdlib.h>
#include <assert.h>
#include "CubitUtilConfigure.h"

class MemoryBlock;

class CUBIT_UTIL_EXPORT MemoryManager
{
private:
  
  const bool     useNew;
    //- Disable memory management (pass all calls to new/delete)
    //- Unused if compiled with -DNDEBUG
 
  char*          objectName;
    //- name of the class for which this memory manager is declared
  
  MemoryBlock*   memBlockStack;
    //- memBlockStack is the head of a stack of MemoryBlock objects
    //- that contain allocated memory chunks
  
  char*          headOfFreeList;
    //- headOfFreeList is the head of a stack that contains free
    //- elements that have not yet been constructed by the allocating
    //- object
  
  size_t         objectSize;
    //- objectSize is the size of the allocating object
  
  int            memAllocatnSize;
    //- memAllocatnSize is the size of the MemoryBlocks that will be
    //- allocated by the allocating object
  
  MemoryManager* next;
    //- next object pointer is used to save memory managers on a
    //- static stack.
  
  int            staticManager;
    //- static flag to indicate that this memory manager is a static
    //- object.  This parameter is used by the destructor when the
    //- program terminates and static object destructors are called.
    //- The destructor calls destroy_memory() which asserts if
    //- the used objects attached to this memory manager have not been
    //- deleted (generally true for static memory managers at program
    //- termination).  Forcing an assert is desirable during normal
    //- program operation but not during a normal program termination
    //- sequence (Typing QUIT at the CUBIT command line prompt and
    //- observing an assert in MemoryManager is not comforting).
  
  static MemoryManager* memoryManagerListHead;
    //- stack head containing all allocated memory managers.
  
  MemoryManager();
  MemoryManager(const MemoryManager&);
    //- do not allow default constructor (must assign objectSize)
  
public:
  
  MemoryManager(const char* name, size_t size, int mem_size = 0,
                int static_flag = 0);
  ~MemoryManager();
    //- constructor / destructor
  
  void   set_memory_allocation_increment(int mem_size = 0);
  int    get_memory_allocation_increment() const;
    //- set / get memory block allocation increment
  
  size_t get_object_size() const;
    //- get object size
  
  void   destroy_memory(int static_flag = 0);
    //- destroy allocated memory
  
  int    compress_memory();
    //- compresses memory for the object by removing unused MemoryBlocks
  
  int    get_allocated_objects();
    //- returns number of object storage locations that have been
    //- allocated
  
  int    get_free_objects();
    //- returns number of objects that are not in use but have been
    //- allocated
  
  int    get_used_objects();
    //- returns number of objects that are currently in use from those
    //- that have been allocated
  
  void*  operator_new(size_t size);
    //- generic operator new
  
  void   operator_delete(void *deadObject, size_t size);
    //- generic operator delete

#ifdef BOYD15
  static int total_allocated_memory();
  //- returns total memory allocated for all memory managers
#endif

  static void show_object_memory(const char* s);
  static void show_all_object_memory();
    //- prints allocation information to the command line.  The
    //- functions are declared static to enable their use throughout
    //- the code.  The character string argument is interpreted as the
    //- class name provided as input when the memory manager was
    //- constructed.
  
  static int compress_object_memory(const char* s);
  static int compress_all_object_memory();
    //- recaptures unused memory blocks and returns them to free store
    //- and returns the amount of reclaimed memory.  The functions are
    //- declared static to enable their use throughout the code.  The
    //- character string argument is interpreted as the class name
    //- provided as input when the memory manager was constructed.
};

#ifdef USE_DYNAMIC_MEMORY_ALLOCATION

#define SetDynamicMemoryAllocation(memManager)                               \
                                                                             \
  void*  operator new(size_t size)                                           \
         {return memManager.operator_new(size);}                             \
  void   operator delete(void *deadObject, size_t size)                      \
         {memManager.operator_delete(deadObject, size);}                     \
  /*  overloaded new and delete operators */                                    \
                                                                             \

#else

#define SetDynamicMemoryAllocation(memManager)                               \
                                                                             \
  /* DO NOT overload new and delete operators ... use global allocation */   \
                                                                             \

#endif

#endif // MEMORY_MANAGER_HPP  

