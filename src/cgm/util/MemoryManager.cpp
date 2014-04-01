//-       Class: MemoryManager
//-       Owner: Jim Hipp
//- Description: MemoryManager provides object information and stack
//-              storage pointer management used in association with
//-              MemoryBlock and MemoryAllocation classes.
//- Checked By: 
//-    Version:
//  Modified 5/7/96 R W Ostensen:  Fixed operator_new function as per recommendation
//  from J Hipp.  Eliminated assignment of headOfFreeList to wild pointer loaction when
//  memAllocatnSize=1.


#include "MemoryBlock.hpp"
#include "MemoryManager.hpp"
#include "CubitMessage.hpp"
#include "ArrayBasedContainer.hpp"
#include "AppUtil.hpp"
#include <assert.h>
#include <string.h>
#ifndef NT
#include <unistd.h>
#endif

MemoryManager* MemoryManager::memoryManagerListHead = NULL;

// constructors
MemoryManager::MemoryManager() : useNew(false)
{
  // do not allow default constructor (must assign objectSize)
  assert(0);
}

MemoryManager::MemoryManager(const MemoryManager&) : useNew(false)
{
  // do not allow copy
  assert(0);
}

MemoryManager::MemoryManager(const char* name, size_t size, int mem_size,
                             int static_flag)
#ifndef NDEBUG
  : useNew(getenv("CUBIT_DISABLE_MEMORY_MANAGERS") != NULL)
#else
  : useNew(false)
#endif
{
    // allocate space for name and copy into buffer
  if (name)
  {
    objectName = new char[strlen(name) + 1];
    strcpy(objectName, name);
  }
  else
    objectName = NULL;
  
    // initialize manager parameters
  staticManager   = static_flag;
  memBlockStack   = NULL;
  headOfFreeList  = NULL;
  
  assert(size >= sizeof(char*));
  objectSize = size;
  
  set_memory_allocation_increment(mem_size);
  
    // attach new manager to static stack
  next = memoryManagerListHead;
  memoryManagerListHead = this;
}

// destructor
MemoryManager::~MemoryManager()
{
  // delete objectName space and destroy block memory
  if (objectName)
    delete [] objectName;
  destroy_memory(staticManager);
  
  // find 'prev' manager pointer
  MemoryManager* prev        = NULL;
  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager != this)
  {
    prev        = mem_manager;
    mem_manager = mem_manager->next;
  }

  // remove memory manager from static list

  if (!next)
  {
    if (!prev)
    {
      // no memory managers left

      memoryManagerListHead = NULL;
    }
    else
    {
      // remove tail node

      prev->next = NULL;
    }
  }
  else
  {
    if (!prev)
    {
      // remove head node

      memoryManagerListHead = next;
    }
    else
    {
      // remove intermediate node

      prev->next = next;
    }
  }
}

// set memory block allocation increment
void MemoryManager::set_memory_allocation_increment(int mem_size)
{
  // set memory allocation increment

  if (mem_size <= 0)
  {
     memAllocatnSize = DEFAULT_MEMORY_ALLOC_SIZE;
  }
  else
  {
    memAllocatnSize = mem_size;
  }
}

// destroy allocated block memory
void MemoryManager::destroy_memory(int static_flag)
{
  // assert if this is not a static memory manager and some of it's objects
  // are still in use
  if (!static_flag) {
    assert (!get_used_objects());
  }

  // else delete block memory if any was allocated
  if (memBlockStack) delete memBlockStack;

  // reset stack pointers to empty
  headOfFreeList = NULL;
  memBlockStack  = NULL;
}

// return number of objects allocated
int MemoryManager::get_allocated_objects()
{
  // return total number of objects allocated for this memory manager

  if (memBlockStack)
  {
    return (memBlockStack->get_memory_allocation() / objectSize);
  }
  else
  {
    return 0;
  }
}

// return free objects
int MemoryManager::get_free_objects()
{
  // return total number of free objects (allocated but not used) for this
  // memory manager

  if (headOfFreeList)
  {
    int i = 0;
    char* list_ptr = headOfFreeList;
    while (list_ptr)
    {
      i++;
      list_ptr = *((char**) list_ptr);
    }

    return i;
  }
  else
  {
    return 0;
  }
}

// return used objects
int MemoryManager::get_used_objects()
{
  // return total number of used objects for this memory manager

  return (get_allocated_objects() - get_free_objects());
}

// print allocation information to the command line
void MemoryManager::show_object_memory(const char* name)
{
  int instance = 0;

  // find all instances of memory managers for object: name

  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, name))
    {
      // if found then print pertinent memory allocation information

      instance++;
      if (instance > 1)
      {
        PRINT_INFO("\nObject Name(%d): %s\n\n", instance, name);
      }
      else
      {
        PRINT_INFO("\nObject Name: %s\n\n", name);
      }

      int a_obj = mem_manager->get_allocated_objects();
      int f_obj = mem_manager->get_free_objects();
      int u_obj = a_obj - f_obj;

      PRINT_INFO("  Object Size: %d     Allocation Increment: %d\n\n",
		 mem_manager->objectSize, mem_manager->memAllocatnSize);
      PRINT_INFO("  Allocated Objects: %d  (bytes) %d\n", a_obj,
		 a_obj * mem_manager->objectSize);
      if (a_obj)
      {
        PRINT_INFO("       Free Objects: %d  (bytes) %d (%d%%)\n", f_obj,
	  	   f_obj * mem_manager->objectSize,
		   (100*f_obj)/a_obj);
        PRINT_INFO("       Used Objects: %d  (bytes) %d (%d%%)\n", u_obj,
		   u_obj * mem_manager->objectSize,
		   (100*u_obj)/a_obj);
      }
    }

    // get next memory manager

    mem_manager = mem_manager->next;
  }

  // if none were found then announce

  if (!instance)
  {
    PRINT_INFO("\nObject: %s was not found ...\n",  name);
  }
}


// show allocation for all memory manager objects
void MemoryManager::show_all_object_memory()
{
  long int a_obj_byte_total = 0;
  long int f_obj_byte_total = 0;
  long int u_obj_byte_total = 0;

  // loop over all memory managers

  PRINT_INFO("\nDynamic Memory Allocation per Object\n\n");

  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager) {
    long int a_obj = mem_manager->get_allocated_objects();
    long int f_obj = mem_manager->get_free_objects();
    long int u_obj = a_obj - f_obj;
    if (a_obj) {
      // sum total allocated memory parameters (bytes)
      a_obj_byte_total += a_obj * mem_manager->objectSize;
      f_obj_byte_total += f_obj * mem_manager->objectSize;
      u_obj_byte_total += u_obj * mem_manager->objectSize;
    }
    mem_manager = mem_manager->next;
  }

  mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    // print pertinent memory allocation information 

    long int a_obj = mem_manager->get_allocated_objects();
    long int f_obj = mem_manager->get_free_objects();
    long int u_obj = a_obj - f_obj;

    if (a_obj)
    {
      PRINT_INFO("\nObject Name: %s\n\n", mem_manager->objectName);

      PRINT_INFO("  Object Size: %d     Allocation Increment: %d\n\n",
	         mem_manager->objectSize, mem_manager->memAllocatnSize);
      if (a_obj_byte_total != 0)
      {
         PRINT_INFO("  Allocated Objects: %ld  (bytes) %ld (%d%% of Total)\n",
                    a_obj, a_obj * mem_manager->objectSize,
                    int((100.0*a_obj * mem_manager->objectSize)/a_obj_byte_total));
         PRINT_INFO("       Free Objects: %ld  (bytes) %ld (%d%%)\n", f_obj,
                    f_obj * mem_manager->objectSize,
                    int((100.0*f_obj)/a_obj));
         PRINT_INFO("       Used Objects: %ld  (bytes) %ld (%d%%)\n", u_obj,
                    u_obj * mem_manager->objectSize,
                    int((100.0*u_obj)/a_obj));
      }
      else
      {
         PRINT_INFO("  Allocated Objects: %ld  (bytes) %ld (100%% of Total)\n",
                    a_obj, a_obj * mem_manager->objectSize);
         PRINT_INFO("       Free Objects: %ld  (bytes) %ld (100%%)\n", f_obj,
                    f_obj * mem_manager->objectSize);
         PRINT_INFO("       Used Objects: %ld  (bytes) %ld (100%%)\n", u_obj,
                    u_obj * mem_manager->objectSize);
      }
    }

    mem_manager = mem_manager->next;
  }

  // print total memory allocation information

  char sizechar;

  PRINT_INFO("\nTotal Memory Allocation Information\n\n");
  int divisor;
  if (a_obj_byte_total > 10000000) {
    sizechar = 'M';
    divisor = 1000000;
  }
  else {
    sizechar = 'K';
    divisor = 1000;
  }
  
  PRINT_INFO("  Allocated Memory: %ld%c (%ld bytes)\n", a_obj_byte_total/divisor, 
             sizechar, a_obj_byte_total);

  if (a_obj_byte_total)
  {
    PRINT_INFO("       Free Memory: %ld%c (%d%%)\n", f_obj_byte_total/divisor, sizechar,
	       int((100.0*f_obj_byte_total)/a_obj_byte_total));
    PRINT_INFO("       Used Memory: %ld%c (%d%%)\n", u_obj_byte_total/divisor, 
               sizechar,
	       int((100.0*u_obj_byte_total)/a_obj_byte_total));
  }

#ifndef JANUS
#ifndef NT
  struct rusage r_usage;
  AppUtil::instance()->apputil_getrusage(r_usage);
  PRINT_INFO("       (System reports %ld%c used, incl. executable)\n", 
              r_usage.ru_maxrss*getpagesize()/divisor, sizechar);
#else
  PRINT_INFO("\n");
#endif // NT
#endif // JANUS

  // print DLList non-Pool allocation information
  PRINT_INFO("\nTotal non-pool ArrayBasedContainer memory allocation  = %u%c\n"
	       "Maximum non-pool ArrayBasedContainer memory allocated = %u%c\n",
             ArrayBasedContainer::current_allocated_memory()/divisor, sizechar, 
             ArrayBasedContainer::maximum_allocated_memory()/divisor, sizechar);
#if 0
  // print HOOPS memory usage
  long allocated = 0;
  long in_use = 0;
    //  DrawingTool::instance()->show_memory(allocated, in_use);
  if (allocated != 0)
     PRINT_INFO("\nGraphics subsystem memory: Allocated = %u%c (%d bytes)\n"
                "                           In-Use    = %u%c (%d%%)\n",
                allocated/divisor, sizechar, allocated, 
                in_use/divisor, sizechar, int((100.0*in_use)/allocated));
  else
     PRINT_INFO("\nGraphics subsystem memory: Allocated = %u%c (%d bytes)\n"
                "                           In-Use    = %u%c (100%%)\n",
                allocated/divisor, sizechar, allocated, 
                in_use/divisor, sizechar);
#endif
}

// compress memory for the requested object
int MemoryManager::compress_object_memory(const char* name)
{
  // find all instances of memory managers for object: name

  int found        = 0;
  int saved_memory = 0;
  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, name))
    {
      // if found then compress memory

      saved_memory += mem_manager->compress_memory();
      found         = 1;
    }

    mem_manager = mem_manager->next;
  }

  return found;
}

// compress all object memory
int MemoryManager::compress_all_object_memory()
{
  // find all instances of memory managers
  int saved_memory = 0;
  MemoryManager* mem_manager = memoryManagerListHead;
  MemoryManager* block_manager = NULL;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, "MemoryBlock"))
    {
      // save block_manager until end

      block_manager = mem_manager;
    }
    else
    {
      // compress memory

      saved_memory += mem_manager->compress_memory();
    }

    mem_manager = mem_manager->next;
  }

  if (block_manager)
  {
    saved_memory += block_manager->compress_memory();
  }

  return saved_memory;
}
      
// generic operator new call
void* MemoryManager::operator_new(size_t size)
{
#ifndef NDEBUG
  if (useNew) return malloc(size);
#endif

  // send requests of "wrong" size to ::new

  if (size != objectSize) return ::new char[size];

  // get new element from head of free list

  char* p = headOfFreeList;

  if(!p)
    {
    // allocate new block

    int block_size  = memAllocatnSize * size;
    char* new_block = ::new char[block_size];
    if (!new_block) return (void*) NULL;

    // link new elements to form the free list

    int fill_limit = (memAllocatnSize - 1) * size;
    for (int j = 0; j < fill_limit; j += size)
    {
      *((char**) &new_block[j]) = &new_block[j + size];
    }
    *((char**) &new_block[fill_limit]) = (char*) NULL;

    // assign new element

    p = new_block;

    // save new block to memory block stack

    memBlockStack = new MemoryBlock(memBlockStack, new_block, block_size);
  }
  //assign head of free list and return p

  headOfFreeList = *((char**) p);
  return (void*) p;
}

// generic operator delete call
void MemoryManager::operator_delete(void *deadObject, size_t size)
{
#ifndef NDEBUG
  if (useNew) 
  {
    free(deadObject);
    return;
  }
#endif

  // requests of "wrong" size to ::delete

  if (size != objectSize)
  {
    ::delete [] ((char*) deadObject);
    return;
  }

  // attach dead element to head of free list

  char* delete_object = (char*) deadObject;
  *((char**) delete_object) = headOfFreeList;
  headOfFreeList = delete_object;
}

// compress memory blocks
int MemoryManager::compress_memory()
{
  // if free objects exist then begin compression algorithm

  if (headOfFreeList)
  {
    // find total number of memory blocks attached to stack
 
    int n_blocks = 0;
    MemoryBlock* mem_block = memBlockStack;
    while (mem_block)
    {
      n_blocks++;
      mem_block = mem_block->next_block();
    }

    if (n_blocks == 0)
    {
       // no available memory to free ... return 0
       // this is here for safety ... n_blocks should never be zero if
       // headOfFreeList is not Null

       return 0;
    }
    else
    {
      // first determine if all objects are free ... if so then perform
      // the easy compression routine

      if (!get_used_objects())
      {
        // all objects are free ... delete all blocks

        int n_bytes = memBlockStack->get_memory_allocation();
        destroy_memory(staticManager);

        // return freed memory

        return n_bytes;
      }

      // else perform the complex routine to remove those memory blocks that
      // have all free elements

      // begin by constructing an integer array for each memory block to
      // tally the number of free objects that each block contains

        // if there are a lot of blocks, we can save a huge amount of
        // time by looking in the last few blocks that contained an
        // element.
      const int use_cache = n_blocks > 8;
      int i, j, k;
      i = j = k = 0;
      const int cache_size = 4;
      MemoryBlock* mem_block_sav[cache_size];
      int i_sav[cache_size];
      for ( i = cache_size; i--; )
      {
        mem_block_sav[i] = NULL;
        i_sav[i] = 0;
      }
      int found = 0;
      
      mem_block = NULL;
      char* list_ptr = NULL;
      
      unsigned int* free_tally = new unsigned int [n_blocks];
      for (i = 0; i < n_blocks; i++) free_tally[i] = 0;

      // loop through free list tallying free elements

      list_ptr = headOfFreeList;
      while (list_ptr)
      {
        // find memory block that owns this element
        
          // check last few blocks for speed
        found = CUBIT_FALSE;
        if ( use_cache )
        {
          for ( i = 0; i < cache_size; i++ )
          {
            mem_block = mem_block_sav[i];
            if ( mem_block &&
                 list_ptr >= mem_block->get_block() &&
                 list_ptr < (mem_block->get_block() + 
                             mem_block->block_size()) )
            {
              k = i_sav[i];
              free_tally[k]++;
              found = CUBIT_TRUE;
              break;
            }
          }
        }
        if ( !found )
        {
            // search through all blocks
          mem_block = memBlockStack;
          for (i = 0; i < n_blocks; i++)
          {
            if ((list_ptr >= mem_block->get_block()) &&
                (list_ptr < (mem_block->get_block() + mem_block->block_size())))
            {
                // increment tally and exit
              
              free_tally[i]++;
                // save
              if ( use_cache && mem_block_sav[j] != mem_block )
              {
                mem_block_sav[j] = mem_block;
                i_sav[j] = i;
                if ( ++j >= cache_size )
                  j = 0;
              }
              break;
            }
            
              //  get next memory block
            mem_block = mem_block->next_block();
          }
        }

        // get next element
        list_ptr = *((char**) list_ptr);
      }

      // zero tally for memory blocks that cannot be removed ... those that
      // have some used elements

      int all_blocks = 0;
      mem_block = memBlockStack;
      for (i = 0; i < n_blocks; i++)
      {
        if (free_tally[i] != (mem_block->block_size() / objectSize))
	{
          free_tally[i] = 0;
          all_blocks++;
        }

        mem_block = mem_block->next_block();
      }

      if (all_blocks == n_blocks)
      {
        // no memory can be saved ... all blocks have some used elements
        // return 0

        delete [] free_tally;
        return 0;
      }

      // adjust free list pointers to remove those that belong to
      // memory blocks that can be deleted
      char* prev_ptr = NULL;
      list_ptr = headOfFreeList;
      while (list_ptr)
      {
        // find memory block that owns this element
          // check last few blocks for speed
        found = CUBIT_FALSE;
        if ( use_cache )
        {
          for ( i = 0; i < cache_size; i++ )
          {
            mem_block = mem_block_sav[i];
            if ( mem_block &&
                 list_ptr >= mem_block->get_block() &&
                 list_ptr < (mem_block->get_block() + 
                             mem_block->block_size()) )
            {
              k = i_sav[i];
              found = CUBIT_TRUE;
              break;
            }
          }
        }
        if ( !found )
        {
          mem_block = memBlockStack;
          for (i = 0; i < n_blocks; i++)
          {
            if ((list_ptr >= mem_block->get_block()) &&
                (list_ptr < (mem_block->get_block() + mem_block->block_size())))
            {
              k = i;
                // save
              if ( use_cache && mem_block_sav[j] != mem_block )
              {
                mem_block_sav[j] = mem_block;
                i_sav[j] = i;
                if ( ++j >= cache_size )
                  j = 0;
              }
              break;
            }
              // get next memory block
            mem_block = mem_block->next_block();
          }
        }
        
        if (free_tally[k])
        {
            // remove element
          
          if (prev_ptr)
          {
            *((char**) prev_ptr) = *((char**) list_ptr);
          }
          else
          {
            headOfFreeList = *((char**) list_ptr);
          }
        }
        else
        {
            // advance prev_ptr  
          prev_ptr = list_ptr;
        }

        // get next element
        list_ptr = *((char**) list_ptr);
      }

      // delete all memory blocks that have free_tally[i] > 0

      i = 0;
      int save_bytes = 0;
      MemoryBlock* prev_block = NULL;
      mem_block               = memBlockStack;
      while (mem_block)
      {
        if (free_tally[i])
        {
          // set previous MemoryBlocks next pointer to skip this block

          if (prev_block)
          {
            prev_block->next_block(mem_block->next_block());
          }
          else
          {
            memBlockStack = mem_block->next_block();
          }

          // set MemoryBlock next pointer to NULL to avoid recusive delete
          // update saved memory and delete mem_block

          mem_block->next_block((MemoryBlock*) NULL);
          save_bytes += mem_block->block_size();
          delete mem_block;

          // update mem_block to point to new current MemoryBlock

          if (prev_block)
          {
            mem_block = prev_block->next_block();
          }
          else
          {
            mem_block = memBlockStack;
          }
        }
        else
        {
          // if block wasn't removed then update previous and current blocks

          prev_block = mem_block;
          mem_block = mem_block->next_block();
        }

        // increment to next block (used by free_tally array)

        ++i;
      }

      // return freed memory (bytes)

      delete [] free_tally;
      return save_bytes;
    }
  }
  else
  {
    // no memory allocated ... return 0

    return 0;
  }
}
