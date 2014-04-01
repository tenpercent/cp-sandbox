//-       Class: MemoryBlock
//-       Owner: Jim Hipp
//- Description: MemoryBlock provides stack storage for block or "Chunk"
//-              memory allocated with overloaded 'new' operators in classes
//-              utilizing a MemoryManager object.
//- Checked By: 
//-    Version:

#ifndef MEMORY_BLOCK_HPP
#define MEMORY_BLOCK_HPP

#include <stdlib.h>
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT MemoryBlock
{
  private:

    char*        block;
    int          size;
    MemoryBlock* next;
    //- block is the pointer to the beginning of the memory block of objects
    //- of type char, size is the size of the memory block, and next is a
    //- pointer to the next MemoryBlock object.

  public:

    MemoryBlock (MemoryBlock* Head, char* Block, int Size)
    : block(Block), size(Size), next(Head) {}
    //- Constructor: performs initialization assignment

   ~MemoryBlock();
    //- Destructor: performs recursive delete of MemoryBlock stack

    int get_memory_allocation();
    //- returns total amount of allocated memory (bytes)

    MemoryBlock* next_block() const;
    void         next_block(MemoryBlock* next_block);
    //- retreive / set next pointer

    int block_size() const;
    //- get block_size

    char* get_block() const;
    //- get block
};

inline MemoryBlock* MemoryBlock::next_block() const
{
  return next;
}

inline void         MemoryBlock::next_block(MemoryBlock* nxt_blck)
{
  next = nxt_blck;
}

inline int          MemoryBlock::block_size() const
{
  return size;
}

inline char*        MemoryBlock::get_block() const
{
  return block;
}

#endif // MEMORY_BLOCK_HPP

