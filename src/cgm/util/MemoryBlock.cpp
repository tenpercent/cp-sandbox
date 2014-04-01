//-       Class: MemoryBlock
//-       Owner: Jim Hipp
//- Description: MemoryBlock provides stack storage for block or "Chunk"
//-              memory allocated with overloaded 'new' operators in classes
//-              utilizing a MemoryManager object.
//- Checked By: 
//-    Version:

#include"MemoryBlock.hpp"


// destructor: recusively deletes stack of allocated memory
MemoryBlock::~MemoryBlock()
{
  delete [] block;
  if (next) 
  {
    delete next;
    next = NULL;
  }
}

// retreive total amount of memory allocated (bytes) including all blocks
// beneath 'this' block on the stack
int MemoryBlock::get_memory_allocation()
{
  MemoryBlock* block_ptr = this;

  // sum all allocated memory

  int amount = 0;
  while (block_ptr)
  {
    amount += block_ptr->size;
    block_ptr = block_ptr->next;
  }

  return amount;
}



