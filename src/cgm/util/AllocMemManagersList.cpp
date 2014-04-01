#include "DLList.hpp"
#include "DynamicArray.hpp"
#include "Queue.hpp"
#include "SDLList.hpp"
#include "Tree.hpp"

const int ARRAYMEMORY_ALLOC_SIZE    = 8192;
const int DLLIST_ALLOC_SIZE         = 1024;
const int SDLLIST_ALLOC_SIZE        = 128;
const int DYNAMIC_ARRAY_ALLOC_SIZE  = 4096;
const int QUEUE_NODE_ALLOC_SIZE     = 512;
const int TREE_ALLOC_SIZE           = 128;

MemoryManager ArrayMemory::memoryManager("ArrayMemory", sizeof(ArrayMemory),
                                         ARRAYMEMORY_ALLOC_SIZE,
                                         STATIC_MEMORY_MANAGER);
MemoryManager DLList::memoryManager("DLList", sizeof(DLList),
                                    DLLIST_ALLOC_SIZE,
                                    STATIC_MEMORY_MANAGER);
MemoryManager SDLList::memoryManager("SDLList", sizeof(SDLList),
                                     SDLLIST_ALLOC_SIZE,
                                     STATIC_MEMORY_MANAGER);
MemoryManager DynamicArray::memoryManager("DynamicArray", sizeof(DynamicArray),
                                          DYNAMIC_ARRAY_ALLOC_SIZE,
                                          STATIC_MEMORY_MANAGER);
MemoryManager QueueNode::memoryManager("QueueNode", sizeof(QueueNode),
                                       QUEUE_NODE_ALLOC_SIZE,
                                       STATIC_MEMORY_MANAGER);
MemoryManager Tree::memoryManager("Tree", sizeof(Tree),
                                  TREE_ALLOC_SIZE,
                                  STATIC_MEMORY_MANAGER);
				  
#include "TDCellIndex.hpp"
const int CELL_INDEX_ALLOC_SIZE     = 1024;
MemoryManager TDCellIndex::memoryManager("TDCellIndex", sizeof(TDCellIndex),
                                         CELL_INDEX_ALLOC_SIZE,
                                         STATIC_MEMORY_MANAGER);
