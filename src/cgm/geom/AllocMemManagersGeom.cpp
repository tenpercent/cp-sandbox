#include "TDCAGE.hpp"
#include "TDUniqueId.hpp"


const int TDCAGE_ALLOC_SIZE         = 1024;
const int TDUNIQUEID_ALLOC_SIZE = 128;

MemoryManager TDUniqueId::memoryManager("TDUniqueId",
                                       sizeof(TDUniqueId),
                                       TDUNIQUEID_ALLOC_SIZE, 
                                       STATIC_MEMORY_MANAGER);


MemoryManager TDCAGE::memoryManager("TDCAGE",
                                    sizeof(TDCAGE),
                                           TDCAGE_ALLOC_SIZE,
                                           STATIC_MEMORY_MANAGER);


