//- Class: TDCellIndex
//- Owner: Ray Ostensen
//- Description: This tool data contains the cellIndex used for sorting
//- nodes in GridSearch
//- Checked By: 
//- Version:

#ifndef TD_CELLINDEX_HPP
#define TD_CELLINDEX_HPP

#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "CastTo.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT TDCellIndex : public virtual ToolData
{
private:

  static MemoryManager memoryManager;
    //- memory management object

  int cellIndex;

  
public:

    TDCellIndex(int cell_index = -1 );
    virtual ~TDCellIndex() {}
    //-constructor and destructor

    static int is_cell_index(const ToolData* td)
     {return (CAST_TO(const_cast<ToolData*>(td), TDCellIndex) != NULL);}
  
    void cell_index ( int cell_index ) { cellIndex = cell_index; }
    int cell_index () { return cellIndex; }
    //- Heading: get/set member data

    SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators

    static void set_memory_allocation_increment(int increment = 0)
                {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment

    static void destroy_memory()
                {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

};

inline TDCellIndex::TDCellIndex(int cell_index) : ToolData(), cellIndex(cell_index)
{}



#endif // TD_CellIndex_HPP




