//- Class: TDExport
//- Owner: Joel Kopp
//- Description: 
//- Checked By: 
//- Version:

#ifndef TD_EXPORT_HPP
#define TD_EXPORT_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "ToolData.hpp"

class TDExport : public ToolData
{
  private:

    static MemoryManager memoryManager;
    //- memory management object

    ElementType elementType;
    //- List of boundary cards containing this node
  
  public:

    TDExport();
    //- constructor

    static int is_export(const ToolData* td) 
      {return (CAST_TO(const_cast<ToolData*>(td), TDExport) != NULL);}
    
    ElementType get_element_type();
    void set_element_type( ElementType type );  

    SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators

    static void set_memory_allocation_increment(int increment = 0)
                {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment

    static void destroy_memory()
                {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object
};

inline TDExport::TDExport() : ToolData() { }

inline ElementType TDExport::get_element_type()
                      { return elementType; }

inline void TDExport::set_element_type ( ElementType type ) {elementType = type;}

#endif // TD_BOUNDARY_KNOWING_HPP
