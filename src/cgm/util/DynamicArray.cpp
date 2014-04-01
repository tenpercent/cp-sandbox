//- Class: DynamicArray
//- Owner: Paul Kinney
//- Checked by:
//- Version: $Id: 

#include "DynamicArray.hpp"
#include <string.h>

#if defined(NO_MEMMOVE)
extern "C" void *CUBIT_memmove(void *s1, const void *s2, size_t n);
#define memmove CUBIT_memmove
#endif

// Method: Standard constructor
// ==========================================================================
DynamicArray::DynamicArray ( int size ) : ArrayBasedContainer ( size )
{}


// Method: Copy constructor
// ==========================================================================
DynamicArray::DynamicArray (const DynamicArray& copy_from ) : 
  ArrayBasedContainer ( copy_from )
{}

DynamicArray::DynamicArray (const ArrayBasedContainer& copy_from ) :
  ArrayBasedContainer ( copy_from )
{}

// Method: Destructor
// ==========================================================================
DynamicArray::~DynamicArray()
{}


// Method: set_size
// Set the accessed space to the specified size. There is no effect to the
// allocated space.
// ==========================================================================

void DynamicArray::set_size(int size)
{
    itemCount = size;
}



// Method: operator[]
// Return the nth position of the array
// ==========================================================================
void*& DynamicArray::operator[](int index)
{
    if(itemCount < index+1) itemCount = index+1;

    if(itemCount > listLength) lengthen_list();

    return listArray[index];
}


// Method: find_item
// Return the index for a particular item in the array
// ==========================================================================
int DynamicArray::find_item(void* item) const
{
    if(!listArray) return -1;

    for(int index = 0; index < itemCount; ++index)
    {
        if(listArray[index] == item)
        {
            return index;
        }
    }
    return -1;
}


// Method: append
// Add the item to the end of the used portion of the array
// ==========================================================================
void DynamicArray::append_item(void* item)
{
//+// Add item to list without knowing how long it is.
    if ( itemCount == listLength )
        lengthen_list();

    listArray[itemCount++] = item;
}



