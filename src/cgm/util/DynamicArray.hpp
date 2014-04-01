//- Class: DynamicArray
//- Description: DynamicArray is an array of generic pointers that is 
//-              accessible by an index and can be grown as items are added.
//-              The type of the pointers is controlled by a macro that
//-              substitutes the appropriate pointer definition. The array
//-              has an initial size (with default) and is grown when indexing
//-              goes beyond the current size. Initial values for all
//-              positions will be zero.
//-
//- Assumptions: The list does not assume that items are stored consecutively
//-              in the array. It is up to the caller to determine if the
//-              data in an array position is valid. This is important if
//-              a conversion is done to a DLList. Data will not be shifted
//-              in the array without explicit action. Removing a pointer
//-              from the list has no effect on the entity pointed to.
//-
//- Owner: Paul Kinney
//- Checked by:
//- Version: $Id: 

#ifndef DYNAMICARRAY_HPP
#define DYNAMICARRAY_HPP

#include "MemoryManager.hpp"
#include "ArrayBasedContainer.hpp"
#include <string.h>
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT DynamicArray : public ArrayBasedContainer
{

public:

    friend class DynamicArrayIterator;

    DynamicArray ( int size );
    //- Standard constructor: Create an array of size {size}. The list will
    //- be grown by {size} each time an index beyond the current end is used.

    DynamicArray (const DynamicArray& copy_from );
    DynamicArray (const ArrayBasedContainer& copy_from );
    //- Create from another DynamicArray

    ~DynamicArray();
    //- Destructor: free all resources

    void set_size(int size);
    //- Change the array size to be {size}. Items after are lost.

    void set_increment(int k);
    //- Change array growth increment

#ifdef BOYD15
    void shift(int from, int to, int count);
    //- Shift count items in either direction.
#endif

  SetDynamicMemoryAllocation(memoryManager)
  //- overloaded new and delete operators

protected:
    void*& operator[](int new_index);
    //- Array dereferencer, on either side of assignment.

    int find_item(void*) const;
    //- Position to item and return that index.

    void append_item(void*);
    //- Insert the item after the last used position.

private:

  static const char* type() { return "DynArr"; }

  static MemoryManager memoryManager;
  //- memory management static objects for performing block allocation
};



//===================================================================
class CUBIT_UTIL_EXPORT DynamicArrayIterator
{
public:
  DynamicArrayIterator() { dynamicArray = NULL; }
    // constructor
  virtual ~DynamicArrayIterator() 
    { dynamicArray = NULL; }
  
    // destructor
  
  int size() const
    { 
      if ( dynamicArray ) 
        return dynamicArray->size();
      return 0; 
    }
  
    // see find and watch below.
  
protected: 

    // see operator [] in derived classes below.
  void* get( int i ) const
    {
      if (dynamicArray && (i < dynamicArray->itemCount) && (i >= 0) )
        return dynamicArray->listArray[i];
      return NULL;
    }

    // data
  const DynamicArray *dynamicArray;
};

//- Define a DynamicArray of a particular type. This becomes a derived class.

#define DynamicArrayDeclare(name,typePtr)                                   \
class name##Iterator;                                                       \
class name : public DynamicArray                                            \
{                                                                           \
  friend class name##Iterator;                                              \
public:                                                                     \
  name ( int size = 0) : DynamicArray (size) {}                             \
  name ( const name& copy_from ) : DynamicArray (copy_from) {}              \
  name ( const ArrayBasedContainer& copy_from ) : DynamicArray (copy_from) {}\
  name& operator=(const name& copy_from)                                    \
               { return (name&)    DynamicArray::operator=(copy_from); }    \
  name& operator=(const name##Iterator& from_iter);                         \
  typePtr& operator[](int new_index)                                        \
               { return (typePtr&) DynamicArray::operator[](new_index); }   \
  int find( typePtr itemPtr ) const { return find_item ( (void*) itemPtr); }\
  void append(typePtr itemPtr) {DynamicArray::append_item((void*)itemPtr);} \
};                                                                          \
                                                                            \
                                                                            \
class name##Iterator : public DynamicArrayIterator                          \
{                                                                           \
  friend class name;                                                        \
public:                                                                     \
  void watch( const name *dynamic_array )                                   \
    { dynamicArray = dynamic_array; }                                       \
  typePtr operator[](int i) const                                     \
     { return (typePtr) get( i ); }                                         \
  int find( typePtr itemPtr ) const                                         \
    { if (dynamicArray)                                                     \
        return ((name*)dynamicArray)->find(itemPtr);                        \
      return -1; }                                                          \
};                                                                          \
inline name& name::operator=(const name##Iterator& from_iter)               \
{ if (from_iter.dynamicArray) *this = *((name*) from_iter.dynamicArray);    \
  else clean_out();                                                         \
  return *this; }                                                           


// Use the dynamicArray iterator like so:
//
// DynamicArrayDeclare(DEdgeUseArray,EdgeUse*);
// DEdgeUseArray edge_use_array;
// DEdgeUseArrayIterator edge_use_iterator;
// edge_use_iterator.watch( edge_use_array );
// for (int i = edge_user_iterator.size(); i--; )
// {
//    CubitNode *node = node_iterator[i];
//    ...
// }
//


#endif //- DYNAMICARRAY_HPP

