//- Class: ArrayBasedContainer
//- Owner: Paul Kinney
//- Checked by:
//- Version: $Id: 

#include "ArrayBasedContainer.hpp"
#include "CubitDefines.h"
#include <string.h>
#include <assert.h>

unsigned int ArrayBasedContainer::currentAllocatedMemory = 0;
unsigned int ArrayBasedContainer::maxAllocatedMemory = 0;

//- Constructor: Create a list with initial size {increment}. The list
//- will be grown by {increment} each time it is filled. Memory for the
//- list is not allocated until the first element is inserted using
//- {insertLink}. 
//- If {increment} is zero, the default increment ({INCREMENT}) will be used
//- From an efficiency standpoint, it is very important that the 
//- increment be set large enough to reduce the number of list 
//- growths, but small enough to not waste memory.
//- It is more efficient to sligthly overestimate the increment than 
//- to underestimate the increment.
ArrayBasedContainer::ArrayBasedContainer ( int increment )
{
  if ( increment > 0)
    {
      listIncrement = increment;
    }
  else
    {
      listIncrement = INCREMENT;
    }
  itemCount  = 0;
  listLength = 0;
  listArray  = 0;
  poolAllocated = CUBIT_FALSE;
  listIsSorted  = CUBIT_FALSE;
}

//- Copy Constructor
//- Uses memory from ArrayMemory if size is default
ArrayBasedContainer::ArrayBasedContainer(const ArrayBasedContainer& from)
{
  listLength = from.listLength;
  listIncrement = from.listIncrement;
  poolAllocated = from.poolAllocated;
  listIsSorted  = from.listIsSorted;
  listArray = NULL;

  if (listLength)
    {
      if (listLength == INCREMENT && (int)listIncrement == INCREMENT) {
	listArray = (void **)new ArrayMemory;
	poolAllocated = CUBIT_TRUE;
      }
      else {
	listArray = new void* [listLength];
	poolAllocated = CUBIT_FALSE;
	currentAllocatedMemory += listLength;
	if (currentAllocatedMemory > maxAllocatedMemory)
	  maxAllocatedMemory = currentAllocatedMemory;
      }
      assert(listArray != 0);
    }

  itemCount     = from.itemCount;

  if (itemCount)
    memcpy (listArray, from.listArray, itemCount*sizeof(void*));
}

ArrayBasedContainer::~ArrayBasedContainer()
{
    if ( listLength )
    { 
      if (poolAllocated) {
	delete ((ArrayMemory*) listArray);
      }
      else {
        delete [] listArray;
	currentAllocatedMemory -= listLength;
      }
      listArray = NULL;
    }
}


void ArrayBasedContainer::lengthen_list()
{
  void **tempArray = 0;
  if (listLength == 0 && (int)listIncrement == INCREMENT)
  {
    poolAllocated = CUBIT_TRUE;
    listArray = (void **)new ArrayMemory;
    assert(listArray != 0);
    listLength = listIncrement;
    return;
  }
  
  // Normal allocation from here to end.

  tempArray = new void* [listLength + listIncrement];
  currentAllocatedMemory += (listLength + listIncrement);
  if (currentAllocatedMemory > maxAllocatedMemory)
    maxAllocatedMemory = currentAllocatedMemory;
  assert(tempArray != 0);
  
  // copy the old stuff into the new array
  if (listLength)
  {
    memcpy ( tempArray, listArray, listLength*sizeof(void*) );
    // delete the old array
    if (poolAllocated)
    {
      delete ((ArrayMemory*) listArray);
    }
    else {
      delete [] listArray;
      currentAllocatedMemory -= listLength;
    }
  }
  
  // set the new equal to the old
  listArray = tempArray;
  
  // update data
  poolAllocated = CUBIT_FALSE;
  listLength += listIncrement;
  listIncrement *= 2;
  listIncrement = listIncrement > 1000000 ? 1000000 : listIncrement;
}

ArrayBasedContainer&
ArrayBasedContainer::operator=(const ArrayBasedContainer& from)
{
  if (this == &from)
    return *this;
  
  // See if the existing listArray length is equal to the from listArray
  // length. If it is, we can reuse it without deleting and newing.
  if (listLength != from.listLength) {
    if (listLength)
      if (poolAllocated) {
	delete ((ArrayMemory*) listArray);
      } else {
	delete [] listArray;
	currentAllocatedMemory -= listLength;
      }
    
    listLength = from.listLength;
    if (listLength)
      {
	if (listLength == INCREMENT && (int)listIncrement == INCREMENT)
	  {
	listArray = (void **)new ArrayMemory;
	poolAllocated = CUBIT_TRUE;
	  }
	else
	  {
	    listArray = new void* [listLength];
	    poolAllocated = CUBIT_FALSE;
	    currentAllocatedMemory += listLength;
	    if (currentAllocatedMemory > maxAllocatedMemory)
	      maxAllocatedMemory = currentAllocatedMemory;
	  }
	assert(listArray != 0);
      }
  }
  itemCount = from.itemCount;
  listIncrement = from.listIncrement;
  listIsSorted  = from.listIsSorted;
  
  if (itemCount)
    memcpy (listArray, from.listArray, itemCount*sizeof(void*));
  
  return *this;
}

ArrayBasedContainer&
              ArrayBasedContainer::operator+=(const ArrayBasedContainer& from)
{
  int tmp_itemCount = itemCount + from.itemCount;
  if (tmp_itemCount >= listLength)
    {
      if ((listLength + (int)listIncrement) < tmp_itemCount)
	{
	  listIncrement = tmp_itemCount - listLength;
	}
      lengthen_list();
    }
  if (from.itemCount == 1)
    {
      listArray[itemCount] = from.listArray[0];
    }
  else
    {
      if (from.itemCount)
	memcpy (&listArray[itemCount], from.listArray,
		from.itemCount*sizeof(void*));
    }
  itemCount += from.itemCount;
  listIsSorted  = listIsSorted && from.listIsSorted;
  return *this;
}

ArrayBasedContainer& 
ArrayBasedContainer::operator-=(const ArrayBasedContainer& from)
{
  if ( !itemCount )
    return *this;
  int i, j;
  void *item;
  for ( i = 0; i < itemCount; i++ )
  {
    item = listArray[i];
    for ( j = from.itemCount; j--; )
      if ( from.listArray[j] == item ) 
      {
        listArray[i] = NULL;
        break;
      }
  }
    // compress, virtual
  compress();
  return *this;
}

void ArrayBasedContainer::compress()
{
  int i, j;
  void *item;
  for ( i = j = 0; i < itemCount; i++ )
    if ((item = listArray[i]) != NULL) // = is ok
      listArray[j++] = item;
}

unsigned int ArrayBasedContainer::current_allocated_memory()
{
  return currentAllocatedMemory;
}

unsigned int ArrayBasedContainer::maximum_allocated_memory()
{
  return maxAllocatedMemory;
}

