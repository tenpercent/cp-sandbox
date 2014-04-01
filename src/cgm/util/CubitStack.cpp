//- Class:        CubitStack
//- Description:  A generic base class to create a stack of pointers
//-               to any type of object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#include <string.h>
#include "CubitStack.hpp"


CubitStack::CubitStack()
{
  initialize();
}

CubitStack::CubitStack(int increment)
{                
  (increment < 1) ? initialize(increment): initialize();
}

void CubitStack::initialize(int increment)
{
  stackSizeIncrement = increment;
  stackArray         = 0;                      
  numberItems        = 0;
  stackSize          = 0;
}

CubitStack::~CubitStack()
{
  if (stackSize)
  {
    delete [] stackArray;
  }
}

void CubitStack::lengthen_stack()
{
  void** temp_array = new void* [stackSize + stackSizeIncrement];
  if (stackSize)
  {
    memcpy( temp_array, stackArray, stackSize*sizeof(void*) );
    delete [] stackArray;
  }
  stackArray  = temp_array;
  stackSize  += stackSizeIncrement;
}

