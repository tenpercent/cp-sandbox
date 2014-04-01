//- Class:        CubitStack
//- Description:  A generic base class to create a stack of pointers
//-               to any type of object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#ifndef   CUBIT_STACK_HPP
#define   CUBIT_STACK_HPP

#include "CubitContainer.hpp"
#include "CubitUtilConfigure.h"

const int STACK_SIZE_INCREMENT = 100;


class CUBIT_UTIL_EXPORT CubitStack: public CubitContainer
{
  public:

    CubitStack();
    CubitStack(int increment);
    virtual ~CubitStack();
   
#ifdef BOYD15
    void reset();
#endif

  protected:

    void  base_push(void* data);
    void* base_pop();

  private:

    void initialize(int increment=STACK_SIZE_INCREMENT);
    void lengthen_stack();

    void** stackArray;

    int stackSizeIncrement;
    int stackSize;          //- the total number of items the stack can
                            //- hold after a lengthen_stack() operation
};


#define Stackdeclare(name, typePtr)                                        \
     class name : public CubitStack                                        \
{                                                                          \
public:                                                                    \
  void    push(typePtr objPtr) { base_push( (void*)objPtr );  }            \
  typePtr pop()                { return (typePtr) base_pop(); }            \
}


#endif // CUBIT_STACK_HP

