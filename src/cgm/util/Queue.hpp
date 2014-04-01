//- Class: Queue
//- Description: Queue is a queue class that accepts generic
//-              pointers to any object.  It is controlled by a macro that
//-              substitutes the specific pointer definitions in the class
//-              declaration. 
//+
//-              The macro {Queuedeclare(name,typePtr)} is defined to
//-              create specific instances of the Queue class.
//-              {name} is the name of the Queue class created, and
//-              {typePtr} is the type of elements stored in the list.
//+
//-              Functions supported are pop, append, and a built in
//-              iterator to step through the items of the list (read only).
//+
//-              The items in the Queue are kept in "QueueNode" class
//-              objects, which are completely hidden except to Queue.
//- Assumptions: 
//-
//- Owner: Scott Mitchell
//- Checked by: 
//- Version: $Id: 

#ifndef QUEUE_HPP
#define QUEUE_HPP

#include "MemoryManager.hpp"
#include "CubitUtilConfigure.h"

class Queue;

class CUBIT_UTIL_EXPORT QueueNode
{
  private:

    static MemoryManager memoryManager;
    //- memory management object

  protected:

    QueueNode *next;
    //- next item in the containing Queue, NULL terminates Queue.

    void *data;  
    //- the node data
  
    QueueNode(void* set_data);
    //- constructor, can only be called by Queue
  
    SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators, using memoryManager
  
  public:

    ~QueueNode( ) {}

    friend class Queue;  
};

inline
QueueNode::QueueNode(void* set_data)
{
  next = NULL;
  data = set_data;
}

class CUBIT_UTIL_EXPORT Queue
{
  private:

    QueueNode *head;
    QueueNode *tail;
    QueueNode *iterator;
    //- QueueNode data

  public:

    Queue();
    //- constructor

    virtual ~Queue();
    //- delete this and all QueueNodes in the Queue. Calls
    //- QueueNode::memoryManager to compress *all* Queues in existence.
  
    void append_queue(void *add_data);
    //- add to end of queue

    void *pop_queue();
    //- remove from front of queue

    void first();
    void step();
    void *queue_get();
    void *queue_get_and_step();
    //- built in iterator for stepping through the queue in sequence
};

inline 
Queue::Queue()
{ 
  iterator = head = tail = NULL; 
}

inline
void* Queue::queue_get() 
{ 
  return iterator ? iterator->data : NULL; 
}

inline
void* Queue::queue_get_and_step()
{
  if (iterator)
  {
    QueueNode* iterator_old = iterator;

    iterator = iterator->next;
    return iterator_old->data;
  }
  else
  {
    return NULL;
  }
}

#define Queuedeclare(name, typePtr)                                          \
                                                                             \
class name : public Queue                                                    \
{                                                                            \
  public:                                                                    \
  name() : Queue()                  {}                                       \
  typePtr  pop()                    {return (typePtr) pop_queue();}          \
  void     append(typePtr add_data) {append_queue((void*) add_data);}        \
  typePtr  get()                    {return (typePtr) queue_get();}          \
  typePtr  get_and_step()           {return (typePtr) queue_get_and_step();} \
}                                                                            \

#endif //- QUEUE_HPP


