// Classes: Queue, QueueNode

#include "Queue.hpp"


Queue::~Queue() 
{
  // delete all QueueNodes in this Queue
  QueueNode *old_head;
  while (head) {
    old_head = head;
    head = head->next;
    delete old_head;
  }
  // free up memory
  QueueNode::memoryManager.compress_memory();
}

void Queue::append_queue( void *add_data ) 
{
  QueueNode *node = new QueueNode(add_data);
  // if the last element was popped, then
  // tail will be non-null but bogus, but head will be null.
  if (head) {    
    tail->next = node;
    tail = node;
  }
  else {
    head = tail = node;
  }  
}

void *Queue::pop_queue() 
{
  if (head) {
    void *return_data = head->data;
    QueueNode *old_head = head;
    head = head->next; //may be null, tail not set for efficiency
    delete old_head;
    return return_data;
  }
  else
    return NULL;
}
