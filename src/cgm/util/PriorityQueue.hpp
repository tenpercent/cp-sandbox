#ifndef PRIORITYQUEUE_HPP
#define PRIORITYQUEUE_HPP
#include <vector>

template <class Object>
class PriorityQueue
{
public:
  typedef bool (*CompareFunc)(Object &a, Object &b);
  PriorityQueue(CompareFunc compare);
  int size() const;
  bool empty() const;

  const Object top() const;
  void push(Object x);
  void pop();

  int where_is_item(Object &a, int start_index=1, bool queue_is_valid=false);
  bool update_item(Object &a, bool queue_is_valid=false);
  bool update_item(int object_index);

  bool validate(int index=1);
  
private:
  std::vector<Object> theItems;
  CompareFunc lessThanFunc;
};

#if defined(TEMPLATE_DEFS_INCLUDED)
  #include "PriorityQueue.cpp"
#endif


#endif

