#include "assert.h"
#include "PriorityQueue.hpp"
#include "CubitMessage.hpp"


#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

template <class Object> MY_INLINE
PriorityQueue<Object>::PriorityQueue(typename PriorityQueue<Object>::CompareFunc compare)
    : theItems(1)
{
  lessThanFunc = compare;
}
template <class Object> MY_INLINE
int PriorityQueue<Object>::size() const
{
  return theItems.size() - 1;
}
template <class Object> MY_INLINE
bool PriorityQueue<Object>::empty() const
{
  if ( size() == 0 )
    return true;
  else
    return false;
}
template <class Object> MY_INLINE
const Object PriorityQueue<Object>::top() const
{
  if ( empty() )
  {
    PRINT_ERROR("Empty priority queue had top request.\n");
    assert(!empty());
    return static_cast<Object>(0);
  }
  return theItems[1];
}
template <class Object> MY_INLINE
void PriorityQueue<Object>::push(Object x)
{
  theItems.push_back(x);
  theItems[0] = x;//stash the item here...
    //Percolate up.
  int hole = size();
    //While x is less than the parrent, do (1): then go up the tree (hole /= 2 is the parent...).
  for ( ; lessThanFunc(x, theItems[hole/2]); hole /= 2 )
  {
      //(1) set the replace the parent with the child.
    theItems[hole] = theItems[hole/2];
  }
    //Now where the item is greater than the parent, set this value.
  theItems[hole] = x;
}
template <class Object> MY_INLINE
void PriorityQueue<Object>::pop()
{
  if ( empty() )
  {
    PRINT_ERROR("Trying to delete min from empty priority queue.\n");
    assert(!empty());
    return;
  }

  int hole = 1;
  int child;

  Object tmp = theItems.back();
  theItems.pop_back();
  int the_size = size();
    //while the current isn't at the end, do (1): then set hole to child.
  for ( ; hole*2 <= the_size; hole = child )
  {
    child = hole*2;
    if ( child != the_size &&
         lessThanFunc(theItems[child+1], theItems[child]) )
      child++;
    if ( lessThanFunc(theItems[child], tmp) )
      theItems[hole] = theItems[child];
    else
      break;
  }
  if ( !empty() )
    theItems[hole] = tmp;
}

template <class Object> MY_INLINE
int PriorityQueue<Object>::where_is_item(Object &a, int start_index, bool queue_is_valid)
{
  int i;

  if (!queue_is_valid)
  {
    for (i=0; i < theItems.size(); ++i)
    {
      if (theItems[i] == a) {return i;}
    }
    return -1;
  }
  
  if (a == theItems[start_index]) {return start_index;}
  if (!lessThanFunc(a, theItems[start_index]) && 2*start_index+1 <= size())
  {
    int result;
    result = where_is_item(a, 2*start_index);
    if (result > 0) {return result;}
    result = where_is_item(a, 2*start_index+1);
    if (result > 0) {return result;}
  }
  
  return -1;
}

template <class Object> MY_INLINE
bool PriorityQueue<Object>::update_item(Object &a, bool queue_is_valid)
{
  int index = where_is_item(a, queue_is_valid);
  if (index == -1) {return false;}
  else {return update_item(index);}
}

template <class Object> MY_INLINE
bool PriorityQueue<Object>::update_item(int object_index)
{
  if (object_index == 0) {return true;}
  
  if (lessThanFunc(theItems[object_index], theItems[object_index/2]))
  {
    Object temp = theItems[object_index];
    theItems[object_index] = theItems[object_index/2];
    theItems[object_index/2] = temp;
    update_item(object_index/2);
    return true;
  }

  int smaller_child = 2*object_index;
  if (2*object_index > size()) {return true;}
  if (2*object_index+1 <= size())
  {
    if (lessThanFunc(theItems[2*object_index+1],theItems[2*object_index])) {++smaller_child;}
  }
  
  if (lessThanFunc(theItems[smaller_child], theItems[object_index]))
  {
    Object temp = theItems[object_index];
    theItems[object_index] = theItems[smaller_child];
    theItems[smaller_child] = temp;
    update_item(smaller_child);
  }

  return true;
}

template <class Object> MY_INLINE
bool PriorityQueue<Object>::validate(int index)
{
  bool fine = true, kid_one = false, kid_two = false;
  if  (2*index <= size())
  {
    kid_one = true;
    if ((lessThanFunc(theItems[2*index], theItems[index])))
    {
      fine = false;
    }
  }
  
  if  (2*index+1 <= size())
  {
    kid_two = true;
    if ((lessThanFunc(theItems[2*index+1], theItems[index])))
    {
      fine = false;
    }
  }
  
  if (fine)
  {
    if (kid_one && !kid_two) {return validate(2*index);}
    if (!kid_one && kid_two) {return validate(2*index+1);}
    if (kid_one && kid_two) {return (validate(2*index) && validate(2*index+1));}
    else {return true;}
  }
  
  return false;
  
}
