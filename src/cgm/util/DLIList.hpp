//- Class: DLIList
//- Description: DLIList is a doubly linked intrinsic list class that accepts
//-              any type of object.  It uses a template to allow the same
//-              list class definition for any data type.  Because the list
//-              makes a copy of each item added to the list, this type of
//-              list should only be used for simple data types such as
//-              integers, doubles, and enumerated types.  Other data types
//-              should use the DLList, which keeps list of pointers.
//-              The list is implemented as an array that is grown by a
//-              specified amount when the list becomes full.
//-              Insertions and deletions at any point other than the end
//-              of the list are handled inefficiently at the current time
//-              since all data from that point to the end is bubbled
//-              down/up to fill/create the void. Operators {+} and {+=} 
//-              are provided as an efficient means of assigning one list
//-              to another or appending one list onto another.
//-
//- Assumptions: All data are stored contiguously in the array with 
//-              empty slots at the end of the array.
//-
//- Owner: Darryl Melander

#ifndef DLILIST_HPP
#define DLILIST_HPP

#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include <string.h>
#include <assert.h>

#define DLI_COUNT_INCREMENT 8
#define DLI_COUNT_FACTOR 1.5
template<class X> class DLIListIterator;

template<class X> class DLIList
{
public:
  friend class DLIListIterator<X>;
  
  explicit DLIList (int list_size = 0);
    //- Constructor: Create a list with initial size 'size'. Memory for the
    //- list is not allocated until the first element is added to the list. 
  
  DLIList(const DLIList<X>& from);  //- Copy Constructor
  
  ~DLIList();
    //- Destructor: Free all resources used by this list.
  
  void step();
    //- Move the pointer to the next element in the list. If pointer reaches
    //- the end of the list, wrap around to the beginning
  
  void step(int n);
    //- Move the pointer {n} positions forward in the list wrapping 
    //- around to the beginning of the list if the end is reached
    //- If {n} is less than zero, move backward.
  
  void back();
    //- Move the pointer to the previous element in the list. If reaches
    //- the beginning of the list, wrap around to the end
  
  void back(int n);
    //- Move the pointer {n} positions backward in the list wrapping 
    //- around to the end of the list if the beginning is reached.
    //- If {n} is less than zero, move forward.
  
  void reset();     //- Set the pointer to the beginning of the list.   
  void last();      //- Set the pointer to the end of the list.
  
  void clean_out();
    //- Delete all elements in the list, reset the pointer to zero.
    //- Does not release memory already allocated for the list. This
    //- call is more efficient than creating a new list repeatedly
    //- within a loop.
  
  void shrink(int k);
    //- itemCount -= k;  Doesn't change array elements.
    //- Sets the index to the new end of the list if the old_index >= 
    //- new_itemCount
  
  CubitBoolean is_at_beginning() const;
    //- returns CUBIT_TRUE if list is at beginning
  
  CubitBoolean is_at_end() const;
    //- returns CUBIT_TRUE if list is at end

  void remove_all_with_value(X val);
    //- Remove all items with value 'val'. This is similar to the DLList
    //- function pair nullify() and compress().  First, change the value
    //- of any elements you want to remove to 'val'.  Then, call this
    //- function to compress the list.
  
  void reverse();
    //- Reverse the items in the list.
  
  DLIList<X>& operator=(const DLIList<X>& from);
  DLIList<X>& operator=(const DLIListIterator<X>&);
    //- Create a copy of a list.  Most efficient way to do this.
  
  DLIList<X>& operator+=(const DLIList<X>& from);
  DLIList<X>& operator+=(void* from);
    //- Append one list to another list.  Most efficient way to do this.

  DLIList<X>& operator-=(const DLIList<X>& from);
    //- Subtract one list from another list.  Retains order.

    //Added by K. Walton 3/1/01
    //Needed to make changes from DLList
  int operator==(const DLIList<X>& from);
  int operator!=(const DLIList<X>& from);

  X & operator[](int my_index) const;
  void merge_unique(const DLIList<X>& merge_list, 
                    bool merge_list_unique = false);
    //- Merges the contents of the list, merge_list, with those of "this"
    //- list, ensuring that items that are being added do not already
    //- exist in "this" list.  If the items are known to appear in the
    //- merge_list only once, then it can be done faster.

  void intersect(const DLIList<X>& merge_list);
  void intersect(void* merge_list);
  void intersect_unordered(const DLIList<X>& merge_list);
  //- The regular version is O(n^2)
  //- The unordered version is O(n + sort_time).

  CubitBoolean append_unique(X new_item);
    //- Appends the new item to the list, if it doesn't already exist
    //- in the list. In either case, index is not changed.
    //- Return CUBIT_TRUE if the item was added, else CUBIT_FALSE.
  
  void uniquify_unordered();
    //- Ensure each element only appears once.
    //- Sorts first for speed, resets index to 0.
    
  void uniquify_ordered();
    //- Ensures each element appears only once.
    //- Preserves order.  Resets index to 0.
  
  X remove ();
    //- Removes the current item from the list.
  
  X remove (X item);
    //- Removes item with value 'item' from the list.
  
  CubitBoolean omit (X old_val);
    //- Finds instance of item by matching value and deleting next instance
    //- of it from the list, wrapping if necessary. The current position of
    //- the list is not changed.
    //- Returns CUBIT_TRUE if the item was found, CUBIT_FALSE if not.
  
  X get () const;
    //- Returns value of current item
  
  X next () const;
    //- Returns value at next position in the list
  
  X next (int n) const;
    //- Returns value at 'n' positions forward in list
  
  X prev () const;
    //- Returns value at previous position in list
  
  X prev (int n) const;
    //- Returns value 'n' positions back in list
  
  X get_and_step ();
    //- Returns value, then steps
  
  X get_and_back ();
    //- Returns value, then goes back
  
  X step_and_get ();
    //- Steps then returns value
  
//  CubitBoolean move_between (X item1, X item2);
  CubitBoolean move_to(X item);
    //- Moves to the next instance of this value, wrapping if necessary
  CubitBoolean is_in_list(X item) const;
    //- Return True if item is in the list, false otherwise.

  X pop();
    //- Returns and removes last value in list
  X push(X value);
    //- Moves to end of list and appends and returns value
  void insert (X new_item);
    //- put new item in list after current item and make it current
    //- Inefficient due to bubbling.
  void insert_first(X new_item);
    //- Adds value at start of list.  Inefficient due to bubbling.
  void append(X new_item);
    //- Adds value to end of list
  X extract();
    //- Removes current value and puts last value
    //- in the list in its place.  If list order isn't important, this is
    //- much more efficient than remove().
  
  X change_to(X new_value);
    //- Changes value of current Item.
  
  void sort();
    //- Orders the list elements from lowest to highest.
    //- Use reverse after this call if you want to go the
    //- other direction.
  
  typedef int (*SortFunction)(X& a, X& b);
  void sort(SortFunction f);
    //- Sorts based on a function you pass in.
    //- SortFunction should return a negative number if 'a' should
    //- be before 'b', a positive number if 'b' comes before 'a', and
    //- zero if they are equal, or relative order doesn't matter.
  
  int size() const
    { return itemCount; }
    //- Returns the number of items in the list

  int get_index()
     { return index; }
    //- Returns current index value
  
  int where_is_item(X item) const;
    //- Return index of item in list, or -1 if item is not in the list.

   int memory_use(CubitBoolean verbose_boolean);
   //- return memory allocated in bytes
   
  void copy_to(X *other_array);
    //- copy this list's listArray into other_array

  void copy_from(X *other_array, const int other_size);
    //- copy items from other_array to this list

  CubitBoolean move_to_nearby(const X objPtr);
    //- searches from the current item both forward and backward.

  int distance_to_nearby(const X objPtr);
    //- Sets indes to the item in the list that matches {body}

  CubitBoolean move_between(X objPtr1, X objPtr2);
    //- Sets the index to the place matching one of the items,
    //- where the previous item is the other item

private:
  void lengthen_list(int by_how_much = DLI_COUNT_INCREMENT,
                     double by_what_factor = DLI_COUNT_FACTOR );
    //- Makes the array bigger. Multiply current size 
    //- "by_what_factor" and add 'by_how_much'.
  
  int  index;      // Index of current item in {listArray}
  X   *listArray;  // Array of values
  int  listLength; // Number of allocated spaces in array
  int  itemCount;  // Number of items stored in array
  
};

// *** inline function definitions *******************************************
template <class X> inline 
DLIList<X>& DLIList<X>::operator=(const DLIListIterator<X>& from_iter)
{ 
  if ( from_iter.dlIList ) 
    *this = *(from_iter.dlIList); 
  else
    clean_out();
  return *this; 
}

template <class X> inline
X & DLIList<X>::operator[](int my_index) const
{
    assert( my_index >= 0 && my_index < itemCount );
    return listArray[my_index];
}

template <class X> inline void DLIList<X>::step()
{
  if (itemCount)
    index++; 
  if (index >= itemCount)
    index = 0;
}

template <class X> inline void DLIList<X>::step(int n)
{
  if (itemCount)
  {
    index += n;
    if (index < 0)
      index = itemCount - (-index)%itemCount;
    if (index >= itemCount)
      index %= itemCount;
  }
}

// Chop off the last (k) items in the list.
template <class X> inline void DLIList<X>::shrink(int k)
{
  if (k > 0) 
  {
    if (itemCount > k)
      itemCount -= k;
    else
      itemCount = 0;
  }
  if (index >= itemCount)
  {
    if (itemCount)
      index = itemCount-1;
    else
      index = 0;
  }
}

template <class X> inline void DLIList<X>::back()
{
  if (itemCount)
  {
    index--;
    if (index < 0)
      index = itemCount - 1;
  }
}

template <class X> inline void DLIList<X>::back(int n)
{
  step(-n);
}

template <class X> inline void DLIList<X>::reset()
{
  index = 0;
}

template <class X> inline void DLIList<X>::last()
{
  if (itemCount)
    index = itemCount - 1;
}

template <class X> inline CubitBoolean DLIList<X>::is_at_beginning() const
{
  if (!index)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

template <class X> inline CubitBoolean DLIList<X>::is_at_end() const
{
  if (itemCount == 0 || index == itemCount-1)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

template <class X> inline void DLIList<X>::clean_out()
{
  itemCount = 0;
  index = 0;
}

template <class X> inline CubitBoolean DLIList<X>::is_in_list(X item) const
{
  return where_is_item(item) >= 0 ? CUBIT_TRUE : CUBIT_FALSE;
}

//- Add item to end of list, do not change current index value.
//- This function used to reduce time required to do an insert 
//- followed by a move_to back to where we were.
template <class X> inline void DLIList<X>::append(X new_item)
{
    // see if the list must be lengthened
  if (itemCount == listLength)
    lengthen_list();
  
  listArray[itemCount++] = new_item;
}

template <class X> inline X DLIList<X>::get() const
{
  if ( !itemCount )
  {
    PRINT_WARNING("Attempted get of empty DLIList\n");
    return X(0);
  }
  else
  {
    return listArray[index];
  }
}

template <class X> inline X DLIList<X>::get_and_step()
{
#ifndef GDSJAAR
  if ( !itemCount )
  {
    PRINT_WARNING("Attempted get_and_step from empty DLIList\n");
    return X(0);
  }
  else
#endif
  {
    X temp = listArray[index++];
    if (index == itemCount)
      index=0;
    return temp;
  }
}

template <class X> inline X DLIList<X>::get_and_back ()
{
   if ( !itemCount )
   {
      PRINT_WARNING("Attempted get_and_back from empty DLIList\n");
      return X(0);
   }
   X temp = listArray[index--];
   if (index < 0)
      index=itemCount-1;
   return temp;
}

template <class X> inline X DLIList<X>::next() const
{
  if (!itemCount)
  {
    PRINT_WARNING("Attempted next of empty DLIList\n");
    return X(0);
  }
  else if (index == itemCount-1)
    return listArray[0];
  else
    return listArray[index+1];
}

template <class X> inline X DLIList<X>::next(int n) const
{
  if ( !itemCount )
  {
    PRINT_WARNING("Attempted next of empty DLIList\n");
    return X(0);
  }
  else
  {
      // return the proper index
      // beware of negative n leading to negative new_index
    int new_index = index+n;
    while (new_index < 0)
      new_index += itemCount;
    if (new_index >= itemCount)
      new_index %= itemCount;
    
    return listArray[new_index];
  }
}

template <class X> inline X DLIList<X>::prev() const
{
  return this->next(-1);
}

template <class X> inline X DLIList<X>::prev(int n) const
{
  return this->next(-n);
}

//- put new item in list at index 0 and make it current
template <class X> inline void DLIList<X>::insert_first(X new_item)
{
  // set index to -1 ... index will be set to 0 in insert
  index = -1;
  insert(new_item);
}

template <class X> inline CubitBoolean DLIList<X>::move_to (X item)
{
  int item_index = where_is_item(item);
  CubitBoolean rv = CUBIT_FALSE;
  if (item_index >= 0)
  {
    index = item_index;
    rv = CUBIT_TRUE;
  }
  return rv;
}

template <class X> inline X DLIList<X>::change_to (X new_value)
{
  X temp = X(0);
  if ( !itemCount )
  {
    PRINT_WARNING("DLIList: Attempted update of empty list\n");
  }
  else
  {
    temp = listArray[index];
    listArray[index] = new_value;
  }
  return temp;
}

// Removes every instance of 'val' in list.
// Set to beginning of list after this call.
template <class X> inline void DLIList<X>::remove_all_with_value(X val)
{
  int j = 0;
  int i = 0;
  
  for ( ; i < itemCount; i++)
    if (listArray[i] != val && j++ != i)
      listArray[j-1] = listArray[i];
  
  itemCount = j;
  index = 0;
}

// Searches for item and removes the next instance from the list.
// If the item was found and removed it is returned, else 0 is returned.
template <class X> inline X DLIList<X>::remove (X item)
{
  X temp = X(0);
  if (move_to(item))
    temp = remove();
  return temp;
}

template <class X> inline int DLIList<X>::where_is_item (X item) const
{
  if (itemCount == 0)
    return -1;
  
    // loop through list searching for item ...
    // if found return index
  
    // Search from current index to end of array
  int i;
  for (i=index; i < itemCount; i++)
    if (listArray[i] == item)
      return i;
  
    // Now search from beginning of array to index...
  for (i = 0; i < index && i < itemCount; i++)
    if (listArray[i] == item)
      return i;
  
    // item is not in array, return -1
  return -1;
}
//- Append this new_item to the list if it doesn't already exist in it.
//- Make sure the index is unchanged.
template <class X> inline CubitBoolean DLIList<X>::append_unique(X new_item)
{
    // Append new_item, if it isn't already there.
  if( where_is_item(new_item) < 0 ) 
  {
    append (new_item);
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}
// Added by CAT for NT port
#if defined(TEMPLATE_DEFS_INCLUDED)
  #include "DLIList.cpp"
#endif


template <class X> class DLIListIterator
{
  friend class DLIList<X>;
public:
    // constructor
  DLIListIterator() :
      mIndex(0),
      dlIList(NULL)
    {}
  
    // destructor
  virtual ~DLIListIterator()
    {}
  
  void watch( const DLIList <X> *dl_i_list )
    {
      dlIList = dl_i_list;
      mIndex = dlIList->index;
    }
  
  void reset()
    { mIndex = 0; }
  
  void step( int n = 1 ) 
    {
      if (!size())
        mIndex = 0;
      else
      {
        mIndex += n;
        while (mIndex < 0)
          mIndex += dlIList->itemCount;
        mIndex %= dlIList->itemCount;
      }
    }
  
  int size() const
    { return dlIList ? dlIList->size() : 0; }
  
  X get() const
    { return next(0); }
  X next( int n = 1 ) const;
  X get_and_step( int n = 1 ) 
    {
      X temp = get();
      step( n );
      return temp;
    }
  
    //- Moves to the next instance of this value, wrapping if necessary
  CubitBoolean move_to(X item)
    { 
      mIndex = where_is_item( item );
      return mIndex >= 0 ? CUBIT_TRUE : CUBIT_FALSE; 
    }
  
    //- Return True if item is in the list, false otherwise.
  CubitBoolean is_in_list(X item) const
    { return where_is_item(item) >= 0 ? CUBIT_TRUE : CUBIT_FALSE; }
  
    // Returns CUBIT_TRUE if we're pointing at the last item in the list.
  CubitBoolean is_at_end() const
    {
      if (size() == 0 || mIndex == size() - 1)
        return CUBIT_TRUE;
      return CUBIT_FALSE;
    }
  
    //- returns CUBIT_TRUE if we're pointing at the firt item in the list.
  CubitBoolean is_at_beginning() const
    {
      if (!mIndex || mIndex == size())
        return CUBIT_TRUE;
      return CUBIT_FALSE;
    }
  
protected:
    // utility functions
  int where_is_item( X item ) const
    { return dlIList ? dlIList->where_is_item( item ) : -1; }
  
    // data
  int mIndex;
  const DLIList<X> *dlIList;
};

template <class X> 
inline X DLIListIterator<X>::next( int n ) const
{
  int size_loc = size();
  if ( size_loc == 0 )
    return NULL;
  int index_loc = mIndex + n;
  while ( index_loc < 0 )
    index_loc += size_loc;
  while ( index_loc >= size_loc )
    index_loc -= size_loc;
  
    // if dlIList == NULL, then size == 0 and returned already 
  return dlIList->listArray[index_loc];
}

#endif //- DLILIST_HPP

