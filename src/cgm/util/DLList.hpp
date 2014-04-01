//- Class: DLList
//- Description: DLList is a doubly linked list class that accepts generic
//-              pointers to any object.  It is controlled by a macro that
//-              substitutes the specific pointer definitions in the class
//-              declaration. The list is implemented as an array that is 
//-              grown by a specified amount when the list becomes full.
//-              Insertions and deletions at any point other than the end
//-              of the list are handled inefficiently at the current time
//-              since all data from that point to the end is bubbled
//-              down/up to fill/create the void. Operators {+} and {+=} 
//-              are provided as an efficient means of assigning one list
//-              to another or appending one list onto another.
//-
//-              The macro {DLListdeclare(name,typePtr)} is defined to
//-              create specific instances of the DLList class.
//-              {name} is the name of the DLList class created, and
//-              {typePtr} is the type of elements stored in the list.
//-              The macro also defines the functions {push(item)} and 
//-              {pop()} which allow {DLLists} to perform like stacks.
//-
//- Assumptions: All data are stored contiguously in the array with 
//-              empty slots at the end of the array.
//-
//- Owner: Greg Sjaardema
//- Checked by: Jim Hipp, 3/28/94
//- Version: $Id: 

#ifndef DLLIST_HPP
#define DLLIST_HPP

#include "ArrayBasedContainer.hpp"
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include <string.h>
#include <assert.h>
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT DLList : public ArrayBasedContainer
{
public:

  DLList (int list_size);
  //- Constructor: Create a list with initial size {list_size}. The list
  //- will be grown by {list_size} each time it is filled. Memory for the
  //- list is not allocated until the first element is inserted using
  //- {insertLink}. 
  //- If {list_size} is zero, the default increment ({INCREMENT}) will be used
  //- From an efficiency standpoint, it is very important that the 
  //- increment be set large enough to reduce the number of list 
  //- growths, but small enough to not waste memory.
  //- It is more efficient to sligthly overestimate the size than 
  //- to underestimate the size.
  
  DLList(const DLList& copy_from);  //- Copy Constructor
  
  virtual ~DLList();
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
  //- call is more efficient creating a new list repeatadly within
  //- a loop.
  
  void shrink(int k);
  //- itemCount -= k;  Doesn't change array elements.
  //- Sets the index to the new end of the list if the old_index >= 
  //- new_itemCount
  
  CubitBoolean is_at_beginning() const;
  //- returns CUBIT_TRUE if list is at beginning
  
  CubitBoolean is_at_end() const;
  //- returns CUBIT_TRUE if list is at end
  
  int get_index();
  //- Returns the current index value.  This function and set_index are
  //- included for efficiency when one wants to return to a previous location 
  //- in the list when it is known that the list has not changed.  This
  //- is the only recommended usage for these functions.
  
  void set_index(int new_index);
  //- Sets the list to the given index.  If the requested index is less than
  //- zero or out of the scope of the number of items in the list, an assert 
  //- occurs.  This function and get_index are included for efficiency when 
  //- one wants to return to a previous location in the list when it is known 
  //- that the list has not changed.  This is the only recommended usage for 
  //- these functions.
  
  void compress();
  //- Remove all null (0) pointers from list. Must be called after
  //- {nullify_link()} is called before calling any list modification 
  //- function except {move_to()} or {nullify_link()}
  
  void reverse();
  //- Reverse the items in the list.
  
  DLList& operator=(const DLList&);
  //- Create a copy of a list.
  
  virtual void merge_unique(DLList& merge_list, 
                            int merge_list_unique = CUBIT_FALSE);
  //- Merges the contents of the list, merge_list, with those of "this"
  //- list, ensuring that items that are being added do not already
  //- exist in "this" list. If the items are known to appear in the
  //- merge_list only once, then it can be done faster.

  void intersect(DLList& merge_list);
  //- Merges the contents of the list merge_list, with those of "this"
  //- list, but only the items that are common to both lists are kept
  //- in "this" list.  This is a rather costly process I think.
  //- if merge_list is the same as the calling list, it just returns
  
  int append_unique(void* new_item);
  //- Appends the new item to the list, if it doesn't already exist
  //- in the list. In either case, index is not changed.
  //- Return CUBIT_TRUE if the item was added, else CUBIT_FALSE.

  SetDynamicMemoryAllocation(memoryManager)
  //- overloaded new and delete operators
	
	//Added by J. Kraftcheck, 5/25/99
	//These are not very efficient (O(n^2))
	//In order to maintain the correct functionality of these
	//operators when either list contains duplicates, it is
	//also necessary to create a temporary list during the check.
	int operator< ( const DLList& list ) const; //subset
	int operator<=( const DLList& list ) const; //subset or equivalent
	int operator> ( const DLList& list ) const; //superset
	int operator>=( const DLList& list ) const; //superset or equivalent
	int operator==( const DLList& list ) const; //non-order-dependend compare
	int operator!=( const DLList& list ) const; //non-order-dependend compare

protected:
  
  void insert_link(void* newBody);
  //- put new item in list after current item and make it current
  
  void insert_link_first(void* new_item);
  //- put new item in list at index 0 and make it current
  
  void append_link(void* newBody);  
  //- put new item at end of list and keep pointer where it is.
  //- This call should be used to add items to a list if order 
  //- of items is not important.
  
  void *nullify_link();
  //- Change the current item to a null pointer and return a pointer
  //- to the item (before nulled). See the discussion for {nullItem}.
  
  void *cut_link();
  //- remove the item at the current location and return a pointer to it.
  //- The next node becomes the current node
  //- Returns {NULL} if there are no items in list
  
  void* extract_link();
  //- remove the item at the current location and return a pointer to it.
  //- used for efficiency in cases where preservation of list order is not
  //- important.  moves last list item (itemCount - 1) to current index and
  //- decrements itemCount.  eliminates the need to perform the list bubble
  //- down (i.e. cut_link) but sacrifices list order in the process.  this
  //- function should not be used when up-stream order from the removed node
  //- is important.  when processing a list using this function the user
  //- should reset the list to the head (index = 0) before beginning to
  //- ensure all list nodes are processed properly.
  //- Returns {NULL} if there are no items in list
  
  CubitBoolean omit_link(void*);
  //- Finds instance of item by matching pointer and deleting all instances
  //- of it from the list. The current position of the list is not changed.
  //- Returns CUBIT_TRUE if the item was found and deleted, CUBIT_FALSE
  //- if not.
  
  void *next_item(int n=1) const;
  //- Return a pointer to the item {n} items after the current item.
  //- Index is not changed.
  //- Returns {NULL} if there are no items in list
  
  void *prev_item(int n=1) const;
  //- Return a pointer to the item {n} items before the current item.
  //- Index is not changed.
  //- Returns {NULL} if there are no items in list
  
  void *get_item() const;
  //- Return a void pointer to the current item.
  //- Index is not changed.
  //- Returns {NULL} if there are no items in list

  void *change_item_to(void* newBody);

  void *move_to_and_remove_item(void* item);
  //- search for item in list and remove if found.  Function uses
  //- move_to_item to find item and cut_link to remove it.  If item is found
  //- it is returned else NULL is returned.

  void *get_item_and_step();
  //- Return a pointer to the current item and increment the index by 1.
  //- Returns {NULL} if there are no items in list
  
  void *step_and_get_item();
  //- Increment the index by 1, then return a pointer to the current item.
  //- Returns {NULL} if there are no items in list
  
  void *get_item_and_back();
  //- Return a pointer to the current item and decrement the index by 1.
  //- Returns {NULL} if there are no items in list
  
  CubitBoolean move_to_item(void* body);
  CubitBoolean move_to_nearby_item(void* body);
  int distance_to_nearby_item(void* body);
  //- Sets the index to the item in the list that matches {body}
  //- Returns {CUBIT_FALSE} if there are no items in list or if 
  //- the requested item is not in the list.
  //- Returns {CUBIT_TRUE} if the item is Found.
    //- move_to_item searches from the current item forward
    //- move_to_nearby_item searches from the current item both
    //- forward and backward.
  
  CubitBoolean move_between_items(void *body1, void *body2);
  //- Sets the index to the place matching one of the items, where the
  //- previous item is the other item (i.e. no order assumed on the
  //- parameters)
  //- Returns {CUBIT_FALSE} if there are no items in list or the items
  //- are not consecutive in the list.
  //- Returns {CUBIT_TRUE} if the consecutive items are found.
  
  //void *change_item_to(void* newBody); //commented out 8-12-98
  //- Sets the item at the current index value to {newBody} and
  //- returns the old value at that location.
  //- Returns {NULL} if there are no items in list
  
  void *nullItem;
  //- Item removed in last {nullify_link} operation. It is saved here
  //- only for execution speed reasons. The primary use of {nullify_link}
  //- is in deleting a partial mesh. After {nullify_link()} is called,
  //- a destructor for that item is called which does a {move_to}.
  //- Since we just nulled that entry, we can't {move_to()} it, so
  //- {move_to()} checks to see if the item matches
  //- {nullItem}. It it does, we return, otherwise we do a normal 
  //- {move_to()}.
  
  int where_is_item( void *item ) const;
  //- Return index of item in list, or -1 if item is not in the list.
  //- This information should never be passed outside the class.

  int  index;                  //- index of current item in {listArray}
  
private:
  
  // int  listIsSorted;  //- This is now a bit in ArrayBasedContainer.hpp
  //- flag used to indicate if an ordered list is no longer ordered
  
  static const char* type() { return "DLList"; }
  //- declare which type of ArrayBasedContainer this is
  
  static MemoryManager memoryManager;
  //- memory management static objects for performing block allocation
};

// *** inline function definitions *******************************************

inline void DLList::step()      { if (itemCount) index++; 
                                  if (index >= itemCount) index = 0;}

inline void DLList::step(int n)  {
  if (itemCount) {
    index += n;
    if (index < 0)   {
      index = itemCount - (-index)%itemCount;
    }
    if (index >= itemCount) {
      index %= itemCount;
    }
  }
}


inline void DLList::back()      { if (itemCount) index--;
		                  if (index < 0) index = itemCount - 1; }
  
inline void DLList::back(int n) { step(-n); }
  
inline void DLList::reset()     { index = 0; }
  
inline void DLList::last()      { if (itemCount) index = itemCount - 1; }
  
inline CubitBoolean DLList::is_at_beginning() const
{if (!index) return CUBIT_TRUE;
 else return CUBIT_FALSE;}

inline CubitBoolean DLList::is_at_end() const
{
  if (itemCount == 0 || index == itemCount-1)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}
 
//- Add item to end of list, do not change current index value.
//- This function used to reduce time required to do an insert 
//- followed by a move_to back to where we were.
 inline void DLList::append_link ( void* new_item )
{
    assert(new_item != NULL);
    // see if the list must be lengthened

    if ( itemCount == listLength )
        lengthen_list();

    listArray[itemCount++] = new_item;
}

inline void *DLList::get_item () const
{
    if ( !itemCount ) {
        PRINT_WARNING("Attempted get of empty DLList\n");
        return NULL;
    }
    else {
      assert(listArray[index] != NULL);
      return listArray[index];
    }
}

inline void *DLList::get_item_and_step () 
{
#ifndef GDSJAAR
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted get_and_step from empty DLList\n");
        return NULL;
    }
    else
#endif
    {
        void *temp = listArray[index++];
        if (index == itemCount) index=0;
        assert(temp != NULL);
        return temp;
    }
}

inline void *DLList::next_item ( int n )  const
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted next of empty DLList\n");
        return NULL;
    }
    else
    {
        // return the proper index
        // beware of negative n leading to negative new_index
        int new_index = index+n;
	while (new_index < 0)
	  new_index += itemCount;

	new_index = new_index%itemCount;

        assert(listArray[new_index] != NULL);
        return listArray[new_index];
    }
}

inline void DLList::clean_out ()
{
    itemCount = 0;
    index = 0;
}


class DLListIterator
{
public:
  DLListIterator() { index = 0; }
    // constructor
  virtual ~DLListIterator() {}
    // destructor
  
  void step( int n = 1 ) 
    { index += n; }

protected:
  
    // data
  int index;
};


// *** macro definitions *****************************************************

// Three separate derived classes are possible from DLList all of which are
// defined by one of the following macros:
// 
//    DLListdeclare(name, typePtr)
//
// The first is the standard DLList declaration macro, the second provides
// sorting/searching functionality, and the third provides sorting/searching
// functionallity for simple intrinsic data lists.  All three are described in
// more detail below.  Two support macros
//
//    CommonDefine(typePtr, notSorted)
//    CommonSortedDefine(name, typePtr)
//
// are also defined but are not meant to be utilized externally.  They simply
// collect class definitions that are common to the three list defintions
// macros described above

// The following macro defintion provides the core functionality for the
// standard DLListdeclare as-well-as the sorted SDLListdeclare macros.
// typePtr is still the type pointer of the data stored by the list.
// notSorted is simply a semi-colon (;) for the standard list (if white space
// is used for the parameter the CPP outputs a warning message that a NULL
// string is being substituted for the parameter ... the proper outcome but
// the excessive warning messages are a real nuisance) and a call to the
// set_sort_flag function (with an argument of False) for the sorted lists.
#define CommonDefine(typePtr, notSorted)                                     \
                                                                             \
    typePtr remove()         {return (typePtr) cut_link();}                  \
    int omit(typePtr objPtr) {return omit_link((void*) objPtr);}             \
    typePtr get() const      {return (typePtr) get_item();}                  \
    CubitBoolean set(typePtr objPtr) {return move_to_item((void*)objPtr);}   \
    typePtr next() const     {return (typePtr) next_item();}                 \
    typePtr next( int n ) const {return (typePtr) next_item(n);}             \
    typePtr prev() const     {return (typePtr) prev_item();}                 \
    typePtr prev(int n) const {return (typePtr) prev_item(n);}               \
    typePtr get_and_step()   {return (typePtr) get_item_and_step();}         \
    typePtr get_and_back()   {return (typePtr) get_item_and_back();}         \
    typePtr step_and_get()   {return (typePtr) step_and_get_item();}         \
                                                                             \
    CubitBoolean move_between(typePtr objPtr1, typePtr objPtr2)              \
    {                                                                        \
        if (nullItem && (((typePtr) nullItem == objPtr1) ||                  \
	                 ((typePtr) nullItem == objPtr2) ))                \
  	  return CUBIT_FALSE;                                              \
        else                                                                 \
	  return (move_between_items(((void*) objPtr1), ((void*) objPtr2))); \
    }                                                                        \
    typePtr pop()                {last(); return remove();}                  \
    typePtr push(typePtr objPtr)                                             \
        { notSorted last(); append_link((void*)objPtr); return objPtr; }     \
                                                                             \
    void insert_first(typePtr objPtr)                                        \
        { notSorted insert_link_first((void*) objPtr); }                     \
    void append(typePtr objPtr)                                              \
        { notSorted append_link((void*) objPtr); }                           \
    typePtr extract()                                                        \
        { notSorted return (typePtr) extract_link(); }                       \
    typePtr nullify()                                                        \
        { notSorted return (typePtr) nullify_link(); }                       \
    typePtr change_to(typePtr objPtr)                                        \
        { notSorted return (typePtr) change_item_to((void*) objPtr); }       \



// The standard DLList declaration macro creates a list class that is derived
// from DLList.  The class (name) stores pointers to data of type typePtr.
// No list ordering capability is provided by the standard declare macro
#define DLListdeclare(name, typePtr)                                         \
                                                                             \
class name : public DLList                                                   \
{                                                                            \
  friend class name##Iterator;                                               \
                                                                             \
  public:                                                                    \
                                                                             \
    name(int list_size=0) : DLList (list_size)                               \
    {                                                                        \
      assert(sizeof(typePtr) == sizeof(void*));                              \
    }                                                                        \
    void insert(typePtr objPtr) { insert_link((void*) objPtr); }             \
    typePtr remove(typePtr objPtr)                                           \
    {                                                                        \
      return (typePtr) move_to_and_remove_item((void*) objPtr);              \
    }                                                                        \
    CubitBoolean move_to(const typePtr objPtr)                               \
    {                                                                        \
      if (nullItem && (typePtr)nullItem == objPtr)                           \
	return CUBIT_FALSE;                                                \
      else                                                                   \
	return (move_to_item((void*) objPtr));                             \
    }                                                                        \
    CubitBoolean move_to_nearby(const typePtr objPtr)                        \
    {                                                                        \
      if (nullItem && (typePtr)nullItem == objPtr)                           \
	return CUBIT_FALSE;                                                \
      else                                                                   \
	return (move_to_nearby_item((void*) objPtr));                      \
    }                                                                        \
    int distance_to_nearby(const typePtr objPtr)                             \
    {                                                                        \
      if (nullItem && (typePtr)nullItem == objPtr)                           \
	return CUBIT_FALSE;                                                \
      else                                                                   \
	return (distance_to_nearby_item((void*) objPtr));                  \
    }                                                                        \
    CubitBoolean is_in_list(const typePtr objPtr) const                      \
    {                                                                        \
      return where_is_item((void*) objPtr) >= 0 ? CUBIT_TRUE : CUBIT_FALSE;  \
    }                                                                        \
    name& operator=(const name##Iterator & iter);                            \
    CommonDefine(typePtr, ;)                                                 \
};                                                                           \
                                                                             \
class name##Iterator : public DLListIterator                                 \
{                                                                            \
  friend class name;                                                         \
public:                                                                      \
  name##Iterator() : DLListIterator() { dlList = NULL; }                     \
  void watch( const name *dl_list )                                          \
    { dlList = dl_list; }                                                    \
  int size() const { return dlList ? dlList->size() : 0; }                   \
  typePtr get() const                                                        \
    {return (typePtr) get_item();}                                           \
  typePtr get_and_step( int n = 1 )                                          \
    {                                                                        \
      typePtr temp = (typePtr) get_item();                                   \
      step( n );                                                             \
      return temp;                                                           \
    }                                                                        \
  typePtr next( int n = 1 ) const                                            \
    {return (typePtr) next_item(n);}                                         \
  CubitBoolean is_in_list( const typePtr objPtr ) const                      \
    { return where_is_item( (void*) objPtr ) >= 0 ? CUBIT_TRUE :             \
      CUBIT_FALSE; }                                                         \
  CubitBoolean move_to(const typePtr objPtr )                                \
    { index = where_is_item( (void*) objPtr );                               \
      return index > 0 ? CUBIT_TRUE : CUBIT_FALSE; }                         \
private:                                                                     \
  const name *dlList;                                                        \
  inline void *next_item( int n = 1 ) const;                                 \
  void *get_item() const  { return next_item(0); }                           \
  int where_is_item( void *item ) const                                      \
    { return dlList ? dlList->where_is_item( item ) : -1; }                  \
};                                                                           \
                                                                             \
inline void *name##Iterator::next_item( int n ) const                        \
{                                                                            \
  int size_loc = size();                                                     \
  if ( size_loc == 0 )                                                       \
    return NULL;                                                             \
  int index_loc = index + n;                                                 \
  while ( index_loc < 0 )                                                    \
    index_loc += size_loc;                                                   \
  while ( index_loc >= size_loc )                                            \
    index_loc -= size_loc;                                                   \
 return dlList->listArray[index_loc];                                        \
}                                                                            \
inline name& name::operator=(const name##Iterator& iter)                     \
{ if ( iter.dlList ) *this = *(iter.dlList);                                 \
  else clean_out();                                                          \
  return *this; }                                                            \

// Use the list iterator like so:
//
// DLListDeclare( DLCubitNodeList, CubitNode * );
// DLCubitNodeList node_list;
// DLCubitNodeListIterator node_iterator;
// node_iterator.watch( node_list );
// node_iterator.reset();
// for (int i = node_iterator.size(); i--; )
// {
//    CubitNode *node = node_iterator.get_and_step();
//    ...
// }
//

#endif //- DLLIST_HPP


