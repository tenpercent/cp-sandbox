//- Class: SDLList
//- Description: SDLList is a doubly linked list class that accepts generic
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
//-              The macro {SDLListdeclare(name,typePtr)} is defined to
//-              create specific instances of the SDLList class.
//-              {name} is the name of the SDLList class created, and
//-              {typePtr} is the type of elements stored in the list.
//-              The macro also defines the functions {push(item)} and 
//-              {pop()} which allow {SDLLists} to perform like stacks.
//-
//- Assumptions: All data are stored contiguously in the array with 
//-              empty slots at the end of the array.
//-
//- Owner: Greg Sjaardema
//- Checked by: Jim Hipp, 3/28/94
//- Version: $Id: 

#ifndef SDLLIST_HPP
#define SDLLIST_HPP

#include "DLList.hpp"
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include <string.h>
#include <assert.h>
#include "CubitUtilConfigure.h"

const int ORDER_ASCENDING  = 0;
const int ORDER_DESCENDING = 1;
//- sort direction specifiers for the SSDLListDeclare Constructor

class CUBIT_UTIL_EXPORT SDLList : public DLList
{
public:
  
  SDLList (int size);
  //- Constructor: Create a list with initial size {size}. The list
  //- will be grown by {size} each time it is filled. Memory for the
  //- list is not allocated until the first element is inserted using
  //- {insertLink}. 
  //- If {size} is zero, the default increment ({INCREMENT}) will be used
  //- From an efficiency standpoint, it is very important that the 
  //- increment be set large enough to reduce the number of list 
  //- growths, but small enough to not waste memory.
  //- It is more efficient to sligthly overestimate the size than 
  //- to underestimate the size.
  
  SDLList(const SDLList& copy_from);  //- Copy Constructor
  SDLList(const DLList& copy_from);  //- Copy Constructor
  
  virtual ~SDLList();
  //- Destructor: Free all resources used by this list.
  
#ifdef BOYD15
  void reverse();
  //- Reverse the items in the list.
#endif
  
  SDLList& operator=(const SDLList&);
  //- Create a copy of a list.
  
  virtual void merge_unique(DLList& merge_list, int merge_list_unique );
  //- Merges the contents of the list, merge_list, with those of "this"
  //- list, ensuring that items that are being added do not already
  //- exist in "this" list. If the items are known to appear in the
  //- merge_list only once, then it can be done faster.
  
#ifdef BOYD15
  int append_unique(void* new_item);
  //- Appends the new item to the list, if it doesn't already exist
  //- in the list. In either case, index is not changed.
#endif

  SetDynamicMemoryAllocation(memoryManager)
  //- overloaded new and delete operators
  
protected:
  
  void insert_link(void* newBody);
  //- put new item in list after current item and make it current
  
  void insert_link_first(void* new_item);
  //- put new item in list at index 0 and make it current
  
  void append_link(void* newBody);  
  //- put new item at end of list and keep pointer where it is.
  //- This call should be used to add items to a list if order 
  //- of items is not important.
  
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
  
  void *change_item_to(void* newBody);
  //- Sets the item at the current index value to {newBody} and
  //- returns the old value at that location.
  //- Returns {NULL} if there are no items in list
  
  void sort_list();
  //- sort_list sorts SDLLists in an ascending or descending fashion
  //- (depending on the assignment of compare_order function pointer).  The
  //- list is sorted using a standard Heap Sort algorithm.
  
  int move_to_item_sorted(void* value);
  //- move_to_item_sorted performs a binary search on the list to locate
  //- the item (object) that has functionType functionName() = value
  //- (see SSDLListdeclare macro for description of functionType and
  //- functionName).  If a matching item is found then the current
  //- index is set to the matching item index and TRUE is returned.  If the
  //- object is not found then index is unchanged and FALSE is returned.
  
  void insert_item_sorted(void* value, void* new_item);
  //- insert_item_sorted uses compare function (compare_order) to perform a
  //- binary search to find the insert index in the list. The function then
  //- calls insert_link after setting the index to the insert position - 1.
  
  void* move_to_and_remove_item_sorted(void* value);
  //- move_to_and_remove_item_sorted uses move_to_item_sorted function to
  //- perform a binary search to locate the item's index which has
  //- functionType functionName() = value.  If found, remove_link is called
  //- to remove and return the item from the list.  Otherwise, NULL is
  //- returned.
  
  int (*compare_order_obj)(void* object_1, void* object_2);
  int (*compare_order)(void* object_1, void* object_2);
  int (*compare_equal)(void* object_1, void* object_2);
  //- compare function pointers are assigned within the macro definition
  //- SSDLListDeclare
  
  void set_sorted_flag(const int sorted_flag);
  //- set listIsSorted flag

  CubitBoolean where_is_item_sorted( void *value, int &item_index );
  //- Return True if item with given value is in the list, else false.
  //- Return in item_index the index of item of that value in list, 
  //- or -1 if an item of that value is not in the list.
  //- This information should never be passed outside the class.

private:
  
  void binary_search(int* min_location, int* max_location, void* item) const;
  //- binary_search performs a standard array based binary search (essentially
  //- bisection) to bracket/locate an insert/search position.
  
  static MemoryManager memoryManager;
  //- memory management static objects for performing block allocation
};

// *** inline function definitions *******************************************

//- Add item to end of list, do not change current index value.
//- This function used to reduce time required to do an insert 
//- followed by a move_to back to where we were.
 inline void SDLList::append_link ( void* new_item )
{
    assert(new_item != NULL);
    // see if the list must be lengthened

    if ( itemCount == listLength )
        lengthen_list();

    listArray[itemCount++] = new_item;
}

inline void SDLList::set_sorted_flag(const int sorted_flag)
{
    listIsSorted = sorted_flag;
}


// *** macro definitions *****************************************************

// Three separate derived classes are possible from SDLList all of which are
// defined by one of the following macros:
// 
//    SDLListdeclare(name, typePtr)
//    SSDLListdeclare(name, typePtr, functionName, functionType)
//    SSDLListIntrinsicdeclare(name, typePtr)
//
// The first is the standard SDLList declaration macro, the second provides
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
// standard SDLListdeclare as-well-as the sorted SSDLListdeclare macros.
// typePtr is still the type pointer of the data stored by the list.
// notSorted is simply a semi-colon (;) for the standard list (if white space
// is used for the parameter the CPP outputs a warning message that a NULL
// string is being substituted for the parameter ... the proper outcome but
// the excessive warning messages are a real nuisance) and a call to the
// set_sort_flag function (with an argument of False) for the sorted lists.
#define CommonDefine2(typePtr, notSorted)                                    \
                                                                             \
    typePtr remove()         {return (typePtr) cut_link();}                  \
    int omit(typePtr objPtr) {return omit_link((void*) objPtr);}             \
    typePtr get() const      {return (typePtr) get_item();}                  \
    typePtr next() const     {return (typePtr) next_item();}                 \
    typePtr next( int n ) const {return (typePtr) next_item(n);}             \
    typePtr prev() const     {return (typePtr) prev_item();}                 \
    typePtr prev(int n) const {return (typePtr) prev_item(n);}                \
    typePtr get_and_step()   {return (typePtr) get_item_and_step();}         \
    typePtr get_and_back()   {return (typePtr) get_item_and_back();}         \
    typePtr step_and_get()   {return (typePtr) step_and_get_item();}         \
                                                                             \
    int move_between(typePtr objPtr1, typePtr objPtr2)                       \
    {                                                                        \
        if (nullItem && (((typePtr) nullItem == objPtr1) ||                  \
	                 ((typePtr) nullItem == objPtr2) ))                  \
  	  return CUBIT_FALSE;                                                \
        else                                                                 \
	  return (move_between_items(((void*) objPtr1), ((void*) objPtr2))); \
    }                                                                        \
    typePtr pop()                {last(); return remove();}                  \
    typePtr push(typePtr objPtr)                                             \
        { last(); append(objPtr); return objPtr; }                           \
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
                                                                             \


// The macro CommonSortedDefine contains common public functions used by
// SSDLListdeclare and SSDLListIntrinsicdeclare
#define CommonSortedDefine(name, typePtr)                                    \
                                                                             \
    name(int size = 0, int direction = ORDER_ASCENDING) : SDLList(size)       \
    {                                                                        \
      assert(sizeof(typePtr) == sizeof(void*));                              \
      if (direction == ORDER_DESCENDING)                                     \
      {                                                                      \
        compare_order = &name::compare_descend;                              \
        compare_order_obj = &name::compare_descend_obj;                      \
      }                                                                      \
      else                                                                   \
      {                                                                      \
        compare_order = &name::compare_ascend;                               \
        compare_order_obj = &name::compare_ascend_obj;                       \
      }                                                                      \
      compare_equal = &name::compare_equate;                                 \
      set_sorted_flag(CUBIT_TRUE);                                           \
    }                                                                        \
                                                                             \
    void sort()                      { sort_list(); }                        \
    void  set_list_order_ascend()                                            \
    {                                                                        \
      compare_order_obj = &name::compare_ascend_obj;                         \
      compare_order = &name::compare_ascend; sort_list();                    \
    }                                                                        \
    void  set_list_order_descend()                                           \
    {                                                                        \
      compare_order_obj = &name::compare_descend_obj;                        \
      compare_order = &name::compare_descend; sort_list();                   \
    }                                                                        \

// The sorted SDLList declaration macro creates an ordered list class that
// is derived from SDLList.  The class (name) stores pointers to data of
// type typePtr.  List ordering is provided by sorting the objects based
// on the returned value of one of the objects functions (functionName()).
// The return type of functionName is functionType.  Searching functionality
// is provided (with move_to) to locate the object (typePtr objPtr) or a 
// value (functionType value) of the object utilizing an binary search.
#define SDLListdeclare(name, typePtr, functionName, functionType)            \
                                                                             \
class name : public SDLList                                                   \
{                                                                            \
                                                                             \
  public:                                                                    \
                                                                             \
    CommonSortedDefine(name, typePtr)                                        \
    CommonDefine(typePtr, set_sorted_flag(CUBIT_FALSE);)                     \
                                                                             \
    void insert(typePtr objPtr)                                              \
    {                                                                        \
      functionType value = objPtr->functionName();                           \
      insert_item_sorted((void*) &value, (void*) objPtr);                    \
    }                                                                        \
    typePtr remove(typePtr objPtr)                                           \
    {                                                                        \
      assert(objPtr != NULL);                                                \
      functionType value = objPtr->functionName();                           \
      return (typePtr) move_to_and_remove_item_sorted((void*) &value);       \
    }                                                                        \
    typePtr remove(functionType value)                                       \
    {                                                                        \
      return (typePtr) move_to_and_remove_item_sorted((void*) &value);       \
    }                                                                        \
    int move_to(typePtr objPtr)                                              \
    {                                                                        \
      if (nullItem && ((typePtr) nullItem == objPtr))                        \
        return CUBIT_FALSE;                                                  \
      else                                                                   \
      {                                                                      \
        assert(objPtr != NULL);                                              \
        functionType value = objPtr->functionName();                         \
        return move_to_item_sorted((void*) &value);                          \
      }                                                                      \
    }                                                                        \
    int move_to(functionType value)                                          \
    {                                                                        \
      return move_to_item_sorted((void*) &value);                            \
    }                                                                        \
    int is_in_list(functionType value)                                       \
    {                                                                        \
      int item_index;                                                        \
      int return_value =                                                     \
        where_is_item_sorted((void*) &value, item_index) ?                   \
        item_index + 1 : 0;                                                  \
      return return_value;                                                   \
    }                                                                        \
    int is_in_list(typePtr objPtr)                                           \
    {                                                                        \
      if (nullItem && ((typePtr) nullItem == objPtr))                        \
        return CUBIT_FALSE;                                                  \
      else                                                                   \
      {                                                                      \
        assert(objPtr != NULL);                                              \
        functionType value = objPtr->functionName();                         \
        return is_in_list( value );                                          \
      }                                                                      \
    }                                                                        \
                                                                             \
  private:                                                                   \
                                                                             \
    static int compare_ascend_obj(void* object_1, void* object_2)            \
    {                                                                        \
      return (((typePtr) object_1)->functionName() >=                        \
	      ((typePtr) object_2)->functionName());                       \
    }                                                                        \
    static int compare_descend_obj(void* object_1, void* object_2)           \
    {                                                                        \
      return (((typePtr) object_1)->functionName() <=                        \
    	      ((typePtr) object_2)->functionName());                       \
    }                                                                        \
                                                                             \
    static int compare_ascend(void* value, void* object)                     \
    {                                                                        \
      return (*((functionType*) value) >=                                    \
	      ((typePtr) object)->functionName());                         \
    }                                                                        \
    static int compare_descend(void* value, void* object)                    \
    {                                                                        \
      return (*((functionType*) value) <=                                    \
    	      ((typePtr) object)->functionName());                         \
    }                                                                        \
    static int compare_equate(void* value, void* object)                     \
    {                                                                        \
      return (*((functionType*) value) ==                                    \
	      ((typePtr) object)->functionName());                         \
    }                                                                        \
};                                                                           \


// The sorted Intrinsic SDLList declaration macro creates an ordered intrinsic
// list class that is derived from SDLList.  The class (name) stores values of
// typePtr through casting to void* pointers.  This macro only works for lists
// of integers or floats.  No capability is currently supported for lists of
// doubles.
#define SSDLListIntrinsicdeclare(name, typePtr)                              \
                                                                             \
class name : public SDLList                                                  \
{                                                                            \
                                                                             \
  public:                                                                    \
                                                                             \
    CommonSortedDefine(name, typePtr)                                        \
    CommonDefine(typePtr, set_sorted_flag(CUBIT_FALSE);)                     \
                                                                             \
    void insert(typePtr objPtr)                                              \
    {                                                                        \
      insert_item_sorted((void*) objPtr, (void*) objPtr);                    \
    }                                                                        \
    typePtr remove(typePtr objPtr)                                           \
    {                                                                        \
      return (typePtr) move_to_and_remove_item_sorted((void*) objPtr);       \
    }                                                                        \
    int move_to(typePtr objPtr)                                              \
    {                                                                        \
      return move_to_item_sorted((void*) objPtr);                            \
    }                                                                        \
    int is_in_list(typePtr objPtr)                                           \
    {                                                                        \
      int item_index;                                                        \
      int return_value =                                                     \
        where_is_item_sorted((void*) &objPtr, item_index) ?                  \
        item_index + 1 : 0;                                                  \
      return return_value;                                                   \
    }                                                                        \
                                                                             \
  private:                                                                   \
                                                                             \
    static int compare_ascend_obj(void* object_1, void* object_2)            \
    {                                                                        \
      return (((typePtr) object_1) >= ((typePtr) object_2));                 \
    }                                                                        \
    static int compare_descend_obj(void* object_1, void* object_2)           \
    {                                                                        \
      return (((typePtr) object_1) <= ((typePtr) object_2));                 \
    }                                                                        \
                                                                             \
    static int compare_ascend(void* object_1, void* object_2)                \
    {                                                                        \
      return (((typePtr) object_1) >= ((typePtr) object_2));                 \
    }                                                                        \
    static int compare_descend(void* object_1, void* object_2)               \
    {                                                                        \
      return (((typePtr) object_1) <= ((typePtr) object_2));                 \
    }                                                                        \
    static int compare_equate(void* object_1, void* object_2)                \
    {                                                                        \
      return (((typePtr) object_1) == ((typePtr) object_2));                 \
    }                                                                        \
};                                                                           \

#endif //- SDLLIST_HPP

