//- Class: SDLList
//- Owner: Greg Sjaardema
//- Checked by:
//- Version: $Id: 

#include "SDLList.hpp"


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
SDLList::SDLList ( int increment ): DLList( increment )
{}

//- Copy Constructor
SDLList::SDLList(const SDLList& from) : DLList(from)
{
  listIsSorted = from.listIsSorted;
}

//- Copy Constructor
SDLList::SDLList(const DLList& from) : DLList(from)
{
  listIsSorted = CUBIT_FALSE;
}

SDLList::~SDLList()
{
}

//- put new item in list after current item and make it current
void SDLList::insert_link ( void* new_item )
{
  DLList::insert_link(new_item);
  listIsSorted = CUBIT_FALSE;
}

//- merge the input list, merge_list, with the current list, making
//- sure not to add duplicate items into the current list
void SDLList::merge_unique ( DLList& merge_list, int merge_list_unique )
{
  DLList::merge_unique(merge_list, merge_list_unique);
  listIsSorted = CUBIT_FALSE;
}

//- put new item in list at index 0 and make it current
void SDLList::insert_link_first(void* new_item)
{
  DLList::insert_link_first(new_item);
  listIsSorted = CUBIT_FALSE;
}

//- remove the item at the current location and return a pointer to it.
//- used for efficiency in cases where preservation of list order is not
//- important.  moves last list item (itemCount - 1) to current index and
//- decrements itemCount.  eliminates the need to perform the list bubble
//- down (i.e. cut_link) but sacrifices list order in the process.  this
//- function should not be used when up-stream order from the removed node is
//- important.  when processing a list using this function the user should
//- reset the list to the head (index = 0) before beginning to ensure all
//- list nodes are processed properly.
//- Returns {NULL} if there are no items in list
void* SDLList::extract_link ()
{
  listIsSorted = CUBIT_FALSE;
  return DLList::extract_link();
}


void* SDLList::change_item_to ( void *new_item )
{
  listIsSorted = CUBIT_FALSE;
  return DLList::change_item_to(new_item);
}

// sort_list sorts ordered lists in an ascending or descending fashion
// (depending on the assignment of compare_order function pointer).  The
// list is sorted using a standard Heap Sort algorithm ("Algorithms in C++",
// Robert Sedgewick).
void SDLList::sort_list()
{
  listIsSorted = CUBIT_TRUE;
  if (itemCount > 1)
  {
    int mid = (itemCount >> 1) + 1;
    int ir = itemCount;
    void* temp_element = NULL;

    while(CUBIT_TRUE)
    {
      if (mid > 1)
      {
	mid--;
	temp_element = listArray[mid - 1];
      }
      else
      {
	ir--;
	temp_element = listArray[ir];
	listArray[ir] = listArray[0];
	if (ir == 1)
	{
	  listArray[0] = temp_element;
	  return;
	}
      }

      int i = mid;
      int j = mid + mid;

      while (j <= ir)
      {
	if (j < ir)
	{
	   if ((*compare_order_obj)(listArray[j], listArray[j - 1])) j++;
	}
	if ((*compare_order_obj)(listArray[j - 1], temp_element))
	{
	  listArray[i - 1] = listArray[j - 1];
	  i = j;
	  j += j;
	}
	else
	{
	  j = ir + 1;
	}
      }

      listArray[i - 1] = temp_element;
    }
  }
}

// binary_search performs a standard array based binary search (essentially
// bisection) to bracket/locate an insert/search position.
void SDLList::binary_search(int* min_location, int* max_location, void* item) 
  const
{
  // sort the list if necessary and initialize binary search parameters
  assert( listIsSorted );

  *min_location = 0;
  *max_location = itemCount - 1;
  int bracket_width = *max_location;
  int compare_location = itemCount >> 1;

    // bracket the location (binary search)

  while(bracket_width > 1)
  {
    if ((*compare_order)(item, listArray[compare_location]))
    {
      *min_location     = compare_location;
      bracket_width     = *max_location - *min_location;
      compare_location += (bracket_width >> 1);
    }
    else
    {
      *max_location     = compare_location;
      bracket_width     = *max_location - *min_location;
      compare_location -= (bracket_width >> 1);
    }
  }

    // if there are duplicate entries in the list, bracket them with min_location
    // and max_location
  while (  *min_location > 0 && (*compare_equal)(item, listArray[*min_location-1]))
    (*min_location)--;

    // (max_location already backets the item)
  
}

// move_to_item_sorted finds the index of the given value, then sets the 
// current index there
int SDLList::move_to_item_sorted(void* value)
{
  int item_index;
  CubitBoolean item_exists = where_is_item_sorted( value, item_index );
  index = item_index; // always
  if ( item_exists ) {
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}


// where_is_item_sorted performs a binary search on the list to locate value.
// If the item with value is found then the current index is set to the item's
// index and TRUE is returned.  If item is not found then index is set to 
// just before where the item would be, and FALSE is returned.
CubitBoolean SDLList::where_is_item_sorted( void *value, int &insert_index ) 
{
  // if entries exist then find insertion index

  if ( 0 == itemCount ) {
    insert_index = 0;
    return CUBIT_FALSE;
  }
  
  if (!listIsSorted)
  {
    sort_list();
  }

  // perform the binary search
  
  int min_location, max_location;
  binary_search(&min_location, &max_location, value);
  
  // if item occupies a place in the list then min_location and/or
  // max_location specifies the location ... otherwise item was not found
  
  if ((*compare_equal)(value, listArray[min_location]))
    {
      insert_index = min_location;
      return CUBIT_TRUE;
    }
  else if ((*compare_equal)(value, listArray[max_location]))
    {
      insert_index = max_location;
      return CUBIT_TRUE;
    }  
  // else not found, 
  insert_index = min_location;
  // is it beyond the end?
  if (itemCount == max_location ) {
    if ((*compare_order)(value, listArray[max_location]))
      insert_index = itemCount;
  }
  return CUBIT_FALSE;
}


// insert_item_sorted performs a binary search to find the insert
// index in the list for new_item. The function then calls insert_link
// after setting the index to the insert position - 1.
void SDLList::insert_item_sorted(void* value, void* new_item)
{
  assert(new_item != NULL);

  // if entries exist then find insertion index

  if ( itemCount > 0)
  {
    // perform the binary search and set index to the min_location

    int min_location, max_location;
    if (!listIsSorted)
    {
      sort_list();
    }
    binary_search(&min_location, &max_location, value);
    index = min_location;

    // test for special cases (insert first or last)

    if (index == 0)
    {
      // if new_items value < first (ascending) or new_items value > first
      // (descending) then insert as first item (index = -1)

      if (!(*compare_order)(value, listArray[index]))
      {
        index--;
      }
    }
    if (index == (itemCount - 2))
    {
      // if new_items value > last (ascending) or new_items value < last
      // (descending) then insert as last item (index = itemCount - 1)

      if ((*compare_order)(value, listArray[index + 1]))
      {
        index++;
      }
    }
  }

  // insert new_item
  insert_link(new_item);
  listIsSorted = CUBIT_TRUE;
}

// move_to_and_remove_item_sorted uses move_to_item_sorted function to perform
// a binary search to locate the objects's index for which the objects
// functionName() = value.  If found, the objects index is set as the current
// index. Then cut_link is called to remove the object from the list.
void* SDLList::move_to_and_remove_item_sorted(void* value)
{
  if (move_to_item_sorted(value))
    return cut_link();
  else
    return (void*) NULL;
}


SDLList& SDLList::operator=(const SDLList& from)
{
  if (this != &from) {
    DLList::operator=(from);
  }
  return *this;
}

