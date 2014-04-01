//- Class: DLList
//- Owner: Greg Sjaardema
//- Checked by:
//- Version: $Id: 

#include "DLList.hpp"


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
DLList::DLList ( int increment ): ArrayBasedContainer( increment )
{
  index      = 0;
  nullItem   = 0;
}

//- Copy Constructor
DLList::DLList(const DLList& from) : ArrayBasedContainer ( from )
{
  index = from.index;
  listIsSorted = CUBIT_FALSE;
}

DLList::~DLList()
{
}

//- put new item in list after current item and make it current
void DLList::insert_link ( void* new_item )
{
  assert(new_item != NULL);
  // see if the list must be lengthened
  if ( itemCount == listLength )
  {
      lengthen_list();
  }
  
  // set new index
  
  if ( itemCount )
  {
      index++;
  }
  else
  {
      index = 0;
  }
  
  // the item must be put in the current index, all higher
  // indexes must be bubbled up one spot
  
  for ( int i = itemCount; i > index; i-- )
  {
      listArray[i] = listArray[i-1];
  }
  
  listArray[index] = new_item;
  itemCount++;
}

//- merge the input list, merge_list, with the current list, making
//- sure not to add duplicate items into the current list
void DLList::merge_unique ( DLList& merge_list, int merge_list_unique )
{
   // MJP Note:
   // This procedure could be much more efficient if sorted lists
   // are used. However, I need this procedure at this time to merge
   // DLLists that already exist. These were not created as sorted
   // lists (SDLLists) and it would be painful to convert them to
   // SDLLists. It would be a lot easier if one could simply sort 
   // a DLList based on the numeric values of its items.
   
   // Save the current index of the merge_list
   int current_index = merge_list.index;
   int old_size = size();   
   int i, j, check_index;
   void* new_item = NULL;   

   for (i = 0; i < merge_list.size(); i++)
   {
      // Get the item from the merge_list and insert it into "this"
      // list if it doesn't already exist there.
      new_item = merge_list.get_item_and_step();
      check_index = merge_list_unique ? old_size : size();
      
      for ( j = 0; j < check_index; j++ )
      {
        if ( listArray[j] == new_item )
        {
          check_index = -1;
          break;
        }
      }
      
      if ( check_index != -1 )
        append_link(new_item);
   }
   
   // Restore the original index of the merge_list
   merge_list.index = current_index;
}


//- put new item in list at index 0 and make it current
void DLList::insert_link_first(void* new_item)
{
  // set index to -1 ... index will be set to 0 in insert_link

  index = -1;
  insert_link(new_item);
}

//- remove the item at the current location and return a pointer to it.
//- The next node becomes the current node
//- Returns {NULL} if there are no items in list
void* DLList::cut_link ()
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted link removal from empty DLList\n");
        return NULL;
    }

    // save the current value

    void *temp = listArray[index];

    // compress memory

    for ( int i = index; i < itemCount-1; i++ )
    {
        listArray[i] = listArray[i+1];
    }

    itemCount--;
    if ( index >= itemCount )
    {
        index = 0;
    }

    return temp;
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
void* DLList::extract_link ()
{
    if ( !itemCount ) {
      PRINT_WARNING("Attempted link removal from empty DLList\n");
      return NULL;
    }
    
    // save the current value
    void *temp = listArray[index];
    
    // assign last node to the current index
    listArray[index] = listArray[itemCount - 1];
    
    // decrement itemCount
    itemCount--;
    if ( index == itemCount && index != 0) {
      // The choices here are index at beginning or end.
      // End seems to work better when 'extract' is being used
      // with calls to 'move_to_item'. 
      index--;
    }
    
    return temp;
}


//+//Added so list removals don't disturb current position. PRK 05-23-94
//+//Corrected for omitting the last item in the list. PRK 09-16-94
//+//Corrected for omitting before the current position. PRK 10-07-94
//- Finds instance of item by matching pointer and deleting all instances
//- of it from the list. The current position of the list is not changed.
//- Returns CUBIT_TRUE if the item was found and deleted, CUBIT_FALSE if not.
CubitBoolean DLList::omit_link(void *oldObjPtr)
{
    int scan_index;
    int squeeze_index = 0;
    CubitBoolean found = CUBIT_FALSE;
    for(scan_index = 0; scan_index < itemCount; ++scan_index)
    {
        if(listArray[scan_index] == oldObjPtr)
        {
            found = CUBIT_TRUE;
            if(index == scan_index) index = squeeze_index - 1;
        }
        else 
        {
            if(scan_index != squeeze_index)
            {
                listArray[squeeze_index] = listArray[scan_index];
                if(index == scan_index) index = squeeze_index;
            }
            ++squeeze_index;
        }
    }

    if(found)
    {
        itemCount = squeeze_index;
//+//   If the first item was deleted and index pointed to it, make an
//+//   adjustment here. If itemCount is zero, don't assign -1 again.
        if (index < 0)
	  index = (itemCount) ? (itemCount - 1) : 0 ;
	if (index >= itemCount)
	  index = (itemCount) ? (itemCount - 1) : 0 ;
    }

    return found;
}


//- Change the current item to a null pointer and return a pointer
//- to the item (before nulled). See the discussion for {nullItem}.
void* DLList::nullify_link ()
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted link nullify from empty DLList\n");
        return NULL;
    }

    // save the current value
      nullItem = listArray[index];
      listArray[index] = 0;
    return nullItem;
}

void *DLList::prev_item ( int n )  const
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted prev of empty DLList\n");
        return NULL;
    }

    // return the proper index
    // beware of negative n
    int new_index = index - n;

    while (new_index < 0)
      new_index += itemCount;
    // beware of negative n leading to new_index >itemCount
    new_index = new_index%itemCount;

    assert(listArray[new_index] != NULL);
    return listArray[new_index];
}

void *DLList::step_and_get_item () 
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted step_and_get from empty DLList\n");
        return NULL;
    }

    if (++index == itemCount)
      index=0;
    void *temp = listArray[index];
    assert(temp != NULL);
    return temp;
}

void *DLList::get_item_and_back ()
{
    if ( !itemCount )
    {
        PRINT_WARNING("Attempted get_and_back from empty DLList\n");
        return NULL;
    }
    void *temp = listArray[index--];
    if (index < 0) index=itemCount-1;
    assert(temp != NULL);
    return temp;
}

// move_to_and_remove_item searches for item and removes it from the list.
// If the item was found and removed it is returned, else NULL is returned.
void* DLList::move_to_and_remove_item(void* item)
{
  if (move_to_item(item))
    return cut_link();
  else
    return (void*) NULL;
}

int DLList::distance_to_nearby_item(void* body)
{
  int old_index = index;
  move_to_nearby_item( body );
  int distance = abs(index - old_index);
  if ( distance > itemCount / 2 )
    distance = itemCount - distance;
  index = old_index;
  return distance;
}


CubitBoolean DLList::move_to_nearby_item ( void* item )
{
      // empty list
  if (itemCount == 0) 
    return CUBIT_FALSE;
  
  // Search current item first...
  void **ptr_up = &(listArray[index]);
  if (*ptr_up == item) 
    return CUBIT_TRUE;

  int i_up, i_down;
  i_up = i_down = index;
  void **ptr_down = ptr_up;
  while (1) 
  {
      // check forward in the list
      // increment
    if ( ++i_up < itemCount )
      ptr_up++;
    else
    {
      i_up = 0;
      ptr_up = listArray;
    }
      // check
    if ( *ptr_up == item )
    {
      index = i_up;
      return CUBIT_TRUE;
    }
    if ( i_up == i_down )
    {
      return CUBIT_FALSE;
    }
    
      // check backward in the list
      // increment
    if ( --i_down >= 0 )
      ptr_down--;
    else
    {
      i_down = itemCount-1;
      ptr_down = &(listArray[i_down]);
    }
      // check
    if ( *ptr_down == item )
    {
      index = i_down;
      return CUBIT_TRUE;
    }
    if ( i_up == i_down )
    {
      return CUBIT_FALSE;
    }
    
  }
}

CubitBoolean DLList::move_to_item ( void* item )
{
  int item_position = where_is_item( item );
  if ( item_position >= 0 ) {
    index = item_position;
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

int DLList::where_is_item( void *item ) const
{
  // test for null input item
  assert(item != NULL);

  if (itemCount == 0) return -1;
  
  // loop through list searching for item ...
  // if found set index and return true

  // Search current item first...
  void **ptr = &(listArray[index]);
  if (*ptr == item) return index;

  // Now, search from current index to end of array
  int i;
  for (i=index+1; i < itemCount; i++) {
    ptr++;
    if ( *ptr == item) {
        // check if move_to_nearby would have been better
        //if ( i - index > 1000 )
        // PRINT_INFO(" Found, i = %d, index = %d, itemCount = %d.\n",
        //           i, index, itemCount);
      return i;
    }
  }
  
  // Now search from beginning of array to index...
  ptr = listArray;
  for (i = 0; i < index; i++) {
    if (*ptr == item) {
        // check if move_to_nearby would have been better
        // if ( i + itemCount - index > 1000 )
        // PRINT_INFO(" Found, i = %d, index = %d, itemCount = %d.\n",
        //           i, index, itemCount);
      return i;
    }
    ptr++;
  }
  
  // item is not in array, return false
  return -1;
}

// Special care for size two lists, so that wrap around is used only when
// needed.  This works faster than the original (commented out below) for
// long lists
CubitBoolean DLList::move_between_items(void *item1, void *item2)
{
  
  assert(item1 != NULL && item2 != NULL);
  if ( !itemCount )
      return CUBIT_FALSE;
  
  for ( int i = 0; i < itemCount; i++ )
    {
      if ( listArray[i] == item1 )
        {
	  //dance around looking for item2
	  if ( i+1 < itemCount ) {
	    if ( listArray[i+1] == item2 ) {
	      index = i+1;
	      return CUBIT_TRUE;
	    }
	  }
	  if ( i > 0 ) {
	    if ( listArray[i-1] == item2 ) {
	      index = i;
	      return CUBIT_TRUE;
	    }
	  }
	  if ( ( i+1 == itemCount && listArray[0] == item2 ) 
	       || ( i == 0 && listArray[itemCount-1] == item2 ) )
	    {
	      index = 0;
	      return CUBIT_TRUE;
	    }
	}
    }
  return CUBIT_FALSE;
}

 
//// following is somewhat slower for long lists, but if you wanted to search
//// for long patterns a modification of it would be good
//int DLList::move_between_items(void *item1, void *item2)
//{
//if ( !itemCount )
//return CUBIT_FALSE;
//
//int match = 0; //which item did the previous index match?
//
//// search for item1 and item2 consecutive
//for ( int i = itemCount+1; i > 0; i--, step())
//{
//    if ( listArray[index] == item1 )
//      if (match == 2)
//	return CUBIT_TRUE;
//      else 
//	match = 1;
//    else if ( listArray[index] == item2 )
//      if (match == 1)
//	return CUBIT_TRUE;
//      else 
//	match = 2;
//    else 
//      match = 0;
//} 
//
//// not found consecutive
//return CUBIT_FALSE;
//}


void* DLList::change_item_to ( void *new_item )
{
  assert(new_item != NULL);
  if ( !itemCount )
    {
      PRINT_WARNING("DLList: Attempted update of empty list\n");
      return NULL;
    }
  void *temp = listArray[index];
    listArray[index] = new_item;
  return temp;
}

DLList& DLList::operator=(const DLList& from)
{
  if (this != &from) {
    ArrayBasedContainer::operator=(from);
    index = from.index;
  }
  return *this;
}

void DLList::compress()
{ 
   int j = 0; 
   int i;
   int new_index = index;

   for ( i = 0; i < itemCount; i++ )
      if (listArray[i] != NULL)
      {
         listArray[j++] = listArray[i];
         if (i < index)
            new_index--;
      }
   
   itemCount = j; 
   index = new_index;
   if (index >= itemCount)
      index = 0;
   nullItem = 0;
}

void DLList::reverse()
{
  int front = 0; 
  int tail  = itemCount-1;
  void *temp;
  
  while (front < tail) {
    temp             = listArray[front];
    listArray[front] = listArray[tail];
    listArray[tail]  = temp;
    tail--;
    front++;
  }
}




int DLList::operator<(  const DLList& list ) const
{
	return (itemCount < list.itemCount) && (list >= *this);
}

int DLList::operator>(  const DLList& list ) const
{
	return (itemCount > list.itemCount) && (*this >= list);
}

int DLList::operator<=( const DLList& list ) const
{
	return list >= *this;
}

int DLList::operator>=( const DLList& list ) const
{
	if( itemCount < list.itemCount ) return 0;
	
	DLList temp_list = list;
	for( int i = 0; i < itemCount; i++ )
		if( temp_list.move_to_item( listArray[i] ) ) temp_list.extract_link();
	
	return temp_list.itemCount == 0;  
}

int DLList::operator==( const DLList& list ) const
{
	return (itemCount == list.itemCount) && (list >= *this);
}

int DLList::operator!=( const DLList& list ) const
{
	return (itemCount != list.itemCount) || !(list >= *this);
}
