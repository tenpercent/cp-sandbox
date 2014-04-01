//- Class: DLIList
//- Owner: Darryl Melander

#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

#include <set>
#include "DLIList.hpp"

template <class X> MY_INLINE X DLIList<X>::push(X value)
{
   last();
   append(value);
   return value;
}

template <class X> MY_INLINE X DLIList<X>::pop()
{
  last();
  return remove();
}

//- Constructor: Create a list with initial size 0. The list
//- will be grown by {increment} each time it is filled. Memory for the
//- list is not allocated until the first element is inserted using
//- {insertLink}. 
template <class X> MY_INLINE DLIList<X>::DLIList (int list_size)
{
   index      = 0;
   itemCount  = 0;
   listLength = 0;
   listArray  = NULL;
   if (list_size)
   {
      listArray = new X [list_size];
      listLength = list_size;
   }
}

//- Copy Constructor
template <class X> MY_INLINE DLIList<X>::DLIList(const DLIList<X>& from)
{
   if (&from != this)
   {
        // Setup the variables
      index = from.index;
      itemCount = from.itemCount;
      if (itemCount)
      {
         listArray  = new X [from.listLength];
         listLength = from.listLength;
           // Now copy the data
         memcpy (listArray, from.listArray, listLength*sizeof(X));
           // Maybe the last one should be itemCount instead of listLength?
      }
      else
      {
         listArray = NULL;
         listLength = 0;
      }
   }
}

// Destructor
template <class X> MY_INLINE DLIList<X>::~DLIList()
{
   if (listArray != NULL)
      delete [] listArray;
}

template <class X> MY_INLINE void DLIList<X>::lengthen_list(int by_how_much,
                                                  double by_what_factor)
{
    // Make a new array
  int new_size = 
    static_cast<int>( static_cast<double>(listLength) * by_what_factor )
    + by_how_much;
  X* temp_list = new X [new_size];

    // If the list wasn't empty, copy over the old stuff
  if (listLength)
  {
    memcpy (temp_list, listArray, listLength*sizeof(X));
    delete [] listArray;
  }
  listLength = new_size;
  listArray = temp_list;
}

//- put new item in list after current item and make it current
template <class X> MY_INLINE void DLIList<X>::insert(X new_item)
{
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
      listArray[i] = listArray[i-1];

   listArray[index] = new_item;
   itemCount++;
}

//- merge the input list, merge_list, with the current list, making
//- sure not to add duplicate items into the current list
template <class X> MY_INLINE 
void DLIList<X>::merge_unique ( const DLIList<X>& merge_list, 
                                bool  merge_list_unique )
{
     // MJP Note:
     // This procedure could be much more efficient if sorted lists
     // are used. However, I need this procedure at this time to merge
     // DLILists that already exist. These were not created as sorted
     // lists (SDLILists) and it would be painful to convert them to
     // SDLILists. It would be a lot easier if one could simply sort 
     // a DLIList based on the numeric values of its items.

     // Save the current index of the merge_list
   int old_size = size();   
   int i, j, check_index;

   X new_item;
   
   for ( i = 0; i < merge_list.size(); i++)
   {
        // Get the item from the merge_list and insert it into "this"
        // list if it doesn't already exist there.
      new_item = merge_list.next(i);
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
        append(new_item);
//      append_unique(new_item);
   }
}

template <class X> MY_INLINE void DLIList<X>::intersect_unordered( 
  const DLIList<X>& merge_list )
{
  if ( &merge_list == this )
     return;

  DLIList <X> intersect_list(merge_list);   // copy input list so can sort
  sort();
  intersect_list.sort();
  
  X* iter1 = listArray;                      // iterator for this array
  X* end1 = iter1 + itemCount;               // end of this array
  X* iter2 = intersect_list.listArray;       // iterstor for other array
  X* end2 = iter2 + intersect_list.itemCount;// end of other array
  X* last_insert = iter1 - 1;                // location of last insert
  
  for ( ; iter1 < end1; ++iter1 )
  {
    while (iter2 < end2 && *iter2 < *iter1)
      ++iter2;
    
    if (iter2 == end2)
      break;
    
    if ((*iter2 == *iter1) &&        // items are the same and ...
        (last_insert < listArray ||  // is the first item or ...
         *iter1 != *last_insert))    // is not the same as the previous item
      *++last_insert = *iter1;
  }
  
  itemCount = last_insert - listArray + 1;
  reset();
}    

template <class X> MY_INLINE void DLIList<X>::intersect ( const DLIList<X>& merge_list )
{
  if ( &merge_list == this )
     return;

  CubitBoolean removed_something = CUBIT_FALSE;

  X new_item = NULL;
  for ( int i = size(); i--; )
  {
    new_item = get();
    if ( !(merge_list.is_in_list(new_item)) )
    {
      change_to(NULL);
      removed_something = CUBIT_TRUE;
    }
    step();
  }

  if ( removed_something )
     remove_all_with_value(NULL);
}

template <class X> MY_INLINE void DLIList<X>::intersect ( void* merge_list )
{
  intersect( *static_cast< DLIList<X>* >(merge_list) );
}


//- remove the item at the current location and return a pointer to it.
//- The next node becomes the current node
//- Returns 0 if there are no items in list
template <class X> MY_INLINE X DLIList<X>::remove ()
{
   if ( !itemCount )
   {
      PRINT_WARNING("Attempted link removal from empty DLIList\n");
      return X(0);
   }

     // save the current value
   X temp = listArray[index];

     // compress memory
   for ( int i = index; i < itemCount-1; i++ )
   {
      listArray[i] = listArray[i+1];
   }

   itemCount--;
   if ( index == itemCount )
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
//- Returns 0 if there are no items in list
template <class X> MY_INLINE X DLIList<X>::extract ()
{
   if ( !itemCount )
   {
      PRINT_WARNING("Attempted link removal from empty DLIList\n");
      return X(0);
   }

     // save the current value
   X temp = listArray[index];

     // assign last node to the current index
   listArray[index] = listArray[itemCount - 1];

     // decrement itemCount
   itemCount--;
   if ( index == itemCount && index )
        // The choices here are index at beginning or end.
        // End seems to work better when 'extract' is being used
        // with calls to 'move_to_item'. 
      index--;

   return temp;
}

//+//Added so list removals don't disturb current position. PRK 05-23-94
//+//Corrected for omitting the last item in the list. PRK 09-16-94
//+//Corrected for omitting before the current position. PRK 10-07-94
//- Finds instance of item by matching value and delets first instance
//- of it from the list. The current position of the list is not changed.
//- Returns CUBIT_TRUE if the item was found and deleted, CUBIT_FALSE if not.
template <class X> MY_INLINE CubitBoolean DLIList<X>::omit(X old_val)
{
   int scan_index;
   int squeeze_index = 0;
   CubitBoolean found = CUBIT_FALSE;
   for(scan_index = 0; scan_index < itemCount; ++scan_index)
   {
      if(listArray[scan_index] == old_val)
      {
         found = CUBIT_TRUE;
         if(index == scan_index)
            index = squeeze_index - 1;
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
      if(index < 0)
         index = (itemCount) ? (itemCount - 1) : 0 ;
   }

   return found;
}

template <class X> MY_INLINE X DLIList<X>::step_and_get () 
{
   if ( !itemCount )
   {
      PRINT_WARNING("Attempted step_and_get from empty DLIList\n");
      return X(0);
   }

   if (++index == itemCount)
      index=0;
   X temp = listArray[index];
   return temp;
}

template <class X> MY_INLINE DLIList<X>& DLIList<X>::
                       operator=(const DLIList<X>& from)
{
   if (this != &from)
   {
      index = from.index;
      itemCount = from.itemCount;
      if (listArray)
         delete [] listArray;
      if (itemCount)
      {
         listArray  = new X [from.listLength];
         listLength = from.listLength;
           // Now copy the data
         memcpy (listArray, from.listArray, listLength*sizeof(X));
           // Maybe the last one should be itemCount instead of listLength?
      }
      else
      {
         listArray = NULL;
         listLength = 0;
      }
   }
   return *this;
}

template <class X> MY_INLINE DLIList<X>& DLIList<X>::
                       operator+=(const DLIList<X>& from)
{
     // Don't do anything if the list being appended is empty.
   if (from.itemCount == 0)
      return *this;

     // Make sure the array is big enough
   int tmp_itemCount = itemCount + from.itemCount;
   if (tmp_itemCount >= listLength)
     lengthen_list(tmp_itemCount - listLength, 2.0 ); 
     // factor of 1.0 can cause huge inefficiencies

     // Now add the 'from' items to the list
   if (from.itemCount == 1)
   {
      listArray[itemCount] = from.listArray[0];
   }
   else
   {
      memcpy (&listArray[itemCount], from.listArray,
              from.itemCount*sizeof(X));
   }

     // Increase the itemCount
   itemCount += from.itemCount;
   return *this;
}
template <class X> MY_INLINE DLIList<X>& DLIList<X>::
                       operator+=(void* from)
{
  DLIList<X>* cast_from = static_cast< DLIList<X>* >(from);
  
     // Don't do anything if the list being appended is empty.
   if (cast_from->itemCount == 0)
      return *this;

     // Make sure the array is big enough
   int tmp_itemCount = itemCount + cast_from->itemCount;
   if (tmp_itemCount >= listLength)
      lengthen_list(tmp_itemCount - listLength, 2.0 );
     // factor of 1.0 can cause huge inefficiencies

     // Now add the 'from' items to the list
   if (cast_from->itemCount == 1)
   {
      listArray[itemCount] = cast_from->listArray[0];
   }
   else
   {
      memcpy (&listArray[itemCount], cast_from->listArray,
              cast_from->itemCount*sizeof(X));
   }

     // Increase the itemCount
   itemCount += cast_from->itemCount;
   return *this;
}
template <class X> MY_INLINE DLIList<X>& DLIList<X>::
                       operator-=(const DLIList<X>& from)
{
    // step through items in from list, removing them from this list.
   X val;
   for (int i = from.itemCount; i--; )
   {
     // quit early if this list is empty.
     if (itemCount == 0)
       break;
     val = from.listArray[i];
     remove_all_with_value(val);
   }   
   return *this;
}

template <class X> MY_INLINE int DLIList<X>::operator==(const DLIList<X> &from)
{
  if(itemCount != from.itemCount)
     return CUBIT_FALSE;
  DLIList<X> temp_list = from;
  for( int i = 0; i < itemCount; i++)
     if(temp_list.move_to(listArray[i]))
        temp_list.remove();
  return temp_list.itemCount == 0;
}
template <class X> MY_INLINE int DLIList<X>::operator!=(const DLIList<X> &from)
{  
  return !( *this == from );
}

// Sorts list from low to high value, according to the > operator
// of the list type.
// List is sorted using a standard Heap Sort algorithm ("Algorithms in C++",
// Robert Sedgewick).
template <class X> MY_INLINE void DLIList<X>::sort(typename DLIList<X>::SortFunction f)
{
    // Only sort it if there is more than one
    // item in the list
  if (itemCount > 1)
  {
    int mid = (itemCount >> 1) + 1;
    int ir = itemCount;
    X temp_element;
    
      // You loop until ir is 1 (ir = iterations remaining)
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
          if (f(listArray[j], listArray[j - 1]) >= 0)
            j++;
        }
        if (f(listArray[j - 1], temp_element) >= 0)
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

template <class X> MY_INLINE void DLIList<X>::sort()
{
    // Only sort it if there is more than one
    // item in the list
  if (itemCount > 1)
  {
    int mid = (itemCount >> 1) + 1;
    int ir = itemCount;
    X temp_element;
    
      // You loop until ir is 1 (ir = iterations remaining)
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
          if (listArray[j] >= listArray[j - 1])
            j++;
        }
        if (listArray[j - 1] >= temp_element)
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

template <class X> MY_INLINE void DLIList<X>::reverse()
{
  int front = 0; 
  int tail  = itemCount-1;
  X temp;

  while (front < tail)
  {
     temp             = listArray[front];
     listArray[front] = listArray[tail];
     listArray[tail]  = temp;
     tail--;
     front++;
  }
}

template <class X> MY_INLINE int DLIList<X>::memory_use(CubitBoolean verbose_boolean)
{
   // report amount of memory allocated

   int my_size = listLength * sizeof(X);

   if (verbose_boolean)
   {
      PRINT_INFO("      DLIList: %d bytes\n", my_size);
   }

   return my_size;
}

template <class X> MY_INLINE void DLIList<X>::copy_to(X *other_array)
    //- copy this list's listArray into other_array
{
  assert(other_array != 0);
  if (itemCount)
    memcpy (other_array, listArray, itemCount*sizeof(X));
}
  
template <class X> MY_INLINE void DLIList<X>::copy_from(X *other_array, const int other_size)
  //- copy other_array into listArray
{
  assert(other_array);
  if (other_size > listLength) {
    lengthen_list(other_size - listLength, 1.0);
  }
  
  memcpy (listArray, other_array, other_size*sizeof(X));
  itemCount = other_size;
  
}

template <class X> MY_INLINE CubitBoolean DLIList<X>::move_to_nearby(const X item)
{
//  if (nullItem && (X)nullItem == item)
//     return CUBIT_FALSE;
//  else 
  if (itemCount == 0)
     return CUBIT_FALSE;
  else
     
  {
    X* ptr_up = &(listArray[index]);
    if (*ptr_up == item)
       return CUBIT_TRUE;

    int i_up, i_down;
    i_up = i_down = index;
    X *ptr_down = ptr_up;
    while (1)
    {
        // check forward in the list increment
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
        //  increment
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
}

template <class X> MY_INLINE int DLIList<X>::distance_to_nearby(const X body)
{
//  if (nullItem && (X)nullItem == body)
//     return CUBIT_FALSE;
//  else
  {
    int old_index = index;
    move_to_nearby( body );
    int distance = abs(index - old_index);
    if (distance > itemCount / 2)
       distance = itemCount - distance;
    index= old_index;
    return distance;
  }
}

template <class X> MY_INLINE CubitBoolean DLIList<X>::move_between(X item1, X item2)
{
//  if (nullItem && (((X)nullItem == item1 ||
//                    ((X)nullItem == item2))))
//     return CUBIT_FALSE;
//  else
  {
    assert(item1 != NULL && item2 != NULL);
    if ( !itemCount )
       return CUBIT_FALSE;

    for ( int i = 0; i < itemCount; i++ )
    {
      if ( listArray[i] == item1 )
      {
          //  dance around looking for item2
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
}

// Removes every instance of 'val' in list.
// Set to beginning of list after this call.
template <class X> MY_INLINE void DLIList<X>::uniquify_unordered()
{
  if ( itemCount < 2 )
    return;
  
  sort();
  
  int j = 1;
  int i = 0;
  
  while (j < itemCount)
  {
    if (listArray[i] != listArray[j])
      listArray[++i] = listArray[j++];
    else
      j++;
  }
  
  itemCount = i + 1;
  index = 0;
}

template <class X> MY_INLINE void DLIList<X>::uniquify_ordered()
{
  std::set<X> encountered;
  X *read, *write, *const end = listArray + itemCount;
  
    // To avoid copying the array onto itself in the case
    // where the list contains no duplicates, loop twice.
    
    // Find first duplicate entry.
  for ( read = listArray; read != end; ++read )
    if ( ! encountered.insert(*read).second )
      break;

    // Now compact array, removing duplicates.  If there
    // are no duplicates, this loop will not run (read == end).
  for ( write = read; read != end; ++read )
    if ( encountered.insert(*read).second )
      *(write++) = *read;
  
  itemCount = write - listArray;
  index = 0;
}


  
