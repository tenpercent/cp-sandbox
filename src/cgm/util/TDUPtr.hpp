//-------------------------------------------------------------------------
// Filename      : TDUPtr.hpp
//
// Purpose       : Smart pointer that becomes NULL when pointed-to object
//                 is destroted.
//
// Special Notes : Uses ToolData mechanism.  Pointed-to object must be
//                 a ToolDataUser
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/03
//-------------------------------------------------------------------------

/* Example Usage:
  
  TDUPtr<Battleship> my_ptr = new Battleship;
  
  // ...
  // do some stuff
  // ...
  
  if ( my_ptr == NULL ) 
    PRINT_INFO("You sunk my Battleship!\n");
  else
    PRINT_INFO("Location: %d, %d\n", my_ptr->x(), my_ptr->y() );
*/


#ifndef TDU_PTR_HPP
#define TDU_PTR_HPP

#include "ToolData.hpp"
#include "CubitMessage.hpp"
#include "CubitUtilConfigure.h"

  // Base class for TDUPtr - used because TDUPtr is a template,
  // so a common base type for all the template types is required.
class TDUPtrListNode
{
  public:
  
  inline TDUPtrListNode() : nextInTD(0) {}
  
  private:
  
  friend class TDPtr;

  virtual void nullify() = 0;
  TDUPtrListNode* nextInTD;
};

  // The actuall pointer class -- this is the one you
  // should be using.
template <class T> class TDUPtr : public TDUPtrListNode
{
public:

  inline TDUPtr() : ptrTo(0) {}
  
  inline TDUPtr( T* ptr ) : ptrTo(0)
    { change_ptr(ptr); }

  inline TDUPtr( const TDUPtr<T>& copy ) : TDUPtrListNode(), ptrTo(0)
    { change_ptr(copy.ptrTo); }
  
  inline virtual ~TDUPtr() 
    { change_ptr(0); }
  
  inline T& operator*() const   { return *ptrTo; }
  inline T* operator->() const  { return ptrTo; }
  inline operator bool() const { return ptrTo != 0; }
  
  inline TDUPtr<T>& operator=( const TDUPtr<T>& copy )
    { change_ptr(copy.ptrTo); return *this; }
    
  inline TDUPtr<T>& operator=( T* ptr )
    { change_ptr(ptr); return *this; }
  
  inline bool operator==( const TDUPtr<T>& other ) const
    { return ptrTo == other.ptrTo; }
  inline bool operator!=( const TDUPtr<T>& other ) const
    { return ptrTo != other.ptrTo; }
  inline bool operator==( const T* ptr ) const
    { return ptrTo == ptr; }
  inline bool operator!=( const T* ptr ) const
    { return ptrTo != ptr; }

  inline void change_ptr( T* ptr );
  
  inline T* ptr() const
    { return ptrTo; }

private:

  virtual void nullify()
    { ptrTo = 0; }

  T* ptrTo;

};

  // The ToolData used to associate pointer objects with TDU.
class CUBIT_UTIL_EXPORT TDPtr : public ToolData
{
public:
  static void add_to_TD( TDUPtrListNode* add, ToolDataUser* to );
  static void remove_from_TD( TDUPtrListNode* remove, ToolDataUser* from );
  
  virtual ~TDPtr();
  
private:

  inline TDPtr() : listHead(0) { }
  
  static int is_ptr( const ToolData* td )
    { return dynamic_cast<const TDPtr*>(td) != 0; }
  
  TDUPtrListNode* listHead;
};

template <class T> inline void TDUPtr<T>::change_ptr( T* ptr )
{
  if( ptr != ptrTo ) 
  {
    TDPtr::remove_from_TD(this, ptrTo);
    ptrTo = ptr;
    TDPtr::add_to_TD(this, ptrTo);
  }
}


#endif

