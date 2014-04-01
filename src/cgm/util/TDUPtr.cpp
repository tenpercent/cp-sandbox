#include "ToolDataUser.hpp"
#include "TDUPtr.hpp"


void TDPtr::remove_from_TD( TDUPtrListNode* ptr, ToolDataUser* tdu )
{
  if ( !tdu )
    return;
  
  TDPtr* td = (TDPtr*)(tdu->get_TD(&is_ptr));
  assert(td != NULL);
  
  if (ptr == td->listHead) {
    td->listHead = ptr->nextInTD;
    if (!td->listHead ) {
      delete tdu->remove_TD(&is_ptr);
    }
  }
  else {
    TDUPtrListNode* prev = td->listHead;
    while (prev->nextInTD != ptr) {
      prev = prev->nextInTD;
      assert(prev != NULL);
    }
    prev->nextInTD = ptr->nextInTD;
  }
  
  ptr->nextInTD = 0;
}

void TDPtr::add_to_TD( TDUPtrListNode* ptr, ToolDataUser* tdu )
{
  if( !tdu )
    return;
  
  TDPtr* td = (TDPtr*)(tdu->get_TD(&is_ptr));
  if( ! td ) {
    tdu->add_TD(td = new TDPtr);
  }
  
  assert(ptr && !ptr->nextInTD);
  ptr->nextInTD = td->listHead;
  td->listHead = ptr;
}

TDPtr::~TDPtr()
{
  while (listHead) {
    TDUPtrListNode* dead = listHead;
    listHead = dead->nextInTD;
    dead->nextInTD = 0;
    dead->nullify();
  }
}


    
