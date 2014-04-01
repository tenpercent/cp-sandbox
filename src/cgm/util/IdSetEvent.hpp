//- Class: IdSetEvent
//- Description: Class describing an event when an observable's id has 
//-              been changed. Contains the old and new ids.  
//- Owner:  

#ifndef IDSET_EVENT_HPP
#define IDSET_EVENT_HPP

#include "CubitEvent.hpp"

class IdSetEvent: public CubitEvent
{
public:
  IdSetEvent(int old_id, int new_id ); 
  
 int get_old_id() const;
 int get_new_id() const;
  
private:
  int oldId;
  int newId;
};

inline
IdSetEvent::IdSetEvent( int old_id, int new_id ) 
    : CubitEvent( ID_SET )
{
  oldId = old_id; 
  newId = new_id; 
}

inline
int IdSetEvent::get_old_id() const
{
  return oldId;
}

inline 
int IdSetEvent::get_new_id() const
{
  return newId;
}

#endif 

