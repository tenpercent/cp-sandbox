
#ifndef UNMERGE_EVENT_HPP
#define UNMERGE_EVENT_HPP

#include "CubitEvent.hpp"

class CubitObservable;

class UnMergeEvent : public CubitEvent
{
  public:
  
    UnMergeEvent( CubitObservable *old_ptr, CubitObservable *new_ptr )
      : CubitEvent( NEW_ENTITY_UNMERGED ),
        old_entity(old_ptr), 
        new_entity(new_ptr)
    {}
    
    CubitObservable* const old_entity;
   
    CubitObservable* const new_entity;
};

#endif
  
