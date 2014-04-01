#ifndef COMPOSITE_COMBINE_EVENT_HPP
#define COMPOSITE_COMBINE_EVENT_HPP

#include "CubitEvent.hpp"

class CompositeCombineEvent : public CubitEvent
{
  public:
    CompositeCombineEvent(const int event_type)
      : CubitEvent( event_type ),
        keep_entity(NULL), 
        delete_entity(NULL)
    {}
  
    CompositeCombineEvent(const int        event_type,
                          CubitObservable *keep_ptr,
                          CubitObservable *delete_ptr )
      : CubitEvent( event_type ),
        keep_entity(keep_ptr), 
        delete_entity(delete_ptr)
    {}
    
    CubitObservable* keep_entity;
   
    CubitObservable* delete_entity;
};

#endif
  
