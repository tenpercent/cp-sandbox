//- Class: CubitEvent
//- Description: Class describing an event. This base class is sufficient
//-              to describe basic unary events, like destroying a model 
//-              entity. Derived classes describe events requiring extra
//-              data, such as merging two model entities.
//- Owner:  

#ifndef CUBITEVENT_HPP
#define CUBITEVENT_HPP

#include "CubitDefines.h"
#include "CubitEventDefines.h"

class CubitEvent
{
public:
  
  typedef int Type;
  typedef int SubType;
  
  CubitEvent();
    // Constructor.
    // Is there a need for no-argument constructor? There usually is...
  
  CubitEvent( Type event_type );
    // Constructor.
    // Takes a Type, so CubitObserver functions taking a CubitEvent argument 
    // can be passed an (EventType) instead as a shorthand.
  
  virtual ~CubitEvent();
    // destructor

    // Get/Set event type.
  Type get_event_type() const;
  void set_event_type( Type event_type );
  
    // Subtype stuff.  This probably isn't the best way to 
    // implement this, but...
  CubitEvent( Type event_type, SubType event_subtype );

  SubType get_event_sub_type() const;
  void set_event_sub_type(SubType sub_event_type);
  
private:

  Type eventType;
    //- takes values in EventType defined in CubitDefines.hpp
  SubType eventSubType;

};

// inline functions
inline CubitEvent::CubitEvent() 
 : eventType(INVALID_EVENT_TYPE), eventSubType( INVALID_EVENT_TYPE )
{}

inline CubitEvent::CubitEvent( CubitEvent::Type event_type )
 : eventType(event_type), eventSubType( INVALID_EVENT_TYPE )
{}

inline CubitEvent::~CubitEvent() 
{}

inline CubitEvent::Type CubitEvent::get_event_type() const
{ return eventType; }

inline void CubitEvent::set_event_type( CubitEvent::Type event_type )
{ eventType = event_type; }


inline CubitEvent::CubitEvent( CubitEvent::Type event_type,
                               CubitEvent::SubType event_subtype )
        : eventType(event_type),
          eventSubType(event_subtype) 
{}

inline CubitEvent::SubType CubitEvent::get_event_sub_type() const
{ return eventSubType; }

inline void CubitEvent::set_event_sub_type(CubitEvent::SubType sub_event_type)
{ eventSubType = sub_event_type; }

#endif

