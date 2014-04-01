//- Class: CubitObservable
//- Description: Observable class - *very* rudimentary implementation
//- Owner:  Tim Tautges

#ifndef CUBITOBSERVABLE_HPP
#define CUBITOBSERVABLE_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "DLIList.hpp" 
#include "CubitUtilConfigure.h"

class CubitEvent;
class CubitObserver;

class CUBIT_UTIL_EXPORT CubitObservable
{

public:

  friend class CubitObserver;

  CubitObservable();
  
  virtual ~CubitObservable();

  void get_observer_list(DLIList<CubitObserver*> &observers);
    //- copy the observer list; return just the non-null observers.

  CubitStatus notify_observers( const CubitEvent &event );
    //R CubitStatus
    //R- An AND of the result returned by each observer for the event.
    //- Notify all observers of an event.
    //- Process will proceed even if an observer returns failure for
    //- a notification.
    
  CubitStatus notify_all_observers( const CubitEvent& event );
    //- Same as notify_observers, but also notifies static observers.

protected:
  void remove_from_observers();
    //- unregister this observable from all observers		


    //- IMPORTANT: keep add_observer and remove_observer private,
    //- so that they can only be accessed from CubitObserver; all register/
    //- unregister actions go through CubitObserver!!!
  virtual CubitStatus add_observer(CubitObserver *observer);
    //- add the observer to this observable; returns CUBIT_FAILURE
    //- if this observer contained this observable already

  virtual CubitStatus remove_observer(CubitObserver *observer);
    //- remove the observer from this observable; returns CUBIT_SUCCESS
    //- if this observable contained this observer

  DLIList <CubitObserver*> *observerList;
    //- pointer to list of observers of this entity
    //- This list may contain NULL pointers. This is used to make sure
    //- each observer is notified once and only once.

    CubitObservable( const CubitObservable& );
    void operator=( const CubitObservable&);
}; // 

#endif

