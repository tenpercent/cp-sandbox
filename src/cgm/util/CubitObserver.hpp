//- Class: CubitObserver
//- Description: Observer class - *very* rudimentary implementation
//- Owner:  Tim Tautges

#ifndef CUBITOBSERVER_HPP
#define CUBITOBSERVER_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "CubitEvent.hpp"
#include "DLIList.hpp"
#include "CubitUtilConfigure.h"

class CubitObservable;

class CUBIT_UTIL_EXPORT CubitObserver
{
public:

  CubitObserver();
  
  virtual ~CubitObserver();

  static CubitStatus notify_static_observers(CubitObservable *observable,
                                             const CubitEvent &observer_event,
                                             CubitBoolean from_observable = CUBIT_FALSE);
    //- notify function for static observers; calls notify_observer for all static observers
  
  static void init_static_observers();
  static void term_static_observers();
    //- Access function to create and destroy the list of observed.
  
  CubitStatus watch_observable(CubitObservable *observable)
    { return register_observable(observable); }
    //- Cause this observer to watch the specified observable.
  
    //- Request that this observer stop watching the specified observable
  CubitStatus stop_watching_observable(CubitObservable *observable)
    { return unregister_observable(observable, CUBIT_FALSE); }
  
   
  virtual CubitStatus notify_observer(
     CubitObservable *observable, const CubitEvent &observer_event,
     CubitBoolean from_observable = CUBIT_FALSE) = 0;
    // This should be called through the
    // CubitObservable::notify_observers generic notify function.
    // NOTE: VERY IMPORTANT: when implementing this function, leaf
    // classes should call register_observable() or
    // unregister_observable() if the event is one of construction or
    // destruction, resp.; this makes sure the observable knows we're
    // observing it, or unhooks us from the observable.  A count of
    // observables being observed is also kept, to prevent stale
    // pointers on observables

  CubitStatus register_observable(CubitObservable *observable);
    //- Cause this observer to watch the specified observable.
  
  CubitStatus unregister_observable(CubitObservable *observable,
                                    CubitBoolean from_observable = CUBIT_FALSE);
    //- remove this observer from the observable's list, and if successful,
    //- decrement the observable count  
  CubitStatus register_static_observer(CubitObserver *observer);
    //- add this observer to the static observer list

  CubitStatus unregister_static_observer(CubitObserver *observer);
    //- remove this observer from the static observer list

private:

  int observableCount;
    //- number of observables to which this observer is registered
  static DLIList<CubitObserver*> *staticObservers;
    //- a static list of observers that observe everybody

  
}; // 

#endif

