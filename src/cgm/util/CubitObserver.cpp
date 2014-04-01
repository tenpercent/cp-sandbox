#include "CubitObserver.hpp"
#include "CubitObservable.hpp"
#include "CubitEvent.hpp"

DLIList <CubitObserver*> *CubitObserver::staticObservers = NULL;
    //- a static list of observers that observe everybody

CubitObserver::CubitObserver()
{
  if (staticObservers == NULL)
  {
    PRINT_ERROR("Local staticObservers need to be initiallized first in\n"
                "your program.\n");
    assert(staticObservers != NULL);
  }
  
  observableCount = 0;
}

CubitObserver::~CubitObserver()
{
  assert(observableCount == 0);
}

CubitStatus CubitObserver::register_observable(CubitObservable *observable)
{
    //- add this observer to the observable's list, and if successful,
    //- increment the observable count
  CubitStatus success = observable->add_observer(this);
  
  if (success == CUBIT_SUCCESS)
    observableCount++;
  
  return success;
}

CubitStatus CubitObserver::unregister_observable(CubitObservable *observable,
                                                 CubitBoolean from_observable)
{
    //- remove this observer from the observable's list, and if successful,
    //- decrement the observable count
  CubitStatus success = CUBIT_SUCCESS;

    //- only call the remove function on the observable if we're not being called
    //- from the observable
  if (from_observable == CUBIT_FALSE)
    success = observable->remove_observer(this);

  if (success == CUBIT_SUCCESS)
    observableCount--;
  
  return success;
}


    //- pass this event onto any static observers, if there are any
CubitStatus CubitObserver::notify_static_observers(CubitObservable *observable,
                                                   const CubitEvent &observer_event,
                                                   CubitBoolean from_observable)
{
    //- pass this event onto any static observers, if there are any
  CubitStatus success = CUBIT_SUCCESS;

  if (NULL == staticObservers) return CUBIT_SUCCESS;

  for (int i = staticObservers->size(); i > 0; i--) {
    if ( staticObservers->get_and_step()->
         notify_observer(observable, observer_event, from_observable) 
         == CUBIT_FAILURE )
      success = CUBIT_FAILURE;
  }

  return success;
}

CubitStatus CubitObserver::register_static_observer(CubitObserver *observer)
{
    //- add this observer to the static observer list
  staticObservers->append_unique(observer);
  return CUBIT_SUCCESS;
}

CubitStatus CubitObserver::unregister_static_observer(CubitObserver *observer)
{
    //- remove this observer from the static observer list
  if (staticObservers->remove(observer))
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}
void CubitObserver::init_static_observers()
{
  if ( !staticObservers )
  {
    staticObservers = new DLIList<CubitObserver*>;
  }
}
void CubitObserver::term_static_observers()
{
  if ( staticObservers )
    delete staticObservers;
  staticObservers = NULL;
}
    //- Access function to create and destroy the list of observed.

