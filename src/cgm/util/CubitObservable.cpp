#include "CubitObservable.hpp"
#include "CubitObserver.hpp"
#include "CubitDefines.h"

CubitObservable::CubitObservable()
{
  observerList = NULL;
}
  
CubitObservable::~CubitObservable()
{
  if (observerList) 
  {
    assert( (observerList->remove_all_with_value( NULL ),
             observerList->size() == 0) );
    delete observerList;
    observerList = NULL;
  }
}
  
CubitStatus CubitObservable::add_observer(CubitObserver *observer)
{
  if ( !observerList )
     observerList = new DLIList <CubitObserver*>;
  if( observerList->append_unique( observer ) == CUBIT_TRUE )
  {
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

CubitStatus CubitObservable::remove_observer(CubitObserver *observer)
{
  CubitStatus success = CUBIT_FAILURE;
  if( observerList != NULL )
  {
    int i;
    for ( i = 0; i < observerList->size(); i++ )
    {
      if ( observerList->get() == observer )
      {
        success = CUBIT_SUCCESS;
        observerList->change_to( NULL );
        break;
      }
      observerList->step();
    }
  }
  return success;
}

void CubitObservable::remove_from_observers()
{
  notify_observers( ENTITY_DESTRUCTED );
  if( observerList != NULL )
     observerList->clean_out();
}
  
void CubitObservable::get_observer_list( DLIList <CubitObserver*> &observers)
{
  if (observerList)
  {
    for( int i = observerList->size(); i > 0; i-- )
    {
      CubitObserver* observer = observerList->get_and_step();
      if( observer )
         observers.append( observer );
    }
  }
}
  
CubitStatus CubitObservable::notify_observers( const CubitEvent &event )
{
  CubitStatus result = CUBIT_SUCCESS;
  if ( ! observerList || observerList->size() == 0 )
     return result;

  CubitObserver *observer = NULL;
  for ( int i = 0; i < observerList->size(); i++ )
  {
    observerList->reset();
    observerList->step(i);
    observer = observerList->get();
    if ( observer != NULL )
       if( !observer->notify_observer(this, event) )
          result = CUBIT_FAILURE;
  }
  observerList->remove_all_with_value( NULL );

  return result;
}

CubitStatus CubitObservable::notify_all_observers( const CubitEvent &event )
{
  CubitStatus d_stat = notify_observers( event );
  CubitStatus s_stat = CubitObserver::notify_static_observers( this, event );
  return (s_stat && d_stat) ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

