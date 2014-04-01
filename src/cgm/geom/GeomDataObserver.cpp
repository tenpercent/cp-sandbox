//------------------------------------------------------------------------
// Class GeomDataObserver
// Description:  Observer class that stores/caches specific geometric
//               information (for example, the area for a surface.
//
// Author: David White
// Creation Date: 9/7/2003
//------------------------------------------------------------------------
#include "GeomDataObserver.hpp"
#include "RefEntity.hpp"
GeomDataObserver::GeomDataObserver(RefEntity* watched)
    : myRefEntity(watched)
{
  measureSet = CUBIT_FALSE;
    //initialize to something weird.
  myMeasure = -CUBIT_DBL_MAX;
}

GeomDataObserver::~GeomDataObserver()
{ 
  unregister_observable( myRefEntity ); 
}

GeomDataObserver* GeomDataObserver::get( RefEntity* on_this )
{
   DLIList<CubitObserver*> list;
   GeomDataObserver* eo = NULL;

   on_this->get_observer_list(list);
   for (int i = list.size(); i--; )
   {
     if ( (eo = dynamic_cast<GeomDataObserver*>(list.step_and_get()) ))
        break;
   }
   return eo;
}

GeomDataObserver* GeomDataObserver::create( RefEntity* on_this ) {
   GeomDataObserver* eo = get(on_this);
   if (eo)
     return eo;

   eo = new GeomDataObserver(on_this);
   eo->register_observable(on_this);
   return eo;
}

CubitStatus GeomDataObserver::notify_observer(
                               CubitObservable* watched,
                               const CubitEvent& event,
                               CubitBoolean )
{
   assert(watched == myRefEntity);

   switch (event.get_event_type())
   {
     case GEOMETRY_TOPOLOGY_MODIFIED:
     case TOPOLOGY_MODIFIED:
     case GEOMETRY_MODIFIED:
     case ENTITY_DESTRUCTED:
     case MODEL_ENTITY_DESTRUCTED:
       break;
     default:
       return CUBIT_SUCCESS;
   }

   /* don't call virtual functions or access
      class data after this! */
   delete this;
   return CUBIT_SUCCESS;
}

