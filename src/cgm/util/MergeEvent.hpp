//- Class: MergeEvent
//- Description: Class describing a merge event, two observables being merged
//-              together, one being kept, and the other lost. 
//- Owner:  

#ifndef MERGE_EVENT_HPP
#define MERGE_EVENT_HPP

#include "CubitEvent.hpp"

class MergeEvent: public CubitEvent
{
public:
  MergeEvent( CubitObservable *lost_entity, CubitObservable *kept_entity );
  
 CubitObservable *get_lost_entity() const;
 CubitObservable *get_kept_entity() const;
  
private:
  CubitObservable *lostEntity;
  CubitObservable *keptEntity;
};

inline
MergeEvent::MergeEvent( CubitObservable *lost_entity, 
                        CubitObservable *kept_entity ) 
    : CubitEvent( ENTITIES_MERGED )
{
  lostEntity = lost_entity;
  keptEntity = kept_entity;
}

inline
CubitObservable *MergeEvent::get_lost_entity() const
{
  return lostEntity;
}

inline 
CubitObservable *MergeEvent::get_kept_entity() const
{
  return keptEntity;
}

#endif 

