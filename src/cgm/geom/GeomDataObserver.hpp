//-------------------------------------------------------------------------
// Filename      : GeomDataObserver.hpp
//
// Special Notes :
//
// Creator       : David R White.
//
// Creation Date : 10/20/02
//
// Owner         : David R. White
//-------------------------------------------------------------------------
#ifndef GEOMDATA_OBSERVER_HPP
#define GEOMDATA_OBSERVER_HPP

#include "CubitObserver.hpp"
#include "CubitGeomConfigure.h"
class RefEntity;

class CUBIT_GEOM_EXPORT GeomDataObserver : public CubitObserver
{
public:
    ///
    /// The following are specific functions for making this a dynamic
    /// observer.
    ///
  static GeomDataObserver* get( RefEntity* on_this );
  static GeomDataObserver* create( RefEntity* on_this );

  virtual CubitStatus notify_observer(
    CubitObservable* watched,
    const CubitEvent& event,
    CubitBoolean );

  virtual ~GeomDataObserver();

    ///
     /// Now the reset are more specific to the GeomDataObserver.
    ///

  double get_measure()
    {return myMeasure;}
    ///
    /// Returns the cached measurment for the RefEntity.
    ///

  void set_measure(double val)
    { measureSet = CUBIT_TRUE;
    myMeasure = val;}
    ///
    /// Sets the measure value.
    ///
    
  CubitBoolean is_measure_set()
    {return measureSet;}

private:

  GeomDataObserver( RefEntity* watched );
  RefEntity* myRefEntity;
  double myMeasure;
  CubitBoolean measureSet;
  
};

#endif 

