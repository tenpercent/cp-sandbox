//---------------------------------------------
// Class: GeomPoint
// Description: Simple point class.
// Created by: David R. White
// Date: 7/9/02
//---------------------------------------------

#ifndef GEOMPOINT_HPP
#define GEOMPOINT_HPP
#include "CubitVector.hpp"

class RefEntity;

class GeomPoint
{
private:
  RefEntity *myOwner;
  CubitVector myPosition;
public:
  GeomPoint( const double x, const double y, const double z,
             RefEntity *owner);
  GeomPoint( const CubitVector &pos, RefEntity *owner );
  ~GeomPoint();
  
  CubitVector& coordinates();
  
  void owner(RefEntity *owner);
  
  RefEntity* owner();
  
};
inline GeomPoint::GeomPoint( const CubitVector &pos, RefEntity *owner )
{
  myOwner = owner;
  myPosition.set(pos);
}
inline GeomPoint::GeomPoint( const double x, const double y, const double z,
                             RefEntity *owner)
{
  myOwner = owner;
  myPosition.set(x,y,z);
}
inline GeomPoint::~GeomPoint()
{}
inline CubitVector& GeomPoint::coordinates()
{return myPosition;}
inline void GeomPoint::owner(RefEntity *owner)
{myOwner = owner;}
inline RefEntity* GeomPoint::owner()
{return myOwner;}
#endif

