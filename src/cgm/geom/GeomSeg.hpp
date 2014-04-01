//---------------------------------------------
// Class: GeomSeg
// Description: Simple segment class.
// Created by: David R. White
// Date: 7/9/02
//---------------------------------------------

#ifndef GEOMSEG_HPP
#define GEOMSEG_HPP

#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "GeomPoint.hpp"
class RefEntity;

class GeomSeg
{
private:
  RefEntity *myOwner;
  GeomSeg *prevSeg, *nextSeg;
  GeomPoint *myStart, *myEnd;
  CubitBox boundingBox;
public:
  GeomSeg( GeomPoint *start, GeomPoint *end, RefEntity *owner );
  ~GeomSeg();
  CubitBox& bounding_box();
  void owner(RefEntity *owner);
  RefEntity* owner();
  GeomPoint* get_start();
  GeomPoint* get_end();
  void set_next( GeomSeg *tmp_next);
  GeomSeg* get_next();
  void set_prev( GeomSeg *tmp_prev);
  GeomSeg* get_prev();
  
};
inline GeomSeg::GeomSeg( GeomPoint *start,
                         GeomPoint *end,
                         RefEntity *owner )
{
  myOwner = owner;
  myStart = start;
  myEnd = end;
  prevSeg = NULL;
  nextSeg = NULL;
  CubitVector min, max;
  CubitVector p1 = start->coordinates();
  CubitVector p2 = end->coordinates();
  if ( p1.x() < p2.x() ) {
    min.x(p1.x());
    max.x(p2.x());
  } else {
    min.x(p2.x());
    max.x(p1.x());
  }
  if ( p1.y() < p2.y() ) {
    min.y(p1.y());
    max.y(p2.y());
  } else {
    min.y(p2.y());
    max.y(p1.y());
  }
  if ( p1.z() < p2.z() ) {
    min.z(p1.z());
    max.z(p2.z());
  } else {
    min.z(p2.z());
    max.z(p1.z());
  }
  boundingBox = CubitBox(min, max);
}
inline GeomSeg::~GeomSeg()
{}
inline CubitBox& GeomSeg::bounding_box()
{return boundingBox;}
inline void GeomSeg::owner(RefEntity *owner)
{myOwner = owner;}
inline RefEntity* GeomSeg::owner()
{return myOwner;}
inline GeomPoint* GeomSeg::get_start()
{return myStart;}
inline GeomPoint* GeomSeg::get_end()
{return myEnd;}
inline void GeomSeg::set_next( GeomSeg *tmp_next)
{nextSeg = tmp_next;}
inline GeomSeg* GeomSeg::get_next()
{return nextSeg;}
inline void GeomSeg::set_prev( GeomSeg *tmp_prev)
{prevSeg = tmp_prev;}
inline GeomSeg* GeomSeg::get_prev()
{return prevSeg;}



#endif

