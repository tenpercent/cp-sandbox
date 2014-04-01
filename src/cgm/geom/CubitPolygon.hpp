//- Class: CubitPolygon
//-
//- Description: This file defines the CubitPolygon class which is basically
//- an ordered list of 2D points and their 2D bounding box. 
//-
//- Owner: Steve Storm

#ifndef CUBITPOLYGON_HPP
#define CUBITPOLYGON_HPP

#include "CubitDefines.h"
#include "Cubit2DPoint.hpp"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

class CUBIT_GEOM_EXPORT CubitPolygon
{
public:
  
  //- Heading: Constructors and Destructor
  CubitPolygon();  //- Default constructor.
  
  CubitPolygon( DLIList<Cubit2DPoint*> &point_list );
  //- Constructor: create polygon from array of points.  Points are copied.

  CubitPolygon( Cubit2DPoint &point );
  //- Constructor: create polygon with initial point.  Point is copied.

  ~CubitPolygon();
  //- Destructor

  void add_point( Cubit2DPoint &point );
  //- Add a point to the end of the point list.  Point is copied.

#ifdef BOYD14
  void reset();
  //- Clean out all of the points
#endif
  
  CubitPointContainment pnt_containment( Cubit2DPoint &point,
                                         double tol = 1e-10 );
  //- Determine if the input point is inside, outside or on the boundary
  //- of the polygon.

  CubitStatus centroid_area( Cubit2DPoint &centroid, double &area );
  //- Get geometrical centroid and area of polygon

#ifdef BOYD14
  void list_coords();
  void list_bbox();
  //- Debug functions
#endif

private:
  
  DLIList<Cubit2DPoint*> pointList;
  Cubit2DPoint minCoord;  // Minimum coordinate of polygon
  Cubit2DPoint maxCoord;  // Maximum coordinate of polygon

};

#endif


