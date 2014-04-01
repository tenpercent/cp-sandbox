//- Class: CubitPolygon
//- Description: This file defines the CubitPolygon class.
//- Owner: Steve Storm
//- Checked by:

#include <math.h>
#include "CubitPolygon.hpp"
#include "CubitMessage.hpp"

CubitPolygon::CubitPolygon()
{
}

CubitPolygon::CubitPolygon( Cubit2DPoint &start_pnt)
{
   add_point( start_pnt );
}

CubitPolygon::CubitPolygon( DLIList<Cubit2DPoint*> &point_list )
{
   int i;
   point_list.reset();
   for ( i = point_list.size(); i--; )
   {
      Cubit2DPoint *pnt = point_list.get_and_step();         
      add_point( *pnt );
   }
}

CubitPolygon::~CubitPolygon()
{
   int i;
   for ( i = pointList.size(); i--; )
      delete pointList.get_and_step();
}

void 
CubitPolygon::add_point( Cubit2DPoint &new_pnt )
{

   // Add the point
   pointList.append( new Cubit2DPoint( new_pnt ) );
   
   // Update bounding box
   if( pointList.size() == 1 )
   {
      minCoord = new_pnt;
      maxCoord = new_pnt;
   }
   else
      new_pnt.update_min_max( minCoord, maxCoord );
}

CubitPointContainment 
CubitPolygon::pnt_containment( Cubit2DPoint &pnt,
                               double tol )
{
   // Method is to fire a ray in negative x-direction and check the number
   // of intersections with the polygon.  Odd - in the loop, Even - outside
   // of the loop.
   
   // Returns: CubitPointContainment ( CUBIT_PNT_OUTSIDE = 0, 
   //                                  CUBIT_PNT_INSIDE = 1, 
   //                                  CUBIT_PNT_BOUNDARY = 2 )
   
   int i,
      c = 0; // Check variable
   
   // First check to see if the point is outside of the bounding box
   // of the polygon - if so, return OUT.
   
   if( pnt.x() < minCoord.x()-tol || pnt.y() < minCoord.y()-tol ||
      pnt.x() > maxCoord.x()+tol || pnt.y() > maxCoord.y()+tol ) 
   {
      return CUBIT_PNT_OUTSIDE;
   }
   
   // Loop on line segments of this polygon
   Cubit2DPoint *start_pnt, *end_pnt;
   pointList.reset();
   for ( i=pointList.size(); i--; )
   {
     start_pnt = pointList.get_and_step();
     end_pnt = pointList.get();

     // Check if point lies on the line segment.  This is necessary because
     // the "fire-ray" code which checks if the point is inside or outside 
     // the polygon would return outside if the point is on the boundary.
     if( pnt.is_on_line_segment( start_pnt, end_pnt, tol ) )
       return CUBIT_PNT_BOUNDARY;

     // Check if point is in bounds of y of current segment (is it a 
	   // candidate to cross segment?)
     if( start_pnt->y()<=pnt.y() && pnt.y()<end_pnt->y() ||
          end_pnt->y()<=pnt.y() && pnt.y()<start_pnt->y() ) 
     {
       // It's a candidate
       // Check if ray fired in negative x-direction crosses segment
       
       if( pnt.x() < (end_pnt->x()-start_pnt->x()) * (pnt.y()-start_pnt->y()) / 
           (end_pnt->y()-start_pnt->y()) + start_pnt->x() ) 
       {
         
         // Keeps track of even or odd number of crossings.
         //  0-even
         //  1-odd
         
         c = !c;
         
       }
     }
   }
   
   if( c )
      // Odd number of crossings
      return CUBIT_PNT_INSIDE;
   else
      // Even number of crossings
      return CUBIT_PNT_OUTSIDE;
}

CubitStatus
CubitPolygon::centroid_area( Cubit2DPoint &centroid, double &area )
{
  // Algorithm from Graphics Gems
  int n = pointList.size();

  if( pointList.size() < 3 )
    return CUBIT_FAILURE;

  int i, j;
  double ai, atmp = 0.0, xtmp = 0.0, ytmp = 0.0;

  Cubit2DPoint *p1, *p2;
  pointList.reset();
  for( i=n-1, j=0; j<n; i=j, j++ )
  {
    p1 = pointList.get_and_step();
    p2 = pointList.get();
    ai = p1->x() * p2->y() - p2->x() * p1->y();
    //ai = x[i] * y[j] - x[j] * y[i];
    atmp += ai;
    xtmp += ( p2->x() + p1->x() ) * ai;
    //xtmp += (x[j] + x[i]) * ai;
    ytmp += ( p2->y() + p1->y() ) * ai;
    //ytmp += (y[j] + y[i]) * ai;
  }

  area = atmp / 2.0;
  if (atmp != 0.0)
  {
    centroid.set( xtmp / (3.0 * atmp), ytmp / (3.0 * atmp) );
    return CUBIT_SUCCESS;
  }

  return CUBIT_FAILURE;
}

