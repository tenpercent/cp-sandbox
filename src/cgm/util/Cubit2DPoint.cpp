//- Class: Cubit2DPoint
//- Description: This file defines the Cubit2DPoint class.
//- Owner: Steve Storm
//- Checked by:

#include <math.h>
#include "Cubit2DPoint.hpp"
#include "CubitMessage.hpp"

#include "DLIList.hpp"

void
Cubit2DPoint::min_max( const Cubit2DPoint &pnt2, 
                       double &xmin, double &xmax,
                       double &ymin, double &ymax ) const
{
   if( xVal < pnt2.x() )
   {
      xmin = xVal;
      xmax = pnt2.x();
   }
   else
   {
      xmin = pnt2.x();
      xmax = xVal;
   }

   if( yVal < pnt2.y() )
   {
      ymin = yVal;
      ymax = pnt2.y();
   }
   else
   {
      ymin = pnt2.y();
      ymax = yVal;
   }
}

void 
Cubit2DPoint::update_min_max( Cubit2DPoint &min,
                              Cubit2DPoint &max ) const
{
   min.x( CUBIT_MIN( xVal, min.x() ) );
   min.y( CUBIT_MIN( yVal, min.y() ) );

   max.x( CUBIT_MAX( xVal, max.x() ) );
   max.y( CUBIT_MAX( yVal, max.y() ) );
}  

CubitBoolean 
Cubit2DPoint::is_on_line_segment( const Cubit2DPoint &end1,
                                  const Cubit2DPoint &end2,
                                  double tol ) const
{
   // Check bounding box
   double xmin, xmax, ymin, ymax;
   end1.min_max( end2, xmin, xmax, ymin, ymax );
   if( xVal < xmin-tol || yVal < ymin-tol ||
       xVal > xmax+tol || yVal > ymax+tol ) 
   {
      return CUBIT_FALSE;
   }

   // Check area of triangle (formula in if is twice the area)
   
	if( fabs( end1.x()*end2.y() + end1.y()*xVal +
	      end2.x()*yVal - end2.y()*xVal -
	      end1.x()*yVal - end1.y()*end2.x() ) < tol )
      return CUBIT_TRUE;
   
   return CUBIT_FALSE;
}
