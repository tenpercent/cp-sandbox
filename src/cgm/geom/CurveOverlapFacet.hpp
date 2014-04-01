//- Class: CurveOverlapFacet
//- Description: Facet definition class for efficient processing
//-              for CurveOverlapTool.  
//- Owner: Corey Ernst 
//- Created: August 11, 2005

#ifndef CurveOverlapFacet_HPP
#define CurveOverlapFacet_HPP

#include "CubitDefines.h"
#include "GMem.hpp"
#include "CubitBox.hpp"
#include "CubitGeomConfigure.h"
#include "GeometryDefines.h"

class CUBIT_GEOM_EXPORT CurveOverlapFacet
{

public:
  
  friend class SurfaceOverlapTool;

  CurveOverlapFacet( GPoint pnt[2] );
  ~CurveOverlapFacet();
  
  double distance_overlapping( CurveOverlapFacet *other_facet );
  double angle( CurveOverlapFacet *other_facet );
  double length();

  bool bbox_overlap( double tol, CurveOverlapFacet *other_facet ) 
  { return boundingBox.overlap( tol, other_facet->boundingBox ); }

  double facet_to_facet_distance( CurveOverlapFacet *other_facet );

protected:
   
private:

  CubitVector p0;
  CubitVector p1;
  CubitBox boundingBox;
  double facetLength;
};

#endif



