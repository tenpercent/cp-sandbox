//- Class: SurfaceOverlapFacet
//- Description: Facet definition class for efficient processing
//-              for SurfaceOverlapTool.  
//- Owner: Steve Storm
//- Created: January 26, 2003

#ifndef SurfaceOverlapFacet_HPP
#define SurfaceOverlapFacet_HPP

#include "CubitDefines.h"
#include "AnalyticGeometryTool.hpp"
#include "GMem.hpp"
#include "CubitBox.hpp"
#include "CubitGeomConfigure.h"

#ifndef CUBIT_MIN_3
#define CUBIT_MIN_3(a,b,c)             (( (a) < (b) ? (a) : (b) ) < \
                                        ( (c) ) ? \
                                        ( (a) < (b) ? (a) : (b) ) : \
                                        ( (c) ))
#endif
#ifndef CUBIT_MAX_3
#define CUBIT_MAX_3(a,b,c)             (( (a) > (b) ? (a) : (b) ) > \
                                        ( (c) ) ? \
                                        ( (a) > (b) ? (a) : (b) ) : \
                                        ( (c) ))
#endif

class CUBIT_GEOM_EXPORT SurfaceOverlapFacet
{

public:
  
  SurfaceOverlapFacet( GPoint pnt[3] );
  ~SurfaceOverlapFacet();
  
  double distance( SurfaceOverlapFacet &other_facet );
  double angle( SurfaceOverlapFacet &other_facet );
  double projected_overlap( SurfaceOverlapFacet &other_facet );

  bool bbox_overlap( double tol, SurfaceOverlapFacet &other_facet ) 
  { return boundingBox.overlap( tol, other_facet.boundingBox ); }

  CubitBox bounding_box() { return boundingBox; }
  void draw( int color ); 
  
  CubitVector centroid();

  CubitVector smallest_edge_midpoint();

protected:
   
private:

  Triangle3 t;
  CubitBox boundingBox;
  static AnalyticGeometryTool *agt;

};

#endif



