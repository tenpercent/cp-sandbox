//- Class: CurveOverlapFacet
//- Description: Facet definition class for efficient processing
//-              for OverlapTool.  
//- Owner: Corey Ernst 
//- Created: August 11, 2005

#include "CurveOverlapFacet.hpp"
#include "IntersectionTool.hpp"

//AnalyticGeometryTool* CurveOverlapFacet::agt = AnalyticGeometryTool::instance();

// Constructor
CurveOverlapFacet::CurveOverlapFacet( GPoint p[2] )
{
  p0.set( p[0].x, p[0].y, p[0].z);
  p1.set( p[1].x, p[1].y, p[1].z);

  CubitVector min;
  min.set( CUBIT_MIN( p[0].x, p[1].x ),
           CUBIT_MIN( p[0].y, p[1].y ),
           CUBIT_MIN( p[0].z, p[1].z ) );
  
  CubitVector max;
  max.set( CUBIT_MAX( p[0].x, p[1].x ),
           CUBIT_MAX( p[0].y, p[1].y ),
           CUBIT_MAX( p[0].z, p[1].z ) );
  
  boundingBox.reset( min, max );

  facetLength = 0.0;
}

// Destructor
CurveOverlapFacet::~CurveOverlapFacet()
{
}

double CurveOverlapFacet::distance_overlapping( CurveOverlapFacet *other_facet )
{
  //determine which is the smaller facet
  CurveOverlapFacet *short_facet, *long_facet;
  if( this->length() < other_facet->length() )
  {
    short_facet = this;
    long_facet = other_facet;
  }
  else
  {
    short_facet = other_facet;
    long_facet = this;
  }

  CubitVector long_facet_vector(long_facet->p1 - long_facet->p0);
  double long_facet_length = long_facet_vector.length();
  long_facet_vector.normalize();
  
  CubitVector l_p0_to_s_p0( short_facet->p0 - long_facet->p0 );   
  CubitVector l_p0_to_s_p1( short_facet->p1 - long_facet->p0 );   

  double dot_to_s_p0 = l_p0_to_s_p0 % long_facet_vector; 
  double dot_to_s_p1 = l_p0_to_s_p1 % long_facet_vector; 

  //no overlap
  if( ( dot_to_s_p0 <= 0.0 && dot_to_s_p1 <= 0.0) || 
      ( dot_to_s_p0 >= long_facet_length  && 
        dot_to_s_p1 >= long_facet_length ) )
  {
    return 0.0;
  }

  //overlap like this:
  // o-----------o
  //        o---------------------o           
  if(  dot_to_s_p0 > 0.0 && dot_to_s_p1 < 0.0 ) 
    return dot_to_s_p0;
  if( dot_to_s_p1 > 0.0 && dot_to_s_p0 < 0.0 ) 
    return dot_to_s_p1;
  if( dot_to_s_p0 > long_facet_length &&
      dot_to_s_p1 < long_facet_length )
    return long_facet_length - dot_to_s_p1;
  if( dot_to_s_p1 > long_facet_length &&
      dot_to_s_p0 < long_facet_length )
    return long_facet_length - dot_to_s_p0;
   

  //overlap like this 
  // o-----------o
  // o---------------------o           
  // or this:
  //    o-----------o
  // o---------------------o           
  // or this:
  // o---------------------o           
  // o---------------------o           
  if( (dot_to_s_p0 <= long_facet_length &&
       dot_to_s_p1 < long_facet_length ) || 
      (dot_to_s_p1 <= long_facet_length &&
       dot_to_s_p0 < long_facet_length ) ) 
    return fabs( dot_to_s_p0 - dot_to_s_p1 );

  return 0.0;

/*
  //create vectors
  CubitVector v0( long_facet->p1 - long_facet->p0 ); 
  CubitVector v1( short_facet->p0 - long_facet->p0 ); 
  CubitVector v2( short_facet->p1 - long_facet->p0 ); 
  CubitVector normalized_v0 = v0;
  normalized_v0.normalize();
  
  //project v1 and v2 onto v0
  double dot1 = normalized_v0 % v1; 
  double dot2 = normalized_v0 % v2; 
  v1 = normalized_v0 * dot1;
  v2 = normalized_v0 * dot2;
  double length0 = v0.length();
  double length1 = v1.length();
  double length2 = v2.length();

  //no overlap at all
  if( (length1 > long_facet->length() && length2 > long_facet->length() ) ||
      (dot1 <= 0 && dot2 <= 0 ) )
  {
    return 0.0;
  }
  else if( dot1 <= 0 && dot2 > 0 )
  {
    return length2; 
  }
  else if( length2 <= 0 && length1 > 0 )
  {
    return length1; 
  }
  else if( length1 > length0 )
  {
    return length0 - length2;
  }
  else if( length2 > length0 )
  {
    return length0 - length1;
  }
  //else
    return fabs( length2 - length1 );
*/
/*
  CubitVector tmp_vec1( p0.x, p0.y, p0.z);
  CubitVector tmp_vec2( p1.x, p1.y, p1.z);

  CubitVector tmp_vec3( other_facet->p0.x, other_facet->p0.y, other_facet->p0.z);
  CubitVector tmp_vec4( other_facet->p1.x, other_facet->p1.y, other_facet->p1.z);


  IntersectionTool int_tool( GEOMETRY_RESABS );
 
  CubitVector closest_point_seg_1;
  CubitVector closest_point_seg_2;
  double sc, tc;
  CubitStatus stat = int_tool.closest_points_on_segments( p0, p1, 
                            other_facet->p0, other_facet->p1,
                            closest_point_seg_1, closest_point_seg_2,
                            sc, tc );
  //Make sure the closest points aren't the end points.  If they are
  //and we are within tolerance, it may be that the tolerance is too big
  //cause we shouldn't be at a cross if the closest point is an end point...
  if ( sc > 0. && sc < 1. && tc > 0. && tc < 1. &&
       closest_point_seg_1.within_tolerance(closest_point_seg_2, 0.001) )
  {
    


  }
  return 0.0;
  */
}

double CurveOverlapFacet::angle( CurveOverlapFacet *other_facet )
{
  CubitVector vec1 = p1 - p0;
  CubitVector vec2 = other_facet->p1 - other_facet->p0;

  double angle = vec1.interior_angle( vec2 );
  return angle;
}

double CurveOverlapFacet::length()
{
  if( facetLength == 0.0 )
  {
    CubitVector tmp_vec = p1 - p0;
    facetLength = tmp_vec.length();
  }
  return facetLength;
}

double CurveOverlapFacet::facet_to_facet_distance( CurveOverlapFacet *other_facet )
{
  CubitVector u = p1 - p0;
  CubitVector v = other_facet->p1 - other_facet->p0;
  CubitVector w = p0 - other_facet->p0;
  double a = u%u;
  double b = u%v;
  double c = v%v;
  double d = u%w;
  double e = v%w;
  double D = a*c - b*b;
  double sc, sN, sD = D;
  double tc, tN, tD = D;

  // compute the line parameters of the two closest points
  if( D < GEOMETRY_RESABS )  ///the lines are almost parallel
  {
    sN = 0.0;        // force using point P0 on segment S1
    sD = 1.0;        // to prevent possible division by 0.0 later
    tN = e;
    tD = c;
  }
  else // get the closest points on the infinite lines
  {
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0) 
    {       // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = e;
      tD = c;
    }
    else if (sN > sD) 
    {  // sc > 1 => the s=1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < 0.0) 
  {           // tc < 0 => the t=0 edge is visible
    tN = 0.0;
    // recompute sc for this edge
    if (-d < 0.0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else 
    {
      sN = -d;
      sD = a;
    }
  }
  else if (tN > tD) 
  {      // tc > 1 => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-d + b) < 0.0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else 
    {
      sN = (-d + b);
      sD = a;
    }
  }

  // finally do the division to get sc and tc
  sc = (fabs(sN) < GEOMETRY_RESABS ? 0.0 : sN / sD);
  tc = (fabs(tN) < GEOMETRY_RESABS ? 0.0 : tN / tD);

  // get the difference of the two closest points
  CubitVector   dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)

  return dP.length();   // return the closest distance
}
