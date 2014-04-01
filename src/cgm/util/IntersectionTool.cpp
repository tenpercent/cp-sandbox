// file: IntersectionTool.cpp
// author: Michael Stephenson
//

#include <math.h>
#include "IntersectionTool.hpp"
#include "CubitVector.hpp"
#include "CubitMessage.hpp"
#include "GeometryDefines.h"
#include "DLIList.hpp"


  
// double IntersectionTool::distance_point_line(const double point[3], 
//                                              const double start[3],
//                                              const double end[3], 
//                                              double &t)
// {
//   int i;
//   double distSq = 0.0;
//   for (i = 0; i < 3; i++)
//     distSq += (point[i] - end[i]) * (point[i] - end[i]);
//   if (distSq < mTolerance) {
//     t = 1.0;
//     return 0.0;
//   }

//   distSq = 0.0;
//   for (i = 0; i < 3; i++)
//     distSq += (point[i] - start[i]) * (point[i] - start[i]);
//   if (distSq < mTolerance) {
//     t = 0.0;
//     return 0.0;
//   }

//   double div = 0.0;
//   for (i = 0; i < 3; i++)
//     div += (end[i] - start[i]) * (end[i] - start[i]);
//   if (div < mTolerance) {
//     return -1.0;
//   }

//   t = sqrt(distSq)/div;
  
//   double pt[3];
//   for (i = 0; i < 3; i++)
//     pt[i] = (end[i] - start[i]) * t + start[i];

//   distSq = 0;
//   for (i = 0; i < 3; i++)
//     distSq += (point[i] - pt[i]) * (point[i] - pt[i]);
  
//   return sqrt(distSq);
// }

double IntersectionTool::parametric_position(const double node[3],
                                           const double pt1[3],
                                           const double pt2[3])
{
  int i;
  double dist_sq, t;
  
    // check for end-points
  double p13[3];
  for (i = 0; i < 3; i++)
    p13[i] = node[i] - pt1[i];
  dist_sq = p13[0] * p13[0] +  p13[1] * p13[1] +  p13[2] * p13[2];
  if (dist_sq < mToleranceSquared) {
    t = 0.0;
    return t;
  }

  double p23[3];
  for (i = 0; i < 3; i++)
    p23[i] = node[i] - pt2[i];
  dist_sq = p23[0] * p23[0] +  p23[1] * p23[1] +  p23[2] * p23[2];
  if (dist_sq < mToleranceSquared) {
    t = 1.0;
    return t;
  }
  
    // t is parametric distance along vector p12
  t = 0.0;
  double p12[3];
  for (i = 0; i < 3; i++)
    p12[i] = pt2[i] - pt1[i];
  
    // point1 and point2 are coincident if dot1 is zero
  double dot1 = p12[0] * p12[0] + p12[1] * p12[1] + p12[2] * p12[2];
  if (dot1 > -mToleranceSquared && dot1 < mToleranceSquared)
    return CUBIT_DBL_MAX;
  t = (p13[0] * p12[0] + p13[1] * p12[1] + p13[2] * p12[2]) / dot1;
  if (t > -mTolerance && t <  mTolerance)
    t = 0.0;
  else if ((t - 1.0) > -mTolerance &&
           (t - 1.0) <  mTolerance)
    t = 1.0;
  return t;
}

int IntersectionTool::point_on_polyline(CubitVector& pt, DLIList<CubitVector*> &pt_list)
{
  int i, ret;
  double t, distance;
  double pt_coords[3];
  double line_coords1[3];
  double line_coords2[3];

  ret = 0;

  pt_coords[0] = pt.x();
  pt_coords[1] = pt.y();
  pt_coords[2] = pt.z();

  pt_list.reset();
  CubitVector *last_pt = pt_list.get_and_step();
  for(i=pt_list.size(); i>1; i--)
  {
    CubitVector *next_pt = pt_list.get_and_step();

    line_coords1[0] = last_pt->x();
    line_coords1[1] = last_pt->y();
    line_coords1[2] = last_pt->z();
    line_coords2[0] = next_pt->x();
    line_coords2[1] = next_pt->y();
    line_coords2[2] = next_pt->z();

    distance = distance_point_line(pt_coords, line_coords1, line_coords2, t);

    if(distance > -GEOMETRY_RESABS && distance < GEOMETRY_RESABS)
    {
      i = 1;
      ret = 1;
    }
    else
      last_pt = next_pt;
  }
  return ret;
}

double IntersectionTool::distance_point_line(const double node[3], 
                                             const double pt1[3],
                                             const double pt2[3], 
                                             double &t)
{
  double dist_sq;
  int i;
  double p12[3];
  for (i = 0; i < 3; i++)
    p12[i] = pt2[i] - pt1[i];
  
  t = parametric_position(node, pt1, pt2);
  if ( t == CUBIT_DBL_MAX )
    return -1.;
  //if ( t== 0.0 || t == 1.0)
  //  return 0.0;
    // is t on vector p12 or its infinite extension
  if (t < 0.0 || t > 1.0)
    return -1.0;
    // calculate point on p12
  double p[3];
  for (i = 0; i < 3; i++)
    p[i] = pt1[i] + t * p12[i];
    
    // return distance from node to point on p12
  dist_sq = 0.0;
  for (i = 0; i < 3; i++)
    dist_sq += (p[i] - node[i]) * (p[i] - node[i]);

  return sqrt(dist_sq);
}

CubitBoolean IntersectionTool::ray_tri_test(const double start[3],
                                            const double dir[3],
                                            const double vert0[3],
                                            const double vert1[3], 
                                            const double vert2[3],
                                            double &t, double &u, double &v)
{
  int i;
  
    // find vectors for two edges sharing vert0
  double edge1[3];
  double edge2[3];
  for (i = 0; i < 3; i++) {
    edge1[i] = vert1[i] - vert0[i];
    edge2[i] = vert2[i] - vert0[i];
  }
  
    // calculate determinate
  double pvec[3];
  pvec[0] = dir[1] * edge2[2] - dir[2] * edge2[1];
  pvec[1] = dir[2] * edge2[0] - dir[0] * edge2[2];
  pvec[2] = dir[0] * edge2[1] - dir[1] * edge2[0];
  double det =
    edge1[0] * pvec[0] + edge1[1] * pvec[1] + edge1[2] * pvec[2];
  
    // if determinate is near zero, the ray is in plane of triangle
  if (det > -mTolerance && det < mTolerance) {
    return CUBIT_FALSE;
  }
  
    // calculate distance from vert0 to ray origin
  double tvec[3];
  for (i = 0; i < 3; i++)
    tvec[i] = start[i] - vert0[i];
  
    // calculate U parameter and test bounds
  double inv_det = 1.0/det;
  u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2]) * inv_det;
  if (u > -mTolerance && u < mTolerance) u = 0.0;
  else if ((u - 1.0) > -mTolerance &&
           (u - 1.0) <  mTolerance) u = 1.0;
  if (u < 0.0 || u > 1.0) {
    return CUBIT_FALSE;
  }

    // calculate V parameter and test bounds
  double qvec[3];
  qvec[0] = tvec[1] * edge1[2] - tvec[2] * edge1[1];
  qvec[1] = tvec[2] * edge1[0] - tvec[0] * edge1[2];
  qvec[2] = tvec[0] * edge1[1] - tvec[1] * edge1[0];
  v = (dir[0] * qvec[0] + dir[1] * qvec[1] + dir[2] * qvec[2]) * inv_det;
  if (v > -mTolerance && v < mTolerance) v = 0.0;
  else if ((v - 1.0) > -mTolerance &&
           (v - 1.0) <  mTolerance) v = 1.0;
  if (v < 0.0 || (u + v - 1.0) > mTolerance) {    
    return CUBIT_FALSE;
  }

    // calculate T, ray intersects triangle
  t = (edge2[0] * qvec[0] + edge2[1] * qvec[1] + 
       edge2[2] * qvec[2]) * inv_det;
  if (t > -mTolerance && t < mTolerance) t = 0.0;
  else if ((t - 1.0) > -mTolerance &&
           (t - 1.0) <  mTolerance) t = 1.0;

  return CUBIT_TRUE;
}


CubitBoolean IntersectionTool::skew_line_test(const double start1[3], 
                                              const double end1[3],
                                              const double start2[3], 
                                              const double end2[3],
                                              double &t, double &u)
{
  t = -1.0;
  u = -1.0;
  
  int i;
  double p13[3], p43[3];
  for (i = 0; i < 3; i++) {
    p13[i] = start1[i] - start2[i];
    p43[i] = end2[i]   - start2[i];
  }
  double len_sq1 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
  if (len_sq1 < mToleranceSquared)
    return CUBIT_FALSE;

  double p21[3];
  for (i = 0; i < 3; i++) {
    p21[i] = end1[i] - start1[i];
  }
  double len_sq2 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];
  if (len_sq2 < mToleranceSquared)
    return CUBIT_FALSE;

  double fact;
  if (len_sq2 < len_sq1) fact = 10.0/sqrt(len_sq1);
  else                   fact = 10.0/sqrt(len_sq2);
  for (i = 0; i < 3; i++) {
    p13[i] *= fact;
    p43[i] *= fact;
    p21[i] *= fact;
  }

  double d1343, d4321, d1321, d4343, d2121;
  d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
  d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
  d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
  d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
  d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

  double denom = d2121 * d4343 - d4321 * d4321;
  if (denom > -mTolerance && denom < mTolerance)
    return CUBIT_FALSE;
  double numer = d1343 * d4321 - d1321 * d4343;

  t = numer / denom;
  if (t > -mTolerance && t < mTolerance) t = 0.0;
  else if ((t - 1.0) > -mTolerance &&
           (t - 1.0) <  mTolerance) t = 1.0;
  if (t < 0.0 || t > 1.0)
    return CUBIT_FALSE;
   
  u = (d1343 + d4321 * t) / d4343;
  if (u > -mTolerance && u < mTolerance) u = 0.0;
  else if ((u - 1.0) > -mTolerance &&
           (u - 1.0) <  mTolerance) u = 1.0;
  if (u < 0.0 || u > 1.0)
    return CUBIT_FALSE;

  return CUBIT_TRUE;
}

CubitStatus IntersectionTool::closest_points_on_segments( CubitVector &p0,
                                                          CubitVector &p1,
                                                          CubitVector &p2,
                                                          CubitVector &p3,
                                                          CubitVector &point_1,
                                                          CubitVector &point_2,
                                                          double &sc, double &tc)
{
  CubitVector   u = p1 - p0;
  CubitVector   v = p3 - p2;
  CubitVector   w = p0 - p2;
  double    a = u%u;     //|u|  always >= 0
  double    b = u%v;
  double    c = v%v;     //|v| always >= 0
  double    d = u%w;
  double    e = v%w;
  double    D = a*c - b*b;       // always >= 0
  double    sN, sD = D;      // sc = sN / sD, default sD = D >= 0
  double    tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
  if (D < GEOMETRY_RESABS) { // the lines are almost parallel
    sN = 0.0;
    tN = e;
    tD = c;
  }
  else {                // get the closest points on the infinite lines
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0) {       // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = e;
      tD = c;
    }
    else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }
  
  if (tN < 0) {           // tc < 0 => the t=0 edge is visible
    tN = 0.0;
      // recompute sc for this edge
    if (-d < 0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  }
  else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
    tN = tD;
      // recompute sc for this edge
    if ((-d + b) < 0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else {
      sN = (-d + b);
      sD = a;
    }
  }
  sc = CUBIT_DBL_MAX;
  tc = CUBIT_DBL_MAX;
    //If these are going to be zero then do it...
  if ( sN < CUBIT_RESABS && sN > -CUBIT_RESABS )
    sc = 0.0;
  if ( tN < CUBIT_RESABS && tN > -CUBIT_RESABS )
    tc = 0.0;

    // finally do the division to get sc and tc
  if ( sD < CUBIT_RESABS && sD > -CUBIT_RESABS && sc != 0.0 )
  {
    PRINT_ERROR("About to divide by zero in closest_points_on_segments.\n");
    return CUBIT_FAILURE;
  }
  if ( tD < CUBIT_RESABS && tD > -CUBIT_RESABS && tc != 0.0 )
  {
    PRINT_ERROR("About to divide by zero in closest_points_on_segments.\n");
    return CUBIT_FAILURE;
  }
  if ( sc != 0.0 )
    sc = sN / sD;
  if ( tc != 0.0 )
    tc = tN / tD;
  
  point_1 = p0 + sc*u;
  point_2 = p2 + tc*v;
  return CUBIT_SUCCESS;
}
