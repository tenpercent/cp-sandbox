// file: GeomIntersectionTool.hpp
// author: Michael Stephenson
//

#ifndef INTERSECTION_TOOL_HPP
#define INTERSECTION_TOOL_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

class CubitVector;
template <class X> class DLIList;

class CUBIT_UTIL_EXPORT IntersectionTool
{
public:
  IntersectionTool(const double &tol = 0.0001);
  virtual ~IntersectionTool() {}

  virtual CubitStatus closest_points_on_segments(CubitVector &p0,
                                                 CubitVector &p1,
                                                 CubitVector &p2,
                                                 CubitVector &p3,
                                                 CubitVector &point_1,
                                                 CubitVector &point_2,
                                                 double &sc, double &tc);
    //- Finds the two closest points on the line segments {p0,p1} and {p2,p3},
    //- and returns them in point_1 and point_2.  Will return failure
    //- if a divide by zero is attempted.

  double distance_point_line(const double point[3], const double start[3], 
                             const double end[3], double &t);
  int point_on_polyline(CubitVector& pt, DLIList<CubitVector*> &pt_list);

  double parametric_position(const double node[3],
                           const double pt1[3],
                           const double pt2[3]);
  
    ///
    /// Returns the parametric position of the node on the
    /// line segement between the points pt1 and pt2.  This
    /// value is between 0 and 1.  If the node is closest
    /// to pt1 it will return 0, if it is closest to pt2, it
    /// will return 1.  If the node is off the line, it
    /// returns the closest parametric value.
    /// If there is an error, like pt1 and pt2 are equal,
    /// then the function returns CUBIT_DBL_MAX.
    ///

  virtual CubitStatus initialize(){return CUBIT_SUCCESS;}

protected:
  double mTolerance;
  double mToleranceSquared;

#ifdef BOYD15
  CubitBoolean does_intersect_2D(const double start1[3], 
                                 const double end1[3],
                                 const double start2[3], 
                                 const double end2[3]);

  CubitBoolean does_intersect_3D(const double start1[3], 
                                 const double end1[3],
                                 const double start2[3], 
                                 const double end2[3]);

  CubitBoolean does_intersect_tri(const double start1[3], 
                                  const double end1[3],
                                  const double tri1[3], 
                                  const double tri2[3], 
                                  const double tri3[3]);
  
  CubitBoolean point_line_test(const double point[3], 
                               const double start[3], 
                               const double end[3], double &t);
#endif

  CubitBoolean ray_tri_test(const double start[3], const double dir[3],
                            const double tri1[3], const double tri2[3], 
                            const double tri3[3],
                            double &t, double &alpha, double &beta);

  CubitBoolean skew_line_test(const double start1[3], const double end1[3],
                              const double start2[3], const double end2[3],
                              double &t, double &u);

  void tolerance(const double &tol);
  double tolerance(void) const;
  
private:
};

inline IntersectionTool::IntersectionTool(const double &tol) 
    : mTolerance(tol)
{
  mToleranceSquared = tol * tol;
}

#endif // INTERSECTION_TOOL_HPP

