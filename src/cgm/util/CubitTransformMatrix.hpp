//       Class: CubitTransformMatrix
//
// Description: A 4-Dimensional Matrix.  Essentially the same as
//              a generic CubitMatrix, except that it has some
//              extra 3D transformation functions.
//
//              All transformations are pre-multiplications,
//              meaning that M*V will transform a point V
//              in the same order transformations are applied to M.
//
//      Owner: Darryl Melander

#ifndef CUBIT_MATRIX_4D_HPP
#define CUBIT_MATRIX_4D_HPP

#include "CubitMatrix.hpp"
#include "CubitVector.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitTransformMatrix : public CubitMatrix
{
public:
  CubitTransformMatrix();
  CubitTransformMatrix(const CubitTransformMatrix& from);
  ~CubitTransformMatrix();
  
  CubitTransformMatrix& translate(const CubitVector& v);
  CubitTransformMatrix& translate (double x, double y, double z);
  
  CubitTransformMatrix& rotate(double degrees, const CubitVector& vector);
  CubitTransformMatrix& rotate(double degrees, char axis);
  CubitTransformMatrix& rotate(double degrees, const CubitVector& axis_from,
                        const CubitVector& axis_to);
  
  CubitTransformMatrix& reflect(const CubitVector& vector);

  CubitTransformMatrix& scale_about_origin (const CubitVector& scale);
  CubitTransformMatrix& scale_about_origin (double x, double y, double z);
  CubitTransformMatrix& scale_about_origin (double scale);
#ifdef BOYD15
  CubitTransformMatrix& scale_about_point  (const CubitVector& scale,
                                     const CubitVector& about_point);
#endif

#ifdef BOYD15
  CubitTransformMatrix& transpose();
#endif
  CubitTransformMatrix& inverse();
  
  CubitVector operator* (const CubitVector& point) const;
  friend CubitVector operator* (const CubitVector& point,
                                const CubitTransformMatrix& matrix);
  
  CubitTransformMatrix operator*(const CubitTransformMatrix& matrix) const;
  CubitMatrix operator*(const CubitMatrix& matrix) const;
  CubitTransformMatrix operator*(double val) const;
    
  //! return the origin of this system
  CubitVector origin() const;
  //! return the x-axis
  CubitVector x_axis() const;
  //! return the y-axis
  CubitVector y_axis() const;
  //! return the z-axis
  CubitVector z_axis() const;
  
  // convenience helper for making transform matrices
  static CubitTransformMatrix construct_matrix(const CubitVector& origin,
                                               const CubitVector& x_axis,
                                               const CubitVector& y_axis);

  void print_me() const;
};

#endif

