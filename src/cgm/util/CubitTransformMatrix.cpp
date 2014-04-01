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

#include "CubitTransformMatrix.hpp"
#include "CubitMessage.hpp"

CubitTransformMatrix::CubitTransformMatrix() 
    : CubitMatrix(4)
{
    // Just creates a 4x4 matrix set to identity
}

CubitTransformMatrix::CubitTransformMatrix(const CubitTransformMatrix& from)
    : CubitMatrix(from)
{
    // Just copies it
}


CubitTransformMatrix::~CubitTransformMatrix()
{
    // Just deletes it
}



CubitTransformMatrix& CubitTransformMatrix::translate(const CubitVector& v)
{
  return translate(v.x(), v.y(), v.z());
}


CubitTransformMatrix& CubitTransformMatrix::translate (double x, double y, double z)
{
  set (0, 3, get(0, 3) + x);
  set (1, 3, get(1, 3) + y);
  set (2, 3, get(2, 3) + z);
  
  return *this;
}

CubitTransformMatrix& CubitTransformMatrix::rotate(double degrees,
                                                   const CubitVector& vector)
{
  double angle, ct, st;
    // Make a copy of vector so we don't have to change it
  CubitVector axis = vector;
  
    // Convert degrees to radians
  angle = -degrees * CUBIT_PI/180.0;
  
    // Take sin and cos
  ct = cos(angle);
  st = sin(angle);
  
    // Normalize the axis vector
  axis.normalize();
  
    // Make an identity matrix
  CubitTransformMatrix mat;
  
    // Setup some calculations that occur repeatedly
  double one_minus_cos = 1.0 - ct;
  double dx_ct = axis.x() * one_minus_cos;
  double dy_ct = axis.y() * one_minus_cos;
  double dz_ct = axis.z() * one_minus_cos;
  double dx_st = axis.x() * st;
  double dy_st = axis.y() * st;
  double dz_st = axis.z() * st;
  
    // Set the values in the matrix
  mat.set(0, 0,     ct + axis.x() * dx_ct);
  mat.set(1, 0, -dz_st + axis.x() * dy_ct);
  mat.set(2, 0,  dy_st + axis.x() * dz_ct);
  
  mat.set(0, 1,  dz_st + axis.y() * dx_ct);
  mat.set(1, 1,     ct + axis.y() * dy_ct);
  mat.set(2, 1, -dx_st + axis.y() * dz_ct);
  
  mat.set(0, 2, -dy_st + axis.z() * dx_ct);
  mat.set(1, 2,  dx_st + axis.z() * dy_ct);
  mat.set(2, 2,     ct + axis.z() * dz_ct);
  
    // Premultiply the matrix
  *this = mat * *this;
    // Return
  return *this;
}

CubitTransformMatrix& CubitTransformMatrix::rotate(double degrees, char axis)
{
  assert (axis == 'x' || axis == 'y' || axis == 'z');
  
    // Convert to Radians, Get the sine and cosine
  double angle = degrees * CUBIT_PI/180.;
  double s, c;
  s = sin(angle);
  c = cos(angle);
  
    // Make an Identity matrix
  CubitTransformMatrix mat;
    // Place values in appropriate places
  switch (axis)
  {
    case 'x':
      mat.set(1, 1, c);
      mat.set(2, 2, c);
      mat.set(1, 2, -s);
      mat.set(2, 1, s);
      break;
    case 'y':
      mat.set(0, 0, c);
      mat.set(0, 2, s);
      mat.set(2, 0, -s);
      mat.set(2, 2, c);
      break;
    case 'z':
      mat.set(0, 0, c);
      mat.set(1, 0, s);
      mat.set(0, 1, -s);
      mat.set(1, 1, c);
      break;
  }
  
    // Perform Pre-Multiplication
  *this = mat * *this;
  
  return *this;
}

CubitTransformMatrix& CubitTransformMatrix::reflect(const CubitVector& vector)
{
  double a, b, c, d;
  CubitVector axis = vector;
  axis.normalize();

  a = axis.x();
  b = axis.y();
  c = axis.z();
  
  d = sqrt(b*b + c*c);
  
    // Make an Identity matrix
  CubitTransformMatrix mat;
  if(d)
  {
      // Place values in appropriate places for negative rotate about x
    mat.set(1, 1, c/d);
    mat.set(1, 2, -b/d);
    mat.set(2, 1, b/d);
    mat.set(2, 2, c/d);
  
      // Perform Pre-Multiplication
    *this = mat * *this;
    mat.set_to_identity();
  }
  
    // Place values in appropriate places for negative rotate about y
  mat.set(0, 0, d);
  mat.set(0, 2, -a);
  mat.set(2, 0, a);
  mat.set(2, 2, d);

  
    // Perform Pre-Multiplication
  *this = mat * *this;
  mat.set_to_identity();
  
    // Place values in appropriate places for reflect across z
  mat.set(2, 2, -1);

  
    // Perform Pre-Multiplication
  *this = mat * *this;
  mat.set_to_identity();
  
    // Place values in appropriate places for rotate about y
  mat.set(0, 0, d);
  mat.set(0, 2, a);
  mat.set(2, 0, -a);
  mat.set(2, 2, d);

  
    // Perform Pre-Multiplication
  *this = mat * *this;
  mat.set_to_identity();

  if(d)
  {
      // Place values in appropriate places for rotate about x
    mat.set(1, 1, c/d);
    mat.set(1, 2, b/d);
    mat.set(2, 1, -b/d);
    mat.set(2, 2, c/d);
  
      // Perform Pre-Multiplication
    *this = mat * *this;
    mat.set_to_identity();
  }
    
  return *this;

}

CubitTransformMatrix& CubitTransformMatrix::rotate(double degrees,
                                     const CubitVector& axis_from,
                                     const CubitVector& axis_to)
{
    // Translate so that axis_from is at origin
  translate (-axis_from); 
  
    // Rotate about specified axis
  rotate (degrees, axis_to - axis_from);
  
    // Translate back
  translate (axis_from);
  
  return *this;
}

CubitTransformMatrix& CubitTransformMatrix::scale_about_origin (const CubitVector& scale)
{
  return scale_about_origin(scale.x(), scale.y(), scale.z());
}

CubitTransformMatrix& CubitTransformMatrix::scale_about_origin (double x, double y, double z)
{
  CubitTransformMatrix mat;
  mat.set(0, 0, x);
  mat.set(1, 1, y);
  mat.set(2, 2, z);
  
    // Perform Pre-Multiplication
  *this = mat * *this;
  
  return *this;
}

CubitTransformMatrix& CubitTransformMatrix::scale_about_origin (double scale)
{
  return scale_about_origin(scale, scale, scale);
}

CubitTransformMatrix& CubitTransformMatrix::inverse()
{
  CubitMatrix matrix = *this;
  matrix = matrix.inverse();
  for( int ii = 0; ii < 4; ii++ )
  {
    for( int jj = 0; jj < 4; jj++ )
    {
      set( ii, jj, matrix.get(ii,jj) );
    }
  }
  return *this;
}


  // Post-multiplication of a point (M*V)
CubitVector CubitTransformMatrix::operator* (const CubitVector& point) const
{
    // Make a sub-matrix, multiply the point by it.
  CubitVector vec = (this->sub_matrix(3, 3))*point;
    // Handle the fourth column here
  vec.x(vec.x() + get(0, 3));
  vec.y(vec.y() + get(1, 3));
  vec.z(vec.z() + get(2, 3));
  
  return vec;
}


  // point * matrix
CubitVector operator* (const CubitVector& point,
                       const CubitTransformMatrix& matrix)
{
    // Make a 1x4 matrix, multiply matrix by it.
  CubitMatrix m1(1,4);
  m1.set(0, 0, point.x());
  m1.set(0, 1, point.y());
  m1.set(0, 2, point.z());
  m1.set(0, 3, 1);
  
    // Perform the multiplication
  m1 = m1 * matrix;
    // The result is a 1x4
  
    // Put the results into a vector (dividing by w), and return
  double w = m1.get(0,3);
  return CubitVector(m1.get(0,0)/w, m1.get(0,1)/w, m1.get(0,2)/w);
}

CubitTransformMatrix CubitTransformMatrix::operator*(
  const CubitTransformMatrix& matrix) const
{
  CubitTransformMatrix rv;
  CubitMatrix temp(4);
  temp = CubitMatrix::operator*(matrix);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      rv.set(i, j, temp.get(i,j));
  return rv;
}

CubitMatrix CubitTransformMatrix::operator*(const CubitMatrix& matrix) const
{
  return CubitMatrix::operator*(matrix);
}

CubitTransformMatrix CubitTransformMatrix::operator*(double val) const
{
  CubitTransformMatrix rv;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      rv.set(i, j, this->get(i,j)*val);
  return rv;
}

void CubitTransformMatrix::print_me() const
{
  PRINT_INFO("%8.4f  %8.4f  %8.4f  %8.4f\n",
    get(0,0), get(0,1), get(0,2), get(0,3));
  PRINT_INFO("%8.4f  %8.4f  %8.4f  %8.4f\n",
    get(1,0), get(1,1), get(1,2), get(1,3));
  PRINT_INFO("%8.4f  %8.4f  %8.4f  %8.4f\n",
    get(2,0), get(2,1), get(2,2), get(2,3));
  PRINT_INFO("%8.4f  %8.4f  %8.4f  %8.4f\n",
    get(3,0), get(3,1), get(3,2), get(3,3));
}

//! return the origin of this system
CubitVector CubitTransformMatrix::origin() const
{
  return CubitVector(this->get(0,3), this->get(1,3), this->get(2,3));
}

//! return the x-axis
CubitVector CubitTransformMatrix::x_axis() const
{
  CubitMatrix tmp1(4,1);
  tmp1.set(0,0,1);
  tmp1.set(1,0,0);
  tmp1.set(2,0,0);
  tmp1.set(3,0,0);
  CubitMatrix tmp2 = (*this) * tmp1;
  return CubitVector(tmp2.get(0,0), tmp2.get(1,0), tmp2.get(2,0));
}

//! return the y-axis
CubitVector CubitTransformMatrix::y_axis() const
{
  CubitMatrix tmp1(4,1);
  tmp1.set(0,0,0);
  tmp1.set(1,0,1);
  tmp1.set(2,0,0);
  tmp1.set(3,0,0);
  CubitMatrix tmp2 = (*this) * tmp1;
  return CubitVector(tmp2.get(0,0), tmp2.get(1,0), tmp2.get(2,0));
}

//! return the z-axis
CubitVector CubitTransformMatrix::z_axis() const
{
  CubitMatrix tmp1(4,1);
  tmp1.set(0,0,0);
  tmp1.set(1,0,0);
  tmp1.set(2,0,1);
  tmp1.set(3,0,0);
  CubitMatrix tmp2 = (*this) * tmp1;
  return CubitVector(tmp2.get(0,0), tmp2.get(1,0), tmp2.get(2,0));
}


CubitTransformMatrix CubitTransformMatrix::construct_matrix(const CubitVector& origin,
                                             const CubitVector& x_axis,
                                             const CubitVector& y_axis)
{
  CubitTransformMatrix mat;
  double angle = x_axis.interior_angle(CubitVector(1,0,0));
  CubitVector axis = angle == 180 ? CubitVector(0,1,0) : CubitVector(1,0,0) * x_axis;
  mat.rotate(angle, axis);
  CubitVector y = mat * CubitVector(0,1,0);
  angle = y_axis.interior_angle(y);
  axis = angle == 180 ? CubitVector(0,0,1) : y * y_axis;
  mat.rotate(angle, axis);
  mat.translate(origin);
  return mat;
}

