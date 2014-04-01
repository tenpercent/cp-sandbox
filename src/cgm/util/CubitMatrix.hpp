//- Class: CubitMatrix
//-
//- Description: This file defines the CubitMatrix class which is a
//- standard NxM Matrix. All relevant arithmetic operators are
//- overloaded so CubitMatrix's can be used similar to built-in types.
//-
//- Owner: Dan Goodrich
//- Checked by:


#ifndef CUBITMATRIX_HPP
#define CUBITMATRIX_HPP

#include "CubitVector.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitMatrix
{
public:
  
    //- Heading: Constructors and Destructor
  CubitMatrix();
    //- Default constructor. 3x3 init to zeros.
  
  CubitMatrix( const int num_rows, const int num_cols);
    //- Initializes NxM matrix to zeros.
  
  CubitMatrix( const CubitVector& vec1, const CubitVector& vec2,
               const CubitVector& vec3);
    //- Constructor: create 3x3 Matrix from three vector components
    //  where each vector forms one column of the matrix.
  
  CubitMatrix( const CubitVector& vec1, const CubitVector& vec2 );
    //- Constructor: create 3x3 Matrix from two vectors (otimes)
    //  ( vec1 otimes vec2 )_ij = vec1_i * vec2_j
  
  CubitMatrix( const int n );
    //- Constructor: create n x n Identity matrix

  CubitMatrix( const CubitMatrix& matrix);
    //- Copy Constructor
  
  virtual ~CubitMatrix();
    //- destructor
  
  
    //- Heading: Set and Inquire Functions
  int num_rows() const { return numRows; }
  int num_cols() const { return numCols; }
  
  void print_matrix() const;
    //- Prints matrix

  void set_to_identity();
  
    //- Change Matrix component row n, column m to val.
  void set(const int n, const int m, const double val)
    {
      assert (n >= 0 && n < numRows);
      assert (m >= 0 && m < numCols);
      matrixPtr[n][m] = val;
    }

    //- add the value to the component row n col m
  void add(const int n, const int m, const double val)
    {
      matrixPtr[n][m] += val;
    }
  
    //- Gets the values of the matrix at position (n,m)
  double get( int n, int m ) const
    {
      assert (n >= 0 && n < numRows);
      assert (m >= 0 && m < numCols);
      return matrixPtr[n][m];
    }
  
    //- Heading: Matrix Algebra functions 
    //  (assert if incompatible element sizes)
  CubitMatrix operator= (const CubitMatrix& matrix);
  CubitMatrix operator* (const CubitMatrix& matrix) const;
  CubitVector operator* (const CubitVector& vector) const;
  CubitMatrix operator* (double val ) const;
  CubitMatrix operator/ (double val ) const;
  CubitMatrix operator+ (const CubitMatrix& matrix) const;
  CubitMatrix operator- (const CubitMatrix& matrix) const;
  
  CubitMatrix& operator+=(const CubitMatrix &matrix);
  CubitMatrix& operator*=(const double multiplier);
    //- Scale the matrix by a linear factor: {this = this * multiplier},
  
  CubitMatrix inverse(); // inverts 1x1, 2x2, or 3x3 matrix
    // asserts if singular
  CubitMatrix symm() const; // returns matrix^transpose * matrix
 
#ifdef BOYD15
  CubitBoolean diagonally_dominant() const;  // returns true if diagonally dominant
#endif
  CubitBoolean positive_definite() const;  // returns true if matrix is positive definite
  
  double determinant () const;
  double inf_norm () const; // infinity norm
  double frobenius_norm_squared() const; // square of frobenius norm of A
  double frobenius_norm_squared_symm() const; // square of frobenius norm of A^t A
  double frobenius_norm_squared_adj() const; // square of frobenius norm of adjoint A
  double frobenius_norm_squared_inv() const; // square of frobenius norm of A-inverse
  double condition() const;  // condition number of A using frobenius norm

  double cofactor (const int row, const int col) const;
#ifdef BOYD15
  double trace () const;
#endif
  CubitMatrix adjoint() const;
  CubitMatrix transpose() const;

     // row and col indicate the row and column to leave out of the
    // resulting matrix. 
  CubitMatrix sub_matrix( const int row, const int col ) const;
  
    // routines to perform Gaussian elimination with pivoting and scaling
    // on 3x3 matrix
  int  gauss_elim (CubitVector &b);
  int  factor (CubitVector &pivot);
  void solve (CubitVector &b, const CubitVector& pivot);

    // routines to solve NxN system (from Numerical Recipes in C)
  CubitStatus solveNxN( CubitMatrix& rhs, CubitMatrix& coef );
  CubitStatus ludcmp( double *indx, double& d );
  CubitStatus lubksb( double *indx, double *b );
  
private:
  double **matrixPtr;
  int numRows;
  int numCols;
};

#endif



