//- Class: CubitMatrix
//- Description: This file defines the CubitMatrix class.
//- Owner: Dan Goodrich
//- Checked by:

#include <assert.h>

#include "CubitMatrix.hpp"
#include "CubitMessage.hpp"
#include "CubitVector.hpp"
#include "CubitDefines.h"

CubitMatrix::CubitMatrix()
{
  numRows = 3;
  numCols = 3;
  
  matrixPtr = new double *[3];
  
  int ii;
  for(  ii = 0; ii < 3; ii++ )
  {
    matrixPtr[ii] = new double [3];
  }
  
    // Initialize matrix to zeros.     
  for( ii = 0; ii < 3; ii++ )
     for( int jj = 0 ; jj < 3; jj++ )
        matrixPtr[ii][jj] = 0.0;
}

CubitMatrix::CubitMatrix( const int n, const int m )
{
  numRows = n;
  numCols = m;
  
  matrixPtr = new double *[n];
  
  int ii;
  for(  ii = 0; ii < n; ii++ )
  {
    matrixPtr[ii] = new double [m];
  }
  
    // Initialize matrix to zeros.     
  for( ii = 0; ii < n; ii++ )
     for( int jj = 0 ; jj < m; jj++ )
        matrixPtr[ii][jj] = 0.0;
}


CubitMatrix::CubitMatrix( const int n )
{
  numRows = n;
  numCols = n;
  
  matrixPtr = new double *[n];
  
  for( int ii = 0; ii < n; ii++ )
    matrixPtr[ii] = new double [n];
  
    // Initialize matrix to identity.
  set_to_identity();
}

CubitMatrix::CubitMatrix (const CubitVector& vec1,
                          const CubitVector& vec2,
                          const CubitVector& vec3 )
{
  numRows = 3;
  numCols = 3;
  
  matrixPtr = new double *[3];
  
  for( int ii = 0; ii < 3; ii++ )
  {
    matrixPtr[ii] = new double [3];
  }
  
    // Initialize the matrix columns to the three vectors
  matrixPtr[0][0] = vec1.x();
  matrixPtr[1][0] = vec1.y();
  matrixPtr[2][0] = vec1.z();
  matrixPtr[0][1] = vec2.x();
  matrixPtr[1][1] = vec2.y();
  matrixPtr[2][1] = vec2.z();
  matrixPtr[0][2] = vec3.x();
  matrixPtr[1][2] = vec3.y();
  matrixPtr[2][2] = vec3.z();
}


CubitMatrix::CubitMatrix(const CubitVector& vec1,
                         const CubitVector& vec2 )
{             
  numRows = 3;
  numCols = 3;
  
  matrixPtr = new double *[3];
  
  for( int ii = 0; ii < 3; ii++ )
  {
    matrixPtr[ii] = new double [3];
  }
  
    // Initialize the matrix elements using otimes (outer product)
  matrixPtr[0][0] = vec1.x() * vec2.x();
  matrixPtr[1][0] = vec1.y() * vec2.x();
  matrixPtr[2][0] = vec1.z() * vec2.x();
  matrixPtr[0][1] = vec1.x() * vec2.y();
  matrixPtr[1][1] = vec1.y() * vec2.y();
  matrixPtr[2][1] = vec1.z() * vec2.y();
  matrixPtr[0][2] = vec1.x() * vec2.z();
  matrixPtr[1][2] = vec1.y() * vec2.z();
  matrixPtr[2][2] = vec1.z() * vec2.z();
}

CubitMatrix::CubitMatrix( const CubitMatrix &matrix )
{
  numRows = matrix.num_rows();
  numCols = matrix.num_cols();
  
  matrixPtr = new double *[numRows];
  
  int ii;
  for(  ii = 0; ii < numRows; ii++ )
  {
    matrixPtr[ii] = new double [numCols];
  }
  
  for( ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
       matrixPtr[ii][jj] = matrix.get( ii, jj );
  }
}



CubitMatrix::~CubitMatrix()
{
  for( int ii = 0; ii < numRows; ii++ )
  {
    delete [] ( matrixPtr[ii] );
  }
  delete [] ( matrixPtr );
}


void CubitMatrix::print_matrix() const
{
  for( int row = 0; row < numRows; row++ )
  {
    for( int col = 0; col < numCols; col++ )
       PRINT_INFO("%8.3f", matrixPtr[row][col]);
    PRINT_INFO("\n");
  }
}

// Sets this matrix equal to 'matrix'.  'this' is
// redimensioned if needed.
CubitMatrix CubitMatrix::operator=(const CubitMatrix& matrix)
{
  int i, j;
  
  if (numRows != matrix.num_rows() ||
      numCols != matrix.num_cols())
  {
      // note (vvyas, 3/2006): corrected array deletions below
    for (i = 0; i < numRows; i++)
    {
        //delete matrixPtr[i];
      delete [] matrixPtr[i];
    }
      //delete matrixPtr;
    delete [] matrixPtr;
    
    numRows = matrix.num_rows();
    numCols = matrix.num_cols();
    matrixPtr = new double*[numRows];
    for (i = 0; i < numRows; i++)
      matrixPtr[i] = new double[numCols];
  }
  
  for(i = 0; i < numRows; i++ )
    for(j = 0; j < numCols; j++)
      matrixPtr[i][j] = matrix.get(i, j);
  
  return *this;
}


// Multiply this ( size NxM ) with the input matrix ( size MxL ).
// return matrix of size NxL
CubitMatrix CubitMatrix::operator*(const CubitMatrix& matrix ) const
{
    // Check that we can multiply them.
  assert( numCols == matrix.num_rows() );
  
  CubitMatrix return_matrix( numRows, matrix.num_cols() );
  
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < matrix.num_cols(); jj++ )
    {
      double temp = 0.0;
      for( int kk = 0; kk < numCols; kk++ )
      {
          //temp += matrixPtr[ii][kk] * matrix.get( kk, jj );
        temp += get( ii, kk ) * matrix.get( kk, jj );
      }
      return_matrix.set( ii, jj, temp );
    }
  }
  return return_matrix;
}



// multiply this times the input vector
CubitVector CubitMatrix::operator* (const CubitVector& vector ) const
{
    // Check that we can multiply them.
  assert( numCols == 3 );
  
  double vec1[3];
  double vec2[3];
  
  vec2[0] = vector.x();
  vec2[1] = vector.y();
  vec2[2] = vector.z();
  
  for( int row = 0; row < numRows; row++ )
  {
    vec1[row] = 0.0;
    for( int col = 0; col < numCols; col++ )
    {
      vec1[row] += ( matrixPtr[row][col] * vec2[col] );
    }
  }
  
  return CubitVector( vec1[0],  vec1[1], vec1[2] );
}


// multiply this times the input scalar
CubitMatrix CubitMatrix::operator*( double val ) const
{
  CubitMatrix matrix( numRows, numCols );
  
  for( int row = 0; row < numRows; row++ )
  {
    for( int col = 0; col < numCols; col++ )
    {
      matrix.set( row, col,( matrixPtr[row][col] * val ) );
    }
  }
  return matrix;
}

// multiply this times the input scalar
CubitMatrix CubitMatrix::operator/( double val ) const
{
  assert( val != 0 );
  CubitMatrix matrix( numRows, numCols );
  
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
    {
      matrix.set( ii, jj,( matrixPtr[ii][jj] / val ) );
    }
  }
  return matrix;
}


// subtract this ( size NxM ) with the input matrix ( size NxM ).
// return matrix of size NxM
CubitMatrix CubitMatrix::operator-(const CubitMatrix& matrix) const
{
  CubitMatrix return_matrix( numRows, numCols );
  
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
    {
      return_matrix.set( ii, jj, matrixPtr[ii][jj] -
                         matrix.get( ii, jj ));
    }
  }
  
  return return_matrix;
}

// add this ( size NxM ) with the input matrix ( size NxM ).
// return matrix of size NxM
CubitMatrix CubitMatrix::operator+(const CubitMatrix& matrix ) const
{
   CubitMatrix return_matrix( numRows, numCols );
   
   for( int ii = 0; ii < numRows; ii++ )
   {
      for( int jj = 0; jj < numCols; jj++ )
      {
         return_matrix.set( ii, jj, matrixPtr[ii][jj] +
                            matrix.get( ii, jj ));
      }
   }

   return return_matrix;
}


CubitMatrix& CubitMatrix::operator+=( const CubitMatrix &matrix )
{
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
    {
      matrixPtr[ii][jj] += matrix.get( ii, jj );
    }
  }
  return *this;
}

CubitMatrix& CubitMatrix::operator*=(const double multiplier)
{
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
    {
      matrixPtr[ii][jj] *= multiplier;
    }
  }
  return *this;
}

// Sets the matrix to all zeros except along diagonal.
// Matrix doesn't have to be square.
void CubitMatrix::set_to_identity()
{
  for (int i = numRows; i--; )
    for (int j = numCols; j--; )
    {
      if (i == j)
        matrixPtr[i][j] = 1;
      else
        matrixPtr[i][j] = 0;
    }
}

/*
// Inverts this matrix, if it is of size NxN, and a 3x3 or
// smaller.
CubitMatrix CubitMatrix::inverse()
{
  CubitMatrix adj_matrix( numRows, numCols );
  double   det;
  
  if( numRows > 4 )
  {
//       rval = recipie_inverse();
//       return rval == CUBIT_TRUE ? CUBIT_TRUE : CUBIT_FALSE;
    PRINT_INFO("Can't handle matrice's greater than 3x3 yet.\n");
  }
  
  det = determinant();
  assert( fabs(det) > CUBIT_DBL_MIN );
  
  adj_matrix = adjoint();
  return adj_matrix * ( 1.0/det );
}
*/

// Inverts this matrix, if it is size 4x4 or bigger
// uses ludcmp and lubksb from numerical recipes.
CubitMatrix CubitMatrix::inverse()
{
  // can't invert a non-square matrix
  assert(numRows == numCols);

  CubitMatrix matrix_inverse( numRows, numCols );
  
  if (numRows <4)
  {
    double   det;
    det = determinant();
    assert( fabs(det) > CUBIT_DBL_MIN );
    double det_inv = 1./det;

    if ( numRows == 1 ) {
      det = determinant();
      assert( fabs(det) > CUBIT_DBL_MIN );
    
      matrix_inverse.set(0,0, matrixPtr[0][0]);
      matrix_inverse *= det_inv;
    }

    if ( numRows == 2 ) {
      matrix_inverse.set(0,0, matrixPtr[1][1]);
      matrix_inverse.set(1,0,-matrixPtr[1][0]);
      matrix_inverse.set(0,1,-matrixPtr[0][1]);
      matrix_inverse.set(1,1, matrixPtr[0][0]);
      matrix_inverse *= det_inv;
    }

    if ( numRows == 3 ) {
      matrix_inverse.set(0,0, matrixPtr[1][1] * matrixPtr[2][2] - matrixPtr[1][2] * matrixPtr[2][1] );
      matrix_inverse.set(1,0, matrixPtr[2][0] * matrixPtr[1][2] - matrixPtr[1][0] * matrixPtr[2][2] );
      matrix_inverse.set(2,0, matrixPtr[1][0] * matrixPtr[2][1] - matrixPtr[1][1] * matrixPtr[2][0] );
      matrix_inverse.set(0,1, matrixPtr[2][1] * matrixPtr[0][2] - matrixPtr[0][1] * matrixPtr[2][2] );
      matrix_inverse.set(1,1, matrixPtr[0][0] * matrixPtr[2][2] - matrixPtr[0][2] * matrixPtr[2][0] );
      matrix_inverse.set(2,1, matrixPtr[0][1] * matrixPtr[2][0] - matrixPtr[0][0] * matrixPtr[2][1] );
      matrix_inverse.set(0,2, matrixPtr[0][1] * matrixPtr[1][2] - matrixPtr[0][2] * matrixPtr[1][1] );
      matrix_inverse.set(1,2, matrixPtr[1][0] * matrixPtr[0][2] - matrixPtr[0][0] * matrixPtr[1][2] );
      matrix_inverse.set(2,2, matrixPtr[0][0] * matrixPtr[1][1] - matrixPtr[1][0] * matrixPtr[0][1] );
    
    }
    matrix_inverse *= det_inv;
  }
  else
  {

    // use numerical recipes Inverse of a Matrix

    int i, j;
    double d;
    double *indx = new double [numRows];
    double *col = new double [numRows];
    CubitMatrix save_matrix = *this;

    CubitStatus rv = ludcmp(indx, d);
    assert(rv == CUBIT_SUCCESS);
    for (j=0; j<numRows; j++)
    {
      for(i=0; i<numRows; i++) 
      {
        col[i] = 0.0;
      }
      col[j] = 1.0;
      rv = lubksb(indx, col);
      assert(rv == CUBIT_SUCCESS);
      for (i=0; i<numRows; i++) 
      {
        matrix_inverse.set(i,j,col[i]);
      }
    }
    delete [] indx;
    delete [] col;
    *this = save_matrix;
  }
  
  return matrix_inverse;
}

CubitBoolean CubitMatrix::positive_definite() const
{

  if ( matrixPtr[0][0] <= 0. ) { return CUBIT_FALSE; }

  double det2x2 = matrixPtr[0][0] * matrixPtr[1][1] - matrixPtr[1][0] * matrixPtr[0][1];
  if ( det2x2 <= 0. ) { return CUBIT_FALSE; }

  if ( determinant() <= 0. ) { return CUBIT_FALSE; }
  
  return CUBIT_TRUE;
}

double CubitMatrix::determinant() const
{
  double det = 0.0;
  
  if( numRows == 1 )
     det = matrixPtr[0][0];
  else if( numRows == 2 )
     det = matrixPtr[0][0] * matrixPtr[1][1] - matrixPtr[0][1]
        * matrixPtr[1][0];
  else if (numRows == 3)
     det = matrixPtr[0][0] * matrixPtr[1][1] * matrixPtr[2][2] +
         matrixPtr[0][1] * matrixPtr[1][2] * matrixPtr[2][0] +
         matrixPtr[0][2] * matrixPtr[1][0] * matrixPtr[2][1] -
         matrixPtr[2][0] * matrixPtr[1][1] * matrixPtr[0][2] -
         matrixPtr[2][1] * matrixPtr[1][2] * matrixPtr[0][0] -
         matrixPtr[2][2] * matrixPtr[1][0] * matrixPtr[0][1];  
  else
  {
    for( int jj = 0; jj < numRows; jj++ )
    {
      det += ( matrixPtr[0][jj] * cofactor( 0, jj ) );
    }
  }
  return det;
}

double CubitMatrix::cofactor( const int row, const int col ) const
{
  double c = 0.0;
  CubitMatrix matrix_sub( numRows - 1, numCols -1 );
  
  matrix_sub = sub_matrix( row, col );
  
  c = matrix_sub.determinant();
  c = (row+col)%2 ? -1*c : c;
  
  return c;
}

CubitMatrix CubitMatrix::adjoint() const
{
  CubitMatrix matrix( numRows, numRows );
  
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numRows; jj++ )
    {
      matrix.set( ii, jj, cofactor( ii, jj ) );
    }
  }
  return matrix.transpose();
}

CubitMatrix CubitMatrix::transpose() const
{
  CubitMatrix return_matrix( numCols, numRows );
  
  for( int ii = 0; ii < numRows; ii++ )
  {
    for( int jj = 0; jj < numCols; jj++ )
    {
      return_matrix.set( jj, ii, matrixPtr[ii][jj] );
    }
  }
  
  return return_matrix;
}

// Creates and returns a matrix that is a copy of 'this',
// except that the indicated row and column are left out.
CubitMatrix CubitMatrix::sub_matrix( const int row, const int col ) const
{
  CubitMatrix matrix (numRows - 1, numCols - 1);
  
  int copy_row = 0;
  for (int source_row = 0; source_row < numRows; source_row++)
  {
    if (source_row != row)
    {
      int copy_col = 0;
      for (int source_col = 0; source_col < numCols; source_col++)
      {
        if (source_col != col)
        {
          matrix.set (copy_row, copy_col, matrixPtr[source_row][source_col]);
          copy_col++;
        }
      }
      copy_row++;
    }
  }
  
  return matrix;
}

double CubitMatrix::inf_norm()  const
{ 
    // infinity norm  = max_i sum_j  | A_ij |
  double matrix_norm = 0., row_norm, v;
  for ( int ii = 0; ii < numRows; ii++ ) {
    row_norm = 0.;
    for( int jj = 0; jj < numCols; jj++ )
    {
      v = fabs( get( ii, jj ) );
      row_norm += v;
    }
    if ( row_norm > matrix_norm )
       matrix_norm = row_norm;
  }
  return matrix_norm;
}

double CubitMatrix::frobenius_norm_squared() const 
{ 
    // frobenius norm-squared  = trace( M^T M )
  
  double matrix_norm=0;
  for ( int ii = 0; ii < numRows; ii++ ) {

    for( int jj = 0; jj < numCols; jj++ )
    {
      matrix_norm += matrixPtr[ii][jj] * matrixPtr[ii][jj];
    }

  }

  return matrix_norm;

}

double CubitMatrix::frobenius_norm_squared_symm()  const
{ 
    // frobenius norm-squared 2 = trace[( M^T M )( M^T M )]
  
  double matrix_norm=0;
  for ( int ii = 0; ii < numRows; ii++ ) {

    for( int jj = 0; jj < numCols; jj++ )
    {
      double b=0;
      for ( int kk = 0; kk < numRows; kk++ ) {
	b += matrixPtr[kk][ii] * matrixPtr[kk][jj];
      }
      matrix_norm += b*b;
    }

  }

  return matrix_norm;

}

double CubitMatrix::frobenius_norm_squared_adj()  const
{ 
    // square of frobenius norm of adjoint
  
  double norm=0;

  if ( numRows == 1 ) { norm=1; }

  if ( numRows == 2 ) {
    norm = this->frobenius_norm_squared();
  }

  if ( numRows == 3 ) {
    norm = 0.5 * ( pow( this->frobenius_norm_squared(), 2 ) - this->frobenius_norm_squared_symm() );
  }

  if ( numRows > 3 ) {
    CubitMatrix adj = this->adjoint();
    norm = adj.frobenius_norm_squared();
  }

  return norm;

}

double CubitMatrix::frobenius_norm_squared_inv()  const
{ 
    // square of frobenius norm of A-inverse
  
  double det = this->determinant();
  assert( det != 0 );

  double norm=this->frobenius_norm_squared_adj()/pow(det,2);

  return norm;

}

double CubitMatrix::condition()  const
{ 
    // condition number of A using frobenius norm 
  
  double norm = ( this->frobenius_norm_squared() ) * (this->frobenius_norm_squared_inv() );

  return sqrt( norm );

}

int CubitMatrix::gauss_elim( CubitVector &b )
{
    CubitVector pivot;
    int ierr = factor( pivot );
    if ( ierr == 0 ) {  solve( b, pivot ); }
    return ierr;
}

int CubitMatrix::factor( CubitVector &pivot )
{
    double pvt[3];

    const int n=3;
    double s[3], tmp;

    int i,j;
    for ( i=0; i<n; i++ ) {
       s[i] = 0.0;
       for ( j=0; j<n; j++ ) {
          tmp = fabs( matrixPtr[i][j] );
          if ( tmp > s[i] ) {
             s[i] = tmp;
	  }
       }
       
       if ( s[i] == 0.0 ) { return(1); }
    
    }

    for ( int k=0; k<n-1; k++ ) {
       double ck = 0.0;
       int i0 = -1;
       for ( i=k; i<n; i++ ) {
          tmp = fabs( matrixPtr[i][k] / s[i] );
          if ( tmp > ck ) {
             ck = tmp;
             i0 = i;
	  }
       }

       pvt[k] = i0;
       if ( ck == 0.0 ) { return(1); }

       if ( i0 != k ) {
          for ( j=k; j<n; j++ ) {
             double swap = matrixPtr[i0][j];
             matrixPtr[i0][j] = matrixPtr[k][j];
             matrixPtr[k][j] = swap;
	  }
       }

       for ( i=k+1; i<n; i++ ) {
          double r = matrixPtr[i][k] / matrixPtr[k][k];
          matrixPtr[i][k] = r;
          for ( j=k+1; j<n; j++ ) {
             matrixPtr[i][j] -= r * matrixPtr[k][j];
	  }
       }

    }

    pivot.set( pvt[0], pvt[1], pvt[2] );
    return(0);
}

void CubitMatrix::solve( CubitVector &b, const CubitVector& pivot )
{
  double rhs[3];
  rhs[0] = b.x();
  rhs[1] = b.y();
  rhs[2] = b.z();

  double pvt[3];
  pvt[0] = pivot.x();
  pvt[1] = pivot.y();
  pvt[2] = pivot.z();

  int j;
  const int n=3;
  for ( int k=0; k<n-1; k++ ) {
     j=(int)pvt[k];
     if ( j != k ) {
        double swap = rhs[k];
        rhs[k] = rhs[j];
        rhs[j] = swap;
     }

     for ( int i=k+1; i<n; i++ ) {
        rhs[i] -= matrixPtr[i][k] * rhs[k];
     }

  }

  rhs[n-1] /= matrixPtr[n-1][n-1];

  for ( int i=n-2; i>-1; i-- ) {
     double sum=0.;
     for ( j=i+1; j<n; j++ ) {
        sum += matrixPtr[i][j] * rhs[j];
     }
     rhs[i] = ( rhs[i] - sum ) / matrixPtr[i][i];
  }

  b.set( rhs[0], rhs[1], rhs[2] );

}


// Here is the recipe for inverting a NxM matrice.
// I did not spend the time trying to convert it to Cubit style.
// Matrix is a double**
// Vector is a double*
// Scalar is a double
// int mxiRecipieInverse(Matrix M1, Matrix M2, int N)
// {
//   Matrix           M1_loc, M2_loc, M3_loc;
//   Vector           col, copycol;
//   Scalar           d;
//   int              i, j, *indx;

//   indx = ((int*)malloc((unsigned long)N*sizeof(int)))-1;

//   M1_loc  = mxInitMatrixR(1, N, 1, N);
//   M2_loc  = mxInitMatrixR(1, N, 1, N);
//   M3_loc  = mxInitMatrixR(1, N, 1, N);
//   col     = mxInitVectorR(1, N);
//   copycol = mxInitVectorR(1, N);
//   if (M1_loc == NULL || M2_loc == NULL || col == NULL || indx == NULL) 
//     return 0;
//   if (M3_loc == NULL || copycol == NULL)
//     printf("\n\nCannot use Improve function\n\n");

//             /* copy the input matrix */
//   for( i = 1; i <= N; i++ )
//     for( j = 1; j <= N; j++ ) {
//       M1_loc[i][j] = M1[i-1][j-1];
//       if (M3_loc != NULL)
//         M3_loc[i][j] = M1[i-1][j-1];
//       M2_loc[i][j] = 0.0;
//     }

//   if (!mxiLudcmp(M1_loc, N, indx, &d)) return 0;
//   for (j=1; j<=N; j++) {
//     for (i=1; i<=N; i++) {
//       col[i]=0.0;
//       if (copycol != NULL)
//         copycol[i] = 0.0;
//     }
//     col[j] = 1.0;
//     if (copycol != NULL) copycol[j] = 1.0;
//     mxiLubksb(M1_loc, N, indx, col);
//     if (copycol != NULL && M3_loc != NULL)
//       if (!mxiImprove(M3_loc, M1_loc, N, indx, copycol, col))
//         return 0;
//     for (i=1; i<=N; i++) M2_loc[i][j]=col[i];
//   }
//             /* copy the inverted matrix */
//   for( i = 1; i <= N; i++ )
//     for( j = 1; j <= N; j++ )
//       M2[i-1][j-1] = M2_loc[i][j];

//   mxFreeMatrixR(M1_loc, 1, N, 1, N);
//   mxFreeMatrixR(M2_loc, 1, N, 1, N);
//   mxFreeMatrixR(M3_loc, 1, N, 1, N);
//   mxFreeVectorR(col, 1, N);
//   mxFreeVectorR(copycol, 1, N);
//   free(indx+1);
//   return 1;
// } /* mxiRecipieInverse */

CubitStatus CubitMatrix::solveNxN( CubitMatrix& rhs, CubitMatrix& coef )
{
  if (numRows != rhs.num_rows() ||
    numRows != numCols) {
    return CUBIT_FAILURE;
  }
  int i;
  double d;
  double *indx = new double [numRows];
  double *b = new double [numRows];
  if (!indx) {
    return CUBIT_FAILURE;
  }
  CubitStatus status = ludcmp(indx, d);
  if (status == CUBIT_SUCCESS) {
    for(i=0; i<numRows; i++) {
      b[i] = rhs.get(i,0);
    }
    status = lubksb(indx, b);
    for (i=0; i<numRows; i++) {
      coef.set(i,0,b[i]);
    }
  }
  delete [] indx;
  delete [] b;
  return status;
}

// from numerical recipies in C: Decompose a NxN matrix into
// Upper and Lower trianglar (in place)
CubitStatus CubitMatrix::ludcmp( double *indx, double& d )
{
  int i, j, k, imax = -1;
  double big, tmp, sum;
  double *vv = new double [numRows];
  if (!vv) {
    return CUBIT_FAILURE;
  }

  d = 1.0; // no row interchanges yet

  // loop over rows to get implicit scale info

  for (i=0; i<numRows; ++i){
    big = 0.0;
    for (j=0; j<numRows; ++j)
      if ((tmp = fabs(matrixPtr[i][j])) > big)
        big = tmp;
    if (big == 0.0) {
        // note (vvyas, 3/2006): corrected array deletion
        // delete vv;
      delete [] vv;
      return CUBIT_FAILURE;
    }
    vv[i] = 1.0/big;
  }

  // loop over columns-Crout's method

  for (j=0; j<numRows; ++j){
    for (i=0; i<j; ++i){
      sum = matrixPtr[i][j];
      for (k=0; k<i; ++k)
        sum -= matrixPtr[i][k]*matrixPtr[k][j];
      matrixPtr[i][j] = sum;
    }
    big = 0.0;                         // initialize pivot search
    for (i=j; i<numRows; ++i){
      sum = matrixPtr[i][j];
      for (k=0; k<j; ++k)
        sum -= matrixPtr[i][k]*matrixPtr[k][j];
      matrixPtr[i][j] = sum;
      if ((tmp = vv[i]*fabs(sum)) > big) {
        big = tmp;
        imax = i;
      }
    }
    if (j != imax) {                   // do we need to change rows
      for (k=0; k<numRows; ++k) {
        tmp = matrixPtr[imax][k];
        matrixPtr[imax][k] = matrixPtr[j][k];
        matrixPtr[j][k] = tmp;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (matrixPtr[j][j] == 0.0) matrixPtr[0][0] = 1.0e-20;
    if (j != numRows-1) {             // divide by the pivot element
      tmp = 1.0/matrixPtr[j][j];
      for (i=j+1; i<numRows; ++i)
        matrixPtr[i][j] *= tmp;
    }
  }                                   // go back for next column

    // note (vvyas 3/2006): corrected array deletion
    // delete vv;
  delete [] vv;
  
  return CUBIT_SUCCESS;
}

// from numerical recipies in C: solve [mat]{x} = {b} by back
// substitution (mat = LU of mat)
CubitStatus CubitMatrix::lubksb( double *indx, double *b )
{
  int i, j, ii, ip;
  double sum;

  // do the forward substitution

  ii = -1;
  for (i=0; i<numRows; ++i){
    ip = (int)indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0)
      for (j=ii; j<=i-1; ++j)
        sum -= matrixPtr[i][j]*b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }

  // do the back substitution

  for (i=numRows-1; i>=0; --i){
    sum = b[i];
    for (j=i+1; j<numRows; ++j)
      sum -= matrixPtr[i][j]*b[j];
    b[i] = sum/matrixPtr[i][i];  // store a component of solution
  }
  return CUBIT_SUCCESS;
}


