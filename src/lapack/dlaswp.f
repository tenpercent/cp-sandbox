      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
c -- LAPACK auxiliary routine (version 3.1) --
c    Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c    November 2006
c
c    .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
c    .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
c    ..
c
c  Purpose
c  =======
c
c  <a name="DLASWP.18"></a><a href="http://www.netlib.org/lapack/explore-html/dlaswp.f.html#DLASWP.1">DLASWP</a> performs a series of row interchanges on the matrix A.
c  One row interchange is initiated for each of rows K1 through K2 of A.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.
c
c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On entry, the matrix of column dimension N to which the row
c          interchanges will be applied.
c          On exit, the permuted matrix.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.
c
c  K1      (input) INTEGER
c          The first element of IPIV for which a row interchange will
c          be done.
c
c  K2      (input) INTEGER
c          The last element of IPIV for which a row interchange will
c          be done.
c
c  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
c          The vector of pivot indices.  Only the elements in positions
c          K1 through K2 of IPIV are accessed.
c          IPIV(K) = L implies rows K and L are to be interchanged.
c
c  INCX    (input) INTEGER
c          The increment between successive values of IPIV.  If IPIV
c          is negative, the pivots are applied in reverse order.
c
c  Further Details
c  ===============
c
c  Modified by
c   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
c
c =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
c     ..
c     .. Executable Statements ..

c     Interchange row I with row IPIV(I) for each of rows K1 through K2.
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
      RETURN
      END
