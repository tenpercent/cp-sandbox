C ======================================================================
      Subroutine FEM3Dtet(
C ======================================================================
     &      XY1, XY2, XY3, XY4,
     &      operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &      label, D, DATA, iSYS, order, 
     &      LDA, A, nRow, nCol)
C ======================================================================
      include 'fem3Dtet.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes an elemental matrix for the bilinear form 
C
C  (1)            <D OpA(u), OpB(v)>               
C
C  where D is a tensor, OpA and OpB are linear operators, and  u and v 
C  are finite element functions, u in space "A", and v in space "B". 
C  Solution of a non-linear problem may involve a Newton-like iterative 
C  loop. In this case D may depend on a discrete function (e.g. 
C  approximation from the previous iterative step). Evaluation of D may 
C  be quite complex and may require additional data. A data array is 
C  passed inside subroutines (see formula (3)). 
C 
C  In order to compute the right hand side, we can use the following 
C  trick:
C
C (2)             f(v) = < D FEM_P0, v > 
C
C  where action of D is given by function f.
C
C   
C  *** BILINER FORMS ***
C  The finite element space "A" is assumed to be richer that the
C  finite element space "B". The possible choices for these spaces are:
C
C    FEM_P0        - piecewise constant
C    FEM_P1        - continuous piecewise linear 
C    FEM_P2        - continuous piecewise quadratic
C    FEM_P3        - continuous piecewise cubic
C    FEM_P1vector  - vector piecewise linear. The unknowns are ordered
C                    first by vertices and then by space directions (x,y,z)
C    FEM_P2vector  - vector piecewise quadratic. 
C    FEM_RT0       - lower order Raviart-Thomas finite elements
C    FEM_ND0       - lower order Nedelec finite elements
C
C    FEM_CR1       - Crouzeix-Raviart finite elements
C
C  The available operators are:
C    IDEN          - identity operator 
C    GRAD          - gradient operator
C    DIV           - divergence operator
C    CURL          - curl operator
C    DUDX          - partial derivative d/dx
C
C  A quadrature formula is chosen as follows:
C    order = 1     - quadrature formula with  1 center point
C    order = 2     - quadrature formula with  4 points inside tetrahedron
C    order = 3     - quadrature formula with  8 points inside tetrahedron
C    order = 5     - quadrature formula with 15 points inside tetrahedron
C    order = 6     - quadrature formula with 24 points inside tetrahedron
C
C  The matrix A is defined by bilinear form (1). The following rules 
C  are applied for numbering basis functions:
C      A) First, basis function associated with vertices are enumerated
C         in the same order as the vertices XYi, i = 1,2,3, and 4
C
C      B) Second, basis function associated with edges are enumerated
C         in the same order as egdes 12,13,14,23,24 and 34.
C
C         If there are more than one degree of freedom per edge, they
C         are ordered by the distance from the lowest end-point.
C
C      C) Third, basis function associated with faces are enumerated
C         in the same order as faces 123,234,341 and 412.
C      
C      D) The vector basis functions with 3 degrees of freedom per 
C         a mesh object (vertex, edge, face) are enumerated first by the
C         corresponding mesh objects and then by space coordinates, x,y,z.
C  
C     LDA   - leading dimention of A
C     nRow  - number of rows of A
C     nCol  - number of columns of A
C
C
C *** DESCRIPTION OF THE TENSOR ***  
C  The external function D has the following STANDARD format
C
C  (3)   INTEGER FUNCTION D(x,y,z, label, DATA, iSYS, Coef)
C
C    (x, y, z) - Real*8 Cartesian coordinates of a 3D point where
C                tensor Diff should be evaluated
C
C    label     - identificator of either element or face 
C
C    DATA(*)   - Real*8 user given data (an array)
C
C    iSYS(21)  - system buffer for information exchange:
C                On Input:
C                     iSYS( 3)    <- tetrahedrom number
C                     iSYS( 4:7)  <- numbers of vertices
C                     iSYS( 8:13) <- numbers of edges
C                     iSYS(14:17) <- numbers of faces
C                     iSYS(18)    <- the number of points
C                     iSYS(19)    <- the number of edges
C                     iSYS(20)    <- the number of faces
C                     iSYS(21)    <- the number of tetrahedra
C                On Output:
C                     iD  = iSYS(1) -> number of rows in tensor Diff
C                     jD  = iSYS(2) -> number of columns in tensor Diff
C
C    Coef(9, 9) - Real*8 matrix. The number of actually used rows and 
C                 columns is iD and jD, respectively. In the case of 
C                 vector  finite elements U = (u, v, w) the following  
C                 ordering should be used:
C                 u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z.
C                  
C                 Examples: 
C                   A) isotropic diffusion problem:  
C                      iD = jD = 1
C                      Coef = diffusion value at the point (x,y,z)
C
C                   B) anisotropic diffusion problem:
C                      iD = jD = 3
C                      Coef(i, j) = diffusion value at the point (x,y,z)
C
C                   C) convection problem:
C                      iD = 1, jD = 3
C                      Coef(1, j) = velocity value at the point (x,y,z)
C
C                   D) linear elasticity problem:
C                      iD = jD = 9
C
C                               |2                |     |1       1       1| 
C                               |  1   1          |     |                 |
C                               |    1       1    |     |                 |
C                             M |  1   1          |     |                 |
C                      Coef = - |        2        | + L |1       1       1|
C                             2 |          1   1  |     |                 |
C                               |    1       1    |     |                 |
C                               |          1   1  |     |                 |
C                               |                2|     |1       1       1|
C              
C                      where M and L are Lame coefficient values at 
C                      the point (x,y,z)
C
C
C  To simplify the problem of writing subroutine D(), we provide a few 
C  examples with constant coefficients. In all examples, parameters
C  x, y and z are obsolete and the coefficients are passed inside the
C  subroutine through array DATA. Examples are placed in file forlibfem.f
C                     
C  (a)  INTEGER FUNCTION ANI_D1x1_one(x, y, z, DATA, iSYS, Coef)
C
C  (b)  INTEGER FUNCTION ANI_D3x3_one(x, y, z, DATA, iSYS, Coef)
C
C  (c)  INTEGER FUNCTION ANI_Dconv_const(x, y, z, DATA, iSYS, Coef)
C
C          DATA(1:3) = constant velocity vector
C
C  (d)  INTEGER FUNCTION ANI_Delas_const(x, y, z, DATA, iSYS, Coef)
C
C          DATA(1) = constant scalar 1st Lame coefficient
C          DATA(2) = constant scalar 2nd Lame coefficient
C
C
C ======================================================================
C *** ADMISSIBLE (Y) and PLANNED (P) OPERATIONS ***
C ======================================================================
C
C                    IDEN  GRAD  DIV   CURL  DUDX  DUDY
C                  -------------------------------------
C     FEM_P0          Yes   Yes   NO     NO   Yes    P
C     FEM_P1          Yes   Yes   NO     NO   Yes    P
C     FEM_P2          Yes   Yes   NO     NO    P     P
C     FEM_P3          Yes   Yes   NO     NO    P     P
C     FEM_P1vector    Yes   Yes  Yes    Yes    P     P
C     FEM_P2vector    Yes   Yes  Yes     P     P     P
C     FEM_P2reduced    P    P     P      P     P     P
C     FEM_RT0         Yes   Yes  Yes    Yes    P     P
C     FEM_RT1          P    P     P      P     P     P
C     FEM_ND0         Yes   P    Yes    Yes    P     P
C     FEM_RT1          P    P     P      P     P     P
C     FEM_BDM1         P    P     P      P     P     P
C     FEM_CR1         Yes   Yes   NO     NO    P     P
C     FEM_CR1vector    P    Yes  Yes     P     P     P
C
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), XY4(*)

      Real*8   DATA(*)
      Integer  label, iSYS(*), D

      EXTERNAL D

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  order, LDA

      Real*8   A(LDA, *)
C ======================================================================
c Local variables
      Real*8  XYP(3, 4)
      Real*8  det, cdet, vol, s
   
      Real*8  PSI(3, 3)
      Real*8  U(9, MaxSize, MaxPointGauss), V(9, MaxSize, MaxPointGauss)
      Real*8  Diff(9, 9, MaxPointGauss),   DU(9, MaxSize, MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(3, MaxPointGauss)
      Real*8  XYL(4, AllPointGauss)

      Integer tensor
      Logical FEMtype, operator
C ======================================================================
      DATA XYL/4 * T1A, 
c ... 4 points (order 2)
     &         T2B,T2A,T2A,T2A,  T2A,T2B,T2A,T2A,
     &         T2A,T2A,T2B,T2A,  T2A,T2A,T2A,T2B,
c ... 8 points (order 3)
     &     T3B,T3A,T3A,T3A,  T3A,T3B,T3A,T3A,  T3A,T3A,T3B,T3A,  
     &     T3A,T3A,T3A,T3B,  T3D,T3C,T3C,T3C,  T3C,T3D,T3C,T3C,
     &         T3C,T3C,T3D,T3C,  T3C,T3C,T3C,T3D,
c ... 15 points (order 5)
     &         4 * T5A,
     &     T5C,T5B,T5B,T5B,  T5B,T5C,T5B,T5B,  T5B,T5B,T5C,T5B,  
     &     T5B,T5B,T5B,T5C,  T5E,T5D,T5D,T5D,  T5D,T5E,T5D,T5D,
     &     T5D,T5D,T5E,T5D,  T5D,T5D,T5D,T5E,  T5F,T5F,T5G,T5G,  
     &     T5G,T5G,T5F,T5F,  T5F,T5G,T5G,T5F,  T5G,T5F,T5F,T5G,
     &     T5F,T5G,T5F,T5G,  T5G,T5F,T5G,T5F,
c ... 24 point (order 6)
     &     T6A,T6B,T6B,T6B,  T6B,T6A,T6B,T6B,  T6B,T6B,T6A,T6B,
     &     T6B,T6B,T6B,T6A,  T6C,T6D,T6D,T6D,  T6D,T6C,T6D,T6D,
     &     T6D,T6D,T6C,T6D,  T6D,T6D,T6D,T6C,  T6E,T6F,T6F,T6F,
     &     T6F,T6E,T6F,T6F,  T6F,T6F,T6E,T6F,  T6F,T6F,T6F,T6E,
     &     T6G,T6H,T6I,T6I,  T6G,T6I,T6H,T6I,  T6G,T6I,T6I,T6H,
     &     T6H,T6G,T6I,T6I,  T6I,T6G,T6H,T6I,  T6I,T6G,T6I,T6H,
     &     T6H,T6I,T6G,T6I,  T6I,T6H,T6G,T6I,  T6I,T6I,T6G,T6H,
     &     T6H,T6I,T6I,T6G,  T6I,T6H,T6I,T6G,  T6I,T6I,T6H,T6G/
C ================================================================
      If(order.LE.0 .OR. order.EQ.4 .OR. order.GT.6) 
     &   Call errMesFEM(2001, 'fem3Dtet', 
     &        'input data are incorrect: order = 1,2,3,5, or 6')

      FEMtype = FEMtypeA.EQ.FEMtypeB
      operator = operatorA.EQ.operatorB

c ... transformation of variables y = PSI * (x - x_0)
      Do i = 1, 3
         XYP(i, 1) = 0D0
         XYP(i, 2) = XY2(i) - XY1(i)
         XYP(i, 3) = XY3(i) - XY1(i)
         XYP(i, 4) = XY4(i) - XY1(i)
      End do

      Call solve3x3(XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(1, 1),
     &              XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(1, 2),
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(1, 3), cdet)

      Call solve3x3(XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(2, 2),
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(2, 3), det)

      Call solve3x3(XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(3, 1), 
     &              XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(3, 2),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(3, 3), det)

c ... weights and points
      vol = dabs(det) / 6
      Call WeightsPoints(XY1, XY2, XY3, XY4, vol, order, XYL,
     &                   XYG, w, iGauss, iL)

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)
         If(nfa.EQ.0) 
     &   Call applyGRA2(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else If(operatorA.EQ.DUDX) Then
         Call applyDUDX(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)

      Else If(operatorA.EQ.DIV) Then
         Call applyDIV( iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)
      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)
      Else If(operatorA.EQ.CURL) Then
         Call applyCURL(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else
         Call errMesFEM(2001, 'fem3Dtet', 'operatorA is not supported')
      End if

c ... error messages
      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem3Dtet',
     &     'the local matrix leading dimension is too small')

      If(idim.LE.0) Call errMesFEM(2001, 'fem3Dtet', 
     &     'operatorA is not supported for this finite element')


c ... compute operatorB * FEMtypeB
      nfb = nfa
      jdim = idim
      if(operator .AND. FEMtype) Goto 100

      If(operatorB.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeB, nfb, jdim, V)
         If(nfb.EQ.0) 
     &   Call applyGRA2(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, cdet)

      Else If(operatorB.EQ.DUDX) Then
         Call applyDUDX(iGauss, XYL(1, iL), PSI, FEMtypeB, nfb, jdim, V)

      Else If(operatorB.EQ.DIV) Then
         Call applyDIV( iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, cdet)
      Else If(operatorB.EQ.IDEN) Then
         Call applyIDEN(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, cdet)
      Else If(operatorB.EQ.CURL) Then
         Call applyCURL(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, cdet)

      Else
         Call errMesFEM(2001, 'fem3Dtet', 'operatorB is not supported')
      End if

c ... error messages
      If(jdim.LE.0) Call errMesFEM(2001, 'fem3Dtet', 
     &     'operatorB is not supported for this finite element')


c ... compute D * U
 100  iD = jdim
      jD = idim

      tensor = TENSOR_NULL
      Do n = 1, iGauss
         tensor = D(XYG(1,n), XYG(2,n), XYG(3,n), 
     &                        label, DATA, iSYS, Diff(1, 1, n))
         If(tensor.EQ.TENSOR_NULL) Goto 200
      End do


 200  Continue
      If(tensor.EQ.TENSOR_NULL) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 
     &        'fem3Dtet', 'the operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do
      Else If(tensor.EQ.TENSOR_SCALAR) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem3Dtet', 
     &        'the operators A and B are not compatible')

         Do n = 1, iGauss
            s = Diff(1, 1, n) * w(n) 
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * s
               End do
            End do
         End do
      Else If(tensor.EQ.TENSOR_SYMMETRIC .OR. 
     &        tensor.EQ.TENSOR_GENERAL) Then
         iD = iSYS(1)
         jD = iSYS(2)

         If(jD.NE.idim .OR. iD.NE.jdim) Call errMesFEM(2001, 
     &        'fem3Dtet', 'the operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, iD
               Do k = 1, nfa
                  s = 0D0
                  Do j = 1, jD
                     s = s + Diff(i, j, n) * U(j, k, n)
                  End do
                  DU(i, k, n) = s * w(n)
               End do
            End do
         End do
      Else
         Call errMesFEM(2001, 'fem3Dtet', 
     &        'the given tensor is not supported') 
      End if


c ... compute <D U, V>
      If(operator .AND. FEMtype .AND. tensor.NE.TENSOR_GENERAL) Then
         Do i = 1, nfa
            Do j = 1, i - 1
               A(i, j) = A(j, i)
            End do

            Do j = i, nfa
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * U(k, j, n)
                  End do
                  A(i, j) = s
               End do
            End do
         End do
      Else If(operator .AND. FEMtype) Then
         Do i = 1, nfa
            Do j = 1, nfa
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * U(k, j, n)
                  End do
                  A(j, i) = s
               End do
            End do
         End do
      Else 
         Do i = 1, nfa
            Do j = 1, nfb
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * V(k, j, n)
                  End do
                  A(j, i) = s
               End do
            End do
         End do
      End if

      nRow = nfb
      nCol = nfa

      Return
      End



c ======================================================================
      Subroutine WeightsPoints(XY1, XY2, XY3, XY4, vol, order, XYL,
     &                         XYG, w, iGauss, iLref)
c ======================================================================
      include 'fem3Dtet.fd'
c ======================================================================
C The procedure is used for effective computing points for 
C numerical integration, b/c of a symmetry.
C
C Remark: A 1-to-1 coresspondance between XYL and XYG should be hold.
c ======================================================================
      Real*8  XY1(3), XY2(3), XY3(3), XY4(3), vol, XYL(4, *)
      Integer order

      Real*8  XYG(3, *), w(*)
      Integer iGauss, iLref

      Real*8  s1, s2, s3, s4
c ======================================================================
      If(order.EQ.1) Then
         iGauss = LDG1
         iLref = 1
         w(1) = W1A * vol

         Do i = 1, 3
            XYG(i, 1) = T1A * (XY1(i) + XY2(i) + XY3(i) + XY4(i))
         End do
      Else If(order.EQ.2) Then
         iGauss = LDG2
         iLref = LDG1 + 1
         Do i = 1, iGauss
            w(i) = W2A * vol
         End do

         Do i = 1, 3
            XYG(i, 1) = T2B * XY1(i) + T2A * (XY2(i) + XY3(i) + XY4(i)) 
            XYG(i, 2) = T2B * XY2(i) + T2A * (XY1(i) + XY3(i) + XY4(i)) 
            XYG(i, 3) = T2B * XY3(i) + T2A * (XY1(i) + XY2(i) + XY4(i)) 
            XYG(i, 4) = T2B * XY4(i) + T2A * (XY1(i) + XY2(i) + XY3(i)) 
         End do
      Else If(order.EQ.3) Then
         iGauss = LDG3
         iLref = LDG1 + LDG2 + 1
         Do i = 1, 4
            w(i) = W3A * vol
            w(i + 4) = W3B * vol
         End do

         Do i = 1, 3
            XYG(i, 1) = XY1(i)
            XYG(i, 2) = XY2(i)
            XYG(i, 3) = XY3(i)
            XYG(i, 4) = XY4(i)

            XYG(i, 5) = T3D * XY1(i) + T3C * (XY2(i) + XY3(i) + XY4(i))
            XYG(i, 6) = T3D * XY2(i) + T3C * (XY1(i) + XY3(i) + XY4(i)) 
            XYG(i, 7) = T3D * XY3(i) + T3C * (XY1(i) + XY2(i) + XY4(i)) 
            XYG(i, 8) = T3D * XY4(i) + T3C * (XY1(i) + XY2(i) + XY3(i)) 
         End do
      Else If(order.EQ.5) Then
         iGauss = LDG5
         iLref = LDG1 + LDG2 + LDG3 + 1
         w(1) = W5A * vol
         Do i = 1, 4
            w(i + 1) = W5B * vol
            w(i + 5) = W5C * vol
         End do
         Do i = 1, 6
            w(i + 9) = W5D * vol
         End do

         Do i = 1, 3
            XYG(i, 1) = T5A * (XY1(i) + XY2(i) + XY3(i) + XY4(i))

            XYG(i, 2) = T5C * XY1(i) + T5B * (XY2(i) + XY3(i) + XY4(i)) 
            XYG(i, 3) = T5C * XY2(i) + T5B * (XY1(i) + XY3(i) + XY4(i)) 
            XYG(i, 4) = T5C * XY3(i) + T5B * (XY1(i) + XY2(i) + XY4(i)) 
            XYG(i, 5) = T5C * XY4(i) + T5B * (XY1(i) + XY2(i) + XY3(i)) 

            XYG(i, 6) = T5E * XY1(i) + T5D * (XY2(i) + XY3(i) + XY4(i)) 
            XYG(i, 7) = T5E * XY2(i) + T5D * (XY1(i) + XY3(i) + XY4(i)) 
            XYG(i, 8) = T5E * XY3(i) + T5D * (XY1(i) + XY2(i) + XY4(i)) 
            XYG(i, 9) = T5E * XY4(i) + T5D * (XY1(i) + XY2(i) + XY3(i)) 

            XYG(i, 10) = T5F*(XY1(i) + XY2(i)) + T5G*(XY3(i) + XY4(i)) 
            XYG(i, 11) = T5G*(XY1(i) + XY2(i)) + T5F*(XY3(i) + XY4(i)) 
            XYG(i, 12) = T5F*(XY1(i) + XY4(i)) + T5G*(XY2(i) + XY3(i)) 
            XYG(i, 13) = T5G*(XY1(i) + XY4(i)) + T5F*(XY2(i) + XY3(i)) 
            XYG(i, 14) = T5F*(XY1(i) + XY3(i)) + T5G*(XY2(i) + XY4(i)) 
            XYG(i, 15) = T5G*(XY1(i) + XY3(i)) + T5F*(XY2(i) + XY4(i)) 
         End do
      Else If(order.EQ.6) Then
         iGauss = LDG6
         iLref = LDG1 + LDG2 + LDG3 + LDG5 + 1
         Do i = 1, 4
            w(i)     = W6A * vol
            w(i + 4) = W6B * vol
            w(i + 8) = W6C * vol
         End do
         Do i = 13, 24
            w(i) = W6D * vol
         End do

         Do n = 1, LDG6
            s1 = XYL(1, iLref + n - 1) 
            s2 = XYL(2, iLref + n - 1) 
            s3 = XYL(3, iLref + n - 1) 
            s4 = XYL(4, iLref + n - 1) 

            Do i = 1, 3
               XYG(i, n) = s1 * XY1(i) + s2 * XY2(i) 
     &                   + s3 * XY3(i) + s4 * XY4(i) 
            End do
         End do
      End if
      Return
      End



c ======================================================================
      Real*8 Function Lfun(i, x, y, z)
c ======================================================================
      Real*8 x, y, z

      If(i.EQ.1) Then
         Lfun = 1D0 - x - y - z
      Else If(i.EQ.2) Then
         Lfun = x
      Else If(i.EQ.3) Then
         Lfun = y
      Else If(i.EQ.4) Then
         Lfun = z
      End if
      Return
      End


c ======================================================================
      Subroutine solve3x3(
c ======================================================================
     &     a11, a12, a13, a,
     &     a21, a22, a23, b,
     &     a31, a32, a33, c, det)
c ======================================================================
      Real*8  a11, a12, a13, a
      Real*8  a21, a22, a23, b
      Real*8  a31, a32, a33, c, det

c Local variables
      Real*8  da, db, dc, s
c ======================================================================
      da =  a22 * a33 - a32 * a23
      db =-(a21 * a33 - a31 * a23) 
      dc =  a21 * a32 - a31 * a22

      det = a11 * da + a12 * db + a13 * dc
      s = 1 / det

      a = da * s
      b = db * s
      c = dc * s
      Return
      End



C ======================================================================
      Subroutine calNormalVec(xy1, xy2, xy3, xyn)
C ======================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xyn(3)

C Local variables
      Real*8 ax, ay, az, bx, by, bz

      ax = xy2(1) - xy1(1)
      ay = xy2(2) - xy1(2)
      az = xy2(3) - xy1(3)

      bx = xy3(1) - xy1(1)
      by = xy3(2) - xy1(2)
      bz = xy3(3) - xy1(3)

      xyn(1) = ay * bz - az * by
      xyn(2) = az * bx - ax * bz 
      xyn(3) = ax * by - ay * bx
     
      Return
      End



C ======================================================================
      Subroutine copyGauss(iGauss, nfa, idim, U)
C ======================================================================
      include 'fem3Dtet.fd'
c ======================================================================
      Real*8  U(9, MaxSize, *)

      Do n = 2, iGauss
         Do i = 1, nfa
            Do k = 1, idim
               U(k, i, n) = U(k, i, 1)
            End do
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine clearU(iGauss, nfa, idim, U)
C ======================================================================
      include 'fem3Dtet.fd'
c ======================================================================
      Real*8  U(9, MaxSize, *)

      Do n = 1, iGauss
         Do i = 1, nfa
            Do k = 1, idim
               U(k, i, n) = 0D0
            End do
         End do
      End do

      Return
      End

