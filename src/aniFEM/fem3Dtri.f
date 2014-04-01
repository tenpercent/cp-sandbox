C ======================================================================
      Subroutine fem3Dtri(
C ======================================================================
     &      XY1D3, XY2D3, XY3D3, 
     &      operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &      label, D, DATAFEM, iSYS, order,  
     &      LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
      include 'fem3Dtri.fd'
      include 'assemble.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes elemental matrix for a bilinear form 
C
C  (1)            <D OpA(u), OpB(v)>               
C
C  where D is a tensor, OpA and OpB are linear operators, and u and v 
C  are finite element functions, u in "A", and v in "B". The matrix is
C  computed for a triangular face of a tetrahedral element. 
C ======================================================================
      Real*8   XY1D3(*), XY2D3(*), XY3D3(*)

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  label, order, LDA, D, nRow, nCol

      Real*8   DATAFEM(*)
      Integer  iSYS(*)
      EXTERNAL D

      Real*8   A(LDA, *)

C ======================================================================
c Local variables
      Real*8  XY1(2), XY2(2), XY3(2), Rm(3, 3)

      Real*8  XYP(2, 3), XYN(2, 3)
      Real*8  det, vol, s
   
      Real*8  PSI(2, 2)
      Real*8  U(4, 12, MaxPnt2DGauss),    V(4, 12, MaxPnt2DGauss)
      Real*8  Diff(4, 4, MaxPnt2DGauss), DU(4, 12, MaxPnt2DGauss)

      Real*8  w(MaxPnt2DGauss), XYG(3, MaxPnt2DGauss)
      Real*8  XYL(3, AllPnt2DGauss)

      Integer i,j,k,n, nfa,nfb, idim,jdim, iD,jD, tensor, iGauss, iL
      Logical FEMtype, operator

C ======================================================================
      DATA XYL/3 * U1A, 
c ... 3 points (order 2)
     &         U2A,U2B,U2B,  U2B,U2A,U2B,  U2B,U2B,U2A,
c ... 7 points (order 5)
     &         U5A,U5A,U5A,  U5B,U5C,U5C,  U5C,U5B,U5C,  U5C,U5C,U5B,
     &                       U5D,U5E,U5E,  U5E,U5D,U5E,  U5E,U5E,U5D,
c ... 12 points (order 6)
     &         U6A,U6B,U6B,  U6B,U6A,U6B,  U6B,U6B,U6A,  
     &         U6C,U6D,U6D,  U6D,U6C,U6D,  U6D,U6D,U6C,  
     &         U6E,U6F,U6G,  U6G,U6E,U6F,  U6F,U6G,U6E,  
     &         U6E,U6G,U6F,  U6G,U6F,U6E,  U6F,U6E,U6G/

C ================================================================
      Call TransformMatrix(XY1D3, XY2D3, XY3D3, Rm)

      Do i = 1, 2
         XY1(i) = 0D0
         XY2(i) = 0D0
         XY3(i) = 0D0
      End do

      Do i = 1, 2
         Do j = 1, 3
            XY1(i) = XY1(i) + XY1D3(j) * Rm(j, i)
            XY2(i) = XY2(i) + XY2D3(j) * Rm(j, i)
            XY3(i) = XY3(i) + XY3D3(j) * Rm(j, i)
         End do
      End do

      If(order.LE.0 .OR. order.GT.6) 
     &  Call errMesFEM(2001, 'fem3Dtri', 'quadrature order is wrong')

      FEMtype = FEMtypeA.EQ.FEMtypeB
      operator = operatorA.EQ.operatorB

c ... transformation of variables y = PSI * (x - x_0)
      Do i = 1, 2
         XYP(i, 1) = 0D0
         XYP(i, 2) = XY2(i) - XY1(i)
         XYP(i, 3) = XY3(i) - XY1(i)
      End do

      Call solve2x2(XYP(1, 2), XYP(2, 2), PSI(1, 1),
     &              XYP(1, 3), XYP(2, 3), PSI(1, 2), det)

      Call solve2x2(XYP(1, 3), XYP(2, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), PSI(2, 2), det)

c ... weights and points
      vol = dabs(det) / 2
      Call WeightsPnt2D(XY1D3, XY2D3, XY3D3, vol, order, 
     &                  XYG, w, iGauss, iL)

c ... exterior normal vectors
      If(FEMtypeA.EQ.FEM_P2reduced .OR. FEMtypeB.EQ.FEM_P2reduced) Then
         Call calNormalExt(XY1, XY2, XY3, XYN(1, 1))
         Call calNormalExt(XY2, XY3, XY1, XYN(1, 2))
         Call calNormalExt(XY3, XY1, XY2, XYN(1, 3))
      End if

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD2D(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                    nfa, idim, U, XYN)

      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN2D(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                    nfa, idim, U, XYP, XYN, det)

      Else
         Call errMesFEM(2001, 'fem3Dtri', 'operatorA is not supported')
      End if

      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem3Dtri',
     &     'the local matrix leading dimension, LDA, is too small')


c ... compute operatorB * FEMtypeB
      nfb = nfa
      jdim = idim
      if(operator .AND. FEMtype) Goto 100

      If(operatorB.EQ.GRAD) Then
         Call applyGRAD2D(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                    nfb, jdim, V, XYN)

      Else If(operatorB.EQ.IDEN) Then
         Call applyIDEN2D(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                    nfb, jdim, V, XYP, XYN, det)

      Else
         Call errMesFEM(2001, 'fem3Dtri', 'operatorB is not supported')
      End if

      If(nfb.GT.LDA) Call errMesFEM(2001, 'fem3Dtri',
     &     'the local matrix second dimension, LDA, is too small')
 

c ... compute D * U
 100  iD = jdim
      jD = idim

      Do n = 1, iGauss
         tensor = D(XYG(1,n), XYG(2,n), XYG(3,n), label, 
     &              DATAFEM, iSYS, Diff(1, 1, n))
         If(tensor.EQ.TENSOR_NULL) Goto 200
      End do
      If(tensor.GE.BC_DIRICHLET) tensor = TENSOR_SCALAR


 200  Continue
      If(tensor.EQ.TENSOR_NULL) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 
     &        'fem3Dtri', 'Operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do
      Else If(tensor.EQ.TENSOR_SCALAR) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem3Dtri', 
     &        'Operators A and B are not compatible')

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
     &        'fem3Dtri', 'the operators A and B are not compatible')

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
         Call errMesFEM(2001, 'fem3Dtri', 
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
            Do j = 1, nfb
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
      Subroutine TransformMatrix(x, y, z, Rm)
c ======================================================================
      Implicit none
      Real*8   x(3), y(3), z(3), Rm(3,3)

c Local variables
      Real*8   a(3), b(3), dnrm, dsqrt
      Integer i

c ======================================================================
      Do i = 1, 3
         a(i) = y(i) - x(i)
         b(i) = z(i) - x(i)
      End do

      dnrm = dsqrt(a(1)**2 + a(2)**2 + a(3)**2)

      Do i = 1, 3
         a(i) = a(i) / dnrm
      End do

      dnrm = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

      Do i = 1, 3
         b(i) = b(i) - a(i) * dnrm
      End do

      dnrm = dsqrt(b(1)**2 + b(2)**2 + b(3)**2)

      Do i = 1, 3
         Rm(i, 1) = a(i)
         Rm(i, 2) = b(i) / dnrm
      End do

c ... debug starts
c     Rm(1,3) = Rm(2,1) * Rm(3,2) - Rm(2,2) * Rm(3,1)
c     Rm(2,3) = Rm(1,2) * Rm(3,1) - Rm(1,1) * Rm(3,2)
c     Rm(3,3) = Rm(1,1) * Rm(2,2) - Rm(1,2) * Rm(2,1)

      Return
      End



c ======================================================================
      Subroutine WeightsPnt2D(XY1, XY2, XY3, vol, order, 
     &                        XYG, w, iGauss, iLref)
c ======================================================================
      include 'fem3Dtri.fd'
c ======================================================================
C The procedure is used for effective computing points for numerical
C integration over triangle in 3D space.
C
C Remark: A 1-to-1 coresspondance between XYL and XYG should be hold.
c ======================================================================
      Real*8  XY1(3), XY2(3), XY3(3), vol
      Integer order

      Real*8  XYG(3, *), w(*)
      Integer iGauss, iLref

c ======================================================================
      If(order.EQ.1) Then
         iGauss = KDG1
         iLref = 1
         w(1) = S1A * vol

         Do i = 1, 3
            XYG(i, 1) = U1A * (XY1(i) + XY2(i) + XY3(i))
         End do

      Else If(order.EQ.2) Then
         iGauss = KDG2
         iLref = KDG1 + 1
         Do i = 1, iGauss
            w(i) = S2A * vol
         End do

         Do i = 1, 3
            XYG(i, 1) = U2B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = U2B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = U2B * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.5) Then
         iGauss = KDG5
         iLref = KDG1 + KDG2 + 1

         w(1) = S5A * vol
         Do i = 1, 3
            XYG(i, 1) = U5A * (XY1(i) + XY2(i) + XY3(i))
         End do

         Do i = 2, 4
            w(i) = S5B * vol
         End do
         Do i = 1, 3
            XYG(i, 2) = U5B * XY1(i) + U5C * (XY2(i) + XY3(i)) 
            XYG(i, 3) = U5B * XY2(i) + U5C * (XY1(i) + XY3(i)) 
            XYG(i, 4) = U5B * XY3(i) + U5C * (XY1(i) + XY2(i)) 
         End do

         Do i = 5, 7
            w(i) = S5C * vol
         End do
         Do i = 1, 3
            XYG(i, 5) = U5D * XY1(i) + U5E * (XY2(i) + XY3(i)) 
            XYG(i, 6) = U5D * XY2(i) + U5E * (XY1(i) + XY3(i)) 
            XYG(i, 7) = U5D * XY3(i) + U5E * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.6) Then
         iGauss = KDG6
         iLref = KDG1 + KDG2 + KDG5 + 1

         Do i = 1, 3
            w(i) = S6A * vol
         End do
         Do i = 1, 3
            XYG(i, 1) = U6A * XY1(i) + U6B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = U6A * XY2(i) + U6B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = U6A * XY3(i) + U6B * (XY1(i) + XY2(i)) 
         End do

         Do i = 4, 6
            w(i) = S6B * vol
         End do
         Do i = 1, 3
            XYG(i, 4) = U6C * XY1(i) + U6D * (XY2(i) + XY3(i)) 
            XYG(i, 5) = U6C * XY2(i) + U6D * (XY1(i) + XY3(i)) 
            XYG(i, 6) = U6C * XY3(i) + U6D * (XY1(i) + XY2(i)) 
         End do

         Do i = 7, 12
            w(i) = S6C * vol
         End do
         Do i = 1, 3
            XYG(i, 7) = U6E * XY1(i) + U6F * XY2(i) + U6G * XY3(i)
            XYG(i, 8) = U6G * XY1(i) + U6E * XY2(i) + U6F * XY3(i)
            XYG(i, 9) = U6F * XY1(i) + U6G * XY2(i) + U6E * XY3(i)

            XYG(i,10) = U6E * XY1(i) + U6G * XY2(i) + U6F * XY3(i)
            XYG(i,11) = U6G * XY1(i) + U6F * XY2(i) + U6E * XY3(i)
            XYG(i,12) = U6F * XY1(i) + U6E * XY2(i) + U6G * XY3(i)
         End do
      End if
      Return
      End



c ======================================================================
      Real*8 Function Lfun2D(i, x, y)
c ======================================================================
      Integer  i
      Real*8   x, y

      If(i.EQ.1) Then
         Lfun2D = 1D0 - x - y
      Else If(i.EQ.2) Then
         Lfun2D = x
      Else If(i.EQ.3) Then
         Lfun2D = y
      End if
      Return
      End



c ======================================================================
      Subroutine solve2x2(
c ======================================================================
     &     a11, a12, a,
     &     a21, a22, b, det)
c ======================================================================
      Real*8  a11, a12, a
      Real*8  a21, a22, b, det

c Local variables
      Real*8  s

c ======================================================================
      det = a11 * a22 - a21 * a12
      s = 1 / det

      a = a22 * s
      b =-a21 * s
      Return
      End



c ======================================================================
      Subroutine invert3x3(A, B, det)
c ======================================================================
      Real*8  A(3, 3), B(3, 3), det

c Local variables
      Real*8  s11, s12, s13

c ======================================================================
      s11 = A(2,2) * A(3,3) - A(3,2) * A(2,3)
      s12 = A(2,1) * A(3,3) - A(3,1) * A(2,3) 
      s13 = A(2,1) * A(3,2) - A(3,1) * A(2,2)

      det = A(1,1) * s11 - A(1,2) * s12 + A(1,3) * s13

      B(1,1) = s11 / det
      B(2,1) =-s12 / det
      B(3,1) = s13 / det

      B(1,2) =-(A(1,2) * A(3,3) - A(3,2) * A(1,3)) / det
      B(2,2) = (A(1,1) * A(3,3) - A(3,1) * A(1,3)) / det
      B(3,2) =-(A(1,1) * A(3,2) - A(3,1) * A(1,2)) / det

      B(1,3) = (A(1,2) * A(2,3) - A(2,2) * A(1,3)) / det
      B(2,3) =-(A(1,1) * A(2,3) - A(2,1) * A(1,3)) / det
      B(3,3) = (A(1,1) * A(2,2) - A(2,1) * A(1,2)) / det

      Return
      End



C ======================================================================
      Subroutine calNormalExt(xy1, xy2, xy3, xyn)
C ======================================================================
C Routines compute external normal vectors to the edge {xy1, xy2}
C of triangle {xy1, xy2, xy3}. This is the twin sister of a similar
C routine from package aniMBA.
C ======================================================================
      Real*8 xy1(2), xy2(2), xy3(2), xyn(2)
      Real*8 x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = dsqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

c ... orientation
      x = xy3(1) - xy1(1)
      y = xy3(2) - xy1(2)

      If( x*xyn(1) + y*xyn(2).GT.0D0 ) Then
         xyn(1) = -xyn(1)
         xyn(2) = -xyn(2)
      End if 

      Return
      End



C ======================================================================
      Subroutine copyGauss2D(iGauss, nfa, idim, U)
C ======================================================================
      Real*8  U(4, 12, *)

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
      Subroutine clearU2D(iGauss, nfa, idim, U)
C ======================================================================
      Real*8  U(4, 12, *)

      Do n = 1, iGauss
         Do i = 1, nfa
            Do k = 1, idim
               U(k, i, n) = 0D0
            End do
         End do
      End do

      Return
      End



c ======================================================================
      Subroutine applyGRAD2D(iGauss, XYL, PSI, FEMtype,
     &                       nfa, dim, U, XYN)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
      include 'fem3Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), U(4, 12, *), XYN(2, 3)

c Data for the reference triangle
      Real*8  XYL(3, *)
      Real*8  GRAD_P1(2, 3), GRAD_P2(2,  6, MaxPnt2DGauss), GRAD_P2r(2)
      Real*8                 GRAD_P3(2, 10, MaxPnt2DGauss), GRAD_P3r(2)

      Integer i, j, k, n, l, mfa
      Real*8  x, y, Lfun2D, s1, s2, s3

      Integer iref(4)

      DATA    GRAD_P1/-1,-1, 1,0, 0,1/
      DATA    iref /1,2,3,1/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 0D0
         End do

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P1 .OR. FEMtype.EQ.FEM_P1DG) Then
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            Do k = 1, dim
               U(k, i, 1) = 0D0
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
            End do
         End do
         Call copyGauss2D(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
c  ...  use formula psi_l = 1 - 2 phi_i
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            l = iref(i + 1)
            Do k = 1, dim
               U(k, l, 1) = 0D0
               Do j = 1, 2
                  U(k, l, 1) = U(k, l, 1) - 2*PSI(j, k) * GRAD_P1(j, i)
               End do
            End do
         End do
         Call copyGauss2D(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         nfa = 6
         dim = 4
         Call clearU2D(1, nfa, dim, U)

         Do i = 1, 3
            l = iref(i + 1)
            Do k = 1, 2
               U(k, l, 1) = 0D0
               Do j = 1, 2
                  U(k, l, 1) = U(k, l, 1) - 2*PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, l + 3, 1) = U(k, l, 1)
            End do
         End do

         Call copyGauss2D(iGauss, nfa, dim, U)

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 2
         If(FEMtype.EQ.FEM_P2vector)  Call clearU2D(iGauss, 12, 4, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               Do k = 1, dim
                  GRAD_P2(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun2D(i, x, y) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, dim
                  GRAD_P2(k, mfa, n) = 
     &               4 * (Lfun2D(i, x, y) * GRAD_P1(k, j) +
     &                    Lfun2D(j, x, y) * GRAD_P1(k, i))
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 2
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, k) * GRAD_P2(j, i, n)
                  End do
               End do
            End do
         End do

         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 12
            dim = 4

            Do n = 1, iGauss
               Do i = 1, 6
                  Do k = 1, 2
                     U(k + 2, i + 6, n) = U(k, i, n)
                  End do
               End do
            End do
         End if

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 4 
         Call clearU2D(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 3, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss2D(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         nfa = 8
         dim = 4 
         Call clearU2D(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 4, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss2D(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            s1 = Lfun2D(1, x, y)
            s2 = Lfun2D(2, x, y)
            s3 = Lfun2D(3, x, y)

            Do k = 1, 2
               GRAD_P3r(k) = 27 * (GRAD_P1(k, 1) * s2 * s3 +
     &                             GRAD_P1(k, 2) * s1 * s3 +
     &                             GRAD_P1(k, 3) * s1 * s2)
            End do

            Do k = 1, 2
               U(k, 4, n) = PSI(1, k) * GRAD_P3r(1) 
     &                    + PSI(2, k) * GRAD_P3r(2)
               U(k + 2, 8, n) = U(k, 4, n)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2reduced) Then
         nfa = 9
         dim = 4 
         Call clearU2D(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 6, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss2D(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               Do k = 1, 2
                  GRAD_P2r(k) = 4 * (Lfun2D(i, x, y) * GRAD_P1(k, j) +
     &                               Lfun2D(j, x, y) * GRAD_P1(k, i))
               End do

               Do k = 1, 2
                  s1 = PSI(1, k) * GRAD_P2r(1) + PSI(2, k) * GRAD_P2r(2)
                  U(k,     i + 3, n) = XYN(1, i) * s1
                  U(k + 2, i + 3, n) = XYN(2, i) * s1
               End do
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10
         dim = 2

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               s1 = Lfun2D(i, x, y)
               s2 = Lfun2D(j, x, y)

               Do k = 1, dim
                  GRAD_P3(k, i, n) = 
     &                GRAD_P1(k, i) * (13.5D0 * s1 * (s1 - 1D0) + 1D0)

                  GRAD_P3(k, i + 3, n) = 
     &                GRAD_P1(k, i) * (6D0 * s1 - 1D0) * s2 +
     &                GRAD_P1(k, j) * (3D0 * s1 - 1D0) * s1

                  GRAD_P3(k, i + 6, n) = 
     &                GRAD_P1(k, j) * (6D0 * s2 - 1D0) * s1 +
     &                GRAD_P1(k, i) * (3D0 * s2 - 1D0) * s2
               End do
            End do

            s1 = Lfun2D(1, x, y)
            s2 = Lfun2D(2, x, y)
            s3 = Lfun2D(3, x, y)

            Do k = 1, dim
               GRAD_P3(k, 10, n) = 27 * 
     &            (GRAD_P1(k, 1) * s2 * s3 +
     &             GRAD_P1(k, 2) * s1 * s3 +
     &             GRAD_P1(k, 3) * s1 * s2)
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 2
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, k) * GRAD_P3(j, i, n)
                  End do
               End do
            End do
         End do

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'applyGRAD2D', 
     &        'unsupported operation for this element type')
      End if
      Return
      End



c ======================================================================
      Subroutine applyIDEN2D(iGauss, XYL, PSI, FEMtype, 
     &                       nfa, dim, U, XYP, XYN, det)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
      include 'fem3Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), XYL(3, *), U(4, 12, *)
      Real*8  XYP(2, *), XYN(2, *), det

C LOCAL VARIABLEs
      Real*8  vol, edge, calEdge
      Real*8  x, y, Lfun2D, s, s1, s2
      Integer i, j, k, l, n, mfa, iP1, iP2, iP3

c Data for the reference triangle
      Integer IPF(3, 3), iref(4)

      DATA    IPF/1,2,3, 2,3,1, 1,3,2/ 
      DATA    iref /1,2,3,1/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 1D0
         End do

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P1 .OR. FEMtype.EQ.FEM_P1DG) Then
         nfa = 3
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = XYL(i, n)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
c  ...  use formula psi_l = 1 - 2 phi_i
         nfa = 3
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               l = iref(i + 1)
               U(1, l, n) = 1 - 2 * XYL(i, n)  
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         nfa = 6
         dim = 2
         Call clearU2D(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            Do i = 1, nfa
               l = iref(i + 1)
               U(1, l,     n) = 1 - 2 * XYL(i, n)  
               U(2, l + 3, n) = 1 - 2 * XYL(i, n)  
            End do
         End do

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 1
         If(FEMtype.EQ.FEM_P2vector) Call clearU2D(iGauss, 12, 2, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               s = Lfun2D(i, x, y)
               U(1, i, n) = (2 * s - 1D0) * s
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               U(1, mfa, n) = 4 * Lfun2D(i, x, y) * Lfun2D(j, x, y)
            End do
         End do

         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 12
            dim = 2

            Do n = 1, iGauss
               Do i = 1, 6
                  U(2, i + 6, n) = U(1, i, n)
               End do
            End do
         End if

c ... next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10 
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               s1 = Lfun2D(i, x, y)
               s2 = Lfun2D(j, x, y)

               U(1, i,     n) = s1 * (3*s1 - 1D0) * (3*s1 - 2D0) / 2
               U(1, i + 3, n) = s1 * (3*s1 - 1D0) * s2 * 4.5D0
               U(1, i + 6, n) = s1 * (3*s2 - 1D0) * s2 * 4.5D0
            End do

            U(1, 10, n) = 27 * Lfun2D(1, x, y) 
     &                       * Lfun2D(2, x, y) * Lfun2D(3, x, y)
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 2

         Call clearU2D(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 3, n) = XYL(i, n)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         nfa = 8
         dim = 2

         Call clearU2D(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 4, n) = XYL(i, n)
            End do
         End do

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            U(1, 4, n) = 27 * Lfun2D(1, x, y) 
     &                      * Lfun2D(2, x, y) * Lfun2D(3, x, y)
            U(2, 8, n) = U(1, 4, n) 
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2reduced) Then
         nfa = 9
         dim = 2

         Call clearU2D(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 6, n) = XYL(i, n)
            End do

            Do i = 1, 3
               j = iref(i + 1)

               s1 = 4 * Lfun2D(i, x, y) * Lfun2D(j, x, y)
               U(1, i + 3, n) = s1 * XYN(1, i)
               U(2, i + 3, n) = s1 * XYN(2, i)
            End do
         End do

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_RT0 .OR. FEMtype.EQ.FEM_BDM1) Then
         nfa = 3
         dim = 2
         vol = dabs(det)

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            iP3 = IPF(3, i)
            
            edge = calEdge(XYP(1, iP1), XYP(1, iP2)) 
            edge = edge / vol
            Do n = 1, iGauss
               Do k = 1, dim
                  U(k, i, n) = edge *
     &                 (XYL(iP1, n) * (XYP(k, iP1) - XYP(k, iP3))
     &                + XYL(iP2, n) * (XYP(k, iP2) - XYP(k, iP3)))
               End do

               If(FEMtype.EQ.FEM_BDM1) Then
                  Do k = 1, dim
                     U(k, i+nfa, n) = edge * 6 *
     &                   (- XYL(iP1, n) * (XYP(k, iP1) - XYP(k, iP3))
     &                    + XYL(iP2, n) * (XYP(k, iP2) - XYP(k, iP3)))
                  End do
               End if
            End do
         End do

         If(FEMtype.EQ.FEM_BDM1) nfa = 6

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'applyIDEN2D', 
     &        'unsupported operation for the given element type')
      End if
      Return
      End

