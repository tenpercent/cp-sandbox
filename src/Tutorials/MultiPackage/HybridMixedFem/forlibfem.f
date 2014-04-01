C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'
C ======================================================================
C  Diffusion tensor is passed through DATA array. DATA is filled 
C  in main.f using subroutine Drotate at the end of this file.
C ======================================================================
      Real*8  x, y, z, DATA(3,*), Coef(9, 9)
      Integer label, iSYS(*)
C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 3

c HMFEM needs inverted diffusion tensor; we use Lapack for this
      Do i = 1, 3
         Do j = i, 3
            Coef(i, j) = DATA(i, j)
         End do
      End do

      Call invertSPDmatrix(3, Coef, 9)

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Boundary condition
C ======================================================================
      Integer Function Dbc(x, y, z, label, DATA, iSYS, eBC)
      Include 'assemble.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), eBC(*)
      Integer label, iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      If (label.ge.1.and.label.le.6) then
          eBC(1) = 0D0
          Dbc = BC_DIRICHLET
      Else if (label.ge.7.and.label.le.12) then
          eBC(1) = 2D0
          Dbc = BC_DIRICHLET
      Else
          stop 'wrong labelF'
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, z, label, DATA, iSYS, F)
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), F(*)
      Integer label, iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      F(1) = 0D0
      Drhs = TENSOR_SCALAR

      Return
      End



C ======================================================================
c Templated routine for elemental matrix
C ======================================================================
      Subroutine FEM3Dext(
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
      Implicit none
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)

      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DATA(*)
      Integer iSYS(*)

      Integer LDA, nRow, nCol
      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc
      Real*8   calSqr  
      External calSqr  

      Real*8   Q(4, 4), B(4), L(4), G(1), det, alpha, s
      Real*8   x, y, z, eBC(1)
      Integer  i, j, k, ir,ic, INFO
C ======================================================================
      nRow = 4
      nCol = 4

c ... set up templates (unknowns are located on faces)
      Do i = 1, 4
         templateR(i) = Fdof
         templateC(i) = Fdof
      End do

c ... compute the mass matrix (Q)
c     pass a copy of iSYS which is changed inside fem3Dtet
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_RT0, IDEN, FEM_RT0,
     &              lbE, Ddiff, DATA, iSYS, 2,
     &              4, Q, ir, ic)

c ... compute right hand side
c     pass a copy of iSYS which is changed inside fem3Dtet
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              lbE, Drhs, DATA, iSYS, 2,
     &              1, G, ir, ic)

c ... compute the constraint matrix
      B(1) = calSqr(XY1, XY2, XY3)
      B(2) = calSqr(XY2, XY3, XY4)
      B(3) = calSqr(XY3, XY4, XY1)
      B(4) = calSqr(XY4, XY1, XY2)

c ... elliminate fluxes: invert the mass matrix first
      call invertSPDmatrix(4,Q,4)

      Do i = 1, 4
         Do j = 1, 4
            Q(i, j) = Q(i, j) * B(j)
         End do
      End do

      Do i = 1, 4
         L(i) = Q(1, i) * B(1) + Q(2, i) * B(2) 
     &        + Q(3, i) * B(3) + Q(4, i) * B(4)
      End do

c ... elliminate central pressure
      alpha = L(1) + L(2) + L(3) + L(4)

      Do i = 1, 4
         Do j = i, 4
            A(i, j) = B(i) * Q(i, j) - L(i) * L(j) / alpha
            A(j, i) = A(i, j)
         End do
         F(i) = L(i) * G(1) / alpha
      End do

c ... impose boundary conditions (assume nRow = nCol)
      Do k = 1, 4
         If(lbF(k).GT.0) Then
            If(k.EQ.1) Then
               x = (XY1(1) + XY2(1) + XY3(1)) / 3
               y = (XY1(2) + XY2(2) + XY3(2)) / 3
               z = (XY1(3) + XY2(3) + XY3(3)) / 3
            ElseIf(k.EQ.2) Then
               x = (XY2(1) + XY3(1) + XY4(1)) / 3
               y = (XY2(2) + XY3(2) + XY4(2)) / 3
               z = (XY2(3) + XY3(3) + XY4(3)) / 3
            ElseIf(k.EQ.3) Then
               x = (XY3(1) + XY4(1) + XY1(1)) / 3
               y = (XY3(2) + XY4(2) + XY1(2)) / 3
               z = (XY3(3) + XY4(3) + XY1(3)) / 3
            ElseIf(k.EQ.4) Then
               x = (XY4(1) + XY1(1) + XY2(1)) / 3
               y = (XY4(2) + XY1(2) + XY2(2)) / 3
               z = (XY4(3) + XY1(3) + XY2(3)) / 3
            End if

c           pass of copy of iSYS which is changed inside Dbc
            i = Dbc(x, y, z, lbF(k), DATA, iSYS, eBC)

            Do i = 1, nRow
               F(i) = F(i) - A(i, k) * eBC(1)
               A(k, i) = 0
               A(i, k) = 0
            End do

            A(k, k) = 1
            F(k) = eBC(1)
         End if
      End do

      Return
      End



C ======================================================================
C  The routine takes an input matrix D and rotates it about i-th
C  coordinate axis by angle theta.
C ======================================================================
      Subroutine Drotate(D, i, theta)
C ======================================================================
      Real*8  D(3, 3), theta, c, s, a, b
C ======================================================================
      j = i + 1
      If( j.GT.3 ) j = 1

      c = dcos(theta)
      s = dsin(theta)

      Do k = 1, 3
         a = D(k, i)
         b = D(k, j)

         D(k, i) = a * c - b * s
         D(k, j) = a * s + b * c
      End do

      Do k = 1, 3
         a = D(i, k)
         b = D(j, k)

         D(i, k) = a * c - b * s
         D(j, k) = a * s + b * c
      End do

      Return
      End

