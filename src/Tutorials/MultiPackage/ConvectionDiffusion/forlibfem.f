C ================================================================
      Subroutine FEM3Dext(
C ================================================================
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'
C ================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)

      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DATA(*)
      Integer iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc
      External Ddiff, Dconv, Drhs, Dbc

      Real*8   B(4, 4)
      Real*8   x, y, z, eBC(1)

C ================================================================
      nRow = 4
      nCol = 4

c ... set up templates
      Do i = 1, nRow
         templateR(i) = Vdof   
      End do

      Do k = 1, nCol
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix 
      label = lbE

      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DATA, iSYS, 1,
     &              LDA, A, ir, ic)

      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P1,
     &              label, Dconv, DATA, iSYS, 1,
     &              4, B, ir, ic)

      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do

c ... compute right hand side
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, DATA, iSYS, 1,
     &              LDA, F, ir, ic)

c ... impose boundary conditions (assume nRow = nCol) at nodes of triangle
      Do k = 1, 4
         If(lbP(k).ne.0) Then
            If(k.EQ.1) Then
               x = XY1(1)
               y = XY1(2)
               z = XY1(3)
            ElseIf(k.EQ.2) Then
               x = XY2(1)
               y = XY2(2)
               z = XY2(3)
            ElseIf(k.EQ.3) Then
               x = XY3(1)
               y = XY3(2)
               z = XY3(3)
            ElseIf(k.EQ.4) Then
               x = XY4(1)
               y = XY4(2)
               z = XY4(3)
            End if

c  ...  Dbc may change iSYS
            ibc = Dbc(x, y, z, lbP(k), DATA, iSYS, eBC)

            If(ibc.EQ.BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End



C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer label, iSYS(*)
C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = 1D-2
      End do

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C  Convection tensor             
C ======================================================================
      Integer Function Dconv(x, y, z, label, DATA, iSYS, Conv)
      include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Conv(9, 9)
      Integer label, iSYS(*)
C ================================================================
      iSYS(1) = 1
      iSYS(2) = 3

      Conv(1,1) = DATA(1)
      Conv(1,2) = DATA(2)
      Conv(1,3) = DATA(3)

      Dconv = TENSOR_GENERAL  

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

      Dbc = BC_DIRICHLET
      eBC(1)  = 0d0

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

      F(1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


