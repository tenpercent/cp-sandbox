C ======================================================================
c Templated routine for elemental matrix
C ======================================================================
      Subroutine fem3Dext( 
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      
      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DATA(*)
      Integer iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc, ANI_Dnull
      External Ddiff, Dconv, Drhs, Dbc, ANI_Dnull

      Real*8   B(30, 30), XYP(3, 4)
      Real*8   x, y, z, eBC(3)

      Integer  iE, iP1, iP2, iP3, iP4, iR1, iR2, iR3, iR4, iR5, iR6
      Integer  iF1, iF2, iF3, iF4, nP, nR, nF, nE
      Real*8   DATANULL(1), DATACONV(42)

C ======================================================================
      Do i = 1, 3
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
         XYP(i, 4) = XY4(i)
      End do

c ... size of the elemental problem (30 velocities + 4 pressures)
      nRow = 34 
      nCol = 34

c ... set up templates 
      Do i = 1, 4
         templateR(i)      = Vdof   !u_x
         templateR(i + 10) = Vdof   !u_y
         templateR(i + 20) = Vdof   !u_z

         templateR(i + 30) = Vdof   !p
      End do

      Do i = 1, 6
         templateR(i + 4)  = Rdof
         templateR(i + 14) = Rdof
         templateR(i + 24) = Rdof
      End do

      Do k = 1, nRow
         templateC(k) = templateR(k)
      End do

c ... preprocess data arrays
      Call decodeISYS(iE, iP1, iP2, iP3, iP4,
     &                iR1, iR2, iR3, iR4, iR5, iR6, 
     &                iF1, iF2, iF3, iF4,  
     &                nP, nR, nF, nE, iSYS)

      k = 0
      Do i = 1, 3
         Do j = 1, 4
            k = k + 1
            DATACONV(k) = XYP(i, j)
         End do
      End do

      Do i = 1, 3
         iiP = (i - 1) * nP
         iiR = 4 * nP + (i - 1) * nR

         DATACONV(k + 1) = DATA(iiP + iP1)
         DATACONV(k + 2) = DATA(iiP + iP2)
         DATACONV(k + 3) = DATA(iiP + iP3)
         DATACONV(k + 4) = DATA(iiP + iP4)

         DATACONV(k + 5) = DATA(iiR + iR1)
         DATACONV(k + 6) = DATA(iiR + iR2)
         DATACONV(k + 7) = DATA(iiR + iR3)
         DATACONV(k + 8) = DATA(iiR + iR4)
         DATACONV(k + 9) = DATA(iiR + iR5)
         DATACONV(k +10) = DATA(iiR + iR6)

         k = k + 10
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:30,1:30) is local vector Laplacian matrix for three velocity components
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Ddiff, DATANULL, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:30,1:30) is local vector convection matrix for three velocity components
c     convection is passed through DATA(*)
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P2vector, IDEN, FEM_P2vector,
     &              label, Dconv, DATACONV, iSYS, 2,
     &              30, B, ir, ic)

c ... summing diffusion and convection matrices
      Do i = 1, 30
         Do j = 1, 30
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do


c     B(1:4,1:30) is local divergence matrix for velocity
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              DIV, FEM_P2vector, IDEN, FEM_P1,
     &              label, ANI_Dnull, DATANULL, iSYS, 2,
     &              30, B, ir, ic)

      Do i = 1, 30
         Do j = 1, 4
            A(i, j + 30) = -B(j, i)
            A(j + 30, i) = -B(j, i)
         End do
      End do

      Do i = 1, 4
         Do j = 1, 4
            A(i + 30, j + 30) = 0
         End do
      End do

c ... compute the right hand side
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P2vector, 
     &              lbE, Drhs, DATANULL, iSYS, 2,
     &              LDA, F, ir, ic)

      Do i = 1, 4
         F(i + 30) = 0
      End do

c ... impose boundary conditions at vertices of tetrahedron
      Do k = 1, 4
         If(lbP(k).ne.0) Then
            x = XYP(1, k)
            y = XYP(2, k)
            z = XYP(3, k)

c  ...  Dbc will change iSYS(1:2)
            ibc = Dbc(x, y, z, lbP(k), DATA, iSYS, eBC)

            If(ibc.EQ.BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k,    eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+10, eBC(2))
               Call applyDIR(LDA, nRow, A, F, k+20, eBC(3))
            End if
         End if
      End do

c ... impose boundary conditions (assume nRow = nCol) at mid-points of edges
      k = 0
      Do i = 1, 3
         Do j = i + 1, 4  
            k = k + 1

            If(lbP(i).GT.0 .AND. lbP(j).GT.0) Then
               x = (XYP(1, i) + XYP(1, j)) / 2
               y = (XYP(2, i) + XYP(2, j)) / 2
               z = (XYP(3, i) + XYP(3, j)) / 2

c  ...  Dbc will change iSYS(1:2)
               ibc = Dbc(x, y, z, lbP(i), DATA, iSYS, eBC)

               If(ibc.EQ.BC_DIRICHLET) Then
                  Call applyDIR(LDA, nRow, A, F, k+4,  eBC(1))
                  Call applyDIR(LDA, nRow, A, F, k+14, eBC(2))
                  Call applyDIR(LDA, nRow, A, F, k+24, eBC(3))
               End if
            End if
         End do
      End do

      Return
      End



C ======================================================================
C Identity diffusion tensor
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
C ======================================================================
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D-1
      Ddiff = TENSOR_SCALAR
      Return
      End



C ======================================================================
C Convection tensor        
      Integer Function Dconv(x, y, z, label, DATA, iSYS, Conv)
C ======================================================================
      include 'fem3Dtet.fd'

      Real*8   x, y, z, DATA(*), Conv(9, 9)
      Integer  iSYS(*)

      Real*8   calVol
      External calVol

      Real*8   xy(3), xy1(3), xy2(3), xy3(3), xy4(3)
      Real*8   s1,s2,s3,s4,s12,s13,s14,s23,s24,s34
      Real*8   d1,d2,d3,d4,d,lam1,lam2,lam3,lam4 
      Real*8   g12,g13,g14,g23,g24,g34

      Integer  i,j

C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 9

      Do i = 1, 3
         Do j = 1, 9
            Conv(i, j) = 0
         End do
      End do

      xy(1)  = x
      xy(2)  = y
      xy(3)  = z
      xy1(1) = DATA( 1)
      xy2(1) = DATA( 2)
      xy3(1) = DATA( 3)
      xy4(1) = DATA( 4)
      xy1(2) = DATA( 5)
      xy2(2) = DATA( 6)
      xy3(2) = DATA( 7)
      xy4(2) = DATA( 8)
      xy1(3) = DATA( 9)
      xy2(3) = DATA(10)
      xy3(3) = DATA(11)
      xy4(3) = DATA(12)

      d1 = calVol(xy, xy2, xy3, xy4)
      d2 = calVol(xy, xy1, xy4, xy3)
      d3 = calVol(xy, xy1, xy2, xy4)
      d4 = calVol(xy, xy1, xy3, xy2)
      d = d1 + d2 + d3 + d4
c     barycentric coordinates of x,y,z
      lam1 = d1 / d
      lam2 = d2 / d
      lam3 = d3 / d
      lam4 = d4 / d

c quadratic interpolation of ux
      s1 = DATA(13)
      s2 = DATA(14)
      s3 = DATA(15)
      s4 = DATA(16)
      s12= DATA(17)
      s13= DATA(18)
      s14= DATA(19)
      s23= DATA(20)
      s24= DATA(21)
      s34= DATA(22)

      g12 = 4 * ( s12 - (s1+s2)/2 )
      g13 = 4 * ( s13 - (s1+s3)/2 )
      g14 = 4 * ( s14 - (s1+s4)/2 )
      g23 = 4 * ( s23 - (s2+s3)/2 )
      g24 = 4 * ( s24 - (s2+s4)/2 )
      g34 = 4 * ( s34 - (s3+s4)/2 )

      Conv(1,1) = s1*lam1 + s2*lam2 + s3*lam3 + s4*lam4
     &          + g12*lam1*lam2 + g13*lam1*lam3 + g14*lam1*lam4
     &          + g23*lam2*lam3 + g24*lam2*lam4 + g34*lam3*lam4
      Conv(2,4) = Conv(1,1)
      Conv(3,7) = Conv(1,1)

c quadratic interpolation of uy
      s1 = DATA(23)
      s2 = DATA(24)
      s3 = DATA(25)
      s4 = DATA(26)
      s12= DATA(27)
      s13= DATA(28)
      s14= DATA(29)
      s23= DATA(30)
      s24= DATA(31)
      s34= DATA(32)

      g12 = 4 * ( s12 - (s1+s2)/2 )
      g13 = 4 * ( s13 - (s1+s3)/2 )
      g14 = 4 * ( s14 - (s1+s4)/2 )
      g23 = 4 * ( s23 - (s2+s3)/2 )
      g24 = 4 * ( s24 - (s2+s4)/2 )
      g34 = 4 * ( s34 - (s3+s4)/2 )

      Conv(1,2) = s1*lam1 + s2*lam2 + s3*lam3 + s4*lam4
     &          + g12*lam1*lam2 + g13*lam1*lam3 + g14*lam1*lam4
     &          + g23*lam2*lam3 + g24*lam2*lam4 + g34*lam3*lam4
      Conv(2,5) = Conv(1,2)
      Conv(3,8) = Conv(1,2)

c quadratic interpolation of uz
      s1 = DATA(33)
      s2 = DATA(34)
      s3 = DATA(35)
      s4 = DATA(36)
      s12= DATA(37)
      s13= DATA(38)
      s14= DATA(39)
      s23= DATA(40)
      s24= DATA(41)
      s34= DATA(42)

      g12 = 4 * ( s12 - (s1+s2)/2 )
      g13 = 4 * ( s13 - (s1+s3)/2 )
      g14 = 4 * ( s14 - (s1+s4)/2 )
      g23 = 4 * ( s23 - (s2+s3)/2 )
      g24 = 4 * ( s24 - (s2+s4)/2 )
      g34 = 4 * ( s34 - (s3+s4)/2 )

      Conv(1,3) = s1*lam1 + s2*lam2 + s3*lam3 + s4*lam4
     &          + g12*lam1*lam2 + g13*lam1*lam3 + g14*lam1*lam4
     &          + g23*lam2*lam3 + g24*lam2*lam4 + g34*lam3*lam4
      Conv(2,6) = Conv(1,3)
      Conv(3,9) = Conv(1,3)

      Dconv = TENSOR_GENERAL

      Return
      End



C ======================================================================
C boundary conditions
      Integer Function Dbc(x, y, z, label, DATA, iSYS, eBC)
C ======================================================================
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), eBC(*)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 3

      If(label.eq.3) Then
         eBC(1) = (y - 0.5) * z * (1 - y) * (1 - z) * 64
         eBC(2) = 0
         eBC(3) = 0
         Dbc = BC_DIRICHLET
      Else If(label.eq.5) Then
         eBC(1) = 0
         eBC(2) = 0
         eBC(3) = 0
         Dbc = BC_NEUMANN
      Else
         eBC(1) = 0
         eBC(2) = 0
         eBC(3) = 0
         Dbc = BC_DIRICHLET
      End if
      Return
      End


C ======================================================================
C Right hand side
      Integer Function Drhs(x, y, z, label, DATA, iSYS, Coef)
C ======================================================================
      Include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 1

      Coef(1) = 0D0
      Coef(2) = 0D0
      Coef(3) = 0D0
      Drhs = TENSOR_GENERAL

      Return
      End 



