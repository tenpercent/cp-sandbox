C
C   ALL ROUTINES IN THIS FILE WILL BE REMOVED IN VERSION 2.4
C
C   THE TEMPLATE.F WILL BECOME OUR BASIS ASSEMBLING ROUTINE    

C ======================================================================
      Subroutine BilinearFormVolume(
C ======================================================================
     &           nP, nE, XYP, IPE, lbE,
     &           operatorA, FEMtypeA, operatorB, FEMtypeB,
     &           D, DATA, order,
     &           assembleStatus, MaxIA, MaxA,
     &                           IA, JA, DA, A, nRow, nCol,
     &           MaxWi, iW)
C ======================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'
C ======================================================================
C  The routine computes a stiffness matrix for bilinear form (2).
C  THIS ROUTINE IS DEPRECATED.
C ======================================================================
C *** GENERAL ***
C  Here is the decription of input parameters for all routines in
C  this file.
C
C  The order of calls of these routines is important to save on 
C  arithmetical operations. BilinearFormVolume( ... ) should be called 
C  first, and BoundaryConditions( ... ) should be called last.
C
C
C *** INTRODUCTION ***
C  The routines assemble the local matrices for the bilinear form
C
C  (2)            <D OpA(u), OpB(v)>
C
C  where D is a tensor, OpA and OpB are linear operators.
C
C  In order to compute the right hand side, we can use the following
C  trick:
C
C  (3)             f(v) = < D v, FEM_P0 >
C
C  where action of D is given by function f, read also comments
C  about formula (4) below.
C
C  A general variational problem may consists of a few bilinear
C  forms. We supply a library for operating with sparce matrices
C  given in a few formats (see files  algebra_*.f).
C
C
C *** GRID ***
C   The grid is the union of elements (tetrahedra, cubes, etc.)
C   At the moment only tetrahedra are allowed.
C
C     nP  - the number of points (P)
C     nF  - the number of faces (F)
C     nE  - the number of elements (E)
C
C     XYP(3, nP) - The Cartesian coordinates of mesh points
C
C     IPF(3, nF) - connectivity list of boundary faces
C
C     IPE(4, nE) - connectivity list of element. On output,
C                  each column of it is ordered by increasing.
C
C     lbF(nF)    - boundary identificator for imposing essential
C                  boundary conditions,
C                    (example: unit cube has 6 boundaries which
C                     may have or not different identificators)
C     lbE(nE)    - element indentificator (a positive number)
C
C
C *** BILINEAR & LINEAR FORMS ***
C   The bilinear form is defined as follows: <D OpA(u), OpB(v)>
C   where D is a tensor, OpA and OpB are linear operators, and
C   u and v are finite element functions.
C
C     operatorA - the available operators are:
C                 IDEN  - identity operator
C                 GRAD  - gradient operator
C                 DIV   - divergence operator
C                 CURL   - rotor operator
C
C                 REMARK: It is assumed that operatorA is richer than
C                          operatorB.
C
C     FEMtypeA - The possible choices of finite elements are available:
C
C                FEM_P0        - piecewise constant
C                FEM_P1        - piecewise linear
C                FEM_P2        - picewise quadratic
C                FEM_P1vector  - vector piecewise linear. The
C                                unknowns are ordered first by
C                                vertices and then by the space
C                                directions (x, y and z)
C                FEM_8P1vector - piecewise linear on the uniform
C                                partition onto 8 tetrahedra
C                FEM_RT0       - lower order Raviart-Thomas finite elements
C                FEM_ND0       - lower order Nedelec finite elements
C                FEM_CR1       - Crouzeix-Raviart finite elements
C
C     operatorB - see operatorA
C
C     FEMtypeB  - see FEMtypeA
C
C
C *** DESCRIPTION OF TENSORS ***
C     D         - the external function D has the following format
C
C     (4)  INTEGER FUNCTION D(x,y,z, label, DATA, iSYS, Diff)
C
C     where
C              (x, y, z) - Real*8 Cartesian coordinates of a 3D point where
C                          tensor Diff should be evaluated
C
C              label     - identificator of either element or face
C
C              DATA      - Real*8 user given data (a number or an array)
C
C              iSYS(21)  - system buffer for information exchange:
C                   On Input:
C                     iSYS( 3)    <- tetrahedrom number
C                     iSYS( 4:7)  <- numbers of vertices
C                     iSYS( 8:17) <- supported ONLY in template.f
C                     iSYS(18)    <- the number of points
C                     iSYS(19)    <- the number of edges
C                     iSYS(20)    <- the number of faces
C                     iSYS(21)    <- the number of tetrahedra
C                   On Output:
C                     iD  = iSYS(1) -> number of rows in tensor Diff
C                     jD  = iSYS(2) -> number of columns in tensor Diff
C
C              Diff(9, jD) - Real*8 matrix. In the case of vector finite
C                            elements U = (u, v, w) the following ordering
C                            should be used:
C                            u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z.
C
C                            REMARK 1: the leading dimension of Diff is 9.
C
C              There are three type of external functions:
C              1. Bilinear form tensors:
C                 Parameter iSys(4) is not used here.
C
C              2. Boundary conditions:
C                 Parameter Diff means the value of the essential boundary
C                 condition at point (x,y,z) which is usually the middle
C                 point of a mesh object (point, face or edge).
C                 At the moment iSys(1) and iSys(2) are not used.
C
C              3. Right hand sides:
C                 At the moment iSys(1) and iSys(2) are not used.
C
C              Examples of finction D may be found in file user.f.
C
C
C     DATA    - Real*8 user given data (a number or an array)
C
C     tensor  - a short description of the tensor. The available tensors are
C               TENSOR_NULL      - identity tensor.
C                                  Example: isotropic diffusion problem in
C                                  homogeneous domains.
C
C               TENSOR_SCALAR    - piecewise constant scalar tensor
C                                  Example: isotropic diffusion problem in
C                                  heterogeneous domains.
C
C               TENSOR_SYMMETRIC - symmetric tensor
C                                  Example: anisotropic diffusion problem in
C                                  heterogeneous/homegeneous domains.
C
C               TENSOR_GENERAL   - general tensor
C
C     order - the order of a quadrature formula:
C             order = 1 - quadrature formula with one center point
C             order = 2 - quadrature formula with 4 points inside tetrahedron
C
C
C *** SPARSE MATRIX ***
C   By default, our sparse matrix is assembled in one of the row formats
C   (IA, JA, DA, A) specified by assembleStatus. Other formats are supported
C   thorough the library of format transformations.
C
C     assembleStatus - some a priory information about the matrix A:
C
C                      MATRIX_SYMMETRIC - symmetric matrix
C                      MATRIX_GENERAL   - general matrix
C
C                      FORMAT_AMG       - format used in AMG
C                      FORMAT_CSR       - diagonal of A is saved only in DA
C
C                      REMARK: Any logical AND can be used to defined
C                              the variable, possible contradictions will
C                              be checked by the code.
C
C     MaxA      - the maximal size of array A - the maximal
C                 number of nonzero entries
C
C     IA, JA, DA, A - sparcity structure of matrix A:
C
C                     IA(nRow + 1) - IA(k + 1) equals to the number of
C                                    nonzero entries in first k rows
C                                    plus 1
C
C                     JA(M)        - list of column indexes of non-zero
C                                    entries ordered by rows,
C                                    M = IA(nRow + 1) - 1
C
C                     A(M)         - non-zero entries ordered as in JA
C
C                     DA(nRow)     - main diagonal of A
C
C     nRow      - the number of rows in A
C     nCol      - the number of columns in A, nCol = nRow for
C                 symmetric bilinear forms.
C
C
C *** WORKING MEMORY ***
C     MaxWi  - the size of the wirking integer array
C
C     iW(MaxWi)  - the integer working array. On output it contains:
C                  L          = iW(1) - the leading dimension of IAE,
C                  IAE(L, nE) = iW(2) - connectivity list between tetrahedra
C                                       and degrees of freedom (points,
C                                       faces, or edges).
C
C                  REMARK: L * nE memory cells are reserved for other
C                          routines assuming that they are called
C                          after computing the stiffness matrix.
C
C ======================================================================
C  Note:
C       Input parameters:  nP, nE, XYP,
C                          operatorA, FEMtypeA, operatorB, FEMtypeB,
C                          D, DATA, order,
C                          assembleStatus, MaxIA, MaxA,
C                          MaxWi, iW
C
C       Input / output:    IPE, iW
C
C       Output parameters: IA, JA, A, nRow, nCol
C
C ======================================================================
C
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Yu. Vassilevski (yuri.vassilevs@gmail.com)
C *** Date:    2000 - 2009 
C *** Complains & Comments: lipnikov@gmail.com
C *** External routines:  AMG1R5, D
C
C ======================================================================
      Real*8   XYP(3, *)
      Integer  IPE(4, *), lbE(*)

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  order, D

      Real*8   DATA(*)
      EXTERNAL D

      Integer  IA(*), JA(*), assembleStatus
      Real*8   A(*), DA(*)

      Integer  iW(*)

C ======================================================================
C Local variables
      Real*8  Aloc(MaxSize, MaxSize)

      Integer label, iSYS(MAXiSYS)
      Logical FEMtype, flagMS, flagMG, flagFA, flagFR

      Integer findY
      Character*80 message

C ======================================================================
      If(FEMtypeA.EQ.FEM_P2 .OR. FEMtypeB.EQ.FEM_P2)
     &   Call errMesFEM(2001, 'BilinearFormVolume',
     &               'FEM_P2 is not supported')

c ... default values for flags
      flagMG = .TRUE.

      flagMS = IAND(assembleStatus, MATRIX_SYMMETRIC).NE.0
      flagMG = IAND(assembleStatus, MATRIX_GENERAL).NE.0

      flagFA = IAND(assembleStatus, FORMAT_AMG).NE.0
      flagFR = IAND(assembleStatus, FORMAT_CSR).NE.0

      If((flagMS .AND. flagMG) .OR.
     &    flagFA .AND. flagFR)
     &     Call errMesFEM(2001, 'BilinearFormVolume',
     &                 'assembleStatus is wrong')

      Call order2D(4, nE, IPE)

      FEMtype = FEMtypeA .EQ. FEMtypeB


c ... memory distribution
      iL = 1
      iIAE = iL + 1
      inEP = iIAE + 6 * nE

      If(FEMtypeA.EQ.FEM_RT0) Then
         iIEF = inEP
         inEF = iIEF + 4 * nE
         inEP = inEF + 4 * nE
      End if

      iIEP = inEP + MaxIA + 1
      iIP1 = iIEP + 6 * nE
      iEnd = iIP1 + MaxIA + 1

      If(iEnd.GT.MaxWi) Then
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iEnd
         Call errMesFEM(1001, 'BilinearFormVolume', message)
      End if



      If(FEMtypeA.EQ.FEM_P1) Then
         L = 4
         Call copy2D(L, nE, IPE, iW(iIAE))
         nRow = nP
      Else If(FEMtypeA.EQ.FEM_RT0) Then
         L = 4
         Call listE2F(nP, nF, nE, IPE, iW(iIAE), iW(inEP), iW(iIEP))
         nRow = nF

         Call backReferences(nF, nE, 4,4, iW(iIAE), iW(inEF), iW(iIEF))
      Else If(FEMtypeA.EQ.FEM_ND0) Then
         L = 6
         Call listE2R(nP, nR, nE, IPE, iW(iIAE), iW(inEP), iW(iIEP))
         nRow = nR
      End if

      If(nRow.GT.MaxIA) Call errMesFEM(1013, 'BilinearFormVolume',
     &                             'local parameter MaxIA is small')

      Call backReferences(nRow, nE, L, L, iW(iIAE), iW(inEP), iW(iIEP))


      If(FEMtype) Then
         nCol = nRow
         iIBE = iIAE
         M = L
      Else
         Call errMesFEM(6001, 'BilinearFormVolume',
     &              'rectangular matrices are not supported')
      End if


c ... shifting some of the memory pointers
      iIAE = iIAE - 1
      iIBE = iIBE - 1
      inEP = inEP - 1
      iIEP = iIEP - 1


c ... the symbolic assembling of the matrix
      Do n = 1, nRow
         iW(iIP1 + n) = 0
      End do

      IA(1) = 1

      i2 = 0
      m2 = 0
      Do n = 1, nRow
         m1 = m2

         i1 = i2 + 1
         i2 = iW(inEP + n)
         If(m2 + L * (i2 - i1) + 1.GT.MaxA) Call errMesFEM(1014,
     &        'BilinearFormVolume', 'local parameter MaxA is small')

         If(flagFA .AND. n.LE.nCol) Then
            m2 = m2 + 1
            JA(m2) = n
            iW(iIP1 + n) = 1
         End if

         Do i = i1, i2
            iE = iW(iIEP + i)
            Do 100 j = 1, L
               iP = iW(iIAE + (iE - 1) * L + j)

               If(flagMS .AND. iP.LT.n) Goto 100

               If(iW(iIP1 + iP).EQ.0) Then
                  iW(iIP1 + iP) = 1

                  m2 = m2 + 1
                  JA(m2) = iP
               End if
 100        Continue
         End do
         IA(n + 1) = m2 + 1

         Do i = m1 + 1, m2
            iW(iIP1 + JA(i)) = 0
         End do
      End do


C ... fill in the sparcity structure
      Do n = 1, IA(nRow + 1)
         A(n) = 0D0
      End do

      Do n = 1, nRow
         DA(n) = 0D0
      End do

      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)
 
         label = lbE(n)
         Call encodeISYSfull(n, iP1,iP2,iP3,iP4, 0,0,0,0,0,0, 
     &                       0,0,0,0, nP, nR, nF, nE, iSYS)

         Call FEM3Dtet(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3), XYP(1,iP4),
     &                 operatorA, FEMtypeA, operatorB, FEMtypeB,
     &                 label, D, DATA, iSYS, order,
     &                 MaxSize, Aloc, ir, ic)

c  ...  correct orientation of basis functions
         If(FEMtypeA.EQ.FEM_RT0) Then
            Do i = 1, 4
               iFt = iW(iIAE + 4 * (n - 1) + i)

               iEt = findY(iFt, iW(iIEF), iW(inEF), n)

               If(n.LT.iEt) Then
                  Do j = 1, ic
                     Aloc(i, j) = -Aloc(i, j)
                     Aloc(j, i) = -Aloc(j, i)
                  End do
               End if
            End do

         Else If(FEMtypeA.EQ.FEM_ND0) Then
            i = 0
            Do i1 = 1, 3
               Do i2 = i1 + 1, 4
                  i = i + 1

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)

                  If(iP1.LT.iP2) Then
                     Do j = 1, ic
                        Aloc(i, j) = -Aloc(i, j)
                        Aloc(j, i) = -Aloc(j, i)
                     End do
                  End if
               End do
            End do
         End if

c  ...  assemble
         Do i = 1, ir
            iP = iW(iIAE + (n - 1) * L + i)
            DA(iP) = DA(iP) + Aloc(i, i)

            js = 1
            If(flagFR .AND. flagMS) js = i + 1

            Do 200 j = js, ic
               iQ = iW(iIBE + (n - 1) * M + j)
               Do k = IA(iP), IA(iP + 1) - 1
                  If(JA(k).EQ.iQ) Then
                     A(k) = A(k) + Aloc(i, j)
                     Goto 200
                  End if
               End do
 200        Continue
         End do
      End do


c ... output of the connectivity list (only for symmetric problems)
      iW(1) = L

      Return
      End



C ======================================================================
      Subroutine LinearFormVolume(
C ======================================================================
     &           nP, nE, XYP, IPE, lbE,
     &           FEMtypeA, 
     &           D, DATA, order,
     &           F, nRow,
     &           MaxWi, iW)
C ======================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'
C ======================================================================
C  The routine computes the right hand side for linear form (3).
C  THIS ROUTINE IS DEPRECATED.
C ======================================================================
      Real*8   XYP(3, *)
      Integer  IPE(4, *), lbE(*)

      Integer  FEMtypeA, order, D

      Real*8   DATA(*)
      EXTERNAL D

      Real*8   F(*)

      Integer  iW(*)

C ======================================================================
C Local variables
      Real*8  Floc(MaxSize, MaxSize)
      Integer label, iSYS(MAXiSYS), findY

C ======================================================================
c ... memory distribution 
      L  = iW(1)

      iL = 1
      iIAE = iL + 1
      iEnd = iIAE + 6 * nE

      If(FEMtypeA.EQ.FEM_RT0) Then
         iIEF = iEnd
         inEF = iIEF + 4 * nE
         iEnd = inEF + 4 * nE
      End if


      Do n = 1, nRow
         F(n) = 0D0
      End do


      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         label = lbE(n)
         Call encodeISYSfull(n, iP1,iP2,iP3,iP4, 0,0,0,0,0,0, 
     &                       0,0,0,0, nP, nR, nF, nE, iSYS)

         Call FEM3Dtet(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3), XYP(1,iP4),
     &                 IDEN, FEMtypeA, IDEN, FEM_P0,
     &                 label, D, DATA, iSYS, order,
     &                 MaxSize, Floc, ir, ic)

c  ...  correct orientation of basis functions
         If(FEMtypeA.EQ.FEM_RT0) Then
            Do i = 1, 4
               iFt = iW(iIAE + 4 * (n - 1) + i - 1)

               iEt = findY(iFt, iW(iIEF), iW(inEF), n)

               If(n.LT.iEt) Floc(1, i) = -Floc(1, i)
            End do
         Else If(FEMtypeA.EQ.FEM_ND0) Then
            i = 0
            Do i1 = 1, 3
               Do i2 = i1 + 1, 4
                  i = i + 1

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)

                  If(iP1.LT.iP2) Floc(1, i) = -Floc(1, i)
               End do
            End do
         End if

c  ...  assemble
         Do i = 1, ic
            m = iW(iIAE + L * (n - 1) + i - 1)
            F(m) = F(m) + Floc(1, i)
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine BoundaryConditions(
C ======================================================================
     &           nP, nF, nE, XYP, IPF, IPE, lbF,
     &           FEMtypeA,
     &           D, DATA,
     &           IA, JA, DA, A, F, nRow,
     &           MaxWi, iW)
C ======================================================================
      Include 'assemble.fd'
      Include 'fem3Dtet.fd'
C ======================================================================
C  The routine imposes essential boundary conditions.
C ======================================================================
      Real*8   XYP(3, *)
      Integer  IPF(3, *), IPE(4, *), lbF(*)

      Integer  FEMtypeA, D

      Real*8   DATA(*)
      EXTERNAL D

      Integer  IA(*), JA(*)
      Real*8   A(*), DA(*), F(*)

      Integer  iW(*)

C ======================================================================
C Local variables
      Integer dof, label, ipr(2, 3), IAA(3), iSYS(MAXiSYS), p
      Real*8  XYA(3, 3), eBC(3)
      Character*80 message

      DATA ipr /1,2, 1,3, 2,3/

C ======================================================================
c ... memory distribution
      iL = 1
      L  = iW(1)

      iIAE = iL + 1
      iIFE = iIAE + L * nE
      inEP = iIFE + 4 * nE
      iIEP = inEP + nRow
      iEnd = iIEP + L * nE

      If(iEnd.GT.MaxWi) Then
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iEnd
         Call errMesFEM(1001, 'BoundaryCondition', message)
      End if


      Call makMb(nP, nF, nE, IPF, IPE, iW(iIFE), iW(inEP), iW(iIEP))


      Do n = 1, nE
         Do 100 i1 = 1, 4
            iF = iW(iIFE + 4 * (n - 1) + i1 - 1)
            If(iF.EQ.0) Goto 100

            If(FEMtypeA.EQ.FEM_P1) Then
               dof = 3
               Do k = 1, dof
                  Do i = 1, 3
                     XYA(i, k) = XYP(i, IPF(k, iF))
                  End do

                  IAA(k) = IPF(k, iF)
               End do

            Else If(FEMtypeA.EQ.FEM_RT0) Then
               dof = 1
               Do i = 1, 3
                  XYA(i, k) = (XYP(i, IPF(1, iF))
     &                       + XYP(i, IPF(2, iF))
     &                       + XYP(i, IPF(3, iF))) / 3
               End do

               IAA(1) = iW(iIAE + L * (n - 1) + i1 - 1)

            Else If(FEMtypeA.EQ.FEM_ND0) Then
               dof = 3
               Do k = 1, dof
                  Do i = 1, 3
                     XYA(i, k) = (XYP(i, IPF(ipr(1, k), n))
     &                          + XYP(i, IPF(ipr(2, k), n))) / 2
                  End do

                  ir = ipr(1, k) + ipr(2, k) - 1
                  IAA(k) = iW(iIAE + L * (n - 1) + ir - 1)
               End do
            End if


c  ...  changing matrix and right-hand side
            Do k = 1, dof
               label = lbF(iF)
               iBC = D(XYA(1,k), XYA(2,k), XYA(3,k), 
     &                           label, DATA, iSYS, eBC)

               If(iBC.NE.BC_DIRICHLET) Goto 100
           
               m = IAA(k)
               F(m) = DA(m) * eBC(1)

               Do 50 i = IA(m), IA(m + 1) - 1
                  p = JA(i)

                  If(p.GT.0 .AND. p.NE.m) Then
                     A(i) = 0D0

                     Do j = IA(p), IA(p + 1) - 1
                        If(JA(j).EQ.m) Then
                           F(p) = F(p) - A(j) * eBC(1)
                           A(j) = 0D0

                           Goto 50
                        End if
                     End do
                  End if
 50            Continue
            End do
 100     Continue
      End do

      Return
      End



C ======================================================================
      Subroutine makMb(nP, nF, nE, IPF, IPE, IFE, nEP, IEP)
C ======================================================================
C  The routine computes connectivity lists for BOUNDARY mesh faces
C ======================================================================
      Integer IPF(3, *), IPE(4, *), IFE(4, *)
      Integer IEP(*), nEP(*)

C ======================================================================
C group (Local variables)
      Integer ipb(3, 4)
      Logical cmpE

      DATA ipb /1,2,3, 1,2,4, 1,3,4, 2,3,4/

C ======================================================================
      Call backReferences(nP, nF, 3, 3, IPF, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 4
            IFE(i, n) = 0

            ip1 = IPE(ipb(1, i), n)
            ip2 = IPE(ipb(2, i), n)
            ip3 = IPE(ipb(3, i), n)

            If(cmpE(ip1, ip2, ip3, IEP, nEP, 0, iF)) Then
               IFE(i, n) = iF
            End if
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine order2D(L, nE, IPE)
C ======================================================================
C  The routine orders each column of 2D array
C ======================================================================
      Integer IPE(L, *)

      Do n = 1, nE
         Do i = 1, L
            im = i
            mv = IPE(i, n)
            Do j = i + 1, L
               If(mv.GT.IPE(j, n)) Then
                  im = j
                  mv = IPE(j, n)
               End if
            End do
            IPE(im, n) = IPE(i,  n)
            IPE(i,  n) = mv
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine copy2D(L, nE, IPE1, IPE2)
C ======================================================================
C  The routine copies IPE1 into IPE2
C ======================================================================
      Integer IPE1(L, *), IPE2(L, *)

      Do n = 1, nE
         Do i = 1, L
            IPE2(i, n) = IPE1(i, n)
         End do
      End do

      Return
      End
