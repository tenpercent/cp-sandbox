c ======================================================================
      Program  mainBC
c ======================================================================
c This program generates a finite element system for the diffusion problem
c
c  -div K grad u + A u = F     in  Omega
c                    u = U_0   on  Gamma_D
c        K du/dn       = G_0   on  Gamma_N
c        K du/dn + S u = G_1   on  Gamma_R
c
c where Omega = [0,1]^3. The user-defined coefficients are  
c    K(x)   - positive definite tensor
c    A(x)   - non-negative reaction 
c    F(x)   - right-hand side
c    U_0(x) - essential (Dirichlet) boundary condition
c    G_0(x) - Neumann boundary condition
c    G_1(x) - Robin boundary condition
c    S(x)   - Robin boundary coefficient
c
c These coefficients are implemented in the following routines: 
c    K->Ddiff,  A->Dreact,  F->Drhs,  {U_0,G_0,G_1}->Dbc,  S->DbcRobCoef
c
c ======================================================================
      implicit none

      integer nvmax,ntmax,nbmax,namax
c ... nvmax - maximum number of mesh nodes (v)
c ... ntmax - maximum number of mesh triangles (t)
c ... nbmax - maximum number of boundary edges (b)
      parameter(nvmax = 15 000, ntmax = 2*nvmax, nbmax = 10 000)

c ... namax - maximum number of non-zero matrix entries
      parameter(namax = 90 000)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)

c ======================================================================
c Mesh definition
c ======================================================================
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb

c ... array of coordinates of mesh nodes
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for triangles and triangle labels
      Integer  tet(4,ntmax), material(ntmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(3,nbmax), labelF(nbmax)

c     IPV the array of fixed (never touched) points.
c     IFV the array of fixed (never touched) edges (faces)
c     IEV the array of fixed (never touched) elements (triangles)
      Integer   MaxPv,     MaxEv,     MaxFv
      Parameter(MaxPv = 9, MaxEv = 1, MaxFv = 1)

      Integer   nPv, IPV(MaxPv), nFv, nEv, IFV(MaxFv), IEV(MaxEv)

c ======================================================================
c For library aniFEM
c ======================================================================
      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax), status, nRow, nCol
      Real*8    A(namax), RHS(nvmax)
      Real*8   DATAFEM(1)

      EXTERNAL FEM3Dext
      Integer  FEM3Dext


c ... local variables
      Integer  i, n, iv

c =====================================================================
c ... load the initial mesh with 4 boundary edges. 
      Call loadMani(
     &      nvmax, nbmax, ntmax,
     &      nv, nb, nt,
     &      vrt, bnd, tet, labelF, material,
     &      nPv, nFv, nEv, IPV, IFV, IEV,
     &      iW, iW, "../data/cube.ani")


c ... save a GMV image of the mesh
      Call saveMgmv(nv, nb, nt, 
     &              vrt, bnd, tet, labelF, material, 
     &              'mesh.gmv', iW)


c ... Assemble the stifness matrix
c     no extra data is provided for the user subroutines Dxxxx
      DATAFEM(1) = 0

c     mark the Dirichlet points on boundary faces 1 & 2
      Do n = 1, nv
         labelP(n) = 0
      End do

      Do n = 1, nb
         If(labelF(n).EQ.1 .OR. labelF(n).EQ.2) then
            Do i = 1, 3
               iv = bnd(i, n)
               labelP(iv) = 1
            End do
         End if
      End do



c     general sparse matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)
      status = IOR(MATRIX_GENERAL, FORMAT_AMG)

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3Dext, DATAFEM, status,
     &     nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, iW)


c ... save the matrix
      Open(1,file='CSRsystem')
        Write(1,*)  nRow, IA(nRow+1)-1
        Write(1,*) (IA(i), i=1,nRow+1)
        Write(1,*) (JA(i), i=1,IA(nRow+1)-1)
        Write(1,*) ( A(i), i=1,IA(nRow+1)-1)
         Write(1,*) (RHS(i), i=1,nRow)
      Close(1)

      Stop 
      End



C ======================================================================
C  The user defined routines required above
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
     &                    lbE, lbF, lbR, lbP, DATAFEM, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      
      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DATAFEM(*)
      Integer iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C LOCAL VARIABLEs
      Integer  Ddiff, Dreact, Drhs, Dbc, DbcRobCoef
      External Ddiff, Dreact, Drhs, Dbc, DbcRobCoef

      Real*8   B(4, 4), C(3, 3), G(3), XYP(3, 4)
      Real*8   x, y, z, eBC(1)

      Integer  i,j,k,l,m, ir, ic, label, ibc
  
      Integer  iref(5), ip(4)
      DATA     iref/1,2,3,4,1/
 
C ======================================================================
      nRow = 4
      nCol = 4

c ... set up templated degrees of freedom for rows and columns. 
c     used convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
      Do i = 1, 4
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:4,1:4) is elemental vector elliptic operator;
c     in other words, for the bilinear form <grad(P1), grad(P1)>
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DATAFEM, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:4,1:4) is elemental mass matrix;
c     in other words, for the bilinear form <P1, P1>
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, Dreact, DATAFEM, iSYS, 5,
     &              4, B, ir, ic)

      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + B(i, j) 
         End do
      End do


c ... compute the right hand side vector using external function Drhs
c     in other words, the linear form <Drhs(x), P1> 
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1, 
     &              lbE, Drhs, DATAFEM, iSYS, 5,
     &              LDA, F, ir, ic)


c ... impose the Neumann and Robin boundary conditions
      Do i = 1, 3
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
         XYP(i, 4) = XY4(i)
      End do

      Do k = 1, 4
         If(lbF(k).GT.0) Then
            l = iref(k + 1)
            m = iref(l + 1)

            x = (XYP(1, k) + XYP(1, l) + XYP(1, m)) / 3
            y = (XYP(2, k) + XYP(2, l) + XYP(2, m)) / 3
            z = (XYP(3, k) + XYP(3, l) + XYP(3, m)) / 3

            ibc = Dbc(x, y, z, lbF(k), DATAFEM, iSYS, eBC)

            If(ibc.EQ.BC_NEUMANN .OR. ibc.EQ.BC_ROBIN) Then
               label = lbF(k)

               Call fem3Dtri(XYP(1, k), XYP(1, l), XYP(1, m),
     &                       IDEN, FEM_P0, IDEN, FEM_P1, 
     &                       label, Dbc, DATAFEM, iSYS, 4, 
     &                       3, G, ir, ic)

               F(k) = F(k) + G(1)
               F(l) = F(l) + G(2)
               F(m) = F(m) + G(3)
            End if

            If(ibc.EQ.BC_ROBIN) Then
               label = lbF(k)

               Call fem3Dtri(XYP(1, k), XYP(1, l), XYP(1, m),
     &                       IDEN, FEM_P1, IDEN, FEM_P1,
     &                       label, DbcRobCoef, DATAFEM, iSYS, 4,
     &                       3, C, ir, ic)

               ip(1) = k 
               ip(2) = l
               ip(3) = m
               Do i = 1, 3
                  Do j = 1, 3 
                     A(ip(i), ip(j)) = A(ip(i), ip(j)) + C(i, j)
                  End do
               End do
            End if

         End if
      End do

c ... impose Dirichlet boundary conditions at triangle nodes
c     this condition should go the last one
      Do k = 1, 3
         If(lbP(k).NE.0) Then
            x = XYP(1, k)
            y = XYP(2, k)

            ibc = Dbc(x, y, z, lbP(k), DATAFEM, iSYS, eBC)

            If(ibc.EQ.BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End



C ======================================================================
c 2x2 diffusion tensor K
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, *)
      Integer label, iSYS(*)

      Integer i, j

      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
      End do

      Coef(1, 1) = 1D0
      Coef(2, 2) = 1D1
      Coef(3, 3) = 1D1

      Coef(1, 2) =-1D0
      Coef(2, 1) =-1D0

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
c Reaction coefficient A
C ======================================================================
      Integer Function Dreact(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 1D0
      Dreact = TENSOR_SCALAR

      Return
      End



C ======================================================================
c Boundary conditions 
C ======================================================================
      Integer Function Dbc(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.LE.2) Then
         Coef(1) = x + y + z  ! non-homogeneous Dirichlet on boundaries 1 & 2
         Dbc = BC_DIRICHLET
      Else If(label.EQ.3) Then
         Coef(1) = x - y 
         Dbc = BC_NEUMANN     ! non-homogeneous Neumann on edge 3
      Else If(label.EQ.4) Then
         Coef(1) = x + z
         Dbc = BC_ROBIN
      Else
         Dbc = BC_NULL
      End if

      Return
      End



C ======================================================================
c Coefficient in Robin boundary condition
C ======================================================================
      Integer Function DbcRobCoef(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.EQ.4) Then
         Coef(1) = 1D0
         DbcRobCoef = BC_ROBIN_COEF
      Else 
         DbcRobCoef = BC_NULL
      End if

      Return
      End



C ======================================================================
c Right hand side F
C ======================================================================
      Integer Function Drhs(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      Include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


