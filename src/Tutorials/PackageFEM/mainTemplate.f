c ======================================================================
      Program  main
c ======================================================================
c This program generates a finite element system for the Stokes problem
c
c  -div grad u  + grad p = 0   in Omega
c        div u           = 0   in Omega
c
c                      u = u_0 on dOmega_1
c                      u = 0   on dOmega_2
c              du/dn - p = 0   on dOmega_3
c
c where Omega is a prism with triangular base, dOmega_1 is the side at z=0, 
c dOmega_3 is the side at z=1, and dOmega_2 is the rest of the boundary. 
c The non-homogeneous boundary condition is 
c
c    u_0 = { 0, 0, 4*y*(1-y) }.
c
c We use the P2 finite elements for the velocity u and P1 finite elements 
c for the pressure p, which are known to be the stable pair. 
c The discretization method results in a symmetric indefinite matrix.
c
c Step 1: we generate a quasi-uniform mesh using library libmba3D.a
c step 2: we generate an algebraic system using library libfem3D.a
c ======================================================================
      implicit none

      integer nvmax,ntmax,nbmax,niamax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary edges
c ... niamax- maximum number of matrix rows
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax  = 15 000,  ntmax = 2*nvmax, nbmax = 10 000)
      parameter(niamax = nvmax*6, namax = 30*niamax)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 1 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)

c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tetrahedra, and boundary facess
      Integer  nv, nt, nb

c ... coordinates of mesh points, labels of points
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for triangles and triangle labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(3,nbmax),labelF(nbmax)


c ... additional mesh arrays are needed to call mbaAnalytic() from libmba3D.a 
c     IPV the array of fixed (never touched) points
c     IFV the array of fixed (never touched) faces
c     IEV the array of fixed (never touched) elements (tets)
      Integer   nPv, IPV(6), nFv, nEv, IFV(1), IEV(1)

c ... library MBA
      Integer   MaxSkipE, MaxQItr, nEStar
      Parameter(MaxSkipE = 100, MaxQItr = 50 000)

      Real*8    Quality
      Parameter(Quality = 4D-1)

      Real*8    rQuality
      Logical   flagAuto
      Integer   iPrint, iERR

      Integer   ANI_Metric_Eucl
      External  ANI_Metric_Eucl


c ... library FEM
      include 'fem3Dtet.fd'
      include 'assemble.fd'
      Integer  IA(niamax), JA(namax)
      Real*8    A(namax), RHS(niamax)
      Integer  status, nRow, nCol

      Real*8   DATAFEM(1)
      Integer  Dbc
      EXTERNAL fem3Dext, Dbc
  

c ... local variables
      Integer   i, iv1,iv2,iv3, ibc, iSYS(MAXiSYS)
      Real*8    x,y,z, eBC(3)

c ======================================================================

c Step 1: generate a mesh
c ... Define a simple mesh for the prism consisting of 3 tets
      nv = 6
      vrt(1,1) = 0d0
      vrt(2,1) = 0d0
      vrt(3,1) = 0d0

      vrt(1,2) = 1d0
      vrt(2,2) = 0d0
      vrt(3,2) = 0d0

      vrt(1,3) = 0d0
      vrt(2,3) = 1d0
      vrt(3,3) = 0d0

      vrt(1,4) = 0d0
      vrt(2,4) = 0d0
      vrt(3,4) = 1d0

      vrt(1,5) = 1d0
      vrt(2,5) = 0d0
      vrt(3,5) = 1d0

      vrt(1,6) = 0d0
      vrt(2,6) = 1d0
      vrt(3,6) = 1d0

      nt = 3
      tet(1,1) = 1
      tet(2,1) = 2
      tet(3,1) = 3
      tet(4,1) = 4
      material(1) = 1 ! material 1
      tet(1,2) = 2
      tet(2,2) = 4
      tet(3,2) = 5
      tet(4,2) = 6
      material(2) = 1
      tet(1,3) = 2
      tet(2,3) = 3
      tet(3,3) = 4
      tet(4,3) = 6
      material(3) = 1

      nb = 8
      bnd(1,1) = 1
      bnd(2,1) = 2
      bnd(3,1) = 3
      labelF(1) = 1   ! Dirichlet face (z=0) - label 1
      bnd(1,2) = 4
      bnd(2,2) = 5
      bnd(3,2) = 6
      labelF(2) = 1   ! Dirichlet face (z=1) - label 2
      bnd(1,3) = 1
      bnd(2,3) = 2
      bnd(3,3) = 4
      labelF(3) = 3   ! Dirichlet face (y=0) - label 3
      bnd(1,4) = 2
      bnd(2,4) = 4
      bnd(3,4) = 5
      labelF(4) = 3   ! Dirichlet face (y=0) - label 3
      bnd(1,5) = 1
      bnd(2,5) = 3
      bnd(3,5) = 4
      labelF(5) = 4   ! Dirichlet face (x=0) - label 4
      bnd(1,6) = 3
      bnd(2,6) = 4
      bnd(3,6) = 6
      labelF(6) = 4   ! Dirichlet face (x=0) - label 4
      bnd(1,7) = 2
      bnd(2,7) = 5
      bnd(3,7) = 6
      labelF(7) = 5   ! Dirichlet face (x+y=1) - label 5
      bnd(1,8) = 2
      bnd(2,8) = 3
      bnd(3,8) = 6
      labelF(8) = 5   ! Dirichlet face (x+y=1) - label 5

c ... six fixed points
      nPv = 6
      IPV(1)=1
      IPV(2)=2
      IPV(3)=3
      IPV(4)=4
      IPV(5)=5
      IPV(6)=6

c ... no fixed faces and elements
      nFv = 0
      nEv = 0


c Generate a quasi-uniform mesh with 4000 tets starting from the simple mesh 
      nEStar   = 400     ! number of tets in generated mesh
      flagAuto = .TRUE.  ! default mesh generation options
c     status   = 1       ! forbid boundary tets (see aniMBA/status.fd)
      status   = 0       ! no options
      iPrint   = 1       ! low level of output information

      Call mbaAnalytic(
c group (M)
     &     nv, nvmax, nb, nbmax, nt, ntmax,
     &     vrt, bnd, tet, labelF, material,
     &     nEStar, 
c group (Dev)
     &     nPv, nFv, nEv, IPV, IFV, IEV, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     ANI_Metric_Eucl, Quality, rQuality,
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)

      If(iERR.GT.1000) Call errMes(iERR, 'main',
     &                     'unspecified error if mbaAnalytic')


c Step 2: generate the finite element system with P2-P1 elements
c     no data is provided for the user subroutine fem3Dext
      DATAFEM(1) = 0

c ... mark the Dirichlet points via labelP
      Do i = 1, nv 
         labelP(i) = 0
      End do

      Do i = 1, nb
         iv1 = bnd(1,i)
         iv2 = bnd(2,i)
         iv3 = bnd(3,i)

         x = (vrt(1,iv1) + vrt(1,iv2) + vrt(1,iv3)) / 3
         y = (vrt(2,iv1) + vrt(2,iv2) + vrt(2,iv3)) / 3
         z = (vrt(3,iv1) + vrt(3,iv2) + vrt(3,iv3)) / 3
         ibc = Dbc(x, y, z, labelF(i), DATAFEM, iSYS, eBC)

         If(ibc.EQ.BC_DIRICHLET) Then
            labelP(iv1) = 1
            labelP(iv2) = 1
            labelP(iv3) = 1
         End if
      End do

c ... general sparse matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)
      status = IOR(MATRIX_GENERAL, FORMAT_CSR)

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3Dext, DATAFEM, status,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, iW)


c ... save the matrix
      open(1,file='CSRsystem',status='UNKNOWN')
         write(1,*) nRow, IA(nRow+1)-1
         write(1,*) (IA(i),i=1,nRow+1)  
         write(1,*) (JA(i),i=IA(1),IA(nRow+1)-1)  
         write(1,*) ( A(i),i=IA(1),IA(nRow+1)-1)  
         write(1,*) (RHS(i),i=1,nRow)  
      close(1)

      Stop 
      End



C ======================================================================
C  Here are the user defined routines
C ======================================================================
c Templated routine for elemental matrix
      Subroutine fem3Dext( 
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DATAFEM, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      
      Integer lbE, lbF(3), lbR(6), lbP(4)
      Real*8  DATAFEM(1)
      Integer iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc, ANI_Dnull
      External Ddiff, Drhs, Dbc, ANI_Dnull

      Real*8   B(30, 30), XYP(3, 4)
      Real*8   x, y, z, eBC(3)

C ======================================================================
      Do i = 1, 3
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
         XYP(i, 4) = XY4(i)
      End do

c ... the size of elemental matrix (30 velocities + 4 pressures)
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

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:30,1:30) is local vector Laplacian matrix for three velocity components
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Ddiff, DATAFEM, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:4,1:30) is local divergence matrix for  velocity
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              DIV, FEM_P2vector, IDEN, FEM_P1,
     &              label, ANI_Dnull, DATAFEM, iSYS, 2,
     &              30, B, ir, ic)

      Do i = 1, 30
         Do j = 1, 4
            A(i, j + 30) = B(j, i) 
            A(j + 30, i) = B(j, i)
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
     &              lbE, Drhs, DATAFEM, iSYS, 2,
     &              LDA, F, ir, ic)

      Do i = 1, 4
         F(i + 30) = 0
      End do

c ... impose boundary conditions (assume nRow = nCol) at nodes of tetrahedron
      Do k = 1, 4
         If(lbP(k).ne.0) Then
            x = XYP(1, k)
            y = XYP(2, k)
            z = XYP(3, k)

c  ...  Dbc may change iSYS
            ibc = Dbc(x, y, z, lbP(k), DATAFEM, iSYS, eBC)

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

c  ...  Dbc may change iSYS
               ibc = Dbc(x, y, z, lbP(i), DATAFEM, iSYS, eBC)

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
C Identity tensor
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
C ======================================================================
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1
      Ddiff = TENSOR_SCALAR
      Return
      End



C ======================================================================
C Boundary conditions
      Integer Function Dbc(x, y, z, label, DATA, iSYS, Coef)
C ======================================================================
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 3
      If(label.eq.1) Then
         Coef(1) = 0
         Coef(2) = 0
         Coef(3) = y * (1 - y) * 4
         Dbc = BC_DIRICHLET
      Else If(label.eq.2) Then
         Coef(1) = 0
         Coef(2) = 0
         Coef(3) = 0
         Dbc = BC_NEUMANN
      Else
         Coef(1) = 0
         Coef(2) = 0
         Coef(3) = 0
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



