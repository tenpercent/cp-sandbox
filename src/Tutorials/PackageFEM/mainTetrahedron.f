c ======================================================================
      Program  mainTetrahedron
c ======================================================================
c Generation of elemental mass and stiffness matrices.
c ======================================================================
      implicit none

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, tetrahedra, boundary edges, and curved edges
      Integer  nv, nt, nb

c ... coordinates of mesh points 
      Real*8   vrt(3,4)

c ... connectivity table for tetrahedra and material labels
      Integer  tet(4,1), material(1)

c ... connectivity table for boundary faces
      Integer  bnd(3,4)


c ... library FEM
      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Real*8   A(100, 100), DATAFEM(1)
      Integer  nRow, nCol, order, LDA, iSYS(MAXiSYS)

      EXTERNAL Ddiff
      Integer  Ddiff


c LOCAL VARIABLEs
      Integer  i, j

c ======================================================================
c ... Define a simple mesh consisting of 1 tetrahedron
c     used convetion: 'v'=vertex, 't'=tetrahedron, 'b'=boundary face
      nv = 3
      vrt(1,1) = 0d0
      vrt(2,1) = 0d0
      vrt(3,1) = 0d0

      vrt(1,2) = 1d0
      vrt(2,2) = 0d0
      vrt(3,2) = 0d0

      vrt(1,3) = 0.1d0
      vrt(2,3) = 1d0
      vrt(3,3) = 0d0

      vrt(1,4) = 0d0
      vrt(2,4) = 0.1d0
      vrt(3,4) = 1d0

      nt = 1
      tet(1,1) = 1
      tet(2,1) = 2
      tet(3,1) = 3
      tet(4,1) = 4
      material(1) = 1 ! x-y<0 - material 1

      Write(*,'(A)') 'Tetrahedron with the following vertices:'
      Do i = 1, 4
         Write(*,'(I3,3F8.3)') i, vrt(1, i), vrt(2, i), vrt(3, i)
      End do

      nb = 4
      bnd(1,1) = 1
      bnd(2,1) = 2
      bnd(3,1) = 3

      bnd(1,2) = 2
      bnd(2,2) = 3
      bnd(3,2) = 4

      bnd(1,3) = 3
      bnd(2,3) = 4
      bnd(3,3) = 1

      bnd(1,4) = 4
      bnd(2,4) = 1
      bnd(3,4) = 2


c     no extra data is provided for the user subroutine Ddiff
      DATAFEM(1) = 0

c     order of the quadrature rule 
      order = 2

c     leading dimension of the local matrix A
      LDA = 100

c ... L2-product of Nedelec functions. No need to populate iSYS
      Call fem3Dtet(
     &     vrt(1,1), vrt(1,2), vrt(1,3), vrt(1,4),
     &     IDEN, FEM_ND0, IDEN, FEM_ND0, 
     &     material(1), Ddiff, DATAFEM, iSYS, order,  
     &     LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 1: Bilinear form <ND0, ND0>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100F10.6)') (A(i, j), j = 1, nCol)
      End do


c ... integral of div(RT0), No need to populate iSYS
      Call fem3Dtet(
     &     vrt(1,1), vrt(1,2), vrt(1,3), vrt(1,4),
     &     DIV, FEM_RT0, IDEN, FEM_P0, 
     &     material(1), Ddiff, DATAFEM, iSYS, order,  
     &     LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 2: Linear form <div(RT0), 1>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100F10.6)') (A(i, j), j = 1, nCol)
      End do


c ... L2-product of CURL(P1^3) functions. No need to populate iSYS
      Call fem3Dtet(
     &     vrt(1,1), vrt(1,2), vrt(1,3), vrt(1,4),
     &     CURL, FEM_P1vector, CURL, FEM_P1vector, 
     &     material(1), Ddiff, DATAFEM, iSYS, order,  
     &     LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 3: Bilinear form <curl(P1^3), curl(P1^3)>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100F9.5)') (A(i, j), j = 1, nCol)
      End do


c ... L2-product of GRAD(P1) functions. No need to populate iSYS
      order = 5
      Call fem3Dtet(
     &     vrt(1,1), vrt(1,2), vrt(1,3), vrt(1,4),
     &     GRAD, FEM_P1, GRAD, FEM_P1, 
     &     material(1), Ddiff, DATAFEM, iSYS, order,  
     &     LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 4: Bilinear form <GRAD P1, GRAD P1>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100F12.8)') (A(i, j), j = 1, nCol)
      End do

      Stop 
      End



C ======================================================================
C  The user defined routines required above
C ======================================================================
c ... 2x2 diffusion tensor Diff
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, DATA(*), Coef(9, *)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0
      Ddiff = TENSOR_SCALAR

      Return
      End



