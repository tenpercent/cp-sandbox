c ======================================================================
      Program  main
c ======================================================================
c THIS ROUTINE WILL BE REMOVED IN VERSION 2.4
c ======================================================================
c This program generates a finite element system for the diffusion problem
c
c -div D grad u = 1  in Omega
c             u = 0  on dOmega
c
c where Omega is a prism with triangular base and
c
c    D = 1 + x^2.
c
c Step 1: we generate a uniform mesh by the uniform refinement of
c         a trivial mesh using library libmba3D.a
c Step 2: we generate the algebraic system using library libfem3D.a
c ======================================================================
      implicit none

      integer nvmax,ntmax,nbmax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax = 15 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 90 000)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)

c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)

c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax),labelF(nbmax)

c     Data for mesh refinement routines
c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)



c ... library FEM
      include 'fem3Dtet.fd'
      include 'assemble.fd'

c ... memory for the global matetx
      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), DA(nvmax), RHS(nvmax)

c ... data for generating elemental matetces
      Real*8   DATAFEM(1)

c ... definition of coefficients
      Integer  Ddiff, Dbc, Drhs
      EXTERNAL Ddiff, Dbc, Drhs

      Integer  status, nRow, nCol


c ... local variables
      Integer  i

c ==========================================================

c Step 1: generate a mesh
c ... Define a simple mesh for the prism consisting of 3 tets
      nv = 6          ! mesh nodes
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

      nt = 3          ! mesh tets
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

      nb = 8          ! boundary edges
      bnd(1,1) = 1
      bnd(2,1) = 2
      bnd(3,1) = 3
      labelF(1) = 1   ! Dirichlet face (z=0) - label 1
      bnd(1,2) = 4
      bnd(2,2) = 5
      bnd(3,2) = 6
      labelF(2) = 1   ! Dirichlet face (z=1) - label 1
      bnd(1,3) = 1
      bnd(2,3) = 2
      bnd(3,3) = 4
      labelF(3) = 2   ! Dirichlet face (y=0) - label 2
      bnd(1,4) = 2
      bnd(2,4) = 4
      bnd(3,4) = 5
      labelF(4) = 2   ! Dirichlet face (y=0) - label 2
      bnd(1,5) = 1
      bnd(2,5) = 3
      bnd(3,5) = 4
      labelF(5) = 3   ! Dirichlet face (x=0) - label 3
      bnd(1,6) = 3
      bnd(2,6) = 4
      bnd(3,6) = 6
      labelF(6) = 3   ! Dirichlet face (x=0) - label 3
      bnd(1,7) = 2
      bnd(2,7) = 5
      bnd(3,7) = 6
      labelF(7) = 4   ! Dirichlet face (x=0) - label 4
      bnd(1,8) = 2
      bnd(2,8) = 3
      bnd(3,8) = 6
      labelF(8) = 4   ! Dirichlet face (x=0) - label 4

c Initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(
     &      nv, nt, vrt, tet,
     &      MapMtr, Ref2MapMtr)

c Refine the mesh uniformly 3 times
      Do i = 1, 3
         Call uniformRefinement(
     &      nv, nvmax, nb, nbmax, nt, ntmax,
     &      vrt, bnd, tet, labelF, material,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
      End do
      Write(*,*)
      Write(*,'(A,I6,A)') 'The mesh has ',nt,' tetrahedra'


c Create a GMV-file with the mesh. The name must have extension .gmv
      Write(*,'(A)') 'Saving mesh image in mesh.gmv'
      Call saveMgmv(nv, nb, nt, 
     &              vrt, bnd, tet, labelF, material, 
     &              'mesh.gmv', iW)


c Step 2: generate the finite element system for continuous P1 elements
c ... Assemble the stifness matetx
c     no extra data is provided for the user subroutine Ddiff
      DATAFEM(1) = 0

c     general sparce matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)
      status = IOR(MATRIX_GENERAL, FORMAT_AMG)

      Call BilinearFormVolume(
     &     nv, nt, vrt, tet, material,
     &     GRAD, FEM_P1, GRAD, FEM_P1,
     &     Ddiff, DATAFEM, 1,
     &     status, nvmax, namax,
     &     IA, JA, DA, A, nRow, nCol,
     &     MaxWi, iW)

      Write(*,*)
      Write(*,'(A,I8)') 'Number of non-zero entries:',  IA(nRow+1) - 1

c ... assemble the right-hand side
      Call LinearFormVolume(
     &     nv, nt, vrt, tet, material,
     &     FEM_P1,
     &     Drhs, DATAFEM, 2,
     &     RHS, nRow,
     &     MaxWi, iW)

      Write(*,'(A,I22)') 'Problem size:', nRow

c ... set up boundary conditions
      Call BoundaryConditions(
     &     nv, nb, nt, vrt, bnd, tet, labelF, 
     &     FEM_P1,
     &     Dbc, DATAFEM,
     &     IA, JA, DA, A, RHS, nRow,
     &     MaxWi, iW)

 
c ... save the matrix
      Open(1,file='CSRsystem')
         Write(1,*)  nRow, IA(nRow+1)-1
         Write(1,*) (IA(i), i=1,nRow+1)
         Write(1,*) (JA(i), i=1,IA(nRow+1)-1)
         Write(1,*) ( A(i), i=1,IA(nRow+1)-1)
         Write(1,*) (RHS(i),i=1,nRow)
      Close(1)

      Stop 
      End



C ================================================================
C  The user defined routines required above
C ================================================================

c ... Diffusion tensor
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9,*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1,1) = 1d0 + x ** 2
      Ddiff = TENSOR_SCALAR

      Return
      End


c ... Boundary condition
      Integer Function Dbc(x, y, z, label, DATA, iSYS, Coef)
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If (label.eq.2) Then
         Coef(1) = 0d0
         Dbc = BC_NEUMANN
      Else
         Coef(1) = 0d0
         Dbc = BC_DIRICHLET
      End if

      Return
      End


c ... Right hand side
      Integer Function Drhs(x, y, z, label, DATA, iSYS, Coef)
      Include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 1D0
      Drhs = TENSOR_SCALAR
 
      Return
      End



