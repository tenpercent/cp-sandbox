c ==========================================================
      Program  main
c ==========================================================
c This program defines a trivial mesh in a prism consisting 
c of three tetrahedra, refines it uniformly 3 times, and then
c adapts this mesh. The adapted mesh is quasi-uniform in the 
c metric defined via functions F,G,H (at the bottom of this 
c file) and contains about 4000 elements.
c
c The program uses libraries libmba3D.a and libview3D.a
c ==========================================================
      implicit none

      integer nvmax,ntmax,nbmax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
      parameter(nvmax = 15 000, ntmax = 2*nvmax, nbmax = 10 000)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)

c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax),labelF(nbmax)

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(6), nbfix, ibfix(1), ntfix, itfix(1) 

c     Data for refinement routines
c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)


c ... library MBA
c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr, nEStar
      Parameter(MaxSkipE = 300, MaxQItr = 50 000)

c Desired final mesh quality 
      Real*8    Quality
      Parameter(Quality = 4D-1)

      Real*8    rQuality
      Logical   flagAuto
      Integer   iPrint, status, iERR

c The routine which defines the nodal metric
      Integer   MetricFunction_user
      External  MetricFunction_user


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

      nvfix = 6        ! mesh fixed points (not to be touched in adaptation)
      Do i = 1, 6      ! actually these points may be detected automatically since 
         ivfix(i) = i  ! they separate boundary edges with different colors
      End do

      nbfix = 0
      ntfix = 0


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
      Call saveMgmv(nv, nb, nt, 
     &              vrt, bnd, tet, labelF, material, 
     &              "prism.gmv", iW)



c Step 2: generate a quasi-uniform mesh ith nEStar = 10000 tetrahedra
c which is quasi-unifrom in the metric defined in subroutine 
c MetricFunction_user (see below).
      nEStar   =  2000   ! number of tetrahedra in generated mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status = 1         ! forbid boundary elements
      iPrint = 1         ! average level of output information


      Call mbaAnalytic(
c group (M)
     &     nv, nvmax, nb, nbmax, nt, ntmax,
     &     vrt, bnd, tet, labelF, material,
     &     nEStar, 
c group (Dev)
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     MetricFunction_user, Quality, rQuality, 
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)


c Create a GMV-file with the adaptive mesh. File name must have extension .gmv
      Call saveMgmv(nv, nb, nt,
     &              vrt, bnd, tet, labelF, material, 
     &              "save.gmv", iW)

      Stop 
      End



C ================================================================
      Integer Function MetricFunction_user(x, y, z, Metric)
C ================================================================
C  This routine creates a metric at the given point (x,y z). The
C  metric is a 3x3 positive definite symmetric tensor:
C                M11   M12   M13
C      Metric =  M12   M22   M23
C                M13   M23   M33
C
C  Only the upper triangular part of array Metric must be defined.
C ================================================================
      Real*8  x, y, z, Metric(3, 3)

      Metric(1,1) = 1D2 
      Metric(2,2) = 1D2 
      Metric(3,3) = 1D0

      Metric(1,2) = 0D0
      Metric(1,3) = 0D0
      Metric(2,3) = 0D0

      MetricFunction_user = 0

      Return
      End

