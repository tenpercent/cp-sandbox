c ==========================================================
      Program  main
c ==========================================================
c This program loads a trivial mesh (data/cub6.out) consisting 
c of six tetrahedra, refines it uniformly 3 times, refines it
c locally 4 times and then  saves the final mesh in GMV format.
c ==========================================================
      implicit none

      integer nvmax,ntmax,nbmax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 100 000)

c ... work memory
      Integer   MaxWi
      Parameter(MaxWi = 5 000 000)

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

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(6), nbfix, ibfix(1), ntfix, itfix(1) 

c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)
c ... array of flags for splitting elements
      Logical SplitFlag(ntmax)

c Local variables
      Integer   i, j, k, ii
      Real*8    xc,yc,zc
c ==========================================================
c ... Load the mesh in aft-format
      Call loadMaft(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/cub6.out")

      Write(*,*)
      Write(*,'(A,I6,A)') 'The loaded mesh has ',nt,' tetrahedra'

c ... Initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(
     &      nv, nt, vrt, tet,
     &      MapMtr, Ref2MapMtr)

c ... Refine the mesh uniformly 3 times
      Do i = 1, 3
         Call uniformRefinement(
     &      nv, nvmax, nb, nbmax, nt, ntmax,
     &      vrt, bnd, tet, labelF, material,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
      End do

      Write(*,*)
      Write(*,'(A,I6,A)')
     &      'The uniformly refined mesh has ',nt,' tetrahedra'

c ... Refine mesh locally 4 times according SplitFlag. The domain of refinement is shrinked.
c     If the domain of refinement is the same, the quality reduces.
      Do ii = 1, 4
         Do i = 1, nv
            xc = vrt(1,i)
            yc = vrt(2,i)
            zc = vrt(3,i)
            If (xc.le.0.5d0-(ii-1)*0.125d0.and.
     &          yc.le.0.5d0-(ii-1)*0.125d0.and.  
     &          zc.le.0.5d0-(ii-1)*0.125d0) Then
                iW(i) = 0
            Else
                iW(i) = 1
            End if
         End do
         Do i = 1, nt
            SplitFlag(i) = .false.  
            k = 0
            Do j = 1, 4
               k = k + iW(tet(j,i))
            End do
            If (k.eq.0) SplitFlag(i) = .true.
         End do
c ... Refine locally: unformly refined elements (Ref2MapMtr(i)>0) 
c                     with SplitFlag(i) = .true. will be split uniformly;
c                     the others may be split non-uniformly
         Call localRefinement(
     &      nv, nvmax, nb, nbmax, nt, ntmax,
     &      vrt, bnd, tet, labelF, material,
     &      MapMtr, Ref2MapMtr, SplitFlag,
     &      iW, MaxWi)
         Write(*,'(A,I6,A)') 
     &         'The locally refined mesh has ',nt,' tetrahedra'
      End do

 
c Create a GMV-file with the mesh. The name must have extension .gmv
      Call saveMgmv(nv, nb, nt, 
     &              vrt, bnd, tet, labelF, material, 
     &              "save.gmv", iW)

      Stop
      End


