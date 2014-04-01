      Program  main
c ==========================================================
c The program demonstrates the use of local refinement/coarsening library librcb3D.a based
c on the Marked Edge Bisection.  It starts from a simple coarse tetrahedral meshes,
c (the mesh may be arbitrary conformal)  
c refines it nlevel times according to the rule given in RefineRule,
c and coarse the refined mesh 4 times according to the rule given in CoarseRule.
c Example of RefineRule, CoarseRule are at the file end.
c Prior the refinement/coarsening the routine  InitializeRCB has to be called (once).
c The rules are designed so that the mesh is refined towards a point source moving along
c the edge of the computational polyhedral domain.
c
c Caution! Curve-linear boundary faces are processed as straight:
c          no crv data are on input/output.
c The program uses  routines [GMVmesh] from libview3D.a and [error] from libmba3D.a.
c
c The output  is  three gmv-pictures of the initial mesh (ini.gmv) and the mesh refined towards
c the initial position of the source (fst.gmv) and the final position of the source (lst.gmv)
c ==========================================================
      implicit none

c ... user defined procedures 
      external  RefineRule,  CoarseRule 

c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
      Integer   nvmax,ntmax,nbmax
      Parameter (nvmax = 200 000, ntmax = 3*nvmax, nbmax = 100 000)

c ... standard mesh arrays 
c ... number of points, tetrahedra and boundary faces
      Integer   nv, nt, nb

c ... coordinates of mesh points 
      Double precision   vrt(3,nvmax)

c ... connectivity table for tetrahedra, and their label
      Integer   tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer   bnd(3,nbmax),labelF(nbmax)

c ... maxlevel - maximum number of levels for refinement and coarsening
      Integer   maxlevel
      Parameter (maxlevel = 60)

c ... work memory (at least 11*ntmax+7)?
      Integer   MaxWi
      Parameter (MaxWi =  12 500 000)
      Integer   iW(MaxWi), WiW

c ... history of bisection
      Logical history(4*maxlevel*ntmax)

c ... list of tetrahedra
      Integer   tetpmax
      Parameter (tetpmax = 50)
      Integer   listtet((tetpmax+1)*ntmax)

c ... number of levels of refinement/coarsening
      integer  nlevel

c ... user data to be passed to CoarseRule (dummy)
      Double precision CoarseRuleData(1)

c ... user data to be passed to RefineRule 
      Double precision RefineRuleData(3)

c ... local variables
      Integer  i,i1,i2, ilevel, iERR, iPrint
      character*7 fname

c ==========================================================

c Step 1: generate initial mesh

c ... Define a simple mesh consisting of 2 tetrahedra
c ... points

      nv = 5

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

      vrt(1,5) = 0d0
      vrt(2,5) = 0d0
      vrt(3,5) = -1d0

c ... tetrahedra

      nt = 2
      tet(1,1) = 1
      tet(2,1) = 2
      tet(3,1) = 3
      tet(4,1) = 4
      material(1) = 1 ! z>0 - material 1

      tet(1,2) = 1
      tet(2,2) = 2
      tet(3,2) = 3
      tet(4,2) = 5
      material(2) = 2 ! z<0 - material 2

c ... boundary faces

      nb = 6

      bnd(1,1) = 1
      bnd(2,1) = 2
      bnd(3,1) = 4
      labelF(1) = 1

      bnd(1,2) = 1
      bnd(2,2) = 3
      bnd(3,2) = 4
      labelF(2) = 2

      bnd(1,3) = 2
      bnd(2,3) = 3
      bnd(3,3) = 4
      labelF(3) = 3

      bnd(1,4) = 1
      bnd(2,4) = 2
      bnd(3,4) = 5
      labelF(4) = 4

      bnd(1,5) = 1
      bnd(2,5) = 3
      bnd(3,5) = 5
      labelF(5) = 5

      bnd(1,6) = 2
      bnd(2,6) = 3
      bnd(3,6) = 5
      labelF(6) = 6


      Write(*,'(A,2I7)')
     &        'Initial mesh:   numbers of nodes and tets:',nv, nt

      fname = 'ini.gmv'
      call GMVmesh(fname, 99,nv,vrt, nt,tet, nb,bnd,labelF)

c Step 2: initialize data structure (work memory is at least 11*ntmax+7)

      iERR     = 0
      call InitializeRCB (nt, ntmax, nv, nvmax, vrt, tet, 
     &      MaxWi, iW, listtet, tetpmax, iERR)
      If(iERR.GT.0) Call errMes(iERR, 'main', 'error in InitializeRCB')
      iPrint   = 1       ! low level of output information

c Step 3: refine initial mesh nlevel times by local bisection. 
c         The rule is defined in RefineRule

      Do i = 1, 15

c coordinates of moving point source
        RefineRuleData(1) = 0
        RefineRuleData(2) = 0
        RefineRuleData(3) = 0.75 - i * 0.1
        
        nlevel = 15
      
        Do ilevel = 1, nlevel
          call LocalRefine (
     &          nv, nvmax, nb, nbmax, nt, ntmax,
     &          vrt, tet, bnd, material,labelF,
     &          RefineRule, ilevel,
     &          maxlevel, history, 
     &          listtet, tetpmax, RefineRuleData,
     &          MaxWi, iW, 
     &          iPrint, iERR)
           If(iERR.GT.0) Call errMes(iERR, 'main',
     &                 'error in LocalRefine')

        End do


        Write(*,'(A,2I7)')
     &        'Refined mesh:   numbers of nodes and tets:',nv, nt
c
        if (i.eq.1) fname = 'fst.gmv'
        if (i.eq.15) fname = 'lst.gmv'
        if (i.eq.1.or.i.eq.15) then
           call GMVmesh(fname, 99,nv,vrt, nt,tet, nb,bnd,labelF)
        end if
      
c Step 4: coarse the refined mesh 2 times by local coarsening. 
c         The rule is defined in CoarseRule
c
        Do ilevel = nlevel, nlevel-3, -1
          call LocalCoarse (
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, tet, bnd, material,labelF,
     &        CoarseRule, ilevel,
     &        maxlevel, history,
     &        listtet, tetpmax, CoarseRuleData,
     &        MaxWi, iW,
     &        iPrint, iERR)
           If(iERR.GT.0) Call errMes(iERR, 'main',
     &                     'error in LocalCoarse')

        End do
c      
        Write(*,'(A,2I7)')
     &        'Coarsened mesh: numbers of nodes and tets:',nv, nt
      End do

c
c
      Stop 
      End

C ================================================================
c User routine  defines the rule for local refinement depending on 
c current level.
c On input: nt  current number of elements
c           tet current connectivity table
c           vrt current coordinates
c           ilevel current level of refinement
c           RefineRuleData  data array for refining
c On output: verf marker for refinement of each element
c            (0 - no need to refine, 1 - refine by single bisection,
c             2 - refine by two   levels of bisection,
c             3 - refine by three levels of bisection preserving the shape)
C ================================================================
      Subroutine RefineRule (nt, tet, vrt, verf, ilevel, xyzp)
C ================================================================
      implicit none
      
      Integer                nt
      Integer                tet(4,*)
      Double precision       vrt(3,*)
      Integer                verf(*)
      Integer                ilevel
      Double precision       xyzp(3)

      Integer                i, j
      Double precision       xy, xyc(3)
C ================================================================

c ... refine towards the point xyzp
      Do i = 1, nt

        Do j = 1, 3
          xyc(j) = (vrt(j,tet(1,i)) + 
     &              vrt(j,tet(2,i)) +
     &              vrt(j,tet(3,i)) +
     &              vrt(j,tet(4,i))) / 4
        End do

        Call Dist3D (xyc(1), xyzp(1),xy)       

        If (xy .le. 11.0/((ilevel-1)**2)) Then 
           verf(i) =  1! one level of bisection
          Else
           verf(i) =  0 ! no need to refine
        End if
      End do  

      Return
      End
C ================================================================
c User routine  defines the rule for local coarsening depending on 
c current level.
c On input: nt  current number of elements
c           tet current connectivity table
c           vrt current coordinates
c           ilevel current level of refinement
c           CoarseRuleData  data array for coarsening
c On output: verf marker for coarsening of each element
c            (0 - no need to coarse, 1 - coarse by single merging,
c             2 - coarse by two   levels of merging,
c             3 - coarse by three levels of merging   preserving the shape)
C ================================================================
      Subroutine CoarseRule (nt, tet, vrt, verf, ilevel, CoarseRuleData)
C ================================================================
      implicit none
      
      Integer                nt
      Integer                tet(4,*)
      Double precision       vrt(3,*)
      Integer                verf(*)
      Integer                ilevel
      Double precision       CoarseRuleData(*) ! dummy in this case

      Integer                i
      Double precision       xy, xy1, xy2, xy3
C ================================================================

c Uniform coarsening with moderate rate of coarsening
c
c Here we exploit the feature of coarsening procedure:
c coarsening of a coarse tetrahedron causes multiple coarsenings of
c finer tetrahedra due to mesh conformity
      Do i = 1, nt
          verf(i) =  1 ! one level  of merging
      End do


      Return
      End

