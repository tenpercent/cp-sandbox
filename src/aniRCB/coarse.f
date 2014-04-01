      Subroutine  LocalCoarse  (
C ================================================================    
     &        nP, nPmax, nF, nFmax, nE, nEmax,
     &        XYP, IPE, IPF, lbE, labelF,
     &        CoarseRule, ilevel,
     &        maxlevel, history,
     &        listtet, tetpmax, CoarseRuleData,
     &        MaxWi, iW,
     &        iPrinE, iERR)
C ================================================================ 

      Implicit none

c ... user procedure
      external  CoarseRule

c ... standard mesh arrays 
c ... number of points, tetrahedra and boundary edges
      Integer  nP, nPmax, nE, nEmax, nF, nFmax

c ... maxlevel - maximum number of bisection levels
      Integer maxlevel

c ... coordinates of mesh points 
      Double precision   XYP(3,nPmax)

c ... connectivity table for tetrahedra and tetrahedra labels
      Integer  IPE(4,nEmax), lbE(nEmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  IPF(3,nFmax), labelF(nFmax)

c ... level of bisection
      Integer  ilevel

c ... history of bisection
      Logical history(4,maxlevel, nEmax)
      
c ... list of tetrahedra
      Integer   tetpmax
      Integer   listtet(tetpmax+1,*)

c ... user data to be passed to CoarseRule
      Double precision CoarseRuleData(*)

c ... work memory
      Integer   MaxWi
      Integer   iW(MaxWi)
      
      Integer   iPrinE, iERR

c ... local point
      Integer   iP_pbisE, iP_flag, iP_level 
      Integer   iP_fpe, ip_edp 
      Integer   ip_IEF, ip_verf, ip_pconfP
      Integer   ip_realP, ip_realE
      Integer   ip_liW, MaxWil
C ================================================================    


c ... ip_IEF, ip_verf, ip_pbisE, ip_pconfP
      iP_pbisE   =                   1
      iP_flag    = iP_pbisE        + 3*nEmax
      iP_level   = iP_flag         + nEmax
      iP_fpe     = iP_level        + nEmax
      ip_edp     = iP_fpe          + 24
      ip_IEF     = iP_edp          + 16
      ip_verf    = ip_IEF          + 4*nEmax
      ip_pconfP  = ip_verf         + nEmax
      ip_realE   = ip_pconfP       + 4*nEmax
      ip_realP   = ip_realE        + nEmax
      ip_liW     = ip_realP        + nPmax
      
      If (MaxWi .le. 15*nEmax+nPmax+41) Then
          iERR = 1010
          call errMes(iERR, 'LocalRefine','MaxWi  is too small')
      End if
      MaxWil = MaxWi - 15*nEmax-nPmax-41

c ... define ip_IEF, ip_pbisE
      call OrientandMakeNeighbors (
     &           nP, nPmax, nF, nFmax, nE, nEmax,
     &           IPE, iW(ip_IEF), IPF, labelF,
     &           listtet, tetpmax,
     &           iERR)
     
c ... define ip_verf
      call CoarseRule (nE, IPE, XYP, iW(ip_verf), ilevel,CoarseRuleData)

c ... bisection
      call Untetr (nP, nPmax, nE, nEmax, IPE, iW(ip_IEF), XYP, 
     &      iW(ip_pbisE), iW(ip_flag), lbE, iW(ip_verf), 
     &      iW(ip_pconfP), iW(ip_level), maxlevel, history, 
     &      iW(ip_realE), iW(ip_realP),
     &      iW(ip_fpe), iW(ip_edp), listtet, tetpmax,
     &      iERR)

c ... define IPF
      call RestoreBndEdges (nF, nFmax, nE, IPE, iW(ip_IEF),
     &                      IPF, labelF, iERR)
                  
      Return
      End

C  ================================================================ 
c  UNTETR  produces (local) coarsening of a given set of tetrahedra.
c  The algorithm is a modified version of method by
c  L. Arnold, A. Mukherjee and L. Pouly 
C  ================================================================ 
       Subroutine Untetr (
C  ================================================================ 
     &           nP, nPmax, nE, nEmax,
     &           IPE, IEF, XYP, 
     &           pbisE, flag, lbE, verf, pconfP,
     &           level, maxlevel, history,
     &           realE, realP,
     &           fpe, edp, listtet, tetpmax, iERR)
C  ================================================================ 
c
c                          l3
c                         /\ \ 
c                        /  \   \ 
c                       /    \     \
c                      /      \       \
c                     /        \        \  l4
c                    /          \   .   / 
c                   /          . \     /
c                  /      .       \   /
c                 /   .            \ /
c              l1 ------------------ l2
c
c  XYP(l, i)    : l-th coordinate of the i-th node
c  verf(k)      : flag whether to bisect tetrahedron k or not
c  IPE(l, k)    : No. of l-th vertex of tetrahedron k (anticlockwise)
c  IEF(l, k)    : neighbour of tetrahedron k opposite to vertex l
c  pconfP(l, k) : flag whether there's a "conflict" on point l
c  pbisE(l,k)   : No. of the refinement and marked edges
c                 of tetrahedron k
C  ================================================================ 
      implicit none
c
      Double precision  XYP(3,*)
      Integer           nP, nPmax, nE, nEmax
      Integer           maxlevel
      Integer           IPE(4,*), IEF(4,*)
      Integer           pbisE(3,*), flag(*)
      Integer           lbE(*)
      Integer           verf(*)
      Integer           pconfP(4,*)
      Integer           level(*)
      Integer           realE(*)
      Integer           realP(*)
      Integer           tetpmax
      Integer           listtet(tetpmax+1,*)
      Integer           iERR
      Logical           history(4,maxlevel, *)
                      
      Integer           fpe(4,6), edp(4,4)
      Integer           i, j, k, m, maxlevelmesh
C  ================================================================ 
c ... initialization:

      Do i = 1, nE
        realE(i) = 1
      End do  
      
      Do i = 1, nP
        realP(i) = 1
      End do  

      Do i = 1, nE
        Do j = 1, 4
          pconfP(j, i) = 0
        End do
      End do

      Do i = 1, nE
       IF (verf(i) .ge. 1) Goto 1
      End do

      Return

c ... coarse elements until there is  no conflict:
1     Continue

      maxlevelmesh=0
      Do i = 1, nE
       IF ((realE(i) .eq. 1) .and.  
     &      (verf(i) .ge. 1) .and.
     &      (level(i) .gt. maxlevelmesh))
     &   maxlevelmesh = level(i)
      End do


      Do i = 1, nE
        If ((verf(i) .ge. 1) 
     &      .and. (realE(i) .eq. 1) 
     &      .and. (level(i) .eq. maxlevelmesh)) 
     &    call Unteil(i, nP, nPmax, nE, nEmax, IPE, IEF, XYP,
     &           fpe, edp, pbisE, flag, lbE, verf, pconfP,
     &           level, maxlevel, history, 
     &           realE, realP, listtet, tetpmax, iERR)
      End do
         
      Do i = 1, nE
        If ((verf(i) .ge. 1) .and. (realE(i) .eq. 1))  Go to 1
      End do

c ... delete point
      k = 0
      Do i = 1, nP 
        If (realP(i) .eq. 1) Then
          k = k + 1
          realP(i) = k
          Do j = 1, 3
            XYP(j,k) = XYP(j,i)
          End do 
        End if
      End do
      nP = k

c ... delete tetrahedra 
      k = 0
      Do i = 1, nE
        If (realE(i) .eq. 1) Then
          k = k + 1
          realE(i) = k
        End if
      End do

      Do i = 1, nE
       If (realE(i) .ge. 1) Then
         Do j = 1, 4
           IPE(j, realE(i)) = realP(IPE(j, i))
           If (IEF(j, i).ge.1) Then
             IEF(j,realE(i)) = realE(IEF(j, i))
            Else
             IEF(j,realE(i)) = IEF(j,i)
           End if  
         End do

         Do m = 1, 3
           pbisE(m, realE(i)) = pbisE(m, i)
         End do  
         flag(realE(i)) = flag(i)
         lbE(realE(i)) = lbE(i)

         level(realE(i)) = level(i)
         Do j=1, level(i)
           Do m = 1, 4
             history(m, j, realE(i)) = history(m, j, i)
           End do  
         End do  
       End if
      End do 
      nE = k


      Return
      End
C  ================================================================ 

C  ================================================================ 
c  UNTEIL coarsen the tetrahedron i from i and neighbor
c
c  INPUT:      i - no. of tetrahedron
C  ================================================================ 
      Subroutine Unteil(
C  ================================================================ 
     &           i, nP, nPmax, nE, nEmax, IPE, IEF, XYP, 
     &           fpe, edp, pbisE, flag, lbE, verf, pconfP,
     &           level, maxlevel, history, realE, realP,
     &           listtet, tetpmax, iERR)
C  ================================================================ 
      implicit none
c
      External          twoedges_on_face
      Logical           twoedges_on_face
c
      Integer           i
      Integer           nP, nPmax, nE, nEmax
      Integer           maxlevel
      Integer           IPE(4,*), IEF(4,*)
      Double precision  XYP(3,*) 
      Integer           fpe(4,6), edp(4,4)
      Integer           pbisE(3,*), flag(*)
      Integer           lbE(*)
      Integer           verf(*)
      Integer           pconfP(4,*)
      Integer           level(*)
      Integer           realE(*)
      Integer           realP(*)
      Integer           tetpmax
      Integer           listtet(tetpmax+1,*)
      Integer           iERR
      Logical           history(4,maxlevel, *)
     
      Integer           pt1, pt2, pt3, pt4
      Integer           ii, j, jj,  k, kk, l, m, v
      Logical           flagcoarse
      Integer           deltet, keeptet, pt1lc, pt2lc
C  ================================================================ 
c
      If (level(i) .eq. 0) Then
         verf(i) = 0
         Return
      End if   

c ... define pt1 for tetrahedron i
      call  Uncode (history(1,level(i),i),
     &              history(2,level(i),i), pt1lc)
      
      call  Uncode (history(3,level(i),i),
     &              history(4,level(i),i), pt2lc)
      

c ... define tetrahedron ii
      ii = IEF(pt1lc,i)

      l  = edp(pt1lc, pt2lc)
      pt3 = fpe(1, 7 - l)
      pt4 = fpe(2, 7 - l)

      verf (ii) = max(verf(ii), 1)

c ... define conflict nodes
      v = IPE(pt2lc, i)
     
      flagcoarse = .false.
      deltet = 0 
      keeptet = 0


c ... control of correspondence of level, pbisE and history
c     of tetrahedra i and ii 
      If (level(i) .ne. level(ii)) Then
         If (pconfP(pt2lc, i) .eq. 1) Return
         Go to 20  
      End if
    
c ... define flag
      If (pconfP(pt2lc, i)  .eq. 0) Then
        flagcoarse = .true.
      End if
      
      If (realP(v)  .eq. 1) Then
        realP(v) = 0
      End if
      
c ... define deltetl, keeptet
      If (i .gt. ii) Then
        deltet = i
        keeptet = ii
        pt1 = pt2lc
        pt2 = pt1lc
       Else
        deltet = ii
        keeptet = i
        pt1 = pt1lc
        pt2 = pt2lc
      End if

c ... delete tet deltet
      realE(deltet) = 0     

c ... define IPE
      IPE(pt2, keeptet) = IPE(pt2, deltet)

c ... define IEF
      IEF(pt1, keeptet) = IEF(pt1, deltet)
       
      j = IEF(pt1, deltet)
      If (j .gt. 0) Then
        Do kk = 1, 4
          If (IEF(kk, j) .eq. deltet) Then
           IEF(kk, j) = keeptet
           Go to 1
          End if  
        End do
1       Continue
      End if

      j = IEF(pt3, deltet)
      If (j .gt. 0) Then
        Do kk = 1, 4
          If (IEF(kk, j) .eq. deltet) Then
           IEF(kk, j) = keeptet
           Go to 2
          End if  
        End do
2      Continue
      End if 

      j = IEF(pt4, deltet)
      If (j .gt. 0) Then
        Do kk = 1, 4
          If (IEF(kk, j) .eq. deltet) Then
           IEF(kk, j) = keeptet
           Go to 3
          End if  
        End do
3       Continue
      End if 
      

c ... define pbisE
      flag(keeptet) = 0
      If (twoedges_on_face(pbisE(1, keeptet),
     &                     pbisE(2, keeptet),fpe) .and.
     &    twoedges_on_face(pbisE(1, keeptet),
     &                     pbisE(3, keeptet),fpe) .and. 
     &    .not.(twoedges_on_face(pbisE(2, keeptet),
     &          pbisE(3, keeptet),fpe)))
     &      flag(keeptet) = 1
      
      pbisE(2, keeptet) = pbisE(1, keeptet)
      pbisE(3, keeptet) = pbisE(1, deltet)
      pbisE(1, keeptet) = edp(pt1, pt2)
      call Order_marked_edges (pbisE(1,keeptet), fpe)

c ... define pconfP
      pconfP(pt2, keeptet) = pconfP(pt2, deltet)
      pconfP(pt3, keeptet) = pconfP(pt3, deltet)
      pconfP(pt4, keeptet) = pconfP(pt4, deltet)

c ... define verf
      verf(keeptet) = max(max(verf(keeptet), verf(deltet)) - 1,
     &              pconfP(pt1,keeptet), 
     &              pconfP(pt2,keeptet), 
     &              pconfP(pt3,keeptet), 
     &              pconfP(pt4,keeptet)) 

c ... define level 
      level(keeptet) = level(keeptet) - 1

c ... modify  listtet     
      k = IPE(pt2,keeptet)
      m = listtet(1, k)
      Do l = 1, m
        If (listtet(l + 1, k) .eq. deltet) Then
          listtet(l + 1, k) = keeptet
          Go to 10
        End if
      End do
      iERR = 1006  !??
      call errMes(iERR, 'Unteil','Error in list of tets')
10    Continue

      If (.not. flagcoarse) Return

20    Continue

      m = listtet(1, v)
      Do l = 1, m
        jj = listtet(l + 1, v)
        If ((jj .eq. deltet) .or.
     &      (jj .eq. keeptet)) Go to 30
        If (realE(jj) .eq. 0)  Go to 30

        Do kk = 1, 4
          If (IPE(kk, jj) .eq. v) Then
             verf (jj) = max(verf(jj), 1)
             pconfP(kk, jj) = 1
             Go to 30
          End if
        End do
        iERR = 1006  !??
        call errMes(iERR, 'Unteil','Error in list of tets')

30      Continue


      End do 

      Return
      End
     
