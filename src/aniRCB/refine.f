      Subroutine  LocalRefine  (
C ================================================================    
     &        nP, nPmax, nF, nFmax, nE, nEmax,
     &        XYP, IPE, IPF, lbE, labelF,
     &        RefineRule, ilevel, maxlevel, history,
     &        listtet, tetpmax, RefineRuleData,
     &        MaxWi, iW,
     &        iPrinE, iERR)
C ================================================================    

      Implicit none

c ... user procedure
      external  RefineRule

c ... standard mesh arrays 
c ... number of points, tetrahedra and boundary edges
      Integer   nP, nPmax, nE, nEmax, nF, nFmax

c ... maxlevel - maximum number of level bisection
      Integer   maxlevel

c ... coordinates of mesh points 
      Double precision   XYP(3,nPmax)

c ... connectivity table for tetrahedra and tet labels
      Integer   IPE(4,nEmax), lbE(nEmax)

c ... connectivity table for boundary edges, and edge labels
      Integer   IPF(3,nFmax), labelF(nFmax)

c ... level of bisection
      Integer   ilevel

c ... history of bisection
      Logical   history(4,maxlevel, nEmax)

c ... list of tetrahedra
      Integer   tetpmax
      Integer   listtet(tetpmax+1,*)
      
c ... user data to be passed to RefineRule
      Double precision RefineRuleData(*)

c ... work memory
      Integer   MaxWi
      Integer   iW(MaxWi)
     
      Integer   iPrinE, iERR

c ... local point
      Integer   iP_pbisE, iP_flag, iP_level 
      Integer   iP_fpe, ip_edp
      Integer   ip_IEF, ip_verf, ip_fconfE
      Integer   ip_liW, MaxWil
C ================================================================    


c ... ip_IEF, ip_verf, ip_pbisE, ip_fconfE
      iP_pbisE   =                   1
      iP_flag    = iP_pbisE        + 3*nEmax
      iP_level   = iP_flag         + nEmax
      iP_fpe     = iP_level        + nEmax
      ip_edp     = iP_fpe          + 24
      ip_IEF     = iP_edp          + 16
      ip_verf    = ip_IEF          + 4*nEmax
      ip_fconfE  = ip_verf         + nEmax
      iP_liW     = ip_fconfE       + 6*nEmax
      
      If (MaxWi .le. 15*nEmax+41) Then
          iERR = 1010
          call errMes(iERR, 'LocalRefine','MaxWi  is too small')
      End if
      MaxWil = MaxWi - 15*nEmax-41


c ... define ip_IEF, ip_pbisE
      call OrientandMakeNeighbors (
     &           nP, nPmax, nF, nFmax, nE, nEmax,
     &           IPE, iW(ip_IEF), IPF, labelF,
     &           listtet, tetpmax,
     &           iERR)
     
c ... define ip_verf
      call RefineRule(nE, IPE, XYP, iW(ip_verf), ilevel, RefineRuleData)
      
c ... bisection
      call Tetr (nP, nPmax, nE, nEmax, IPE, iW(ip_IEF), XYP, 
     &      iW(ip_pbisE),  lbE, iW(ip_verf), 
     &      iW(ip_fconfE), iW(ip_flag),
     &      iW(ip_level), maxlevel, history,
     &      listtet, tetpmax,
     &      iW(ip_fpe), iW(ip_edp), iERR)

c ... define IPF
      call RestoreBndEdges (nF, nFmax, nE, IPE, iW(ip_IEF),
     &                      IPF, labelF, iERR)
                  
      Return
      End

C  ================================================================ 
c  TETR produces (local) refinement of a given set of tetrahedra 
c  the algorithm is a version of method by L. Arnold, 
c  A. Mukherjee and L. Pouly 
C  ================================================================ 
       Subroutine Tetr (
C  ================================================================ 
     &           nP, nPmax, nE, nEmax,
     &           IPE, IEF, XYP, 
     &           pbisE,  lbE, verf, fconfE, flag, 
     &           level, maxlevel, history,
     &           listtet, tetpmax,
     &           fpe, edp, iERR)
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
c
c  XYP(l, i)    : l-th coordinate of the i-th node
c  verf(k)      : flag whether to bisect tetrahedron k or not
c  IPE(l, k)    : No. of l-th vertex of tet k 
c  IEF(l, k)    : neighbour of tet k opposite to vertex l
c  fconfE(l, k) : flag whether there's a "conflict" on edge l
c  pbisE(k)     : No. of the refinement and marked edges
c                 of tetrahedra k
C  ================================================================ 
      implicit none

      Double precision XYP(3,*)
      Integer          nP, nPmax, nE, nEmax
      Integer          maxlevel
      Integer          IPE(4,*), IEF(4,*)
      Integer          pbisE(3,*)
      Integer          lbE(*)
      Integer          verf(*)
      Integer          fconfE(6,*)
      Integer          flag(*)
      Integer          level(*)
      Integer          iERR
      Logical          history(4,maxlevel,*)
      Integer          tetpmax
      Integer          listtet(tetpmax + 1, *)

      Integer          minlevel
                     
      Integer          fpe(4,6), edp(4,4)
      Integer          i, j
C  ================================================================ 
c ... initialization:

      Do i = 1, nE
        Do j = 1, 6
          fconfE(j, i) = 0
        End do
      End do

c ... bisect elemenes until there is  no conflict:
1     Continue

      Do i = 1, nE
        If (verf(i) .ne. 0)
     &    call Teil_tet(i, nP, nPmax, nE, nEmax, IPE, IEF, XYP,
     &           fpe, edp,  pbisE,  lbE, verf, fconfE, flag,
     &           level, maxlevel, history,
     &           listtet, tetpmax,
     &           iERR)
      End do

      Do i = 1, nE
        If (verf(i) .ne. 0)  Go to 1
      End do

c
      Return
      End
C  ================================================================ 

C  ================================================================ 
c  TEIL_TET divides the tetrahedron i into two new ones
c
c  INPUT:      i - no. of tetrahedron
C  ================================================================ 
      Subroutine Teil_tet(
C  ================================================================ 
     &           i, nP, nPmax, nE, nEmax, IPE, IEF, XYP, 
     &           fpe, edp, pbisE,  lbE, verf, fconfE, flag,         
     &           level, maxlevel, history,
     &           listtet, tetpmax,
     &           iERR)
C  ================================================================ 
      implicit none

      Integer           i
      Integer           nP, nPmax, nE, nEmax
      Integer           maxlevel
      Integer           IPE(4,*), IEF(4,*)
      Double precision  XYP(3,*) 
      Integer           fpe(4,6), edp(4,4)
      Integer           pbisE(3,*)
      Integer           flag(*)
      Integer           lbE(*)
      Integer           verf(*)
      Integer           fconfE(6,*)
      Integer           level(*)
      Integer           tetpmax
      Integer           listtet(tetpmax + 1, *)
      Integer           iERR
      Logical           history(4,maxlevel, *)
     

      Integer           ii, l, j
      Integer           pt1, pt2, pt3, pt4
      Integer           pbisEfater(3), flagfater
      Integer           flagconf
      Integer           k, npt,  newverf
C  ================================================================ 

c ... l - refinement edge
      l = pbisE(1,i)
      flagconf = fconfE(l, i)
      
c ... one more tetrahedron: ii - No. of new tetrahedron
      If (nE .ge. nEmax) Then
        iERR = 1006
        call errMes(iERR, 'teil','nEmax is too small')
      End if  

c ... introduce new tetrahedron
      nE = nE+1
      lbE(nE) = lbE(i)
      ii = nE

c ... define level 
      level(i)  = level(i) + 1
      level(ii) = level(i) 
      

c ... define the value of  some variables for the
c     new tetrahedron and change values for the old one

c ... points tetrahedron
      pt1 = fpe(1, l)
      pt2 = fpe(2, l)
      pt3 = fpe(1, 7 - l)
      pt4 = fpe(2, 7 - l)

c ... confl tetrahedron
      fconfE(edp(pt2, pt3), ii) = fconfE(edp(pt2, pt3), i)
      fconfE(edp(pt2, pt4), ii) = fconfE(edp(pt2, pt4), i)
      fconfE(edp(pt3, pt4), ii) = fconfE(edp(pt3, pt4), i)

      fconfE(edp(pt1, pt2), ii) = 0 
      fconfE(edp(pt1, pt3), ii) = 0
      fconfE(edp(pt1 ,pt4), ii) = 0

      fconfE(edp(pt2, pt1), i) = 0 
      fconfE(edp(pt2, pt3), i) = 0 
      fconfE(edp(pt2, pt4), i) = 0

c ... IEF
      IEF(pt1, ii) = IEF(pt1, i)
      IEF(pt1, i)  = 1
      IEF(pt2, ii) = 1
      IEF(pt3, ii) = IEF(pt3, i)
      IEF(pt4, ii) = IEF(pt4, i)

c ... pbisE
      pbisEfater(1) = pbisE(1,i)
      pbisEfater(2) = pbisE(2,i)
      pbisEfater(3) = pbisE(3,i)
      flagfater     = flag(i)

      call marked_edges (
     &      pbisEfater,  flagfater,
     &      pbisE(1,i), flag(i),
     &      pbisE(1,ii), flag(ii), fpe)


c ...is there a "conflict" for i or ii:
      newverf = max(verf(i) - 1, 0)
      
      verf (i) = newverf
      If (fconfE(edp(pt1, pt3), i) .ne.0)
     &        verf (i) = max(verf(i), 1)
      If (fconfE(edp(pt1, pt4), i) .ne.0)
     &        verf (i) = max(verf(i), 1)
      If (fconfE(edp(pt3, pt4), i) .ne.0)
     &        verf (i) = max(verf(i), 1)

      verf (ii) = newverf
      If (fconfE(edp(pt2, pt3), ii) .ne.0)
     &        verf (ii) = max(verf(ii), 1)
      If (fconfE(edp(pt2, pt4), ii) .ne.0)
     &        verf (ii) = max(verf(ii), 1)
      If (fconfE(edp(pt3, pt4), ii) .ne.0)
     &        verf (ii) = max(verf(ii), 1)

c ... IPE
      IPE(pt2, ii) = IPE(pt2, i) 
      IPE(pt3, ii) = IPE(pt3, i) 
      IPE(pt4, ii) = IPE(pt4, i) 

c ... listtet
      k = IPE(pt2, i)
      Do j = 1, listtet(1, k)
        If (listtet(j+1, k) .eq. i) Then
          listtet(j+1, k) = ii
          Go to 2
        End if  
      End do
2     Continue

      k = IPE(pt3, i)
      listtet(1, k) = listtet(1, k) + 1
      j = listtet(1, k)
      listtet(j + 1, k) = ii

      k = IPE(pt4, i)
      listtet(1, k) = listtet(1, k) + 1
      j = listtet(1, k)
      listtet(j + 1, k) = ii

c ... adjust all variables in the elements meeting
c ... at the refinement edge:

      If (flagconf .ne. 0) Then
        npt =  flagconf
        listtet(1, npt) = listtet(1, npt) + 2
        If (listtet(1, npt) .ge. tetpmax) Then
          iERR = 1007 
          call errMes(iERR, 'Teil_tet','tetpmax is too small')
        End if  
        j = listtet(1, npt) - 1
        listtet(j + 1, npt) = i
        listtet(j + 2, npt) = ii
       Else
        call randk (IPE(pt1,i),IPE(pt2,i),
     &       nP,nPmax,XYP, iERR)
        npt = nP
        call conf (nE, i, ii, IPE(pt1,i),IPE(pt2,i), npt,
     &        IPE, fconfE, verf, edp, listtet, tetpmax)
        listtet(1, npt) = 2
        listtet(2, npt) = i
        listtet(3, npt) = ii
      End if

      IPE(pt2, i)  = npt
      IPE(pt1, ii) = npt

c ... define history
      call  Code (pt1, history(1,level(i), i),
     &                 history(2,level(i), i)) 
      call  Code (pt2, history(3,level(i), i),
     &                 history(4,level(i), i)) 
        
      call  Code (pt2, history(1,level(ii), ii),
     &                 history(2,level(ii), ii))
      call  Code (pt1, history(3,level(ii), ii),
     &                 history(4,level(ii), ii))


      Return
      End

C  ================================================================ 
c  RANDK the refinement adge is now a boundary edge
C  ================================================================ 
      Subroutine  randk (
C  ================================================================ 
     &       pl1, pl2, nP, nPmax, XYP, iERR)
C  ================================================================ 
      Implicit none

      Integer           pl1, pl2
      Integer           nP,nPmax
      Double precision  XYP(3,*)
      Integer           iERR

C  ================================================================
c ...add a new node:
      If (nP.ge.nPmax) Then
         iERR = 1003
         call errMes(iERR, 'randk','nPmax is too small')
      End if  

      nP = nP+1
      XYP(1, nP) = (XYP(1, pl1) +  XYP(1, pl2))/2.d0
      XYP(2, nP) = (XYP(2, pl1) +  XYP(2, pl2))/2.d0
      XYP(3, nP) = (XYP(3, pl1) +  XYP(3, pl2))/2.d0

      Return
      End

C  ================================================================ 
c  CONF  poses "conflict" on tetrahedra sharing the edge,
c       on which tetrahedron i is bisected 
C  ================================================================ 
      Subroutine  conf  (
C  ================================================================ 
     &       nE, i, ii, pt1, pt2, npl,
     &       IPE, fconfE, verf, edp, listtet, tetpmax)
C  ================================================================ 
      Implicit none

      Integer           nE, i, ii
      Integer           pt1, pt2, npl
      Integer           IPE(4,*), fconfE(6,*), verf(*)
      Integer           edp(4,4)
      Integer           tetpmax
      Integer           listtet(tetpmax + 1, *)

      Integer           j1, j2, k1, k2
      Integer           kk, kk1, kk2
C  ================================================================

      Do j1 = 1, listtet(1, pt1)

        k1 = listtet(j1 + 1, pt1)
        If ((k1 .eq. i) .or. (k1 .eq. ii)) Go to 1

        Do j2 = 1, listtet(1, pt2) 

          k2 = listtet(j2 + 1, pt2)

          If (k1 .eq. k2) Then
            Do  kk = 1, 4
              If (IPE(kk, k1) .eq. pt1)  kk1 = kk
              If (IPE(kk, k1) .eq. pt2)  kk2 = kk
            End do  
      
            verf(k1) = max(verf(k1), 1)
            fconfE(edp(kk1, kk2), k1) = npl
          End if
        
        End do
1       Continue
      End do

         
      Return 
      End

C  ================================================================ 
c  Marked of edges 
C  ================================================================ 
        Subroutine marked_edges (
C  ================================================================ 
     &      markedg,  flag,
     &      markedg1, flag1,
     &      markedg2, flag2,
     &      fpe)
C  ================================================================ 
        implicit none
c
        External          twoedges_on_face
        Logical           twoedges_on_face
c
        Integer           markedg(3),   flag
        Integer           markedg1(3),  flag1
        Integer           markedg2(3),  flag2
        Integer           fpe(4,6)
        
        Integer           i, j
C  ================================================================ 

         markedg1(1)  =  markedg(2)
         markedg2(1) =  markedg(3)

         flag1 = 0
         flag2 = 0

 
         If (flag .eq. 0) Then

           If (twoedges_on_face(markedg(1),markedg(2),fpe) .and.
     &         twoedges_on_face(markedg(1),markedg(3),fpe) .and.
     &         twoedges_on_face(markedg(2),markedg(3),fpe)) Then
              flag1 = 1
              flag2 = 1
           End if

           Call Edges_opozite_vertex (markedg1(1), fpe(2,markedg(1)),
     &                                markedg1(2), markedg1(3), fpe)     

           Call Edges_opozite_vertex (markedg2(1), fpe(1,markedg(1)),
     &                                markedg2(2), markedg2(3), fpe)     

         Else
           
           Call Edge_with_vertex (fpe(1,markedg(1)), markedg(1),
     &                                markedg(2),  markedg1(2), fpe)      
           markedg1(3) = 7 - markedg1(2)

           Call Edge_with_vertex (fpe(2,markedg(1)), markedg(1),
     &                                markedg(3),  markedg2(2), fpe)      

           markedg2(3) = 7 - markedg2(2)
        End if

        Call Order_marked_edges (markedg1, fpe)
        Call Order_marked_edges (markedg2, fpe)
     
     
        Return
        End

C  ================================================================ 
c  chick two edges on face
C  ================================================================ 
        Logical function twoedges_on_face (
C  ================================================================ 
     &       l, k, fpe)
C  ================================================================ 
        implicit none
c
        Integer           l, k, fpe(4,6)
        
        Integer           i, j
C  ================================================================ 

        twoedges_on_face = .false.
     
        Do i = 3, 4
          Do j = 3, 4
            If (fpe(i,k) .eq. fpe(j,l)) Then
               twoedges_on_face = .true.
               Return
            End if
          End do
        End do  
     
     
        Return
        End

C  ================================================================ 
c  define  two edges opozite the vertex and 
c  different the given edge
C  ================================================================
        Subroutine Edges_opozite_vertex (
C  ================================================================ 
     &       l1, pt, l2, l3, fpe)
C  ================================================================ 
        implicit none
c
        Integer           l1, pt, fpe(4,6)
        Integer           l2, l3
        
        Integer           fl, i
C  ================================================================ 

        fl = 0
        Do i = 1, 6
          If (i. eq. l1) go to 1
          If (fpe(1,i). eq. pt) go to 1
          If (fpe(2,i). eq. pt) go to 1
          If (fl. eq. 0) Then
            l2 = i
            fl = 1 
           Else  
            l3 = i
          End if
1         Continue
        End do
          
        Return
        End


C  ================================================================ 
c  define  two edge near the vertex and different the given edges
C  ================================================================
        Subroutine Edge_with_vertex ( 
C  ================================================================ 
     &        pt, l1, l2, l3, fpe)
C  ================================================================ 
        implicit none
c
        Integer           pt, l1, l2,fpe(4,6)
        Integer           l3
        
        Integer           i
C  ================================================================ 

        Do i = 1, 6
          If (i. eq. l1) go to 1
          If (i. eq. l2) go to 1
          If ((fpe(1,i). eq. pt) .or. (fpe(2,i). eq. pt)) Then
            l3 = i
            Return
          End if
1         Continue
        End do
          
        Return
        End

C  ================================================================ 
c  order the marked edges
C  ================================================================
        Subroutine Order_marked_edges ( 
C  ================================================================ 
     &        markedg, fpe)
C  ================================================================ 
        implicit none
c
        Integer           markedg(3), fpe(4,6)
        
        Integer           pt, k 
C  ================================================================ 

        pt = fpe(1, markedg(1))

        If ((fpe(1, markedg(3)) .eq. pt) .or. 
     &      (fpe(2, markedg(3)) .eq. pt)) Then
           k = markedg(3)
           markedg(3) = markedg(2)
           markedg(2) = k
        End if

        Return
        End
 
