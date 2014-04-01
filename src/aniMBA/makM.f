C ================================================================
      Subroutine makM(
C ================================================================
c group (M)
     &      nP, nF, MaxF, nE, MaxE, 
     &      XYP, IPF, IPE, 
     &      ICP, IPP, IEP, IFE, IEE,
     &      IHolP, IHolF, IHolE,
     &      IEPw, nEPw,
     &      status,
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (ERR)
     &      iERR)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'status.fd'
C ================================================================
C Routine analyses the initial geometry and creates auxiliary 
C cross-refrences.
C
C Pre-conditions:  1. connectivity structure {IPE(5, *), XYP(3, *)}
C
C Post-conditions: 1. additional mesh structures, IEP, IFE, IEE
C                  2. coloring mesh point according to given
C                     surface and volume colors, and lists of
C                     fixed points, triangles, and tetrahedra.
C ================================================================
C group (M)
C     Integer MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(5, *), IPF(4, *)

      Integer ICP(*), IPP(*), IEP(*)
      Integer IFE(4, *), IEE(4, *)
      Integer IHolP(*), IHolF(*), IHolE(*)
      Integer IEPw(*), nEPw(*)

      Integer status

c group (Dev)
      Integer nPv, nFv, nEv, IPV(*), IFV(*), IEV(*)

C ================================================================
C group (Functions)
      Logical cmpE, cmpR, cmpP, crossPoint

C group (Local variables)
      Integer ip(5), IPFs(4)
      Logical flagFBF, flagBND, ifXnode
      Real*8  calVol

      Integer Mlist(2, MaxS-iMface+1)

C ================================================================
      iERR = 0
 
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


      IHolP(1) = 0
      IHolF(1) = 0
      IHolE(1) = 0

      Do n = 1, nP
         ICP(n) = 0
      End do


C ... create an auxiliary structure
      Call backReferences(nP, nE, 4, 5, IPE, nEPw, IEPw)

C ... create IEE & IEP
      Do n = 1, MaxE
         Do i = 1, 4
            IEE(i, n) = 0
            IFE(i, n) = 0
         End do
      End do

      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            IEP(ip1) = n

            If(cmpE(ip1, ip2, ip3, IEPw, nEPw, n, iE2)) Then
               IEE(i1, n) = iE2
           End if
         End do
      End do


c ... create an auxiliary structure
      Call backReferences(nP, nF, 3, 4, IPF, nEPw, IEPw)


c ... create IFE: basic, fictitious, material
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            If(cmpE(ip1, ip2, ip3, IEPw, nEPw, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do


c ... create missing boundaries
      Do n = 1, nF
         If(IPF(4, n).GE.iVface) Call errMes(1011, 'makM',
     &                  'reserved boundary identificator is used')
      End do

      flagFBF = ifXnode(status, ANIFixBoundaryFaces)
      Do n = 1, nE
         Do i1 = 1, 4
            If(IEE(i1, n).EQ.0 .AND. IFE(i1, n).EQ.0) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               IPFs(1) = IPE(i1, n)
               IPFs(2) = IPE(i2, n)
               IPFs(3) = IPE(i3, n)
               IPFs(4) = iVface

               Call facAdd(iF, nF, MaxF, IHolF)
               Call facUpd(1, IPF, iF, IPFs)

               IFE(i1, n) = iF

               If(flagFBF) Then
                  Do i = 1, 3
                    iPt = IPFs(i)
                    Call addXnode(ICP(iPt), jTnode)
                  End do
               End if
            End if
         End do
      End do


c ... create material boundaries
      k = 0
      kMax = MaxS - iMface
      Do n = 1, nE
         Do i1 = 1, 4
            iFt = IFE(i1, n)
            iEt = IEE(i1, n)
            If(iFt.EQ.0 .AND. iEt.NE.0) Then
               mat1 = IPE(5, n)
               mat2 = IPE(5, iEt)
               If(mat1.NE.mat2) Then
c  ...  search for this pair in the list of material interfaces
                  Do i = 1, k
                     If((Mlist(1, i).EQ.mat1 .AND.
     &                   Mlist(2, i).EQ.mat2) .OR.
     &                  (Mlist(1, i).EQ.mat2 .AND.
     &                   Mlist(2, i).EQ.mat1)) Then
                        iCface = iMface + i - 1
                        Goto 1
                     End if
                  End do

c  ...  make a new material interface
                  k = k + 1
                  If(k.GT.kMax) Call errMes(1010,
     &              'makM', 'not enough memory for material faces')
                  Mlist(1, k) = mat1
                  Mlist(2, k) = mat2
                  iCface = iMface + k - 1

  1               i2 = ip(i1 + 1)
                  i3 = ip(i2 + 1)

                  IPFs(1) = IPE(i1, n)
                  IPFs(2) = IPE(i2, n)
                  IPFs(3) = IPE(i3, n)
                  IPFs(4) = iCface

                  Call facAdd(iF, nF, MaxF, IHolF)
                  Call facUpd(1, IPF, iF, IPFs)

                  IFE(i1, n) = iF

                  Do i = 1, 4
                     If(IEE(i, iEt).EQ.n) IFE(i, iEt) = iF
                  End do
               End if
            End if
         End do
      End do


c ... color the points (boundary points)
c ... new auxiliary structure accumulates all faces
      Call backReferences(nP, nF, 3, 4, IPF, nEPw, IEPw)

      Do n = 1, nE
         Do i1 = 1, 4
            iF = IFE(i1, n)
            iE = IEE(i1, n)
            If(iF.NE.0) Then
               i2 = i1

               Do k = 1, 3
                  iP1 = IPE(i2, n)

                  Call addXnode(ICP(iP1), jSnode)
                  If(iE.EQ.0) Call addXnode(ICP(iP1), jBnode)

                  i2 = ip(i2 + 1)
               End do
            End if
         End do
      End do


c ... coloring the points (edge points)
      Do n = 1, nP
         If(.NOT.ifXnode(ICP(n), jBnode)) Call addXnode(ICP(n), jInode)
         If(cmpP(n, IEPw, nEPw, IPF))     Call addXnode(ICP(n), jRnode)
         If(cmpR(n, IEPw, nEPw, XYP, IPF)) 
     &                                    Call addXnode(ICP(n), jRnode)
      End do


c ... coloring the points (cross points)
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nEPw(n)

         If(ifXnode(ICP(n), jSnode)) Then
            If(crossPoint(XYP, i2 - i1 + 1, IEPw(i1), n, IPF)) Then
               Call addXnode(ICP(n), jVnode)
            End if
         End if
      End do


c ... coloring the points (fix vertices)
      Do n = 1, nPv
         iP1 = IPV(n)
         Call addXnode(ICP(iP1), jVnode)
      End do


c ... color the points of fix faces and elements
      Do n = 1, nFv
         iF = IFV(n)
         Do i = 1, 3
            iP1 = IPF(i, iF)
            Call addXnode(ICP(iP1), jTnode)
         End do
      End do

      Do n = 1, nEv
         iE = IEv(n)
         Do i = 1, 4
            iP1 = IPE(i, iE)
            Call addXnode(ICP(iP1), jTnode)
         End do
      End do


c ... color the points (fixed boundary points)
      If(ifXnode(status, ANIFixSurfacePoints)) Then
         Do n = 1, nF
            Do i = 1, 3
               iP1 = IPF(i, n)
               Call addXnode(ICP(iP1), jVnode)
            End do
         End do
      End if


c ... marking points with a multi-connected superelement
      Do n = 1, nP
         nEPw(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, 4
            i1 = IPE(i, n)
            nEPw(i1) = nEPw(i1) + 1
         End do
      End do

      Do n = 1, nP
         If(ifXnode(ICP(n), jSnode)) Then
            Call makSP(n, IEP, IPE, IEE, MaxE, lE, IEPw)
            If(nEPw(n).NE.lE) Then
               Call addXnode(status, ANIMultiConnectedGeometry)
               Call addXnode(ICP(n), jVnode)
            End if
         End if
      End do


c ... coloring the T-points if any (case of inner boundaries)
c ... array IPF(4, *) is overloaded to mark only boundary edges
      Do n = 1, nE
         Do i = 1, 4
            iFt = IFE(i, n)
            If(iFt.GT.0) IPF(4, iFt) = -IPF(4, iFt)
         End do
      End do

      Do n = 1, nP
         IEPw(n) = 0
      End do

      Do 100 n = 1, nF
         If(IPF(4, n).GT.0) Goto 100
         IPF(4, n) = -IPF(4, n)

         flagBND = .FALSE.
         Do i = 1, 3
            If(IPP(IPF(i, n)).EQ.0) flagBND = .TRUE.
         End do

         If(flagBND) Then
            Do i = 1, 3
               IEPw(IPF(i, n)) = 1
            End do
         End if
 100  Continue

      Do n = 1, nP
         If(IPP(n).GT.0) Then 
            Call addXnode(ICP(n), jTnode)
            If(IEPw(n).EQ.0) Then 
               Call addXnode(ICP(n), jInode)
               Call delXnode(ICP(n), jBnode)
            End if
         End if
      End do


c ... delete isolated points
      Do n = 1, nP
         IEPw(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, 4
            IEPw(IPE(i, n)) = 1
         End do
      End do

      Do n = 1, nP
         If(IEPw(n).EQ.0) Then
            Call pntDel(n, nP, ICP, IHolP)
         End if
      End do

c ... delete isolated faces
      Do n = 1, nF
         IEPw(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, 4
            iFt = IFE(i, n)
            If(iFt.GT.0) IEPw(iFt) = 1
         End do
      End do

      mF = nF
      Do n = 1, mF
         If(IEPw(n).EQ.0) Then
            Call facDel(n, nF, IPF, IHolF)

            Do k = 1, nFv
               If(IFV(k).EQ.n) Then
                  IFV(k) = IFV(nFv)
                  nFv = nFv - 1
                  Goto 500
               End if
            End do
 500     Continue
         End if
      End do


c ... order boundary faces clockwise looking from inside domain
      Do n = 1, nE
         Do i1 = 1, 4
            iF = IFE(i1, n)
            If(iF.GT.0) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)
               i4 = ip(i3 + 1)

               ip4 = IPE(i4, n)

               ip1 = IPF(1, iF)
               ip2 = IPF(2, iF)
               ip3 = IPF(3, iF)

               If(calVol(XYP(1, ip1), XYP(1, ip2), 
     &                   XYP(1, ip3), XYP(1, ip4)).GT.0D0) Then
                  Call swapii(IPF(1, iF), IPF(2, iF))   
               End if
            End if
         End do
      End do

      Return
      End



C ================================================================
      Subroutine updM(
C ================================================================
c group (M)
     &           nP, nF, nE, 
     &           XYP, IPF, IPE,
     &           ICP, IPP, IFE, IEE,
     &           IHolP, IHolF, IHolE,
     &           status,
c group (Dev)
     &           nPv, nFv, nEv, IPV, IFV, IEV,
c group (Q)
     &           HesP, qE,
c group (W)
     &           IPw)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'status.fd'
C ================================================================
C Routine removes holes (deleted elements) from the mesh data
C structures. The optimal complexity algorithms are implemented.
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(5, *)
      Real*8  XYP(3, *)

      Integer ICP(*), IPP(*), IFE(4, *), IEE(4, *)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer status

c group (Dev)
      Integer nPv, nFv, nEv, IPV(*), IFV(*), IEV(*)

C group (Q)
      Real*8  HesP(6, *), qE(*)

c group (W)
      Integer IPw(*)
C ================================================================
C group (Functions)
      Logical ifXnode, flagDTF

C ================================================================
      flagDTF = ifXnode(status, ANIDeleteTemporaryFaces)


c ... recovering V-, VB- and I-points from TV-, TVB, and TB-points
      Do n = 1, nP
        Call delXnode(ICP(n), jTnode)
      End do
 

c ... delete references to material or fictitious faces
      lE = nE + IHolE(1)
      Do n = 1, lE
         Do i = 1, 4
            iF = IFE(i, n)
            If(iF.GT.0) Then
               If(IPF(4, iF).EQ.iVface .AND. flagDTF .OR.
     &            IPF(4, iF).GE.iMface) IFE(i, n) = iFface
            End if
         End do
      End do


c ... delete all material faces
      lF = nF + IHolF(1)
      Do n = 1, lF
         If(IPF(4, n).GE.iMface) Then
            Call facDel(n, nF, IPF, IHolF) 
        End if
      End do


      lF = nF + IHolF(1)
      If(flagDTF) Then
c  ...  delete all fictitious faces
        Do n = 1, lF
           If(IPF(4, n).EQ.iVface) Then
              Call facDel(n, nF, IPF, IHolF)
           End if
        End do
      Else
c  ...  change the color of fictitious faces if it's possible
        icfree = iVface

        Do 1 ic = 1, iVface - 1
           Do n = 1, lF
              If(IPF(4, n).EQ.ic) Goto 1
           End do
           icfree = ic
           Goto 2
 1      Continue

 2      Continue
        Do n = 1, lF
           If(IPF(4, n).EQ.iVface) IPF(4, n) = icfree
        End do
      End if


      nP = nP + IHolP(1)
      nF = nF + IHolF(1)
      nE = nE + IHolE(1)


c ... fill-in holes in the list of mesh points
      Do n = 1, nP
         IPw(n) = 0
      End do

      lP = IHolP(1)
      Do n = 1, lP
         iP = IHolP(n + 1)
         IPw(iP) = -1
      End do

      mP = 0
      Do n = 1, nP
         If(IPw(n).EQ.0) Then
            mP = mP + 1
            IPw(n) = mP

            Do i = 1, 3
               XYP(i, mP) = XYP(i, n)
            End do
         
            ICP(mP) = ICP(n)
            IPP(mP) = IPP(n)

            Do i = 1, 6
               HesP(i, mP) = HesP(i, n)
            End do
         End if
      End do

      Do n = mP + 1, nP
         ICP(n) = 0
         IPP(n) = 0
      End do

      Do n = 1, nPv
         IPV(n) = IPw(IPV(n))
      End do

      Do n = 1, nF
         If(IPF(1, n).GT.0) Then
            Do i = 1, 3
               IPF(i, n) = IPw(IPF(i, n)) 
            End do
         End if
      End do

      Do n = 1, nE
         If(IPE(1, n).GT.0) Then
            Do i = 1, 4
               IPE(i, n) = IPw(IPE(i, n))
            End do
         End if
      End do

      nP = mP


c ... fill-in holes in the list of mesh faces
      lF = IHolF(1)
      Do 200 n = 1, lF
         iF = IHolF(n + 1)

         Do m = nF, iF + 1, -1
            If(IPF(1, m).GT.0) Then
               kF = m
               Goto 20
            End if
         End do
         Goto 200

 20      Do i = 1, 4
            IPF(i, iF) = IPF(i, kF)
         End do


         Do k = 1, nFv
            If(IFV(k).EQ.kF) IFV(k) = iF
         End do


C  ...  auxiliary structures
         Do k = 1, nE
            Do i = 1, 4
               If(IFE(i, k).EQ.kF) IFE(i, k) = iF
            End do
         End do

         IPF(1, kF) = 0
 200     nF = nF - 1


c ... fill-in holes in the list of mesh elements
      lE = IHolE(1)
      Do 300 n = 1, lE
         iE = IHolE(n + 1)

         Do m = nE, iE + 1, -1
            If(IPE(1, m).GT.0) Then
               kE = m
               Goto 30
            End if
         End do
         Goto 300

 30      Do i = 1, 5
            IPE(i, iE) = IPE(i, kE)
         End do

         qE(iE) = qE(kE)

         Do k = 1, nEv
            If(IEV(k).EQ.kE) IEV(k) = iE
         End do


C  ...  auxiliary structures
         Do i = 1, 4
            IFE(i, iE) = IFE(i, kE)
            IEE(i, iE) = IEE(i, kE)
            iEt = IEE(i, iE)
            If(iEt.GT.0) Then
              Do j = 1, 4
                 If(IEE(j, iEt).EQ.kE) IEE(j, iEt) = iE
              End do
            End if
         End do

         IPE(1, kE) = 0
 300     nE = nE - 1

      Return
      End



C ================================================================
      Logical Function cmpE(i1, i2, i3, IEP, nEP, iE1, iE2)
C ================================================================
C cmpE = TRUE if iE2 != iE1 and iE2 = {i1, i2, i3, *}  
C ================================================================
      Integer IEP(*), nEP(*)

C group (Local variables)
      Integer ib(3), ie(3), ip(3)

      ip(1) = i1
      ip(2) = i2
      ip(3) = i3
      Do i = 1, 3
         If(ip(i).EQ.1) Then
            ib(i) = 1
         Else
            ib(i) = nEP(ip(i) - 1) + 1
         End if
         ie(i) = nEP(ip(i))
      End do

      Do 10 i = ib(1), ie(1)
         iE2 = IEP(i)
         If(iE2.EQ.iE1) Goto 10
         Do j = ib(2), ie(2)
            If(iE2.EQ.IEP(j)) Then
               Do k = ib(3), ie(3)
                  If(iE2.EQ.IEP(k)) Then
                     cmpE = .TRUE.
                     Goto 1000
                  End if
               End do
               Goto 10
            End if
         End do
 10   Continue

      cmpE = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function cmpF(i1, i2, IFP, nFP, iF1, iF2)
C ================================================================
      Integer IFP(*), nFP(*)

C group (Local variables)
      Integer ib(2), ie(2), ip(2)

      ip(1) = i1
      ip(2) = i2
      Do i = 1, 2
         If(ip(i).EQ.1) Then
            ib(i) = 1
         Else
            ib(i) = nFP(ip(i) - 1) + 1
         End if
         ie(i) = nFP(ip(i))
      End do

      Do 10 i = ib(1), ie(1)
         iF2 = IFP(i)
         If(iF2.EQ.iF1) Goto 10
         Do j = ib(2), ie(2)
            If(iF2.EQ.IFP(j)) Then
               cmpF = .TRUE.
               Goto 1000
            End if
         End do
 10   Continue

      cmpF = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function cmpP(iP, IFP, nFP, IPF)
C ================================================================
C cmpP = TRUE if point iP belongs to a common edge of two faces 
C with different color. Otherwise cmpP = FALSE.
C
C Remark: the routine doesn't say anything about cross points.
C ================================================================
      Integer IFP(*), nFP(*), IPF(4, *)

      If(iP.EQ.1) Then
         ib = 1
      Else
         ib = nFP(iP - 1) + 1
      End if
      ie = nFP(iP)

      Do i = ib, ie
         ICF1 = IPF(4, IFP(i))
         Do j = i + 1, ie
            ICF2 = IPF(4, IFP(j))
            If(ICF1.NE.ICF2) Then
               cmpP = .TRUE.
               Goto 1000
            End if
         End do
      End do

      cmpP = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function cmpR(iP, IFP, nFP, XYP, IPF)
C ================================================================
C cmpR = TRUE if point iP belongs to a common edge of 2 flat faces
C and the angle between these faces is smaller than 120 degrees. 
C ================================================================
      Real*8  XYP(3, *)
      Integer IFP(*), nFP(*), IPF(4, *)

      Integer iref(4)
      Real*8  c, angle2Faces
      Logical check22

C ================================================================
      DATA    iref/1, 2, 3, 1/
C ================================================================

      If(iP.EQ.1) Then
         ib = 1
      Else
         ib = nFP(iP - 1) + 1
      End if
      ie = nFP(iP)

      Do 20 i = ib, ie
         iF = IFP(i)

         Do 10 j = i + 1, ie
            jF = IFP(j)

            Do i1 = 1, 3
               i2 = iref(i1 + 1)

               iP1 = IPF(i1, iF)
               iP2 = IPF(i2, iF)

               Do j1 = 1, 3
                  j2 = iref(j1 + 1)

                  jP1 = IPF(j1, jF)
                  jP2 = IPF(j2, jF)


                  If(check22(iP1, iP2, jP1, jP2)) Then
                     i3 = iref(i2 + 1)
                     j3 = iref(j2 + 1)

                     iP3 = IPF(i3, iF)
                     jP3 = IPF(j3, jF)

                     c = angle2Faces(XYP(1, iP1), XYP(1, iP2), 
     &                               XYP(1, iP3), XYP(1, jP3))
c  ...               cos(120) = -0.5
                     If(c.GT.-0.5D0) Then
                        cmpR = .TRUE.
                        Goto 1000
                     End if
                  End if
               End do
            End do 
 10        Continue
 20   Continue

      cmpR = .FALSE.
 1000 Return
      End



C ================================================================
      Logical function crossPoint(XYP, nF, IFP, iP, IPF)
C ================================================================
C The number of edges with end point iP is evaluated.
C The array IFP is destroyed in out algorithm. 
C ================================================================
      Real*8  XYP(3, *)
      Integer IFP(*), IPF(4, *)

      Integer iRs(3)
      Real*8  angle2Edges, ang   
      Logical check13

C ================================================================
      crossPoint = .FALSE.
      
c ... analyze pairs of different faces
      lR = 0
      Do n = 1, nF
         iF1 = IFP(n)
         Do 20 m = n + 1, nF
            iF2 = IFP(m)
            If(IPF(4, iF1).EQ.IPF(4, iF2)) Goto 20

            Do 10 i = 1, 3
               iPt = IPF(i, iF1)
               If(iPt.EQ.iP) Goto 10

               If(check13(iPt, IPF(1, iF2), 
     &                         IPF(2, iF2), IPF(3, iF2))) Then
                  Call findSE(lR, iRs, iPt, nPt)
                  If(nPt.EQ.0) Then
                     lR = lR + 1
                     If(lR.GE.3) Then
                        crossPoint = .TRUE.
                        Goto 1000
                     End if

                     iRs(lR) = iPt
                  End if
               End if
 10         Continue
 20      Continue
      End do

      If(lR.LE.1) Then
         crossPoint = .FALSE.
      Else
         ang = angle2Edges(XYP(1, iP), XYP(1, iRs(1)), XYP(1, iRs(2)))
         If(ang.GT.0D0) crossPoint = .TRUE.
      End if

 1000 Return
      End



C ================================================================
      Integer Function countColors(N, ICE)
C ================================================================
C Compute the amount of different colors in array ICE(N).
C The array is destroyed in our algorithm. The different colors
C are gathered in the first entrices of the array.
C
C Remark : only positive colors are counted
C ================================================================
      Integer ICE(*)

      countColors = 0
      If(N.EQ.0) Return

      M = N
  1   countColors = countColors + 1
      k = countColors
      ic = ICE(k)
      Do i = countColors + 1, M
         If(ICE(i).NE.ic) Then
            k = k + 1
            ICE(k) = ICE(i)
         End if
      End do

      M = k
      If(countColors.LT.M) Goto 1

c ... delete negative colors
      countColors = 0
      Do i = 1, M
         If(ICE(i).GT.0) Then
            countColors = countColors + 1
            ICE(countColors) = ICE(i)
         End if
      End do
      Return
      End



C ================================================================
      Subroutine scale2Cube(nP, XYP, flag)
C ================================================================
C Routine scales the model to the square [0.1, 0.9]^2. We allow
C 10% freedom for curved edges.
C ================================================================
      Real*8  XYP(3, *)
      Logical flag

C ================================================================
      Real*8  minXYP(3), maxXYP(3), scale, size

      Real*8  refXYP(3), scaXYP(3)
      Common /aniCRV/refXYP, scaXYP

C ================================================================
      If(flag) Then
         Do i = 1, 3
            minXYP(i) = XYP(i, 1)
            maxXYP(i) = XYP(i, 1)
         End do

         Do n = 2, nP
            Do i = 1, 3
               minXYP(i) = min(minXYP(i), XYP(i, n))
               maxXYP(i) = max(maxXYP(i), XYP(i, n))
            End do
         End do

c  ...  add 5% for the bounding box
c        Do i = 1, 3
c           size = (maxXYP(i) - minXYP(i)) / 20 
c           minXYP(i) = minXYP(i) - size
c           maxXYP(i) = maxXYP(i) + size
c        End do


         Do i = 1, 3
            refXYP(i) = minXYP(i)
            scaXYP(i) = 1D0 / (maxXYP(i) - minXYP(i))
         End do

         scale = min(scaXYP(1), scaXYP(2))
         scale = min(scale,     scaXYP(3))

         Do i = 1, 3
            scaXYP(i) = scale
         End do

         Do n = 1, nP
            Do i = 1, 3
               XYP(i, n) = (XYP(i, n) - refXYP(i)) * scaXYP(i)
            End do
         End do
      Else
         Do i = 1, 3
            scaXYP(i) = 1D0 / scaXYP(i)
         End do

         Do n = 1, nP
            Do i = 1, 3
               XYP(i, n) = refXYP(i) + XYP(i, n) * scaXYP(i)
            End do
         End do
      End if
      Return
      End



C ================================================================
      Subroutine scaleBack(XYPi, XYPo)
C ================================================================
C  Routine computes physical coordinates of point XYPi 
C ================================================================
      Real*8   XYPi(3), XYPo(3)

      Real*8   refXYP(3), scaXYP(3)
      Common  /aniCRV/refXYP, scaXYP

C ================================================================
      Do i = 1, 3
         XYPo(i) = refXYP(i) + XYPi(i) / scaXYP(i)
      End do

      Return
      End



C ==============================================================
      Subroutine RandR(XY1, XY2, XY3, XY4, rOut, rIn)
C ==============================================================
C Computes curcumscribed and inscribed radii for the tetrahedron
C given by forth vertices. 
C ==============================================================
      Real*8  XY1(3), XY2(3), XY3(3), XY4(3)
      Real*8  rOut, rIn, calVol, calSqr, calEdge

      Real*8  v(3, 4), C(3), F(3), vol, sqr
C ==============================================================
      Do i = 1, 3
         v(1, i) = XY1(i) - XY4(i)
         v(2, i) = XY2(i) - XY4(i)
         v(3, i) = XY3(i) - XY4(i)
         v(i, 4) = 0D0
      End do

      Do i = 1, 3
         F(i) = 0D0
         Do j = 1, 3
            F(i) = F(i) + v(i, j) * (XY1(j) + XY4(j)) 
         End do
         F(i) = F(i) / 2
      End do

      vol = calVol(v(1, 1), v(1, 2), v(1, 3), v(1, 4))

      C(1) = calVol(F, v(1, 2), v(1, 3), v(1, 4)) / vol
      C(2) = calVol(v(1, 1), F, v(1, 3), v(1, 4)) / vol
      C(3) = calVol(v(1, 1), v(1, 2), F, v(1, 4)) / vol
     
      rOut = calEdge(XY1, C)

      sqr = calSqr(XY1, XY2, XY3) + calSqr(XY2, XY3, XY4)
     &    + calSqr(XY3, XY4, XY1) + calSqr(XY4, XY1, XY2)

      rIn = 6 * dabs(vol) / sqr
      Return
      End



C ==============================================================
      Subroutine copyMeshData(nP, nE, XYP,  HesP,  IPE, 
     &                                XYPw, HesPw, IPEw)
C ==============================================================
      Integer IPE(4, *), IPEw(4, *)
      Real*8  XYP(3, *), XYPw(3, *), HesP(6, *), HesPw(6, *)
C ==============================================================
      Do n = 1, nP
         Do i = 1, 3
             XYPw(i, n) = XYP(i, n)
          End do
       End do

       Do n = 1, nP
         Do i = 1, 6
            HesPw(i, n) = HesP(i, n)
         End do
      End do

      Do n = 1, nE
         Do i = 1, 4
            IPEw(i, n) = IPE(i, n)
         End do
      End do

      Return
      End
 

C ==============================================================
      Real*8  Function surfaceArea(nF, XYP, IPF, ic)
C ==============================================================
C The routine computes area of surface maked as ic. If ic <= 0,
C area of the total surface is computed.
C ==============================================================
      Real*8  XYP(3, *)
      Integer IPF(4, *)

c (Local variables)
      Real*8  calSqr, s

      s = 0D0
      Do 10 n = 1, nF
         If(IPF(1, n).LE.0) Goto 10

         If(IPF(4, n).EQ.ic .OR. ic.LE.0) Then
            iP1 = IPF(1, n)
            iP2 = IPF(2, n)
            iP3 = IPF(3, n)

            s = s + calSqr(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
         End if
 10   Continue

      surfaceArea = s
      Return
      End


C ==============================================================
      Real*8  Function fixedArea(nFv, XYP, IPF, IFV)
C ==============================================================
C The routine computes area of surface maked as ic. 
C ==============================================================
      Real*8  XYP(3, *)
      Integer IPF(4, *), IFV(*)

c (Local variables)
      Real*8  calSqr, s

      s = 0D0
      Do 10 n = 1, nFv
         iF = IFV(n)

         iP1 = IPF(1, iF)
         iP2 = IPF(2, iF)
         iP3 = IPF(3, iF)

         s = s + calSqr(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
 10   Continue

      fixedArea = s
      Return
      End


C ==============================================================
      Real*8  Function domainVolume(nE, XYP, IPE, ic)
C ==============================================================
C The routine computes volume of subdomain maked as ic. 
C If ic = 0, volume the whole domain is computed.
C ==============================================================
      Real*8  XYP(3, *)
      Integer IPE(5, *)

c (Local variables)
      Real*8  calVol, s

      s = 0D0
      Do 10 n = 1, nE
         If(IPE(1, n).LE.0) Goto 10

         If(IPE(5, n).EQ.ic .OR. ic.LE.0) Then
            iP1 = IPE(1, n)
            iP2 = IPE(2, n)
            iP3 = IPE(3, n)
            iP4 = IPE(4, n)

            s = s + dabs(calVol(XYP(1, iP1), XYP(1, iP2), 
     &                          XYP(1, iP3), XYP(1, iP4)))
         End if
 10   Continue

      domainVolume = s 
      Return
      End


C ==============================================================
      Real*8  Function fixedVolume(nEv, XYP, IPE, IEV)
C ==============================================================
C The routine computes volume of a fixed domain.
C ==============================================================
      Real*8  XYP(3, *)
      Integer IPE(5, *), IEV(*)

c (Local variables)
      Real*8  calVol, s
                                                            
      s = 0D0
      Do 10 n = 1, nEv
         iE = IEV(n)
                                                                     
         iP1 = IPE(1, iE)
         iP2 = IPE(2, iE)
         iP3 = IPE(3, iE)
         iP4 = IPE(4, iE)
                                                     
         s = s + dabs(calVol(XYP(1, iP1), XYP(1, iP2),
     &                       XYP(1, iP3), XYP(1, iP4)))
 10   Continue
                                                                   
      fixedVolume = s
      Return
      End
 
