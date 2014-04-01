C ================================================================
      Subroutine chkM(
C ================================================================
c group (M)
     &            nP, nF, nE,
     &            XYP, IPF, IPE,
     &            ICP, IFE, IEE,
     &            rR, status,
c group (W)
     &            iCheck, nEPw, IEPw)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'status.fd'
C ================================================================
C Routine checks topology of the input and output meshes.
C
C iCheck : 0 - check everything
C          1 - don't check for isolated faces and points marked as
C              destroyed structures
C
C nEPw  :  working array of size nP
C IEPw  :  working array of size 4*nE
C
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(5, *)
      Real*8  XYP(3, *)

      Integer ICP(*), IFE(4, *), IEE(4, *)

      Real*8  rR
      Integer status

C group (W)
      Integer nEPw(*), IEPw(*)

C ================================================================
C group (Local function)
      Logical ifXnode
      Integer countColors

C group (Local variables)
      Integer ip(5)
      Real*8  calVol, v1, vv, rOut, rIn, mutualOrientation
      Logical check33, flagFBE, flag

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      rR = 0D0
      crvPREC = 2D-8
      
      flagFBE = ifXnode(status, ANIForbidBoundaryElements)


C ... check for face markers
      Do n = 1, nF
         iCLRf = IPF(4, n)
         If(iCLRf.LE.0 .AND. iCheck.EQ.0) 
     &      Call errMes(5001, 'chkM', 'wrong face ID')
         If(iCLRf.GT.MaxS)
     &      Call errMes(5002, 'chkM', 'face ID is out of limits')
      End do


      Do n = 1, nP
         If(iCheck.EQ.0 .AND. ICP(n).EQ.0)
     &      Call errMes(5003, 'chkM', 'wrong point color')

         If(ifXnode(ICP(n), jBnode)) Then
            If(.NOT.ifXnode(ICP(n), jSnode)) 
     &         Call errMes(5103, 'chkM', 'wrong point color')

            If(ifXnode(ICP(n), jInode))
     &         Call errMes(5103, 'chkM', 'wrong point color')
         End if
      End do


      Do n = 1, nF
         nFac = n

         Do i = 1, 4
            If(iCheck.EQ.0 .AND. IPF(i, n).LE.0) Then
               iERR = 5012
               Goto 400
            End if
         End do

         If(iCheck.EQ.0 .AND. IPF(4, n).GE.iVface) 
     &      Call errMes(4103, 'chkM',
     &                 'reserved boundary identificator is used')
      End do


      Do n = 1, nP
         nEPw(n) = 0
      End do

      Do n = 1, nE
         nTet = n

         flag = flagFBE
         Do i = 1, 4
            iPt = IPE(i, n)
            If(iPt.LE.0)
     &         Call errMes(5006, 'chkM', 'wrong connectivity table')
            nEPw(iPt) = nEPw(iPt) + 1

            If(.NOT.ifXnode(ICP(iPt), jBnode)) flag = .FALSE.
         End do

         If(iCheck.EQ.0 .AND. flag) 
     &      Call errMes(5022, 'chkM', 'boundary element')


         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)
         If(iP1.EQ.iP2 .OR. iP1.EQ.iP3 .OR. iP1.EQ.iP4 .OR.
     &      iP2.EQ.iP3 .OR. iP2.EQ.iP4 .OR. iP3.EQ.iP4) Then
            iERR = 5013
            Goto 500
         End if


         If(IPE(5, n).LE.0)
     &      Call errMes(5023, 'chkM', 'non-positive element label')


         iF1 = IFE(1, n)
         iF2 = IFE(2, n)
         iF3 = IFE(3, n)
         iF4 = IFE(4, n)

         If(iF1.EQ.iF2 .AND. iF1.GT.0 .OR.
     &      iF1.EQ.iF3 .AND. iF1.GT.0 .OR.
     &      iF1.EQ.iF4 .AND. iF1.GT.0 .OR.
     &      iF2.EQ.iF3 .AND. iF2.GT.0 .OR.
     &      iF2.EQ.iF4 .AND. iF2.GT.0 .OR.
     &      iF3.EQ.iF4 .AND. iF3.GT.0) Then
            iERR = 5014
            Goto 500
         End if

         iE1 = IEE(1, n)
         iE2 = IEE(2, n)
         iE3 = IEE(3, n)
         iE4 = IEE(4, n)

         If(iE1.EQ.iE2 .AND. iE1.NE.0 .OR.
     &      iE1.EQ.iE3 .AND. iE1.NE.0 .OR.
     &      iE1.EQ.iE4 .AND. iE1.NE.0 .OR.
     &      iE2.EQ.iE3 .AND. iE2.NE.0 .OR.
     &      iE2.EQ.iE4 .AND. iE2.NE.0 .OR.
     &      iE3.EQ.iE4 .AND. iE3.NE.0) Then
            iERR = 5015
            Goto 500
         End if


         Do 20 i1 = 1, 4
            iF = IFE(i1, n)
            iE = IEE(i1, n)
            If(iF.EQ.0 .AND. iE.EQ.0) Then
               iERR = 5016
               Goto 500
            End if

            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)
            iP3 = IPE(i3, n)

            If(iF.NE.0 .AND. iF.NE.iFface) Then
               jP1 = IPF(1, iF)
               jP2 = IPF(2, iF)
               jP3 = IPF(3, iF)

               If(.NOT.check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  iERR = 5017
                  Goto 500
               End if
            End if

            If(iE.NE.0) Then
               If(iF.NE.0) Then
                  Do j1 = 1, 4
                     If(IFE(j1, iE).EQ.iF) Goto 10
                  End do

                  iERR = 5018
                  Goto 500
               End if

  10           Do j1 = 1, 4
                  j2 = ip(j1 + 1)
                  j3 = ip(j2 + 1)

                  jP1 = IPE(j1, iE)
                  jP2 = IPE(j2, iE)
                  jP3 = IPE(j3, iE)

                  If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                     If(IEE(j1, iE).NE.n) Then
                        iERR = 5019
                        Goto 500
                     End if

                     i4  = ip(i3 + 1)
                     iP4 = IPE(i4, n)

                     j4  = ip(j3 + 1)
                     jP4 = IPE(j4, iE)

c                    v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, iP4))
c                    v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, jP4))

                     vv = mutualOrientation(
     &                    XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                    XYP(1, iP4), XYP(1, jP4))

                     If(vv.GE.0D0 .AND. 
     &                  .NOT.ifXnode(status, ANITangledMesh)) Then
                        iERR = 5020
                        Goto 500
                     End if
                     Goto 20
                  End if
               End do

               iERR = 5021
               Goto 500
            End if
 20      Continue
      End do

      If (iCheck.eq.0) then
        Do n = 1, nP
           If(nEPw(n).EQ.0) Call errMes(5007, 'chkM', 'isolated point')
        End do
      End if


C ... check color of edge points
      Call backReferences(nP, nE, 4, 5, IPE, nEPw, IEPw)
      Do n = 1, 4 * nE
         IEPw(n) = IPE(5, IEPw(n))
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nEPw(n)

         ic = countColors(i2 - i1 + 1, IEPw(i1))
         If(ic.GE.3 .AND.
     &     .NOT.(ifXnode(ICP(n), jRnode) .OR. 
     &           ifXnode(ICP(n), jVnode))) Then
            Write(*,*) ic, ICP(n), n
            Call errMes(5008, 'chkM', 'wrong color of edge point')
         End if
      End do


      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call RandR(XYP(1, iP1), XYP(1, iP2),
     &              XYP(1, iP3), XYP(1, iP4), rOut, rIn)

         rR = max(rR, rOut / rIn)
      End do
      Return

 400  Write(*, 5002) nFac, iERR, (IPF(i, nFac), i = 1, 4),
     &                           (ICP(IPF(i, nFac)), i = 1, 3)

      Call errMes(iERR, 'chkM', 'faces are wrong')

 500  Continue
      Write(*, 4000) 
      Write(*, 5000) nTet, v1, iERR, (IPE(i, nTet), i = 1, 4),
     &                               (IFE(i, nTet), i = 1, 4),
     &                               (IEE(i, nTet), i = 1, 4),
     &                               (ICP(IPE(i, nTet)), i = 1, 4)
      Do k = 1, 4
         iE = IEE(k, nTet)
         If(iE.GT.0) Then
            iP1 = IPE(1, iE)
            iP2 = IPE(2, iE)
            iP3 = IPE(3, iE)
            iP4 = IPE(4, iE)

            v1 = calVol(XYP(1, iP1), XYP(1, iP2), 
     &                  XYP(1, iP3), XYP(1, iP4))

            Write(*, 5000) iE, v1, iERR, (IPE(i, iE), i = 1, 4),
     &                                   (IFE(i, iE), i = 1, 4),
     &                                   (IEE(i, iE), i = 1, 4),
     &                                   (ICP(IPE(i, iE)), i = 1, 4)
         End if
      End do
      Do k = 1, 4
         iF = IFE(k, nTet)
         If(iF.GT.0) Then
            Write(*, 5004) iF, nF, (IPF(i, iF), i = 1, 3),
     &                             (ICP(IPF(i, max(1, iF))), i = 1, 3)
         End if
      End do

c     Do n = 1, nE
c        IEPw(n) = 1
c     End do
c     Call saveMgmv(nP, nF, nE,  
c    &              XYP, IPF, IPE, IEPw, IEPw, 
c    &              'error.gmv', IEPw(nE+1))

      Call errMes(iERR, 'chkM', 'tetrahedra are wrong')

 4000 Format(/,'=============== ERROR details =================')

 5000 Format('Error in checking tet =', I7, '  vol =', E16.9, 
     &       '   iERR=', I4, /,
     &       'Points =', 4I7, /,
     &       'Faces  =', 4I7, /,
     &       'Tetras =', 4I7, /,
     &       'colors =', 4I7, /)

c5001 Format('Error =', E12.4, /,
c    &       'Curvilinear face=', I4, ' attributs=', 2I4)

 5002 Format('Error in checking face =', I6, '   iERR=', I5, /,
     &       'Points =', 4I6, /,
     &       'colors =', 3I6)

 5004 Format('Face ', I7, '  (', I6, ')  of the bad tetrahedron', /,
     &       'Points =', 3I7, /,
     &       'colors =', 3I7, /)

      End

