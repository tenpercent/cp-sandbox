C ================================================================
      Subroutine check_mesh(
C ================================================================
     &            nP, nF, nE,
     &            XYP, IPF, IPE, lbF, lbE,
     &            nEP, IEP, IFE, IEE)
C ================================================================
C Routine checks topology of the input and output meshes.
C
C nEP :  working array of size nP
C IEP :  working array of size 4*nE
C IFE :  working array of size 4*nE
c IEE :  working array of size 4*nE
C
C ================================================================
C group (M)
      Integer IPF(3, *), IPE(4, *), lbF(*), lbE(*)
      Real*8  XYP(3, *)

C group (W)
      Integer IFE(4, *), IEE(4, *)
      Integer nEP(*), IEP(*)

C ================================================================
C group (Local variables)
      Integer ip(5)
      Real*8  calVol, calEdge, v1, v2
      Logical cmpE, check33

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


C ... check faces
      Do n = 1, nF
         nFac = n

         Do i = 1, 3
            If(IPF(i, n).LE.0) Goto 400
         End do

         If(lbF(n).LE.0) Goto 400
      End do


C ... create an auxiliary structure
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)


C ... create IEE
      Do n = 1, nE
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

            If(cmpE(ip1, ip2, ip3, IEP, nEP, n, iE2)) Then
               IEE(i1, n) = iE2
           End if
         End do
      End do


c ... create an auxiliary structure
      Call backReferences(nP, nF, 3, 3, IPF, nEP, IEP)


c ... create IFE: basic, fictitious, material
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            If(cmpE(ip1, ip2, ip3, IEP, nEP, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do


C ... check for faces
      Do n = 1, nP
         nEP(n) = 0
      End do

      Do n = 1, nE
         nTet = n

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)
         If(iP1.EQ.iP2 .OR. iP1.EQ.iP3 .OR. iP1.EQ.iP4 .OR.
     &      iP2.EQ.iP3 .OR. iP2.EQ.iP4 .OR. iP3.EQ.iP4) Then
            iERR = 5013
            Goto 500
         End if


         If(lbE(n).LE.0)
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

                     v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                           XYP(1, iP3), XYP(1, iP4))
                     v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                           XYP(1, iP3), XYP(1, jP4))

                     If(v1 * v2.GE.0D0) Then
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

      Return


 400  Continue
      Write(*, 5002) nFac, (IPF(i, nFac), i = 1, 3), lbF(nFac)
      Return


 500  Write(*, 5000) nTet, iERR, (IPE(i, nTet), i = 1, 4),
     &                           (IFE(i, nTet), i = 1, 4),
     &                           (IEE(i, nTet), i = 1, 4)
      Do k = 1, 4
         iE = IEE(k, nTet)
         If(iE.GT.0) Then
            Write(*, 5000) iE, iERR, (IPE(i, iE), i = 1, 4),
     &                               (IFE(i, iE), i = 1, 4),
     &                               (IEE(i, iE), i = 1, 4)
         End if
      End do
      Call errMes(iERR, 'chkM', 'tetrahedra are wrong')

 5000 Format('Error in checking tetrahedron =', I7, '   iERR=', I5, /,
     &       'Points =', 4I7, /,
     &       'Faces  =', 4I7, /,
     &       'Tetras =', 4I7, /)

 5002 Format('Error in face =', I6, ',  Points =', 3I6, ',  label=',I4)

 5004 Format('Face ', I7, '  (', I6, ')  of the bad tetrahedron', /,
     &       'Points =', 3I7, /)

      End



 
