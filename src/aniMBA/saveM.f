C ================================================================
C Saves mesh in ANI format
C ================================================================
      Subroutine saveMani(
C ================================================================
c group (M)
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (I)
     &      flagI, IPP, IFF,
     &      fName)
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

c group (Dev)
      Integer IPV(*), IFV(*), IEV(*)

c group (I)
      Logical flagI
      Integer IPP(*), IFF(*)

      Character*(*) fName

C ================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do

      Write(*,'(A,A)') 'Saving mesh ', fName(1:i+3)

      kps = 11
      kfs = kps + nP + 2
      kes = kfs + nF + 2
      kcs = kes + nE + 2

      kfp = kcs + 2
      kff = kfp + nPv + 2
      kfe = kff + nFv + 2

      kis = kfe + nEv + 2

      Open(10, file=fName(1:i+3), status='UNKNOWN')
      Write(10, '(3(A,I8),A)') 'T points:        ', 
     &          nP, ' (lines ', kps, ' - ', kps + nP - 1, ')'
      Write(10, '(3(A,I8),A)') 'T faces:         ',
     &          nF, ' (lines ', kfs, ' - ', kfs + nF - 1, ')'
      Write(10, '(3(A,I8),A)') 'T elements:      ',
     &          nE, ' (lines ', kes, ' - ', kes + nE - 1, ')'

      Write(10, '(A)') 'T curved faces:     0'

      If(nPv.NE.0) Then
         Write(10, '(3(A,I8),A)') 'T fixed points:  ',
     &             nPv, ' (lines ', kfp, ' - ', kfp + nPv - 1, ')'
      Else
         Write(10, '(A,I8)') 'T fixed points:  ', nPv
      End if

      If(nFv.NE.0) Then
         Write(10, '(3(A,I8),A)') 'T fixed faces:   ', 
     &             nFv, ' (lines ', kff, ' - ', kff + nFv - 1, ')'
      Else
         Write(10, '(A,I8)') 'T fixed faces:   ', nFv
      End if

      If(nEv.NE.0) Then
         Write(10, '(3(A,I8),A)') 'T fixed elements:', 
     &             nEv, ' (lines ', kfe, ' - ', kfe + nEv - 1, ')'
      Else
         Write(10, '(A,I8)') 'T fixed elements:', nEv 
      End if

      Write(10, '(L1,A)') flagI, ' interface data:'

      Write(10,*)
      Write(10,*) nP, ' # of nodes'
      Do n = 1, nP
         Write(10,'(3E23.15)') (XYP(j, n), j = 1, 3)
      End do

      Write(10,*)
      Write(10,*) nF, ' # of faces'
      Do n = 1, nF
         Write(10,*) (IPF(j, n), j = 1, 3), lbF(n)
      End do


      Write(10,*)
      Write(10,*) nE, ' # of elements'
      Do n = 1, nE
         Write(10,*) (IPE(j, n), j = 1, 4), lbE(n)
      End do


c ... saving curvilinear faces
      Write(10,*)
      Write(10,*) '0 # of curvilinear faces'


c ... saving the fixed mesh points, faces and elements
      Write(10,*)
      Write(10,*) nPv, ' # number of fixed points'
      Do n = 1, nPv
         Write(10,*) IPV(n)
      End do


      Write(10,*)
      Write(10,*) nFv, ' # number of fixed faces'
      Do n = 1, nFv
         Write(10,*) IFV(n)
      End do


      Write(10,*)
      Write(10,*) nEv, ' # number of fixed elements'
      Do n = 1, nEv
         Write(10,*) IEV(n)
      End do

      If(flagI) Then
         Write(10,*)
         Write(10,*) nP, nF, ' interface data'
         Do n = 1, nP
            Write(10,*) IPP(n)
         End do
         Do n = 1, nF
            Write(10,*) IFF(n)
         End do
      End if
      Close(10)

 5000 Format(2I6,6F15.11)
      Return
      End



C ================================================================
      Subroutine saveS(nP, Sol, fName)
C ================================================================
      Real*8  Sol(*)
      Character*(*) fName

C ================================================================
      Write(*,'(A,A)') 'Saving solution ', fName

      Open(10, file=fName, status='UNKNOWN')
      Write(10,*) nP
      Write(10,*)

      Do n = 1, nP
         Write(10,*) Sol(n)
      End do
      Close(10)

      Return
      End



C =====================================================================
      Subroutine saveMgmv(
C =====================================================================
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      fName, iW)
C =====================================================================
C Routine saves mesh in the GMV file. (New version)
C
C *** Remarks:
C        1. The size of the working memory is nE
C =====================================================================
C group (M)
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*), iW(*)
      Character*(*) fName

C group (Local variables)
      Integer countColors
      Real*8  calVol, v

C =====================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.gmv')
         i = i + 1
      End do

      Write(*,'(A,A)') 'Saving GMV image ', fName(1:i+3)

      Open(10, file=fName(1:i+3), status='UNKNOWN')

      Write(10, '(A)') 'gmvinput ascii'
      Write(10, *)
      Write(10, *) 'nodev ', nP
      Do n = 1, nP
         Write(10, *) (XYP(i, n), i = 1, 3)
      End do

c ... save cells
      Write(10, *)
      Write(10 ,*) 'cells ', nE

      Do n = 1, nE
         iP1 = IPE(1, n) 
         iP2 = IPE(2, n) 
         iP3 = IPE(3, n) 
         iP4 = IPE(4, n) 

         v = calVol(XYP(1, iP1), XYP(1, iP2), 
     &              XYP(1, iP3), XYP(1, iP4))
         If(v.GT.0D0) Then
            Write(10, *) ' tet 4 ', iP1, iP2, iP3, iP4
         Else
            Write(10, *) ' tet 4 ', iP1, iP3, iP2, iP4
         End if
      End do
      Write(10, *)
      Write(10, *)

c ... save materials
      Do n = 1, nE
         iW(n) = lbE(n)
      End do
      ic = countColors(nE, iW)

      Write(10, *) 'material ', ic, ' 0'
      Do i = 1, ic
         If(iW(i).LT.10) Then
            Write(10,'(A,I1)') 'mat', iW(i)
         Else If(iW(i).LT.100) Then
            Write(10,'(A,I2)') 'mat', iW(i)
         Else If(iW(i).LT.1000) Then
            Write(10,'(A,I3)') 'mat', iW(i)
         Else If(iW(i).LT.10000) Then
            Write(10,'(A,I4)') 'mat', iW(i)
         Else 
            Call errMes(6004, 'saveMgmv', 'missing code')
         End if
      End do

      Write(10, *)  (lbE(i), i = 1, nE)

c ... save faces
      Write(10, '(A)') 'polygons'
      Do n = 1, nF
         Write(10, *) lbF(n), ' 3 '
         Do k = 1, 3
            Write(10, *) (XYP(k, IPF(i, n)), i = 1, 3)
         End do
      End do

      Write(10, '(A)') 'endpoly'
      Write(10, *)
      Write(10, '(A)') 'endgmv'

      Close(10)

      Return
      End



C ================================================================
      Subroutine saveMgeo(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, 
     &      XYP, IPF, IPE, lbP, lbF, lbE,
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
     &      fName)
C ================================================================
      include 'color.fd'
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbP(*), lbF(*), lbE(*)

c group (Dev)
      Integer IPV(*), IFV(*), IEV(*)

      Character*(*) fName

C group (Local variables)
      Logical flag, ifXnode
      Character*70 fNameExt

C ================================================================
c ... counting the boundary faces (array lbF is overloaded)
      nB = 0
      Do n = 1, nF
         flag = .TRUE.

         Do i = 1, 3
            If(.NOT.ifXnode(lbP(IPF(i, n)), jBnode)) flag = .FALSE.
         End do

         If(flag) Then
            nB = nB + 1
            lbF(n) = -lbF(n)
         End if
      End do


      fNameExt = fName // '.geo'
      Write(*,'(A,A)') 'Saving mesh ', fNameExt


      Open(10, file=fNameExt, status='UNKNOWN')
      Write(10,*) nP, nE, nB

      Write(10,*)
      Do n = 1, nP
         Write(10,*) n, (XYP(j, n), j = 1, 3)
      End do

      Write(10,*)
      Do n = 1, nE
         Write(10,*) n, (IPE(j, n), j = 1, 4), lbE(n)
      End do

      Write(10,*)
      k = 0
      Do n = 1, nF
         If(lbF(n).LT.0) Then
            lbF(n) = -lbF(n)
            k = k + 1
            Write(10,*) k, (IPF(j, n), j = 1, 3), lbF(n)
         End if
      End do
      Close(10)

      Return
      End



C ================================================================
C Saves mesh in AFT format
C ================================================================
      Subroutine saveMaft(
C ================================================================
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      fName)
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

      Character*(*) fName

C ================================================================
      i = 1
      Do while(fName(i:i+3) .NE. '.out')
         i = i + 1
      End do

      Write(*,'(A,A)') 'Saving aft-mesh ', fName(1:i+3)

      Open(10, file=fName(1:i+3), status='UNKNOWN', ERR=1000)

c ... write points
      Write(10, *) nP
      Do n = 1, nP
         Write(10, *) (XYP(i, n), i = 1, 3)
      End do 

c ... write elements 
      Write(10, *) nE
      Do n = 1, nE
         Write(10, *) (IPE(i, n), i = 1, 4), lbE(n)
      End do

c ... write faces    
      Write(10, *) nF
      Do n = 1, nF
         Write(10, *) (IPF(i, n), i = 1, 3), lbF(n)
      End do
      
      Close(10)
      Return
 1000 Continue
      Call errMesIO(4001, 'saveMaft', 'Input file name is wrong')

      Return
      End
