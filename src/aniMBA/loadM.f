C ================================================================
      Subroutine loadMani(
C ================================================================
c group (M)
     &      MaxP, MaxF, MaxE,   
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (I)
     &      IPP, IFF,
     &      fName)
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

c group (Dev)
      Integer IPV(*), IFV(*), IEV(*)

c group (I)
      Integer IPP(*), IFF(*)

      Character*(*) fName

c group (Local variables)
      Integer iCrv(2)
      Real*8  rCrv(6)

      Logical flagP, flagF, flagE, flagC, flagPV, flagFV, flagEV
      logical flagI

C ================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do

      Write(*,'(A,A)') 'Loading mesh ', fName(1:i+3)

      Open(10, file=fName(1:i+3), status='OLD', ERR=1000)

      Read(10,*) flagP 
      Read(10,*) flagF
      Read(10,*) flagE
      Read(10,*) flagC
      Read(10,*) flagPV
      Read(10,*) flagFV
      Read(10,*) flagEV
      Read(10,*) flagI

      Read(10,*)
      Read(10,*) nP
      If(nP.GT.MaxP) Call errMesIO(1003, 'loadMani',
     &                   'local parameter MaxP is small')
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j = 1, 3)
      End do


      Read(10,*)
      Read(10,*) nF
      If(nF.GT.MaxF) Call errMesIO(1004, 'loadMani',
     &                   'local parameter MaxF is small')
      Do n = 1, nF
         Read(10,*) (IPF(j, n), j = 1, 3), lbF(n)
      End do


      Read(10,*)
      Read(10,*) nE
      If(nE.GT.MaxE) Call errMesIO(1006, 'loadMani',
     &                   'local parameter MaxE is small')
      Do n = 1, nE
         Read(10,*) (IPE(j, n), j = 1, 4), lbE(n)
      End do


c ... reading curvilinear faces
      Read(10,*)
      Read(10,*) nC
      Do n = 1, nC
         Read(10, *) (iCrv(i), i = 1, 2), (rCrv(i), i = 1, 6)
      End do


c ... reading the fixed mesh points, faces and elements
      Read(10,*)
      Read(10,*) nPv
      Do n = 1, nPv
         Read(10,*) IPV(n)
      End do


      Read(10,*)
      Read(10,*) nFv
      Do n = 1, nFv
         Read(10,*) IFV(n)
      End do


      Read(10,*)
      Read(10,*) nEv
      Do n = 1, nEv
         Read(10,*) IEV(n)
      End do

c ... reading interface data
      If(flagI) Then
         Read(10,*)
         Read(10,*) nPo, nFo
         Do n = 1, nP
            IPP(n) = 0
         End do
         Do n = 1, nF
            IFF(n) = 0
         End do
         Do n = 1, min(nP, nPo)
            Read(10,*) IPP(n)
         End do
         Do n = 1, min(nF, nFo)
            Read(10,*) IFF(n)
         End do
      End if
      Close(10)

      Return

 1000 Continue
      Call errMesIO(4001, 'loadMani', 'Input file name is wrong ')

      Return
      End



C ================================================================
      Subroutine loadMani_header(nP, nF, nE, nPv, nFv, nEv, fName)
C ================================================================
C Routines reads the header from the input file
C ================================================================
      Integer       nP, nF, nE, nC, nPv, nFv, nEv
      Character*(*) fName

      Logical       flag
      Character*30  fNameExt, text, text2

C ================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'ani'

      Write(*,'(A,A)') 'Reading header of mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)

      Read(10,*) flag, text, nP
      Read(10,*) flag, text, nF
      Read(10,*) flag, text, nE
      Read(10,*) flag, text, text2, nC
      Read(10,*) flag, text, text2, nPv
      Read(10,*) flag, text, text2, nFv
      Read(10,*) flag, text, text2, nEv

      Close(10)
      Return

 1000 Continue
      Call errMes(4001, 'loadMani_header', 'File name is wrong')

      Return
      End



C ================================================================
      Subroutine loadS(nP, Sol, fName)
C ================================================================
      Real*8  Sol(*)
      Character*(*) fName

C group (Local variables)
      Character*70 fNameExt

C ================================================================
      Write(*,'(A,A)') 'Loading solution ', fName

c ... read the solution associated to the mesh
      fNameExt = fName

      Open(10, file=fNameExt, status='OLD', ERR = 1000)
      Read(10,*) nPw
      If(nPw.NE.nP) Call errMesIO(4001, 'loadS',
     &             'the lines number differs from the points number')

      Read(10,*)
      Do n = 1, nP
         Read(10,*) Sol(n)
      End do
      Close(10)

      Return

 1000 Continue
      Call errMesIO(4001, 'loadS', 'Input file name is wrong ')

      Return
      End



C ================================================================
      Subroutine loadMaft(
C ================================================================
c group (M)
     &      MaxP, MaxF, MaxE,   
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (I)
     &      IPP, IFF,
     &      fName)
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

c group (Dev)
      Integer IPV(*), IFV(*), IEV(*)

c group (I)
      Integer IPP(*), IFF(*)
      Character*(*) fName

C ================================================================
      i = 1
      Do while(fName(i:i+3) .NE. '.out')
         i = i + 1
      End do

      Write(*,'(A,A)') 'Loading mesh ', fName(1:i+3)

      Open(10, file=fName(1:i+3), status='OLD', ERR=1000)


c ... read points
      Read(10, *) nP
      If(nP.GT.MaxP) Call errMesIO(1005, 'loadMaft',
     &                   'local parameter MaxP is small')
      Do n = 1, nP
         Read(10, *) (XYP(i, n), i = 1, 3)
      End do 

c ... read elements 
      Read(10, *) nE
      If(nE.GT.MaxE) Call errMesIO(1006, 'loadMaft',
     &                   'local parameter MaxE is small')
      Do n = 1, nE
         Read(10, *) (IPE(i, n), i = 1, 4), lbE(n)
      End do

c ... read faces    
      Read(10, *) nF
      If(nF.GT.MaxF) Call errMesIO(1007, 'loadMaft',
     &                   'local parameter MaxF is small')
      Do n = 1, nF
         Read(10, *) (IPF(i, n), i = 1, 3), lbF(n)
      End do
      
      Close(10)


c ... set up the other parameters to default values
      nPv = 0
      nFv = 0
      nEv = 0

      Return

 1000 Continue
      Call errMesIO(4001, 'loadMaft', 'Input file name is wrong')

      Return
      End



C ================================================================
      Subroutine loadMgmv(
C ================================================================
c group (M)
     &      MaxP, MaxF, MaxE,   
     &      nP, nF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (I)
     &      IPP, IFF,
     &      fName)
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

c group (Dev)
      Integer IPV(*), IFV(*), IEV(*)

c group (I)
      Integer IPP(*), IFF(*)

      Character*(*) fName

c group (Local variables
      Character*30 text

C ================================================================
      Write(*,'(A,A)') 'Loading mesh ', fName

      Open(10, file=fName, status='OLD', ERR=1000)


c ... read header
      Do i = 1, 2
         Read(10, *)
      End do
  

c ... read points
      Read(10, *) text, nP
      Do n = 1, nP
         Read(10, *) (XYP(i, n), i = 1, 3)
      End do 


c ... read elements 
      Read(10, *) text, nE
      Do n = 1, nE
         Read(10, *) text, k, (IPE(i, n), i = 1, 4)
      End do

      Read(10, *) text, nMat
      Do n = 1, nMat
         Read(10, *) text
      End do
      
      Read(10, *) (lbE(n), n = 1, nE)
      Close(10)


c ... set up the other parameters to default values
      nPv = 0
      nF  = 0
      nFv = 0
      nC  = 0
      nEv = 0

      Return

 1000 Continue
      Call errMesIO(4001, 'loadMgmv', 'Input file name is wrong')

      Return
      End




