C ======================================================================
      Subroutine BilinearFormTemplate(
C ======================================================================
     &           nP, nFb, nE, XYP, lbP, IPF, lbF, IPE, lbE, 
     &           fem3Dext, DATA, assembleStatus, 
     &           MaxIA, MaxA, IA, JA, A, F, nRow, nCol,
     &           MaxWi, iW)
C ======================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'
C ======================================================================
C  The routine assembles a stiffness matrix for a bilinear form B.
C  Here we describe new way to assemble elemental matrices. See 
C  comments in file assemble.f for the way to create elemental
C  matrices.
C ======================================================================
C *** PROBLEM INFORMATION ***
C  The routines assemble elemental matrices for the bilinear form
C
C              B(u, v)
C
C  specified by the user. Each elemental matrix is a combination of 
C  fem3Dtet calls because fem3Dext may consists of a few simple bilinear 
C  forms when u and v are vectors of unknowns. 
C  The elemental matrix is described by array 'template' and its 
C  contruction may use array DATA if neccessary. 
C
C
C *** SUBROUTINE fem3Dext ***
C  
C     SUBROUTINE fem3Dext(XY1, XY2, XY3, XY4,
C    &                    lbEloc, lbFloc, lbRloc, lbPloc, DATA, iSYS,
C    &                    LDA, A, F, nRow, nCol,
C    &                    templateR, templateC)
C
C     Parameters:
C              XY1(3) - Real*8 Cartesian coordinates of the first
C                         vertex of the triangle.
C              XY2(3) - second vertex of the triangle
C              XY3(3) - third vertex of the triangle
C              XY4(3) - fourth vertex of the triangle
C
C              lbEloc    - label of the tetrahedron (material label)
C              lbFloc(4) - labels of its faces (boundary labels),
C                          lbFloc(i) = 0 for internal faces
C              lbRloc(6) - labels of its edges (copied from global lbR)
C              lbPloc(4) - labels of points (copied from global lbP)
C
C              DATA   - Real*8 user given data (a number or an array)
C
C              iSYS   - system buffer for information exchange:
C                     iSYS( 3)    <- tetrahedrom number
C                     iSYS( 4:7)  <- numbers of vertices
C                     iSYS( 8:13) <- numbers of edges
C                     iSYS(14:17) <- numbers of faces
C                     iSYS(18)    <- the number of points
C                     iSYS(19)    <- the number of edges
C                     iSYS(20)    <- the number of faces
C                     iSYS(21)    <- the number of tetrahedra
C
C              LDA       - leading dimention of matrix A
C              A(LDA, *) - elemental matrix, the order of degrees of freedom
C                          should be in an agreement with array template
C              F(nRow)   - the vector of the elemental right-hand side     
C 
C              nRow   - the number of rows in A
C              nCol   - the number of columns in A
C
C              templateR(nRow) - array of degrees of freedom for rows;
C                                the D.O.F. must follows in groups, e.g. 
C                                four for vertices, six for edges, etc.
C
C              templateC(nCol) - array of degrees of freedom for columns. 
C
C
C *** SPARSE MATRIX ***
C
C   By default, our sparse matrix is assembled in one of the row formats
C   (IA, JA, A) specified by assembleStatus. Other formats are supported
C   thorough the library of format transformations.
C
C     assembleStatus - some a priory information about the matrix A:
C
C              MATRIX_SYMMETRIC - symmetric matrix
C              MATRIX_GENERAL   - general matrix
C
C              FORMAT_AMG       - format used in AMG
C              FORMAT_CSR       - compressed row format (diagonal is not saved)
C
C                REMARK: Any logical AND can be used to defined
C                        this variable, possible contradictions will
C                        be checked by the code.
C
C     MaxA   - the maximal size of array A - the maximal
C              number of nonzero entries
C
C     IA, JA, A - row sparcity structure of the matrix A:
C
C                IA(k + 1) - IA(k) equals to the number of
C                            nonzero entries in the k-th row
C
C                JA(M) - list of column indexes of non-zero
C                        entries ordered by rows; M = IA(nRow + 1) - 1
C
C                A(M)  - non-zero entries ordered as in JA
C
C     nRow      - the number of rows in A
C     nCol      - the number of columns in A, nCol = nRow for
C                 symmetric bilinear forms.
C
C
C  The following rules are applied for numbering of unknowns:
C      A) 1st, unknowns associated with vertices are numerated. 
C
C      B) 2nd, unknowns associated with faces are numerated.
C
C      C) 3rd, unknowns associated with edges are numerated.
C
C      D) 4th, unknowns associated with elements are numerated.
C
C      E) The unknowns corresponding to vector functions (i.e. 
C         with 3 degrees of freedom per a mesh object (vertex, face,
C         edge, or element) are enumerated first by the corresponding 
C         mesh objects and then by the space coordinates: x, y, and z.
C
C
C *** WORKING MEMORY ***
C     MaxWi  - the size of the working integer array
C
C     iW(MaxWi)  - the integer working array. 
C
C ======================================================================
C *** Note:
C       Input parameters:  nP, nE, XYP, IPE, lbE,
C                          fem3Dext, DATA, assembleStatus,
C                          MaxIA, MaxA, MaxWi
C
C       Input / output:    iSYS, iW
C
C       Output parameters: IA, JA, A, F, nRow, nCol
C
C ======================================================================
C
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Yu. Vassilevski (yuri.vassilevski@gmail.ru)
C *** Date:   2005 - 2007
C *** Complains & Comments: lipnikov@gmail.com
C
C ================================================================
      Real*8   XYP(3, *)
      Integer  lbP(*), IPF(3, *), lbF(*), IPE(4, *), lbE(*)

      Integer  iSYS(MAXiSYS)
      Real*8   DATA(*)
      EXTERNAL fem3Dext

      Integer  IA(*), JA(*), assembleStatus
      Real*8   A(*), F(*)

      Integer  iW(*)

C ======================================================================
C Local variables
      Integer templateR(MaxSize), templateC(MaxSize)
      Real*8  Aloc(MaxSize, MaxSize), Floc(MaxSize)
      Integer indR(MaxSize), indC(MaxSize)

      Integer lbFloc(4), lbRloc(6), lbPloc(4), iIEE(1), p
      Logical flagMS, flagMG, flagAMG, flagCSR

      Logical flagP2P, flagP2R, flagP2F, flagP2E 
      Logical flagR2P, flagR2R, flagR2F, flagR2E 
      Logical flagF2P, flagF2R, flagF2F, flagF2E 
      Logical flag,    flagE2R, flagE2F

      Integer findY, non_zero, iusage, rusage

      Character*80 message

C ======================================================================
      nF = 0
      MaxF = 4 * nE

      nR = 0
      MaxR = 6 * nE

c ... default values for flags
      flagMG = .TRUE.

      flagMS = IAND(assembleStatus, MATRIX_SYMMETRIC).NE.0
      flagMG = IAND(assembleStatus, MATRIX_GENERAL).NE.0
      If(flagMG) Then
         flagMS = .FALSE.
         Write(*,'(A)') 'FEM: matrix with full structure'
      Else
         Write(*,'(A)') 'FEM: upper tiangular part of matrix'
      End if

      flagAMG = IAND(assembleStatus, FORMAT_AMG).NE.0
      flagCSR = IAND(assembleStatus, FORMAT_CSR).NE.0
      If(flagCSR) Then
         Write(*,'(A)') '     sparse compressed row format'
      Else 
         Write(*,'(A)') '     sparse compressed AMG row format'
      End if


      If(flagMS .AND. flagMG  .OR.  flagAMG .AND. flagCSR)
     &   Call errMesFEM(2001, 'BilinearFormTemplate', 
     &                     'assembleStatus is wrong')


c ... order indices of vertices in the connectivity table E->P
      Call order2D(4, nE, IPE)


c ... get template array
      iP1 = IPE(1, 1)
      iP2 = IPE(2, 1)
      iP3 = IPE(3, 1)
      iP4 = IPE(4, 1)

      Do i = 1, 4
         lbFloc(i) = 0
         lbPloc(i) = 0
      End do

      Do i = 1, 6
         lbRloc(i) = 0
      End do

      Call encodeISYSfull(1, iP1,iP2,iP3,iP4, 0,0,0,0,0,0, 
     &                    0,0,0,0, nP, nR, nF, nE, iSYS)

      Call fem3Dext(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, iP4),
     &              lbE(1), lbFloc, lbRloc, lbPloc, DATA, iSYS,
     &              MaxSize, Aloc, Floc, ir, ic, 
     &              templateR, templateC)
 
      If(ir.GT.MaxSize .OR. ic.GT.MaxSize) 
     &  Call errMesFEM(2001, 'BilinearFormTemplate', 'MaxSize is small')


c ... analyze rows
      irp = 0
      irf = 0
      irr = 0
      ire = 0
      Do i = 1, ir
         If(templateR(i).EQ.Vdof) Then
            irp = irp + 1
         Else If(templateR(i).EQ.Fdof .OR.
     &           templateR(i).EQ.FdofOrient) Then
            irf = irf + 1
         Else If(templateR(i).EQ.Rdof .OR.
     &           templateR(i).EQ.RdofOrient) Then
            irr = irr + 1
         Else If(templateR(i).EQ.Edof) Then
            ire = ire + 1
         End if
      End do
      irp = irp / 4 
      irf = irf / 4 
      irr = irr / 6 


c ... analyze columns
      icp = 0
      icf = 0
      icr = 0
      ice = 0
      Do i = 1, ic
         If(templateC(i).EQ.Vdof) Then
            icp = icp + 1
         Else If(templateC(i).EQ.Fdof .OR.
     &           templateC(i).EQ.FdofOrient) Then
            icf = icf + 1
         Else If(templateC(i).EQ.Rdof .OR.
     &           templateC(i).EQ.RdofOrient) Then
            icr = icr + 1
         Else If(templateC(i).EQ.Edof) Then
            ice = ice + 1
         End if
      End do
      icp = icp / 4
      icf = icf / 4
      icr = icr / 6


c ... default pointers to lists X -> Y where X,Y = P,R,F,E
      Call initPointers(irp, irp, icp, icp, 
     &                  inPP, iIPP, flagP2P, i, j, flag)
      Call initPointers(irp, irr, icp, icr, 
     &                  inRP, iIRP, flagP2R, inPR, iIPR, flagR2P)
      Call initPointers(irp, irf, icp, icf, 
     &                  inFP, iIFP, flagP2F, inPF, iIPF, flagF2P)
      Call initPointers(irp, ire, icp, ice, 
     $                  inEP, iIEP, flagP2E, i, j, flag)

      Call initPointers(irr, irr, icr, icr, 
     $                  inRR, iIRR, flagR2R, i, j, flag)
      Call initPointers(irr, irf, icr, icf, 
     $                  inFR, iIFR, flagR2F, inRF, iIRF, flagF2R)
      Call initPointers(irr, ire, icr, ice, 
     $                  inER, iIER, flagR2E, inRE, iIRE, flagE2R)

      Call initPointers(irf, irf, icf, icf, 
     $                  inFF, iIFF, flagF2F, i, j, flag)
      Call initPointers(irf, ire, icf, ice, 
     $                  inEF, iIEF, flagF2E, inFE, iIFE, flagE2F)


c ... resolve dependencies
      If(irp.GT.0) flagP2E = .TRUE.
      If(irr.GT.0) Then
         flagR2E = .TRUE.
         flagE2R = .TRUE.
      End if
      If(irf.GT.0) Then
         flagF2E = .TRUE.
         flagE2F = .TRUE.
      End if

      If(icr.GT.0) flagE2R = .TRUE.
      If(icf.GT.0) flagE2F = .TRUE.

      flagE2F = .TRUE.

      iEnd = 2

c ... populate lists E -> R
      If(flagE2R) Then
         Call memory(iIRE, iEnd, 6*nE)
         Call memory(inEP, iEnd, nP)
         Call memory(iIEP, iEnd, 4*nE)
         If(iEnd.GT.MaxWi) Goto 5000

         Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))
      End if


c ... populate lists E -> F
      If(flagE2F) Then
         Call memory(iIFE, iEnd, 4*nE)
         Call memory(inEP, iEnd, nP)
         Call memory(iIEP, iEnd, 4*nE)
         If(iEnd.GT.MaxWi) Goto 5000

         Call listE2F(nP, nF, nE, IPE, iW(iIFE), iW(inEP), iW(iIEP))
      End if


c ... populate lists P -> E
      If(flagP2E) Then
         Call memory(inEP, iEnd, nP)
         Call memory(iIEP, iEnd, 4*nE)
         If(iEnd.GT.MaxWi) Goto 5000

         Call backReferences(nP, nE, 4,4, IPE, iW(inEP), iW(iIEP))
      End if


c ... populate lists F -> E
      If(flagF2E) Then
         Call memory(inEF, iEnd, nF)
         Call memory(iIEF, iEnd, 4*nE)
         If(iEnd.GT.MaxWi) Goto 5000

         Call backReferences(nF, nE, 4,4, iW(iIFE), iW(inEF), iW(iIEF))
      End if


c ... populate lists R -> E
      If(flagR2E) Then
         Call memory(inER, iEnd, nR)
         Call memory(iIER, iEnd, 6*nE)
         If(iEnd.GT.MaxWi) Goto 5000

         Call backReferences(nR, nE, 6,6, iW(iIRE), iW(inER), iW(iIER))
      End if


c ... populate lists (P->P) = (P->E) x (E->P)
      If(flagP2P) Then
         MaxX= 25 * nP

         inPP = iEnd
         iIPP = inPP + nP
 100     iiW  = iIPP + MaxX
         iEnd = iiW  + nP
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nP, nP, nE, iW(inEP), iW(iIEP), 4, IPE, 
     &                 nX, MaxX, iW(inPP), iW(iIPP), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 5 * nP
            Goto 100
         End if

         iEnd = iIPP + nX 
      End if
 

c ... populate lists (P->R) = (P->E) x (E->R)
      If(flagP2R) Then
         MaxX = 5 * nR

         inRP = iEnd
         iIRP = inRP + nP
 200     iiW  = iIRP + MaxX
         iEnd = iiW  + nR
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nP, nR, nE, iW(inEP), iW(iIEP), 6, iW(iIRE), 
     &                 nX, MaxX, iW(inRP), iW(iIRP), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 2 * nR
            Goto 200
         End if

         iEnd = iIRP + nX 
      End if
 

c ... populate lists (P->F) = (P->E) x (E->F)
      If(flagP2F) Then
         MaxX = 5 * nF

         inFP = iEnd
         iIFP = inFP + nR
 300     iiW  = iIFP + MaxX
         iEnd = iiW  + nF
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nP, nF, nE, iW(inEP), iW(iIEP), 4, iW(iIFE), 
     &                 nX, MaxX, iW(inFP), iW(iIFP), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 2 * nF
            Goto 300
         End if

         iEnd = iIFP + nX 
      End if


c ... populate lists (R->P) = (R->E) x (E->P)
      If(flagR2P) Then
         MaxX = 5 * nR

         inPR = iEnd
         iIPR = inPR + nR
 400     iiW  = iIPR + MaxX
         iEnd = iiW  + nX
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nR, nP, nE, iW(inER), iW(iIER), 4, IPE,
     &                 nX, MaxX, iW(inPR), iW(iIPR), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 2 * nR
            Goto 400
         End if

         iEnd = iIPR + nX 
      End if


c ... populate lists (R->R) = (R->E) x (E->R)
      If(flagR2R) Then
         MaxX = 10 * nR

         inRR = iEnd
         iIRR = inRR + nR
 500     iiW  = iIRR + MaxX
         iEnd = iiW  + nR
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nR, nR, nE, iW(inER), iW(iIER), 6, iW(iIRE),
     &                 nX, MaxX, iW(inRR), iW(iIRR), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 2 * nR
            Goto 500
         End if

         iEnd = iIRR + nX
      End if


c ... populate lists (R->F) = (R->E) x (E->F)
      If(flagR2F) Then
         MaxX = 10 * nR

         inFR = iEnd
         iIFR = inFR + nR
 600     iiW  = iIFR + MaxX
         iEnd = iiW  + nF 
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nR, nF, nE, iW(inER), iW(iIER), 4, iW(iIFE),
     &                 nX, MaxX, iW(inFR), iW(iIFR), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + 2 * nR
            Goto 600
         End if

         iEnd = iIFR + nX
      End if


c ... populate lists (F->P) = (F->E) x (E->P)
      If(flagF2P) Then
         MaxX = 3 * nF

         inPF = iEnd
         iIPF = inPF + nR
 700     iiW  = iIPF + MaxX
         iEnd = iiW  + nP
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nF, nP, nE, iW(inEF), iW(iIEF), 4, IPE,
     &                 nX, MaxX, iW(inPF), iW(iIPF), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + nF
            Goto 700
         End if

         iEnd = iIPF + nX
      End if


c ... populate lists (F->F) = (F->E) x (E->F)
      If(flagF2F) Then
         MaxX = 7 * nF

         inFF = iEnd
         iIFF = inFF + nF
 800     iiW  = iIFF + MaxX
         iEnd = iiW  + nF
         If(iEnd.GT.MaxWi) Goto 5000

         Call listConv(nF, nF, nE, iW(inEF), iW(iIEF), 4, iW(iIFE),
     &                 nX, MaxX, iW(inFF), iW(iIFF), iW(iiW), iERR)

         If(iERR.GT.0) Then
            MaxX = MaxX + nF
            Goto 800
         End if

         iEnd = iIFF + nX
      End if


c ... create map E -> Fb
      iIFbE = iEnd
      inFbP = iIFbE + 4 * nE
      iIFbP = inFbP + nP
      iEnd  = iIFbP + 3 * nFb 
      If(iEnd.GT.MaxWi) Goto 5000
      
      Call listE2Fb(nP, nFb, nE, IPF, IPE, 
     &              iW(iIFbE), iW(inFbP), iW(iIFbP))


      nRow = irp * nP + irf * nF + irr * nR + ire * nE
      nCol = icp * nP + icf * nF + icr * nR + ice * nE

      If(nRow.GT.MaxIA) Call errMesFEM(1013, 'BilinearFormTemplate',
     &                        'local parameter MaxIA is small')

      iusage = (iEnd * 100) / MaxWi
      rusage = 0


c ... the symbolic assemble of the matrix 
      k = 1
      IA(1) = 1

      Do m = 1, irp
         ip1 = 0
         if1 = 0
         ir1 = 0
         ie1 = 0

         Do n = 1, nP 
            ip2 = iW(inPP + n - 1) 
            if2 = iW(inFP + n - 1) 
            ir2 = iW(inRP + n - 1) 
            ie2 = iW(inEP + n - 1) 

            iEnd = IA(k) + icp * (ip2 - ip1) 
     &                   + icr * (ir2 - ir1)
     &                   + icf * (if2 - if1)
     &                   + ice * (ie2 - ie1)
            If(iEnd.GT.MaxA) Goto 6000

            p = IA(k)
            Call fillJA(
     &           p, JA, 
     &           nP, nR, nF, nE, iW(iIPP), iW(iIRP), iW(iIFP), iW(iIEP),
     &           icp, icr, icf, ice, 
     &           ip1, ir1, if1, ie1, ip2, ir2, if2, ie2)

            k = k + 1
            IA(k) = iEnd

            ip1 = ip2
            ir1 = ir2
            if1 = if2
            ie1 = ie2
         End do
      End do

      Do m = 1, irr
         ip1 = 0
         ir1 = 0
         if1 = 0
         ie1 = 0

         Do n = 1, nR 
            ip2 = iW(inPR + n - 1)
            ir2 = iW(inRR + n - 1) 
            if2 = iW(inFR + n - 1) 
            ie2 = iW(inER + n - 1) 

            iEnd = IA(k) + icp * (ip2 - ip1)
     &                   + icr * (ir2 - ir1)
     &                   + icf * (if2 - if1)
     &                   + ice * (ie2 - ie1)
            If(iEnd.GT.MaxA) Goto 6000

            p = IA(k)
            Call fillJA(
     &           p, JA, 
     &           nP, nR, nF, nE, iW(iIPR), iW(iIRR), iW(iIFR), iW(iIER),
     &           icp, icr, icf, ice, 
     &           ip1, ir1, if1, ie1, ip2, ir2, if2, ie2)

            k = k + 1
            IA(k) = iEnd

            ip1 = ip2
            ir1 = ir2
            if1 = if2
            ie1 = ie2
         End do
      End do

      Do m = 1, irf
         ip1 = 0
         ir1 = 0
         if1 = 0
         ie1 = 0

         Do n = 1, nF 
            ip2 = iW(inPF + n - 1)
            ir2 = iW(inRF + n - 1) 
            if2 = iW(inFF + n - 1) 
            ie2 = iW(inEF + n - 1) 

            iEnd = IA(k) + icp * (ip2 - ip1)
     &                   + icr * (ir2 - ir1)
     &                   + icf * (if2 - if1)
     &                   + ice * (ie2 - ie1)
            If(iEnd.GT.MaxA) Goto 6000

            p = IA(k)
            Call fillJA(
     &           p, JA, 
     &           nP, nR, nF, nE, iW(iIPF), iW(iIRF), iW(iIFF), iW(iIEF),
     &           icp, icr, icf, ice, 
     &           ip1, ir1, if1, ie1, ip2, ir2, if2, ie2)

            k = k + 1
            IA(k) = iEnd

            ip1 = ip2
            ir1 = ir2
            if1 = if2
            ie1 = ie2
         End do
      End do

      Do m = 1, ire
         Do n = 1, nE 
            ip1 = 4 * (n - 1)
            ir1 = 6 * (n - 1)
            if1 = 4 * (n - 1)
            iIEE(1) = n

            ip2 = ip1 + 4 
            ir2 = ir1 + 6
            if2 = if1 + 4

            iEnd = IA(k) + icp * 4 + icr * 6 + icf * 4 + ice
            If(iEnd.GT.MaxA) Goto 6000

            p = IA(k)
            Call fillJA(
     &           p, JA, 
     &           nP, nR, nF, nE, IPE, iW(iIRE), iW(iIFE), iIEE,
     &           icp, icr, icf, ice, 
     &           ip1, ir1, if1, 0, ip2, ir2, if2, 1)

            k = k + 1
            IA(k) = iEnd
         End do
      End do

      iEnd = IA(nRow + 1) - 1


C ... comply with the AMG format
      If(flagAMG) Then
         Do 900 n = 1, nRow
            Do j = IA(n), IA(n + 1) - 1
               If(JA(j).EQ.n) Then
                  Call swapii(JA(j), JA(IA(n)))
                  Goto 900
               End if
            End do
 900     Continue
      End if


C ... comply with the AMG format
      If(flagCSR) Then
         Do n = 1, nRow
            j1 = IA(n)
            j2 = IA(n + 1)
            Call order2D(j2 - j1, 1, JA(j1))
         End do 
      End if


C ... fill in the sparsity structure
      isrR = irp * nP
      isrF = isrR + irr * nR
      isrE = isrF + irf * nF

      iscR = icp * nP
      iscF = iscR + icr * nR
      iscE = iscF + icf * nF
      
      Do n = 1, nRow
         F(n) = 0D0
      End do

      Do n = 1, IA(nRow + 1) - 1
         A(n) = 0D0
      End do

      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

c  ...   save edge labels
         Do i = 1, 4
            lbFloc(i) = 0

            iFt = iW(iIFbE + 4 * (n - 1) + i - 1)
            If(iFt.GT.0) lbFloc(i) = lbF(iFt)
         End do

         Do i = 1, 4
            iPt = IPE(i, n)
            lbPloc(i) = lbP(iPt)
         End do

         iiR = iIRE + 6*(n-1) + 1
         iiF = iIFE + 4*(n-1) + 1
         Call encodeISYS(n, iP1,iP2,iP3,iP4, iW(iiR), iW(iiF),
     &                      nP, nR, nF, nE, iSYS)

         Call fem3Dext(
     &        XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, iP4),
     &        lbE(n), lbFloc, lbRloc, lbPloc, DATA, iSYS,
     &        MaxSize, Aloc, Floc, ir, ic, 
     &        templateR, templateC)

c  ...  create arrays of row indices 
         irp = 0
         irr = 0
         irf = 0
         ire = 0

         i = 1
         Do While(i.LE.ir)
            If(templateR(i).EQ.Vdof) Then
               indR(i)     = irp * nP + iP1 
               indR(i + 1) = irp * nP + iP2 
               indR(i + 2) = irp * nP + iP3 
               indR(i + 3) = irp * nP + iP4 

               i = i + 4
               irp = irp + 1
            Else If(templateR(i).EQ.Rdof .OR.
     &              templateR(i).EQ.RdofOrient) Then
               Do k = 0, 5
                  indR(i + k) = isrR + irr * nR + iW(iIRE + 6*(n-1) + k)
               End do

               i = i + 6
               irr = irr + 1
            Else If(templateR(i).EQ.Fdof .OR. 
     &              templateR(i).EQ.FdofOrient) Then
               Do k = 0, 3
                  indR(i + k) = isrF + irf * nF + iW(iIFE + 4*(n-1) + k)
               End do

               i = i + 4
               irf = irf + 1
            Else If(templateR(i).EQ.Edof) Then
               indR(i) = isrE + ire * nE + n 

               i = i + 1
               ire = ire + 1
            End if
         End do

c  ...  create arrays of column indices 
         icp = 0
         icr = 0
         icf = 0
         ice = 0

         i = 1
         Do While(i.LE.ic)
            If(templateC(i).EQ.Vdof) Then
               indC(i)     = icp * nP + iP1 
               indC(i + 1) = icp * nP + iP2 
               indC(i + 2) = icp * nP + iP3 
               indC(i + 3) = icp * nP + iP4 

               i = i + 4
               icp = icp + 1
            Else If(templateC(i).EQ.Rdof .OR.
     &              templateC(i).EQ.RdofOrient) Then
               Do k = 0, 5
                  indC(i + k) = iscR + icr * nR + iW(iIRE + 6*(n-1) + k)
               End do

               i = i + 6
               icr = icr + 1
            Else If(templateC(i).EQ.Fdof .OR.
     &              templateC(i).EQ.FdofOrient) Then
               Do k = 0, 3
                  indC(i + k) = iscF + icf * nF + iW(iIFE + 4*(n-1) + k)
               End do

               i = i + 4
               icf = icf + 1
            Else If(templateC(i).EQ.Edof) Then
               indC(i) = iscE + ice * nE + n 

               i = i + 1
               ice = ice + 1
            End if
         End do


c  ...  correct orientation of the basis functions in the local matrix
c  ...  F-unknowns: the global orientation is from bigger E to smaller E
c  ...  R-unknowns: the global orientation is from bigger P to smaller P
         i = 1
         Do While(i.LE.ir)
            If(templateR(i).EQ.FdofOrient) Then
               Do k = 0, 3
                  iFt = iW(iIFE + 4*(n-1) + k)

                  iEt = findY(iFt, iW(iIEF), iW(inEF), n)

                  If(n.LT.iEt) Then 
                     Do j = 1, ic
                        Aloc(i, j) = -Aloc(i, j)
                     End do 
                     Floc(i) = -Floc(i)
                  End if

                  i = i + 1
               End do
            Else If(templateR(i).EQ.RdofOrient) Then
               Do i1 = 1, 3
                  Do i2 = i1 + 1, 4
                     iP1 = IPE(i1, n)
                     iP2 = IPE(i2, n)

                     If(iP1.LT.iP2) Then 
                        Do j = 1, ic
                           Aloc(i, j) = -Aloc(i, j)
                        End do 
                        Floc(i) = -Floc(i)
                     End if

                     i = i + 1
                  End do
               End do
            Else
               i = i + 1
            End if
         End do 
   

c  ...  correct orientation of the basis functions in the local matrix
c  ...  F-unknowns: the global orientation is from bigger E to smaller E
c  ...  R-unknowns: the global orientation is from bigger P to smaller P
         i = 1
         Do While(i.LE.ic)
            If(templateC(i).EQ.FdofOrient) Then
               Do k = 0, 3
                  iFt = iW(iIFE + 4*(n-1) + k)

                  iEt = findY(iFt, iW(iIEF), iW(inEF), n)

                  If(n.LT.iEt) Then 
                     Do j = 1, ir
                        Aloc(j, i) = -Aloc(j, i)
                     End do 
                  End if

                  i = i + 1
               End do
            Else If(templateC(i).EQ.RdofOrient) Then
               Do i1 = 1, 3
                  Do i2 = i1 + 1, 4
                     iP1 = IPE(i1, n)
                     iP2 = IPE(i2, n)

                     If(iP1.LT.iP2) Then 
                        Do j = 1, ir
                           Aloc(j, i) = -Aloc(j, i)
                        End do 
                     End if

                     i = i + 1
                  End do
               End do
            Else
               i = i + 1
            End if
         End do 
   

c  ...  assemble the right hand side
         Do i = 1, ir
            iP = indR(i)
            F(iP) = F(iP) + Floc(i)
         End do


c  ...  assemble the elemental matrix
         Do i = 1, ir
            iP = indR(i)

            Do 1000 j = 1, ic
               iQ = indC(j)

               Do k = IA(iP), IA(iP + 1) - 1
                  If(JA(k).EQ.iQ) Then
                     A(k) = A(k) + Aloc(i, j)
                     Goto 1000
                  End if
               End do
 1000       Continue
         End do
      End do
 

c ... print statistics
      non_zero = IA(nRow + 1) - 1 

      Write(*,'(5X,A,I7,A,I7,A,I9,A)') 'matrix size:', nRow,' x', nCol, 
     &                            '  (', non_zero, ' non-zero entries)'

      Write(*,'(A,3(I2,A))') '     memory usage: ', 
     &      iusage, '% (Integer) and ', rusage, '% (Real*8)'


      Return

c ... error messages
 5000 Continue
      Write(message,'(A,I10)')
     &    'The approximate size of iW is ', iEnd
      Call errMesFEM(1001, 'BilinearFormTemplate', message)

 6000 Continue 
      Call errMesFEM(1014,
     &    'BilinearFormTemplate', 'local parameter MaxA is small')

      Return
      End



C ======================================================================
      Subroutine fillJA(p, JA, 
     &                  nP,  nR,  nF, nE,  IPX, IRX, IFX, IEX,
     &                  icp, icr, icf, ice, 
     &                  ip1, ir1, if1, ie1, ip2, ir2, if2, ie2)
C ======================================================================
C  Technical subroutine for incorporating a part of the local 
C  sparsity structure into the global one.
C ======================================================================
      Integer p, JA(*)
      Integer IPX(*), IRX(*), IFX(*), IEX(*)

C ====================================================================== 
      isR = icp * nP
      isF = isR + icr * nR
      isE = isF + icf * nF

      Do l = 1, icp
         ishift = (l - 1) * nP

         Do j = ip1 + 1, ip2
            JA(p) = ishift + IPX(j)
            p = p + 1
         End do
      End do

      Do l = 1, icr
         ishift = isR + (l - 1) * nR

         Do j = ir1 + 1, ir2
            JA(p) = ishift + IRX(j)
            p = p + 1
        End do
      End do

      Do l = 1, icf
         ishift = isF + (l - 1) * nF

         Do j = if1 + 1, if2
            JA(p) = ishift + IFX(j)
            p = p + 1
        End do
      End do

      Do l = 1, ice
         ishift = isE + (l - 1) * nE

         Do j = ie1 + 1, ie2
            JA(p) = ishift + IEX(j)
            p = p + 1
         End do
      End do

      Return
      End



C ====================================================================== 
      Subroutine initPointers(irX, irY, icX, icY, 
     &                        inXY, iIXY, flagXY, inYX, iIYX, flagYX)
C ====================================================================== 
C We need to initialize pointers to zero to avoid passing wrong
C pointers even if arrays to which they point are not used.
C ====================================================================== 
      Logical flagXY, flagYX

      inXY = 1
      iIXY = 1

      flagXY = .FALSE.
      If(irX.GT.0 .AND. icY.GT.0)  flagXY = .TRUE.

      inYX = 1
      iIYX = 1

      flagYX = .FALSE.
      If(irY.GT.0 .AND. icX.GT.0)  flagYX = .TRUE.

      Return 
      End



C ====================================================================== 
      Subroutine memory(ipointer, ifree, isize)
C ====================================================================== 
C Some pointers were initialized to one to point to a valid memory
C Therefore the free memory starts at 2.
C ====================================================================== 
      If(ipointer.LE.1) Then
         ipointer = ifree
         ifree = ifree + isize
      End if

      Return
      End



C ======================================================================
      Integer Function findY(iX, IYX, nYX, iY)
C ======================================================================
C findY returns pointer to object Y, other than iY, in the sublist 
c of X(iX) -> {iY1, iY2, ..., iYk} of map X -> Y.
C If not found, findY = 0.
C ======================================================================
      Integer IYX(*), nYX(*)

      If(iX.EQ.1) Then
         ib = 1
      Else
         ib = nYX(iX - 1) + 1
      End if
      ie = nYX(iX)

      Do i = ib, ie
         iYt = IYX(i)
         If(iYt.NE.iY) Then
            findY = iYt
            Goto 1000
         End if
      End do

      findY = 0
 1000 Return
      End


