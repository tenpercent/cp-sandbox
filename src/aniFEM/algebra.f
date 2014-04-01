c ======================================================================
      Subroutine mulAgen(N, IA, JA, A, U, AU)
c ======================================================================
c  Computing  AU = A * U 
c ======================================================================
      Integer N, IA(*), JA(*)
      Real*8  A(*), U(*), AU(*)
c ======================================================================
      Do k = 1, N
         AU(k) = 0D0
      End do

      Do k = 1, N
         Do i = IA(k), IA(k + 1) - 1
            AU(k) = AU(k) + A(i) * U(JA(i))
         End do
      End do
      Return
      End



c ======================================================================
      Subroutine AMG2CSC(NA, IA, JA, A, NB, IB, JB, B)
c ======================================================================
c  Converting the AMG format (CSR) to the compressed sparse 
C  column format (CSC).
C
C  *** Note: 1. size(IB) = max(JA) + 1
C            2. size(JB) = size(JA)
C            3. size(B) = size(A) 
c ======================================================================
      Integer NA, IA(*), JA(*)
      Integer NB, IB(*), JB(*)
      Real*8  A(*), B(*)

      M = IA(NA + 1) - 1

c ... compute the number of rows
      NB = 0
      Do i = 1, M
         NB = max(NB, JA(i))
      End do


c ... compute the number of non-zero elements in each row
      Do i = 1, NB
         IB(i) = 0
      End do

      Do i = 1, M
         j = JA(i)
         IB(j) = IB(j) + 1 
      End do 

      Do i = 1, NB
         IB(i + 1) = IB(i + 1) + IB(i)
      End do
      IB(NB + 1) = IB(NB) + 1


c ... copy matrix A to matrix B
      Do i = NA, 1, -1
         Do j = IA(i), IA(i + 1) - 1
            k = JA(j)

            l = IB(k)
            IB(k) = l - 1

            JB(l) = i
             B(l) = A(j)
         End do
      End do

      Do i = 1, NB
         IB(i) = IB(i) + 1
      End do

      Return
      End



c ======================================================================
      Subroutine ZeroBasedFormat(N, IA, JA)
c ======================================================================
C Routine converts sparcity format to the 0-based format.
c ======================================================================
      Integer N, IA(*), JA(*)

      Do i = 1, N + 1
         IA(i) = IA(i) - 1
      End do

      Do i = 1, IA(N + 1)
         JA(i) = JA(i) - 1
      End do

      Return
      End



c ======================================================================
      Subroutine CSC2AMG(NA, IA, JA, A, NB, IB, JB, B)
c ======================================================================
c  Converting the compressed 0-based column format (CSC) to the
C  AMG (row-wice) format.
C
C  *** Note: 1. size(IB) = max(JA) + 1
C            2. size(JB) = size(JA)
C            3. size(B) = size(A) 
c ======================================================================
      Integer NA, IA(*), JA(*)
      Integer NB, IB(*), JB(*)
      Real*8  A(*), B(*)

      Call AMG2CSC(NA, IA, JA, A, NB, IB, JB, B)

      Do 10 i = 1, NB
         j1 = IB(i)
         j2 = IB(i + 1) - 1

         Do j = j1, j2
            k = JB(j)
            If(k.EQ.i) Then
               Call swapii(JB(j), JB(j1)) 
               Call swapdd( B(j),  B(j1))
               Goto 10
            End if
         End do
 10   Continue

      Return
      End



C ======================================================================
      Subroutine comABgen(nA, IA, JA, A, 
     &                    nB, IB, JB, B, 
     &                    nC, IC, JC, C, MaxC, iW)
C ======================================================================
c  Computing  C = A * B where all matrices are in AMG sparse row format
C ======================================================================
      Integer IA(*), JA(*), IB(*), JB(*), IC(*), JC(*)
      Real*8  A(*), B(*), C(*)

      Integer iW(*)
C ======================================================================
      nC = nA

      mC = 0
      iC(1) = 1

c ... the number of columns of B
      nColB = 0
      Do n = 1, IB(nB+1) - 1
         nColB = max(nColB, JB(n)) 
      End do
      Do n = 1, nColB
         iW(n) = 0
      End do
      
      Do i = 1, nA
         If(mC + nColB.GT.MaxC) Call errMesFEM(1004, 
     &        'algebra', 'Not enough memory for the matrix')

         mC = mC + 1
         JC(mC) = i
         C(mC)  = 0D0
         iW(i)  = mC

         Do n = IA(i), IA(i + 1) - 1
            k = JA(n)
            Do l = IB(k), IB(k + 1) - 1
               j = JB(l)
               m = iW(j)
               If(m.NE.0) Then 
                  C(m) = C(m) + A(k) * B(l)
               Else
                  mC = mC +  1
                  JC(mC) = j
                  C(mC)  = A(k) * B(l)
                  iW(j)  = mC
               End if
            End do
         End do
         IC(i + 1) = mC + 1

         Do k = IC(i), mC 
            iW(JC(k)) = 0
         End do
      End do

      Return
      End



c ======================================================================
      Subroutine addBgen(nA, IA, JA, A, nB, IB, JB, B)
c ======================================================================
c  Computing  A = A + B. The sparcity structure of A must be
c  bigger than that of B.
c ======================================================================
      Integer IA(*), JA(*), IB(*), JB(*)
      Real*8  A(*), B(*)

c ======================================================================
      Do i = 1, nB
         k1 = IA(i)   
         k2 = IA(i + 1) + 1

         j1 = IB(i)
         j2 = IB(i + 1) - 1

         Do j = j1, j2
            Call findSE(k2 - k1, JA(k1), JB(j), k)

            If(k.EQ.0) Call errMesFEM(2011, 'algebra',
     &         'Structure of B bigger than structure of A.')

            k = k1 + k - 1
            A(k) = A(k) + B(j)
         End do
      End do
 
      Return
      End



c ======================================================================
      Subroutine printA(nA, IA, JA, A, F)
c ======================================================================
      Integer IA(*), JA(*)
      Real*8  A(*), F(*), asum, amin, amax
c ======================================================================
      Do i = 1, nA
         j1 = IA(i)   
         j2 = IA(i + 1) - 1

         amin = A(j1)
         amax = A(j1)
         asum = 0d0
         Do j = j1, j2
            amin = min(amin, A(j))
            amax = max(amax, A(j))
            asum = asum + A(j)
         End do
         Write(*, '(A,I3,A,200I9)') 'Row[', i, ']:', (JA(j), j=j1,j2)
         Write(*, '(A,200E9.1)')    '         ', (A(j), j=j1,j2)
         Write(*, '(A,3E12.4,A,E12.4)') 
     &      'Row sum, min/max =', asum, amin, amax, '  f=', F(i)
      End do
 
      Return
      End



c ======================================================================
      Subroutine printMap(n, nXY, IXY)
c ======================================================================
      Integer nXY(*), IXY(*)

c ======================================================================
      j2 = 0
      Do i = 1, n
         j1 = j2 + 1
         j2 = nXY(i)

         Write(*, '(A,I3,A,100I6)') 'Row[', i, ']:', (IXY(j), j=j1,j2)
      End do
 
      Return
      End



C ======================================================================
      Subroutine invertSPDmatrix(N, A, LDA)
C ======================================================================
C  Routine inverts symmetric positive definite matrix A of order N.
C  It is a driver for two LAPACK routines.
C ======================================================================
      Integer  N, LDA
      Real*8   A(LDA, *)
C ======================================================================
      INFO = 0

      call DPOTRF( 'U', N, A, LDA, INFO )
      If(INFO.NE.0) Stop 'DPOTRF in invertSPDmatrix failed'

      call DPOTRI( 'U', N, A, LDA, INFO )
      If(INFO.NE.0) Stop 'DPOTRI in invertSPDmatrix failed'

      Do i = 1, N
         Do j = i+1, N
            A(j, i) = A(i, j)
         End do
      End do

      Return
      End


