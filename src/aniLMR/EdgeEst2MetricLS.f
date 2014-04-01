C ======================================================================
C Routine generates nodal tensor metric from edge-based error estimates
C using the Least Square solution of the following system
C
C    (M(a_i) e_k,e_k) = \eta_k
C
C for all edges k incident to a mesh node a_i.
C ======================================================================
      Subroutine EdgeEst2MetricLS(nP, nE, XYP, IPE,   
     &                            Error, Metric,       
     &                            MaxWr, MaxWi, rW, iW) 
C ======================================================================
      include 'makS.fd'
      include 'magic.fd'
C ======================================================================
c Input:
      Integer nP, nE       ! number of nodes and elements
      Real*8  Error(6, *)  ! edge error estimates data
      Real*8  XYP(3, *)    ! coordinates of mesh nodes
      Integer IPE(4, *)    ! connectivity table for tets

c Output:
      Real*8  Metric(6, *) ! node-based metric 

c Working arrays:
      Integer MaxWr, MaxWi
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)

c Local variables
      Real*8  ebuf(MaxS), vbuf(3, MaxS), det

C ======================================================================
c ... form list edge -> elements
      iIRE = 1
      inEP = iIRE + 6 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 4 * nE
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))


c ... collect errors on edges
      Do i = 1, nE
         Do j = 1, 6
            iR = iW(iIRE + (i-1)*6 + j-1)
            rW(iR) = error(j, i)
         End do
      End do


c ... form list edges -> points
      MaxR = nR
      iIPR = 1
      inEP = iIPR + 2 * MaxR
      iIEP = inEP + nP
      iEnd = iIEP + 4 * nE + nP
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call listR2P(nP, nR, nE, MaxR, IPE, iW(iIPR), iW(inEP), iW(iIEP))


c ... form list points -> edges 
      inRP = inEP
      iIRP = inRP + nP
      iEnd = iIRP + 2 * nR
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call backReferences(nP, nR, 2,2, iW(iIPR), iW(inRP), iW(iIRP))


c ... compute nodal values of metric
      i2 = 0
      Do i = 1, nP
         i1 = i2 + 1
         i2 = iW(inRP + i - 1)

         m = 0
         Do j = i1, i2
            iR  = iW(iIRP + j - 1)

            iP1 = iW(iIPR + (iR-1) * 2)
            iP2 = iW(iIPR + (iR-1) * 2 + 1)

            m = m + 1
            ebuf(m) = rW(iR)
            Do k = 1, 3
               vbuf(k, m) = XYP(k, iP1) - XYP(k, iP2)
            End do
         End do    

         Call NodalMetric(m, ebuf, vbuf, Metric(1, i))

c ... take spectral modulo of metric, |M|
         Call SpectralModule(Metric(1, i), det)
      End do

      Return
      End
      


C ======================================================================
      Subroutine NodalMetric(k, values, xy, metric)
C ======================================================================
C  This routine uses least square linear solution to the system
C  (Metric xy_i, xy_i) = values_i, i=1,\dots,k
C ======================================================================
      implicit none

      Integer  k
      Real*8   xy(3,*), values(*), metric(6)

c (local variables)
      Real*8   A(6, 6), S(6), work(30)
      Integer  i, j, ipiv(6), info

C ======================================================================
      Do i = 1, 6
         Do j = 1, 6
            A(i, j) = 0D0
         End do
         S(i) = 0D0
      End do

c ... generate the least squares matrix
      Do i = 1, k
         A(1,1) = A(1,1) + xy(1, i)**4 
         A(1,2) = A(1,2) + xy(1, i)**2 * xy(2, i)**2
         A(1,3) = A(1,3) + xy(1, i)**2 * xy(3, i)**2
         A(1,4) = A(1,4) + 2 * xy(1, i)**3 * xy(2, i)
         A(1,5) = A(1,5) + 2 * xy(2, i)**2 * xy(2, i) * xy(3, i)
         A(1,6) = A(1,6) + 2 * xy(1, i)**3 * xy(3, i)

         A(2,2) = A(2,2) + xy(2, i)**4 
         A(2,3) = A(2,3) + xy(2, i)**2 * xy(3, i)**2
         A(2,4) = A(2,4) + 2 * xy(2, i)**3 * xy(1, i)
         A(2,5) = A(2,5) + 2 * xy(2, i)**3 * xy(3, i)
         A(2,6) = A(2,6) + 2 * xy(2, i)**2 * xy(1, i) * xy(3, i)

         A(3,3) = A(3,3) + xy(3, i)**4 
         A(3,4) = A(3,4) + 2 * xy(3, i)**2 * xy(1, i) * xy(2, i) 
         A(3,5) = A(3,5) + 2 * xy(3, i)**3 * xy(2, i) 
         A(3,6) = A(3,6) + 2 * xy(3, i)**3 * xy(1, i) 

         A(4,4) = A(4,4) + 4 * xy(1, i)**2 * xy(2, i)**2
         A(4,5) = A(4,5) + 4 * xy(1, i) * xy(2, i)**2 * xy(3, i)
         A(4,6) = A(4,6) + 4 * xy(2, i) * xy(1, i)**2 * xy(3, i)

         A(5,5) = A(5,5) + 4 * xy(2, i)**2 * xy(3, i)**2
         A(5,6) = A(5,6) + 4 * xy(1, i) * xy(3, i)**2 * xy(2, i)

         A(6,6) = A(6,6) + 4 * xy(1, i)**2 * xy(3, i)**2
      End do

c ... generate the RHS
      Do i = 1, k
         S(1) = S(1) + xy(1, i)**2 * values(i)
         S(2) = S(2) + xy(2, i)**2 * values(i)
         S(3) = S(3) + xy(3, i)**2 * values(i)

         S(4) = S(4) + 2 * xy(1, i) * xy(2, i) * values(i)
         S(5) = S(5) + 2 * xy(2, i) * xy(3, i) * values(i)
         S(6) = S(6) + 2 * xy(1, i) * xy(3, i) * values(i)
      End do

c ... fix the matrix (a bug???)
      Do i = 1, 6
         Do j = i + 1, 6
            A(j, i) = A(i, j)
         End do
      End do

      Do i = 1, 6
         A(i, i) = A(i, i) + 1e-16
      End do

      Call dsysv('U', 6, 1, A, 6, ipiv, S, 6, work, 30, info)

      If(info.NE.0) Call errMesLMR(1001, 
     &                  'NodalMetric', 'Lapack routine dsysv')

c ... write the metric
      Do i = 1, 6
         metric(i) = S(i)
      End do

      Return
      End

