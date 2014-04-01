C ======================================================================
c Routine generates isotropic metric from element-based error estimates.
c The nodal metric is recovered with the ZZ interpolation algorithm. 
C ======================================================================
      Subroutine CellEst2MetricZZ(nP, nF, nE, XYP, IPE, IPF,  
     &                            Error, Metric,       
     &                            MaxWr, MaxWi, rW, iW) 
C ======================================================================
c Input:
      Integer nP, nF, nE   ! numbers of nodes, element, boundary edges
      Real*8  XYP(3, *)    ! coordinates of mesh nodes
      Integer IPE(4, *)    ! connectivity table for elements
      Integer IPF(3, *)    ! boundary edges data
      Real*8  Error(*)     ! element-based error estimates

c Output:
      Real*8  Metric(6, *) ! node-based metric

c Working arrays:
      Integer MaxWr, MaxWi
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)

C ================================================================
c ... memory test
      iEnd = 3 * nP + 4 * nE
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 'CellEst2MetricZZ', 
     &                      'increase size of work memory MaxWi')

c ... ZZ-interpolation of element-based metric into nodes
      MaxWrWork = MaxWr - nP

      Call P02P1(nP, nF, nE, XYP, IPF, IPE,
     &           Error, rW,
     &           MaxWrWork, MaxWi, rW(nP+1), iW)


c ... create the metric
      Do i = 1, nP
         Metric(1, i) = rW(i)
         Metric(2, i) = Metric(1, i)
         Metric(3, i) = Metric(1, i)
         Metric(4, i) = 0
         Metric(5, i) = 0
         Metric(6, i) = 0
      End do

      Return
      End

