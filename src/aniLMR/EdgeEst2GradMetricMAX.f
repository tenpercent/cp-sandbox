C ======================================================================
      Subroutine EdgeEst2GradMetricMAX(Error,  
     &                                 nP, nE, XYP, IPE, 
     &                                 metric, MaxWr, rW)
C ======================================================================
C  Routines calculates a nodal metric using the method of shifts. 
C  The input data are error estimates prescribed to mesh edges.
C ======================================================================
C  Input:
C     nP  - the number of nodes
C     nE  - the number of tetrahedra
C
C     Error(6, nE) - error estimates given on edges
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table for tets
C
C  Output: 
C     Metric(6, nE) - nodal tensor metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr, MaxWr >= nP
C ======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(3, *)
      Integer  IPE(4, *)

C Function and metric 
      Real*8   Error(6, *),  Metric(6, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, n, iP1, iP2, iP3, iP4
      Real*8   det, H(6)

C ======================================================================
c ... memory check
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'EdgeEst2GradMetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do

c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call TetGradMetric(XYP(1, iP1), XYP(1, iP2), 
     &                      XYP(1, iP3), XYP(1, iP4), Error, H)

         det = H(1) * (H(2) * H(3) - H(5) ** 2)
     &       - H(4) * (H(4) * H(3) - H(5) * H(6))
     &       + H(6) * (H(4) * H(5) - H(6) * H(2))

c  ...  update the maximum metric
         Do i = 1, 4
            iP1 = IPE(i, n)

            If(det.GT.rW(iP1)) Then
                rW(iP1) = det
                Do j = 1, 6
                   Metric(j, iP1) = H(j)
                End do
            End if
         End do
      End do

      Return
      End



C=======================================================================
      Subroutine EdgeEst2GradMetricCell(Error,  
     &                                  nP, nE, XYP, IPE, 
     &                                  Metric, MaxWr, rW)
C=======================================================================
C  Routine calculates the cell-based metric using the edge-based
C  error estimates Error.
C=======================================================================
C  Input:
C     nP  - the number of nodes
C     nE  - the number of tetrahedra
C
C     Error(6, nE)  - given error estimates on edges
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table for tets
C
C  Output: 
C     Metric(6, nE) - cell metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr, MaxWr >= nP
C=======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(3, *)
      Integer  IPE(4, *)

C Function and metric 
      Real*8   Error(6, *), Metric(6, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, n, iP1, iP2, iP3, iP4
      Real*8   H(6)

C=======================================================================
c ... memory test
      If(nP.GT.MaxWr) Call errMesLMR(1002, 'EdgeEst2GradMetricCell', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do

c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call TetGradMetric(XYP(1, iP1), XYP(1, iP2), 
     &                      XYP(1, iP3), XYP(1, iP4), Error, H)

         Do i = 1, 6
            Metric(i, n) = H(i)
         End do
      End do

      Return
      End



C=======================================================================
      Subroutine TetGradMetric(xy1, xy2, xy3, xy4, Error, H)
C=======================================================================
C The routine computes the gradient metric H for a tetrahedron given by
C its vertices. 
C=======================================================================
      implicit none
      include 'fem3Dtet.fd'
      include 'assemble.fd'

C=======================================================================
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3)
      Real*8   Error(6), H(6)

c Local variables
      Integer  i, j, label
      Real*8   B(10, 10), s, p, DATAFEM

      Integer  ANI_Dnull, iSYS(MAXiSYS)
      EXTERNAL ANI_Dnull

C=======================================================================
c ... calculate matrix B = \grad b_i . \grad b_j where b_i is the bubble
      label = 1
      Call  FEM3Dtet(
     &      xy1, xy2, xy3, xy4,
     &      GRAD, FEM_P2, GRAD, FEM_P2, 
     &      label, ANI_Dnull, DATAFEM, iSYS, 2, 
     &      10, B, i, j)

c ... compute alpha's using the given discrete function. 
      p = 0D0
      s = 0D0
      Do i = 1, 6
         Do j = 1, 6
            s = s + B(4+i, 4+j) * Error(i) * Error(j)
         End do
         p = p + Error(i)
      End do

      s = s / p

c ... gamma's (H) are scaled by alphas
      Do i = 1, 6
         H(i) = s * Error(i) 
      End do

c ... find the SPD matrix H such that (H e_i, e_i) = alpha(i), 
      Call updateH(xy1, xy2, xy3, xy4, H)

      Return
      End


