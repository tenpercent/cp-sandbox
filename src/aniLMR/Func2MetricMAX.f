C=======================================================================
      Subroutine Func2MetricMAX(Func,  
     &                          nP, nE, XYP, IPE, 
     &                          Metric, MaxWr, rW)
C=======================================================================
C  Input:
C     Func       - Real*8 Function f(xy), where xy(3)
C
C     nP         - the number of nodes
C     nE         - the number of tetrahedra
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table for tets
C
C  Output: 
C     Metric(6, nE) - nodal metric 
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
      Real*8   Func, Metric(6, *)
      EXTERNAL Func

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, k, n, i1, i2, iP1, iP2, iP3, iP4
      Real*8   det, Up1, Up2, Ur1, xyt(3), H(6), Error(6)

C=======================================================================
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'Func2MetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do

      Do n = 1, nE
         k = 0
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               Up1 = Func(XYP(1, iP1))
               Up2 = Func(XYP(1, iP2))

               Do i = 1, 3
                  xyt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
               End do

               Ur1 = Func(xyt)

               k = k + 1
               Error(k) = 4 * (Ur1 - (Up1 + Up2) / 2)
            End do
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call TetMetric(XYP(1,iP1), XYP(1,iP2), 
     &                  XYP(1,iP3), XYP(1,iP4), Error, H)

         det = H(1) * (H(2) * H(3) - H(5) ** 2) 
     &       - H(4) * (H(4) * H(3) - H(5) * H(6)) 
     &       + H(6) * (H(4) * H(5) - H(6) * H(2))

         Do i = 1, 4
            iP1 = IPE(i, n)

c  ...  update the maximum metric
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
      Subroutine Func2MetricCell(Func,  
     &                           nP, nE, XYP, IPE, 
     &                           Metric, MaxWr, rW)
C=======================================================================
C  Input:
C     Func       - Real*8 Function f(xy), where xy(3)
C
C     nP         - the number of nodes
C     nE         - the number of tetrahedra
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table for tets
C
C  Output: 
C     Metric(6, nE) - cell metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr
C=======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(3, *)
      Integer  IPE(4, *)

C Function and metric 
      Real*8   Func, Metric(6, *)
      EXTERNAL Func

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, k, n, i1, i2, iP1, iP2, iP3, iP4
      Real*8   Up1, Up2, Ur1, xyt(3), H(6), Error(6)

C=======================================================================
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'Func2MetricCell', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do

      Do n = 1, nE
         k = 0
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
 
               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               Up1 = Func(XYP(1, iP1))
               Up2 = Func(XYP(1, iP2))

               Do i = 1, 3
                  xyt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
               End do

               Ur1 = Func(xyt)

               k = k + 1
               Error(k) = 4 * (Ur1 - (Up1 + Up2) / 2)
            End do
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call TetMetric(XYP(1,iP1), XYP(1,iP2), 
     &                  XYP(1,iP3), XYP(1,iP4), Error, H)

         Do i = 1, 6
            Metric(i, n) = H(i)
         End do
      End do

      Return
      End


