C ======================================================================
      Subroutine EdgeEst2MetricMAX(Error,
     &                             nP, nE, XYP, IPE, 
     &                             Metric, MaxWr, rW)
C ======================================================================
C  Routine uses method of shifts to calculate the nodal metric 
C  corresponding to edge-based error estimates in array Error.
C ======================================================================
C  Input:
C     Error(6, nE) - error estimates on edges
C
C     nP - the number of nodes
C     nE - the number of triangles
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table
C
C  Output: 
C     Metric(6, nP) - nodal metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr, MaxWr > nP
C ======================================================================
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
      Integer  i, j, n, iP1, iP2, iP3, iP4
      Real*8   det, H(6)

C ======================================================================
c ... memory check
      If(nP.GT.MaxWr) Call errMesLMR(1002, 'EdgeEst2MetricMAX', 
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

         Call TetMetric(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3), XYP(1,iP4),
     &                  Error(1, n), H)

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


C ======================================================================
      Subroutine EdgeEst2MetricCell(Error,  
     &                              nP, nE, XYP, IPE, 
     &                              Metric, MaxWr, rW)
C ======================================================================
C  Input:
C     nP - the number of nodes
C     nE - the number of triangles
C
C     Error(6, nE) - given function on edges
C
C     XYP(3, nP) - Cartesian coordinates of the nodes
C     IPE(4, nE) - the connecticity table
C
C  Output: 
C     Metric(6, nE) - cell metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr
C ======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(3, *)
      Integer  IPE(4, *)

C Function and metric 
      Real*8   Error(6, *)
      Real*8   Metric(6, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  n, iP1, iP2, iP3, iP4

C ======================================================================
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'EdgeEst2MetricCell', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         iP4 = IPE(4, n)

         Call TetMetric(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3), XYP(1,iP4),
     &                  Error(1, n), Metric(1, n))
      End do

      Return
      End



C ======================================================================
      Subroutine TetMetric(xy1, xy2, xy3, xy4, Error, H)
C ======================================================================
C Routine computes the metric H for a tet given by its vertices.
C ======================================================================
      implicit none

      include 'magic.fd'
C ======================================================================
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3)
      Real*8   Error(6), H(6)

c Local variables
      Integer  i, j
      Real*8   B(6, 6), s, p, t

C ======================================================================
c ... calculate matrix B = b_i . b_j where b_i is the bubble function 
      s = 1D0 / 1260
      t = 1D0 / 2520

      Do i = 1, 6
         B(i, i) = s
         Do j = i + 1, 6
            B(i, j) = t 
            B(j, i) = t 
         End do
      End do

c ... define edge errors
      s = 0D0
      p = 0D0
      Do i = 1, 6
          Do j = i, 6
            s = s + B(i, j) * Error(i) * Error(j)
          End do
         p = p + Error(i)
      End do

      s = s / p

      Do i = 1, 6
         H(i) = Error(i) * s
      End do

c ... find the SPD matrix H such that (H e_i, e_i) = alpha(i), 
      Call updateH(xy1, xy2, xy3, xy4, H)

      Return
      End



C=======================================================================
      Subroutine updateH(xy1, xy2, xy3, xy4, H)
C=======================================================================
C  Routine finds the SPD matrix H such that (M e_i, e_i) = H(i)
C  and conditions of Lemma 1 are satisfied. The result is returned in H.
C
C  If H is close to singular, we modify it as described in Lemma 1. This
C  must be verified. If H is indefinite, we take its spectral module. 
C=======================================================================
      implicit none

      include 'magic.fd'
      include 'fem3Dtet.fd'

      Real*8   xy1(*), xy2(*), xy3(*), xy4(*), H(6)

      Integer  i,k, iLoop, info, IPIV(6)
      Real*8   ax(3, 6), B(6, 6), H0(6), s, det

C=======================================================================
      Do i = 1, 3
         ax(i, 1) = xy1(i) - xy2(i)
         ax(i, 2) = xy1(i) - xy3(i)
         ax(i, 3) = xy1(i) - xy4(i)
         ax(i, 4) = xy2(i) - xy3(i)
         ax(i, 5) = xy2(i) - xy4(i)
         ax(i, 6) = xy3(i) - xy4(i)
      End do

      Do i = 1, 6
         H0(i) = H(i)
      End do

      Do iLoop = 1, 2
         Do i = 1, 6
            B(i, 1) = ax(1, i) ** 2
            B(i, 2) = ax(2, i) ** 2
            B(i, 3) = ax(3, i) ** 2
            B(i, 4) = 2 * ax(1, i) * ax(2, i)
            B(i, 5) = 2 * ax(2, i) * ax(3, i)
            B(i, 6) = 2 * ax(1, i) * ax(3, i)

            H(i) = H0(i)
         End do

         Call dgesv(6, 1, B, 6, IPIV, H, 6, info)
         If(info.NE.0) Call errMesLMR(3011, 'updateH',
     &                      'Error in the LAPACK routine dgesv')

         s   = (H(1) + H(2) + H(3)) ** 3

         det = H(1) * (H(2) * H(3) - H(5) ** 2) 
     &       - H(4) * (H(4) * H(3) - H(5) * H(6)) 
     &       + H(6) * (H(4) * H(5) - H(6) * H(2))

         If(dabs(det).LT.1D-8*dabs(s)) Then
            If(iLoop.EQ.2) Call errMesLMR(6007, 'updateH',
     &                          'System error: H must be nonsingular')

            k = 1
            Do i = 2, 6
               If(H0(i).GT.H0(k)) k = i
            End do
 
            If(H0(k).GT.0D0) Then
               H0(k) = H0(k) * 1.001
            Else
               H(1) = AniEigenvalue
               H(2) = AniEigenvalue
               H(3) = AniEigenvalue
               H(4) = 0d0
               H(5) = 0d0
               H(6) = 0d0

               Goto 100
            End if
         Else
            Call SpectralModule(H, det)
            Goto 100
         End if
      End do

 100  Continue

      Return
      End

