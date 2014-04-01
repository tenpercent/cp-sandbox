C ======================================================================
      Subroutine Lp_norm(nP, Lp, Metric)
C ======================================================================
C  Routine computes the metric for L_p norm using the metric
C  generated for the maximum norm.
C
C     Lp    - norm in which the metric is optimal:
C             Lp > 0  means  L_p     norm
C             Lp = 0  means  maximum norm (L_infinity)
C             Lp < 0  means  H_1     norm (not implemented yet)
C ======================================================================
      Real*8  Metric(6, *), Lp, det
C ======================================================================
      If(Lp.EQ.0d0) Return
      If(Lp.LT.0d0) Stop 'Lp_gradnorm: negative Lp'

      Do n = 1, nP
         det = Metric(1, n) * (Metric(2, n) * Metric(3, n) 
     &                       - Metric(5, n) ** 2) -
     &         Metric(4, n) * (Metric(4, n) * Metric(3, n) 
     &                       - Metric(5, n) * Metric(6, n)) +
     &         Metric(6, n) * (Metric(4, n) * Metric(5, n) 
     &                      - Metric(6, n) * Metric(2, n))

         det = det ** (-1D0 / (2 * Lp + 2))

         Do i = 1, 6
            Metric(i, n) = Metric(i, n) * det
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine Lp_gradnorm(nP, Lp, Metric)
C ======================================================================
C  Routine computes the metric for L_p norm of the gradient
C  using the metric generated for the maximum norm.
C
C     Lp    - norm for which the metric is to be adjusted:
C             Lp > 0  means  L_p     norm
C             Lp = 0  means  maximum norm (L_infinity)
C ======================================================================
      Real*8  Metric(6, *), Lp, det
C ======================================================================
      If(Lp.EQ.0d0) Return
      If(Lp.LT.0d0) Stop 'Lp_gradnorm: negative Lp'

      Do n = 1, nP
         det = Metric(1, n) * (Metric(2, n) * Metric(3, n) 
     &                       - Metric(5, n) ** 2) -
     &         Metric(4, n) * (Metric(4, n) * Metric(3, n) 
     &                       - Metric(5, n) * Metric(6, n)) +
     &         Metric(6, n) * (Metric(4, n) * Metric(5, n) 
     &                       - Metric(6, n) * Metric(2, n))

         det = det ** (-1D0 / (2 + Lp))

         Do i = 1, 6
            Metric(i, n) = Metric(i, n) * det
         End do
      End do

      Return
      End

