C ======================================================================
      Integer Function ANI_Metric_Eucl(x, y, z, Metric)
C ======================================================================
      include 'metric.fd'
C ======================================================================
C  This routine creates a metric at point (x,y,z). The
C  metric is a 3x3 positive definite symmetric tensor:
C
C                M11   M12   M13
C      Metric =  M12   M22   M23
C                M13   M23   M33
C
C  Only the upper triangular part of Metric must be defined.
C ======================================================================
      Real*8  x, y, z, Metric(3, 3)

      Metric(1,1) = 1D0
      Metric(2,2) = 1D0
      Metric(3,3) = 1D0

      Metric(1,2) = 0D0
      Metric(1,3) = 0D0
      Metric(2,3) = 0D0

      ANI_Metric_Eucl = ANI_METRIC_TENSOR

      Return
      End

