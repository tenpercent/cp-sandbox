C ================================================================
      Subroutine makQ(
C ================================================================
     &     nLoop,
c group (M)
     &     nP, nE, XYP, IPE, nEv, IEV,
     &     nEStar, hStar,
c group (Q)
     &     HesP, detG, qE)
C ================================================================
      include 'magic.fd'
      include 'cubature.fd'
C ================================================================
C Quality computation for mesh elements.
C
C Pre-conditions:  1. connectivity structure {IPE(4, *), XYP(3, *)}
C                  2. tensor metric field HesP(6, *)
C
C Post-conditions: 1. determinant of the tensor metric, detG(*)
C                  2. quality of mesh elements, qE(*) 
C
C Remark: nLoop is used in parallel codes.
C         It should be equal to 1 in serial nodal/analytic codes.
C         It should be equal to 0 in serial fixquality code.
C ================================================================
C group (M)
      Real*8  XYP(3, *)
      Integer IPE(4, *), IEV(*)

      Real*8  hStar

C group (Q)
      Real*8  HesP(6, *)
      Real*8  detG(*), qE(*)

C group (Local variables)
      Real*8  d1, d2, d3, d4, dsum, calVol, HesAvg(6), detAvg
      Real*8  VStar, Vn

C ================================================================
      Do n = 1, nP
         Call calDet(HesP(1, n), detG(n))
      End do

      If(nLoop.EQ.1) Then
         VStar = 0D0
         Do n = 1, nE
            i1 = IPE(1, n)
            i2 = IPE(2, n)
            i3 = IPE(3, n)
            i4 = IPE(4, n)

c  ...  1-point quadrature
            Do i = 1, 6
               HesAvg(i) = (HesP(i, i1) + HesP(i, i2) +
     &                      HesP(i, i3) + HesP(i, i4)) / 4
            End do
            Call calDet(HesAvg, detAvg)
 
            Vn = calVol(XYP(1, i1), XYP(1, i2), XYP(1, i3), XYP(1, i4))
            VStar = VStar + dabs(Vn) * dsqrt(dabs(detAvg))

c  ...  4-point quadrature
c           d1 = detG(i1)
c           d2 = detG(i2)
c           d3 = detG(i3)
c           d4 = detG(i4)
c          
c           dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) 
c    &           + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) 
c    &           + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) 
c    &           + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) 

c           Vn = calVol(XYP(1, i1), XYP(1, i2), XYP(1, i3), XYP(1, i4))
c           VStar = VStar + dabs(Vn) * dsum / 4
         End do
         hStar = (VStar / nEStar * magicNumber *
     &                    12D0 / dsqrt(2D0)) ** 3.333D-1
      End if

      Do n = 1, nE
         i1 = IPE(1, n)
         i2 = IPE(2, n)
         i3 = IPE(3, n)
         i4 = IPE(4, n)

         Call calQE(
     &        HesP(1, i1), XYP(1, i1),
     &        HesP(1, i2), XYP(1, i2),
     &        HesP(1, i3), XYP(1, i3),
     &        HesP(1, i4), XYP(1, i4),
     &        hStar, qE(n), Vn)
      End do


c ... set quality of fixed elements to 1
      Do n = 1, nEv
         qE(IEV(n)) = 1D0
      End do 

      Return
      End



C ================================================================
      Subroutine updQa(n, XYP, IPE, IEE, qE)
C ================================================================
C Initial quality modification for tangled elements and 
C their closest (face-) neighboors.
C ================================================================
      Real*8  XYP(3, *), qE(*)
      Integer IPE(5, *), IEE(4, *)

C (Local variables)
      Integer ip(5)
      Real*8  calVol, v1, v2
      Logical check33

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         iE = IEE(i1, n)
         If(iE.LE.0) Goto 100

         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPE(i1, n)
         iP2 = IPE(i2, n)
         iP3 = IPE(i3, n)

         Do j1 = 1, 4
            j2 = ip(j1 + 1)
            j3 = ip(j2 + 1)

            jP1 = IPE(j1, iE)
            jP2 = IPE(j2, iE)
            jP3 = IPE(j3, iE)

            If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
               i4  = ip(i3 + 1)
               iP4 = IPE(i4, n)

               j4  = ip(j3 + 1)
               jP4 = IPE(j4, iE)

               v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                     XYP(1, iP3), XYP(1, iP4))
               v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                     XYP(1, iP3), XYP(1, jP4))

               If(v1 * v2.GE.0D0) Then
                  qE(n)  = -dabs(qE(n))
                  qE(iE) = -dabs(qE(iE))
               End if
               Goto 100
            End if
         End do
 100  Continue

      Return
      End


C ================================================================
      Subroutine updQb(nEs, lE, iEs, XYP, IPEs, qEs)
C ================================================================
C Dynamic quality modification for tangled elements inside
C a super-element.
C
C Remark: non-efficient, time-consuming, but robust algorithm.
C ================================================================
      Real*8  XYP(3, *), qEs(*)
      Integer iEs(*), IPEs(5, *)

C group (Local variables)
      Integer ip(5)
      Real*8  calVol, v1, v2
      Logical check33

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)
         iP3 = IPEs(i3, nEs)

         Do 20 k = 1, lE
            If(iEs(k).LT.0)  Goto 20
            If(k.EQ.nEs)     Goto 20

            Do j1 = 1, 4
               j2 = ip(j1 + 1)
               j3 = ip(j2 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)
               jP3 = IPEs(j3, k)

               If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  i4  = ip(i3 + 1)
                  iP4 = IPEs(i4, nEs)

                  j4  = ip(j3 + 1)
                  jP4 = IPEs(j4, k)

                  v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                        XYP(1, iP3), XYP(1, iP4))
                  v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                        XYP(1, iP3), XYP(1, jP4))

                  If(v1 * v2.GE.0D0) Then
                     qEs(nEs) = -dabs(qEs(nEs))
                     qEs(k)   = -dabs(qEs(k))
                  End if

                  Goto 100
               End if
            End do
 20      Continue
 100  Continue

      Return
      End


C ================================================================
      Real*8 Function calVol(xy1, xy2, xy3, xy4)
C ================================================================
C The oriented volume of the tetrahedron given by 4 vertices
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3)

C group (Local variables)
      Real*8 v1(3), v2(3), v3(3)

      Do i = 1, 3
         v1(i) = xy1(i) - xy4(i)
         v2(i) = xy2(i) - xy4(i)
         v3(i) = xy3(i) - xy4(i)
      End do

      calVol = (v1(1) * (v2(2) * v3(3) - v2(3) * v3(2)) +
     &          v1(2) * (v2(3) * v3(1) - v2(1) * v3(3)) +
     &          v1(3) * (v2(1) * v3(2) - v2(2) * v3(1))) /
     &          6D0
      Return
      End


C ================================================================
      Real*8 Function mutualOrientation(xy1, xy2, xy3, xy4, xy5)
C ================================================================
C Mutual orientation of points p4 an p5 with respect to plane p1-p2-p3.
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3), xy5(3)

C group (Local variables)
      Real*8 v2(3), v3(3), v4(3), v5(3), ax, ay, az, vol1, vol2

      Do i = 1, 3
         v2(i) = xy2(i) - xy1(i)
         v3(i) = xy3(i) - xy1(i)
         v4(i) = xy4(i) - xy1(i)
         v5(i) = xy5(i) - xy1(i)
      End do

      ax = v2(2) * v3(3) - v2(3) * v3(2)
      ay = v2(3) * v3(1) - v2(1) * v3(3)
      az = v2(1) * v3(2) - v2(2) * v3(1)

      vol1 = ax * v4(1) + ay * v4(2) + az * v4(3)
      vol2 = ax * v5(1) + ay * v5(2) + az * v5(3)

      mutualOrientation = vol1 * vol2

      Return
      End


C ================================================================
      Real*8 Function calSqr(xy1, xy2, xy3)
C ================================================================
C  Routine calculates area of a triangle formed by three points
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3)

C Local variables
      Real*8 ax, ay, az, bx, by, bz

      ax = xy1(1) - xy3(1)
      ay = xy1(2) - xy3(2)
      az = xy1(3) - xy3(3)

      bx = xy2(1) - xy3(1)
      by = xy2(2) - xy3(2)
      bz = xy2(3) - xy3(3)

      calSqr = 5D-1 * dsqrt((ay * bz - az * by) ** 2 +
     &                      (ax * bz - az * bx) ** 2 +
     &                      (ax * by - ay * bx) ** 2)
      Return
      End


C ======================================================================
      Subroutine calDet(HesP, detG)
C ======================================================================
C Routine computes the determinant of the metric HesP.
C
C Remark: robustness of the overall code was increased by
C         replacing errMes() with the computation of |H|.
C ======================================================================
      Real*8 HesP(6), detG

C ======================================================================
      detG = HesP(1) * (HesP(2) * HesP(3) - HesP(5) ** 2) -
     &       HesP(4) * (HesP(4) * HesP(3) - HesP(5) * HesP(6)) +
     &       HesP(6) * (HesP(4) * HesP(5) - HesP(6) * HesP(2))

      If(detG.LE.0D0) Then
         Call SpectralModule(HesP, detG)
      End if
      
      Return
      End



C ======================================================================
      Subroutine SpectralModule(HesP, detG)
C ======================================================================
      include 'magic.fd'
C ======================================================================
C Routines sets the minimal eigenvalue of HesP to some constant
C ======================================================================
      Real*8  HesP(6), detG
    
C Local arrays for LAPACK
      Real*8  A(3, 3), E(3), rW(15)
      Integer info

C ======================================================================
        A(1, 1) = HesP(1)
        A(2, 2) = HesP(2)
        A(3, 3) = HesP(3)

        A(1, 2) = HesP(4)
        A(2, 3) = HesP(5)
        A(1, 3) = HesP(6)

        Call dsyev('V', 'U', 3, A, 3, E, rW, 15, info)
        If(info.NE.0) Call errMes(3011, 'calDet',
     &                    'Error in LAPACK routine dsyev')

        E(1) = dabs(E(1))
        E(2) = dabs(E(2))
        E(3) = dabs(E(3))

        E(1) = max( E(1), E(3) * AniRatio )
      If(E(1).EQ.0D0) Then
         Do i = 1, 3
            E(i) = AniEigenvalue
         End do
      End if

        Do i = 1, 6
           HesP(i) = 0D0
        End do

        Do i = 1, 3
          HesP(1) = HesP(1) + E(i) * A(1, i) ** 2 
          HesP(2) = HesP(2) + E(i) * A(2, i) ** 2
          HesP(3) = HesP(3) + E(i) * A(3, i) ** 2

          HesP(4) = HesP(4) + E(i) * A(1, i) * A(2, i) 
          HesP(5) = HesP(5) + E(i) * A(2, i) * A(3, i) 
          HesP(6) = HesP(6) + E(i) * A(1, i) * A(3, i) 
        End do

      detG = E(1) * E(2) * E(3)

      Return
      End



C ================================================================
      Subroutine calQE(
     &           Hes1, xy1, Hes2, xy2, Hes3, xy3, Hes4, xy4,
     &      hStar, qE, volume)
C ================================================================
C Computing quality of tetrahedron {xy1, ..., xy4} assuming that
C the metric field is linear.
C ================================================================
      include 'cubature.fd'
C ================================================================
      Real*8 Hes1(6), xy1(3)
      Real*8 Hes2(6), xy2(3)
      Real*8 Hes3(6), xy3(3)
      Real*8 Hes4(6), xy4(3)
      Real*8 hStar, qE, volume
C ================================================================
c  local variables
      Real*8 HesAvg(6), d1, d2, d3, d4, dsum, detAvg
      Real*8 calVol, Pk, Vk
      Real*8 F, x1, y1, z1

C ================================================================
      Pk = 0D0
      Do n = 1, 6
         If(n.EQ.1) Then
            x1 = xy1(1) - xy4(1)
            y1 = xy1(2) - xy4(2)
            z1 = xy1(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.2) Then
            x1 = xy2(1) - xy4(1)
            y1 = xy2(2) - xy4(2)
            z1 = xy2(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.3) Then
            x1 = xy3(1) - xy4(1)
            y1 = xy3(2) - xy4(2)
            z1 = xy3(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes3(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.4) Then
            x1 = xy1(1) - xy3(1)
            y1 = xy1(2) - xy3(2)
            z1 = xy1(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.5) Then
            x1 = xy2(1) - xy3(1)
            y1 = xy2(2) - xy3(2)
            z1 = xy2(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.6) Then
            x1 = xy1(1) - xy2(1)
            y1 = xy1(2) - xy2(2)
            z1 = xy1(3) - xy2(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes2(i)) / 2
            End do
         End if
         Pk = Pk + dsqrt(HesAvg(1) * x1 * x1 +
     &                   HesAvg(2) * y1 * y1 +
     &                   HesAvg(3) * z1 * z1 +
     &               2 * HesAvg(4) * x1 * y1 +
     &               2 * HesAvg(5) * y1 * z1 +
     &               2 * HesAvg(6) * x1 * z1)
      End do

c ... 1-point quadrature rule
      Do i = 1, 6
         HesAvg(i) = (Hes1(i) + Hes2(i) + Hes3(i) + Hes4(i)) / 4
      End do
      Call calDet(HesAvg, detAvg)

      volume = calVol(xy1, xy2, xy3, xy4)
      Vk = volume * dsqrt(dabs(detAvg)) 

c ... 4-point quadrature rule
c     d1 = det1
c     d2 = det2 
c     d3 = det3 
c     d4 = det4 
c          
c     dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) 
c    &     + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) 
c    &     + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) 
c    &     + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) 

c     volume = calVol(xy1, xy2, xy3, xy4)
c     Vk = volume * dsum / 4

          qE = 1832.8208D0 * dabs(Vk) / (Pk ** 3) 

      If(hStar.GT.0D0) Then
          x1 = Pk / (6 * hStar)
          x1 = min(x1, 1D0 / x1)

         F  = (x1 * (2D0 - x1)) ** 5
         qE = qE * F
      End if

      Return
      End



C ================================================================
      Subroutine calQF(
     &      Hes1, xy1, Hes2, xy2, Hes3, xy3, Hes4, xy4,
     &      hStar, iF, iR)
C ================================================================
      Real*8  Hes1(6), xy1(3)
      Real*8  Hes2(6), xy2(3)
      Real*8  Hes3(6), xy3(3)
      Real*8  Hes4(6), xy4(3)
      Real*8  hStar

      Integer iF(4), iR(6)

C group (Local variables)
      Real*8  HesAvg(6)
      Real*8  qF(4), qR(6)
      Real*8  F, x1, y1, z1

C group (Function)
      F(x1) = (x1 * (2D0 - x1)) ** 5

C ================================================================
      Do n = 1, 6
         If(n.EQ.1) Then
            x1 = xy1(1) - xy2(1)
            y1 = xy1(2) - xy2(2)
            z1 = xy1(3) - xy2(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes2(i)) / 2
            End do
         Else If(n.EQ.2) Then
            x1 = xy1(1) - xy3(1)
            y1 = xy1(2) - xy3(2)
            z1 = xy1(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.3) Then
            x1 = xy1(1) - xy4(1)
            y1 = xy1(2) - xy4(2)
            z1 = xy1(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.4) Then
            x1 = xy2(1) - xy3(1)
            y1 = xy2(2) - xy3(2)
            z1 = xy2(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.5) Then
            x1 = xy2(1) - xy4(1)
            y1 = xy2(2) - xy4(2)
            z1 = xy2(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.6) Then
            x1 = xy3(1) - xy4(1)
            y1 = xy3(2) - xy4(2)
            z1 = xy3(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes3(i) + Hes4(i)) / 2
            End do
         End if

         x1 = dsqrt(HesAvg(1) * x1 * x1 +
     &              HesAvg(2) * y1 * y1 +
     &              HesAvg(3) * z1 * z1 +
     &          2 * HesAvg(4) * x1 * y1 +
     &          2 * HesAvg(5) * y1 * z1 +
     &          2 * HesAvg(6) * x1 * z1) / hStar
         x1 = min(x1, 1D0 / x1)

         iR(n) = n
         qR(n) = F(x1)
      End do


      qF(1) = min(qR(1), qR(2))
      qF(1) = min(qF(1), qR(4))

      qF(2) = min(qR(4), qR(5))
      qF(2) = min(qF(2), qR(6))

      qF(3) = min(qR(2), qR(3))
      qF(3) = min(qF(3), qR(6))

      qF(4) = min(qR(1), qR(3))
      qF(4) = min(qF(4), qR(5))

      Do i = 1, 4
         iF(i) = i
      End do


      Do i = 1, 3
         iMin = i
         Do j = i + 1, 4
            If(qF(j).LT.qF(iMin)) iMin = j
         End do

         x1 = qF(i)
         qF(i) = qF(iMin)
         qF(iMin) = x1

         k = iF(i)
         iF(i) = iF(iMin)
         iF(iMin) = k
      End do


      Do i = 1, 5
         iMin = i
         Do j = i + 1, 6
            If(qR(j).LT.qR(iMin)) iMin = j
         End do

         x1 = qR(i)
         qR(i) = qR(iMin)
         qR(iMin) = x1

         k = iR(i)
         iR(i) = iR(iMin)
         iR(iMin) = k
      End do

      Return
      End



C ================================================================
      Subroutine HesBnd(lP, iPs, ICP, HesP, HesPs)
C ================================================================
      include 'color.fd'
C ================================================================
      Integer iPs(*), ICP(*)
      Real*8  HesP(6, *), HesPs(6)

C (Local variables)
      Real*8  hesB(6), hesI(6)
      Logical ifXnode

C ================================================================
      nI = 0
      nB = 0

      Do i = 1, 6
         hesI(i) = 0D0
         hesB(i) = 0D0
      End do

      Do n = 1, lP
         iP = iPs(n)
         If(ifXnode(ICP(iP), jInode)) Then
            nI = nI + 1
            Do i = 1, 6
               hesI(i) = hesI(i) + HesP(i, iP)
            End do
         Else
            nB = nB + 1
            Do i = 1, 6
               hesB(i) = hesB(i) + HesP(i, iP)
            End do
         End if
      End do

      If(nI.GT.0) Then
         Do i = 1, 6
            HesPs(i) = hesI(i) / nI
         End do
      Else If(nB.GT.0) Then
         Do i = 1, 6
            HesPs(i) = hesB(i) / nB
         End do
      Else
         Call errMes(6001, 'HesBnd', 'system error')
      End if

      Return
      End



C ======================================================================
      Subroutine iniQ_analytic(nP, XYP, MetricFunction, HesP)
C ======================================================================
C  Three Fortran routines below create a metric field which
C  is 3x3 variable positive definite diagonal tensor HesP,
C
C             M11   M12   M13 
C      HesP = M12   M22   M23 
C             M13   M23   M33
C
C  where Mij = Mij(x, y, z).
C
C  The tensor element are enumerated in the following order:
C  HesP_{11}, HesP_{22}, HesP_{33}, HesP_{12}, HesP_{23}, HesP_{13}
C
C ======================================================================
      include  'metric.fd'

      Real*8   XYP(3, *)
      Real*8   HesP(6, *)

      Integer  MetricFunction
      EXTERNAL MetricFunction

      Real*8   x, y, z, Metric(3, 3)
C =====================================================================
      Do n = 1, nP
         x = XYP(1, n)
         y = XYP(2, n)
         z = XYP(3, n)

         i = MetricFunction(x, y, z, Metric)

         HesP(1, n) = Metric(1,1)
         HesP(2, n) = Metric(2,2)
         HesP(3, n) = Metric(3,3)
         HesP(4, n) = Metric(1,2)
         HesP(5, n) = Metric(2,3)
         HesP(6, n) = Metric(1,3)
      End do

      Return
      End


C ================================================================
      Real*8 Function avgQ(nE, qE, L1E, L2E)
C ================================================================
      Real*8  qE(*)
      Integer L2E(*), L1E(2, *) 

      avgQ = 0D0

      iE = L2E(1)
      Do n = 1, nE
         avgQ = avgQ + qE(iE)
         iE = L1E(2, iE)
      End do

      avgQ = avgQ / nE

      Return
      End

