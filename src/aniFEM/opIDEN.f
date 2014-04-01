c ======================================================================
      Subroutine applyIDEN(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, det)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), XYL(4, *), U(9, MaxSize, *)
      Real*8  XYP(3, *), det

c Local variables 
      Integer i, j, k, n, l, iP1,iP2,iP3,iP4, mfa
      Real*8  vol, sqr, calSqr, calEdge, xyn(3), xym(3)
      Real*8  x, y, z, Lfun, s, s1, s2, s3
c ======================================================================
c Data for the reference tetrahedron
      Integer IPF(4, 4), IPR(4, 6)

      DATA    IPF/1,2,3,4, 2,3,4,1, 1,3,4,2, 1,2,4,3/ 
      DATA    IPR/1,2,3,4, 1,3,4,2, 1,4,2,3, 2,3,1,4, 2,4,3,1, 3,4,1,2/ 
c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 1D0
         End do

      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 4
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = XYL(i, n)
            End do
         End do

      Else If(FEMtype.EQ.FEM_CR1) Then
         nfa = 4
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               l = IPF(4, i)
               U(1, i, n) = 1 - 3 * XYL(l, n)
            End do
         End do

      Else If(FEMtype.EQ.FEM_P2 .OR.
     &        FEMtype.EQ.FEM_P2vector) Then
         nfa = 10
         dim = 1

         If(FEMtype.EQ.FEM_P2vector) Call clearU(iGauss, 30, 3, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            z = XYL(4, n)
            Do i = 1, 4
               s = Lfun(i, x, y, z)
               U(1, i, n) = s * (2 * s - 1D0)
            End do

            mfa = 4
            Do i = 1, 4
               Do j = i + 1, 4
                  mfa = mfa + 1
                  U(1, mfa, n) = 4 * Lfun(i, x, y, z) * Lfun(j, x, y, z)
               End do
            End do
         End do

         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 30
            dim = 3

            Do n = 1, iGauss
               Do i = 1, 10
                  U(2, i + 10, n) = U(1, i, n) 
                  U(3, i + 20, n) = U(1, i, n) 
               End do
            End do
         End if

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 12
         dim = 3

         Call clearU(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            Do i = 1, 4
               U(1, i,     n) = XYL(i, n)
               U(2, i + 4, n) = XYL(i, n)
               U(3, i + 8, n) = XYL(i, n)
            End do
         End do

      Else If(FEMtype.EQ.FEM_P0vector) Then
         nfa = 3
         dim = 3

         Call clearU(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            U(1, 1, n) = 1D0
            U(2, 2, n) = 1D0
            U(3, 3, n) = 1D0
         End do

      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 20
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            z = XYL(4, n)

c  ...  basis functions at vertices
            Do i = 1, 4
               s1 = Lfun(i, x, y, z)

               U(1, i, n) = s1 * (3*s1 - 1D0) * (3*s1 - 2D0) / 2
            End do

c  ...  basis functions on edges
            mfa = 4
            Do i = 1, 3
               Do j = i + 1, 4
                  s1 = Lfun(i, x, y, z)
                  s2 = Lfun(j, x, y, z)

                  mfa = mfa + 1
                  U(1, mfa,     n) = s1 * (3*s1 - 1D0) * s2 * 4.5D0
                  U(1, mfa + 6, n) = s1 * (3*s2 - 1D0) * s2 * 4.5D0
               End do
            End do

c  ...  basis functions on edges
            mfa = 16
            Do i = 1, 4
               mfa = mfa + 1

               iP1 = IPF(1, i)
               iP2 = IPF(2, i)
               iP3 = IPF(3, i)

               s1 = Lfun(iP1, x, y, z)
               s2 = Lfun(iP2, x, y, z)
               s3 = Lfun(iP3, x, y, z)

               U(1, mfa, n) = 27 * s1 * s2 * s3
            End do
         End do

      Else If(FEMtype.EQ.FEM_RT0) Then
         nfa = 4
         dim = 3
         vol = dabs(det) / 2

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            iP3 = IPF(3, i)
            iP4 = IPF(4, i)
            
            sqr = calSqr(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3)) 
            sqr = sqr / vol
            Do n = 1, iGauss
               Do k = 1, dim
                  U(k, i, n) = sqr *
     &                 (XYL(iP1, n) * (XYP(k, iP1) - XYP(k, iP4))
     &                + XYL(iP2, n) * (XYP(k, iP2) - XYP(k, iP4))
     &                + XYL(iP3, n) * (XYP(k, iP3) - XYP(k, iP4)))
               End do
            End do
         End do

      Else If(FEMtype.EQ.FEM_ND0) Then
         nfa = 6
         dim = 3
         vol = det

         Do i = 1, nfa
            iP1 = IPR(1, i)
            iP2 = IPR(2, i)
            iP3 = IPR(3, i)
            iP4 = IPR(4, i)

            Call calNormalVec(XYP(1,iP1), XYP(1,iP3), XYP(1,iP4), xyn)
            Call calNormalVec(XYP(1,iP2), XYP(1,iP3), XYP(1,iP4), xym)

            sqr = calEdge(XYP(1, iP1), XYP(1, iP2))
            sqr = sqr / vol
            Do n = 1, iGauss
               Do k = 1, dim
                  U(k, i, n) = sqr * (XYL(iP1, n) * xyn(k) 
     &                              + XYL(iP2, n) * xym(k))
               End do
            End do
         End do

      Else
         nfa = 0
         dim = -FEMtype
      End if
      Return
      End

