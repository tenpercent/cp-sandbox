c ======================================================================
      Subroutine applyGRAD(iGauss, XYL, PSI, FEMtype, nfa, dim, U)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), U(9, MaxSize, *)

c Data for the reference tetrahedron
      Integer IPF(4, 4)

      Real*8  XYL(4, *)
      Real*8  GRAD_P1(3, 4), GRAD_P2(3, 10, MaxPointGauss)
      Real*8                 GRAD_P3(3, 20, MaxPointGauss)

      Integer i, j, k, n, l, mfa, iP1,iP2,iP3
      Real*8  x, y, z, Lfun, s1, s2, s3

      DATA    IPF/1,2,3,4, 2,3,4,1, 1,3,4,2, 1,2,4,3/ 

      DATA    GRAD_P1/-1,-1,-1, 1,0,0, 0,1,0, 0,0,1/
c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 0D0
         End do

      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 4
         dim = 3 
         Do i = 1, nfa
            Do k = 1, dim
               U(k, i, 1) = 0D0
               Do j = 1, 3
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
         nfa = 4
         dim = 3 
         Do i = 1, nfa
            l = IPF(4, i)
            Do k = 1, dim
               U(k, i, 1) = 0D0
               Do j = 1, 3
                  U(k, i, 1) = U(k, i, 1) - 3*PSI(j, k) * GRAD_P1(j, l)
               End do
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         nfa = 12
         dim = 9
         Call clearU(1, nfa, dim, U)

         Do i = 1, 4
            l = IPF(4, i)
            Do k = 1, 3
               U(k, l, 1) = 0D0
               Do j = 1, 2
                  U(k, l, 1) = U(k, l, 1) - 3*PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 3, l + 4, 1) = U(k, l, 1)
            End do
         End do

         Call copyGauss(iGauss, nfa, dim, U)

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR.
     &        FEMtype.EQ.FEM_P2vector) Then
         nfa = 10
         dim = 3

         If(FEMtype.EQ.FEM_P2vector) Call clearU(iGauss, 30, 9, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            z = XYL(4, n)

            Do i = 1, 4
               Do k = 1, dim
                  GRAD_P2(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, y, z) - 1D0) 
               End do
            End do

            mfa = 4
            Do i = 1, 4
               Do j = i + 1, 4
                  mfa = mfa + 1
                  Do k = 1, dim
                     GRAD_P2(k, mfa, n) = 
     &                    4 * (Lfun(i, x, y, z) * GRAD_P1(k, j) +
     &                         Lfun(j, x, y, z) * GRAD_P1(k, i))
                  End do
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 3
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, k) * GRAD_P2(j, i, n)
                  End do
               End do
            End do
         End do
         
         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 30
            dim = 9

            Do n = 1, iGauss
               Do i = 1, 10
                  Do k = 1, 3
                     U(k + 3, i + 10, n) = U(k, i, n)
                     U(k + 6, i + 20, n) = U(k, i, n)
                  End do
               End do
            End do
         End if

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 12
         dim = 9 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 4
            Do k = 1, 3
               Do j = 1, 3
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 3, i + 4, 1) = U(k, i, 1)
               U(k + 6, i + 8, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)
   
c ... next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 20
         dim = 3

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            z = XYL(4, n)

c  ...  basis functions at vertices
            Do i = 1, 4
               s1 = Lfun(i, x, y, z)

               Do k = 1, dim
                  GRAD_P3(k, i, n) = 
     &                GRAD_P1(k, i) * (13.5D0 * s1 * (s1 - 1D0) + 1D0)
               End do
            End do

c  ...  basis functions on edges
            mfa = 4
            Do i = 1, 3
               Do j = i + 1, 4
                  s1 = Lfun(i, x, y, z)
                  s2 = Lfun(j, x, y, z)

                  mfa = mfa + 1
                  Do k = 1, dim
                     GRAD_P3(k, mfa, n) = 
     &                   GRAD_P1(k, i) * (6D0 * s1 - 1D0) * s2 +
     &                   GRAD_P1(k, j) * (3D0 * s1 - 1D0) * s1

                     GRAD_P3(k, mfa + 6, n) = 
     &                   GRAD_P1(k, j) * (6D0 * s2 - 1D0) * s1 +
     &                   GRAD_P1(k, i) * (3D0 * s2 - 1D0) * s2
                  End do
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

               Do k = 1, dim
                  GRAD_P3(k, mfa, n) = 27 * 
     &               (GRAD_P1(k, iP1) * s2 * s3 +
     &                GRAD_P1(k, iP2) * s1 * s3 +
     &                GRAD_P1(k, iP3) * s1 * s2)
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 3
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, k) * GRAD_P3(j, i, n)
                  End do
               End do
            End do
         End do

      Else
         nfa = 0
         dim = -FEMtype
      End if
      Return
      End



c ======================================================================
      Subroutine applyGRA2(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, det)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), XYL(4, *), U(9, MaxSize, *)
      Real*8  XYP(3, 4), det

c Local variables 
      Integer i, j, k, l, iP1,iP2,iP3,iP4, nfc
      Real*8  V(9, MaxSize, MaxPointGauss)
      Real*8  vol, sqr, calSqr
c ======================================================================
c Data for the reference tetrahedron
      Integer IPF(4, 4)

      DATA    IPF/1,2,3,4, 2,3,4,1, 1,3,4,2, 1,2,4,3/ 
c ======================================================================
      If(FEMtype.EQ.FEM_RT0) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V)

         nfa = 4
         dim = 9
         vol = det / 2

         Call clearU(1, nfa, dim, U)

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            iP3 = IPF(3, i)
            iP4 = IPF(4, i)
            
            sqr = calSqr(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3)) 
            sqr = sqr / vol
           
            l = 0
            Do j = 1, 3
               Do k = 1, 3
                 l = l + 1
                 U(l, i, 1) = 
     &                (V(k, iP1, 1) * (XYP(j, iP1) - XYP(j, iP4))
     &               + V(k, iP2, 1) * (XYP(j, iP2) - XYP(j, iP4))
     &               + V(k, iP3, 1) * (XYP(j, iP3) - XYP(j, iP4)))
               End do
            End do

            vol = -vol
         End do

         Call copyGauss(iGauss, nfa, dim, U)
      Else
         nfa = 0
         dim = -FEMtype
      End if
      Return
      End




