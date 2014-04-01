c ======================================================================
      Subroutine applyDIV(iGauss, XYL, PSI, FEMtype, 
     &                    nfa, dim, U, XYP, det)
c ======================================================================
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), XYL(4, *), U(9, MaxSize, *)
      Real*8  XYP(3, 4), det
c ======================================================================
c Local variables 
      Real*8  V(9, MaxSize, MaxPointGauss)
      Real*8  vol, sqr, calSqr, calEdge, xyn(3), xym(3), s
c ======================================================================
c Data for the reference tetrahedron
      Integer IPF(4, 4), IPR(4, 6)

      DATA    IPF/1,2,3,4, 2,3,4,1, 1,3,4,2, 1,2,4,3/ 
      DATA    IPR/1,2,3,4, 1,3,4,2, 1,4,2,3, 2,3,1,4, 2,4,3,1, 3,4,1,2/ 
c ======================================================================
      If(FEMtype.EQ.FEM_RT0) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V)

         nfa = 4
         dim = 1
         vol = det / 2

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            iP3 = IPF(3, i)
            iP4 = IPF(4, i)
            
            sqr = calSqr(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3)) 
            sqr = sqr / vol
           
            Do n = 1, iGauss
               s = 0D0
               Do j = 1, 3
                 s = s + (V(j, iP1, 1) * (XYP(j, iP1) - XYP(j, iP4))
     &                  + V(j, iP2, 1) * (XYP(j, iP2) - XYP(j, iP4))
     &                  + V(j, iP3, 1) * (XYP(j, iP3) - XYP(j, iP4)))
               End do
               U(1, i, n) = sqr * s 
            End do

            vol = -vol
         End do

      Else If(FEMtype.EQ.FEM_ND0) Then
c  ...  check that FEM_ND0 is divergence-free
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V)

         nfa = 6
         dim = 1
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

            s = 0D0
            Do j = 1, 3
               s = s + V(j, iP1, 1) * xyn(j) + V(j, iP2, 1) * xym(j)
            End do
            U(1, i, 1) = sqr * s
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V)

         nfa = 12
         dim = 1 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 4
            U(1, i,     1) = V(1, i, 1)
            U(1, i + 4, 1) = V(2, i, 1)
            U(1, i + 8, 1) = V(3, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         Call applyGRAD(1, XYL, PSI, FEM_CR1, nfc, dim, V)

         nfa =12 
         dim = 1
         Call clearU(1, nfa, dim, U)

         Do i = 1, 4
            U(1, i,     1) = V(1, i, 1)
            U(1, i + 4, 1) = V(2, i, 1)
            U(1, i + 8, 1) = V(3, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2vector) Then
         Call applyGRAD(iGauss, XYL, PSI, FEM_P2, nfc, dim, V)

         nfa = 30
         dim = 1
         Call clearU(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            Do i = 1, 10
               U(1, i,      n) = V(1, i, n)
               U(1, i + 10, n) = V(2, i, n)
               U(1, i + 20, n) = V(3, i, n)
            End do
         End do

      Else
         nfa = 0
         dim = -FEMtype
      End if
      Return
      End

