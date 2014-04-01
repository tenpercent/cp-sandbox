c ======================================================================
      Subroutine applyCURL(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, det)
c ======================================================================
      implicit none
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), XYL(4, *), U(9, MaxSize, *)
      Real*8  XYP(3, *), det

c Local variables 
      Integer i, j, k, n, iP1,iP2,iP3,iP4, nfc
      Real*8  V(9, MaxSize, MaxPointGauss)
      Real*8  vol, sqr, calEdge, xyn(3), xym(3), grad_xyz(3)

c Data for the reference tetrahedron
      Integer IPR(4, 6)
      Real*8  GRAD_P1(3, 4)

      DATA    IPR/1,2,3,4, 1,3,4,2, 1,4,2,3, 2,3,1,4, 2,4,3,1, 3,4,1,2/ 
      DATA    GRAD_P1/-1,-1,-1, 1,0,0, 0,1,0, 0,0,1/

c ======================================================================
      If(FEMtype.EQ.FEM_RT0) Then
         nfa = 4
         dim = 1

         Do i = 1, nfa
            Do n = 1, iGauss
               U(1, i, n) = 0D0
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_ND0) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V)

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

            U(1, i, 1) = sqr
     &           * (V(2, iP1, 1) * xyn(3) + V(2, iP2, 1) * xym(3)
     &            - V(3, iP1, 1) * xyn(2) - V(3, iP2, 1) * xym(2))

            U(2, i, 1) = sqr
     &           * (V(3, iP1, 1) * xyn(1) + V(3, iP2, 1) * xym(1)
     &            - V(1, iP1, 1) * xyn(3) - V(1, iP2, 1) * xym(3))

            U(3, i, 1) = sqr
     &           * (V(1, iP1, 1) * xyn(2) + V(1, iP2, 1) * xym(2)
     &            - V(2, iP1, 1) * xyn(1) - V(2, iP2, 1) * xym(1))
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 12
         dim = 3
         Call clearU(1, nfa, dim, U)

         Do i = 1, 4
            Do k = 1, 3
               grad_xyz(k) = 0D0
               Do j = 1, 3
                  grad_xyz(k) = grad_xyz(k) + PSI(j, k) * GRAD_P1(j, i)
               End do

               U(2, i, 1) =  grad_xyz(3)
               U(3, i, 1) = -grad_xyz(2)

               U(1, i + 4, 1) = -grad_xyz(3)
               U(3, i + 4, 1) =  grad_xyz(1)

               U(1, i + 8, 1) =  grad_xyz(2)
               U(2, i + 8, 1) = -grad_xyz(1)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

      Else
         nfa = 0
         dim = -FEMtype
      End if
      Return
      End

