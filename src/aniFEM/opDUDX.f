c ======================================================================
      Subroutine applyDUDX(iGauss, XYL, PSI, FEMtype, nfa, dim, U)
c ======================================================================
      include 'fem3Dtet.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(3, 3), U(9, MaxSize, *)

c Data for the reference triangle
      Real*8  XYL(4, *)
      Real*8  GRAD_P1(3, 4)

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
         dim = 1 

c  ...   choose ixy=1 for DuDx and ixy=2 for DuDy and ixy=3 for DuDz
         ixy = 1
         Do i = 1, nfa
            U(1, i, 1) = 0D0
            Do j = 1, 3
               U(1, i, 1) = U(1, i, 1) + PSI(j, ixy) * GRAD_P1(j, i)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

      Else
         nfa = 0
         dim = -FEMtype
      End if

      Return
      End




