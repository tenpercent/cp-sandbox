C =========================================================
      subroutine tetL8error(U, nv, vrt, nt, tet, error)
C =========================================================
      implicit none

      Real*8  U(*), vrt(3,*), error(*)
      Integer tet(4, *), nv, nt

      Real*8   Func
      EXTERNAL Func

      Integer   nStep
      Parameter(nStep = 20)

      Integer i, j, k, l, n, iv1, iv2, iv3, iv4
      Real*8  xyz(3), uint, diff, da, d1, d2, d3, d4

C =========================================================
      da = 1d0 / nStep

      Do n = 1, nt
         error(n) = 0d0

         iv1 = tet(1, n)
         iv2 = tet(2, n)
         iv3 = tet(3, n)
         iv4 = tet(4, n)

         Do i = 0, nStep
         Do j = 0, nStep - i
         Do k = 0, nStep - i - j
c ... baricentric coordinates
            d1 = da * i
            d2 = da * j
            d3 = da * k
            d4 = max(0d0, 1d0 - (d1 + d2 + d3))

            Do l = 1, 3
               xyz(l) = d1 * vrt(l, iv1) + d2 * vrt(l, iv2) 
     &                + d3 * vrt(l, iv3) + d4 * vrt(l, iv4) 
            End do             

            uint = d1 * U(iv1) + d2 * U(iv2) 
     &           + d3 * U(iv3) + d4 * U(iv4) 

            diff = dabs(Func(xyz) - uint)
            error(n) = max(error(n), diff)
         End do
         End do
         End do
      End do

      Return
      End

         

