C ======================================================================
      Subroutine triIntL8(U, nv, vrt, nt, tet, error)
C ======================================================================
C  Routine calculates the maximum norm of the interpolation error in
C  tetrahedra by spliting them in smaller tetrahedra.
C ======================================================================
      Integer   NX
      Parameter(NX = 20)
C ======================================================================
      Real*8   U(*), vrt(3, *), error(*)
      Integer  tet(4, *)

      Real*8   xyt(3), l1, l2, l3, l4, da, uint, diff
      
      Real*8   Func
      EXTERNAL Func

C ======================================================================
      da = 1D0 / NX

      Do n = 1, nt
         error(n) = 0D0

         iP1 = tet(1, n)
         iP2 = tet(2, n)
         iP3 = tet(3, n)
         iP4 = tet(4, n)

         Do i = 0, NX
         Do j = 0, max(0, NX-i)
         Do k = 0, max(0, NX-i-j)
            l1 = i * da
            l2 = j * da
            l3 = k * da
            l4 = 1D0 - l1 - l2 - l3

            Do m = 1, 3
               xyt(m) = vrt(m, iP1) * l1 + vrt(m, iP2) * l2  
     &                + vrt(m, iP3) * l3 + vrt(m, iP4) * l4 
            End do

            uint = U(iP1) * l1 + U(iP2) * l2 + U(iP3) * l3 + U(iP4) * l4
            diff = dabs(Func(xyt) - uint)
            error(n) = max(error(n), diff)
         End do
         End do
         End do
      End do

      Return
      End


         
C ======================================================================
      Subroutine edgeIntL8(U, nv, vrt, nt, tet, error)
C ======================================================================
C  Routine calculates the maximum norm of the interpolation error on
C  mesh edges by spltting edges into small subedges.
C ======================================================================
      Integer   NX
      Parameter(NX = 20)
C =========================================================
      Real*8   U(*), vrt(3, *), error(6, *)
      Integer  tet(4, *)

      Real*8   xyt(3), l1, l2, da, uint, diff

      Real*8   Func
      EXTERNAL Func

C =========================================================
      da = 1D0 / NX

      Do n = 1, nt
         k = 0
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
               k = k + 1
               error(k, n) = 0D0

               iP1 = tet(i1, n)
               iP2 = tet(i2, n)

               Do i = 0, NX
                  l1 = i * da
                  l2 = 1D0 - l1
                  Do m = 1, 3
                     xyt(m) = vrt(m, iP1) * l1 + vrt(m, iP2) * l2
                  End do

                  uint = U(iP1) * l1 + U(iP2) * l2
                  diff = dabs(Func(xyt) - uint)
                  error(k, n) = max(error(k, n), diff)
               End do
            End do
         End do
      End do

      Return
      End



