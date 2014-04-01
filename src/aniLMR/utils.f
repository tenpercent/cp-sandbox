c ======================================================================
      Subroutine GradU(xy1, xy2, xy3, xy4, u1, u2, u3, u4, GRAD)
c ======================================================================
c Routine computes gradient of a linear function given
c by value u_i at points xy_i. 
c ======================================================================
      Real*8  xy1(*), xy2(*), xy3(*), xy4(*), GRAD(*)
      Real*8  u1, u2, u3, u4

      Real*8  A(3, 3)
      Integer ipiv(3), info
c ======================================================================
      Do i = 1, 3
         A(i, 1) = xy1(i) - xy4(i)
         A(i, 2) = xy2(i) - xy4(i)
         A(i, 3) = xy3(i) - xy4(i)
      End do

      GRAD(1) = u1 - u4
      GRAD(2) = u2 - u4
      GRAD(3) = u3 - u4

      Call dgesv(3, 1, A, 3, ipiv, GRAD, 3, info)

      If(info.NE.0) Then
         Call errMesLMR(3011, "GradU", "Error in LAPACK library") 
      End if

      Return
      End

