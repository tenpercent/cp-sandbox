      subroutine getLInterpCoef_Detailed (x1, x2, x3, x4, 
     &                                    y1, y2, y3, y4, 
     &                                    z1, z2, z3, z4, 
     &                                    a,  b,  c,  d)

c ===  | A x1 + B y1 + C z1 + D = 0
c ===  | A x2 + B y2 + C z2 + D = 0
c ===  | A x3 + B y3 + C z3 + D = 0
c ===  | A x4 + B y4 + C z4 + D = 1

        real*8 :: x1, x2, x3, x4, 
     &            y1, y2, y3, y4, 
     &            z1, z2, z3, z4, 
     &            a, b, c, d

        real*8 :: det, det_a, det_b, det_c, det_d

        real*8 :: det4
        real*8, parameter :: mindouble = 1d-8
c =====================================================================
        det = det4(x1, y1, z1, 1d0,
     &             x2, y2, z2, 1d0,
     &             x3, y3, z3, 1d0,
     &             x4, y4, z4, 1d0)

        det_a = det4(0d0, y1, z1, 1d0,
     &               0d0, y2, z2, 1d0,
     &               0d0, y3, z3, 1d0,
     &               1d0, y4, z4, 1d0)

        det_b = det4(x1, 0d0, z1, 1d0,
     &               x2, 0d0, z2, 1d0,
     &               x3, 0d0, z3, 1d0,
     &               x4, 1d0, z4, 1d0)

        det_c = det4(x1, y1, 0d0, 1d0,
     &               x2, y2, 0d0, 1d0,
     &               x3, y3, 0d0, 1d0,
     &               x4, y4, 1d0, 1d0)

        det_d = det4(x1, y1, z1, 0d0,
     &               x2, y2, z2, 0d0,
     &               x3, y3, z3, 0d0,
     &               x4, y4, z4, 1d0)

        if (abs(det) < mindouble) then
          write (*,*) 'no solutions in interpolation'
          return
        end if

        a = det_a / det
        b = det_b / det
        c = det_c / det
        d = det_d / det

        return
      end subroutine getLInterpCoef_Detailed
  

      subroutine getLInterpCoef (p1, p2, p3, p4, A, B, C, D)
c ===   F(p1) = F(p2) = F(p3) = 0; F(p4) = 1
c ===   F(x, y, z) = Ax + By + Cz + D
        real*8, dimension(1:3) :: p1, p2, p3, p4
        real*8 :: A, B, C, D

        call getLInterpCoef_Detailed(p1(1), p2(1), p3(1), p4(1),
     &                               p1(2), p2(2), p3(2), p4(2),
     &                               p1(3), p2(3), p3(3), p4(3),
     &                               A,     B,     C,     D)
        return
      end subroutine getLInterpCoef


      function LInterpValue (p1, p2, p3, p4, point)
c ===   F(p1) = F(p2) = F(p3) = 0; F(p4) = 1
c ===   F(x, y, z) = Ax + By + Cz + D
c ===   return F(point)
        real*8 :: LInterpValue
        real*8, dimension(1:3) :: p1, p2, p3, p4, point
c ===   local variables
        real*8 :: a, b, c, d
c ============
        call getLInterpCoef(p1, p2, p3, p4, a, b, c, d)

        LInterpValue = a * point(1) + b * point(2) + c * point(3) + d

        return
      end function LInterpValue
