c === det (2 x 2)
      function det2 (x11, x12, x21, x22)
          implicit none

          real*8 :: x11, x12, x21, x22, det2

          det2 = x11 * x22 - x21 * x12

      return
      end function det2


c === det (3 x 3)
      function det3 (x11, x12, x13, x21, x22, x23, x31, x32, x33)
          implicit none
          real*8 :: x11, x12, x13, x21, x22, x23, x31, x32, x33, det3

          det3 =   x11 * x22 * x33
     &           + x12 * x23 * x31
     &           + x13 * x21 * x32
     &           - x12 * x21 * x33
     &           - x13 * x22 * x31
     &           - x11 * x23 * x32
          return
      end function det3


c === det (4 x 4)
      function det4 (x1, x2, x3, x4,
     &               y1, y2, y3, y4,
     &               z1, z2, z3, z4,
     &               t1, t2, t3, t4)
          implicit none

          real*8 ::  x1, x2, x3, x4,
     &               y1, y2, y3, y4,
     &               z1, z2, z3, z4,
     &               t1, t2, t3, t4

          real*8 :: det4

          det4 = - t4 * x3 * y2 * z1 
     &           + t3 * x4 * y2 * z1 
     &           + t4 * x2 * y3 * z1 
     &           - t2 * x4 * y3 * z1 
     &           - t3 * x2 * y4 * z1 
     &           + t2 * x3 * y4 * z1 
     &           + t4 * x3 * y1 * z2 
     &           - t3 * x4 * y1 * z2 
     &           - t4 * x1 * y3 * z2 
     &           + t1 * x4 * y3 * z2 
     &           + t3 * x1 * y4 * z2 
     &           - t1 * x3 * y4 * z2 
     &           - t4 * x2 * y1 * z3 
     &           + t2 * x4 * y1 * z3 
     &           + t4 * x1 * y2 * z3 
     &           - t1 * x4 * y2 * z3 
     &           - t2 * x1 * y4 * z3 
     &           + t1 * x2 * y4 * z3 
     &           + t3 * x2 * y1 * z4 
     &           - t2 * x3 * y1 * z4 
     &           - t3 * x1 * y2 * z4 
     &           + t1 * x3 * y2 * z4 
     &           + t2 * x1 * y3 * z4 
     &           - t1 * x2 * y3 * z4

          return
      end function det4
      