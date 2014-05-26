      function PointInTet (p1, p2, p3, p4, point)
          implicit none
c ===     return type
          logical :: PointInTet
c ===     input variables
          real*8, dimension(1:3) :: p1, p2, p3, p4, point
c ===     local variables
          real*8 :: vol_1, vol_2, vol_3, vol_4, vol
c ===     ./dirac.f
          real*8 :: tetrahedron_volume

          real*8, parameter :: mindouble = 1d-8
c === === === === === === ===
          vol_1 = tetrahedron_volume(p2, p3, p4, point)
          vol_2 = tetrahedron_volume(p1, p3, p4, point)
          vol_3 = tetrahedron_volume(p1, p2, p4, point)
          vol_4 = tetrahedron_volume(p1, p2, p3, point)

          vol = tetrahedron_volume(p1, p2, p3, p4)

          PointInTet = (abs(vol - vol_1 - vol_2 - vol_3 - vol_4) 
     &                  < mindouble)
          return
      end function PointInTet
