      function det2 (x11, x12, x21, x22)
      implicit none

          real*8 :: x11, x12, x21, x22, det2

          det2 = x11 * x22 - x21 * x12

      return
      end function det2
c====================================================================
      function det4 (x11, x12, x13, x14,
     &               x21, x22, x23, x24,
     &               x31, x32, x33, x34,
     &               x41, x42, x43, x44)
      implicit none

          real*8 :: x11, x12, x13, x14,
     &              x21, x22, x23, x24,
     &              x31, x32, x33, x34,
     &              x41, x42, x43, x44

          real*8 :: det4, det2

          det4 = det2 (x11, x12, x21, x22) * det2 (x33, x34, x43, x44) -
     &           det2 (x13, x14, x23, x24) * det2 (x31, x32, x41, x42)

      return
      end function det4
c====================================================================
      function tetrahedron_volume (p1, p2, p3, p4)
      implicit none

          real*8, dimension(1:3) :: p1, p2, p3, p4
          real*8 :: tetrahedron_volume, det4

          tetrahedron_volume = abs (det4 (1d0, p1(1), p1(2), p1(3),
     &                               1d0, p2(1), p2(2), p2(3),
     &                               1d0, p3(1), p3(2), p3(3),
     &                               1d0, p4(1), p4(2), p4(3)) / 6d0)

      return
      end function tetrahedron_volume
c====================================================================
c==========================================================================================
      function euclid_distance (p1, p2)
      implicit none

        real*8, dimension (1:3) :: p1, p2
        real*8 :: euclid_distance

        euclid_distance = sqrt((p1(1) - p2(1)) ** 2 + 
     &                          (p1(2) - p2(2)) ** 2 + 
     &                          (p1(3) - p2(3)) ** 2)
      return
      end function euclid_distance
c==========================================================================================
      function search_nearest (node_coordinates, 
     &                         total_mesh_nodes, 
     &                         point_coordinates)
      implicit none
        integer :: total_mesh_nodes
        real*8, dimension (1:3, 1:total_mesh_nodes) :: node_coordinates
        real*8, dimension (1:3) :: point_coordinates
        integer :: search_nearest
        integer :: i

        real*8 :: euclid_distance

        real*8 :: min_distance, current_distance

        min_distance = euclid_distance (
     &      point_coordinates, 
     &      node_coordinates(1 : 3, 1))
        current_distance = min_distance

        do i = 2, total_mesh_nodes
            current_distance = euclid_distance (
     &          point_coordinates, 
     &          node_coordinates(1 : 3, i))
            if (current_distance < min_distance) then
                min_distance = current_distance
                search_nearest = i
            end if
        end do

      return
      end function search_nearest
c====================================================================
      function delta_dirac(tetrahedra_nodes, 
     &                     total_tetrahedra, 
     &                     node_coordinates, 
     &                     total_mesh_nodes, 
     &                     point_coordinates,
     &                     nearest_node_coordinates)
      implicit none

        integer :: total_tetrahedra, total_mesh_nodes
        integer, dimension (1 : 4, 1 : total_tetrahedra) ::
     &      tetrahedra_nodes
        real*8, dimension (1 : 3, 1 : total_mesh_nodes) ::
     &      node_coordinates
        real*8, dimension (1 : 3) :: point_coordinates, 
     &      nearest_node_coordinates
        real*8 :: delta_dirac

        integer :: nearest_node, search_nearest

        external backReferences

        integer, parameter :: L = 4, M = 4

        integer, dimension (1 : total_mesh_nodes) :: nEP
        integer, dimension (1 : 4 * total_tetrahedra) :: IEP

        real*8 :: volume = 0d0
        real*8 :: calVol

        integer :: i

        nearest_node = search_nearest(node_coordinates, 
     &                         total_mesh_nodes, 
     &                         point_coordinates)
        nearest_node_coordinates = 
     &    node_coordinates(1 : 3, nearest_node)

        call backReferences (total_mesh_nodes, 
     &                       total_tetrahedra, 
     &                       L, M, tetrahedra_nodes, nEP, IEP) 

        do i = nEP(nearest_node - 1) + 1, nEP(nearest_node)
            volume = volume + abs(
     &        calVol (
     &          node_coordinates(1 : 3, tetrahedra_nodes(1, IEP(i))), 
     &          node_coordinates(1 : 3, tetrahedra_nodes(2, IEP(i))),
     &          node_coordinates(1 : 3, tetrahedra_nodes(3, IEP(i))), 
     &          node_coordinates(1 : 3, tetrahedra_nodes(4, IEP(i)))))
        end do

        if (volume .eq. 0) then
          write (*,*) "volume eq 0"
          delta_dirac = 1
        else
          delta_dirac = 4d0 / volume
        end if
      return
      end function delta_dirac
