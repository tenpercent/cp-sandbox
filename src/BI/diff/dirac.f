      program test_dirac
        implicit none
        integer, parameter :: max_tetrahedra = 400000,                                    ! "ntmax"
     &                        max_mesh_nodes = 70000,                                     ! "nvmax"
     &                        max_boundary_faces = 5000,                                  ! "nbmax"
     &                        max_nonzero_matrix_entries = 1000000                        ! "namax"
        integer, parameter :: max_work_memory_size = 8000000                              ! "MaxWi"
        integer, dimension (1 : max_work_memory_size) :: work_buffer                      ! "iW"

        integer :: total_mesh_nodes,                                                      ! "nv"
     &             total_tetrahedra,                                                      ! "nt"
     &             total_boundary_faces                                                   ! "nb"

        real*8, dimension (1 : 3, 1 : max_mesh_nodes) :: 
     &             node_coordinates                                                       ! "vrt"
        integer :: node_material_labels                                                   ! "labelP"

        integer, dimension (1 : 4, 1 : max_tetrahedra) ::
     &             tetrahedra_nodes                                                       ! "tet"
        integer, dimension (1 : max_tetrahedra) :: 
     &             tetrahedra_material_labels                                             ! "material"

        integer, dimension (1 : 3, 1 : max_boundary_faces) :: face_nodes                  ! "bnd"
        integer, dimension (1 : max_boundary_faces) :: 
     &             face_material_labels                                                   ! "labelF"

        integer, parameter :: max_fixed_points = 9,                                       ! "maxPV"
     &                        max_fixed_faces = 1,                                        ! "maxFV"
     &                        max_fixed_elements = 1                                      ! "maxEV"

        integer :: total_fixed_points,                                                    ! "nPv"
     &             total_fixed_faces,                                                     ! "nFv"
     &             total_fixed_elements                                                   ! "nEv"
        integer, dimension (1 : max_fixed_points)   :: fixed_points                       ! "IPV"
        integer, dimension (1 : max_fixed_faces)    :: fixed_faces                        ! "IFV"
        integer, dimension (1 : max_fixed_elements) :: fixed_elements                     ! "IEV"
C === AniFEM stuff        
        include 'fem3Dtet.fd'
        include 'assemble.fd'
C === "sr" means "sparse row"
        integer, dimension (1 : max_mesh_nodes) :: 
     &             sr_row_nonzero_counter                                                 ! "JA"
        integer, dimension (1 : max_nonzero_matrix_entries) :: 
     &             sr_nonzero_entries_columns                                             ! "IA"
        real*8, dimension (1 : max_nonzero_matrix_entries) :: 
     &             sr_nonzero_entries                                                     ! "A"
C === "sle" means "system of linear equations"
        real*8, dimension (1 : max_mesh_nodes) :: sle_right_hand_side                     ! "RHS"
        real*8, dimension (1 : max_mesh_nodes) :: sle_solution                            ! "SOL"
        integer :: total_columns, total_rows                                              ! "nCol", "nRow"
        integer :: matrix_type                                                            ! "status"
C === seems that these are arrays of coefficients to be changed
        real*8, dimension (1 : 44) :: DATAFEM_R                                           ! "DATAFEMR"
        real*8, dimension (1 : 1) :: DATAFEM_E                                            ! "DATAFEME"
c === external procedures block
        external loadMani

        real*8, dimension (1 : 3) :: point_coordinates
        data point_coordinates /5d-1, 5d-1, 5d-1/
        real*8 :: dirac_value, delta_dirac

        call loadMani(
     &                max_mesh_nodes, 
     &                max_boundary_faces, 
     &                max_tetrahedra,
     &                total_mesh_nodes, 
     &                total_boundary_faces, 
     &                total_tetrahedra,
     &                node_coordinates, 
     &                face_nodes, 
     &                tetrahedra_nodes, 
     &                face_material_labels, 
     &                tetrahedra_material_labels,
     &                total_fixed_points, 
     &                total_fixed_faces, 
     &                total_fixed_elements, 
     &                fixed_points, 
     &                fixed_faces, 
     &                fixed_elements,
     &                work_buffer, 
     &                work_buffer, 
     &                "../../../../data/cube.ani")

        dirac_value = delta_dirac(tetrahedra_nodes, 
     &                     total_tetrahedra, 
     &                     node_coordinates, 
     &                     total_mesh_nodes, 
     &                     point_coordinates)

        write (*,*) dirac_value

      end program test_dirac
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
     &                     point_coordinates)
      implicit none

        integer :: total_tetrahedra, total_mesh_nodes
        integer, dimension (1 : 4, 1 : total_tetrahedra) ::
     &      tetrahedra_nodes
        real*8, dimension (1 : 3, 1 : total_mesh_nodes) ::
     &      node_coordinates
        real*8, dimension (1 : 3) :: point_coordinates
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
c === check if nearest_node is not 0
        write (*, *) nearest_node

        call backReferences (total_mesh_nodes, 
     &                       total_tetrahedra, 
     &                       L, M, tetrahedra_nodes, nEP, IEP) 

c === check if tetrahedra_nodes is not empty
c        do i = 1, total_tetrahedra
c          write (*,*) tetrahedra_nodes(1:4, i)
c        end do
c        write (*,*) nEP(1 : total_mesh_nodes)

        do i = nEP(nearest_node - 1) + 1, nEP(nearest_node)
c === print tetrahedra which contain nearest node; fail
            write (*, *) tetrahedra_nodes(1 : 4, IEP(i))
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
