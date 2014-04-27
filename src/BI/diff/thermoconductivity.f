      program thermoconductivity
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
        external BilinearFormTemplate
        external FEM3Dext
        integer :: LSolver
        external GMVscalarTet

C === local variables block
c === "ls" means "linear solver"
        logical :: ls_parameters_file_exists
        integer :: ls_method, ls_preconditioner
        integer :: ls_status
C =========================
c ===   load solver parameters
        inquire(file = "solver.txt", exist = ls_parameters_file_exists)
        if (ls_parameters_file_exists) then                                 ! solver.txt exists, load it 
          open(12, file = 'solver.txt')
          read(12, *) ls_method, ls_preconditioner
          close(12)
        else                                                                ! otherwise use default settings
          ls_method = 0
          ls_preconditioner = 0
        end if
c ===   load mesh
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
c ===   impose Dirichlet boundary conditions
c ===   but not now

c ===   read parameters
c ===   but not now

c ===   set bit mask for matrix type
        matrix_type = IOR(MATRIX_GENERAL, FORMAT_AMG)
c ===   construct finite elements
        call BilinearFormTemplate(
     &                            total_mesh_nodes, 
     &                            total_boundary_faces, 
     &                            total_tetrahedra, 
     &                            node_coordinates, 
     &                            node_material_labels, 
     &                            face_nodes, 
     &                            face_material_labels, 
     &                            tetrahedra_nodes, 
     &                            tetrahedra_material_labels,
     &                            FEM3Dext, 
     &                            DATAFEM_R, 
     &                            matrix_type,
     &                            max_mesh_nodes, 
     &                            max_nonzero_matrix_entries, 
     &                            sr_nonzero_entries_columns, 
     &                            sr_row_nonzero_counter, 
     &                            sr_nonzero_entries, 
     &                            sle_right_hand_side, 
     &                            total_rows, 
     &                            total_columns,
     &                            max_work_memory_size, 
     &                            work_buffer)
c ===   launch solver
        ls_status = LSolver(
     &                 ls_method,
     &                 ls_preconditioner, 
     &                 total_rows,
     &                 sr_row_nonzero_counter,
     &                 sr_nonzero_entries_columns,
     &                 sr_nonzero_entries,
     &                 sle_right_hand_side,
     &                 sle_solution)
c ===   write results to VTK file
        call GMVscalarTet(sle_solution, 
     &                    "mesh.vtk", 
     &                    10,
     &                    total_mesh_nodes,
     &                    node_coordinates, 
     &                    total_tetrahedra,
     &                    tetrahedra_nodes, 
     &                    total_boundary_faces,
     &                    face_nodes,
     &                    face_material_labels)  
c ===   I may be wrong about that      
      end program thermoconductivity

c      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
c     &                    lbE, lbF, lbR, lbP, DATAFEM, iSYS,
c     &                    LDA, A, F, nRow, nCol,
c     &                    templateR, templateC)
C ======================================================================
c      implicit none
c      Include 'fem3Dtet.fd'
c      Include 'assemble.fd'

c      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
c      
c      Integer lbE, lbF(4), lbR(6), lbP(4)
c      Real*8  DATAFEM(*)
c      Integer iSYS(*), LDA, nRow, nCol

c      Real*8  A(LDA, *), F(*)
c      Integer templateR(*), templateC(*)

c      return
c      end
c====================================================================
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
      function euclid_distance (p1, p2)
      implicit none

        real*8, dimension (1:3) :: p1, p2
        real*8 :: euclid_distance

        euclid_distance = sqrt((p1(1) - p2(1)) ** 2 + 
     &                          (p1(2) - p2(2)) ** 2 + 
     &                          (p1(3) - p2(3)) ** 2)
      return
      end function euclid_distance
c====================================================================
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

        integer, parameter :: L = 4, M = 4, max_adjacent = 1000

        integer, dimension (1 : M, 1 : max_adjacent) :: IPE

        integer, dimension (1 : total_mesh_nodes) :: nEP, IEP

        real*8 :: volume = 0d0
        real*8 :: calVol

        integer :: i

        nearest_node = search_nearest(node_coordinates, 
     &                         total_mesh_nodes, 
     &                         point_coordinates)

        call backReferences (total_mesh_nodes, 
     &                       total_tetrahedra, 
     &                       L, M, IPE, nEP, IPE) 

        do i = nEP(nearest_node) + 1, nEP(nearest_node + 1)
            volume = volume + abs(calVol (IPE(1, i), 
     &                                    IPE(2, i), 
     &                                    IPE(3, i), 
     &                                    IPE(4, i)))
        end do

        delta_dirac = 4d0 / volume
      return
      end function delta_dirac

