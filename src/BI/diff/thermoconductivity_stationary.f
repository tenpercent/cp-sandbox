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
        integer, dimension (1 : max_mesh_nodes) :: node_material_labels                   ! "labelP"

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
C ===   whatever
        real*8, dimension (1 : 4) :: DATAFEM
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
        real*8 :: delta_dirac
        real*8 :: delta_dirac_value
        real*8, dimension (1 : 3) :: source_coordinates
        real*8, dimension (1 : 3) :: nearest_node_coordinates
        data source_coordinates /3 * 5d-1/
        common /dirac/ delta_dirac_value, nearest_node_coordinates

        integer :: i
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
     &                max_mesh_nodes,                ! nvmax
     &                max_boundary_faces,            ! nbmax
     &                max_tetrahedra,                ! ntmax
     &                total_mesh_nodes,              ! nv
     &                total_boundary_faces,          ! nb
     &                total_tetrahedra,              ! nt
     &                node_coordinates,              ! vrt
     &                face_nodes,                    ! bnd
     &                tetrahedra_nodes,              ! tet
     &                face_material_labels,          ! labelF
     &                tetrahedra_material_labels,    ! material
     &                total_fixed_points,            ! nPv
     &                total_fixed_faces,             ! nFv
     &                total_fixed_elements,          ! nEv
     &                fixed_points,                  ! IPV
     &                fixed_faces,                   ! IFV
     &                fixed_elements,                ! IEV
     &                work_buffer,                   ! iW
     &                work_buffer,                   ! iW
     &                "../../../../data/cube.ani")   
        write (*,*) "loadMani has finished!"

        do i = 1, total_boundary_faces
c ===   only neumann boundary conditions    
            face_material_labels(i) = 1
        end do

        do i = 1, total_mesh_nodes
c ===   again, no dirichlet boundary condition
            node_material_labels(i) = 1
        end do            

c === mesh is ready; calculate dirac function
        delta_dirac_value = delta_dirac(tetrahedra_nodes, 
     &                     total_tetrahedra, 
     &                     node_coordinates, 
     &                     total_mesh_nodes, 
     &                     source_coordinates,
     &                     nearest_node_coordinates)
        write (*,*) "Dirac function has been computed!"
        write (*,*) delta_dirac_value
        
c ===   set bit mask for matrix type
        matrix_type = IOR(MATRIX_GENERAL, FORMAT_AMG)
c ===   Assemble the stiffness matrix
c       no extra data is provided for the user subroutines Dxxxx
        DATAFEM(1) = delta_dirac_value
        DATAFEM(2) = nearest_node_coordinates(1)
        DATAFEM(3) = nearest_node_coordinates(2)
        DATAFEM(4) = nearest_node_coordinates(3)
c === test if mesh is initialized

c        write (*,*) "boundary labels before: "
c        do i = 1, total_boundary_faces
c            write (*, *) face_material_labels(i)
c        end do

c        write (*,*) "element labels before: "
c        do i = 1, total_tetrahedra
c            write (*, *) tetrahedra_material_labels(i)
c        end do

        write (*,*) "point labels before: "
        do i = 1, total_mesh_nodes
            write (*, *) node_material_labels(i)
        end do

c ===   construct finite elements
        call BilinearFormTemplate(
     &                            total_mesh_nodes,                     ! nv
     &                            total_boundary_faces,                 ! nb
     &                            total_tetrahedra,                     ! nt
     &                            node_coordinates,                     ! vrt
     &                            node_material_labels,                 ! labelP
     &                            face_nodes,                           ! bnd
     &                            face_material_labels,                 ! labelF
     &                            tetrahedra_nodes,                     ! tet
     &                            tetrahedra_material_labels,           ! material
     &                            FEM3Dext,                             ! 
     &                            DATAFEM,                              ! 
     &                            matrix_type,                          ! status
     &                            max_mesh_nodes,                       ! nvmax
     &                            max_nonzero_matrix_entries,           ! namax
     &                            sr_nonzero_entries_columns,           ! IA
     &                            sr_row_nonzero_counter,               ! JA
     &                            sr_nonzero_entries,                   ! A
     &                            sle_right_hand_side,                  ! RHS
     &                            total_rows,                           ! nRow
     &                            total_columns,                        ! nCol
     &                            max_work_memory_size,                 ! MaxWi
     &                            work_buffer)                          ! iW

        write (*,*) "BilinearFormTemplate has finished!"

c ===   test what's up with mesh
c        write (*,*) "boundary labels after: "
c        do i = 1, total_boundary_faces
c            write (*, *) face_material_labels(i)
c        end do

c        write (*,*) "element labels after: "
c        do i = 1, total_tetrahedra
c            write (*, *) tetrahedra_material_labels(i)
c        end do

        write (*,*) "point labels after: "
        do i = 1, total_mesh_nodes
            write (*, *) node_material_labels(i)
        end do

c ===   launch solver
        ls_status = LSolver(
     &                 ls_method,                                        ! isolv
     &                 ls_preconditioner,                                ! iprec
     &                 total_rows,                                       ! nRow
     &                 sr_nonzero_entries_columns,                       ! IA
     &                 sr_row_nonzero_counter,                           ! JA
     &                 sr_nonzero_entries,                               ! A
     &                 sle_right_hand_side,                              ! b
     &                 sle_solution)                                     ! x

      write (*,*) "LSolver has finished!"

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

      write (*,*) "GMVscalarTet has finished!"

c ===   I may be wrong about that      
      end program thermoconductivity

C ======================================================================
C  The user defined routines required above
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
c ======================================================================
c This program generates a finite element system for the diffusion problem
c
c  -div K grad u + A u = F     in  Omega
c                    u = U_0   on  Gamma_D
c        K du/dn       = G_0   on  Gamma_N
c        K du/dn + S u = G_1   on  Gamma_R
c
c where Omega = [0,1]^3. The user-defined coefficients are  
c    K(x)   - positive definite tensor
c    A(x)   - non-negative reaction 
c    F(x)   - right-hand side
c    U_0(x) - essential (Dirichlet) boundary condition
c    G_0(x) - Neumann boundary condition
c    G_1(x) - Robin boundary condition
c    S(x)   - Robin boundary coefficient
c
c These coefficients are implemented in the following routines: 
c    K->Ddiff,  A->Dreact,  F->Drhs,  {U_0,G_0,G_1}->Dbc,  S->DbcRobCoef
c
c ======================================================================
      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
     &                    lbE, lbF, lbR, lbP, DATAFEM, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      
      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DATAFEM(*)
      Integer iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C LOCAL VARIABLEs
      Integer  Ddiff, Dreact, Drhs, Dbc
      External Ddiff, Dreact, Drhs, Dbc

      Real*8   B(4, 4), C(3, 3), G(3), XYP(3, 4)
      Real*8   x, y, z, eBC(1)

      Integer  i,j,k,l,m, ir, ic, label, ibc
  
      Integer  iref(5), ip(4)
      DATA     iref /1,2,3,4,1/
      
      Logical print_mesh_labels
      parameter (print_mesh_labels = .true.) 

C ======================================================================
      nRow = 4
      nCol = 4

c ... set up templated degrees of freedom for rows and columns. 
c     used convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
      Do i = 1, 4
         templateC(i) = Vdof
         templateR(i) = Vdof
      End do

      if (print_mesh_labels) then
          write (*,'(a, 1i2)') "element label: ", lbE
          write (*,'(a, 4i2)') "face labels: ", lbF
          write (*,'(a, 6i2)') "edge labels: ", lbR
          write (*,'(a, 4i2)') "point labels: ", lbP
      end if

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:4,1:4) is elemental vector elliptic operator;
c     in other words, for the bilinear form <grad(P1), grad(P1)>
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DATAFEM, iSYS, 2,
     &              LDA, A, ir, ic)
c === === === === === === === === === === === === === === ===

c     B(1:4,1:4) is elemental mass matrix;
c     in other words, for the bilinear form <P1, P1>
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, Dreact, DATAFEM, iSYS, 5,
     &              4, B, ir, ic)

      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + B(i, j) 
         End do
      End do

c ... compute the right hand side vector using external function Drhs
c     in other words, the linear form <Drhs(x), P1> 
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1, 
     &              lbE, Drhs, DATAFEM, iSYS, 5,
     &              LDA, F, ir, ic)

c ... impose the Neumann boundary conditions
      Do i = 1, 3
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
         XYP(i, 4) = XY4(i)
      End do

      Do k = 1, 4
c === iterate through faces
        if (lbF(k) .eq. 1) then      
          l = iref(k + 1)
          m = iref(l + 1)

          x = (XYP(1, k) + XYP(1, l) + XYP(1, m)) / 3
          y = (XYP(2, k) + XYP(2, l) + XYP(2, m)) / 3
          z = (XYP(3, k) + XYP(3, l) + XYP(3, m)) / 3

          ibc = Dbc(x, y, z, lbF(k), DATAFEM, iSYS, eBC)
           
          label = lbF(k)

          Call fem3Dtri(XYP(1, k), XYP(1, l), XYP(1, m),
     &                    IDEN, FEM_P0, IDEN, FEM_P1, 
     &                    label, Dbc, DATAFEM, iSYS, 4, 
     &                    3, G, ir, ic)

          F(k) = F(k) + G(1)
          F(l) = F(l) + G(2)
          F(m) = F(m) + G(3)
        else
          if (print_mesh_labels) then
            write (*,'(a, i1, a)') "labelF eq ", 
     &                             lbF(k), 
     &                             " neq 1; something's wrong"
          end if
        end if
      End do

      Return
      End

C ======================================================================
c diffusion tensor K
C ======================================================================

      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, *)
      Integer label, iSYS(*)

      Integer i, j

      Integer thermo_coefficient
      parameter (thermo_coefficient = 1d0)
c === should fix it later      

      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
      End do

      Coef(1, 1) = thermo_coefficient
      Coef(2, 2) = thermo_coefficient
      Coef(3, 3) = thermo_coefficient

      Ddiff = TENSOR_SYMMETRIC

      Return
      End

C ======================================================================
c Reaction coefficient A
C ======================================================================

      Integer Function Dreact(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1
c     A == 0
      Coef(1) = 0D0
      Dreact = TENSOR_SCALAR

      Return
      End

C ======================================================================
c Boundary conditions 
C ======================================================================

      Integer Function Dbc(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      Include 'assemble.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Dbc = BC_NEUMANN

c === \int\limits_{\partial\Omega} g_N d S = 1.
c === where \Omega = [0, 1]^3 
c === and g_N = const

      Coef(1) = 1d0 / 6d0

      Return
      End

C ======================================================================
c Right hand side F
C ======================================================================

      Integer Function Drhs(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      Include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(*)
      Integer label, iSYS(*)

      real*8 dirac_value
      real*8 dirac_nz_coordinates (3) 

      real*8 min_dist
      parameter (min_dist = 1d-6)

      dirac_value = DATA(1)
      dirac_nz_coordinates(1 : 3) = DATA(2 : 4)

      iSYS(1) = 1
      iSYS(2) = 1

      if    (abs(x - dirac_nz_coordinates(1)) < min_dist
     & .and. abs(y - dirac_nz_coordinates(2)) < min_dist
     & .and. abs(z - dirac_nz_coordinates(3)) < min_dist
     &    ) then 
        coef(1) = dirac_value
      else
        coef(1) = 0d0
      end if

      Drhs = TENSOR_SCALAR

      Return
      End