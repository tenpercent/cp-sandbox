      function get_parallelogram_area(vertice_a, vertice_b, vertice_c)
        implicit none  

        real*8, dimension(1 : 3) :: vertice_a, vertice_b, vertice_c
        real*8, dimension (1 : 3) :: vector_ab, vector_ac
        real*8, dimension(1 : 3) :: cross_product
        real*8 :: get_parallelogram_area

        interface
          function get_cross_product (vector_ab, vector_ac)
            implicit none
            real*8, dimension(1 : 3) :: get_cross_product
            real*8, dimension (1 : 3) :: vector_ab, vector_ac
          end function get_cross_product

          function get_3d_vector (vertice_a, vertice_b)
            implicit none
            real*8, dimension(1 : 3) :: get_3d_vector
            real*8, dimension(1 : 3) :: vertice_a, vertice_b
          end function get_3d_vector

          function get_vector_length (vector)
            implicit none
            real*8 :: get_vector_length
            real*8, dimension (1 : 3) :: vector
          end function get_vector_length
        end interface

        vector_ab = get_3d_vector (vertice_a, vertice_b)
        vector_ac = get_3d_vector (vertice_a, vertice_c)

        cross_product = get_cross_product (vector_ab, vector_ac)

        get_parallelogram_area = get_vector_length (cross_product)

      end function get_parallelogram_area

      function get_vector_length (vector)
        implicit none
        real*8, dimension(1 : 3) :: vector
        real*8 :: get_vector_length

        get_vector_length = 
     &    dSqrt (vector(1) ** 2 + vector(2) ** 2 + vector(3) ** 2)
      end function get_vector_length

      function get_cross_product (vector_ab, vector_ac)
        implicit none
        real*8, dimension(1 : 3) :: get_cross_product
        real*8, dimension (1 : 3) :: vector_ab, vector_ac

        get_cross_product(1) = vector_ab(2) * vector_ac(3) - 
     &                         vector_ab(3) * vector_ac(2)
        get_cross_product(2) = vector_ab(3) * vector_ac(1) - 
     &                         vector_ab(1) * vector_ac(3)
        get_cross_product(3) = vector_ab(1) * vector_ac(2) - 
     &                         vector_ab(2) * vector_ac(1)
        
      end function get_cross_product

      function get_3d_vector (vertice_a, vertice_b)
        implicit none
        real*8, dimension(1 : 3) :: get_3d_vector
        real*8, dimension(1 : 3) :: vertice_a, vertice_b
        integer :: i
        do i = 1, 3
          get_3d_vector(i) = vertice_b(i) - vertice_a(i)
        end do
      end function get_3d_vector

      function get_normal (vertice_a, vertice_b, vertice_c)
        implicit none
        real*8, dimension (1 : 3) :: get_normal, 
     &                               vertice_a, vertice_b, vertice_c, 
     &                               vector_ab, vector_ac
        real*8  :: norm_coef
        integer :: i
        interface
          function get_vector_length (vector)
            implicit none
            real*8 :: get_vector_length
            real*8, dimension (1 : 3) :: vector
          end function get_vector_length

          function get_3d_vector (vertice_a, vertice_b)
            implicit none
            real*8, dimension(1 : 3) :: get_3d_vector
            real*8, dimension(1 : 3) :: vertice_a, vertice_b
          end function get_3d_vector

          function get_cross_product (vector_ab, vector_ac)
            implicit none
            real*8, dimension(1 : 3) :: get_cross_product
            real*8, dimension(1 : 3) :: vector_ab, vector_ac
          end function get_cross_product
        end interface
        vector_ab = get_3d_vector(vertice_a, vertice_b)
        vector_ac = get_3d_vector(vertice_a, vertice_c)
        get_normal = get_cross_product(vector_ab, vector_ac)
        norm_coef = 1 / get_vector_length(get_normal)
        do i = 1, 3
          get_normal(i) = get_normal(i) * norm_coef
        end do
      end function get_normal  

      program test_get_area
        implicit none
        real*8, dimension (1 : 3) :: vertice_a, vertice_b, vertice_c
        real*8 :: get_parallelogram_area
        real*8 :: area = 0
c ===   read vertices
        write (*,*) "input vertice_a coordinates:"
        read (*,*) vertice_a(1), vertice_a(2), vertice_a(3)

        write (*,*) "input vertice_b coordinates:"
        read (*,*) vertice_b(1), vertice_b(2), vertice_b(3)

        write (*,*) "input vertice_c coordinates:"
        read (*,*) vertice_c(1), vertice_c(2), vertice_c(3)
c ===   calculate
        area = get_parallelogram_area(vertice_a, vertice_b, vertice_c)
c ===   print parallelogram area
        write(*,'(A, 1F12.6)') "parallelogram area: ", area
      end program test_get_area
