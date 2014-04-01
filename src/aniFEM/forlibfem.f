C ======================================================================
      Integer Function ANI_Dnull(x,y,z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(*)
      Integer iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 1D0
      ANI_Dnull = TENSOR_NULL

      Return
      End



C ======================================================================
C Example of identity tensor
C ======================================================================
      Integer Function ANI_D3x3_one(x,y,z, label, DATA, iSYS, Coef)
      include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer iSYS(*)
C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = 1D0
      End do

      ANI_D3x3_one = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Example of Dirichlet boundary condition
C ======================================================================
      Integer Function ANI_D1x1_zero(x,y,z, label, DATA, iSYS, Coef)
      Include 'assemble.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(*)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 0D0
      ANI_D1x1_zero = BC_DIRICHLET

      Return
      End



C ======================================================================
C Example of source a function
C ======================================================================
      Integer Function ANI_D1x1_one(x,y,z, label, DATA, iSYS, Coef)
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(*)
      Integer iSYS(*)

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1) = 1D0
      ANI_D1x1_one = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Example of a constant convection
C ======================================================================
      Integer Function ANI_Dconv_const(x,y,z, label, DATA, iSYS, Coef)
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 3

      Do i = 1, 3
         Coef(1, i) = DATA(i)
      End do
      ANI_Dconv_const = TENSOR_GENERAL

      Return
      End



C ======================================================================
C Example of a constant elasticity
C ======================================================================
      Integer Function ANI_Delas_const(x,y,z, label, DATA, iSYS, Coef)
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DATA(*), Coef(9, 9)
      Integer iSYS(*)
C ======================================================================
C Local variables
      Real*8  LameM, LameL
      Integer im(9, 9)

      DATA    im/2, 0, 0, 0, 0, 0, 0, 0, 0,
     &           0, 1, 0, 1, 0, 0, 0, 0, 0,
     &           0, 0, 1, 0, 0, 0, 1, 0, 0,
     &           0, 1, 0, 1, 0, 0, 0, 0, 0,
     &           0, 0, 0, 0, 2, 0, 0, 0, 0,
     &           0, 0, 0, 0, 0, 1, 0, 1, 0,
     &           0, 0, 1, 0, 0, 0, 1, 0, 0,
     &           0, 0, 0, 0, 0, 1, 0, 1, 0,
     &           0, 0, 0, 0, 0, 0, 0, 0, 2/
C ======================================================================
      iSYS(1) = 9
      iSYS(2) = 9

      LameM = DATA(1)
      LameL = DATA(2)

      Do i = 1, 9
         Do j = i, 9
            Coef(i, j) = LameM * im(i, j) 
            Coef(j, i) = Coef(i, j)
         End do
      End do

      Do i = 1, 9, 4
         Do j = 1, 9, 4
            Coef(i, j) = Coef(i, j) + LameL
         End do
      End do
      ANI_Delas_const = TENSOR_SYMMETRIC

      Return
      End

