C ======================================================================
C  The user defined routines required above
C ======================================================================
      Subroutine FEM3Dext (XY1, XY2, XY3, XY4,
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
      Real*8  DATAFEM1(6)
      Real*8  DATAFEM2(6)
      Real*8  DATAFEM3(6)

      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      External matrRotate

      Real*8   C(3, 3), G(3), XYP(3, 4)
      Real*8   x, y, z, eBC(1)

      Integer  i, j, k, l, m, ir, ic, label, ibc

      Integer  iref(5), ip(4)
      DATA     iref /1, 2, 3, 4, 1/


C ======================================================================
      nRow = 4
      nCol = 4

c ... set up templated degrees of freedom for rows and columns.
c     used convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
      Do i = 1, 4
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c anisotropy of material numbers 1,2,3
      if (label .EQ. 1 .or. 
     &    label .EQ. 2 .or.
     &    label .EQ. 3) then

        Call matrRotate(XY1,XY2,XY3,XY4,label,
     &                  DATAFEM,DATAFEM1)

c     A(1:4,1:4) is elemental vector elliptic operator;
c     in other words, for the bilinear form <grad(P1), grad(P1)>

        Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DATAFEM1, iSYS, 2,
     &              LDA, A, ir, ic)

c ... compute the right hand side vector using external function Drhs
c     in other words, the linear form <Drhs(x), P1>
        Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, DATAFEM1, iSYS, 5,
     &              LDA, F, ir, ic)


      Else

        Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DATAFEM, iSYS, 2,
     &              LDA, A, ir, ic)

        Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, DATAFEM, iSYS, 5,
     &              LDA, F, ir, ic)

      End if


c ... impose the Neumann boundary conditions
      Do i = 1, 3
        XYP(i, 1) = XY1(i)
        XYP(i, 2) = XY2(i)
        XYP(i, 3) = XY3(i)
        XYP(i, 4) = XY4(i)
      End do

      Do k = 1, 4
        If(lbF(k) .GT. 0) Then
            l = iref(k + 1)
            m = iref(l + 1)

            x = (XYP(1, k) + XYP(1, l) + XYP(1, m)) / 3
            y = (XYP(2, k) + XYP(2, l) + XYP(2, m)) / 3
            z = (XYP(3, k) + XYP(3, l) + XYP(3, m)) / 3

            ibc = Dbc(x, y, z, lbF(k), DATAFEM, iSYS, eBC)

            If(ibc .EQ. BC_NEUMANN) Then
               label = lbF(k)

               Call fem3Dtri(XYP(1, k), XYP(1, l), XYP(1, m),
     &                       IDEN, FEM_P0, IDEN, FEM_P1,
     &                       label, Dbc, DATAFEM, iSYS, 4,
     &                       3, G, ir, ic)

               F(k) = F(k) + G(1)
               F(l) = F(l) + G(2)
               F(m) = F(m) + G(3)
            End if


        End if
      End do




c ... impose Dirichlet boundary conditions at triangle nodes
c     this condition should go the last one
      Do k = 1, 4
         If(lbP(k) .NE. 0) Then
            x = XYP(1, k)
            y = XYP(2, k)
            z = XYP(3, k)!
            ibc = Dbc(x, y, z, lbP(k), DATAFEM, iSYS, eBC)

            If(ibc .EQ. BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End

C ======================================================================
c 3x3 diffusion tensor K
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, *)
      Integer label, iSYS(*)

      Integer i, j

      iSYS(1) = 3
      iSYS(2) = 3

      if (label .EQ. 1 .or. 
     &    label .EQ. 2 .or.
     &    label .EQ. 3) then

        Coef(1, 1) = DATA(1)
        Coef(2, 2) = DATA(4)
        Coef(3, 3) = DATA(6)

        Coef(1, 2) = DATA(2)
        Coef(2, 1) = DATA(2)

        Coef(1, 3) = DATA(3)
        Coef(3, 1) = DATA(3)

        Coef(3, 2) = DATA(5)
        Coef(2, 3) = DATA(5)

      Else

        Coef(1, 1) = DATA(label+7)
        Coef(2, 2) = DATA(label+7)
        Coef(3, 3) = DATA(label+7)

        Coef(1, 2) = 0D0
        Coef(2, 1) = 0D0

        Coef(1, 3) = 0D0
        Coef(3, 1) = 0D0

        Coef(3, 2) = 0D0
        Coef(2, 3) = 0D0
      End if

      Ddiff = TENSOR_SYMMETRIC

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

        If(label .EQ. 111) Then
           Coef(1) = 0! Dirichlet at point
           Dbc = BC_DIRICHLET
        Else If (label .EQ. 1) Then
           Coef(1) = 0
           Dbc = BC_NEUMANN     ! potok na poverhnosti tela
        Else If (label .EQ. DATA(1)) Then
           Coef(1) = DATA(3)!1D0
           Dbc = BC_NEUMANN
        Else If (label.EQ.DATA(2)) Then
           Coef(1) = -DATA(4)!-1D0
           Dbc = BC_NEUMANN
        Else If (label.GT.1 .AND. label.LT.30) Then
           Coef(1) = 0
           Dbc = BC_NEUMANN     ! potok na poverhnosti tela
        Else
           Dbc = BC_NULL
        End if

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

        iSYS(1) = 1
        iSYS(2) = 1

        Coef(1) = 0D0
        Drhs = TENSOR_SCALAR

      Return
      End

c==========================================================================
c povorot
c==========================================================================
      Subroutine matrRotate(XY1,XY2,XY3,XY4,label,DATA,DATAR)
      implicit none


      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      Real*8  DATA(*), DATAR(*)

C LOCAL VARIABLEs

      Integer i, j, k, l
      Integer label

      Real*8 e1, e2, e3, ed ! e1 (this is a comment)
      Real*8 m1, m2, m3, md ! e2
      Real*8 n1, n2, n3, nd ! e3
      Real*8 f1, f2, f3, fd !

      Real*8 KB(3,3), T(3,3), KF(3,3)
c====================================================
      ed = sqrt ( (XY2(1) - XY1(1)) ** 2 + 
     &            (XY2(2) - XY1(2)) ** 2 +
     &            (XY2(3) - XY1(3)) ** 2)
      e1 = (XY2(1) - XY1(1)) / ed
      e2 = (XY2(2) - XY1(2)) / ed
      e3 = (XY2(3) - XY1(3)) / ed


      fd = sqrt( (XY3(1) - XY1(1)) ** 2 + 
     &           (XY3(2) - XY1(2)) ** 2 +
     &           (XY3(3) - XY1(3)) ** 2)
      f1 = (XY3(1) - XY1(1)) / fd
      f2 = (XY3(2) - XY1(2)) / fd
      f3 = (XY3(3) - XY1(3)) / fd

      md = sqrt( (e2 * f3 - e3 * f2) ** 2 + 
     &           (e3 * f1 - e1 * f3) ** 2 +
     &           (e1 * f2 - e2 * f1) ** 2)
      m1 = (e2 * f3 - e3 * f2) / md
      m2 = (e3 * f1 - e1 * f3) / md
      m3 = (e1 * f2 - e2 * f1) / md

      nd = sqrt( (e2 * m3 - m2 * e3) ** 2 + 
     &           (e3 * m1 - e1 * m3) ** 2 +
     &           (e1 * m2 - e2 * m1) ** 2)
      n1 = (e2 * m3 - m2 * e3) / nd
      n2 = (e3 * m1 - e1 * m3) / nd
      n3 = (e1 * m2 - e2 * m1) / nd

c provodimost posle povorota k1
c zapolnyeam matrix, posle povorota
c matrix povorota T
      T(1, 1) = e1
      T(2, 1) = e2
      T(3, 1) = e3
      T(1, 2) = m1
      T(2, 2) = m2
      T(3, 2) = m3
      T(1, 3) = n1
      T(2, 3) = n2
      T(3, 3) = n3

      Do i = 1, 3
      	Do j = 1, 3
          KF(i, j) = 0D0
        End do
      End do

      if (label .EQ. 1) then
        KB(1, 1) = DATA(5)
        KB(2, 1) = 0D0     ! DATA(2)
        KB(3, 1) = 0D0     ! DATA(3)
        KB(1, 2) = 0D0     ! DATA(2)
        KB(2, 2) = DATA(6) ! DATA(4)
        KB(3, 2) = 0D0     ! DATA(5)
        KB(1, 3) = 0D0     ! DATA(3)
        KB(2, 3) = 0D0     ! DATA(5)
        KB(3, 3) = DATA(6)

      Else if (label .EQ. 2) then
        KB(1, 1) = DATA(7)
        KB(2, 1) = 0D0
        KB(3, 1) = 0D0
        KB(1, 2) = 0D0
        KB(2, 2) = DATA(8)
        KB(3, 2) = 0D0
        KB(1, 3) = 0D0
        KB(2, 3) = 0D0
        KB(3, 3) = DATA(8)

      Else if (label .EQ. 3) then
        KB(1, 1) = DATA(9)
        KB(2, 1) = 0D0
        KB(3, 1) = 0D0
        KB(1, 2) = 0D0
        KB(2, 2) = DATA(10)
        KB(3, 2) = 0D0
        KB(1, 3) = 0D0
        KB(2, 3) = 0D0
        KB(3, 3) = DATA(10)

      End if

c ::: KF += T * KB * T.transposed
      Do i = 1, 3
      	Do j = 1, 3
	        Do k = 1, 3
	          Do l = 1, 3
	            KF(i, j) = KF(i, j) + T(i, k) * T(j, l) * KB(k, l) 
	          End do
	        End do
       	End do
      End do

      DATAR(1) = KF(1, 1)
      DATAR(2) = KF(2, 1)
      DATAR(3) = KF(3, 1)
      DATAR(4) = KF(2, 2)
      DATAR(5) = KF(3, 2)
      DATAR(6) = KF(3, 3)

      Return
      End

c ==============================================================
      Integer Function Ddiff1(x, y, z, label, DATA, iSYS, Coef)
      implicit none
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DATA(*), Coef(9, *)
      Integer label, iSYS(*)

      Integer i, j

      iSYS(1) = 3
      iSYS(2) = 3

      Coef(1, 1) = 1D0
      Coef(2, 2) = 1D0
      Coef(3, 3) = 1D0

      Coef(1, 2) = 0D0
      Coef(2, 1) = 0D0

      Coef(1, 3) = 0D0
      Coef(3, 1) = 0D0

      Coef(3, 2) = 0D0
      Coef(2, 3) = 0D0

      Ddiff1 = TENSOR_SYMMETRIC

      Return
      End

c 3x3 diffusion tensor K for results
c ======================================================================
      Subroutine Ddiff2(label, DATA, Coef)
      implicit none


      Real*8  DATA(*), Coef(3, 3)
      Integer label

      Integer i, j

      if (label .EQ. 1 .or.
     &    label .EQ. 2 .or.
     &    label .EQ. 3) then

        Coef(1, 1) = DATA(1)
        Coef(2, 2) = DATA(4)
        Coef(3, 3) = DATA(6)

        Coef(1, 2) = DATA(2)
        Coef(2, 1) = DATA(2)

        Coef(1, 3) = DATA(3)
        Coef(3, 1) = DATA(3)

        Coef(3, 2) = DATA(5)
        Coef(2, 3) = DATA(5)

      Else

        Coef(1, 1) = DATA(label + 7)
        Coef(2, 2) = DATA(label + 7)
        Coef(3, 3) = DATA(label + 7)

        Coef(1, 2) = 0D0
        Coef(2, 1) = 0D0

        Coef(1, 3) = 0D0
        Coef(3, 1) = 0D0

        Coef(3, 2) = 0D0
        Coef(2, 3) = 0D0

      End if

      Return
      End

c=======================================================
c vicheslenie chuvstvitelnosti
c=======================================================
      Subroutine Sens(n, E11, E12, E21, E22, SnR, SnI)
      implicit none

      Integer n
      Real*8 E11(3, *), E12(3, *), E21(3, *), E22(3, *)
      Real*8 SnI, SnR

      SnR = E11(1, n) * E21(1, n) + 
     &      E11(2, n) * E21(2, n) + 
     &      E11(3, n) * E21(3, n) - 
     &      E12(1, n) * E22(1, n) - 
     &      E12(2, n) * E22(2, n) - 
     &      E12(3, n) * E22(3, n)

      SnI = E12(1, n) * E21(1, n) + 
     &      E12(2, n) * E21(2, n) + 
     &      E12(3, n) * E21(3, n) +
     &      E11(1, n) * E22(1, n) + 
     &      E11(2, n) * E22(2, n) + 
     &      E11(3, n) * E22(3, n)

      Return
      End

c==========================================================
c vycheslenie plotnosti toka
c==========================================================
      Subroutine Jeval(XY1,XY2,XY3,XY4,label,DATAFEMR,DATAFEMI,
     & RV1,RV2,J1,J2,J3,J4,J5,J6)

      include 'fem3Dtet.fd'

      Real*8  RV1(*),RV2(*)! napryajennost
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      Real*8  DATAFEMR(*),DATAFEMI(*)

      Real*8  DATAFEM1(6), DATAFEM2(6)


C LOCAL VARIABLEs

      Integer label

      Real*8 J1,J2,J3,J4,J5,J6

      Real*8 CoeffE(3,3), CoeffI(3,3), RV3(3),RV4(3)

      External Ddiff2
      External matrRotate

c =====================================

c ... plotnost toka v zavisimosti ot nomera materiala

      if (label.EQ.1 .or. label.EQ.2 .or.
     & label.EQ.3) then

      Call matrRotate(XY1,XY2,XY3,XY4,label,
     & DATAFEMR,DATAFEM1)

      Call Ddiff2(label,DATAFEM1,CoeffE)

      Call matrRotate(XY1,XY2,XY3,XY4,label,
     & DATAFEMI,DATAFEM2)

      Call Ddiff2(label,DATAFEM2,CoeffI)

c... real part J

      Call dGemv('N',3,3,1D0,CoeffE,3,RV1,1,0d0,RV3,1)

      Call dGemv('N',3,3,1D0,CoeffI,3,RV2,1,0d0,RV4,1)

       J1=RV3(1)-RV4(1)
       J2=RV3(2)-RV4(2)
       J3=RV3(3)-RV4(3)

c... imaginary part J

      Call dGemv('N',3,3,1D0,CoeffI,3,RV1,1,0d0,RV3,1)

      Call dGemv('N',3,3,1D0,CoeffE,3,RV2,1,0d0,RV4,1)

      J4=RV3(1)+RV4(1)
      J5=RV3(2)+RV4(2)
      J6=Rv3(3)+RV4(3)


      Else

      Call Ddiff2(label,DATAFEMR,CoeffE)

      Call Ddiff2(label,DATAFEMI,CoeffI)


c... real part J

      Call dGemv('N',3,3,1D0,CoeffE,3,RV1,1,0d0,RV3,1)

      Call dGemv('N',3,3,1D0,CoeffI,3,RV2,1,0d0,RV4,1)

      J1=RV3(1)-RV4(1)
      J2=RV3(2)-RV4(2)
      J3=RV3(3)-RV4(3)

c      Jr(1,n)=RV3(1)-RV4(1)
c      Jr(2,n)=RV3(2)-RV4(2)
c      Jr(3,n)=RV3(3)-RV4(3)

c... imaginary part J

      Call dGemv('N',3,3,1D0,CoeffI,3,RV1,1,0d0,RV3,1)

      Call dGemv('N',3,3,1D0,CoeffE,3,RV2,1,0d0,RV4,1)

      J4=RV3(1)+RV4(1)
      J5=RV3(2)+RV4(2)
      J6=RV3(3)+RV4(3)

      End if

      Return
      End

c=======================================================
c sostavlenie blochnoi matrix
c=======================================================
      Subroutine Mblock(IA1,IA2,A1,A2,JA1,JA2,
     & RHS1,RHS2,nrow1,IA,A,JA,RHS)

      implicit none
c==============================================
c Local variables
c==============================================
      Integer IA1(*),IA2(*),JA1(*),JA2(*)
      Integer IA(*),JA(*)
      Real*8 RHS1(*), RHS2(*),RHS(*)
      Real*8 A1(*),A2(*),A(*)

      Integer nrow1
      Integer s,t,h,l,i,k

      include 'mmax.fd'
      Integer nvmax1
      parameter(nvmax1=2*nvmax)

      Integer  C2(nvmax1)
c=============================================

      IA(1)=1
      C2(1)=0
      s=1
      t=1
      h=1
      l=1

      Do i=2, (nRow1+1)
      IA(i)=IA(i-1)+IA1(i)-IA1(i-1)+IA2(i)-IA2(i-1)
      C2(i)=IA1(i)-IA1(i-1)
      RHS(i-1)=RHS1(i-1)

       Do k=(IA(i-1)+C2(i)-1), IA(i-1),-1
          if ((IA1(i)-s). GT. 0) then
              A(k)=A1(IA1(i)-s)
              JA(k) =JA1(IA1(i)-s)
          End if
             s=s+1
       End do
           s=1

       Do k=(IA(i)-1),(IA(i-1)+C2(i)),-1
          if ((IA2(i)-t). GT. 0) then
             A(k)=-A2(IA2(i)-t)
             JA(k)=JA2(IA2(i)-t)+nRow1
          End if
             t=t+1
          End do
             t=1
      End do

      Do i=(nRow1+2), (2*nRow1+1)
        IA(i)=IA(i-1)+IA1(i-nRow1)-IA1(i-1-nRow1)+
     &  IA2(i-nRow1)-IA2(i-1-nRow1)
        C2(i)=IA2(i-nRow1)-IA2(i-1-nRow1)
        RHS(i-1)=RHS2(i-1-nrow1)
        Do k=(IA(i-1)+C2(i)-1), IA(i-1),-1

         if ((IA2(i-nRow1)-h). GT. 0) then
              A(k)=A2(IA2(i-nRow1)-h)
              JA(k) =JA2(IA2(i-nRow1)-h)
         End if
             h=h+1
        End do
           h=1

        Do k=(IA(i)-1),(IA(i-1)+C2(i)),-1
          if ((IA1(i-nRow1)-l). GT. 0) then
               A(k)=A1(IA1(i-nRow1)-l)
               JA(k)=JA1(IA1(i-nRow1)-l)+nRow1
          End if
             l=l+1
         End do
             l=1
      End do

      Return
      End



c=======================================================
c sostavlenie blochnoi matrix, when imaginary part first
c=======================================================
      Subroutine Mblock_Im(IA1,IA2,A1,A2,JA1,JA2,
     & RHS1,RHS2,nrow1,IA,A,JA,RHS)

      implicit none
c==============================================
c Local variables
c==============================================
      Integer IA1(*),IA2(*),JA1(*),JA2(*)
      Integer IA(*),JA(*)
      Real*8 RHS1(*), RHS2(*),RHS(*)
      Real*8 A1(*),A2(*),A(*)

      Integer nrow1
      Integer s,t,h,l,i,k

      include 'mmax.fd'
      Integer nvmax1
      parameter(nvmax1=2*nvmax)

      Integer  C2(nvmax1)
c=============================================

      IA(1)=1
      C2(1)=0
      s=1
      t=1
      h=1
      l=1

      Do i=2, (nRow1+1)
      IA(i)=IA(i-1)+IA2(i)-IA2(i-1)+IA1(i)-IA1(i-1)
      C2(i)=IA2(i)-IA2(i-1)
      RHS(i-1)=-RHS2(i-1)

       Do k=(IA(i-1)+C2(i)-1), IA(i-1),-1
          if ((IA2(i)-s). GT. 0) then
              A(k)=A2(IA2(i)-s)
              JA(k) =JA2(IA2(i)-s)
          End if
             s=s+1
       End do
           s=1

       Do k=(IA(i)-1),(IA(i-1)+C2(i)),-1
          if ((IA2(i)-t). GT. 0) then
             A(k)=-A1(IA1(i)-t)
             JA(k)=JA1(IA1(i)-t)+nRow1
          End if
             t=t+1
          End do
             t=1
      End do

      Do i=(nRow1+2), (2*nRow1+1)
        IA(i)=IA(i-1)+IA2(i-nRow1)-IA2(i-1-nRow1)+
     &  IA1(i-nRow1)-IA1(i-1-nRow1)
        C2(i)=IA1(i-nRow1)-IA1(i-1-nRow1)
        RHS(i-1)=RHS1(i-1-nrow1)
        Do k=(IA(i-1)+C2(i)-1), IA(i-1),-1

         if ((IA1(i-nRow1)-h). GT. 0) then
              A(k)=A1(IA1(i-nRow1)-h)
              JA(k) =JA1(IA1(i-nRow1)-h)
         End if
             h=h+1
        End do
           h=1

        Do k=(IA(i)-1),(IA(i-1)+C2(i)),-1
          if ((IA2(i-nRow1)-l). GT. 0) then
               A(k)=A2(IA2(i-nRow1)-l)
               JA(k)=JA2(IA2(i-nRow1)-l)+nRow1
          End if
             l=l+1
         End do
             l=1
      End do

      Return
      End


c============================================
c vycheslenie opredelitelya
c============================================
      Subroutine Det(qw1,qw2,qw3,qw4,qw5,qw6,qw7,qw8,qw9,qw)
      implicit none

      Real*8 qw1, qw2, qw3, qw4, qw5, qw6
      Real*8 qw7,qw8,qw9,qw

      qw=(qw1*qw5*qw9+qw7*qw2*qw6+qw3*qw4*qw8-
     & qw3*qw5*qw7-qw2*qw4*qw9-qw1*qw8*qw6)

      Return
      End


c============================================
c nahojdenie constant dlya ploskostei
c============================================
      Subroutine ConsPlane(XY1,XY2,XY3,A,B,C,in)
      implicit none

      Real*8 XY1(*), XY2(*), XY3(*)
      Real*8 A, B, C
      Real*8 D, D1, D2, D3
      Real*8 t1, t2, t3
      Integer in


      External Det

      t1 = XY1(1) ** 2 + XY1(2) ** 2 + XY1(3) ** 2
      t2 = XY2(1) ** 2 + XY2(2) ** 2 + XY2(3) ** 2
      t1 = XY3(1) ** 2 + XY3(2) ** 2 + XY3(3) ** 2

      if (t1.NE.0 .and. t2.NE.0 .and. t3.NE.0) then
      Call Det(XY1(1),XY1(2),XY1(3),XY2(1),XY2(2),
     & XY2(3),XY3(1),XY3(2),XY3(3),D)

      Call Det(-1D0,XY1(2),XY1(3),-1D0,XY2(2),
     & XY2(3),-1D0,XY3(2),XY3(3),D1)

      Call Det(XY1(1),-1D0,XY1(3),XY2(1),-1D0,
     & XY2(3),XY3(1),-1D0,XY3(3),D2)

      Call Det(XY1(1),XY1(2),-1D0,XY2(1),XY2(2),
     & -1D0,XY3(1),XY3(2),-1D0,D3)

      if (D.NE.0) then
        A = D1 / D
        B = D2 / D
        C = D3 / D
        in = 1
      else
        stop 'determinant=0'
      end if

      Else if (t1.EQ.0) then
      D=XY2(1)*XY3(2)-XY2(2)*XY3(1)
      D1=-(XY2(3)*XY3(2)-XY2(2)*XY3(3))
      D2=-(XY2(1)*XY3(3)-XY2(3)*XY3(1))

      if (D.NE.0) then
      A=D1/D
      B=D2/D
      C=1d0
      in=0
      else
      stop 'determinant=0'
      end if

      Else if (t2.EQ.0) then
      D=XY1(1)*XY3(2)-XY1(2)*XY3(1)
      D1=-(XY1(3)*XY3(2)-XY1(2)*XY3(3))
      D2=-(XY1(1)*XY3(3)-XY1(3)*XY3(1))

      if (D.NE.0) then
      A=D1/D
      B=D2/D
      C=1d0
      in=0
      else
      stop 'determinant=0'
      end if

      Else if (t3.EQ.0) then
      D=XY1(1)*XY2(2)-XY1(2)*XY2(1)
      D1=-(XY1(3)*XY2(2)-XY1(2)*XY2(3))
      D2=-(XY1(1)*XY2(3)-XY1(3)*XY2(1))

      if (D.NE.0) then
      A=D1/D
      B=D2/D
      C=1
      in=0
      else
      stop 'determinant=0'
      end if


      end if

      Return
      End
