c ======================================================================
      Program bioimpedance
c ======================================================================
c This program generates a finite element system for the diffusion problem
c
c  -div K grad u = 0     in  Omega
c              u = U_0   on  Gamma_D
c        K du/dn = G_1   on  Gamma_N=Gamma/GammaD
c
c where Omega = rassmatrivaemaya oblast
c    K(x)   - positive definite tensor
c    F(x)   - right-hand side
c    U_0(x) - essential (Dirichlet) boundary condition
c    G_0(x) - Neumann boundary condition
c
c These coefficients are implemented in the following routines:
c    K->Ddiff, F->Drhs,  {U_0,G_0,G_1}->Dbc
c
c ======================================================================
      implicit none

      include 'mmax.fd'
      integer namax1,nvmax1 ! dlya A1,JA1
      parameter(namax1=4*namax,nvmax1=2*nvmax)
c ... work memory
      Integer   MaxWi
c ... old value:
c ... Parameter(MaxWi = 800000000)
c ... new value:
      Parameter(MaxWi = 8000000)

      Integer  iW(MaxWi)
c ======================================================================
c Mesh definition
c ======================================================================
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb

c ... array of coordinates of mesh nodes
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for triangles and triangle labels
      Integer  tet(4,ntmax), material(ntmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(3,nbmax), labelF(nbmax)

c ======================================================================
c For library aniFEM
c ======================================================================
      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA1(nvmax), JA1(namax)
      Integer  IA2(nvmax), JA2(namax)

      Integer  status, nRow1, nCol1
      Integer  nRow2, nCol2
      Integer  nRow, nCol

      Real*8   A1(namax), RHS1(nvmax)
      Real*8   A2(namax), RHS2(nvmax)

      Real*8   DATAFEMR(44),DATAFEMI(44)
      Real*8   DATAFEME(1) ! dlya napryjonnosti E identical matrix

      Integer  IA(nvmax1), JA(namax1) !? dimension
      Real*8   A(namax1), RHS(nvmax1) !? dimension
      Integer  C2(nvmax1) !? dimesions
      Real*8   SOL(nvmax1), RES(nvmax1) ! for solving of system

c ... Integer  FEM3Dext
      EXTERNAL FEM3Dext

c ... local variables
      Integer  i, n,k,h,l,s,t
      Integer  iv,iv1,iv2,iv3
      Real*8   q1,q2,q3,q4,q5! dlini storon triangle, poluperimetr,area
      Real*8   Stotal(12),Stotal1, Stotal2 ! total area pod elektrodami
      Real*8   MesR(12),MesI(12) !srednii potencial
      Real*8   SilaToka, pr,pi !silatoka, pr,pi - naprajeniya (rel,imaginary parts)
      Real*8   ImpR, ImpI ! impedance real, imaginary parts

c ... dlya gradienta
      Real*8 GV(4,4), SV(4), RV1(3),SV2(4),GV2(4,4),RV2(3)! bilinear form, uzlovoe reshenie, gradient
      Real*8 Er1(3,ntmax), Ei1(3,ntmax) ! gradinet po vsem elemantam
      Real*8 Jr1(3,ntmax), Ji1(3,ntmax) ! plotnost toka dlya vseh elemtov
      Real*8 Er2(3,ntmax), Ei2(3,ntmax) ! gradinet po vsem elemantam
      Real*8 Jr2(3,ntmax), Ji2(3,ntmax) ! plotnost toka dlya vseh elemtov
      Real*8 XY1(3), XY2(3), XY3(3), XY4(3)! coordinati uzlov elementa
      Real*8 SnsI(ntmax), SnsR(ntmax) ! sensitivity function
      Real*8 SnsM(ntmax) ! remapped sensitivity


      External Ddiff2
      External Ddiff
      External Ddiff1
      External matrRotate
      External Sens
      External Jeval
      External Mblock

      Integer LSolver

      Integer Ddiff,nr1,nc1,LDA,Ddiff1
      Integer av,bv,cv,dv! nomera uzlov tetraedra
      Integer iSYS(5)

      Real*8 qw1,qw2,qw3,qw4,qw5,qw6,qw7,qw8,qw9,qw
      Real*8 SenI,SenR
      Real*8 J1,J2,J3,J4,J5,J6, J7,J8
      Real*8 Escl(ntmax),Jscl(ntmax) !scalar U,E,J,Sens

      Integer Num, NI1,NI2,NU1,NU2 ! chislo datchikov
      Integer chE,chJ,chGmv,chVtk

C      Real*8 SNX(47,75), SNY(21,75) ! for mesh xOz
      Integer  INFO
      Integer  j
      Real*8   st, sm(44)

      Integer LSmeth, LSprec ! linear solver method and preconditioner
      Logical LSFEX          ! for file with specific linear solver params

      Integer DPid
      Real*8  DPdist, Pdist
      Real*8  PT(3)


c =====================================================================
c ... check and load linear solver parameters
      inquire(file="solver.txt",exist=LSFEX) ! check the presence of solver.txt

      if (LSFEX) then                        ! solver.txt exists, load it 
	open(12, file='solver.txt')
	read(12,*) LSmeth,LSprec
	close(12)
      else                                   ! otherwise use default settings
        LSmeth = 0
        LSprec = 0
      end if

c =====================================================================
c ... load the initial mesh
      Call loadbin(
     &      nvmax, nbmax, ntmax,
     &      nv, nb, nt,
     &      vrt, bnd, tet, labelF, material,
     &      "mesh.ani*")



c     mark the Dirichlet point near target point PT

      PT(1) = 0D0
      PT(2) = 0D0
      PT(3) = 0D0
      Do n = 1, nv
         labelP(n) = 0
         PT(1) = PT(1) + vrt(1,n)
         PT(2) = PT(2) + vrt(2,n)
         PT(3) = PT(3) + vrt(3,n)
      End do
      PT(1) = PT(1) / nv
      PT(2) = PT(2) / nv
      PT(3) = PT(3) / nv
      DPid = 0
      DPdist = 0D0

      Do n = 1, nv
         Pdist = sqrt((vrt(1,n)-PT(1))**2+(vrt(2,n)-PT(2))**2+
     &                (vrt(3,n)-PT(3))**2)
         If(DPid.EQ.0 .OR. Pdist.LT.DPdist) then
           DPid = n
           DPdist = Pdist
         End if
      End do
c      labelP(DPid) = 111
c      write(*,*) 'Fix point: ', DPid, ', dist = ', DPdist
c      write(*,*) vrt(1,DPid), vrt(2,DPid), vrt(3,DPid)


c     general sparse matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)
100      status = IOR(MATRIX_GENERAL, FORMAT_AMG)


c     user data
      open(12, file='param.txt')
      read(12,*) Num, NI1,NI2,NU1,NU2,J7,J8
      read(12,*) (DATAFEMR(i),i=5,40)
      read(12,*) (DATAFEMI(i),i=5,40)
c      read(12,*) J7,J8
      read(12,*) chE,chJ,chGmv,chVtk
      close(12)

      DATAFEMR(1)=NI1+1
      DATAFEMR(2)=NI2+1
      DATAFEMI(1)=NI1+1
      DATAFEMI(2)=NI2+1

      DATAFEME(1)=0D0

c for Neuman conditions vicheslenie ploshadi dlya zond.toka
      Stotal1=0D0
      Stotal2=0D0
c nado, kone4no, vinesti v otdelnuu proceduru dlya ploshadi
      Do n = 1, nb
      If(labelF(n).EQ.DATAFEMR(1)) then
            iv1=bnd(1,n)
            iv2=bnd(2,n)
            iv3=bnd(3,n)

        q1=(vrt(1,iv1)-vrt(1,iv2))**2+
     &     (vrt(2,iv1)-vrt(2,iv2))**2+
     &     (vrt(3,iv1)-vrt(3,iv2))**2

        q2=(vrt(1,iv2)-vrt(1,iv3))**2+
     &     (vrt(2,iv2)-vrt(2,iv3))**2+
     &     (vrt(3,iv2)-vrt(3,iv3))**2

        q3=(vrt(1,iv3)-vrt(1,iv1))**2+
     &     (vrt(2,iv3)-vrt(2,iv1))**2+
     &     (vrt(3,iv3)-vrt(3,iv1))**2

        q4=(sqrt(q1)+sqrt(q2)+sqrt(q3))/2
        q5=sqrt(q4*(q4-sqrt(q1))*(q4-sqrt(q2))*(q4-sqrt(q3)))
        Stotal1=Stotal1+q5
     		
		else If(labelF(n).EQ.DATAFEMR(2)) then
            iv1=bnd(1,n)
            iv2=bnd(2,n)
            iv3=bnd(3,n)

        q1=(vrt(1,iv1)-vrt(1,iv2))**2+
     &     (vrt(2,iv1)-vrt(2,iv2))**2+
     &     (vrt(3,iv1)-vrt(3,iv2))**2

        q2=(vrt(1,iv2)-vrt(1,iv3))**2+
     &     (vrt(2,iv2)-vrt(2,iv3))**2+
     &     (vrt(3,iv2)-vrt(3,iv3))**2

        q3=(vrt(1,iv3)-vrt(1,iv1))**2+
     &     (vrt(2,iv3)-vrt(2,iv1))**2+
     &     (vrt(3,iv3)-vrt(3,iv1))**2

        q4=(sqrt(q1)+sqrt(q2)+sqrt(q3))/2
        q5=sqrt(q4*(q4-sqrt(q1))*(q4-sqrt(q2))*(q4-sqrt(q3)))
        Stotal2=Stotal2+q5
		End if		

	 End do
    
      DATAFEMR(3)=J7/Stotal1
      DATAFEMI(3)=J8/Stotal1
     
      DATAFEMR(4)=J7/Stotal2
      DATAFEMI(4)=J8/Stotal2

c       for  real part

         Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     FEM3Dext, DATAFEMR, status,
     &     nvmax, namax, IA1, JA1, A1, RHS1, nRow1, nCol1,
     &     MaxWi, iW)

c       for  imaginary part

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     FEM3Dext, DATAFEMI, status,
     &     nvmax, namax, IA2, JA2, A2, RHS2, nRow2, nCol2,
     &     MaxWi, iW)


c ============================================
c zapolnyaem blokami matricu. proverit gde kakie bloki
c po idee nrow1=nrow2
c ============================================
       write(*,*) 'constitute the matrix'

       call Mblock(IA1,IA2,A1,A2,JA1,JA2,
     & RHS1,RHS2,nRow1,IA,A,JA,RHS)

c      nCol=nRow

C ================================================================
C  STAGE 2: solve linear system
C ================================================================
      nRow=2*nRow1

      INFO = LSolver(LSmeth,LSprec, nRow,IA,JA,A,RHS,SOL)

C      if (INFO.ne.0) stop 'Linear solver failed'

c      call GMVscalarVrtMat(SOL,"mesh.gmv", 10,
c     &  nv, vrt, nt, tet, material, nb, bnd, labelF)

      Do i=1,Num
       Stotal(i)=0D0
       MesR(i)=0D0
       MesI(i)=0D0
      End do

      Do n = 1, nb

      If(labelF(n).GE.2) then
            i=labelF(n)-1

            iv1=bnd(1,n)
            iv2=bnd(2,n)
            iv3=bnd(3,n)

        q1=(vrt(1,iv1)-vrt(1,iv2))**2+
     &     (vrt(2,iv1)-vrt(2,iv2))**2+
     &     (vrt(3,iv1)-vrt(3,iv2))**2

        q2=(vrt(1,iv2)-vrt(1,iv3))**2+
     &     (vrt(2,iv2)-vrt(2,iv3))**2+
     &     (vrt(3,iv2)-vrt(3,iv3))**2

        q3=(vrt(1,iv3)-vrt(1,iv1))**2+
     &     (vrt(2,iv3)-vrt(2,iv1))**2+
     &     (vrt(3,iv3)-vrt(3,iv1))**2

        q4=(sqrt(q1)+sqrt(q2)+sqrt(q3))/2
        q5=sqrt(q4*(q4-sqrt(q1))*(q4-sqrt(q2))*(q4-sqrt(q3)))
        Stotal(i)=Stotal(i)+q5

        MesR(i)=MesR(i)+(SOL(iv1)+SOL(iv2)+SOL(iv3))*q5/3
        MesI(i)=MesI(i)+(SOL(iv1+nRow1)+SOL(iv2+nRow1)+
     &        SOL(iv3+nRow1))*q5/3

      End if

      End do

      Do i=1,Num

      MesR(i)=MesR(i)/Stotal(i)
      MesI(i)=MesI(i)/Stotal(i)

      End do


c ... raschet complex napejenii(pr,pi)=raznica srednih potencialov pod elktrodami
c ... impedancov (real and imaginary part)

      Open(1,file='napryajeniya.txt')
      Open(2,file='impedance.txt')
      SilaToka= (J7**2+J8**2)
      Do n=1, Num

      Do i=1, Num


      if (n.LT.i) then
      pr=MesR(n)-MesR(i)
      pi=MesI(n)-MesI(i)

c ... save the naprejenia

        Write(1,*)  n, i, pr, pi


c ... impedance/absolut values

      ImpR=dabs((pr*J7+pi*J8))/SilaToka
      ImpI=dabs((pi*J7-pr*J8))/SilaToka

        Write(2,*) n, i, ImpR, ImpI
      End if

      End do

      End do

      Close(1)
      Close(2)

      write(*,*) 'file nanapryajeniya.txt is completed'
      write(*,*) 'file impedance.txt is completed'
      write(*,*) 'E and J evaluating'
c ... vycheslenie gradienta

      Do n = 1,nt

      av=tet(1,n)
      bv=tet(2,n)
      cv=tet(3,n)
      dv=tet(4,n)

      XY1(1)=vrt(1,av)
      XY1(2)=vrt(2,av)
      XY1(3)=vrt(3,av)

      XY2(1)=vrt(1,bv)
      XY2(2)=vrt(2,bv)
      XY2(3)=vrt(3,bv)

      XY3(1)=vrt(1,cv)
      XY3(2)=vrt(2,cv)
      XY3(3)=vrt(3,cv)

      XY4(1)=vrt(1,dv)
      XY4(2)=vrt(2,dv)
      XY4(3)=vrt(3,dv)

      qw1=XY1(1)-XY4(1)
      qw2=XY1(2)-XY4(2)
      qw3=XY1(3)-XY4(3)

      qw4=XY2(1)-XY4(1)
      qw5=XY2(2)-XY4(2)
      qw6=XY2(3)-XY4(3)

      qw7=XY3(1)-XY4(1)
      qw8=XY3(2)-XY4(2)
      qw9=XY3(3)-XY4(3)

      qw=dabs(qw1*qw5*qw9+qw7*qw2*qw6+qw3*qw4*qw8-
     & qw3*qw5*qw7-qw2*qw4*qw9-qw1*qw8*qw6)/6

c... real part E

      SV(1)=SOL(av)/qw
      SV(2)=SOL(bv)/qw
      SV(3)=SOL(cv)/qw
      SV(4)=SOL(dv)/qw

      Call fem3dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P0vector,
     &              material, Ddiff1, DATAFEME, iSYS, 2,
     &              4, GV, nr1, nc1)


      Call dGemv('N',nr1,nc1,1D0,GV,4,SV,1,0d0,RV1,1)

      Er1(1,n)=RV1(1)
      Er1(2,n)=RV1(2)
      Er1(3,n)=RV1(3)

c... imaginary part E

      SV2(1)=SOL(av+nrow1)/qw
      SV2(2)=SOL(bv+nrow1)/qw
      SV2(3)=SOL(cv+nrow1)/qw
      SV2(4)=SOL(dv+nrow1)/qw

      Call dGemv('N',nr1,nc1,1D0,GV,4,SV2,1,0d0,RV2,1)

      Ei1(1,n)=RV2(1)
      Ei1(2,n)=RV2(2)
      Ei1(3,n)=RV2(3)

      Escl(n)=sqrt(Er1(1,n)**2+Er1(2,n)**2+Er1(3,n)**2+
     & Ei1(1,n)**2+Ei1(2,n)**2+Ei1(3,n)**2)

c ... plotnost toka v zavisimosti ot nomera materiala

      Call Jeval(XY1,XY2,XY3,XY4,material(n),DATAFEMR,DATAFEMI,
     & RV1,RV2,J1,J2,J3,J4,J5,J6)

      Jr1(1,n)=J1
      Jr1(2,n)=J2
      Jr1(3,n)=J3

      Ji1(1,n)=J4
      Ji1(2,n)=J5
      Ji1(3,n)=J6

      Jscl(n)=sqrt(Jr1(1,n)**2+Jr1(2,n)**2+Jr1(3,n)**2+
     & Ji1(1,n)**2+Ji1(2,n)**2+Ji1(3,n)**2)

      End do


c ... save the matrix

      if (chE.EQ.1 .and. chJ.EQ.1) then

          Open(1,file='E.txt')
          Open(2,file='J.txt')
          Do n=1,nt
              Write(1,*)  n
              Write(1,*) (Er1(i,n), i=1,3)
              Write(1,*) (Ei1(i,n), i=1,3)
              Write(1,*)  Escl(n)
              write(2,*)  n
              Write(2,*) (Jr1(i,n), i=1,3)
              Write(2,*) (Ji1(i,n), i=1,3)
              Write(2,*)   Jscl(n)
          End do
          Close(1)
          Close(2) 
          write(*,*) 'files E.txt and J.txt are completed'

      Else if (chE.EQ.1 .and. chJ.EQ.0) then
          Open(1,file='E.txt')
          Do n=1,nt
              Write(1,*)  n
              Write(1,*) (Er1(i,n), i=1,3)
              Write(1,*) (Ei1(i,n), i=1,3)
              Write(1,*)  Escl(n)
          End do
          Close(1)
          write(*,*) 'files E.txt is completed'
      Else if (chE.EQ.0 .and. chJ.EQ.1) then
          Open(2,file='J.txt')
          Do n=1,nt
              write(2,*)  n
              Write(2,*) (Jr1(i,n), i=1,3)
              Write(2,*) (Ji1(i,n), i=1,3)
              Write(2,*)   Jscl(n)
          End do
          Close(2) 
          write(*,*) 'files J.txt is completed'
      End if




      if (chGmv.EQ.1) then

          call GMVvectorTet(Er1,"er.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

          call GMVvectorTet(Ei1,"ei.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

          call GMVvectorTet(Jr1,"jr.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

          call GMVvectorTet(Ji1,"ji.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

          call GMVscalarTet(Escl,"escl.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

          call GMVscalarTet(Jscl,"jscl.gmv", 10,
     &  nv, vrt, nt, tet, nb, bnd, labelF)

      End if

      call savebin(6, "data1:ure.uim.ere.eim.jre.jim.", nv, 1, SOL,
     & nv, 1, SOL(nv), nt, 3, Er1, nt, 3, Ei1, nt, 3, Jr1, nt, 3, Ji1)

      if (chVtk.EQ.1) then
c     save VTK file for ParaView
        call VTKSpecial("mesh.vtk", 10,  "output",
     &  nv, vrt, nt, tet, material, nb, bnd, labelF,
     &  SOL, "U",
     &  SOL, "U_Re",
     &  SOL, "U_Im",
     &  Er1,  "E_Re",
     &  Ei1,  "E_Im",
     &  Jr1,  "J_Re",
     &  Ji1,  "J_Im",
     &  SOL, "X_Re",
     &  SOL, "X_Im")

      End if


c for potential electrodes

      DATAFEMR(1)=NU1+1
      DATAFEMR(2)=NU2+1
      DATAFEMI(1)=NU1+1
      DATAFEMI(2)=NU2+1

c for Neuman conditions vicheslenie ploshadi dlya zond.toka
      Stotal1=0D0
      Stotal2=0D0

      Do n = 1, nb
      If(labelF(n).EQ.DATAFEMR(1)) then
            iv1=bnd(1,n)
            iv2=bnd(2,n)
            iv3=bnd(3,n)

        q1=(vrt(1,iv1)-vrt(1,iv2))**2+
     &     (vrt(2,iv1)-vrt(2,iv2))**2+
     &     (vrt(3,iv1)-vrt(3,iv2))**2

        q2=(vrt(1,iv2)-vrt(1,iv3))**2+
     &     (vrt(2,iv2)-vrt(2,iv3))**2+
     &     (vrt(3,iv2)-vrt(3,iv3))**2

        q3=(vrt(1,iv3)-vrt(1,iv1))**2+
     &     (vrt(2,iv3)-vrt(2,iv1))**2+
     &     (vrt(3,iv3)-vrt(3,iv1))**2

        q4=(sqrt(q1)+sqrt(q2)+sqrt(q3))/2
        q5=sqrt(q4*(q4-sqrt(q1))*(q4-sqrt(q2))*(q4-sqrt(q3)))
        Stotal1=Stotal1+q5

      else If(labelF(n).EQ.DATAFEMR(2)) then
            iv1=bnd(1,n)
            iv2=bnd(2,n)
            iv3=bnd(3,n)

        q1=(vrt(1,iv1)-vrt(1,iv2))**2+
     &     (vrt(2,iv1)-vrt(2,iv2))**2+
     &     (vrt(3,iv1)-vrt(3,iv2))**2

        q2=(vrt(1,iv2)-vrt(1,iv3))**2+
     &     (vrt(2,iv2)-vrt(2,iv3))**2+
     &     (vrt(3,iv2)-vrt(3,iv3))**2

        q3=(vrt(1,iv3)-vrt(1,iv1))**2+
     &     (vrt(2,iv3)-vrt(2,iv1))**2+
     &     (vrt(3,iv3)-vrt(3,iv1))**2

        q4=(sqrt(q1)+sqrt(q2)+sqrt(q3))/2
        q5=sqrt(q4*(q4-sqrt(q1))*(q4-sqrt(q2))*(q4-sqrt(q3)))
        Stotal2=Stotal2+q5
                End if
      End do


      DATAFEMR(3)=J7/Stotal1
      DATAFEMI(3)=J8/Stotal1

      DATAFEMR(4)=J7/Stotal2
      DATAFEMI(4)=J8/Stotal2

      status = IOR(MATRIX_GENERAL, FORMAT_AMG)
c       for  real part

         Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     FEM3Dext, DATAFEMR, status,
     &     nvmax, namax, IA1, JA1, A1, RHS1, nRow1, nCol1,
     &     MaxWi, iW)

c       for  imaginary part

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     FEM3Dext, DATAFEMI, status,
     &     nvmax, namax, IA2, JA2, A2, RHS2, nRow2, nCol2,
     &     MaxWi, iW)


c ============================================
c zapolnyaem blokami matricu
c ============================================

       call Mblock(IA1,IA2,A1,A2,JA1,JA2,
     & RHS1,RHS2,nRow1,IA,A,JA,RHS)

C ================================================================
C  STAGE 2: solve linear system
C ================================================================
      nRow=2*nRow1

      INFO = LSolver(LSmeth,LSprec, nRow,IA,JA,A,RHS,SOL)

C      if (INFO.ne.0) stop 'Linear solver failed'

c      call GMVscalarVrtMat(SOL,"mesh.gmv", 10,
c     &  nv, vrt, nt, tet, material, nb, bnd, labelF)




c ... vycheslenie gradienta

      Do n = 1,nt

      av=tet(1,n)
      bv=tet(2,n)
      cv=tet(3,n)
      dv=tet(4,n)

      XY1(1)=vrt(1,av)
      XY1(2)=vrt(2,av)
      XY1(3)=vrt(3,av)

      XY2(1)=vrt(1,bv)
      XY2(2)=vrt(2,bv)
      XY2(3)=vrt(3,bv)

      XY3(1)=vrt(1,cv)
      XY3(2)=vrt(2,cv)
      XY3(3)=vrt(3,cv)

      XY4(1)=vrt(1,dv)
      XY4(2)=vrt(2,dv)
      XY4(3)=vrt(3,dv)

      qw1=XY1(1)-XY4(1)
      qw2=XY1(2)-XY4(2)
      qw3=XY1(3)-XY4(3)

      qw4=XY2(1)-XY4(1)
      qw5=XY2(2)-XY4(2)
      qw6=XY2(3)-XY4(3)

      qw7=XY3(1)-XY4(1)
      qw8=XY3(2)-XY4(2)
      qw9=XY3(3)-XY4(3)

      qw=dabs(qw1*qw5*qw9+qw7*qw2*qw6+qw3*qw4*qw8-
     & qw3*qw5*qw7-qw2*qw4*qw9-qw1*qw8*qw6)/6

c... real part E

      SV(1)=SOL(av)/qw
      SV(2)=SOL(bv)/qw
      SV(3)=SOL(cv)/qw
      SV(4)=SOL(dv)/qw

      Call fem3dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P0vector,
     &              material, Ddiff1, DATAFEME, iSYS, 2,
     &              4, GV, nr1, nc1)

      Call dGemv('N',nr1,nc1,1D0,GV,4,SV,1,0d0,RV1,1)

      Er2(1,n)=RV1(1)
      Er2(2,n)=RV1(2)
      Er2(3,n)=RV1(3)

c... imaginary part E

      SV2(1)=SOL(av+nrow1)/qw
      SV2(2)=SOL(bv+nrow1)/qw
      SV2(3)=SOL(cv+nrow1)/qw
      SV2(4)=SOL(dv+nrow1)/qw

      Call fem3dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P0vector,
     &              material, Ddiff1, DATAFEME, iSYS, 2,
     &              4, GV, nr1, nc1)

      Call dGemv('N',nr1,nc1,1D0,GV,4,SV2,1,0d0,RV2,1)

      Ei2(1,n)=RV2(1)
      Ei2(2,n)=RV2(2)
      Ei2(3,n)=RV2(3)

c ... plotnost toka v zavisimosti ot nomera materiala

      Call Jeval(XY1,XY2,XY3,XY4,material(n),DATAFEMR,DATAFEMI,
     & RV1,RV2,J1,J2,J3,J4,J5,J6)

      Jr2(1,n)=J1
      Jr2(2,n)=J2
      Jr2(3,n)=J3

      Ji2(1,n)=J4
      Ji2(2,n)=J5
      Ji2(3,n)=J6

      End do

c....funkciya chuvstvitelnosti

C      Do n=1,47
C      Do i=1,75
C      SNX(n,i)=0d0
C      End do
C      End do

C      Do n=1,41
C      Do i=1,75
C      SNY(n,i)=0d0
C      End do
C      End do

      st = 0d0
      Do n=1,44
        sm(n) = 0d0
      End do

C     save sensitivity for further use
c      Open(1,file='sensitivity.txt')

      Do n=1,nt

      Call Sens(n,Jr1,Ji1,Jr2,Ji2,SenR,SenI)

      SnsR(n)=SenR*1e12/J7/J7
      SnsI(n)=SenI*1e12/J7/J7

c      if (material(n).NE.2 .and. material(n).NE.3) then
      av=tet(1,n)
      bv=tet(2,n)
      cv=tet(3,n)
      dv=tet(4,n)

      XY1(1)=vrt(1,av)
      XY1(2)=vrt(2,av)
      XY1(3)=vrt(3,av)

      XY2(1)=vrt(1,bv)
      XY2(2)=vrt(2,bv)
      XY2(3)=vrt(3,bv)

      XY3(1)=vrt(1,cv)
      XY3(2)=vrt(2,cv)
      XY3(3)=vrt(3,cv)

      XY4(1)=vrt(1,dv)
      XY4(2)=vrt(2,dv)
      XY4(3)=vrt(3,dv)

      qw1=XY1(1)-XY4(1)
      qw2=XY1(2)-XY4(2)
      qw3=XY1(3)-XY4(3)

      qw4=XY2(1)-XY4(1)
      qw5=XY2(2)-XY4(2)
      qw6=XY2(3)-XY4(3)

      qw7=XY3(1)-XY4(1)
      qw8=XY3(2)-XY4(2)
      qw9=XY3(3)-XY4(3)

      qw=dabs(qw1*qw5*qw9+qw7*qw2*qw6+qw3*qw4*qw8-
     & qw3*qw5*qw7-qw2*qw4*qw9-qw1*qw8*qw6)/6

      qw1=(XY1(1)+XY2(1)+XY3(1)+XY4(1))/4
      qw2=(XY1(2)+XY2(2)+XY3(2)+XY4(2))/4
      qw3=(XY1(3)+XY2(3)+XY3(3)+XY4(3))/4

C      h=(qw1+5)/10+24
C      s=(qw2+5)/10+11
C      k=(qw3+5)/10+38
      
C      SNX(h,k)=SNX(h,k)+qw*SnsR(n)
C      SNY(s,k)=SNY(s,k)+qw*SnsR(n)
      st = st + qw*SnsR(n)
      sm(material(n)) = sm(material(n)) + qw*SnsR(n)
      call addsens(SnsR(n), material(n), qw, qw1, qw2, qw3)
c      Write(1,*) SnsR(n), material(n), qw, qw1, qw2, qw3

c      End if

      End do
c      Close(1)

      Open(1,file='sensitivity.txt')
      Write(1,*) 'Total sens = ', st
      Do n=1, 32
      Write(1,*)'Label',n,' = ',sm(n),'(',sm(n)/st*100d0,'%)'
      End do
      Close(1)
      call sensremap(nt, SnsR, SnsM)

C      Open(1,file='sens-xz.vtk')
C      Write(1,'(A)') '# vtk DataFile Version 2.0'
C      Write(1,'(A)') 'Sensitivity, plane XZ'
C      Write(1,'(A)') 'ASCII'
C      Write(1,'(A)') 'DATASET RECTILINEAR_GRID'
C      Write(1,'(A)') 'DIMENSIONS 47 1 75'
C      Write(1,'(A)') 'X_COORDINATES 47 float'
C      Write(1,*) (i*10,i=-23,23)
C      Write(1,'(A)') 'Y_COORDINATES 1 float'
C      Write(1,'(A)') '0'
C      Write(1,'(A)') 'Z_COORDINATES 75 float'
C      Write(1,*) (i*10,i=-37,37)
C      Write(1,'(A)') 'POINT_DATA 3525'
C      Write(1,'(A)') 'SCALARS scalars float'
C      Write(1,'(A)') 'LOOKUP_TABLE default'
C      Do i=1,75
C        Write(1,*)  (SNX(n,i),n=1,47)
C      End do
C      Close(1)
C      
C      Open(1,file='sens-yz.vtk')
C      Write(1,'(A)') '# vtk DataFile Version 2.0'
C      Write(1,'(A)') 'Sensitivity, plane XZ'
C      Write(1,'(A)') 'ASCII'
C      Write(1,'(A)') 'DATASET RECTILINEAR_GRID'
C      Write(1,'(A)') 'DIMENSIONS 1 21 75'
C      Write(1,'(A)') 'X_COORDINATES 1 float'
C      Write(1,'(A)') '0'
C      Write(1,'(A)') 'Y_COORDINATES 21 float'
C      Write(1,*) (i*10,i=-10,10)
C      Write(1,'(A)') 'Z_COORDINATES 75 float'
C      Write(1,*) (i*10,i=-37,37)
C      Write(1,'(A)') 'POINT_DATA 1575'
C      Write(1,'(A)') 'SCALARS scalars float'
C      Write(1,'(A)') 'LOOKUP_TABLE default'
C      Do i=1,75
C         Write(1,*)  (SNY(n,i),n=1,21)
C      End do
C      Close(1)
      

      call savebin(6, "data2:ure.uim.ere.eim.jre.jim.", nv, 1, SOL,
     & nv, 1, SOL(nv), nt, 3, Er2, nt, 3, Ei2, nt, 3, Jr2, nt, 3, Ji2)
      call savebin(3, "sens:sre.sim.smp.",
     & nt, 1, SnsR, nt, 1, SnsI, nt, 1, SnsM)
      if (chVtk.EQ.1) then
c     save VTK file for ParaView
        call VTKSpecial("mesh2.vtk", 10,  "output",
     &  nv, vrt, nt, tet, material, nb, bnd, labelF,
     &  SOL, "U",
     &  SOL, "U_Re",
     &  SOL, "U_Im",
     &  Er2,  "E_Re",
     &  Ei2,  "E_Im",
     &  Jr2,  "J_Re",
     &  Ji2,  "J_Im",
     &  SnsR,  "S_Re",
     &  SnsI,  "S_Im")

      End if

      Stop
      End

