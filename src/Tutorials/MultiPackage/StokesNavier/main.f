c ======================================================================
      Program  main
c ======================================================================
c The program generates and solves a nonlinear finite element system for the
c Navier-Stokes problem:
c 
c  -div grad u  + (u.grad) u + grad p = 0   in Omega
c        div u                        = 0   in Omega
c 
c                                   u = u_0 on dOmega_1
c                                   u = 0   on dOmega_2
c                           du/dn - p = 0   on dOmega_3
c
c where Omega is a backstep domain,
c dOmega_1 is the side at x=0, dOmega_3 is the side at x=1,
c and dOmega_2 is the rest of the boundary. The non-homogeneous
c boundary condition is
c
c    u_0 = { 64*(y-0.5)*(1-y)*z*(1-z), 0, 0 }.
c 
c We use the P2 finite elements for the velocity u and P1
c finite elements for the pressure p, which are known to
c be stable. The discretization method results in a system of
c nonlinear algebraic equations with the saddle point structure.
c
c *** Contents ***
c    main.f      - Main program with calls routines from the libraries
c 
c    forlibfem.f - The user-prepared file for the library libfem3D (matrix
c                  and right-hand side generators). It provides routines
c                  for computing the Stokes matrix, the nonlinear function
c                  of Navier-Stokes residual  and the right-hand side.
c 
c    forlibinb.f - The user-prepared file for the library libinb (nonlinear
c                  solver based on Inexact Newton Backtracking). It provides
c                  routines for nonlinear function evaluation and the
c                  evaluation of preconditioner for the jacobian.
c
c Step 1: we generate a quasi-uniform mesh using library libmba3D.a
c Step 2: we generate the algebraic Stokes system using library libfem3D.a
c Step 3: we initialize the preconditioner by factorizing the Stokes matrix with library libLU.a
c Step 4: we solve the nonlinear Navier-Stokes system using
c         the Inexact Newton-Krylov Backtracking with library libINB.a                
c ======================================================================
      implicit none

C=====GLOBAL PARAMETERS======
      integer nvmax,ntmax,nbmax,niamax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary edges
c ... niamax- maximum number of dofs (equations)
c ... namax - maximum number of non-zero Stokes matrix entries
      parameter(nvmax  =  5 500,   ntmax = 6*nvmax, nbmax = 10 000)
      parameter(niamax = nvmax*34, namax = 40*niamax)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 10 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, tetrahedra, and boundary facess
      Integer  nv, nt, nb

c ... coordinates of mesh points, labels of points
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for triangles and triangle labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(3,nbmax),labelF(nbmax)


c ... additional mesh arrays are needed to call aniMetric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(12), nbfix, ibfix(1), ntfix, itfix(1) 


C=====For LIBRARY LIBMBA======

      Integer   MaxSkipE, MaxQItr, nEStar
      Parameter(MaxSkipE = 100, MaxQItr = 50 000)

      Real*8    Quality
      Parameter(Quality = 4D-1)

      Real*8    rQuality
      Logical   flagAuto
      Integer   iPrint, iERR

      Integer   ANI_Metric_Eucl
      External  ANI_Metric_Eucl


C=====For LIBRARY LIBFEM======

      include 'fem3Dtet.fd'
      include 'assemble.fd'
      Integer  IA(niamax), JA(namax)
      Real*8    A(namax), RHS(niamax), SOL(niamax), RES(niamax)
      Integer  status, nRow, nCol

      Integer  Dbc, DATAFEM(1)
      EXTERNAL fem3Dext, Dbc
  

C=====For LIBRARY LIBLU======

      Integer  symbolic(2), numeric(2), sys
      Real*8   control(20), infoUMF(90)

      Integer  IB(niamax), JB(namax)
      Real*8    B(namax)

C=====For LIBRARY LIBINB======

      Integer  IPREVEC,INFO(5)
      Real*8   RESID, STPTOL
      External fnlin, prevec


C=====LOCAL VARIABLES======

      Integer   iSYS(5) 
      Integer   i,j,k,n, iux,iuy,iuz, ip, iv1,iv2,iv3,lbbnd, ibc
      Real*8    x,y,z, eBC(3), rmax
      Integer   iIRE, inEP, iIEP, iEnd, ipR, LenWork, nnz, nr

c ======================================================================
c Step 1: load a mesh
      Call loadMani(
     &      nvmax, nbmax, ntmax,
     &      nv, nb, nt,
     &      vrt, bnd, tet, labelF, material,
     &      nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &      iW, iW, "../data/ramp.ani")


c ... generate a quasi-uniform mesh with 2000 tets starting from ramp.ani
      nEStar   = 2000    ! number of tets in generated mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status   = 1       ! forbid boundary tets (see aniMBA/status.fd)
      iPrint   = 1       ! low level of output information


      Call mbaAnalytic(
c group (M)
     &     nv, nvmax, nb, nbmax, nt, ntmax,
     &     vrt, bnd, tet, labelF, material,
     &     nEStar,
c group (Dev)
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     ANI_Metric_Eucl, Quality, rQuality,
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)

      If(iERR.GT.1000) Call errMes(iERR, 'main',
     &                     'unspecified error if aniAnalytic')


c Step 2: generate the finite element Stokes system with P2-P1 elements

c ... mark the Dirichlet points via labelP
      Do i = 1, nv
         labelP(i) = 0
      End do
c ... set labels for all Dirichlet faces
      Do n = 1, nb
         If(labelF(n).NE.5) Then
            Do i = 1, 3
               iv1 = bnd(i,n)
               labelP(iv1) = labelF(n)
            End do
         End if
      End do
c ... overwrite labels for inhomogeneous Dirichlet faces
      Do n = 1, nb
         If(labelF(n).EQ.3) Then
            Do i = 1, 3
               iv1 = bnd(i,n)
               labelP(iv1) = labelF(n)
            End do
         End if
      End do


c     data provided for the user subroutine Dconv
      Do i = 1, niamax
         rW(i) = 0D0
      End Do

c ... general sparce matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)
      status = IOR(MATRIX_GENERAL, FORMAT_CSR)

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3Dext, rW, status,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, iW)


c Step 3: solve the Stokes finite element system using the LU decomposition
      Call AMG2CSC(nRow,IA,JA,A, nCol,IB,JB,B)

c ... converting to the 0-based format
      Call ZeroBasedFormat(nCol, IB, JB)

c ... set up default control parameters & print only error messages
      Call umf4def(control)
      control(1) = 1

c ... pre-order and symbolic analysis
      Call umf4sym(nCol, nCol, IB,JB,B, symbolic,control,infoUMF)
      If(infoUMF(1) .lt. 0) Then
         Write(*,*) 'Error occurred in umf4sym: ', infoUMF(1)
         Stop
      End if

c ... numeric factorization
      Call umf4num(IB,JB,B, symbolic,numeric,control,infoUMF)
      If(infoUMF(1) .lt. 0) Then
         Write(*,*) 'Error occurred in umf4num: ', infoUMF(1)
         Stop
      End if

c ... free the symbolic analysis data
      Call umf4fsym(symbolic)

c ... solve Ax=b, without iterative refinement
      sys = 0
      Call umf4sol(sys, SOL, RHS, numeric, control, infoUMF)
      If(infoUMF(1) .lt. 0) Then
         Write(*,*) 'Error occurred in umf4sol: ', infoUMF(1)
         Stop
      End if
c  ... numeric factorization data will be used in the preconditioning in INB solver


c ... check the residual for Stokes
      Call mulAgen(nRow, IA, JA, A, SOL, RES)
      rmax = 0
      Do i = 1, nRow
         rmax = max(rmax, RES(i) - RHS(i))
      End do
      Write(*,'(A,E12.6)') 'LU:   Maximal norm of residual: ', rmax

c ... SOL is the solution of the Stokes problem

c ... prepare real*8 and integer data to pass to the INB solver

c ... IRE(6, nt) is a map from elements to edges; the number of edge is nr
c ... nEP and IEP are auxiliary arrays
      iIRE = 6
      inEP = iIRE + 6 * nt
      iIEP = inEP + nv
      iEnd = iIEP + 4 * nt

      If(iEnd.GT.MaxWi) Then
         Write(*,*) 'Size of the integer array iW is small'
         Stop
      End if

      Call listE2R(nv, nr, nt, tet, iW(iIRE), iW(inEP), iW(iIEP))

      k = iIRE + 6 * nt

c ... pass general mesh parameters to the internal function fnlin
      nnz = IA(nRow+1)-1
      iW(1) = nv
      iW(2) = nt
      iW(3) = nb
      iW(4) = nnz+1
      iW(5) = MaxWi - 5 - 11*nt - nv - 4*nb - (nRow+1) - (nnz+1)

      If(iW(5).le.0) Then
         Write(*,*) 'Size of the integer array iW is small'
         Stop
      End if

c ... pass labels of points
      Do i = 1, nv
         iW(k) = labelP(i)
         k = k + 1
      End do

c ... pass boundary faces
      Do i = 1, nb
         Do j = 1, 3
            iW(k) = bnd(j,i)
            k = k + 1
         End do
      End do

c ... pass labels of faces 
      Do i = 1, nb
         iW(k) = labelF(i)
         k = k + 1
      End do

c ... pass connectivity table for tetrahedra
      Do i = 1, nt
         Do j = 1, 4
            iW(k) = tet(j,i)
            k = k + 1
         End do
      End do

c ... pass labels of tetrahedra
      Do i = 1, nt
         iW(k) = material(i)
         k = k + 1
      End do

      LenWork = 11*nRow  ! work memory for INB
      ipR = LenWork + 1

      Call dcopy(3*nv,vrt,1,rW(ipR),1)

      If(LenWork + 3*nv + nRow + nnz .gt. MaxWr) then
         Write(*,*) 'Size of the Real*8 array rW is small'
         Stop
      End if

c ... call the driver for the solution by inexact Newton-Krylov backtracking
      RESID = 1d-6
      STPTOL = 1d-7
      IPREVEC = nRow
      INFO(1) = 0   ! initializing successful termination flag for Newton
      INFO(2) = 10  ! maximal number of linear iterations
      INFO(3) = 10  ! maximal number of nonlinear iterations
      INFO(4) = 10  ! maximal number of backtracks
      INFO(5) = 1   ! print level 1 (0 mum, 1 nonlinear residuals, 2 for linear residual)

      call slInexactNewton(
     &      prevec,IPREVEC,numeric,control,
     &      fnlin, rW(ipR), iW,
     &      nRow,SOL,
     &      RESID, STPTOL,
     &      rW,LenWork,
     &      INFO)

      If(INFO(1).ne.0) then
         Write(*,*) 'failed to solve the problem with',
     &              ' desired accuracy, INFO=',INFO(1)
         Stop
      End if

      Write(*,'(A,I3)') '     Number of linear iterations:   ',INFO(2)
      Write(*,'(A,I3)') '     Number of nonlinear iterations:',INFO(3)
      Write(*,'(A,I3)') '     Number of backtracks:          ',INFO(4)
      Write(*,'(A,I3)') '     Number of function evaluations:',INFO(5)

c ... free the numeric factorization data
      Call umf4fnum (numeric)


c Step 4: draw the solution in GMV file
c ... unknowns in SOL: nodal x-component of velocity starts at SOL(iux),
c                      nodal y-component of velocity starts at SOL(iuy),
c                      nodal z-component of velocity starts at SOL(iuz),
c                      nodal pressure starts at SOL(ip)

      iux = 1
      iuy = iux + nv
      iuz = iuy + nv
      ip  = iuz + nv

      Call GMVscalarVrt(SOL(iux), "velocity_x.gmv", 10,
     &                nv, vrt, nt, tet, nb, bnd, labelF)

      Call GMVscalarVrt(SOL(ip), "pressure.gmv", 10,
     &                nv, vrt, nt, tet, nb, bnd, labelF)

      Stop
      End



