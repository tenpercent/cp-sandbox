c ==========================================================
      Program  main
c ==========================================================
c This program generates a finite element system for the Stokes problem
c
c  -div grad u  + grad p = 0   in Omega
c        div u           = 0   in Omega
c
c                      u = u_0 on dOmega_1
c                      u = 0   on dOmega_2
c              du/dn - p = 0   on dOmega_3
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
c be the stable pair. The discretization method results in 
c a symmetric indefinite matrix.
c
c Step 1: we generate a quasi-uniform mesh using library libmba3D.a
c Step 2: we generate an algebraic system and matrices for preconditioner 
c         using library libfem3D.a
c Step 3: we solve the finite element system using library libILU.a
c ==========================================================
      implicit none


C=====GLOBAL PARAMETERS======

      integer nvmax,ntmax,nbmax,niamax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary edges
c ... niamax- maximum number of matrix rows
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax  =  5 500,   ntmax = 6*nvmax, nbmax = 10 000)
      parameter(niamax = nvmax*34, namax = 40*niamax)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 10 000 000, MaxWi = 10 000 000)

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

      Real*8   DATAFEM(1)
      Integer  Dbc
      EXTERNAL fem3Dext, Dbc
      EXTERNAL fem3DextMassP1, fem3DextLaplP2
  

C=====For LIBRARY LIBILU======

      EXTERNAL matvec, prevec
      Integer  imatvec(1), iprevec(8)
      Real*8   tau1, tau2, partlur, partlurout
      Integer  verb, UsedWr, UsedWi, ipBCG
      Integer  iter, info, nunit
      Real*8   resid

      EXTERNAL ddot
      Real*8   ddot


C=====LOCAL VARIABLES======

      Integer   i,n, iux,iuy,iuz, ip, iv1, nLapl, iLoop
      Real*8    x,y,z, eBC(3), rmax
      integer  ipLaplI,ipLaplR,ipMassI,ipMassR,ipFreeR,ipFreeI

c ==========================================================
c Step 1: load a mesh
      Call loadMani(
     &      nvmax, nbmax, ntmax,
     &      nv, nb, nt,
     &      vrt, bnd, tet, labelF, material,
     &      nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &      iW, iW, "../data/ramp.ani")

c ... generate a quasi-uniform mesh with 1000 tets starting from ramp.ani
      nEStar   = 7000    ! number of tets in generated mesh
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



c Step 2: generate the finite element system with P2-P1 elements
c     no data is provided for the user subroutine fem3Dext
      DATAFEM(1) = 0

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


c ... general sparse matrices in the AMG format (modified CSR format
c     where the diagonal element goes first)
      status = IOR(MATRIX_GENERAL, FORMAT_CSR)

c ... generate Laplacian for P2-FEM for a velocity component
      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3DextLaplP2, DATAFEM, status,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, iW)
      nLapl = nRow
      write(*,'(A)') '     P2 Laplacian is ready.'

c ... initialization of the preconditioner for P2 Laplacian
      verb    = 0     ! verbose no
      tau1    = 1d-2  ! absolute threshold for L,U
      tau2    = 1d-3  ! absolute threshold for T,R
      partlur = 0.5   ! even partition of memory between LU and R
      iERR    = 0     ! error code

      ipLaplI = 1
      ipLaplR = 1

      Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &           rW(ipLaplR), iW(ipLaplI), MaxWr, MaxWi, 
     &           partlur, partlurout,
     &           UsedWr, UsedWi, iERR)

      If(iERR.NE.0) Then
         Write(*,'(A,I4)')
     &      'Initialization1 of iluoo has failed, iERR=', iERR
         Stop
      End if
      write(*,'(A)') 'ILU: preconditioner for P2 Laplacian is ready.'

c ... generate mass matrix for P1 FEM
      ipMassI = UsedWi + 1
      ipMassR = UsedWr + 1

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3DextMassP1, DATAFEM, status,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi-ipMassI, iW(ipMassI))
      write(*,'(A)') '     P1 mass matrix is ready.'

c ... generate preconditioner for P1 mass matrix
      Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &           rW(ipMassR), iW(ipMassI), 
     &           MaxWr-ipMassR, MaxWi-ipMassI, 
     &           partlur, partlurout,
     &           UsedWr, UsedWi, iERR)

      If(iERR.NE.0) Then
         Write(*,'(A,I4)')
     &      'Initialization2 of iluoo has failed, iERR=', iERR
         Stop
      End if
      write(*,'(A)') 'ILU: preconditioner for P1 mass matrix is ready.'

      ipFreeI = ipMassI + UsedWi + 1
      ipFreeR = ipMassR + UsedWr + 1

c ... generate the Stokes P2-P1 matrix
      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &     fem3Dext, DATAFEM, status,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi-ipFreeI, iW(ipFreeI))

      If(ipFreeR + 8*nRow + nLapl.GT.MaxWr) Then
         Write(*,'(A,I7)') 'Increase MaxWr to ', ipFreeR + 8*nRow
         Stop
      End if


c Step 3: solve the finite element system  using libILU library
      ipBCG = ipFreeR +  nLapl            ! room for BCG work vectors
      iter  = 20000                       ! max number of iterations
      info  = 0                           ! no troubles on input
      nunit = 6                           ! output to display
      call dcopy(nRow, 0d0, 0, SOL, 1)    ! initial guess
      resid = ddot( nRow, RHS, 1, RHS, 1 )
      resid = 1d-7*dsqrt( resid )         ! final residual

      iprevec(1) = nRow
      iprevec(2) = ipLaplI
      iprevec(3) = ipLaplR
      iprevec(4) = ipMassI
      iprevec(5) = ipMassR
      iprevec(6) = ipFreeR
      iprevec(7) = nLapl            
      iprevec(8) = nv

      imatvec(1) = nRow
      Call slpbcgs(prevec,  iprevec, iW, rW,
     &             matvec,  imatvec, IA, JA, A,
     &             rW(ipBCG), nRow, 8,
     &             nRow, RHS, SOL,
     &             iter, resid, info, nunit)

      If(info.NE.0) Stop 'BiCGStab had failed'

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
     &                  nv, vrt, nt, tet, nb, bnd, labelF)

      Call GMVscalarVrt(SOL(ip), "pressure.gmv", 10,
     &                  nv, vrt, nt, tet, nb, bnd, labelF)


      End
