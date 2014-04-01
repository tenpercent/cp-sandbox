c ==========================================================
      Program main
c ==========================================================
c The adaptive iterative solution of the following boudary value 
c problem (bvp):

c -div D grad u = 1 in Omega
c             u = 0 on dOmega_1
c         du/dn = 0 on dOmega_2
c
c where Omega is the ramp (0,1)^3 \ [0,.5]^2 x (0,1)
c dOmega_2 = {z=1}or{z=0}, dOmega_1 = dOmega \ dOmega_2 and 
c
c    D = diag{1,  1, 100}
c
c The program generates optimal mesh with nEStar elements.
c ==========================================================
      implicit none

C=====GLOBAL PARAMETERS======

      integer nvmax,ntmax,nbmax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax = 30 000, ntmax = 6*nvmax, nbmax = 10 000)
      parameter(namax = 300 000)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 2 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


C=====FOR LIBRARY LIBFEM======

      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), DA(nvmax), RHS(nvmax), SOL(nvmax), RES(nvmax)

      Real*8   DATAFEM(1)

      EXTERNAL Ddiff, Dbc, Drhs
      Integer  Ddiff, Dbc, Drhs
      Integer  status, nRow, nCol


C=====MESH ALLOCATION======

c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)

c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax),labelF(nbmax)

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(12), nbfix, ibfix(1), ntfix, itfix(1) 
      DATA      nvfix/0/,  nbfix/0/, ntfix/0/


C=====FOR LIBRARY LIBMLR======

      Real*8   Lp
      Real*8   Metric(6, nvmax)


C=====FOR LIBRARY LIBMBA======

c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr
      Parameter(MaxSkipE = 500, MaxQItr = 100 000)

c Desired final mesh quality 
      Real*8    Quality
      Parameter(Quality = 4D-1)

c Number of adaptive loops (metric generation -> mesh adaptation)
      Integer   nLOOPs
      Parameter(nLOOPs = 7)

c Other parameters, see comments before the library call
      Real*8   rQuality
      Logical  flagAuto
      Integer  nEStar, iPrint, iERR


C=====For LIBRARY LIBILU======

      Real*8   tau1, tau2, partlur, partlurout
      Integer  verb, UsedWr, UsedWi, ipBCG

      Integer  iter, info, nunit
      Real*8   resid

      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)


C=====LOCAL VARIABLES======

      Integer   i, iLoop
      Real*8    rmax

c ==========================================================
C STEP1: Load a regular mesh occupying part of a unit cube.

      Call loadMani(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/ramp.ani")

      Write(*,*)
      Write(*,'(A,I6,A)') 'The loaded mesh has ',nt,' tetrahedra'


C STEP2: Start adaptive iterations

      Do iLoop = 1, nLOOPs
         Write(*,'(/,2(A,I2))') '===> LOOP: ', iLoop, ' out of ', nLOOPs

c  ===   assemble the stifness matrix
c        no data is provided for the user subroutine Ddiff
         DATAFEM(1) = 0

c        symmetric sparse matrix in a 0-based CSC format used in UMFPACK
c        status = IOR(MATRIX_SYMMETRIC, FORMAT_CSC)

c        general sparse matrix in the AMG format (modified CSR format
c        where the diagonal element goes first)
         status = IOR(MATRIX_GENERAL, FORMAT_AMG)

         Call BilinearFormVolume(
     &        nv, nt, vrt, tet, material,
     &        GRAD, FEM_P1, GRAD, FEM_P1,
     &        Ddiff, DATAFEM, 1,
     &        status, nvmax, namax,
     &        IA, JA, DA, A, nRow, nCol,
     &        MaxWi, iW)

         Write(*,'(A,I8)') 'Number of non-zero entries:', IA(nRow+1)-1

c  ...   assemble the right-hand side
         Call LinearFormVolume(
     &        nv, nt, vrt, tet, material,
     &        FEM_P1,
     &        Drhs, DATAFEM, 2,
     &        RHS, nRow,
     &        MaxWi, iW)

         Write(*,'(A,I22)') 'Problem size:', nRow

c  ...   set up boundary conditions
         Call BoundaryConditions(
     &        nv, nb, nt, vrt, bnd, tet, labelF,
     &        FEM_P1,
     &        Dbc, DATAFEM,
     &        IA, JA, DA, A, RHS, nRow,
     &        MaxWi, iW)


c  ===  call the driver for LU factorization and solution
c  ...  initialization of the preconditioner
         verb    = 0     ! verbose no
         tau1    = 1d-2  ! absolute threshold for L,U
         tau2    = 1d-3  ! absolute threshold for T,R
         partlur = 0.5   ! even partition of memory between LU and R
         iERR    = 0     ! error code

         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,'(A,I4)')
     &         'Initialization of iluoo has failed, iERR=', iERR
            Stop
         End if


c  ...   iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop
         End if

         ipBCG = UsedWr + 1

         iter  = 20000  ! max number of iterations
         resid = 1d-12  ! final residual
         info  = 0      ! no troubles on input
         nunit = 6      ! output to display 
         Do i = 1, nRow ! initial guess
            SOL(i) = 0d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL,
     &                iter, resid, info, nunit)

         If(info.NE.0) Stop 'BiCGStab had failed'


c  ...  check the residual
         Call mulAgen(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') '   Maximal norm of residual: ', rmax


C  ===   draw the mesh and the solution
         Call GMVscalarVrt(SOL,"mesh.gmv", 10,
     &                     nv, vrt, nt, tet, nb, bnd, labelF)

         If(iLoop.eq.nLOOPs) Stop


c  ===  generate metric (from SOL) optimal for the L_p norm
c        Lp = 0             ! maximum norm
         Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL,
     &                    vrt, nv, tet, nt, bnd, nb, Metric,
     &                    MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c  ===  generate the adaptive mesh
         nEStar   = 10000   ! the final number of elements
         flagAuto = .TRUE.  ! default mesh generation options
         status   = 1       ! forbid boundary triangles (see aniMBA/status.fd)
c        status   = 0       ! no options
         iPrint   = 1       ! average level of output information


         Call mbaNodal(
c group (M)
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, bnd, tet, labelF, material,
     &        nEStar, 
c group (Dev)
     &        nvfix, nbfix, ntfix, ivfix, ibfix, itfix, 
     &        flagAuto, status,
c group (Q)
     &        MaxSkipE, MaxQItr,
     &        Metric, Quality, rQuality,
c group (W)
     &        MaxWr, MaxWi, rW, iW,
     &        iPrint, iERR)

         If(iERR.GT.1000) Call errMes(iERR, 'main',
     &                        'unspecified error if mbaNodal')

      End do

      End





