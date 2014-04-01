c ==========================================================
      Program HybridMixedFEM
c ==========================================================
c  The adaptive solution for the following bvp:
c
c
c   -div u = 0         in  Omega    (mass conservation eqn)
c        u = D grad p  in  Omega    (constitutive equation)
c        p = 0         on  dOmega_1
c        p = 2         on  dOmega_2
c
c
c  where Omega is the unit cube (dOmega_1) with a cubic hole 
c  (dOmega_2) in the center. The full anisotropic diffusion 
c  tensor D is defined by three rotations of a diagonal 
c  anisotropic tensor (see routine Drotate at the end of file 
c  forlibfem.f).
c
c  We use the mixed-hybrid finite element method with the lowest 
c  order Raviart-Thomas elements. The method results in a 
c  problem with a symmetric positive definite stifness matrix 
c  for the Lagrange multipliers.
c ==========================================================
      implicit none

C=====GLOBAL PARAMETERS======

      integer nvmax,ntmax,nbmax,niamax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... niamax- maximum number of rows
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax = 10 000, ntmax = 6*nvmax, nbmax = 10 000)
      parameter(niamax = 11*nvmax)
      parameter(namax  =  7*niamax)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 10 000 000, MaxWi = 20 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


C=====DATA FOR LIBRARY LIBFEM======

      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA(niamax), JA(namax)
      Real*8    A(namax)

      Real*8   RHS(niamax), SOL(niamax), RES(niamax)

      Real*8   SOLatTets(ntmax), VECatTets(3,ntmax), SOLatNodes(nvmax)

c  allocate data (none) for user defined routines
      Real*8   DATAFEM(3, 3)

      EXTERNAL FEM3Dext, Ddiff, Dbc, Drhs
      Integer  Ddiff, Dbc, Drhs
      Integer  status, nRow, nCol


C=====MESH ALLOCATION======

c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax),labelF(nbmax)

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(16), nbfix, ibfix(1), ntfix, itfix(1) 

c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)


C=====DATA FOR LIBRARY LIBILU======

      Real*8   tau1, tau2, partlur, partlurout
      Integer  verb, UsedWr, UsedWi, ipBCG

      Integer  iter, info, nunit
      Real*8   resid

      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)


C=====DATA FOR LIBRARY LIBMLR======

      Real*8   Lp
      Real*8   Metric(6, nvmax)


C=====DATA FOR LIBRARY LIBMBA======

c Basket capacity and the maximum number of local modifications
      Integer   MaxSkipE, MaxQItr
      Parameter(MaxSkipE = 300, MaxQItr = 50 000)

c Desired final mesh quality
      Real*8    Quality
      Parameter(Quality = 4D-1)

c Number of adaptive loops (metric generation -> mesh adaptation)
      Integer   nLOOPs
      Parameter(nLOOPs = 4)

c Other parameters, see comments before the library call
      Real*8   rQuality
      Logical  flagAuto
      Integer  nEStar, iPrint


C=====LOCAL VARIABLES======

      Integer   i, j, iERR, iLoop
      Real*8    rmax, theta, PI

c ==========================================================
C STEP1: Load a regular mesh and refine it uniformly (optional)

c ... Load the mesh in aft-format
      Call loadMaft(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/cube_hole.out")

      Write(*,*)
      Write(*,'(A,I6,A)') 'The loaded mesh has ',nt,' tetrahedra'

c ... initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(nv, nt, vrt, tet, MapMtr, Ref2MapMtr)

c ... refine the mesh uniformly 
      Do i = 1, 1
         Call uniformRefinement(
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, bnd, tet, labelF, material,
     &        MapMtr, Ref2MapMtr,
     &        iW, MaxWi)
      End do
      Write(*,'(A,I6,A)') 'The refined mesh has ',nt,' tetrahedra'



C STEP2: Define the diffusion tensor
c     define the full diffusion tensor via rotations about coordinate
c     axes of the diagonal tensor

      Do i = 1, 3
         Do j = 1, 3
            DATAFEM(i, j) = 0D0
         End do
      End do

      DATAFEM(1, 1) = 1D2
      DATAFEM(2, 2) = 1D1
      DATAFEM(3, 3) = 1D0

      PI = 4 * datan(1D0)

      theta = PI / 3
      Call Drotate(DATAFEM, 1, theta)
      
      theta = PI / 4
      Call Drotate(DATAFEM, 2, theta)

      theta = PI / 6
      Call Drotate(DATAFEM, 3, theta)


C STEP3: Start adaptive iterations
      Do iLoop = 1, nLOOPs
         Write(*,'(/,2(A,I2))') '===> LOOP: ', iLoop, ' out of ', nLOOPs


C STEP3a: Generate FEM system
c  ===   assemble the stifness matrix

c        fill labels of mesh points with zeros
         Do i = 1, nv
            labelP(i) = 0
         End do

c        general sparce matrix in the AMG format (modifed CSR format
c        where the diagonal element goes first)
         status = IOR(MATRIX_GENERAL, FORMAT_AMG)

c        matrix is assembled for Lagrange multipliers defined on faces
         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &        FEM3Dext, DATAFEM, status,
     &        niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, iW)


C STEP3b: solve the system
c        initialization of the preconditioner
         verb    = 0     ! verbose no
         tau1    = 1d-2  ! absolute threshold for L,U
         tau2    = 3d-4  ! absolute threshold for T,R
         partlur = 0.5   ! even partition of memory between LU and R
         iERR    = 0     ! error code

         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,'(A,I4)')
     &      'Initialization of iluoo has failed, iERR = ', iERR
            Stop
         End if

c        iterative solution
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
            SOL(i) = 1d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL,
     &                iter, resid, info, nunit)

         If(info.NE.0) Stop 'BiCGStab failed'

c        check the residual
         Call mulAgen(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') '   Maximal norm of residual: ', rmax


c STEP3c: compute cell-centered velocities VECatTets
c             and cell-centered pressures  SOLatTets
         Call HMFEMrecoverUP(nv, vrt, nt, tet, material, DATAFEM,
     &                       SOL, SOLatTets, VECatTets, iW, MaxWi)

C  ===   draw the mesh and the solution recovered at cells
         Call GMVscalarTet(SOLatTets,"pressure.gmv", 10,
     &                     nv, vrt, nt, tet, nb, bnd, labelF)
C  ===   draw the mesh and the vector field recovered at cells
         Call GMVvectorTet(VECatTets,"velocity.gmv", 10,
     &                     nv, vrt, nt, tet, nb, bnd, labelF)

         If (iLoop.eq.nLOOPs) Stop


c STEP3d: interpolate cell-centered piecewise constant solution 
c        to mesh vertices using the ZZ interpolation method
         Call P02P1(nv, nb, nt, vrt, bnd, tet,
     &              SOLatTets, SOLatNodes, MaxWr, MaxWi, rW, iW)


c STEP3e: generate metric (from SOLatNodes) optimal for the L_p norm
c        we choose L_1 norm since to get a smoother mesh aroung singularities
c        Lp = 0             ! maximum norm
         Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOLatNodes,
     &                    vrt, nv, tet, nt, bnd, nb, Metric,
     &                    MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)

c STEP3f: generate an adaptive mesh with nEStar tetrahedra which is 
c        quasi-unifrom in metric Metric
         nEStar   = 30000   ! number of tets in generated mesh
         flagAuto = .TRUE.  ! default mesh generation options
         status = 1         ! forbid boundary triangles
         iPrint = 1         ! average level of output information

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
     &                        'unspecified error in mbaNodal')
      End do


      Stop
      End



c ==========================================================
      Subroutine HMFEMrecoverUP(nv, vrt, nt, tet, material, DATAFEM,
     &                          SOL, SOLatTets, VECatTets, iW, MaxWi)
c ==========================================================
c  Recovers pressures SOLatTets at cell centers and velocities 
c  VECatTets at cells centers from Lagrange multipliers SOL 
c  defined on mesh faces
c ==========================================================
      Implicit none
      include 'fem3Dtet.fd'

      Integer  nv, nt, MaxWi
      Real*8   vrt(3,*), DATAFEM(*)
      Integer  tet(4,*), material(*), iW(*)
      Real*8   SOL(*), SOLatTets(*), VECatTets(3,*)

      EXTERNAL FEM3Dext, Ddiff, Drhs
      Integer  Ddiff, Drhs

      Real*8    calSqr, calVol
      External  calSqr, calVol

C=====LOCAL VARIABLES======
      Integer   n,i,j,k, i1,i2,i3,i4, ir,ic, ifc, info
      Integer   nf, iIFE, inEP, iIEP, iiW, iSYS(5)

      Real*8    localA(5, 5), localF(5), xyzc(3), rW(100), s
      Real*8    face_area(4), vol

      Integer   Face2OpPnt(4)
      Data      Face2OpPnt / 4, 1, 2, 3 /

c ==========================================================
c     make a list E -> F
      iIFE = 1
      inEP = iIFE + 4 * nt
      iIEP = inEP + nv
      iiW  = iIEP + 4 * nt
      If(iiW .gt.MaxWi) Stop 'Please increase MaxWi'

      Call listE2F(nv, nf, nt, tet, iW(iIFE), iW(inEP), iW(iIEP))

c     loop over mesh cells
      Do n = 1, nt
         i1 = tet(1, n)
         i2 = tet(2, n)
         i3 = tet(3, n)
         i4 = tet(4, n)

c        compute mass matrix using material properties
c        we use the user-specified routine Ddiff()
         iSYS(1) = n
         iSYS(2) = i1
         iSYS(3) = i2
         iSYS(4) = i3
         iSYS(5) = i4

         Call fem3Dtet(
     &        vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4),
     &        IDEN, FEM_RT0, IDEN, FEM_RT0,
     &        material(n), Ddiff, DATAFEM, iSYS, 2,
     &        5, localA, ir, ic)

c        add mass conservation equation (constraint matrix)
c        we use the divergent theorem here as we did in forlibfem.f
         localA(1, 5) = calSqr(vrt(1,i1), vrt(1,i2), vrt(1,i3))
         localA(2, 5) = calSqr(vrt(1,i2), vrt(1,i3), vrt(1,i4))
         localA(3, 5) = calSqr(vrt(1,i3), vrt(1,i4), vrt(1,i1))
         localA(4, 5) = calSqr(vrt(1,i4), vrt(1,i1), vrt(1,i2))

         localA(5, 5) = 0D0

         Do i = 1, 4
            face_area(i) = localA(i, 5)
         End do

c        compute right hand side (Lagrange muptipliers are in SOL)
         Do i = 1, 4
            ifc = iW(iIFE + 4*(n-1) + i-1)
            localF(i) = localA(i, 5) * SOL(ifc)
         End do

c        we use the user-specified routine Drhs()
         Call fem3Dtet(
     &        vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4),
     &        IDEN, FEM_P0, IDEN, FEM_P0,
     &        material(n), Drhs, DATAFEM, iSYS, 2,
     &        1, localF(5), ir, ic)

c        solve the 5x5 system localA x = localF. The first four 
c        entries in the solution x are velocities on tetrahedron faces.
c        The 5-th entry is the solution p at the mass center of the
c        tetrahedron. The user may use the LAPACK routine dsysv.

         Call dsysv('U', 5, 1, localA, 5, iW, localF, 5, rW, 100, info)
         If (info.ne.0) stop 'error in dsysv'

c        recover pressure
         SOLatTets(n) = localF(5)

c        recover velocity vector 
         vol = calVol(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))
         vol = dabs(vol)

         Do i = 1, 3
            xyzc(i) = (vrt(i,i1) + vrt(i,i2) + vrt(i,i3) + vrt(i,i4))/4
            VECatTets(i,n) = 0d0
         End do

         Do j = 1, 4 ! loop over faces
            k = Face2OpPnt(j)

            s = face_area(k) / (3*vol) * localF(j)

            i1 = tet(k, n)
            Do i = 1, 3
               VECatTets(i,n) = VECatTets(i,n) + s * (xyzc(i)-vrt(i,i1))
            End do
         End do
      End do

      Return
      End

