c ======================================================================
      Program ConvectionDiffusion
c ======================================================================
c The adaptive iterative solution of the boundary value problem:
c
c  div M grad u + (M+L) grad div u = 0 in Omega
c                                u = u_0 on Gamma1
c                                u = u_0 on Gamma2
c
c where Omega is the union of two parallepipeds rotated about z-axis,
c M and L are Lame coefficients, and pull boxes away form each other.
c
c  u_0 = (0, 0, z).
c 
c  GMV file contains only z-component of the displacement u.
c
c ======================================================================
      implicit none

C=====GLOBAL PARAMETERS======

      integer nvmax,ntmax,nbmax,nfmax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... nfmax - maximum number of columns in matrix
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax = 10 000,  ntmax = 6*nvmax, nbmax = 10 000)
      parameter(nfmax = 3*nvmax, namax = 40*nfmax)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 7 000 000, MaxWi = 7 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


C=====FOR LIBRARY LIBFEM======

      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), SOL(nfmax), RES(nfmax)

      Real*8   DATAFEM(1)

      EXTERNAL Dbc,  FEM3Dext
      Integer  Dbc
      Integer  status, nRow, nCol, iSYS(5)

      
C=====MESH ALLOCATION======

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)
      Integer  labelP(nvmax)

c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax), material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax), labelF(nbmax)

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(8), nbfix, ibfix(1), ntfix, itfix(1) 
      DATA      nvfix/0/,  nbfix/0/, ntfix/0/


C=====FOR LIBRARY LIBMLR======

      Real*8   Lp
      Real*8   Metric(6, nvmax)


C=====FOR LIBRARY LIBMBA======

c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr
      Parameter(MaxSkipE = 400, MaxQItr = 40 000)

c Desired final mesh quality 
      Real*8    Quality
      Parameter(Quality = 4D-1)

c Number of adaptive loops (metric generation -> mesh adaptation)
      Integer   nLOOPs
      Parameter(nLOOPs = 5)

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

      Integer  i, iLoop, ibc, iv1, iv2, iv3
      Real*8   x, y, z, rmax, eBC(1)

c ======================================================================
c Load a mesh for a prism
      Call loadMani(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/twoboxes.ani")

      Write(*,*)
      Write(*,'(A,2(I6,A))') 'The loaded mesh has ',
     &     nv,' points and ', nt,' tetrahedra'


c BEGIN ADAPTIVE ITERATIONS

      Do iLoop = 1, nLOOPs
c  ===  assemble the stifness matrix
         Write(*,'(/,2(A,I2))') '===> LOOP: ', iLoop, ' out of ', nLOOPs

c  ...  mark the Dirichlet points via labelP. 
         Do i = 1, nv
            labelP(i) = 0
         End do

         Do i = 1, nb
            If(labelF(i).EQ.1 .OR. labelF(i).EQ.18) Then
               iv1 = bnd(1, i)
               iv2 = bnd(2, i)
               iv3 = bnd(3, i)

               x = (vrt(1, iv1) + vrt(1, iv2) + vrt(1, iv2)) / 3
               y = (vrt(2, iv1) + vrt(2, iv2) + vrt(2, iv2)) / 3
               z = (vrt(3, iv1) + vrt(3, iv2) + vrt(3, iv2)) / 3

               ibc = Dbc(x, y, z, labelF(i), DATAFEM, iSYS, eBC)

               If(ibc.EQ.BC_DIRICHLET) Then
                  labelP(iv1) = max(labelP(iv1), labelF(i))
                  labelP(iv2) = max(labelP(iv2), labelF(i))
                  labelP(iv3) = max(labelP(iv3), labelF(i))
               End if
            End if
         End do

c  ...  no data provided for the user subroutine Dconv used in FEM3Dext
         DATAFEM(1) = 0


c  ...  general sparce matrix in the AMG format (modifed CSR format
c       where the diagonal element goes first)
         status = IOR(MATRIX_GENERAL, FORMAT_AMG)

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelP, bnd, labelF, tet, material,
     &        FEM3Dext, DATAFEM, status,
     &        nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, iW)


c  ===  call the driver for ILU factorization and solution
c  ...  initialization of the preconditioner
         verb    = 0     ! verbose no
         tau1    = 4d-3  ! absolute threshold for L,U
         tau2    = 4d-4  ! absolute threshold for T,R
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


c  ...  iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop
         End if

         ipBCG = UsedWr + 1

         iter  = 2000   ! max number of iterations
         resid = 1d-10  ! final residual
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
         Write(*,'(A,E12.6)') '     maximal norm of residual: ', rmax


C  ===   draw the adapted mesh and the solution
c        Call GMVscalarVrt(SOL, "mesh.gmv", 10,               ! x-displacements
c    &                   nv, vrt, nt, tet, nb, bnd, labelF)
c        Call GMVscalarVrt(SOL(nv+1), "mesh.gmv", 10,         ! y-displacements
c    &                   nv, vrt, nt, tet, nb, bnd, labelF)
         Call GMVscalarVrt(SOL(2*nv+1), "mesh.gmv", 10,       ! z-displacements
     &                   nv, vrt, nt, tet, nb, bnd, labelF)

         If(iLoop.eq.nLOOPs) Stop


c  ===  generate metric (from SOL) optimal for the L_p norm
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm

c        Call Nodal2MetricVAR(SOL,                            ! adapt to x-displacements
         Call Nodal2MetricVAR(SOL(2*nv+1),                    ! adapt to z-displacements
     &                        vrt, nv, tet, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c  ===  generate the adaptive mesh
         nEStar   = 30000   ! final number of elements
         flagAuto = .TRUE.  ! default mesh generation options
c        status   = 1       ! forbid boundary elements (see aniMBA/status.fd)
         status   = 0       ! no options
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





