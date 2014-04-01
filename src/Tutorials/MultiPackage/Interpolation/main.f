C =========================================================
c The program illustrates organization of an adaptive loop
c for building a mesh minimizing the interpolation error 
c for a given function (Func). The function is taken from 
c D'Azevedo's paper on optimal triangulations and adapted
c it to 3D.
C =========================================================


C =========================================================
      Real*8 Function Func(xyz)
C =========================================================
      implicit none
      Real*8 xyz(3)
      Real*8 x0, y0, z0, a, dx, dy, dz

      x0 = 0.5d0
      y0 = 0.5d0
      z0 =-0.2d0
  
      a = dsqrt(1.0d1)

      dx = (xyz(1) - x0) ** 2
      dy = (xyz(2) - y0) ** 2
      dz = (a*xyz(3) - z0) ** 2

      Func = (dz + dy - dz) / (dx + dy + dz) ** 2

      Return
      End


C =========================================================
      program AdaptiveInterpolation
C =========================================================
      implicit none

C=====GLOBAL PARMAETERS======

      integer nvmax,ntmax,nbmax,namax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... namax - maximum number of non-zero matrix entries
      parameter(nvmax = 15 000, ntmax = 6*nvmax, nbmax = 10 000)
      parameter(namax = 90 000)

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


C=====MESH ALLOCATION======

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
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
      Integer   nvfix, ivfix(8), nbfix, ibfix(1), ntfix, itfix(1) 
      DATA      nvfix/0/,  nbfix/0/, ntfix/0/


C=====FOR LIBRARY LIBMLR======

      Real*8   Lp
      Real*8   Metric(6, nvmax)


C=====For Library LIBMBA======

c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr
      Parameter(MaxSkipE = 100, MaxQItr = 50 000)

c Desired final mesh quality 
      Real*8    Quality
      Parameter(Quality = 4D-1)

c Number of adaptive loops (metric generation -> mesh adaptation)
      Integer   nLOOPs
      Parameter(nLOOPs = 5)

c Other parameters, see comments before the library call
      Real*8   rQuality
      Logical  flagAuto
      Integer  nEStar, status, iPrint, iERR


C=====LOCAL VARIABLES======

c function and nodal function U
      Real*8   Func, U(nvmax)
      EXTERNAL Func

c errors: maximum (L8) and gradient (H1)
      Real*8   error(ntmax)
      Real*8   errMin, errMax, errMean

c other
      integer i,  iLoop
      
C =========================================================
C STEP1: Load a regular mesh for a unit cube.

      Call loadMani(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/cube.ani")

      Write(*,*)
      Write(*,'(A,I6,A)') 'The loaded mesh has ',nt,' tetrahedra'


C STEP2: Start adaptive iterations. first, we estimate the 
c        interpolation error, then adapt the mesh to the
c        function U given by Func(xyz)
      Do iLoop = 1, nLOOPs
         Write(*,'(/,2(A,I2))') '===> LOOP: ', iLoop, ' out of ', nLOOPs

c define the function at mesh nodes
         Do i = 1, nv
            U(i) = Func(vrt(1, i))
         End do


c compute errors in maximum norm
         call tetL8error(U, nv, vrt, nt, tet, error)

         errMin  = error(1)
         errMax  = error(1)
         errMean = 0
         do i = 1, nt
            errMin  = min(errMin, error(i))
            errMax  = max(errMax, error(i))
            errMean = errMean + error(i)
         end do
         errMean = errMean / nt
         write(*,5003) 'L8', errMin, errMax, errMean

c Create a GMV-file with the mesh. The name must have extension .gmv
         Call saveMgmv(nv, nb, nt, 
     &                 vrt, bnd, tet, labelF, material, 
     &                 "save.gmv", iW)

         If (iLoop.eq.nLOOPs) Stop

c generate metric (from U) optimal for the L_p norm
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(U,
     &                    vrt, nv, tet, nt, bnd, nb, Metric,
     &                    MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c adapt the mesh to nodal function U by generating
c a quasi-uniform mesh ith nEStar = 10000 tetrahedra
c which is quasi-unifrom in the metric defined in subroutine 
c MetricPoint (see below).
      nEStar   = 10000   ! number of triangles in generated mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status = 1         ! forbid boundary triangles
      iPrint = 1         ! average level of output information


      Call mbaNodal(
c group (M)
     &     nv, nvmax, nb, nbmax, nt, ntmax,
     &     vrt, bnd, tet, labelF, material,
     &     nEStar, 
c group (Dev)
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     Metric, Quality, rQuality,
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)


         If(iERR.GT.1000) Call errMes(iERR, 'main',
     &                      'unspecified error in mbaNodal')


      End do


 5003 Format(A, ' errors (min/max/mean):', 3E10.3)

      Stop
      End


