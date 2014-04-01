C ======================================================================
c The program demonstrates metric generation from edge or cell based
c data representing the local error. In the program, the interpolation 
C error is calculated using a given function Func. The error calculation 
C may be replaced by any a posteriori error estimates.
C ======================================================================


C ======================================================================
      Real*8 Function Func(xy)
C ======================================================================
      implicit none
      Real*8 xy(3), a

      a = 2.00D0
      Func = exp(2*(dabs(xy(1))**a + dabs(xy(2))**a)) + xy(3)**2

      Return
      End



C ======================================================================
      Subroutine FuncGrd(xy, fGrd)
C ======================================================================
      implicit none
      Real*8 xy(3), fGrd(3), a, u

      a = 2.00D0
      u = exp(2*(dabs(xy(1))**a + dabs(xy(2))**a))

      fGrd(1) = u * 2*a * xy(1)**(a-1)
      fGrd(2) = u * 2*a * xy(2)**(a-1)
      fGrd(3) = 2D0

      Return
      End



C ======================================================================
      Program Est2MetricRecovery
C ======================================================================
      implicit none
      integer nvmax,ntmax,nbmax
c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)

c working memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr =  1 000 000, MaxWi = 5 000 000)

      Integer iW(MaxWi)
      Real*8  rW(MaxWr)


c ... for library aniMBA
      Integer    MaxPv,     MaxEv,     MaxFv
      Parameter(MaxPv = 1, MaxEv = 1, MaxFv = 1)

      Integer  nv, nt, nb, nPv, nFv, nEv
      Real*8   vrt(3,nvmax)
      Integer  material(ntmax), tet(4,ntmax), bnd(3,nbmax), lbF(nbmax)
      Integer  IPV(MaxPv), IFV(MaxFv), IEV(MaxEv)

c Function
      Real*8   metric_h(6,nvmax), Func
      EXTERNAL Func

c Errors
      Real*8   Error(6,ntmax), U(nvmax)


c local variables
      Integer  i, n
      
C ======================================================================
c Step 1.  Load and refine the mesh
C ======================================================================
c ... load the initial mesh
      Call loadMani(
     &      nvmax, nbmax, ntmax,
     &      nv, nb, nt,
     &      vrt, bnd, tet, lbF, material,
     &      nPv, nFv, nEv, IPV, IFV, IEV,
     &      iW, iW, "../data/cube.ani")

c     Call saveMgmv(nv, nb, nt, 
c    &              vrt, bnd, tet, lbF, material, 
c    &              "cube.gmv", iW)

      Write(*,'(A,I6,A)') 'The input mesh containes ',nt,' tetrahedra'


C ======================================================================
c Step 2. Compute Errors
C ======================================================================
         Do i = 1, nv
            U(i) = Func(vrt(1,i))
         End do

c ... edge-based error estimates 
      Call edgeIntL8(U, nv, vrt, nt, tet, Error)


C ======================================================================
c Step 3. Compute Metric
C ======================================================================
c ... generate nodal tensor metric from edge-based errors by LeastSquares
         Call EdgeEst2MetricLS(nv, nt, vrt, tet, 
     &                         Error, metric_h,       
     &                         MaxWr, MaxWi, rW, iW) 


c ... generate nodal tensor metric from edge-based errors by method of shifts
c        Call EdgeEst2MetricMAX(Error,
c    &                          nv, nt, vrt, tet, 
c    &                          metric_h, MaxWr, rW)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nv
           Write(10, '(6E18.10)') (metric_h(i, n), i=1,6)
        End do
      Close(10)

      Write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'

      Stop 
      End


