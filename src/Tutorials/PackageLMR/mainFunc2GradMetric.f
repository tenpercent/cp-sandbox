C ======================================================================
C The program demonstrates building the optimal metric for the gradient
C of the P1 interpolation error. The metric depends on Lp norm of error
C the user wants to minimize. The interpolation error is calculated
c using a given function Func. The interpolation error may be replaced
C with any error estimates.
C ======================================================================


C ======================================================================
      Real*8 Function Func(xy)
C ======================================================================
      implicit none
      Real*8  xy(3)

      Func = xy(1)**2 * xy(2)
     &     + xy(2)**3 + dtanh(10*(dsin(5*xy(2)) - 2*xy(1))) + xy(3)**2

      Return
      End



C ======================================================================
      Subroutine FuncGrd(xy, fGrd)
C ======================================================================
      implicit none
      Real*8   xy(3), fGrd(3), u

      u = (dcosh(10*(dsin(5*xy(2)) - 2*xy(1))))**2
      fGrd(1) = 2*xy(1)*xy(2) - 20d0/u
      fGrd(2) = xy(1)**2 + 3*xy(2)**2 + 50*dcos(5*xy(2)) / u
      fGrd(3) = 2D0

      Return
      End



C ======================================================================
      Program mainFunc2Metric
C ======================================================================
      implicit none

c Maximum number of nodes, elements and boundary edges
      Integer   MaxP,  MaxE,  MaxF, MaxPV, MaxEV, MaxFV
      Parameter(MaxP = 20000, MaxF = 10000, MaxE = 2 * MaxP)
      Parameter(MaxPV = MaxP, MaxFV = MaxF, MaxEV = MaxE)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 2 000 000)

C ======================================================================
C group (M)
      Real*8  XYP(3, MaxP)
      Integer IPE(4, MaxE), IPF(3, MaxF), IPV(MaxPV)
      Integer nP, nF, nE

C group (Dev)
      Integer IFV(MaxFV), IEV(MaxEV), lbF(MaxF), lbE(MaxE)
      Integer nPv, nFv, nEv

C group (Q)
      Real*8  Metric(6, MaxP)

C group (W)
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)

C Local variables
      Integer  i, n

      Real*8   Func, Lp
      EXTERNAL Func

C ======================================================================
c ... load the initial mesh
      Call loadMani(
     &      MaxP, MaxF, MaxE,
     &      nP, nF, nE,
     &      XYP, IPF, IPE, lbF, lbE,
     &      nPv, nFv, nEv, IPV, IFV, IEV,
     &      iW, iW, "../data/cube.ani")


c ... generate metric (from SOL) optimal for the L_p norm
      Lp = 0             ! maximum norm
c     Lp = 1             ! L_1 norm

c the recovery of continuous metric based on methods of shifts for
c minimizing gradient of the interpolation error
      Call Func2GradMetricMAX(Func,  
     &                        nP, nE, XYP, IPE, 
     &                        Metric, MaxWr, rW)


      If(Lp.GT.0) Call Lp_gradnorm(nP, Lp, Metric)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nP
           Write(10, '(6E18.10)') (Metric(i, n), i=1,6)
        End do
      Close(10)

      Write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'

      Stop
      End


