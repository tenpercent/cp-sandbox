C ==========================================================
      Program main
C ==========================================================
C Program loads a mesh (ramp.ani), defines a mesh function 
C (SOL) and then generates a metric based on the Hession of 
C discrete function SOL. 
C ==========================================================
c Maximum number of nodes, elements and boundary edges
      Integer   MaxP, MaxE, MaxF
      Parameter(MaxP = 30000, MaxE = 4 * MaxP, MaxF = 14000)

c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr 
      Parameter(MaxSkipE = 500, MaxQItr = 20 000) 
 
c Available memory 
      Integer   MaxWr, MaxWi 
      Parameter(MaxWr = 4 500 000, MaxWi = 6 000 000) 
 
C ==========================================================
C group (M)
      Real*8  XYP(3, MaxP)
      Integer IPE(4, MaxE), IPF(3, MaxF)
      Integer IPV(MaxP), lbF(MaxF), lbE(MaxE)

C Group (Dev)
      Integer IFV(MaxF), IEV(MaxE) 

C group (Q)
      Real*8  Sol(MaxP)
      Real*8  Metric(6, MaxP)

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)

C ==========================================================
c ... load the initial mesh
      Call loadMani(
     &      MaxP, MaxF, MaxE,
     &      nP, nF, nE,
     &      XYP, IPF, IPE, lbF, lbE,
     &      nPv, nFv, nEv, IPV, IFV, IEV,
     &      iW, iW, "../data/ramp.ani")


c ... define a discrete function at mesh points. The function
c ... has a moderate boundary layer around x=1.
      Do n = 1, nP
         x = XYP(1, n)
         y = XYP(2, n)
         z = XYP(3, n)

         Sol(n) = (x**8)*(1-x) * y*(1-y) * z*(1-z)
      End do


c ... generate metric (from SOL) optimal for the L_p norm
      Lp = 0             ! maximum norm
c     Lp = 1             ! L_1 norm
      Call Nodal2MetricVAR(Sol,
     &                     XYP, nP, IPE, nE, IPF, nF, Metric,
     &                     MaxWr, rW, MaxWi, iW)
c     Call Nodal2MetricZZ( Sol,
c    &                     XYP, nP, IPE, nE, IPF, nF, Metric,
c    &                     MaxWr, rW, MaxWi, iW)


      If(Lp.GT.0) Call Lp_norm(nP, Lp, Metric)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nP
           Write(10, '(6E18.10)') (Metric(i, n), i=1,6)
        End do
      Close(10)

      Write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'

      Stop
      End


