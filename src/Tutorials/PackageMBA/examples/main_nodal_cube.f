C ==========================================================
      Program main
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
 
c Desired number of triangles in the final mesh 
      Integer   nEStar 
      Parameter(nEStar = 10 000) 
 
c Desired final mesh quality 
      Real*8    Quality 
      Parameter(Quality = 4D-1) 
 
C ==========================================================
C group (M)
      Real*8  XYP(3, MaxP)
      Integer IPE(4, MaxE), IPF(3, MaxF)
      Integer IPV(MaxP), lbF(MaxF), lbE(MaxE)

C Group (Dev)
      Integer IFV(MaxF), IEV(MaxE) 
      Logical flagAuto
      Integer status

C group (Q)
      Real*8  Sol(MaxP), rQuality, Lp
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
     &      iW, iW, "cube.ani")


C ... draw the initial mesh
      Call saveMgmv(
     &     nP, nF, nE, XYP, IPF, IPE, lbF, lbE, 
     &     "cube.gmv", iW)


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
      Call metricNodal(Sol,
     &                 XYP, nP, IPE, nE, IPF, nF, Metric,
     &                 MaxWr, rW, MaxWi, iW)

      If(Lp.GT.0) Call Lp_norm(nP, Lp, Metric)


c ... generate adaptive mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status = 1         ! forbid boundary triangles
      iPrint = 1         ! average level of output information


      Call mbaNodal(
c group (M)
     &     nP, MaxP, nF, MaxF, nE, MaxE,
     &     XYP, IPF, IPE, lbF, lbE,
     &     nEStar, 
c group (Dev)
     &     nPv, nFv, nEv, IPV, IFV, IEV, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     Metric, Quality, rQuality,
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)


c ... draw the final mesh
      Call saveMgmv(
     &     nP, nF, nE, XYP, IPF, IPE, lbF, lbE, 
     &     "save.gmv", iW)


c ... save the final mesh
      Call saveMani(
     &     nP, nF, nE,  
     &     XYP, IPF, IPE, lbF, lbE,
     &     nPv, nFv, nEv, IPV, IFV, IEV,
     &     .FALSE., iW, iW, "save.ani")

      Stop
      End





