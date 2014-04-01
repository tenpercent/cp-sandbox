C ==========================================================
C ANI3D Version 1.4
C The file contains:
C    (1) main program
C    (2) user routine  iniQ
C    (3) user routines tkrFncX
C ==========================================================
      Program Main_Metric
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
      Real*8  rQuality

      Integer  ANI_Metric_Eucl
      EXTERNAL ANI_Metric_Eucl

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
     &      iW, iW, "ramp.ani")


C ... draw the initial mesh
      Call saveMgmv(
     &     nP, nF, nE, XYP, IPF, IPE, lbF, lbE, 
     &     "ramp.gmv", iW)


c ... generate adaptive mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status = 1         ! forbid boundary triangles
      iPrint = 1         ! average level of output information


      Call mbaAnalytic(
c group (M)
     &     nP, MaxP, nF, MaxF, nE, MaxE,
     &     XYP, IPF, IPE, lbF, lbE,
     &     nEStar, 
c group (Dev)
     &     nPv, nFv, nEv, IPV, IFV, IEV, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     ANI_Metric_Eucl, Quality, rQuality,
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

