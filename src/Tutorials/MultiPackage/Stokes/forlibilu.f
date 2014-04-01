c Preconditioner evaluation for the Stokes matrix
      SUBROUTINE prevec( iprevec, dummy, x, y, iwork, dwork)
      INTEGER   iprevec(*), dummy, iwork(*), n
      REAL*8    x(*),y(*),dwork(*)
      Integer   nv,nLapl,ne 
      Integer   ipLaplI,ipLaplR,ipMassI,ipMassR,ipFreeR,iX,iY,iZ

      n        =  iprevec(1)
c     call dcopy (n,x,1,y,1)

      ipLaplI  =  iprevec(2) - 1
      ipLaplR  =  iprevec(3) - 1
      ipMassI  =  iprevec(4) - 1
      ipMassR  =  iprevec(5) - 1
      ipFreeR  =  iprevec(6) - 1
      nLapl    =  iprevec(7)
      nv       =  iprevec(8)

      ne = nLapl - nv

c ... block diagonal preconditioner

c ... precondition Ux
      iX = 1
      iZ = ipFreeR + 1
      call dcopy (nv,x(iX),1,dwork(iZ),1)
      iX = 4*nv + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,x(iX),1,dwork(iZ),1)
      iZ = ipFreeR + 1
      call iluoo_solve (nLapl, dwork(ipLaplR + 3*nLapl+1),
     &                         iwork(ipLaplI + nLapl+2), 
     &                         iwork(ipLaplI + 2*nLapl+3), 
     &                         iwork(ipLaplI + 4*nLapl+3), dwork(iZ))
      iY = 1
      iZ = ipFreeR + 1
      call dcopy (nv,dwork(iZ),1,y(iY),1)
      iY = 4*nv + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,dwork(iZ),1,y(iY),1)

c ... precondition Uy
      iX = nv + 1
      iZ = ipFreeR + 1
      call dcopy (nv,x(iX),1,dwork(iZ),1)
      iX = 4*nv + ne + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,x(iX),1,dwork(iZ),1)
      iZ = ipFreeR + 1
      call iluoo_solve (nLapl, dwork(ipLaplR + 3*nLapl+1),
     &                         iwork(ipLaplI + nLapl+2), 
     &                         iwork(ipLaplI + 2*nLapl+3), 
     &                         iwork(ipLaplI + 4*nLapl+3), dwork(iZ))
      iY = nv + 1
      iZ = ipFreeR + 1
      call dcopy (nv,dwork(iZ),1,y(iY),1)
      iY = 4*nv + ne + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,dwork(iZ),1,y(iY),1)

c ... precondition Uz
      iX = 2*nv + 1
      iZ = ipFreeR + 1
      call dcopy (nv,x(iX),1,dwork(iZ),1)
      iX = 4*nv + 2*ne + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,x(iX),1,dwork(iZ),1)
      iZ = ipFreeR + 1
      call iluoo_solve (nLapl, dwork(ipLaplR + 3*nLapl+1),
     &                         iwork(ipLaplI + nLapl+2), 
     &                         iwork(ipLaplI + 2*nLapl+3), 
     &                         iwork(ipLaplI + 4*nLapl+3), dwork(iZ))
      iY = 2*nv + 1
      iZ = ipFreeR + 1
      call dcopy (nv,dwork(iZ),1,y(iY),1)
      iY = 4*nv + 2*ne + 1
      iZ = ipFreeR + nv + 1
      call dcopy (ne,dwork(iZ),1,y(iY),1)

c ... precondition P (scaling by -0.1 for better convergence)
      iX = 3*nv + 1
      iY = 3*nv + 1
      call dcopy (nv,x(iX),1,y(iY),1)
      call dscal (nv,-1d-1,y(iY),1)
      call iluoo_solve (nv, dwork(ipMassR + 3*nv+1),
     &                      iwork(ipMassI + nv+2), 
     &                      iwork(ipMassI + 2*nv+3), 
     &                      iwork(ipMassI + 4*nv+3), y(iY))

      return
      end

