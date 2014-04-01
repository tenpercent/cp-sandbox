C ========================================================================
c This is slpbcgsnewton, the BiCGSTAB routine for determining (trial) inexact
C Newton steps.
C ========================================================================
      
      subroutine slpbcgsnewton (
     $   prevec, IPREVEC, iW,rW,   
     $   f, rpar, ipar,
     $   n, SOL, fcur,  step, 
     $   ITER_LIN, nfe, rsnrm,  
     $   r,rwork,
     $   INFO,iPrint)

      implicit none 
c---->-------------------------------------------------------------------
c   Argument types:
c

      EXTERNAL f, ddot, prevec
      INTEGER ipar(*),IPREVEC(*),iW(*)
      INTEGER N, ITER_LIN, INFO, nfe, iPrint
      DOUBLE PRECISION SOL(n), fcur(n), step(n),  rpar(*), rwork(n,6),
     $  r(n), rW(*) 
      DOUBLE PRECISION  rsnrm, ddot

c-------------------------------------------------------------------------
c   Local parameters
c
      double precision abstol, alpha, beta, omega, 
     $     rho, rhomns, tau, temp, TMP
      integer i, istb, itask, itrmjv

      double precision dlamch
      external dlamch
c
      double precision sfmin
      data sfmin /0.0d0/

      integer rtil,p,phat,v,t,rwork1
      rwork1    = 1
      rtil      = 2
      p         = 3
      phat      = 4
      v         = 5
      t         = 6

c ------------------------------------------------------------------------
c  Initialize sfmin only on first entry.
c
      if ( sfmin .eq. 0.0d0 ) sfmin = dlamch( 's' )
c ------------------------------------------------------------------------
c Set the stopping tolerance, initialize the step, etc. 
c ------------------------------------------------------------------------
      abstol = rsnrm
      do 10 i = 1, n
         step(i) = 0.d0
 10   continue
      istb = 0
c ------------------------------------------------------------------------ 
c Set up r and rtil. 
c ------------------------------------------------------------------------
      call dcopy(n,fcur,1,r,1)
      temp = -1.d0
      call dscal(n,temp,r,1)
      call dcopy(n,r,1,rwork(1,rtil),1)
c ------------------------------------------------------------------------
c Top of the iteration loop. 
c ------------------------------------------------------------------------
 100  continue
      istb = istb + 1
c ------------------------------------------------------------------------
c Perform the first "half-iteration". 
c ------------------------------------------------------------------------
      rho = ddot(n,rwork(1,rtil),1,r,1)
      if (istb .eq. 1) then 
         call dcopy(n,r,1,rwork(1,p),1)
      else
         if ( abs(rhomns) .lt. sfmin*abs(rho) ) then
            INFO = 4
            goto 900
         else
            beta = (rho/rhomns)*(alpha/omega)
            call daxpy(n,-omega,rwork(1,v),1, rwork(1,p), 1)
            call dscal(n,beta,rwork(1,p),1)
            call daxpy(n,1.d0,r,1,rwork(1,p),1)
         endif
      endif
      itask = 2
      call jvnewton(n, SOL, fcur, f, rpar, ipar,        
     $          itask, nfe, rwork(1,p),rwork(1,phat), 
     $          prevec, IPREVEC, iW,rW,
     $          rwork(1,rwork1), itrmjv)
         if (itrmjv .gt. 0) then 
            INFO = 2
            go to 900
         endif
      itask = 0
      call jvnewton(n, SOL, fcur, f, rpar, ipar,        
     $     itask, nfe, rwork(1,phat), rwork(1,v), 
     $     prevec, IPREVEC, iW,rW,
     $     rwork(1,rwork1), itrmjv)
      if (itrmjv .gt. 0) then 
         INFO = 1
         go to 900
      endif
      tau = ddot(n,rwork(1,rtil),1,rwork(1,v),1)
      if ( abs(tau) .lt. sfmin*abs(rho) ) then
         INFO = 4
         goto 900
      else
         alpha = rho/tau
      endif
      call daxpy(n,-alpha,rwork(1,v),1,r,1)
      call daxpy(n,alpha,rwork(1,phat),1,step,1)
      TMP = ddot(n,r,1,r,1)
      rsnrm = dsqrt(TMP)
c ------------------------------------------------------------------------
c Test for termination. 
c ------------------------------------------------------------------------
      if (rsnrm .le. abstol) then 
         if (iPrint.ge.2) write(*,*)'BCG iterate ',istb,' ||r||=',rsnrm
         INFO = 0
         go to 900
      endif
c ------------------------------------------------------------------------
c Perform the second "half-iteration". 
c ------------------------------------------------------------------------
      itask = 2
      call jvnewton(n, SOL, fcur, f, rpar, ipar,        
     $         itask, nfe, r, rwork(1,phat), 
     $         prevec, IPREVEC, iW,rW,
     $         rwork(1,rwork1), itrmjv)
         if (itrmjv .gt. 0) then 
            INFO = 2
            go to 900
         endif
      itask = 0
      call jvnewton(n, SOL, fcur, f, rpar, ipar,        
     $          itask, nfe, rwork(1,phat), rwork(1,t), 
     $          prevec, IPREVEC,iW,rW,
     $          rwork(1,rwork1), itrmjv)
      if (itrmjv .gt. 0) then 
         INFO = 1
         go to 900
      endif
      TMP = ddot(n,rwork(1,t),1,rwork(1,t),1)
      tau = dsqrt(TMP)
      tau = tau*tau
      temp = ddot(n,rwork(1,t),1,r,1)
      if ( tau .le. sfmin*abs(temp) ) then
         INFO = 4
         goto 900
      else
         omega = temp/tau
      endif
      if ( abs(omega) .lt. sfmin*abs(alpha) ) then 
         INFO = 4
         go to 900
      endif
      call daxpy(n,-omega,rwork(1,t),1,r,1)
      call daxpy(n,omega,rwork(1,phat),1,step,1)
      TMP = ddot(n,r,1,r,1)
      rsnrm = dsqrt(TMP)
c ------------------------------------------------------------------------
c Test for termination. 
c ------------------------------------------------------------------------
      if (iPrint.ge.2) write(*,*)'BCG iterate ',istb,' ||r||=',rsnrm
      if (rsnrm .le. abstol) then 
         INFO = 0
         go to 900
      endif
      if (istb .ge. ITER_LIN) then 
         INFO = 3
         go to 900
      endif
c ------------------------------------------------------------------------
c If continuing, update and return to the top of the iteration loop. 
c ------------------------------------------------------------------------
      rhomns = rho
      go to 100
c ------------------------------------------------------------------------
c Bottom of the iteration loop. 
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
c ------------------------------------------------------------------------
c Returning number of solver iteration
      ITER_LIN = istb  
      return
      end
