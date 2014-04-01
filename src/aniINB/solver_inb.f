C =======================================================
C       Subroutine slInexactNewton sets up parameters of Newton
C       iterative method and calls driver subroutine drvnewton.
C =======================================================
        subroutine slInexactNewton(
     $       prevec, IPREVEC,iW,rW,
     $       f, rpar, ipar,      
     $       n, SOL, 
     $       RESID, STPTOL,
     $       WORK, LENWORK,
     $       INFO)

      implicit none  

      include 'inbdflts.h'
      include 'inbparam.h'
      integer n, ipar(*), IPREVEC(*), INFO(*), LENWORK 
      double precision SOL(n), RESID, STPTOL, WORK(*), rpar(*) 
      double precision ddot, rW(*)
      integer iW(*)
      external f, ddot, prevec
      integer  ITER_LIN, ITER_BT, ITER, iPrint 


       if (LENWORK.lt.11*n) then
           INFO(1) = 7
           return
       end if
       ITER_LIN = INFO(2)
       ITER     = INFO(3)
       ITER_BT  = INFO(4)
       iPrint   = INFO(5)
       etamax = DFLT_ETA_MAX
       thmin = DFLT_THMIN
       thmax = DFLT_THMAX
       
      call drvnewton(
     $       prevec, IPREVEC, iW,rW,   
     $       f, rpar, ipar,
     $       n, SOL,
     $       RESID, STPTOL, ITER, 
     $       ITER_LIN, ITER_BT,          
     $       WORK,
     $       INFO, iPrint )
     
        INFO(2) = ITER_LIN
        INFO(3) = ITER
        INFO(4) = ITER_BT
     
      return
      end

C ===============================================================
C This is subroutine drvnewton, the driver routine for the Newton
C iterative method.
C ===============================================================

        subroutine drvnewton(
     $     prevec, IPREVEC, iW,rW,   
     $     f, rpar, ipar,
     $     n, SOL, 
     $     RESID, STPTOL, ITER, 
     $     ITER_LIN, ITER_BT,         
     $     WORK,
     $     INFO, iPrint)

      implicit none  
c---->------------------------------------------------------------
c       Argument types:
c

      EXTERNAL f, ddot, prevec
      INTEGER n, ITER, ITER_LIN, ITER_BT
      INTEGER ipar(*), IPREVEC(*), iW(*), INFO(*), iPrint
      double precision SOL(n), WORK(n,6), rpar(*), rW(*)
      DOUBLE PRECISION RESID, STPTOL, ddot

      include 'inbparam.h'
      include 'inbdflts.h'  
c -----------------------------------------------------------------
c Explanation: 
c  prevec   : extern : Precondition-vector routine
c
c  IPREVEC, iW,rW  Configuration data for 'prevec'
c
c  f  : extern :  name of user-supplied subroutine for evaluating the function 
c            the zero of which is sought; this routine has the form
c
c                 subroutine f(n, SOL, fcur, rpar, ipar, itrmf)
c
c  rpar    = real parameter/work array passed to the f and jacv routines. 
c
c  ipar    = integer parameter/work array passed to the f and jacv routines.
c
c  n       = dimension of the problem.
c
c  SOL    = vector of length n, initial guess on input and final 
c            approximate solution on output. 
c
c  WORK   = real work vector for use by the Krylov solver.  
c
c  RESID    = stopping tolerance on the f-norm.
c
c  STPTOL  = stopping tolerance on the steplength.
c
c  ITER  = maximum allowable number of nonlinear iterations.
c
c  ITER_LIN  = maximum allowable number of linear iterations.
c
c  ITER_BT  = maximum allowable number of back-tracking iterations.
c
c  iPrint  = print level (0 mum, 1 nonlinear residuals, 2 for linear residuals)
c
c ------------------------------------------------------------------------
c       Local variables:
c

      double precision alpha, epsmach, eta, etamin, fcnrm, 
     $     flmnrm, fpnrm, oftjs, oftlm, redfac, rsnrm, 
     $     stpnrm, temp 
      integer ibt, itrmbt, itrmf, itrmks, ibtmax,iksmax,
     $      ncall, nni, nli, nfe

      double precision dlamch, TMP
      external dlamch

      integer fcur, xpls, fpls, step, r, rwork

c ------------------------------------------------------------------------
      data ncall / 0 /
      save ncall, alpha, epsmach 
c ------------------------------------------------------------------------
c Initialize.
c ------------------------------------------------------------------------

       fcur  = 1
       xpls  = 2
       fpls  = 3
       step  = 4
       r     = 5
       rwork = 6

      if ( ncall .eq. 0 ) epsmach = 2.0d0*dlamch( 'e' )
      ncall = 1
      nni = 0
      nli = 0
      ibt = 0
      nfe = 0
      iksmax= ITER_LIN 
      ibtmax = ITER_BT
      alpha = DFLT_CHOICE1_EXP 
      eta = .5d0
      
      call f(n, SOL, WORK(1,fcur), rpar, ipar, itrmf)
      nfe = nfe+1
      if ( itrmf .ne. 0 ) then
         INFO(1) = 2
         go to 900
      endif
      fcnrm = dsqrt(ddot(n, WORK(1,fcur), 1, WORK(1,fcur), 1)) 
c ------------------------------------------------------------------------ 
c Nonlinear iteration loop.
c ------------------------------------------------------------------------
 100  continue
c ------------------------------------------------------------------------ 
c Test for stopping. 
c------------------------------------------------------------------------
      if (iPrint.ge.1) write(*,'(A,I2,A,E15.8)') 
     &                 'INB: Newton step ', nni, '  ||f|| =',fcnrm
      if (fcnrm .le. RESID) then 
         INFO(1) = 0
         go to 900
      endif
      if (nni .gt. 0 .and. stpnrm .le. STPTOL .and. itrmks .eq. 0) then 
         INFO(1) = 0
         go to 900
      endif
      if (nni .ge. ITER) then 
         INFO(1) = 1
         go to 900
      endif
c ------------------------------------------------------------------------
c Compute the (trial) inexact Newton step with the Krylov solver. 
c ------------------------------------------------------------------------
         ITER_LIN = iksmax
         rsnrm = fcnrm* eta
         
         call slpbcgsnewton (
     $    prevec, IPREVEC, iW, rW,   
     $    f, rpar, ipar,
     $    n, SOL, WORK(1,fcur), WORK(1,step),  
     $    ITER_LIN, nfe, rsnrm,
     $    WORK(1,r),WORK(1,rwork), 
     $    itrmks,iPrint)
     
        nli = nli + ITER_LIN
c------------------------------------------------------------------------
      if (itrmks .eq. 1 .or. itrmks .eq. 2) then
         INFO(1) = itrmks + 2
         go to 900
      endif
      if (itrmks .ge. 3) then 
         if (rsnrm/fcnrm .gt. 1.d0) then 
            INFO(1) = 5
            go to 900
         else
            temp = dlog(RESID/fcnrm)/
     $           dlog(rsnrm/((1.d0 + 10.d0*epsmach)*fcnrm))
            if (temp .gt. 1000.d0*dfloat(ITER - nni)) then 
               INFO(1) = 5
               go to 900
            endif
         endif
      endif
c ------------------------------------------------------------------------
c Compute the original value of f(transpose)*Js for backtracking; the 
c original value of f(transpose)*(linear model) is also computed for 
c later use.   
c -----------------------------------------------------------------------
      oftlm = -ddot(n, WORK(1,fcur), 1, WORK(1,r), 1)
      oftjs = oftlm - fcnrm**2
c ------------------------------------------------------------------------
c Determine an acceptable step via backtracking. 
c ------------------------------------------------------------------------
        ITER_BT = ibtmax
      call btnewton(n, SOL, fcnrm, WORK(1,step), eta, WORK(1,xpls), 
     $     WORK(1,fpls), fpnrm, oftjs, 
     $     redfac,  ITER_BT, nfe, f, rpar, ipar, itrmbt)
        ibt = ibt + ITER_BT
      if (itrmbt .eq. 1) then
         INFO(1) = 6
         go to 900
      else if (itrmbt .eq. 2) then
         INFO(1) = 2
         go to 900
      endif
c ------------------------------------------------------------------------
c Set eta for next iteration. 
c ------------------------------------------------------------------------
         etamin = eta**alpha
         temp = 1.d0 - redfac
         flmnrm = dsqrt((temp*fcnrm)**2 + 2.d0*redfac*temp*oftlm + 
     $        (redfac*rsnrm)**2)
         eta = dabs(fpnrm - flmnrm)/fcnrm 
         if (etamin .le. eta_cutoff) etamin = 0.d0
         if (eta .lt. etamin) eta = etamin 
         if (eta .gt. etamax) eta = etamax 
         if (eta*fpnrm .le. 2.d0*RESID) eta = (.8d0*RESID)/fpnrm
c ------------------------------------------------------------------------
c Update SOL, fcur, fcnrm, stpnrm, nni for next iteration.
c ------------------------------------------------------------------------
      call dcopy(n, WORK(1,xpls), 1, SOL, 1)
      call dcopy(n, WORK(1,fpls), 1, WORK(1,fcur), 1)
      fcnrm = fpnrm
      TMP = ddot(n,WORK(1,step),1,WORK(1,step),1)
      stpnrm = dsqrt(TMP)
      nni = nni + 1
c ------------------------------------------------------------------------
c Return to top of loop for next iteration.
c ------------------------------------------------------------------------
      go to 100
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
        ITER = nni
        ITER_LIN = nli
        ITER_BT = ibt
        INFO(5) = nfe
        
      RESID = fcnrm
c ------------------------------------------------------------------------ 
      return
      end


