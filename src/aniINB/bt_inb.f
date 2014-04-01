C=======================================================================
C       This is btnewton, the backtracking routine for the (inexact)
C       Newton iterations.
C=======================================================================

       subroutine btnewton(n, SOL, fcnrm, step, eta, xpls, fpls, fpnrm, 
     $     oftjs, redfac, ITER_BT ,nfe, 
     $     f, rpar, ipar, 
     $     itrmbt)

      implicit none 

      integer n,  ipar(*), itrmbt, ITER_BT,nfe
      double precision SOL(n), fcnrm, step(n), eta, ddot,
     $     xpls(n), fpls(n), fpnrm, oftjs, redfac, rpar(*) 
      external f, ddot

c
      double precision t, theta, TMP
      integer i, itrmf, ibtmax, ibt
      include 'inbparam.h'


c Initialize.
c ------------------------------------------------------------------------
      ibtmax = ITER_BT  
      t = 1.d-4
      ibt  = 0
      redfac = 1.d0
c ------------------------------------------------------------------------
c Backtracking loop. 
c ------------------------------------------------------------------------
 100  continue
      do 110 i = 1, n
         xpls(i) = SOL(i) + step(i) 
 110  continue
      call f(n, xpls, fpls, rpar, ipar, itrmf)
      nfe = nfe+1
      if (itrmf .ne. 0) then
         itrmbt = 2
         go to 900
      endif
      TMP = ddot(n,fpls,1,fpls,1)
      fpnrm = dsqrt(TMP)
      ibt = ibt + 1
c ------------------------------------------------------------------------
c If t-condition is met or backtracking is turned off, return. 
c ------------------------------------------------------------------------
      if (fpnrm .le. (1.d0 - t*(1.d0-eta))*fcnrm .or. ibtmax .eq. -1) 
     $ then 
         itrmbt = 0
         go to 900 
      endif
c ------------------------------------------------------------------------
c Otherwise, ... 
c ------------------------------------------------------------------------
      if (ibt .gt. ibtmax) then 
         itrmbt = 1
         go to 900
      endif
c ------------------------------------------------------------------------
c ... choose theta ...
c ------------------------------------------------------------------------
      theta = -(oftjs*redfac)/(fpnrm**2 - fcnrm**2 - 2.d0*oftjs*redfac)
      if(theta .lt. thmin) theta = thmin
      if(theta .gt. thmax) theta = thmax
c ------------------------------------------------------------------------
c ... then reduce the step, increase eta, update redfac ... 
c ------------------------------------------------------------------------
      call dscal(n, theta, step, 1)
      eta = 1.d0 - theta*(1.d0 - eta)
      redfac = theta*redfac
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
      go to 100
c ------------------------------------------------------------------------
c All returns made here.
c -----------------------------------------------------------------------
 900  continue
       ITER_BT = ibt
c ------------------------------------------------------------------------ 
      return
      end
