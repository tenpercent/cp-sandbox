      subroutine jvnewton(n, xcur, fcur, 
     $     f, rpar, ipar,        
     $     itask, nfe, v, z,
     $     prevec, IPREVEC,iW,rW, 
     $     rwork1, 
     $     itrmjv)

      implicit none  

      integer n, ipar(*), itask,nfe, itrmjv, IPREVEC(*),iW(*)
      double precision xcur(n), fcur(n), rpar(*), v(n), z(n), 
     $     rwork1(n),   ddot , rW(*)
      external f,  ddot, prevec

c ------------------------------------------------------------------------
c
c This is jvnewton, the routine for controlling evaluation of products 
c J*v or J*P(inverse)*v or P(inverse)*v, where J is the Jacobian of f 
c and P is a right preconditioning operator. 
c
c ------------------------------------------------------------------------
c 
c Explanation: 
c
c  n       = dimension of the problem.
c
c  xcur    = vector of length n, initial guess on input and final 
c            approximate solution on output. 
c
c  fcur    = vector of length n, value of f at xcur. 
c
c  f       = name of user-supplied subroutine for evaluating the function 
c            the zero of which is sought; this routine has the form 
c
c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
c
c            where xcur is the array containing the current x value, fcur 
c            is f(xcur) on output, rpar and ipar are, respectively, real 
c            and integer parameter/work arrays for use by the subroutine,
c            and itrmf is an integer termination flag.  The meaning of
c            itrmf is as follows:
c              0 => normal termination; desired function value calculated.
c              1 => failure to produce f(xcur).
c
c  rpar    = real parameter/work array passed to the f routine. 
c
c  ipar    = integer parameter/work array passed to the f routine. 
c
c  nfe     = number of function evaluations
c
c  itask   = flag for determining which product is produced.
c              0 => z = J*v
c              1 => z = J*P(inverse)*v 
c              2 => z = P(inverse)*v 
c
c  v       = vector to be multiplied. 
c
c  z       = desired product. 
c
c  rwork1  = vector of length n, work array. 
c
c  itrmjv  = termination flag; values have the following meanings: 
c              0 => normal termination; desired product evaluated. 
c              1 => failure to produce J*v.
c              2 => failure to produce P(inverse)*v. 
c
c
c ------------------------------------------------------------------------
c If z = J*v is desired (itask = 0), then copy v into rwork1; if 
c z = J*P(inverse)*v or z = P(inverse)*v is desired (itask = 1,2), 
c then compute P(inverse)*v in rwork1. 
c ------------------------------------------------------------------------
      if (itask .eq. 0) then 
         call dcopy(n, v, 1, rwork1, 1)
      else
        call prevec(IPREVEC,0,v,rwork1,iW,rW)
      endif
c ------------------------------------------------------------------------
c If only z = P(inverse)*v is desired (itask = 2), then copy rwork1 into 
c z and exit.
c ------------------------------------------------------------------------
      if (itask .eq. 2) then 
         call dcopy(n, rwork1, 1, z, 1)
         itrmjv = 0
         go to 900
      endif
c ------------------------------------------------------------------------
c If z = J*v or z = J*P(inverse)*v is desired (itask = 0, 1), then 
c compute J*rwork1 in z by either analytic evaluation (ijacv = 1) or 
c finite-differences (ijacv = 0, -1). 
c ------------------------------------------------------------------------
      call fdnewton(n, xcur, fcur, f, rpar, ipar, nfe,
     $          rwork1, z, itrmjv)
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      return
      end
      
      subroutine fdnewton(n, xcur, fcur, f, rpar, ipar, nfe, 
     $     v, z, 
     $     itrmjv)

      implicit none  

      integer n, ipar(*), itrmjv,nfe
      double precision xcur(n), fcur(n), rpar(*), v(n), z(n), 
     $      ddot
      external f, ddot

c 
      double precision eps, epsmach, temp, TMP
      integer i, itrmf, ncall
c
      double precision dlamch
      external dlamch
c
c ------------------------------------------------------------------------
      data ncall / 0 /
      save ncall, epsmach 
c ------------------------------------------------------------------------
c Set epsmach (machine epsilon) on first call. 
c ------------------------------------------------------------------------
      if (ncall .eq. 0) epsmach = 2.0d0*dlamch( 'e' )
      ncall = 1
c ------------------------------------------------------------------------
c Compute z = J*v by finite-differences: First, set eps = ||v||for later 
c use in computing the difference step; then evaluate the difference 
c formula according to ijacv and ifdord. 
c ------------------------------------------------------------------------
      TMP = ddot(n,v,1,v,1)
      eps = dsqrt(TMP)
      if (eps .eq. 0.d0) then 
         itrmjv = 1
         go to 900
      endif
        TMP = ddot(n,xcur,1,xcur,1)
        temp = dsqrt(TMP)
         eps = dsqrt((1.d0 + temp)*epsmach)/eps
         do 100 i = 1, n
            v(i) = xcur(i) + eps*v(i)
 100     continue
         call f(n, v, z, rpar, ipar, itrmf)
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         do 110 i = 1, n
            z(i) = (z(i) - fcur(i))/eps
 110     continue
         itrmjv = 0 
         go to 900
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      return
      end
