c ======================================================================
c The program solves a simple system
c        2 x^2 + y^2 + z^2 = 1
c        x^2 + 2 y^2 + z^2 = 1
c        x^2 + y^2 + 2 z^2 = 1
c using inexact Newton method with inner linear BiCGstab iterations.
c We assume that the Jacobian is unknown. 
c The precondtioner is identity matrix.
c ======================================================================
        PROGRAM  INB_SIMPLE
        IMPLICIT NONE
        
c     Maximum order of problem
        integer   maxn
        Parameter(maxn = 3)

c    Work array
        integer   MaxWr
        Parameter(MaxWr = 11*maxn)
        double precision rW(MaxWr)
c 11*maxn - real space for inexact Newton method with inner BiCGstab solver

c     problem data
        integer          N
     
c     InexactNewton data
        EXTERNAL fSimple, ddot, prevecIden

c  fSimple    is the vector function of the nonlinear system
c  prevecIden is the preconditioner evaluation for the jacobian (identity here)

        INTEGER          IPREVEC, iParSimple(1), iDum(1), INFO(5)
        DOUBLE PRECISION RESID, STPTOL, SOL(maxn), rParSimple(1),rDum(1)

c
c  INFO(*) is the control array for the driver routine slInexactNewton() of INB package
c
c  Explanation of INFO(5): 
c       on  input, its components  are as follows:
c           INFO(1) <= initial value for successful termination flag (=0)
c           INFO(2) <= maximal number of linear iterations per Newton step 
c           INFO(3) <= maximal number of nonlinear iterations 
c           INFO(4) <= maximal number of backtracks   
c           INFO(5) <= printing level (0 keeps mum,1 nonlinear residuals,2 linear residuals)
c
c       on output, its components  are as follows:
c
c       INFO(1) - termination flag; its value has the folowing meanings:
c              0 => normal termination: ||F||.le.RESID or ||step||.le.STPTOL.
c              1 => maximum nonlinear iterations reached without success.
c              2 => failure to evaluate F.
c              3 => in jvnewton, J*v failure.
c              4 => in jvnewton, P(inverse)*v failure.
c              5 => in drvnewton, insufficient initial model norm reduction
c                   for adequate progress. NOTE: This can occur for several
c                   reasons; examine itrmks on return from the Krylov
c                   solver for further information. 
c              6 => in btnewton, failure to reach an acceptable step through
c                   backtracking.
c              7 => insufficient work memory for slInexactNewton (should be at least 11*N)
c       INFO(2) - number of linear iterations
c       INFO(3) - number of nonlinear iterations
c       INFO(4) - number of backtracks
c       INFO(5) - number of function evaluations

c .. order of the problem
        N = 3


        IPREVEC = N
        RESID   = 1d-10
        STPTOL  = 1d-7

        INFO(1) = 0   ! initializing successful termination flag for Newton
        INFO(2) = 1000! maximal number of linear iterations 
        INFO(3) = 100 ! maximal number of nonlinear iterations 
        INFO(4) = 10  ! maximal number of backtracks   
        INFO(5) = 1   ! print level 

c Initial guess for nonlinear iterations
        SOL(1) = 0.d0
        SOL(2) = 0.d0
        SOL(3) = 0.d0
        
        call slInexactNewton(prevecIden, IPREVEC, iDum, rDum,   
     $                       fSimple, rParSimple, iParSimple,
     $                       N, SOL,
     $                       RESID, STPTOL, 
     $                       rW, MaxWr,
     $                       INFO)
        
        if(INFO(1).ne.0) then
           write(*,*) 'Failed to solve the problem INFO=',INFO(1)
        else
           write(*,*)
           write(*,*) 'Solution: x =',SOL(1),' y =',SOL(2),' z =',SOL(3)
           write(*,*) '||F(SOL)|| = ',RESID
           write(*,*) 'Number of linear iterations: ',INFO(2)
           write(*,*) 'Number of nonlinear iterations: ',INFO(3)
           write(*,*) 'Number of backtracks: ',INFO(4)
           write(*,*) 'Number of function evaluations: ',INFO(5)
        end if
        
        Stop
        End



c --------------------------------------------------------------------
c This user routine evaluates a nonlinear residual. 
c INPUT:
c    n       dimension of vectors
c    xcur    current vector 
c    rpar    double precision user-supplied parameters
c    ipar    integer user-supplied parameters
c
c OUTPUT:
c    fcur    nonlinear residual vector (zero for the solution)
c    itrmf   flag for successful termination of the routine
c --------------------------------------------------------------------
      subroutine fSimple(n, xcur, fcur, rpar, ipar, itrmf)
c --------------------------------------------------------------------
      implicit none
      integer          itrmf, n, ipar(*)
      double precision xcur(n), fcur(n), rpar(*)
c --------------------------------------------------------------------

      if (n.ne.3) then
          itrmf = 1
          return
      end if

      fcur(1) = 2*xcur(1)**2 + xcur(2)**2 + xcur(3)**2 - 1d0
      fcur(2) = xcur(1)**2 + 2*xcur(2)**2 + xcur(3)**2 - 1d0
      fcur(3) = xcur(1)**2 + xcur(2)**2 + 2*xcur(3)**2 - 1d0

c  Set  termination flag for success.
      itrmf = 0

      Return
      End



c --------------------------------------------------------------------
      Subroutine prevecIden( iprevec, dummy, x, y, iwork, dwork)
c --------------------------------------------------------------------
c This is the identity preconditioner. 
c --------------------------------------------------------------------
      Integer   iprevec(*), dummy, iwork(*)
      Real*8    x(*),y(*),dwork(*)
      integer   n
c --------------------------------------------------------------------
      n = iprevec(1)

      Call dcopy(n,x,1,y,1)

      Return
      End

