c ======================================================================
c The program reads data for Bratu problem, initializes ILU2 preconditioner,
c and solves Bratu problem using inexact Newton method with inner linear 
c BiCGstab iterations (assuming that the jacobian is not known). The precondtioner
c for BiCGstab is ILU2 preconditioner for the grid Laplacian operator (aniILU package)
c ======================================================================
c  Delta u + d  u'_x + lambda exp(u) = 0   in domain
c                                 u  = 0   on boundary
c       Parameters:
c
c       n       - dimension of the discrete problem 
c       x       - vector of the initial guess and the solution
c
c       nx,ny,nz- number of grid nodes in X,Y and Z directions, n = nx*ny*nz
c       lamdba  - Bratu problem parameter (e.g.1)
c       d       - Bratu problem parameter (e.g.1)
c ======================================================================

      PROGRAM  INB_BRATU
      Implicit None
        
c     Maximum order of problem
      Integer   maxn
      Parameter(maxn = 250000)

c    Work arrays
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 11*maxn + 20*maxn, MaxWi = 30*maxn)

      Double precision rW(MaxWr)
      Integer          iW(MaxWi)
c 11*maxn - real space for inexact Newton method with inner BiCGstab solver
c 20*maxn - real    space for preconditioner and matrix
c 30*maxn - integer space for preconditioner and matrix

c     Bratu problem data
      Integer          nx, ny, nz
      Double precision h, cl,cr, h2l, lambda, d
     
c     Pointers for matrix kept in CSR format
      Integer  iPia, iPja, rPa


c     InexactNewton data
      EXTERNAL fBratu, ddot, prevec2, laplace_csr 
c  fbratu is the vector function of the nonlinear system
c  prevec2 is the preconditioner evaluation for the jacobian
c  make_csr generates a matrix for which the preconditioner is initialized

      Integer          IPREVEC, N, iParBratu(3), INFO(5), rPwork,LenrW
      Double Precision RESID, STPTOL, SOL(maxn), rParBratu(3) 

c
c  INFO(*) is the control array for the driver routine slInexactNewton() of INB package
c
c  Explanation of INFO(5): 
c       on  input, its components  are as follows:
c           INFO(1) <= initial value for successful termination flag (=0)
c           INFO(2) <= maximal number of linear iterations per Newton step 
c           INFO(3) <= maximal number of nonlinear iterations 
c           INFO(4) <= maximal number of backtracks   
c           INFO(5) <= printing level (0 keeps mum)
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
 

c     ILU data
      Double Precision tau1, tau2, partlur, partlurout
      Integer          verb, ierr, UsedWr, UsedWi, iPilu, rPilu
      Integer          MaxWrILU, MaxWiILU

c  Local variables
      Integer          i


C =======================================================
C       STAGE 1:  set parameters of Bratu problem
C =======================================================
      nx = 10
      d  = 1
      lambda = 1
      ny = nx
      nz = nx
      write(*,'(A)') 'Parameters for the Bratu problem:'
      write(*,'(A,I3,2(A,F7.4))') 
     &     'nx = ny = nz =', nx, ', d =', d, ', lambda =', lambda
        
      iParBratu(1) = nx
      iParBratu(2) = ny
      iParBratu(3) = nz
      N = nx*ny*nz
      h = 1.d0/dfloat(nx + 1) 

      cl  = 1.d0 - h*d/2.d0 
      cr  = 1.d0 + h*d/2.d0 
      h2l = h*h*lambda

      rParBratu(1) = cl
      rParBratu(2) = cr
      rParBratu(3) = h2l

      If(N.gt.maxn) Stop 'increase maxn'


c======================================================
c     STAGE 2: Initializing in CSR format the matrix of 
c              discrete Laplacian to be used in 
c              defining the preconditioner        
c =====================================================
      iPia = 1
      iPja = iPia + N+1
      rPa  = 1

      call laplace_csr(nx,ny,nz, iW(iPia), iW(iPja), rW(rPa))


C ======================================================
C     STAGE 3: Initialization of ILU2 preconditioner
C ======================================================
      verb = 0
      tau1 = 5d-2
      tau2 = 5d-3
      partlur = 0.5d0
      ierr = 0

      iPilu = iPja + 7*N
      rPilu = rPa + 7*N

      MaxWrILU = MaxWr - rPilu
      MaxWiILU = MaxWi - iPilu

      Call iluoo(n, iW(iPia), iW(iPja), rW(rPa), tau1, tau2, verb,
     &           rW(rPilu), iW(iPilu), MaxWrILU, MaxWiILU,
     &           partlur,partlurout, UsedWr, UsedWi, ierr)

      If(ierr.ne.0) then
         write(*,*) 'Initialization of iluoo has failed, ierr=', ierr
         stop
      End if

c Shift the preconditioner data to the beginning of work array 
      Do i = 1, UsedWi
         iW(i) = iW(i + iPilu - 1)
      End do

      Do i = 1, UsedWr
         rW(i) = rW(i + rPilu - 1)
      End do


C ======================================================
C     STAGE 4: iterative solution by slInexactNewton
C ======================================================
      IPREVEC = N
      RESID   = 1d-10
      STPTOL  = 1d-7

      INFO(1) = 0    ! initializing successful termination flag for Newton
      INFO(2) = 1000 ! maximal number of linear iterations 
      INFO(3) = 100  ! maximal number of nonlinear iterations 
      INFO(4) = 10   ! maximal number of backtracks   
      INFO(5) = 0    ! print level (0 keeps mum)

      rPwork = rPa + UsedWr
      LenrW  = MaxWr - rPwork + 1

c Initial guess for nonlinear iterations
      Do i = 1, N
         SOL(i) = 0.d0
      End do
        
      Call slInexactNewton(prevec2, IPREVEC, iW, rW,   
     &                     fBratu, rParBratu, iParBratu,
     &                     N, SOL,
     &                     RESID, STPTOL, 
     &                     rW(rPwork), LenrW,
     &                     INFO)
        
      If(INFO(1).ne.0) Then
         Write(*,*) 'Failed to solve the problem INFO=', INFO(1)
      Else
         Write(*,*)
         Write(*,*) '||F(SOL)|| = ',RESID
         Write(*,*) 'Number of linear iterations: ',INFO(2)
         Write(*,*) 'Number of nonlinear iterations: ',INFO(3)
         Write(*,*) 'Number of backtracks: ',INFO(4)
         Write(*,*) 'Number of function evaluations: ',INFO(5)
      End if
        
      Stop
      End



c --------------------------------------------------------------------
c This user routine evaluates a nonlinear residual for the Bratu 
c test problem. 
c
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
      Subroutine fBratu(n, xcur, fcur, rpar, ipar, itrmf)
c --------------------------------------------------------------------
      Implicit none
      Integer          i, itrmf, j, k, j1, j2, n, ipar(*), nx, ny, nz 
      Double precision cl, cr, h2l, xcur(n), fcur(n), rpar(*)
c --------------------------------------------------------------------

c  Set local variables from user-supplied parameters.
      nx  = ipar(1)
      ny  = ipar(2)
      nz  = ipar(3)
      cl  = rpar(1) 
      cr  = rpar(2)
      h2l = rpar(3)

c  Evaluate nonlinear residual.
      Do 100 k = 1, nz 
        Do 110 j = 1, ny 
          j1 = (k - 1)*nx*ny + (j - 1)*nx + 2 
          j2 = (k - 1)*nx*ny + j*nx - 1

          Do 120 i = j1, j2
            fcur(i) = cr*xcur(i+1) + cl*xcur(i-1) - 6.d0*xcur(i) 
     $              + h2l*dexp(xcur(i)) 
 120      Continue

          j1 = j1 - 1
          fcur(j1) = cr*xcur(j1+1) - 6.d0*xcur(j1) + h2l*dexp(xcur(j1)) 

          j2 = j2 + 1 
          fcur(j2) = cl*xcur(j2-1) - 6.d0*xcur(j2) + h2l*dexp(xcur(j2)) 

          If(j.ne.1) then 
            Do 130 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i-nx)
 130        Continue
          Endif

          If(j.ne.ny) then 
            Do 140 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i+nx)
 140        Continue
          Endif

          If(k.ne.1) then 
            Do 150 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i-nx*ny)
 150        Continue
          Endif

          If(k.ne.nz) then 
            Do 160 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i+nx*ny)
 160        Continue
          Endif
 110    Continue 
 100  Continue 

c  Set  termination flag for success.

      itrmf = 0

      Return
      End



c --------------------------------------------------------------------
      Subroutine laplace_csr(nx,ny,nz,ia,ja,a)
c --------------------------------------------------------------------
c This routine generates the 7-point Laplacian matrix and stores it
c in the compressed sparce row format.
c --------------------------------------------------------------------
      Implicit none
      Integer   nx, ny,nz
      Real *8   a(*)
      Integer   ia(*), ja(*)
      Integer   i,j,k,m,l
c --------------------------------------------------------------------
      
      ia(1) = 1
      m = 1
      Do k = 1, nz 
        Do j = 1, ny 
          Do i = 1, nx
             m = m + 1
             ia(m) = ia(m-1)
             ia(m) = ia(m)+1
             if (i.gt.1)  ia(m) = ia(m)+1
             if (i.lt.nx) ia(m) = ia(m)+1
             if (j.gt.1)  ia(m) = ia(m)+1
             if (j.lt.ny) ia(m) = ia(m)+1
             if (k.gt.1)  ia(m) = ia(m)+1
             if (k.lt.nz) ia(m) = ia(m)+1
          End do
        End do
      End do
              
      m = 0
      Do k = 1, nz 
        Do j = 1, ny 
          Do i = 1, nx
             m = m + 1
             l = ia(m)
             a(l) = 6d0
             ja(l) = m  
             if (i.gt.1) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m-1
             end if
             if (i.lt.nx) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m+1
             end if
             if (j.gt.1) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m-nx
             end if
             if (j.lt.ny) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m+nx
             end if
             if (k.gt.1) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m-nx*ny
             end if
             if (k.lt.nz) then
                 l = l+1
                 a(l) = -1d0
                 ja(l) = m+nx*ny
             end if
          End do
        End do
      End do
             
      Return
      End
                                                     

