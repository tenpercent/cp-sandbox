C ================================================================
C The program reads a system in CSR format, initialize ILU2 preconditioner,
C and solve the system by preconditioned BiCGstab method
C ================================================================
      Program  Test
      implicit none

C Maximum order and maximum number of non-zeros
      Integer   maxn,maxnz 
      Parameter(maxn=100000, maxnz=1000000)

C Arrays for matrix  kept in CSR format
      Integer ia(maxn+1), ja(maxnz)
      Real*8   a(maxnz), f(maxn), u(maxn)

C Work arrays
      Integer   MaxWr,MaxWi
      Parameter(MaxWr=5*maxnz, MaxWi=6*maxnz)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)

C BiCGStab data
      External  matvec, prevec2
      Integer   ITER, INFO, NUNIT
      Real*8    RESID

C ILU data
      Real*8   tau1,tau2,partlur,partlurout
      Integer  verb, ierr, UsedWr, UsedWi

C External routines
      Real*8   ddot
      External ddot,dcopy    ! from blas library

C Local variables
      Integer  imatvec(1), iprevec(1), ipBCG
      Real*8   resinit
      Integer  n,i,j,nz

C ================================================================
C  STAGE 1: read system in CSR format
C ================================================================

      open(10,file='CSRsystem',status='OLD',ERR=1000)
        read(10,*,ERR=1001) n
        if (n.gt.maxn) then
           write(*,*)'increase maxn to', n
           stop
        end if
        read(10,*,ERR=1001) (ia(i),i=1,n+1)
        nz=ia(n+1)-1
        if (nz.gt.maxnz) then
           write(*,*)'increase maxnz to', nz
           stop
        end if
        read(10,*,ERR=1001) (ja(j),j=1,nz)
        read(10,*,ERR=1001) (a(j),j=1,nz)
        read(10,*,ERR=1001) (f(i),i=1,n)
      close(10)


C ================================================================
C  STAGE 2: initialization of the preconditioner
C ================================================================
      verb=0 ! verbose no
      tau1 =1d-2
      tau2 =1d-3
      partlur = 0.5
      ierr = 0

      call iluoo (n, ia, ja, a, tau1, tau2, verb,
     &   rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &   UsedWr, UsedWi, ierr)
      if (ierr.ne.0) then
         write(*,*)'initialization of  iluoo failed, ierr=',ierr
         stop
      end if
      write (*,*) ' Recommendation for ILU2 + BiCGs:'
      write (*,*) ' partlur=',partlurout
      write (*,*) ' MaxWr=',UsedWr+8*n
      write (*,*) ' MaxWi=',UsedWi

      if (UsedWr+8*n.gt.MaxWr) then
          write(*,*) 'Increase MaxWr to ',UsedWr+8*n
          stop
      end if
      ipBCG  = UsedWr + 1


C ================================================================
C  STAGE 3: set initial guess and compute initial residual
C ================================================================
c  set initial guess to 0
      call dcopy(n,0d0,0,u,1)
c  compute initial residual norm
      resinit=ddot(n,f,1,f,1)
      resinit=dsqrt(resinit)
      if (resinit.eq.0d0) stop 'rhs=0, nothing to solve!'

C ================================================================
C  STAGE 4: iterative solution
C ================================================================
      ITER = 1000             ! max number of iterations
      RESID = 1d-8 * resinit  ! threshold for \|RESID\|
      INFO  = 0               ! no troubles on imput
      NUNIT = 6               ! output channel
      iprevec(1) = n          ! single entry required: system size 
      imatvec(1) = n          ! single entry required: system size 


      call slpbcgs(
     >  prevec2, iprevec, iW,rW,
     >  matvec,  imatvec, ia,ja,a,
     >  rW(ipBCG), n, 8,
     >  n, f, u,
     >  ITER, RESID,
     >  INFO, NUNIT )

      if (INFO.ne.0) stop 'BiCGStab failed'

      Stop

 1000 Continue
      Stop 'Cannot open file CSRsystem'
 1001 Continue
      Stop 'Corrupted data in CSRsystem'

      End

