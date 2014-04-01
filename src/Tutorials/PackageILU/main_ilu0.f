C ================================================================
C The program reads a system in CSR format, initialize ILU0 preconditioner,
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

C Work arrays keep ILU factors and 8 BCG vectors 
      Integer   MaxWr,MaxWi
      Parameter(MaxWr=maxnz+8*maxn, MaxWi=maxnz+2*maxn+1)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)

C BiCGStab data
      External  matvec, prevec0
      Integer   ITER, INFO, NUNIT
      Real*8    RESID

C ILU0 data
      Integer   ierr, ipaLU, ipjLU, ipjU, ipiw

C External routines
      Real*8   ddot
      External ddot,dcopy    ! from blas library

C Local variables
      Integer  imatvec(1), iprevec(1), ipBCG
      Real*8   resinit, zero
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
      ipaLU=1
      ipBCG=ipaLU+nz
      ipjU =1
      ipjLU=ipjU+n+1
      ipiw =ipjLU+nz !work array of length n
      call ilu0(n,a,ja,ia, rW(ipaLU),iW(ipjLU),iW(ipjU),iW(ipiw),ierr)
      if (ierr.ne.0) then
         write(*,*)'initialization of  ilu0 failed, zero pivot=',ierr
         stop
      end if
c Keep data in rW and iW up to rW(nz) and iW(nz+n+1) !



C ================================================================
C  STAGE 3: set initial guess and compute initial residual
C ================================================================
c  set initial guess to 0
      zero = 0d0
      call dcopy(n,zero,0,u,1)
c  compute initial residual norm
      resinit=ddot(n,f,1,f,1)
      resinit=dsqrt(resinit)
      if (resinit.eq.0d0) stop 'rhs=0, nothing to solve!'

C ================================================================
C  STAGE 4: iterative solution
C ================================================================
      ITER = 1000             ! max number of iterations
      RESID = 1d-8 * resinit  ! threshold for \|RESID\|
      INFO  = 0               ! no troubles on input
      NUNIT = 6               ! output channel
      iprevec(1) = n          ! single entry required: system size 
      imatvec(1) = n          ! single entry required: system size 
     

      call slpbcgs(
     >  prevec0, iprevec, iW,rW,
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

