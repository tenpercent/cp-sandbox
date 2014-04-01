c=======================================================
c Solve linear system
c=======================================================
      Integer Function LSolver(isolv,iprec,N,IA,JA,A,b,x)
      implicit none
c ======================================================================
c method for linear solver and preconditioner
c solv: 1-BiCGstab, 2-GMRES
c prec: 1-ilu(0), 2-ilu(1), 3-ilu(2),
c       4-ilut(10,0.01), 5-iluoo(0.01,0.001)
c ======================================================================
      Integer  isolv, iprec
      Integer  solv, prec
      Character*12 methname,precname
c ======================================================================
c linear system in CSR format
c ======================================================================
      Integer  N
      Integer  IA(*), JA(*) 
      Real*8   A(*), b(*), x(*)
c ======================================================================
c work memory
c ======================================================================
      include 'mmax.fd'
      Integer   MaxWr, MaxWi
c ... Parameter(MaxWr = 29 * (namax + 8 * nvmax) * 1)
c ... old value:
c ... Parameter(MaxWr = 2100000000)
c ... new value:
      Parameter(MaxWr = 21000000)

c ... Parameter(MaxWi = 18 * (namax + 2 * nvmax + 1) * 2)
c ... old value:
c ... Parameter(MaxWr = 2100000000)
c ... new value:
      Parameter(MaxWi = 21000000)
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)
c ### Gmres
      INTEGER   IRESTART
      PARAMETER (IRESTART = 50)
      INTEGER   MW, NW, MH, NH
c ### MW maximal length of solution vector
      PARAMETER (MW = 4 * nvmax, NW = (IRESTART + 1) + 2)
      PARAMETER (MH = IRESTART + 1, NH = IRESTART + 6 )
      REAL*8    H (MH, NH)
c ### Work memory
      real*8    work (MW, NW)
      integer   iworkmax
c ... old value:
c ... PARAMETER (iworkmax = 2147000000)
c ... new value:
      PARAMETER (iworkmax = 21470000)
      integer   iwork (iworkmax)
c External routines
      Real*8   ddot
      External ddot,dcopy    ! from blas library
      External  matvec,matveck,prevec,prevec0,prevec1,prevec2,gdsum
      External  prevec0alu,prevec2alu
      Integer   igdsum
      Integer   ITER, INFO, NUNIT
      Real*8    RESID
c ILU0 data
      Integer   ierr, ipaLU, ipjLU, ipjU, ipiw
c L and U factors
      integer   nzlumax
      PARAMETER (nzlumax=4*NAMAX*5)
      integer   jlu(nzlumax),ju(MW+1)
      real*8    alu(nzlumax)
C ILU data
      Real*8   tau1,tau2,partlur,partlurout
      Integer  verb, UsedWr, UsedWi

c Local variables
      Integer  imatvec(1), iprevec(1), ipBCG
      Real*8   resinit, zero
      real*8    droptol
c Permutation -------------------------------------------------------
      integer   perm(MW),inperm(MW)
c Local--------------------------------------------------------------
      integer   i,j
      integer   nwk,lfil
      real      etime
      real      t(2)
      real      t1,t2
      Logical LSFEX          ! for file with specific linear solver params

      solv = isolv
      prec = iprec
      if (solv.eq.0)  solv = 2
      if (prec.eq.0)  prec = 5
      methname = ''
      precname = ''
c Check of memory
      if (4*MW+nzlumax.ge.iworkmax) then
           write (*,*) 'increase iworkmax to', 4*MW+nzlumax
           stop
      end if     

c Permutations (identical)
      do i=1,N
           perm(i) = i
      end do
      do i=1,n
       inperm(perm(i))=i
      end do 


C ================================================================
C  STAGE 1: initialization of the preconditioner
C ================================================================
      t1 = etime(t)
      goto (1,2,3,4,5) prec
c     
c Init ILU(0)
c preprocess the matrix
1     continue
      precname='ilu(0)'
      call iluk0 (n, a, ja, ia, alu, jlu, ju,
     >           perm,inperm,work,iwork,ierr)
      if (ierr.eq.0) then
         write(*,*)'Initialization of ILU(0) is ok'
      else
         write(*,*)'Initialization of ILU(0) is failed'
         stop
      end if
      goto 9
c Init ILU(1)
c preprocess the matrix
2     continue
      lfil = 1
      nwk=nzlumax
      call iluk(n,a,ja,ia,lfil,alu,jlu,ju,nwk,
     >          perm,inperm,work,iwork,ierr)
      write(*,*) ' nnz for a =', ia(n+1) - ia(1)
      write(*,*) ' nnz for ilu =', jlu(n+1) -jlu(1) + n
      if (ierr.eq.0) then
         write(*,'(a,i2,a)')'Initialization of ILU(',lfil,') is ok'
      else
         write(*,'(a,i2,a)')'Initialization of ILU(',lfil,') is failed'
         stop
      end if
      goto 9
c Init ILU(2)
c preprocess the matrix
3     continue
      lfil = 2
      nwk=nzlumax
      call iluk(n,a,ja,ia,lfil,alu,jlu,ju,nwk,
     >          perm,inperm,work,iwork,ierr)
      write(*,*) ' nnz for a =', ia(n+1) - ia(1)
      write(*,*) ' nnz for ilu =', jlu(n+1) -jlu(1) + n
      if (ierr.eq.0) then
         write(*,'(a,i2,a)')'Initialization of ILU(',lfil,') is ok'
      else
         write(*,'(a,i2,a)')'Initialization of ILU(',lfil,') is failed'
         stop
      end if
      goto 9
c Init ILUt
c preprocess the matrix
4     continue
c ... check and load linear solver parameters
      inquire(file="solver.txt",exist=LSFEX) ! check the presence of solver.txt
      if (LSFEX) then                        ! solver.txt exists, load it 
	open(12, file='solver.txt')
	read(12,*) ! method, preconditioner
	read(12,*) lfil, droptol
	close(12)
      else                                   ! otherwise use default settings
	lfil = 10
	droptol = 0.01
      end if
      nwk=nzlumax
      call ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,nwk,
     &          perm,inperm,work,iwork,ierr)
      write(*,*) ' nnz for a =', ia(n+1) - ia(1)
      write(*,*) ' nnz for ilu =', jlu(n+1) -jlu(1) + n
      if (ierr.eq.0) then
         write(*,'(a,i2,a)')'Initialization of ILUt is ok'
      else
         write(*,'(a,i2,a)')'Initialization of ILUt is failed'
         stop
      end if
      goto 9
      
c Init ILUoo
c preprocess the matrix
5     continue
c ... check and load linear solver parameters
      inquire(file="solver.txt",exist=LSFEX) ! check the presence of solver.txt
      if (LSFEX) then                        ! solver.txt exists, load it 
	open(12, file='solver.txt')
	read(12,*) ! method, preconditioner
	read(12,*) ! lfil, droptol for ilut
	read(12,*) tau1, tau2
	close(12)
      else                                   ! otherwise use default settings
	tau1 =1d-2
	tau2 =1d-3
      end if
      verb=0 ! verbose no
      partlur = 0.5
      ierr = 0
      call iluoo (n, ia, ja, a, tau1, tau2, verb,
     &   rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &   UsedWr, UsedWi, ierr)
      if (ierr.ne.0) then
         write(*,*)'initialization of  iluoo failed, ierr=',ierr
         stop
      end if
      write (*,*) ' Recommendation for ILU2:'
      write (*,*) ' partlur=',partlurout
      write (*,*) ' MaxWr=',UsedWr
      write (*,*) ' MaxWi=',UsedWi

      if (UsedWr.gt.MaxWr) then
          write(*,*) 'Increase MaxWr to ',UsedWr
          stop
      end if
      ipBCG  = UsedWr + 1
      goto 9
      
9     continue
      t1 = etime(t) - t1
      write(*,*)'Preconditioner initialization took ',t1,' s'
C ================================================================
C  STAGE 2: set initial guess and compute initial residual
C ================================================================
c  set initial guess to 0
      zero = 0d0
      call dcopy(N,zero,0,x,1)
c  compute initial residual norm
      resinit=ddot(N,b,1,b,1)
      resinit=dsqrt(resinit)
      if (resinit.eq.0d0) then
	 LSolver = 0
	 write(*,*)'rhs=0, nothing to solve!'
	 return
      end if

C ================================================================
C  STAGE 3: iterative solution
C ================================================================
      write(*,*) 'solving system'
      ITER = 5000             ! max number of iterations
      RESID = 1d-12 * resinit  ! threshold for \|RESID\|
      INFO  = 0               ! no troubles on input
      NUNIT = 6               ! output channel
      iprevec(1) = N          ! single entry required: system size
      imatvec(1) = N          ! single entry required: system size

      t2 = etime(t)
      if (solv.eq.1) then ! BiCGstab
	if (prec.eq.5) then
	  call slpbcgs(
     >      prevec2, iprevec, iW,rW,
     >      matvec,  imatvec, IA,JA,A,
     >      rW(ipBCG), N, 8,
     >      N, b, x,
     >      ITER, RESID,
     >      INFO, NUNIT )
	else
	  write(*,*)'Unsupported preconditioner. Stop'
	  stop
	end if
      else if (solv.eq.2) then ! GMRES
	if (1.le.prec.and.prec.le.3) then
	  call slgmres(
     >    gdsum,   igdsum,
     >    prevec,  iprevec, alu,jlu,ju,perm,inperm,
     >    matveck, imatvec, a,ia,ja,
     >    WORK, MW, NW, H, MH, NH,
     >    N, b, x,
     >    ITER, RESID,
     >    INFO, NUNIT )
	else  if (prec.eq.4) then
	  call slgmres(
     >    gdsum,  igdsum,
     >    prevec1, iprevec, alu,jlu,ju,perm,inperm,
     >    matveck, IMATVEC, a,ia,ja,
     >    WORK, MW, NW, H, MH, NH,
     >    N, b, x,
     >    ITER, RESID,
     >    INFO, NUNIT )
	else  if (prec.eq.5) then
	  call slgmres(
     >    gdsum,  igdsum,
     >    prevec2alu, iprevec, rW,iW,ju,perm,inperm,
     >    matveck, IMATVEC, a,ia,ja,
     >    WORK, MW, NW, H, MH, NH,
     >    N, b, x,
     >    ITER, RESID,
     >    INFO, NUNIT )
	else
	  write(*,*)'Unknown preconditioner. Stop'
	  stop
	end if
      else
         write(*,*) 'Unknown solver. Stop'
         stop
      end if
      t2 = etime(t) - t2

      if (INFO.ne.0) then
	 LSolver = INFO
	 write(*,*)'Solver failed, info = ', INFO
      else
         LSolver = 0
         write(*,*)'Linear solver terminated successfully'
      end if
      write(*,*)'Iterations ',t2,' s'
      write(*,*)'Total time ',t1+t2,' s'

      Return
      End



      Integer Function LSolve(N,IA,JA,A,b,x)
      implicit none
c ======================================================================
c linear system in CSR format
c ======================================================================
      Integer  N
      Integer  IA(*), JA(*) 
      Real*8   A(*), b(*), x(*)
      integer LSolver

      LSolve = LSolver(0,0,N,IA,JA,A,b,x)
      End


c ======================================================================
c wrappers
c ======================================================================
      SUBROUTINE prevec2alu( iprevec, dummy, x, y, dwork, iwork,
     &                       dummyju, dummyperm, dummyinperm)
      INTEGER   iprevec(*), dummy, iwork(*), n
      REAL*8    x(*),y(*),dwork(*)
      INTEGER   dummyju(*), dummyperm(*), dummyinperm(*)

      call prevec2(iprevec, dummy, x, y, iwork, dwork)
      return
      end

      SUBROUTINE prevec0alu( iprevec, dummy, x, y, dwork, iwork,
     &                       dummyju, dummyperm, dummyinperm)
      INTEGER   iprevec(*), dummy, iwork(*), n
      REAL*8    x(*),y(*),dwork(*)
      INTEGER   dummyju(*), dummyperm(*), dummyinperm(*)

      call prevec0(iprevec, dummy, x, y, iwork, dwork)
      return

      end
