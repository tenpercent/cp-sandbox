c-------------------------------------------------------------
c    This is the simplest test program which reads a system in CSR format
c    and solves it by UMFPACK package and tests the residual.
c    Since UMFPACK uses sorted CSC format, the matrix is transformed
c    from CSR to sorted CSC by the use of routines from aux_umf.f
c    (Other useful operations may be learned from UMFPACK/Demo/umf4hb.f)
c
c    link:  -liblu.a -lblas 
c    where:  liblu.a = UMFPACK/Demo/umf4_f77wrapper.o + libumfpack.a + libamd.a
c
c  authors: Yuri Vassilevski & Konstantin Lipnikov
c  for questions: yuri.vassilevski@gmail.com
c-------------------------------------------------------------
      PROGRAM TESTUMFPACK
      IMPLICIT NONE

c from: umf4_f77wrapper.f
      EXTERNAL  umf4def, umf4sym, umf4num, umf4sol, umf4fsym, umf4fnum
c from: aux_umf.f
      EXTERNAL  IDSORT, TRANSP

      INTEGER   NMAX,  NZMAX
c maximum order of the system
      PARAMETER ( NMAX = 50000 )
c maximum number of non-zero entries in the matrix
      PARAMETER ( NZMAX= NMAX * 5 )
c system data
      integer   IA(NMAX+1),JA(NZMAX), iwk(NZMAX) ! iwk is needed ONLY to transform CSR to CSC
      real*8    A(NZMAX), RHS(NMAX), SOL(NMAX), RES(NMAX)
c umfpack data
      integer symbolic(2), numeric(2), sys
      real*8  control(20), info(90)
c local
      integer i,j,N,NZ,ierr
c read CSR system
      open(10,file='input_main')
        read(10,*) n,nz
        if (n.gt.NMAX) then
           write(*,*)'increase NMAX to', n
           stop
        end if
        if (nz.gt.NZMAX) then
           write(*,*)'increase NZMAX to', nz
           stop
        end if
        read(10,*) (a(j),j=1,nz)
        read(10,*) (ia(j),j=1,n+1)
        read(10,*) (ja(j),j=1,nz)
        read(10,*) (RHS(j),j=1,n)
      close(10)
      if (nz.ne.ia(n+1)-1) stop 'Inconsistency of data'

c order within each row and transform system to CSC  and make it 0-bazed
      do i=1,n
       j=ia(i+1)-ia(i)
       call IDSORT (JA(ia(i)), A(ia(i)), j, 2)
      end do
      call TRANSP (n,n,a,ja,ia,iwk,ierr)
      if (ierr.ne.0) stop 'transp failed'
      do i = 1, n+1
          IA(i) = IA(i) - 1
      end do
      do j = 1, nz
          JA (j) = JA (j) - 1
      end do
c now the matrix is in 0-bazed CSC format with sorted row indexes.

c set  the default control parameters in the Control array
      call umf4def( control )
c print error messages only
      control (1) = 1

c pre-order and symbolic analysis
      call umf4sym( n,n, ia,ja,a, symbolic,control,info)
c check umf4sym error condition
      if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4sym: ', info (1)
            stop
      endif
c print time of symbolic analysis
      print 80,  info (16)
80    format ('time of symbolic analysis:', e10.2, ' (sec)')

c numeric factorization
      call umf4num( ia,ja,a, symbolic,numeric,control,info)
c check umf4num error condition
      if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4num: ', info (1)
            stop
      endif

c print statistics for the numeric factorization
c call umf4pinf (control, info) could also be done.
      print 90, info (66),
     $    (info (41) * info (4)) / 2**20,
     $    (info (42) * info (4)) / 2**20,
     $    info (43), info (44), info (45)
90    format ('numeric factorization:',/,
     $    '   time:    ', e10.2, /,
     $    '   actual numeric LU statistics:', /,
     $    '   size of LU:    ', f10.2, ' (MB)', /,
     $    '   memory needed: ', f10.2, ' (MB)', /,
     $    '   flop count:    ', e10.2, /
     $    '   nnz (L):       ', f10.0, /
     $    '   nnz (U):       ', f10.0)

c free the symbolic analysis data
      call umf4fsym (symbolic)

c solve Ax=b, without iterative refinement
      sys = 0
      call umf4sol (sys, SOL, RHS, numeric, control, info)
      if (info (1) .lt. 0) then
          print *, 'Error occurred in umf4sol: ', info (1)
          stop
      endif


c free the numeric factorization data
      call umf4fnum (numeric)

c check the residual, the matrix is in zero-based CSC format.  
      call residCSC (n, IA, JA, A,  SOL, RHS, RES )

      end

      
c-------------------------------------------------------------------
c Compute the residual, r = Ax-b, its max-norm, and print the max-norm
C Note that A is zero-based CSC format.
      subroutine residCSC (n,  Ap, Ai, Ax, x, b, r)
      integer          n,  Ap(n+1), Ai(*), j, i, p
      double precision Ax(*), x(n), b(n), r(n), rmax, aij
      do i = 1, n
          r(i) = -b(i)
      end do

      do j = 1,n
          do p = Ap(j) + 1, Ap(j+1)
              i = Ai(p) + 1
              aij = Ax(p)
              r(i) = r(i) + aij * x(j)
          end do
      end do
      rmax = 0
      do i = 1, n
          rmax = max (rmax, r(i))
      end do
      print *, 'norm (A*x-b): ', rmax
      return
      end
