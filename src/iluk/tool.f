      subroutine check_resid(n,a,ia,ja,rhs,sol,work)
      implicit none
c
c Input
c
      REAL*8   a(*),rhs(*),sol(*), work(*)
      integer  ia(*),ja(*),n
c
c Local
c
      integer  i
      REAL*8   resid
c
c Multiply A*sol
c
      call matveck (n, 1d0, sol, 0d0, work, a,ja,ia)
c
c Calculation if resid
c
      resid=0d0
      do i=1,n
        resid = resid +(work(i)-rhs(i))**2
      end do
      write (*,'(a,E16.9)') ' True residual of solution = ',dsqrt(resid)

      end
c-----------------------------------------------------------------------
      subroutine dvperm (n, x, perm)
      integer n, perm(n)
      real*8 x(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of a real vector x
c according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c       x(perm(j)) :== x(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n     = length of vector x.
c perm  = integer array of length n containing the permutation  array.
c x     = input vector
c
c on return:
c----------
c x     = vector x permuted according to x(perm(*)) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      real*8 tmp, tmp1
c
      init      = 1
      tmp       = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1        = x(ii)
      x(ii)     = tmp
      next        = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp       = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
c
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue
c
      return
c-------------------end-of-dvperm---------------------------------------
c-----------------------------------------------------------------------
      end


c-----------------------------------------------------------------------
      logical function idnan(d)
      double precision d, x
      integer          ix(2), is, inf
      parameter        (inf = '7FF00000'x)
      equivalence      (x,ix)

      x=d
      is=iand (ix(2), '7FFFFFFF'x)
      if (is .ne. inf) then
         idnan = (is .gt. inf)
         return
      endif
      idnan = (ix(1) .ne. 0)

      return
      end

c-----------------------------------------------------------------------
      logical function idinf (d)
      double precision d, x
      integer          ix(2), is, inf
      parameter        (inf = '7FF00000'x)
      equivalence      (x,ix)

      x=d
      is=iand (ix(2), '7FFFFFFF'x)
      idinf = (is .eq. inf .and. ix(1) .eq. 0)

      return
      end


      subroutine loadsys(n,ia,ja,a,rhs,MW,NZMAX)
      integer n,ia(*),ja(*),MW,NZMAX
      double precision a(*), rhs(*)
      open(10,file='a.dat')
        read(10,*)n
        if (n.gt.MW) then
           write(*,*)'increase MW to', n
           stop
        end if
        read(10,*) (ia(j),j=1,n+1)
        nz=ia(n+1)-1
        if (nz.gt.NZMAX) then
           write(*,*)'increase NZMAX to', nz
           stop
        end if
        read(10,*) (ja(j),j=1,nz)
        read(10,*) (a(j),j=1,nz)
        read(10,*)(RHS(i),i=1,n)
      close(10)
      end

