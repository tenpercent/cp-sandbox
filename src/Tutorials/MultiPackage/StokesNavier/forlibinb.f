c ======================================================================
c  This user routine evaluates the preconditioner-vector product for
c  the unknown Jacobian. The precondtioner is the factorized Stokes system.
c
c  INPUT:
c    iprevec - control parameters for the routine (not used)
c    dummy   - different dummy allow to use different preconditioners
c    x       - input vector
c    iwork   - integer data arrays for the preconditioner
c    dwork   - real*8 data arrays for the preconditioner
c
c  OUTPUT:
c    y       - preconditioner-vector product
c ======================================================================
      Subroutine prevec(iprevec, dummy, x, y, iwork, dwork)
c ======================================================================
      Integer    iprevec(*), dummy, iwork(*)
      Real*8     x(*), y(*), dwork(*), infoUMF(90)

      Integer    sys
c ======================================================================
      sys = 0

      Call umf4sol(sys, y, x, iwork, dwork, infoUMF)

      If(infoUMF(1) .lt. 0) Then
         Write(*,*) 'Error occurred in umf4sol in prevec: ', infoUMF(1)
         Stop
      End if

      Return
      End



c ======================================================================
c This user routine evaluates the nonlinear residual for the FEM problem.
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
c ======================================================================
      Subroutine fnlin(n, xcur, fcur, rpar, ipar, itrmf)
c ======================================================================
      Implicit none

      Integer  itrmf, n, ipar(*)
      Real*8   xcur(n), fcur(n), rpar(*)

      Integer  nv,nt,nb,MaxWi,na,nRow,nCol
      Integer  ipvrt,ipRHS,ipA,iplabelP,ipbnd,iptet,iplbb
      Integer  ipmaterial,ipIRE,ipIA,ipJA,ipiW

      Integer  Dbc, status
      EXTERNAL Dbc, fem3Dext
c ======================================================================
      include 'assemble.fd'
c ======================================================================
c  Set local variables from user-supplied parameters.
      nv    = ipar(1)
      nt    = ipar(2)
      nb    = ipar(3)
      na    = ipar(4)
      MaxWi = ipar(5)

c  structure of rpar
      ipvrt     = 1
      ipRHS     = ipvrt + 3*nv
      ipA       = ipRHS + n

c  structure of ipar
      ipIRE    = 6
      iplabelP = ipIRE + 6*nt
      ipbnd    = iplabelP + nv
      iplbb    = ipbnd + 3*nb
      iptet    = iplbb + nb
      ipmaterial=iptet + 4*nt
      ipIA     = ipmaterial + nt
      ipJA     = ipIA + n+1
      ipiW     = ipJA + na

c     symmetric sparse matrix in a 0-based CSC format used in UMFPACK
c     status = IOR(MATRIX_SYMMETRIC, FORMAT_CSC)
c     general sparse matrix in the AMG format (modifed CSR format
c     where the diagonal element goes first)

      status = IOR(MATRIX_GENERAL, FORMAT_AMG)

      Call BilinearFormTemplate(
     &     nv, nb, nt, rpar(ipvrt), ipar(iplabelP), ipar(ipbnd),
     &                 ipar(iplbb), ipar(iptet), ipar(ipmaterial),
     &     fem3Dext, xcur, status,
     &     n,na,ipar(ipIA),ipar(ipJA),rpar(ipA),rpar(ipRHS),nRow,nCol,
     &     MaxWi, ipar(ipiW))


c  Evaluate nonlinear residual.
      Call mulAgen(nRow, ipar(ipIA),ipar(ipJA),rpar(ipA), xcur,fcur)
      Call daxpy(nRow,-1d0,rpar(ipRHS),1,fcur,1)

c  Set  termination flag for success.

      itrmf = 0

      Return
      End



