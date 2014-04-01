C ======================================================================
      Subroutine encodeISYS(iE, iP1, iP2, iP3, iP4, iiR, iiF,
     &                          nP, nR, nF, nE, iSYS)
C ======================================================================
      implicit none
      include 'assemble.fd'

      Integer iE, iP1, iP2, iP3, iP4, iiR(6), iiF(4)
      Integer     nP, nR, nF, nE, i
      Integer iSYS(MAXiSYS)
C ======================================================================
c ... tet number
      iSYS(3) = iE

c ... vertices of this tet
      iSYS(4) = iP1
      iSYS(5) = iP2
      iSYS(6) = iP3
      iSYS(7) = iP4

c ... edges of this tet 
      Do i = 1, 6
         iSYS(i + 7) = iiR(i)
      End do

c ... faces of this tet
      Do i = 1, 4
         iSYS(i + 13) = iiF(i)
      End do

c ... total number of mesh objects
      iSYS(18) = nP
      iSYS(19) = nR
      iSYS(20) = nF
      iSYS(21) = nE

      Return
      End



C ======================================================================
      Subroutine encodeISYSfull(iE, iP1, iP2, iP3, iP4,
     &                          iR1, iR2, iR3, iR4, iR5, iR6, 
     &                          iF1, iF2, iF3, iF4,  
     &                          nP, nR, nF, nE, iSYS)
C ======================================================================
      implicit none
      include 'assemble.fd'

      Integer iE, iP1, iP2, iP3, iP4, iR1, iR2, iR3, iR4, iR5, iR6
      Integer     iF1, iF2, iF3, iF4, nP, nR, nF, nE
      Integer iSYS(MAXiSYS)
C ======================================================================
c ... tet number
      iSYS(3) = iE

c ... vertices of this tet
      iSYS(4) = iP1
      iSYS(5) = iP2
      iSYS(6) = iP3
      iSYS(7) = iP4

c ... edges of this tet 
      iSYS( 8) = iR1
      iSYS( 9) = iR2
      iSYS(10) = iR3
      iSYS(11) = iR4
      iSYS(12) = iR5
      iSYS(13) = iR6

c ... faces of this tet
      iSYS(14) = iF1
      iSYS(15) = iF2
      iSYS(16) = iF3
      iSYS(17) = iF4

c ... total number of mesh objects
      iSYS(18) = nP
      iSYS(19) = nR
      iSYS(20) = nF
      iSYS(21) = nE

      Return
      End

C ======================================================================
      Subroutine decodeISYS(iE, iP1, iP2, iP3, iP4,
     &                          iR1, iR2, iR3, iR4, iR5, iR6, 
     &                          iF1, iF2, iF3, iF4,  
     &                          nP, nR, nF, nE, iSYS)
C ======================================================================
      implicit none
      include 'assemble.fd'

      Integer iE, iP1, iP2, iP3, iP4, iR1, iR2, iR3, iR4, iR5, iR6
      Integer     iF1, iF2, iF3, iF4, nP, nR, nF, nE
      Integer iSYS(MAXiSYS)
C ======================================================================
c ... tetrahedron number
      iE = iSYS(3)

c ... vertices of this tetrahedron
      iP1 = iSYS(4)
      iP2 = iSYS(5)
      iP3 = iSYS(6)
      iP4 = iSYS(7)

c ... edges of this tetrahedron
      iR1 = iSYS( 8)
      iR2 = iSYS( 9)
      iR3 = iSYS(10)
      iR4 = iSYS(11)
      iR5 = iSYS(12)
      iR6 = iSYS(13)

c ... faces of this tetrahedron
      iF1 = iSYS(14)
      iF2 = iSYS(15)
      iF3 = iSYS(16)
      iF4 = iSYS(17)

c ... total number of mesh objects
      nP = iSYS(18)
      nR = iSYS(19)
      nF = iSYS(20)
      nE = iSYS(21)

      Return
      End


