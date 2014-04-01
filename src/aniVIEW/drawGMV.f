C =====================================================================
      Subroutine GMVmesh(fname, Chan,
     &                   Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh in a GMV file.
C =====================================================================
      Implicit None

      Real*8        vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ',nvrt
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)')(vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ',ntet

      Do i = 1, ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(2,i3), vrt(2,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3E24.15)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVscalarVrt(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and nodal solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1, Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3E24.15)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'variable'
      Write(Chan,'(A)') 'solution 1'
      Do i = 1, Nvrt
         Write(Chan,'(E24.15)') U(i)
      End do

      Write(Chan,'(A)') 'endvars'
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVscalarVrtMat(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet,lbe, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and nodal solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*), lbe(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

      Integer       cmap(100)

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1, 100
         cmap(i) = 0
      End do

      Do i = 1, Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if

         j = max(1, lbe(i))
         j = min(lbe(i), 100)
         cmap(j) = cmap(j) + 1
      End do

      j = 0
      Do i = 1, 100
         If (cmap(i).GT.0) Then
            j = j + 1
            cmap(i) = j
         End if
      End do

      Write(Chan,*)
      Write(Chan,*) 'material ', j, ' 0'
      Do i = 1, 100
         If (cmap(i).GT.0) Then
            If(i.LT.10) Then
               Write(10,'(A,I1)') 'mat', i
            Else If(i.LT.100) Then
               Write(10,'(A,I2)') 'mat', i
            Else If(i.LT.1000) Then
               Write(10,'(A,I3)') 'mat', i
            End if
         End if
      End do
      Write(Chan,*)  (cmap(min(max(1,lbE(i)),100)), i = 1, Ntet)

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3E24.15)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'variable'
      Write(Chan,'(A)') 'solution 1'
      Do i = 1, Nvrt
         Write(Chan,'(E24.15)') U(i)
      End do

      Write(Chan,'(A)') 'endvars'
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVscalarTet(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and cellwise solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1, ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, Nbnd
         Write(Chan,'(2I6,3E24.15)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'variable'
      Write(Chan,'(A)') 'solution 0'
      Do i = 1, Ntet
         Write(Chan,'(E24.15)') U(i)
      End do

      Write(Chan,'(A)') 'endvars'
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVvectorTet(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and cellwise vector solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(3,*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1,Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3E24.15)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3E24.15)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'velocity 0'
      Write(Chan,'(10000000E24.15)') (U(1,i),i=1,Ntet)
      Write(Chan,'(10000000E24.15)') (U(2,i),i=1,Ntet)
      Write(Chan,'(10000000E24.15)') (U(3,i),i=1,Ntet)

      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End


      
C =====================================================================
      Real*8 Function DETofTET(v1, v2, v3, v4)
C =====================================================================
      Real*8  v1(3), v2(3), v3(3), v4(3)
    
      Real*8  dx12,dy12,dz12, dx13,dy13,dz13, dx14,dy14,dz14
C =====================================================================

      dx12 = v2(1) - v1(1)
      dy12 = v2(2) - v1(2)
      dz12 = v2(3) - v1(3)
      dx13 = v3(1) - v1(1)
      dy13 = v3(2) - v1(2)
      dz13 = v3(3) - v1(3)
      dx14 = v4(1) - v1(1)
      dy14 = v4(2) - v1(2)
      dz14 = v4(3) - v1(3)

      DETofTET = dx12*dy13*dz14 + dy12*dz13*dx14 + dz12*dx13*dy14
     &         - dx14*dy13*dz12 - dy14*dz13*dx12 - dz14*dx13*dy12

      Return
      End



      
C =====================================================================
      Subroutine VTKscalarVrtMat(U, fname, Chan,
     &                        header, scalarname,
     &                        Nvrt,vrt, Ntet,tet,lbe, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and nodal solution U in a VTK file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*), lbe(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName, header, scalarname

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') '# vtk DataFile Version 3.0'
      Write(Chan,'(A)') header
      Write(Chan,'(A)') 'ASCII'
      Write(Chan,'(A)') 'DATASET UNSTRUCTURED_GRID'

      Write(Chan,*) 'POINTS ', Nvrt, ' double'
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'CELLS ', Ntet+nbnd, 5*Ntet+4*nbnd

      Do i = 1, Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         If(det.GT.0) Then
            Write(Chan,*) 4, i1-1, i2-1, i4-1, i3-1
         Else
            Write(Chan,*) 4, i1-1, i2-1, i3-1, i4-1
         End if
      End do

      Do i = 1, nbnd
         Write(Chan,*) 3, (bnd(j,i)-1, j=1,3)
      End do

      Write(Chan,*) 'CELL_TYPES ', Ntet+nbnd

      Write(Chan,*) (10, i=1, Ntet)
      Write(Chan,*) (5, i=1, nbnd)

      Write(Chan,*) 'CELL_DATA ', Ntet+nbnd
      Write(Chan,*) 'SCALARS material int 1'
      Write(Chan,*) 'LOOKUP_TABLE default'

      Write(Chan,*) (lbe(i), i=1, Ntet)
      Write(Chan,*) (lbbnd(i), i=1, nbnd)

      Write(Chan,*) 'POINT_DATA ', Nvrt
      Write(Chan,*) 'SCALARS ',scalarname,' double 1'
      Write(Chan,*) 'LOOKUP_TABLE default'

      Do i = 1, Nvrt
         Write(Chan,'(E24.15)') U(i)
      End do

      Close(Chan)

      Return
      End


C =====================================================================
      Subroutine VTKSpecial(fname, Chan,
     &                      header,
     &                      Nvrt,vrt, Ntet,tet,lbe, Nbnd,bnd,lbbnd,
     &                      U, Uname,
     &                      Ur, Urname,
     &                      Ui, Uiname,
     &                      Er, Ername,
     &                      Ei, Einame,
     &                      Jr, Jrname,
     &                      Ji, Jiname,
     &                      Sr, Srname,
     &                      Si, Siname)
C =====================================================================
C Routine saves mesh and nodal solution U in a VTK file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*), lbe(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName, header, Uname, Urname, Uiname
      Character*(*) Ername, Einame, Jrname, Jiname, Srname, Siname
      Real*8        Ur(*), Ui(*), Er(3,*), Ei(3,*), Jr(3,*), Ji(3,*)
      Real*8        Sr(*), Si(*)

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') '# vtk DataFile Version 3.0'
      Write(Chan,'(A)') header
      Write(Chan,'(A)') 'ASCII'
      Write(Chan,'(A)') 'DATASET UNSTRUCTURED_GRID'

      Write(Chan,*) 'POINTS ', Nvrt, ' double'
      Do i = 1, nvrt
         Write(Chan,'(3E24.15)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'CELLS ', Ntet, 5*Ntet

      Do i = 1, Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))

         If(det.GT.0) Then
            Write(Chan,*) 4, i1-1, i2-1, i4-1, i3-1
         Else
            Write(Chan,*) 4, i1-1, i2-1, i3-1, i4-1
         End if
      End do

c      Do i = 1, nbnd
c         Write(Chan,*) 3, (bnd(j,i)-1, j=1,3)
c      End do

      Write(Chan,*) 'CELL_TYPES ', Ntet

      Write(Chan,*) (10, i=1, Ntet)
c      Write(Chan,*) (5, i=1, nbnd)

      Write(Chan,*) 'CELL_DATA ', Ntet
      Write(Chan,*) 'SCALARS material int 1'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Write(Chan,*) (lbe(i), i=1, Ntet)
c      Write(Chan,*) (lbbnd(i), i=1, nbnd)


      Write(Chan,*) 'VECTORS ', Ername,' double'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') (Er(j,i), j=1,3)
      End do
      Write(Chan,*) 'VECTORS ', Einame,' double'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') (Ei(j,i), j=1,3)
      End do
      Write(Chan,*) 'VECTORS ', Jrname,' double'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') (Jr(j,i), j=1,3)
      End do
      Write(Chan,*) 'VECTORS ', Jiname,' double'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') (Ji(j,i), j=1,3)
      End do
      Write(Chan,*) 'SCALARS ', Srname,' double 1'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') Sr(i)
      End do
      Write(Chan,*) 'SCALARS ', Siname,' double 1'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Do i = 1, Ntet
         Write(Chan,'(3E24.15)') Si(i)
      End do

      Write(Chan,*) 'POINT_DATA ', Nvrt
      Write(Chan,*) 'SCALARS ', Uname,' double 2'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Do i = 1, Nvrt
         Write(Chan,'(2E24.15)') U(i), U(i+Nvrt)
      End do

      Write(Chan,*) 'SCALARS ', Urname,' double 1'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Do i = 1, Nvrt
         Write(Chan,'(E24.15)') Ur(i)
      End do

      Write(Chan,*) 'SCALARS ', Uiname,' double 1'
      Write(Chan,*) 'LOOKUP_TABLE default'
      Do i = 1, Nvrt
         Write(Chan,'(E24.15)') Ui(i+Nvrt)
      End do


      Close(Chan)

      Return
      End
