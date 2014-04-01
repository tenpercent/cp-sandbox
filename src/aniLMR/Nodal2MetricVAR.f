C ======================================================================
C  The routine generates the metric "metric' on the basis of given
C  function "u" and the mesh 'tet, vrt, bnd'
C
C  Limitation to the mesh: each boundary node must have at least
C                          one adjacent interior node or one adjacent 
C                          boundary node adjacent to an interior one.
C ======================================================================
      Subroutine Nodal2MetricVAR(U,
     &                           Vrt,Nvrt, Tet,Ntet, Bnd,Nbnd, metric,
     &                           Nrmem,rmem, Nimem,imem)
C ======================================================================
C  Input: 
C     U -    function defined at mesh nodes
C     Nvrt - the number of nodes
C     Vrt  - coords of the nodes
C     Ntet - the number of tetrahedra
C     Tet  - the connectivity table
C     Nbnd - the number of boundary triangles
C     Bnd  - the list of boundary triangles
C
C  Output: 
C     metric - metric to be defined. Is nothing else but the discrete
C              Hessian reduced to an elliptic form
C
C  Work arrays: rmem - d.p., of length  Nrmem
C               imem - integer, of length Nimem
C
C ======================================================================
      implicit None

      double precision  U(*)

C mesh
      integer           Nvrt,Ntet,Nbnd
      double precision  Vrt(3,*)
      integer           Tet(4,*), Bnd(3,*)

C metric to be defined
      double precision  metric(6,*)

C work arrays
      integer           Nrmem,Nimem
      integer           imem(*)
      double precision  rmem(*)

C Local variables
      double precision  DX12,DY12,DZ12,DX13,DY13,DZ13,DX14,DY14,DZ14,DET
      double precision  Dxu,Dyu,Dzu,DxNi,DyNi,DzNi
      integer           i,j,k,m,it,iv,kv,L1,L2,L3,L4, Ni(4), ki,km
      integer           lbuf,kbuf(1000)
      double precision  meas(6)
      double precision  BAS(4,3)
      double precision  A11,A12,A13,A21,A22,A23,A31,A32,A33,B12,B13,B23
      logical           flag

C ======================================================================
c ... memory check
      if (Nrmem.lt.Nvrt.or.Nimem.lt.2*Nvrt+4*Ntet) then
         write(*,*) 'Increase Nrmem to', Nvrt
         write(*,*) 'Increase Nimem to', 2*Nvrt+4*Ntet
         Stop
      end if

c ... gradients of nodal basis functions
      BAS(1,1) = -1.
      BAS(1,2) = -1.
      BAS(1,3) = -1.

      BAS(2,1) = 1.
      BAS(2,2) = 0.
      BAS(2,3) = 0.

      BAS(3,1) = 0.
      BAS(3,2) = 1.
      BAS(3,3) = 0.

      BAS(4,1) = 0.
      BAS(4,2) = 0.
      BAS(4,3) = 1.

c ... initialize BadNode,NhostTet,HostTet
      do iv = 1, Nvrt
c        NhostTet(iv) = 0
         imem(iv+Nvrt) = 0
      end do

c ... generate NhostTet
      do it = 1, Ntet
         do i = 1, 4
c            NhostTet(tet(i,it)) = NhostTet(tet(i,it)) + 1
            iv  =  tet(i,it)
            imem(iv+Nvrt) = imem(iv+Nvrt) + 1
         end do
      end do

c ... pointers to HostTet
      imem( 1 ) = 2*Nvrt + 1
      do iv = 2, Nvrt
         imem( iv ) = imem( iv-1 ) + imem( Nvrt + iv-1 )
      end do

      do iv = 1, Nvrt
         imem( Nvrt + iv ) = 0
      end do

      do it = 1, Ntet
         do i = 1, 4
            iv = tet(i,it)
            k = imem( iv ) + imem( Nvrt + iv )
            imem( k ) = it
            imem( Nvrt + iv ) = imem( Nvrt + iv ) + 1
         end do
      end do

c ... restore BadNode, NhostTet
      do iv = 1, Nvrt
c        BadNode(iv) = 1
c        NhostTet(iv) = 0
         imem(iv) = 0
         imem(iv+Nvrt) = 0
      end do

      do it = 1, Ntet
         do i = 1, 4
c           NhostTet(tet(i,it)) = NhostTet(tet(i,it)) + 1
            iv  =  tet(i,it)
            imem(iv+Nvrt) = imem(iv+Nvrt) + 1
         end do
      end do


c ... mark boundary nodes
      do i = 1, Nbnd
c        NhostTet(bnd(1,i)) = -abs(NhostTet(bnd(1,i)))
         imem(bnd(1,i)+Nvrt) = -abs(imem(bnd(1,i)+Nvrt))
         imem(bnd(1,i)) = 1
         imem(bnd(2,i)+Nvrt) = -abs(imem(bnd(2,i)+Nvrt))
         imem(bnd(2,i)) = 1
         imem(bnd(3,i)+Nvrt) = -abs(imem(bnd(3,i)+Nvrt))
         imem(bnd(3,i)) = 1
      end do

c ... initialize Mass and metric
      do iv = 1, Nvrt
c        Mass(iv) = 0d0
         rmem(iv) = 0d0
         do i = 1, 6
            metric(i,iv) = 0d0
         end do
      end do


c ... recover the Hessian: H_ij  = M^{-1} integral_Om du/dx_i dNi/dx_j d Om
      ki = 1
      do iv = 1, Nvrt
         do i  = 1, iabs(imem(Nvrt+iv))
            if (imem(Nvrt+iv).lt.0) then
               ki = ki + 1
            else
               it = imem(2*Nvrt+ki)
               ki = ki + 1

               L1 = Tet(1,it)
               L2 = Tet(2,it)
               L3 = Tet(3,it)
               L4 = Tet(4,it)

c  ...  nodal basis function
               Ni(1) = 0
               Ni(2) = 0
               Ni(3) = 0
               Ni(4) = 0
               if (iv.eq.L1) Ni(1) = 1d0
               if (iv.eq.L2) Ni(2) = 1d0
               if (iv.eq.L3) Ni(3) = 1d0
               if (iv.eq.L4) Ni(4) = 1d0

               DX12 = VRT(1,L1) - VRT(1,L2)
               DX13 = VRT(1,L1) - VRT(1,L3)
               DX14 = VRT(1,L1) - VRT(1,L4)
               DY12 = VRT(2,L1) - VRT(2,L2)
               DY13 = VRT(2,L1) - VRT(2,L3)
               DY14 = VRT(2,L1) - VRT(2,L4)
               DZ12 = VRT(3,L1) - VRT(3,L2)
               DZ13 = VRT(3,L1) - VRT(3,L3)
               DZ14 = VRT(3,L1) - VRT(3,L4)

               DET  = DX12*DY13*DZ14+DY12*DZ13*DX14+DZ12*DX13*DY14-
     &                DX14*DY13*DZ12-DY14*DZ13*DX12-DZ14*DX13*DY12

c              Mass(iv) = Mass(iv) + DABS(DET)/24d0
               rmem(iv) = rmem(iv) + DABS(DET)/24d0

               A11 =  (DY13*DZ14-DY14*DZ13)/DET
               A12 = -(DX13*DZ14-DX14*DZ13)/DET
               A13 =  (DX13*DY14-DX14*DY13)/DET

               A21 = -(DY12*DZ14-DY14*DZ12)/DET
               A22 =  (DX12*DZ14-DX14*DZ12)/DET
               A23 = -(DX12*DY14-DX14*DY12)/DET

               A31 =  (DY12*DZ13-DY13*DZ12)/DET
               A32 = -(DX12*DZ13-DX13*DZ12)/DET
               A33 =  (DX12*DY13-DX13*DY12)/DET

               B12 = A12
               B13 = A13
               B23 = A23
               A12 = A21
               A13 = A31
               A23 = A32
               A21 = B12
               A31 = B13
               A32 = B23

               Dxu = U(L1)*(A11*BAS(1,1) + A12*BAS(1,2) + A13*BAS(1,3))
     &             + U(L2)*(A11*BAS(2,1) + A12*BAS(2,2) + A13*BAS(2,3))
     &             + U(L3)*(A11*BAS(3,1) + A12*BAS(3,2) + A13*BAS(3,3))
     &             + U(L4)*(A11*BAS(4,1) + A12*BAS(4,2) + A13*BAS(4,3))

               Dyu = U(L1)*(A21*BAS(1,1) + A22*BAS(1,2) + A23*BAS(1,3))
     &             + U(L2)*(A21*BAS(2,1) + A22*BAS(2,2) + A23*BAS(2,3))
     &             + U(L3)*(A21*BAS(3,1) + A22*BAS(3,2) + A23*BAS(3,3))
     &             + U(L4)*(A21*BAS(4,1) + A22*BAS(4,2) + A23*BAS(4,3))

               Dzu = U(L1)*(A31*BAS(1,1) + A32*BAS(1,2) + A33*BAS(1,3))
     &             + U(L2)*(A31*BAS(2,1) + A32*BAS(2,2) + A33*BAS(2,3))
     &             + U(L3)*(A31*BAS(3,1) + A32*BAS(3,2) + A33*BAS(3,3))
     &             + U(L4)*(A31*BAS(4,1) + A32*BAS(4,2) + A33*BAS(4,3))

               DxNi = Ni(1)*(A11*BAS(1,1) + A12*BAS(1,2) + A13*BAS(1,3))
     $              + Ni(2)*(A11*BAS(2,1) + A12*BAS(2,2) + A13*BAS(2,3))
     $              + Ni(3)*(A11*BAS(3,1) + A12*BAS(3,2) + A13*BAS(3,3))
     $              + Ni(4)*(A11*BAS(4,1) + A12*BAS(4,2) + A13*BAS(4,3))

               DyNi = Ni(1)*(A21*BAS(1,1) + A22*BAS(1,2) + A23*BAS(1,3))
     $              + Ni(2)*(A21*BAS(2,1) + A22*BAS(2,2) + A23*BAS(2,3))
     $              + Ni(3)*(A21*BAS(3,1) + A22*BAS(3,2) + A23*BAS(3,3))
     $              + Ni(4)*(A21*BAS(4,1) + A22*BAS(4,2) + A23*BAS(4,3))

               DzNi = Ni(1)*(A31*BAS(1,1) + A32*BAS(1,2) + A33*BAS(1,3))
     $              + Ni(2)*(A31*BAS(2,1) + A32*BAS(2,2) + A33*BAS(2,3))
     $              + Ni(3)*(A31*BAS(3,1) + A32*BAS(3,2) + A33*BAS(3,3))
     $              + Ni(4)*(A31*BAS(4,1) + A32*BAS(4,2) + A33*BAS(4,3))

               metric(1,iv) = metric(1,iv) - DABS(DET)*Dxu*DxNi/6d0
               metric(2,iv) = metric(2,iv) - DABS(DET)*Dyu*DyNi/6d0
               metric(3,iv) = metric(3,iv) - DABS(DET)*Dzu*DzNi/6d0
               metric(4,iv) = metric(4,iv) 
     &                       - DABS(DET)*(Dxu*DyNi+Dyu*DxNi)/12d0
               metric(5,iv) = metric(5,iv)
     &                       - DABS(DET)*(Dyu*DzNi+Dzu*DyNi)/12d0
               metric(6,iv) = metric(6,iv)
     &                       - DABS(DET)*(Dxu*DzNi+Dzu*DxNi)/12d0
            end if
         end do
      end do


c ... inner nodes
      do iv = 1, Nvrt
         if (imem(Nvrt+iv).gt.0) then
            do j = 1, 6
               metric(j,iv) = metric(j,iv) / rmem(iv)
            end do
         end if
      end do


c .. extrapolate the metric from interior nodes to boundary and average
c    extrapolations. Exclude respective nodes from the set of bad ones.
      ki = 1
      do iv = 1, Nvrt
         if (imem(Nvrt+iv).gt.0) then
            do i = 1, imem(Nvrt+iv)
               ki = ki + 1
            end do
         else
            lbuf = 0
            do i  = 1, -imem(Nvrt+iv)
               it = imem(2*Nvrt+ki)
               ki = ki + 1

               do j = 1, 4
                  kv = tet(j,it)
                  flag = .true.

                  do k = 1, lbuf
                     if (kbuf(k).eq.kv) flag = .false.
                  end do

                  if (flag.and.imem(kv).eq.0) then
                     lbuf = lbuf + 1
                     kbuf(lbuf) = kv
                  end if
               end do
            end do

            do k = 1, lbuf
               kv = kbuf(k)
               do j = 1, 6
                  metric(j,iv) = metric(j,iv)
     &                          + metric(j,kv)/lbuf
               end do
               imem(Nvrt+iv) = iabs(imem(Nvrt+iv))
            end do
         end if
      end do

      do iv = 1, Nvrt
         if (imem(iv).ne.0.and.imem(Nvrt+iv).gt.0) then
            imem(Nvrt+iv) = - imem(Nvrt+iv)
            imem(iv)      = 0
         end if
      end do


c ... treat the BadNodes by a secondary extrapolation
      km = 1
      do iv = 1, Nvrt
c        if (BadNode(iv).eq.1) then
         if (imem(iv).eq.1) then
            i = 0
            do k = 1, 6
               meas(k) = 0d0
            end do

c           do m = 1, -NhostTet(iv)
            do m = 1, -imem(Nvrt+iv)
c              it = HostTet(m,iv)
               it = imem(2*Nvrt+km)
               km = km + 1
               do j = 1, 4
c                 if (BadNode(tet(j,it)).eq.0) then
                  if (imem(tet(j,it)).eq.0) then
                     i = i + 1
                     do k = 1, 6
                        meas(k) = meas(k) + metric(k,tet(j,it))
                     end do
                  end if
               end do
            end do

            if (i.eq.0) then
               write(*,*)'two arms rule failed on this mesh... '
               do k = 1, 3
                  meas(k) = 1d0
               end do
c              Stop 
            else
               do k = 1, 6
                  metric(k,iv) = meas(k) / i
               end do
            end if
         else
            do m = 1, iabs(imem(Nvrt+iv))
               km = km + 1
            end do
         end if
      end do


c ... make the Hessian to be elliptic using SpectralModule from libanimba.a
      do iv = 1, Nvrt
         Call SpectralModule(metric(1,iv), det)
      end do

      return
      end

