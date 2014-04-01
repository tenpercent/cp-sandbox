C=========================================================================
C  The routine generates the metric from a given continuous piecewise 
C  linear function u and the mesh (tet, vrt, bnd).
C  The algorithm is based on the ZZ recovery method for interior nodes.
C=========================================================================
      Subroutine Nodal2MetricZZ(u, 
     &                          vrt,Nvrt, tet,Ntet, bnd,Nbnd, metric, 
     &                          Nrmem,rmem, Nimem,imem)
C=========================================================================
C  Input: 
C     u    - function defined at mesh vertices
C     Nvrt - the number of vertices
C     vrt  - coordinates of these vertices
C     Ntet - the number of tetrahedra
C     tet  - the connecticity table of tetrahedra
C     Nbnd - the number of boundary faces
C     Bnd  - the list of boundary faces
C
C  Output: 
C     metric(6,*) - tensor metric
C
C  Work arrays: rmem - real*8 array of length Nrmem
C               imem - integer array of length Nimem
C=========================================================================
      implicit None

      Real*8    u(*)

C Mesh
      Integer   Nvrt, Ntet, Nbnd
      Real*8    vrt(3, *)
      Integer   tet(4, *), bnd(3, *)

C Metric
      Real*8    metric(6, *)

C Work arrays
      integer   Nrmem, Nimem
      integer   imem(*)
      Real*8    rmem(*)

C Local variables
      Integer   n, ip1, ip2, ip3, ip4
      Integer   igx, igy, igz, ivx, ivy, ivz, iEnd
      Real*8    grad(3), det

C=========================================================================
c ... memory check
      If(Nimem.lt.3*Nvrt + 4*Ntet) Then
         Call errMesLMR(1001, "Nodal2MetricZZ", 
     &                         "Not enough memory for Integer arrays")
      End if

      igx  = Nvrt
      igy  = igx + Ntet 
      igz  = igy + Ntet
      ivx  = igz + Ntet
      ivy  = ivx + Nvrt
      ivz  = ivy + Nvrt
      iEnd = ivz + Nvrt

      If(Nrmem.lt.iEnd) Then
         Call errMesLMR(1002, "Nodal2MetricZZ", 
     &                         "Not enough memory for Real*8 arrays")
      End if


c ... recover the constant gradients g = grad(u)
      Do n = 1, Ntet
         ip1 = tet(1, n)
         ip2 = tet(2, n)
         ip3 = tet(3, n)
         ip4 = tet(4, n)

         Call GradU(vrt(1,ip1), vrt(1,ip2), vrt(1,ip3), vrt(1,ip4),
     &              u(ip1), u(ip2), u(ip3), u(ip4), grad)

         rmem(igx + n) = grad(1)
         rmem(igy + n) = grad(2)
         rmem(igz + n) = grad(3)
      End do


c ... recover the linear gradients (gx,gy,gz) -> (vx,vy,vz)
      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igx+1), rmem(ivx+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igy+1), rmem(ivy+1), Nvrt, Nimem, rmem, imem)


      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igz+1), rmem(ivz+1), Nvrt, Nimem, rmem, imem)


c ... recover the constant gradients g = grad(vx)
      Do n = 1, Ntet
         ip1 = tet(1, n)
         ip2 = tet(2, n)
         ip3 = tet(3, n)
         ip4 = tet(4, n)

         Call GradU(vrt(1,ip1), vrt(1,ip2), vrt(1,ip3), vrt(1,ip4),
     &              rmem(ivx+ip1), rmem(ivx+ip2), 
     &              rmem(ivx+ip3), rmem(ivx+ip4), grad)

         rmem(igx + n) = grad(1)
         rmem(igy + n) = grad(2)
         rmem(igz + n) = grad(3)
      End do


c ... recover the linear gradients (gx,gy,gz) -> (uxx,uxy,uxz)
      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igx+1), rmem(ivx+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igy+1), rmem(ivy+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igz+1), rmem(ivz+1), Nvrt, Nimem, rmem, imem)

      Do n = 1, Nvrt
         metric(1, n) = rmem(ivx + n)
         metric(4, n) = rmem(ivy + n) / 2
         metric(6, n) = rmem(ivz + n) / 2
      End do


c ... recover the constant gradients g = grad(vy)
      Do n = 1, Ntet
         ip1 = tet(1, n)
         ip2 = tet(2, n)
         ip3 = tet(3, n)
         ip4 = tet(4, n)

         Call GradU(vrt(1,ip1), vrt(1,ip2), vrt(1,ip3), vrt(1,ip4),
     &              rmem(ivy+ip1), rmem(ivy+ip2), 
     &              rmem(ivy+ip3), rmem(ivy+ip4), grad)

         rmem(igx + n) = grad(1)
         rmem(igy + n) = grad(2)
         rmem(igz + n) = grad(3)
      End do


c ... recover the linear gradients (gx,gy,gz) -> (uyx,uyy,uyz)
      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igx+1), rmem(ivx+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igy+1), rmem(ivy+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igz+1), rmem(ivz+1), Nvrt, Nimem, rmem, imem)

      Do n = 1, Nvrt
         metric(4, n) = metric(4, n) + rmem(ivx + n) / 2 
         metric(2, n) = rmem(ivy + n)
         metric(5, n) = rmem(ivz + n) / 2
      End do


c ... recover the constant gradients g = grad(vz)
      Do n = 1, Ntet
         ip1 = tet(1, n)
         ip2 = tet(2, n)
         ip3 = tet(3, n)
         ip4 = tet(4, n)

         Call GradU(vrt(1,ip1), vrt(1,ip2), vrt(1,ip3), vrt(1,ip4),
     &              rmem(ivz+ip1), rmem(ivz+ip2), 
     &              rmem(ivz+ip3), rmem(ivz+ip4), grad)

         rmem(igx + n) = grad(1)
         rmem(igy + n) = grad(2)
         rmem(igz + n) = grad(3)
      End do


c ... recover the linear gradients (gx,gy,gz) -> (uzx,uzy,uzz)
      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igx+1), rmem(ivx+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igy+1), rmem(ivy+1), Nvrt, Nimem, rmem, imem)

      Call P02P1(Nvrt, Nbnd, Ntet, vrt, bnd, tet, 
     &           rmem(igz+1), rmem(ivz+1), Nvrt, Nimem, rmem, imem)

      Do n = 1, Nvrt
         metric(6, n) = metric(6, n) + rmem(ivx + n) / 2 
         metric(5, n) = metric(5, n) + rmem(ivy + n) / 2
         metric(3, n) = rmem(ivz + n)
      End do


c ... Make the Hessian elliptic using SpectralModule() from libanimba.a 
      Do n = 1, Nvrt
         Call SpectralModule(metric(1,n), det)
      End do

      Return
      End 


