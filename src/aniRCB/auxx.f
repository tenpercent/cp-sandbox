c ==============================================================
c Initialization routine
c ==============================================================
      Subroutine InitializeRCB (
c ==============================================================
     &           nE, nEmax, nP, nPmax, XYP, IPE,
     &           MaxWi, iW, listtet, tetpmax, iERR)
      implicit none
c ==============================================================
      Integer             nE, nEmax, nP, nPmax, MaxWi, iERR
      Double precision    XYP(3,*)
      Integer             IPE(4,*)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      Integer             iW(*) 

      Integer             iP_pbisE, iP_flag
      Integer             iP_level 
      Integer             iP_fpe, ip_edp
c ==============================================================
      
      If (MaxWi.lt.15*nEmax+nPmax+41) Then
          iERR=1 ! LocalRefine needs at least 15*nEmax+nPmax+41
          return
      end if

      iP_pbisE   =              1
      iP_flag    = iP_pbisE   + 3*nEmax
      iP_level   = iP_flag    + nEmax
      iP_fpe     = iP_level   + nEmax
      ip_edp     = iP_fpe     + 24
  
      call InitializeMeshData (nE, nP,  XYP, IPE, iW(iP_pbisE),
     &        iW(iP_flag), iW(iP_level), iW(iP_fpe), iW(iP_edp),
     &        listtet, tetpmax, iERR)

      Return
      End
c ==============================================================
c ==============================================================
      Subroutine InitializeMeshData (
c ==============================================================
     &           nE, nP, XYP, IPE, pbisE, flag, level, fpe, edp,
     &           listtet, tetpmax, iERR)
      implicit none
c ==============================================================
      Integer             nE, nP
      Double precision    XYP(3,*)
      Integer             IPE(4,*)
      Integer             pbisE(3,*)
      Integer             flag(*)
      Integer             level(*)
      Integer             fpe(4,6), edp(4,4)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      Integer             iERR
      
      Integer             i, i1, i2, j, j1, k, m, l
      Integer             pt1, pt2, pt3, pt4
      Integer             tt1, tt2, tt3
      Integer             k1, k2, k3, kk, kkt1, kkt2, kk1, kk2, fl
      Double precision    s, s1(3)
c ==============================================================
c
c ... define fpe
c
      fpe(1,1) = 1
      fpe(2,1) = 2
      fpe(3,1) = 1
      fpe(4,1) = 2 

      fpe(1,2) = 1
      fpe(2,2) = 3
      fpe(3,2) = 1
      fpe(4,2) = 3 

      fpe(1,3) = 1
      fpe(2,3) = 4
      fpe(3,3) = 2
      fpe(4,3) = 3 

      fpe(1,4) = 2
      fpe(2,4) = 3
      fpe(3,4) = 1
      fpe(4,4) = 4 

      fpe(1,5) = 2
      fpe(2,5) = 4
      fpe(3,5) = 2
      fpe(4,5) = 4 

      fpe(1,6) = 3
      fpe(2,6) = 4
      fpe(3,6) = 3
      fpe(4,6) = 4 
c
c ... define  
c
      edp(1,1) = 0
      edp(2,1) = 1
      edp(3,1) = 2
      edp(4,1) = 3

      edp(1,2) = 1
      edp(2,2) = 0
      edp(3,2) = 4
      edp(4,2) = 5

      edp(1,3) = 2
      edp(2,3) = 4 
      edp(3,3) = 0
      edp(4,3) = 6

      edp(1,4) = 3
      edp(2,4) = 5
      edp(3,4) = 6
      edp(4,4) = 0

c
c ... define level
c
      Do i = 1, nE
        level(i)=0
      End do
c
c ... generate list of tet
c
      Call Initializelisttet (nE, nP, IPE, listtet, tetpmax)

c
c ... define pbisE
c
      Do i = 1, nE
        Do j = 1, 3
          pbisE(j, i) = 0
        End do
        flag(i) = 0
      End do

1     Continue

      fl = 0
      s = 0
      i1 = 0
      j1 = 0

      Do i = 1, nE
        If (pbisE(1,i) .eq. 0) Then
          fl = 1
          Do j = 1, 6
            Call Dist3D (XYP(1, IPE(fpe(1, j), i)),
     &                   XYP(1, IPE(fpe(2, j), i)), s1(1))       
            If (s1(1) .gt. s) Then
               s = s1(1)
               i1 = i
               j1 = j
            End if
          End do
        End if
      End do

      If (fl .eq. 0)  Go to 3

      pt1 = IPE(fpe(1, j1), i1)
      pt2 = IPE(fpe(2, j1), i1)
      Do k1 = 1, listtet(1, pt1)
        tt1 = listtet(k1 + 1, pt1)
        Do k2 = 1, listtet(1, pt2)
          tt2 = listtet(k2 + 1, pt2)
          If (tt1 .ne. tt2) Go to 2
          Do kk = 1, 4
            If (IPE(kk, tt1) .eq. pt1) kk1 = kk 
            If (IPE(kk, tt1) .eq. pt2) kk2 = kk 
          End do
          
          If (pbisE(1, tt1) .eq. 0) Then
            pbisE(1, tt1) = edp (kk1, kk2)
           Else
            pt3 = IPE(fpe(1, pbisE(1, tt1)), tt1)
            pt4 = IPE(fpe(2, pbisE(1, tt1)), tt1)
            If ((pt3 .ne. pt1) .and. (pt3 .ne. pt2)) 
     &            pbisE(3, tt1) = edp (kk1, kk2) 
            If ((pt4 .ne. pt1) .and. (pt4 .ne. pt2)) 
     &            pbisE(2, tt1) = edp (kk1, kk2) 
          End if
2         Continue
        End do
      End do
      
      Go to 1

3     Continue

      Do i = 1, nE
c
        If (pbisE(2,i) .ne. 0) Go to 4

        l = pbisE(1, i)
        pt1 = IPE(fpe(2, l), i)
        pt2 = IPE(fpe(1, l), i)
        pt3 = IPE(fpe(1, 7 - l), i)
        pt4 = IPE(fpe(2, 7 - l), i)
        Call Dist3D (XYP(1, pt2), XYP(1, pt3), s1(1))        
        Call Dist3D (XYP(1, pt2), XYP(1, pt4), s1(2))        
        Call Dist3D (XYP(1, pt3), XYP(1, pt4), s1(3))        

        s = s1(1)
        j = 1
        Do k = 2, 3
          If (s1(k) .gt. s) Then
            s = s1(k)
            j = k
          End if
        End do

        If (j . eq. 1) Then
          kk1 = fpe(1, l)
          kk2 = fpe(1, 7 -l)
         Else    
          If (j . eq. 2) Then
            kk1 = fpe(1, l)
            kk2 = fpe(2, 7 -l)
           Else 
            kk1 = fpe(1, 7 -l)
            kk2 = fpe(2, 7 -l)
          End if  
        End if

        pbise(2,i) = edp(kk1, kk2)
        
        Do k1 = 1, listtet(1, pt2)
          tt1 = listtet(k1 + 1, pt2)
          If (tt1 .ne. i) Then 
            Do k2 = 1, listtet(1, pt3)
              tt2 = listtet(k2 + 1, pt3)
              If (tt2 .ne. i) Then
              
                Do k3 = 1, listtet(1, pt4)
                  tt3 = listtet(k3 + 1, pt4)
                  If ((tt1 .eq. tt2) .and. (tt1 .eq. tt3)) Then
                    Do k = 1, 4
                      If ((IPE(k, tt1) .ne. pt2) .and.
     &                    (IPE(k, tt1) .ne. pt3) .and.                
     &                    (IPE(k, tt1) .ne. pt4)) kk = k
                      If (IPE(k, tt1) .eq. IPE(kk1, i))  kkt1 = k
                      If (IPE(k, tt1) .eq. IPE(kk2, i))  kkt2 = k
                    End do
                    If (fpe(1,pbisE(1,tt1)) .eq. kk) Then
                      pbisE(3,tt1) = edp(kkt1, kkt2)
                      Go to 4
                    End if
                    If (fpe(2,pbisE(1,tt1)) .eq. kk) Then
                      pbisE(2,tt1) = edp(kkt1, kkt2)
                      Go to 4
                    End if
                    iERR = 1011
                     call errMes(iERR, 'InitializeMeshData',
     &                'Mistake in initial data')
                  End if
                End do
              End if
            End do
          End if
        End do  
c
4       Continue
        If (pbisE(3,i) .ne. 0) Go to 5
        l = pbisE(1, i)
        pt1 = IPE(fpe(1, l), i)
        pt2 = IPE(fpe(2, l), i)
        pt3 = IPE(fpe(1, 7 - l), i)
        pt4 = IPE(fpe(2, 7 - l), i)
        Call Dist3D (XYP(1, pt2), XYP(1, pt3), s1(1))        
        Call Dist3D (XYP(1, pt2), XYP(1, pt4), s1(2))        
        Call Dist3D (XYP(1, pt3), XYP(1, pt4), s1(3))        

        s = s1(1)
        j = 1
        Do k = 2, 3
          If (s1(k) .gt. s) Then
            s = s1(k)
            j = k
          End if
        End do

        If (j . eq. 1) Then
          kk1 = fpe(2, l)
          kk2 = fpe(1, 7 -l)
         Else    
          If (j . eq. 2) Then
            kk1 = fpe(2, l)
            kk2 = fpe(2, 7 -l)
           Else 
            kk1 = fpe(1, 7 -l)
            kk2 = fpe(2, 7 -l)
          End if  
        End if

        pbise(3,i) = edp(kk1, kk2)
        
        Do k1 = 1, listtet(1, pt2)
          tt1 = listtet(k1 + 1, pt2)
          If (tt1 .ne. i) Then 
            Do k2 = 1, listtet(1, pt3)
              tt2 = listtet(k2 + 1, pt3)
              If (tt2 .ne. i) Then
              
                Do k3 = 1, listtet(1, pt4)
                  tt3 = listtet(k3 + 1, pt4)
                  If ((tt1 .eq. tt2) .and. (tt1 .eq. tt3)) Then
                    Do k = 1, 4
                      If ((IPE(k, tt1) .ne. pt2) .and.
     &                    (IPE(k, tt1) .ne. pt3) .and.                
     &                    (IPE(k, tt1) .ne. pt4)) kk = k
                      If (IPE(k, tt1) .eq. IPE(kk1, i))  kkt1 = k
                      If (IPE(k, tt1) .eq. IPE(kk2, i))  kkt2 = k
                    End do
                    If (fpe(1,pbisE(1,tt1)) .eq. kk) Then
                      pbisE(3,tt1) = edp(kkt1, kkt2)
                      Go to 5
                    End if
                    If (fpe(2,pbisE(1,tt1)) .eq. kk) Then
                      pbisE(2,tt1) = edp(kkt1, kkt2)
                      Go to 5
                    End if
                    iERR = 1011
                     call errMes(iERR, 'InitializeMeshData',
     &                'Mistake in initial data')
                  End if
                End do
              End if
            End do
          End if
        End do  
c
5       Continue
      End do

      Return
      End
c ==============================================================
      subroutine Ordertriple (a,b,iERR)
c ==============================================================
c ==============================================================
      implicit none
      
      Integer      a(3), b(3), iERR
      Integer      c(3), d(2), i, j, k
c ==============================================================

       If (a(1) .eq. a(2)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if    
       If (a(1) .eq. a(3)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if    
       If (a(2) .eq. a(3)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if 

       Do i = 1, 3 
         c(i) = a(i)
       End do

       k = 1
       Do i = 2, 3
         If (c(i) .lt. c(k)) k = i
       End do

       j = 0
       Do i = 1, 3
         If (i .ne. k) Then
             j = j + 1
             d(j) = i
         End if
       End do

       b(1) = c(k)
       If (c(d(1)) .lt. c(d(2))) Then
          b(2) = c(d(1))  
          b(3) = c(d(2))
        Else
          b(2) = c(d(2))  
          b(3) = c(d(1))
       End if

      Return
      End
c ==============================================================
c Addition auxiliary procedure.
c ==============================================================
      subroutine OrientandMakeNeighbors (
c ==============================================================
     &           nP, nPmax, nF, nFmax, nE, nEmax, 
     &           IPE, IEF, IPF, labelF, listtet, tetpmax,
     &           iERR)
      implicit none
c ==============================================================
      Integer             nP, nPmax, nF, nFmax, nE, nEmax
      Integer             IPE(4,*),IEF(4,*),IPF(3,*)
      Integer             labelF(*)
      Integer             tetpmax
      Integer             listtet((tetpmax + 1), *)
      Integer             iERR

      Integer             pf(3)
      Integer             i,j,i1,i2,i3
      Integer             k, m, maxlabelF
c ==============================================================
c
c  Check of the orientation of boundary faces 
c

       iERR = 0
       Do i = 1, nF
         call Ordertriple (IPF(1,i),IPF(1,i),iERR)
       end do 
   
c
c ... generate list of tetetrahedra
c
       Call Initializelisttet (nE, nP, IPE, listtet, tetpmax)

c        
c ... form IEF
c
       maxlabelF=0
       Do i = 1, nF
         j = labelF(i)
         If (j .gt. maxlabelF) maxlabelF = j
       End do
       maxlabelF = maxlabelF +1

       Do i = 1, nE
         Do j = 1, 4
           IEF(j,i)=-maxlabelF
         End do
       End do

       Do i = 1, nE
         Do j = 1, 4
           If (IEF(j,i) .eq. -maxlabelF) Then

             m = 0
             Do k = 1, 4
               If (k .ne. j) Then
                  m = m + 1
                  pf(m) =  IPE(k,i)
               End if
             End do

             Do i1 = 1, listtet(1, pf(1))
               If (listtet(i1 + 1, pf(1)) .ne. i) Then
                 Do i2 = 1, listtet(1, pf(2))
                   If (listtet(i2 + 1, pf(2)) .ne. i) Then
                     Do i3 = 1, listtet(1, pf(3))
                       If ((listtet(i1 + 1, pf(1)) .eq.
     &                      listtet(i2 + 1, pf(2))) .and.
     &                     (listtet(i1 + 1, pf(1)) .eq.
     &                      listtet(i3 + 1, pf(3)))) Then
                             IEF(j,i) = listtet(i1 + 1, pf(1))
                             Go to 10
                       End if                      
                     End do
                   End if    
                 End do
               End if
             End do
           
             call Ordertriple (pf,pf,iERR)
             
             Do k = 1, nF
               If ((pf(1) .eq. IPF(1,k)) .and.
     &             (pf(2) .eq. IPF(2,k)) .and.
     &             (pf(3) .eq. IPF(3,k))) Then 
                    IEF(j,i) = - labelF(k)
                    Go to 10
               End if
             End do
             iERR = 1011
             call errMes(iERR, 'Init_mesh_gen',
     &                'Mistake in the bound date')
           End if
10         Continue

         End do
       End do 
c


      Return
      End
      
c ==============================================================
c Recover boundary edges.
c ==============================================================
      Subroutine RestoreBndEdges (
     &           nF, nFmax, nE, 
     &           IPE, IEF, IPF, labelF, iERR)
      implicit none
c ==============================================================
      Integer             nF, nFmax, nE
      Integer             IPE(4,*), IEF(4,*)
      Integer             IPF(3,*), labelF(*),  iERR

      Integer             i,j,k,m,m1
c ==============================================================

      nF = 0

      Do i = 1, nE
        Do j =1, 4
        
          k = IEF(j,i)
          If (k .le. 0) Then
            nF = nF + 1
            If (nF .gt. nFmax) Then 
             iERR = 1004
             call errMes(iERR, 'Post_procesing','nFmax is too small')
            End if
            
            m1 = 0
            Do m = 1, 4
              If (m .ne. j) Then
                 m1 = m1 + 1
                 IPF (m1, nF) = IPE(m, i)
              End if
            End do
            labelF(nF) = -k            
            call Ordertriple (IPF(1,nF),IPF(1,nF),iERR)
          End if        
          
        End do
      End do
        
      Return
      End
c ==============================================================
c Codering the number of the point tet
c ==============================================================
      Subroutine Code (k, fl1, fl2)
c ==============================================================
      implicit none
      
      Integer        k
      Logical        fl1, fl2
c ==============================================================
      
      If (k .ge. 3) Then
        fl1 = .true.
       Else
        fl1 = .false.
      End if

      If ((k .eq. 2) .or. (k .eq. 4)) Then
        fl2 = .true.
       Else
        fl2 = .false.
      End if

      Return
      End
c      
c ==============================================================
c Uncodering the number of the point tet
c ==============================================================
      Subroutine Uncode (fl1, fl2, k)
c ==============================================================
      implicit none
      
      Logical        fl1, fl2
      Integer        k
c ==============================================================
     
      If ((.not. fl1) .and. (.not. fl2)) k = 1
      If ((.not. fl1) .and. (      fl2)) k = 2
      If ((      fl1) .and. (.not. fl2)) k = 3
      If ((      fl1) .and. (      fl2)) k = 4

      Return
      End
c      
c ==============================================================
c ==============================================================
      Subroutine Initializelisttet (
c ==============================================================
     &           nE, nP, IPE, listtet, tetpmax)
      implicit none
c ==============================================================
      Integer             nE, nP
      Integer             IPE(4,*)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      
      Integer             i, j, k, m
c ==============================================================
       Do i = 1, nP
         listtet(1, i) = 0
       End do

       Do i = 1, nE
         Do j = 1, 4
           k = IPE(j, i)
           listtet(1, k) = listtet(1, k) + 1
           m = listtet(1, k)
           listtet(m + 1, k) = i
         End do
       End do


      Return
      End
c ==============================================================
      subroutine Dist3D (pt1,pt2,res)
c ==============================================================
c ==============================================================
      implicit none
      
      Double precision  pt1(3), pt2(3), res
c ==============================================================

       res = dsqrt((pt1(1) - pt2(1)) **2 +
     &             (pt1(2) - pt2(2)) **2 +
     &             (pt1(3) - pt2(3)) **2)
                         
      Return
      End
      

