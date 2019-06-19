      program body_ax

!==========================================
! FDLIB, CFDLAB, BEMLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!==========================================

!---------------------------------------------------
! Axisymmetric streaming (uniform)
! potential flow past a stationary body.
!
! This program solves
! an integral equation of the second kind
! for the disturbance harmonic potential
! over the body contour
! and computes streamlines
!
! The contour of the body in the xy upper half-plane
! consists of a Number of SeGments (NSG)
!
! The segments may be straight lines or circular arcs.
!
! Each segment is discretized into a number of elements
!
!
! Symbols:
! --------
!
! NSG: Number of segments defining the body contour
!
! NE(i): Number of elements on the ith segment
!
! RT(i): Stretch ratio of elements on ith segment
!
! Itp(i): Index for shape of the ith segment:
!         1 for a straight segment
!         2 for a circular arc
!
! (Xe, Ye):  end-nodes of elements on a segment
! (Xm, Ym):  mid-nodes of elements on a segment
! (Xw, Yw):  end-nodes of elements on all segments
!
! (X0, Y0):  coordinates of collocation points
!
! T0(i): angle subtended from the center of a circular element
!        at ith collocation point
!
! arel(i):  axisymmetric surface area of ith element
!
! phi: disturbance potential at collocation points
!
! dphidn0: normal derivative of disturbance potential
!                 at collocation points
!
! dphids0: tangential derivative of potential
!                 at collocation points
!
! cp:      pressure coefficient at collocation points
!
! NGL: Number of Gaussian points for integration
!      over each element
!
! Icross: Index for stopping the computation of the streamlines
!         Default value is 0
!         If Icross = 1, computation stops when a streamline
!         crosses the yz plane
!
! Notes:
! ------
!
! Normal vector points into the flow
!
!----------------------------
      use precision_mod
      use parameters_mod
      use gauss_legendre, only : gauss_legendre_init
      use common_mod01
      use common_mod02

      Implicit Double Precision (a-h,o-z)

      real(FP_KIND), dimension(MAX_SEGMENTS, MAX_ELEMS) :: phi
      real(FP_KIND), dimension(MAX_DIM)                 :: velt, veln, cp
      real(FP_KIND), dimension(MAX_DIM, MAX_DIM)        :: al       ! for the linear system
      real(FP_KIND), dimension(MAX_DIM)                 :: bl, sol  ! ditto

      real(FP_KIND), dimension (3600)                   :: xstr, ystr ! for streamlines
      character(LEN=64) :: arg

      real(kind=8) t1, t2, t11, t22
      logical :: detailed_iso = .false.

!----------
! constants
!----------

      Null = 0
      None = 1
      Ntwo = 2
      iflow = 1
      ifname = ''
      call cpu_time(t1)

!------
! input
!------

      if (iargc() >= 1) then
          call getarg(1, arg)
          read (arg,*) iflow
          if (iargc() >= 2) then
            call getarg(2, ifname)
            if (iargc() >= 3) then
              detailed_iso = .true.
            end if
          end if
      else
        flow_sel: do while (iflow /= 50 .and. Iflow /= 51)

          write (*,*)
          write (*,*) " Enter: "
          write (*,*)
          write (*,*) " 0 to quit"
          write (*,*) " 50 for flow past a sphere"
          write (*,*) " 51 for flow past a triangular torus"
          write (*,*) "-------------------------------------"

          read (*,*) iflow

          if (iflow == 0) stop

          if (iflow /= 50 .and. Iflow /= 51) then
            write (6,*) "Invalid selection; please try again"
          end if
        end do flow_sel
      end if

!-----------
! initialize
!-----------

!----------------------
!  Compute:
!
!  the element distribution
!  the collocation points
!  the normal derivative of the
!      disturbance potential
!      at the collocation points
!
!----------------------

      call body_ax_geo (Ncl)

!---------------------
! open streamline file
!---------------------

      open (9,file="body_ax.str")
      open (12, file='body_ax_.asy')
      write (12, *) 'import graph;'
      write (12, *) 'size(350,350,IgnoreAspect);'
      write (12, *) 'pair[][] toru = {'

!-------------------------
! Printing session
!
! Print boundary geometry
! to accompany streamlines
!-------------------------

      Ntot = 0
      do j=1,NSG
       Ntot = Ntot+NE(j)+1
      end do
      write (9,102) Ntot

      do j=1,NSG
!      write (9,102) NE(j)+1
!      write (6,102) NE(j)+1
       write (12, *) '  {'
       do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
        write (12, *) '(', XW(j,i), ',' ,YW(j,i), '),'
!       write (6,102) i,XW(j,i),YW(j,i)
       end do
       write (12, *) '  },'
      end do
      write (12, *) '};'

!---
! print the point vortex
!---

      write (9,102) None
      write (9,102) None,Xlvr,Ylvr

!---------------------------
! Display collocation points
! and disturbance velocity
!---------------------------

      write (6,*)
      write (6,*) Ncl," Collocation points "

!     write (6,*)
!     write (6,*) "Collocation points"
!     write (6,*) "and disturbance normal velocity"
!     write (6,*)
!
!     Do i=1,Ncl
!       write (6,102) i,X0(i),Y0(i),dphidn0(i)
!     End Do

!--------
! Prepare
!--------

      call gauss_legendre_init (NGL)

!---------------------------------------------
! Generate the linear system
! for the potential at the collocation points
!
! Generate the influence matrix
! consisting of integrals of the
! single-layer potential
!
! Compute the rhs by evaluating the dlp
!-------------------------------------------------------

      call cpu_time(t11)
      matrix_rhs: do, i=1,Ncl    ! loop over collocation points

       BL(i) = ZERO    ! initialize the right-hand side

!      write (6,*)        " Collocation point :",i

       j = 0              ! counter

       do, k=1,NSG         ! loop over segments

        if (Itp(k) == 2) then
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        end if

        do, L=1,NE(k)       ! loop over elements

         X1 = XW(K,L)
         Y1 = YW(K,L)
         T1 = TW(K,L)

         X2 = XW(K,L+1)
         Y2 = YW(K,L+1)
         T2 = TW(K,L+1)

         j = j+1
         Ising = 0
         if (i == j) Ising = 1

         call body_ax_sdlp &
 
           (X0(i),Y0(i),T0(i) &
           ,X1,Y1,T1 &
           ,X2,Y2,T2 &
           ,NGL &
           ,Ising &
           ,Itp(k) &
           ,rad,xcnt,ycnt &
           ,QQQ &
           ,WWW &
           )

         AL(i,j) = WWW

         BL(i) = BL(i) + QQQ*dphidn0(j)

         end do

        end do

        AL(i,i) = AL(i,i) - HALF

      end do matrix_rhs

      call cpu_time(t22)
      print '(a,f0.3,a)', 'stiffness matrix and rhs construction took: ', t22-t11, ' s.'

!------------------------
! Solve the linear system
!------------------------

      write (6,*) " body_ax: Size of the linear system:",Ncl

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call cpu_time(t11)
      call gel &
 
         (Ncl,AL,BL,SOL &
         ,Isym_g &
         ,Iwlpvt &
         ,Deter &
         ,Istop &
         )

      write (6,*) " body_ax: Linear system solved"
      call cpu_time(t22)
      print '(a,f0.3,a)', 'linear system solution took: ', t22-t11, ' s.'

!-------------------
! Display the system
!-------------------

!     Do I=1,Ncl
!      write (6,101) (AL(i,j),j=1,Ncl),BL(I),SOL(I)
!     End Do

!------------------------
! Distribute the solution
!------------------------

      k = 0        ! counter

      Do i=1,NSG
       Do j=1,NE(i)
        k = k+1
        phi(i,j) = SOL(k)
       End Do
      End Do

!------------------------------------------------
! Compute the tangential disturbance velocity
!             by taking tangential derivative of phi
!
! Compute the pressure coefficient cp
!             drag force
!------------------------------------------------

!     write (6,*)
!     write (6,*) "Tang and normal tot vel; press coeff"
!     write (6,*)

      Forcex = ZERO

      k = 0       ! counter

      Do i=1,NSG
        Do j=1,NE(i)

         k = k+1

!--------------------------------
! tangential disturbance velocity
! computed by differentiating phi
!--------------------------------

         if (j == 1) then           ! forward differencing

             dphids0 = (phi(i,j+1)-phi(i,j))/(s0(k+1)-s0(k))

         else if(j == NE(i)) then  ! backward differencing

             dphids0 = (phi(i,j)-phi(i,j-1))/(s0(k)-s0(k-1))

         else                      ! second-order differencing

             g1 = phi(i,j-1)
             g2 = phi(i,j)
             g3 = phi(i,j+1)
             h1 = s0(k-1)
             h2 = s0(k)
             h3 = s0(k+1)
             aa = ((g3-g2)/(h3-h2)-(g1-g2)/(h1-h2))/(h3-h1)
             bb = (g3-g2)/(h3-h2)-aa*(h3-h2)
             dphids0 = bb

         end if

!---------------
! total velocity
!---------------

         velx = dphidn0(k)*vnx0(k)+dphids0*tnx0(k)
         vely = dphidn0(k)*vny0(k)+dphids0*tny0(k)

         velx = velx + Vx        ! add incident flow

!---
! add line vortex ring
!---

         Iopt = 1

         call lvr_fs &
 
            (Iopt &
            ,X0(k),Y0(k) &
            ,Xlvr,Ylvr &
            ,ulvr,vlvr &
            ,psi &
            )

         velx = velx + cr*ulvr
         vely = vely + cr*vlvr

!-------------------------------
! tangential and normal velocity
! and pressure coefficient
!-------------------------------

         velt(k) = velx*tnx0(k)+vely*tny0(k)
         veln(k) = velx*vnx0(k)+vely*vny0(k)

!------------
! axial force
!------------

         cp(k) = 1.0D0-velt(k)**2/Vx**2

         Forcex = Forcex + cp(k)*vnx0(k)*arel(k)

!        write (6,102) k,velt(k),veln,cp(k)

        end do
      end do

      write (6,*)
      write (6,888) Forcex

!-------------------
! print the solution
!-------------------

      open (3,file="body_ax.out")

!     write (6,*)
!     write (6,*) "x, y, arc length, phi,"
!     write (6,*) "tangential and normal velocity, cp"
!     write (6,*)

      k = 0

      do, i=1,NSG

!       write (6,102) NE(i)
        write (3,102) NE(i)

        do, j=1,NE(i)
         k = k+1
         write (3,102) j,X0(k),Y0(k),s0(k),phi(i,j) &
                       ,velt(k),veln(k),cp(k)
!        write (6,102) j,X0(k),Y0(k),s0(k),phi(i,j)
!    +                 ,velt(k),veln(k),cp(k)
        end do

      end do

      write (3,102) Null

      close (3)

!-----------------------------------------------------------

!     write (6,*)
!     write (6,*) "          MENU"
!     write (6,*)
!     write (6,*) " Enter 0 to quit"
!     write (6,*) "       1 to draw streamlines"
!     write (6,*) " ---------------------------"
!     read  (5,*) menu

!     If(menu == 0) goto 99

!     write (6,*)
!     write (6,*) " Enter maximum number of points before "
!     write (6,*) "                  stopping or inquiring"
!     write (6,*) " --------------------------------------"
!     read  (5,*) Mstr

!     write (6,*)
!     write (6,*) " Enter 0 to stop when number is exceeded"
!     write (6,*) "       1 to inquire for continuation"
!     write (6,*) " --------------------------------------"
!     read  (5,*) Isc

!     write (6,*)
!     write (6,*) "Integration method"
!     write (6,*)
!     write (6,*) " Enter 0 to quit"
!     write (6,*) "       2 for the second-order Runge-Kutta"
!     write (6,*) "       4 for the fourth-order Runge-Kutta"
!     write (6,*) " ----------------------------------------"
!     read  (5,*) IRK

!     If(IRK == 0) goto 99

!     write (6,*)
!     write (6,*) " Enter the spatial step"
!     write (6,*) " ----------------------"
!     read  (5,*) Dl
!
!-----------------------

      menu  = 1
      isc   = 0
      mstr  = 1200
      if (detailed_iso) then
          Dl    = 0.05_FP_KIND * 0.125_FP_KIND
          irk   = 4   ! Runge-Kutta 4th order
      else
          Dl    = 0.05_FP_KIND
          irk   = 2
      end if

!--------------------------
! prepare for crossing test
!--------------------------

      if (iflow == 50) then          ! sphere
          icross = 1
      else if (iflow == 51) then     ! torus
          icross = 0
      end if

!--------------------------
! begin drawing streamlines
!--------------------------

      write (6,*)
      read  (4,*)
      write (12, *) 'pair[][] stream_lines = {'

  22  Continue

      read (4,*) X00,Y00

!------------------------------
! will stop if X00=99 or Y00=99
!------------------------------

      if ((abs(X00-99) < 0.0000001) .or. (abs(Y00-99) < 0.0000001)) then
        close (4)
        goto 97
      end if

      Xcross = X00   ! to be used for crossing check

!------

      L = 1     ! local counter for inquiry
      K = 1     ! total counter

  20  Continue

      xstr(L) = X00
      ystr(L) = Y00

!---------------
! integrate ODEs
!---------------

      call velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux1,Uy1 &
         )

!     write (6,104) L,X00,Y00,Ux1,Uy1
!     print *, 'L: ', L, ' X00: ', X00, ' Y00: ', Y00, ' Ux1: ', Ux1, ' Uy1: ', Uy1

      step = Dl/sqrt(Ux1**2+Uy1**2)     ! set the frozen-time step

      Xsv = X00  ! save
      Ysv = Y00  ! save

!----------------------
      If(IRK == 2) then
!----------------------

        steph = 0.5*step
        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1

        call velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux2,Uy2 &
         )

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)

!---------------------------
      Else If(IRK == 4) then
!---------------------------

        steph = 0.5*step
        step6 = step/6.0

        X00 = Xsv + steph * Ux1
        Y00 = Ysv + steph * Uy1

        call velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux2,Uy2 &
         )

        X00 = Xsv + steph * Ux2
        Y00 = Ysv + steph * Uy2

        call velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux3,Uy3 &
         )

        X00 = Xsv + step * Ux3
        Y00 = Ysv + step * Uy3

        call velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux4,Uy4 &
         )

        X00 = Xsv + step/6.0D0 * (Ux1+2.0*Ux2+2.0*Ux3+Ux4)
        Y00 = Ysv + step/6.0D0 * (Uy1+2.0*Uy2+2.0*Uy3+Uy4)

!-----------
      end if
!-----------

      K = K+1
      L = L+1

!---------------------------
! test for x=0 plane crossing
!---------------------------

      If(Icross == 1) then
       test = Xcross*X00
       If(test.lt.0) then
        write (6,*) " Crossed the x=0 plane: I will stop"
        goto 21
       End If
      End If

!-------------------------
! test for sphere crossing
!-------------------------

      If(Iflow == 50) then

      crosss = sqrt((X00-xcenter)**2  + (Y00-ycenter)**2)

      if (crosss < actis(1)) then
        L = L-1
        goto 23
      end if

      End If

!-----------------------
! window crossing checks
!-----------------------

      If(X00.lt.xwmin) goto 21
      If(X00.gt.xwmax) goto 21
      If(Y00.lt.ywmin) goto 21
      If(Y00.gt.ywmax) goto 21

!---------------------
! point capacity check
!---------------------

      If(K < mstr) goto 20    ! continue the streamline

      If(Isc == 0)  goto 21    ! abandon this streamline

      K = 1    ! reset local counter

      write (6,*)
      write (6,*) "Continue this streamline ?"
      write (6,*)
      write (6,*) "Enter 0 for no, 1 for yes"
      write (6,*) "---------------------------"
      read  (5,*) Icon

      if (Icon == 1) goto 20

!--------------------
! End of a streamline
!--------------------

  21  Continue

      xstr(L) = X00
      ystr(L) = Y00

  23  Continue

!-------------------------
! Print out the streamline
!-------------------------

      write (6,*) " One streamline with ",L," points completed"

      write (9,*) L  ! number of points

      write (12, *) '  {'
      do I=1,L
        write (9,104) I,xstr(I),ystr(I)
        write (12, *) '(', xstr(i), ',' ,ystr(i), '),'
      end do
      write (12, *) '  },'

      goto 22

!--------------------------
! End of streamline drawing
!--------------------------

  97  Continue
      write (12, *) '};'
      close(12)

      write (9,104) Null
      close (9)

!-----
! Done
!-----

      call cpu_time(t2)
      print '(a,f0.3,a)', 'calculation took: ', t2-t1, ' s.'

  99  Continue

 101  Format (20(1x,f7.3))
 102  Format (1x,i5,20(1x,f8.5))
 104  Format (1x,i5,20(1x,f9.5))
 150  Format (1X,10(1X,f10.5))
 888  Format (1X," Axial Force: ",f10.5)

      Stop
      End
