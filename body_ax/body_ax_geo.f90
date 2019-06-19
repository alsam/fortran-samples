      subroutine body_ax_geo (ncl)

!==========================================
! FDLIB, CFDLAB, BEMLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!==========================================

!---------------------------------------------
! Element distribution
! and computation of collocation points
! and disturbance normal velocity
!
! SYMBOLS
! -------
!
! NSG: Number of segments
!
! NE(i): Number of elements on ith segment
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
! CAPACITY
! --------
!
!  10 segments
! 128 elements per segment
!
!-------------------------------------------------

      use precision_mod
      use parameters_mod
      use common_mod01
      use common_mod02
      implicit none

      ! parameters
      integer, intent(out) :: ncl

      ! internal variables
      real(FP_KIND), dimension(MAX_ELEMS + 1) :: xe, ye, te, se
      real(FP_KIND), dimension(MAX_ELEMS)     :: xm, ym, sm
      real(FP_KIND), dimension(MAX_SEGMENTS)  :: rt

      integer       :: i, ic, iopt, isym, itry
      real(FP_KIND) :: angle, cr_new, ddl, ddx, ddy, dth, psi, rad, sinit, ulvr, vlvr
      real(FP_KIND) :: xfirst, yfirst, xsecond, ysecond, xthird, ythird, ycenter_new

!--------------
! common blocks
!--------------

      common/gr1/cr_new,ycenter_new,Itry   ! graphics

!--------
!  SPHERE
!--------

      if (iflow == 50) then

      if (ifname /= '') then
        open (4,file=ifname)
      else
        open (4,file='sphere.dat')
      end if

      read (4,*) NGL
      read (4,*) rad
      read (4,*) xcenter
      read (4,*) Vx       ! velocity of indident flow
      read (4,*) cr       ! line vortex ring strength
      read (4,*)
      read (4,*) NE(1)    ! one segment consisting of arc elms
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

!----------
! over-ride
!----------

      if (Itry > 1) then               ! graphics
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      end if                           ! graphics

!--------------------------------
! place the lvr inside the center
!--------------------------------

      Xlvr = xcenter
      Ylvr = HALF*rad

!---
! one semi-circular contour
! with evenly distributed arcs
!---

      ycenter = ZERO        ! sphere center is on the x axis

      dth = pi/NE(1)

      Do i=1,NE(1)+1       ! points arranged in the countecl sense
        angle = (i-ONE)*dth
        te(i) = angle
        se(i) = angle*rad
        xe(i) = xcenter+rad*cos(angle)
        ye(i) =         rad*sin(angle)
      End Do

!---
! prepare
!---

      NSG = 1           ! one circular segment
      actis(1) = rad
      xcntr(1) = xcenter
      ycntr(1) = ycenter

!---
! semicircular contour
!---

      Itp(1) = 2    ! circular arcs

      Do i=1,NE(1)+1
       tw(1,i) = te(i)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
      End Do

!---
! Collocation points
!---

      ncl = NE(1)

      Do i=1,NE(1)        ! collocation points

        t0(i) = HALF*(te(i)+te(i+1))
        x0(i) = xcenter+rad*Dcos(t0(i))
        y0(i) = ycenter+rad*Dsin(t0(i))
        s0(i) = HALF*(se(i)+se(i+1))

        arel(i) = dth*rad*pi2*y0(i)
        tnx0(i) =-sin(t0(i))
        tny0(i) = cos(t0(i))
        vnx0(i) = tny0(i)
        vny0(i) =-tnx0(i)

        dphidn0(i) = -Vx*vnx0(i)

        Iopt = 1

        call lvr_fs &
          (Iopt &
          ,x0(i),y0(i) &
          ,Xlvr,Ylvr &
          ,ulvr,vlvr &
          ,psi &
          )

        dphidn0(i) = dphidn0(i)-cr*(ulvr*vnx0(i)+vlvr*vny0(i))

      end do

!------------------------------------
! Flow past a triangular torus
!
! Important:
!
! Segments should be arranged in the
! counterclockwise sense
!
! Normal vector points into the fluid
!------------------------------------

      else if (iflow == 51) then

      if (ifname /= '') then
        open (4,file=ifname)
      else
        open (4,file='torus_trgl.dat')
      end if

      read (4,*) NGL
      read (4,*) xfirst,yfirst     ! first vertex
      read (4,*) xsecond,ysecond   ! second vertex
      read (4,*) xthird,ythird     ! third vertex
      read (4,*) Vx
      read (4,*) cr                ! line vortex ring strength
      read (4,*)
      read (4,*) NE(1),RT(1)
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

!----------
! over-ride
!----------

      if (Itry > 1) then
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      end if

!------------------------------
! place the lvr at the centroid
! of the triangle
!------------------------------

      Xlvr = (xfirst+xsecond+xthird)/THREE
      Ylvr = (yfirst+ysecond+ythird)/THREE

!-------------
! preparations
!-------------

      NSG   = 3
      Ic    = 0       ! collocation point counter
      sinit = ZERO    ! initialize arc length

!---
! side # 1
!---

      Itp(1) = 1    ! straight segment
      Isym   = 1

      call elm_line &
 
         (NE(1) &
         ,RT(1) &
         ,xfirst,yfirst &
         ,xsecond,ysecond &
         ,sinit &
         ,Isym &
         ,Xe,Ye,se &
         ,Xm,Ym,sm &
         )

      Do i=1,NE(1)+1
       XW(1,i) = Xe(i)
       YW(1,i) = Ye(i)
      End Do

!---
! collocation points
!---

      do i=1,NE(1)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx**2+ddy**2)

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        Iopt = 1

        call lvr_fs &
 
           (Iopt &
           ,x0(Ic),y0(Ic) &
           ,Xlvr,Ylvr &
           ,ulvr,vlvr &
           ,psi &
           )

        dphidn0(Ic) = dphidn0(Ic) &
                    -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))

      end do

      sinit = se(NE(1)+1)

!---
! side # 2
!---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line &
 
         (NE(2) &
         ,RT(2) &
         ,xsecond,ysecond &
         ,xthird,ythird &
         ,sinit &
         ,Isym &
         ,Xe,Ye,se &
         ,Xm,Ym,sm &
         )

      do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
      end do

!---
! collocation points
!---

      do i=1,NE(2)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx**2+ddy**2)

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        call lvr_fs &
 
          (Iopt &
          ,x0(Ic),y0(Ic) &
          ,Xlvr,Ylvr &
          ,ulvr,vlvr &
          ,psi &
          )

        dphidn0(Ic) = dphidn0(Ic) &
                    -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))

      end do

      sinit = se(NE(2)+1)

!---
! side # 3
!---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line &
 
         (NE(3) &
         ,RT(3) &
         ,xthird,ythird &
         ,xfirst,yfirst &
         ,sinit &
         ,Isym &
         ,Xe,Ye,se &
         ,Xm,Ym,sm &
         )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
      End Do

!---
! collocation points
!---

      Do i=1,NE(3)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx**2+ddy**2)

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        call lvr_fs &
 
         (Iopt &
         ,x0(Ic),y0(Ic) &
         ,Xlvr,Ylvr &
         ,ulvr,vlvr &
         ,psi &
         )

        dphidn0(Ic) = dphidn0(Ic) &
                    -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))

      End Do

      ncl = Ic          ! number of collocation points

!-----------
      End If       ! End of geometry module
!-----------

!-----
! Done
!-----

      Return
      End
