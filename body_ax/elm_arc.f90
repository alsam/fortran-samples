      subroutine elm_arc &
 
        (N &
        ,ratio &
        ,Xcnt,Ycnt &
        ,radius &
        ,angle1,angle2 &
        ,sinit &
        ,Isym &
        ,Xe,Ye,Te,se &
        ,Xm,Ym,Tm,sm &
        )

!-----------------------------------------
! FDLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!----------------------------------------

!----------------------------------------------------------------
!  Disretization of a circular arc into N elements
!
!  SYMBOLS:
!  -------
!
!  Xcnt, Ycnt:      center of the arc
!  ratio:      radius of the arc
!
!  ratio: If Isym = 0, ratio of length of LAST to FIRST element
!         If Isym = 1, ratio of length of MID  to FIRST element
!
!  Xe, Ye, Te:      coordinates and polar angle for the end-nodes
!  Xm, Ym, Tm:      coordinates and polar angle for the mid-nodes
!
!  sinit: assigned arc length at the first point
!
!----------------------------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      ! parameters
      integer,                                 intent(in)  :: n, isym
      real(FP_KIND),                           intent(in)  :: ratio, xcnt, ycnt, &
                                                              radius, angle1, angle2, sinit
      real(FP_KIND), dimension(MAX_ELEMS + 1), intent(out) :: xe, ye, te, se
      real(FP_KIND), dimension(MAX_ELEMS),     intent(out) :: xm, ym, tm, sm

      ! internal variables
      integer       :: n1, nh, nh1, i
      real(FP_KIND) :: alpha, deltat, texp, tmp1, tmp2, angleh, factor

!------------
! one element
!------------

      if (N.eq.1) then

       xe(1) = xcnt ! FIXME was X1
       ye(1) = ycnt ! FIXME was Y1
       te(1) = angle1
       se(1) = sinit

       xe(2) = xcnt ! FIXME was X2
       ye(2) = ycnt ! FIXME was Y2
       te(2) = angle2 ! was angle12 - looks like typo
       se(2) = se(1) ! FIXME was '+sqrt( (X2-X1)**2+(Y2-Y1)**2)'

       Go to 99

      end if

!--------
! prepare
!--------

      N1 = N+1

!---------------------------
! non-symmetric distribution
!---------------------------

      If(Isym.eq.0) then

      texp  = 1.0D0/(N-1.0D0)
      alpha = ratio**texp

      if (abs(alpha-1.0D0).gt.0.0000001) then
       factor = (1.0D0-alpha)/(1.0D0-alpha**N)
      else
       factor = 1.0D0/(N-1.0D0+1.0D0)
      end if

      deltat = (angle2-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*cos(te(1))
      Ye(1) = Ycnt + radius*sin(te(1))
      se(1) = sinit

      Do i=2,N1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*cos(Te(i))
        Ye(i) = Ycnt + radius*sin(Te(i))
        se(i) = se(i-1)+radius*abs(deltat)
        deltat = deltat*alpha
      End Do

      Go to 99

      End If

!-----------------------
! symmetric distribution
!----------------------

!---
! even number of elements
!---

      If(mod(N,2).eq.0) then

      angleh = 0.5D0*(angle1+angle2)      ! mid-point

      Nh  = N/2
      Nh1 = Nh+1

      texp  = 1.0D0/(Nh-1.0D0)
      alpha = ratio**texp

      If(abs(alpha-1.0D0).gt.0.0000001) then
       factor = (1.0D0-alpha)/(1.0D0-alpha**Nh)
      Else
       factor = 1.0D0/(Nh-1.0D0+1.0D0)
      End If

      deltat = (angleh-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*cos(te(1))
      Ye(1) = Ycnt + radius*sin(te(1))
      se(1) = sinit

      Do i=2,Nh1                          ! up to mid-point
        te(i) = te(i-1)+deltat
        Xe(i) = Xcnt + radius*cos(te(i))
        Ye(i) = Ycnt + radius*sin(te(i))
        se(i) = se(i-1)+radius*abs(deltat)
        deltat = deltat*alpha
      End Do

      deltat = deltat/alpha

      Do i=Nh1+1,N1                        ! reflect
        te(i) = te(i-1)+deltat
        Xe(i) = Xcnt + radius*cos(te(i))
        Ye(i) = Ycnt + radius*sin(te(i))
        se(i) = se(i-1)+radius*abs(deltat)
        deltat = deltat/alpha
      End Do

      Go to 99

      End If

!---
! odd number of elements
!---

      texp  = 2.0D0/(N-1.0D0)
      alpha = ratio**texp

      if (abs(alpha-1.0D0) > 0.0000001) then
       tmp1   = (N+1.0)/2.0D0
       tmp2   = (N-1.0)/2.0D0
       factor = (1.0D0-alpha)/(2.0D0-alpha**tmp1-alpha**tmp2)
      else
       factor = 1.0/N
      end if

      deltat = (angle2-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*cos(te(1))
      Ye(1) = Ycnt + radius*sin(te(1))
      se(1) = sinit

      do i=2,(N+3)/2                      ! up to mid point + 1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*cos(te(i))
        Ye(i) = Ycnt + radius*sin(te(i))
        se(i) = se(i-1)+radius*abs(deltat)
        deltat = deltat*alpha
      end do

      deltat = deltat/(alpha**2)

      Do i=(N+5)/2,N1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*cos(te(i))
        Ye(i) = Ycnt + radius*sin(te(i))
        se(i) = se(i-1)+radius*abs(deltat)
        deltat = deltat/alpha
      End Do

!-----------------------------------

  99  Continue

!-----------------------
! compute the mid-points
!-----------------------

      Do i=1,N
       tm(i) = 0.5D0*(te(i)+te(i+1))
       sm(i) = 0.5D0*(se(i)+se(i+1))
       Xm(i) = Xcnt + radius*cos(tm(i))
       Ym(i) = Ycnt + radius*sin(tm(i))
      End Do

!-----
! Done
!-----

      Return
      End
