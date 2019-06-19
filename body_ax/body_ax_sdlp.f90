      subroutine body_ax_sdlp &
 
         (X0,Y0,T0 &
         ,X1,Y1,T1 &
         ,X2,Y2,T2 &
         ,NGL &
         ,Ising &
         ,Itype &
         ,Rad,xcnt,ycnt &
         ,QQQ &
         ,WWW &
         )

!=========================================
! FDLIB, BEMLIB, CFDLAB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!=========================================

!----------------------------------------------------------
! Compute the single-layer and double-layer potential over
! a straight segment or circular arc
!
! LEGEND:
! -------
!
! QQQ: single-layer potential
! WWW: double-layer potential
!----------------------------------------------------------

      use precision_mod
      use parameters_mod
      use gauss_legendre, only : zz, ww

      implicit none

      ! parameters
      integer,       intent(in)  :: ngl ! number of gauss-legendre quadrature points
      integer,       intent(in)  :: ising, itype
      real(FP_KIND), intent(in)  :: x0, y0, t0, x1, y1, t1, x2, y2, t2, rad, xcnt, ycnt
      real(FP_KIND), intent(out) :: qqq, www

      ! internal variables
      integer        :: i, iopt
      real(FP_KIND)  :: cs, sn, dgdx, dgdy, dists, dr, g, ornt, t, td, tm, vnx, vny, x, xm, xd, y, ym, yd

!-----------
! initialize
!-----------

      Iopt = 2   ! for the Green's function

      QQQ = ZERO
      WWW = ZERO

!---------------------------
! prepare for the quadrature
!---------------------------

      if (Itype == 1) then     ! straight segments

        XM = HALF*(X2+X1)
        XD = HALF*(X2-X1)
        YM = HALF*(Y2+Y1)
        YD = HALF*(Y2-Y1)
        DR = sqrt(XD*XD+YD*YD)

        vnx =  YD/DR           ! unit normal vector
        vny = -XD/DR

      else                     ! circular arcs

        TM = HALF*(T2+T1)
        TD = HALF*(T2-T1)
        DR = Rad*abs(TD)
        ornt = ONE             ! orientation index
        if (TD < ZERO) ornt = -ONE

      end if

!---
! loop over Gaussian points
!---

      gauss_pt: do, i=1,NGL

        if (itype == 1) then

          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)

        else

          T  = TM + TD*ZZ(i)
          cs = cos(t)
          sn = sin(t)
          X  = xcnt + Rad*cs
          Y  = ycnt + Rad*sn
          vnx = cs * ornt  ! unit normal vector
          vny = sn * ornt  ! when arc is counter-clockwise,
                           ! normal vector points away from center
        end if

        call lgf_ax_fs &
 
          (Iopt &
          ,X,Y &
          ,X0,Y0 &
          ,G &
          ,dGdx &
          ,dGdy &
          )

!--------------------------------------------------
!  treat the slp singularity
!
!  Subtract off
!  the logarithmic singularity corresponding to the
!  free-space Green's function
!
!  NOTE: The double-layer singularity is not treated
!
!--------------------------------------------------

        if (Ising == 1) then
          if (Itype == 1) Dists = (X-X0)**2+(Y-Y0)**2
          if (Itype == 2) Dists = ( Rad*(T0-T) )**2
          G  = G + log(Dists)/PI4
        end if
  
        QQQ = QQQ + G * Y * WW(i)
  
        WWW = WWW + (dGdx*vnx+dGdy*vny) * Y * WW(i)

      end do gauss_pt

!-------------------------
! finish up the quadrature
!-------------------------

      QQQ = QQQ * DR
      WWW = WWW * DR

!------------------------------------
! add slp singularity back to the slp
!------------------------------------

      if (Ising == 1) then
        QQQ = QQQ - TWO * DR * (LOG(DR) - ONE) / PI2
      end if

!-----
! Done
!-----

      return
      end
