      subroutine lgf_ax_fs &
 
         (Iopt &
         ,x,s &
         ,x0,s0 &
         ,G &
         ,Gx,Gs &
         )

!=========================================
! FDLIB, CFDLAB, BEMLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!=========================================

!-----------------------------------------
! Free-space axisymmetric Green's function.
! of Laplace's equation
!
!  Iopt =  1 compute only the Green's function
!       ne 1 compute Green's function and gradient
!-------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      ! parameters
      integer,                                 intent(in)   :: iopt
      real(FP_KIND),                           intent(in)   :: x, s, x0, s0
      real(FP_KIND),                           intent(out)  :: g, gx, gs

      ! internal variables
      real(FP_KIND) :: dx, dxs, ss0s, rks, rksc, f, e, den, RJ10, RI10, RJ30, RI30, RJ31, RI31, cf

!--------
! prepare
!--------

      Dx = x-x0
      Dxs = Dx*Dx
      ss0s = (s+s0)**2

      rks = FOUR*s*s0/(Dxs+ss0s)

      call ell_int (rks,F,E)

!-----------------
! Green's function
!-----------------

      RJ10 = F
      den = sqrt(Dxs+ss0s)

      RI10 = FOUR*RJ10/den

      G = RI10/pi4

      If(Iopt == 1) Go to 99

!---------------------
! compute: I30 and I31
!---------------------

      rksc = ONE-rks
      RJ30 = E/rksc
      RJ31 = (-TWO*F+(TWO-rks)*E/rksc)/rks
      cf   = FOUR/den**3
      RI30 = cf*RJ30
      RI31 = cf*RJ31

!---------
! gradient
!---------

      Gx = - dx * RI30
      Gs = - s*RI30+s0*RI31

      Gx = Gx/pi4
      Gs = Gs/pi4

!-----
! Done
!-----

  99  Continue

      Return
      End
