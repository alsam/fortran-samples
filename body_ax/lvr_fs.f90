      subroutine lvr_fs &
                        (Iopt &
                        ,x,s &
                        ,x0,s0 &
                        ,u,v &
                        ,psi &
                        )

!-----------------------------------------
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!----------------------------------------

!--------------------------------------------------
! Velocity components (u, v)
! and Stokes streamfunction (psi)
! due to a line vortex ring
!
! If Iopt.eq.1 the subroutine computed the velocity
! If Iopt.ne.1 the subroutine computed the velocity
!              and the Stokes stream function
!--------------------------------------------------
 
      use precision_mod
      use parameters_mod
      implicit none

      ! parameters
      integer,                                 intent(in)   :: iopt
      real(FP_KIND),                           intent(in)   :: x, s, x0, s0
      real(FP_KIND),                           intent(out)  :: u, v, psi

      ! internal variables
      real(FP_KIND) :: dx, dxs, rks, rk, f, e, RJ30, RJ31, RI30, RI31, RJ11, RI11, cf

!---
! prepare
!---

      Dx  = x-x0
      Dxs = Dx**2

      rks = FOUR*s*s0/(Dxs+(s+s0)**2)

      call ell_int (rks,F,E)

      RJ30 = E/(ONE-rks)
      RJ31 = (-TWO*F + E*(TWO-rks)/(ONE-rks))/rks

      cf = FOUR/sqrt((Dxs+(s+s0)**2)**3)

      RI30 = cf * RJ30
      RI31 = cf * RJ31

      u = (-s*RI31+s0*RI30)/pi4
      v =  Dx*RI31/pi4

      If(Iopt == 1) Go to 99
!---
! compute the Stokes stream function
!---

      rk   = sqrt(rks)
      RJ11 = ((TWO-rks)*F-TWO*E)/rks

      cf = FOUR/sqrt((Dxs+(s+s0)**2))

      RI11 = cf * RJ11
      psi  = rk*s*s0*RI11/pi4

!-----
! Done
!-----

 99   Continue

      Return
      End
