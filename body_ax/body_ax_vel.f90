      subroutine velocity &
 
         (phi &
         ,DphiDn0 &
         ,X00,Y00 &
         ,Ux,Uy &
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

!---------------------------------------------------
! Compute the velocity at the point: X00, Y00
! using the boundary-integral
! representation
!
! Velocity is computed by numerical differentiation
! of the potential
!-----------------------------------------------------

      use precision_mod
      use parameters_mod
      use common_mod01
      implicit none
 
      ! parameters
      real(FP_KIND),                                      intent(in)  :: x00, y00
      real(FP_KIND),                                      intent(out) :: ux,  uy
      real(FP_KIND), dimension(MAX_SEGMENTS, MAX_ELEMS),  intent(in)  :: phi
      real(FP_KIND), dimension(MAX_SEGMENTS * MAX_ELEMS), intent(in)  :: dphidn0

      ! internal variables
      integer       :: iopt, ising, j, k, l
      real(FP_KIND) :: eps2, phi1, phi2, phi3, phi4, psi, qqq, www, rad, t1, t2, ulvr, vlvr, unused
      real(FP_KIND) :: xcnt, ycnt, x01, x02, y01, y02, x1, y1, x2, y2

!--------
! prepare
!--------

      eps2 = 2.0*eps   ! for centered differences

      Ising = 0

!---
! initialize
!---

      phi1 = ZERO 
      phi2 = ZERO 
      phi3 = ZERO 
      phi4 = ZERO 

      unused = ZERO

      x01 = X00-eps
      x02 = X00+eps

      y01 = Y00-eps
      y02 = Y00+eps

!---------------------------------
! Boundary integral representation
!---------------------------------

      j = 0         ! collocation point counter

      do K=1,NSG

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        do L=1,NE(K)

        x1 = XW(K,L)
        y1 = YW(K,L)
        t1 = TW(K,L)

        x2 = XW(K,L+1)
        y2 = YW(K,L+1)
        t2 = TW(K,L+1)

        j = j+1

         call body_ax_sdlp &
 
           (X01,Y00,unused &
           ,X1,Y1,T1 &
           ,X2,Y2,T2 &
           ,NGL &
           ,Ising &
           ,Itp(k) &
           ,rad,xcnt,ycnt &
           ,QQQ &
           ,WWW &
           )

        phi1 = phi1 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_ax_sdlp &
 
            (X02,Y00,unused &
            ,X1,Y1,T1 &
            ,X2,Y2,T2 &
            ,NGL &
            ,Ising &
            ,Itp(k) &
            ,rad,xcnt,ycnt &
            ,QQQ &
            ,WWW &
            )

        phi2 = phi2 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_ax_sdlp &
 
            (X00,Y01,unused &
            ,X1,Y1,T1 &
            ,X2,Y2,T2 &
            ,NGL &
            ,Ising &
            ,Itp(k) &
            ,rad,xcnt,ycnt &
            ,QQQ &
            ,WWW &
            )

        phi3 = phi3 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_ax_sdlp &
 
            (X00,Y02,unused &
            ,X1,Y1,T1 &
            ,X2,Y2,T2 &
            ,NGL &
            ,Ising &
            ,Itp(k) &
            ,rad,xcnt,ycnt &
            ,QQQ &
            ,WWW &
            )

        phi4 = phi4 - dphidn0(j)*QQQ + phi(k,l)*WWW

        end do

      end do

!--------------------------
! numerical differentiation
!--------------------------

      Ux = (phi2-phi1)/eps2
      Uy = (phi4-phi3)/eps2

!--------------------------------
! add the incident streaming flow
!--------------------------------

      Ux = Ux + Vx

!---------------------------------------------
! add the velocity due to the line vortex ring
!---------------------------------------------

      Iopt = 1

      call lvr_fs &
 
             (Iopt &
             ,X00,Y00 &
             ,Xlvr,Ylvr &
             ,ulvr,vlvr &
             ,psi &
             )

      Ux = Ux + cr*ulvr
      Uy = Uy + cr*vlvr

!-----
! Done
!-----

      return
      end
