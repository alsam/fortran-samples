      subroutine ell_int (RK2,F,E)

!-----------------------------------------
! FDLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!----------------------------------------

!------------------------------------------------
! This program accompanies the book:
!
! C. Pozrikidis
! ``Numerical Computation in Science and Engineering''
!  Oxford University Press, 1998
!------------------------------------------------

!-------------------------------------------
!  Evaluation of complete elliptic integrals
!  of the first and second kind
!  using a recursion formula
!
!  SYMBOLS:
!  --------
!
!  F:      first kind
!  E:      second kind
!  acc: specified accuracy
!-------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      real(FP_KIND), intent(in)  :: RK2
      real(FP_KIND), intent(out) :: F, E
      real(FP_KIND)              :: RK, G, B, C, D

!----------
! launching
!----------

      RK = sqrt(RK2)
      F  = pih
      E  = ONE
      G  = ONE
      B  = RK

      do
        C = sqrt(ONE - B**2)
        B = (ONE - C) / (ONE + C)
        D = F*B
        F = F+D
        G = HALF*G*B
        E = E+G
  
        if (abs(D) <= acc) exit

      end do

      E = F * (ONE - HALF*RK2*E)

!-----
! Done
!-----

      return
      end
